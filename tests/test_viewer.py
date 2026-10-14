"""Interactive views (crevice.viewer): mesh extraction, scene description and HTML.

Pure-Python checks: no browser is started. Surfaces must be the display
surfaces the native PyMOL/VMD/ChimeraX scenes draw from the same fields.
"""
from __future__ import annotations

import base64
import gzip
import json
import re
import sys

import numpy as np
import pytest

pytest.importorskip("skimage")

from crevice import pore_profile, void_cast
from crevice.cast_segments import SEGMENT_HEX, lining_residues, segment_casts
from crevice.viewer import (CreviceView, Layer, cast_mesh, field_mesh, read_dx_field, view)

from test_cast_segments import caged_channel
from test_capped_channels import OPTIONS, _publish_capped


def _decode(spec):
    raw = base64.b64decode(spec["data"])
    if spec["encoding"] == "gzip+base64":
        raw = gzip.decompress(raw)
    return np.frombuffer(raw, dtype=np.uint32 if spec["dtype"] == "uint32" else np.float32)


@pytest.fixture(scope="module")
def caged():
    frame = caged_channel()
    profile = pore_profile(frame, **OPTIONS)
    cast = void_cast(frame, profile=profile, spacing=0.5, min_radius=0.0, max_grid_points=400000)
    return frame, profile, cast


@pytest.fixture(scope="module")
def published(tmp_path_factory):
    return _publish_capped(tmp_path_factory.mktemp("view"))


def test_sphere_mesh_lies_on_the_isosurface_with_outward_normals():
    axis = np.arange(-6, 7) * 0.5
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    field = 2.0 - np.sqrt(x**2 + y**2 + z**2)  # 0.5 at radius 1.5
    rotation = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))[0]
    if np.linalg.det(rotation) < 0:
        rotation[0] *= -1
    origin = np.array([10.0, -4.0, 2.0])
    deltas = 0.5 * rotation
    mesh = field_mesh(field, origin - np.full(3, 6) @ deltas, deltas)
    radii = np.linalg.norm(mesh["vertices"] - origin, axis=1)
    assert mesh["vertices"].dtype == np.float32 and mesh["faces"].dtype == np.uint32
    assert np.allclose(radii, 1.5, atol=0.08)
    outward = np.einsum("ij,ij->i", mesh["normals"], mesh["vertices"] - origin)
    assert (outward > 0).all() and np.allclose(np.linalg.norm(mesh["normals"], axis=1), 1, atol=1e-5)
    # Faces wind counter-clockwise seen from outside (consistent with the normals).
    v = mesh["vertices"][mesh["faces"]]
    face_normal = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    proper = np.linalg.norm(face_normal, axis=1) > 1e-6  # marching cubes can emit zero-area faces
    assert (np.einsum("ij,ij->i", face_normal, v.mean(axis=1) - origin)[proper] > 0).all()
    assert field_mesh(np.zeros((3, 3, 3)), (0, 0, 0), np.eye(3)) is None


def test_dx_meshes_reproduce_the_native_display_surfaces(caged, tmp_path):
    from scipy.spatial import cKDTree
    from crevice.presentation import write_viewer_structure
    from crevice.volume_export import write_volume_viewer_bundle
    frame, profile, cast = caged
    segments = segment_casts(frame, profile, cast)
    viewer = write_viewer_structure(frame, tmp_path / "v.pdb")
    files = write_volume_viewer_bundle(cast, structure_path=viewer, output_dir=tmp_path, prefix="t", frame=frame,
                                       profile=profile, segments=segments)
    expected = np.load(tmp_path / "t_geometry_reference.npz")["surface_vertices"]
    tree = cKDTree(expected)
    for name in ["lumen"] + [s.metadata["segment"] for s in segments]:
        field, origin, deltas = read_dx_field(files[f"segment_{name}_display_dx"])
        mesh = field_mesh(field, origin, deltas)
        distance, _ = tree.query(mesh["vertices"])
        assert distance.max() < 2e-3, name  # DX values carry 6 significant digits
    # The in-memory path computes the same segment surface as the file path.
    from dataclasses import replace
    part = replace(cast, components=(segments[0],), points=tuple(segments[0].points))
    distance, _ = tree.query(cast_mesh(part)["vertices"])
    assert distance.max() < 2e-3


def test_publish_directory_scene_has_toggleable_segments_lining_and_centre_line(published):
    v = view(published)
    names = v.layer_names
    assert names[0] == "protein" and names[1:3] == ["lumen", "exit_1"]
    assert {"lining_lumen", "lining_exit_1", "centre_line"} <= set(names)
    assert v.layer("lumen").color == SEGMENT_HEX["lumen"] and v.layer("exit_1").color == SEGMENT_HEX["exit"]
    assert not v.layer("lining_lumen").visible and not v.layer("centre_line").visible
    assert v.layer("lining_exit_1").color == SEGMENT_HEX["exit"]
    v.toggle("lining")
    assert v.layer("lining_lumen").visible and v.layer("lining_exit_1").visible
    v.toggle("exit_1", False).toggle("protein")
    assert not v.layer("exit_1").visible and not v.layer("protein").visible
    with pytest.raises(KeyError, match="layers are"):
        v.toggle("nonexistent")
    scene = v.scene()
    assert scene["schema"] == "crevice.interactive_view.v1" and scene["structure"]["format"] == "pdb"
    lumen = next(layer for layer in scene["layers"] if layer["name"] == "lumen")
    vertices, faces = _decode(lumen["vertices"]), _decode(lumen["faces"])
    assert vertices.size == 3 * lumen["vertex_count"] and faces.max() < lumen["vertex_count"]
    np.testing.assert_array_equal(vertices.reshape(-1, 3), v.layer("lumen").data["vertices"])
    sticks = next(layer for layer in scene["layers"] if layer["name"] == "lining_lumen")
    assert sticks["residues"] and all(r["chain"] == "A" for r in sticks["residues"])
    q = np.asarray(scene["camera"]["quaternion"])
    assert np.isclose(np.linalg.norm(q), 1)


def test_quaternion_matches_the_native_camera_rotation():
    from crevice.viewer import _rotation_quaternion
    rotation = np.linalg.qr(np.random.default_rng(3).normal(size=(3, 3)))[0]
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    x, y, z, w = _rotation_quaternion(rotation)
    rebuilt = np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
                        [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
                        [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])
    np.testing.assert_allclose(rebuilt, rotation, atol=1e-9)


def test_standalone_html_is_self_contained(published, tmp_path):
    v = view(published)
    path = v.save_html(tmp_path / "view.html")
    text = path.read_text()
    assert text.startswith("<!DOCTYPE html>") and "3Dmol.js 2.5.5" in text and "BSD-3-Clause" in text
    assert "Copyright (c) 2014, University of Pittsburgh" in text
    assert not re.search(r"<script[^>]+src=", text) and "<link" not in text
    scene = json.loads(re.search(r"var creviceScene = (\{.*?\});\n", text, re.S).group(1).replace("<\\/", "</"))
    assert [layer["name"] for layer in scene["layers"]] == v.layer_names
    cdn = v.to_html(inline_js=False)
    assert 'src="https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"' in cdn and len(cdn) < len(text)


def test_cli_view_writes_the_html_with_requested_layers(published, tmp_path, capsys):
    from crevice.cli import main
    out = tmp_path / "v.html"
    assert main(["view", str(published), "-o", str(out), "--show", "lining_lumen", "--hide", "protein"]) == 0
    assert f"view file={out}" in capsys.readouterr().out
    scene = json.loads(re.search(r"var creviceScene = (\{.*?\});\n", out.read_text(), re.S).group(1)
                       .replace("<\\/", "</"))
    state = {layer["name"]: layer["visible"] for layer in scene["layers"]}
    assert state["lining_lumen"] and not state["protein"] and state["lumen"]
    assert main(["view", str(published)]) == 0
    assert (published / "capped_view.html").is_file()
    assert main(["view", str(tmp_path / "missing")]) == 2


def test_region_scene_shows_surfaces_member_waters_and_snapshots(tmp_path):
    from crevice.water_density import write_scalar_dx
    axis = np.arange(10) - 4.5
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    write_scalar_dx((x**2 + y**2 + z**2 < 9).astype(float), np.zeros(3), 1.0, tmp_path / "r_display.dx")
    pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n"
    for i in (0, 5):
        (tmp_path / f"r_analysis_frame_{i:06d}.pdb").write_text(pdb)
    (tmp_path / "r_w_000005.pdb").write_text(
        "HETATM    1  O   HOH W   1       4.000   4.000   4.000  1.00  0.00           O\nEND\n")
    scene = {"prefix": "r", "protein_pdb": "r_analysis_frame_000000.pdb",
             "surfaces": [{"key": "cavity", "dx": "r_display.dx", "level": .5, "color_rgb": [.69, .73, .73],
                           "opacity": .18, "meaning": "reference cavity display"}],
             "residues": [{"residue": "A:ALA1", "viewer_chain": "A", "viewer_resid": 1}],
             "display_residues": ["A:ALA1"],
             "snapshots": [{"pdb": "r_analysis_frame_000000.pdb", "source_frame_index": 0},
                           {"pdb": "r_analysis_frame_000005.pdb", "source_frame_index": 5}],
             "member_waters": {"color_rgb": [1., .5, 0.], "sphere_radius_A": .7, "definition": "geometric",
                               "snapshots": [{"source_frame_index": 0, "pdb": None},
                                             {"source_frame_index": 5, "pdb": "r_w_000005.pdb"}]}}
    (tmp_path / "r_analysis_scene.json").write_text(json.dumps(scene))
    first = view(tmp_path)
    assert first.layer_names == ["protein", "region", "lining_region"]
    last = view(tmp_path, snapshot=-1)
    assert last.layer_names == ["protein", "region", "lining_region", "waters"]
    assert last.layer("waters").color == "#ff8000" and last.layer("waters").data["points"] == [[4.0, 4.0, 4.0]]
    assert "frame 5" in last.title
    with pytest.raises(ValueError, match="snapshot"):
        view(tmp_path, snapshot=7)
    with pytest.raises(FileNotFoundError, match="no CREVICE"):
        view(tmp_path, kind="cast")


def test_in_memory_results(caged):
    frame, profile, cast = caged
    segments = segment_casts(frame, profile, cast)
    lining = lining_residues(frame, cast, segments)
    v = view(frame=frame, cast=cast, segments=segments, profile=profile, lining=lining, waters=[[0., 0., 0.]])
    assert v.layer_names[:3] == ["protein", "lumen", "exit_1"] or v.layer_names[:3] == ["protein", "entry_1", "lumen"]
    assert "centre_line" in v.layer_names and v.layer("waters").visible
    plain = view(frame=frame, cast=cast)
    assert plain.layer("cast").color == SEGMENT_HEX["lumen"]  # a channel cast stays teal
    with pytest.raises(ValueError, match="in-memory"):
        view()


def test_jupyter_display_needs_py3dmol_with_a_clear_message(published, monkeypatch):
    v = CreviceView([Layer("a", "A", "spheres", "#ff8000", data={"points": [[0, 0, 0]]})])
    monkeypatch.setitem(sys.modules, "py3Dmol", None)
    with pytest.raises(ImportError, match=r"pip install 'crevice\[view\]'"):
        v.jupyter_html()
    assert "CreviceView(a)" == repr(v)


def test_jupyter_fragment_drives_the_py3dmol_viewer(published):
    pytest.importorskip("py3Dmol")
    fragment = view(published).jupyter_html(width=400, height=300)
    uid = re.search(r'id="3dmolviewer_(\d+)"', fragment).group(1)
    assert f"CreviceView.build(viewer_{uid}," in fragment and f"crevice_controls_{uid}" in fragment

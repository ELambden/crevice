"""ENTRY / LUMEN / EXIT cast segments and lining residues (synthetic structures only)."""
from __future__ import annotations

import csv
import json
import math

import numpy as np
import pytest

from crevice import Atom, StructureFrame, pore_profile, void_cast
from crevice.cast_segments import (SEGMENT_HEX, lining_residues, segment_casts, segment_roles, segment_summary)
from crevice.channel_exits import BulkEnvelope

from test_capped_channels import OPTIONS, _publish_capped, capped_channel


def caged_channel():
    """Carbon lattice block (2 A, x, y, z in [-8, 8]) pierced by a straight lumen (r < 3 A).

    Four 4-atom pillars at radius 5.5 A stand on the top face (z = 10-16 A). The
    4.4 A gaps between them let the enclosure probe escape sideways (the
    upper mouth is open) but keep the 6 A bulk probe out of the cage, so the
    space inside the cage is a vestibule beyond the upper mouth. The bottom
    face opens straight into bulk (no vestibule).
    """
    coords = [(2.0*i, 2.0*j, 2.0*k) for i in range(-4, 5) for j in range(-4, 5) for k in range(-4, 5)
              if math.hypot(2.0*i, 2.0*j) >= 3.0]
    for angle in (0, 90, 180, 270):
        x, y = 5.5*math.cos(math.radians(angle)), 5.5*math.sin(math.radians(angle))
        coords += [(round(x, 3), round(y, 3), z) for z in (10.0, 12.0, 14.0, 16.0)]
    coords.sort(key=lambda p: (p != (8.0, 0.0, -8.0), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


@pytest.fixture(scope="module")
def caged():
    frame = caged_channel()
    profile = pore_profile(frame, **OPTIONS)
    cast = void_cast(frame, profile=profile, spacing=0.5, min_radius=0.0, max_grid_points=400000)
    return frame, profile, cast


def _keys(points):
    return {tuple(round(v, 6) for v in p.position) for p in points}


def test_vestibule_segment_is_disjoint_beyond_the_mouth_and_outside_bulk(caged):
    frame, profile, cast = caged
    before = _keys(p for c in cast.components for p in c.points)
    segments = segment_casts(frame, profile, cast)
    assert [(c.metadata["segment"], c.kind, c.metadata["end"]) for c in segments] == [
        ("exit_1", "axial_vestibule", "upper")]
    exit_1 = segments[0]
    assert exit_1.volume == pytest.approx(exit_1.point_count*0.5**3) and exit_1.volume > 50
    # The lumen is untouched and the segment never overlaps it.
    assert _keys(p for c in cast.components for p in c.points) == before
    assert before.isdisjoint(_keys(exit_1.points))
    assert all(p.t > cast.metadata["mouth_planes_A"][1] for p in exit_1.points)
    # The vestibule stops at the bulk boundary (outer boundary, not a width limit).
    envelope = BulkEnvelope(frame.atoms, profile.axis_origin, cast.metadata["grid_basis"], profile.points[-1].t, 1,
                            spacing=1.0, bulk_radius=6.0)
    assert all(envelope.depth(p.position) > 0 for p in exit_1.points)
    # Space-filling: the vestibule is wider than the lumen it continues.
    assert max(p.raw_clearance for p in exit_1.points) > max(p.raw_clearance for c in cast.components
                                                             for p in c.points)
    rows = segment_summary(cast, segments)
    assert [r["segment"] for r in rows] == ["lumen", "exit_1"]
    assert rows[0]["volume_A3"] == cast.total_volume
    assert rows[1]["volume_A3"] == exit_1.volume


def test_labels_follow_the_axis_and_can_be_overridden(caged):
    frame, profile, cast = caged
    assert segment_roles(profile) == ({"lower": "entry", "upper": "exit"}, "axis_start_entry")
    flipped = segment_casts(frame, profile, cast, entry_end="end")
    assert [(c.metadata["segment"], c.metadata["end"], c.metadata["label_rule"]) for c in flipped] == [
        ("entry_1", "upper", "override_end")]
    with pytest.raises(ValueError, match="entry_end"):
        segment_roles(profile, "middle")


def test_capped_end_legs_are_exits_and_the_open_end_is_the_entry():
    frame = capped_channel(second_portal=True)
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    cast = void_cast(frame, profile=profile, spacing=0.5, min_radius=0.0, max_grid_points=400000)
    roles, rule = segment_roles(profile)
    assert roles == {"lower": "entry", "upper": "exit"} and rule == "capped_end_exits"
    segments = segment_casts(frame, profile, cast)
    assert [(c.metadata["segment"], c.kind) for c in segments] == [("exit_1", "lateral_exit"),
                                                                   ("exit_2", "lateral_exit")]
    assert _keys(segments[0].points).isdisjoint(_keys(segments[1].points))
    override = segment_casts(frame, profile, cast, entry_end="end")
    assert [c.metadata["segment"] for c in override] == ["entry_1", "entry_2"]
    # Leg volumes are those of the lateral exit casts.
    from crevice.channel_exits import lateral_exit_casts
    assert [c.volume for c in segments] == [c.volume for c in lateral_exit_casts(frame, profile, cast)]


def test_lining_residues_match_a_brute_force_distance_search(caged):
    frame, profile, cast = caged
    segments = segment_casts(frame, profile, cast)
    rows = lining_residues(frame, cast, segments, cutoff=3.3)
    groups = {"lumen": [p.position for c in cast.components for p in c.points],
              "exit_1": [p.position for p in segments[0].points]}
    for name, points in groups.items():
        points = np.asarray(points)
        expected = {}
        for atom in frame.atoms:
            d = float(np.min(np.linalg.norm(points-np.asarray(atom.coord), axis=1)))
            if d <= 3.3:
                expected[atom.residue_key.label] = d
        found = {r["residue"]: r["min_distance_A"] for r in rows if r["segment"] == name}
        assert set(found) == set(expected) and found
        assert all(found[k] == pytest.approx(expected[k], abs=1e-4) for k in found)
        assert all(r["role"] == f"{name}_lining" for r in rows if r["segment"] == name)
    # The upper pillar atoms line the vestibule, not the lumen (which reaches z = 10.3 A).
    pillars = {a.residue_key.label for a in frame.atoms if a.z >= 14}
    assert pillars & {r["residue"] for r in rows if r["segment"] == "exit_1"}
    assert not pillars & {r["residue"] for r in rows if r["segment"] == "lumen"}
    with pytest.raises(ValueError, match="cutoff"):
        lining_residues(frame, cast, segments, cutoff=0)


def test_bundle_draws_each_segment_and_lining_group_as_its_own_object(caged, tmp_path):
    from crevice.presentation import write_viewer_structure
    from crevice.volume_export import write_volume_viewer_bundle
    frame, profile, cast = caged
    segments = segment_casts(frame, profile, cast)
    lining = lining_residues(frame, cast, segments)
    viewer = write_viewer_structure(frame, tmp_path/"v.pdb")
    files = write_volume_viewer_bundle(cast, structure_path=viewer, output_dir=tmp_path, prefix="t", frame=frame,
                                       profile=profile, segments=segments, lining=lining)
    pml = (tmp_path/"t_volume.pml").read_text()
    assert "crevice_surface('crevice_lumen', 'lumen')" in pml and "crevice_surface('crevice_exit_1', 'exit_1')" in pml
    assert 'crevice_surface("crevice_volume"' not in pml
    assert "cmd.disable('crevice_lumen_lining')" in pml and "cmd.disable('crevice_exit_1_lining')" in pml
    assert 'cmd.extend("crevice_lining", crevice_lining)' in pml and "cmd.label" not in pml
    tcl = (tmp_path/"t_volume.tcl").read_text()
    assert "mol rename $crevice_volume crevice_lumen" in tcl and "mol rename $crevice_segment crevice_exit_1" in tcl
    assert "mol showrep $crevice_protein $crevice_lining_rep(lumen) off" in tcl
    cxc = (tmp_path/"t_volume.cxc").read_text()
    assert "rename #2 crevice_lumen" in cxc and "crevice_exit_1" in cxc and f"{SEGMENT_HEX['exit']}" in cxc
    assert "name frozen crevice_lumen_lining #1/A:" in cxc
    scene = json.loads((tmp_path/"t_scene.json").read_text())
    assert [s["pymol_object"] for s in scene["cast_segments"]] == ["crevice_lumen", "crevice_exit_1"]
    metadata = json.loads((tmp_path/"t_volume_metadata.json").read_text())
    # The display objects partition the displayed cast exactly once.
    assert sum(s["display_points"] for s in metadata["cast_segments"]) == metadata["display_points"]
    assert metadata["measured_volume_A3"] == cast.total_volume
    assert metadata["cast_segments"][1]["measured_volume_A3"] == segments[0].volume
    assert "segment_exit_1_dx" in files and (tmp_path/"t_exit_1_cast.dx").exists()
    with pytest.raises(ValueError, match="segments"):
        write_volume_viewer_bundle(cast, structure_path=viewer, output_dir=tmp_path/"r", prefix="t", frame=frame,
                                   segments=segments, regions={"a": cast.points})


def test_exit_casts_alone_still_get_their_own_objects(tmp_path):
    from crevice.channel_exits import lateral_exit_casts
    from crevice.presentation import write_viewer_structure
    from crevice.volume_export import write_volume_viewer_bundle
    frame = capped_channel()
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    cast = void_cast(frame, profile=profile, spacing=0.5, min_radius=0.0, max_grid_points=400000)
    viewer = write_viewer_structure(frame, tmp_path/"v.pdb")
    write_volume_viewer_bundle(cast, structure_path=viewer, output_dir=tmp_path, prefix="t", frame=frame,
                               profile=profile, exit_casts=lateral_exit_casts(frame, profile, cast))
    pml = (tmp_path/"t_volume.pml").read_text()
    assert "crevice_surface('crevice_exit_1', 'exit_1')" in pml and "crevice_exit_1_color" in pml


def test_publish_reports_segments_and_lining(tmp_path):
    out = _publish_capped(tmp_path, "--lining-cutoff", "3.0")
    cast = json.loads((out/"capped_void_cast.json").read_text())
    meta = cast["metadata"]
    names = [r["segment"] for r in meta["cast_segments"]]
    assert names == ["lumen", "exit_1"] and meta["segment_label_rule"]["rule"] == "capped_end_exits"
    assert meta["cast_segments"][0]["volume_A3"] == cast["total_volume"]
    assert meta["segments_total_volume_A3"] == pytest.approx(sum(r["volume_A3"] for r in meta["cast_segments"]))
    assert meta["lining_cutoff_A"] == 3.0
    manifest = json.loads((out/"capped_manifest.json").read_text())
    assert set(manifest["profile_status"]["cast_segment_volumes_A3"]) == {"lumen", "exit_1"}
    rows = list(csv.DictReader((out/"capped_cast_segments.csv").open()))
    assert [r["segment"] for r in rows] == ["lumen", "exit_1"]
    lining = list(csv.DictReader((out/"capped_lining_residues.csv").open()))
    assert {r["segment"] for r in lining} == {"lumen", "exit_1"}
    assert all(float(r["min_distance_A"]) <= 3.0 for r in lining)
    contacts = list(csv.DictReader((out/"capped_residue_contacts.csv").open()))
    lumen = {r["residue"] for r in lining if r["segment"] == "lumen"}
    assert all((r["lumen_lining"] == "True") == (r["residue"] in lumen) for r in contacts)
    (tmp_path/"flip").mkdir()
    flipped = _publish_capped(tmp_path/"flip", "--entry-end", "end")
    meta = json.loads((flipped/"capped_void_cast.json").read_text())["metadata"]
    assert [r["segment"] for r in meta["cast_segments"]] == ["entry_1", "lumen"]

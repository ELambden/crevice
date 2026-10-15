"""Curved and branched pathways (crevice.paths) on synthetic carved blocks.

Software tests on analytic geometry only: a curved tube of known centre line
and radius, a Y-branched tube with a known junction, and a closed cavity with
no route to bulk. No biological result is implied.
"""
import csv
import json
import math

import numpy as np
import pytest

from crevice.models import Atom, StructureFrame
from crevice.paths import find_paths, write_paths_outputs

STEP = 1.5
CARVE = 4.0


def _block(xr, yr, zr, carve):
    """Carbon atoms on a 1.5 Å lattice filling a box, minus those whose centre ``carve`` returns True for."""
    xs, ys, zs = (np.arange(a, b + 1e-9, STEP) for a, b in (xr, yr, zr))
    grid = np.array(np.meshgrid(xs, ys, zs, indexing="ij")).reshape(3, -1).T
    keep = grid[~carve(grid)]
    atoms = tuple(Atom(i + 1, "C", "ALA", "A", i // 10 + 1, *map(float, p), "C") for i, p in enumerate(keep))
    return StructureFrame(atoms), keep


def _distance_to_polyline(points, line):
    line = np.asarray(line, float)
    best = np.full(len(points), np.inf)
    for a, b in zip(line, line[1:]):
        ab = b - a
        t = np.clip(((points - a) @ ab) / (ab @ ab), 0, 1)
        best = np.minimum(best, np.linalg.norm(points - (a + t[:, None] * ab), axis=1))
    return best


# Curved tube: circular arc of radius 18 Å centred at (0, 0, -14) in the xz plane,
# from the bottom face (z = -10) on one side, over the top (0, 0, 4), to the bottom face on the other.
ARC_R, ARC_C = 18.0, np.array([0.0, 0.0, -14.0])
ARC_THETA = math.acos(4.0 / ARC_R)          # half-angle from the top to z = -10
ARC = np.array([[ARC_R * math.sin(t), 0.0, -14.0 + ARC_R * math.cos(t)]
                for t in np.linspace(-ARC_THETA - 0.3, ARC_THETA + 0.3, 400)])


@pytest.fixture(scope="module")
def curved():
    frame, atoms = _block((-21, 21), (-9, 9), (-10, 10), lambda p: _distance_to_polyline(p, ARC) < CARVE)
    return frame, atoms, find_paths(frame, seed=(0.0, 0.0, 4.0), axis="x")


def test_curved_tube_length_radius_and_centre_line(curved):
    frame, atoms, result = curved
    assert result.status == "resolved"
    assert len(result.branches) == 2
    assert sorted(b.side for b in result.branches) == ["+", "-"]
    expected_half = ARC_R * ARC_THETA            # arc length from the top to the bottom face
    start = np.asarray(result.start["node_A"])
    offset = ARC_R * math.atan2(start[0], start[2] + 14.0)   # signed arc offset of the snapped start
    for branch in result.branches:
        pts = np.asarray(branch.points)
        # Centre line follows the analytic arc.
        radial = np.abs(np.linalg.norm(pts[:, [0, 2]] - ARC_C[[0, 2]], axis=1) - ARC_R)
        inside = pts[:, 2] > -8.0
        assert radial[inside].max() < 0.75
        # Length: arc from the start to the bottom face, within the mouth tolerance.
        sign = 1 if branch.side == "+" else -1
        assert abs(branch.length - (expected_half - sign * offset)) < 3.0
        # Radius: exact free-sphere radius on the analytic centre line.
        on_arc = ARC_C + ARC_R * np.column_stack([np.sin(np.linspace(0, sign * ARC_THETA, 200)), np.zeros(200),
                                                  np.cos(np.linspace(0, sign * ARC_THETA, 200))])
        true = np.linalg.norm(on_arc[:, None, :] - atoms[None, :, :], axis=2).min(axis=1) - 1.7
        interior = np.asarray(branch.radii)[inside]
        assert abs(np.median(interior) - np.median(true[on_arc[:, 2] > -8])) < 0.3
        assert branch.bottleneck_radius == pytest.approx(min(branch.radii))
        assert branch.bottleneck_radius > 0.8
        assert all(m == "relaxed" for m in branch.metadata["centre_line"].values())
    # The two routes separate at the start.
    assert len(result.junctions) == 1 and result.junctions[0]["is_start"]


Y_JUNCTION = np.array([0.0, 0.0, 0.0])
Y_LINES = [np.array([[0, 0, -8.0], [0, 0, 0.0]]), np.array([[0, 0, 0.0], [12, 0, 12.0]]),
           np.array([[0, 0, 0.0], [-12, 0, 12.0]])]


@pytest.fixture(scope="module")
def ybranch():
    carve = lambda p: np.min([_distance_to_polyline(p, line) for line in Y_LINES], axis=0) < CARVE
    frame, atoms = _block((-18, 18), (-9, 9), (-16, 10), carve)
    return frame, find_paths(frame, seed=(0.0, 0.0, -6.0))


def test_y_branch_two_branches_and_junction(ybranch):
    frame, result = ybranch
    assert result.status == "resolved"
    assert len(result.branches) == 2
    exits = sorted(b.exit_point[0] for b in result.branches)
    assert exits[0] < -6 and exits[1] > 6
    inner = [j for j in result.junctions if not j["is_start"]]
    assert len(inner) == 1
    assert np.linalg.norm(np.asarray(inner[0]["position_A"]) - Y_JUNCTION) < 3.0
    assert sorted(inner[0]["branches"]) == [1, 2]
    # Both branches share the stem segment and pass the junction at the same arc length.
    first, second = result.branches
    assert first.segments[0] == second.segments[0]
    assert first.junctions[0]["arc_length_A"] == pytest.approx(second.junctions[0]["arc_length_A"], abs=1e-6)
    # Each branch has a lining residue list.
    assert first.lining and second.lining


def test_closed_cavity_is_unresolved():
    centre = np.zeros(3)
    frame, _ = _block((-12, 12), (-12, 12), (-12, 12), lambda p: np.linalg.norm(p - centre, axis=1) < 5.0)
    result = find_paths(frame, seed=(0.0, 0.0, 0.0))
    assert result.status == "unresolved_no_bulk_path"
    assert result.branches == []
    assert result.unresolved and "reaches bulk" in result.unresolved[0]["reason"]


def test_seed_blocked_and_start_options(ybranch):
    frame, _ = ybranch
    blocked = find_paths(frame, seed=(15.0, 0.0, -14.0), seed_tolerance=1.0)
    assert blocked.status == "unresolved_seed_blocked"
    with pytest.raises(ValueError):
        find_paths(frame, seed=(0, 0, 0), seed_residues="1")


def test_outputs_written(tmp_path, ybranch):
    frame, result = ybranch
    files = write_paths_outputs(result, tmp_path, prefix="y", frame=frame)
    names = {p.name for p in files}
    assert {"y_paths.json", "y_paths.csv", "y_path_points.csv", "y_path_lining.csv", "y_path_junctions.csv",
            "y_paths.png", "y_paths.pml", "y_paths.cxc", "y_paths_vmd.tcl"} <= names
    data = json.loads((tmp_path / "y_paths.json").read_text())
    assert data["status"] == "resolved" and len(data["branches"]) == 2
    rows = list(csv.DictReader((tmp_path / "y_path_points.csv").open()))
    assert {"branch", "arc_length_A", "x_A", "y_A", "z_A", "radius_A"} <= set(rows[0])
    pml = (tmp_path / "y_paths.pml").read_text()
    assert "path_1" in pml and "path_2" in pml and "label" not in pml.lower().replace("labels", "")


def test_unresolved_outputs(tmp_path):
    frame, _ = _block((-12, 12), (-12, 12), (-12, 12), lambda p: np.linalg.norm(p, axis=1) < 5.0)
    result = find_paths(frame, seed=(0.0, 0.0, 0.0))
    files = write_paths_outputs(result, tmp_path, prefix="c", frame=frame)
    rows = list(csv.DictReader((tmp_path / "c_paths.csv").open()))
    assert rows[0]["status"] == "unresolved_no_bulk_path"
    assert (tmp_path / "c_paths.png").exists()


def test_view_hook_reads_paths_directory(tmp_path, ybranch):
    from crevice.viewer import view
    frame, result = ybranch
    write_paths_outputs(result, tmp_path, prefix="y", frame=frame, figure=False)
    v = view(tmp_path)
    assert v.layer_names[:1] == ["protein"] and {"path_1", "path_2"} <= set(v.layer_names)
    scene = v.scene()
    assert all(layer["kind"] == "line" for layer in scene["layers"] if layer["name"].startswith("path_"))
    assert view(tmp_path / "y_paths.json").layer_names == v.layer_names

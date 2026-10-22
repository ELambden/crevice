"""MDAnalysis Universe / AtomGroup input gives the same results as the path route.

Software tests on a synthetic channel written to PDB + XTC (and, when the
local GLUT1 files are present, a few frames of that real trajectory). They
check equivalence of the two input routes, not any scientific result.
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")

from crevice import analyze_trajectory, inspect_trajectory, load_structure, load_trajectory
from crevice.parser import load_structure_report
from crevice.trajectory import load_trajectory_report

ROOT = Path(__file__).resolve().parents[1]


def _channel_files(tmp_path, frames=6):
    """A ring-wall channel along z (one ALA per ring atom), a water inside, PDB + XTC."""
    lines, serial = [], 0
    for ring, z in enumerate(np.arange(-9.0, 9.1, 1.5)):
        for k in range(10):
            serial += 1
            angle = 2 * math.pi * k / 10
            x, y = 4.2 * math.cos(angle), 4.2 * math.sin(angle)
            resid = ring * 10 + k + 1
            lines.append(f"ATOM  {serial:5d}  CA  ALA A{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    serial += 1
    lines.append(f"HETATM{serial:5d}  OH2 TIP3W 900    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           O")
    pdb = tmp_path / "channel.pdb"
    pdb.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8", newline="\n")
    universe = mda.Universe(str(pdb))
    base = universe.atoms.positions.copy()
    rng = np.random.default_rng(7)
    xtc = tmp_path / "channel.xtc"
    with mda.Writer(str(xtc), len(universe.atoms), dt=10) as writer:
        for i in range(frames):
            universe.atoms.positions = base + rng.normal(scale=0.15, size=base.shape)
            universe.dimensions = [60, 60, 60, 90, 90, 90]
            universe.trajectory.ts.frame = i
            universe.trajectory.ts.time = 10.0 * i
            writer.write(universe.atoms)
    return pdb, xtc


def test_load_structure_reads_the_current_frame_of_a_universe(tmp_path):
    pdb, xtc = _channel_files(tmp_path)
    universe = mda.Universe(str(pdb), str(xtc))
    universe.trajectory[3]
    frame, report = load_structure_report(universe)
    by_path = load_trajectory(str(pdb), str(xtc), start=3, stop=4).frames[0]
    assert frame == by_path
    assert report["input_format"] == "mdanalysis" and report["input_object"] == "MDAnalysis.Universe"
    assert universe.trajectory.ts.frame == 3
    group = universe.select_atoms("resid 1:20")
    assert load_structure(group) == load_trajectory(str(pdb), str(xtc), selection="resid 1:20", start=3, stop=4).frames[0]
    with pytest.raises(ValueError, match="model_index"):
        load_structure(universe, model_index=2)


def test_gro_path_and_universe_give_identical_frames():
    gro = ROOT / "prod.gro"
    if not gro.is_file():
        pytest.skip("local GLUT1 prod.gro is not present")
    assert load_structure(str(gro)) == load_structure(mda.Universe(str(gro)))


def test_trajectory_loader_iterates_a_universe_and_restores_its_frame(tmp_path):
    pdb, xtc = _channel_files(tmp_path)
    universe = mda.Universe(str(pdb), str(xtc))
    universe.trajectory[2]
    by_object, report = load_trajectory_report(universe, start=1, stride=2)
    by_path, path_report = load_trajectory_report(str(pdb), str(xtc), start=1, stride=2)
    assert list(by_object.frames) == list(by_path.frames)
    assert by_object.times() == by_path.times() and by_object.source == by_path.source
    assert {k: v for k, v in report.items() if k != "input_object"} == path_report
    assert universe.trajectory.ts.frame == 2
    group = universe.select_atoms("resid 5:30")
    grouped, group_report = load_trajectory_report(group)
    assert list(grouped.frames) == list(load_trajectory(str(pdb), str(xtc), selection="resid 5:30").frames)
    assert group_report["selection"] == "AtomGroup"
    assert inspect_trajectory(group)["selected_atom_count"] == len(group)
    with pytest.raises(ValueError, match="not both"):
        load_trajectory(universe, str(xtc))


def test_analyze_trajectory_results_match_between_routes(tmp_path):
    pdb, xtc = _channel_files(tmp_path, frames=4)
    options = dict(analyses=("profile", "features"), profile_kwargs={"axis": "z", "samples": 31})
    by_path = analyze_trajectory(load_trajectory(str(pdb), str(xtc)), **options)
    by_object = analyze_trajectory(load_trajectory(mda.Universe(str(pdb), str(xtc))), **options)
    assert json.dumps(by_path.to_dict(), sort_keys=True, default=str) == \
        json.dumps(by_object.to_dict(), sort_keys=True, default=str)
    assert any(frame.profile is not None for frame in by_path.frames)


def test_hydration_trajectory_accepts_a_universe(tmp_path):
    from crevice.hydration_trajectory import analyze_hydration_trajectory
    pdb, xtc = _channel_files(tmp_path, frames=3)
    options = dict(residue_ids=["A:ALA1", "A:ALA2"], sasa_samples=32, max_frames=3, selection="protein")
    by_path = analyze_hydration_trajectory(str(pdb), str(xtc), **options)
    universe = mda.Universe(str(pdb), str(xtc))
    by_object = analyze_hydration_trajectory(universe, **options)
    for key, value in by_path["arrays"].items():
        assert np.array_equal(np.asarray(value), np.asarray(by_object["arrays"][key]), equal_nan=np.asarray(value).dtype.kind == "f"), key
    by_group = analyze_hydration_trajectory(universe.select_atoms("protein"), **{k: v for k, v in options.items() if k != "selection"})
    assert np.array_equal(by_group["arrays"]["water_count"], by_path["arrays"]["water_count"])
    assert universe.trajectory.ts.frame == 0


def test_residue_interactions_and_source_classes_accept_md_objects(tmp_path):
    from crevice.residue_interactions import trajectory_residue_interactions
    from crevice.water_membership import WaterSource
    pdb, xtc = _channel_files(tmp_path, frames=3)
    options = dict(residue_ids=["A:ALA1", "A:ALA11"], selection="protein")
    found_path, times_path = trajectory_residue_interactions(str(pdb), str(xtc), **options)
    universe = mda.Universe(str(pdb), str(xtc))
    found_object, times_object = trajectory_residue_interactions(universe, **options)
    assert np.array_equal(times_path, times_object)
    assert json.dumps(found_path, sort_keys=True, default=str) == json.dumps(found_object, sort_keys=True, default=str)
    trajectory = load_trajectory(universe)
    source = WaterSource(universe.select_atoms("protein"), None, protein_selection="ignored for a group",
                         reference_frame=trajectory.frames[0], frame_times=dict(enumerate(trajectory.times())))
    assert source.topology is universe and len(source.protein_indices) == len(trajectory.frames[0].atoms)

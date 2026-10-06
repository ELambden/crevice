"""Trajectory reader contracts: selection, frame ranges, elements and reports.

The synthetic fixtures are written and read back through MDAnalysis, so these
exercise real GROMACS readers rather than a mock. They do not validate periodic
boundary handling, and they make no claim about any biological system.
"""

from __future__ import annotations

import json

import MDAnalysis as mda  # required runtime dependency: fail, do not skip, if absent
import numpy as np
import pytest

from crevice import trajectory as traj_module
from crevice.trajectory import (
    DEFAULT_MAX_FRAMES,
    inspect_trajectory,
    load_trajectory,
    load_trajectory_report,
)

def build_universe(n_residues=4, per_residue=3, names=("N", "CA", "CB"), segid="A"):
    count = n_residues * per_residue
    universe = mda.Universe.empty(
        count, n_residues=n_residues,
        atom_resindex=np.repeat(np.arange(n_residues), per_residue),
        residue_segindex=np.zeros(n_residues, dtype=int), trajectory=True,
    )
    universe.add_TopologyAttr("names", list(names) * n_residues)
    universe.add_TopologyAttr("resnames", ["ALA"] * n_residues)
    universe.add_TopologyAttr("resids", list(range(1, n_residues + 1)))
    universe.add_TopologyAttr("segids", [segid])
    universe.dimensions = [50, 50, 50, 90, 90, 90]
    universe.atoms.positions = np.arange(count * 3, dtype=np.float32).reshape(count, 3)
    return universe


@pytest.fixture
def gromacs(tmp_path):
    """A five-frame GROMACS pair whose coordinates encode the frame index."""

    universe = build_universe()
    count = len(universe.atoms)
    topology = tmp_path / "system.gro"
    trajectory = tmp_path / "run.xtc"
    universe.atoms.write(str(topology))
    base = np.arange(count * 3, dtype=np.float32).reshape(count, 3)
    with mda.Writer(str(trajectory), n_atoms=count) as writer:
        for index in range(5):
            universe.atoms.positions = base + index * 10.0
            writer.write(universe.atoms)
    return topology, trajectory


# --- inspection -------------------------------------------------------------

def test_inspect_reports_counts_without_loading_frames(gromacs):
    topology, trajectory = gromacs
    info = inspect_trajectory(topology, trajectory)
    assert info["frame_count"] == 5
    assert info["selected_atom_count"] == 12
    assert info["selected_residue_count"] == 4
    assert info["total_atom_count"] == 12
    # GRO carries no segment ID, so the reader assigns its default.
    assert info["segment_ids"] == ["SYSTEM"]
    assert info["box_dimensions"][:3] == [50.0, 50.0, 50.0]


def test_inspect_honours_the_selection(gromacs):
    topology, trajectory = gromacs
    info = inspect_trajectory(topology, trajectory, selection="name CA")
    assert info["selected_atom_count"] == 4
    assert info["total_atom_count"] == 12


def test_inspect_is_json_serializable(gromacs):
    topology, trajectory = gromacs
    assert json.loads(json.dumps(inspect_trajectory(topology, trajectory)))["frame_count"] == 5


# --- frame ranges -----------------------------------------------------------

def test_all_frames_load_by_default(gromacs):
    topology, trajectory = gromacs
    assert len(load_trajectory(topology, trajectory).frames) == 5


def test_stride_subsamples_at_read_time(gromacs):
    topology, trajectory = gromacs
    loaded, report = load_trajectory_report(topology, trajectory, stride=2)
    assert len(loaded.frames) == 3
    assert [frame.frame_index for frame in loaded.frames] == [0, 2, 4]
    assert report["available_frame_count"] == 5
    assert report["loaded_frame_count"] == 3


def test_start_and_stop_bound_the_range(gromacs):
    topology, trajectory = gromacs
    loaded = load_trajectory(topology, trajectory, start=1, stop=4)
    assert [frame.frame_index for frame in loaded.frames] == [1, 2, 3]


def test_stop_beyond_the_end_is_clamped(gromacs):
    topology, trajectory = gromacs
    assert len(load_trajectory(topology, trajectory, stop=99).frames) == 5


def test_frames_carry_their_source_index_not_a_counter(gromacs):
    """A strided frame must report where it came from in the input."""
    topology, trajectory = gromacs
    loaded = load_trajectory(topology, trajectory, start=1, stride=3)
    assert [frame.frame_index for frame in loaded.frames] == [1, 4]


@pytest.mark.parametrize("kwargs,message", [
    ({"stride": 0}, "stride must be at least 1"),
    ({"start": -1}, "start must be non-negative"),
    ({"start": 99}, "beyond the last frame"),
    ({"start": 3, "stop": 2}, "must be greater than start"),
])
def test_invalid_frame_ranges_are_rejected(gromacs, kwargs, message):
    topology, trajectory = gromacs
    with pytest.raises(ValueError, match=message):
        load_trajectory(topology, trajectory, **kwargs)


def test_max_frames_guard_names_the_remedy(gromacs):
    topology, trajectory = gromacs
    with pytest.raises(ValueError, match="raise stride"):
        load_trajectory(topology, trajectory, max_frames=2)


def test_max_frames_can_be_disabled(gromacs):
    topology, trajectory = gromacs
    assert len(load_trajectory(topology, trajectory, max_frames=None).frames) == 5


def test_default_guard_is_documented_and_generous(gromacs):
    topology, trajectory = gromacs
    assert DEFAULT_MAX_FRAMES >= 100
    assert len(load_trajectory(topology, trajectory).frames) == 5


# --- selection --------------------------------------------------------------

def test_selection_restricts_the_loaded_atoms(gromacs):
    topology, trajectory = gromacs
    loaded, report = load_trajectory_report(topology, trajectory, selection="name CA")
    assert len(loaded.frames[0].atoms) == 4
    assert report["selection"] == "name CA"
    assert report["selected_atom_count"] == 4


def test_empty_selection_is_rejected(gromacs):
    topology, trajectory = gromacs
    with pytest.raises(ValueError, match="matched no atoms"):
        load_trajectory(topology, trajectory, selection="name ZZ")


def test_invalid_selection_syntax_names_the_selection(gromacs):
    topology, trajectory = gromacs
    with pytest.raises(ValueError, match="rejected selection"):
        load_trajectory(topology, trajectory, selection="not a valid ((")


# --- missing and malformed inputs -------------------------------------------

def test_missing_topology_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="Topology"):
        load_trajectory(tmp_path / "absent.gro")


def test_missing_trajectory_is_reported(gromacs, tmp_path):
    topology, _ = gromacs
    with pytest.raises(FileNotFoundError, match="Trajectory"):
        load_trajectory(topology, tmp_path / "absent.xtc")


def test_unreadable_topology_names_the_file(tmp_path):
    broken = tmp_path / "broken.gro"
    broken.write_text("not a gro file\n")
    with pytest.raises(ValueError, match="MDAnalysis could not read"):
        load_trajectory(broken)


# --- atom fidelity ----------------------------------------------------------

def test_coordinates_track_the_frame(gromacs):
    topology, trajectory = gromacs
    loaded = load_trajectory(topology, trajectory)
    first = loaded.frames[0].atoms[0]
    third = loaded.frames[2].atoms[0]
    assert third.x == pytest.approx(first.x + 20.0, abs=1e-2)


def test_elements_use_the_shared_inference_rule(gromacs):
    """A reader supplying no elements must not fall back to the first letter."""
    topology, trajectory = gromacs
    loaded, report = load_trajectory_report(topology, trajectory)
    assert report["element_source"] == "name inference"
    assert [atom.element for atom in loaded.frames[0].atoms[:3]] == ["N", "C", "C"]


def test_charmm_hydrogen_names_are_recognised(tmp_path):
    """CHARMM HT1/HT2 style names are hydrogen, not something beginning with H."""
    universe = build_universe(n_residues=2, per_residue=3, names=("N", "HT1", "HT2"))
    topology = tmp_path / "h.gro"
    universe.atoms.write(str(topology))
    loaded = load_trajectory(topology)
    assert [atom.element for atom in loaded.frames[0].atoms[:3]] == ["N", "H", "H"]
    assert len(loaded.frames[0].selected_atoms(include_hydrogen=False)) == 2


def test_segment_ids_stand_in_for_absent_chain_ids(gromacs):
    """GRO has no chain or segment column, so every atom shares one label."""
    topology, trajectory = gromacs
    loaded, report = load_trajectory_report(topology, trajectory)
    assert report["chain_id_source"] == "segid"
    assert report["chain_ids"] == ["SYSTEM"]
    assert all(atom.chain_id == "SYSTEM" for atom in loaded.frames[0].atoms)


def test_residue_identity_is_preserved(gromacs):
    topology, trajectory = gromacs
    loaded = load_trajectory(topology, trajectory)
    residues = loaded.frames[0].residues()
    assert len(residues) == 4
    assert [residue.key.resid for residue in residues] == [1, 2, 3, 4]


# --- times ------------------------------------------------------------------

def test_times_account_for_stride_and_start(gromacs):
    topology, trajectory = gromacs
    loaded, report = load_trajectory_report(topology, trajectory, start=1, stride=2)
    step = report["source_time_step_ps"]
    if step is None:
        pytest.skip("this writer recorded no time step")
    assert report["loaded_time_step_ps"] == step * 2
    assert loaded.times()[0] == pytest.approx(step * 1)
    assert loaded.times()[1] == pytest.approx(step * 3)


# --- report -----------------------------------------------------------------

def test_report_records_reader_provenance_and_limits(gromacs):
    topology, trajectory = gromacs
    _, report = load_trajectory_report(topology, trajectory)
    assert report["reader"] == "mdanalysis"
    assert report["reader_version"] == mda.__version__
    assert report["periodic_boundary_handling"] == "check"
    assert report["periodic_check_scope"].startswith("known bonds")
    assert report["completeness_validated"] is False
    assert json.loads(json.dumps(report))["loaded_frame_count"] == 5


def test_missing_mdanalysis_names_the_required_dependency(gromacs, monkeypatch):
    import builtins

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "MDAnalysis":
            raise ImportError("no MDAnalysis")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match=r"MDAnalysis is a required CREVICE dependency"):
        load_trajectory(*gromacs)

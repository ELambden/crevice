"""Synthetic regression: diagnostics for cast results that can flip with the grid.

The dominant rolling-probe cast keeps at most ``max_components`` regions and
groups cores whose probe-radius supports touch. Both decisions are made on a
grid, so a near tie between regions or a neck whose clearance is within about
one grid tolerance of the probe can change the reported volume
discontinuously with spacing or grid offset (1GRM.pdb with its D-residues
dropped: 181.9 A^3 at 0.5 A, 433.3 A^3 at 0.25 A). These checks only flag such
runs; no result changes. PDB files also record polymer residues such as
D-amino acids as HETATM, which the default "protein" obstacle selection drops;
that is flagged too. Synthetic geometry only; not a real-system validation.
"""
import math

import numpy as np

from crevice.cavity_selection import _grouping_sensitivity
from crevice.models import Atom, StructureFrame
from crevice.rolling import _selection_stability, rolling_probe_cast, select_cast_atoms


def _shell(center, n=250, radius=7.0, start=1):
    golden = math.pi * (3 - math.sqrt(5))
    atoms = []
    for i in range(n):
        y = 1 - (2 * i + 1) / n
        r = math.sqrt(1 - y * y)
        atoms.append(Atom(start + i, "C", "ALA", "A", start + i, center[0] + radius * math.cos(golden * i) * r,
                          center[1] + radius * y, center[2] + radius * math.sin(golden * i) * r, "C"))
    return atoms


def test_two_equal_chambers_with_one_region_requested_are_a_near_tie():
    frame = StructureFrame(tuple(_shell((0, 0, 0)) + _shell((30, 0, 0), start=251)))
    cast = rolling_probe_cast(frame, spacing=0.5, enclosure_fraction=0.9, min_component_volume=50,
                              max_components=1, selection_mode="dominant")
    stability = cast.metadata["selection_stability"]
    assert len(cast.components) == 1
    assert stability["status"] == "sensitive" and stability["near_tie"] is not None
    assert "dropped by the region limit" in stability["message"]
    both = rolling_probe_cast(frame, spacing=0.5, enclosure_fraction=0.9, min_component_volume=50,
                              max_components=2, selection_mode="dominant")
    assert both.metadata["selection_stability"]["status"] == "stable"
    single = rolling_probe_cast(StructureFrame(tuple(_shell((0, 0, 0)))), spacing=0.5, enclosure_fraction=0.9,
                                min_component_volume=50, max_components=1, selection_mode="dominant")
    assert single.metadata["selection_stability"] == {**single.metadata["selection_stability"],
                                                      "status": "stable", "near_tie": None,
                                                      "grid_sensitive_grouping": []}


def _neck_grid(neck_clearance):
    clearance = np.full((24, 3, 3), 1.5)
    clearance[10:14] = neck_clearance
    depth = np.full(clearance.shape, 10.0)
    rays = np.full(clearance.shape, 26)
    groups = {(1, 0): [(1, np.array([[2, 1, 1]]))], (2, 0): [(2, np.array([[20, 1, 1]]))]}
    return groups, clearance, depth, rays


def test_a_neck_within_one_grid_tolerance_of_the_probe_is_reported():
    # probe 0.8 A at 0.5 A spacing: tolerance 0.433 A, so a 0.75 A neck may open or close with the grid.
    report = _grouping_sensitivity(*_neck_grid(0.75), 26, 0.8, 0.5)
    assert report["tolerance_A"] == 0.5 * math.sqrt(3) / 2
    assert report["joined_groups"] == [[1, 2]]
    # A 0.2 A neck is far below the probe and stays closed.
    assert _grouping_sensitivity(*_neck_grid(0.2), 26, 0.8, 0.5)["joined_groups"] == []
    # The tolerance shrinks with the spacing: at 0.05 A the 0.75 A neck is resolved as closed.
    assert _grouping_sensitivity(*_neck_grid(0.75), 26, 0.8, 0.05)["joined_groups"] == []


def test_stability_report_only_flags_groups_touching_retained_regions():
    dominance = {"regions": [{"region_id": 1, "core_ids": [3]}, {"region_id": 2, "core_ids": [5]}],
                 "grouping_sensitivity": {"joined_groups": [[5, 7]]}}
    decisions = [{"source_component_id": 1, "volume_A3": 300.0, "decision": "retained"},
                 {"source_component_id": 2, "volume_A3": 100.0, "decision": "component_limit"}]
    assert _selection_stability(decisions, [1], dominance, 0.5)["status"] == "stable"
    report = _selection_stability(decisions, [1, 2], dominance, 0.5)
    assert report["status"] == "sensitive" and report["grid_sensitive_grouping"] == [[5, 7]]


def _residue(chain, number, name, hetero, x0, start):
    names = ("N", "CA", "C", "O")
    offsets = ((0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (2.5, 1.1))
    return [Atom(start + k, atom, name, chain, number, x0 + dx, dy, 0.0, atom[0], hetero=hetero)
            for k, (atom, (dx, dy)) in enumerate(zip(names, offsets))]


def test_peptide_bonded_hetatm_residues_are_reported_when_dropped():
    # ALA - DLE (HETATM) - ALA with 1.33 A peptide bonds (C at x0+2.0, next N at x0+3.33), plus a water.
    atoms = (_residue("A", 1, "ALA", False, 0.0, 1) + _residue("A", 2, "DLE", True, 3.33, 5)
             + _residue("A", 3, "ALA", False, 6.66, 9)
             + [Atom(13, "O", "HOH", "A", 100, 20.0, 0.0, 0.0, "O", hetero=True)])
    frame = StructureFrame(tuple(atoms))
    _, report = select_cast_atoms(frame, "protein")
    assert report["excluded_residue_atom_counts"] == {"DLE": 4, "HOH": 1}
    assert report["polymer_hetero_residues_excluded"] == {"DLE": 1}
    _, report = select_cast_atoms(frame, "all")
    assert report["polymer_hetero_residues_excluded"] == {}

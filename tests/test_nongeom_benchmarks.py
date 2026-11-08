"""Fast checks of the non-geometric validation helpers in ``benchmarks/nongeom`` (software tests, synthetic data).

The long studies (coverage, ensembles, real trajectories) are scripts; these
tests keep their building blocks honest: the independent reference counter
agrees event by event with CREVICE on random walkers, the exact ring MSD has
the right limits, the scripted catalogue matches its hand-written answers,
and the MD driver runs end to end on a small synthetic system.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")
NONGEOM = Path(__file__).resolve().parents[1] / "benchmarks" / "nongeom"
sys.path.insert(0, str(NONGEOM))
common = pytest.importorskip("nongeom_common")

from crevice.permeation import CylinderChannel, analyze_permeation  # noqa: E402
from crevice.uncertainty import blocking_standard_error  # noqa: E402


def test_reference_counter_and_collective_coordinate_match_crevice_on_random_walkers():
    rng = np.random.default_rng(4)
    box_z, n, frames = 20.0, 30, 3000
    z = common.ring_walkers(n, frames, D=0.125, dt=1.0, box_z=box_z, rng=rng)
    pos = np.zeros((frames, n, 3))
    pos[:, :, 2] = z
    u = common.memory_universe(["TIP3"] * n, pos, box=(30.0, 30.0, box_z), dt=1.0)
    ch = CylinderChannel("pore", (0, 0, 1), 5.0, 15.0, point=(15.0, 15.0, 10.0), radius=3.0)
    result = analyze_permeation(u, species=["TIP3"], channels=[ch], align=False, msd_blocks=5)
    s, img, lat, _ = common.reference_axial_lateral(u, np.arange(6, 6 + n), (15.0, 15.0, 10.0))
    states = common.reference_states(img, lat, 5.0, 15.0, 3.0)
    ref = common.reference_events(states, s, 5.0, 15.0, s_image=img, convention="compartment")
    assert len(ref) > 50
    assert common.crevice_event_keys(result, 6, list(range(frames))) == common.event_keys(ref)
    n_ref = common.reference_collective_coordinate(s, 5.0, 15.0, box_z)
    assert np.max(np.abs(np.asarray(result.frames["series"]["pore/TIP3"]["n"]) - n_ref)) < 1e-9


def test_exact_ring_msd_limits_and_single_file_relabelling():
    N, D, L, P = 50, 0.3, 10.0, 25.0
    short = common.ring_n_msd([1e-4], N, D, L, P)[0] / 2e-4
    long = common.ring_n_msd([1e7], N, D, L, P)[0] / 2e7
    assert short == pytest.approx(N * D / (L * P), rel=1e-3)
    assert long == pytest.approx(N * D / P ** 2, rel=1e-3)
    x = common.ring_walkers(8, 500, D=0.5, dt=1.0, box_z=20.0, rng=np.random.default_rng(2))
    y = common.single_file_relabel(x, 20.0)
    assert np.allclose(y.sum(axis=1), x.sum(axis=1))
    assert np.all(np.diff(y, axis=1) > 0) and np.all(y[:, -1] < y[:, 0] + 20.0)
    assert np.abs(np.diff(y, axis=0)).max() < 10.0


def test_scripted_permeation_catalogue_matches_its_expected_counts():
    import permeation_counting
    rows = permeation_counting.run_cases()
    known = {"n_entry_through_opposite_periodic_face"}   # documented n(t) limitation, see the README
    bad = [r["case"] for r in rows if r["case"] not in known and not r["crevice_matches_expected"]]
    assert not bad
    assert all(r.get("reference_matches_expected", True) for r in rows)


def test_jonsson_transcription_equals_crevice_blocking_on_power_of_two_series():
    import uncertainty_reference_worker as worker
    for phi in (0.0, 0.8):
        x = common.ar1(2048, phi, np.random.default_rng(int(phi * 10)))
        assert blocking_standard_error(x)["standard_error"] == pytest.approx(worker.jonsson_block(x.copy()), rel=1e-12)


def test_md_channel_analysis_runs_on_a_synthetic_system(tmp_path):
    import md_channel_analysis
    frames = 40
    box = (30.0, 30.0, 40.0)
    zs_up = np.clip(np.linspace(-12, 12, frames), -12, 12)
    paths = [zs_up, zs_up[::-1], np.full(frames, -15.0)]
    resnames = ["TIP3", "SOD", "CLA"]
    pos = np.zeros((frames, 3, 3))
    for j, p in enumerate(paths):
        pos[:, j, 2] = p
    u = common.memory_universe(resnames, pos, box=box, dt=100.0)
    top, traj = common.write_md_files(u, tmp_path)
    systems = {"toy": {"species": ["TIP3", "SOD", "CLA"], "alignment_selection": "name CA", "efield_V_per_nm": 0.05,
                       "pores": [{"name": "pore", "point_selection": "name CA", "axis": [0, 0, 1],
                                  "extent_relative_A": [-5.0, 5.0], "radius_A": 3.0}],
                       "sensitivity": {"radius_A": [2.0], "extent_delta_A": [2.0]}}}
    (tmp_path / "systems.json").write_text(json.dumps(systems))
    lit = tmp_path / "lit.csv"
    lit.write_text("quantity,system,value,uncertainty,units,conditions,method,source,doi,verified,compare_with\n"
                   "conductance_per_trimer_nS,toy,1.0,0.1,nS,,,,,,toy\n")
    summary = md_channel_analysis.main(["toy", "--topology", str(top), "--trajectory", str(traj), "--out-dir",
                                        str(tmp_path / "out"), "--systems", str(tmp_path / "systems.json"),
                                        "--literature", str(lit), "--windows", "1"])
    rows = list(__import__("csv").DictReader(line for line in (tmp_path / "out" / "toy_permeation_definitions.csv")
                                             .read_text().splitlines() if not line.startswith("#")))
    for convention in ("compartment", "strict"):
        base = {r["species"]: r for r in rows if r["definition"] == "base" and r["convention"] == convention}
        assert (base["TIP3"]["crevice_up"], base["TIP3"]["reference_up"]) == ("1", "1")
        assert (base["SOD"]["crevice_down"], base["SOD"]["reference_down"]) == ("1", "1")
    assert all(r["events_identical_to_reference_all_species"] == "True" for r in rows
               if r.get("events_identical_to_reference_all_species"))
    cur = summary["currents"]
    # One Na+ moves 24 A down in a 40 A box over 3.9 ns: displacement current -24/40 e / 3900 ps.
    assert cur["displacement_current_pA"] == pytest.approx(-24 / 40 * 1.602176634e5 / 3900, rel=1e-3)
    assert cur["crevice_count_current_sum_over_pores_pA"] == pytest.approx(-1.602176634e5 / 3900, rel=1e-6)
    assert math.isfinite(summary["literature_comparison"][0]["ratio"])


def test_md_channel_reference_ignores_periodic_face_crossings_and_counts_all_outside_lumen_transits(tmp_path):
    # Regression (OmpF trimer, 28 Oct 2026): the reference pass flagged every ion that crossed the periodic
    # face in the bulk as an unresolved crossing (unwrapped step not minimum-imaged), and listed outside-lumen
    # transits only for molecules that had been in that lumen, so no event list matched CREVICE's.
    import csv

    import md_channel_analysis
    frames = 40
    box = (30.0, 30.0, 40.0)
    through = np.linspace(-12, 12, frames)
    pos = np.zeros((frames, 3, 3))
    pos[:, 0, 2] = through                      # K+ through the pore (up)
    pos[:, 1, 2] = through                      # K+ through the slab 10 A beside the lumen (outside_lumen)
    pos[:, 1, 0] = 10.0
    pos[:, 2, 2] = np.linspace(-7.0, -23.0, frames)    # Cl- from near the mouth across the periodic face (no event)
    pos[:, 2, 0] = 10.0
    u = common.memory_universe(["POT", "POT", "CLA"], pos, box=box, dt=100.0)
    top, traj = common.write_md_files(u, tmp_path)
    systems = {"toy": {"species": ["POT", "CLA"], "alignment_selection": "name CA", "efield_V_per_nm": 0.05,
                       "pores": [{"name": "pore", "point_selection": "name CA", "axis": [0, 0, 1],
                                  "extent_relative_A": [-5.0, 5.0], "radius_A": 3.0}]}}
    (tmp_path / "systems.json").write_text(json.dumps(systems))
    lit = tmp_path / "lit.csv"
    lit.write_text("quantity,system,value,uncertainty,units,conditions,method,source,doi,verified,compare_with\n")
    md_channel_analysis.main(["toy", "--topology", str(top), "--trajectory", str(traj), "--out-dir",
                              str(tmp_path / "out"), "--systems", str(tmp_path / "systems.json"),
                              "--literature", str(lit), "--windows", "1"])
    rows = list(csv.DictReader(line for line in (tmp_path / "out" / "toy_permeation_definitions.csv")
                               .read_text().splitlines() if not line.startswith("#")))
    base = {(r["species"], r["convention"]): r for r in rows if r["definition"] == "base"}
    k, cl = base[("POT", "compartment")], base[("CLA", "compartment")]
    assert (k["crevice_up"], k["reference_up"]) == ("1", "1")
    assert (k["crevice_outside_lumen"], k["reference_outside_lumen"]) == ("1", "1")
    assert (cl["crevice_unresolved"], cl["reference_unresolved"]) == ("0", "0")
    assert all(r["events_identical_to_reference_all_species"] == "True" for r in rows
               if r.get("events_identical_to_reference_all_species"))

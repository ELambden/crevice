"""Automatic probe of crevice paths (``probe_radius="auto"``) on synthetic carved blocks.

Software tests on analytic geometry only: a narrow tube that only small probes
pass, a wide tube with a narrow side slit that only small probes see, and a
tube that is narrow on one side of the start and wide on the other. No
biological result is implied.
"""
import csv
import json

import numpy as np
import pytest

from crevice.paths import (PROBE_LADDER, FreeSpace, _exposure_levels, _probe_arg, _ray_escapes, find_paths,
                           write_paths_outputs)
from test_paths import ARC, CARVE, _block, _distance_to_polyline

NARROW = np.array([[0.75, 0.75, -30.0], [0.75, 0.75, 30.0]])      # nearest atoms 2.37 Å away: clearance ~0.6 Å
WIDE = np.array([[0.0, 0.0, -30.0], [0.0, 0.0, 30.0]])             # clearance ~2.3 Å
SLIT = np.array([[0.0, 0.75, 0.75], [30.0, 0.75, 0.75]])           # narrow side slit along +x


def _strip_timing(record):
    """to_dict() without wall-clock fields (they differ between runs)."""
    sel = record["settings"].get("probe_selection")
    if sel:
        sel.pop("timing_s", None)
        for row in sel["candidates"]:
            row.pop("search_s", None)
    return record


@pytest.fixture(scope="module")
def narrow():
    frame, _ = _block((-9, 9), (-9, 9), (-9, 9), lambda p: _distance_to_polyline(p, NARROW) < 2.0)
    return frame


@pytest.fixture(scope="module")
def slit():
    carve = lambda p: (_distance_to_polyline(p, WIDE) < CARVE) | (_distance_to_polyline(p, SLIT) < 2.0)
    frame, _ = _block((-15, 15), (-9, 9), (-12, 12), carve)
    return frame


@pytest.fixture(scope="module")
def two_sided():
    upper = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 30.0]])
    lower = np.array([[0.75, 0.75, -30.0], [0.75, 0.75, 0.0]])
    slit = np.array([[0.0, 0.75, 5.25], [30.0, 0.75, 5.25]])
    carve = lambda p: ((_distance_to_polyline(p, upper) < CARVE) | (_distance_to_polyline(p, lower) < 2.0)
                       | (_distance_to_polyline(p, slit) < 2.0))
    frame, _ = _block((-15, 15), (-9, 9), (-12, 15), carve)
    return frame


def test_exposure_levels_match_single_probe_rays():
    frame, _ = _block((-21, 21), (-9, 9), (-10, 10), lambda p: _distance_to_polyline(p, ARC) < CARVE)
    atoms = frame.selected_atoms(include_hydrogen=False, include_hetero=False)
    space = FreeSpace(atoms, np.asarray([a.coord for a in atoms]), 0.5, 6.0, 20_000_000)
    probes = [p for p in PROBE_LADDER if p < 6.0]
    for required in (20, 1, 26):
        levels = _exposure_levels(space.clearance, space.in_bulk, probes, required)
        for j, p in enumerate(probes):
            blocked = space.clearance <= p
            single = (_ray_escapes(blocked, space.in_bulk) > 26 - required) & ~blocked & ~space.in_bulk
            assert np.array_equal(single, (levels > j) & (space.clearance > p) & ~space.in_bulk)


def test_narrow_tube_needs_small_probe(narrow):
    fixed = find_paths(narrow, seed=(0.75, 0.75, 0.0), axis="z", probe_radius=0.8)
    assert fixed.status == "unresolved_seed_blocked"
    auto = find_paths(narrow, seed=(0.75, 0.75, 0.0), axis="z", workers=1)
    assert auto.status == "resolved"
    assert sorted(b.side for b in auto.branches) == ["+", "-"]
    selection = auto.settings["probe_selection"]
    assert selection["mode"] == "auto"
    assert {v["probe_A"] for v in selection["chosen_by_side"].values()} == {0.5}
    assert auto.settings["probe_radius_A"] == 0.5
    assert all(0.55 < b.bottleneck_radius < 0.7 for b in auto.branches)
    # Probes at or above the widest route width are pruned, not searched.
    statuses = {(r["probe_A"], r["side"]): r["status"] for r in selection["candidates"]}
    assert statuses[(0.5, "+")] == "resolved" and statuses[(0.8, "+")] in {"pruned", "unresolved_seed_blocked"}
    assert selection["timing_s"]["searches_pruned"] > 0


def test_wide_tube_ignores_small_probe_slit(slit):
    small = find_paths(slit, seed=(0.0, 0.0, 0.0), probe_radius=0.5)
    assert len(small.branches) == 3                       # main tube both ways + the slit
    auto = find_paths(slit, seed=(0.0, 0.0, 0.0), workers=1)
    chosen = auto.settings["probe_selection"]["chosen_by_side"]["all"]
    assert chosen["probe_A"] == 0.65 and chosen["routes"] == 2
    assert len(auto.branches) == 2
    assert all(b.bottleneck_radius > 2.0 for b in auto.branches)
    rows = {r["probe_A"]: r for r in auto.settings["probe_selection"]["candidates"]}
    assert rows[0.5]["routes"] == 3 and rows[0.5]["group"] != rows[0.65]["group"]


def test_auto_equals_explicit_run_at_the_chosen_probe(slit):
    auto = find_paths(slit, seed=(0.0, 0.0, 0.0), workers=1).to_dict()
    fixed = find_paths(slit, seed=(0.0, 0.0, 0.0), probe_radius=0.65).to_dict()
    for branch in auto["branches"]:
        assert branch.pop("probe_radius_A") == 0.65
    assert auto["branches"] == fixed["branches"]
    assert auto["segments"] == fixed["segments"] and auto["junctions"] == fixed["junctions"]
    assert auto["start"] == fixed["start"]
    assert "probe_selection" not in fixed["settings"] and fixed["settings"]["probe_radius_A"] == 0.65


def test_probe_chosen_per_side(two_sided, tmp_path):
    result = find_paths(two_sided, seed=(0.0, 0.0, 3.75), axis="z", workers=1)
    chosen = result.settings["probe_selection"]["chosen_by_side"]
    assert chosen["-"]["probe_A"] == 0.5
    assert chosen["+"]["probe_A"] == 0.65 and chosen["+"]["routes"] == 1   # the slit is seen only below 0.65 Å
    by_side = {b.side: b.metadata["probe_radius_A"] for b in result.branches}
    assert by_side == {"+": chosen["+"]["probe_A"], "-": 0.5}
    files = {p.name for p in write_paths_outputs(result, tmp_path, prefix="t", figure=True, scenes=False)}
    assert "t_path_probes.csv" in files
    rows = list(csv.DictReader((tmp_path / "t_paths.csv").open()))
    assert {r["side"]: float(r["probe_radius_A"]) for r in rows} == by_side
    probes = list(csv.DictReader((tmp_path / "t_path_probes.csv").open()))
    assert sum(r["chosen"] == "True" for r in probes) == 2


def test_parallel_equals_serial(two_sided):
    serial = _strip_timing(find_paths(two_sided, seed=(0.0, 0.0, 3.75), axis="z", workers=1).to_dict())
    parallel = _strip_timing(find_paths(two_sided, seed=(0.0, 0.0, 3.75), axis="z", workers=3).to_dict())
    assert json.dumps(serial, sort_keys=True, default=float) == json.dumps(parallel, sort_keys=True, default=float)


def test_no_route_reports_probes_tried():
    frame, _ = _block((-12, 12), (-12, 12), (-12, 12), lambda p: np.linalg.norm(p, axis=1) < 5.0)
    result = find_paths(frame, seed=(0.0, 0.0, 0.0), workers=1)
    assert result.status == "unresolved_no_bulk_path"
    assert "any candidate probe" in result.start["reason"]
    selection = result.settings["probe_selection"]
    assert selection["chosen_probe_A"] is None
    assert selection["candidates"] and selection["chosen_by_side"]["all"]["probe_A"] is None
    blocked = find_paths(frame, seed=(11.0, 11.0, 11.0), seed_tolerance=0.5, workers=1)
    assert blocked.status == "unresolved_seed_blocked" and "every candidate probe" in blocked.start["reason"]


def test_probe_argument_and_validation(narrow):
    assert _probe_arg("auto") == "auto" and _probe_arg(" AUTO ") == "auto" and _probe_arg("0.5") == 0.5
    import argparse
    for bad in ("x", "-1", "nan"):
        with pytest.raises(argparse.ArgumentTypeError):
            _probe_arg(bad)
    with pytest.raises(ValueError):
        find_paths(narrow, seed=(0.75, 0.75, 0.0), probe_radius="big")
    with pytest.raises(ValueError):
        find_paths(narrow, seed=(0.75, 0.75, 0.0), workers=0)
    with pytest.raises(ValueError):
        find_paths(narrow, seed=(0.75, 0.75, 0.0), bulk_radius=0.4)


def test_cli_auto_prints_choice(narrow, tmp_path, capsys):
    from crevice.cli import build_parser
    from crevice.viewer import _frame_pdb
    pdb = tmp_path / "narrow.pdb"
    pdb.write_text(_frame_pdb(narrow))
    args = build_parser().parse_args(["paths", str(pdb), "-o", str(tmp_path / "out"), "--seed", "0.75,0.75,0",
                                      "--axis", "z", "--workers", "1", "--no-scenes"])
    assert args.probe_radius == "auto"
    args.func(args)
    out = capsys.readouterr().out
    assert "probe=auto chosen: + 0.5 Å, - 0.5 Å" in out
    data = json.loads((tmp_path / "out" / "narrow_paths.json").read_text())
    assert data["settings"]["probe_selection"]["chosen_by_side"]["+"]["probe_A"] == 0.5
    assert (tmp_path / "out" / "narrow_path_probes.csv").exists()

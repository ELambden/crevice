"""Pure helpers of ``benchmarks/timing_trajectory.py`` (no comparator or trajectory needed)."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1] / "benchmarks"
sys.path.insert(0, str(BENCH))
timing = pytest.importorskip("timing_trajectory")


def test_parse_cpus_ranges_and_lists():
    assert timing.parse_cpus(None) is None
    assert timing.parse_cpus("") is None
    assert timing.parse_cpus("0") == [0]
    assert timing.parse_cpus("3,0-1,1") == [0, 1, 3]


def test_descendants_follow_the_process_tree_only():
    parents = {10: 1, 11: 10, 12: 11, 13: 10, 20: 1, 21: 20}
    assert timing.descendants(parents, 10) == {10, 11, 12, 13}
    assert timing.descendants(parents, 21) == {21}


def _region():
    return {"geometry": {"spacing_A": 0.2, "probe_radius_A": 1.4, "region_margin_A": 2.0, "profile_padding_A": 4.0,
                         "axis": [0.0, 0.0, 1.0], "grid_phase": [0.0, 0.0, 0.0], "max_grid_points": 8000000},
            "system": {"protein_selection": "protein and not name H*", "obstacle_selection": "not resname TIP3",
                       "pbc": "check", "topology": {"path": "top.gro"}, "trajectory": {"path": "traj.xtc"}},
            "reference": {"asset": {"path": "ref.dx"}}}


def test_region_spacing_args_keep_the_definition_but_the_spacing():
    args = timing.region_spacing_args(_region(), Path("t.gro"), Path("x.xtc"), Path("r.dx"), 0.35)
    assert args[:4] == ["trajectory", "--cavities", "t.gro", "x.xtc"]
    pairs = dict(zip(args[4::2], args[5::2]))
    assert pairs["--spacing"] == "0.35"
    assert pairs["--geometry-mode"] == "reference-region"
    assert pairs["--probe-radius"] == "1.4" and pairs["--region-margin"] == "2.0"
    assert pairs["--obstacle-selection"] == "not resname TIP3"
    assert "--axis=0.0,0.0,1.0" in args and "--grid-phase=0.0,0.0,0.0" in args


def test_load_system_resolves_paths_and_region_defaults(tmp_path):
    work = tmp_path / "work"
    (work / "ex").mkdir(parents=True)
    (work / "ex" / "region.json").write_text(json.dumps(_region()))
    system = tmp_path / "sys.json"
    system.write_text(json.dumps({"label": "demo", "region": "{work_dir}/ex/region.json"}))
    data = timing.load_system(system, work)
    assert data["label"] == "demo"
    assert data["topology"] == (work / "ex" / "top.gro").resolve()
    assert data["trajectory"] == (work / "ex" / "traj.xtc").resolve()
    assert data["reference_dx"] == (work / "ex" / "ref.dx").resolve()
    assert data["obstacle_selection"] == "not resname TIP3"
    over = timing.load_system(system, work, {"trajectory": "other.xtc", "label": None})
    assert over["trajectory"] == (tmp_path / "other.xtc").resolve() and over["label"] == "demo"


def test_build_cases_skip_missing_tools_and_follow_selection(tmp_path, monkeypatch):
    monkeypatch.setattr(timing, "find_tool", lambda name: None)
    monkeypatch.setattr(timing, "hole_python", lambda: None)
    system = {"selection": "protein", "topology": Path("t.gro"), "trajectory": Path("x.xtc"),
              "region": Path("r.json"), "region_data": _region(), "reference_dx": Path("r.dx"),
              "obstacle_selection": "not resname TIP3",
              "profile": {"axis": [0, 0, 1], "origin": [1, 2, 3], "search_radius": 12}}
    cases = timing.build_cases(system, tmp_path, frames=50, spacings=[0.5, 0.7], workers=2, selected=set(timing.CASES))
    names = [c["case"] for c in cases]
    assert names == ["crevice_profile", "hole", "crevice_region", "crevice_region_0.5", "crevice_region_0.7",
                     "mdpocket_input", "mdpocket"]
    by = {c["case"]: c for c in cases}
    assert by["hole"]["skip"] and by["mdpocket"]["skip"] and not by["crevice_region"]["skip"]
    assert "--stop" in by["crevice_profile"]["cmd"] and "50" in by["crevice_profile"]["cmd"]
    for name in ("crevice_profile", "crevice_region", "crevice_region_0.5"):
        cmd = by[name]["cmd"]  # the CLI refuses more than 500 frames unless --max-frames allows them
        assert cmd[cmd.index("--max-frames") + 1] == "50"
    cmd = by["crevice_region_0.5"]["cmd"]
    assert cmd[cmd.index("--workers") + 1] == "2" and cmd[cmd.index("--spacing") + 1] == "0.5"
    only = timing.build_cases(system, tmp_path, frames=None, spacings=[0.5], workers=1, selected={"crevice_region"})
    assert [c["case"] for c in only] == ["crevice_region"] and "--stop" not in only[0]["cmd"]


def test_frames_written_counts_results(tmp_path):
    (tmp_path / "x_frames.csv").write_text("# note\nframe_index,profile_status\n0,resolved\n1,unresolved\n2,resolved\n")
    assert timing.frames_written("crevice_profile", tmp_path) == (3, 2)
    (tmp_path / "r_cavity_frames.csv").write_text("frame_index,status,volume_A3\n0,regional_measured,800\n1,unresolved,\n")
    assert timing.frames_written("crevice_region_0.7", tmp_path) == (2, 1)
    (tmp_path / "mdpout_descriptors.txt").write_text("snapshot pock_volume\n1 500\n2 510\n")
    assert timing.frames_written("mdpocket", tmp_path) == (2, 2)
    (tmp_path / "frames_written.txt").write_text("21\n")
    assert timing.frames_written("mdpocket_input", tmp_path) == (21, 21)
    assert timing.frames_written("hole", tmp_path / "missing") == (None, None)


def _row(case, wall, cpu, rss, code=0, frames=100):
    return {"case": case, "tool": case.upper(), "frames": frames, "spacing_A": 0.7, "grid": "volume grid", "cpus": 1,
            "workers": 1, "exit_code": code, "wall_s": wall, "cpu_user_s": cpu, "cpu_sys_s": 0.0, "peak_rss_mb": rss,
            "peak_tree_rss_mb": rss + 1, "frames_written": frames, "frames_with_result": frames,
            "load_1min_before": 0.5, "load_1min_after": 0.7}


def test_summarise_medians_rates_and_markdown():
    rows = [_row("a", 10.0, 9.0, 100), _row("a", 20.0, 19.0, 120), _row("a", 30.0, 29.0, 110),
            _row("b", 5.0, 5.0, 50, code=1),
            {"case": "c", "tool": "C", "status": "skipped", "note": "tool not found"}]
    summary = {r["case"]: r for r in timing.summarise(rows)}
    a = summary["a"]
    assert a["median_wall_s"] == 20.0 and a["median_cpu_s"] == 19.0
    assert a["frames_per_s"] == pytest.approx(5.0) and a["s_per_frame"] == pytest.approx(0.2)
    assert a["max_peak_rss_mb"] == 120 and a["max_peak_tree_rss_mb"] == 121
    assert a["status"] == "ok" and a["failed_runs"] == 0
    assert summary["b"]["status"] == "failed" and summary["b"]["failed_runs"] == 1
    assert summary["b"]["frames_per_s"] is None and summary["b"]["s_per_frame"] is None
    assert summary["c"]["status"] == "skipped"
    table = timing.markdown_table(list(summary.values()))
    lines = table.strip().splitlines()
    assert len(lines) == 5 and all(line.count("|") == 10 for line in lines)
    assert "skipped: tool not found" in table and "(1 failed)" in table


@pytest.mark.skipif(not hasattr(os, "wait4"), reason="os.wait4 (POSIX) needed")
def test_measure_records_time_and_memory(tmp_path):
    res = timing.measure([sys.executable, "-c", "import time; x = bytearray(30 * 2**20); time.sleep(0.5)"],
                         log=tmp_path / "log.txt", cpus=None)
    assert res["exit_code"] == 0
    assert res["wall_s"] >= 0.5
    assert res["peak_rss_mb"] >= 30
    if Path("/proc/self/statm").exists():
        assert res["peak_tree_rss_mb"] >= 30
    assert (tmp_path / "log.txt").read_text().startswith("$ ")


def test_cpu_max_mhz_is_empty_without_cpus_and_never_raises():
    assert timing.cpu_max_mhz(None) == ""
    assert timing.cpu_max_mhz([10 ** 6]) == ""
    value = timing.cpu_max_mhz([0])
    assert value == "" or all(part.isdigit() for part in value.split(";"))


def test_portable_replaces_longest_paths_first():
    text = "/data/ws/bench/x.json /data/ws/a.gro /usr/bin/python"
    out = timing.portable(text, {"/data/ws": "<work-dir>", "/data/ws/bench": "<work>", "/usr/bin/python": "python",
                                 "": "ignored"})
    assert out == "<work>/x.json <work-dir>/a.gro python"

"""Time CREVICE and the comparators on the same benchmark inputs.

Usage (after the comparison scripts have prepared their inputs)::

    python benchmarks/timing_vs_comparators.py [--repeat 3] [--work-dir DIR]

Each command runs as a separate process, as a user would run it, measured as in
``scripts/timing.py`` (wall time, user and system CPU time and peak resident
memory of the process tree's largest process, via ``os.wait4``; POSIX only).
The machine, CPU count and the 1/5/15-minute load average before and after
every case are recorded, because timings on a shared machine depend on load.
Single-threaded BLAS everywhere; CREVICE path searches use ``--workers 1`` so
that every tool runs on one core (CAVER is a Java program and may use more
than one thread; its CPU time is reported).

Cases (when the tool and inputs are available):

* ``profile``: HOLE 2.3.1 vs ``crevice profile`` with the automatic probe and
  with the chosen probe fixed, on 1GRM, 1OED and 3UKM (matched inputs, HOLE
  simple.rad; 3UKM with Bondi radii and lateral exits, as resolved).
* ``tunnels``: CAVER 3.0.3 vs ``crevice paths`` (automatic probe and 0.9 Å),
  on 1CQW, 1MJ5 and 1EA5.
* ``volumes``: pyKVFinder and POVME vs ``crevice cavities --cast`` on 1L83 and 2NNQ.

Output: ``results/timing/timing.csv`` (every run) and
``results/timing/summary.csv`` (median per case).
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BENCH, RESULTS, crevice_env, default_work_dir, find_tool, machine_record,  # noqa: E402
                    python_for_crevice, write_csv)


def measure(cmd, *, cwd=None, env=None, stdin=None, log: Path) -> dict:
    with open(log, "w") as out:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT,
                                stdin=open(stdin) if stdin else None)
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - start
    return {"exit_code": os.waitstatus_to_exitcode(status), "wall_s": wall, "cpu_user_s": usage.ru_utime,
            "cpu_sys_s": usage.ru_stime, "peak_rss_mb": usage.ru_maxrss / 1024}


def portable(text: str, work: Path) -> str:
    """Command line with local paths replaced by placeholders (results are committed)."""
    for path, label in ((str(work.resolve()), "<work-dir>"), (str(work), "<work-dir>"),
                        (str(find_tool("hole") or "\0"), "hole"), (str(find_tool("java") or "\0"), "java"),
                        (str(find_tool("caver") or "\0"), "caver.jar"),
                        (str(find_tool("caver").parent) if find_tool("caver") else "\0", "caver-home"), (python_for_crevice(), "python"),
                        (os.environ.get("CREVICE_KVFINDER_PYTHON", "\0"), "python-kvfinder"),
                        (os.environ.get("CREVICE_POVME_PYTHON", "\0"), "python2.7"),
                        (os.environ.get("CREVICE_POVME_SCRIPT", "\0"), "POVME3.py"), (str(BENCH), "benchmarks")):
        text = text.replace(path, label)
    return text


def cases(work: Path):
    py = python_for_crevice()
    env = crevice_env()
    hole = find_tool("hole")
    for acc, radius_set, extra in (("1GRM", "simple", []), ("1OED", "simple", []), ("3UKM", "bondi", ["--lateral-exits"])):
        cur = work / "channels" / acc
        comp = cur / f"comparison_{radius_set}.json"
        if not comp.exists():
            continue
        axis = json.loads((cur / "axis.json").read_text())
        summ = json.loads(comp.read_text())["summary"]
        rad = cur / f"hole_{radius_set}" / "m.rad"   # the radius file HOLE read (copied there by the comparison)
        if hole is not None and (cur / f"hole_{radius_set}" / "hole.inp").exists():
            yield (f"profile_{acc}", "HOLE 2.3.1", [str(hole)], cur / f"hole_{radius_set}",
                   str(cur / f"hole_{radius_set}" / "hole.inp"), None)
        base = [py, "-m", "crevice", "profile", str(cur / "curated.pdb"), "--offline", "-q",
                "--axis=" + ",".join(f"{v:.10f}" for v in axis["axis_unit"]),
                "--origin=" + ",".join(f"{v:.6f}" for v in axis["origin"]), "--section-spacing", "0.25",
                "--radii", str(rad), *extra, "-o", "/tmp/crevice_timing_profile.json"]
        yield (f"profile_{acc}", "CREVICE profile (auto probe)", base, None, None, env)
        if summ.get("crevice_chosen_probe_A"):
            yield (f"profile_{acc}", f"CREVICE profile (probe {summ['crevice_chosen_probe_A']:g} A fixed)",
                   base + ["--enclosure-radius", f"{summ['crevice_chosen_probe_A']:g}"], None, None, env)
    java, jar = find_tool("java"), find_tool("caver")
    for acc in ("1CQW", "1MJ5", "1EA5"):
        cur = work / "tunnels" / acc
        if not (cur / "curation.json").exists():
            continue
        seed = json.loads((cur / "curation.json").read_text())["seed"]
        pdb = cur / "pdb" / f"{acc}.pdb"
        if java and jar and (cur / "caver" / "config.txt").exists():
            yield (f"tunnels_{acc}", "CAVER 3.0.3", [str(java), "-Xmx4g", "-jar", str(jar), "-home", str(jar.parent),
                                                     "-pdb", str(pdb.parent), "-conf", str(cur / "caver" / "config.txt"),
                                                     "-out", "/tmp/crevice_timing_caver"], None, None, None)
        for probe in ("auto", "0.9"):
            yield (f"tunnels_{acc}", f"CREVICE paths (probe {probe}, 1 worker)",
                   [py, "-m", "crevice", "paths", str(pdb), "--seed=" + ",".join(f"{v:.3f}" for v in seed), "--radii",
                    "bondi", "--no-scenes", "--probe-radius", probe, "--workers", "1", "-q", "-o",
                    "/tmp/crevice_timing_paths"], None, None, env)
    kv = os.environ.get("CREVICE_KVFINDER_PYTHON")
    pv, ps = os.environ.get("CREVICE_POVME_PYTHON"), os.environ.get("CREVICE_POVME_SCRIPT")
    for acc in ("1L83", "2NNQ"):
        cur = work / "pockets" / acc
        if not (cur / "curation.json").exists():
            continue
        seed = json.loads((cur / "curation.json").read_text())["seed"]
        pdb = cur / f"{acc}.pdb"
        s = ",".join(f"{v:.3f}" for v in seed)
        yield (f"volumes_{acc}", "CREVICE cast (defaults)",
               [py, "-m", "crevice", "cavities", "--cast", str(pdb), f"--seed={s}", "--radii", "bondi", "--skip-analysis-views", "-q",
                "--out-dir", "/tmp/crevice_timing_cast"], None, None, env)
        if kv:
            yield (f"volumes_{acc}", "pyKVFinder 0.9.5 (defaults)",
                   [kv, str(BENCH / "kvfinder_worker.py"), str(pdb), s, "/tmp/crevice_timing_kv.json"], None, None, None)
        if pv and ps and (cur / "povme" / "povme.ini").exists():
            penv = dict(os.environ, PYTHONPATH=str(Path(ps).resolve().parent.parent))
            yield (f"volumes_{acc}", "POVME 3.0.34", [pv, ps, "povme.ini"], cur / "povme", None, penv)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--repeat", type=int, default=3)
    args = ap.parse_args(argv)
    machine = machine_record()
    logs = args.work_dir / "timing_logs"
    logs.mkdir(parents=True, exist_ok=True)
    rows = []
    for n, (case, tool, cmd, cwd, stdin, env) in enumerate(cases(args.work_dir)):
        for rep in range(args.repeat):
            load_before = list(os.getloadavg())
            res = measure(cmd, cwd=cwd, env=env, stdin=stdin, log=logs / f"{n:02d}_{rep}.log")
            rows.append({"case": case, "tool": tool, "repeat": rep, **res,
                         "load_1min_before": load_before[0], "load_1min_after": os.getloadavg()[0],
                         "command": portable(" ".join(cmd), args.work_dir)})
            print(case, tool, rep, f"{res['wall_s']:.2f} s", res["exit_code"])
    summary = []
    for key in dict.fromkeys((r["case"], r["tool"]) for r in rows):
        rs = [r for r in rows if (r["case"], r["tool"]) == key]
        summary.append({"case": key[0], "tool": key[1], "repeats": len(rs),
                        "median_wall_s": statistics.median(r["wall_s"] for r in rs),
                        "median_cpu_s": statistics.median(r["cpu_user_s"] + r["cpu_sys_s"] for r in rs),
                        "max_peak_rss_mb": max(r["peak_rss_mb"] for r in rs),
                        "failed_runs": sum(r["exit_code"] != 0 for r in rs),
                        "load_1min_range": f"{min(r['load_1min_before'] for r in rs):.1f}-"
                                           f"{max(r['load_1min_after'] for r in rs):.1f}"})
    note = (f"Machine: {machine.get('cpu_model', '?')}, {machine['cpu_count']} logical CPUs, {machine['platform']}, "
            f"Python {machine['python']}; load average at start {machine.get('load_average_1_5_15')}. "
            "Shared machine: timings depend on load (see load columns).")
    write_csv(args.results / "timing" / "timing.csv", rows, comment=note)
    write_csv(args.results / "timing" / "summary.csv", summary, comment=note)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

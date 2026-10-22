"""Time a fixed set of CREVICE commands and write the results to a CSV file.

Each case is run as a separate ``python -m crevice`` process so that start-up,
input parsing and output writing are included, as a user would see them. For
every run the script records the wall time, the CPU time and the peak resident
memory of the command (POSIX only; the largest single process, which for
commands with worker processes is usually the main process), together with the
CREVICE version, the Python version, the platform, the CPU count and the
settings used.

Cases
-----
``publish_1grm``
    The canonical 1GRM ``publish`` command (samples 101, 0.25 Å sections and
    cast, gemmi parser, automatic axis).
``publish_3ukm``
    ``publish`` on the 3UKM biological assembly with lateral exits and a
    0.55 Å enclosure probe, as in the 3UKM worked example.
``profile_1grm``, ``profile_3ukm``
    ``profile`` with the automatic enclosure probe (3UKM: assembly 1, no
    heteroatoms, lateral exits).
``region_trajectory_glut1``
    ``region-trajectory`` of a prepared region definition over the first
    ``--frames`` frames, with ``--water-membership``.
``cavity_trajectory_glut1``
    ``cavity-trajectory`` (default rolling mode, 0.7 Å lattice so the
    whole-protein grid stays under the default 2 million point guard) of the
    same frames, using the region definition's reference map.
``paths_glut1``
    ``paths`` from the region's landmarks on frame 0, lipids as walls.
``paths_3ukm``
    ``paths`` on the 3UKM assembly from a seed in the selectivity filter.

Inputs
------
The structures are read from ``--cache-dir`` (``1GRM.cif`` and
``3UKM-assembly1.cif``); ``--fetch`` downloads missing ones with
``crevice fetch``. The trajectory cases need a prepared region definition
(``--region``): its topology, trajectory and reference map are taken from the
definition. By default the definition of the GLUT1 excerpt from
``crevice fetch --example glut1-excerpt`` is used when it is present. Cases whose
inputs are missing are reported as skipped.

Examples
--------
Full set, three repeats, outputs kept for inspection::

    python scripts/timing.py --repeat 3 --out timing.csv --keep-outputs runs/

Quick subset::

    python scripts/timing.py --quick --out timing_quick.csv

The timings depend on the machine and on its load; report the median of
several repeats and the CPU count recorded in the file.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

# 3UKM selectivity-filter point and pore axis (assembly 1 coordinates).
UKM_SEED = "20.423964217778888,-33.313221063824905,56.64516420783265"
UKM_AXIS = "-0.0965476084,-0.8077609078,0.5815504064"

CASES = ("publish_1grm", "publish_3ukm", "profile_1grm", "profile_3ukm",
         "region_trajectory_glut1", "cavity_trajectory_glut1", "paths_glut1", "paths_3ukm")
QUICK_CASES = ("publish_1grm", "profile_1grm", "region_trajectory_glut1", "paths_3ukm")
QUICK_FRAMES = 5
DEFAULT_FRAMES = 20

FIELDS = ("case", "repeat", "status", "wall_s", "cpu_user_s", "cpu_sys_s", "peak_rss_mb", "exit_code",
          "crevice_version", "python", "platform", "machine", "cpu_count", "workers", "threads",
          "frames", "command")


@dataclass
class Case:
    name: str
    args: list[str] = field(default_factory=list)
    missing: str = ""


def default_cache_dir() -> Path:
    """The structure cache: ``$CREVICE_WORKSPACE_DIR/pdb`` when set, else ``.crevice/pdb``."""
    workspace = os.environ.get("CREVICE_WORKSPACE_DIR", "")
    for item in workspace.split(os.pathsep):
        if item and (Path(item) / "pdb").is_dir():
            return Path(item) / "pdb"
    return Path(".crevice") / "pdb"


def default_region() -> Path:
    root = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(root) / "crevice" / "examples" / "glut1-excerpt" / "glut1_site_candidate_region.json"


def crevice_version(env: dict[str, str]) -> str:
    code = "import crevice; print(crevice.__version__)"
    result = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
    return result.stdout.strip() or "unknown"


def region_inputs(region: Path) -> dict[str, Path] | None:
    """Topology, trajectory and reference map named by a region definition (paths relative to it)."""
    if not region.is_file():
        return None
    data = json.loads(region.read_text())
    base = region.resolve().parent
    found = {
        "topology": base / data["system"]["topology"]["path"],
        "trajectory": base / data["system"]["trajectory"]["path"],
        "reference": base / data["reference"]["asset"]["path"],
    }
    return found if all(path.is_file() for path in found.values()) else None


def build_cases(names, cache: Path, region: Path, frames: int, workers: int) -> list[Case]:
    grm = cache / "1GRM.cif"
    ukm = cache / "3UKM-assembly1.cif"
    glut = region_inputs(region)
    w = ["--workers", str(workers)]
    cases = []
    for name in names:
        case = Case(name)
        if name == "publish_1grm":
            case.args = ["publish", str(grm), "--offline", "--model", "0", "--parser", "gemmi", "--axis", "auto",
                         "--samples", "101", "--section-spacing", "0.25", "--cast-spacing", "0.25",
                         "--dpi", "100", "--out-dir", "{out}", "--prefix", "1GRM"]
        elif name == "publish_3ukm":
            case.args = ["publish", str(ukm), "--offline", "--exclude-hetero", "--lateral-exits",
                         "--enclosure-radius", "0.55", "--dpi", "150", "--out-dir", "{out}", "--prefix", "3UKM"]
        elif name == "profile_1grm":
            case.args = ["profile", str(grm), "--offline", "--axis", "auto", "-o", "{out}/1GRM_profile.json",
                         "--csv", "{out}/1GRM_profile.csv"]
        elif name == "profile_3ukm":
            case.args = ["profile", str(ukm), "--offline", "--exclude-hetero", "--lateral-exits",
                         "-o", "{out}/3UKM_profile.json", "--csv", "{out}/3UKM_profile.csv"]
        elif name == "paths_3ukm":
            case.args = ["paths", str(ukm), f"--seed={UKM_SEED}", f"--axis={UKM_AXIS}", *w,
                         "-o", "{out}", "--prefix", "3UKM"]
        elif name.endswith("_glut1"):
            if glut is None:
                case.missing = f"region definition or its inputs not found: {region}"
            elif name == "region_trajectory_glut1":
                case.args = ["region-trajectory", str(region), "--out-dir", "{out}", "--stop", str(frames),
                             "--water-membership", *w]
            elif name == "cavity_trajectory_glut1":
                case.args = ["cavity-trajectory", str(glut["topology"]), str(glut["trajectory"]),
                             "--reference-volume-dx", str(glut["reference"]), "--spacing", "0.7", "--stop", str(frames),
                             "--out-dir", "{out}", "--prefix", "glut1", *w]
            elif name == "paths_glut1":
                case.args = ["paths", str(glut["topology"]), "--trajectory", str(glut["trajectory"]),
                             "--frames", "0", "--md-selection", "not resname TIP3", "--include-hetero",
                             "--pbc", "none", "--region", str(region), "--axis", "z", *w,
                             "-o", "{out}", "--prefix", "glut1"]
        else:
            raise SystemExit(f"unknown case {name!r}; choose from {', '.join(CASES)}")
        if not case.missing:
            for path in (grm,) if "1grm" in name else (ukm,) if "3ukm" in name else ():
                if not path.is_file():
                    case.missing = f"{path} not found (use --fetch or --cache-dir)"
        cases.append(case)
    return cases


def run_once(args: list[str], out: Path, env: dict[str, str], log: Path) -> tuple[int, float, float, float, float | None]:
    """Run one command; return exit code, wall, user and system CPU seconds and peak RSS in MiB."""
    command = [sys.executable, "-m", "crevice", *[a.replace("{out}", str(out)) for a in args]]
    out.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as handle:
        start = time.perf_counter()
        process = subprocess.Popen(command, env=env, stdout=handle, stderr=subprocess.STDOUT)
        if hasattr(os, "wait4"):
            _, status, usage = os.wait4(process.pid, 0)
            wall = time.perf_counter() - start
            process.returncode = os.waitstatus_to_exitcode(status)
            # ru_maxrss is KiB on Linux and bytes on macOS.
            scale = 1 / (1024 * 1024) if sys.platform == "darwin" else 1 / 1024
            return process.returncode, wall, usage.ru_utime, usage.ru_stime, usage.ru_maxrss * scale
        code = process.wait()
        return code, time.perf_counter() - start, float("nan"), float("nan"), None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, default=Path("crevice_timing.csv"), help="CSV file to write")
    parser.add_argument("--quick", action="store_true",
                        help=f"Run only {', '.join(QUICK_CASES)} with {QUICK_FRAMES} frames")
    parser.add_argument("--cases", help="Comma-separated subset of: " + ", ".join(CASES))
    parser.add_argument("--repeat", type=int, default=1, help="Runs per case (report the median)")
    parser.add_argument("--frames", type=int, help=f"Trajectory frames (default {DEFAULT_FRAMES}, quick {QUICK_FRAMES})")
    parser.add_argument("--workers", type=int, default=1, help="--workers for the commands that take it")
    parser.add_argument("--threads", type=int, default=1,
                        help="OMP/OpenBLAS/MKL thread count for the child processes")
    parser.add_argument("--cache-dir", type=Path, default=None, help="Directory with 1GRM.cif and 3UKM-assembly1.cif")
    parser.add_argument("--fetch", action="store_true", help="Download missing structures with crevice fetch")
    parser.add_argument("--region", type=Path, default=None,
                        help="Prepared region definition for the GLUT1 cases (default: the glut1-excerpt example)")
    parser.add_argument("--keep-outputs", type=Path, default=None,
                        help="Keep each run's outputs in DIR/<case>/run<N> (default: temporary, deleted)")
    options = parser.parse_args(argv)

    cache = options.cache_dir or default_cache_dir()
    region = options.region or default_region()
    frames = options.frames or (QUICK_FRAMES if options.quick else DEFAULT_FRAMES)
    names = options.cases.split(",") if options.cases else (QUICK_CASES if options.quick else CASES)
    env = dict(os.environ)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[var] = str(options.threads)
    env.setdefault("CREVICE_PROGRESS", "0")

    if options.fetch:
        cache.mkdir(parents=True, exist_ok=True)
        wanted = [("1GRM", None, cache / "1GRM.cif"), ("3UKM", "1", cache / "3UKM-assembly1.cif")]
        for accession, assembly, path in wanted:
            if not path.is_file():
                extra = ["--assembly", assembly] if assembly else []
                subprocess.run([sys.executable, "-m", "crevice", "fetch", accession, "--format", "cif",
                                "--cache-dir", str(cache), *extra], env=env, check=False)

    cases = build_cases(names, cache, region, frames, options.workers)
    common = {
        "crevice_version": crevice_version(env),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "workers": options.workers,
        "threads": options.threads,
    }
    keep = options.keep_outputs
    scratch = Path(tempfile.mkdtemp(prefix="crevice-timing-")) if keep is None else keep
    rows = []
    try:
        for case in cases:
            case_frames = frames if case.name.endswith("trajectory_glut1") else ""
            if case.missing:
                print(f"{case.name}: skipped ({case.missing})", file=sys.stderr)
                rows.append({**common, "case": case.name, "repeat": "", "status": "skipped",
                             "frames": case_frames, "command": case.missing})
                continue
            for repeat in range(1, options.repeat + 1):
                out = scratch / case.name / f"run{repeat}"
                if out.exists():
                    shutil.rmtree(out)
                log = scratch / case.name / f"run{repeat}.log"
                code, wall, user, system, rss = run_once(case.args, out, env, log)
                print(f"{case.name} run {repeat}: {wall:.2f} s, exit {code}", file=sys.stderr)
                rows.append({**common, "case": case.name, "repeat": repeat,
                             "status": "ok" if code == 0 else "failed",
                             "wall_s": f"{wall:.3f}", "cpu_user_s": f"{user:.3f}", "cpu_sys_s": f"{system:.3f}",
                             "peak_rss_mb": "" if rss is None else f"{rss:.1f}", "exit_code": code,
                             "frames": case_frames, "command": "crevice " + " ".join(case.args)})
                if keep is None:
                    shutil.rmtree(out, ignore_errors=True)
    finally:
        if keep is None:
            shutil.rmtree(scratch, ignore_errors=True)

    options.out.parent.mkdir(parents=True, exist_ok=True)
    with open(options.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{'case':26s} {'median wall s':>14s} {'peak RSS MiB':>13s}")
    for case in cases:
        walls = [float(r["wall_s"]) for r in rows if r["case"] == case.name and r["status"] == "ok"]
        peaks = [float(r["peak_rss_mb"]) for r in rows if r["case"] == case.name and r.get("peak_rss_mb")]
        if walls:
            print(f"{case.name:26s} {statistics.median(walls):14.2f} {max(peaks) if peaks else float('nan'):13.1f}")
        else:
            print(f"{case.name:26s} {'skipped' if case.missing else 'failed':>14s}")
    print(f"wrote {options.out}")
    return 1 if any(r["status"] == "failed" for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())

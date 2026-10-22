"""Time CREVICE trajectory analyses, HOLE and MDpocket on the same MD trajectory.

Usage::

    python benchmarks/timing_trajectory.py SYSTEM.json [--frames N] [--repeat 3] [--cpus 0]
        [--spacings 0.35,0.7] [--cases LIST] [--dry-run]

``SYSTEM.json`` describes one trajectory, so a new system needs only a new
file. Paths are relative to the file; ``{work_dir}`` is replaced by the work
directory (``--work-dir``)::

    {
      "label": "glut1-excerpt",
      "description": "Public GLUT1 excerpt (21 frames)",
      "topology": "{work_dir}/examples/glut1-excerpt/glut1_excerpt.gro",
      "trajectory": "{work_dir}/examples/glut1-excerpt/glut1_excerpt.xtc",
      "selection": "protein",
      "region": "{work_dir}/examples/glut1-excerpt/glut1_site_candidate_region.json",
      "profile": {"axis": [0, 0, 1], "origin": [x, y, z], "search_radius": 12}
    }

``region`` (a prepared region definition) enables the volume cases and
``profile`` (channel axis, a point inside the channel in first-frame
coordinates and optional ``search_radius``, ``enclosure_radius`` and
``section_spacing``) the channel-profile cases; ``topology`` and
``trajectory`` default to the region definition's own.

Cases (``--cases``; those whose inputs or tools are missing are skipped):

* ``crevice_profile``: ``crevice trajectory --topology ... --analysis profile``
  (per-frame channel profile; one process, no workers).
* ``hole``: HOLE 2.3.1 on every frame through mdahole2
  (``trajectory_worker.py hole``), same atoms, frames fitted on the same CA
  atoms, HOLE's ``simple.rad``, the same channel point and direction.
* ``crevice_region``: ``crevice trajectory --region DEFINITION`` at the
  definition's own grid spacing.
* ``crevice_region_<s>``: the same measurement at grid spacing ``s`` Å
  (``crevice trajectory --cavities --geometry-mode reference-region`` with the
  definition's selections, probe, margin and axis; one case per
  ``--spacings`` value).
* ``mdpocket_input``: MDAnalysis fit of the frames on the same CA atoms and
  writing of the non-water atoms (``trajectory_worker.py mdpocket-input``);
  MDpocket needs superimposed frames, so this is part of its cost.
* ``mdpocket``: MDpocket (fpocket 4.1.5) pocket characterisation of the
  reference region (``--selected_pocket``, the region's grid points on a 1 Å
  lattice, as in ``mdpocket_vs_region.py``) on those frames.

Every command runs as a separate process (``--cpus`` pins it with
``taskset``). For each run the script records the wall time, the user and
system CPU time of the process and its children, the peak resident memory of
the largest single process (``os.wait4``) and of the whole process tree
(summed resident memory sampled from ``/proc`` every 0.2 s, Linux; shared
pages count once per process, so with workers it is an upper bound), the CPUs
allowed and their maximum clock, the workers,
the grid spacing, the frames requested and the frames written by the tool,
and the 1-minute load average before and after. CREVICE commands use
``--workers`` (default 1) where they take it and single-threaded BLAS.

Output (``--results``, default ``results/timing_trajectory/<label>/``):
``timing.csv`` (every run), ``summary.csv`` (median per case) and
``summary.md`` (the table). POSIX only.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BENCH, RESULTS, crevice_env, default_work_dir, find_tool, machine_record,  # noqa: E402
                    python_for_crevice, workspace_dir, write_csv)

CASES = ("crevice_profile", "hole", "crevice_region", "crevice_region_spacings", "mdpocket_input", "mdpocket")
DEFAULT_SPACINGS = (0.35, 0.7)
QUIET = ["--skip-analysis-views", "-q"]


# ----------------------------------------------------------------------------- system description
def load_system(path: Path, work_dir: Path, overrides: dict | None = None) -> dict:
    """Read a system file; resolve its paths (relative to the file, ``{work_dir}`` expanded)."""
    path = Path(path)
    data = json.loads(path.read_text()) if path.is_file() else {}
    data.update({k: v for k, v in (overrides or {}).items() if v is not None})
    base = path.resolve().parent if path.is_file() else Path.cwd()

    def resolve(value):
        text = str(value).replace("{work_dir}", str(work_dir))
        p = Path(text)
        return (p if p.is_absolute() else base / p).resolve()

    for key in ("topology", "trajectory", "region"):
        if data.get(key):
            data[key] = resolve(data[key])
    data.setdefault("label", path.stem)
    data.setdefault("selection", "protein")
    if data.get("region"):
        region = json.loads(Path(data["region"]).read_text())
        rbase = Path(data["region"]).resolve().parent
        data["region_data"] = region
        data.setdefault("topology", rbase / region["system"]["topology"]["path"])
        data.setdefault("trajectory", rbase / region["system"]["trajectory"]["path"])
        data["reference_dx"] = rbase / region["reference"]["asset"]["path"]
        data.setdefault("obstacle_selection", region["system"].get("obstacle_selection", "not resname TIP3"))
    data.setdefault("obstacle_selection", "not resname TIP3")
    return data


def region_spacing_args(region: dict, topology: Path, trajectory: Path, reference_dx: Path, spacing: float) -> list[str]:
    """``crevice trajectory --cavities`` arguments equivalent to ``--region`` at another grid spacing.

    Mirrors :func:`crevice.region_definition.region_trajectory_args`: the
    definition's protein and obstacle selections, PBC handling,
    reference-region geometry, probe radius, margin, profile padding, axis and
    grid phase, with ``spacing`` (Å) in place of the definition's.
    """
    g, s = region["geometry"], region["system"]
    return ["trajectory", "--cavities", str(topology), str(trajectory), "--reference-volume-dx", str(reference_dx),
            "--selection", s["protein_selection"], "--obstacle-selection", s["obstacle_selection"], "--pbc", s["pbc"],
            "--geometry-mode", "reference-region", "--spacing", f"{spacing:g}",
            "--probe-radius", str(g["probe_radius_A"]), "--region-margin", str(g["region_margin_A"]),
            "--profile-padding", str(g["profile_padding_A"]),
            "--axis=" + ",".join(map(str, g["axis"])), "--grid-phase=" + ",".join(map(str, g["grid_phase"]))]


def hole_python() -> str | None:
    """Interpreter with mdahole2.

    ``$CREVICE_MDAHOLE2_PYTHON``, else ``<workspace>/tools/mdahole2-venv``
    (``install_comparators.sh``), else this interpreter if it has mdahole2.
    """
    value = os.environ.get("CREVICE_MDAHOLE2_PYTHON")
    if value and Path(value).exists():
        return value
    venv = workspace_dir() / "tools" / "mdahole2-venv" / "bin" / "python"
    if venv.exists():
        return str(venv)
    try:
        import mdahole2  # noqa: F401
        return sys.executable
    except ImportError:
        return None


# ----------------------------------------------------------------------------- cases
def build_cases(system: dict, work: Path, *, frames: int | None, spacings, workers: int, selected) -> list[dict]:
    """The commands to time, as dicts (``case``, ``tool``, ``cmd``, ``cwd``, ``env``, ``spacing_A``, ``skip``)."""
    py = python_for_crevice()
    env = crevice_env()
    env.setdefault("CREVICE_PROGRESS", "0")
    # --max-frames: the CLI refuses more than 500 frames by default.
    stop = ["--stop", str(frames), "--max-frames", str(frames)] if frames else []
    top, traj, sel = system.get("topology"), system.get("trajectory"), system["selection"]
    cases = []

    def add(case, tool, cmd=None, *, cwd=None, env_=None, spacing="", grid="", skip="", out=None):
        if case.split("_spacing")[0] in selected or case in selected or (
                case.startswith("crevice_region_") and "crevice_region_spacings" in selected):
            cases.append({"case": case, "tool": tool, "cmd": cmd, "cwd": cwd, "env": env_, "spacing_A": spacing,
                          "grid": grid, "skip": skip, "out": out})

    prof = system.get("profile")
    if prof:
        axis = ",".join(f"{v:g}" for v in prof["axis"])
        origin = ",".join(f"{v:g}" for v in prof["origin"])
        out = work / "crevice_profile"
        cmd = [py, "-m", "crevice", "trajectory", str(traj), "--topology", str(top), "--selection", sel,
               "--analysis", "profile", f"--axis={axis}", f"--origin={origin}", *stop, "-o", str(out / "profile.json"),
               *QUIET]
        for key, flag in (("search_radius", "--search-radius"), ("enclosure_radius", "--enclosure-radius"),
                          ("section_spacing", "--section-spacing"), ("radii", "--radii")):
            if prof.get(key) is not None:
                cmd += [flag, str(prof[key])]
        add("crevice_profile", "CREVICE trajectory (channel profile)", cmd, env_=env,
            spacing=prof.get("section_spacing", 0.5), grid="section grid", out=out)
        hole, hpy = find_tool("hole"), hole_python()
        rad = None
        if hole is not None:
            rad = next((c for c in (hole.parent.parent / "share/hole2/rad/simple.rad", hole.parent.parent / "rad/simple.rad")
                        if c.exists()), None)
        skip = "" if hole and hpy else "HOLE or an environment with mdahole2 not found (CREVICE_HOLE, CREVICE_MDAHOLE2_PYTHON)"
        out_h = work / "hole"
        cmd = None
        if not skip:
            cmd = [hpy, str(BENCH / "trajectory_worker.py"), "hole", str(top), str(traj), "--select", sel,
                   f"--cpoint={origin}", f"--cvect={axis}", "--hole", str(hole), "--tmp", str(out_h),
                   "--out", str(out_h / "hole_frames.csv")]
            if frames:
                cmd += ["--stop", str(frames)]
            if rad:
                cmd += ["--radii", str(rad)]
        add("hole", "HOLE 2.3.1 (mdahole2, per frame)", cmd, env_=dict(os.environ, OMP_NUM_THREADS="1"), spacing=0.2,
            grid="plane spacing", skip=skip, out=out_h)
    region = system.get("region_data")
    if region:
        w = ["--workers", str(workers)]
        out = work / "crevice_region"
        add("crevice_region", "CREVICE trajectory --region",
            [py, "-m", "crevice", "trajectory", "--region", str(system["region"]), "--out-dir", str(out), *stop, *w,
             "--skip-water-density", "--skip-interactions", *QUIET],
            env_=env, spacing=region["geometry"]["spacing_A"], grid="volume grid", out=out)
        for s in spacings:
            out = work / f"crevice_region_{s:g}"
            add(f"crevice_region_{s:g}", f"CREVICE trajectory --cavities (reference region at {s:g} A)",
                [py, "-m", "crevice", *region_spacing_args(region, top, traj, system["reference_dx"], s),
                 "--out-dir", str(out), "--prefix", "region", *stop, *w, "--skip-water-density", "--skip-interactions",
                 *QUIET], env_=env, spacing=s, grid="volume grid", out=out)
        mdp = work / "mdpocket"
        add("mdpocket_input", "MDpocket input (MDAnalysis fit and write)",
            [py, str(BENCH / "trajectory_worker.py"), "mdpocket-input", str(top), str(traj), "--select", sel,
             "--obstacles", system["obstacle_selection"], "--out-dir", str(mdp), *(["--stop", str(frames)] if frames else [])],
            env_=env, out=mdp)
        tool = find_tool("mdpocket")
        add("mdpocket", "MDpocket 4.1.5 (selected pocket)",
            ["./mdpocket", "--trajectory_file", "traj.xtc", "--trajectory_format", "xtc", "-f", "top.pdb",
             "--selected_pocket", "region.pdb", "-o", "mdpout"],
            cwd=mdp, spacing=1.0, grid="MDpocket grid",
            skip="" if tool else "mdpocket not found (CREVICE_MDPOCKET or install_comparators.sh)", out=mdp)
    return cases


def prepare_mdpocket(system: dict, mdp: Path) -> None:
    """Selected-pocket file (reference region on a 1 Å lattice) and a short-path mdpocket link in ``mdp``."""
    import numpy as np

    sys.path.insert(0, str(BENCH.parent / "src"))
    from crevice.residue_evidence import read_binary_dx

    mdp.mkdir(parents=True, exist_ok=True)
    grid, origin, deltas = read_binary_dx(system["reference_dx"])
    pts = origin + np.argwhere(grid > 0.5) @ deltas
    inside = []
    for p in np.unique(np.round(pts), axis=0):
        k = np.round(np.linalg.solve(deltas.T, p - origin)).astype(int)
        if np.all(k >= 0) and np.all(k < grid.shape) and grid[tuple(k)] > 0.5:
            inside.append(p)
    lines = [f"ATOM  {i % 100000:5d}  C   STP     1    {x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00"
             for i, (x, y, z) in enumerate(inside, 1)]
    (mdp / "region.pdb").write_text("\n".join(lines) + "\nEND\n")
    link = mdp / "mdpocket"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(find_tool("mdpocket"))   # mdpocket aborts when argv[0] or its arguments are long paths


# ----------------------------------------------------------------------------- measurement
def descendants(parents: dict[int, int], root: int) -> set[int]:
    """``root`` and every process below it, from a ``{pid: parent pid}`` map."""
    children: dict[int, list[int]] = {}
    for pid, ppid in parents.items():
        children.setdefault(ppid, []).append(pid)
    found, stack = set(), [root]
    while stack:
        pid = stack.pop()
        if pid not in found:
            found.add(pid)
            stack.extend(children.get(pid, []))
    return found


def _process_table() -> dict[int, int]:
    table = {}
    for entry in os.scandir("/proc"):
        if entry.name.isdigit():
            try:
                stat = Path(entry.path, "stat").read_text()
                table[int(entry.name)] = int(stat[stat.rindex(")") + 2:].split()[1])
            except (OSError, ValueError):
                pass
    return table


def _rss_mb(pid: int) -> float:
    try:
        return int(Path(f"/proc/{pid}/statm").read_text().split()[1]) * os.sysconf("SC_PAGE_SIZE") / 2 ** 20
    except (OSError, ValueError, IndexError):
        return 0.0


class TreeMemorySampler(threading.Thread):
    """Samples the summed resident memory of a process tree (Linux ``/proc``); ``peak_mb`` is the maximum."""

    def __init__(self, pid: int, interval: float = 0.2):
        super().__init__(daemon=True)
        self.pid, self.interval, self.peak_mb = pid, interval, 0.0
        self._stop_event = threading.Event()
        self.available = Path("/proc/self/statm").exists()

    def run(self) -> None:
        while self.available and not self._stop_event.is_set():
            total = sum(_rss_mb(p) for p in descendants(_process_table(), self.pid))
            self.peak_mb = max(self.peak_mb, total)
            self._stop_event.wait(self.interval)

    def stop(self) -> float | None:
        self._stop_event.set()
        self.join()
        return self.peak_mb if self.available else None


def measure(cmd: list[str], *, cwd=None, env=None, log: Path, cpus: list[int] | None) -> dict:
    """Run one command (pinned to ``cpus`` with taskset when given); wall, CPU, peak RSS."""
    full = (["taskset", "-c", ",".join(map(str, cpus))] if cpus else []) + list(cmd)
    with open(log, "w") as out:
        out.write("$ " + " ".join(shlex.quote(str(c)) for c in full) + "\n")
        out.flush()
        start = time.perf_counter()
        proc = subprocess.Popen(full, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT)
        sampler = TreeMemorySampler(proc.pid)
        sampler.start()
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - start
        tree = sampler.stop()
    scale = 1 / 2 ** 20 if sys.platform == "darwin" else 1 / 1024   # ru_maxrss: bytes on macOS, KiB on Linux
    return {"exit_code": os.waitstatus_to_exitcode(status), "wall_s": wall, "cpu_user_s": usage.ru_utime,
            "cpu_sys_s": usage.ru_stime, "peak_rss_mb": usage.ru_maxrss * scale, "peak_tree_rss_mb": tree}


def frames_written(case: str, out: Path | None) -> tuple[int | None, int | None]:
    """Frames in a run's per-frame output and how many of them gave a result.

    A check that every requested frame was analysed. "Gave a result":
    ``profile_status`` resolved (CREVICE profile), a non-empty HOLE profile, a
    measured region volume (``status`` not unresolved/failed), a row of
    MDpocket's descriptors, or a frame written for MDpocket. ``(None, None)`` when the output is not found.
    """
    if out is None or not out.exists():
        return None, None
    if case == "mdpocket_input":
        count = out / "frames_written.txt"
        n = int(count.read_text()) if count.exists() else None
        return n, n
    patterns = {"crevice_profile": "*_frames.csv", "hole": "hole_frames.csv", "mdpocket": "mdpout_descriptors.txt"}
    pattern = "*_cavity_frames.csv" if case.startswith("crevice_region") else patterns.get(case)
    found = sorted(out.glob(pattern)) if pattern else []
    if not found:
        return None, None
    if case == "mdpocket":
        n = max(0, sum(1 for line in found[0].read_text().splitlines() if line.strip()) - 1)
        return n, n
    with open(found[0]) as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    if case == "crevice_profile":
        ok = sum(r.get("profile_status") == "resolved" for r in rows)
    elif case == "hole":
        ok = sum(int(r.get("points") or 0) > 0 for r in rows)
    else:
        ok = sum(bool(r.get("volume_A3")) and not str(r.get("status", "")).startswith(("unresolved", "fail"))
                 for r in rows)
    return len(rows), ok


def summarise(rows: list[dict]) -> list[dict]:
    """Median wall and CPU time, maximum peak memory and frame rate per case."""
    out = []
    for key in dict.fromkeys(r["case"] for r in rows):
        rs = [r for r in rows if r["case"] == key]
        ok = [r for r in rs if r.get("exit_code") == 0]
        first = rs[0]
        if first.get("status") == "skipped":
            out.append({"case": key, "tool": first["tool"], "status": "skipped", "note": first.get("note", "")})
            continue
        use = ok or rs
        wall = statistics.median(r["wall_s"] for r in use)
        cpu = statistics.median(r["cpu_user_s"] + r["cpu_sys_s"] for r in use)
        counts = [r.get("frames_written") for r in use if r.get("frames_written")]
        frames = first.get("frames") or (int(statistics.median(counts)) if counts else None)
        tree = [r["peak_tree_rss_mb"] for r in use if r.get("peak_tree_rss_mb") is not None]
        written = {r.get("frames_written") for r in use}
        resolved = {r.get("frames_with_result") for r in use}
        out.append({"case": key, "tool": first["tool"], "status": "ok" if ok else "failed", "repeats": len(rs),
                    "failed_runs": len(rs) - len(ok), "frames": frames,
                    "frames_written": ";".join(str(w) for w in sorted(written, key=str)),
                    "frames_with_result": ";".join(str(w) for w in sorted(resolved, key=str)),
                    "spacing_A": first.get("spacing_A", ""), "grid": first.get("grid", ""),
                    "cpus": first.get("cpus", ""), "cpu_max_mhz": first.get("cpu_max_mhz", ""),
                    "workers": first.get("workers", ""),
                    "median_wall_s": wall, "median_cpu_s": cpu, "cpu_per_wall": cpu / wall if wall else None,
                    # No frame rate when every repeat failed: nothing was measured per frame.
                    "frames_per_s": frames / wall if ok and frames and wall else None,
                    "s_per_frame": wall / frames if ok and frames else None,
                    "max_peak_rss_mb": max(r["peak_rss_mb"] for r in use),
                    "max_peak_tree_rss_mb": max(tree) if tree else None,
                    "load_1min_range": f"{min(r['load_1min_before'] for r in rs):.1f}-"
                                       f"{max(r['load_1min_after'] for r in rs):.1f}"})
    return out


def markdown_table(summary: list[dict]) -> str:
    """The summary as a Markdown table (wall/CPU seconds, frames per second, memory)."""
    head = ("| Case | Tool | grid (Å) | CPUs / workers | frames (with result) | median wall (s) | median CPU (s) | frames/s "
            "| peak RSS, largest process / tree (MB) |")
    lines = [head, "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]

    def num(v, digits=3):
        if v in (None, ""):
            return "–"
        v = float(v)
        return f"{v:.{digits}g}" if abs(v) < 10 ** digits else f"{v:.0f}"

    for r in summary:
        if r.get("status") == "skipped":
            lines.append(f"| {r['case']} | {r['tool']} | | | | skipped: {r.get('note', '')} | | | |")
            continue
        grid = f"{r['spacing_A']:g} ({r['grid']})" if isinstance(r.get("spacing_A"), (int, float)) else (r.get("grid") or "–")
        fail = f" ({r['failed_runs']} failed)" if r.get("failed_runs") else ""
        lines.append(f"| {r['case']} | {r['tool']}{fail} | {grid} | {r.get('cpus', '')} / {r.get('workers') or '–'} "
                     f"| {r.get('frames') or 'all'} ({r.get('frames_with_result') or '–'}) | {num(r['median_wall_s'])} | {num(r['median_cpu_s'])} "
                     f"| {num(r.get('frames_per_s'))} | {num(r['max_peak_rss_mb'], 4)} / {num(r.get('max_peak_tree_rss_mb'), 4)} |")
    return "\n".join(lines) + "\n"


def cpu_max_mhz(cpus) -> str:
    """Maximum clock (MHz) of each CPU in ``cpus`` from Linux sysfs, e.g. ``4900`` or ``3700;4900``.

    Hybrid processors mix fast and efficient cores, so the pinned CPU's
    maximum clock is recorded with every timing. Empty when unavailable.
    """
    found = set()
    for cpu in cpus or []:
        try:
            khz = Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/cpuinfo_max_freq").read_text().strip()
            found.add(int(khz) // 1000)
        except (OSError, ValueError):
            pass
    return ";".join(str(v) for v in sorted(found))


def portable(text: str, replacements: dict[str, str]) -> str:
    """``text`` with local paths replaced by labels (the result tables are committed).

    Longer paths are replaced first, so a directory inside another keeps its
    own label.
    """
    for path in sorted((p for p in replacements if p), key=len, reverse=True):
        text = text.replace(path, replacements[path])
    return text


def parse_cpus(text: str | None) -> list[int] | None:
    """``0,2-3`` -> ``[0, 2, 3]``; ``None`` or empty -> no pinning."""
    if not text:
        return None
    cpus = []
    for part in text.split(","):
        lo, _, hi = part.partition("-")
        cpus.extend(range(int(lo), int(hi or lo) + 1))
    return sorted(set(cpus))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("system", type=Path, help="system description (JSON); see the module docstring")
    ap.add_argument("--topology", help="override the system file's topology")
    ap.add_argument("--trajectory", help="override the system file's trajectory")
    ap.add_argument("--region", help="override the system file's prepared region definition")
    ap.add_argument("--label", help="override the system label (results subdirectory)")
    ap.add_argument("--frames", type=int, default=0, help="analyse the first N frames (default: all)")
    ap.add_argument("--repeat", type=int, default=3, help="runs per case; the summary reports the median")
    ap.add_argument("--cpus", help="CPU list for taskset, e.g. 0 or 0-1 (default: no pinning)")
    ap.add_argument("--workers", type=int, default=1, help="--workers for the CREVICE commands that take it")
    ap.add_argument("--spacings", default=",".join(f"{s:g}" for s in DEFAULT_SPACINGS),
                    help="extra grid spacings (Å) for the region volume (default %(default)s)")
    ap.add_argument("--cases", default=",".join(CASES), help="comma-separated subset of: " + ", ".join(CASES))
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=None, help="default: results/timing_trajectory/<label>")
    ap.add_argument("--keep-outputs", action="store_true", help="keep the last run's outputs in the work directory")
    ap.add_argument("--dry-run", action="store_true", help="print the commands and exit")
    args = ap.parse_args(argv)

    system = load_system(args.system, args.work_dir, {"topology": args.topology, "trajectory": args.trajectory,
                                                       "region": args.region, "label": args.label})
    label = system["label"]
    work = args.work_dir / "timing_trajectory" / label
    results = args.results or RESULTS / "timing_trajectory" / label
    cpus = parse_cpus(args.cpus)
    spacings = [float(s) for s in args.spacings.split(",") if s]
    selected = set(args.cases.split(","))
    unknown = selected - set(CASES)
    if unknown:
        raise SystemExit(f"unknown cases {sorted(unknown)}; choose from {', '.join(CASES)}")
    cases = build_cases(system, work, frames=args.frames or None, spacings=spacings, workers=args.workers,
                        selected=selected)
    if args.dry_run:
        for c in cases:
            print(f"# {c['case']}: {c['tool']}" + (f" (skipped: {c['skip']})" if c["skip"] else ""))
            if c["cmd"]:
                print(("cd " + shlex.quote(str(c["cwd"])) + " && " if c["cwd"] else "")
                      + " ".join(shlex.quote(str(x)) for x in c["cmd"]))
        return 0

    machine = machine_record()
    allowed = (cpus or sorted(os.sched_getaffinity(0))) if hasattr(os, "sched_getaffinity") else cpus
    places = {str(work.resolve()): "<work>", str(args.work_dir.resolve()): "<work-dir>",
              str(args.system.resolve().parent): "<system-dir>", str(BENCH.parent): "<repo>",
              python_for_crevice(): "python", str(hole_python() or ""): "python-mdahole2",
              str(find_tool("hole") or ""): "hole", str(Path.home()): "~"}
    for key in ("topology", "trajectory", "region"):
        if system.get(key):
            places[str(Path(system[key]).parent)] = f"<{key}-dir>"
    logs = work / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    rows = []
    for n, c in enumerate(cases):
        base = {"case": c["case"], "tool": c["tool"], "frames": args.frames or "", "spacing_A": c["spacing_A"],
                "grid": c["grid"], "cpus": len(allowed) if allowed else "", "cpu_max_mhz": cpu_max_mhz(allowed),
                "workers": args.workers if c["case"].startswith("crevice_region") else 1}
        if c["skip"]:
            print(c["case"], "skipped:", c["skip"])
            rows.append({**base, "status": "skipped", "note": c["skip"]})
            continue
        if c["case"] == "mdpocket" and not (work / "mdpocket" / "traj.xtc").exists():
            print("mdpocket skipped: run the mdpocket_input case first")
            rows.append({**base, "status": "skipped", "note": "no MDpocket input (mdpocket_input case)"})
            continue
        for rep in range(args.repeat):
            out = c["out"]
            if out is not None and c["case"] not in ("mdpocket_input", "mdpocket") and out.exists():
                shutil.rmtree(out)
            if c["case"] == "mdpocket":
                prepare_mdpocket(system, out)
                for old in out.glob("mdpout*"):
                    old.unlink()
            elif c["case"] == "hole":
                out.mkdir(parents=True, exist_ok=True)
            load_before = os.getloadavg()[0]
            res = measure(c["cmd"], cwd=c["cwd"], env=c["env"], log=logs / f"{n:02d}_{c['case']}_{rep}.log", cpus=cpus)
            written, with_result = frames_written(c["case"], out)
            row = {**base, "repeat": rep, "status": "ok" if res["exit_code"] == 0 else "failed", **res,
                   "frames_written": written, "frames_with_result": with_result, "load_1min_before": load_before,
                   "load_1min_after": os.getloadavg()[0],
                   "command": portable(" ".join(shlex.quote(str(x)) for x in c["cmd"]), places)}
            rows.append(row)
            print(c["case"], rep, f"{res['wall_s']:.1f} s", f"exit {res['exit_code']}",
                  f"frames {written} (with result {with_result})", flush=True)
        if not args.keep_outputs and c["out"] is not None and c["case"] not in ("mdpocket_input",) and c["out"].exists():
            if c["case"] == "mdpocket":
                for old in c["out"].glob("mdpout*"):
                    old.unlink()
            else:
                shutil.rmtree(c["out"], ignore_errors=True)
    if not args.keep_outputs and "mdpocket_input" in selected:
        shutil.rmtree(work / "mdpocket", ignore_errors=True)
    summary = summarise(rows)
    note = (f"System: {label}. {system.get('description', '')}\n"
            f"Machine: {machine.get('cpu_model', '?')}, {machine['cpu_count']} logical CPUs, {machine['platform']}, "
            f"Python {machine['python']}; load average at start {machine.get('load_average_1_5_15')}. "
            f"CPUs allowed per run: {','.join(map(str, allowed)) if allowed else 'all'} "
            f"(max MHz {cpu_max_mhz(allowed) or '?'}); frames: {args.frames or 'all'}; "
            f"repeats: {args.repeat}. Timings depend on machine load (see the load columns).")
    write_csv(results / "timing.csv", rows, comment=note)
    write_csv(results / "summary.csv", summary, comment=note)
    (results / "summary.md").write_text("<!-- " + note.replace("\n", " ") + " -->\n\n" + markdown_table(summary))
    print(markdown_table(summary))
    return 1 if any(r.get("status") == "failed" for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())

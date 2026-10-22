"""Shared helpers for the CREVICE benchmark scripts.

Paths
-----
Inputs and comparator outputs are large and are never committed. They live in
a *work directory*, by default ``<workspace>/benchmarks`` where ``<workspace>``
is the first directory in ``$CREVICE_WORKSPACE_DIR`` or ``.crevice`` in the
repository (the same rule as ``scripts/workspace_paths.py``). Every script
takes ``--work-dir`` to override it. Small result tables are written to
``benchmarks/results/`` (``--results``).

Comparators are located, in order, from an environment variable, from
``<workspace>/tools`` and from ``PATH`` (see :func:`find_tool`); none is
vendored. ``install_comparators.sh`` installs them into ``<workspace>/tools``.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent
REPO = BENCH.parent
INPUTS = BENCH / "inputs"
RESULTS = BENCH / "results"


def workspace_dir() -> Path:
    """The first ``$CREVICE_WORKSPACE_DIR`` entry, else ``<repo>/.crevice``."""
    for item in os.environ.get("CREVICE_WORKSPACE_DIR", "").split(os.pathsep):
        if item:
            return Path(item).resolve()
    return REPO / ".crevice"


def default_work_dir() -> Path:
    return workspace_dir() / "benchmarks"


def python_for_crevice() -> str:
    """Interpreter used for ``python -m crevice`` (this one, or ``$CREVICE_PYTHON``)."""
    return os.environ.get("CREVICE_PYTHON", sys.executable)


def crevice_env() -> dict[str, str]:
    """Environment for CREVICE subprocesses: single-threaded BLAS and ``src`` importable."""
    env = dict(os.environ)
    src = str(REPO / "src")
    env["PYTHONPATH"] = src + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    return env


TOOL_CANDIDATES = {
    # name: (environment variable, paths under <workspace>/tools, names on PATH)
    "hole": ("CREVICE_HOLE", ["hole2/bin/hole"], ["hole"]),
    "java": ("CREVICE_JAVA", ["openjdk/bin/java"], ["java"]),
    "caver": ("CREVICE_CAVER_JAR", ["caver/caver_3.0.3/caver/caver.jar", "caver/caver.jar"], []),
    "fpocket": ("CREVICE_FPOCKET", ["fpocket/bin/fpocket"], ["fpocket"]),
    "mdpocket": ("CREVICE_MDPOCKET", ["fpocket/bin/mdpocket"], ["mdpocket"]),
    "kvfinder_python": ("CREVICE_KVFINDER_PYTHON", [], []),
    "povme_python": ("CREVICE_POVME_PYTHON", [], []),
}


def find_tool(name: str) -> Path | None:
    """Locate a comparator executable or file, or return ``None``."""
    env_var, rel, names = TOOL_CANDIDATES[name]
    value = os.environ.get(env_var)
    if value and Path(value).exists():
        return Path(value)
    for base in [workspace_dir() / "tools"]:
        for item in rel:
            if (base / item).exists():
                return base / item
    for item in names:
        found = shutil.which(item)
        if found:
            return Path(found)
    return None


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(line for line in handle if not line.startswith("#"))]


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None, comment: str | None = None) -> Path:
    """Write rows as CSV (fields from the first row unless given); floats keep full precision."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        if comment:
            for line in comment.splitlines():
                handle.write(f"# {line}\n")
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k)) for k in fields})
    return path


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}" if abs(value) >= 1e5 or (value != 0 and abs(value) < 1e-3) else f"{value:.4f}"
    if isinstance(value, (list, tuple)):
        return ";".join(str(fmt(v)) for v in value)
    return value


def run_logged(cmd: list[str], log: Path, *, cwd: Path | None = None, env: dict | None = None,
               timeout: float | None = None) -> tuple[int, float]:
    """Run a command, append the command line and its output to ``log``; return (exit code, wall s)."""
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
        code, out = proc.returncode, proc.stdout + proc.stderr
    except subprocess.TimeoutExpired as exc:
        code, out = -9, f"TIMEOUT after {timeout} s\n{exc.stdout or ''}{exc.stderr or ''}"
    wall = time.perf_counter() - start
    with open(log, "a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(shlex.quote(str(c)) for c in cmd) + "\n")
        handle.write(out if isinstance(out, str) else out.decode(errors="replace"))
        handle.write(f"[exit {code}, {wall:.2f} s]\n")
    return code, wall


def machine_record() -> dict:
    """Machine, Python and load at the time of the call (for timing tables)."""
    record = {"platform": platform.platform(), "machine": platform.machine(), "python": platform.python_version(),
              "cpu_count": os.cpu_count()}
    try:
        record["load_average_1_5_15"] = list(os.getloadavg())
    except OSError:
        pass
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                record["cpu_model"] = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return record


def write_json(path: Path, data) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path

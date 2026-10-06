"""Install a built wheel in a fresh virtual environment and check it.

The wheel is installed with its declared runtime dependencies (NumPy, SciPy,
MDAnalysis, Matplotlib and scikit-image) but without extras, so
the check exercises the minimum supported installation. It verifies metadata
(version single-sourcing and unconditional required dependencies), imports the
installed copy rather than the source tree, runs the console scripts on an
analytic control, compares the SciPy-accelerated and pure-Python spatial
paths, and checks that optional mmCIF intake names its install extra.

Use ``--no-index --find-links DIR`` to install from a local wheelhouse offline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import venv
import zipfile


REQUIRED = ("numpy", "scipy", "mdanalysis", "matplotlib", "scikit-image")

INTAKE_CHECKS = """
import json, pathlib
from crevice import pdb

root = pathlib.Path.cwd()
cache = root / "cache"
cache.mkdir(exist_ok=True)
structure = cache / "1GRM.pdb"
structure.write_bytes((root / "wall.pdb").read_bytes())

record = pdb.cached_record("1GRM", structure, "pdb")
pdb.write_provenance(record)
provenance_verified = pdb.verify_cached_structure(structure).sha256 == record.sha256

structure.write_bytes(structure.read_bytes() + b"ATOM\\n")
try:
    pdb.verify_cached_structure(structure)
    tamper_detected = False
except RuntimeError:
    tamper_detected = True

grammar_ok = (
    pdb.normalize_accession("1grm") == "1GRM"
    and not pdb.is_accession("nope")
    and pdb.structure_url("1grm", fmt="cif") == "https://files.rcsb.org/download/1GRM.cif"
    and pdb.structure_url("1grm", fmt="cif", assembly=2)
        == "https://files.rcsb.org/download/1GRM-assembly2.cif"
)

local = pdb.resolve_structure(structure, cache_dir=cache)
local_file_wins = local.kind == "local_file"
try:
    pdb.resolve_structure("9ZZZ", cache_dir=cache, offline=True)
    offline_refused = False
except RuntimeError:
    offline_refused = True

try:
    import gemmi  # noqa: F401
    gemmi_installed = True
except ImportError:
    gemmi_installed = False
try:
    from crevice.mmcif import load_mmcif
    load_mmcif(structure)
    message = ""
except RuntimeError as exc:
    message = str(exc)

print(json.dumps({
    "provenance_verified": provenance_verified,
    "tamper_detected": tamper_detected,
    "accession_grammar_ok": grammar_ok,
    "local_file_wins": local_file_wins,
    "offline_refused": offline_refused,
    "gemmi_installed": gemmi_installed,
    "gemmi_message": message,
}))
"""

SEGMENT_CHECKS = """
import json
from crevice import Atom
from crevice.spatial import SpatialIndex
from crevice.grid import GridConnectivity
s = SpatialIndex((Atom(1, 'CA', 'ALA', 'A', 1, 0, 0, 0, element='C'),))
uses_scipy = s._tree is not None
g = GridConnectivity(s, (-3, 0, 0), 6, 0.2)
tree = {'crossing': s.segment_clearance((-3, 0, 0), (3, 0, 0)),
        'offset': s.segment_clearance((-3, 2, 0), (3, 2, 0)),
        'blocked': not g.allows((0, 0, 0), (1, 0, 0), 1.3, 1.3)}
# The pure-Python fallback path is retained and must agree with the tree path.
s._tree = None
fallback = {'crossing': s.segment_clearance((-3, 0, 0), (3, 0, 0)),
            'offset': s.segment_clearance((-3, 2, 0), (3, 2, 0))}
print(json.dumps({'uses_scipy': uses_scipy, 'tree': tree, 'fallback': fallback}))
"""

VERSIONS = """
import json, importlib.metadata as m, crevice
print(json.dumps({'source': crevice.__file__, 'version': crevice.__version__,
                  'distribution_version': m.version('crevice'),
                  'dependencies': {n: m.version(n) for n in
                                   ('numpy', 'scipy', 'MDAnalysis', 'matplotlib',
                                    'scikit-image')}}))
"""


def _requirement_name(requirement: str) -> str:
    return re.split(r"[\s;<>=!~\[(]", requirement.strip(), maxsplit=1)[0].lower()


def wheel_metadata(wheel: Path) -> dict:
    """Return the wheel's Version and its unconditional/extra requirements."""
    with zipfile.ZipFile(wheel) as archive:
        name = next(n for n in archive.namelist() if n.endswith(".dist-info/METADATA"))
        text = archive.read(name).decode()
    version = re.search(r"^Version: (.+)$", text, re.M).group(1).strip()
    requires = re.findall(r"^Requires-Dist: (.+)$", text, re.M)
    unconditional = sorted(r for r in requires if "extra ==" not in r)
    return {"version": version, "requires_dist_unconditional": unconditional,
            "requires_dist_count": len(requires)}


def main():
    """Install the wheel in a fresh virtual environment, run the checks, write the report.

    Options: ``wheel`` (path of the built ``.whl``), ``--report`` (JSON report
    path), ``--find-links DIR`` (local wheelhouse for pip; repeatable),
    ``--no-index`` (install only from ``--find-links``) and ``--timeout``
    (seconds per subprocess, default 900).

    Raises
    ------
    AssertionError
        If the wheel's version or required dependencies are wrong, or a check
        in the installed environment fails.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--find-links", action="append", default=[], type=Path,
                        help="local wheelhouse passed to pip (repeatable)")
    parser.add_argument("--no-index", action="store_true",
                        help="install only from --find-links (offline)")
    parser.add_argument("--timeout", type=float, default=900.0,
                        help="seconds allowed for each subprocess (default 900)")
    args = parser.parse_args()
    wheel = args.wheel.resolve(strict=True)

    metadata = wheel_metadata(wheel)
    filename_version = wheel.name.split("-")[1]
    assert metadata["version"] == filename_version, (metadata, wheel.name)
    declared = {_requirement_name(r) for r in metadata["requires_dist_unconditional"]}
    missing_required = sorted(set(REQUIRED) - declared)
    assert not missing_required, f"wheel does not require {missing_required}"

    env = {k: v for k, v in os.environ.items() if k not in {"PYTHONPATH", "PYTHONHOME"}}
    with tempfile.TemporaryDirectory(prefix="crevice-wheel-") as tmp:
        root = Path(tmp)
        venv.EnvBuilder(with_pip=True).create(root / "venv")
        python = root / "venv/bin/python"

        def run(*command):
            return subprocess.run([str(a) for a in command], cwd=root, env=env, text=True,
                                  capture_output=True, timeout=args.timeout, check=True).stdout

        pip = [python, "-m", "pip", "install", "--disable-pip-version-check"]
        if args.no_index:
            pip.append("--no-index")
        for link in args.find_links:
            pip += ["--find-links", link.resolve()]
        run(*pip, wheel)
        installed = json.loads(run(python, "-c", VERSIONS))
        assert str(root / "venv") in installed["source"], "Imported source tree instead of installed wheel"
        assert installed["version"] == installed["distribution_version"] == \
            metadata["version"], installed

        source = root / "wall.pdb"
        lines = []
        for serial, (x, y, z) in enumerate([(5, 0, -2), (-5, 0, -2), (0, 5, -2), (0, -5, -2),
                                            (5, 0, 2), (-5, 0, 2), (0, 5, 2), (0, -5, 2)], 1):
            lines.append(f"ATOM  {serial:5d}  CA  ALA A{serial:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}           C")
        source.write_text("\n".join(lines) + "\nEND\n")
        run(root / "venv/bin/crevice", "--help")
        run(root / "venv/bin/crevice", "profile", source, "--axis", "z", "--samples", "3",
            "--search-radius", "0", "-o", root / "profile.json")
        profile = json.loads((root / "profile.json").read_text())
        assert profile["point_count"] == 3 and abs(profile["min_radius"] - 3.3) < 1e-8

        segments = json.loads(run(python, "-c", SEGMENT_CHECKS))
        for path in ("tree", "fallback"):
            assert abs(segments[path]["crossing"] + 1.7) < 1e-12, segments
            assert abs(segments[path]["offset"] - 0.3) < 1e-12, segments
        assert segments["tree"]["blocked"] and segments["uses_scipy"], segments

        checks_script = root / "intake_checks.py"
        checks_script.write_text(INTAKE_CHECKS)
        intake = json.loads(run(python, checks_script))
        assert intake["provenance_verified"] and intake["tamper_detected"]
        assert intake["accession_grammar_ok"], intake
        assert intake["local_file_wins"] and intake["offline_refused"]
        if not intake["gemmi_installed"]:
            assert "crevice[structures]" in intake["gemmi_message"], intake["gemmi_message"]
        freeze = run(python, "-m", "pip", "freeze", "--disable-pip-version-check").split()

    report = {"status": "passed", "wheel": str(wheel),
              "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
              "builder_python": sys.version.split()[0],
              "wheel_metadata": metadata,
              "runtime_dependencies_installed": True, "optional_extras_installed": False,
              "installed_package": installed, "installed_environment": freeze,
              "checks": ["wheel filename/METADATA/__version__ versions agree",
                         "NumPy, SciPy, MDAnalysis, Matplotlib, scikit-image are unconditional requirements",
                         "isolated import of the installed wheel (not the source tree)",
                         "crevice --help console script",
                         "analytic 3.3 A clearance through crevice profile",
                         "SciPy segment clearance equals pure-Python fallback; blocked edge",
                         "cache provenance record and tamper detection",
                         "accession grammar and download URL construction",
                         "input resolution: local file wins, offline refuses",
                         "optional mmCIF intake names its install extra when Gemmi is absent"],
              "profile_min_radius_A": profile["min_radius"], "segment_checks": segments,
              "intake_checks": intake}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

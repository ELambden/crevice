"""Canonical 1GRM/4PYP publish numbers under the 0.3 compact defaults (real structures; optional).

Runs the reference command (samples 101, gemmi parser, 0.25 Å for 1GRM and
0.5 Å for 4PYP) with the default outputs and checks the reference numbers:
1GRM minimum radius 1.3372879041047168 Å and cast 344.65625 Å³ (0.5.0: per-plane
enclosure probe, so the mouths and samples move and each cast plane is swept
by its own probe; 0.4.0 gave 1.337158628192804 Å and 292.40625 Å³, 0.2 to 0.3.3
1.330188807251824 Å); 4PYP unresolved (no profile files), fallback cast
2439.625 Å³ (unchanged since 0.2).
Skipped when the cached mmCIF files are absent. The axis is assumed, not
curated: these are regression numbers, not validated biological results.
"""
import json
import sys
from pathlib import Path

import pytest

from crevice.cli import main

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from workspace_paths import pdb_cache_dir  # noqa: E402


@pytest.mark.parametrize("name,spacing,expected", [
    ("1GRM", "0.25", {"status": "resolved", "min_radius": 1.3372879041047168, "volume": 344.65625}),
    ("4PYP", "0.5", {"status": "unresolved", "volume": 2439.625}),
])
def test_canonical_publish_numbers_are_unchanged_by_the_compact_defaults(name, spacing, expected, tmp_path):
    source = pdb_cache_dir() / f"{name}.cif"
    if not source.is_file():
        pytest.skip(f"{source} not cached")
    out = tmp_path / name
    assert main(["report", str(source), "--offline", "--model", "0", "--parser", "gemmi", "--axis", "auto",
                 "--samples", "101", "--section-spacing", spacing, "--cast-spacing", spacing,
                 "--dpi", "100", "--out-dir", str(out), "--prefix", name]) == 0
    manifest = json.loads((out / f"{name}_manifest.json").read_text(encoding="utf-8"))
    assert manifest["void_cast"]["total_volume"] == expected["volume"]
    status = manifest["profile_status"]["status"]
    assert status == expected["status"]
    if status == "resolved":
        assert json.loads((out / f"{name}_profile.json").read_text(encoding="utf-8"))["min_radius"] == expected["min_radius"]
    else:
        assert not list(out.glob("*_profile.json")) and not list(out.glob("*_profile.csv"))
    assert not list(out.glob("*hydration*")) and not list(out.glob("*.html"))

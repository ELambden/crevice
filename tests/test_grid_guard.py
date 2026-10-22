"""Grid-spacing guard: a warning recorded in the manifest, never a changed setting (software tests)."""
import json
import math

import pytest

from crevice.grid_guard import RECOMMENDED_SPACING_A, grid_resolution_check


def test_rule_and_recommendation():
    check = grid_resolution_check(0.5, 0.8, clearance_source="rolling probe radius", context="cast")
    assert check["status"] == "warning" and check["threshold_A"] == 1.0
    assert "0.35" in check["message"] and check["recommended_spacing_A"] == RECOMMENDED_SPACING_A
    assert grid_resolution_check(0.35, 0.8, clearance_source="x")["status"] == "ok"
    assert grid_resolution_check(0.25, None, clearance_source="none")["status"] == "not_assessed"
    # The 1GRM profile minimum radius (1.34 A) is resolved at 0.5 A but not at 0.7 A.
    assert grid_resolution_check(0.5, 1.337158628192804, clearance_source="profile")["status"] == "ok"
    assert grid_resolution_check(0.7, 1.337158628192804, clearance_source="profile")["status"] == "warning"


def _pocket_pdb(path):
    lines, serial = [], 0
    for i in range(-3, 4):
        for j in range(-3, 4):
            for k in range(-3, 4):
                if max(abs(i), abs(j), abs(k)) == 3:
                    serial += 1
                    lines.append(f"ATOM  {serial:5d}  CA  ALA A{serial:4d}    {2.0*i:8.3f}{2.0*j:8.3f}{2.0*k:8.3f}  1.00  0.00           C")
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8", newline="\n")
    return path


@pytest.mark.parametrize("spacing", ["0.5", "0.35"])
def test_cast_records_an_unassessed_check_and_never_warns(tmp_path, capsys, spacing):
    """A cast measures no narrowest width, so the guard is recorded as not assessed and stays silent."""
    from crevice.cli import main
    pdb = _pocket_pdb(tmp_path / "box.pdb")
    out = tmp_path / spacing
    assert main(["cast", str(pdb), "--out-dir", str(out), "--spacing", spacing, "--min-component-volume", "1",
                 "--outputs", "none", "--quiet"]) == 0
    manifest = json.loads((out / "box_manifest.json").read_text(encoding="utf-8"))
    check = manifest["grid_resolution"]
    assert check["status"] == "not_assessed" and check["spacing_A"] == float(spacing)
    assert check["narrowest_clearance_A"] is None and "not measured" in check["clearance_source"]
    assert "crevice: warning:" not in capsys.readouterr().err

"""--sensitivity: one-at-a-time reruns over margin, probe-ladder neighbours and grid (software tests)."""
import csv
import json

import pytest

from crevice.sensitivity import (SENSITIVITY_COLUMNS, ladder_neighbours, run_sensitivity, standard_variants,
                                 write_sensitivity)


def test_variants_follow_the_ladder_and_skip_the_baseline_grid():
    assert ladder_neighbours(0.8) == (0.75, 0.9)
    assert ladder_neighbours(3.0)[-1] < 3.0
    assert standard_variants(grid=0.25) == [("grid_spacing_A", 0.5)]
    assert standard_variants(margin=0.5) == [("margin_A", 1.5)]  # no negative margin


def test_run_sensitivity_spread_and_failed_variants(tmp_path):
    def compute(overrides):
        if overrides.get("probe_radius_A") == 0.9:
            raise ValueError("no channel")
        return {"volume_A3": 100.0 + 10 * overrides.get("grid_spacing_A", 0.0)}
    result = run_sensitivity(compute, {"probe_radius_A": 0.8, "grid_spacing_A": 0.35},
                             standard_variants(probe=0.8, grid=0.35), command="test")
    assert [r["variant"] for r in result["runs"]] == ["baseline", "probe_radius_A=0.75", "probe_radius_A=0.9",
                                                       "grid_spacing_A=0.25", "grid_spacing_A=0.5"]
    assert result["runs"][2]["status"].startswith("failed: ValueError")
    spread = result["spread"]["volume_A3"]
    assert spread["minimum"] == 100.0 and spread["maximum"] == 105.0 and spread["runs"] == 4
    assert spread["relative_range_percent"] == pytest.approx(5.0)
    written = write_sensitivity(result, tmp_path / "s.csv", tmp_path / "s.json")
    rows = list(csv.DictReader((tmp_path / "s.csv").open(encoding="utf-8")))
    assert tuple(rows[0]) == SENSITIVITY_COLUMNS and len(rows) == 5
    assert json.loads(written[1].read_text(encoding="utf-8"))["command"] == "test"


def _box(path):
    lines, serial = [], 0
    for i in range(-3, 4):
        for j in range(-3, 4):
            for k in range(-3, 4):
                if max(abs(i), abs(j), abs(k)) == 3:
                    serial += 1
                    lines.append(f"ATOM  {serial:5d}  CA  ALA A{serial:4d}    {2.0*i:8.3f}{2.0*j:8.3f}{2.0*k:8.3f}  1.00  0.00           C")
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8", newline="\n")
    return path


def test_cast_sensitivity_is_opt_in_and_recorded(tmp_path, capsys):
    from crevice.cli import main
    pdb = _box(tmp_path / "box.pdb")
    plain = tmp_path / "plain"
    assert main(["cast", str(pdb), "--out-dir", str(plain), "--min-component-volume", "1", "--outputs", "none", "--quiet"]) == 0
    assert not (plain / "box_sensitivity.csv").exists()
    out = tmp_path / "sens"
    assert main(["cast", str(pdb), "--out-dir", str(out), "--min-component-volume", "1", "--outputs", "none",
                 "--quiet", "--sensitivity"]) == 0
    rows = list(csv.DictReader((out / "box_sensitivity.csv").open(encoding="utf-8")))
    assert {r["variant"] for r in rows} == {"baseline", "probe_radius_A=0.75", "probe_radius_A=0.9", "grid_spacing_A=0.25"}
    baseline = [r for r in rows if r["variant"] == "baseline" and r["quantity"] == "cast_volume_A3"][0]
    manifest = json.loads((out / "box_manifest.json").read_text(encoding="utf-8"))
    assert float(baseline["result"]) == pytest.approx(manifest["void_cast"]["total_volume"])
    assert manifest["files"]["sensitivity_csv"].endswith("box_sensitivity.csv")
    assert "sensitivity cast_volume_A3" in capsys.readouterr().out


def test_cavity_trajectory_sensitivity_reruns_a_frame_subset(tmp_path):
    from crevice.cli import main
    from test_hydration import md_fixture, reference_dx
    source, trajectory = md_fixture(tmp_path)
    dx = reference_dx(tmp_path / "ref.dx")
    out = tmp_path / "out"
    assert main(["cavity-trajectory", str(source), str(trajectory), "--reference-volume-dx", str(dx),
                 "--geometry-mode", "reference-region", "--max-frames", "3", "--out-dir", str(out), "--dpi", "60",
                 "--bootstrap-replicates", "200", "--outputs", "none", "--quiet", "--sensitivity",
                 "--sensitivity-frames", "2"]) == 0
    result = json.loads((out / "crevice_sensitivity.json").read_text(encoding="utf-8"))
    assert result["command"] == "cavity-trajectory" and "2 evenly spaced frames" in result["note"]
    parameters = {run["parameter"] for run in result["runs"]}
    assert {"margin_A", "probe_radius_A", "grid_spacing_A"} <= parameters
    manifest = json.loads((out / "crevice_cavity_manifest.json").read_text(encoding="utf-8"))
    assert manifest["grid_resolution"]["spacing_A"] == 0.5
    assert manifest["files"]["sensitivity_json"].endswith("crevice_sensitivity.json")

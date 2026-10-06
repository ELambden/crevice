"""End-to-end output and batch error regressions."""
import json
from pathlib import Path

import pytest

from crevice.cli import main
from crevice.static_suite import run_static_benchmark_suite
from test_crevice import _pdb_fixture


def test_cli_distribution_artifacts(tmp_path):
    from test_ensemble_cast_network import ring
    from test_core_api import _pdb_atom

    structure = tmp_path / "input.pdb"
    structure.write_text("\n".join(_pdb_atom(a.serial, a.name, a.resname, a.chain_id,
                          a.resid, a.x, a.y, a.z, a.element) for a in ring().atoms) + "\nEND\n")
    assert main(["trajectory", str(structure), str(structure), "-o", str(tmp_path / "trajectory.json"),
                 "--axis", "z", "--samples", "5", "--search-radius", "0", "--align",
                 "--png", str(tmp_path / "band.png"), "--distribution-csv", str(tmp_path / "band.csv")]) == 0
    band = json.loads((tmp_path / "band.json").read_text())
    assert band["frame_count"] == 2
    assert all(r["upper_radius"] == r["lower_radius"] for r in band["rows"])


def test_publish_map_when_cavity_analysis_skipped(tmp_path):
    structure = tmp_path / "input.pdb"
    structure.write_text(_pdb_fixture())
    assert main(["publish", str(structure), "--out-dir", str(tmp_path / "out"),
                 "--axis", "z", "--search-radius", "0", "--samples", "5",
                 "--cast-spacing", "1", "--skip-cavities", "--skip-tunnels", "--skip-network", "--dpi", "60"]) == 0
    manifest = json.loads((tmp_path / "out/input_manifest.json").read_text())
    for key in ("volume_dx", "volume_pml", "volume_tcl", "volume_cxc"):
        assert Path(manifest["files"][key]).is_file()


def test_batch_reports_bad_file_and_continues(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "1GRM.pdb").write_text("invalid structure")
    (cache / "2CHB.pdb").write_text(_pdb_fixture())
    result = run_static_benchmark_suite(cache_dir=cache, output_dir=tmp_path / "out",
                 pdb_ids=("1GRM", "2CHB"), include_network=False, include_tunnels=False,
                 include_cavities=False, cast_spacing=1, profile_samples=5, search_radius=0,
                 probe_radii=(0,), dpi=60)
    assert result["systems"][0]["status"] == "error"
    assert "ValueError" in result["systems"][0]["error"]
    assert result["systems"][1]["status"] != "error"
    manifest = json.loads(Path(result["manifests"]["2CHB"]).read_text())
    for key in ("volume_dx", "volume_pml", "volume_tcl", "volume_cxc"):
        assert Path(manifest["files"][key]).is_file()


def test_suite_rows_carry_curation_status(tmp_path):
    """Every summary row states whether the system was curated, error or not."""
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "1GRM.pdb").write_text("invalid structure")
    (cache / "2CHB.pdb").write_text(_pdb_fixture())
    result = run_static_benchmark_suite(cache_dir=cache, output_dir=tmp_path / "out",
                 pdb_ids=("1GRM", "2CHB", "4PYP"), include_network=False, include_tunnels=False,
                 include_cavities=False, cast_spacing=1, profile_samples=5, search_radius=0,
                 probe_radii=(0,), dpi=60)
    statuses = {row["pdb_id"]: row["curation_status"] for row in result["systems"]}
    assert statuses == {"1GRM": "uncurated", "2CHB": "set_aside", "4PYP": "uncurated"}
    header = Path(result["summary_csv"]).read_text().splitlines()[0]
    assert "curation_status" in header.split(",")


def test_vmd_script_paths_parse_as_tcl(tmp_path):
    tkinter = pytest.importorskip("tkinter")
    from crevice import void_cast, write_volume_viewer_bundle
    from test_ensemble_cast_network import ring
    from crevice import pore_profile
    frame = ring()
    cast = void_cast(frame, mode="channel", spacing=1,
                     profile=pore_profile(frame, samples=13, search_radius=0))
    structure = tmp_path / 'wall $literal [command] "quoted".pdb'
    structure.write_text(_pdb_fixture())
    files = write_volume_viewer_bundle(cast, structure_path=structure, output_dir=tmp_path)
    interpreter = tkinter.Tcl()
    interpreter.eval("proc mol {args} {return 0}")
    for command in ("color", "display", "axes", "material", "molinfo", "graphics", "scale"):
        interpreter.eval(f"proc {command} {{args}} {{}}")
    # Source ensures Tcl's [info script] resolves sibling files as in VMD.
    interpreter.call("source", files["volume_tcl"])
    assert Path(str(interpreter.getvar("crevice_structure"))) == structure


def test_native_vmd_volume_render(tmp_path):
    import os
    import shutil
    import subprocess
    import sys
    from crevice import pore_profile, void_cast, write_volume_viewer_bundle
    from test_ensemble_cast_network import ring

    if shutil.which("vmd") is None:
        pytest.skip("Native VMD executable not installed")
    root = Path(__file__).resolve().parents[1]
    frame = ring()
    cast = void_cast(frame, mode="channel", spacing=1,
                     profile=pore_profile(frame, samples=13, search_radius=0))
    structure = tmp_path / "structure.pdb"
    structure.write_text(_pdb_fixture())
    files = write_volume_viewer_bundle(cast, structure_path=structure, output_dir=tmp_path)
    env = {**os.environ, "PYTHONPATH": str(root / "src")}
    result = subprocess.run([sys.executable, str(root / "scripts/check_vmd.py"), files["volume_tcl"],
                             "--output", str(tmp_path / "render")], env=env,
                            capture_output=True, text=True, timeout=150)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads((tmp_path / "render/vmd_verification.json").read_text())
    assert report["status"] == "passed"
    assert report["foreground_fraction"] > 0.01
    preview = tmp_path / "render/vmd_native.png"
    rotated = preview.read_bytes()
    command = [sys.executable, str(root / "scripts/check_vmd.py"), files["volume_tcl"],
               "--output", str(tmp_path / "render"), "--rotation", "0", "0", "0"]
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=150)
    assert result.returncode == 0, result.stdout + result.stderr
    assert rotated != preview.read_bytes(), "VMD did not apply the requested camera rotation"

    broken_scene = tmp_path / "broken.tcl"
    broken_scene.write_text("error {intentional render-check failure}\n")
    command[2] = str(broken_scene)
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=150)
    assert result.returncode != 0
    report = json.loads((tmp_path / "render/vmd_verification.json").read_text())
    assert report["status"] == "failed", "An earlier image must not mask a failed render"

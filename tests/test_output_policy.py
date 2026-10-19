"""Output policy (crevice.outputs) and typed settings (crevice.config): software checks on synthetic inputs.

* exact file inventories of the compact default, selected (``--outputs figures``,
  ``--outputs pymol``) and ``--all`` modes (``tests/data/output_inventory_0.3.json``);
* ``--legacy-outputs`` reproduces the CREVICE 0.2.0 file sets exactly
  (``tests/data/legacy_inventory_0.2.0.json``) with unchanged numbers;
* no viewer script or HTML unless its file kind is requested (renderer guard);
* CLI and Python give the same numbers; ``--dry-run`` writes nothing;
  settings round-trip through ``--print-config``/``--config``;
* the typed residue interactions run by default with the numbers of the
  hydration workflow, without water terms.

No case is a biological result.
"""
import json
from pathlib import Path

import pytest

from crevice.cli import main
from crevice.outputs import OutputSelectionError, file_kind, resolve_outputs
from inventory_cases import CASES, channel_pdb, inventory, run_case

DATA = Path(__file__).parent / "data"
LEGACY = json.loads((DATA / "legacy_inventory_0.2.0.json").read_text(encoding="utf-8"))["cases"]
EXPECTED = json.loads((DATA / "output_inventory_0.3.json").read_text(encoding="utf-8"))["cases"]
SCRIPT_SUFFIXES = (".pml", ".vmd", ".tcl", ".cxc", ".bild", "_chimerax.py", ".html")


@pytest.mark.parametrize("name", sorted(CASES))
def test_legacy_outputs_reproduce_the_020_file_set_and_numbers(name, tmp_path):
    from crevice.outputs import COMMANDS
    argv0 = CASES[name](tmp_path / "probe", tmp_path / "probe_out")[0][0]
    extra = ["--legacy-outputs"] if argv0 in COMMANDS else []
    with pytest.warns(DeprecationWarning) if extra else _nothing():
        result = run_case(name, tmp_path / "run", extra)
    assert result["exit"] == 0
    assert result["files"] == LEGACY[name]["files"]
    assert result["numbers"] == LEGACY[name]["numbers"]


class _nothing:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.mark.parametrize("mode", sorted(EXPECTED))
def test_exact_inventories_of_default_selected_and_all_modes(mode, tmp_path):
    case, *extra = mode.split(" ")
    result = run_case(case, tmp_path, extra)
    assert result["exit"] == 0
    assert result["files"] == EXPECTED[mode]["files"]
    # The measurements never depend on the output selection.
    legacy = LEGACY[case]["numbers"]
    assert {k: v for k, v in result["numbers"].items() if k in legacy} == \
        {k: v for k, v in legacy.items() if k in result["numbers"]}


def test_renderer_guard_writes_no_scene_or_html_unless_requested(tmp_path):
    for case in ("publish", "cast", "residue-evidence", "hydration-static"):
        files = run_case(case, tmp_path / case, ["--outputs", "figures"])["files"]
        assert not [f for f in files if f.endswith(SCRIPT_SUFFIXES)], case
    files = run_case("publish", tmp_path / "pymol", ["--outputs", "pymol"])["files"]
    assert [f for f in files if f.endswith(SCRIPT_SUFFIXES)] == ["c_volume.pml"]
    assert not [f for f in files if f.endswith(".png")]


def test_default_publish_has_three_viewers_key_figures_and_no_water_files(tmp_path):
    files = run_case("publish", tmp_path, [])["files"]
    assert {"c_volume.pml", "c_volume.vmd", "c_volume.tcl", "c_volume.cxc"} <= set(files)
    assert [f for f in files if f.endswith(".png")] == ["c_profile_radius.png", "c_residue_contacts.png"]
    assert "c_interactions.json" in files
    assert not [f for f in files if "hydration" in f or "water" in f or f.endswith(".html")]
    assert not [f for f in files if f.startswith(("c_cavities", "c_tunnels", "c_network"))]
    manifest = json.loads((tmp_path / "out" / "c_manifest.json").read_text(encoding="utf-8"))
    for path in manifest["files"].values():
        assert Path(path).is_file(), path
    run = manifest["run"]
    assert run["outputs"]["mode"] == "default" and run["outputs"]["optional_analyses_run"] == ["interactions"]
    assert run["settings"]["geometry"]["samples"] == 21 and run["non_default_settings"]["samples"] == "cli"
    assert run["inputs"][0]["role"] == "structure" and len(run["inputs"][0]["sha256"]) == 64


def test_residue_interactions_equal_the_hydration_workflow_without_water_terms(tmp_path):
    from crevice.parser import load_structure
    from crevice.hydration_workflow import static_hydration_bundle
    from crevice.residue_interactions import static_residue_interactions, write_residue_interactions
    from test_hydration import md_fixture
    source, _ = md_fixture(tmp_path)
    frame = load_structure(source)
    residues = ["A:ALA1", "A:ALA2", "A:ALA3"]
    full = static_hydration_bundle(frame, tmp_path / "full", prefix="x", residue_ids=residues, sasa_samples=32, dpi=40)
    alone = write_residue_interactions(static_residue_interactions(frame, residue_ids=residues), tmp_path / "alone", "x")
    a = json.loads(Path(full["interaction_json"]).read_text(encoding="utf-8"))
    b = json.loads(Path(alone["interaction_json"]).read_text(encoding="utf-8"))
    assert [e for e in a["edges"] if e["kind"] != "water_bridge_hbond"] == b["edges"]
    assert b["metadata"]["water_terms"].startswith("not measured")
    strip = {"mean_hydrogen_bonded_waters", "chi1_wet_endpoint_intervals", "chi1_dry_endpoint_intervals",
             "mean_chi1_step_wet_endpoints_degrees", "mean_chi1_step_dry_endpoints_degrees"}
    assert [{k: v for k, v in r.items() if k not in strip} for r in a["residues"]] == \
        [{k: v for k, v in r.items() if k not in strip} for r in b["residues"]]


def test_dry_run_writes_nothing_and_prints_the_plan(tmp_path, capsys):
    out = tmp_path / "out"
    assert main(["publish", str(channel_pdb(tmp_path)), "--offline", "--out-dir", str(out), "--all", "--dry-run"]) == 0
    assert not out.exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["channel.identities.json", "channel.pdb"]
    err = capsys.readouterr().err
    assert "crevice publish plan (--all)" in err and "nothing was read or written" in err


def test_cli_and_python_give_identical_numbers_and_files(tmp_path):
    from crevice.config import run, settings
    from crevice.figures import write_static_publication_bundle
    from crevice.parser import load_structure
    pdb = channel_pdb(tmp_path)
    args = dict(samples=21, search_radius=8, dpi=40, prefix="c")   # ints become floats as on the command line
    assert main(["publish", str(pdb), "--offline", "--samples", "21", "--search-radius", "8", "--dpi", "40",
                 "--prefix", "c", "--out-dir", str(tmp_path / "cli")]) == 0
    record = run(settings("publish", str(pdb), offline=True, out_dir=str(tmp_path / "py"), **args))
    assert record["plan"]["mode"] == "default"
    files = write_static_publication_bundle(load_structure(pdb), structure_path=pdb, output_dir=tmp_path / "api",
                                            prefix="c", samples=21, search_radius=8.0, dpi=40, outputs="default")
    names = {d: inventory(tmp_path / d) for d in ("cli", "py", "api")}
    assert names["cli"] == names["py"]
    # The bundle writer has no CLI wrapper extras (input report, residue interactions).
    assert set(names["api"]) <= set(names["cli"])
    for d in ("cli", "py", "api"):
        assert (tmp_path / d / "c_profile.json").read_text(encoding="utf-8") == (tmp_path / "cli" / "c_profile.json").read_text(encoding="utf-8")
        assert (tmp_path / d / "c_void_cast.json").read_text(encoding="utf-8") == (tmp_path / "cli" / "c_void_cast.json").read_text(encoding="utf-8")
    assert Path(files["volume_pml"]).is_file() and "volume_vmd" in files


def test_compute_writes_nothing_and_export_matches_the_wrapper(tmp_path):
    from crevice.parser import load_structure
    from crevice.static_workflow import compute_static_workflow, export_static_workflow
    pdb = channel_pdb(tmp_path)
    before = sorted(p.name for p in tmp_path.iterdir())
    result = compute_static_workflow(load_structure(pdb), samples=21, search_radius=8)
    assert sorted(p.name for p in tmp_path.iterdir()) == before
    assert result.status == "resolved" and result.cavities is not None and result.cast.total_volume > 0
    export_static_workflow(result, structure_path=pdb, output_dir=tmp_path / "out", prefix="c", dpi=40,
                           outputs="figures")
    assert not [f for f in inventory(tmp_path / "out") if f.endswith(SCRIPT_SUFFIXES)]


def test_settings_round_trip_through_print_config_and_config(tmp_path, capsys):
    pdb = channel_pdb(tmp_path)
    argv = ["publish", str(pdb), "--out-dir", str(tmp_path / "o"), "--samples", "31", "--outputs", "figures",
            "--cast-spacing", "0.75", "--focus", "1", "2", "3", "--origin", "0,0,1"]
    assert main(argv + ["--print-config"]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["schema"] == "crevice.settings/1" and printed["geometry"]["samples"] == 31
    assert printed["outputs"]["outputs"] == "figures" and printed["geometry"]["origin"] == [0.0, 0.0, 1.0]
    config = tmp_path / "run.json"
    config.write_text(json.dumps(printed), encoding="utf-8", newline="\n")
    assert main(["publish", str(pdb), "--config", str(config), "--print-config"]) == 0
    assert json.loads(capsys.readouterr().out) == printed
    # Command-line flags override the file.
    assert main(["publish", str(pdb), "--config", str(config), "--samples", "41", "--print-config"]) == 0
    assert json.loads(capsys.readouterr().out)["geometry"]["samples"] == 41
    from crevice.config import RunSettings
    assert RunSettings.from_dict(printed).to_dict() == printed


def test_settings_are_validated_before_any_input_is_read(tmp_path, capsys):
    missing = str(tmp_path / "absent.pdb")
    assert main(["publish", missing, "--out-dir", str(tmp_path / "o"), "--section-spacing", "-1"]) == 2
    assert "--section-spacing must be positive" in capsys.readouterr().err
    assert main(["residue-evidence", missing, "--volume-dx", str(tmp_path / "no.dx"), "--out-dir", str(tmp_path)]) == 2
    assert "no.dx does not exist" in capsys.readouterr().err
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"command": "publish", "geometry": {"no_such_option": 1}}), encoding="utf-8", newline="\n")
    with pytest.raises(SystemExit):
        main(["publish", missing, "--out-dir", str(tmp_path / "o"), "--config", str(bad)])
    assert not (tmp_path / "o").exists()


def test_output_vocabulary_resolution_rules():
    plan = resolve_outputs("publish", "default,hydration")
    assert {"pymol", "vmd", "chimerax", "key-figures"} <= plan.formats and plan.runs("hydration") and plan.runs("interactions")
    plan = resolve_outputs("publish", "none")
    assert plan.formats == {"data"} and plan.analyses == {"interactions"}
    assert resolve_outputs("publish", "viewers").viewers == ("pymol", "vmd", "chimerax")
    assert not resolve_outputs("publish", skip=["interactions"]).runs("interactions")
    assert resolve_outputs("cavity-trajectory", default_analyses=("hydration", "interactions")).runs("hydration")
    legacy = resolve_outputs("publish", legacy=True, skip=["hydration"])
    assert not legacy.runs("interactions")   # 0.2 ran them only inside hydration
    for bad in [dict(outputs="nonsense"), dict(outputs="figures", all_outputs=True),
                dict(outputs="hydration", skip=["hydration"])]:
        with pytest.raises(OutputSelectionError):
            resolve_outputs("publish", **bad)
    with pytest.raises(OutputSelectionError):
        resolve_outputs("analyze", "cavities")   # analyze selects cavities with --analysis
    assert file_kind("x_end_view.cxc") == ("extra-views", "chimerax")
    assert file_kind("x_analysis_frame_000001.pdb") == ("scenes", None)
    assert file_kind("x_structure.cif") == ("structures", None)


def test_all_records_missing_prerequisites_instead_of_inventing_outputs(tmp_path):
    result = run_case("publish-unresolved", tmp_path, ["--all"])
    assert result["exit"] == 0 and result["numbers"]["profile_status"] == "unresolved"
    assert not [f for f in result["files"] if "_profile." in f or "fingerprint" in f]
    manifest = json.loads((tmp_path / "out" / "p_manifest.json").read_text(encoding="utf-8"))
    assert "no through-channel profile" in manifest["run"]["outputs"]["notes"]["fingerprint"]

"""CREVICE 0.4 command layout: every 0.3 command name is a deprecated alias of its new form.

Software checks on synthetic inputs (no case is a biological result):

* the layout covers all 25 CREVICE 0.3 commands, and each 0.4 form translates
  to exactly the 0.3 parser (identical parsed options and defaults);
* ``crevice --help`` lists only the nine public commands;
* a 0.3 name prints one deprecation line on standard error, a 0.4 form none;
* every 0.3 command and its 0.4 form, run on the same synthetic inputs into the
  same output path, write byte-identical files and the same standard output.
  The only allowed differences are the recorded command line (``run.argv`` in
  manifests) and wall-clock timings, which differ between any two runs.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import shutil
from pathlib import Path

import pytest

from crevice.cli import (COMMAND_LAYOUT, DEPRECATED_COMMANDS, build_parser, canonical_command, command_form, main,
                         translate_argv)
from inventory_cases import HYDRATION_ARGS, PROFILE_ARGS, channel_pdb, md_inputs, pocket_pdb, region_definition

#: The 25 commands of CREVICE 0.3.3.
COMMANDS_03 = ("region-init", "region-prepare", "region-compare", "region-trajectory", "hydration",
               "cavity-trajectory", "cast", "profile", "residues", "cavities", "tunnels", "network",
               "residue-evidence", "features", "trajectory", "publish", "fingerprint", "fetch", "fetch-example",
               "view", "analyze", "benchmark", "static-suite", "paths", "track-cavities")


def new_form(argv: list[str]) -> list[str]:
    """The 0.4 spelling of a 0.3 command line (the value of --example/--region follows the flag)."""
    words = command_form(argv[0]).split()[1:]
    return [*words, *argv[1:]]


def _parsed(argv):
    return vars(build_parser().parse_args(argv))


# --- Layout and translation ---------------------------------------------------------------


def test_layout_covers_every_03_command_once():
    analyses = [mode.analysis for modes in COMMAND_LAYOUT.values() for mode in modes]
    assert sorted(analyses) == sorted(COMMANDS_03) and len(set(analyses)) == 25
    assert list(COMMAND_LAYOUT) == ["fetch", "profile", "cavities", "paths", "residues", "report", "trajectory",
                                    "region", "view"]
    assert set(DEPRECATED_COMMANDS) == set(COMMANDS_03) - set(COMMAND_LAYOUT)


@pytest.mark.parametrize("analysis", COMMANDS_03)
def test_new_form_translates_to_the_03_parser_with_identical_options(analysis):
    positional = {"fetch": [], "fetch-example": ["glut1-excerpt"], "benchmark": [], "static-suite": [],
                  "cavity-trajectory": ["top.gro", "traj.xtc"], "track-cavities": ["top.gro"],
                  "region-compare": ["l.json", "r.json"], "view": ["out"]}.get(analysis, ["input.pdb"])
    required = {"region-init": ["--region-id", "r", "--output", "d.json"], "cast": ["--out-dir", "o"],
                "region-prepare": ["--out-dir", "o"], "region-compare": ["--out-dir", "o"],
                "region-trajectory": ["--out-dir", "o"], "hydration": ["--out-dir", "o"],
                "cavity-trajectory": ["--reference-volume-dx", "r.dx", "--out-dir", "o"],
                "residue-evidence": ["--volume-dx", "v.dx", "--out-dir", "o"], "publish": ["--out-dir", "o"],
                "fingerprint": ["--out-dir", "o"], "analyze": ["--out-dir", "o"], "paths": ["-o", "o", "--seed", "0,0,0"],
                "track-cavities": ["-o", "o"]}.get(analysis, [] if analysis in ("fetch", "fetch-example", "benchmark",
                                                                                "static-suite", "view", "trajectory")
                                                   else ["-o", "x.json"])
    old = [analysis, *positional, *required]
    new = new_form(old)
    assert translate_argv(new).argv == old
    assert translate_argv(old).deprecated == (analysis in DEPRECATED_COMMANDS)
    assert _parsed(new) == _parsed(old)
    assert canonical_command(" ".join(command_form(analysis).split()[1:])) == analysis


def test_mode_flags_never_shadow_an_option_of_the_default_mode():
    sub = next(a for a in build_parser()._actions if isinstance(a, argparse._SubParsersAction))
    shared = []
    for command, modes in COMMAND_LAYOUT.items():
        flags = {m.flag for m in modes if m.flag and m.flag.startswith("--")}
        for mode in modes:
            options = set(sub.choices[mode.analysis]._option_string_actions)
            for flag in flags & options:
                shared.append((mode.analysis, flag))
                assert mode.flag is not None   # never the default mode
    # The one documented overlap: cavity tracking's own --region (a definition to centre on); --track wins.
    assert shared == [("track-cavities", "--region")]
    assert translate_argv(["trajectory", "--region", "d.json", "--track", "t.gro", "-o", "o"]).argv == \
        ["track-cavities", "--region", "d.json", "t.gro", "-o", "o"]


def test_value_modes_and_conflicting_modes():
    assert translate_argv(["fetch", "--example=glut1-excerpt", "--overwrite"]).argv == \
        ["fetch-example", "glut1-excerpt", "--overwrite"]
    assert translate_argv(["fetch", "--example", "--list"]).argv == ["fetch-example", "--list"]
    assert translate_argv(["-q", "trajectory", "--out-dir", "o", "--region", "d.json"]).argv == \
        ["-q", "region-trajectory", "d.json", "--out-dir", "o"]
    assert translate_argv(["region", "compare", "a", "b", "--out-dir", "o"]).argv == \
        ["region-compare", "a", "b", "--out-dir", "o"]
    assert translate_argv(["profile", "x.pdb", "-o", "o.json", "--", "--fingerprint"]).argv[0] == "profile"
    for argv in (["cavities", "--cast", "--tunnels", "x.pdb"], ["residues", "x", "--network", "--features"],
                 ["report", "--suite", "--benchmark"]):
        with pytest.raises(SystemExit):
            build_parser().parse_args(argv)


def test_help_lists_only_the_nine_public_commands(capsys):
    with pytest.raises(SystemExit):
        main(["--help"])
    text = capsys.readouterr().out
    for command in COMMAND_LAYOUT:
        assert f"\n    {command} " in text
    for old in DEPRECATED_COMMANDS:
        assert f"\n    {old} " not in text
    with pytest.raises(SystemExit):
        main(["cavities", "--cast", "--help"])
    text = capsys.readouterr().out
    assert text.startswith("usage: crevice cavities --cast") and "crevice cavities --tunnels" in text
    with pytest.raises(SystemExit):
        main(["region", "--help"])
    assert "init" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(["region", "bogus"])


def test_old_names_print_one_deprecation_line(capsys):
    assert main(["fetch-example", "--list"]) == 0
    old = capsys.readouterr()
    assert old.err.count("\n") == 1 and "'crevice fetch-example' is deprecated; use 'crevice fetch --example'" in old.err
    assert main(["fetch", "--example", "--list"]) == 0
    new = capsys.readouterr()
    assert new.err == "" and new.out == old.out


# --- Byte-identical outputs -----------------------------------------------------------------


def _is_timing(key: str) -> bool:
    return key == "seconds" or key.endswith("_seconds")


def _mask_timings(value):
    if isinstance(value, dict):
        return {k: ("<time>" if _is_timing(k) else _mask_timings(v)) for k, v in value.items()}
    if isinstance(value, list):
        return [_mask_timings(v) for v in value]
    return value


def _snapshot(root: Path) -> dict:
    """Relative path -> bytes, with the recorded command line (``run.argv``) masked in JSON.

    Wall-clock timings (``*_seconds`` keys and ``seconds`` columns, written by
    the benchmark suite) differ between any two runs and are masked too.
    """
    import csv
    files = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        data = path.read_bytes()
        if path.suffix == ".json":
            try:
                record = json.loads(data)
            except ValueError:
                record = None
            if isinstance(record, dict):
                if isinstance(record.get("run"), dict) and "argv" in record["run"]:
                    record["run"]["argv"] = "<command line>"
                data = json.dumps(_mask_timings(record), sort_keys=True).encode()
        elif path.suffix == ".csv" and b"seconds" in data.split(b"\n", 1)[0]:
            rows = list(csv.reader(io.StringIO(data.decode("utf-8"))))
            timing = [i for i, name in enumerate(rows[0]) if _is_timing(name)]
            data = repr([[("<time>" if i in timing and n else c) for i, c in enumerate(r)]
                         for n, r in enumerate(rows)]).encode()
        files[str(path.relative_to(root))] = data
    return files


def _run(argv, work: Path):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue(), err.getvalue()


def _case(name: str, inputs: Path, out: Path) -> list[str]:
    """0.3 command line of case ``name``; only its own inputs are built; outputs go under ``out``."""
    pdb = str(channel_pdb(inputs))
    pocket = str(pocket_pdb(inputs)) if name == "cast" else ""
    md = name in ("residue-evidence", "cavity-trajectory", "hydration", "track-cavities", "region-init")
    source, trajectory, dx = (str(p) for p in md_inputs(inputs / "md")) if md else ("", "", "")
    defn = str(region_definition(inputs / "region")) if name == "region-trajectory" else ""
    landmarks = inputs / "landmarks.json"
    landmarks.write_text("[]", encoding="utf-8")
    cache = inputs / "cache"
    cache.mkdir(exist_ok=True)
    shutil.copy(pdb, cache / "1GRM.pdb")
    raw_definition = ""
    if name == "region-prepare":
        from test_region_definition import definition
        raw_definition = definition(inputs / "raw-region")[0]
    short = ["--max-frames", "3", "--dpi", "40", "--bootstrap-replicates", "200", *HYDRATION_ARGS]
    cases = {
        "fetch": ["fetch", "--list"],
        "fetch-example": ["fetch-example", "--list"],
        "profile": ["profile", pdb, *PROFILE_ARGS, "-o", f"{out}/p.json"],
        "fingerprint": ["fingerprint", pdb, *PROFILE_ARGS, "--dpi", "40", "--out-dir", str(out), "--prefix", "c"],
        "cavities": ["cavities", pdb, "--offline", "-o", f"{out}/cav.json"],
        "cast": ["cast", pocket, "--offline", "--min-component-volume", "30", "--dpi", "40", *HYDRATION_ARGS,
                 "--out-dir", str(out), "--prefix", "p"],
        "tunnels": ["tunnels", pdb, "--offline", "-o", f"{out}/tun.json"],
        "paths": ["paths", pdb, "--offline", "-o", str(out), "--seed", "0,0,0", "--probe-radius", "1.0"],
        "residues": ["residues", pdb, *PROFILE_ARGS, "-o", f"{out}/res.csv"],
        "network": ["network", pdb, *PROFILE_ARGS, "-o", f"{out}/net.json"],
        "residue-evidence": ["residue-evidence", source, "--offline", "--volume-dx", dx, "--dpi", "40",
                             *HYDRATION_ARGS, "--out-dir", str(out)],
        "features": ["features", pdb, *PROFILE_ARGS, "-o", f"{out}/feat.json"],
        "publish": ["publish", pdb, *PROFILE_ARGS, "--cast-spacing", ".5", "--dpi", "40", *HYDRATION_ARGS,
                    "--out-dir", str(out), "--prefix", "c"],
        "analyze": ["analyze", pdb, *PROFILE_ARGS, "--dpi", "40", *HYDRATION_ARGS, "--out-dir", str(out),
                    "--prefix", "c"],
        "static-suite": ["static-suite", "1GRM", "--cache-dir", str(cache), "--format", "pdb", "--out-dir", str(out),
                         "--samples", "11", "--search-radius", "8", "--probe-radius", "0", "--cast-spacing", "1",
                         "--skip-cavities", "--skip-tunnels", "--skip-network", "--dpi", "40"],
        "benchmark": ["benchmark", "1GRM", "--cache-dir", str(cache), "--format", "pdb", "--samples", "11",
                      "--search-radius", "8", "--output", f"{out}/bench.json"],
        "trajectory": ["trajectory", pdb, pdb, "-o", f"{out}/t.json", "--samples", "21", "--search-radius", "8",
                       "--no-align", "--confidence", "0", "--dpi", "40", *HYDRATION_ARGS],
        "cavity-trajectory": ["cavity-trajectory", source, trajectory, "--reference-volume-dx", dx,
                              "--geometry-mode", "reference-region", "--workers", "1", *short, "--out-dir", str(out)],
        "hydration": ["hydration", source, "--trajectory", trajectory, *short, "--out-dir", str(out)],
        "region-trajectory": ["region-trajectory", defn, "--max-frames", "2", "--workers", "1", "--dpi", "40",
                              "--bootstrap-replicates", "200", *HYDRATION_ARGS, "--out-dir", str(out)],
        "track-cavities": ["track-cavities", source, trajectory, "-o", str(out), "--centre", "0,0,0",
                           "--radius", "20", "--workers", "1", "--dpi", "40"],
        "region-init": ["region-init", source, "--trajectory", trajectory, "--region-id", "site",
                        "--landmarks-json", str(landmarks), "--reference-volume-dx", dx,
                        "--output", f"{out}/site_v1.json"],
        "region-prepare": ["region-prepare", str(raw_definition), "--out-dir", str(out)],
        "view": ["view", str(inputs / "published"), "-o", f"{out}/v.html"],
    }
    return cases[name]


def _published(inputs: Path) -> None:
    """A small report bundle for ``view``."""
    target = inputs / "published"
    if not target.exists():
        assert main(["report", str(channel_pdb(inputs)), *PROFILE_ARGS, "--cast-spacing", "1", "--dpi", "40",
                     "--outputs", "scenes", "--out-dir", str(target), "--prefix", "c"]) == 0


@pytest.mark.parametrize("name", COMMANDS_03)
def test_old_command_and_new_form_write_identical_outputs(name, tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    work = tmp_path / "work"
    out = work / "out"
    if name == "view":
        _published(inputs)
    if name == "region-compare":
        _compare_case(inputs, tmp_path)
        return
    old = _case(name, inputs, out)
    results = {}
    for label, argv in (("old", old), ("new", new_form(old))):
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        code, stdout, stderr = _run(argv, work)
        assert code == 0, (label, stderr)
        notes = [line for line in stderr.splitlines() if "is deprecated; use" in line]
        assert len(notes) == (1 if label == "old" and name in DEPRECATED_COMMANDS else 0)
        results[label] = (stdout, [l for l in stderr.splitlines() if l not in notes], _snapshot(out))
    assert results["old"][0] == results["new"][0]
    assert results["old"][1] == results["new"][1]
    old_files, new_files = results["old"][2], results["new"][2]
    assert sorted(old_files) == sorted(new_files)
    different = [path for path in old_files if old_files[path] != new_files[path]]
    assert not different, different


def _compare_case(inputs: Path, tmp_path: Path) -> None:
    """region-compare needs two region-trajectory runs on the same frames."""
    from test_region_definition import definition
    from crevice.region_definition import prepare_region, write_json
    path, data, _ = definition(inputs)
    stats = []
    for region in ("left", "right"):
        data["region_id"] = region
        write_json(path, data)
        prepared = prepare_region(path, inputs / f"{region}-review")
        run = inputs / region
        assert main(["trajectory", "--region", prepared["definition_json"], "--out-dir", str(run), "--skip-hydration",
                     "--stop", "2", "--dpi", "40"]) == 0
        stats.append(str(run / f"{region}_cavity_statistics.json"))
    out = tmp_path / "work" / "out"
    snapshots = {}
    for label, argv in (("old", ["region-compare", *stats, "--out-dir", str(out)]),
                        ("new", ["region", "compare", *stats, "--out-dir", str(out)])):
        if out.exists():
            shutil.rmtree(out)
        code, stdout, _ = _run(argv, tmp_path)
        assert code == 0
        snapshots[label] = (stdout, _snapshot(out))
    assert snapshots["old"] == snapshots["new"]

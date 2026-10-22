"""Recorded paths are relative to the recording file, with ``/`` on every platform (0.4.0).

Software tests on synthetic inputs: an output folder is written, moved, and
read back by every reader; Windows behaviour is exercised with Windows-style
text (``windows=True``, ``PureWindowsPath``) because CI on Windows is separate;
a manifest written by CREVICE 0.3.3 (stored in ``tests/data``) is still read.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path, PureWindowsPath

import pytest

from crevice._portable_paths import (is_absolute_text, portable_json, posix_text, record_path, resolve_files,
                                     resolve_recorded)
from crevice.cli import main
from test_crevice import _pdb_fixture

DATA = Path(__file__).parent / "data" / "portable_paths"
PUBLISH_ARGS = ["--axis", "z", "--search-radius", "0", "--samples", "5", "--cast-spacing", "1",
                "--skip-cavities", "--skip-tunnels", "--skip-network", "--dpi", "60", "--offline"]


# --- the helper ---------------------------------------------------------------

def test_windows_paths_are_recorded_relative_with_forward_slashes():
    base = r"C:\Users\me\run\out"
    assert record_path(r"C:\Users\me\run\out\sub\x_volume.dx", base, windows=True) == "sub/x_volume.dx"
    assert record_path(r"C:\Users\me\run\out", base, windows=True) == "."
    # Outside the folder (or on another drive): absolute, with "/".
    assert record_path(r"C:\Users\me\data\in.pdb", base, windows=True) == "C:/Users/me/data/in.pdb"
    assert record_path(r"D:\md\prod.xtc", base, windows=True) == "D:/md/prod.xtc"
    assert posix_text(r"inputs\prot.pdb", windows=True) == "inputs/prot.pdb"
    assert posix_text("odd\\name.pdb", windows=False) == "odd\\name.pdb"   # a legal POSIX file name


def test_posix_paths_are_recorded_relative():
    assert record_path("/run/out/a/b.json", "/run/out", windows=False) == "a/b.json"
    assert record_path("/run/outside.json", "/run/out", windows=False) == "/run/outside.json"
    assert record_path("/run/out2/x.json", "/run/out", windows=False) == "/run/out2/x.json"   # not a prefix match


def test_portable_json_converts_paths_and_keeps_everything_else():
    base = r"C:\work\out"
    record = {"files": {"volume_dx": r"C:\work\out\p_volume.dx", "scene": r"C:\work\out\views\p_scene.json"},
              "inputs": [{"argument": r"..\data\prot.pdb", "path": r"C:\work\data\prot.pdb", "sha256": "ab" * 32}],
              "argv": ["publish", r"..\data\prot.pdb", "--out-dir", r"C:\work\out"],
              "selection": "protein and not name H*", "units": "oxygen/Å³", "volume_A3": 292.40625,
              "rows": [{"path": r"C:\work\out\a.csv"}, ("x", 1)]}
    out = portable_json(record, base, windows=True)
    assert out["files"] == {"volume_dx": "p_volume.dx", "scene": "views/p_scene.json"}
    assert out["inputs"][0] == {"argument": "../data/prot.pdb", "path": "C:/work/data/prot.pdb", "sha256": "ab" * 32}
    assert out["argv"] == ["publish", "../data/prot.pdb", "--out-dir", "C:/work/out"]   # as typed, "/" only
    assert out["selection"] == record["selection"] and out["units"] == record["units"]
    assert out["volume_A3"] == 292.40625
    assert out["rows"] == [{"path": "a.csv"}, ["x", 1]]
    assert record["files"]["volume_dx"] == r"C:\work\out\p_volume.dx"   # the input is not modified
    assert "\\" not in json.dumps(out, ensure_ascii=False)


def test_relative_output_directory_and_working_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "x.csv").write_text("a\n", encoding="utf-8")
    (tmp_path / "input.pdb").write_text("END\n", encoding="utf-8")
    out = portable_json({"csv": os.path.join("results", "x.csv"), "input": "input.pdb", "word": "results",
                         "prefix": "x", "missing": "results-not-a-file.txt"}, "results")
    assert out["csv"] == "x.csv"
    assert out["input"] == (tmp_path / "input.pdb").as_posix()   # outside the folder: absolute
    assert out["word"] == "results" and out["prefix"] == "x"   # not paths
    assert out["missing"] == "results-not-a-file.txt"


def test_resolve_recorded_reads_new_and_0_3_records(tmp_path):
    (tmp_path / "sub").mkdir()
    target = tmp_path / "sub" / "p_volume.dx"
    target.write_text("x", encoding="utf-8")
    assert resolve_recorded("sub/p_volume.dx", tmp_path) == target
    assert resolve_recorded(str(target), tmp_path) == target
    # 0.3.x: absolute paths of another machine, in either style; subdirectories are kept.
    assert resolve_recorded("/srv/other/run/out/sub/p_volume.dx", tmp_path) == target
    assert resolve_recorded(r"C:\Users\other\out\sub\p_volume.dx", tmp_path) == target
    assert resolve_recorded(r"C:\Users\other\out\p_volume.dx", tmp_path) is not None   # by name, as 0.3.3 did
    assert resolve_recorded("/elsewhere/none.dx", tmp_path, must_exist=True) is None
    assert resolve_recorded(None, tmp_path) is None
    assert is_absolute_text(r"C:\x") and is_absolute_text("C:/x") and is_absolute_text("/x")
    assert is_absolute_text(r"\\server\share\x") and not is_absolute_text("a/b")
    assert PureWindowsPath(r"C:\a\b").as_posix() == "C:/a/b"


# --- whole output folders --------------------------------------------------------

def _texts(root):
    for path in sorted(Path(root).rglob("*")):
        if path.is_file() and path.suffix in {".json", ".csv"}:
            yield path, path.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def publish_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("portable")
    structure = root / "inputs" / "input.pdb"
    structure.parent.mkdir()
    structure.write_text(_pdb_fixture(), encoding="utf-8", newline="\n")
    out = root / "out"
    assert main(["publish", str(structure), "--out-dir", str(out), *PUBLISH_ARGS]) == 0
    return root, structure, out


def test_data_records_hold_no_path_of_the_output_folder(publish_run):
    root, structure, out = publish_run
    manifest = json.loads((out / "input_manifest.json").read_text(encoding="utf-8"))
    assert manifest["files"]["volume_dx"] == "input_volume.dx"
    assert manifest["files"]["structure_file"] == "input_viewer.pdb"
    assert all(not is_absolute_text(v) for v in manifest["files"].values())
    # The input lives outside the folder: absolute (with "/") next to its SHA-256; argv is as typed.
    record = manifest["run"]["inputs"][0]
    assert record["path"] == structure.resolve().as_posix() and len(record["sha256"]) == 64
    for path, text in _texts(out):
        if path.name == "input_manifest.json":
            run = json.loads(text).pop("run")
            assert str(out) in json.dumps(run["argv"])   # the command line as typed
            text = json.dumps({k: v for k, v in json.loads(text).items() if k != "run"})
        assert str(out) not in text, path.name


def test_moved_folder_is_read_by_every_reader(publish_run, tmp_path):
    from crevice.viewer import view
    _, _, out = publish_run
    moved = tmp_path / "copied to another machine"
    shutil.copytree(out, moved)
    manifest = json.loads((moved / "input_manifest.json").read_text(encoding="utf-8"))
    for key, path in resolve_files(manifest["files"], moved).items():
        assert Path(path).is_file() and Path(path).parent == moved, key
    shown = view(moved)
    assert {"protein", "cast"} <= {layer.name for layer in shown.layers}
    assert shown.source == str(moved / "input_manifest.json")
    metadata = json.loads((moved / "input_volume_metadata.json").read_text(encoding="utf-8"))
    assert resolve_recorded(metadata["unsmoothed_reference_dx"], moved) == moved / "input_volume.dx"
    page = (moved / "input_view.html")
    if page.is_file():   # viewer-html is optional in the default output set
        assert str(out) not in page.read_text(encoding="utf-8")


def test_numbers_are_unchanged_only_paths_differ(publish_run):
    """Compared with the 0.3.3 manifest of the same command, only path strings differ."""
    _, _, out = publish_run
    new = json.loads((out / "input_manifest.json").read_text(encoding="utf-8"))
    old = json.loads((DATA / "crevice-0.3.3_publish_manifest.json").read_text(encoding="utf-8"))
    for record in (new, old):
        record.pop("run")
    assert new["void_cast"] == old["void_cast"]
    assert {k: Path(v).name for k, v in new["files"].items()} == {k: Path(v).name for k, v in old["files"].items()}


@pytest.mark.parametrize("style", ["posix", "windows"])
def test_0_3_manifest_with_absolute_paths_still_opens(publish_run, tmp_path, style):
    """A stored 0.3.3 manifest (absolute paths of the writing machine) in a moved folder.

    The stored manifest was written by CREVICE 0.3.3 (``crevice publish`` on this
    fixture); only its folder prefix was replaced by ``/srv/someone/crevice-0.3.3``.
    """
    from crevice.viewer import view
    _, _, out = publish_run
    folder = tmp_path / "old bundle"
    shutil.copytree(out, folder)
    text = (DATA / "crevice-0.3.3_publish_manifest.json").read_text(encoding="utf-8")
    if style == "windows":
        old = json.loads(text)
        old["files"] = {k: "C:\\Users\\someone\\run\\out\\" + Path(v).name for k, v in old["files"].items()}
        text = json.dumps(old)
    (folder / "input_manifest.json").write_text(text, encoding="utf-8")
    manifest = json.loads(text)
    assert all(is_absolute_text(v) for v in manifest["files"].values())
    for key, path in resolve_files(manifest["files"], folder).items():
        assert Path(path).is_file() and Path(path).parent == folder, key
    assert {"protein", "cast"} <= {layer.name for layer in view(folder).layers}


def test_python_run_returns_usable_paths(tmp_path, monkeypatch):
    from crevice.config import run, settings
    structure = tmp_path / "input.pdb"
    structure.write_text(_pdb_fixture(), encoding="utf-8", newline="\n")
    monkeypatch.chdir(tmp_path)
    result = run(settings("publish", "input.pdb", out_dir="rel-out", axis="z", search_radius=0,
                          samples=5, cast_spacing=1, skip_cavities=True, skip_tunnels=True,
                          skip_network=True, dpi=60, offline=True))
    monkeypatch.chdir(tmp_path.parent)
    assert all(Path(p).is_file() for p in result["files"].values())
    saved = json.loads((tmp_path / "rel-out" / "input_manifest.json").read_text(encoding="utf-8"))
    assert saved["files"]["volume_dx"] == "input_volume.dx"
    assert saved["run"]["inputs"][0]["argument"] == "input.pdb"


def test_region_comparison_reads_moved_runs(tmp_path):
    """region-trajectory outputs (reader, statistics, profile) are compared after the folders move."""
    from crevice.region_comparison import compare_regions
    from crevice.region_definition import prepare_region, write_json
    from test_region_definition import definition
    path, data, _ = definition(tmp_path)
    moved = []
    for name in ("left", "right"):
        data["region_id"] = name
        write_json(path, data)
        prepared = prepare_region(path, tmp_path / (name + "-review"))
        out = tmp_path / name
        assert main(["region-trajectory", prepared["definition_json"], "--out-dir", str(out), "--skip-hydration",
                     "--stop", "2", "--dpi", "60"]) == 0
        reader = json.loads((out / "reader.json").read_text(encoding="utf-8")) if (out / "reader.json").exists() \
            else json.loads((out / f"{name}_reader.json").read_text(encoding="utf-8"))
        assert is_absolute_text(reader["topology"])   # the input is outside the folder
        statistics = json.loads((out / f"{name}_cavity_statistics.json").read_text(encoding="utf-8"))
        leaked = [line for line in json.dumps(statistics, indent=0).splitlines() if str(out) + os.sep in line or line.endswith(str(out) + '",')]
        assert not leaked, leaked[:5]
        target = tmp_path / "moved" / name
        target.parent.mkdir(exist_ok=True)
        shutil.move(str(out), target)
        moved.append(target / f"{name}_cavity_statistics.json")
    result = compare_regions(moved[0], moved[1], tmp_path / "comparison")
    assert result["frame_count"] == 2 and result["paired_volume_difference_A3"]["mean"] == 0

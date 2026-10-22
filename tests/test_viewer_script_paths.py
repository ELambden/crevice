"""Generated viewer scripts hold no absolute path of the writing machine (0.4.0).

Software tests on synthetic inputs. Every PyMOL/VMD/ChimeraX script (and
HTML page, BILD file) written by ``crevice report --all`` and by trajectory
runs is scanned for the output folder's absolute path and the home
directory. Scripts name their data files relative to their own folder and
look for them next to the script, then in the working directory, with an
error naming a missing file. The lookup snippets are executed here (Python
with a stand-in for PyMOL's ``cmd``; Tcl through ``tkinter.Tcl``); the native
viewers themselves are run separately (``scripts/check_*.py``).
"""
from __future__ import annotations

import re
import shutil
import types
from pathlib import Path

import pytest

from crevice._script_paths import PYMOL_INPUT, script_relative, tcl_input_proc
from inventory_cases import run_case

SCRIPT_SUFFIXES = {".pml", ".py", ".cxc", ".vmd", ".tcl", ".html", ".bild"}


@pytest.fixture(scope="module", params=["publish", "cavity-trajectory", "region-trajectory", "trajectory"])
def written(request, tmp_path_factory):
    """``report --all`` and three trajectory runs (``--all``) on synthetic inputs."""
    scratch = tmp_path_factory.mktemp(request.param)
    result = run_case(request.param, scratch, ["--all"])
    assert result["exit"] == 0
    return scratch, scratch / "out"


def _forbidden(*folders):
    texts = set()
    for folder in folders:
        folder = Path(folder)
        for form in (folder, folder.resolve()):
            texts |= {str(form), form.as_posix()}
    home = Path.home()
    if len(home.parts) > 1:
        texts |= {str(home), home.as_posix()}
    return texts


def test_scripts_hold_no_absolute_path_of_the_writing_machine(written):
    scratch, out = written
    scripts = [p for p in sorted(out.rglob("*")) if p.is_file() and p.suffix in SCRIPT_SUFFIXES]
    assert {".pml", ".cxc", ".vmd"} <= {p.suffix for p in scripts}
    forbidden = _forbidden(out, scratch)
    leaks = {str(p.relative_to(out)): sorted(f for f in forbidden if f in p.read_text(encoding="utf-8", errors="replace"))
             for p in scripts}
    assert not {k: v for k, v in leaks.items() if v}


_PML_INPUT = re.compile(r"""crevice_input\((['"])(.+?)\1\)""")
_TCL_INPUT = re.compile(r'crevice_input "([^"]+)"')
_CXC_OPEN = re.compile(r'^open ("[^"]+"|\S+)', re.M)


def test_moved_folder_names_only_files_inside_it(written, tmp_path):
    """After the folder is copied elsewhere, every file a script names is next to that script."""
    _, out = written
    moved = tmp_path / "copied to another machine"
    shutil.copytree(out, moved)
    checked = 0
    for script in sorted(moved.rglob("*")):
        text = script.read_text(encoding="utf-8", errors="replace") if script.suffix in SCRIPT_SUFFIXES else ""
        names = []
        if script.suffix == ".pml":
            names = [m[1] for m in _PML_INPUT.findall(text)]
        elif script.suffix in {".vmd", ".tcl"}:
            names = _TCL_INPUT.findall(text)
        elif script.suffix == ".cxc":
            names = [m.strip('"') for m in _CXC_OPEN.findall(text)]
        for name in names:
            assert (script.parent / name).is_file(), (script.name, name)
            checked += 1
    assert checked


def test_profile_pymol_script_names_an_outside_input_relatively(tmp_path):
    """``profile --pymol``: the input structure lies outside the folder and is named by a relative path."""
    from crevice.cli import main
    from inventory_cases import PROFILE_ARGS, channel_pdb
    pdb = channel_pdb(tmp_path / "inputs")
    out = tmp_path / "out"
    out.mkdir()
    assert main(["profile", str(pdb), *PROFILE_ARGS, "-o", str(out / "p.json"), "--pdb", str(out / "p_dummy.pdb"),
                 "--pymol", str(out / "p.pml"), "--pymol-png", str(out / "p.png")]) == 0
    text = (out / "p.pml").read_text(encoding="utf-8")
    assert "crevice_input('../inputs/channel.pdb')" in text and "crevice_input('p_dummy.pdb')" in text
    assert "_crevice_output('p.png')" in text
    assert not any(f in text for f in _forbidden(tmp_path))


def test_script_relative():
    assert script_relative("/run/out/x_viewer.pdb", "/run/out") == "x_viewer.pdb"
    assert script_relative("/run/out/sub/x.dx", "/run/out") == "sub/x.dx"
    assert script_relative("/run/inputs/in.pdb", "/run/out") == "../inputs/in.pdb"
    assert script_relative("/data/in.pdb", "/run/out") == "in.pdb"   # would climb to the root: name only


def _pymol_namespace(script):
    cmd = types.SimpleNamespace(_pymol=types.SimpleNamespace(__script__=str(script)))
    namespace = {"cmd": cmd}
    exec("from pathlib import Path\nimport inspect\n" + PYMOL_INPUT, namespace)
    return namespace


def test_pymol_lookup_script_folder_then_working_directory(tmp_path, monkeypatch):
    folder, elsewhere = tmp_path / "moved bundle", tmp_path / "cwd"
    folder.mkdir()
    elsewhere.mkdir()
    (folder / "x_viewer.pdb").write_text("END\n", encoding="utf-8")
    (elsewhere / "only_here.npz").write_text("x", encoding="utf-8")
    monkeypatch.chdir(elsewhere)
    crevice_input = _pymol_namespace(folder / "x_volume.pml")["crevice_input"]
    assert crevice_input("x_viewer.pdb") == str(folder / "x_viewer.pdb")
    assert crevice_input("only_here.npz") == str(elsewhere / "only_here.npz")
    assert crevice_input("../somewhere/x_viewer.pdb") == str(folder / "x_viewer.pdb")   # bare name next to the script
    with pytest.raises(FileNotFoundError, match="x_scene.json.*next to this script"):
        crevice_input("x_scene.json")


def test_tcl_lookup_script_folder_then_working_directory(tmp_path, monkeypatch):
    tkinter = pytest.importorskip("tkinter")
    try:
        tcl = tkinter.Tcl()
    except tkinter.TclError as error:   # pragma: no cover - Tcl library missing
        pytest.skip(f"Tcl unavailable: {error}")
    folder, elsewhere = tmp_path / "moved bundle", tmp_path / "cwd"
    folder.mkdir()
    elsewhere.mkdir()
    (folder / "x_viewer.pdb").write_text("END\n", encoding="utf-8")
    (elsewhere / "only_here.dx").write_text("x", encoding="utf-8")
    monkeypatch.chdir(elsewhere)
    tcl.eval("cd {" + elsewhere.as_posix() + "}")
    tcl.eval("set crevice_root {" + folder.as_posix() + "}")
    tcl.eval(tcl_input_proc())
    assert Path(tcl.eval('crevice_input "x_viewer.pdb"')) == folder / "x_viewer.pdb"
    assert Path(tcl.eval('crevice_input "only_here.dx"')) == elsewhere / "only_here.dx"
    # Evaluated outside "source" (as under "vmd -e"), the error also says how to run the script.
    with pytest.raises(tkinter.TclError, match="x_volume.dx was not found next to this script.*VMD -e"):
        tcl.eval('crevice_input "x_volume.dx"')

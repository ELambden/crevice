"""How generated viewer scripts find their data files (0.4.0).

Generated PyMOL, VMD and ChimeraX scripts never contain an absolute path of
the machine that wrote them. Each data file is named by its path relative to
the script's own folder (normally just the file name, because the files sit
together in the output folder). When a script runs it looks for the file
next to itself, then in the working directory, and stops with an error naming
the file when neither has it. An output folder therefore keeps working when
it is moved, renamed or copied to another machine.

The snippets here are pasted into the generated text; the functions they
define (``crevice_input`` in Python and Tcl) are part of the scripts, not of
the CREVICE package.
"""
from __future__ import annotations

import os
from pathlib import Path

__all__ = ["script_relative", "pymol_input_lines", "PYMOL_INPUT", "tcl_input_proc"]


def script_relative(path: str | os.PathLike, script_dir: str | os.PathLike) -> str:
    """Path of ``path`` relative to ``script_dir``, with ``/``, for a generated script.

    Files inside the folder give their relative path (``"x_viewer.pdb"``,
    ``"sub/x.dx"``); a file elsewhere gives a ``../`` path, or only its name
    when the two share no folder below the file-system root or drive (a
    relative path climbing to the root would spell out the absolute path).
    The result never holds an absolute path.
    """
    target = Path(os.path.abspath(path))
    base = Path(os.path.abspath(script_dir))
    try:
        common = Path(os.path.commonpath([target, base]))
    except ValueError:   # another drive (Windows)
        return target.name
    if common == Path(common.anchor):
        return target.name
    return Path(os.path.relpath(target, base)).as_posix()


def pymol_input_lines(indent: str = "") -> list[str]:
    """Python lines (inside a PyMOL ``python`` block) defining ``crevice_input(name)``.

    ``crevice_input`` returns the file ``name`` (relative to the script) found
    next to the script, else in the working directory (each tried with the
    relative path, then with the bare file name), and raises
    ``FileNotFoundError`` naming the file otherwise. ``_crevice_script_root``
    is the script's folder.

    Parameters
    ----------
    indent : str, default ""
        Prefix added to every line (for a block nested in other code).
    """
    lines = [
        "_crevice_script_root = Path(getattr(cmd._pymol, \"__script__\",",
        "    inspect.currentframe().f_code.co_filename)).resolve().parent",
        "def crevice_input(name):",
        "    for directory in (_crevice_script_root, Path.cwd()):",
        "        for candidate in dict.fromkeys((name, Path(name).name)):",
        "            path = directory / candidate",
        "            if path.is_file():",
        "                return str(path)",
        "    raise FileNotFoundError(\"CREVICE: data file %r was not found next to this script (%s) or in the \"",
        "                            \"working directory (%s); keep the script with the other files of its \"",
        "                            \"output folder\" % (name, _crevice_script_root, Path.cwd()))",
    ]
    return [indent + line for line in lines]


#: :func:`pymol_input_lines` as one block of text.
PYMOL_INPUT = "\n".join(pymol_input_lines())


def tcl_input_proc(root_variable: str = "crevice_root") -> str:
    """Tcl ``proc crevice_input {name}`` for VMD scripts.

    ``root_variable`` names the global holding the script folder
    (``[file dirname [file normalize [info script]]]``). The procedure
    returns the file next to the script, else in the working directory, and
    raises a Tcl error naming the file otherwise.
    """
    return ("proc crevice_input {name} {\n"
            f"    global {root_variable}\n"
            f"    set crevice_hint [expr {{[info script] eq {{}} ? {{ (VMD -e does not tell a script where it is: start VMD in the folder of the script, or type: source PATH)}} : {{}}}}]\n"
            "    set candidates [list $name]\n"
            "    if {[file tail $name] ne $name} {lappend candidates [file tail $name]}\n"
            f"    foreach directory [list ${root_variable} [pwd]] {{\n"
            "        foreach candidate $candidates {\n"
            "            set path [file join $directory $candidate]\n"
            "            if {[file exists $path]} {return $path}\n"
            "        }\n"
            "    }\n"
            f"    error \"CREVICE: data file $name was not found next to this script (${root_variable}) or in the "
            "working directory ([pwd]); keep the script with the other files of its output folder$crevice_hint\"\n"
            "}\n")

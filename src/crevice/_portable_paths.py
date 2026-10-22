"""How CREVICE records file paths in its outputs and finds them again.

From 0.4.0 every path that an output file records (manifests, run records,
scene JSON, review files, suite reports) is written

- **relative to the directory of the file that records it**, with ``/``
  separators on every platform, when it lies inside that directory (output
  files, copied structures, inputs that were placed in the output folder);
- as an **absolute path with ``/`` separators** when it lies outside (an input
  structure or trajectory elsewhere on disk, a file on another Windows drive).
  These cannot move with the folder; they are kept as provenance, next to
  their SHA-256 where the record has one.

So an output folder can be moved or copied to another machine or operating
system. Readers use :func:`resolve_recorded`, which also reads 0.3.x folders
whose manifests hold absolute paths of the writing machine (in either
separator style): a path that does not exist is looked for inside the folder by
its trailing components and then by file name.

The functions take a ``windows`` flag (default: the running platform) so that
Windows behaviour can be exercised on any system with Windows-style text.
"""
from __future__ import annotations

import ntpath
import os
import posixpath
from pathlib import Path, PurePosixPath, PureWindowsPath

__all__ = ["posix_text", "record_path", "portable_json", "dumps_record", "resolve_recorded", "resolve_files",
           "is_absolute_text"]

_ON_WINDOWS = os.name == "nt"
#: Values recorded as typed (only separators become ``/`` on Windows): the
#: command line, the input arguments of the run record and its settings, which
#: are the options as given (``--config`` accepts them back, relative to the
#: working directory like any option), and embedded copies of a region
#: definition (its relative paths are relative to the definition file).
#: ``parent.key`` entries match a key inside that parent only.
VERBATIM_KEYS = frozenset({"argv", "argument", "run.settings", "run.non_default_settings", "definition"})


def _module(windows):
    return ntpath if (_ON_WINDOWS if windows is None else windows) else posixpath


def posix_text(text, windows=None) -> str:
    """``text`` with ``/`` separators: backslashes become ``/`` on Windows; POSIX text is unchanged.

    On POSIX a backslash is a legal file-name character, so it is kept.

    Parameters
    ----------
    text : str or os.PathLike
        Path text.
    windows : bool, optional
        Treat the text as a Windows path (default: the running platform).

    Returns
    -------
    str
    """
    text = os.fspath(text)
    return text.replace("\\", "/") if (_ON_WINDOWS if windows is None else windows) else text


def is_absolute_text(text) -> bool:
    """True for an absolute path written on any platform (``/x``, ``C:\\x``, ``C:/x``, ``\\\\server\\share``).

    Parameters
    ----------
    text : str or os.PathLike
        Path text, in either separator style.

    Returns
    -------
    bool
    """
    text = os.fspath(text)
    return PurePosixPath(text).is_absolute() or PureWindowsPath(text).is_absolute() or text.startswith("\\\\")


def _inside(path, base, module):
    """``path`` relative to ``base`` (lexically) or None when it is not inside ``base``."""
    try:
        relative = module.relpath(path, base)
    except ValueError:                      # Windows: another drive
        return None
    if relative == module.pardir or relative.startswith(module.pardir + module.sep) or module.isabs(relative):
        return None
    return relative


def record_path(path, base, windows=None) -> str:
    """The text recorded for ``path`` in a file that lives in directory ``base``.

    Parameters
    ----------
    path : str or os.PathLike
        File or directory being recorded (absolute, or relative to the
        working directory).
    base : str or os.PathLike
        Directory of the recording file (for a manifest, its own directory).
    windows : bool, optional
        Treat both as Windows paths (default: the running platform).

    Returns
    -------
    str
        ``/``-separated path relative to ``base`` when ``path`` is inside it
        (``"."`` for ``base`` itself), otherwise the absolute path with ``/``
        separators.
    """
    module = _module(windows)
    native = module is os.path
    text, base_text = os.fspath(path), os.fspath(base)
    absolute = module.normpath(module.abspath(text) if native else text)
    base_abs = module.normpath(module.abspath(base_text) if native else base_text)
    relative = _inside(absolute, base_abs, module)
    if relative is None and native:         # one side may be written through a symbolic link
        relative = _inside(os.path.realpath(absolute), os.path.realpath(base_abs), module)
    if relative is not None:
        return relative.replace("\\", "/") if module is ntpath else relative
    return absolute.replace("\\", "/") if module is ntpath else absolute


def _portable_text(text, base, base_text, windows):
    module = _module(windows)
    native = module is os.path
    if module.isabs(text) or (module is ntpath and is_absolute_text(text) and not text.startswith("/")):
        relative = record_path(text, base, windows)
        return relative if not is_absolute_text(relative) else posix_text(text, windows)
    if base_text not in ("", ".") and not is_absolute_text(base_text):
        for sep in {module.sep, "/"}:
            prefix = base_text.rstrip("/\\") + sep
            if text.startswith(prefix) and len(text) > len(prefix):
                return posix_text(text[len(prefix):], windows)
    if (native and (os.sep in text or "/" in text or "." in text) and os.path.exists(text)
            and not os.path.exists(os.path.join(base_text, text))):    # already relative to base: keep
        recorded = record_path(text, base, windows)
        if recorded != ".":
            return recorded     # inside: relative to base; outside: absolute (a relative text means base-relative)
    return text


def portable_json(value, base, windows=None, verbatim_keys=VERBATIM_KEYS):
    """A copy of a JSON-serialisable ``value`` with its recorded paths made portable.

    Every string that is a path of a file or directory inside ``base`` is
    rewritten relative to ``base`` with ``/`` separators; other absolute paths
    keep their location and get ``/`` separators. A relative string that names
    an existing file outside ``base`` (relative to the working directory)
    becomes absolute, because relative text in a record always means relative
    to the record's directory. Strings that are not paths are unchanged: a
    relative string is rewritten only when it starts with the directory as
    given, or names an existing file (it must contain ``/``, the separator or
    a dot).

    Parameters
    ----------
    value : object
        Data about to be written (dicts, lists, tuples, strings, numbers).
    base : str or os.PathLike
        Directory of the file being written.
    windows : bool, optional
        Treat the text as Windows paths (default: the running platform).
    verbatim_keys : set of str
        Dictionary keys (or ``parent.key`` pairs) whose values are recorded as
        typed, with ``/`` separators on Windows (``VERBATIM_KEYS``).

    Returns
    -------
    object
        The converted copy; ``value`` itself is not modified.
    """
    base_text = os.fspath(base)

    def verbatim(item):
        if isinstance(item, (str, os.PathLike)):
            return posix_text(item, windows)
        if isinstance(item, dict):
            return {key: verbatim(inner) for key, inner in item.items()}
        if isinstance(item, (list, tuple)):
            return [verbatim(inner) for inner in item]
        return item

    def walk(item, name=None):
        if isinstance(item, str):
            return _portable_text(item, base, base_text, windows)
        if isinstance(item, os.PathLike):
            return _portable_text(os.fspath(item), base, base_text, windows)
        if isinstance(item, dict):
            return {key: (verbatim(inner) if key in verbatim_keys or f"{name}.{key}" in verbatim_keys
                          else walk(inner, key)) for key, inner in item.items()}
        if isinstance(item, (list, tuple)):
            return [walk(inner, name) for inner in item]
        return item

    return walk(value)


def resolve_recorded(text, base, must_exist=False):
    """Find the file that an output record names.

    Parameters
    ----------
    text : str or os.PathLike or None
        The recorded path: relative to ``base`` (0.4.0 and later) or absolute
        (0.3.x, or a file outside the output folder), in either separator style.
    base : str or os.PathLike
        Directory of the record (the manifest's directory).
    must_exist : bool
        Return None instead of a best guess when nothing is found.

    Returns
    -------
    Path or None
        The first existing candidate: the path itself (absolute) or ``base``
        joined with it (relative); then, inside ``base``, the trailing
        components of the recorded path (longest first, so a moved 0.3.x
        folder with subdirectories is found) and finally its file name.
        Without a match, the primary candidate, or None when ``must_exist``.
    """
    if text is None or text == "":
        return None
    text, base = os.fspath(text), Path(base)
    windows_style = "\\" in text or PureWindowsPath(text).drive
    parts = [p for p in (PureWindowsPath(text) if windows_style else PurePosixPath(text)).parts
             if p not in ("/", "\\") and not p.endswith((":", ":\\", ":/")) and not p.startswith("\\\\")]
    if is_absolute_text(text):
        primary = Path(text)
        parts = [p for p in parts if p not in ("..", ".")]
        candidates = [primary] + [base.joinpath(*parts[i:]) for i in range(len(parts))]
    else:
        primary = base.joinpath(*parts) if parts else base
        candidates = [primary] + [base.joinpath(*parts[i:]) for i in range(1, len(parts))
                                  if ".." not in parts[i:]]
    for candidate in candidates:
        try:
            if candidate.exists():
                return candidate
        except OSError:                     # e.g. a Windows drive path on POSIX
            continue
    return None if must_exist else primary


def dumps_record(value, target, **kwargs) -> str:
    """``json.dumps`` of ``value`` with its paths made portable for a file written at ``target``.

    Parameters
    ----------
    value : object
        Data to write.
    target : str or os.PathLike
        The JSON file it will be written to (paths become relative to its directory).
    **kwargs
        Passed to :func:`json.dumps`.

    Returns
    -------
    str
    """
    import json
    return json.dumps(portable_json(value, os.path.dirname(os.fspath(target)) or "."), **kwargs)


def resolve_files(files, base):
    """A manifest's ``files`` index with every path resolved against ``base``.

    Parameters
    ----------
    files : dict
        ``{key: recorded path}`` as read from a manifest (other values are kept).
    base : str or os.PathLike
        The manifest's directory.

    Returns
    -------
    dict
        ``{key: str}`` paths usable from any working directory
        (see :func:`resolve_recorded`).
    """
    if not isinstance(files, dict):
        return files
    return {key: str(resolve_recorded(value, base)) if isinstance(value, str) and value else value
            for key, value in files.items()}

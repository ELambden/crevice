"""Locate the local structure cache and viewer toolchains used by the check scripts.

The package default is ``.crevice/pdb`` for downloaded structures (``--cache-dir``).
The native-viewer check scripts additionally look for locally unpacked viewers
under ``<workspace>/tools`` (for example ``tools/bin/chimerax``). Candidate
workspace directories are searched per file, in this order:

1. each directory listed in ``$CREVICE_WORKSPACE_DIR`` (separated by
   ``os.pathsep``), when set;
2. ``<repository>/.crevice``.

The first candidate containing the requested path wins. If none does, the
``<repository>/.crevice`` path is returned so error messages name the default.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_WORKSPACE_DIR = ".crevice"


def workspace_dirs(root: Path = REPO) -> list[Path]:
    """Return candidate workspace directories in lookup order.

    Parameters
    ----------
    root : Path, default the repository root
        Repository whose ``.crevice`` directory is the last candidate.

    Returns
    -------
    list of Path
        The ``$CREVICE_WORKSPACE_DIR`` entries, then ``root/.crevice``.
    """
    configured = os.environ.get("CREVICE_WORKSPACE_DIR", "")
    candidates = [Path(item) for item in configured.split(os.pathsep) if item]
    candidates.append(root / DEFAULT_WORKSPACE_DIR)
    return candidates


def workspace_path(relative: str, root: Path = REPO) -> Path:
    """Return the first existing ``<workspace>/<relative>``, else the default one.

    Parameters
    ----------
    relative : str
        Path inside a workspace directory, for example ``"pdb"``.
    root : Path, default the repository root
        Repository root (see :func:`workspace_dirs`).

    Returns
    -------
    Path
        The first existing candidate, or ``root/.crevice/<relative>``.
    """
    candidates = workspace_dirs(root)
    for base in candidates:
        if (base / relative).exists():
            return base / relative
    return candidates[-1] / relative


def pdb_cache_dir(root: Path = REPO) -> Path:
    """Local structure cache used by the check scripts.

    Parameters
    ----------
    root : Path, default the repository root
        Repository root.

    Returns
    -------
    Path
        ``<workspace>/pdb``.
    """
    return workspace_path("pdb", root)


def local_chimerax(root: Path = REPO) -> Path:
    """Locally unpacked ChimeraX launcher (may not exist).

    Parameters
    ----------
    root : Path, default the repository root
        Repository root.

    Returns
    -------
    Path
        ``<workspace>/tools/bin/chimerax``.
    """
    return workspace_path("tools/bin/chimerax", root)

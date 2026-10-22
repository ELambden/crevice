"""MDAnalysis ``Universe`` and ``AtomGroup`` objects as CREVICE inputs.

Every public function that reads a structure or a trajectory from a path
(:func:`crevice.parser.load_structure`, :func:`crevice.trajectory.load_trajectory`
and the MD workflows built on them) also accepts an
:class:`MDAnalysis.core.universe.Universe` or an
:class:`MDAnalysis.core.groups.AtomGroup` in place of the topology path:

* a ``Universe`` is used as loaded (its topology, its trajectory and any
  transformations attached to it); the function's own ``selection`` is
  applied to it exactly as it would be to the files;
* an ``AtomGroup`` *is* the selection: its atoms (in index order, each atom
  once) replace the ``selection`` argument, which is then ignored. Water,
  obstacle and environment selections of the hydration and cavity workflows
  are still made from the group's whole ``Universe``. An updating group is
  frozen to the atoms it holds when the function is called.

Static loaders read the object's *current* frame; trajectory loaders iterate
its trajectory with the usual ``start``/``stop``/``stride``. The loaders return
the ``Universe`` to the frame it was on. Results are identical to the path route
on the same files (``tests/test_md_native_input.py``). Reader reports record
the file names MDAnalysis holds for the Universe, plus ``input_object``
(``"MDAnalysis.Universe"`` or ``"MDAnalysis.AtomGroup"``); for an AtomGroup
``selection`` is recorded as ``"AtomGroup"`` by the loaders.

The helpers here are internal glue; the public behaviour is documented on
each entry point.
"""

from __future__ import annotations

import sys
from contextlib import contextmanager

__all__ = ["is_md_object", "resolve_md_input", "md_source_label", "md_input_paths", "preserve_frame",
           "accepts_md_objects"]


def _mda():
    """The imported MDAnalysis package, or ``None`` if it was never imported.

    An object can only be a Universe or AtomGroup if MDAnalysis was already
    imported by the caller, so this check never imports it.
    """
    return sys.modules.get("MDAnalysis")


def is_universe(value) -> bool:
    """Whether ``value`` is an MDAnalysis ``Universe``."""
    mda = _mda()
    return mda is not None and isinstance(value, mda.Universe)


def is_atomgroup(value) -> bool:
    """Whether ``value`` is an MDAnalysis ``AtomGroup`` (including updating groups)."""
    mda = _mda()
    return mda is not None and isinstance(value, mda.core.groups.AtomGroup)


def is_md_object(value) -> bool:
    """Whether ``value`` is an MDAnalysis ``Universe`` or ``AtomGroup``.

    Examples
    --------
    >>> from crevice.mdinput import is_md_object
    >>> is_md_object("protein.gro")
    False
    """
    return is_universe(value) or is_atomgroup(value)


def universe_of(value):
    """The ``Universe`` of a Universe or AtomGroup.

    Parameters
    ----------
    value : Universe or AtomGroup
        The MDAnalysis object.

    Returns
    -------
    Universe
    """
    return value.universe if is_atomgroup(value) else value


def atomgroup_selection(group) -> str:
    """An MDAnalysis selection string matching exactly the atoms of ``group``.

    Consecutive atom indices are written as inclusive ``index a:b`` ranges, so
    the string stays short for the usual contiguous protein selections.

    Raises
    ------
    ValueError
        If the group is empty.
    """
    import numpy as np
    indices = np.unique(np.asarray(group.indices, dtype=int))
    if not len(indices):
        raise ValueError("The AtomGroup is empty")
    breaks = np.flatnonzero(np.diff(indices) != 1)
    starts = np.concatenate(([indices[0]], indices[breaks + 1]))
    stops = np.concatenate((indices[breaks], [indices[-1]]))
    return " or ".join(f"index {a}:{b}" for a, b in zip(starts.tolist(), stops.tolist()))


def resolve_md_input(topology, trajectory=None, selection="protein"):
    """Normalise a ``(topology, trajectory, selection)`` triple for the MD readers.

    Parameters
    ----------
    topology : str, Path, Universe or AtomGroup
        A path is returned unchanged. A Universe or AtomGroup is returned as
        its Universe.
    trajectory : str or Path, optional
        Must be ``None`` with a Universe or AtomGroup (the object already
        carries its trajectory).
    selection : str
        Returned unchanged, except for an AtomGroup, whose atoms replace it
        (``atomgroup_selection``).

    Returns
    -------
    tuple
        ``(topology_or_universe, trajectory, selection)``.

    Raises
    ------
    ValueError
        If a trajectory path is given together with an MDAnalysis object.
    """
    if not is_md_object(topology):
        return topology, trajectory, selection
    if trajectory is not None:
        raise ValueError("Pass either an MDAnalysis Universe/AtomGroup or topology/trajectory paths, not both; "
                         "load the trajectory into the Universe instead")
    if is_atomgroup(topology):
        return topology.universe, None, atomgroup_selection(topology)
    return topology, None, selection


def _filename(value):
    name = getattr(value, "filename", None)
    return name if isinstance(name, str) and name else None


def md_source_label(topology, trajectory=None):
    """The source label recorded on frames: the trajectory path, else the topology path.

    For a Universe or AtomGroup this is the trajectory (else topology) file
    name MDAnalysis recorded, so frames read through a Universe carry the same
    ``source`` as frames read from the same paths; ``"MDAnalysis.Universe (in
    memory)"`` when no file name is known.
    """
    if not is_md_object(topology):
        return str(trajectory or topology)
    top, traj = md_input_paths(topology)
    return str(traj or top or "MDAnalysis.Universe (in memory)")


def md_input_paths(topology, trajectory=None):
    """File paths behind an input, for provenance checks.

    Returns
    -------
    tuple
        ``(topology_path, trajectory_path)`` as given for paths; for a
        Universe the topology and trajectory file names MDAnalysis recorded
        (``None`` where none is known, for example an in-memory trajectory).
    """
    if not is_md_object(topology):
        return topology, trajectory
    universe = universe_of(topology)
    top, traj = _filename(universe), _filename(universe.trajectory)
    return top, (traj if traj != top else None)


@contextmanager
def preserve_frame(topology):
    """Return a Universe to its current frame after the block (no-op for paths).

    Parameters
    ----------
    topology : str, Path, Universe or AtomGroup
        The input; only an MDAnalysis object is restored.
    """
    if not is_md_object(topology):
        yield
        return
    trajectory = universe_of(topology).trajectory
    current = trajectory.ts.frame
    try:
        yield
    finally:
        if trajectory.ts.frame != current:
            trajectory[current]


def accepts_md_objects(selection: str = "selection"):
    """Decorator: let ``topology`` be a Universe or AtomGroup (``trajectory`` then omitted).

    The wrapped function receives the Universe as ``topology``,
    ``trajectory=None`` and, for an AtomGroup, an index selection of the
    group's atoms as its ``selection`` parameter (named by ``selection``). The
    Universe is returned to its current frame afterwards. Path arguments pass
    through untouched.
    """
    import functools
    import inspect

    def wrap(func):
        signature = inspect.signature(func)

        @functools.wraps(func)
        def inner(*args, **kwargs):
            bound = signature.bind_partial(*args, **kwargs)
            topology = bound.arguments.get("topology")
            if not is_md_object(topology):
                return func(*args, **kwargs)
            chosen = bound.arguments.get(selection, signature.parameters[selection].default)
            universe, trajectory, chosen = resolve_md_input(topology, bound.arguments.get("trajectory"), chosen)
            bound.arguments.update(topology=universe, trajectory=trajectory)
            bound.arguments[selection] = chosen
            with preserve_frame(universe):
                return func(*bound.args, **bound.kwargs)

        return inner

    return wrap


def stored_source(topology, trajectory):
    """What a reopenable source keeps: the Universe itself, or the two paths as strings.

    Parameters
    ----------
    topology : str, Path, Universe or AtomGroup
        Topology path or MDAnalysis object.
    trajectory : str or Path or None
        Trajectory path (``None`` with an MDAnalysis object).

    Returns
    -------
    tuple
        ``(Universe, None)`` or ``(str(topology), str(trajectory))``.
    """
    return (universe_of(topology), None) if is_md_object(topology) else (str(topology), str(trajectory))

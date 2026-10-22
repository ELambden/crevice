"""Output policy: which analyses run and which files a command writes.

Since CREVICE 0.3 the bundle commands write a compact default set: the
manifest, the core data (JSON, CSV and the measured DX maps that downstream
commands read), a few key figures and one scene set for each native viewer
(PyMOL, VMD and ChimeraX). The water analyses (hydration) are opt-in except
in commands whose purpose they are; the typed residue–residue interactions
still run by default. Everything else is
opt-in with ``--outputs`` (CLI) or ``outputs=`` (Python), ``--all`` writes
every applicable output. (``--legacy-outputs``, which reproduced the CREVICE
0.2 file set during the 0.3 series, was removed in 0.4; use ``--all``.)

The vocabulary has two kinds of token:

*File kinds* decide which kinds of file are written for the analyses that
run. ``data`` (JSON/CSV/NPZ records, measured DX maps, manifests) is always
written.

=================  ==========================================================
``key-figures``    the command's key Matplotlib figures (see
                   :data:`KEY_FIGURES`)
``figures``        every Matplotlib figure (includes ``key-figures``)
``pymol``          PyMOL ``.pml`` scenes
``vmd``            VMD ``.vmd`` and ``.tcl`` scenes
``chimerax``       ChimeraX ``.cxc``, ``_chimerax.py`` and ``.bild`` files
``viewers``        shorthand for ``pymol,vmd,chimerax``
``scenes``         scene data shared by the viewers (viewer PDB, display DX
                   maps, scene JSON, meshes); added automatically with any
                   viewer or ``viewer-html``
``viewer-html``    a self-contained interactive HTML view
                   (:func:`crevice.viewer.view`)
``structures``     copies of the input and dummy-atom PDBs (cast spheres,
                   cavity and tunnel centres)
``extra-views``    alias and variant scenes: ``PREFIX_pore_cast.*`` copies,
                   ``PREFIX_end_view.*`` and ``PREFIX_pore_spheres.pml``
                   (with the matching viewer kind)
=================  ==========================================================

*Analyses* add optional measurements; the command's core analyses always run.

=================  ==========================================================
``hydration``      water analyses: residue water contacts and occupancy,
                   SASA, water bridges, water density and the hydration
                   scenes (implies ``interactions``; ``--skip-water-density``
                   and ``--skip-interactions`` still apply)
``interactions``   typed residue–residue interactions (hydrogen bonds,
                   salt-bridge candidates, aromatic contacts) and residue
                   conformations of the lining/boundary residues; on by
                   default (``--skip-interactions`` turns it off). Water
                   bridges need ``hydration``
``cavities``       separate cavity search (``report``)
``tunnels``        tunnel search (``report``)
``network``        residue/channel contact network (``report``)
``fingerprint``    pore fingerprint of a resolved channel (``report``)
=================  ==========================================================

``default`` stands for the command's default selection and ``none`` for data
only. Listing any file kind replaces the default file kinds; listing any
analysis replaces the default optional analyses; the other half keeps its
default. So ``--outputs figures`` writes data plus every figure (no viewer
scene), ``--outputs hydration`` adds hydration to the default files, and
``--outputs default,vmd`` adds VMD scenes to the defaults.

The defaults are defined in one place, :data:`COMMANDS`. A run's resolved
:class:`OutputPlan` is printed before any expensive work, recorded in the
manifest, and consulted by the writers through :func:`wants` while it is
active (:func:`use_outputs`). With no active plan every writer behaves as in
CREVICE 0.2 (everything is written), so direct calls of the ``write_*``
functions are unchanged unless they are given ``outputs=``.
"""
from __future__ import annotations

import contextlib
import functools
import re
from contextvars import ContextVar
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

#: File kinds, in display order.
FORMATS = ("data", "key-figures", "figures", "pymol", "vmd", "chimerax", "scenes", "viewer-html",
           "structures", "extra-views")
#: Native viewer file kinds.
VIEWERS = ("pymol", "vmd", "chimerax")
#: Optional analyses, in display order.
ANALYSES = ("hydration", "interactions", "cavities", "tunnels", "network", "fingerprint")
#: Shorthand tokens.
ALIASES = {"viewers": VIEWERS}
#: Meta tokens.
META = ("default", "none")

#: Figure file-name endings (after the prefix) that count as key figures.
KEY_FIGURES = ("_profile_radius.png", "_residue_contacts.png", "_residue_evidence.png", "_hydration_residues.png",
               "_cavity_volume_timeseries.png", "_cavity_radius_profile.png", "_cavity_residues.png",
               "_fingerprint.png", "_permeation.png")

#: Descriptions printed by ``--outputs help`` and in the plan.
DESCRIPTIONS = {
    "data": "manifest, JSON/CSV/NPZ records and measured DX maps (always written)",
    "key-figures": "the key figures: " + ", ".join(e.lstrip("_") for e in KEY_FIGURES),
    "figures": "every Matplotlib figure",
    "pymol": "PyMOL scenes (.pml)",
    "vmd": "VMD scenes (.vmd, .tcl)",
    "chimerax": "ChimeraX scenes (.cxc, _chimerax.py, .bild)",
    "viewers": "pymol, vmd and chimerax",
    "scenes": "scene data shared by the viewers (added with any viewer)",
    "viewer-html": "self-contained interactive HTML view (PREFIX_view.html)",
    "structures": "input copy and dummy-atom PDBs",
    "extra-views": "alias and variant scenes (pore_cast.*, end_view.*, pore_spheres.pml)",
    "hydration": "water contacts/occupancy, SASA, water bridges, water density and hydration scenes (implies interactions)",
    "interactions": "typed residue-residue interactions (H-bonds, salt bridges, aromatic) and conformations",
    "cavities": "separate cavity search",
    "tunnels": "tunnel search",
    "network": "residue/channel contact network",
    "fingerprint": "pore fingerprint of a resolved channel",
    "default": "the command's default selection",
    "none": "data only",
}


@dataclass(frozen=True)
class CommandOutputs:
    """Output policy of one command: what it can write and what it writes by default.

    Attributes
    ----------
    analyses : tuple of str
        Optional analyses the command accepts (subset of :data:`ANALYSES`).
    default_analyses : tuple of str
        Optional analyses run by default.
    formats : tuple of str
        File kinds the command can write.
    default_formats : tuple of str
        File kinds written by default (``data`` is implied).
    core : str
        Human-readable list of the analyses that always run.
    """

    analyses: tuple = ()
    default_analyses: tuple = ()
    formats: tuple = FORMATS
    default_formats: tuple = ("key-figures", "pymol", "vmd", "chimerax")
    core: str = ""


_WATER = ("hydration", "interactions")

#: The single place where the per-command defaults are defined.
COMMANDS = {
    "publish": CommandOutputs(analyses=("hydration", "interactions", "cavities", "tunnels", "network", "fingerprint"),
                              default_analyses=("interactions",),
                              core="profile, cast (with ENTRY/LUMEN/EXIT segments), residue contacts"),
    "cast": CommandOutputs(analyses=_WATER, default_analyses=("interactions",), core="rolling-probe cast"),
    "analyze": CommandOutputs(analyses=_WATER, default_analyses=("interactions",), default_formats=(),
                              core="the --analysis selection (default profile, residues)"),
    "residue-evidence": CommandOutputs(analyses=_WATER, default_analyses=("interactions",),
                                       core="boundary residue evidence"),
    "hydration": CommandOutputs(analyses=(), core="hydration (with water density and typed interactions)"),
    "cavity-trajectory": CommandOutputs(analyses=_WATER, default_analyses=("interactions",),
                                        core="per-frame cavity geometry and statistics"),
    "region-trajectory": CommandOutputs(analyses=_WATER, default_analyses=("interactions",),
                                        core="per-frame region geometry and statistics"),
    "trajectory": CommandOutputs(analyses=_WATER, default_analyses=("interactions",), default_formats=(),
                                 formats=("data", "key-figures", "figures", "pymol", "vmd", "chimerax", "scenes"),
                                 core="the --analysis selection (default profile, features); explicit paths "
                                      "(--png, --distribution-csv, --network-json) as given"),
}


class OutputSelectionError(ValueError):
    """An invalid or contradictory ``--outputs``/``--all`` selection."""


@dataclass(frozen=True)
class OutputPlan:
    """The resolved output selection of one command run.

    Attributes
    ----------
    command : str
        Command name (a key of :data:`COMMANDS`).
    mode : {"default", "selected", "all"}
        How the plan was chosen.
    formats : frozenset of str
        File kinds written (``data`` always; ``scenes`` added with any viewer).
    analyses : frozenset of str
        Optional analyses that run.
    requested : tuple of str
        The ``--outputs`` tokens as given (empty for the default).
    skipped : tuple of str
        Optional analyses removed by ``--skip-*`` flags.
    """

    command: str
    mode: str
    formats: frozenset
    analyses: frozenset
    requested: tuple = ()
    skipped: tuple = ()

    @property
    def viewers(self) -> tuple:
        """Native viewers whose scenes are written, in :data:`VIEWERS` order."""
        return tuple(v for v in VIEWERS if v in self.formats)

    def runs(self, analysis: str) -> bool:
        """Whether optional ``analysis`` runs in this plan.

        Parameters
        ----------
        analysis : str
            An analysis token (:data:`ANALYSES`).

        Returns
        -------
        bool
        """
        return analysis in self.analyses

    def writes(self, kind: str) -> bool:
        """Whether files of ``kind`` are written (``key-figures`` is implied by ``figures``).

        Parameters
        ----------
        kind : str
            A file-kind token (:data:`FORMATS`).

        Returns
        -------
        bool
        """
        if kind == "key-figures":
            return bool({"key-figures", "figures"} & self.formats)
        return kind in self.formats

    def wants(self, path: str | Path) -> bool:
        """Whether this plan writes the file ``path`` (classified by :func:`file_kind`).

        Parameters
        ----------
        path : str or Path
            Output file path; only the file name is used.

        Returns
        -------
        bool
        """
        kind, viewer = file_kind(path)
        if kind == "extra-views":
            return "extra-views" in self.formats and (viewer is None or viewer in self.formats)
        if kind == "key-figures":
            return self.writes("key-figures")
        if kind == "scenes":
            return "scenes" in self.formats
        if kind == "structures" and Path(path).name.endswith("_pore_cast.pdb") and "extra-views" in self.formats:
            return True   # the sphere scene (pore_spheres.pml) loads these dummy atoms
        return self.writes(kind)

    def to_dict(self) -> dict:
        """JSON record of the plan (written to the manifest as ``outputs``).

        Returns
        -------
        dict
        """
        policy = COMMANDS[self.command]
        return {"mode": self.mode, "requested": list(self.requested),
                "file_kinds": [f for f in FORMATS if f in self.formats],
                "optional_analyses_run": [a for a in ANALYSES if a in self.analyses],
                "optional_analyses_not_run": [a for a in policy.analyses if a not in self.analyses],
                "skipped_by_flag": list(self.skipped), "core_analyses": policy.core,
                "definition": "CREVICE 0.3 output policy (crevice.outputs); data and the manifest are always written"}

    def describe(self) -> str:
        """Multi-line, human-readable plan (printed before expensive work).

        Returns
        -------
        str
        """
        policy = COMMANDS[self.command]
        label = {"default": "default outputs", "selected": "--outputs " + ",".join(self.requested),
                 "all": "--all"}[self.mode]
        kinds = [f for f in FORMATS if f in self.formats and f in policy.formats]
        not_kinds = [f for f in policy.formats if f not in kinds and not (f == "key-figures" and "figures" in kinds)]
        from .cli import command_form
        lines = [f"{command_form(self.command)} plan ({label}):", f"  core analyses: {policy.core}"]
        if policy.analyses:
            run = [a for a in policy.analyses if a in self.analyses]
            off = [a for a in policy.analyses if a not in self.analyses]
            lines.append("  optional analyses: " + (", ".join(run) if run else "none")
                         + (f"; not run: {', '.join(off)}" if off else ""))
        lines.append("  file kinds: " + ", ".join(kinds) + (f"; not written: {', '.join(not_kinds)}" if not_kinds else ""))
        if self.mode == "default":
            lines.append("  add files or analyses with --outputs (e.g. --outputs default,hydration), "
                         "everything with --all")
        return "\n".join(lines)


def file_kind(path: str | Path) -> tuple:
    """Classify an output file name into a file kind.

    Parameters
    ----------
    path : str or Path
        Output path; only the name is used.

    Returns
    -------
    tuple of (str, str or None)
        The file kind (one of :data:`FORMATS`) and, for viewer scripts, the
        viewer it belongs to (``"pymol"``, ``"vmd"`` or ``"chimerax"``).

    Examples
    --------
    >>> file_kind("out/1GRM_volume.pml")
    ('pymol', 'pymol')
    >>> file_kind("1GRM_pore_cast.vmd")
    ('extra-views', 'vmd')
    >>> file_kind("1GRM_profile_radius.png"), file_kind("1GRM_void_cast.csv")
    (('key-figures', None), ('data', None))
    """
    name = Path(path).name
    viewer = _viewer_of(name)
    if _EXTRA_VIEWS.search(name):
        return "extra-views", viewer
    if name.endswith(".html"):
        return "viewer-html", None
    if viewer is not None:
        return viewer, viewer
    if name.endswith(".png"):
        return ("key-figures" if name.endswith(KEY_FIGURES) else "figures"), None
    if _STRUCTURES.search(name):
        return "structures", None
    if _SCENES.search(name):
        return "scenes", None
    return "data", None


def _viewer_of(name: str):
    if name.endswith(".pml"):
        return "pymol"
    if name.endswith((".vmd", ".tcl")):
        return "vmd"
    if name.endswith((".cxc", "_chimerax.py", ".bild")):
        return "chimerax"
    return None


_EXTRA_VIEWS = re.compile(r"_(pore_cast\.(pml|vmd|tcl|cxc)|pore_spheres\.pml|end_view\.(pml|vmd|tcl|cxc))$")
_STRUCTURES = re.compile(r"_(structure(\.[A-Za-z0-9]+)+|pore_cast\.pdb|cavities\.pdb|tunnels\.pdb)$")
_SCENES = re.compile(
    r"(_scene\.json|_volume_metadata\.json|_display\.dx|_display_local\.dx|_pymol_local\.dx|_pymol_mesh\.npz"
    r"|_geometry_reference\.npz|_viewer\.pdb|_viewer\.identities\.json|_analysis_protein\.pdb"
    r"|_analysis_protein\.identities\.json|_analysis_frame_\d+\.pdb|_analysis_frame_\d+\.identities\.json"
    r"|_analysis_member_waters_\d+\.pdb|_analysis_geometry\.npz|_analysis_reference\.tsv|_analysis_cavity\.dx"
    r"|_analysis_cavity_local\.dx|_analysis_evidence\.csv|_residue_reference\.npz|_residue_reference\.tsv"
    r"|_region_\d+\.dx|_region_\d+_local\.dx|_display_regions\.json|_region_annotations\.json)$")


def _tokens(outputs) -> list:
    if outputs is None:
        return []
    if isinstance(outputs, str):
        outputs = outputs.split(",")
    tokens = []
    for item in outputs:
        tokens += [t.strip().lower() for t in str(item).split(",") if t.strip()]
    return tokens


def resolve_outputs(command: str, outputs=None, *, all_outputs: bool = False,
                    skip: Iterable[str] = (), default_analyses: Iterable[str] | None = None) -> OutputPlan:
    """Resolve an output selection into an :class:`OutputPlan`.

    Parameters
    ----------
    command : str
        Command name, a key of :data:`COMMANDS` (``"publish"``, ``"cast"``,
        ``"analyze"``, ``"residue-evidence"``, ``"hydration"``,
        ``"cavity-trajectory"``, ``"region-trajectory"``, ``"trajectory"``).
    outputs : str or sequence of str, optional
        ``--outputs`` tokens (comma-separated string or list). ``None`` or
        empty: the command's defaults.
    all_outputs : bool, default False
        ``--all``: every file kind and optional analysis of the command.
    skip : iterable of str, default ()
        Optional analyses removed by ``--skip-*`` flags (for example
        ``"hydration"``); naming one of them in ``outputs`` is an error.
    default_analyses : iterable of str, optional
        Overrides the command's default optional analyses (used for
        ``--water-membership``, which turns hydration on).

    Returns
    -------
    OutputPlan

    Raises
    ------
    OutputSelectionError
        Unknown tokens, a token the command does not support, contradictory
        flags, or an analysis both requested and skipped.

    Examples
    --------
    >>> plan = resolve_outputs("publish")
    >>> sorted(plan.formats), sorted(plan.analyses)
    (['chimerax', 'data', 'key-figures', 'pymol', 'scenes', 'vmd'], ['interactions'])
    >>> plan = resolve_outputs("publish", "figures,hydration")
    >>> sorted(plan.formats), sorted(plan.analyses)
    (['data', 'figures'], ['hydration', 'interactions'])
    """
    if command not in COMMANDS:
        raise OutputSelectionError(f"{command} has no output policy")
    policy = COMMANDS[command]
    tokens = _tokens(outputs)
    skip = tuple(dict.fromkeys(s for s in skip if s))
    if tokens and all_outputs:
        raise OutputSelectionError("choose one of --outputs and --all")
    if "legacy" in tokens:
        raise OutputSelectionError("the CREVICE 0.2 file set ('legacy', --legacy-outputs) was removed in 0.4; "
                                   "use --all (Python: outputs='all') or list the outputs you need")
    defaults = tuple(policy.default_analyses if default_analyses is None else default_analyses)
    if all_outputs:
        formats, analyses, mode = set(policy.formats), set(policy.analyses), "all"
    else:
        known = set(FORMATS) | set(ALIASES) | set(ANALYSES) | set(META)
        unknown = [t for t in tokens if t not in known]
        if unknown:
            raise OutputSelectionError(f"unknown --outputs value(s): {', '.join(unknown)}; choose from "
                                       + ", ".join([*FORMATS, *ALIASES, *ANALYSES, *META]))
        expanded = []
        for token in tokens:
            expanded += list(ALIASES.get(token, (token,)))
        unsupported = [t for t in expanded if t not in META and t != "data"
                       and t not in policy.formats and t not in policy.analyses]
        if unsupported:
            raise OutputSelectionError(f"crevice {command} does not write: {', '.join(unsupported)}")
        kinds = [t for t in expanded if t in FORMATS]
        chosen = [t for t in expanded if t in ANALYSES]
        both_skipped = [a for a in chosen if a in skip]
        if both_skipped:
            raise OutputSelectionError(f"{', '.join(both_skipped)} is both requested in --outputs and skipped")
        formats = set(policy.default_formats) if (not kinds and "none" not in expanded) else set(kinds)
        analyses = set(defaults) if not chosen else set(chosen)
        if "default" in expanded:
            formats |= set(policy.default_formats)
            analyses |= set(defaults)
        mode = "selected" if tokens else "default"
    formats.add("data")
    if formats & (set(VIEWERS) | {"viewer-html"}):
        formats.add("scenes")
    if "hydration" in analyses and "interactions" in policy.analyses:
        analyses.add("interactions")   # the hydration analysis includes the typed interactions
    analyses = {a for a in analyses if a in policy.analyses and a not in skip}
    return OutputPlan(command, mode, frozenset(formats), frozenset(analyses), tuple(tokens), skip)


@dataclass
class OutputLog:
    """What happened to optional outputs during a run (recorded in the manifest).

    Attributes
    ----------
    blocked : set of str
        Resolved paths a writer skipped because the plan does not want them.
    notes : dict of str to str
        Optional outputs that were wanted but not produced, with the reason
        (for example a fingerprint without a resolved channel).
    """

    blocked: set = field(default_factory=set)
    notes: dict = field(default_factory=dict)


_ACTIVE: ContextVar = ContextVar("crevice_output_plan", default=None)


@contextlib.contextmanager
def use_outputs(plan: OutputPlan | None) -> Iterator[OutputLog]:
    """Activate ``plan`` for every writer called inside the block.

    Parameters
    ----------
    plan : OutputPlan or None
        The plan; ``None`` writes everything (CREVICE 0.2 behaviour).

    Yields
    ------
    OutputLog
        Collects the skipped paths and notes of the block.
    """
    log = OutputLog()
    token = _ACTIVE.set((plan, log) if plan is not None else None)
    try:
        yield log
    finally:
        _ACTIVE.reset(token)


def active_plan() -> OutputPlan | None:
    """The plan in effect, or ``None`` when every output is written.

    Returns
    -------
    OutputPlan or None
    """
    state = _ACTIVE.get()
    return state[0] if state else None


def active_log() -> OutputLog | None:
    """The :class:`OutputLog` of the active plan, or ``None``.

    Returns
    -------
    OutputLog or None
    """
    state = _ACTIVE.get()
    return state[1] if state else None


def wants(path: str | Path) -> bool:
    """Whether the active plan writes ``path``; records the path when it does not.

    Writers call this before writing an optional file. With no active plan
    the answer is always ``True``.

    Parameters
    ----------
    path : str or Path
        Output file path.

    Returns
    -------
    bool
    """
    state = _ACTIVE.get()
    if state is None:
        return True
    plan, log = state
    if plan.wants(path):
        return True
    log.blocked.add(str(Path(path).resolve()))
    return False


def writes(kind: str) -> bool:
    """Whether the active plan writes files of ``kind`` (``True`` without a plan).

    Parameters
    ----------
    kind : str
        A file kind (:data:`FORMATS`) or ``"viewers"`` (any native viewer).

    Returns
    -------
    bool
    """
    plan = active_plan()
    if plan is None:
        return True
    if kind == "viewers":
        return bool(plan.viewers)
    return plan.writes(kind)


def note(name: str, reason: str) -> None:
    """Record that a wanted optional output ``name`` was not produced, with ``reason``.

    Parameters
    ----------
    name : str
        Output or analysis name.
    reason : str
        Why it was not produced.
    """
    log = active_log()
    if log is not None:
        log.notes[name] = reason


def drop_blocked(files: dict) -> dict:
    """Remove entries whose path was skipped by the active plan (in place).

    Parameters
    ----------
    files : dict of str to str
        Output key to path, as returned by the ``write_*`` functions.

    Returns
    -------
    dict
        The same dictionary.
    """
    log = active_log()
    if log is None or not log.blocked:
        return files
    for key, value in list(files.items()):
        if isinstance(value, (str, Path)) and str(Path(value).resolve()) in log.blocked:
            del files[key]
    return files


def plan_from(command: str, outputs) -> OutputPlan | None:
    """Turn a Python ``outputs=`` argument into a plan for ``command``.

    Parameters
    ----------
    command : str
        Command whose policy applies.
    outputs : None, str, sequence of str or OutputPlan
        ``None`` keeps the active plan (everything when none is active);
        ``"all"`` is :func:`resolve_outputs` with ``all_outputs=True``;
        other strings and sequences are
        ``--outputs`` tokens.

    Returns
    -------
    OutputPlan or None
    """
    if outputs is None:
        return active_plan()
    if isinstance(outputs, OutputPlan):
        return outputs
    if outputs == "all":
        return resolve_outputs(command, all_outputs=True)
    return resolve_outputs(command, outputs)


@contextlib.contextmanager
def outputs_for(command: str, outputs) -> Iterator[OutputPlan | None]:
    """Activate the plan for a Python ``outputs=`` argument, if one is given.

    Parameters
    ----------
    command : str
        Command whose policy applies.
    outputs : None, str, sequence of str or OutputPlan
        As :func:`plan_from`; ``None`` leaves the active plan unchanged.

    Yields
    ------
    OutputPlan or None
        The plan in effect.
    """
    if outputs is None:
        yield active_plan()
        return
    plan = plan_from(command, outputs)
    with use_outputs(plan):
        yield plan


def output_option(command: str):
    """Decorator giving a writer a keyword-only ``outputs`` argument for ``command``'s policy.

    ``outputs=None`` (the default) keeps the active plan, which writes
    everything when no plan is active (CREVICE 0.2 behaviour). Any other value
    (see :func:`plan_from`) applies to that call and every file it writes.

    Parameters
    ----------
    command : str
        Command whose output policy the writer follows (a key of :data:`COMMANDS`).

    Returns
    -------
    callable
        Decorator for writers with a keyword-only ``outputs`` parameter.
    """
    def decorate(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            outputs = kwargs.get("outputs")
            if outputs is None:
                return func(*args, **kwargs)
            with outputs_for(command, outputs):
                return func(*args, **kwargs)
        return wrapper
    return decorate

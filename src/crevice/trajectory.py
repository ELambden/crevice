"""Load MD trajectories and run CREVICE analyses frame by frame.

Loading
    :func:`load_trajectory_report` reads a topology/trajectory pair with
    MDAnalysis (a required dependency), keeps an atom selection
    (``"protein"`` by default) and stores only the requested frames as
    compact coordinate arrays; a :class:`~crevice.models.StructureFrame` is
    built for a frame only when it is accessed. Coordinates are in Å and
    times in ps, as reported by MDAnalysis. Periodic-boundary handling is
    explicit (``pbc="check"``, ``"none"`` or ``"unwrap"``).

Analysis
    :func:`analyze_trajectory` superimposes each frame on the first frame
    (Kabsch fit on CA atoms, :func:`align_frame_report`) and runs the
    requested analyses: channel profile, lining residues, cavities, tunnels,
    residue networks and numeric features. The channel axis and origin are
    fixed from the first frame, so every profile is measured along the same
    aligned axis. Results are :class:`~crevice.models.TrajectoryAnalysis`
    objects.

Main entry points: :func:`load_trajectory`, :func:`load_trajectory_report`,
:func:`inspect_trajectory`, :func:`analyze_trajectory`,
:func:`align_frame_report`, and the per-frame summaries
:func:`profile_timeseries`, :func:`cavity_timeseries` and
:func:`network_timeseries`. Command-line equivalents are
``crevice trajectory`` and ``crevice cavity-trajectory``.

Assumptions and limitations
    Frames must already be whole (molecules not split across periodic
    boundaries); the default check detects split bonds but does not prove
    every fragment is whole. Alignment removes rigid-body motion only.
    Per-frame geometry is subject to all the limitations of the underlying
    analyses, and frame-to-frame variation is not a free-energy or kinetic
    measurement.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

from .analysis import annotate_residues
from .cavities import detect_cavities
from .channels import pore_profile
from .ml import profile_features
from .models import Atom, FrameAnalysis, StructureFrame, TrajectoryAnalysis
from .networks import build_cavity_network, build_residue_network, network_metrics
from .parser import load_structure
from .tunnels import find_tunnels
from .radii import RadiusSet, radii_option


#: Refuse to materialise more than this many frames without an explicit choice.
DEFAULT_MAX_FRAMES = 500


@dataclass(frozen=True)
class Trajectory:
    """An ordered sequence of frames with optional times.

    Attributes
    ----------
    frames : sequence of StructureFrame
        The frames. For trajectories loaded from MD files this is a lazy
        sequence that builds each frame on access.
    source : str, default ""
        Path or description of the source.
    time_step : float or None, default None
        Time between consecutive stored frames (ps for MD input, already
        multiplied by the load stride), or ``None``.
    time_origin : float, default 0.0
        Time of the first stored frame.
    frame_times : tuple of float or None, optional
        Explicit time of every frame (ps for MD input, taken from the reader).
        Must have one entry per frame when given.

    Raises
    ------
    ValueError
        If ``frames`` is empty or ``frame_times`` has the wrong length.
    """
    frames: Sequence[StructureFrame]
    source: str = ""
    time_step: float | None = None
    time_origin: float = 0.0
    frame_times: tuple[float | None, ...] | None = None

    def __post_init__(self) -> None:
        if not self.frames:
            raise ValueError("Trajectory requires at least one frame")
        if self.frame_times is not None and len(self.frame_times)!=len(self.frames):
            raise ValueError("frame_times must match the number of frames")

    def times(self) -> tuple[float | None, ...]:
        """Time of each frame.

        Returns
        -------
        tuple of float or None
            ``frame_times`` when set; otherwise ``time_origin + i * time_step``;
            all ``None`` when neither is known.
        """
        if self.frame_times is not None:
            return self.frame_times
        if self.time_step is None:
            return tuple(None for _ in self.frames)
        return tuple(
            self.time_origin + index * self.time_step for index in range(len(self.frames))
        )


class _CoordinateFrames(Sequence):
    """Compact coordinate snapshots, materializing one atom frame on access."""
    def __init__(self,template,coordinates,indices):
        self._template=template;self._coordinates=tuple(coordinates);self._indices=tuple(indices)
        for xyz in self._coordinates:xyz.setflags(write=False)
    def __len__(self):return len(self._coordinates)
    def __getitem__(self,index):
        if isinstance(index,slice):return tuple(self[i] for i in range(*index.indices(len(self))))
        xyz=self._coordinates[index]
        atoms=tuple(replace(a,x=float(p[0]),y=float(p[1]),z=float(p[2])) for a,p in zip(self._template.atoms,xyz))
        return replace(self._template,atoms=atoms,frame_index=self._indices[index])


def trajectory_from_frames(
    frames: Iterable[StructureFrame],
    *,
    source: str = "",
    time_step: float | None = None,
) -> Trajectory:
    """Wrap already-loaded frames as a :class:`Trajectory`.

    If there are several frames and all have ``frame_index == 0`` (for example
    separately loaded PDB files), they are renumbered 0, 1, 2, ... so that
    results can be told apart.

    Parameters
    ----------
    frames : iterable of StructureFrame
    source : str, default ""
    time_step : float, optional
        Time between frames, in the units you want reported.

    Returns
    -------
    Trajectory

    Examples
    --------
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.trajectory import trajectory_from_frames
    >>> frame = StructureFrame((Atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),))
    >>> trajectory = trajectory_from_frames([frame, frame, frame], time_step=10.0)
    >>> [f.frame_index for f in trajectory.frames], trajectory.times()
    ([0, 1, 2], (0.0, 10.0, 20.0))
    """

    frames=tuple(frames)
    if len(frames)>1 and all(f.frame_index==0 for f in frames):
        frames=tuple(replace(f,frame_index=i) for i,f in enumerate(frames))
    return Trajectory(frames,source=source,time_step=time_step)


def load_static_frames(paths: Sequence[str | Path]) -> Trajectory:
    """Load several structure files as the frames of one trajectory.

    Each file is read with :func:`crevice.parser.load_structure` (default
    options, first model) and given ``frame_index`` equal to its position in
    ``paths``. No times are assigned. The files should contain the same atoms in
    the same order if they will be aligned or compared.

    Parameters
    ----------
    paths : sequence of str or pathlib.Path

    Returns
    -------
    Trajectory
        With ``source="static-file-sequence"``.
    """

    frames = tuple(replace(load_structure(path),frame_index=i) for i,path in enumerate(paths))
    return Trajectory(frames, source="static-file-sequence")


def _open_universe(topology: str | Path, trajectory: str | Path | None = None):
    """Open an MDAnalysis Universe, with explicit errors for a missing package,
    missing files or unreadable input.
    """
    try:
        import MDAnalysis as mda  # type: ignore[import-not-found]
    except ImportError as exc:
        raise ImportError(
            "MDAnalysis is a required CREVICE dependency but could not be imported; "
            "reinstall CREVICE with its dependencies or install MDAnalysis>=2.7"
        ) from exc
    topology_path = Path(topology)
    if not topology_path.exists():
        raise FileNotFoundError(f"Topology {topology_path} does not exist")
    if trajectory is not None and not Path(trajectory).exists():
        raise FileNotFoundError(f"Trajectory {Path(trajectory)} does not exist")
    try:
        return (mda.Universe(str(topology_path), str(trajectory)) if trajectory
                else mda.Universe(str(topology_path)))
    except Exception as exc:  # MDAnalysis raises a wide range of reader errors
        target = f"{topology_path} with {trajectory}" if trajectory else str(topology_path)
        raise ValueError(f"MDAnalysis could not read {target}: {exc}") from exc


def _frame_range(total: int, start: int, stop: int | None, stride: int) -> range:
    """Validated ``range(start, min(stop, total), stride)`` of frame indices."""
    if stride < 1:
        raise ValueError("stride must be at least 1")
    if start < 0:
        raise ValueError("start must be non-negative")
    if start >= total:
        raise ValueError(f"start {start} is beyond the last frame ({total} frames available)")
    end = total if stop is None else min(stop, total)
    if end <= start:
        raise ValueError(f"stop {stop} must be greater than start {start}")
    return range(start, end, stride)


def inspect_trajectory(
    topology: str | Path,
    trajectory: str | Path | None = None,
    *,
    selection: str = "protein",
) -> dict:
    """Describe a trajectory without materialising any frame.

    Cheap enough to run before committing to a load, so frame counts and
    selection sizes can be checked against available memory first.

    Parameters
    ----------
    topology : str or pathlib.Path
        Topology or single-frame structure file readable by MDAnalysis.
    trajectory : str or pathlib.Path, optional
        Trajectory file (XTC, DCD, ...). Omit to inspect the topology alone.
    selection : str, default "protein"
        MDAnalysis selection whose size is reported.

    Returns
    -------
    dict
        ``total_atom_count``, ``selected_atom_count``,
        ``selected_residue_count``, ``frame_count``, ``time_step_ps``,
        ``total_time_ps``, ``box_dimensions`` (Å and degrees, first frame),
        whether elements and chain IDs are present, and the segment IDs.

    Raises
    ------
    ImportError, FileNotFoundError, ValueError
        If MDAnalysis is missing, a file does not exist, or MDAnalysis cannot
        read the input or rejects the selection.
    """

    universe = _open_universe(topology, trajectory)
    try:
        selected = universe.select_atoms(selection)
    except Exception as exc:
        raise ValueError(f"MDAnalysis rejected selection {selection!r}: {exc}") from exc
    step = float(getattr(universe.trajectory, "dt", 0.0) or 0.0) or None
    frame_count = len(universe.trajectory)
    dimensions = getattr(universe, "dimensions", None)
    return {
        "topology": str(topology),
        "trajectory": str(trajectory) if trajectory else None,
        "selection": selection,
        "total_atom_count": len(universe.atoms),
        "selected_atom_count": len(selected),
        "selected_residue_count": len(selected.residues) if len(selected) else 0,
        "frame_count": frame_count,
        "time_step_ps": step,
        "total_time_ps": step * (frame_count - 1) if step and frame_count else None,
        "box_dimensions": [float(value) for value in dimensions] if dimensions is not None else None,
        "has_elements": bool(_optional_attribute(selected, "elements")),
        "has_chain_ids": any(value.strip() for value in _optional_attribute(selected, "chainIDs")),
        "segment_ids": sorted(set(_optional_attribute(selected, "segids"))),
    }


def load_trajectory(
    topology: str | Path,
    trajectory: str | Path | None = None,
    *,
    selection: str = "protein",
    start: int = 0,
    stop: int | None = None,
    stride: int = 1,
    max_frames: int | None = DEFAULT_MAX_FRAMES,
    pbc: str = "check",
) -> Trajectory:
    """Load a frame range from an MD trajectory.

    Same as :func:`load_trajectory_report` without the report; see it for
    the full description of each option.

    Parameters
    ----------
    topology : str or Path
        Topology or coordinate file readable by MDAnalysis (for example GRO).
    trajectory : str or Path, optional
        Trajectory file (for example XTC); ``None`` reads the topology's frames.
    selection : str, default "protein"
        MDAnalysis atom selection.
    start, stop, stride : int
        Frame range (``stop`` exclusive; ``None`` = to the end).
    max_frames : int or None
        Refuse to load more frames than this (``None``: no limit).
    pbc : {"check", "none"}, default "check"
        Periodic-boundary handling of the reader.

    Returns
    -------
    Trajectory

    Examples
    --------
    >>> from crevice.trajectory import load_trajectory
    >>> traj = load_trajectory("system.gro", "production.xtc", stride=10)  # doctest: +SKIP
    >>> len(traj.frames), traj.times()[:2]                                # doctest: +SKIP
    """

    return load_trajectory_report(
        topology, trajectory, selection=selection, start=start, stop=stop,
        stride=stride,max_frames=max_frames,pbc=pbc,
    )[0]


def load_trajectory_report(
    topology: str | Path,
    trajectory: str | Path | None = None,
    *,
    selection: str = "protein",
    start: int = 0,
    stop: int | None = None,
    stride: int = 1,
    max_frames: int | None = DEFAULT_MAX_FRAMES,
    pbc: str = "check",
) -> tuple[Trajectory, dict]:
    """Load a frame range and report what was read and what was skipped.

    Only the requested frames are materialised, so a long trajectory is
    subsampled at read time rather than loaded and then discarded. Elements are
    inferred with the shared parser rule when the topology has none. Default PBC
    checks reject detected split bonds; they do not prove every fragment is
    whole. Explicit unwrap requires a single bonded fragment and never guesses
    missing topology bonds.

    Parameters
    ----------
    topology : str or pathlib.Path
        Topology or structure file readable by MDAnalysis (GRO, PDB, TPR, PSF,
        ...).
    trajectory : str or pathlib.Path, optional
        Trajectory file; omit to load the topology's own coordinates.
    selection : str, default "protein"
        MDAnalysis selection of the atoms to keep. Atoms outside the MDAnalysis
        ``protein`` selection are marked ``hetero=True``.
    start : int, default 0
        First frame index (zero-based).
    stop : int, optional
        One past the last frame index; defaults to the end.
    stride : int, default 1
        Keep every ``stride``-th frame.
    max_frames : int or None, default 500
        Refuse to load more frames than this. Pass ``None`` to remove the
        limit, or raise ``stride``.
    pbc : {"check", "none", "unwrap"}, default "check"
        ``"check"``: in every frame with a box, compare bonded pairs (topology
        bonds within the selection, or else consecutive peptide C-N links) and
        raise if any pair is longer than 3 Å directly but shorter than 2.5 Å
        under the minimum-image convention. ``"none"``: no check.
        ``"unwrap"``: make the selection whole with MDAnalysis (requires
        topology bonds and exactly one bonded fragment).

    Returns
    -------
    trajectory : Trajectory
        Frames with 1-based atom serials within the selection, chain IDs from
        ``chainID`` (else ``segid``), ``frame_index`` set to the source frame
        number, and ``frame_times`` from the reader (ps).
    report : dict
        Reader and version, selection, atom counts (selected and excluded, by
        residue name), frame range, source and loaded time steps (ps), frame
        times, element and chain-ID sources and the PBC handling and checks
        performed. ``completeness_validated`` is always ``False``.

    Raises
    ------
    ImportError, FileNotFoundError
        If MDAnalysis is missing or a file does not exist.
    ValueError
        For an invalid ``pbc`` value, an empty or rejected selection, an invalid
        frame range, too many frames, a split bond (``pbc="check"``) or a
        failed unwrap.
    """

    if pbc not in {"check","none","unwrap"}:
        raise ValueError("pbc must be check, none or unwrap")
    universe = _open_universe(topology, trajectory)
    try:
        selected_atoms = universe.select_atoms(selection)
    except Exception as exc:
        raise ValueError(f"MDAnalysis rejected selection {selection!r}: {exc}") from exc
    if not len(selected_atoms):
        raise ValueError(f"Selection {selection!r} matched no atoms")

    total_frames = len(universe.trajectory)
    wanted = _frame_range(total_frames, start, stop, stride)
    if max_frames is not None and len(wanted) > max_frames:
        raise ValueError(
            f"{len(wanted)} frames exceeds max_frames={max_frames}; raise stride, "
            f"narrow start/stop, or pass max_frames explicitly to load them all"
        )

    from collections import Counter
    protein_indices = set(universe.select_atoms("protein").indices)
    selected_indices = set(selected_atoms.indices)
    excluded = universe.atoms[[i for i in range(len(universe.atoms)) if i not in selected_indices]]
    elements = _element_names(selected_atoms)
    chain_source, chain_ids = _chain_labels(selected_atoms)
    resnames = [str(value) for value in selected_atoms.resnames]
    names = [str(value) for value in selected_atoms.names]
    _require_radii(names, resnames, elements, selection, str(trajectory or topology))
    resids = [int(value) for value in selected_atoms.resids]

    pbc_pairs,pbc_pair_source=_periodic_check_pairs(selected_atoms,chain_ids,resids,names)
    pbc_checked=0
    coordinates=[]
    template=None
    frame_times=[]
    for output_index, frame_index in enumerate(wanted):
        ts=universe.trajectory[frame_index]
        frame_times.append(float(ts.time))
        if pbc=="unwrap":
            try:
                if len(selected_atoms.fragments)!=1:
                    raise ValueError("Automatic unwrap requires one bonded fragment; reimage multimer assemblies explicitly")
                selected_atoms.unwrap(compound="fragments",reference=None)
            except (AttributeError,ValueError) as exc:
                raise ValueError("PBC unwrap requires reliable topology bonds and a valid box; use a bonded topology or preprocessed whole coordinates") from exc
        positions = selected_atoms.positions
        if pbc!="none" and len(pbc_pairs) and ts.dimensions is not None:
            _check_periodic_pairs(positions,pbc_pairs,ts.dimensions,frame_index)
            pbc_checked+=1
        coordinates.append(positions.copy())
        if template is None:
            icodes=_optional_attribute(selected_atoms,"icodes")
            atoms=tuple(Atom(serial=serial,name=names[serial-1],resname=resnames[serial-1],
                    chain_id=chain_ids[serial-1],resid=resids[serial-1],
                    x=float(position[0]),y=float(position[1]),z=float(position[2]),element=elements[serial-1],
                    icode=icodes[serial-1] if icodes else '',
                    hetero=int(selected_atoms.indices[serial-1]) not in protein_indices)
                    for serial,position in enumerate(positions,start=1))
            template=StructureFrame(atoms,source=str(trajectory or topology),model_index=0,frame_index=frame_index)

    step = float(getattr(universe.trajectory, "dt", 0.0) or 0.0) or None
    time_step = step * stride if step else None
    report = {
        "reader": "mdanalysis",
        "reader_version": _mdanalysis_version(),
        "topology": str(topology),
        "trajectory": str(trajectory) if trajectory else None,
        "selection": selection,
        "total_atom_count": len(universe.atoms),
        "excluded_atom_count": len(excluded),
        "excluded_residue_atom_counts": dict(sorted(Counter(str(r) for r in excluded.resnames).items())),
        "selected_residue_atom_counts": dict(sorted(Counter(resnames).items())),
        "hetero_definition": "outside MDAnalysis protein selection",
        "selected_atom_count": len(selected_atoms),
        "selected_residue_count": len(selected_atoms.residues),
        "available_frame_count": total_frames,
        "loaded_frame_count":len(coordinates),
        "coordinate_storage":"compact arrays; atom frames materialized on access",
        "frame_indices": [wanted.start, wanted.stop, wanted.step],
        "source_time_step_ps": step,
        "loaded_time_step_ps": time_step,
        "frame_times_ps":frame_times,
        "time_origin_source":"reader timestamps",
        "element_source": "reader" if _optional_attribute(selected_atoms, "elements") else "name inference",
        "chain_id_source": chain_source,
        "chain_ids": sorted(set(chain_ids)),
        "periodic_boundary_handling":pbc,
        "periodic_check_pair_source":pbc_pair_source,
        "periodic_check_pair_count":len(pbc_pairs),
        "periodic_checked_frame_count":pbc_checked,
        "periodic_check_scope":"known bonds when available; otherwise consecutive peptide C-N links only",
        "completeness_validated": False,
    }
    return Trajectory(
        _CoordinateFrames(template,coordinates,wanted),source=str(trajectory or topology),
        time_step=time_step, time_origin=frame_times[0],frame_times=tuple(frame_times),
    ), report


def _optional_attribute(selected_atoms, name: str) -> list[str]:
    """Read a topology attribute the reader may simply not provide."""

    try:
        values = getattr(selected_atoms, name)
    except (AttributeError, ValueError, KeyError):
        return []
    return [] if values is None else [str(value) for value in values]


def _element_names(selected_atoms) -> list[str]:
    """Use reader elements when present, else the shared inference rule."""

    from .radii import infer_element

    supplied = _optional_attribute(selected_atoms, "elements")
    names = [str(value) for value in selected_atoms.names]
    resnames = [str(value) for value in selected_atoms.resnames]
    return [
        infer_element(names[index], supplied[index] if index < len(supplied) else "",
                      resnames[index])
        for index in range(len(names))
    ]


def _require_radii(names, resnames, elements, selection, source) -> None:
    """Refuse an MD selection that contains atoms without a radius in the set in effect.

    Trajectory frames must stay atom-for-atom identical to the MDAnalysis
    selection (hydration and obstacle code index into it), so such atoms are
    not dropped silently: the error explains how to supply radii or exclude
    them from the selection by residue name.

    Raises
    ------
    ValueError
        If any selected atom has no radius
        (:meth:`crevice.radii.RadiusSet.has_radius`).
    """
    from .radii import partition_atoms_with_radii, unrecognised_element_message
    atoms = [Atom(i + 1, n, r, "", 0, 0.0, 0.0, 0.0, e) for i, (n, r, e) in enumerate(zip(names, resnames, elements))]
    kept, by_element, by_name = partition_atoms_with_radii(atoms)
    if by_element:
        residues = " ".join(sorted({key.split(":", 1)[0] for key in by_name}))
        raise ValueError(unrecognised_element_message(by_element, by_name, source, excluded=False)
                         + f" Or exclude them from the MDAnalysis selection, for example "
                           f"'({selection}) and not resname {residues}'.")


def _chain_labels(selected_atoms) -> tuple[str, list[str]]:
    """Prefer chain IDs; fall back to segment IDs, recording which was used."""

    chain_ids = _optional_attribute(selected_atoms, "chainIDs")
    if any(value.strip() for value in chain_ids):
        return "chainID", chain_ids
    segids = _optional_attribute(selected_atoms, "segids")
    if any(value.strip() for value in segids):
        return "segid", segids
    return "none", ["" for _ in range(len(selected_atoms))]


def _mdanalysis_version() -> str:
    """Installed MDAnalysis version string."""
    import MDAnalysis  # type: ignore[import-not-found]

    return str(getattr(MDAnalysis, "__version__", "unknown"))


@radii_option
def analyze_trajectory(
    trajectory: Trajectory | Sequence[StructureFrame],
    *,
    analyses: Sequence[str] = ("profile",),
    stride: int = 1,
    profile_kwargs: dict | None = None,
    cavity_kwargs: dict | None = None,
    tunnel_kwargs: dict | None = None,
    network_kwargs: dict | None = None,
    align: bool = True,
    alignment_residues: Sequence[str] | None = None,
    on_profile_error: str = "record",
    retain_networks: bool = True,
    lining_residues: Sequence[str] = (),
    network_occupancy: float = .75,
    progress_callback=None,
    radii: RadiusSet | str | None = None,
    probe_fallback: bool | None = None,
) -> TrajectoryAnalysis:
    """Run selected CREVICE analyses over a trajectory or frame sequence.

    Steps:

    1. The first frame is the reference. If a profile, residues or features
       are requested, the channel is resolved once on the reference: in
       connected mode (``search_radius > 0``) its axis direction and middle
       point become the fixed ``axis``/``origin`` for every frame (unless an
       ``origin`` was given); with ``search_radius=0`` the reference axis and
       atom centroid are used.
    2. Every ``stride``-th frame is superimposed on the reference
       (:func:`align_frame_report`) when ``align`` is true.
    3. The requested analyses run on the aligned frame with the pinned
       enclosure probe. If that probe does not resolve the channel in a frame
       and the per-frame fallback is on (the default when the probe was chosen
       automatically), the frame is profiled again with the automatic choice
       made on that frame alone (same axis and origin). The frame records
       which probe it used (``metadata["enclosure_probe"]``).
    4. When a profile exists, residues limiting its constriction are also
       recorded (``metadata["constriction"]``).

    Parameters
    ----------
    trajectory : Trajectory or sequence of StructureFrame
    analyses : sequence of str, default ("profile",)
        Any of ``"profile"``, ``"residues"`` (lining residues along the
        profile), ``"cavities"``, ``"tunnels"``, ``"network"`` and
        ``"features"``. Duplicates are ignored.
    stride : int, default 1
        Analyse every ``stride``-th frame of ``trajectory``.
    profile_kwargs : dict, optional
        Options for :func:`crevice.channels.pore_profile`. ``axis`` and
        ``origin`` are overwritten from the reference frame as described above.
        Without ``enclosure_radius`` (or with ``None``) the enclosure probe is
        chosen automatically once, on the reference frame
        (:func:`crevice.sections.connected_profile`), and that value is then
        passed explicitly for every frame, so frames never re-run the
        automatic choice. A number is used unchanged for every frame. The
        ``crevice trajectory`` CLI does the same (``--enclosure-radius auto``,
        the default, or a number). The choice and its reason are recorded in
        ``metadata["enclosure_probe"]``. A frame that the pinned probe cannot
        resolve may be retried with a per-frame automatic choice; see
        ``probe_fallback``.
    cavity_kwargs : dict, optional
        Options for :func:`crevice.cavities.detect_cavities`.
    tunnel_kwargs : dict, optional
        Options for :func:`crevice.tunnels.find_tunnels`.
    network_kwargs : dict, optional
        Options for :func:`crevice.networks.build_cavity_network` (when a
        profile exists) or :func:`crevice.networks.build_residue_network`.
        These default to the surface-gap metric; ``crevice trajectory`` passes
        ``distance_metric="center"``.
    align : bool, default True
        Superimpose every frame on the first frame before analysis.
    alignment_residues : sequence of str, optional
        Residue labels (``"A:LEU42"``) restricting the alignment atoms.
    on_profile_error : {"record", "raise"}, default "record"
        What to do when a channel cannot be resolved: record the reason in the
        frame metadata and continue, or raise
        :class:`crevice.sections.ChannelResolutionError`. If the reference
        frame fails under ``"record"``, no frame gets a profile.
    retain_networks : bool, default True
        Keep each frame's network object. If ``False`` (and ``"network"`` is
        requested), networks are summarised on the fly by
        :class:`crevice.residue_dynamics.ResidueDynamics` and only the summary
        is kept in ``metadata["residue_dynamics"]``, which saves memory.
    lining_residues : sequence of str, default ()
        Residues passed to the dynamics summary (``retain_networks=False``).
    network_occupancy : float, default 0.75
        Minimum edge occupancy (fraction of frames) for the dynamics summary.
    progress_callback : callable, optional
        Called with each :class:`~crevice.models.FrameAnalysis` as it
        completes.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.
    probe_fallback : bool or None, default None
        Per-frame fallback for the enclosure probe. When the pinned probe
        does not resolve a frame's channel (``ChannelResolutionError``), retry
        that frame with ``enclosure_radius=None``, i.e. the automatic choice
        (:func:`crevice.sections.connected_profile`) made on that frame, with
        the same fixed axis and origin. ``None`` (default) enables the
        fallback only when the probe itself was chosen automatically (no
        ``enclosure_radius`` in ``profile_kwargs``); an explicit probe stays
        strict. ``True`` enables it for an explicit probe as well; ``False``
        disables it. A fallback frame's profile is measured with a different
        wall definition (its own probe); radii inside a channel's resolving
        window barely depend on the probe, but the mouths can move, so the
        fallback frames are listed and counted. Has no effect with
        ``search_radius=0`` (no enclosure probe) or when the reference frame
        is unresolved (no axis exists).

    Returns
    -------
    TrajectoryAnalysis
        One :class:`~crevice.models.FrameAnalysis` per analysed frame. Frame
        metadata records the alignment report, ``profile_status`` (for example
        ``"resolved"``, ``"unresolved"``, ``"unvalidated_axis_scan"``) and the
        reason for a missing profile. Trajectory metadata includes how many
        profiles were resolved, per residue the number of frames in which it
        limited the constriction, and ``enclosure_probe``: ``mode``
        (``"auto"``, ``"explicit"`` or ``"not_used"`` for ``search_radius=0``),
        ``chosen_A`` (the pinned probe; ``None`` if the automatic choice
        failed), ``reason``, in automatic mode the full ``selection`` record
        of the reference frame, and ``fallback``: ``enabled``, ``policy``,
        ``pinned_frames``, ``fallback_frames``, ``unresolved_frames`` (counts),
        ``fallback_frame_indices`` and ``fallback_probes_A`` (source frame
        index to probe). ``None`` when no profile was requested.

        Each frame's ``metadata["enclosure_probe"]`` records ``source``
        (``"pinned"``, ``"fallback"``, ``"none"`` when no probe resolved the
        frame, or ``"not_attempted"``), ``radius_A`` (the effective probe of
        the profile, ``None`` if unresolved), ``pinned_A``, and for a retried
        frame ``pinned_failure`` and the fallback's ``reason`` or
        ``fallback_failure``.

    Raises
    ------
    ValueError
        For an invalid ``stride``, unknown analysis names or an invalid
        ``on_profile_error``; alignment and analysis errors propagate.
    crevice.sections.ChannelResolutionError
        With ``on_profile_error="raise"``.

    Examples
    --------
    Two identical frames of a synthetic cylinder, profiled along z:

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.trajectory import analyze_trajectory
    >>> atoms = [Atom(12 * k + j + 1, f"C{j}", "ALA", "A", k + 1,
    ...               6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6),
    ...               -10.5 + 1.5 * k, "C")
    ...          for k in range(15) for j in range(12)]
    >>> frame = StructureFrame(tuple(atoms))
    >>> result = analyze_trajectory([frame, frame], profile_kwargs={"axis": "z"})
    >>> [round(f.profile.min_radius, 2) for f in result.frames]
    [4.3, 4.3]
    """

    if stride < 1:
        raise ValueError("stride must be at least 1")
    traj = trajectory if isinstance(trajectory, Trajectory) else trajectory_from_frames(trajectory)
    selected = tuple(dict.fromkeys(analyses))
    valid = {"profile", "residues", "cavities", "tunnels", "network", "features"}
    unknown = set(selected) - valid
    if unknown:
        raise ValueError(f"Unknown trajectory analyses: {', '.join(sorted(unknown))}")

    if on_profile_error not in {"record","raise"}:
        raise ValueError("on_profile_error must be record or raise")
    from .sections import ChannelResolutionError
    reference_error=None
    reference_profile=None
    enclosure_probe=None
    accumulator=None
    if "network" in selected and not retain_networks:
        from .residue_dynamics import ResidueDynamics
        accumulator=ResidueDynamics(lining_residues=lining_residues,min_occupancy=network_occupancy)
    profile_kwargs = dict(profile_kwargs or {})
    cavity_kwargs = dict(cavity_kwargs or {})
    tunnel_kwargs = dict(tunnel_kwargs or {})
    network_kwargs = dict(network_kwargs or {})
    frame_results: list[FrameAnalysis] = []
    times = traj.times()
    reference = traj.frames[0]
    if any(item in selected for item in ("profile", "residues", "features")):
        from .geometry import axis_from_name

        atoms = reference.selected_atoms(
            include_hydrogen=profile_kwargs.get("include_hydrogen", False),
            include_hetero=profile_kwargs.get("include_hetero", True))
        if profile_kwargs.get("search_radius", 6.0) > 0:
            automatic_probe = profile_kwargs.get("enclosure_radius") is None
            try:
                reference_profile=pore_profile(reference,**profile_kwargs)
                profile_kwargs["axis"]=reference_profile.axis_direction
                if automatic_probe:
                    # Automatic probe: chosen once on the reference frame, then
                    # pinned so every frame uses the same enclosure definition
                    # (frames never re-run the automatic choice).
                    selection=reference_profile.metadata["enclosure_probe_selection"]
                    profile_kwargs["enclosure_radius"]=selection["chosen_A"]
                    enclosure_probe={"mode":"auto","chosen_A":selection["chosen_A"],"reason":selection["reason"],
                                     "chosen_on":"reference frame (first frame)","applied":"fixed for every frame",
                                     "selection":selection}
                else:
                    enclosure_probe={"mode":"explicit","chosen_A":profile_kwargs["enclosure_radius"],
                                     "reason":"enclosure_radius given explicitly",
                                     "applied":"fixed for every frame"}
                if profile_kwargs.get("origin") is None:
                    profile_kwargs["origin"]=reference_profile.points[len(reference_profile.points)//2].position
            except ChannelResolutionError as exc:
                selection=getattr(exc,"probe_selection",None)
                enclosure_probe=({"mode":"auto","chosen_A":None,"reason":selection["reason"],
                                  "chosen_on":"reference frame (first frame)","applied":"no frame profiled",
                                  "selection":selection} if automatic_probe and selection is not None else
                                 {"mode":"auto" if automatic_probe else "explicit",
                                  "chosen_A":profile_kwargs.get("enclosure_radius"),
                                  "reason":"reference frame unresolved","applied":"no frame profiled"})
                if on_profile_error=="raise":raise
                reference_error=str(exc)
        else:
            enclosure_probe={"mode":"not_used","chosen_A":None,
                             "reason":"search_radius=0 fixed-axis scan does not use an enclosure probe"}
            origin, direction, _ = axis_from_name(profile_kwargs.get("axis", "auto"), atoms)
            profile_kwargs["axis"] = direction
            if profile_kwargs.get("origin") is None:profile_kwargs["origin"]=origin

    fallback_enabled=False
    if enclosure_probe is not None and enclosure_probe.get("mode") in {"auto","explicit"} and reference_error is None:
        fallback_enabled=(enclosure_probe["mode"]=="auto") if probe_fallback is None else bool(probe_fallback)
    fallback_counts={"pinned":0,"fallback":0,"none":0}
    fallback_probes={}

    for output_index, frame_index in enumerate(range(0, len(traj.frames), stride)):
        frame = traj.frames[frame_index]
        probe_use=None
        alignment={"status":"not_requested"}
        if align:
            frame,alignment=align_frame_report(frame,reference,residues=alignment_residues)
        profile_status="not_requested"
        profile_reason=None
        constriction=None
        profile = None
        contacts = ()
        cavities = ()
        tunnels = ()
        network = None
        features = {}
        if "profile" in selected or "residues" in selected or "features" in selected:
            if reference_error is not None:
                profile_status="not_attempted_reference_unresolved"
                profile_reason=reference_error
                probe_use={"source":"not_attempted","radius_A":None,"pinned_A":None}
            else:
                pinned=profile_kwargs.get("enclosure_radius")
                try:
                    profile=pore_profile(frame,**profile_kwargs)
                    profile_status=profile.metadata.get("status","resolved")
                    probe_use={"source":"pinned","radius_A":profile.metadata.get("enclosure_radius_A",pinned),"pinned_A":pinned}
                except ChannelResolutionError as exc:
                    probe_use={"source":"none","radius_A":None,"pinned_A":pinned,"pinned_failure":str(exc)}
                    if fallback_enabled:
                        try:
                            # Per-frame fallback: the automatic choice on this frame alone,
                            # with the reference axis and origin unchanged.
                            profile=pore_profile(frame,**{**profile_kwargs,"enclosure_radius":None})
                            profile_status=profile.metadata.get("status","resolved")
                            selection=profile.metadata.get("enclosure_probe_selection") or {}
                            probe_use.update(source="fallback",radius_A=profile.metadata.get("enclosure_radius_A",selection.get("chosen_A")),
                                             reason=selection.get("reason"))
                        except ChannelResolutionError as again:
                            probe_use["fallback_failure"]=str(again)
                            if on_profile_error=="raise":raise
                            profile_status="unresolved";profile_reason=f"{exc}; per-frame automatic probe also failed: {again}"
                    else:
                        if on_profile_error=="raise":raise
                        probe_use["fallback"]="disabled"
                        profile_status="unresolved";profile_reason=str(exc)
                if probe_use["source"] in fallback_counts:fallback_counts[probe_use["source"]]+=1
                if probe_use["source"]=="fallback":
                    fallback_probes[frame.frame_index if frame.frame_index is not None else frame_index]=probe_use["radius_A"]
        if profile is not None:
            from .residue_evidence import profile_constriction_evidence
            constriction=profile_constriction_evidence(frame,profile)
        if profile is not None and ("residues" in selected or "features" in selected):
            contacts = annotate_residues(frame,profile.points,include_hydrogen=profile_kwargs.get("include_hydrogen",False),
                                         include_hetero=profile_kwargs.get("include_hetero",True))
        if "cavities" in selected:
            cavities = detect_cavities(frame, **cavity_kwargs)
        if "tunnels" in selected:
            tunnels = find_tunnels(frame, **tunnel_kwargs)
        if "network" in selected:
            if profile is not None:
                if not contacts:
                    contacts = annotate_residues(frame, profile.points)
                network = build_cavity_network(frame, profile, contacts=contacts, **network_kwargs)
            else:
                network = build_residue_network(frame, **network_kwargs)
        if "features" in selected:
            features = profile_features(profile, contacts) if profile else {}
            if network is not None:
                metrics = network_metrics(network)
                features.update(
                    {
                        "network_node_count": float(metrics["node_count"]),
                        "network_edge_count": float(metrics["edge_count"]),
                        "network_density": float(metrics["density"]),
                    }
                )
        if accumulator is not None and network is not None:
            accumulator.add(frame,network,frame_index=frame.frame_index,time=times[frame_index])
        frame_results.append(
            FrameAnalysis(
                frame_index=frame.frame_index if frame.frame_index is not None else frame_index,
                time=times[frame_index],
                profile=profile,
                cavities=tuple(cavities),
                tunnels=tuple(tunnels),
                residue_contacts=tuple(contacts),
                network=network if retain_networks else None,
                features=features,
                metadata={"trajectory_index":output_index,"source_index":frame_index,"alignment":alignment,
                          "profile_status":profile_status,"profile_reason":profile_reason,"enclosure_probe":probe_use,
                          "constriction":constriction,
                          "network_summary_only":not retain_networks,
                          "network_node_count":len(network.nodes) if network else 0,
                          "network_edge_count":len(network.edges) if network else 0},
            )
        )
        if progress_callback is not None:
            progress_callback(frame_results[-1])
    if enclosure_probe is not None and enclosure_probe.get("mode") in {"auto","explicit"}:
        enclosure_probe["fallback"]={
            "enabled":fallback_enabled,
            "policy":("per-frame automatic choice when the pinned probe does not resolve a frame"
                      if fallback_enabled else
                      "off: explicit probe kept strict (enable with probe_fallback=True / --probe-fallback)"
                      if enclosure_probe["mode"]=="explicit" else
                      "off" if reference_error is None else "not applicable: reference frame unresolved"),
            "pinned_frames":fallback_counts["pinned"],"fallback_frames":fallback_counts["fallback"],
            "unresolved_frames":fallback_counts["none"],
            "fallback_frame_indices":sorted(fallback_probes),
            "fallback_probes_A":{str(k):v for k,v in sorted(fallback_probes.items())}}
    from collections import Counter
    constraint_frames=[f for f in frame_results if f.metadata.get("constriction") and f.metadata['constriction']['status']=='measured_path_constraints']
    constraints=Counter(row['residue'] for f in constraint_frames for row in f.metadata['constriction']['residues'])
    constraint_rows=[{'residue':residue,'limiting_frames':count,'eligible_profile_frames':len(constraint_frames),
                     'total_frames':len(frame_results),'conditional_frequency':count/len(constraint_frames)}
                    for residue,count in sorted(constraints.items(),key=lambda item:(-item[1],item[0]))]
    return TrajectoryAnalysis(
        frames=tuple(frame_results),
        analyses=selected,
        metadata={"source": traj.source, "stride": stride, "input_frame_count": len(traj.frames),
                  "aligned_to_first_frame": align, "alignment_residues": alignment_residues,
                  "profile_reference":"first_frame_fixed_axis_and_origin",
                  "reference_profile_status":"unresolved" if reference_error else "available",
                  "reference_profile_reason":reference_error,
                  "enclosure_probe":enclosure_probe,
                  "resolved_profile_count":sum(f.profile is not None for f in frame_results),
                  "constriction_residue_evidence":constraint_rows,
                  "residue_dynamics":accumulator.finish(aligned=align) if accumulator is not None else None},
    )


def align_frame_report(frame: StructureFrame, reference: StructureFrame,
                *, residues: Sequence[str] | None = None) -> tuple[StructureFrame, dict]:
    """Rigid least-squares fit with matching CA atoms (heavy atoms if <3 CA).

    Coordinates must already be made whole across periodic boundaries. Residue
    labels may restrict the fit to a stable scaffold.

    The fit atoms are the non-hydrogen (``element`` not ``"H"`` or ``"D"``), non-HETATM atoms
    (optionally only those in ``residues``), matched by residue, atom name and
    alternate location. If at least three CA atoms match they alone are used;
    otherwise all matched atoms. The optimal proper rotation (Kabsch, with a
    reflection correction) and translation are applied to *every* atom of
    ``frame``.

    Parameters
    ----------
    frame : StructureFrame
        Frame to move.
    reference : StructureFrame
        Frame to fit onto.
    residues : sequence of str, optional
        Residue labels to use for the fit; all must exist in ``reference``.

    Returns
    -------
    aligned : StructureFrame
        ``frame`` with transformed coordinates.
    report : dict
        ``atom_count``, ``atom_selection`` (``"CA"`` or
        ``"heavy_atom_fallback"``), ``rmsd_before_A`` and ``rmsd_after_A`` over
        the fit atoms, ``rotation_rows`` and ``translation_A`` with the
        convention ``x @ rotation + translation``.

    Raises
    ------
    ValueError
        If ``residues`` is empty or contains unknown labels, the selected atom
        identities are ambiguous or differ between the frames, fewer than three
        atoms match, or the fit atoms are collinear.

    Examples
    --------
    >>> from dataclasses import replace
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.trajectory import align_frame_report
    >>> xyz = [(0.0, 0.0, 0.0), (3.8, 0.0, 0.0), (3.8, 3.8, 0.0), (3.8, 3.8, 3.8)]
    >>> reference = StructureFrame(tuple(Atom(i + 1, "CA", "ALA", "A", i + 1, *p, "C")
    ...                                  for i, p in enumerate(xyz)))
    >>> moved = StructureFrame(tuple(replace(a, x=10.0 - a.y, y=a.x) for a in reference.atoms))
    >>> aligned, report = align_frame_report(moved, reference)
    >>> report["atom_selection"], round(report["rmsd_after_A"], 6)
    ('CA', 0.0)
    """
    import numpy as np

    def key(atom):
        return atom.residue_key, atom.name, atom.altloc

    def selected(source):
        atoms = [a for a in source.atoms if a.element.upper() not in {"H", "D"} and not a.hetero
                 and (residues is None or a.residue_key.label in residues)]
        mapping = {key(a): a for a in atoms}
        if len(mapping) != len(atoms):
            raise ValueError("Ambiguous atom identities in alignment selection")
        return mapping

    if residues is not None:
        available={a.residue_key.label for a in reference.atoms}
        if not residues or set(residues)-available:
            raise ValueError("Alignment selection contains unknown residues or is empty")
    mobile, target = selected(frame), selected(reference)
    if mobile.keys() != target.keys():
        raise ValueError("Alignment requires identical selected atom identities")
    keys = [k for k in target if k[1] == "CA"]
    if len(keys) < 3:
        keys = list(target)
    if len(keys) < 3:
        raise ValueError("Alignment requires at least three matching non-collinear atoms")
    x, y = (np.asarray([mapping[k].coord for k in keys]) for mapping in (mobile, target))
    xc, yc = x.mean(axis=0), y.mean(axis=0)
    if np.linalg.matrix_rank(x - xc) < 2 or np.linalg.matrix_rank(y - yc) < 2:
        raise ValueError("Alignment atoms must not be collinear")
    u, _, vt = np.linalg.svd((x - xc).T @ (y - yc))
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    coords = (np.asarray([a.coord for a in frame.atoms]) - xc) @ rotation + yc
    aligned=replace(frame,atoms=tuple(replace(a,x=float(p[0]),y=float(p[1]),z=float(p[2]))
                                      for a,p in zip(frame.atoms,coords)))
    fit=(x-xc)@rotation+yc
    report={"status":"fitted","method":"Kabsch proper rigid least-squares",
            "atom_count":len(keys),"atom_selection":"CA" if all(k[1]=="CA" for k in keys) else "heavy_atom_fallback",
            "rmsd_before_A":float(np.sqrt(np.mean(np.sum((x-y)**2,axis=1)))),
            "rmsd_after_A":float(np.sqrt(np.mean(np.sum((fit-y)**2,axis=1)))),
            "rotation_rows":rotation.tolist(),"translation_A":(yc-xc@rotation).tolist(),
            "transform_convention":"row coordinates: x @ rotation + translation",
            "selected_residues":list(residues) if residues is not None else None,
            "periodic_repair_performed":False}
    return aligned,report


def align_frame(frame: StructureFrame,reference: StructureFrame,*,residues=None)->StructureFrame:
    """Fit a frame to matching reference atoms and return only the moved frame.

    Parameters
    ----------
    frame, reference : StructureFrame
    residues : sequence of str, optional
        Residue labels restricting the fit atoms.

    Returns
    -------
    StructureFrame

    See Also
    --------
    align_frame_report : Same fit, with RMSD and transform diagnostics.
    """
    return align_frame_report(frame,reference,residues=residues)[0]


def profile_timeseries(trajectory: Trajectory | Sequence[StructureFrame], **kwargs) -> list[dict[str, float | int | None]]:
    """Per-frame profile radius summaries.

    Runs :func:`analyze_trajectory` with ``analyses=("profile",)``, alignment to
    the first frame, and ``kwargs`` as ``profile_kwargs``.

    Parameters
    ----------
    trajectory : Trajectory or sequence of StructureFrame
    **kwargs
        Options for :func:`crevice.channels.pore_profile`.

    Returns
    -------
    list of dict
        ``frame_index``, ``time``, ``min_radius``, ``mean_radius`` and
        ``max_radius`` (Å; ``None`` when the frame has no profile).

    Examples
    --------
    Fixed-axis scans (``search_radius=0``) are fast:

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.trajectory import profile_timeseries
    >>> atoms = [Atom(12 * k + j + 1, f"C{j}", "ALA", "A", k + 1,
    ...               6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6),
    ...               -10.5 + 1.5 * k, "C")
    ...          for k in range(15) for j in range(12)]
    >>> frame = StructureFrame(tuple(atoms))
    >>> rows = profile_timeseries([frame, frame], axis="z", search_radius=0, samples=5)
    >>> [(row["frame_index"], round(row["min_radius"], 2)) for row in rows]
    [(0, 4.3), (1, 4.3)]
    """

    result = analyze_trajectory(trajectory, analyses=("profile",), profile_kwargs=kwargs)
    return [
        {
            "frame_index": frame.frame_index,
            "time": frame.time,
            "min_radius": frame.profile.min_radius if frame.profile else None,
            "mean_radius": frame.profile.mean_radius if frame.profile else None,
            "max_radius": frame.profile.max_radius if frame.profile else None,
        }
        for frame in result.frames
    ]


def cavity_timeseries(trajectory: Trajectory | Sequence[StructureFrame], **kwargs) -> list[dict[str, float | int | None]]:
    """Per-frame cavity counts and total cavity volume.

    Runs :func:`analyze_trajectory` with ``analyses=("cavities",)``, alignment
    to the first frame, and ``kwargs`` as ``cavity_kwargs``.

    Parameters
    ----------
    trajectory : Trajectory or sequence of StructureFrame
        Frames to analyse.
    **kwargs
        Options for :func:`crevice.cavities.detect_cavities`.

    Returns
    -------
    list of dict
        ``frame_index``, ``time``, ``cavity_count`` and
        ``total_cavity_volume`` (Å³, sum over the returned cavities, which are
        capped at ``max_cavities``).
    """

    result = analyze_trajectory(trajectory, analyses=("cavities",), cavity_kwargs=kwargs)
    return [
        {
            "frame_index": frame.frame_index,
            "time": frame.time,
            "cavity_count": len(frame.cavities),
            "total_cavity_volume": sum(cavity.volume for cavity in frame.cavities),
        }
        for frame in result.frames
    ]


def network_timeseries(trajectory: Trajectory | Sequence[StructureFrame], **kwargs) -> list[dict[str, float | int | None]]:
    """Per-frame residue-network size and density.

    Runs :func:`analyze_trajectory` with ``analyses=("network",)``, alignment
    to the first frame, and ``kwargs`` as ``network_kwargs``. Without a profile
    this builds :func:`crevice.networks.build_residue_network` networks with its
    defaults (surface-gap metric, heteroatoms excluded).

    Parameters
    ----------
    trajectory : Trajectory or sequence of StructureFrame
        Frames to analyse.
    **kwargs
        Options for :func:`crevice.networks.build_residue_network`.

    Returns
    -------
    list of dict
        ``frame_index``, ``time``, ``node_count``, ``edge_count`` and
        ``density`` (see :func:`crevice.networks.network_metrics`).
    """

    result = analyze_trajectory(trajectory, analyses=("network",), network_kwargs=kwargs)
    rows = []
    for frame in result.frames:
        metrics = network_metrics(frame.network) if frame.network else {"node_count": 0, "edge_count": 0, "density": 0.0}
        rows.append(
            {
                "frame_index": frame.frame_index,
                "time": frame.time,
                "node_count": metrics["node_count"],
                "edge_count": metrics["edge_count"],
                "density": metrics["density"],
            }
        )
    return rows


def _periodic_check_pairs(selected_atoms,chain_ids,resids,names):
    """Use supplied bonds, else an explicitly limited peptide-link screen."""
    import numpy as np
    from MDAnalysis.exceptions import NoDataError
    lookup={int(v):i for i,v in enumerate(selected_atoms.indices)}
    try:
        bonds=selected_atoms.universe.bonds.to_indices()
        pairs=[(lookup[int(a)],lookup[int(b)]) for a,b in bonds if int(a) in lookup and int(b) in lookup]
        if pairs:return np.asarray(pairs,dtype=int),'topology_bonds'
    except (NoDataError,AttributeError):pass
    atoms={(chain,resid,name):i for i,(chain,resid,name) in enumerate(zip(chain_ids,resids,names))}
    pairs=[(i,atoms[(chain,resid+1,'N')]) for (chain,resid,name),i in atoms.items()
           if name=='C' and (chain,resid+1,'N') in atoms]
    return np.asarray(pairs,dtype=int).reshape(-1,2),'consecutive_peptide_C_N_only'


def _check_periodic_pairs(positions,pairs,box,frame_index):
    """Raise if any pair looks split across the box: direct distance > 3 Å,
    minimum-image distance < 2.5 Å and at least 0.5 Å shorter than direct.
    """
    import numpy as np
    from MDAnalysis.lib.distances import calc_bonds
    box=np.asarray(box)
    if not np.isfinite(box).all() or np.any(box[:3]<=0):
        raise ValueError('Periodic checks require finite positive box lengths')
    left=positions[pairs[:,0]];right=positions[pairs[:,1]]
    raw=np.linalg.norm(left-right,axis=1);mic=calc_bonds(left,right,box=box)
    split=(raw>3)&(mic<2.5)&(raw>mic+.5)
    if split.any():
        raise ValueError(f'Frame {frame_index}: {int(split.sum())} bonded/peptide links cross the periodic boundary; make the protein whole before alignment')

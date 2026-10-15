"""Cavity identities across trajectory frames: births, deaths, splits and merges.

``crevice track-cavities`` (and :func:`track_cavities`) measures every cavity
component of each aligned trajectory frame on one fixed lattice, optionally
restricted to a neighbourhood (for example a region's landmark anchor), and
links components of consecutive measured frames by voxel overlap. Each linked
chain is a persistent *identity*; the topology changes between frames are
recorded as events.

Definitions (see :doc:`/methods/cavity-tracking`):

* **Frames and lattice.** Each frame is rigidly aligned to the first frame
  (C-alpha fit, :func:`crevice.trajectory.align_frame_report`), then cast with
  :func:`crevice.rolling.rolling_probe_cast` (``selection_mode="all"``: every
  inner-probe component inside the rolling bulk-probe envelope) on a lattice
  with world axes, spacing ``spacing`` and a node at the origin, so voxels of
  different frames coincide exactly.
* **Components.** Cast components with swept volume >= ``min_volume``. With a
  ``centre`` only components having a voxel within ``radius`` of it are kept.
  A component's voxel set is its swept sample set (the cast's measured
  points); its volume is the cast volume, its centroid the voxel mean, and
  its lining residues those with a heavy-atom centre within ``lining_cutoff``
  (default 3.3 Å) of one of its voxels (the cast lining definition of
  ``crevice.cast_segments.LINING_DEFINITION``).
* **Links.** Component *A* of the previous measured frame and *B* of the
  current frame are linked when their voxel sets overlap with Jaccard index
  ``|A & B| / |A | B| >= link_jaccard`` (default 0.2) or containment
  ``|A & B| / min(|A|, |B|) >= link_containment`` (default 0.5; the smaller
  set lies mostly inside the larger, as for a small piece budding off).
* **Continuation.** Links are taken in order of decreasing Jaccard index
  (ties: larger overlap, then component order); a link continues an identity
  when neither end has been used yet (one-to-one, greedy).
* **Events.** For each previous component, no link: ``death``; links to two
  or more current components: ``split`` (parent identity, child identities;
  the continued child keeps the parent's identity). For each current
  component, no link: ``birth`` (new identity); links to two or more previous
  components: ``merge`` (parents, child; the continued parent's identity
  survives and the other parents end by merging). A current component linked
  only to already-continued components gets a new identity started by the
  split. Frames whose cast fails are recorded as unresolved and skipped; the
  next measured frame is compared with the last measured one (``gap`` noted).

Memory stays bounded: only the previous frame's voxel sets are kept, and the
trajectory is read in chunks. Identities are geometric continuity under the
stated probes, lattice and thresholds; they are not biological states, and an
event can be produced by noise in a small or barely connected component.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

from .radii import RadiusSet, radii_fields, radii_option

LINK_DEFINITION = ("components of consecutive measured frames are linked when Jaccard |A&B|/|A|B| >= link_jaccard "
                   "or containment |A&B|/min(|A|,|B|) >= link_containment, on one fixed aligned lattice")
CONTINUATION_DEFINITION = ("links taken in decreasing Jaccard order (ties: overlap, then order); a link continues "
                           "an identity when neither component has been used yet")
EVENT_DEFINITION = ("birth: current component with no link; death: previous component with no link; split: previous "
                    "component linked to >= 2 current components; merge: current component linked to >= 2 previous "
                    "components")
INTERPRETATION = ("geometric continuity of measured cavity components under the stated probes, lattice and "
                  "thresholds; not a biological state, transport event or proof of connectivity")

_OFF = 1 << 20
_BITS = 21


def voxel_keys(points, spacing):
    """Integer keys of lattice voxels at world positions on a lattice with a node at the origin.

    Parameters
    ----------
    points : array_like, shape (n, 3)
        Voxel centres, Å (multiples of ``spacing``).
    spacing : float
        Lattice spacing, Å.

    Returns
    -------
    numpy.ndarray of int64
        Sorted unique keys.
    """
    import numpy as np
    idx = np.rint(np.asarray(points, dtype=float).reshape(-1, 3) / spacing).astype(np.int64) + _OFF
    return np.unique((idx[:, 0] << (2 * _BITS)) | (idx[:, 1] << _BITS) | idx[:, 2])


@dataclass
class TrackedComponent:
    """One measured component of one frame, as the tracker sees it.

    Attributes
    ----------
    keys : numpy.ndarray of int64
        Sorted unique voxel keys (:func:`voxel_keys`).
    volume : float
        Volume, Å³.
    centroid : tuple of float
        Voxel mean, Å.
    lining : list of str
        Lining residue labels.
    source_id : int
        Component id in the frame's cast.
    """
    keys: object
    volume: float
    centroid: tuple
    lining: list = field(default_factory=list)
    source_id: int = 0


class CavityTracker:
    """Assign persistent identities to components of successive frames.

    Feed frames in time order with :meth:`update`; read ``events``,
    ``observations`` and :meth:`lifetimes`. Only the previous measured
    frame's components are held.

    Parameters
    ----------
    link_jaccard : float, default 0.2
        Jaccard threshold of a link.
    link_containment : float, default 0.5
        Containment threshold of a link (``None`` or 0 disables it).

    Attributes
    ----------
    events : list of dict
        ``event`` (``birth``, ``death``, ``split``, ``merge``), ``frame_index``,
        ``time_ps``, ``identities`` (involved), ``parents``, ``children``,
        ``jaccard`` (link values) and ``gap_frames``.
    observations : list of dict
        One row per (identity, measured frame).

    Examples
    --------
    A blob that splits in two:

    >>> import numpy as np
    >>> from crevice.tracking import CavityTracker, TrackedComponent
    >>> def blob(xs):
    ...     keys = np.unique(np.asarray(xs, dtype=np.int64))
    ...     return TrackedComponent(keys, float(len(keys)), (0.0, 0.0, 0.0))
    >>> tracker = CavityTracker()
    >>> _ = tracker.update(0, 0.0, [blob(range(0, 20))])
    >>> [e["event"] for e in tracker.update(1, 1.0, [blob(range(0, 9)), blob(range(11, 20))])]
    ['split']
    """

    def __init__(self, *, link_jaccard: float = 0.2, link_containment: float | None = 0.5):
        if not 0 < link_jaccard <= 1:
            raise ValueError("link_jaccard must lie in (0, 1]")
        if link_containment is not None and not 0 <= link_containment <= 1:
            raise ValueError("link_containment must lie in [0, 1]")
        self.link_jaccard = float(link_jaccard)
        self.link_containment = float(link_containment or 0.0)
        self.events, self.observations = [], []
        self._previous = None      # list of (identity, TrackedComponent)
        self._previous_frame = None
        self._next_id = 1
        self._ident = {}           # identity -> record
        self._measured = 0

    def _new(self, frame_index, time_ps, reason, parents=()):
        identity = self._next_id
        self._next_id += 1
        self._ident[identity] = {"identity": identity, "first_frame": frame_index, "first_time_ps": time_ps,
                                 "start_reason": reason, "parents": list(parents), "frames": 0,
                                 "last_frame": frame_index, "last_time_ps": time_ps, "end_reason": "present_at_end",
                                 "volumes": [], "merged_into": None}
        return identity

    def _links(self, current):
        import numpy as np
        prev = self._previous
        if not prev or not current:
            return []
        keys = np.concatenate([c.keys for _, c in prev])
        labels = np.concatenate([np.full(len(c.keys), i) for i, (_, c) in enumerate(prev)])
        order = np.argsort(keys, kind="stable")
        keys, labels = keys[order], labels[order]
        sizes = [len(c.keys) for _, c in prev]
        links = []
        for j, comp in enumerate(current):
            pos = np.searchsorted(keys, comp.keys)
            pos = np.minimum(pos, len(keys) - 1)
            hit = keys[pos] == comp.keys
            if not hit.any():
                continue
            counts = np.bincount(labels[pos[hit]], minlength=len(prev))
            for i in np.flatnonzero(counts):
                inter = int(counts[i])
                union = sizes[i] + len(comp.keys) - inter
                jac = inter / union
                cont = inter / min(sizes[i], len(comp.keys))
                if jac >= self.link_jaccard or (self.link_containment and cont >= self.link_containment):
                    links.append((int(i), j, jac, inter, cont))
        return links

    def update(self, frame_index, time_ps, components):
        """Add one measured frame; return the events it produced.

        Parameters
        ----------
        frame_index : int
            Source frame number.
        time_ps : float
            Frame time, ps.
        components : sequence of TrackedComponent
            Components of this frame.

        Returns
        -------
        list of dict
            Events of this step (see ``events``).
        """
        current = list(components)
        links = self._links(current)
        gap = 0 if self._previous_frame is None else None
        new_events = []
        assigned = [None] * len(current)
        if self._previous is None:
            for j, _ in enumerate(current):
                assigned[j] = self._new(frame_index, time_ps, "first_frame")
        else:
            prev_ids = [identity for identity, _ in self._previous]
            used_prev, used_cur = set(), set()
            for i, j, jac, inter, _ in sorted(links, key=lambda x: (-x[2], -x[3], x[0], x[1])):
                if i in used_prev or j in used_cur:
                    continue
                used_prev.add(i)
                used_cur.add(j)
                assigned[j] = prev_ids[i]
            by_prev, by_cur = {}, {}
            for i, j, jac, inter, cont in links:
                by_prev.setdefault(i, []).append((j, jac))
                by_cur.setdefault(j, []).append((i, jac))
            for j in range(len(current)):
                if assigned[j] is None:
                    parents = [prev_ids[i] for i, _ in by_cur.get(j, [])]
                    assigned[j] = self._new(frame_index, time_ps, "split" if parents else "birth", parents)
            for j in range(len(current)):
                if j not in by_cur:
                    new_events.append({"event": "birth", "identities": [assigned[j]], "parents": [],
                                       "children": [assigned[j]], "jaccard": []})
            for i, identity in enumerate(prev_ids):
                targets = by_prev.get(i, [])
                if not targets:
                    new_events.append({"event": "death", "identities": [identity], "parents": [identity],
                                       "children": [], "jaccard": []})
                    self._ident[identity]["end_reason"] = "death"
                elif len(targets) >= 2:
                    kids = [assigned[j] for j, _ in targets]
                    new_events.append({"event": "split", "identities": sorted({identity, *kids}),
                                       "parents": [identity], "children": kids,
                                       "jaccard": [round(v, 4) for _, v in targets]})
                if targets and i not in used_prev:
                    rec = self._ident[identity]
                    rec["end_reason"] = "merge"
                    rec["merged_into"] = assigned[max(targets, key=lambda t: t[1])[0]]
            for j, sources in by_cur.items():
                if len(sources) >= 2:
                    parents = [prev_ids[i] for i, _ in sources]
                    new_events.append({"event": "merge", "identities": sorted({assigned[j], *parents}),
                                       "parents": parents, "children": [assigned[j]],
                                       "jaccard": [round(v, 4) for _, v in sources]})
            gap = frame_index - self._previous_frame - 1 if frame_index is not None else 0
        for event in new_events:
            event.update(frame_index=frame_index, time_ps=time_ps, gap_frames=gap)
        self.events.extend(new_events)
        for j, comp in enumerate(current):
            rec = self._ident[assigned[j]]
            rec["frames"] += 1
            rec["last_frame"], rec["last_time_ps"] = frame_index, time_ps
            rec["volumes"].append(comp.volume)
            self.observations.append({"identity": assigned[j], "frame_index": frame_index, "time_ps": time_ps,
                                      "volume_A3": comp.volume, "voxels": int(len(comp.keys)),
                                      "centroid_x_A": comp.centroid[0], "centroid_y_A": comp.centroid[1],
                                      "centroid_z_A": comp.centroid[2], "source_component": comp.source_id,
                                      "lining_residues": ";".join(comp.lining)})
        self._previous = list(zip(assigned, current))
        self._previous_frame = frame_index
        self._measured += 1
        return new_events

    def lifetimes(self) -> list[dict]:
        """Lifetime summary, one row per identity.

        Returns
        -------
        list of dict
            ``identity``, ``first_frame``, ``last_frame``, ``first_time_ps``,
            ``last_time_ps``, ``lifetime_ps`` (last minus first time),
            ``frames_present``, ``start_reason`` (``first_frame``, ``birth``,
            ``split``), ``parents``, ``end_reason`` (``present_at_end``,
            ``death``, ``merge``), ``merged_into``, ``mean_volume_A3`` and
            ``max_volume_A3``.
        """
        rows = []
        for rec in self._ident.values():
            v = rec["volumes"]
            rows.append({"identity": rec["identity"], "first_frame": rec["first_frame"], "last_frame": rec["last_frame"],
                         "first_time_ps": rec["first_time_ps"], "last_time_ps": rec["last_time_ps"],
                         "lifetime_ps": rec["last_time_ps"] - rec["first_time_ps"], "frames_present": rec["frames"],
                         "start_reason": rec["start_reason"], "parents": ";".join(map(str, rec["parents"])),
                         "end_reason": rec["end_reason"], "merged_into": rec["merged_into"] or "",
                         "mean_volume_A3": round(sum(v) / len(v), 3) if v else "",
                         "max_volume_A3": round(max(v), 3) if v else ""})
        return rows


# --------------------------------------------------------------------------
# Per-frame measurement
# --------------------------------------------------------------------------

def _crop_hetero(frame, margin):
    """Keep non-hetero atoms and hetero atoms within ``margin`` of the protein bounding box."""
    from dataclasses import replace
    import numpy as np
    protein = [a for a in frame.atoms if not a.hetero]
    if not protein:
        return frame
    xyz = np.asarray([a.coord for a in protein])
    low, high = xyz.min(axis=0) - margin, xyz.max(axis=0) + margin
    keep = tuple(a for a in frame.atoms
                 if not a.hetero or bool(np.all((np.asarray(a.coord) >= low) & (np.asarray(a.coord) <= high))))
    return replace(frame, atoms=keep)


def measure_frame(frame, reference, *, spacing=0.5, probe_radius=1.4, outer_radius=6.0, min_volume=5.0,
                  centre=None, radius=8.0, alignment_residues=None, include_hetero=False, lining_cutoff=3.3,
                  max_grid_points=16_000_000):
    """Align one frame and return its tracked components.

    Parameters
    ----------
    frame : StructureFrame
        Frame to measure (original coordinates).
    reference : StructureFrame
        Alignment reference (first frame).
    spacing : float, default 0.5
        Lattice spacing, Å.
    probe_radius : float, default 1.4
        Inner (cavity) probe radius, Å.
    outer_radius : float, default 6.0
        Rolling bulk-probe radius, Å.
    min_volume : float, default 5.0
        Smallest component volume kept, Å³.
    centre : Coord, optional
        Keep only components with a voxel within ``radius`` of this point
        (aligned coordinates, Å).
    radius : float, default 8.0
        Neighbourhood radius, Å.
    alignment_residues : sequence of str, optional
        Residues for the rigid fit; default all matching C-alpha atoms.
    include_hetero : bool, default False
        Hetero atoms within ``outer_radius`` + 4 Å of the protein box are
        obstacles (for MD input: lipids and ions).
    lining_cutoff : float, default 3.3
        Atom-centre to voxel distance of a lining residue, Å.
    max_grid_points : int, default 16000000
        Grid guard of the cast.

    Returns
    -------
    dict
        ``frame_index``, ``status`` (``measured``), ``components`` (list of
        :class:`TrackedComponent`), ``component_count_all`` (before the
        neighbourhood filter) and ``alignment`` (fit report).
    """
    import numpy as np
    from scipy.spatial import cKDTree
    from .rolling import rolling_probe_cast
    from .trajectory import align_frame_report
    aligned, fit = align_frame_report(frame, reference, residues=alignment_residues)
    if include_hetero:
        aligned = _crop_hetero(aligned, outer_radius + 4.0)
    cast = rolling_probe_cast(aligned, spacing=spacing, probe_radius=probe_radius, outer_radius=outer_radius,
                              min_component_volume=min_volume, selection_mode="all", grid_basis=np.eye(3),
                              grid_anchor=(0.0, 0.0, 0.0), include_hetero=include_hetero,
                              max_grid_points=max_grid_points)
    atoms = aligned.selected_atoms(include_hydrogen=False, include_hetero=include_hetero)
    xyz = np.asarray([a.coord for a in atoms])
    tree = cKDTree(xyz)
    out = []
    for comp in cast.components:
        pts = np.asarray([p.position for p in comp.points])
        if centre is not None and np.min(np.linalg.norm(pts - np.asarray(centre), axis=1)) > radius:
            continue
        hits = cKDTree(pts).query_ball_tree(tree, lining_cutoff)
        near = sorted({k for row in hits for k in row})
        lining = list(dict.fromkeys(atoms[k].residue_key.label for k in near))
        out.append(TrackedComponent(voxel_keys(pts, spacing), float(comp.volume), tuple(map(float, pts.mean(axis=0))),
                                    lining, int(comp.id)))
    return {"frame_index": frame.frame_index, "status": "measured", "components": out,
            "component_count_all": len(cast.components),
            "alignment": {"rmsd_after_A": fit.get("rmsd_after_A"), "atom_count": fit.get("atom_count")}}


_WORKER = None


def _init_worker(frames, reference, settings, radii=None):
    global _WORKER
    _WORKER = (frames, reference, settings)
    if radii is not None:
        from .radii import _ACTIVE_RADII
        _ACTIVE_RADII.set(radii)


def _measure_worker(index):
    frames, reference, settings = _WORKER
    return index, _safe_measure(frames[index], reference, settings)


def _safe_measure(frame, reference, settings):
    try:
        return measure_frame(frame, reference, **settings)
    except ValueError as exc:
        return {"frame_index": frame.frame_index, "status": "unresolved", "reason": str(exc), "components": [],
                "component_count_all": None, "alignment": None}


@radii_option
def track_cavities(topology, trajectory=None, *, selection="protein", start=0, stop=None, stride=1, chunk=50,
                   spacing=0.5, probe_radius=1.4, outer_radius=6.0, min_volume=5.0, centre=None, radius=8.0,
                   region=None, alignment_residues=None, include_hetero=False, lining_cutoff=3.3,
                   link_jaccard=0.2, link_containment=0.5, workers=1, pbc="check", max_grid_points=16_000_000,
                   progress_callback=None, radii: RadiusSet | str | None = None):
    """Track cavity components through a trajectory (bounded memory).

    Parameters
    ----------
    topology : str or Path
        Topology/structure file (GRO, PDB, ...).
    trajectory : str or Path, optional
        Trajectory file; omit for the topology's frames.
    selection : str, default "protein"
        MDAnalysis selection (for lipids as obstacles, for example
        ``"not resname TIP3"`` with ``include_hetero=True``).
    start, stop, stride : int
        Frame range (zero-based, ``stop`` exclusive, default to the end).
    chunk : int, default 50
        Frames read per chunk; only one chunk is held at a time.
    spacing, probe_radius, outer_radius, min_volume : float
        Cast settings (:func:`measure_frame`).
    centre : Coord, optional
        Neighbourhood centre in the first frame's coordinates, Å.
    radius : float, default 8.0
        Neighbourhood radius, Å.
    region : str, Path or dict, optional
        Region definition: its landmark anchor in the first frame becomes
        ``centre`` (:func:`crevice.paths.region_seed`).
    alignment_residues : sequence of str, optional
        Residues of the rigid fit.
    include_hetero : bool, default False
        Hetero atoms near the protein are obstacles.
    lining_cutoff : float, default 3.3
        Lining distance, Å.
    link_jaccard, link_containment : float
        Link thresholds (:class:`CavityTracker`).
    workers : int, default 1
        Worker processes per chunk (fork where available). Results do not
        depend on it.
    pbc : {"check", "none"}, default "check"
        Periodic-boundary bond check of the reader.
    max_grid_points : int, default 16000000
        Cast grid guard.
    progress_callback : callable, optional
        Called as ``progress_callback(done, frame_record)``.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call (including every worker process) comes from that set,
        and results with a ``metadata`` dict record it as
        ``metadata["radii"]``; see :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict
        ``frames`` (per-frame status rows), ``events``, ``observations``,
        ``lifetimes``, ``settings`` (with definitions and the radius set) and
        ``interpretation``.
    """
    import numpy as np
    from .trajectory import load_trajectory_report, inspect_trajectory
    if workers < 1 or chunk < 1 or stride < 1:
        raise ValueError("workers, chunk and stride must be positive integers")
    first, _ = load_trajectory_report(topology, trajectory, selection=selection, start=start, stop=start + 1,
                                      pbc=pbc)
    reference = first.frames[0]
    region_residues = None
    if region is not None:
        from .paths import region_seed
        centre, region_residues = region_seed(reference, region)
    total = None
    if stop is None:
        info = inspect_trajectory(topology, trajectory, selection=selection)
        total = int(info["frame_count"])
        stop = total
    settings = dict(spacing=spacing, probe_radius=probe_radius, outer_radius=outer_radius, min_volume=min_volume,
                    centre=None if centre is None else tuple(map(float, centre)), radius=radius,
                    alignment_residues=alignment_residues, include_hetero=include_hetero,
                    lining_cutoff=lining_cutoff, max_grid_points=max_grid_points)
    tracker = CavityTracker(link_jaccard=link_jaccard, link_containment=link_containment)
    frame_rows, done = [], 0
    from .radii import active_radii
    for chunk_start in range(start, stop, chunk * stride):
        chunk_stop = min(stop, chunk_start + chunk * stride)
        loaded, _ = load_trajectory_report(topology, trajectory, selection=selection, start=chunk_start,
                                           stop=chunk_stop, stride=stride, pbc=pbc, max_frames=None)
        frames, times = loaded.frames, loaded.times()
        if workers == 1:
            results = [_safe_measure(f, reference, settings) for f in frames]
        else:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor
            context = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
            with ProcessPoolExecutor(max_workers=int(workers), mp_context=context, initializer=_init_worker,
                                     initargs=(frames, reference, settings, active_radii())) as pool:
                results = [None] * len(frames)
                for i, rec in pool.map(_measure_worker, range(len(frames))):
                    results[i] = rec
        for frame, time_ps, rec in zip(frames, times, results):
            row = {"frame_index": frame.frame_index, "time_ps": float(time_ps), "status": rec["status"],
                   "components": len(rec["components"]), "components_all": rec["component_count_all"],
                   "total_volume_A3": round(sum(c.volume for c in rec["components"]), 3),
                   "alignment_rmsd_A": (rec["alignment"] or {}).get("rmsd_after_A"), "reason": rec.get("reason", "")}
            if rec["status"] == "measured":
                events = tracker.update(frame.frame_index, float(time_ps), rec["components"])
                row["events"] = ";".join(e["event"] for e in events)
            else:
                row["events"] = ""
            frame_rows.append(row)
            done += 1
            if progress_callback:
                progress_callback(done, row)
        del loaded, frames, results
    return {"frames": frame_rows, "events": tracker.events, "observations": tracker.observations,
            "lifetimes": tracker.lifetimes(),
            "settings": {**{k: v for k, v in settings.items()}, "selection": selection, "start": start, "stop": stop,
                         "stride": stride, "chunk": chunk, "link_jaccard": link_jaccard,
                         "link_containment": link_containment, "workers": workers, "pbc": pbc,
                         "region_landmarks": region_residues, "reference_frame": reference.frame_index,
                         "definitions": {"link": LINK_DEFINITION, "continuation": CONTINUATION_DEFINITION,
                                         "events": EVENT_DEFINITION,
                                         "cast": "rolling_probe_cast selection_mode='all' on a world-axis lattice "
                                                 "with a node at the origin, after C-alpha alignment to the first frame",
                                         "lining": "residue with a heavy-atom centre within lining_cutoff of a "
                                                   "component voxel"},
                         **radii_fields()},
            "interpretation": INTERPRETATION}


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------

#: Okabe-Ito colours for the most persistent identities (rank order); others grey.
IDENTITY_COLORS = ("#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9")
_OTHER = "#b8bec4"
_EVENT_MARKERS = {"birth": ("^", "Birth"), "death": ("x", "Death"), "split": ("D", "Split"), "merge": ("o", "Merge")}


def _ranked_identities(result):
    """Identities ordered by persistence (frames present x mean volume), most persistent first."""
    rows = result["lifetimes"]
    return [r["identity"] for r in sorted(rows, key=lambda r: (-(r["frames_present"] * float(r["mean_volume_A3"] or 0)),
                                                               r["identity"]))]


def plot_tracking(result, path, *, dpi: int = 200, max_lanes: int = 30):
    """Identity volumes over time (top) and identity lifetimes with events (bottom).

    Parameters
    ----------
    result : dict
        From :func:`track_cavities` (or the same keys).
    path : str or Path
        Output image.
    dpi : int, default 200
        Resolution.
    max_lanes : int, default 30
        Largest number of identity lanes drawn in the lower panel (the most
        persistent; the rest are counted in the figure metadata).

    Returns
    -------
    pathlib.Path

    Notes
    -----
    7.0 x 5.6 inch, shared time axis (ns). Top: volume (Å³, log scale) of each identity,
    lines joining consecutive measured frames; the six most persistent
    identities (frames present x mean volume) in Okabe-Ito colours
    (``IDENTITY_COLORS``), all others thin grey. Bottom: one lane per identity
    (most persistent at the top), a line from first to last frame with a dot
    at every frame present, same colours; event markers in dark ink: birth
    triangle, death cross, split diamond (on the parent lane), merge circle
    (on the surviving lane). Unresolved frames are grey ticks on the time axis.
    With more than 100 measured frames the event markers are smaller, and
    with more than 300 identities only the 300 most persistent are drawn in
    the top panel (recorded in the metadata).
    No titles or value labels unless annotating.
    """
    from pathlib import Path
    from matplotlib.lines import Line2D
    from .figures import _pyplot, _save_figure, _style_axes
    from .presentation import FigureText
    plt = _pyplot()
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    text = FigureText("Cavity identities tracked across trajectory frames by voxel overlap: volume of each identity "
                      "against time (top) and identity lifetimes with birth, death, split and merge events (bottom)")
    ranked = _ranked_identities(result)
    color = {ident: (IDENTITY_COLORS[k] if k < len(IDENTITY_COLORS) else _OTHER) for k, ident in enumerate(ranked)}
    series = {}
    for row in result["observations"]:
        series.setdefault(row["identity"], []).append((row["frame_index"], row["time_ps"] / 1000.0, row["volume_A3"]))
    measured = [r["frame_index"] for r in result["frames"] if r["status"] == "measured"]
    order = {f: k for k, f in enumerate(measured)}
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(7.0, 5.6), sharex=True, constrained_layout=True,
                                      gridspec_kw={"height_ratios": [1.0, 1.25]})
    dense = len(measured) > 100
    shown = ranked if len(ranked) <= 300 else ranked[:300]
    if len(shown) < len(ranked):
        text.note(f"Top panel: the {len(shown)} most persistent of {len(ranked)} identities drawn")
    for ident in reversed(shown):
        pts = series.get(ident, [])
        main = color[ident] != _OTHER
        run = []
        for k, p in enumerate(pts):
            if run and order.get(p[0], -9) != order.get(run[-1][0], -9) + 1:
                top.plot([q[1] for q in run], [q[2] for q in run], color=color[ident], linewidth=1.6 if main else 0.8,
                         zorder=3 if main else 2)
                run = []
            run.append(p)
        if run:
            top.plot([q[1] for q in run], [q[2] for q in run], color=color[ident], linewidth=1.6 if main else 0.8,
                     zorder=3 if main else 2, marker="o" if len(run) == 1 else None, markersize=3)
    lanes = ranked[:max_lanes]
    lane = {ident: len(lanes) - 1 - k for k, ident in enumerate(lanes)}
    for ident in lanes:
        pts = series.get(ident, [])
        if not pts:
            continue
        y = lane[ident]
        bottom.plot([pts[0][1], pts[-1][1]], [y, y], color=color[ident], linewidth=1.2, zorder=2)
        bottom.scatter([p[1] for p in pts], [y] * len(pts), s=9, color=color[ident], zorder=3, linewidths=0)
    used = set()
    for event in result["events"]:
        marker, name = _EVENT_MARKERS[event["event"]]
        on = event["parents"] if event["event"] in ("split", "death") else event["children"]
        for ident in on:
            if ident in lane:
                style = (dict(facecolor="white", edgecolor="#222222") if marker in "^Do" else dict(color="#222222"))
                bottom.scatter([event["time_ps"] / 1000.0], [lane[ident]], marker=marker, s=12 if dense else 34,
                               linewidth=0.6 if dense else 1.0, zorder=5, **style)
                used.add(event["event"])
    for row in result["frames"]:
        if row["status"] != "measured":
            bottom.axvline(row["time_ps"] / 1000.0, color="#9aa0a6", linewidth=0.6, linestyle=":", zorder=1)
    handles = [Line2D([], [], marker=_EVENT_MARKERS[e][0], linestyle="", markerfacecolor="white" if
                      _EVENT_MARKERS[e][0] in "^Do" else "#222222", markeredgecolor="#222222", markersize=6)
               for e in ("birth", "death", "split", "merge") if e in used]
    labels = [_EVENT_MARKERS[e][1] for e in ("birth", "death", "split", "merge") if e in used]
    if handles:
        bottom.legend(handles, labels, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    if not result["observations"]:
        text.notice(top, 0.5, 0.5, "No component measured", transform=top.transAxes, ha="center", va="center")
    bottom.set_yticks([])
    bottom.set_ylim(-0.8, max(len(lanes), 1) - 0.2)
    if result["observations"]:
        top.set_yscale("log")
    top.set_ylabel("Volume (Å³, log scale)")
    bottom.set_ylabel("Identities")
    bottom.set_xlabel("Time (ns)")
    top.grid(True, color="#e3e3e3", linewidth=0.6)
    bottom.grid(True, axis="x", color="#e3e3e3", linewidth=0.6)
    for ax in (top, bottom):
        _style_axes(ax)
    bottom.spines["left"].set_visible(False)
    if len(ranked) > len(lanes):
        text.note(f"{len(ranked) - len(lanes)} less persistent identities not drawn as lanes")
    counts = {e: sum(1 for x in result["events"] if x["event"] == e) for e in _EVENT_MARKERS}
    text.note("Events: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    text.note(f"Identities: {len(ranked)}; coloured: {', '.join(map(str, ranked[:len(IDENTITY_COLORS)]))}")
    text.title(top, "CREVICE cavity tracking")
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


def write_tracking_outputs(result, directory, *, prefix: str = "crevice", dpi: int = 200, figure: bool = True):
    """Write the tracking tables, JSON summary and figure.

    Parameters
    ----------
    result : dict
        From :func:`track_cavities`.
    directory : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File prefix.
    dpi : int, default 200
        Figure resolution.
    figure : bool, default True
        Write ``<prefix>_tracking.png``.

    Returns
    -------
    list of pathlib.Path

    Outputs
    -------
    <prefix>_tracking_events.csv
        ``event``, ``frame_index``, ``time_ps``, ``gap_frames``, ``parents``,
        ``children``, ``jaccard`` (``;``-separated).
    <prefix>_tracking_identities.csv
        One row per identity and measured frame: ``identity``, ``frame_index``,
        ``time_ps``, ``volume_A3``, ``voxels``, centroid ``*_A``,
        ``source_component``, ``lining_residues``.
    <prefix>_tracking_lifetimes.csv
        :meth:`CavityTracker.lifetimes`.
    <prefix>_tracking_frames.csv
        Per-frame status, component counts, total volume and events.
    <prefix>_tracking.json
        Settings, definitions, counts and interpretation.
    <prefix>_tracking.png
        :func:`plot_tracking`.
    """
    import csv
    import json
    from pathlib import Path
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    written = []

    def table(name, rows, columns):
        path = directory / f"{prefix}_{name}.csv"
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow({c: row.get(c, "") for c in columns})
        written.append(path)

    events = [{**e, "parents": ";".join(map(str, e["parents"])), "children": ";".join(map(str, e["children"])),
               "jaccard": ";".join(map(str, e["jaccard"]))} for e in result["events"]]
    table("tracking_events", events, ["event", "frame_index", "time_ps", "gap_frames", "parents", "children", "jaccard"])
    obs = [{**o, **{k: round(o[k], 3) for k in ("centroid_x_A", "centroid_y_A", "centroid_z_A", "volume_A3")}}
           for o in result["observations"]]
    table("tracking_identities", obs, ["identity", "frame_index", "time_ps", "volume_A3", "voxels", "centroid_x_A",
                                        "centroid_y_A", "centroid_z_A", "source_component", "lining_residues"])
    table("tracking_lifetimes", result["lifetimes"],
          ["identity", "first_frame", "last_frame", "first_time_ps", "last_time_ps", "lifetime_ps", "frames_present",
           "start_reason", "parents", "end_reason", "merged_into", "mean_volume_A3", "max_volume_A3"])
    table("tracking_frames", result["frames"], ["frame_index", "time_ps", "status", "components", "components_all",
                                                "total_volume_A3", "events", "alignment_rmsd_A", "reason"])
    counts = {e: sum(1 for x in result["events"] if x["event"] == e) for e in _EVENT_MARKERS}
    summary = {"settings": result["settings"], "interpretation": result["interpretation"],
               "frames": len(result["frames"]),
               "measured_frames": sum(1 for r in result["frames"] if r["status"] == "measured"),
               "identities": len(result["lifetimes"]), "event_counts": counts,
               "most_persistent": [r for r in sorted(result["lifetimes"], key=lambda r: -r["frames_present"])[:10]]}
    path = directory / f"{prefix}_tracking.json"
    path.write_text(json.dumps(summary, indent=1, default=str) + "\n")
    written.append(path)
    if figure:
        written.append(plot_tracking(result, directory / f"{prefix}_tracking.png", dpi=dpi))
    return written


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------

def add_tracking_command(subparsers):
    """Register ``crevice track-cavities`` on the CLI subparsers.

    Parameters
    ----------
    subparsers : object
        The ``crevice`` subcommand registry (the object ``add_subparsers`` returns).

    Returns
    -------
    argparse.ArgumentParser
        The parser (handler :func:`run_tracking_command`).
    """
    from .geometry import parse_coord
    parser = subparsers.add_parser(
        "track-cavities", help="Track cavity components across trajectory frames: identities, births, deaths, "
                               "splits and merges")
    parser.add_argument("topology", help="Topology or structure file (GRO, PDB, ...)")
    parser.add_argument("trajectory", nargs="?", help="Trajectory file (XTC, DCD, ...); omit for the topology frames")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", default=None, help="File prefix (default: trajectory or topology stem)")
    parser.add_argument("--md-selection", default="protein",
                        help="MDAnalysis selection (default protein; use 'not resname TIP3' with --include-hetero "
                             "for lipids as obstacles)")
    parser.add_argument("--include-hetero", action="store_true",
                        help="Hetero atoms (lipids, ions) near the protein are obstacles")
    parser.add_argument("--start", type=int, default=0, help="First frame (zero-based)")
    parser.add_argument("--stop", type=int, help="One past the last frame (default: end)")
    parser.add_argument("--stride", type=int, default=1, help="Frame stride")
    parser.add_argument("--chunk", type=int, default=50, help="Frames read per chunk (bounds memory; default 50)")
    where = parser.add_mutually_exclusive_group()
    where.add_argument("--centre", type=parse_coord,
                       help="Keep components with a voxel within --radius of this point x,y,z (first-frame coordinates, Å)")
    where.add_argument("--region", help="Region definition JSON: its landmark anchor in the first frame is the centre")
    parser.add_argument("--radius", type=float, default=8.0, help="Neighbourhood radius in Å (default 8)")
    parser.add_argument("--spacing", type=float, default=0.5, help="Lattice spacing in Å (default 0.5)")
    parser.add_argument("--probe-radius", type=float, default=1.4, help="Cavity probe radius in Å (default 1.4)")
    parser.add_argument("--outer-radius", type=float, default=6.0, help="Bulk-probe radius in Å (default 6)")
    parser.add_argument("--min-volume", type=float, default=5.0, help="Smallest component volume kept, Å³ (default 5)")
    parser.add_argument("--lining-cutoff", type=float, default=3.3,
                        help="Atom-centre to voxel distance of a lining residue, Å (default 3.3)")
    parser.add_argument("--link-jaccard", type=float, default=0.2,
                        help="Jaccard overlap that links components of consecutive frames (default 0.2)")
    parser.add_argument("--link-containment", type=float, default=0.5,
                        help="Containment |A&B|/min(|A|,|B|) that also links them (default 0.5; 0 disables)")
    parser.add_argument("--workers", type=int, default=1, help="Worker processes (default 1)")
    parser.add_argument("--pbc", choices=["check", "none"], default="check",
                        help="Periodic-boundary bond check (none for wrapped lipids)")
    parser.add_argument("--max-grid-points", type=int, default=16_000_000, help="Cast grid guard")
    parser.add_argument("--dpi", type=int, default=200,
                        help="PNG figure resolution in dots per inch (default 200; use --dpi 400 or more for print/publication)")
    parser.set_defaults(func=run_tracking_command)
    return parser


def run_tracking_command(args):
    """Run ``crevice track-cavities`` from parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed ``track-cavities`` arguments.
    """
    from pathlib import Path
    result = track_cavities(args.topology, args.trajectory, selection=args.md_selection, start=args.start,
                            stop=args.stop, stride=args.stride, chunk=args.chunk, spacing=args.spacing,
                            probe_radius=args.probe_radius, outer_radius=args.outer_radius,
                            min_volume=args.min_volume, centre=args.centre, radius=args.radius, region=args.region,
                            include_hetero=args.include_hetero, lining_cutoff=args.lining_cutoff,
                            link_jaccard=args.link_jaccard, link_containment=args.link_containment,
                            workers=args.workers, pbc=args.pbc, max_grid_points=args.max_grid_points)
    prefix = args.prefix or Path(args.trajectory or args.topology).stem
    write_tracking_outputs(result, args.output_dir, prefix=prefix, dpi=args.dpi)
    counts = {e: sum(1 for x in result["events"] if x["event"] == e) for e in _EVENT_MARKERS}
    measured = sum(1 for r in result["frames"] if r["status"] == "measured")
    print(f"track-cavities frames={len(result['frames'])} measured={measured} identities={len(result['lifetimes'])} "
          + " ".join(f"{k}={v}" for k, v in counts.items()))

"""Curved and branched pathways from a seed to bulk solvent (``crevice paths``).

A *pathway* here is a curved centre line through the free space of a
structure, from a start (a seed point, the centroid of seed residues, or the
landmark anchor of a region definition) to bulk solvent, with the free-sphere
radius measured along its arc length. Several distinct routes from the same
start form a *pathway tree*: they share their first part and separate at
*junctions*. Unlike ``crevice profile`` no axis is assumed, and unlike planar
sections the radius is the radius of the largest atom-free sphere centred on
the line (as in HOLE, CAVER and MOLE).

Definitions (all distances in Å, atoms as van der Waals spheres of the radius
set in effect; see :doc:`/methods/paths`):

* **Obstacles.** Heavy atoms of the structure: non-HETATM atoms by default;
  with ``include_hetero=True`` also HETATM atoms (for MD input: lipids and
  ions outside the protein selection) whose spheres can reach the grid box.
* **Grid.** A cubic lattice (``spacing``, default 0.5 Å, world axes) over the
  bounding box of the protein (non-hetero) atoms, padded by the largest atom
  radius, ``bulk_radius`` and three lattice steps. Each node stores its
  *clearance*, the exact distance to the nearest atom surface.
* **Bulk.** The rolling-probe definition used everywhere in CREVICE
  (:func:`crevice.rolling.rolling_probe_cast`, lateral exits): bulk-probe
  centres are nodes with clearance >= ``bulk_radius`` (default 6 Å)
  face-connected, along atom-clear segments, to the faces of the grid box; a
  point is in bulk when it lies within ``bulk_radius`` of such a centre. With
  membrane obstacles included, the box faces in the bilayer are blocked, so
  only solvent on either side counts as bulk.
* **Free nodes and steps.** A node is free when its clearance exceeds
  ``probe_radius``. Routes move
  between free nodes in the 26 lattice directions. Diagonal steps are not
  tested analytically during the search (the clearance along a step of at
  most 0.87 Å can fall below its endpoints by at most a few hundredths of an
  Å); every reported centre line is measured exactly afterwards.
* **Widest routes.** A maximin (widest-path) Dijkstra from the start labels
  every free node outside bulk with its widest route; equally wide routes
  prefer the shorter one, then lattice order (deterministic). A route's
  *ranking width* is its smallest node clearance outside the start's own
  free sphere (all routes share that sphere, as in the lateral-exit search).
  Bulk nodes are terminals and are never expanded, so each route ends at the
  first bulk node it meets: the **mouth** (exit point).
* **Distinct routes.** Terminals are visited in order of decreasing width,
  then increasing length. A route is kept when, compared with every route
  already kept, (a) its exit point is more than ``separation`` from the
  other's exit point (default ``2 * bulk_radius``: no single bulk-probe
  position touches both) and (b) the parts of the two routes beyond their
  junction (the last node they share) and farther than ``route_separation``
  (default ``bulk_radius``) from it never come within ``route_separation`` of
  each other. For two routes that share only the start this is exactly the
  lateral-exit rule of :func:`crevice.channel_exits.lateral_exit_legs`. At
  most ``max_paths`` routes are kept (default 8).
* **Sides.** With an ``axis``, the search runs once per half-space, as for
  lateral exits: the ``+`` search may only use nodes on the ``+`` side of the
  plane through the start normal to the axis (plus the start's own free
  sphere), the ``-`` search the other side, each with its own distinct routes
  (``max_paths`` per side). Without the half-space a maximin route may reach
  the far side through a wider detour (a crevice around the protein) rather
  than through a narrower direct passage. A side with no route is recorded as
  unresolved. Without an axis, one search covers all free space.
* **Exposed nodes.** Besides bulk, a route also ends at an *exposed* node:
  a free node from which more than ``26 - ceil(26 * min_enclosure)`` (by
  default at least 7) of the 26 straight lattice rays reach a bulk node
  before meeting a blocked node (clearance <= probe); a ray that leaves the
  grid first does not count.
  This stops routes from running along surface grooves outside the protein
  envelope that the 6 Å bulk probe cannot enter. Each branch records whether
  it ends in ``bulk`` or at an ``exposed`` node; ``min_enclosure=0`` gives
  the pure bulk definition.
* **Tree, junctions and segments.** The kept routes are paths of one
  shortest-path tree, so they share prefixes exactly. A *junction* is a node
  where kept routes separate (the start itself when they separate at once).
  The tree is cut at the start, the junctions and the mouths into
  *segments*; a *branch* is the full route from the start to one mouth.
* **Centre line.** Each segment is resampled every ``spacing`` and relaxed
  with fixed ends: a Laplacian (shortening) step, then a centring step along
  the clearance gradient perpendicular to the line; a move is accepted only
  if the interpolated clearance at the new position is at least the smaller
  of its old clearance and the lattice route's bottleneck, so relaxation can
  never narrow the route below its lattice bottleneck. If the exact
  clearance of the relaxed line falls to ``probe_radius`` or below, the
  segment keeps the straightened lattice route instead (recorded).
* **Radius.** Every 0.25 Å along each branch, the exact distance to the
  nearest atom surface (a free-sphere radius, not a planar section radius,
  not reduced by the probe). The **bottleneck** is the sample of smallest
  radius, reported with its position and arc length.
* **Lining residues.** The lining shell of :mod:`crevice.fingerprint`: an
  atom lines a sample when ``|x - c| - r_atom - R <= lining_shell`` (default
  1.4 Å), with ``c`` the sample centre and ``R`` its radius; hydrogens, water
  and recognised monatomic ions are excluded; a residue lines a branch when
  any of its atoms lines a sample. Its anchor is the arc length of its
  smallest gap.
* **Unresolved.** When the start fits no probe within ``seed_tolerance``,
  lies in bulk, or no free route reaches bulk, nothing is interpolated: the
  result's ``status`` says so and no branch is reported. With an ``axis``,
  each side without a route is recorded in ``unresolved``.

* **Automatic probe** (``probe_radius="auto"``, the default). Every probe of
  ``PROBE_LADDER`` (0.5-3.0 Å, the ladder of ``crevice profile``'s automatic
  enclosure probe) below ``bulk_radius`` is searched on one clearance field,
  bulk and exposure computation (``_exposure_levels``), from one start
  (snapped at the smallest probe that resolves it). Probes at or above the
  widest route width found at the smallest probe are pruned: a route at a
  larger probe is also a route at the smallest one, so none can exist.
  Per side, the route sets found are grouped (same number of routes,
  matched one to one as not distinct with ranking widths within
  ``max(0.05 Å, 2 %)``); the group found at the most probes (at least two;
  ties: the smaller probe) wins and its smallest probe is used
  (``AUTO_PROBE_RULE``). Every candidate, the groups, the choice and its
  reason are recorded in ``settings["probe_selection"]``. ``min_enclosure``
  is the same for every candidate.

Pathways are geometric routes for the stated probes and grid. They are not
transport pathways, permeation mechanisms or functional assignments; the
absence of a route on the grid does not prove that none exists.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import heapq
import math

from .radii import RadiusSet, atom_vdw_radius, radii_fields, radii_option

PATH_DEFINITION = (
    "curved centre line of the widest (maximin node clearance, then shortest) lattice route of free nodes "
    "(clearance > probe_radius) from the start to the first node in bulk, relaxed toward the local clearance "
    "maximum without dropping below the lattice bottleneck")
BULK_DEFINITION = (
    "space within bulk_radius of a bulk-probe centre: a node with clearance >= bulk_radius face-connected along "
    "atom-clear segments to the faces of the padded grid box (rolling-probe envelope, as crevice cast and lateral exits)")
DISTINCT_PATH_DEFINITION = (
    "routes ranked by width beyond the start's free sphere, then length; a further route is distinct when its exit "
    "point is more than separation_A from every kept exit point and the parts of the two routes beyond their "
    "junction and farther than route_separation_A from it never come within route_separation_A of each other")
RADIUS_DEFINITION = ("distance from each centre-line sample to the nearest atom van der Waals surface "
                     "(free-sphere radius, HOLE-like), not reduced by the probe; not a planar section radius")
LINING_DEFINITION = ("atom (no H, water or monatomic ions) whose van der Waals surface lies within lining_shell_A "
                     "of a sample's free sphere: |x - c| - r_atom - R <= lining_shell_A")
INTERPRETATION = ("geometric route for the stated probe, bulk probe and grid; not a transport pathway, permeation "
                  "mechanism or functional assignment")

#: Default largest number of distinct routes.
MAX_PATHS = 8
#: Candidate probes (Å) of ``probe_radius="auto"``: the ladder of ``crevice profile``'s automatic enclosure probe.
PROBE_LADDER = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20,
                1.30, 1.40, 1.50, 1.60, 1.80, 2.00, 2.20, 2.40, 2.70, 3.00)
#: Automatic probe rule (recorded in ``settings["probe_selection"]["rule"]``).
AUTO_PROBE_RULE = (
    "per side: every ladder probe below the bulk probe is searched on one shared clearance field from one "
    "start (snapped at the smallest probe that resolves it); probes at or above the widest route width found "
    "at the smallest probe are pruned (they cannot have a route); the route sets found are grouped (same "
    "number of routes, one-to-one matched as not distinct under the distinct-route rule); the group found at "
    "the most probes (at least 2; ties: the smaller probe) wins and its smallest probe is used")
#: Okabe-Ito colours, one per branch in rank order (colour-blind-safe).
BRANCH_COLORS = ("#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9", "#000000", "#999999")
#: Line styles cycled with the colours (secondary encoding in figures).
BRANCH_STYLES = ("-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 1)), (0, (8, 2)))

_STEPS = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)]


@dataclass
class PathBranch:
    """One route of a pathway tree, from the start to a mouth in bulk.

    Attributes
    ----------
    id : int
        Rank, 1 = widest (ranking width, then shortest).
    points : list of tuple
        Centre-line samples (Å), at most 0.25 Å apart, from the start to the mouth.
    radii : list of float
        Free-sphere radius at each sample (Å, ``RADIUS_DEFINITION``).
    arc_lengths : list of float
        Arc length of each sample from the start (Å).
    length : float
        Branch length (Å).
    bottleneck_radius : float
        Smallest sample radius (Å).
    bottleneck_position : tuple
        Sample of smallest radius (Å).
    bottleneck_arc_length : float
        Its arc length from the start (Å).
    exit_point : tuple
        The mouth: first lattice node in bulk (Å).
    ranking_width : float
        Lattice route width beyond the start's free sphere (Å), used to rank.
    side : str or None
        ``"+"`` or ``"-"``: side of the start the mouth lies on along ``axis``
        (``None`` without an axis).
    end_type : str
        ``"bulk"`` (the mouth is in the rolling-probe bulk) or ``"exposed"``
        (the mouth is an exposed node; see :func:`find_paths` ``min_enclosure``).
    segments : list of int
        Segment ids from the start to the mouth.
    junctions : list of dict
        ``junction`` id and ``arc_length_A`` of every junction the branch passes.
    lining : list of dict
        Lining residues (``residue``, ``min_gap_A``, ``anchor_arc_length_A``,
        ``samples``), ordered by anchor.
    bottleneck_residues : list of str
        Residues lining the bottleneck sample.
    metadata : dict
        ``centre_line`` (``relaxed``, ``lattice_straightened`` or
        ``lattice_fallback`` per segment), ``lattice_route_length_A``,
        ``lattice_nodes`` and ``below_probe`` (True when the measured
        bottleneck is at or below the probe radius: the lattice route passed,
        but the continuous centre line touches the probe, a grid effect).
    """
    id: int
    points: list
    radii: list
    arc_lengths: list
    length: float
    bottleneck_radius: float
    bottleneck_position: tuple
    bottleneck_arc_length: float
    exit_point: tuple
    ranking_width: float
    side: str | None = None
    end_type: str = "bulk"
    segments: list = field(default_factory=list)
    junctions: list = field(default_factory=list)
    lining: list = field(default_factory=list)
    bottleneck_residues: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """JSON record of the branch.

        Returns
        -------
        dict
            Every attribute, with points as lists and ``start_point``.
        """
        return {"branch": self.id, "side": self.side, "length_A": self.length,
                "bottleneck_radius_A": self.bottleneck_radius,
                "bottleneck_position_A": list(self.bottleneck_position),
                "bottleneck_arc_length_A": self.bottleneck_arc_length,
                "start_point_A": list(self.points[0]), "exit_point_A": list(self.exit_point),
                "end_type": self.end_type, "ranking_width_A": self.ranking_width,
                "segments": list(self.segments), "junctions": list(self.junctions),
                "bottleneck_residues": list(self.bottleneck_residues), "lining": list(self.lining),
                "points_A": [list(p) for p in self.points], "radii_A": list(self.radii),
                "arc_lengths_A": list(self.arc_lengths), **self.metadata}


@dataclass
class PathNetwork:
    """Result of :func:`crevice.paths.find_paths`: the pathway tree from one start, or why there is none.

    Attributes
    ----------
    status : str
        ``"resolved"`` (at least one branch), ``"unresolved_seed_blocked"``,
        ``"unresolved_seed_in_bulk"`` or ``"unresolved_no_bulk_path"``.
    start : dict
        Start report: requested point, snapped lattice node, its clearance,
        snap distance and how the start was given.
    branches : list of PathBranch
        Widest first.
    junctions : list of dict
        ``junction`` id, ``position_A``, ``branches`` and ``is_start``.
    segments : list of dict
        ``segment`` id, ``from``/``to`` node kinds and ids, ``points_A`` and
        ``centre_line`` (``relaxed`` or ``lattice_fallback``).
    unresolved : list of dict
        Unresolved cases (for example a side without a route).
    settings : dict
        Every parameter, the grid, the radius set and the definitions.
    """
    status: str
    start: dict
    branches: list
    junctions: list
    segments: list
    unresolved: list
    settings: dict

    def to_dict(self) -> dict:
        """JSON record of the whole result.

        Returns
        -------
        dict
            ``status``, ``start``, ``branches`` (:meth:`PathBranch.to_dict`),
            ``junctions``, ``segments``, ``unresolved``, ``settings`` and
            ``interpretation``.
        """
        return {"status": self.status, "start": self.start,
                "branches": [b.to_dict() for b in self.branches], "junctions": self.junctions,
                "segments": self.segments, "unresolved": self.unresolved, "settings": self.settings,
                "interpretation": INTERPRETATION}

    def summary_rows(self) -> list[dict]:
        """One table row per branch (or one ``unresolved`` row), with units in column names.

        Returns
        -------
        list of dict
            ``branch``, ``status``, ``side``, ``length_A``,
            ``bottleneck_radius_A``, ``bottleneck_arc_length_A``,
            bottleneck, start and exit coordinates, ``ranking_width_A``,
            ``probe_radius_A`` (the probe of the branch's search), ``junction_count``, ``lining_residue_count``,
            ``bottleneck_residues`` and ``lining_residues`` (``;``-separated).
        """
        if not self.branches:
            return [{"branch": "", "status": self.status, "side": "", "length_A": "", "bottleneck_radius_A": "",
                     "bottleneck_arc_length_A": "", "bottleneck_x_A": "", "bottleneck_y_A": "",
                     "bottleneck_z_A": "", "start_x_A": self.start.get("node_A", ["", "", ""])[0] if self.start.get("node_A") else "",
                     "start_y_A": self.start["node_A"][1] if self.start.get("node_A") else "",
                     "start_z_A": self.start["node_A"][2] if self.start.get("node_A") else "",
                     "exit_x_A": "", "exit_y_A": "", "exit_z_A": "", "ranking_width_A": "", "probe_radius_A": "",
                     "junction_count": "",
                     "lining_residue_count": "", "end_type": "", "bottleneck_residues": "", "lining_residues": ""}]
        rows = []
        for b in self.branches:
            rows.append({"branch": b.id, "status": self.status, "side": b.side or "",
                         "length_A": round(b.length, 3), "bottleneck_radius_A": round(b.bottleneck_radius, 3),
                         "bottleneck_arc_length_A": round(b.bottleneck_arc_length, 3),
                         "bottleneck_x_A": round(b.bottleneck_position[0], 3),
                         "bottleneck_y_A": round(b.bottleneck_position[1], 3),
                         "bottleneck_z_A": round(b.bottleneck_position[2], 3),
                         "start_x_A": round(b.points[0][0], 3), "start_y_A": round(b.points[0][1], 3),
                         "start_z_A": round(b.points[0][2], 3),
                         "exit_x_A": round(b.exit_point[0], 3), "exit_y_A": round(b.exit_point[1], 3),
                         "exit_z_A": round(b.exit_point[2], 3), "ranking_width_A": round(b.ranking_width, 3),
                         "probe_radius_A": b.metadata.get("probe_radius_A", self.settings.get("probe_radius_A")),
                         "junction_count": len(b.junctions), "lining_residue_count": len(b.lining),
                         "end_type": b.end_type, "bottleneck_residues": ";".join(b.bottleneck_residues),
                         "lining_residues": ";".join(r["residue"] for r in b.lining)})
        return rows


# --------------------------------------------------------------------------
# Free space: clearance grid and bulk
# --------------------------------------------------------------------------

class FreeSpace:
    """Clearance grid and rolling-probe bulk around a structure (world axes).

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms (van der Waals spheres of the radius set in effect).
    core : numpy.ndarray, shape (n, 3)
        Coordinates whose bounding box (padded) defines the grid box, usually
        the protein atoms.
    spacing : float
        Lattice spacing, Å.
    bulk_radius : float
        Rolling bulk-probe radius, Å.
    max_grid_points : int
        Allocation guard.

    Attributes
    ----------
    origin : numpy.ndarray
        World position of node (0, 0, 0), Å.
    dims : tuple of int
        Grid shape.
    clearance : numpy.ndarray of float32
        Distance from each node to the nearest atom surface, Å.
    in_bulk : numpy.ndarray of bool
        Nodes inside the bulk region.
    atoms : tuple of Atom
        Obstacle atoms that can reach the box.

    Raises
    ------
    ValueError
        If the grid would exceed ``max_grid_points``.
    """

    def __init__(self, atoms, core, spacing, bulk_radius, max_grid_points):
        import numpy as np
        from scipy import ndimage
        from .rolling import _SphereQueries, _connected
        xyz = np.asarray([a.coord for a in atoms], dtype=float)
        radii = np.asarray([atom_vdw_radius(a) for a in atoms], dtype=float)
        rmax = float(radii.max())
        pad = rmax + bulk_radius + 3 * spacing
        core = np.asarray(core, dtype=float)
        low, high = core.min(axis=0) - pad, core.max(axis=0) + pad
        dims = tuple(int(v) for v in np.ceil((high - low) / spacing).astype(int) + 1)
        total = math.prod(dims)
        if total > max_grid_points:
            raise ValueError(f"Path grid needs {total:,} points at {spacing:g} Å; max_grid_points={max_grid_points:,}. "
                             "Use a coarser --spacing or raise the limit")
        reach = rmax + bulk_radius + spacing
        keep = np.all((xyz >= low - reach) & (xyz <= low + (np.asarray(dims) - 1) * spacing + reach), axis=1)
        self.atoms = tuple(a for a, k in zip(atoms, keep) if k)
        xyz, radii = xyz[keep], radii[keep]
        queries = _SphereQueries(xyz, radii)
        clearance = np.empty(total, dtype=np.float64)
        for start in range(0, total, 262144):
            stop = min(start + 262144, total)
            ijk = np.array(np.unravel_index(np.arange(start, stop), dims)).T
            clearance[start:stop] = queries.points(low + ijk * spacing)[0]
        clearance = clearance.reshape(dims)
        labels = _connected(clearance >= bulk_radius, clearance, low, spacing, bulk_radius, queries)
        faces = np.unique(np.concatenate([labels[0].ravel(), labels[-1].ravel(), labels[:, 0].ravel(),
                                          labels[:, -1].ravel(), labels[:, :, 0].ravel(), labels[:, :, -1].ravel()]))
        centres = np.isin(labels, faces[faces > 0])
        del labels
        self.in_bulk = ndimage.distance_transform_edt(~centres, sampling=spacing) <= bulk_radius + 1e-9
        self.clearance = clearance.astype(np.float32)
        self.origin, self.spacing, self.dims, self.bulk_radius = low, float(spacing), dims, float(bulk_radius)
        self.grid_points = total

    def position(self, index):
        """World coordinates (Å) of lattice indices.

        Parameters
        ----------
        index : array_like, shape (..., 3)
            Integer lattice indices.

        Returns
        -------
        numpy.ndarray
        """
        import numpy as np
        return self.origin + np.asarray(index, dtype=float) * self.spacing

    def interpolate(self, points):
        """Trilinear clearance (Å) at world points; outside the box the nearest edge value.

        Parameters
        ----------
        points : array_like, shape (n, 3)

        Returns
        -------
        numpy.ndarray, shape (n,)
        """
        import numpy as np
        from scipy.ndimage import map_coordinates
        coords = ((np.asarray(points, dtype=float) - self.origin) / self.spacing).T
        return map_coordinates(self.clearance, coords, order=1, mode="nearest")

    def gradient(self, points, step=None):
        """Central-difference gradient of the interpolated clearance at world points.

        Parameters
        ----------
        points : array_like, shape (n, 3)
        step : float, optional
            Difference step, Å; default half the spacing.

        Returns
        -------
        numpy.ndarray, shape (n, 3)
        """
        import numpy as np
        points = np.asarray(points, dtype=float)
        h = self.spacing / 2 if step is None else step
        out = np.empty_like(points)
        for axis in range(3):
            delta = np.zeros(3)
            delta[axis] = h
            out[:, axis] = (self.interpolate(points + delta) - self.interpolate(points - delta)) / (2 * h)
        return out


# --------------------------------------------------------------------------
# Start selection
# --------------------------------------------------------------------------

def residue_seed(frame, residues):
    """Centroid (Å) of the heavy atoms of the named residues.

    Parameters
    ----------
    frame : StructureFrame
    residues : str or sequence of str
        Residue tokens, comma-separated in a string: a CREVICE label
        (``"A:GLN282"``), a name and number (``"GLN282"``) or a number
        (``"282"``); an optional ``@ATOM+ATOM`` suffix limits the atoms (for
        example ``"GLN282@OE1+NE2"``).

    Returns
    -------
    centroid : tuple of float
        Mean of the per-residue centroids (each residue counts once).
    matched : list of str
        Labels of the matched residues.

    Raises
    ------
    ValueError
        If a token matches no residue.
    """
    import re
    import numpy as np
    tokens = [t.strip() for t in (residues.split(",") if isinstance(residues, str) else residues) if t.strip()]
    if not tokens:
        raise ValueError("no seed residues given")
    centroids, matched = [], []
    for token in tokens:
        spec, _, names = token.partition("@")
        names = {n.strip().upper() for n in names.split("+") if n.strip()}
        chain, _, rest = spec.rpartition(":")
        found = re.fullmatch(r"([A-Za-z0-9]*?)(-?\d+)([A-Za-z]?)", rest)
        if not found:
            raise ValueError(f"cannot read seed residue {token!r}; use e.g. GLN282, 282, A:GLN282 or GLN282@OE1+NE2")
        resname, resid, icode = found.group(1).upper(), int(found.group(2)), found.group(3)
        groups = {}
        for atom in frame.atoms:
            key = atom.residue_key
            if key.resid != resid or (resname and key.resname.upper() != resname) or (icode and key.icode != icode):
                continue
            if chain and chain != "_" and key.chain_id != chain:
                continue
            if atom.element.upper() in {"H", "D"} or (names and atom.name.upper() not in names):
                continue
            groups.setdefault(key, []).append(atom.coord)
        if not groups:
            raise ValueError(f"seed residue {token!r} matches no heavy atoms")
        for key, coords in sorted(groups.items()):
            centroids.append(np.mean(coords, axis=0))
            matched.append(key.label + (f"@{'+'.join(sorted(names))}" if names else ""))
    return tuple(float(v) for v in np.mean(centroids, axis=0)), matched


def region_seed(frame, definition):
    """Landmark anchor (Å) of a CREVICE region definition in this frame.

    The anchor is the equal-weight mean of the landmark centroids, as in
    :func:`crevice.region_definition.propose_landmark_region`.

    Parameters
    ----------
    frame : StructureFrame
        Frame whose atoms carry the landmark residues.
    definition : dict or str or Path
        Region definition (or its JSON path; only ``landmarks`` is read, no
        file hashes are checked here).

    Returns
    -------
    anchor : tuple of float
    landmarks : list of str
        Residue labels of the landmarks.

    Raises
    ------
    ValueError
        If the definition has no landmarks or one does not resolve.
    """
    import json
    from pathlib import Path
    import numpy as np
    from .region_definition import resolve_landmarks
    data = json.loads(Path(definition).read_text()) if not isinstance(definition, dict) else definition
    rows = data.get("landmarks") or []
    if not rows:
        raise ValueError("region definition has no residue landmarks")
    resolved = resolve_landmarks(frame, rows)
    anchor = np.mean([r["centroid_A"] for r in resolved], axis=0)
    return tuple(float(v) for v in anchor), [r["residue"] for r in resolved]


# --------------------------------------------------------------------------
# Search
# --------------------------------------------------------------------------

def _ray_escapes(blocked, bulk):
    """Number of the 26 lattice rays from each node that meet a bulk node before a blocked node.

    A ray that leaves the grid without meeting either does not escape.
    """
    import numpy as np
    own = np.where(bulk, 1, np.where(blocked, 2, 0)).astype(np.int8)
    count = np.zeros(blocked.shape, dtype=np.uint8)
    for d in _STEPS:
        state = own.copy()
        a = next(i for i, v in enumerate(d) if v)
        sign = d[a]
        n = blocked.shape[a]
        for i in (range(n - 2, -1, -1) if sign == 1 else range(1, n)):
            target, source = [slice(None)] * 3, [slice(None)] * 3
            target[a], source[a] = i, i + sign
            for j in range(3):
                if j != a and d[j]:
                    target[j] = slice(None, -1) if d[j] == 1 else slice(1, None)
                    source[j] = slice(1, None) if d[j] == 1 else slice(None, -1)
            here = state[tuple(target)]
            state[tuple(target)] = np.where(here == 0, state[tuple(source)], here)
        count += state == 1
    return count


def _widest_routes(space, roots, exempt, probe, max_nodes, allowed=None, outside=None):
    """Maximin-then-length Dijkstra from root nodes; returns labels, parents and terminals in pop order.

    ``outside`` (bulk or exposed nodes, the route terminals) defaults to
    ``space.outside`` when set, else ``space.in_bulk``.
    """
    import numpy as np
    from scipy import ndimage
    free = space.clearance > probe
    if outside is None:
        outside = space.in_bulk if getattr(space, "outside", None) is None else space.outside
    inner = free & ~outside
    if allowed is not None:
        inner &= allowed
    labels, _ = ndimage.label(inner, structure=np.ones((3, 3, 3), bool))
    root_labels = {int(labels[tuple(r)]) for r in roots}
    domain = np.isin(labels, list(root_labels)) & inner
    del labels
    terminal_mask = ndimage.binary_dilation(domain, structure=np.ones((3, 3, 3), bool)) & free & outside
    if allowed is not None:
        terminal_mask &= allowed
    active = domain | terminal_mask
    # Work inside the bounding box of the active nodes (plus one layer of padding).
    where = np.argwhere(active)
    lo = np.maximum(where.min(axis=0) - 1, 0)
    hi = np.minimum(where.max(axis=0) + 2, np.asarray(space.dims))
    box = tuple(slice(a, b) for a, b in zip(lo, hi))
    sub_active, sub_terminal = active[box], terminal_mask[box]
    shape = sub_active.shape
    flat = np.flatnonzero(sub_active.ravel())
    if len(flat) > max_nodes:
        raise ValueError(f"path search domain has {len(flat):,} nodes (max_nodes={max_nodes:,}); use a coarser spacing")
    compact = np.full(math.prod(shape), -1, dtype=np.int64)
    compact[flat] = np.arange(len(flat))
    ijk = np.array(np.unravel_index(flat, shape)).T
    n = len(flat)
    neighbours = np.full((n, 26), -1, dtype=np.int64)
    lengths = np.empty(26)
    for s, step in enumerate(_STEPS):
        other = ijk + np.asarray(step)
        ok = np.all((other >= 0) & (other < np.asarray(shape)), axis=1)
        idx = np.full(n, -1, dtype=np.int64)
        idx[ok] = compact[np.ravel_multi_index(other[ok].T, shape)]
        neighbours[:, s] = idx
        lengths[s] = space.spacing * math.sqrt(step[0] ** 2 + step[1] ** 2 + step[2] ** 2)
    world_ijk = ijk + lo
    clear = space.clearance[tuple(world_ijk.T)].astype(float)
    is_terminal = sub_terminal.ravel()[flat]
    exempt_mask = exempt(space.position(world_ijk)) if exempt is not None else np.zeros(n, bool)
    width = np.full(n, -np.inf)
    walked = np.full(n, np.inf)
    parent = np.full(n, -1, dtype=np.int64)
    heap = []
    for r in roots:
        c = compact[np.ravel_multi_index(tuple(np.asarray(r) - lo), shape)]
        if c < 0:
            continue
        width[c], walked[c] = np.inf, 0.0
        heap.append((-np.inf, 0.0, int(c)))
    heapq.heapify(heap)
    done = np.zeros(n, bool)
    terminals = []
    while heap:
        key, dist, node = heapq.heappop(heap)
        if done[node]:
            continue
        done[node] = True
        if is_terminal[node]:
            terminals.append(node)
            continue
        nb = neighbours[node]
        valid = nb >= 0
        nbv = nb[valid]
        cur = width[node]
        new_w = np.where(exempt_mask[nbv], cur, np.minimum(cur, clear[nbv]))
        new_d = dist + lengths[valid]
        better = (~done[nbv]) & ((new_w > width[nbv]) | ((new_w == width[nbv]) & (new_d < walked[nbv] - 1e-9)))
        for v, w, d in zip(nbv[better].tolist(), new_w[better].tolist(), new_d[better].tolist()):
            width[v], walked[v], parent[v] = w, d, node
            heapq.heappush(heap, (-w, d, v))
    return {"ijk": world_ijk, "clear": clear, "width": width, "walked": walked, "parent": parent,
            "terminals": terminals, "domain_nodes": int(domain.sum()), "terminal_nodes": int(terminal_mask.sum()),
            "explored": int(done.sum())}


def _route(search, node):
    out = []
    while node >= 0:
        out.append(int(node))
        node = search["parent"][node]
    return out[::-1]


def _min_pair_distance(a, b):
    import numpy as np
    from scipy.spatial import cKDTree
    if not len(a) or not len(b):
        return math.inf
    return float(cKDTree(np.asarray(b)).query(np.asarray(a))[0].min())


def _resample(points, step):
    """Polyline resampled at (almost) uniform arc steps, keeping both ends."""
    import numpy as np
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return points
    seg = np.linalg.norm(np.diff(points, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    total = arc[-1]
    if total <= 0:
        return points[[0, -1]]
    count = max(1, int(math.ceil(total / step)))
    targets = np.linspace(0.0, total, count + 1)
    return np.column_stack([np.interp(targets, arc, points[:, d]) for d in range(3)])


def _relax(space, points, floor, iterations, step):
    """Laplacian + perpendicular centring relaxation with fixed ends (see module docstring)."""
    import numpy as np
    p = _resample(points, step)
    if len(p) < 3:
        return p
    for it in range(iterations):
        old_c = space.interpolate(p[1:-1])
        lap = p[1:-1] + 0.5 * ((p[:-2] + p[2:]) / 2 - p[1:-1])
        tangent = p[2:] - p[:-2]
        tangent /= np.maximum(np.linalg.norm(tangent, axis=1)[:, None], 1e-12)
        g = space.gradient(lap)
        g -= np.sum(g * tangent, axis=1)[:, None] * tangent
        cand = lap + 0.5 * step * g
        new_c = space.interpolate(cand)
        accept = new_c >= np.minimum(old_c, floor) - 1e-6
        p[1:-1] = np.where(accept[:, None], cand, p[1:-1])
        if it % 10 == 9:
            p = _resample(p, step)
    return _resample(p, step)


def _straighten(spatial, nodes):
    """Straight cuts along a lattice route, each no narrower than the stretch it replaces."""
    edges = [spatial.segment_clearance(tuple(a), tuple(b)) for a, b in zip(nodes, nodes[1:])]
    if not edges:
        return [tuple(n) for n in nodes]
    kept, i = [tuple(nodes[0])], 0
    while i < len(nodes) - 1:
        j = len(nodes) - 1
        while j > i + 1 and spatial.segment_clearance(tuple(nodes[i]), tuple(nodes[j])) < min(edges[i:j]) - 1e-9:
            j -= 1
        kept.append(tuple(nodes[j]))
        i = j
    return kept


def _axis_vector(axis):
    import numpy as np
    if axis is None:
        return None
    if isinstance(axis, str):
        named = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}
        if axis.strip().lower() in named:
            axis = named[axis.strip().lower()]
        else:
            axis = tuple(float(v) for v in axis.split(","))
    v = np.asarray(axis, dtype=float)
    if v.shape != (3,) or not np.isfinite(v).all() or np.linalg.norm(v) == 0:
        raise ValueError("axis must be x, y, z or a finite nonzero x,y,z vector")
    return v / np.linalg.norm(v)


def _exposure_levels(clearance, bulk, probes, required):
    """Exposure of every node for a whole probe ladder at once.

    For a ray, the largest probe at which it escapes is the smallest node
    clearance met before its first bulk node (``inf`` from a bulk node,
    ``-inf`` when it leaves the grid first). A node is exposed at probe
    ``p`` when more than ``26 - required`` rays escape, that is when the
    ``(27 - required)``-th largest of its 26 ray values exceeds ``p``. This
    is :func:`_ray_escapes` for every probe in one sweep per direction.

    Parameters
    ----------
    clearance : numpy.ndarray of float32
        Node clearances (Å).
    bulk : numpy.ndarray of bool
        Bulk nodes.
    probes : sequence of float
        Candidate probes (Å), ascending.
    required : int
        ``ceil(26 * min_enclosure)``, 1 to 26.

    Returns
    -------
    numpy.ndarray of uint8
        Level ``L``: at ``probes[j]`` a free, non-bulk node is exposed exactly
        when ``L > j`` (comparisons in float32, as the single-probe rule).
    """
    import numpy as np
    k = 27 - int(required)
    thresholds = np.asarray(probes, dtype=np.float32)
    clear = np.asarray(clearance, dtype=np.float32)
    top = np.zeros((k,) + clear.shape, dtype=np.uint8)
    init = np.where(bulk, np.float32(np.inf), np.float32(-np.inf)).astype(np.float32)
    for d in _STEPS:
        value = init.copy()
        a = next(i for i, v in enumerate(d) if v)
        sign = d[a]
        n = clear.shape[a]
        for i in (range(n - 2, -1, -1) if sign == 1 else range(1, n)):
            target, source = [slice(None)] * 3, [slice(None)] * 3
            target[a], source[a] = i, i + sign
            for j in range(3):
                if j != a and d[j]:
                    target[j] = slice(None, -1) if d[j] == 1 else slice(1, None)
                    source[j] = slice(1, None) if d[j] == 1 else slice(None, -1)
            t, s = tuple(target), tuple(source)
            value[t] = np.where(bulk[t], value[t], np.minimum(clear[t], value[s]))
        x = np.searchsorted(thresholds, value, side="left").astype(np.uint8)
        del value
        for i in range(k):
            high = np.maximum(top[i], x)
            np.minimum(top[i], x, out=x)
            top[i] = high
    return top[k - 1].copy()


def _snap_start(space, spatial, point, seed_tolerance, probe, outside_box, start):
    """Snap the start to the widest free node near ``point`` (updates ``start``).

    Returns ``(root, root_pos, root_clear)`` or ``(status, reason)`` when unresolved.
    """
    import numpy as np
    p = np.asarray(point)
    spacing = space.spacing
    lo_i = np.maximum(np.floor((p - seed_tolerance - space.origin) / spacing).astype(int), 0)
    hi_i = np.minimum(np.ceil((p + seed_tolerance - space.origin) / spacing).astype(int) + 1, np.asarray(space.dims))
    if np.any(hi_i <= lo_i):
        return "unresolved_seed_blocked", "start lies outside the grid box"
    box = tuple(slice(a, b) for a, b in zip(lo_i, hi_i))
    idx = np.argwhere(np.ones(tuple(hi_i - lo_i), bool)) + lo_i
    dist = np.linalg.norm(space.position(idx) - p, axis=1)
    c = space.clearance[box].ravel().astype(float)
    bulk = outside_box(box).ravel()
    near = (dist <= seed_tolerance) & (c > probe)
    start["point_clearance_A"] = float(spatial.nearest_surface(tuple(point))[0])
    if not near.any():
        return ("unresolved_seed_blocked",
                f"no free node (clearance > {probe:g} Å) within {seed_tolerance:g} Å of the start")
    inner = near & ~bulk
    if not inner.any():
        return "unresolved_seed_in_bulk", "every free node near the start is in bulk"
    order = np.lexsort((dist[inner], -c[inner]))
    pick = np.flatnonzero(inner)[order[0]]
    root = idx[pick]
    root_pos = space.position(root)
    root_clear = float(spatial.nearest_surface(tuple(root_pos))[0])
    start.update(node_A=root_pos.tolist(), node_clearance_A=root_clear, snap_distance_A=float(dist[pick]))
    return root, root_pos, root_clear


def _side_allowed(space, root, root_pos, root_clear, axis_vec, sign):
    """Half-space mask (``sign`` side of the plane through the start normal to the axis) plus the start sphere."""
    import numpy as np
    dims, spacing = space.dims, space.spacing
    proj = sum(((space.origin[d] + np.arange(dims[d]) * spacing - root_pos[d]) * axis_vec[d])
               .reshape([-1 if k == d else 1 for k in range(3)]) for d in range(3))
    allowed = sign * proj >= -1e-9
    r = int(math.ceil(root_clear / spacing)) + 1
    lo = np.maximum(root - r, 0)
    hi = np.minimum(root + r + 1, np.asarray(dims))
    local = np.argwhere(np.ones(tuple(hi - lo), bool)) + lo
    inside = np.linalg.norm(space.position(local) - root_pos, axis=1) <= root_clear
    allowed[tuple(local[inside].T)] = True
    return allowed


def _routes_overlap(pos_a, pos_b, ids_a, ids_b, separation, route_separation):
    """Whether two routes are *not* distinct (``DISTINCT_PATH_DEFINITION``)."""
    import numpy as np
    if np.linalg.norm(pos_a[-1] - pos_b[-1]) <= separation:
        return True
    shared = 0
    while shared < min(len(ids_a), len(ids_b)) and ids_a[shared] == ids_b[shared]:
        shared += 1
    junction = pos_a[shared - 1] if shared else None
    one, two = pos_a[shared:], pos_b[shared:]
    if junction is not None:
        one = one[np.linalg.norm(one - junction, axis=1) > route_separation]
        two = two[np.linalg.norm(two - junction, axis=1) > route_separation]
    return not len(one) or not len(two) or _min_pair_distance(one, two) <= route_separation


def _side_routes(space, outside, root, root_pos, root_clear, probe, axis_vec, side, sign, max_nodes,
                 separation, route_separation, max_paths):
    """Widest routes and the distinct ones kept on one side (or everywhere) at one probe."""
    import time
    import numpy as np
    began = time.perf_counter()
    allowed = None if sign is None else _side_allowed(space, root, root_pos, root_clear, axis_vec, sign)
    exempt = lambda q: np.linalg.norm(q - root_pos, axis=1) <= root_clear
    search = _widest_routes(space, [root], exempt, probe, max_nodes, allowed=allowed, outside=outside)
    del allowed
    stats = {"side": side, "search_domain_nodes": search["domain_nodes"], "terminal_nodes": search["terminal_nodes"],
             "nodes_explored": search["explored"], "terminals_reached": len(search["terminals"])}
    widths = search["width"][search["terminals"]] if search["terminals"] else np.zeros(0)
    widest = float(widths.max()) if len(widths) else -math.inf
    positions = space.position(search["ijk"])
    side_kept, truncated, mouths = [], False, []
    for term in search["terminals"]:
        # Same test as _routes_overlap, with the cheap exit-point distance checked first for all kept routes.
        if mouths and np.min(np.linalg.norm(np.asarray(mouths) - positions[term], axis=1)) <= separation:
            continue
        route = _route(search, term)
        if any(_routes_overlap(positions[route], positions[other], route, other, separation, route_separation)
               for other in side_kept):
            continue
        if len(side_kept) >= max_paths:
            truncated = True
            break
        side_kept.append(route)
        mouths.append(positions[term])
    kept = []
    for route in side_kept:
        ijk = search["ijk"][route]
        steps = np.linalg.norm(np.diff(positions[route], axis=0), axis=1)
        width = float(search["width"][route[-1]])
        kept.append({"side": side, "ids": np.ravel_multi_index(ijk.T, space.dims).tolist(),
                     "walked": np.concatenate([[0.0], np.cumsum(steps)]).tolist(),
                     "width": width if math.isfinite(width) else root_clear,
                     "end": "bulk" if space.in_bulk[tuple(ijk[-1])] else "exposed", "probe": probe})
    return {"stats": stats, "kept": kept, "truncated": truncated, "widest_A": widest,
            "seconds": time.perf_counter() - began}


@radii_option
def find_paths(frame, *, seed=None, seed_residues=None, region=None, probe_radius: float | str | None = "auto",
               spacing: float = 0.5, bulk_radius: float = 6.0, seed_tolerance: float = 3.0,
               max_paths: int = MAX_PATHS, separation: float | None = None, route_separation: float | None = None,
               axis=None, include_hetero: bool = False, lining_shell: float = 1.4, relax_iterations: int = 60,
               min_enclosure: float = 0.75, max_grid_points: int = 20_000_000, max_nodes: int = 6_000_000,
               workers: int | None = None, radii: RadiusSet | str | None = None) -> PathNetwork:
    """Curved, possibly branched, widest routes from a start to bulk solvent.

    See the module docstring for every definition (grid, bulk, widest route,
    distinct routes, junctions, centre line, radius, lining residues) and
    ``AUTO_PROBE_RULE`` for the automatic probe.

    Parameters
    ----------
    frame : StructureFrame
        Structure. Non-hetero atoms define the grid box; obstacles are heavy
        atoms (plus hetero atoms with ``include_hetero``).
    seed : Coord, optional
        Start point, Å.
    seed_residues : str or sequence of str, optional
        Start at the centroid of these residues (:func:`residue_seed`).
    region : dict, str or Path, optional
        Region definition; start at its landmark anchor (:func:`region_seed`).
        Exactly one of ``seed``, ``seed_residues`` and ``region`` is required.
    probe_radius : float, "auto" or None, default "auto"
        Required node clearance, Å (threshold only; radii are not reduced).
        ``"auto"`` (or ``None``) tries every probe of ``PROBE_LADDER`` below
        ``bulk_radius`` on one shared clearance field and keeps, per side,
        the probe chosen by ``AUTO_PROBE_RULE``; the choice is recorded in
        ``settings["probe_selection"]``. A number runs that probe only,
        exactly as before.
    spacing : float, default 0.5
        Lattice spacing, Å.
    bulk_radius : float, default 6.0
        Rolling bulk-probe radius, Å; also sets the default separations.
    seed_tolerance : float, default 3.0
        The start snaps to the free node of largest clearance within this
        distance of the requested point (ties: nearest), Å.
    max_paths : int, default 8
        Largest number of distinct routes.
    separation : float, optional
        Minimum distance between distinct exit points, Å; default ``2 * bulk_radius``.
    route_separation : float, optional
        Minimum distance between the parts of distinct routes beyond their
        junction, Å; default ``bulk_radius``.
    axis : str or Coord, optional
        ``x``, ``y``, ``z`` or a vector: search each half-space (``+``/``-``
        side of the plane through the start normal to the axis) separately,
        label each branch with its side and record a side without a route as
        unresolved. Without it one search covers all free space.
    include_hetero : bool, default False
        Treat hetero atoms (ligands; for MD input lipids and ions) as obstacles.
    lining_shell : float, default 1.4
        Lining shell, Å (``LINING_DEFINITION``).
    relax_iterations : int, default 60
        Centre-line relaxation iterations (0 keeps the straightened lattice route).
    min_enclosure : float, default 0.75
        A free node counts as *exposed* (and ends a route like bulk) when more
        than ``26 - ceil(26 * min_enclosure)`` of the 26 lattice rays from it
        reach bulk before a node with clearance <= ``probe_radius``; 0
        disables the exposure rule. Fixed for every candidate probe.
    max_grid_points : int, default 20000000
        Grid allocation guard.
    max_nodes : int, default 6000000
        Search-domain guard.
    workers : int, optional
        Worker processes for the candidate probes of ``"auto"`` (default:
        automatic, ``min(4, CPUs, candidate searches)``, reduced when
        available memory is short; see
        :func:`crevice.cavity_trajectory.auto_workers`). Results do not
        depend on it. Ignored for an explicit probe.
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

    Returns
    -------
    PathNetwork
        ``status`` ``"resolved"`` with branches, junctions and segments, or an
        unresolved status with no branch.

    Raises
    ------
    ValueError
        For invalid parameters, a missing or ambiguous start, or a grid or
        search domain over its guard.

    Examples
    --------
    A straight tunnel (radius about 2.6 Å) through a block, from its middle:

    >>> import numpy as np
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.paths import find_paths
    >>> g = np.arange(-9, 9.1, 1.5)
    >>> atoms = [Atom(i + 1, "C", "ALA", "A", i + 1, x, y, z, "C")
    ...          for i, (x, y, z) in enumerate((x, y, z) for x in g for y in g for z in g if x * x + y * y > 18)]
    >>> frame = StructureFrame(tuple(atoms))
    >>> result = find_paths(frame, seed=(0, 0, 0), probe_radius=0.8, relax_iterations=0)
    >>> result.status, [(round(b.length, 1), round(b.bottleneck_radius, 1)) for b in result.branches]
    ('resolved', [(7.5, 2.6), (9.0, 2.6)])

    With the automatic probe the same two routes are found and the choice is recorded:

    >>> auto = find_paths(frame, seed=(0, 0, 0), relax_iterations=0, workers=1)
    >>> auto.settings["probe_selection"]["chosen_probe_A"], len(auto.branches)
    (0.5, 2)
    """
    import numpy as np
    from .spatial import SpatialIndex
    auto = probe_radius is None or (isinstance(probe_radius, str) and probe_radius.strip().lower() == "auto")
    if not auto:
        if isinstance(probe_radius, str):
            try:
                probe_radius = float(probe_radius)
            except ValueError:
                raise ValueError("probe_radius must be a number or 'auto'") from None
        probe_radius = float(probe_radius)
    for name, value in ((("probe_radius", probe_radius),) if not auto else ()) + (
            ("spacing", spacing), ("bulk_radius", bulk_radius), ("seed_tolerance", seed_tolerance)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not auto and bulk_radius <= probe_radius:
        raise ValueError("bulk_radius must exceed probe_radius")
    probes = [p for p in PROBE_LADDER if p < bulk_radius] if auto else [probe_radius]
    if not probes:
        raise ValueError("bulk_radius must exceed the smallest automatic probe "
                         f"({PROBE_LADDER[0]:g} Å); give an explicit probe_radius")
    if not math.isfinite(min_enclosure) or not 0 <= min_enclosure <= 1:
        raise ValueError("min_enclosure must lie in [0, 1]")
    if int(max_paths) < 1:
        raise ValueError("max_paths must be a positive integer")
    if workers is not None and int(workers) < 1:
        raise ValueError("workers must be a positive integer")
    if sum(x is not None for x in (seed, seed_residues, region)) != 1:
        raise ValueError("give exactly one start: seed, seed_residues or region")
    separation = 2.0 * bulk_radius if separation is None else float(separation)
    route_separation = bulk_radius if route_separation is None else float(route_separation)
    axis_vec = _axis_vector(axis)
    start = {}
    if seed is not None:
        point = tuple(float(v) for v in seed)
        start.update(kind="point")
    elif seed_residues is not None:
        point, matched = residue_seed(frame, seed_residues)
        start.update(kind="residues", residues=matched)
    else:
        point, matched = region_seed(frame, region)
        start.update(kind="region_landmarks", residues=matched)
    if len(point) != 3 or not all(math.isfinite(v) for v in point):
        raise ValueError("seed must be a finite x,y,z coordinate")
    start["requested_A"] = list(point)
    atoms = frame.selected_atoms(include_hydrogen=False, include_hetero=include_hetero)
    core = np.asarray([a.coord for a in atoms if not a.hetero] or [a.coord for a in atoms], dtype=float)
    import time
    clock = time.perf_counter()
    space = FreeSpace(atoms, core, spacing, bulk_radius, max_grid_points)
    spatial = SpatialIndex(space.atoms)
    field_s = time.perf_counter() - clock
    space.outside = space.in_bulk
    required = math.ceil(26 * min_enclosure)
    levels, levels_s = None, 0.0
    if not auto and min_enclosure > 0:
        blocked = space.clearance <= probe_radius
        space.exposed = (_ray_escapes(blocked, space.in_bulk) > 26 - required) & ~blocked & ~space.in_bulk
        space.outside = space.in_bulk | space.exposed
        del blocked
    elif auto and min_enclosure > 0 and required > 0:
        clock = time.perf_counter()
        levels = _exposure_levels(space.clearance, space.in_bulk, probes, required)
        levels_s = time.perf_counter() - clock
    settings = {"probe_radius_A": "auto" if auto else probe_radius, "spacing_A": spacing,
                "bulk_radius_A": bulk_radius,
                "seed_tolerance_A": seed_tolerance, "max_paths": int(max_paths), "separation_A": separation,
                "route_separation_A": route_separation, "axis": axis_vec.tolist() if axis_vec is not None else None,
                "include_hetero": include_hetero, "lining_shell_A": lining_shell, "min_enclosure": min_enclosure,
                "relax_iterations": int(relax_iterations), "grid_origin_A": space.origin.tolist(),
                "grid_dimensions": list(space.dims), "grid_points": space.grid_points,
                "obstacle_atoms": len(space.atoms),
                "definitions": {"path": PATH_DEFINITION, "bulk": BULK_DEFINITION, "distinct": DISTINCT_PATH_DEFINITION,
                                "radius": RADIUS_DEFINITION, "lining": LINING_DEFINITION},
                "diagonal_step_note": "diagonal lattice steps are not tested analytically during the search; "
                                      "reported centre lines are measured exactly",
                **radii_fields()}

    def unresolved(status, reason):
        return PathNetwork(status, {**start, "reason": reason}, [], [], [], [{"status": status, "reason": reason}],
                           settings)

    runs = [(None, None)] if axis_vec is None else [("+", 1.0), ("-", -1.0)]
    if not auto:
        snap = _snap_start(space, spatial, point, seed_tolerance, probe_radius, lambda box: space.outside[box], start)
        if len(snap) == 2:
            return unresolved(*snap)
        root, root_pos, root_clear = snap
        kept, searches, truncated = [], [], False
        for side, sign in runs:
            found = _side_routes(space, space.outside, root, root_pos, root_clear, probe_radius, axis_vec, side, sign,
                                 max_nodes, separation, route_separation, max_paths)
            searches.append(found["stats"])
            kept += found["kept"]
            truncated = truncated or found["truncated"]
        settings.update(exposed_nodes=int(space.exposed.sum()) if min_enclosure > 0 else 0, searches=searches,
                        paths_truncated=truncated,
                        search_mode="one half-space search per side of the start along the axis"
                        if axis_vec is not None else "one search over all free space")
        if not kept:
            reached = sum(x["search_domain_nodes"] for x in searches)
            return unresolved("unresolved_no_bulk_path",
                              f"no route of free nodes (clearance > {probe_radius:g} Å) from the start reaches bulk "
                              f"({reached} free nodes reachable)")
        return _assemble_paths(frame, space, spatial, kept, root_pos, start, settings, axis_vec, lining_shell,
                               include_hetero, relax_iterations)
    selection = _auto_probe_search(space, spatial, levels, probes, point, seed_tolerance, start, runs, axis_vec,
                                   max_nodes, separation, route_separation, max_paths, workers)
    selection["record"].setdefault("timing_s", {}).update(clearance_field_s=round(field_s, 3),
                                                          exposure_levels_s=round(levels_s, 3))
    settings.update(probe_selection=selection["record"], searches=selection["searches"],
                    paths_truncated=selection["truncated"],
                    exposed_nodes=selection["exposed_nodes"],
                    search_mode="one half-space search per side of the start along the axis, probe chosen per side"
                    if axis_vec is not None else "one search over all free space")
    chosen = selection["record"]["chosen_probe_A"]
    if chosen is not None:
        settings["probe_radius_A"] = chosen
    if selection["status"] is not None:
        return unresolved(selection["status"], selection["reason"])
    return _assemble_paths(frame, space, spatial, selection["kept"], selection["root_pos"], start, settings,
                           axis_vec, lining_shell, include_hetero, relax_iterations, auto=True)


def _assemble_paths(frame, space, spatial, kept, root_pos, start, settings, axis_vec, lining_shell, include_hetero,
                    relax_iterations, auto=False):
    """Tree, segments, centre lines, radii and lining of the kept routes (each carries its probe)."""
    import numpy as np
    dims, spacing = space.dims, space.spacing
    kept.sort(key=lambda k: (-k["width"], k["walked"][-1]))
    gpos = lambda n: space.origin + np.asarray(np.unravel_index(int(n), dims), dtype=float) * spacing
    walked_of = {}
    for k in kept:
        for n, w in zip(k["ids"], k["walked"]):
            walked_of.setdefault(n, w)
    # Tree: junctions = nodes with >= 2 distinct children among kept routes.
    children = {}
    for k in kept:
        for a, b in zip(k["ids"], k["ids"][1:]):
            children.setdefault(a, set()).add(b)
    junction_nodes = sorted({n for n, ch in children.items() if len(ch) >= 2}, key=lambda n: (walked_of[n], n))
    junction_ids = {n: i + 1 for i, n in enumerate(junction_nodes)}
    cut = set(junction_nodes) | {kept[0]["ids"][0]}
    segments, segment_ids = [], {}
    branch_segments = []
    for k_route in kept:
        route = k_route["ids"]
        ids, startk = [], 0
        marks = [k for k, n in enumerate(route) if n in cut and k > 0] + [len(route) - 1]
        for k in sorted(set(marks)):
            key = (route[startk], route[k])
            if key not in segment_ids:
                segment_ids[key] = len(segment_ids) + 1
                segments.append({"segment": segment_ids[key], "nodes": route[startk:k + 1],
                                 "probe": k_route["probe"]})
            else:
                seg = segments[segment_ids[key] - 1]
                seg["probe"] = min(seg["probe"], k_route["probe"])
            ids.append(segment_ids[key])
            startk = k
        branch_segments.append(ids)
    # Centre lines.
    for seg in segments:
        node_ids = seg.pop("nodes")
        probe = seg.pop("probe")
        nodes = np.asarray([gpos(n) for n in node_ids])
        lattice = np.asarray(_straighten(spatial, nodes))
        floor = float(min(space.clearance[np.unravel_index(n, dims)] for n in node_ids))
        line, mode = lattice, "lattice_straightened"
        if relax_iterations > 0 and len(nodes) >= 3:
            relaxed = _relax(space, nodes, floor, int(relax_iterations), spacing)
            dense = _resample(relaxed, 0.25)
            if float(spatial.nearest_surface_values([tuple(q) for q in dense]).min()) > probe:
                line, mode = relaxed, "relaxed"
            else:
                mode = "lattice_fallback"
        seg["points_A"] = [list(map(float, q)) for q in line]
        seg["centre_line"] = mode
    seg_by_id = {s["segment"]: s for s in segments}
    for seg in segments:
        first = tuple(seg["points_A"][0])
        last = tuple(seg["points_A"][-1])
        seg["from"] = _node_kind(first, root_pos, junction_nodes, gpos, junction_ids)
        seg["to"] = _node_kind(last, root_pos, junction_nodes, gpos, junction_ids)
    # Branch samples, radii, bottleneck, lining.
    from .fingerprint import _lining
    branches = []
    junction_records = {j: {"junction": junction_ids[j], "position_A": gpos(j).tolist(), "branches": [],
                            "is_start": False} for j in junction_nodes}
    for rank, (k_route, seg_ids) in enumerate(zip(kept, branch_segments), 1):
        route = k_route["ids"]
        poly, ends = [], []
        for sid in seg_ids:
            pts = seg_by_id[sid]["points_A"]
            poly.extend(pts if not poly else pts[1:])
            ends.append(float(np.linalg.norm(np.diff(np.asarray(poly), axis=0), axis=1).sum()))
        dense = _resample(poly, 0.25)
        radii_s = spatial.nearest_surface_values([tuple(q) for q in dense]).tolist()
        arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(dense, axis=0), axis=1))]).tolist()
        worst = int(np.argmin(radii_s))
        rows = [{"position": tuple(q), "raw_clearance_A": r} for q, r in zip(dense.tolist(), radii_s)]
        per_sample, records = _lining(frame, rows, lining_shell, include_hetero)
        lining = sorted(({"residue": key.label, "min_gap_A": round(rec["gap"], 3),
                          "anchor_arc_length_A": round(arc[rec["index"]], 3), "samples": len(rec["samples"])}
                         for key, rec in records.items()), key=lambda r: (r["anchor_arc_length_A"], r["residue"]))
        junctions = []
        for sid, end in zip(seg_ids, ends):
            to = seg_by_id[sid]["to"]
            if to["kind"] == "junction":
                n = junction_nodes[to["junction"] - 1]
                junctions.append({"junction": to["junction"], "arc_length_A": round(end, 3)})
                junction_records[n]["branches"].append(rank)
        if route[0] in junction_ids:
            junctions.insert(0, {"junction": junction_ids[route[0]], "arc_length_A": 0.0})
            junction_records[route[0]]["branches"].append(rank)
        metadata = {"centre_line": {str(s): seg_by_id[s]["centre_line"] for s in seg_ids},
                    "lattice_route_length_A": k_route["walked"][-1], "lattice_nodes": len(route),
                    "below_probe": bool(radii_s[worst] <= k_route["probe"])}
        if auto:
            metadata["probe_radius_A"] = k_route["probe"]
        branches.append(PathBranch(
            rank, [tuple(map(float, q)) for q in dense], radii_s, arc, arc[-1], radii_s[worst],
            tuple(map(float, dense[worst])), arc[worst], tuple(map(float, gpos(route[-1]))),
            k_route["width"], k_route["side"], k_route["end"], seg_ids, junctions, lining,
            sorted(k.label for k in per_sample[worst]), metadata))
    junction_list = list(junction_records.values())
    for rec in junction_list:
        rec["is_start"] = bool(np.allclose(rec["position_A"], root_pos))
    unresolved_rows = []
    if axis_vec is not None:
        for side in ("+", "-"):
            if not any(b.side == side for b in branches):
                reason = f"no distinct route to bulk on the {side} side of the start along the axis"
                if auto:
                    reason += " at any candidate probe"
                unresolved_rows.append({"status": "unresolved_no_bulk_path_on_side", "side": side, "reason": reason})
    return PathNetwork("resolved", start, branches, junction_list, segments, unresolved_rows, settings)


def _node_kind(point, root_pos, junction_nodes, gpos, junction_ids):
    import numpy as np
    if np.allclose(point, root_pos):
        return {"kind": "start"}
    for n in junction_nodes:
        if np.allclose(point, gpos(n)):
            return {"kind": "junction", "junction": junction_ids[n]}
    return {"kind": "mouth"}


# --------------------------------------------------------------------------
# Automatic probe
# --------------------------------------------------------------------------

_AUTO_WORKER = None
#: Prune ladder probes at or above the smallest probe's widest route width (False only for benchmarks).
_PRUNE = True
_SAME_WIDTH_A = 0.05
_SAME_WIDTH_FRACTION = 0.02


def _init_auto_worker(state):
    global _AUTO_WORKER
    _AUTO_WORKER = state


def _auto_task(task):
    return task, _candidate_routes(_AUTO_WORKER, *task)


def _candidate_routes(state, j, side, sign):
    """Routes at ladder probe ``j`` on one side, on the shared clearance field."""
    space, levels, probes = state["space"], state["levels"], state["probes"]
    p = probes[j]
    outside = space.in_bulk if levels is None else space.in_bulk | ((levels > j) & (space.clearance > p))
    found = _side_routes(space, outside, state["root"], state["root_pos"], state["root_clear"], p, state["axis"],
                         side, sign, state["max_nodes"], state["separation"], state["route_separation"],
                         state["max_paths"])
    found["exposed_nodes"] = int((outside & ~space.in_bulk).sum())
    return found


def _same_route_set(a, b, space, separation, route_separation):
    """Whether two kept-route lists (different probes) are the same routes: equal counts and a one-to-one
    matching in which every pair is *not distinct* under ``DISTINCT_PATH_DEFINITION`` and has ranking widths
    within ``max(0.05 Å, 2 %)`` (the tolerance of ``crevice profile``'s automatic probe)."""
    import numpy as np
    if len(a) != len(b):
        return False
    pos = lambda r: space.position(np.array(np.unravel_index(np.asarray(r["ids"]), space.dims)).T)
    pa, pb = [pos(r) for r in a], [pos(r) for r in b]
    same_width = lambda x, y: abs(x - y) <= max(_SAME_WIDTH_A, _SAME_WIDTH_FRACTION * max(x, y))
    ok = [[same_width(a[i]["width"], b[k]["width"]) and
           _routes_overlap(pa[i], pb[k], a[i]["ids"], b[k]["ids"], separation, route_separation)
           for k in range(len(b))] for i in range(len(a))]
    match = [-1] * len(b)

    def augment(i, seen):
        for k in range(len(b)):
            if ok[i][k] and not seen[k]:
                seen[k] = True
                if match[k] < 0 or augment(match[k], seen):
                    match[k] = i
                    return True
        return False
    return all(augment(i, [False] * len(b)) for i in range(len(a)))


def _auto_probe_search(space, spatial, levels, probes, point, seed_tolerance, start, runs, axis_vec, max_nodes,
                       separation, route_separation, max_paths, workers):
    """Run the probe ladder on one clearance field and choose a probe per side (``AUTO_PROBE_RULE``)."""
    import time
    import numpy as np
    f32 = lambda v: float(np.float32(v))

    def outside_box(j):
        if levels is None:
            return lambda box: space.in_bulk[box]
        return lambda box: space.in_bulk[box] | ((levels[box] > j) & (space.clearance[box] > probes[j]))

    rows, results = {}, {}
    record = {"mode": "auto", "rule": AUTO_PROBE_RULE, "ladder_A": list(PROBE_LADDER), "candidates_A": list(probes),
              "min_enclosure_note": "min_enclosure is the same for every candidate probe",
              "start_snap_probe_A": None, "chosen_probe_A": None, "chosen_by_side": {}, "reason": None,
              "candidates": [], "groups": []}
    out = {"record": record, "searches": [], "truncated": False, "exposed_nodes": 0, "status": None, "reason": None,
           "kept": [], "root_pos": None}
    # Start: snapped at the smallest candidate that resolves it, then shared by every candidate.
    snap, j0 = None, None
    for j, p in enumerate(probes):
        trial = dict(start)
        snap = _snap_start(space, spatial, point, seed_tolerance, p, outside_box(j), trial)
        if len(snap) == 3:
            start.clear()
            start.update(trial)
            j0 = j
            break
        for side, _ in runs:
            rows[(j, side)] = {"probe_A": p, "side": side, "status": snap[0], "reason": snap[1]}
    if j0 is None:
        start.update(trial)
        record["candidates"] = [rows[k] for k in sorted(rows, key=lambda k: (k[0], str(k[1])))]
        record["reason"] = f"{snap[1]} at every candidate probe ({probes[0]:g}-{probes[-1]:g} Å)"
        out.update(status=snap[0], reason=record["reason"])
        return out
    root, root_pos, root_clear = snap
    record["start_snap_probe_A"] = probes[j0]
    lattice_clear = float(space.clearance[tuple(root)])
    state = {"space": space, "levels": levels, "probes": probes, "root": root, "root_pos": root_pos,
             "root_clear": root_clear, "axis": axis_vec, "max_nodes": max_nodes, "separation": separation,
             "route_separation": route_separation, "max_paths": max_paths}
    tasks = []
    for j in range(j0, len(probes)):
        for side, sign in runs:
            if not lattice_clear > f32(probes[j]):
                rows[(j, side)] = {"probe_A": probes[j], "side": side, "status": "unresolved_seed_blocked",
                                   "reason": "the start node's clearance does not exceed this probe"}
            else:
                tasks.append((j, side, sign))
    first = {}
    for task in tasks:
        first.setdefault(task[1], task)
    from .cavity_trajectory import auto_workers
    if workers is None:
        workers = auto_workers(max(1, len(tasks)), max_grid_points=space.grid_points)
    workers = max(1, int(workers))

    # A route at a larger probe is also a route at the smallest one, so the smallest probe's widest route
    # width bounds every larger probe: candidates at or above it cannot have a route. They are recorded as
    # pruned whether or not a worker had already started them (such a search finds nothing), so the record
    # does not depend on the number of workers or on timing.
    widest = {}

    def prunable(task):
        j, side, _ = task
        return (_PRUNE and task != first[side] and side in widest
                and not widest[side] > f32(probes[j]))

    def finished(task, found):
        results[task[:2]] = found
        if task == first[task[1]]:
            widest[task[1]] = found["widest_A"]
            timing["smallest_probe_searches_s"] = max(timing.get("smallest_probe_searches_s", 0.0),
                                                      round(time.perf_counter() - clock, 3))

    timing = {"workers": workers}
    clock = time.perf_counter()
    if workers == 1 or len(tasks) <= 1:
        for task in tasks:
            if not prunable(task):
                finished(task, _candidate_routes(state, *task))
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed
        context = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
        with ProcessPoolExecutor(max_workers=min(workers, len(tasks)), mp_context=context,
                                 initializer=_init_auto_worker, initargs=(state,)) as pool:
            # Submitted in ascending probe order, so the smallest probe of each side starts first.
            futures = {pool.submit(_auto_task, task): task for task in tasks}
            for future in as_completed(futures):
                if future.cancelled():
                    continue
                task, found = future.result()
                finished(task, found)
                if task == first[task[1]]:
                    for other, pending in futures.items():
                        if not other.done() and prunable(pending):
                            other.cancel()
    timing["all_probe_searches_s"] = round(time.perf_counter() - clock, 3)
    for task in tasks:
        if prunable(task):
            j, side, _ = task
            results.pop((j, side), None)
            bound = widest[side]
            rows[(j, side)] = {"probe_A": probes[j], "side": side, "status": "pruned",
                               "reason": ("no route at the smallest probe" if not math.isfinite(bound) else
                                          f"at or above the widest route width {bound:.3f} Å found at "
                                          f"{probes[first[side][0]]:g} Å")}
    timing["searches_kept"] = len(results)
    timing["searches_pruned"] = sum(r["status"] == "pruned" for r in rows.values())
    timing["search_time_sum_s"] = round(sum(r["seconds"] for r in results.values()), 3)
    record["timing_s"] = timing
    for (j, side), found in results.items():
        kept = found["kept"]
        rows[(j, side)] = {"probe_A": probes[j], "side": side,
                           "status": "resolved" if kept else "unresolved_no_bulk_path",
                           "routes": len(kept), "widest_route_width_A": max((k["width"] for k in kept), default=None),
                           "route_widths_A": [round(k["width"], 4) for k in kept],
                           "exit_points_A": [[round(v, 3) for v in space.position(
                               np.unravel_index(k["ids"][-1], space.dims))] for k in kept],
                           "end_types": [k["end"] for k in kept],
                           "lattice_lengths_A": [round(k["walked"][-1], 3) for k in kept],
                           "routes_truncated": found["truncated"], "exposed_nodes": found["exposed_nodes"],
                           "search_s": round(found["seconds"], 3),
                           **{key: found["stats"][key] for key in ("search_domain_nodes", "nodes_explored")}}
    # Choose per side.
    chosen = {}
    for side, _ in runs:
        side_groups = []
        for j in range(len(probes)):
            res = results.get((j, side))
            if res is None or not res["kept"]:
                continue
            for group in side_groups:
                if _same_route_set(results[(group["members"][0], side)]["kept"], res["kept"], space, separation,
                                   route_separation):
                    group["members"].append(j)
                    break
            else:
                side_groups.append({"members": [j]})
        for g, group in enumerate(side_groups, 1):
            group["group"] = g
            for j in group["members"]:
                rows[(j, side)]["group"] = g
        pick = _pick_group(side_groups)
        label = "all" if side is None else side
        if pick is None:
            chosen[side] = None
            record["chosen_by_side"][label] = {"probe_A": None, "reason": "no route at any candidate probe"}
        else:
            group, why = pick
            j = group["members"][0]
            chosen[side] = j
            record["chosen_by_side"][label] = {"probe_A": probes[j], "group": group["group"],
                                               "group_probes_A": [probes[m] for m in group["members"]],
                                               "routes": len(results[(j, side)]["kept"]), "reason": why}
        record["groups"] += [{"side": side, "group": g["group"], "probes_A": [probes[m] for m in g["members"]],
                              "routes": len(results[(g["members"][0], side)]["kept"])} for g in side_groups]
    record["candidates"] = [rows[k] for k in sorted(rows, key=lambda k: (k[0], str(k[1])))]
    picked = [probes[j] for j in chosen.values() if j is not None]
    record["chosen_probe_A"] = min(picked) if picked else None
    record["reason"] = "; ".join(f"{s}: {v['reason']}" if len(runs) > 1 else v["reason"]
                                 for s, v in record["chosen_by_side"].items())
    for side, j in chosen.items():
        if j is None:
            continue
        found = results[(j, side)]
        out["searches"].append({**found["stats"], "probe_radius_A": probes[j]})
        out["kept"] += found["kept"]
        out["truncated"] = out["truncated"] or found["truncated"]
        out["exposed_nodes"] = max(out["exposed_nodes"], found["exposed_nodes"])
    out["root_pos"] = root_pos
    if not out["kept"]:
        reached = max((r.get("search_domain_nodes", 0) for r in rows.values()), default=0)
        out.update(status="unresolved_no_bulk_path",
                   reason=f"no route of free nodes from the start reaches bulk at any candidate probe "
                          f"({probes[0]:g}-{probes[-1]:g} Å; {reached} free nodes reachable at the smallest)")
    return out


def _pick_group(groups):
    """``AUTO_PROBE_RULE`` choice among the route-set groups of one side: ``(group, reason)`` or ``None``."""
    if not groups:
        return None
    eligible = [g for g in groups if len(g["members"]) >= 2]
    if eligible:
        best = max(eligible, key=lambda g: (len(g["members"]), -g["members"][0]))
        return best, (f"smallest probe of the most persistent route set (found unchanged at "
                      f"{len(best['members'])} candidate probes)")
    best = groups[0]
    return best, "no route set persisted over two candidate probes; smallest probe with a route used"


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------

def _write_csv(path, rows, columns):
    import csv
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    return path


def plot_paths(result, path, *, dpi: int = 200):
    """Free-sphere radius against arc length, one line per branch, junctions marked.

    Parameters
    ----------
    result : PathNetwork
        From :func:`find_paths`.
    path : str or Path
        Output image (PNG recommended; the description is embedded).
    dpi : int, default 200
        Resolution.

    Returns
    -------
    pathlib.Path
        The written figure.

    Notes
    -----
    6.4 x 3.8 inch. Each branch is a line in its Okabe-Ito colour
    (``BRANCH_COLORS``) and its own line style; junctions are hollow circles
    on each branch that passes them; the bottleneck of each branch is a small
    filled dot; a thin dotted grey line marks the probe radius (one per
    side when the automatic probe chose different probes). Legend entries
    name the branches (with side), the junction marker and the probe; no
    titles or value labels unless annotating (``--annotate``). An unresolved
    result gives an empty panel with a notice of its status.
    """
    from pathlib import Path
    from matplotlib.lines import Line2D
    from .figures import _pyplot, _save_figure, _style_axes
    from .presentation import FigureText
    plt = _pyplot()
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    text = FigureText("Free-sphere radius against arc length from the start for each route of a pathway tree "
                      "to bulk; hollow circles mark junctions where routes separate")
    fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
    probes = sorted({b.metadata.get("probe_radius_A", result.settings.get("probe_radius_A", 0.0))
                     for b in result.branches})
    if not result.branches:
        text.notice(ax, 0.5, 0.5, f"No route to bulk ({result.status})", transform=ax.transAxes, ha="center",
                    va="center", fontsize=10, color="#4a4a4a")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    for b in result.branches:
        color = BRANCH_COLORS[(b.id - 1) % len(BRANCH_COLORS)]
        style = BRANCH_STYLES[(b.id - 1) % len(BRANCH_STYLES)]
        label = f"Branch {b.id}" + (f" ({b.side})" if b.side else "")
        ax.plot(b.arc_lengths, b.radii, color=color, linestyle=style, linewidth=1.8, label=label)
        ax.scatter([b.bottleneck_arc_length], [b.bottleneck_radius], s=22, color=color, zorder=4)
        for j in b.junctions:
            k = min(range(len(b.arc_lengths)), key=lambda i: abs(b.arc_lengths[i] - j["arc_length_A"]))
            ax.scatter([b.arc_lengths[k]], [b.radii[k]], s=46, facecolor="white", edgecolor=color, linewidth=1.4,
                       zorder=5)
        text.note(f"Branch {b.id}: length {b.length:.2f} Å, bottleneck {b.bottleneck_radius:.2f} Å at "
                  f"{b.bottleneck_arc_length:.2f} Å" + (f", side {b.side}" if b.side else ""))
        text.annotate(ax, f"{b.bottleneck_radius:.2f} Å", (b.bottleneck_arc_length, b.bottleneck_radius),
                      xytext=(5, -12), textcoords="offset points", fontsize=8, color="#333333")
    if result.branches:
        for probe in probes:
            ax.axhline(probe, color="#8d969e", linestyle=":", linewidth=1.0)
        handles, labels = ax.get_legend_handles_labels()
        if any(b.junctions for b in result.branches):
            handles.append(Line2D([], [], marker="o", linestyle="", markerfacecolor="white",
                                  markeredgecolor="#333333", markersize=6.5))
            labels.append("Junction")
        handles.append(Line2D([], [], color="#8d969e", linestyle=":", linewidth=1.0))
        labels.append("Probe radius" if len(probes) == 1 else
                      "Probe radii (" + ", ".join(f"{p:g}" for p in probes) + " Å)")
        ax.legend(handles, labels, frameon=False, fontsize=8, loc="best")
        ax.set_ylim(bottom=0.0)
        ax.set_xlim(left=0.0)
    for u in result.unresolved:
        text.note(f"Unresolved: {u.get('reason', u.get('status'))}")
    ax.set_xlabel("Arc length from start (Å)")
    ax.set_ylabel("Free-sphere radius (Å)")
    text.title(ax, "CREVICE pathway radii")
    ax.grid(True, color="#e3e3e3", linewidth=0.6)
    _style_axes(ax)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


def _rgb(hex_color):
    return tuple(int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))


def _scene_lines(result):
    """(name, rgb, polyline) per branch, decimated to about 0.5 Å steps."""
    out = []
    for b in result.branches:
        pts = b.points[::2] + ([b.points[-1]] if len(b.points) % 2 == 0 else [])
        out.append((f"path_{b.id}", _rgb(BRANCH_COLORS[(b.id - 1) % len(BRANCH_COLORS)]), pts))
    return out


def write_path_scenes(result, directory, prefix, structure_name=None, *, tube_radius: float = 0.25):
    """PyMOL, ChimeraX and VMD scene scripts: one object per branch, no labels.

    Parameters
    ----------
    result : PathNetwork
    directory : str or Path
        Output directory; scripts load files from their own directory.
    prefix : str
        File prefix.
    structure_name : str, optional
        File name (in ``directory``) of the viewer PDB to load as a grey cartoon.
    tube_radius : float, default 0.25
        Centre-line tube radius, Å (shows the route, not its width).

    Returns
    -------
    list of pathlib.Path
        ``<prefix>_paths.pml``, ``<prefix>_paths.cxc`` with one
        ``<prefix>_path_<n>.bild`` per branch, and ``<prefix>_paths_vmd.tcl``.

    Notes
    -----
    Branch colours are ``BRANCH_COLORS`` (Okabe-Ito) in rank order. These
    scripts are generated text; they have not been run in every viewer.
    """
    from pathlib import Path
    directory = Path(directory)
    lines = _scene_lines(result)
    written = []
    pml = ["# CREVICE pathway scene (PyMOL): one CGO object per branch, no labels.", "python",
           "from pathlib import Path", "import inspect", "from pymol import cmd", "from pymol.cgo import CYLINDER",
           "_crevice_script_root = Path(getattr(cmd._pymol, '__script__',",
           "    inspect.currentframe().f_code.co_filename)).resolve().parent",
           "def crevice_input(name, fallback):",
           "    for directory in (_crevice_script_root, Path.cwd(), Path(fallback).parent):",
           "        path = directory / name",
           "        if path.is_file():",
           "            return str(path)",
           "    raise FileNotFoundError(name)",
           "cmd.bg_color('white')"]
    if structure_name:
        fallback = str((directory / structure_name).resolve())
        pml += [f"cmd.load(crevice_input({structure_name!r}, {fallback!r}), 'crevice_protein')",
                "cmd.hide('everything', 'crevice_protein')", "cmd.show('cartoon', 'crevice_protein')",
                "cmd.color('grey70', 'crevice_protein')", "cmd.set('cartoon_transparency', 0.25, 'crevice_protein')"]
    for name, rgb, pts in lines:
        body = []
        for a, b in zip(pts, pts[1:]):
            body += ["CYLINDER", *[round(v, 3) for v in a], *[round(v, 3) for v in b], tube_radius, *rgb, *rgb]
        pml.append(f"cmd.load_cgo([{', '.join(str(x) for x in body)}], {name!r})")
    pml.append("cmd.orient('crevice_protein')" if structure_name else "cmd.zoom()")
    pml.append("python end")
    path = directory / f"{prefix}_paths.pml"
    path.write_text("\n".join(pml) + "\n")
    written.append(path)
    cxc = ["# CREVICE pathway scene (ChimeraX): one BILD model per branch, no labels", "set bgColor white"]
    if structure_name:
        cxc += [f"open {structure_name}", "hide atoms", "cartoon", "color #1 grey", "transparency #1 25 target c"]
    for name, rgb, pts in lines:
        bild = directory / f"{prefix}_{name}.bild"
        rows = [f".color {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f}"]
        rows += [f".cylinder {a[0]:.3f} {a[1]:.3f} {a[2]:.3f} {b[0]:.3f} {b[1]:.3f} {b[2]:.3f} {tube_radius}"
                 for a, b in zip(pts, pts[1:])]
        bild.write_text("\n".join(rows) + "\n")
        written.append(bild)
        cxc.append(f"open {bild.name} name {name}")
    cxc.append("view")
    path = directory / f"{prefix}_paths.cxc"
    path.write_text("\n".join(cxc) + "\n")
    written.append(path)
    tcl = ["# CREVICE pathway scene (VMD): one graphics molecule per branch, no labels",
           "set crevice_here [file dirname [file normalize [info script]]]", "color Display Background white"]
    if structure_name:
        tcl += [f"mol new [file join $crevice_here {structure_name}] type pdb", "mol modstyle 0 top NewCartoon",
                "mol modcolor 0 top ColorID 2"]
    for k, (name, rgb, pts) in enumerate(lines):
        cid = 23 + k if 23 + k < 33 else 32
        tcl += [f"color change rgb {cid} {rgb[0]:.3f} {rgb[1]:.3f} {rgb[2]:.3f}", "mol new",
                f"mol rename top {name}", f"graphics top color {cid}"]
        tcl += [f"graphics top cylinder {{{a[0]:.3f} {a[1]:.3f} {a[2]:.3f}}} {{{b[0]:.3f} {b[1]:.3f} {b[2]:.3f}}} "
                f"radius {tube_radius} resolution 12" for a, b in zip(pts, pts[1:])]
    path = directory / f"{prefix}_paths_vmd.tcl"
    path.write_text("\n".join(tcl) + "\n")
    written.append(path)
    return written


def paths_view(result, frame=None, *, title=None):
    """Interactive view (:class:`crevice.viewer.CreviceView`) with one line layer per branch.

    Parameters
    ----------
    result : PathNetwork
    frame : StructureFrame, optional
        Protein for the cartoon layer (hetero atoms are dropped).
    title : str, optional
        Page title.

    Returns
    -------
    crevice.viewer.CreviceView
        Layers ``protein`` and ``path_<n>`` (Okabe-Ito colours, 0.12 Å tubes,
        no text in the canvas).
    """
    from dataclasses import replace
    from .viewer import CreviceView, Layer, PROTEIN_HEX, PROTEIN_OPACITY, _frame_pdb
    layers, structure = [], None
    if frame is not None:
        protein = replace(frame, atoms=tuple(a for a in frame.atoms if not a.hetero) or frame.atoms)
        structure = _frame_pdb(protein)
        layers.append(Layer("protein", "Protein", "cartoon", PROTEIN_HEX, opacity=PROTEIN_OPACITY))
    for b in result.branches:
        layers.append(Layer(f"path_{b.id}", f"Branch {b.id}" + (f" ({b.side})" if b.side else ""), "line",
                            BRANCH_COLORS[(b.id - 1) % len(BRANCH_COLORS)],
                            data={"paths": [[list(p) for p in b.points[::2]]]},
                            description=f"length {b.length:.1f} Å, bottleneck {b.bottleneck_radius:.2f} Å"))
    return CreviceView(layers, structure=structure, title=title or "CREVICE pathways", source="crevice.paths",
                       notes=[INTERPRETATION])


def paths_view_from_json(path, *, title=None):
    """Interactive view of a written ``*_paths.json`` (used by ``crevice view``).

    Parameters
    ----------
    path : str or Path
        A ``<prefix>_paths.json`` from :func:`write_paths_outputs`; the
        protein is read from ``<prefix>_paths_structure.pdb`` beside it, if present.
    title : str, optional
        Page title.

    Returns
    -------
    crevice.viewer.CreviceView
        Layers ``protein`` and ``path_<n>``, as :func:`paths_view`.
    """
    import json
    from pathlib import Path
    from .viewer import CreviceView, Layer, PROTEIN_HEX, PROTEIN_OPACITY
    path = Path(path)
    data = json.loads(path.read_text())
    structure_path = path.with_name(path.name.replace("_paths.json", "_paths_structure.pdb"))
    structure = structure_path.read_text() if structure_path.is_file() else None
    layers = [Layer("protein", "Protein", "cartoon", PROTEIN_HEX, opacity=PROTEIN_OPACITY)] if structure else []
    for b in data.get("branches", []):
        n = b["branch"]
        layers.append(Layer(f"path_{n}", f"Branch {n}" + (f" ({b['side']})" if b.get("side") else ""), "line",
                            BRANCH_COLORS[(n - 1) % len(BRANCH_COLORS)], data={"paths": [b["points_A"][::2]]},
                            description=f"length {b['length_A']:.1f} Å, bottleneck {b['bottleneck_radius_A']:.2f} Å"))
    return CreviceView(layers, structure=structure, title=title or f"{path.stem} - CREVICE view", source=str(path),
                       notes=[INTERPRETATION])


def write_paths_outputs(result, directory, *, prefix: str = "crevice", frame=None, dpi: int = 200,
                        figure: bool = True, scenes: bool = True):
    """Write the JSON record, CSV tables, figure and viewer scenes of a pathway result.

    Parameters
    ----------
    result : PathNetwork
    directory : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File prefix.
    frame : StructureFrame, optional
        Structure for the viewer scenes (protein atoms only are written).
    dpi : int, default 200
        Figure resolution.
    figure : bool, default True
        Write ``<prefix>_paths.png``.
    scenes : bool, default True
        Write the PyMOL/ChimeraX/VMD scenes (and the viewer PDB with ``frame``).

    Returns
    -------
    list of pathlib.Path
        Written files.

    Outputs
    -------
    <prefix>_paths.json
        :meth:`PathNetwork.to_dict`.
    <prefix>_paths.csv
        :meth:`PathNetwork.summary_rows` (one row per branch, or one status row).
    <prefix>_path_points.csv
        ``branch``, ``sample``, ``arc_length_A``, ``x_A``, ``y_A``, ``z_A``, ``radius_A``.
    <prefix>_path_lining.csv
        ``branch``, ``residue``, ``min_gap_A``, ``anchor_arc_length_A``, ``samples``.
    <prefix>_path_junctions.csv
        ``junction``, ``x_A``, ``y_A``, ``z_A``, ``is_start``, ``branch``, ``arc_length_A``.
    <prefix>_path_probes.csv
        Automatic probe only: one row per candidate probe and side
        (``probe_A``, ``side``, ``status``, ``group``, ``chosen``, routes,
        widths, end types, lengths, search sizes, ``reason``).
    <prefix>_paths.png, scenes
        See :func:`plot_paths` and :func:`write_path_scenes`.
    """
    import json
    from dataclasses import replace
    from pathlib import Path
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    path = directory / f"{prefix}_paths.json"
    path.write_text(json.dumps(result.to_dict(), indent=1, default=float) + "\n")
    written.append(path)
    rows = result.summary_rows()
    written.append(_write_csv(directory / f"{prefix}_paths.csv", rows, list(rows[0])))
    points = [{"branch": b.id, "sample": i, "arc_length_A": round(s, 3), "x_A": round(p[0], 3),
               "y_A": round(p[1], 3), "z_A": round(p[2], 3), "radius_A": round(r, 4)}
              for b in result.branches for i, (s, p, r) in enumerate(zip(b.arc_lengths, b.points, b.radii))]
    written.append(_write_csv(directory / f"{prefix}_path_points.csv", points,
                              ["branch", "sample", "arc_length_A", "x_A", "y_A", "z_A", "radius_A"]))
    lining = [{"branch": b.id, **r} for b in result.branches for r in b.lining]
    written.append(_write_csv(directory / f"{prefix}_path_lining.csv", lining,
                              ["branch", "residue", "min_gap_A", "anchor_arc_length_A", "samples"]))
    junctions = []
    for j in result.junctions:
        for b in result.branches:
            for k in b.junctions:
                if k["junction"] == j["junction"]:
                    junctions.append({"junction": j["junction"], "x_A": round(j["position_A"][0], 3),
                                      "y_A": round(j["position_A"][1], 3), "z_A": round(j["position_A"][2], 3),
                                      "is_start": j["is_start"], "branch": b.id, "arc_length_A": k["arc_length_A"]})
    written.append(_write_csv(directory / f"{prefix}_path_junctions.csv", junctions,
                              ["junction", "x_A", "y_A", "z_A", "is_start", "branch", "arc_length_A"]))
    selection = result.settings.get("probe_selection")
    if selection:
        columns = ["probe_A", "side", "status", "group", "chosen", "routes", "widest_route_width_A",
                   "route_widths_A", "end_types", "lattice_lengths_A", "exposed_nodes", "search_domain_nodes",
                   "nodes_explored", "reason"]
        chosen = {side: v["probe_A"] for side, v in selection["chosen_by_side"].items()}
        probe_rows = [{**{c: (";".join(map(str, r[c])) if isinstance(r.get(c), list) else r.get(c, ""))
                          for c in columns},
                       "side": r["side"] or "all",
                       "chosen": chosen.get(r["side"] or "all") == r["probe_A"]} for r in selection["candidates"]]
        written.append(_write_csv(directory / f"{prefix}_path_probes.csv", probe_rows, columns))
    if figure:
        written.append(plot_paths(result, directory / f"{prefix}_paths.png", dpi=dpi))
    if scenes:
        structure_name = None
        if frame is not None:
            from .presentation import write_viewer_structure
            protein = replace(frame, atoms=tuple(a for a in frame.atoms if not a.hetero) or frame.atoms)
            structure = write_viewer_structure(protein, directory / f"{prefix}_paths_structure.pdb")
            structure_name = structure.name
            written += [structure, structure.with_suffix(".identities.json")]
        written += write_path_scenes(result, directory, prefix, structure_name)
    return written


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------

def _probe_arg(value):
    """``--probe-radius`` value: ``auto`` or a positive number (Å)."""
    import argparse
    if value.strip().lower() == "auto":
        return "auto"
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("use 'auto' or a number of Å") from None
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("the probe radius must be finite and positive")
    return number


def add_paths_command(subparsers):
    """Register ``crevice paths`` on the CLI subparsers.

    Parameters
    ----------
    subparsers : object
        The ``crevice`` subcommand registry (the object ``add_subparsers`` returns).

    Returns
    -------
    argparse.ArgumentParser
        The ``paths`` parser (handler :func:`run_paths_command`).
    """
    from .cli import _add_structure_arg
    from .geometry import parse_coord
    parser = subparsers.add_parser(
        "paths", help="Curved, branched widest routes from a seed or region to bulk solvent (pathway tree)")
    _add_structure_arg(parser)
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", default=None, help="File prefix (default: structure file stem)")
    start = parser.add_mutually_exclusive_group(required=True)
    start.add_argument("--seed", type=parse_coord, help="Start point x,y,z (Å)")
    start.add_argument("--seed-residues",
                       help="Start at the centroid of these residues, comma-separated: GLN282, 282, A:GLN282 or "
                            "GLN282@OE1+NE2 (atoms joined by +); each residue counts once")
    start.add_argument("--region", help="Region definition JSON; start at the mean of its landmark centroids")
    parser.add_argument("--probe-radius", type=_probe_arg, default="auto",
                        help="Required clearance (Å) of every route node; a threshold, radii are not reduced. 'auto' "
                             "(default) searches the probe ladder 0.5-3.0 Å on one shared clearance field and keeps, "
                             "per side, the smallest probe of the route set found at the most probes (recorded in "
                             "the JSON); a number runs that probe only")
    parser.add_argument("--workers", type=int, default=None,
                        help="Worker processes for the candidate probes of --probe-radius auto (default: automatic, "
                             "min(4, available CPUs, searches), reduced if available memory is short); results do "
                             "not depend on it")
    parser.add_argument("--spacing", type=float, default=0.5, help="Search lattice spacing in Å (default 0.5)")
    parser.add_argument("--bulk-radius", type=float, default=6.0,
                        help="Rolling probe radius (Å) that defines bulk solvent; also sets the default distinct-route "
                             "separations (default 6)")
    parser.add_argument("--seed-tolerance", type=float, default=3.0,
                        help="The start snaps to the widest free node within this distance (Å, default 3)")
    parser.add_argument("--max-paths", type=int, default=MAX_PATHS, help="Largest number of distinct routes (default 8)")
    parser.add_argument("--separation", type=float,
                        help="Minimum distance (Å) between exit points of distinct routes (default 2 x bulk radius)")
    parser.add_argument("--route-separation", type=float,
                        help="Minimum distance (Å) between distinct routes beyond their junction (default bulk radius)")
    parser.add_argument("--axis", help="x, y, z or x,y,z (use --axis=-1,0,0 for a leading minus): search each side of "
                                       "the plane through the start normal to this axis separately (half-spaces), "
                                       "label branches + or - and record a side without a route as unresolved")
    parser.add_argument("--min-enclosure", type=float, default=0.75,
                        help="A free node from which more than 26 - ceil(26 x this) of its 26 straight lattice rays "
                             "reach bulk before an atom-blocked node counts as exposed and ends a route like bulk; "
                             "0 uses bulk only (default 0.75: at least 7 rays)")
    parser.add_argument("--include-hetero", action="store_true",
                        help="Treat hetero atoms (ligands; for MD input everything outside the protein selection, "
                             "e.g. lipids and ions with --md-selection 'not resname TIP3') as obstacles")
    parser.add_argument("--lining-shell", type=float, default=1.4,
                        help="Lining shell (Å) between an atom surface and a sample's free sphere (default 1.4)")
    parser.add_argument("--relax-iterations", type=int, default=60,
                        help="Centre-line relaxation iterations; 0 keeps the straightened lattice route (default 60)")
    parser.add_argument("--max-grid-points", type=int, default=20_000_000, help="Grid allocation guard")
    parser.add_argument("--trajectory", help="Trajectory file (with a GRO/topology STRUCTURE): analyse --frames")
    parser.add_argument("--frames", default="0",
                        help="With --trajectory: comma-separated zero-based frame indices (default 0)")
    parser.add_argument("--pbc", choices=["check", "none"], default="check",
                        help="With --trajectory: periodic-boundary bond check (none for wrapped lipids)")
    parser.add_argument("--html", action="store_true",
                        help="Also write a self-contained interactive view <prefix>_paths_view.html (one layer per branch)")
    parser.add_argument("--no-scenes", action="store_true", help="Skip the PyMOL/ChimeraX/VMD scenes")
    parser.add_argument("--dpi", type=int, default=200,
                        help="PNG figure resolution in dots per inch (default 200; use --dpi 400 or more for print/publication)")
    parser.set_defaults(func=run_paths_command)
    return parser


def run_paths_command(args):
    """Run ``crevice paths`` from parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed ``paths`` arguments.
    """
    from pathlib import Path
    from .cli import _load_input
    frames = []
    if args.trajectory:
        from .trajectory import load_trajectory_report
        for index in [int(v) for v in args.frames.split(",") if v.strip()]:
            trajectory, _ = load_trajectory_report(args.structure, args.trajectory, selection=args.md_selection,
                                                   start=index, stop=index + 1, pbc=args.pbc)
            frames.append((f"_f{index}", trajectory.frames[0]))
        prefix = args.prefix or Path(args.trajectory).stem
    else:
        frame, source, _ = _load_input(args)
        frames.append(("", frame))
        prefix = args.prefix or Path(str(source.path)).stem
    for suffix, frame in frames:
        result = find_paths(frame, seed=args.seed, seed_residues=args.seed_residues, region=args.region,
                            probe_radius=args.probe_radius, spacing=args.spacing, bulk_radius=args.bulk_radius,
                            seed_tolerance=args.seed_tolerance, max_paths=args.max_paths, separation=args.separation,
                            route_separation=args.route_separation, axis=args.axis,
                            include_hetero=args.include_hetero, lining_shell=args.lining_shell,
                            relax_iterations=args.relax_iterations, min_enclosure=args.min_enclosure,
                            max_grid_points=args.max_grid_points, workers=args.workers)
        result.settings["frame_index"] = frame.frame_index
        write_paths_outputs(result, args.output_dir, prefix=prefix + suffix, frame=frame, dpi=args.dpi,
                            scenes=not args.no_scenes)
        if args.html:
            paths_view(result, frame).save_html(Path(args.output_dir) / f"{prefix}{suffix}_paths_view.html")
        widest = result.branches[0] if result.branches else None
        print(f"paths{suffix} status={result.status} branches={len(result.branches)} junctions="
              f"{len(result.junctions)}" + (f" widest_bottleneck={widest.bottleneck_radius:.3f} length="
                                              f"{widest.length:.2f}" if widest else "")
              + "".join(f" unresolved_side={u['side']}" for u in result.unresolved if u.get("side")))
        selection = result.settings.get("probe_selection")
        if selection:
            choice = ", ".join(f"{side} {v['probe_A']:g} Å" if v["probe_A"] is not None else f"{side} none"
                               for side, v in selection["chosen_by_side"].items())
            print(f"paths{suffix} probe=auto chosen: {choice} ({selection['reason']})")

"""Lateral exits for capped channel ends (opt-in, ``lateral_exits=True``).

A through-channel profile normally requires both ends to continue straight
along the channel axis, unobstructed, beyond the atom envelope. Some channels
are capped: an extracellular or intracellular domain blocks the straight exit,
and ions leave sideways through portals (for example the K2P channel TWIK-1).
With ``lateral_exits=True`` a blocked end is accepted only if a lateral exit
exists, defined here as follows.

Definitions (all distances in Å, atoms as van der Waals spheres):

* **Start.** The last sampled centre at the capped end of the axial profile
  (the refined mouth section, where lateral enclosure is lost).
* **Half-space.** The exit may only use points on the far side of the plane
  through the start, normal to the channel axis. It therefore cannot turn back
  into the enclosed channel it came from.
* **Free path.** A path of straight segments between nodes of a cubic lattice
  (``exit_spacing``) anchored at the start and aligned with the channel frame,
  26-neighbour steps. Every segment must clear the enclosure probe
  (``enclosure_radius``), checked analytically, as for the axial sections.
* **Bulk.** Space reachable by a large rolling probe (``bulk_radius``, default
  6 Å, the same default as the rolling-probe cast's outer probe) that is
  connected to the outer faces of a padded box around the atoms, restricted to
  the same half-space; the far face at the start plane does not count as
  outside. A point is in bulk when it lies inside a sphere swept by such a
  probe. This is the rolling-probe molecular envelope of
  :func:`crevice.rolling.rolling_probe_cast`, computed on a coarser grid
  (``bulk_spacing``).
* **Exits.** Free paths from the start to bulk, found widest first: a path's
  width is its smallest node and segment clearance, and a maximin
  (widest-path) Dijkstra labels every lattice node with its widest route.
  Widths are counted **beyond the mouth sample's inscribed sphere** (the
  sphere of the start's own clearance, already measured by the axial
  profile); otherwise every route would share the mouth's width and all
  exits would tie. Nodes and segments inside that sphere must still clear
  the enclosure probe. Equally wide routes are ranked by length (the shorter
  wins), then by lattice order, so the result is deterministic. The search
  continues after the first bulk node, and every **distinct** exit is
  reported (:func:`lateral_exit_legs`): a further exit counts when its exit
  point is more than ``2 * bulk_radius`` from every exit already kept (no
  single bulk-probe position touches both) and the parts of the two routes
  farther than ``bulk_radius`` from the start never come within
  ``bulk_radius`` of each other (they leave through different openings).
  Symmetric portals, for example the two C2-related portals of a K2P dimer,
  therefore give one leg each. At most ``MAX_EXIT_LEGS`` (8) legs per end.
  Each leg reports its whole-route bottleneck (including the mouth, as
  before), ``ranking_width_raw_clearance_A`` (the lattice-route width beyond
  the mouth sphere used for ranking) and ``beyond_mouth_min_raw_clearance_A``
  (the smallest clearance of the leg's samples outside that sphere).
* **Leg.** The lattice path is shortened by straight cuts, each no narrower
  than the stretch of path it replaces, then resampled every 0.25 Å. Each leg sample reports the
  distance to the nearest atom surface (a free-sphere radius, like a tunnel
  radius). These are **not** planar cross-section radii, so they are kept
  separate from the axial profile samples and its minimum.

If no free path reaches bulk (for example the cap closes the channel) the
channel stays unresolved. Nothing is interpolated.

Outputs and display. :func:`crevice.sections.connected_profile` stores the
legs of a capped end in ``metadata["exits"][end]["legs"]`` (``ExitLeg.to_dict``
plus the rank ``leg``: points, radii, bottleneck and its position and arc
length, exit point, length, lattice vertices), with ``leg_count`` and the
widest leg's ``bottleneck_radius_A`` beside them.
``crevice.presentation.scene_geometry`` copies every leg's vertices into the
scene JSON as ``exit_paths`` (with ``end`` and ``leg``), and the generated
viewer scenes draw each leg as a thin centre-line tube: PyMOL CGO object
``crevice_exits`` (radius 0.12 Å), VMD graphics colour id 23, and ChimeraX
cylinders appended to ``<prefix>_mouths.bild``, all in blue RGB (0.16, 0.47,
0.84), next to the orange mouth rings. The tube shows the route only, not the
portal width. The profile radius plots draw each leg's radii as a
bluish-green (``#009e73``) dash-dot continuation beyond the mouth against path
length along the leg. No text or labels are drawn. None of this is written
unless a lateral exit exists.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import heapq
import math

from .geometry import add, scale, sub, dot, distance
from .radii import atom_vdw_radius

EXIT_DEFINITION = ("widest free path (maximin segment clearance > enclosure probe) from the capped "
                   "end of the axial profile into bulk, restricted to the half-space beyond that end; "
                   "bulk = space swept by a rolling probe of bulk_radius connected to the outside")
LEG_RADIUS_DEFINITION = ("distance from each leg sample to the nearest atom van der Waals surface, "
                         "minus probe_radius; a free-sphere (tunnel-like) radius, not a planar section radius")


@dataclass
class ExitLeg:
    """One lateral exit leg from a capped channel end to bulk.

    Attributes
    ----------
    end : str
        ``"lower"`` or ``"upper"``: the capped end of the axial profile.
    points : tuple of Coord
        Samples every 0.25 Å along the straightened route, Å, from the mouth
        sample's centre to the exit point in bulk.
    clearances : tuple of float
        Distance of each sample to the nearest atom surface, Å (a free-sphere
        radius, like a tunnel radius, not a planar section radius).
    bottleneck_clearance : float
        Smallest straight-segment clearance of the route, Å.
    bottleneck_position : Coord
        Sample of smallest clearance.
    length : float
        Route length, Å.
    lattice_nodes_explored : int
        Search nodes labelled when this leg was found.
    metadata : dict
        Rank ``leg`` (1 = widest), ``bottleneck_arc_length_from_start_A``,
        ``arc_lengths_A``, exit-grid and bulk settings, axial and lateral offsets
        of the exit point (Å), ``segment_clearances_A``, the straight-route
        ``vertices`` and, from :func:`lateral_exit_legs`,
        ``ranking_width_raw_clearance_A``, ``beyond_mouth_min_raw_clearance_A``
        and ``mouth_sphere_radius_A``.

    Notes
    -----
    ``to_dict(probe_radius)`` gives the JSON record stored in
    ``metadata["exits"][end]["legs"]`` of the profile, with radii reduced by
    ``probe_radius``.
    """
    end: str
    points: tuple
    clearances: tuple
    bottleneck_clearance: float
    bottleneck_position: tuple
    length: float
    lattice_nodes_explored: int
    metadata: dict = field(default_factory=dict)

    def to_dict(self, probe_radius: float = 0.0) -> dict:
        """JSON record of the leg, with radii reduced by the measurement probe.

        Parameters
        ----------
        probe_radius : float, default 0.0
            Measurement probe radius, Å, subtracted from every clearance (never
            below zero) for ``radii_A`` and ``bottleneck_radius_A``.

        Returns
        -------
        dict
            ``type`` (``"lateral"``), ``end``, ``points``, ``radii_A``,
            ``raw_clearances_A``, ``bottleneck_radius_A``,
            ``bottleneck_raw_clearance_A``, ``bottleneck_position``, ``exit_point``,
            ``length_A``, ``lattice_nodes_explored``, ``definition``,
            ``radius_definition`` and every ``metadata`` entry.
        """
        return {"type": "lateral", "end": self.end,
                "points": [list(p) for p in self.points],
                "radii_A": [max(0.0, c-probe_radius) for c in self.clearances],
                "raw_clearances_A": list(self.clearances),
                "bottleneck_radius_A": max(0.0, self.bottleneck_clearance-probe_radius),
                "bottleneck_raw_clearance_A": self.bottleneck_clearance,
                "bottleneck_position": list(self.bottleneck_position),
                "exit_point": list(self.points[-1]),
                "length_A": self.length,
                "lattice_nodes_explored": self.lattice_nodes_explored,
                "definition": EXIT_DEFINITION,
                "radius_definition": LEG_RADIUS_DEFINITION,
                **self.metadata}


class BulkEnvelope:
    """Rolling large-probe bulk region in the half-space beyond one channel end.

    A grid (spacing ``spacing``) in the channel frame covers the atoms padded by
    the largest atom radius, ``bulk_radius`` and two grid steps, cut at the start
    plane ``t_limit`` so only the half-space beyond the capped end is included.
    Probe centres with clearance >= ``bulk_radius`` that connect continuously to
    the box faces (except the face on the start plane) are bulk centres; a point
    is in bulk when it is within ``bulk_radius`` of a bulk centre.

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms (van der Waals spheres).
    origin : Coord
        Channel axis origin, Å.
    basis : sequence of Coord
        Rows ``u``, ``v``, ``axis`` of the channel frame.
    t_limit : float
        Axial coordinate of the start plane, Å.
    sign : {-1, 1}
        Side of the start plane that is searched (+1 = increasing ``t``).
    spacing : float, default 1.0
        Grid spacing of the envelope, Å (``bulk_spacing``).
    bulk_radius : float, default 6.0
        Rolling bulk-probe radius, Å.

    Attributes
    ----------
    distance : numpy.ndarray
        Distance (Å) from each grid point to the nearest bulk centre.
    grid_points : int
        Size of the envelope grid.

    Raises
    ------
    ValueError
        If the grid would exceed 20 million points.
    """

    def __init__(self, atoms, origin, basis, t_limit, sign, *, spacing=1.0, bulk_radius=6.0):
        import numpy as np
        from scipy import ndimage
        from .rolling import _SphereQueries, _connected
        self.origin_world = np.asarray(origin, dtype=float)
        self.basis = np.asarray(basis, dtype=float)  # rows: u, v, axis
        xyz = np.asarray([a.coord for a in atoms], dtype=float)
        local = (xyz-self.origin_world) @ self.basis.T
        radii = np.asarray([atom_vdw_radius(a) for a in atoms], dtype=float)
        pad = float(radii.max())+bulk_radius+2*spacing
        low = local.min(axis=0)-pad
        high = local.max(axis=0)+pad
        if sign > 0:
            low[2] = t_limit
        else:
            high[2] = t_limit
        dims = tuple(int(v) for v in np.ceil((high-low)/spacing).astype(int)+1)
        if math.prod(dims) > 20_000_000:
            raise ValueError("Bulk-envelope grid for the lateral exit exceeds 20 million points; increase exit bulk spacing")
        self.low, self.spacing, self.dims, self.sign = low, spacing, dims, sign
        self.bulk_radius = bulk_radius
        queries = _SphereQueries(local, radii)
        ijk = np.indices(dims).reshape(3, -1).T
        clearance = np.empty(len(ijk))
        for start in range(0, len(ijk), 262144):
            clearance[start:start+262144] = queries.points(low+ijk[start:start+262144]*spacing)[0]
        clearance = clearance.reshape(dims)
        labels = _connected(clearance >= bulk_radius, clearance, low, spacing, bulk_radius, queries)
        faces = [labels[0], labels[-1], labels[:, 0], labels[:, -1], labels[:, :, -1 if sign > 0 else 0]]
        outside = np.unique(np.concatenate([f.ravel() for f in faces]))
        outside = outside[outside > 0]
        centres = np.isin(labels, outside)
        self.distance = ndimage.distance_transform_edt(~centres, sampling=spacing)
        self.grid_points = int(math.prod(dims))

    def local(self, point):
        """Coordinates of a point in the channel frame.

        Parameters
        ----------
        point : Coord
            Cartesian point, Å.

        Returns
        -------
        numpy.ndarray, shape (3,)
            ``(u, v, axis)`` components relative to the channel origin, Å.
        """
        import numpy as np
        return (np.asarray(point, dtype=float)-self.origin_world) @ self.basis.T

    def depth(self, point) -> float:
        """Distance beyond the swept bulk region (0 inside bulk); outside the box counts as bulk.

        Parameters
        ----------
        point : Coord
            Cartesian point, Å (rounded to the nearest envelope grid point).

        Returns
        -------
        float
            ``max(0, distance to the nearest bulk centre - bulk_radius)``, Å.
        """
        import numpy as np
        index = np.rint((self.local(point)-self.low)/self.spacing).astype(int)
        if np.any(index < 0) or np.any(index >= np.asarray(self.dims)):
            return 0.0
        return max(0.0, float(self.distance[tuple(index)])-self.bulk_radius)


#: Largest number of distinct lateral exit legs reported per capped end.
MAX_EXIT_LEGS = 8
DISTINCT_EXIT_DEFINITION = (
    "exits are ranked by their width beyond the mouth sample's inscribed sphere, then by route length; "
    "a further exit is distinct when its exit point is more than exit_separation_A from every kept exit "
    "point (default twice the bulk probe radius: no single bulk-probe position touches both) and the "
    "parts of the two routes farther than route_separation_A from the start (default the bulk probe "
    "radius) never come within route_separation_A of each other")

_STEPS = tuple((i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0))


def lateral_exit(geometry, start, sign, *, end, spacing=0.5, bulk_radius=6.0, bulk_spacing=1.0,
                 max_nodes=1_000_000):
    """Return the widest free path from ``start`` (a Section) to bulk, or None.

    ``geometry`` is the channel's :class:`crevice.sections.SectionGeometry`;
    its enclosure probe is the clearance every segment must exceed. This is
    the first (widest) leg of :func:`lateral_exit_legs`, found with
    ``max_legs=1``.

    Parameters
    ----------
    geometry : crevice.sections.SectionGeometry
        The channel's section geometry.
    start : crevice.sections.Section
        The capped end's last profile sample.
    sign : {-1, 1}
        Direction of the capped end along the axis.
    end : {"lower", "upper"}
        Label copied into the leg.
    spacing : float, default 0.5
        Exit-search lattice spacing, Å (``--exit-spacing``).
    bulk_radius : float, default 6.0
        Rolling bulk-probe radius, Å (``--exit-bulk-radius``).
    bulk_spacing : float, default 1.0
        Bulk-envelope grid spacing, Å.
    max_nodes : int, default 1000000
        Lattice nodes the search may label.

    Returns
    -------
    ExitLeg or None
        The widest leg, or ``None`` when no free path reaches bulk.
    """
    legs, _ = lateral_exit_legs(geometry, start, sign, end=end, spacing=spacing, bulk_radius=bulk_radius,
                                bulk_spacing=bulk_spacing, max_nodes=max_nodes, max_legs=1)
    return legs[0] if legs else None


def lateral_exit_legs(geometry, start, sign, *, end, spacing=0.5, bulk_radius=6.0, bulk_spacing=1.0,
                      max_nodes=1_000_000, max_legs=MAX_EXIT_LEGS, separation=None, route_separation=None):
    """All distinct lateral exits from ``start`` to bulk, widest first.

    Parameters
    ----------
    geometry : crevice.sections.SectionGeometry
        The channel's section geometry; ``geometry.probe`` (the enclosure
        probe, Å) is the clearance every lattice node and straight segment
        must exceed, and its frame (u, v, axis) orients the lattice.
    start : crevice.sections.Section
        The capped end's last profile sample; every leg starts at its centre.
    sign : {-1, 1}
        Direction of the capped end along the axis; legs stay in the
        half-space ``sign * (t - start.t) >= 0``.
    end : {"lower", "upper"}
        Label copied into each leg.
    spacing : float
        Lattice spacing (Å) of the maximin search (CLI ``--exit-spacing``).
        Finer lattices find narrower gaps and straighter routes; cost grows
        as ``spacing**-3``.
    bulk_radius : float
        Rolling bulk-probe radius (Å, ``--exit-bulk-radius``) that defines
        bulk (:class:`BulkEnvelope`); also sets the distinct-exit distances.
    bulk_spacing : float
        Grid spacing (Å) of the bulk envelope.
    max_nodes : int
        Lattice nodes the search may label. If it is exceeded before the
        first exit is found, ``ValueError``; after that, the search stops and
        ``info["search_truncated"]`` is True.
    max_legs : int
        Largest number of legs returned (``info["legs_truncated"]`` records
        a cut).
    separation : float, optional
        Minimum distance (Å) between exit points of distinct legs; default
        ``2 * bulk_radius``.
    route_separation : float, optional
        Minimum distance (Å) between the outer parts (farther than this from
        the start) of the routes of distinct legs; default ``bulk_radius``.

    Returns
    -------
    legs : list of ExitLeg
        Widest beyond the mouth first, then shortest. ``legs[0]`` is exactly
        the leg :func:`lateral_exit` returns. Empty if no free path reaches
        bulk.
    info : dict
        ``distinct_exit_definition``, ``exit_separation_A``,
        ``route_separation_A``, ``max_legs``, ``legs_truncated``,
        ``search_truncated``, ``lattice_nodes_explored``, ``terminals_examined``.

    Notes
    -----
    A maximin (widest-path) Dijkstra with a route-length tie-break runs over
    the whole reachable lattice. Widths ignore clearances inside the mouth
    sample's inscribed sphere (all routes share it). Bulk nodes are terminal
    and are popped in order of decreasing width, then increasing route
    length, so each terminal carries its widest (then shortest) route. A
    terminal gives a new leg only if (a) its exit point is more than
    ``separation`` from every accepted exit point (default: the bulk probe's
    diameter, so no single bulk-probe position touches both), and (b) the
    parts of its route and of every accepted route farther than
    ``route_separation`` from the start stay more than ``route_separation``
    apart (default ``bulk_radius``: the routes leave through different
    openings). Symmetric portals, such as the two C2-related portals of a K2P
    dimer, therefore give one leg each instead of a tie-break.
    """
    spatial, probe = geometry.spatial, geometry.probe
    u, v, axis = geometry.u, geometry.v, geometry.direction
    separation = 2.0*bulk_radius if separation is None else separation
    route_separation = bulk_radius if route_separation is None else route_separation
    envelope = BulkEnvelope(geometry.atoms, geometry.origin, (u, v, axis), start.t, sign,
                            spacing=bulk_spacing, bulk_radius=bulk_radius)
    origin = start.center
    w = scale(axis, sign)

    def position(node):
        i, j, k = node
        return add(origin, add(scale(u, i*spacing), add(scale(v, j*spacing), scale(w, k*spacing))))

    clear = {}

    def clearance(node):
        if node not in clear:
            clear[node] = spatial.nearest_surface(position(node))[0]
        return clear[node]

    root = (0, 0, 0)
    # The start sample's own inscribed sphere is the channel mouth, already
    # measured by the axial profile. Clearances inside it are not used to rank
    # exits (otherwise every route would share the mouth's width and tie);
    # every node and segment must still clear the enclosure probe.
    mouth_sphere = clearance(root)

    def in_mouth(node):
        return distance(position(node), origin) <= mouth_sphere

    width = {root: math.inf}
    parent = {root: None}
    counter = 0
    travelled = {root: 0.0}
    # Widest first; among equally wide routes the shorter one (lexicographic
    # maximin-then-length Dijkstra, exact on the lattice).
    heap = [(-round(width[root], 10), 0.0, counter, root)]
    terminals = []
    info = {"distinct_exit_definition": DISTINCT_EXIT_DEFINITION, "exit_separation_A": separation,
            "route_separation_A": route_separation, "max_legs": max_legs, "legs_truncated": False, "search_truncated": False}
    while heap:
        key, walked, _, node = heapq.heappop(heap)
        current = width[node]
        if -key < round(current, 10) or (-key == round(current, 10) and walked > travelled[node] + 1e-9):
            continue
        p = position(node)
        if envelope.depth(p) == 0.0:
            # A bulk node: terminal, never expanded. Popped in order of
            # decreasing width, then increasing route length.
            terminals.append(node)
            if max_legs == 1:
                break
            continue
        if len(width) > max_nodes:
            if not terminals:
                raise ValueError(f"Lateral-exit search exceeded {max_nodes} lattice nodes; increase exit_spacing")
            info["search_truncated"] = True
            break
        for step in _STEPS:
            other = (node[0]+step[0], node[1]+step[1], node[2]+step[2])
            if other[2] < 0:
                continue
            c = clearance(other)
            if c <= probe:
                continue
            length = spacing*math.sqrt(step[0]**2+step[1]**2+step[2]**2)
            inside = in_mouth(other)
            new = current if inside else min(current, c)
            # (clear[node]+c-length)/2 is a lower bound on the segment clearance.
            bound = (clear[node]+c-length)/2
            if (bound <= probe) if inside else (bound < new):
                segment = spatial.segment_clearance(p, position(other))
                if segment <= probe:
                    continue
                if not inside:
                    new = min(new, segment)
            rounded, total = round(new, 10), travelled[node] + length
            known = round(width[other], 10) if other in width else -math.inf
            if rounded > known or (rounded == known and total < travelled[other] - 1e-9):
                width[other] = new
                travelled[other] = total
                parent[other] = node
                counter += 1
                heapq.heappush(heap, (-rounded, total, counter, other))
    info["lattice_nodes_explored"] = len(width)
    info["terminals_examined"] = len(terminals)
    # Distinct exits, in terminal order (widest, then shortest route first).
    legs, accepted = [], []
    for node in terminals:
        p = position(node)
        if any(distance(p, other_exit) <= separation for other_exit, _ in accepted):
            continue
        route = []
        walk = node
        while walk is not None:
            route.append(position(walk))
            walk = parent[walk]
        route.reverse()
        route[0] = origin
        outer = [q for q in route if distance(q, origin) > route_separation]
        if not outer or any(min(distance(q, r) for q in outer for r in other_outer) <= route_separation
                            for _, other_outer in accepted if other_outer):
            continue
        if len(legs) >= max_legs:
            info["legs_truncated"] = True
            break
        leg = _leg_from_route(spatial, route, end, clearance(root), len(width), origin, axis,
                              spacing=spacing, bulk_radius=bulk_radius, bulk_spacing=bulk_spacing,
                              bulk_grid_points=envelope.grid_points, rank=len(legs)+1)
        leg.metadata["ranking_width_raw_clearance_A"] = width[node] if math.isfinite(width[node]) else None
        outer_samples = [c for q, c in zip(leg.points, leg.clearances) if distance(q, origin) > mouth_sphere]
        leg.metadata["beyond_mouth_min_raw_clearance_A"] = min(outer_samples) if outer_samples else None
        leg.metadata["mouth_sphere_radius_A"] = mouth_sphere
        legs.append(leg)
        accepted.append((p, outer))
    return legs, info


def _leg_from_route(spatial, nodes, end, root_clearance, explored, origin, axis, *, spacing, bulk_radius,
                    bulk_spacing, bulk_grid_points, rank):
    """Straighten, resample and measure one lattice route (start ... bulk node)."""
    edges = [spatial.segment_clearance(a, b) for a, b in zip(nodes, nodes[1:])]
    bottleneck = min(edges) if edges else root_clearance
    # Straight cuts, each no narrower than the stretch of route it replaces:
    # fewer lattice kinks, and neither the whole-route bottleneck nor the
    # width beyond the mouth used to rank exits can drop.
    kept, i = [nodes[0]], 0
    while i < len(nodes)-1:
        j = len(nodes)-1
        while j > i+1 and spatial.segment_clearance(nodes[i], nodes[j]) < min(edges[i:j])-1e-9:
            j -= 1
        kept.append(nodes[j])
        i = j
    samples = [kept[0]]
    for a, b in zip(kept, kept[1:]):
        n = max(1, int(math.ceil(distance(a, b)/0.25)))
        samples += [add(a, scale(sub(b, a), s/n)) for s in range(1, n+1)]
    segment_clearances = [spatial.segment_clearance(a, b) for a, b in zip(kept, kept[1:])]
    exact = min(segment_clearances) if segment_clearances else root_clearance
    worst = min(range(len(samples)), key=lambda s: spatial.nearest_surface(samples[s])[0])
    length = sum(distance(a, b) for a, b in zip(kept, kept[1:]))
    arc = [0.0]
    for a, b in zip(samples, samples[1:]):
        arc.append(arc[-1]+distance(a, b))
    return ExitLeg(end, tuple(samples), tuple(spatial.nearest_surface(p)[0] for p in samples), exact,
                   samples[worst], length, explored,
                   metadata={"leg": rank,
                             "bottleneck_arc_length_from_start_A": arc[worst],
                             "arc_lengths_A": arc,
                             "exit_grid_spacing_A": spacing, "bulk_radius_A": bulk_radius,
                             "bulk_grid_spacing_A": bulk_spacing, "bulk_grid_points": bulk_grid_points,
                             "axial_offset_of_exit_point_A": dot(sub(samples[-1], origin), axis),
                             "lateral_offset_of_exit_point_A": math.dist(
                                 samples[-1], add(origin, scale(axis, dot(sub(samples[-1], origin), axis)))),
                             "segment_clearances_A": segment_clearances,
                             "vertices": [list(p) for p in kept]})


EXIT_CAST_DEFINITION = (
    "probe-swept cast of a lateral exit leg: on the channel cast's lattice, beyond the capped mouth "
    "plane and outside bulk (the leg's rolling bulk-probe envelope), core nodes inside the free sphere of "
    "a leg sample (distance to the sample <= its atom-surface clearance) whose own clearance is >= the "
    "enclosure probe, swept by (enclosure probe - min_radius); swept nodes need clearance >= min_radius, "
    "must lie outside bulk, must connect (continuous edge test) to a node on the leg "
    "centre line, and a node shared by several legs belongs to the leg with the nearest centre-line sample. "
    "Reported separately; never part of the axial channel volume")


def lateral_exit_casts(frame, profile, cast):
    """Probe-swept casts of the lateral exit legs of a capped channel end.

    Each lateral exit leg of ``profile`` (``metadata["exits"][end]["legs"]``)
    becomes one component of kind ``"lateral_exit"`` on the same lattice as
    the channel cast, so that the two join at the capped mouth. The channel
    cast, its volume and every profile value are unchanged; the leg casts are
    reported separately.

    Definition (all distances in Å, atoms as van der Waals spheres; see
    ``EXIT_CAST_DEFINITION``):

    * Lattice: the channel cast's ``grid_origin``, ``grid_basis`` and spacing.
      Only lattice planes beyond the capped mouth plane are used (axial index
      below 0 for the lower end, above the last channel plane for the upper
      end), so leg casts never overlap the channel cast.
    * Bulk is excluded: a node in bulk (inside the space swept by the leg's
      rolling bulk probe, ``bulk_radius_A``, :class:`BulkEnvelope`) is never
      part of a leg cast. A leg cast is therefore the portal space between the
      capped mouth and bulk solvent, not the solvent outside the protein.
    * Core: nodes inside the free sphere of at least one leg sample (distance
      to the sample at most that sample's atom-surface clearance) whose own
      clearance is at least the enclosure probe ``q`` (the cast's
      ``enclosure_radius_A``). Where the leg widens towards bulk the spheres
      grow, so the cast follows the opening.
    * Sweep: core nodes are swept by ``q - min_radius``, as for the channel
      cast; swept nodes need clearance ``>= min_radius``.
    * Connectivity: nodes are clustered with the channel cast's continuous
      edge test (:class:`crevice.grid.GridConnectivity`); a cluster is kept
      only if it contains the lattice node nearest a leg sample.
    * Shared nodes: a node reached from several legs belongs to the leg whose
      nearest centre-line sample is closest (lower rank on ties), so leg
      volumes never double count.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile and cast (the cast's own atom selection is
        applied).
    profile : PoreProfile
        Resolved ``lateral_exits`` profile.
    cast : VoidCast
        Channel cast of ``profile`` (``export_mode ==
        "connected_section_fill"``).

    Returns
    -------
    tuple of VoidComponent
        One component per leg that has at least one node, ordered by end
        then leg rank; ids continue after the channel components. Metadata:
        ``end``, ``leg``, ``definition``, ``core_point_count``,
        ``enclosure_radius_A``, ``probe_sweep_radius_A``, grid origin, basis
        and spacing. Empty when the profile has no lateral exit, the cast is
        not a connected channel cast, or the cast is focus-cropped.

    Examples
    --------
    >>> from crevice.models import PoreProfile
    >>> lateral_exit_casts(None, PoreProfile((0, 0, 0), (0, 0, 1), ()), None)
    ()
    """
    import numpy as np
    from .grid import GridConnectivity
    from .models import ChannelPoint, VoidComponent
    from .spatial import SpatialIndex
    from .voids import _cluster_point_indices
    exits = profile.metadata.get("exits") or {}
    legs = [(label, leg) for label in ("lower", "upper") if (exits.get(label) or {}).get("type") == "lateral"
            for leg in exits[label].get("legs", ())]
    if (not legs or cast is None or cast.metadata.get("export_mode") != "connected_section_fill"
            or cast.metadata.get("focus_points")):
        return ()
    meta = cast.metadata
    spacing, min_radius = cast.spacing, cast.min_radius
    q = meta["enclosure_radius_A"]
    sweep = max(0.0, q-min_radius)
    origin = np.asarray(meta["grid_origin"], dtype=float)
    basis = np.asarray(meta["grid_basis"], dtype=float)
    nz = int(meta["grid_dimensions"][2])
    low = meta["mouth_planes_A"][0]
    spatial = SpatialIndex(frame.selected_atoms(include_hydrogen=meta.get("include_hydrogen", False),
                                                include_hetero=meta.get("include_hetero", True)))
    connectivity = GridConnectivity(spatial, tuple(origin), spacing, min_radius, basis=tuple(map(tuple, basis)))
    clear = {}

    def clearance(idx):
        if idx not in clear:
            clear[idx] = spatial.nearest_surface(connectivity.position(idx))
        return clear[idx][0]

    envelopes = {}

    def envelope(label, leg):
        if label not in envelopes:
            point = profile.points[0] if label == "lower" else profile.points[-1]
            envelopes[label] = BulkEnvelope(spatial.atoms, profile.axis_origin, basis, point.t,
                                            -1 if label == "lower" else 1,
                                            spacing=leg.get("bulk_grid_spacing_A", 1.0),
                                            bulk_radius=leg.get("bulk_radius_A", 6.0))
        return envelopes[label]

    def beyond(label, k):
        return k < 0 if label == "lower" else k > nz-1

    reach = int(math.ceil(sweep/spacing))
    offsets = np.asarray([d for d in np.ndindex(2*reach+1, 2*reach+1, 2*reach+1)], dtype=int)-reach
    offsets = offsets[(offsets**2).sum(1)*spacing**2 <= sweep**2+1e-12]
    per_leg, seeds, lines = [], [], []
    for label, leg in legs:
        samples = np.asarray(leg["points"], dtype=float)
        radii = np.asarray(leg["raw_clearances_A"], dtype=float)
        local = (samples-origin) @ basis.T/spacing
        lines.append(samples)
        bulk = envelope(label, leg)
        inside = {}

        def enclosed(idx):
            if idx not in inside:
                inside[idx] = beyond(label, idx[2]) and bulk.depth(connectivity.position(idx)) > 0.0
            return inside[idx]
        core = set()
        for centre, r in zip(local, radii):
            if r <= 0:
                continue
            n = r/spacing
            lo, hi = np.floor(centre-n).astype(int), np.ceil(centre+n).astype(int)
            box = np.stack(np.meshgrid(*[np.arange(a, b+1) for a, b in zip(lo, hi)], indexing="ij"), -1).reshape(-1, 3)
            box = box[((box-centre)**2).sum(1) <= n*n+1e-9]
            for idx in map(tuple, box.tolist()):
                if idx not in core and enclosed(idx) and clearance(idx) >= q:
                    core.add(idx)
        swept = set()
        for idx in core:
            for d in offsets.tolist():
                other = (idx[0]+d[0], idx[1]+d[1], idx[2]+d[2])
                if enclosed(other) and clearance(other) >= min_radius:
                    swept.add(other)
        nearest = {tuple(np.rint(p).astype(int).tolist()) for p in local}
        per_leg.append((core, swept))
        seeds.append({idx for idx in nearest if idx in swept})
    # A node reached from several legs goes to the nearest centre line.
    owner = {}
    for rank, (_, swept) in enumerate(per_leg):
        for idx in swept:
            if idx in owner:
                p = np.asarray(connectivity.position(idx))
                mine = float(np.min(np.linalg.norm(lines[rank]-p, axis=1)))
                theirs = float(np.min(np.linalg.norm(lines[owner[idx]]-p, axis=1)))
                if mine < theirs-1e-9:
                    owner[idx] = rank
            else:
                owner[idx] = rank
    components = []
    next_id = max((c.id for c in cast.components), default=0)+1
    for rank, ((label, leg), (core, _)) in enumerate(zip(legs, per_leg)):
        nodes = {idx for idx, o in owner.items() if o == rank}
        clusters = _cluster_point_indices(nodes, edge_allowed=lambda a, b: connectivity.allows(
            a, b, clearance(a), clearance(b)))
        kept = sorted(idx for cluster in clusters if seeds[rank] & set(cluster) for idx in cluster)
        if not kept:
            continue
        points = []
        for i, idx in enumerate(kept):
            raw, atom = clear[idx]
            points.append(ChannelPoint(i, connectivity.position(idx), raw, raw, low+idx[2]*spacing,
                                       atom.serial, atom.residue_key.label))
        centre = tuple(float(v) for v in np.mean([p.position for p in points], axis=0))
        components.append(VoidComponent(
            next_id, "lateral_exit", centre, len(points)*spacing**3, tuple(points),
            nearest_residues=tuple(sorted({p.nearest_residue for p in points})),
            metadata={"end": label, "leg": leg.get("leg", rank+1), "definition": EXIT_CAST_DEFINITION,
                      "core_point_count": len(core & set(kept)), "enclosure_radius_A": q,
                      "probe_sweep_radius_A": sweep, "grid_origin": tuple(origin.tolist()),
                      "grid_basis": tuple(map(tuple, basis.tolist())), "grid_spacing": spacing,
                      "leg_bottleneck_radius_A": leg.get("bottleneck_radius_A"),
                      "leg_length_A": leg.get("length_A")}))
        next_id += 1
    return tuple(components)


def exit_cast_summary(components) -> list[dict]:
    """JSON rows describing lateral exit casts, one per leg.

    Parameters
    ----------
    components : sequence of VoidComponent
        From :func:`lateral_exit_casts`.

    Returns
    -------
    list of dict
        ``component_id``, ``end``, ``leg``, ``kind`` (``"lateral_exit"``),
        ``volume_A3``, ``point_count``, ``center``, ``min_clearance_A``,
        ``max_clearance_A``, ``core_point_count``, ``nearest_residues`` and
        ``definition``. Volumes are separate from, and never added to, the
        axial channel volume.

    Examples
    --------
    >>> exit_cast_summary(())
    []
    """
    return [{"component_id": c.id, "end": c.metadata["end"], "leg": c.metadata["leg"], "kind": c.kind,
             "volume_A3": c.volume, "point_count": c.point_count, "center": list(c.center),
             "min_clearance_A": c.min_radius, "max_clearance_A": c.max_radius,
             "core_point_count": c.metadata["core_point_count"],
             "nearest_residues": list(c.nearest_residues), "definition": c.metadata["definition"]}
            for c in components]

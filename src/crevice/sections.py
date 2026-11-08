"""Connected transverse free space in an atom-sphere model.

A section is eligible only when its grid component cannot reach the transverse
search boundary. Both nodes and edges exclude probe-inflated atom spheres.
This is a finite-resolution enclosure test, not a continuum proof. A single
axis describes only passages that are monotone along that direction.
"""

from __future__ import annotations

import contextvars
import functools
import math
from dataclasses import dataclass, replace

from .geometry import add, sub, scale, dot, distance, molecular_basis
from .models import Coord, ChannelPoint
from .radii import atom_vdw_radius
from .spatial import SpatialIndex


class ChannelResolutionError(ValueError):
    """A supplied region cannot be resolved as a unique through channel."""


class OriginClearanceError(ValueError):
    """A supplied channel origin lies too close to the atoms of this structure.

    Raised by :func:`connected_profile` when the clearance at the supplied
    ``origin`` is not greater than the enclosure and measurement probes. It is
    a :class:`ValueError` (an invalid input for a single structure); per-frame
    trajectory analysis (:func:`crevice.trajectory.analyze_trajectory`) records
    such a frame as unresolved instead of stopping.

    Attributes
    ----------
    clearance : float
        Clearance at the origin, Å (distance to the nearest atomic surface).
    """

    def __init__(self, message: str, clearance: float):
        """Store the message and the clearance at the origin (Å)."""
        super().__init__(message)
        self.clearance = float(clearance)


#: Section-centre refinement (:meth:`SectionGeometry.refine`) halves its step until it is below this, Å.
REFINEMENT_TOLERANCE_A = 1e-4
#: Guard on the number of improving moves of one section-centre refinement.
REFINEMENT_MAX_MOVES = 10_000
#: Radius, Å, within which a refinement vertex step looks for a three-atom clearance maximum.
REFINEMENT_VERTEX_RADIUS_A = 0.35
#: Guard on vertex / pattern-search rounds of one section-centre refinement.
REFINEMENT_MAX_ROUNDS = 32
_SQRT2 = math.sqrt(2.0)

#: The eight in-plane neighbour steps, in the order the flood fill visits them.
_STEPS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))


@dataclass
class Section:
    """One enclosed transverse section of the channel at axial coordinate ``t``.

    Attributes
    ----------
    t : float
        Axial coordinate of the plane, Å, relative to the axis origin.
    indices : tuple of tuple of int
        In-plane grid indices ``(i, j)`` of the connected free nodes.
    center : Coord
        Section centre, Å: the node of largest 3D atom-surface clearance (after
        :meth:`SectionGeometry.refine`, the refined local maximum).
    clearance : float
        Atom-surface clearance at ``center``, Å.
    window_bridged : bool
        True when the section reaches the plane's search square and is
        enclosed only within the bridging square twice as wide
        (:mod:`crevice.adaptive_channel`).
    """

    t: float
    indices: tuple[tuple[int, int], ...]
    center: Coord
    clearance: float
    window_bridged: bool = False


@functools.lru_cache(maxsize=64)
def _triples(n):
    """Index arrays ``(i, j, k)``, ``i < j < k < n``, in lexicographic order."""
    import itertools
    import numpy as np

    columns = np.array(list(itertools.combinations(range(n), 3)), dtype=np.intp).reshape(-1, 3)
    return tuple(np.ascontiguousarray(columns[:, m]) for m in range(3))


def _surface_values(points, xyz, radii):
    """Clearance (minimum of ``distance - radius`` over the given atoms) of each point.

    Same per-atom formula as :func:`crevice.spatial._point_surface_distances`.
    """
    import numpy as np

    points = np.asarray(points, dtype=float)
    dx = points[:, 0:1] - xyz[:, 0]
    dy = points[:, 1:2] - xyz[:, 1]
    dz = points[:, 2:3] - xyz[:, 2]
    return (np.sqrt(dx * dx + dy * dy + dz * dz) - radii).min(axis=1)


class _FrameState:
    """Probe-independent data of one section frame (atoms, origin, axis, basis).

    Holds the spatial index, the atoms projected on the frame (Python floats,
    as before, and NumPy columns) and memo tables of point and segment
    clearances. Memoised values are the exact values the spatial index
    returns, so sharing them changes no result.
    """

    _MEMO_LIMIT = 400_000

    def __init__(self, atoms, origin, direction, basis=None):
        import numpy as np

        self.u, self.v = basis or molecular_basis(direction, atoms)
        self.spatial = SpatialIndex(atoms)
        self.projected = [
            (
                dot(sub(a.coord, origin), self.u),
                dot(sub(a.coord, origin), self.v),
                dot(sub(a.coord, origin), direction),
                atom_vdw_radius(a),
            )
            for a in atoms
        ]
        columns = np.array(self.projected, dtype=float).reshape(-1, 4)
        self.px, self.py, self.pz, self.pr = (np.ascontiguousarray(columns[:, k]) for k in range(4))
        self._clear, self._segment = {}, {}

    def clearances(self, points):
        """Atom-surface clearances of ``points`` (memoised, batched)."""
        memo = self._clear
        missing = [p for p in dict.fromkeys(points) if p not in memo]
        if missing:
            if len(memo) > self._MEMO_LIMIT:
                # Clearing drops values this call still returns: compute every requested point.
                memo.clear()
                missing = list(dict.fromkeys(points))
            memo.update(zip(missing, self.spatial.nearest_surface_values(missing).tolist()))
        return [memo[p] for p in points]

    def local_atoms(self, center, clearance, reach):
        """Centres and radii of the atoms that can limit clearance within ``reach`` of ``center``.

        These are the atoms whose centres lie within ``clearance + 2*reach +
        r_max`` of ``center`` (``clearance`` is the clearance at ``center``):
        clearance is 1-Lipschitz, so every atom nearest to the surface of a
        point within ``reach`` is among them, and :func:`_surface_values` over
        them gives that point's exact clearance. Used by section refinement,
        which makes many small queries around one point. None without SciPy.

        Returns
        -------
        tuple of numpy.ndarray or None
            ``(xyz, radii)`` of shapes ``(k, 3)`` and ``(k,)``, in atom order.
        """
        import numpy as np

        spatial = self.spatial
        if spatial._tree is None:
            return None
        index = spatial._tree.query_ball_point(center, clearance + 2 * reach + spatial._max_radius + 1e-6)
        index = np.sort(np.asarray(index, dtype=np.intp))
        return spatial._xyz[index], spatial._radius_array[index]

    def segment_clearance(self, start, end):
        """Memoised :meth:`crevice.spatial.SpatialIndex.segment_clearance`."""
        key = (start, end)
        value = self._segment.get(key)
        if value is None:
            if len(self._segment) > self._MEMO_LIMIT:
                self._segment.clear()
            value = self._segment[key] = self.spatial.segment_clearance(start, end)
        return value


#: Frame cache shared by the probe geometries of one profile (see :mod:`crevice.adaptive_channel`).
_ACTIVE_FRAMES = contextvars.ContextVar("crevice_section_frames", default=None)


class _FrameCache:
    """Shares :class:`_FrameState` between geometries on the same atoms and frame."""

    def __init__(self):
        self._frames = {}

    def frame(self, atoms, origin, direction, basis):
        """The shared state for this frame, built on first use."""
        key = (id(atoms), len(atoms), tuple(origin), tuple(direction), basis)
        entry = self._frames.get(key)
        if entry is None or entry[0] is not atoms:
            entry = self._frames[key] = (atoms, _FrameState(atoms, origin, direction, basis))
        return entry[1]


class SectionGeometry:
    """Transverse grids perpendicular to a channel axis, and their enclosed sections.

    Each plane at axial coordinate ``t`` carries a square grid of
    ``2*ceil(extent/spacing)+1`` nodes per side spanned by the in-plane basis
    ``(u, v)``. Atoms are van der Waals spheres inflated by ``probe``; a node is
    free outside every inflated disk of that plane.

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms.
    origin : Coord
        Axis origin, Å.
    direction : Coord
        Unit axis direction.
    spacing : float, default 0.5
        In-plane grid spacing, Å.
    extent : float, default 6.0
        Half-width of the transverse search square, Å (``search_radius``).
    probe : float, default 0.0
        Enclosure-probe radius, Å, added to every atom radius.
    basis : tuple of Coord, optional
        In-plane axes ``(u, v)``; default from
        :func:`crevice.geometry.molecular_basis`.
    shared : object, optional
        Internal ``_FrameCache`` of probe-independent frame data (spatial index,
        projected atoms, node clearances) shared by several geometries on the
        same atoms and frame, as in the automatic enclosure-probe ladder. It
        changes no result.

    Raises
    ------
    ValueError
        If a plane would need more than 250000 grid points.
    """

    def __init__(self, atoms, origin, direction, *, spacing=0.5, extent=6.0, probe=0.0, basis=None, shared=None):
        self.atoms = tuple(atoms)
        self.origin, self.direction = origin, direction
        self.spacing, self.probe = spacing, probe
        self.half = int(math.ceil(extent / spacing))
        if (2 * self.half + 1) ** 2 > 250_000:
            raise ValueError("Section search exceeds 250000 points; reduce search_radius or increase section_spacing")
        # Everything below depends on the atoms and the frame only, never on
        # the probe: a ``_FrameCache`` shares it between the probes of the
        # automatic ladder.
        state = (
            shared.frame(self.atoms, origin, direction, basis)
            if shared is not None
            else _FrameState(self.atoms, origin, direction, basis)
        )
        self._state = state
        self.u, self.v = state.u, state.v
        self.spatial = state.spatial
        self.projected = state.projected
        self._sections = {}
        self._refined = {}

    def position(self, i, j, t):
        """Cartesian position of in-plane node ``(i, j)`` on the plane at ``t``.

        Parameters
        ----------
        i, j : int or float
            In-plane grid indices along ``u`` and ``v`` (0 = on the axis).
        t : float
            Axial coordinate, Å.

        Returns
        -------
        Coord
            ``origin + t*direction + i*spacing*u + j*spacing*v``, Å.
        """
        return add(
            self.origin,
            add(scale(self.direction, t), add(scale(self.u, i * self.spacing), scale(self.v, j * self.spacing))),
        )

    def _plane_masks(self, t):
        """Free nodes and blocked grid edges of one plane (numpy, exact).

        Returns ``(free, blocked)``: ``free`` is an ``(n, n)`` boolean array
        (index ``[i+h, j+h]``) of nodes outside every probe-inflated atom disk;
        ``blocked[k]`` marks, for the k-th neighbour step in ``_STEPS``, edges
        from node ``(i, j)`` whose closed segment touches a disk.

        The arithmetic is the same, operation by operation, as the scalar
        tests it replaced (node: ``(i*ds-x)**2 + (j*ds-y)**2 > rr + 1e-12``;
        edge: closest point of the segment to the disk centre), so results are
        bit-identical. Each disk is only evaluated on the index box that can
        reach it (disk radius plus one diagonal step, with a margin); nodes
        and edges outside that box cannot touch the disk. This is a speed-up
        only: it changes no result.
        """
        import numpy as np

        h, ds = self.half, self.spacing
        n = 2 * h + 1
        grid = np.arange(-h, h + 1) * ds
        index = np.arange(-h, h + 1)
        inside = np.zeros((n, n), dtype=bool)
        blocked = np.zeros((len(_STEPS), n, n), dtype=bool)
        reach = ds * math.sqrt(2.0)
        for x, y, z, r in self.projected:
            if not abs(z - t) <= r + self.probe:
                continue
            # Exact intersection of the inflated sphere with this plane.
            rr = (r + self.probe) ** 2 - (z - t) ** 2
            bound = math.sqrt(max(rr, 0.0) + 1e-12) + 1e-6
            for extent, target in ((bound, "node"), (bound + reach, "edge")):
                i0 = max(-h, int(math.floor((x - extent) / ds)) - 1)
                i1 = min(h, int(math.ceil((x + extent) / ds)) + 1)
                j0 = max(-h, int(math.floor((y - extent) / ds)) - 1)
                j1 = min(h, int(math.ceil((y + extent) / ds)) + 1)
                if i0 > i1 or j0 > j1:
                    continue
                si, sj = slice(i0 + h, i1 + h + 1), slice(j0 + h, j1 + h + 1)
                if target == "node":
                    inside[si, sj] |= ((grid[si] - x) ** 2)[:, None] + ((grid[sj] - y) ** 2)[None, :] <= rr + 1e-12
                    continue
                ii = index[si].astype(float)[:, None]
                jj = index[sj].astype(float)[None, :]
                for k, (di, dj) in enumerate(_STEPS):
                    fraction = np.maximum(
                        0.0, np.minimum(1.0, ((x / ds - ii) * di + (y / ds - jj) * dj) / (di * di + dj * dj))
                    )
                    blocked[k, si, sj] |= ((ii + fraction * di) * ds - x) ** 2 + (
                        (jj + fraction * dj) * ds - y
                    ) ** 2 <= rr + 1e-12
        return ~inside, blocked

    def _plane_graph(self, t):
        """Free nodes and blocked free-to-free edges of one plane (batched, exact).

        Same arithmetic as :meth:`_plane_masks`, element by element, on the
        same per-disk index boxes, but evaluated for all disks of the plane in
        a few NumPy operations instead of a Python loop over atoms. ``free``
        equals ``_plane_masks(t)[0]``; ``blocked[k, i, j]`` equals
        ``_plane_masks(t)[1][k, i, j]`` wherever node ``(i, j)`` and its
        k-th neighbour are both free (the only edges the flood fill in
        :meth:`sections` reads) and is False elsewhere. The disk radius
        ``rr`` keeps Python scalar arithmetic, as before.
        """
        import numpy as np

        h, ds, probe = self.half, self.spacing, self.probe
        n = 2 * h + 1
        state = self._state
        inside = np.zeros(n * n, dtype=bool)
        blocked = np.zeros((len(_STEPS), n * n), dtype=bool)
        near = np.nonzero(np.abs(state.pz - t) <= state.pr + probe)[0]
        if len(near) == 0:
            return ~inside.reshape(n, n), blocked.reshape(len(_STEPS), n, n)
        xs, ys = state.px[near], state.py[near]
        rr = np.array(
            [(r + probe) ** 2 - (z - t) ** 2 for r, z in zip(state.pr[near].tolist(), state.pz[near].tolist())]
        )
        threshold = rr + 1e-12
        bound = np.sqrt(np.maximum(rr, 0.0) + 1e-12) + 1e-6
        grid = np.arange(-h, h + 1) * ds

        def pairs(extent):
            i0 = np.maximum(-h, np.floor((xs - extent) / ds).astype(np.int64) - 1)
            i1 = np.minimum(h, np.ceil((xs + extent) / ds).astype(np.int64) + 1)
            j0 = np.maximum(-h, np.floor((ys - extent) / ds).astype(np.int64) - 1)
            j1 = np.minimum(h, np.ceil((ys + extent) / ds).astype(np.int64) + 1)
            ni, nj = np.maximum(i1 - i0 + 1, 0), np.maximum(j1 - j0 + 1, 0)
            sizes = ni * nj
            owner = np.repeat(np.arange(len(sizes)), sizes)
            local = np.arange(int(sizes.sum())) - np.repeat(np.cumsum(sizes) - sizes, sizes)
            width = nj[owner]
            return owner, i0[owner] + local // width, j0[owner] + local % width

        owner, ii, jj = pairs(bound)
        hit = ((grid[ii + h] - xs[owner]) ** 2 + (grid[jj + h] - ys[owner]) ** 2) <= threshold[owner]
        inside[((ii + h) * n + (jj + h))[hit]] = True
        free = ~inside
        owner, ii, jj = pairs(bound + ds * math.sqrt(2.0))
        start = (ii + h) * n + (jj + h)
        live = free[start]
        owner, ii, jj, start = owner[live], ii[live], jj[live], start[live]
        # All eight steps at once: one row per (pair, step).
        steps = np.array(_STEPS, dtype=np.int64)
        k = np.tile(np.arange(len(_STEPS)), len(ii))
        owner, ii, jj, start = (np.repeat(a, len(_STEPS)) for a in (owner, ii, jj, start))
        di, dj = steps[k, 0], steps[k, 1]
        ti, tj = ii + di, jj + dj
        ok = (np.abs(ti) <= h) & (np.abs(tj) <= h)
        ok[ok] = free[(ti[ok] + h) * n + (tj[ok] + h)]
        owner, ii, jj, start, k = owner[ok], ii[ok], jj[ok], start[ok], k[ok]
        di, dj = di[ok].astype(float), dj[ok].astype(float)
        x, y = xs[owner], ys[owner]
        fi, fj = ii.astype(float), jj.astype(float)
        fraction = np.maximum(0.0, np.minimum(1.0, ((x / ds - fi) * di + (y / ds - fj) * dj) / (di * di + dj * dj)))
        touch = ((fi + fraction * di) * ds - x) ** 2 + ((fj + fraction * dj) * ds - y) ** 2 <= threshold[owner]
        blocked[k[touch], start[touch]] = True
        return free.reshape(n, n), blocked.reshape(len(_STEPS), n, n)

    def sections(self, t):
        """Enclosed sections (connected free components) of the plane at ``t``.

        A node is free when the probe-inflated atom disks of this plane leave
        it uncovered; 8-connected free nodes are joined when the closed edge
        between them also stays outside every disk. Components that reach the
        search boundary are exterior and dropped; each remaining component is
        a :class:`Section` centred on its node of largest 3D clearance.
        """
        if t in self._sections:
            return list(self._sections[t])
        h = self.half
        free_mask, blocked = self._plane_graph(t)
        result = []
        components = _symmetric_components(free_mask, blocked, h)
        if components is not None:
            for component in components:
                values = self._state.clearances([self.position(i, j, t) for i, j in component])
                scored = [(value, i, j) for value, (i, j) in zip(values, component)]
                clearance, i, j = max(
                    scored, key=lambda row: (round(row[0], 10), -(row[1] ** 2) - row[2] ** 2, -row[1], -row[2])
                )
                result.append(Section(t, component, self.position(i, j, t), clearance))
            self._sections[t] = tuple(result)
            return result
        free = {(int(i) - h, int(j) - h) for i, j in zip(*free_mask.nonzero())}
        while free:
            first = min(free)
            free.remove(first)
            stack, component, exterior = [first], [], False
            while stack:
                node = stack.pop()
                component.append(node)
                i, j = node
                exterior |= abs(i) == h or abs(j) == h
                for k, (di, dj) in enumerate(_STEPS):
                    other = (i + di, j + dj)
                    if other not in free:
                        continue
                    # Minimum distance of each disk center to the closed edge.
                    if not blocked[k, i + h, j + h]:
                        free.remove(other)
                        stack.append(other)
            if not exterior:
                # Real 3D clearance, not the distance to planar disk boundaries.
                values = self._state.clearances([self.position(i, j, t) for i, j in component])
                scored = [(value, i, j) for value, (i, j) in zip(values, component)]
                clearance, i, j = max(
                    scored, key=lambda row: (round(row[0], 10), -(row[1] ** 2) - row[2] ** 2, -row[1], -row[2])
                )
                result.append(Section(t, tuple(sorted(component)), self.position(i, j, t), clearance))
        self._sections[t] = tuple(result)
        return result

    def refine(self, section, steps):
        """Converged local clearance maximum inside this component.

        In-plane clearance is the minimum over atoms of ``distance - radius``
        (3D distance, measured in the section plane). Its local maxima are
        vertices: points equidistant from three atom surfaces whose in-plane
        projections enclose them. Starting at the section centre (the node of
        largest clearance), two deterministic searches alternate (at most
        ``REFINEMENT_MAX_ROUNDS`` rounds):

        * a local ascent (``_pattern_search``): at each step length the
          better of the eight in-plane pattern moves and the move that
          maximises the linearised clearance (``_linear_move``, which
          follows the narrow two-atom ridges on which fixed-direction moves
          stall), repeated until neither improves clearance, then the step is
          halved, from half the grid spacing until it falls below
          ``REFINEMENT_TOLERANCE_A`` (and at least ``steps`` halvings);
        * a vertex step (``_vertex_ascent``): the highest vertex within
          ``REFINEMENT_VERTEX_RADIUS_A`` that improves clearance is solved for
          exactly, and the search ends at a vertex with none higher within
          reach. Such a vertex is an exact strict local maximum.

        Every move must stay in the component's grid cells and keep the
        straight in-plane segment clear of the enclosure probe, so the centre
        stays in the same free space. The result is the maximum of the hill
        that holds the grid centre, not the highest maximum of the section: a
        section can hold a higher hill that the grid samples less well, but a
        centre moved there could break the straight joins to neighbouring
        sections. Results are memoised per plane, start node and ``steps``.

        Parameters
        ----------
        section : Section
            Section to refine.
        steps : int
            ``refinement_steps``: 0 disables refinement (the grid-node centre
            is kept); a positive value is the minimum number of step halvings.

        Returns
        -------
        Section
            Same plane and nodes, with the refined centre and its clearance.
        """
        if steps <= 0:
            return section
        key = (section.t, section.center, steps)
        cached = self._refined.get(key)
        if cached is None:
            cached = self._refined[key] = self._climb(
                section.center, section.clearance, steps, self._component_cells(section)
            )
        return replace(section, center=cached[0], clearance=cached[1])

    def _climb(self, best, clearance, steps, inside):
        """Local clearance maximum from one start: vertex steps, pattern search when none applies."""
        step, at_vertex = self.spacing / 2, False
        for _ in range(REFINEMENT_MAX_ROUNDS):
            vertex = self._vertex_ascent(best, clearance, inside)
            if vertex is not None:
                (best, clearance), at_vertex = vertex, True
                continue
            if at_vertex:
                # A valid vertex is an exact local maximum; none higher within reach.
                break
            # No improving vertex within reach of a non-vertex point: pattern
            # search (the whole method without SciPy), then vertices again.
            moved = self._pattern_search(
                best, clearance, step, min(self.spacing / 2**steps, REFINEMENT_TOLERANCE_A), inside
            )
            if moved[1] <= clearance:
                break
            best, clearance, step = moved
        return best, clearance

    def _component_cells(self, section):
        """Test of whether points lie in a grid cell with a corner node in ``section``.

        Returns ``f(points) -> bool array``. Refinement stays in these cells:
        the enclosure test is a grid test, so a free continuous route through
        a gap narrower than the grid can join a component to outside space,
        and an unbounded search could follow it out of the section.
        """
        import numpy as np

        index = np.array(section.indices, dtype=np.intp).reshape(-1, 2)
        h = int(np.abs(index).max(initial=0))  # bridged sections extend past ``half``
        mask = np.zeros((2 * h + 3, 2 * h + 3), dtype=bool)
        mask[index[:, 0] + h + 1, index[:, 1] + h + 1] = True
        base = np.array(add(self.origin, scale(self.direction, section.t)))
        u, v = np.array(self.u), np.array(self.v)

        def inside(points):
            rel = np.asarray(points, dtype=float) - base
            a, b = np.floor(rel @ u / self.spacing), np.floor(rel @ v / self.spacing)
            ok = (a >= -h - 1) & (a <= h) & (b >= -h - 1) & (b <= h)
            a, b = (np.where(ok, w, 0).astype(np.intp) + h + 1 for w in (a, b))
            return ok & (mask[a, b] | mask[a + 1, b] | mask[a, b + 1] | mask[a + 1, b + 1])

        return inside

    def _clear_move(self, start, start_value, end, end_value):
        """Whether the straight move ``start -> end`` clears the enclosure probe."""
        # Clearance is 1-Lipschitz, so on a segment of length L it is at least
        # (start + end - L)/2; only when that bound fails is the exact check run.
        return (
            start_value + end_value - distance(start, end)
        ) / 2 > self.probe + 1e-9 or self._state.segment_clearance(start, end) > self.probe

    def _evaluator(self, center, clearance, reach):
        """Exact clearances near ``center``: ``(f(points) -> list of float, local atoms or None)``.

        ``local`` is :meth:`_FrameState.local_atoms` (``(xyz, radii)``, every
        atom that can limit a point within ``reach``), None without SciPy.
        """
        local = self._state.local_atoms(center, clearance, reach)
        if local is None:
            return (lambda points: self._state.clearances([tuple(p) for p in points.tolist()])), None
        xyz, radii = local
        return (lambda points: _surface_values(points, xyz, radii).tolist()), local

    def _linear_move(self, best, clearance, step, local):
        """Best move of the linearised clearance within ``step`` of ``best`` (see :meth:`refine`), or None.

        Each atom's clearance ``f_i = |p - a_i| - r_i`` is convex, so its
        linearisation at ``best`` (value ``f_i``, in-plane gradient ``g_i``) is
        a lower bound everywhere. The move ``d`` (``|d_u|, |d_v| <= step``)
        maximising ``min_i f_i + g_i.d`` over the atoms that can be nearest in
        that square (``f_i <= clearance + 2*sqrt(2)*step``; the others stay
        above the bound) therefore has true clearance at least that maximum.
        The maximum of this piecewise-linear concave function is at a corner of
        the square, where two planes meet on an edge, or where three planes
        meet; all are enumerated. Returns the point (Cartesian) when the bound
        exceeds ``clearance``.
        """
        import numpy as np

        xyz, radii = local
        p = np.array(best)
        rel = p - xyz
        dist = np.sqrt((rel * rel).sum(axis=1))
        f = dist - radii
        near = np.flatnonzero(f <= clearance + 2 * _SQRT2 * step + 1e-9)
        if len(near) == 0:
            return None
        f, gx, gy = f[near], (rel[near] @ np.array(self.u)) / dist[near], (rel[near] @ np.array(self.v)) / dist[near]
        moves = [np.array([[-step, -step], [-step, step], [step, -step], [step, step]])]
        if len(near) >= 2:
            i, j = np.triu_indices(len(near), 1)
            a, b, c = gx[i] - gx[j], gy[i] - gy[j], f[j] - f[i]  # a*du + b*dv = c
            with np.errstate(invalid="ignore", divide="ignore"):
                for edge in (-step, step):
                    moves.append(np.stack([np.full_like(a, edge), (c - a * edge) / b], axis=1))
                    moves.append(np.stack([(c - b * edge) / a, np.full_like(a, edge)], axis=1))
        if len(near) >= 3:
            i, j, k = _triples(len(near))
            a11, a12, c1 = gx[i] - gx[j], gy[i] - gy[j], f[j] - f[i]
            a21, a22, c2 = gx[i] - gx[k], gy[i] - gy[k], f[k] - f[i]
            det = a11 * a22 - a12 * a21
            with np.errstate(invalid="ignore", divide="ignore"):
                moves.append(np.stack([(c1 * a22 - a12 * c2) / det, (a11 * c2 - c1 * a21) / det], axis=1))
        moves = np.concatenate(moves)
        moves = moves[np.isfinite(moves).all(axis=1) & (np.abs(moves) <= step * (1 + 1e-9)).all(axis=1)]
        moves = np.clip(moves, -step, step)
        bound = (f[None, :] + moves[:, 0:1] * gx[None, :] + moves[:, 1:2] * gy[None, :]).min(axis=1)
        m = int(bound.argmax())
        if not bound[m] > clearance + 1e-12:
            return None
        return p + moves[m, 0] * np.array(self.u) + moves[m, 1] * np.array(self.v)

    def _pattern_search(self, best, clearance, step, stop, inside):
        """Converging local ascent: pattern and linearised moves (see :meth:`refine`).

        At each step length the candidates are the eight in-plane pattern
        moves and the linearised move (:meth:`_linear_move`); the highest
        exact clearance that satisfies ``inside`` (:meth:`_component_cells`)
        and clears the probe is taken, repeated until none improves, then the
        step is halved. Runs while the step is at least ``stop``. Returns
        ``(centre, clearance, next step)``.
        """
        import numpy as np

        moves = 0
        u, v = np.array(self.u), np.array(self.v)
        di, dj = (np.array([d[k] for d in _STEPS], dtype=float)[:, None] for k in (0, 1))
        # Exact clearances from the atoms near the current anchor (rebuilt when
        # the search could leave a neighbourhood of four grid spacings).
        reach = 4 * self.spacing
        anchor, (evaluate, local) = best, self._evaluator(best, clearance, reach)
        while step >= stop and moves < REFINEMENT_MAX_MOVES:
            if distance(anchor, best) + step * _SQRT2 > reach:
                anchor, (evaluate, local) = best, self._evaluator(best, clearance, reach)
            # The operations of geometry.add/scale, element-wise: best + (u*(di*step) + v*(dj*step)).
            candidates = np.array(best) + (u * (di * step) + v * (dj * step))
            linear = None if local is None else self._linear_move(best, clearance, step, local)
            if linear is not None:
                candidates = np.vstack([candidates, linear])
            chosen, score = None, clearance
            allowed = inside(candidates)
            for k, value in enumerate(evaluate(candidates)):
                if (
                    allowed[k]
                    and value > score + 1e-12
                    and self._clear_move(best, clearance, tuple(candidates[k].tolist()), value)
                ):
                    chosen, score = k, value
            if chosen is not None:
                best, clearance, moves = tuple(candidates[chosen].tolist()), score, moves + 1
            else:
                step /= 2
        return best, clearance, step

    def _vertex_ascent(self, best, clearance, inside):
        """Highest three-atom clearance vertex within reach that improves on ``best``, or None.

        For atoms ``i`` with in-plane projections ``(x_i, y_i)`` (relative to
        ``best``), squared out-of-plane offsets ``h_i**2`` and radii ``r_i``,
        a point ``(X, Y)`` with clearance ``c`` to three atoms satisfies
        ``(X-x_i)**2 + (Y-y_i)**2 + h_i**2 = (c + r_i)**2``. Differences of
        the three equations are linear in ``(X, Y, c)``; substituting gives a
        quadratic in ``c``. A root is a local maximum when ``(X, Y)`` lies in
        the triangle of the three projections. Kept: roots within
        ``REFINEMENT_VERTEX_RADIUS_A`` of ``best`` whose exact clearance (all
        atoms) equals the root (no other atom is closer) and exceeds
        ``clearance``, inside the section's cells (``inside``) and reachable
        by a probe-clear straight move; the highest
        wins (ties: first triple in atom order).

        None without SciPy.
        """
        import numpy as np

        radius = REFINEMENT_VERTEX_RADIUS_A
        local = self._state.local_atoms(best, clearance, radius)
        if local is None:
            return None
        xyz, radii = local
        p = np.array(best)
        rel = xyz - p
        u, v, n = np.array(self.u), np.array(self.v), np.array(self.direction)
        excess = np.sqrt((rel * rel).sum(axis=1)) - radii - clearance
        x, y, h2 = rel @ u, rel @ v, (rel @ n) ** 2
        # Any atom limiting a point q within ``radius``: f_i(best) <= f(q) + |q - best| <= clearance + 2*radius.
        found = self._vertices(np.flatnonzero(excess <= 2 * radius + 1e-9), x, y, h2, radii, clearance, radius)
        if not found:
            return None
        X, Y, c = (np.array(w) for w in zip(*found))
        points = p + X[:, None] * u + Y[:, None] * v
        values = _surface_values(points, xyz, radii)
        good = (values >= c - 1e-7) & (values > clearance + 1e-12) & inside(points)
        for m in sorted(np.flatnonzero(good).tolist(), key=lambda m: -values[m]):
            end = tuple(points[m].tolist())
            if self._clear_move(best, clearance, end, float(values[m])):
                return end, float(values[m])
        return None

    @staticmethod
    def _vertices(atoms, x, y, h2, r, clearance, radius):
        """Three-atom clearance maxima among ``atoms`` (see :meth:`_vertex_ascent`).

        Returns ``(X, Y, c)`` of the roots inside their triangle, above
        ``clearance`` and within ``radius`` of the origin, in triple order.
        """
        import numpy as np

        if len(atoms) < 3:
            return []
        i, j, k = (atoms[c] for c in _triples(len(atoms)))
        a11, a12, a21, a22 = 2 * (x[j] - x[i]), 2 * (y[j] - y[i]), 2 * (x[k] - x[i]), 2 * (y[k] - y[i])
        det = a11 * a22 - a12 * a21
        keep = np.abs(det) > 1e-9
        i, j, k, a11, a12, a21, a22, det = (w[keep] for w in (i, j, k, a11, a12, a21, a22, det))
        area = det / 4

        def inside(X, Y):
            with np.errstate(invalid="ignore", divide="ignore"):
                lj = ((X - x[i]) * (y[k] - y[i]) - (Y - y[i]) * (x[k] - x[i])) / area
                lk = ((x[j] - x[i]) * (Y - y[i]) - (y[j] - y[i]) * (X - x[i])) / area
            return (lj >= -1e-9) & (lk >= -1e-9) & (1 - lj - lk >= -1e-9)

        q = x * x + y * y + h2 - r * r
        c1, c2, d1, d2 = q[j] - q[i], q[k] - q[i], 2 * (r[j] - r[i]), 2 * (r[k] - r[i])
        # (X, Y) = (ax, ay) + c*(bx, by)
        ax, ay = (c1 * a22 - a12 * c2) / det, (a11 * c2 - c1 * a21) / det
        bx, by = -(d1 * a22 - a12 * d2) / det, -(a11 * d2 - d1 * a21) / det
        wx, wy = ax - x[i], ay - y[i]
        qa, qb, qc = bx * bx + by * by - 1, 2 * (wx * bx + wy * by - r[i]), wx * wx + wy * wy + h2[i] - r[i] ** 2
        found = []
        rmin = np.minimum(np.minimum(r[i], r[j]), r[k])
        with np.errstate(invalid="ignore", divide="ignore"):
            root = np.sqrt(qb * qb - 4 * qa * qc)
            linear = np.abs(qa) < 1e-12
            for c in (
                np.where(linear, -qc / qb, (-qb + root) / (2 * qa)),
                np.where(linear, np.nan, (-qb - root) / (2 * qa)),
            ):
                X, Y = ax + c * bx, ay + c * by
                ok = (
                    np.isfinite(c)
                    & (c > clearance + 1e-12)
                    & (c + rmin > 0)
                    & (X * X + Y * Y <= radius * radius)
                    & inside(X, Y)
                )
                found.extend(zip(X[ok].tolist(), Y[ok].tolist(), c[ok].tolist()))
        return found


def _symmetric_components(free_mask, blocked, h):
    """Interior components of a plane graph, or None when an edge test is asymmetric.

    The flood fill in :meth:`SectionGeometry.sections` follows an edge from
    the node it visits, testing that node's own segment. When every free edge
    gives the same answer from both ends (the usual case; the two tests
    differ only by rounding), the graph is undirected and the flood fill
    finds exactly its connected components. This labels them with
    :func:`scipy.sparse.csgraph.connected_components` and returns, in the
    flood fill's discovery order (smallest node first), the sorted node
    tuples of the components that never reach the search boundary. Returns
    None (use the flood fill) without SciPy or when any free edge is
    asymmetric, so the result is always the flood fill's.
    """
    try:
        import numpy as np
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
    except ImportError:
        return None
    n = 2 * h + 1
    free = free_mask.ravel()
    flat = blocked.reshape(len(_STEPS), -1)
    rows, cols = [], []
    i_all, j_all = np.divmod(np.arange(n * n), n)
    for k, (di, dj) in enumerate(_STEPS):
        if (di, dj) <= (0, 0):
            continue  # each undirected edge once, from its smaller end
        back = _STEPS.index((-di, -dj))
        ok = free & (i_all + di < n) & (j_all + dj >= 0) & (j_all + dj < n)
        a = np.nonzero(ok)[0]
        b = a + di * n + dj
        keep = free[b]
        a, b = a[keep], b[keep]
        forward, backward = ~flat[k, a], ~flat[back, b]
        if (forward != backward).any():
            return None
        rows.append(a[forward])
        cols.append(b[forward])
    rows = np.concatenate(rows) if rows else np.zeros(0, dtype=np.intp)
    cols = np.concatenate(cols) if cols else np.zeros(0, dtype=np.intp)
    graph = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n * n, n * n)).tocsr()
    _, labels = connected_components(graph, directed=False)
    nodes = np.nonzero(free)[0]
    if len(nodes) == 0:
        return []
    node_labels = labels[nodes]
    boundary = (i_all[nodes] == 0) | (i_all[nodes] == n - 1) | (j_all[nodes] == 0) | (j_all[nodes] == n - 1)
    exterior = np.zeros(labels.max() + 1, dtype=bool)
    exterior[node_labels[boundary]] = True
    order = np.argsort(node_labels, kind="stable")
    grouped, grouped_labels = nodes[order], node_labels[order]
    cuts = np.nonzero(np.diff(grouped_labels))[0] + 1
    starts = np.concatenate(([0], cuts))
    ends = np.concatenate((cuts, [len(grouped)]))
    components = []
    for start, end in zip(starts.tolist(), ends.tolist()):
        if exterior[grouped_labels[start]]:
            continue
        members = grouped[start:end]
        components.append((int(members[0]), tuple(zip((i_all[members] - h).tolist(), (j_all[members] - h).tolist()))))
    components.sort()
    return [members for _, members in components]


#: Bridging window half-width as a multiple of the plane's search square (see :mod:`crevice.adaptive_channel`).
BRIDGE_WINDOW_FACTOR = 2.0
#: Fixed enclosure probe (Å) of the fallback cast written when no channel resolves (``void_cast``).
LEGACY_ENCLOSURE_RADIUS = 0.8


def connected_profile(
    atoms,
    *,
    axis,
    origin,
    samples,
    search_radius,
    refinement_steps,
    probe_radius,
    section_spacing=0.5,
    enclosure_radius=None,
    lateral_exits=False,
    exit_bulk_radius=6.0,
    exit_spacing=0.5,
):
    """Resolve one connected, laterally enclosed channel and sample it.

    Entry point of the channel tracer of :mod:`crevice.adaptive_channel`.
    ``enclosure_radius=None`` (or ``"auto"``) gives each plane its own
    enclosure probe, the smallest of
    :data:`crevice.adaptive_channel.ADAPTIVE_ENCLOSURE_LADDER` at which its
    lumen is enclosed; a number runs that one probe in every plane (a plane
    whose lumen is narrower than it ends the channel there).

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms (already selected).
    axis : str or Coord
        ``"auto"`` (try the principal directions and keep a unique channel),
        ``"x"``, ``"y"``, ``"z"`` or a direction vector.
    origin : Coord or None
        Seed inside the intended channel, Å. ``None`` uses the atom centroid.
    samples : int
        Minimum number of output samples (axial steps stay at most 0.75 Å).
    search_radius : float
        Base half-width of the transverse search square, Å.
    refinement_steps : int
        Minimum step halvings of each in-plane clearance climb.
    probe_radius : float
        Measurement probe radius, Å, subtracted from each reported radius.
    section_spacing : float, default 0.5
        Transverse grid spacing, Å.
    enclosure_radius : float, "auto" or None, default None
        Enclosure probe, Å, or None/``"auto"`` for the per-plane probe.
    lateral_exits : bool, default False
        Accept a capped end that has a lateral exit to bulk.
    exit_bulk_radius : float, default 6.0
        With ``lateral_exits``: rolling bulk-probe radius, Å.
    exit_spacing : float, default 0.5
        With ``lateral_exits``: lattice spacing, Å, of the exit search.

    Returns
    -------
    PoreProfile
        See :func:`crevice.adaptive_channel.adaptive_profile`.

    Raises
    ------
    ChannelResolutionError
        When no unique connected channel is resolved.
    OriginClearanceError
        For a supplied ``origin`` without clearance for the smallest probe.
    """
    from .adaptive_channel import ADAPTIVE_ENCLOSURE_LADDER, adaptive_profile

    fixed = not (
        enclosure_radius is None
        or (isinstance(enclosure_radius, str) and enclosure_radius.lower() in {"auto", "adaptive"})
    )
    return adaptive_profile(
        atoms,
        axis=axis,
        origin=origin,
        samples=samples,
        search_radius=search_radius,
        refinement_steps=refinement_steps,
        probe_radius=probe_radius,
        section_spacing=section_spacing,
        lateral_exits=lateral_exits,
        exit_bulk_radius=exit_bulk_radius,
        exit_spacing=exit_spacing,
        ladder=(float(enclosure_radius),) if fixed else ADAPTIVE_ENCLOSURE_LADDER,
        fixed=fixed,
    )


def _refinement_metadata(refinement_steps):
    """Description of the in-plane centre climb (:mod:`crevice.adaptive_channel`) for profile metadata."""
    return {
        "method": "converging local ascent within the plane from the previous plane's centre: best of eight "
        "in-plane pattern moves and the linearised-clearance move, step halved when none improves; then "
        "the exact three-atom clearance vertex (local maximum) within reach; repeated from the enclosed "
        "component's highest grid node, keeping the higher top when both lie in one lumen",
        "initial_step": "half the section grid spacing",
        "position_tolerance_A": REFINEMENT_TOLERANCE_A,
        "vertex_search_radius_A": REFINEMENT_VERTEX_RADIUS_A,
        "minimum_halvings": max(1, refinement_steps),
        "max_moves": REFINEMENT_MAX_MOVES,
        "limitations": "local maxima; a higher maximum in another passage joined only through a narrower neck "
        "is not taken",
    }


def channel_cast(
    frame,
    profile,
    *,
    spacing,
    min_radius,
    max_grid_points,
    max_components,
    max_export_points,
    min_component_points,
    focus_points,
    focus_radius,
    min_component_volume,
    max_component_volume,
    include_hydrogen,
    include_hetero,
):
    """Probe-swept lumen from all connected enclosed probe-center sections.

    The full non-circular section contributes, rather than just the inscribed
    sphere at its center. Sweeping the enclosure probe recovers near-wall volume;
    nodes still satisfy the requested atom clearance. Axial ends are explicit
    mouth planes, so the sweep never extends into bulk solvent beyond them.

    Algorithm: on a 3D grid aligned with the profile's axis and section basis,
    each plane between the profile mouths contributes its enclosed section
    nearest the interpolated profile centre (if the straight segment to it
    clears the probe); these core nodes are swept by a sphere of radius
    ``max(enclosure probe, min_radius) - min_radius``; swept nodes with
    atom-surface clearance >= ``min_radius`` are clustered with continuous edge
    tests (:class:`crevice.grid.GridConnectivity`), and the final components are
    selected by size. Planes with no eligible section are recorded as
    unresolved (and warned about), never filled.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile.
    profile : PoreProfile
        Resolved ``"axial-connected"`` profile from :func:`connected_profile`;
        its atom selection must match ``include_hydrogen``/``include_hetero``.
    spacing : float
        Cast grid spacing, Å (never coarsened automatically).
    min_radius : float
        Minimum atom-surface clearance of a cast node, Å.
    max_grid_points : int
        Largest allowed grid; larger requests raise ``ValueError``.
    max_components : int or None
        Maximum number of components kept.
    max_export_points : int
        Largest number of display samples; the measurement keeps every node.
    min_component_points : int
        Components with fewer nodes are dropped.
    focus_points : sequence of Coord or None
        Keep only nodes within ``focus_radius`` of one of these points.
    focus_radius : float
        Focus crop radius, Å.
    min_component_volume, max_component_volume : float or None
        Inclusive volume limits (Å³) of the final components.
    include_hydrogen, include_hetero : bool
        Atom selection; must equal the profile's.

    Returns
    -------
    VoidCast
        Mode ``"channel"``; metadata include the grid origin, basis and
        dimensions, enclosure and sweep radii, ``export_mode
        = "connected_section_fill"``, ``mouth_planes_A``,
        ``unresolved_section_coordinates_A`` and ``status`` (``resolved`` or
        ``partial_or_empty``).

    Raises
    ------
    ValueError
        For an atom-selection mismatch or a grid above ``max_grid_points``.
    """
    from bisect import bisect_right
    from itertools import product
    from .grid import GridConnectivity
    from .models import VoidCast, VoidComponent
    from .voids import _cluster_point_indices, _select_final_components, _renumber_points, _sample_points_by_axis

    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    for name, value in (("include_hydrogen", include_hydrogen), ("include_hetero", include_hetero)):
        if profile.metadata.get(name) != value:
            raise ValueError("Channel profile and cast must use the same atom selection")
    q = max(profile.metadata["enclosure_radius_A"], min_radius)
    # Per-plane enclosure probes (automatic mode): each cast plane uses the larger probe of the two
    # profile samples around it, so a plane is never filled with a probe its lumen leaks at.
    plane_probes = profile.metadata.get("enclosure_probes_A")
    extent = max(profile.metadata.get("section_search_radii_A") or [profile.metadata["section_search_radius_A"]])
    geometry = SectionGeometry(
        atoms,
        profile.axis_origin,
        profile.axis_direction,
        spacing=spacing,
        extent=extent,
        probe=q,
        basis=profile.metadata["section_basis"],
    )
    per_probe = {q: geometry}

    def geometry_at(probe):
        if probe not in per_probe:
            state = geometry._state

            class _Same:
                def frame(self, *args):
                    return state

            per_probe[probe] = SectionGeometry(
                atoms,
                profile.axis_origin,
                profile.axis_direction,
                spacing=spacing,
                extent=extent,
                probe=probe,
                basis=(geometry.u, geometry.v),
                shared=_Same(),
            )
        return per_probe[probe]

    low, high = profile.points[0].t, profile.points[-1].t
    nz = int(math.floor((high - low) / spacing + 1e-9)) + 1
    h = geometry.half
    dims = (2 * h + 1, 2 * h + 1, nz)
    count = math.prod(dims)
    if count > max_grid_points:
        raise ValueError(
            f"Channel grid needs {count} points at {spacing:g} A; increase max_grid_points "
            "or explicitly choose a coarser spacing. No automatic coarsening was applied."
        )
    origin = geometry.position(-h, -h, low)
    basis = (geometry.u, geometry.v, geometry.direction)
    connectivity = GridConnectivity(geometry.spatial, origin, spacing, min_radius, basis=basis)
    ts = [p.t for p in profile.points]
    core, gaps, plane_q = {}, [], []
    for k in range(nz):
        t = low + k * spacing
        j = max(1, min(len(ts) - 1, bisect_right(ts, t)))
        a, b = profile.points[j - 1], profile.points[j]
        f = (t - a.t) / (b.t - a.t)
        target = add(a.position, scale(sub(b.position, a.position), f))
        qk = q if plane_probes is None else max(plane_probes[j - 1], plane_probes[j], min_radius)
        # Per-plane mode: a cast plane between two samples can leak at their probes; like the profile,
        # it then takes the next larger probe of the ladder (below the centre line's clearance).
        if plane_probes is None:
            candidates = [qk]
        else:
            from .adaptive_channel import ADAPTIVE_ENCLOSURE_LADDER

            reach = geometry._state.clearances([target])[0]
            candidates = [qk] + [p for p in ADAPTIVE_ENCLOSURE_LADDER if qk < p < reach]
        eligible = []
        for qk in candidates:
            plane = geometry_at(qk)
            eligible = [s for s in plane.sections(t) if plane._state.segment_clearance(target, s.center) >= qk]
            if eligible:
                break
        plane_q.append(qk)
        if not eligible:
            gaps.append(t)
            continue
        section = min(eligible, key=lambda s: distance(target, s.center))
        core.update(((i + h, j + h, k), qk) for i, j in section.indices)
    # Keep only full 3D components that attach to the measured centerline.
    # Separate planes cannot bypass the continuous segment rule.
    radius = (
        max(0.0, max(core.values(), default=q) - min_radius)
        if plane_probes is None
        else max(0.0, max((max(v, spacing) for v in core.values()), default=q) - min_radius)
    )
    offsets_for = {}
    swept = set()
    for idx, qk in core.items():
        # Per-plane mode: a plane probe can be smaller than the cast grid, and a sweep shorter than one
        # spacing reaches no node, dropping the shell between the wall and the probe centres; sweep at
        # least one spacing there. A fixed probe sweeps exactly its radius, as before.
        sweep = max(0.0, (qk if plane_probes is None else max(qk, spacing)) - min_radius)
        if sweep not in offsets_for:
            cells = int(math.ceil(sweep / spacing))
            offsets_for[sweep] = [
                d
                for d in product(range(-cells, cells + 1), repeat=3)
                if spacing**2 * sum(x * x for x in d) <= sweep**2 + 1e-12
            ]
        swept.update(tuple(i + d for i, d in zip(idx, offset)) for offset in offsets_for[sweep])
    swept = {idx for idx in swept if all(0 <= i < n for i, n in zip(idx, dims))}
    nodes = {}
    for idx in sorted(swept):
        position = connectivity.position(idx)
        raw, nearest = geometry.spatial.nearest_surface(position)
        if raw >= min_radius:
            nodes[idx] = ChannelPoint(
                0, position, raw, raw, low + idx[2] * spacing, nearest.serial, nearest.residue_key.label
            )
    if focus_points:
        nodes = {
            idx: p for idx, p in nodes.items() if any(distance(p.position, f) <= focus_radius for f in focus_points)
        }
    clusters = _cluster_point_indices(
        set(nodes), edge_allowed=lambda a, b: connectivity.allows(a, b, nodes[a].raw_clearance, nodes[b].raw_clearance)
    )
    candidates = []
    for idx, cluster in enumerate(clusters, 1):
        if len(cluster) < min_component_points:
            continue
        points = _renumber_points(sorted((nodes[i] for i in cluster), key=lambda p: (p.t, p.position)))
        center = tuple(sum(p.position[j] for p in points) / len(points) for j in range(3))
        candidates.append(
            VoidComponent(
                idx,
                "section_channel",
                center,
                len(points) * spacing**3,
                points,
                nearest_residues=tuple(sorted({p.nearest_residue for p in points})),
                metadata={
                    "source_component_id": idx,
                    "grid_origin": origin,
                    "grid_basis": basis,
                    "grid_spacing": spacing,
                    "source_kind": "probe_enclosed_channel",
                    "selection_cropped": bool(focus_points),
                },
            )
        )
    components, selection = _select_final_components(
        candidates, max_components=max_components, min_volume=min_component_volume, max_volume=max_component_volume
    )
    full = [p for c in components for p in c.points]
    displayed = _renumber_points(_sample_points_by_axis(full, max_export_points))
    if gaps:
        import warnings

        warnings.warn(
            f"Channel cast is unresolved in {len(gaps)} sections at {spacing:g} A; refine the cast grid. "
            "See unresolved_section_coordinates_A.",
            RuntimeWarning,
            stacklevel=2,
        )
    return VoidCast(
        components,
        displayed,
        spacing,
        min_radius,
        "channel",
        metadata={
            "axis": profile.metadata["axis"],
            "axis_origin": profile.axis_origin,
            "axis_direction": profile.axis_direction,
            "grid_origin": origin,
            "grid_basis": basis,
            "grid_dimensions": dims,
            "grid_point_count": count,
            "requested_spacing_A": spacing,
            "effective_spacing_A": spacing,
            "spacing_adjusted": False,
            "enclosure_radius_A": max(plane_q, default=q),
            "probe_sweep_radius_A": radius,
            "plane_enclosure_radii_A": plane_q if plane_probes is not None else None,
            "export_mode": "connected_section_fill",
            "measurement_point_count": len(full),
            "core_point_count": len(core),
            "display_subsampled": len(full) > len(displayed),
            "max_export_points": max_export_points,
            "max_components": max_components,
            "selection": selection,
            "focus_points": focus_points,
            "focus_radius": focus_radius if focus_points else None,
            "min_component_volume": min_component_volume,
            "max_component_volume": max_component_volume,
            "include_hydrogen": include_hydrogen,
            "include_hetero": include_hetero,
            "connectivity": connectivity.to_dict(),
            "unresolved_section_coordinates_A": gaps,
            "status": "resolved" if not gaps and len(components) == 1 else "partial_or_empty",
            "mouth_planes_A": (low, low + (nz - 1) * spacing),
            "channel_mouths": profile.metadata.get("channel_mouths"),
            "selection_note": "Probe-swept connected lumen, capped at mouth planes; finite-resolution atom-sphere geometry",
            "boundary_validation": "voxel boundaries not continuously validated",
            "limitations": "Narrower crevices than the enclosure probe are excluded; grid refinement is required to assess volume convergence",
        },
    )

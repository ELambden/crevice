"""Channel tracing with a per-plane enclosure probe (the automatic profile mode).

A fixed enclosure probe has to be smaller than the narrowest constriction (or
the channel breaks there) and larger than every gap in the channel wall (or
the lumen leaks sideways and is not enclosed). In many real pores no single
probe does both: in KcsA the 0.75 Å selectivity filter needs a probe below
0.75 Å, while the central cavity leaks through gaps between the helices at
every probe up to 0.7 Å. This module lets each plane use its own probe.

Algorithm
---------
Planes are normal to the channel axis. In each plane:

1. **Lumen point.** Starting from the previous plane's centre (the seed for
   the first plane), the in-plane clearance is climbed to its local maximum
   (the same exact refinement as :meth:`crevice.sections.SectionGeometry.refine`,
   without the section constraint). A plane's centre therefore never jumps to
   a separate pocket in the wall that happens to be enclosed.
2. **Per-plane enclosure probe.** (A plane whose clearance is too close to
   every probe for the section grid, less than half a grid diagonal above it,
   is re-gridded locally at half the spacing, up to ``MAX_GRID_REFINEMENT``
   times; whether a constriction resolves then does not depend on where the
   samples fall.) The plane's probe is the smallest probe of
   :data:`ADAPTIVE_ENCLOSURE_LADDER` (below the centre's clearance) at which
   the free component holding the centre is laterally enclosed: it does not
   reach the transverse search square or, failing that, a square
   :data:`crevice.sections.BRIDGE_WINDOW_FACTOR` times as wide. A component
   that is enclosed at one probe stays enclosed at every larger probe (the free
   space only shrinks), so the smallest enclosing probe is found by bisection
   over the ladder. This probe is the size of the widest lateral opening of
   the lumen at that plane, and it is reported for every sample
   (``metadata["enclosure_probes_A"]``).
3. **Widest hill.** A wide lumen can hold several clearance maxima. The
   climb is repeated from the enclosed component's highest grid node, and the
   higher of the two maxima is the plane's centre when it lies in the same
   component and the straight segment between them keeps at least
   ``HILL_RIDGE_FRACTION`` (90 %) of the lower maximum's clearance, so both
   lie in one lumen. A neighbouring passage reached only through a narrower
   neck is never taken, even when the enclosed component at a small probe
   contains it.
4. **Join.** A plane continues the channel when the straight segment from the
   previous centre clears the smaller of the two planes' probes.

The channel ends (its **mouth**) at the first plane, going outwards from the
seed, whose lumen is not enclosed at any probe below its clearance or cannot
be joined, or at the van der Waals envelope of the atoms along the axis
(planes beyond it cut no atom); the mouth is bracketed by bisection to 0.025 Å. An end is an
axial exit when the straight continuation along the axis beyond the atom
envelope clears that end's probe, open (``bulk``) when its sample already lies
in bulk solvent (the rolling bulk-probe envelope beyond the mouth); otherwise it
is capped, and is accepted only with ``lateral_exits`` (:mod:`crevice.channel_exits`,
searched at the end plane's probe).

With ``lateral_exits``, a capped end is moved to the plane where the widest
route to bulk leaves the lumen (``LATERAL_BRANCH_RULE``, ``_lateral_branches``):
when a lateral route from a plane beyond both the seed and the channel's
narrowest plane is clearly wider than the axial continuation beyond that plane
(its planes and the end's own exit), the lumen ends at the plane with the widest
such route (mouth reason ``lateral_branch``) and that end's lateral exits are
reported. Open ends never move. Side portals that
leave a chamber obliquely, which no single plane shows as an opening wider than
the lumen, are found this way.

The measured radius is the clearance at each centre minus ``probe_radius``,
exactly as in the fixed-probe mode; the enclosure probes only decide where the
channel is and where it ends. A fixed ``enclosure_radius`` keeps the
single-probe method of :func:`crevice.sections.connected_profile`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from .geometry import add, distance, dot, principal_axes, scale, sub
from .models import ChannelPoint, Coord, PoreProfile
from . import _kernels
from .spatial import SpatialIndex

#: Probes (Å) a plane may use for its enclosure test, smallest first.
ADAPTIVE_ENCLOSURE_LADDER = (
    0.30,
    0.35,
    0.40,
    0.45,
    0.50,
    0.55,
    0.60,
    0.65,
    0.70,
    0.75,
    0.80,
    0.90,
    1.00,
    1.10,
    1.20,
    1.40,
    1.60,
    1.80,
    2.00,
    2.40,
    3.00,
)
ADAPTIVE_ENCLOSURE_RULE = (
    "per-plane enclosure probe: in each plane the lumen centre is the in-plane clearance maximum climbed from the "
    "previous centre; the plane's probe is the smallest ladder probe below that clearance at which the free component "
    "holding the centre is laterally enclosed (search square, or a square twice as wide); adjacent centres join when "
    "the straight segment clears the smaller of their probes; the channel ends at the first plane, outwards from the "
    "seed, that is not enclosed or cannot be joined"
)
#: Axial step (Å) of the outward discovery trace and the largest output step.
DISCOVERY_STEP_A = 0.75
#: Mouth bracket width (Å).
MOUTH_TOLERANCE_A = 0.025
#: Search-square levels (multiples of ``search_radius``) a plane may use; a wide lumen gets a wider square.
WINDOW_LEVELS = (1, 2, 4)
#: A plane uses the smallest level whose half-width is at least this many times the lumen clearance.
WINDOW_PER_CLEARANCE = 2.0
#: A higher clearance hill replaces the climbed centre only when the straight segment between the two tops
#: keeps at least this fraction of the lower top's clearance (both tops lie in the same lumen).
HILL_RIDGE_FRACTION = 0.9
#: Coarse lattice stride (nodes) of the exact best-node search of a section (``_Tracer._best_node``).
BEST_NODE_STRIDE = 4
#: Largest local grid refinement (spacing divided by this) for planes narrower than the grid resolves.
MAX_GRID_REFINEMENT = 8
#: Climbs start from a point whose clearance may be slightly negative; moves only need to clear this.
_CLIMB_FLOOR_A = -1.0


@dataclass(frozen=True)
class Lumen:
    """The channel's section in one plane.

    Attributes
    ----------
    t : float
        Axial coordinate of the plane, Å.
    center : Coord
        Lumen centre: the local in-plane clearance maximum, Å.
    clearance : float
        Atom-surface clearance at ``center``, Å.
    probe : float
        Smallest ladder probe at which the lumen is laterally enclosed, Å.
    section : crevice.sections.Section
        The enclosed free component at ``probe`` (grid nodes).
    window_bridged : bool
        Enclosed only in the bridging square, twice as wide as the plane's search square.
    search_radius : float
        Half-width of the plane's search square, Å (``search_radius`` times its window level).
    section_spacing : float
        Grid spacing of the plane's sections, Å (finer than ``section_spacing`` in planes too narrow for it).
    enclosed : bool
        False for a resampled plane in a thin slab where the lumen is not enclosed (measured at its climbed
        centre; ``probe`` is then the previous plane's and ``section`` is None).
    joined : bool
        False when the straight segment from the previous sample does not clear the smaller probe.
    """

    t: float
    center: Coord
    clearance: float
    probe: float
    section: object
    window_bridged: bool = False
    search_radius: float = 0.0
    section_spacing: float = 0.0
    enclosed: bool = True
    joined: bool = True


def _smallest_true(values, test, hint=None):
    """``(value, test(value))`` for the smallest of the sorted ``values`` whose test is truthy, or None.

    The test is assumed monotone (false below some value, true from it on). The search starts at ``hint``
    (the largest value not above it), gallops down or up with doubling steps, then bisects, so a hint at the
    answer costs two tests instead of a full bisection.
    """
    if not values:
        return None
    memo = {}

    def check(k):
        if k not in memo:
            memo[k] = test(values[k])
        return memo[k]

    if hint is None:
        k = len(values) // 2
    else:
        k = max([i for i, v in enumerate(values) if v <= hint + 1e-12] or [0])
    if check(k):
        hi, step = k, 1
        lo = k - step
        while lo >= 0 and check(lo):
            hi, step = lo, step * 2
            lo = hi - step
        lo = max(lo, -1)
    else:
        lo, step = k, 1
        hi = k + step
        while hi < len(values) and not check(hi):
            lo, step = hi, step * 2
            hi = lo + step
        if hi >= len(values):
            if not check(len(values) - 1):
                return None
            hi = len(values) - 1
    # Invariant: check(hi) is true; every index <= lo is false (or lo = -1).
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if check(mid):
            hi = mid
        else:
            lo = mid
    return values[hi], memo[hi]


class _Tracer:
    """Per-plane lumen finder on one axis frame (shares one frame cache over all probes)."""

    def __init__(self, atoms, origin, direction, *, spacing, extent, ladder, refinement_steps, shared):
        from .sections import BRIDGE_WINDOW_FACTOR, SectionGeometry

        self.atoms, self.origin, self.direction = atoms, origin, direction
        self.spacing, self.extent, self.ladder = spacing, extent, tuple(ladder)
        self.steps = max(1, int(refinement_steps))
        self.climber = SectionGeometry(
            atoms, origin, direction, spacing=spacing, extent=extent, probe=_CLIMB_FLOOR_A, shared=shared
        )
        basis = (self.climber.u, self.climber.v)
        self.u, self.v = basis
        self.state = self.climber._state
        self._basis, self._shared = basis, shared
        self._levels = {}
        self.geometries = {p: self.geometry(p, 1) for p in self.ladder}
        self.window = BRIDGE_WINDOW_FACTOR * WINDOW_LEVELS[-1] * self.climber.half * spacing
        # Trace limits: the van der Waals envelope of the atoms along the axis. Beyond it a plane cuts no
        # atom, and only probe-inflated spheres could "enclose" it, which would push mouths into bulk.
        self.low = min(z - r for _, _, z, r in self.climber.projected)
        self.high = max(z + r for _, _, z, r in self.climber.projected)
        self._memo, self._graphs = {}, {}
        #: A plane may use a probe only this far (half a grid-cell diagonal) below its clearance: then the grid
        #: node nearest the centre is free and joined to it, and enclosure is monotone in the probe.
        self.grid_margin = spacing * math.sqrt(2.0) / 2

    def geometry(self, probe, level, refine=1):
        """Section geometry at ``probe`` in a square ``level`` times the base half-width, grid ``spacing/refine``."""
        from .sections import SectionGeometry

        key = (probe, level, refine)
        if key not in self._levels:
            self._levels[key] = SectionGeometry(
                self.atoms,
                self.origin,
                self.direction,
                spacing=self.spacing / refine,
                extent=self.extent * level,
                probe=probe,
                basis=self._basis,
                shared=self._shared,
            )
        return self._levels[key]

    def fits(self, level, refine):
        """Whether a square of this level and grid refinement stays within the 250000-node plane limit."""
        half = math.ceil(self.extent * level / (self.spacing / refine))
        return (2 * half + 1) ** 2 <= 250_000

    def level(self, clearance):
        """Window level for a lumen of this clearance: the base square while it is at least twice as wide."""
        for level in WINDOW_LEVELS:
            if WINDOW_PER_CLEARANCE * clearance <= self.extent * level:
                return level
        return WINDOW_LEVELS[-1]

    # -- one plane -----------------------------------------------------------------------------------------
    def plane_point(self, point, t):
        """``point`` moved along the axis onto the plane at ``t``."""
        return add(point, scale(self.direction, t - dot(sub(point, self.origin), self.direction)))

    def climb(self, point, t):
        """Local in-plane clearance maximum from ``point`` projected on the plane at ``t``."""
        import numpy as np

        start = self.plane_point(point, t)
        clearance = self.state.clearances([start])[0]
        base = np.array(add(self.origin, scale(self.direction, t)))
        u, v, w = np.array(self.u), np.array(self.v), self.window

        def inside(points):
            rel = np.asarray(points, dtype=float) - base
            return (np.abs(rel @ u) <= w) & (np.abs(rel @ v) <= w)

        return self.climber._climb(start, clearance, self.steps, inside)

    def _containing(self, geometry, sections, center):
        """The section among ``sections`` whose cells hold ``center`` with a probe-clear step to one of its nodes."""
        for section in sections:
            if not geometry._component_cells(section)([center])[0]:
                continue
            node = min(section.indices, key=lambda ij: distance(geometry.position(ij[0], ij[1], section.t), center))
            if self.state.segment_clearance(center, geometry.position(node[0], node[1], section.t)) > geometry.probe:
                return section
        return None

    def _labels(self, geometry, t):
        """Connected free components of one plane: ``(free, labels, exterior)`` (flat arrays), or None.

        NumPy path (without numba). The graph is that of :meth:`crevice.sections.SectionGeometry._plane_graph`:
        free nodes, and 8-neighbour edges followed when the closed segment clears the probe tested from both
        ends (the two tests differ only by rounding; an edge either end blocks is blocked, as in the compiled
        kernel). None without SciPy; the caller then uses ``sections``.
        """
        key = (geometry.probe, geometry.half, geometry.spacing, t)
        if key in self._graphs:
            return self._graphs[key]
        result = None
        try:
            import numpy as np
            from scipy.sparse import coo_matrix
            from scipy.sparse.csgraph import connected_components
        except ImportError:
            np = None
        if np is not None:
            from .sections import _STEPS

            h = geometry.half
            n = 2 * h + 1
            free_mask, blocked = geometry._plane_graph(t)
            free = free_mask.ravel()
            flat = blocked.reshape(len(_STEPS), -1)
            i_all, j_all = np.divmod(np.arange(n * n), n)
            rows, cols = [], []
            for k, (di, dj) in enumerate(_STEPS):
                if (di, dj) <= (0, 0):
                    continue
                back = _STEPS.index((-di, -dj))
                ok = free & (i_all + di < n) & (j_all + dj >= 0) & (j_all + dj < n)
                a = np.nonzero(ok)[0]
                b = a + di * n + dj
                keep = free[b]
                a, b = a[keep], b[keep]
                both = ~flat[k, a] & ~flat[back, b]
                rows.append(a[both])
                cols.append(b[both])
            rows = np.concatenate(rows) if rows else np.zeros(0, dtype=np.intp)
            cols = np.concatenate(cols) if cols else np.zeros(0, dtype=np.intp)
            graph = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n * n, n * n)).tocsr()
            _, labels = connected_components(graph, directed=False)
            edge = (i_all == 0) | (i_all == n - 1) | (j_all == 0) | (j_all == n - 1)
            exterior = np.zeros(labels.max() + 1, dtype=bool)
            exterior[labels[free & edge]] = True
            result = (free, labels, exterior)
        if len(self._graphs) > 4096:
            self._graphs.clear()
        self._graphs[key] = result
        return result

    def _component(self, geometry, t, i, j):
        """Grid indices ``(ii, jj)`` of the enclosed free component holding node ``(i, j)``, or None.

        Compiled kernel (:mod:`crevice._kernels`) when numba is available, otherwise the NumPy labelling of
        :meth:`_labels`; both give the same component.
        """
        import numpy as np

        h = geometry.half
        if abs(i) > h or abs(j) > h:
            return None
        if _kernels.HAVE_NUMBA:
            state = self.state
            free, blocked = _kernels.plane_graph(
                state.px,
                state.py,
                state.pz,
                state.pr,
                float(t),
                h,
                float(geometry.spacing),
                float(geometry.probe),
                _kernels.STEPS,
            )
            status, nodes = _kernels.component_from(free, blocked, h, i, j, _kernels.STEPS)
            if status != 1:
                return None
            order = np.lexsort((nodes[:, 1], nodes[:, 0]))
            return nodes[order, 0], nodes[order, 1]
        graph = self._labels(geometry, t)
        if graph is None:
            return False
        free, labels, exterior = graph
        n = 2 * h + 1
        if not free[(i + h) * n + (j + h)] or exterior[labels[(i + h) * n + (j + h)]]:
            return None
        members = np.nonzero((labels == labels[(i + h) * n + (j + h)]) & free)[0]
        return members // n - h, members % n - h

    def _section_at(self, geometry, t, center, small=None):
        """The enclosed :class:`~crevice.sections.Section` holding ``center`` at this geometry's probe, or None.

        ``center`` belongs to the component of the nearest corner node of its grid cell (free, and joined to
        ``center`` by a probe-clear straight segment, whenever the probe is at least ``grid_margin`` below the
        centre's clearance; otherwise the result is None). The section
        centre is the component's node of largest clearance (ties as in ``sections``).
        """
        import numpy as np
        from .sections import Section

        ds = geometry.spacing
        rel = sub(center, add(self.origin, scale(self.direction, t)))
        x, y = dot(rel, geometry.u) / ds, dot(rel, geometry.v) / ds
        a, b = math.floor(x), math.floor(y)
        # The nearest corner of the centre's cell; with the grid margin it is free and joined to the centre.
        i, j = min(
            ((i, j) for i in (a, a + 1) for j in (b, b + 1)), key=lambda ij: ((ij[0] - x) ** 2 + (ij[1] - y) ** 2, ij)
        )
        if self.state.segment_clearance(center, geometry.position(i, j, t)) <= geometry.probe:
            return None
        component = self._component(geometry, t, i, j)
        if component is False:
            return self._containing(geometry, geometry.sections(t), center)
        if component is None:
            return None
        ii, jj = component
        best, value = self._best_node(geometry, t, ii, jj)
        indices = tuple(zip(ii.tolist(), jj.tolist()))
        # ``small``: half-width (nodes) of the plane's own square; reaching it marks the section as bridged.
        bridged = small is not None and bool((np.maximum(np.abs(ii), np.abs(jj)) >= small).any())
        return Section(t, indices, geometry.position(int(ii[best]), int(jj[best]), t), value, bridged)

    def _best_node(self, geometry, t, ii, jj):
        """Index and clearance of the component node of largest clearance (ties as in ``sections``).

        Exact without evaluating every node: clearance is 1-Lipschitz, so a node's clearance is at most that
        of the coarse lattice point (every ``BEST_NODE_STRIDE``-th node) nearest to it plus their distance.
        Coarse points are evaluated first; a node is evaluated only when that bound reaches the best value
        found at a coarse point that belongs to the component (a 1e-9 margin keeps every node that could tie).
        """
        import numpy as np

        ds = geometry.spacing
        base = np.array(add(self.origin, scale(self.direction, t)))
        u, v = np.array(geometry.u), np.array(geometry.v)

        def values_at(i, j):
            return self.state.spatial.nearest_surface_values(base + (i * ds)[:, None] * u + (j * ds)[:, None] * v)

        stride = BEST_NODE_STRIDE
        evaluate = np.ones(len(ii), dtype=bool)
        if len(ii) > 4 * stride * stride:
            ci, cj = np.rint(ii / stride).astype(np.int64) * stride, np.rint(jj / stride).astype(np.int64) * stride
            span = int(max(np.abs(ci).max(), np.abs(cj).max())) + 1
            keys, owner = np.unique((ci + span) * (2 * span + 1) + (cj + span), return_inverse=True)
            owner = owner.ravel()
            coarse = values_at(keys // (2 * span + 1) - span, keys % (2 * span + 1) - span)
            member = (ii % stride == 0) & (jj % stride == 0)
            if member.any():
                floor = coarse[owner[member]].max()
                bound = coarse[owner] + ds * np.hypot(ii - ci, jj - cj)
                evaluate = bound >= floor - 1e-9
        chosen = np.flatnonzero(evaluate)
        si, sj = ii[chosen], jj[chosen]
        values = values_at(si, sj)
        best = int(np.lexsort((-sj, -si, -(si * si + sj * sj), np.round(values, 10)))[-1])
        return int(chosen[best]), float(values[best])

    def enclosed_at(self, t, center, probe, level=1, refine=1):
        """``(section, bridged)`` when the component holding ``center`` is enclosed at ``probe``, else None.

        Enclosed in the bridging square, twice the plane's search square (``level`` times
        ``search_radius``), on a grid of ``section_spacing / refine``. That square contains the plane's own
        square node for node, so a component enclosed in the smaller square is enclosed, unchanged, in the
        larger one: testing the larger square alone gives the same answer as testing the smaller one first.
        ``bridged`` says the component reaches the smaller square's edge.
        """
        from .sections import BRIDGE_WINDOW_FACTOR

        wide = int(BRIDGE_WINDOW_FACTOR * level)
        if not _kernels.HAVE_NUMBA or not self.fits(wide, refine):
            # NumPy path: the plane's own square first (four times fewer nodes), then the wider one.
            section = self._section_at(self.geometry(probe, level, refine), t, center)
            if section is not None or not self.fits(wide, refine):
                return (section, False) if section is not None else None
            section = self._section_at(self.geometry(probe, wide, refine), t, center)
            return (section, True) if section is not None else None
        section = self._section_at(
            self.geometry(probe, wide, refine), t, center, small=self.geometry(probe, level, refine).half
        )
        if section is None:
            return None
        return replace(section, window_bridged=False), section.window_bridged

    def lumen(self, t, start, hint=None):
        """The :class:`Lumen` of the plane at ``t`` climbed from ``start``, or ``(None, centre, clearance)``.

        ``hint`` (the previous plane's probe) only orders the search; the result is the smallest enclosing
        ladder probe either way, because enclosure is monotone in the probe.
        """
        key = (t, start)
        if key in self._memo:
            return self._memo[key]
        center, clearance = self.climb(start, t)
        level, refine = self.level(clearance), 1
        usable = [p for p in self.ladder if p < clearance - self.grid_margin - 1e-9]
        # A plane narrower than the grid can resolve (no probe half a grid diagonal below its clearance)
        # is re-gridded locally at half the spacing, up to MAX_GRID_REFINEMENT times.
        while (
            not usable and clearance > self.ladder[0] and refine < MAX_GRID_REFINEMENT and self.fits(level, 2 * refine)
        ):
            refine *= 2
            margin = self.grid_margin / refine
            usable = [p for p in self.ladder if p < clearance - margin - 1e-9]
        found = _smallest_true(usable, lambda p: self.enclosed_at(t, center, p, level, refine), hint)
        if found:
            probe, (section, bridged) = found
            # The plane's widest point in this enclosed component: also climb from the component's best grid
            # node (a wide lumen can hold several clearance hills) and keep the higher hill.
            geometry = self.geometry(probe, level, refine)
            other, other_clearance = self.climb(section.center, t)
            # Only a hill of the same lumen: the straight segment between the two tops stays above
            # HILL_RIDGE_FRACTION of the lower top's clearance (a neighbouring passage joined only through a
            # narrower neck is not taken).
            if (
                other_clearance > clearance + 1e-9
                and geometry._component_cells(section)([other])[0]
                and self.state.segment_clearance(center, other) > HILL_RIDGE_FRACTION * clearance
            ):
                center, clearance = other, other_clearance
            result = Lumen(t, center, clearance, probe, section, bridged, self.extent * level, self.spacing / refine)
        else:
            result = (None, center, clearance)
        self._memo[key] = result
        return result

    def joins(self, a, b):
        """Whether centres ``a`` and ``b`` join (the segment clears the smaller probe)."""
        return self.state.segment_clearance(a.center, b.center) > min(a.probe, b.probe)

    def step(self, previous, t):
        """``(lumen, None)`` continuing ``previous`` into the plane at ``t``, or ``(None, reason)``."""
        found = self.lumen(t, previous.center, previous.probe)
        if not isinstance(found, Lumen):
            return None, {"reason": "not_enclosed", "clearance_A": found[2]}
        if not self.joins(previous, found):
            return None, {
                "reason": "join_blocked",
                "clearance_A": found.clearance,
                "segment_clearance_A": self.state.segment_clearance(previous.center, found.center),
            }
        return found, None

    # -- outward traces --------------------------------------------------------------------------------------
    def outward(self, seed, sign, ts):
        """Follow the lumen from ``seed`` through the planes ``ts`` (ordered outwards); stop at the first failure."""
        path, previous = [], seed
        for t in ts:
            lumen, failure = self.step(previous, t)
            if lumen is None:
                return path, t, failure
            path.append(lumen)
            previous = lumen
        return path, None, None

    def mouth(self, inside, outside_t):
        """Bracket the mouth between the enclosed ``inside`` lumen and the failing plane at ``outside_t``."""
        failure = None
        while abs(outside_t - inside.t) > MOUTH_TOLERANCE_A:
            middle = (outside_t + inside.t) / 2
            lumen, why = self.step(inside, middle)
            if lumen is None:
                outside_t, failure = middle, why
            else:
                inside = lumen
        return inside, {
            "coordinate_A": (inside.t + outside_t) / 2,
            "bracket_A": sorted((inside.t, outside_t)),
            "enclosed_coordinate_A": inside.t,
            "open_coordinate_A": outside_t,
            "axial_bracket_width_A": abs(outside_t - inside.t),
            "reason": (failure or {}).get("reason", "not_enclosed"),
        }

    def discover(self, seed_point):
        """Discovery trace at ``DISCOVERY_STEP_A`` from the seed plane; ``None`` if the seed plane is open.

        Returns ``(lower, lower_mouth, upper, upper_mouth, seed_lumen)``.
        """
        seed = self.lumen(0.0, seed_point)
        if not isinstance(seed, Lumen):
            return None
        ends, sides = [], []
        for sign, limit in ((-1, self.low), (1, self.high)):
            n = int(math.floor(abs(limit) / DISCOVERY_STEP_A))
            ts = [sign * DISCOVERY_STEP_A * k for k in range(1, n + 1)]
            if not ts or abs(ts[-1] - limit) > 1e-9:
                ts.append(limit)
            path, failed, _ = self.outward(seed, sign, ts)
            sides.append(path)
            last = path[-1] if path else seed
            if failed is None:
                ends.append(
                    (
                        last,
                        {
                            "coordinate_A": last.t,
                            "bracket_A": [last.t, last.t],
                            "enclosed_coordinate_A": last.t,
                            "open_coordinate_A": None,
                            "axial_bracket_width_A": 0.0,
                            "reason": "atom_envelope",
                        },
                    )
                )
            else:
                ends.append(self.mouth(last, failed))
        (lower, lower_mouth), (upper, upper_mouth) = ends
        #: Discovery lumens in axial order, from the lower mouth to the upper mouth (used by the lateral-branch
        #: rule, ``_lateral_branches``).
        planes = {}
        for lumen in (lower, *sides[0], seed, *sides[1], upper):
            if lower.t <= lumen.t <= upper.t:
                planes.setdefault(lumen.t, lumen)
        self.discovery_planes = [planes[t] for t in sorted(planes)]
        return lower, lower_mouth, upper, upper_mouth, seed

    def end_margin(self, lumen, sign):
        """Clearance of the straight axial continuation from an end centre beyond the atom envelope.

        The target lies one ångström beyond the probe-inflated atom envelope,
        as in the fixed-probe mode (``crevice.sections._open_ends``).
        """
        projected = self.climber.projected
        if sign < 0:
            target = min(z - r for _, _, z, r in projected) - lumen.probe - 1
        else:
            target = max(z + r for _, _, z, r in projected) + lumen.probe + 1
        end = add(lumen.center, scale(self.direction, target - lumen.t))
        return self.state.segment_clearance(lumen.center, end)


def _resample(tracer, seed_point, lower, upper, samples):
    """Trace evenly spaced planes between the discovery mouths, which never move.

    The extent of the channel is fixed by the discovery trace and the bisected mouths, so it does not depend
    on the number of samples. A resampled plane that falls in a thin slab where the lumen is not enclosed (or
    whose centre cannot be joined to the previous one) is still measured, at its climbed centre, and is
    flagged (``enclosed=False`` / ``joined=False``) instead of truncating the channel.

    Returns the list of :class:`Lumen` samples.
    """
    low, high = lower.t, upper.t
    count = max(samples, int(math.ceil((high - low) / DISCOVERY_STEP_A)) + 1)
    ts = [low + (high - low) * i / (count - 1) for i in range(count)]
    k0 = min(range(count), key=lambda k: (abs(ts[k]), k))

    def sample(t, start, previous):
        found = tracer.lumen(t, start, previous.probe if previous is not None else None)
        if not isinstance(found, Lumen):
            _, center, clearance = found
            probe = previous.probe if previous is not None else tracer.ladder[0]
            found = Lumen(t, center, clearance, probe, None, False, 0.0, 0.0, enclosed=False)
        if previous is not None and not tracer.joins(previous, found):
            found = replace(found, joined=False)
        return found

    seed = sample(ts[k0], seed_point, None)
    path = [seed]
    for direction in (1, -1):
        previous, side = seed, []
        for t in ts[k0 + 1 :] if direction > 0 else ts[:k0][::-1]:
            previous = sample(t, previous.center, previous)
            side.append(previous)
        path = path + side if direction > 0 else side[::-1] + path
    return path


def _trace_direction(
    atoms, origin, direction, seed_point, *, samples, search_radius, section_spacing, refinement_steps, ladder, shared
):
    """Trace one channel along ``direction`` from ``seed_point``; None when the seed plane is not enclosed."""
    tracer = _Tracer(
        atoms,
        origin,
        direction,
        spacing=section_spacing,
        extent=search_radius,
        ladder=ladder,
        refinement_steps=refinement_steps,
        shared=shared,
    )
    found = tracer.discover(seed_point)
    if found is None:
        return None
    lower, lower_mouth, upper, upper_mouth, seed = found
    path = _resample(tracer, seed.center, lower, upper, samples)
    return {
        "tracer": tracer,
        "path": path,
        "mouths": (lower_mouth, upper_mouth),
        "adjustments": [],
        "seed": seed,
        "span": path[-1].t - path[0].t,
    }


def _seed_points(atoms, origin, direction, *, search_radius, section_spacing, ladder, shared):
    """Candidate seeds on the plane through ``origin`` when no origin was supplied: centres of enclosed sections."""
    from .sections import SectionGeometry

    seeds = []
    for probe in ladder:
        geometry = SectionGeometry(
            atoms, origin, direction, spacing=section_spacing, extent=search_radius, probe=probe, shared=shared
        )
        for section in geometry.sections(0.0):
            if all(distance(section.center, s) > 2 * section_spacing for s in seeds):
                seeds.append(section.center)
    return seeds


def adaptive_profile(
    atoms,
    *,
    axis,
    origin,
    samples,
    search_radius,
    refinement_steps,
    probe_radius,
    section_spacing=0.5,
    lateral_exits=False,
    exit_bulk_radius=6.0,
    exit_spacing=0.5,
    ladder=ADAPTIVE_ENCLOSURE_LADDER,
    fixed=False,
):
    """Resolve a channel with a per-plane enclosure probe (see the module docstring).

    Called by :func:`crevice.sections.connected_profile` for the automatic
    probe (``enclosure_radius=None``); use :func:`crevice.pore_profile`.

    Parameters
    ----------
    atoms : sequence of Atom
        Obstacle atoms (already selected).
    axis : str or Coord
        ``"auto"`` (try the principal directions and keep a unique channel),
        ``"x"``, ``"y"``, ``"z"`` or a direction vector.
    origin : Coord or None
        Seed inside the intended channel, Å; its clearance must exceed the
        smallest usable probe. ``None``: seeds are the enclosed sections of
        the plane through the atom centroid.
    samples : int
        Minimum number of output samples; more are used so that axial steps
        stay at most 0.75 Å.
    search_radius : float
        Base half-width of the transverse search square, Å; a plane whose
        lumen is wide uses 2 or 4 times this (``WINDOW_LEVELS``).
    refinement_steps : int
        Minimum step halvings of each in-plane clearance climb (at least 1).
    probe_radius : float
        Measurement probe radius, Å, subtracted from each reported radius.
    section_spacing : float, default 0.5
        Transverse grid spacing, Å. A plane can use a probe only half a grid
        diagonal below its clearance; a plane too narrow for the grid is
        re-gridded locally at up to ``MAX_GRID_REFINEMENT`` times finer, so
        the narrowest traceable clearance is the smallest ladder probe plus
        ``0.71 * section_spacing / MAX_GRID_REFINEMENT`` (0.32 Å at 0.25 Å,
        0.34 Å at 0.5 Å).
    lateral_exits : bool, default False
        Accept a capped end that has a lateral exit to bulk.
    exit_bulk_radius : float, default 6.0
        With ``lateral_exits``: rolling bulk-probe radius, Å.
    exit_spacing : float, default 0.5
        With ``lateral_exits``: lattice spacing, Å, of the exit search.
    ladder : sequence of float, default ADAPTIVE_ENCLOSURE_LADDER
        Probes a plane may use, Å (values below ``probe_radius`` are raised
        to it).
    fixed : bool, default False
        ``ladder`` is one explicitly chosen probe (``enclosure_radius``):
        recorded as ``enclosure_mode = "fixed"`` without a selection record.

    Returns
    -------
    PoreProfile
        Method ``"axial-connected"`` with the metadata of the fixed-probe mode,
        plus ``enclosure_mode = "adaptive"``, ``enclosure_probes_A`` (each
        sample's probe), ``enclosure_probe_range_A`` and
        ``enclosure_probe_selection`` (``mode = "adaptive"``, the ladder and
        the rule). ``enclosure_radius_A`` is the largest per-plane probe.

    Raises
    ------
    crevice.sections.ChannelResolutionError
        When the seed plane is not enclosed at any probe, the axis or region
        is ambiguous, or an end is capped without ``lateral_exits``.
    crevice.sections.OriginClearanceError
        For a supplied ``origin`` whose clearance does not exceed the smallest
        usable probe.
    """
    from .geometry import axis_from_name
    from .sections import (
        _ACTIVE_FRAMES,
        ChannelResolutionError,
        OriginClearanceError,
        _FrameCache,
        _refinement_metadata,
    )

    ladder = tuple(sorted({max(float(probe_radius), float(q)) for q in ladder}))
    inferred_origin, explicit, label = axis_from_name(axis, atoms)
    supplied = origin is not None
    origin = inferred_origin if origin is None else tuple(origin)
    if supplied:
        clearance = SpatialIndex(atoms).nearest_surface(origin)[0]
        if clearance <= ladder[0]:
            raise OriginClearanceError(
                "Channel origin must have clearance greater than the enclosure and measurement "
                f"probes (per-plane enclosure probe: the smallest candidate is {ladder[0]:.2f} A; "
                f"clearance at the origin {clearance:.2f} A)",
                clearance,
            )
    automatic = isinstance(axis, str) and axis.lower() == "auto"
    directions = principal_axes(atoms) if automatic else (explicit,)
    token = _ACTIVE_FRAMES.set(_ACTIVE_FRAMES.get() or _FrameCache())
    try:
        shared = _ACTIVE_FRAMES.get()
        candidates, diagnostics = [], []
        for direction in directions:
            seeds = (
                [origin]
                if supplied
                else _seed_points(
                    atoms,
                    origin,
                    direction,
                    search_radius=search_radius,
                    section_spacing=section_spacing,
                    ladder=ladder,
                    shared=shared,
                )
            )
            traces = [
                r
                for r in (
                    _trace_direction(
                        atoms,
                        origin,
                        direction,
                        s,
                        samples=samples,
                        search_radius=search_radius,
                        section_spacing=section_spacing,
                        refinement_steps=refinement_steps,
                        ladder=ladder,
                        shared=shared,
                    )
                    for s in seeds
                )
                if r is not None and len(r["path"]) >= 3
            ]
            # The same channel reached from two seeds: keep one (longest, then widest, then first).
            distinct = []
            for trace in sorted(
                traces, key=lambda r: (-round(r["span"], 9), -round(min(p.clearance for p in r["path"]), 9))
            ):
                mid = trace["path"][len(trace["path"]) // 2]
                if all(
                    distance(min(o["path"], key=lambda p: abs(p.t - mid.t)).center, mid.center) > 2 * section_spacing
                    for o in distinct
                ):
                    distinct.append(trace)
            best = distinct[0] if distinct else None
            diagnostics.append(
                {
                    "direction": direction,
                    "seeds": len(seeds),
                    "channels": len(distinct),
                    "through_span_A": best["span"] if best else 0.0,
                }
            )
            if best:
                ambiguous = any(o["span"] >= 0.9 * best["span"] for o in distinct[1:])
                candidates.append((best, ambiguous))
        if not candidates:
            raise ChannelResolutionError(
                "No connected channel was resolved: the seed plane's lumen is not laterally enclosed at any "
                f"per-plane enclosure probe ({ladder[0]:.2f}-{ladder[-1]:.2f} A). Supply a channel axis/origin, "
                "increase search_radius or reduce section_spacing; use cavity analysis for closed voids"
            )
        candidates.sort(key=lambda row: -row[0]["span"])
        if len(candidates) > 1 and candidates[1][0]["span"] >= 0.9 * candidates[0][0]["span"]:
            raise ChannelResolutionError(
                "Ambiguous channel axis: multiple directions have comparable connected spans; "
                "supply an explicit axis and origin"
            )
        best, ambiguous = candidates[0]
        if ambiguous:
            raise ChannelResolutionError(
                "Ambiguous channel region: multiple disconnected passages; supply an origin inside the intended channel"
            )
        options = dict(
            samples=samples,
            search_radius=search_radius,
            section_spacing=section_spacing,
            refinement_steps=refinement_steps,
            ladder=ladder,
            lateral_exits=lateral_exits,
            exit_bulk_radius=exit_bulk_radius,
            exit_spacing=exit_spacing,
            refinement=_refinement_metadata(refinement_steps),
            fixed=fixed,
        )
        if lateral_exits:
            branched = _lateral_branches(best, samples=samples, bulk_radius=exit_bulk_radius, spacing=exit_spacing)
            if branched is not best:
                try:
                    return _build_profile(branched, origin, probe_radius, label, automatic, diagnostics, **options)
                except ChannelResolutionError:
                    # The exit search from the new end found no leg (its resampled centre differs from the
                    # discovery centre the rule searched from): keep the channel as traced.
                    pass
        return _build_profile(best, origin, probe_radius, label, automatic, diagnostics, **options)
    except ChannelResolutionError as exc:
        # A refusal of the per-plane probe carries its record, so callers can report why.
        if getattr(exc, "probe_selection", None) is None and not fixed:
            exc.probe_selection = {
                "mode": "adaptive",
                "rule": ADAPTIVE_ENCLOSURE_RULE,
                "ladder_A": list(ladder),
                "probe_radius_A": probe_radius,
                "chosen_A": None,
                "probe_range_A": None,
                "reason": str(exc).split(". ")[0].split(": ")[0],
            }
        raise
    finally:
        _ACTIVE_FRAMES.reset(token)


#: Bulk-envelope grid spacing (Å) of the lateral-exit search (``lateral_exit_legs`` default).
_EXIT_BULK_SPACING_A = 1.0


def _exit_probes(ladder, end_probe):
    """Probes for a capped end's lateral-exit search: the end probe, about half of it, the smallest ladder probe."""
    below = [q for q in ladder if q <= end_probe]
    half = min(below, key=lambda q: (abs(q - end_probe / 2), q))
    return tuple(sorted({end_probe, half, below[0]}, reverse=True))


#: A lateral route from the lumen must be wider than the rest of the axial path by more than this (Å) to end the
#: lumen there (``_lateral_branches``).
BRANCH_MARGIN_A = 0.1
#: ... or by more than this fraction of the continuation's width, whichever is larger.
BRANCH_FRACTION = 0.1
#: The branch searches use a lattice no finer than the probe divided by this (and no finer than the exit spacing).
BRANCH_LATTICE_DIVISOR = 4.0
LATERAL_BRANCH_RULE = (
    "with lateral exits, a capped end moves to the plane from which the widest lateral route to bulk leaves the "
    "lumen, when that route is clearly wider than the axial continuation beyond the plane: among the discovery planes "
    "beyond both the seed plane and the channel's narrowest discovery plane, the continuation's width beyond a plane "
    "is the smallest clearance of the discovery planes beyond it and of the end's own exit (its widest leg; none "
    "without one), all widths capped at the bulk-probe radius; a plane qualifies when the lateral-exit search from "
    "its centre into the half-space beyond it finds a route to bulk clearing that width plus max(BRANCH_MARGIN_A, "
    "BRANCH_FRACTION of it); the end becomes the qualifying plane with the widest route. Open ends never move"
)
BULK_END_DEFINITION = (
    "the end sample already lies in bulk solvent (inside the rolling bulk-probe envelope beyond the mouth): an open "
    "end, although the straight continuation along the axis does not clear the end plane's enclosure probe"
)


def _end_envelope(tracer, lumen, sign, bulk_radius, t_limit=None):
    """Bulk envelope beyond the plane ``t_limit`` (default the end plane) on the ``sign`` side."""
    from .channel_exits import BulkEnvelope

    geometry = tracer.climber
    return BulkEnvelope(
        geometry.atoms,
        geometry.origin,
        (tracer.u, tracer.v, tracer.direction),
        lumen.t if t_limit is None else t_limit,
        sign,
        spacing=_EXIT_BULK_SPACING_A,
        bulk_radius=bulk_radius,
    )


def _end_in_bulk(envelope, lumen):
    """Whether the end sample's centre lies in bulk (``BulkEnvelope.depth == 0``)."""
    return envelope.depth(lumen.center) == 0.0


def _end_exit_search(tracer, lumen, sign, label, *, bulk_radius, spacing, envelope):
    """The usual lateral-exit search of a capped end: ``(legs, info, tried, exit_probe)``."""
    from .channel_exits import lateral_exit_legs

    if lumen.section is not None:
        section = replace(lumen.section, center=lumen.center, clearance=lumen.clearance)
    else:
        from .sections import Section

        section = Section(lumen.t, (), lumen.center, lumen.clearance)
    # The exit must clear a probe: the end plane's own probe first, then smaller ladder probes (an end plane
    # that is enclosed only at a large probe can still drain through a narrower lateral route).
    tried = _exit_probes(tracer.ladder, lumen.probe)
    legs, info = [], {}
    for exit_probe in tried:
        geometry = tracer.geometry(exit_probe, tracer.level(lumen.clearance))
        legs, info = lateral_exit_legs(
            geometry,
            section,
            sign,
            end=label,
            spacing=spacing,
            bulk_radius=bulk_radius,
            bulk_spacing=_EXIT_BULK_SPACING_A,
            envelope=envelope,
        )
        if legs:
            break
    return legs, info, tried, exit_probe


def _classify_end(tracer, lumen, sign, label, *, bulk_radius, spacing, forced_lateral=False, cache=None):
    """``(kind, margin, payload)`` of one end: ``"axial"``, ``"bulk"`` or ``"lateral"`` (payload: the search).

    Results are cached per end lumen (``cache``), so the lateral-branch rule and the exit metadata search once.
    """
    key = (label, lumen.t, lumen.center, forced_lateral)
    if cache is not None and key in cache:
        return cache[key]
    margin = tracer.end_margin(lumen, sign)
    if margin > lumen.probe and not forced_lateral:
        result = ("axial", margin, None)
    else:
        # The bulk envelope does not depend on the probe: one per end, reused by every probe tried.
        envelope = _end_envelope(tracer, lumen, sign, bulk_radius)
        if not forced_lateral and _end_in_bulk(envelope, lumen):
            result = ("bulk", margin, None)
        else:
            result = (
                "lateral",
                margin,
                _end_exit_search(
                    tracer, lumen, sign, label, bulk_radius=bulk_radius, spacing=spacing, envelope=envelope
                ),
            )
    if cache is not None:
        cache[key] = result
    return result


def _exit_metadata(tracer, path, segments, probe_radius, *, bulk_radius, spacing, branches=None, cache=None):
    """Exit classification of both ends, each searched at its own end probe (lateral exits requested)."""
    from .sections import ChannelResolutionError

    if not (math.isfinite(spacing) and spacing > 0):
        raise ValueError("exit_spacing must be finite and positive")
    branches = branches or {}
    exits, widths = {}, [min(segments)]
    for label, lumen, sign in (("lower", path[0], -1), ("upper", path[-1], 1)):
        if not (math.isfinite(bulk_radius) and bulk_radius > lumen.probe):
            raise ValueError("exit_bulk_radius must be finite and larger than the enclosure probe")
        kind, margin, found = _classify_end(
            tracer,
            lumen,
            sign,
            label,
            bulk_radius=bulk_radius,
            spacing=spacing,
            forced_lateral=label in branches,
            cache=cache,
        )
        if kind == "axial":
            exits[label] = {
                "type": "axial",
                "straight_exit_clearance_A": margin,
                "enclosure_probe_A": lumen.probe,
                "definition": "unobstructed straight continuation along the axis beyond the atom envelope",
            }
            continue
        if kind == "bulk":
            exits[label] = {
                "type": "bulk",
                "straight_exit_clearance_A": margin,
                "enclosure_probe_A": lumen.probe,
                "definition": BULK_END_DEFINITION,
            }
            continue
        legs, info, tried, exit_probe = found
        if not legs:
            raise ChannelResolutionError(
                f"The {label} end is capped (straight exit clearance {margin:.2f} A) and no lateral exit "
                f"to bulk clears a probe of {', '.join(f'{q:.2f}' for q in tried)} A; the channel is not resolved"
            )
        info = {**info, "exit_probe_A": exit_probe, "exit_probes_tried_A": list(tried[: tried.index(exit_probe) + 1])}
        records = [leg.to_dict(probe_radius) for leg in legs]
        widest = max(legs, key=lambda leg: leg.bottleneck_clearance)
        exits[label] = {
            "type": "lateral",
            "end": label,
            "straight_exit_clearance_A": margin,
            "enclosure_probe_A": lumen.probe,
            "leg_count": len(records),
            "bottleneck_radius_A": max(0.0, widest.bottleneck_clearance - probe_radius),
            "bottleneck_raw_clearance_A": widest.bottleneck_clearance,
            **info,
            **({"lateral_branch": branches[label]} if label in branches else {}),
            "legs": records,
        }
        widths.append(widest.bottleneck_clearance)
    lateral = [k for k, e in exits.items() if e["type"] == "lateral"]
    return {
        "exits": exits,
        "path_type": "capped_channel_lateral_exit" if lateral else "through_channel",
        "path_bottleneck_radius_A": max(0.0, min(widths) - probe_radius),
        "path_bottleneck_definition": (
            "smallest of the axial continuous bottleneck and, for each capped end, "
            "the bottleneck of its widest lateral-exit leg"
        ),
        "exit_policy": {
            "lateral_exits": True,
            "bulk_radius_A": bulk_radius,
            "exit_grid_spacing_A": spacing,
            "capped_ends": lateral,
        },
    }


def _lateral_branches(trace, *, samples, bulk_radius, spacing):
    """Apply the lateral-branch rule (``LATERAL_BRANCH_RULE``) to both ends of a traced channel.

    Only a capped end (its straight continuation blocked, so it leaves laterally in any case) moves; open ends,
    ends in bulk and through channels are never changed. For a capped end, the discovery planes beyond both the
    seed plane and the channel's narrowest discovery plane (the lumen always keeps its seed and its main
    constriction) up to the end's own sample are examined. ``W_k``, the width of the axial continuation beyond
    plane ``k``, is the smallest of the clearances of the planes beyond it and the end's exit width (its widest
    leg, none for an end without an exit), capped at ``bulk_radius`` (any passage at least as wide
    as the bulk probe is bulk-wide). A lateral route from plane ``k`` is the lateral-exit search from its centre
    into the half-space beyond it (first leg only), at probe ``q_k = W_k + max(BRANCH_MARGIN_A,
    BRANCH_FRACTION * W_k)`` (at least the smallest ladder probe). Every node and segment of such a route clears
    ``q_k``, so it cannot pass any plane beyond ``k`` inside the lumen: it is lateral. Planes are tested where
    their clearance exceeds ``q_k`` by a lattice-cell diagonal (the search's resolution; the lattice spacing is
    the exit spacing, or ``q_k / BRANCH_LATTICE_DIVISOR`` when that is coarser). The end becomes the tested plane
    with the widest lateral route (searched with the probe raised to the best width found so far). Consecutive
    candidate planes whose centres are joined by segments clearing the innermost one's probe are first tested
    together, by one search from the innermost plane (a route from any of them is a route from it).
    Bulk for these searches is the rolling-probe bulk of the half-space beyond the innermost plane examined (the
    outer of the seed and narrowest planes).

    Returns
    -------
    dict
        ``trace`` itself when no end changes; otherwise a new trace (resampled between the new mouths) with
        ``branches`` (per changed end: plane, widths, probe) and the shared ``exit_cache``.
    """
    from types import SimpleNamespace

    from .channel_exits import lateral_exit_legs
    from .sections import Section

    tracer, path = trace["tracer"], trace["path"]
    planes = getattr(tracer, "discovery_planes", None)
    if not planes or len(planes) < 3:
        return trace
    cache = trace.setdefault("exit_cache", {})
    seed_index = next((i for i, plane in enumerate(planes) if plane.t == 0.0), None)
    if seed_index is None:
        return trace
    narrowest = min(range(len(planes)), key=lambda i: (planes[i].clearance, abs(i - seed_index)))
    base = tracer.climber
    new_ends, branches, mouths = {}, {}, list(trace["mouths"])
    for index, (label, sign, end) in enumerate((("lower", -1, path[0]), ("upper", 1, path[-1]))):
        # Only a capped end moves: an open end already leaves along the axis.
        kind, _, found = _classify_end(tracer, end, sign, label, bulk_radius=bulk_radius, spacing=spacing, cache=cache)
        if kind != "lateral":
            continue
        legs = found[0]
        exit_width = max(leg.bottleneck_clearance for leg in legs) if legs else -math.inf
        # The planes beyond both the seed plane and the channel's narrowest discovery plane (excluded): the lumen
        # always keeps its seed and its main constriction.
        inner = min(seed_index, narrowest) if sign < 0 else max(seed_index, narrowest)
        side = planes[:inner][::-1] if sign < 0 else planes[inner + 1 :]
        if len(side) < 2:
            continue
        side = side[:-1] + [end]  # inner to outer; the mouth plane as the profile's own end sample

        # Widths at least as large as the bulk probe are all bulk-wide: they are compared capped at bulk_radius.
        def cap(width):
            return min(width, bulk_radius)

        state = {"envelope": None}
        searched = {}

        def search_spacing(probe):
            # Lattice resolution proportional to the width tested (the exit spacing at most, for narrow routes).
            return max(spacing, probe / BRANCH_LATTICE_DIVISOR)

        def route(k, probe):
            """Width of the widest lateral route from plane ``k`` clearing ``probe``, or None."""
            if (k, probe) not in searched:
                if state["envelope"] is None:
                    state["envelope"] = _end_envelope(tracer, end, sign, bulk_radius, t_limit=planes[inner].t)
                plane = side[k]
                geometry = SimpleNamespace(
                    spatial=base.spatial,
                    probe=probe,
                    u=tracer.u,
                    v=tracer.v,
                    direction=tracer.direction,
                    atoms=base.atoms,
                    origin=base.origin,
                )
                try:
                    legs, _ = lateral_exit_legs(
                        geometry,
                        Section(plane.t, (), plane.center, plane.clearance),
                        sign,
                        end=label,
                        spacing=search_spacing(probe),
                        bulk_radius=bulk_radius,
                        bulk_spacing=_EXIT_BULK_SPACING_A,
                        max_legs=1,
                        envelope=state["envelope"],
                    )
                except ValueError:  # the search exceeded its node budget before reaching bulk
                    legs = []
                searched[(k, probe)] = legs[0].bottleneck_clearance if legs else None
            return searched[(k, probe)]

        # Width of the axial continuation beyond each plane (pure: plane clearances outward and the end's exit),
        # and the probe a lateral route from that plane must clear to be clearly wider.
        following, probes, continuation = cap(min(end.clearance, exit_width)), {}, {}
        for k in range(len(side) - 2, -1, -1):
            continuation[k] = following
            probes[k] = max(following + max(BRANCH_MARGIN_A, BRANCH_FRACTION * following), tracer.ladder[0])
            following = min(cap(side[k].clearance), following)
        # Candidates: planes wider than their probe by a lattice-cell diagonal (the exit search's resolution).
        candidates = [
            k
            for k in probes
            if probes[k] < bulk_radius
            and cap(side[k].clearance) > probes[k] + math.sqrt(3.0) * search_spacing(probes[k])
        ]
        runs = []
        for k in sorted(candidates, reverse=True):  # outer to inner
            if runs and runs[-1][-1] == k + 1:
                runs[-1].append(k)
            else:
                runs.append([k])
        best = None  # (width, plane index)
        for run in runs:
            inner = run[-1]
            # Consecutive centres joined by segments clearing the innermost plane's probe (the smallest in the
            # run): a route from any plane of the run is then also a route from the innermost one (its half-space
            # contains the others'), so one search from it settles a run without any route.
            chained = all(
                tracer.state.segment_clearance(side[i].center, side[i - 1].center) > probes[inner] for i in run[:-1]
            )
            if chained and len(run) > 1 and route(inner, probes[inner]) is None:
                continue
            for k in sorted(run, key=lambda i: (-side[i].clearance, -i)):
                probe = max(probes[k], best[0]) if best is not None else probes[k]
                if cap(side[k].clearance) <= probe:
                    continue
                width = route(k, probe)
                if width is not None:
                    best = (width, k)
        branch = {}
        if best is not None:
            width, k = best
            branch[k] = {
                "plane_A": side[k].t,
                "next_plane_A": side[k + 1].t,
                "lateral_route_width_A": width,
                "continuation_width_A": continuation[k] if math.isfinite(continuation[k]) else None,
                "search_probe_A": probes[k],
                "plane_clearance_A": side[k].clearance,
            }
        if not branch:
            continue
        k = min(branch)
        new_ends[label] = side[k]
        branches[label] = {
            **branch[k],
            "margin_A": BRANCH_MARGIN_A,
            "margin_fraction": BRANCH_FRACTION,
            "inner_plane_A": planes[inner].t,
            "rule": LATERAL_BRANCH_RULE,
        }
        t, outer = side[k].t, side[k + 1].t
        mouths[index] = {
            "coordinate_A": t,
            "bracket_A": sorted((t, outer)),
            "enclosed_coordinate_A": t,
            "open_coordinate_A": outer,
            "axial_bracket_width_A": abs(outer - t),
            "reason": "lateral_branch",
            "traced_mouth": trace["mouths"][index],
        }
    if not new_ends:
        return trace
    lower = new_ends.get("lower", path[0])
    upper = new_ends.get("upper", path[-1])
    new_path = _resample(tracer, trace["seed"].center, lower, upper, samples)
    return {
        **trace,
        "path": new_path,
        "mouths": tuple(mouths),
        "branches": branches,
        "exit_cache": cache,
        "span": new_path[-1].t - new_path[0].t,
    }


def _build_profile(
    trace,
    origin,
    probe_radius,
    label,
    automatic,
    diagnostics,
    *,
    samples,
    search_radius,
    section_spacing,
    refinement_steps,
    ladder,
    lateral_exits,
    exit_bulk_radius,
    exit_spacing,
    refinement,
    fixed=False,
):
    from .sections import BRIDGE_WINDOW_FACTOR, ChannelResolutionError

    tracer, path = trace["tracer"], trace["path"]
    margins = [tracer.end_margin(path[0], -1), tracer.end_margin(path[-1], 1)]
    capped = [m <= lumen.probe for m, lumen in zip(margins, (path[0], path[-1]))]
    ends_in_bulk = []
    if any(capped) and not lateral_exits and math.isfinite(exit_bulk_radius) and exit_bulk_radius > 0:
        # A capped end whose sample already lies in bulk is open (BULK_END_DEFINITION).
        for k, (label, lumen, sign) in enumerate((("lower", path[0], -1), ("upper", path[-1], 1))):
            if capped[k] and _end_in_bulk(_end_envelope(tracer, lumen, sign, exit_bulk_radius), lumen):
                capped[k] = False
                ends_in_bulk.append(label)
    if any(capped) and not lateral_exits:
        ends = " and ".join(name for name, c in zip(("lower", "upper"), capped) if c)
        raise ChannelResolutionError(
            f"The channel's {ends} end is capped: the straight continuation along the axis beyond the mouth does not "
            "clear the end plane's enclosure probe. Add lateral_exits (--lateral-exits) to accept a capped end with a "
            "lateral exit to bulk"
        )
    points = []
    for i, lumen in enumerate(path):
        raw, nearest = tracer.state.spatial.nearest_surface(lumen.center)
        points.append(
            ChannelPoint(
                i, lumen.center, max(0.0, raw - probe_radius), raw, lumen.t, nearest.serial, nearest.residue_key.label
            )
        )
    segments = [tracer.state.segment_clearance(a.position, b.position) for a, b in zip(points, points[1:])]
    exit_metadata = {}
    if lateral_exits:
        exit_metadata = _exit_metadata(
            tracer,
            path,
            segments,
            probe_radius,
            bulk_radius=exit_bulk_radius,
            spacing=exit_spacing,
            branches=trace.get("branches"),
            cache=trace.get("exit_cache"),
        )
        exit_metadata["exit_policy"]["capped_mouth_adjustments"] = trace["adjustments"]
    arc = [0.0]
    for a, b in zip(points, points[1:]):
        arc.append(arc[-1] + distance(a.position, b.position))
    probes = [lumen.probe for lumen in path]
    lower_mouth, upper_mouth = trace["mouths"]
    selection = {
        "mode": "adaptive",
        "rule": ADAPTIVE_ENCLOSURE_RULE,
        "ladder_A": list(ladder),
        "probe_radius_A": probe_radius,
        "chosen_A": None,
        "probe_range_A": [min(probes), max(probes)],
        "reason": (f"per-plane enclosure probe, {min(probes):.2f}-{max(probes):.2f} A over {len(probes)} samples"),
    }
    middle = path[len(path) // 2]
    return PoreProfile(
        origin,
        tracer.direction,
        tuple(points),
        probe_radius,
        method="axial-connected",
        metadata={
            "axis": "auto:connected-sections" if automatic else label,
            "axis_candidates": diagnostics,
            "status": "resolved",
            "seed_position": middle.center,
            "seed_clearance_A": middle.clearance,
            "seed_policy": "interior of the selected connected channel",
            "volume_estimate_definition": "integral of inscribed sphere areas; not the full non-circular lumen volume",
            "section_spacing_A": section_spacing,
            "enclosure_mode": "fixed" if fixed else "adaptive",
            "enclosure_radius_A": max(probes),
            "enclosure_probes_A": probes,
            "enclosure_probe_range_A": [min(probes), max(probes)],
            **({} if fixed else {"enclosure_probe_selection": selection}),
            "section_basis": [tracer.u, tracer.v],
            "section_search_radius_A": search_radius,
            "section_areas_A2": [
                len(l.section.indices) * l.section_spacing**2 if l.section is not None else None for l in path
            ],
            "unenclosed_samples": [i for i, l in enumerate(path) if not l.enclosed],
            "unjoined_samples": [i for i, l in enumerate(path) if not l.joined],
            "sample_flags_definition": (
                "unenclosed_samples: resampled planes between the "
                "mouths where the lumen is not laterally enclosed at any probe (thin open slabs "
                "the 0.75 A discovery planes stepped over), measured at their climbed centre "
                "and given the previous plane's probe; unjoined_samples: samples whose segment "
                "from the previous one does not clear the smaller probe. The mouths come from "
                "the discovery trace and do not depend on the sample count."
            ),
            "section_spacings_A": [l.section_spacing for l in path],
            "section_search_radii_A": [l.search_radius for l in path],
            "window_bridged_samples": [i for i, l in enumerate(path) if l.window_bridged],
            "window_levels": {
                "levels": list(WINDOW_LEVELS),
                "per_clearance": WINDOW_PER_CLEARANCE,
                "rule": (
                    "each plane's search square is the smallest of search_radius "
                    "times 1, 2 or 4 whose half-width is at least "
                    f"{WINDOW_PER_CLEARANCE:g} times the lumen clearance"
                ),
            },
            "window_bridging": (
                "a plane whose lumen component reaches its search square is "
                f"enclosed if it is enclosed in a square {BRIDGE_WINDOW_FACTOR:g} "
                "times as wide; window_bridged_samples lists those samples"
            ),
            "requested_samples": samples,
            "samples": len(points),
            "search_radius": search_radius,
            "refinement_steps": refinement_steps,
            "section_refinement": refinement,
            "reaction_coordinate": "projection onto selected axis; arc lengths recorded separately",
            "arc_lengths_A": arc,
            "segment_clearances_A": segments,
            "continuous_bottleneck_radius_A": min(segments) - probe_radius,
            "end_clearances_A": margins,
            "channel_mouths": {
                "lower": lower_mouth,
                "upper": upper_mouth,
                "definition": (
                    "first plane, outwards from the seed, whose lumen is not laterally "
                    "enclosed at any probe below its clearance or cannot be joined"
                ),
                "axial_tolerance_A": MOUTH_TOLERANCE_A,
                "limitations": "brackets exclude transverse-grid uncertainty",
                **({"lateral_branch_rule": LATERAL_BRANCH_RULE} if trace.get("branches") else {}),
            },
            "mouth_adjustments": trace["adjustments"],
            "selection": "lumen centres climbed plane to plane, each enclosed at its own probe",
            "limitations": "finite grid; principal-direction proposals; monotone unbranched path; not a unique biological assignment",
            **(
                {
                    "ends_in_bulk": {
                        "ends": ends_in_bulk,
                        "bulk_radius_A": exit_bulk_radius,
                        "definition": BULK_END_DEFINITION,
                    }
                }
                if ends_in_bulk
                else {}
            ),
            **exit_metadata,
        },
    )

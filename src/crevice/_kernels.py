"""Compiled kernels for the hot loops of the channel profile (optional numba).

With numba installed (``pip install crevice[fast]``) the in-plane grid of a
channel section, the flood fill from the lumen centre and the lateral-exit
search of a capped end (:func:`exit_search`, :func:`distinct_exits`) run as
compiled loops; without it, the same results come from the NumPy code in
:mod:`crevice.sections`, :mod:`crevice.adaptive_channel` and
:mod:`crevice.channel_exits`. Both paths use the same arithmetic, operation
by operation, so they give the same free nodes, blocked edges, components and
exit legs (``tests/test_kernels.py`` checks this).

Set the environment variable ``CREVICE_DISABLE_NUMBA=1`` to force the NumPy
path. Compiled functions are cached on disk after their first use.
"""

from __future__ import annotations

import math
import os

import numpy as np

#: In-plane neighbour steps, in the order of ``crevice.sections._STEPS``.
STEPS = np.array(((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)), dtype=np.int64)

try:
    if os.environ.get("CREVICE_DISABLE_NUMBA", "").strip().lower() in {"1", "true", "yes", "on"}:
        raise ImportError("disabled by CREVICE_DISABLE_NUMBA")
    from numba import njit, objmode

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - exercised by the NumPy-path tests
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        """Stand-in decorator: returns the function unchanged (pure Python)."""
        if args and callable(args[0]):
            return args[0]
        return lambda function: function

    import contextlib

    def objmode(**types):
        """Stand-in for ``numba.objmode``: a no-op context (the code already runs in Python)."""
        return contextlib.nullcontext()


@njit(cache=True)
def plane_graph(px, py, pz, pr, t, h, ds, probe, steps):
    """Free nodes and blocked free-to-free edges of the section plane at ``t``.

    The same tests as ``crevice.sections.SectionGeometry._plane_graph``:
    a node ``(i, j)`` is covered when ``(i*ds - x)**2 + (j*ds - y)**2 <=
    (r + probe)**2 - (z - t)**2 + 1e-12`` for an atom ``(x, y, z, r)`` (in-plane
    coordinates), and an edge from a free node to a free neighbour is blocked
    when the closest point of the segment to the atom centre lies within that
    disk. Each atom is tested only on the index box that can reach it.

    Parameters
    ----------
    px, py, pz, pr : numpy.ndarray
        Atom coordinates in the section frame (in-plane ``u``, ``v``, axial
        ``z``) and radii, Å (``crevice.sections._FrameState`` columns).
    t : float
        Axial coordinate of the plane, Å.
    h : int
        Half-width of the square in nodes (``2*h + 1`` nodes per side).
    ds : float
        Grid spacing, Å.
    probe : float
        Enclosure probe, Å.
    steps : numpy.ndarray
        The eight neighbour steps, :data:`STEPS`.

    Returns
    -------
    free : numpy.ndarray
        ``(n, n)`` booleans, ``n = 2*h + 1``.
    blocked : numpy.ndarray
        ``(8, n, n)`` booleans; only edges between two free nodes are set.
    """
    n = 2 * h + 1
    inside = np.zeros((n, n), dtype=np.bool_)
    reach = ds * math.sqrt(2.0)
    count = px.shape[0]
    for a in range(count):
        dz = pz[a] - t
        if not abs(dz) <= pr[a] + probe:
            continue
        x, y = px[a], py[a]
        rr = (pr[a] + probe) ** 2 - dz * dz
        bound = math.sqrt(max(rr, 0.0) + 1e-12) + 1e-6
        thr = rr + 1e-12
        i0 = max(-h, int(math.floor((x - bound) / ds)) - 1)
        i1 = min(h, int(math.ceil((x + bound) / ds)) + 1)
        j0 = max(-h, int(math.floor((y - bound) / ds)) - 1)
        j1 = min(h, int(math.ceil((y + bound) / ds)) + 1)
        for i in range(i0, i1 + 1):
            gx = i * ds - x
            for j in range(j0, j1 + 1):
                gy = j * ds - y
                if gx * gx + gy * gy <= thr:
                    inside[i + h, j + h] = True
    free = ~inside
    blocked = np.zeros((8, n, n), dtype=np.bool_)
    for a in range(count):
        dz = pz[a] - t
        if not abs(dz) <= pr[a] + probe:
            continue
        x, y = px[a], py[a]
        rr = (pr[a] + probe) ** 2 - dz * dz
        bound = math.sqrt(max(rr, 0.0) + 1e-12) + 1e-6 + reach
        thr = rr + 1e-12
        i0 = max(-h, int(math.floor((x - bound) / ds)) - 1)
        i1 = min(h, int(math.ceil((x + bound) / ds)) + 1)
        j0 = max(-h, int(math.floor((y - bound) / ds)) - 1)
        j1 = min(h, int(math.ceil((y + bound) / ds)) + 1)
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                if not free[i + h, j + h]:
                    continue
                fi = float(i)
                fj = float(j)
                for k in range(8):
                    di = steps[k, 0]
                    dj = steps[k, 1]
                    ti = i + di
                    tj = j + dj
                    if abs(ti) > h or abs(tj) > h or not free[ti + h, tj + h]:
                        continue
                    f = ((x / ds - fi) * di + (y / ds - fj) * dj) / (di * di + dj * dj)
                    f = max(0.0, min(1.0, f))
                    ex = (fi + f * di) * ds - x
                    ey = (fj + f * dj) * ds - y
                    if ex * ex + ey * ey <= thr:
                        blocked[k, i + h, j + h] = True
    return free, blocked


@njit(cache=True)
def component_from(free, blocked, h, si, sj, steps):
    """Flood fill of the free component holding node ``(si, sj)`` (grid indices from ``-h`` to ``h``).

    An edge is followed when it is not blocked from either end (the graph of
    ``crevice.sections._symmetric_components``). The fill stops as soon as
    it reaches the boundary of the square.

    Parameters
    ----------
    free, blocked : numpy.ndarray
        Output of :func:`plane_graph`.
    h : int
        Half-width of the square in nodes.
    si, sj : int
        Start node.
    steps : numpy.ndarray
        The eight neighbour steps, :data:`STEPS`.

    Returns
    -------
    status : int
        0 the node is not free, 1 the component is enclosed, 2 it reaches the boundary.
    nodes : (m, 2) int array
        Grid indices of the component (enclosed only; empty otherwise).
    """
    n = 2 * h + 1
    empty = np.zeros((0, 2), dtype=np.int64)
    if not free[si + h, sj + h]:
        return 0, empty
    seen = np.zeros((n, n), dtype=np.bool_)
    stack_i = np.empty(n * n, dtype=np.int64)
    stack_j = np.empty(n * n, dtype=np.int64)
    out = np.empty((n * n, 2), dtype=np.int64)
    top = 0
    count = 0
    stack_i[0] = si
    stack_j[0] = sj
    top = 1
    seen[si + h, sj + h] = True
    while top > 0:
        top -= 1
        i = stack_i[top]
        j = stack_j[top]
        if abs(i) == h or abs(j) == h:
            return 2, empty
        out[count, 0] = i
        out[count, 1] = j
        count += 1
        for k in range(8):
            ti = i + steps[k, 0]
            tj = j + steps[k, 1]
            if abs(ti) > h or abs(tj) > h:
                continue
            if seen[ti + h, tj + h] or not free[ti + h, tj + h]:
                continue
            back = 7 - k  # _STEPS is symmetric: the reverse of step k is step 7 - k
            if blocked[k, i + h, j + h] or blocked[back, ti + h, tj + h]:
                continue
            seen[ti + h, tj + h] = True
            stack_i[top] = ti
            stack_j[top] = tj
            top += 1
    return 1, out[:count].copy()


# -- lateral-exit search (crevice.channel_exits) --------------------------------------------------------------
#
# The search of :func:`crevice.channel_exits.lateral_exit_legs` as one compiled loop. The node and segment
# clearances are computed here from the atom spheres instead of the SpatialIndex KD-tree, with the same formulas
# operation by operation (``crevice.spatial._point_surface_distances`` and ``_segment_surface_distances``), over
# a candidate set that always contains every atom that can attain the minimum; the minimum of the same values is
# the same number. ``tests/test_kernels.py`` checks every piece against the NumPy path.

#: Lattice steps of the exit search, in the order of ``crevice.channel_exits._STEPS``.
EXIT_STEPS = np.array(
    [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1) if (i, j, k) != (0, 0, 0)], dtype=np.int64
)
#: Edge (nodes) of the lattice blocks the compiled exit search allocates together.
EXIT_BLOCK = 8
#: Edge (Å) of the atom cells of the compiled clearance queries.
CELL_SIZE_A = 2.0
#: Extra reach (Å) of every candidate search: far above floating-point error, so no candidate is ever missed.
_REACH_MARGIN = 1e-6
#: A bulk-envelope index this close (in grid units) to a rounding boundary is not trusted to the compiled
#: arithmetic (NumPy's matrix product may round differently); when that rounding matters, the NumPy expression
#: decides (``depth_callback``).
_DEPTH_AMBIGUITY = 1e-6

#: ``exit_search`` status codes.
EXIT_OK, EXIT_NO_TERMINAL, EXIT_TRUNCATED, EXIT_FALLBACK = 0, 1, 2, 3


def atom_cells(xyz, radii, cell=CELL_SIZE_A):
    """Atoms sorted into cubic cells for the compiled clearance queries.

    Parameters
    ----------
    xyz : numpy.ndarray, shape (n, 3)
        Atom centres, Å.
    radii : numpy.ndarray, shape (n,)
        Atom radii, Å.
    cell : float
        Cell edge, Å.

    Returns
    -------
    tuple
        ``(xyz, radii, origin, cell, dims, starts, order)``: contiguous float64 centres and radii, the corner of
        cell ``(0, 0, 0)``, the cell edge, the number of cells per axis, and the atoms of flat cell ``c`` as
        ``order[starts[c]:starts[c + 1]]`` (flat index ``(a * dims[1] + b) * dims[2] + c``).
    """
    xyz = np.ascontiguousarray(xyz, dtype=np.float64).reshape(-1, 3)
    radii = np.ascontiguousarray(radii, dtype=np.float64)
    origin = xyz.min(axis=0) - 1e-6
    dims = np.maximum(1, np.floor((xyz.max(axis=0) - origin) / cell).astype(np.int64) + 1)
    index = np.minimum(np.floor((xyz - origin) / cell).astype(np.int64), dims - 1)
    flat = (index[:, 0] * dims[1] + index[:, 1]) * dims[2] + index[:, 2]
    order = np.argsort(flat, kind="stable").astype(np.int64)
    starts = np.zeros(int(dims.prod()) + 1, dtype=np.int64)
    np.add.at(starts, flat + 1, 1)
    return xyz, radii, origin.astype(np.float64), float(cell), dims, np.cumsum(starts), order


@njit(cache=True)
def round10(x):
    """``round(x, 10)`` exactly as Python computes it (correctly rounded, ties to even), for ``|x| < 1e5``.

    Python rounds the exact binary value of ``x`` to ten decimals and converts the decimal back to the nearest
    double. Here ``x * 1e10`` is split exactly into a double and its rounding error (Dekker's product), the
    integer ``n`` nearest that exact value is found (ties to even), and ``n / 1e10`` is the correctly rounded
    quotient (``n`` and ``1e10`` are exact doubles), which is the same double.
    """
    if not math.isfinite(x):
        return x
    a = abs(x)
    p = a * 1e10
    # Dekker's exact product: a * 1e10 == p + err exactly (no fused multiply-add is generated without fastmath).
    split = 134217729.0
    c = split * a
    ah = c - (c - a)
    al = a - ah
    bh = 10000000000.0
    bl = 0.0
    err = ((ah * bh - p) + ah * bl + al * bh) + al * bl
    fl = math.floor(p)
    frac = p - fl
    if frac > 0.5 or (frac == 0.5 and err > 0.0):
        n = fl + 1.0
    elif frac < 0.5 or (frac == 0.5 and err < 0.0):
        n = fl
    else:
        n = fl if fl == 2.0 * math.floor(fl / 2.0) else fl + 1.0
    r = n / 1e10
    return -r if x < 0 else r


@njit(cache=True)
def _cell_range(lo, hi, origin, cell, dims, axis):
    a = int(math.floor((lo - origin[axis]) / cell))
    b = int(math.floor((hi - origin[axis]) / cell))
    return max(a, 0), min(b, dims[axis] - 1)


@njit(cache=True)
def point_clearance(px, py, pz, reach, xyz, radii, origin, cell, dims, starts, order):
    """Smallest ``distance - radius`` over the atoms of every cell within ``reach`` of the point (inf if none).

    With ``reach`` at least the point's clearance plus the largest atom radius this is the exact clearance, the
    value of ``SpatialIndex.nearest_surface_values`` (same formula, same operations).

    Parameters
    ----------
    px, py, pz : float
        Query point, Å.
    reach : float
        Search radius, Å.
    xyz, radii, origin, cell, dims, starts, order
        The atom cells of :func:`atom_cells` (centres and radii in Å, cell corner, edge, cells per axis, and the
        atoms of each cell).
    """
    best = math.inf
    a0, a1 = _cell_range(px - reach, px + reach, origin, cell, dims, 0)
    b0, b1 = _cell_range(py - reach, py + reach, origin, cell, dims, 1)
    rr = reach * reach
    for a in range(a0, a1 + 1):
        x0 = origin[0] + a * cell
        gx = max(x0 - px, 0.0, px - (x0 + cell))
        for b in range(b0, b1 + 1):
            y0 = origin[1] + b * cell
            gy = max(y0 - py, 0.0, py - (y0 + cell))
            rest = rr - gx * gx - gy * gy
            if rest < 0.0:
                continue
            # The cells of this row within reach are contiguous in ``order``: one slice, no per-cell test.
            dz = math.sqrt(rest)
            c0, c1 = _cell_range(pz - dz, pz + dz, origin, cell, dims, 2)
            if c0 > c1:
                continue
            base = (a * dims[1] + b) * dims[2]
            for s in range(starts[base + c0], starts[base + c1 + 1]):
                n = order[s]
                dx = px - xyz[n, 0]
                dy = py - xyz[n, 1]
                dz2 = pz - xyz[n, 2]
                value = math.sqrt(dx * dx + dy * dy + dz2 * dz2) - radii[n]
                if value < best:
                    best = value
    return best


@njit(cache=True)
def segment_clearance_below(sx, sy, sz, ex, ey, ez, limit, probe, rmax, xyz, radii, origin, cell, dims, starts, order):
    """``min(limit, clearance of the closed segment)``, exact, or a value ``<= probe`` as soon as one is found.

    The segment clearance is the value of ``SpatialIndex.segment_clearances`` (the operations of
    ``crevice.spatial._segment_surface_distances``). Only atoms within ``L/2 + limit + rmax`` of the midpoint can
    give a value below ``limit``; atoms farther away give values above it and do not change the result.

    Parameters
    ----------
    sx, sy, sz, ex, ey, ez : float
        Segment start and end, Å (the order of the NumPy call: start = the expanded node).
    limit : float
        Values at or above ``limit`` are not needed; it is returned instead.
    probe : float
        The search stops at the first value ``<= probe``.
    rmax : float
        Largest atom radius, Å.
    xyz, radii, origin, cell, dims, starts, order
        The atom cells of :func:`atom_cells` (centres and radii in Å, cell corner, edge, cells per axis, and the
        atoms of each cell).
    """
    d0 = ex - sx
    d1 = ey - sy
    d2 = ez - sz
    length_squared = d0 * d0 + d1 * d1 + d2 * d2
    zero = length_squared == 0
    safe = 1.0 if zero else length_squared
    mx = sx + 0.5 * d0
    my = sy + 0.5 * d1
    mz = sz + 0.5 * d2
    reach = 0.5 * math.sqrt(length_squared) + limit + rmax + _REACH_MARGIN
    best = limit
    a0, a1 = _cell_range(mx - reach, mx + reach, origin, cell, dims, 0)
    b0, b1 = _cell_range(my - reach, my + reach, origin, cell, dims, 1)
    rr = reach * reach
    for a in range(a0, a1 + 1):
        x0 = origin[0] + a * cell
        gx = max(x0 - mx, 0.0, mx - (x0 + cell))
        for b in range(b0, b1 + 1):
            y0 = origin[1] + b * cell
            gy = max(y0 - my, 0.0, my - (y0 + cell))
            rest = rr - gx * gx - gy * gy
            if rest < 0.0:
                continue
            dz = math.sqrt(rest)
            c0, c1 = _cell_range(mz - dz, mz + dz, origin, cell, dims, 2)
            if c0 > c1:
                continue
            base = (a * dims[1] + b) * dims[2]
            for s in range(starts[base + c0], starts[base + c1 + 1]):
                n = order[s]
                w0 = xyz[n, 0] - sx
                w1 = xyz[n, 1] - sy
                w2 = xyz[n, 2] - sz
                fraction = max(0.0, min(1.0, (w0 * d0 + w1 * d1 + w2 * d2) / safe))
                if zero:
                    fraction = 0.0
                e0 = w0 - d0 * fraction
                e1 = w1 - d1 * fraction
                e2 = w2 - d2 * fraction
                value = math.sqrt(e0 * e0 + e1 * e1 + e2 * e2) - radii[n]
                if value < best:
                    best = value
                    if best <= probe:
                        return best
    return best
    return best


@njit(cache=True)
def _grow_float(a, n):
    out = np.empty(n, dtype=np.float64)
    out[: a.shape[0]] = a
    return out


@njit(cache=True)
def _grow_int(a, n):
    out = np.empty(n, dtype=np.int64)
    out[: a.shape[0]] = a
    return out


@njit(cache=True)
def _grow_byte(a, n):
    out = np.zeros(n, dtype=np.int8)
    out[: a.shape[0]] = a
    return out


@njit(cache=True)
def _heap_less(k1, k2, k3, a, b):
    if k1[a] != k1[b]:
        return k1[a] < k1[b]
    if k2[a] != k2[b]:
        return k2[a] < k2[b]
    return k3[a] < k3[b]


@njit(cache=True)
def _heap_push(k1, k2, k3, kn, size, v1, v2, v3, vn):
    i = size
    k1[i], k2[i], k3[i], kn[i] = v1, v2, v3, vn
    while i > 0:
        parent = (i - 1) // 2
        if _heap_less(k1, k2, k3, i, parent):
            k1[i], k1[parent] = k1[parent], k1[i]
            k2[i], k2[parent] = k2[parent], k2[i]
            k3[i], k3[parent] = k3[parent], k3[i]
            kn[i], kn[parent] = kn[parent], kn[i]
            i = parent
        else:
            break
    return size + 1


@njit(cache=True)
def _heap_pop(k1, k2, k3, kn, size):
    """Move the smallest entry to index ``size - 1``; return the new size."""
    last = size - 1
    k1[0], k1[last] = k1[last], k1[0]
    k2[0], k2[last] = k2[last], k2[0]
    k3[0], k3[last] = k3[last], k3[0]
    kn[0], kn[last] = kn[last], kn[0]
    i = 0
    while True:
        left = 2 * i + 1
        if left >= last:
            break
        child = left
        if left + 1 < last and _heap_less(k1, k2, k3, left + 1, left):
            child = left + 1
        if _heap_less(k1, k2, k3, child, i):
            k1[i], k1[child] = k1[child], k1[i]
            k2[i], k2[child] = k2[child], k2[i]
            k3[i], k3[child] = k3[child], k3[i]
            kn[i], kn[child] = kn[child], kn[i]
            i = child
        else:
            break
    return last


@njit(cache=True)
def _position(origin, u, v, w, spacing, i, j, k, out):
    # crevice.channel_exits position(): add(origin, add(scale(u, i*s), add(scale(v, j*s), scale(w, k*s))))
    fi = i * spacing
    fj = j * spacing
    fk = k * spacing
    for c in range(3):
        out[c] = origin[c] + (u[c] * fi + (v[c] * fj + w[c] * fk))


@njit(cache=True)
def _bulk_state(p, env_origin, env_basis, env_low, env_spacing, env_dims, env_distance, bulk_radius):
    """1 when ``BulkEnvelope.depth(p) == 0.0``, 0 when not, 2 when only the NumPy rounding can tell.

    An index whose value is within ``_DEPTH_AMBIGUITY`` of a rounding boundary may round either way in the NumPy
    path; every combination of the two candidates of such indices is tried, and the answer is certain when all
    combinations agree.
    """
    f = np.empty(3, dtype=np.float64)
    ambiguous = np.zeros(3, dtype=np.bool_)
    for a in range(3):
        loc = (
            (p[0] - env_origin[0]) * env_basis[a, 0]
            + (p[1] - env_origin[1]) * env_basis[a, 1]
            + (p[2] - env_origin[2]) * env_basis[a, 2]
        )
        f[a] = (loc - env_low[a]) / env_spacing
        r = f[a] - math.floor(f[a])
        ambiguous[a] = abs(r - 0.5) < _DEPTH_AMBIGUITY
    result = -1
    index = np.empty(3, dtype=np.int64)
    for combo in range(8):
        skip = False
        for a in range(3):
            up = (combo >> a) & 1
            if up and not ambiguous[a]:
                skip = True
            if ambiguous[a]:
                index[a] = int(math.floor(f[a])) + up
            else:
                index[a] = int(math.floor(f[a] + 0.5))
        if skip:
            continue
        bulk = 0
        for a in range(3):
            if index[a] < 0 or index[a] >= env_dims[a]:
                bulk = 1
        if bulk == 0:
            depth = env_distance[(index[0] * env_dims[1] + index[1]) * env_dims[2] + index[2]] - bulk_radius
            bulk = 1 if max(0.0, depth) == 0.0 else 0
        if result == -1:
            result = bulk
        elif result != bulk:
            return 2
    return result


def depth_callback(i, j, k):
    """Exact ``BulkEnvelope.depth(...) == 0.0`` of lattice node ``(i, j, k)`` (set by ``channel_exits``)."""
    return 1 if _DEPTH_NODE[0](i, j, k) else 0


#: The NumPy-path bulk test of the running exit search (``channel_exits`` sets it around each call).
_DEPTH_NODE = [None]


@njit(cache=True)
def exit_search(
    xyz,
    radii,
    rmax,
    cell_origin,
    cell,
    cell_dims,
    starts,
    order,
    origin,
    u,
    v,
    w,
    spacing,
    probe,
    root_clearance,
    env_origin,
    env_basis,
    env_low,
    env_spacing,
    env_dims,
    env_distance,
    bulk_radius,
    box_low,
    box_dims,
    max_nodes,
    first_terminal_only,
):
    """Widest-path (maximin) search of ``crevice.channel_exits.lateral_exit_legs``, compiled.

    The same lexicographic Dijkstra as the NumPy path: entries ``(-round(width, 10), route length, counter)``
    popped in order, stale entries skipped, bulk nodes terminal, the 26 lattice steps in the same order, the
    same mouth-sphere rule, lower bound and tie-breaks. Lattice nodes ``(i, j, k)`` (``k >= 0``) are stored in
    blocks of ``EXIT_BLOCK**3`` inside the index box ``box_low + [0, box_dims)``.

    Parameters
    ----------
    xyz, radii, rmax, cell_origin, cell, cell_dims, starts, order
        Atom centres and radii (Å), the largest radius and the atom cells of :func:`atom_cells`.
    origin, u, v, w : numpy.ndarray
        Start of the search (the mouth sample's centre) and the lattice directions (``w`` points away from
        the channel).
    spacing : float
        Lattice spacing, Å.
    probe : float
        Enclosure probe every node and segment must exceed, Å.
    root_clearance : float
        Clearance of the start (the mouth sphere radius), Å, as the NumPy path computes it.
    env_origin, env_basis, env_low, env_spacing, env_dims, env_distance, bulk_radius
        The :class:`crevice.channel_exits.BulkEnvelope` (frame origin and rows, grid corner, spacing, shape,
        flattened distance to the nearest bulk centre, bulk-probe radius).
    box_low, box_dims : numpy.ndarray
        Lattice index box of the search.
    max_nodes : int
        Lattice nodes the search may label.
    first_terminal_only : bool
        Stop at the first bulk terminal (``max_legs == 1``).

    Returns
    -------
    status : int
        ``EXIT_OK``; ``EXIT_NO_TERMINAL`` (``max_nodes`` exceeded before any terminal); ``EXIT_TRUNCATED``
        (exceeded after one); ``EXIT_FALLBACK`` (a node outside the box or a width beyond ``round10``'s range:
        the caller reruns the NumPy path).
    labelled : int
        Lattice nodes labelled (``len(width)`` of the NumPy path).
    terminals : numpy.ndarray
        Node ids of the bulk terminals in the order they were popped.
    parent, width : numpy.ndarray
        Parent node id (-1 at the start) and width of every node id, int64 and float64.
    block_ijk : numpy.ndarray
        Lattice index (int64, shape (m, 3)) of the first node of each allocated block (node id ``slot * EXIT_BLOCK**3 + local``).
    """
    eb = EXIT_BLOCK
    bb = eb * eb * eb
    nbi = (box_dims[0] + eb - 1) // eb
    nbj = (box_dims[1] + eb - 1) // eb
    nbk = (box_dims[2] + eb - 1) // eb
    block_index = np.full(nbi * nbj * nbk, -1, dtype=np.int64)
    cap = 64
    block_ijk = np.empty((cap, 3), dtype=np.int64)
    clear = np.empty(cap * bb, dtype=np.float64)
    known = np.zeros(cap * bb, dtype=np.int8)  # 1 clearance computed, 2 also labelled
    width = np.empty(cap * bb, dtype=np.float64)
    wround = np.empty(cap * bb, dtype=np.float64)
    travelled = np.empty(cap * bb, dtype=np.float64)
    parent = np.empty(cap * bb, dtype=np.int64)
    blocks = 0
    hcap = 4096
    h1 = np.empty(hcap, dtype=np.float64)
    h2 = np.empty(hcap, dtype=np.float64)
    h3 = np.empty(hcap, dtype=np.int64)
    hn = np.empty(hcap, dtype=np.int64)
    hsize = 0
    terminals = np.empty(1024, dtype=np.int64)
    nterm = 0
    steps = EXIT_STEPS
    lengths = np.empty(26, dtype=np.float64)
    for s in range(26):
        lengths[s] = spacing * math.sqrt(steps[s, 0] ** 2 + steps[s, 1] ** 2 + steps[s, 2] ** 2)
    p = np.empty(3, dtype=np.float64)
    q = np.empty(3, dtype=np.float64)
    status = EXIT_OK

    # The root (0, 0, 0).
    ri, rj, rk = -box_low[0], -box_low[1], -box_low[2]
    if ri < 0 or rj < 0 or rk < 0 or ri >= box_dims[0] or rj >= box_dims[1] or rk >= box_dims[2]:
        return EXIT_FALLBACK, 0, terminals[:0], parent[:0], width[:0], block_ijk[:0]
    block_index[((ri // eb) * nbj + rj // eb) * nbk + rk // eb] = 0
    block_ijk[0, 0] = (ri // eb) * eb + box_low[0]
    block_ijk[0, 1] = (rj // eb) * eb + box_low[1]
    block_ijk[0, 2] = (rk // eb) * eb + box_low[2]
    blocks = 1
    root = ((ri % eb) * eb + rj % eb) * eb + rk % eb
    clear[root] = root_clearance
    known[root] = 2
    width[root] = math.inf
    wround[root] = math.inf
    travelled[root] = 0.0
    parent[root] = -1
    labelled = 1
    mouth_sphere = root_clearance
    counter = 0
    hsize = _heap_push(h1, h2, h3, hn, hsize, -math.inf, 0.0, 0, root)
    while hsize > 0:
        hsize = _heap_pop(h1, h2, h3, hn, hsize)
        key, walked, node = h1[hsize], h2[hsize], hn[hsize]
        current = width[node]
        if -key < wround[node] or (-key == wround[node] and walked > travelled[node] + 1e-9):
            continue
        slot = node // bb
        local = node % bb
        ni = block_ijk[slot, 0] + local // (eb * eb)
        nj = block_ijk[slot, 1] + (local // eb) % eb
        nk = block_ijk[slot, 2] + local % eb
        _position(origin, u, v, w, spacing, ni, nj, nk, p)
        # BulkEnvelope.depth(p) == 0.0. The envelope index is rint of a NumPy matrix product, whose last bits
        # depend on the BLAS; where an index is at a rounding boundary and the candidate cells disagree, the
        # NumPy path itself decides (``depth_callback``, the exact expression of BulkEnvelope.depth).
        state = _bulk_state(p, env_origin, env_basis, env_low, env_spacing, env_dims, env_distance, bulk_radius)
        if state == 2:
            with objmode(state="int64"):
                state = depth_callback(ni, nj, nk)
        bulk = state == 1
        if bulk:
            if nterm == terminals.shape[0]:
                terminals = _grow_int(terminals, 2 * nterm)
            terminals[nterm] = node
            nterm += 1
            if first_terminal_only:
                break
            continue
        if labelled > max_nodes:
            status = EXIT_NO_TERMINAL if nterm == 0 else EXIT_TRUNCATED
            break
        cn = clear[node]
        tn = travelled[node]
        for s in range(26):
            oi = ni + steps[s, 0]
            oj = nj + steps[s, 1]
            ok = nk + steps[s, 2]
            if ok < 0:
                continue
            bi, bj, bk = oi - box_low[0], oj - box_low[1], ok - box_low[2]
            if bi < 0 or bj < 0 or bk < 0 or bi >= box_dims[0] or bj >= box_dims[1] or bk >= box_dims[2]:
                return EXIT_FALLBACK, labelled, terminals[:0], parent[:0], width[:0], block_ijk[:0]
            flat = ((bi // eb) * nbj + bj // eb) * nbk + bk // eb
            oslot = block_index[flat]
            if oslot < 0:
                if blocks == cap:
                    cap *= 2
                    grown = np.empty((cap, 3), dtype=np.int64)
                    grown[:blocks] = block_ijk[:blocks]
                    block_ijk = grown
                    clear = _grow_float(clear, cap * bb)
                    known = _grow_byte(known, cap * bb)
                    width = _grow_float(width, cap * bb)
                    wround = _grow_float(wround, cap * bb)
                    travelled = _grow_float(travelled, cap * bb)
                    parent = _grow_int(parent, cap * bb)
                oslot = blocks
                block_index[flat] = oslot
                block_ijk[oslot, 0] = (bi // eb) * eb + box_low[0]
                block_ijk[oslot, 1] = (bj // eb) * eb + box_low[1]
                block_ijk[oslot, 2] = (bk // eb) * eb + box_low[2]
                blocks += 1
            other = oslot * bb + ((bi % eb) * eb + bj % eb) * eb + bk % eb
            length = lengths[s]
            _position(origin, u, v, w, spacing, oi, oj, ok, q)
            if known[other] == 0:
                # Clearance is 1-Lipschitz: the neighbour's is at most cn + length, so its nearest atom lies
                # within cn + length + rmax.
                clear[other] = point_clearance(
                    q[0],
                    q[1],
                    q[2],
                    cn + length + rmax + _REACH_MARGIN,
                    xyz,
                    radii,
                    cell_origin,
                    cell,
                    cell_dims,
                    starts,
                    order,
                )
                known[other] = 1
            c = clear[other]
            if c <= probe:
                continue
            dx = q[0] - origin[0]
            dy = q[1] - origin[1]
            dz = q[2] - origin[2]
            inside = math.sqrt(dx * dx + dy * dy + dz * dz) <= mouth_sphere
            new = current if inside else min(current, c)
            bound = (cn + c - length) / 2
            if inside:
                if bound <= probe:
                    # Only whether the segment clears the probe matters: any limit above the probe gives that.
                    segment = segment_clearance_below(
                        p[0],
                        p[1],
                        p[2],
                        q[0],
                        q[1],
                        q[2],
                        probe + 1e-6,
                        probe,
                        rmax,
                        xyz,
                        radii,
                        cell_origin,
                        cell,
                        cell_dims,
                        starts,
                        order,
                    )
                    if segment <= probe:
                        continue
            elif bound < new:
                segment = segment_clearance_below(
                    p[0],
                    p[1],
                    p[2],
                    q[0],
                    q[1],
                    q[2],
                    new,
                    probe,
                    rmax,
                    xyz,
                    radii,
                    cell_origin,
                    cell,
                    cell_dims,
                    starts,
                    order,
                )
                if segment <= probe:
                    continue
                new = segment
            if abs(new) >= 1e5 and math.isfinite(new):
                return EXIT_FALLBACK, labelled, terminals[:0], parent[:0], width[:0], block_ijk[:0]
            rounded = round10(new)
            total = tn + length
            known_width = wround[other] if known[other] == 2 else -math.inf
            if rounded > known_width or (rounded == known_width and total < travelled[other] - 1e-9):
                if known[other] != 2:
                    known[other] = 2
                    labelled += 1
                width[other] = new
                wround[other] = rounded
                travelled[other] = total
                parent[other] = node
                counter += 1
                if hsize == h1.shape[0]:
                    h1 = _grow_float(h1, 2 * hsize)
                    h2 = _grow_float(h2, 2 * hsize)
                    h3 = _grow_int(h3, 2 * hsize)
                    hn = _grow_int(hn, 2 * hsize)
                hsize = _heap_push(h1, h2, h3, hn, hsize, -rounded, total, counter, other)
    return status, labelled, terminals[:nterm], parent, width, block_ijk[:blocks]


@njit(cache=True)
def _node_ijk(node, block_ijk):
    eb = EXIT_BLOCK
    bb = eb * eb * eb
    slot = node // bb
    local = node % bb
    return (
        block_ijk[slot, 0] + local // (eb * eb),
        block_ijk[slot, 1] + (local // eb) % eb,
        block_ijk[slot, 2] + local % eb,
    )


@njit(cache=True)
def distinct_exits(terminals, parent, block_ijk, origin, u, v, w, spacing, separation, route_separation, max_legs):
    """The distinct-exit selection of ``crevice.channel_exits.lateral_exit_legs`` over the search terminals.

    A terminal is kept when its exit point is more than ``separation`` from every kept exit point and the parts
    of its route and of every kept route farther than ``route_separation`` from the start stay more than
    ``route_separation`` apart (the route starts at ``origin``; same distances as the NumPy path).

    Parameters
    ----------
    terminals, parent, block_ijk
        Output of :func:`exit_search`.
    origin, u, v, w, spacing
        The search lattice (as for :func:`exit_search`).
    separation, route_separation : float
        Distinct-exit distances, Å.
    max_legs : int
        Largest number of kept exits.

    Returns
    -------
    kept : numpy.ndarray
        Node ids of the kept terminals, in order.
    truncated : bool
        A further distinct exit was found after ``max_legs``.
    """
    kept = np.empty(max_legs, dtype=np.int64)
    exits = np.empty((max_legs, 3), dtype=np.float64)
    outer_start = np.zeros(max_legs + 1, dtype=np.int64)
    outer = np.empty((256, 3), dtype=np.float64)
    nkept = 0
    truncated = False
    p = np.empty(3, dtype=np.float64)
    route = np.empty((256, 3), dtype=np.float64)
    for t in range(terminals.shape[0]):
        node = terminals[t]
        i, j, k = _node_ijk(node, block_ijk)
        _position(origin, u, v, w, spacing, i, j, k, p)
        near = False
        for a in range(nkept):
            dx = p[0] - exits[a, 0]
            dy = p[1] - exits[a, 1]
            dz = p[2] - exits[a, 2]
            if math.sqrt(dx * dx + dy * dy + dz * dz) <= separation:
                near = True
                break
        if near:
            continue
        # The route, start (``origin``) first.
        n = 0
        walk = node
        while walk >= 0:
            if n == route.shape[0]:
                grown = np.empty((2 * n, 3), dtype=np.float64)
                grown[:n] = route[:n]
                route = grown
            i, j, k = _node_ijk(walk, block_ijk)
            _position(origin, u, v, w, spacing, i, j, k, route[n])
            n += 1
            walk = parent[walk]
        route[n - 1, 0] = origin[0]
        route[n - 1, 1] = origin[1]
        route[n - 1, 2] = origin[2]
        base = outer_start[nkept]
        m = base
        for r in range(n - 1, -1, -1):
            dx = route[r, 0] - origin[0]
            dy = route[r, 1] - origin[1]
            dz = route[r, 2] - origin[2]
            if math.sqrt(dx * dx + dy * dy + dz * dz) > route_separation:
                if m == outer.shape[0]:
                    grown = np.empty((2 * m, 3), dtype=np.float64)
                    grown[:m] = outer[:m]
                    outer = grown
                outer[m] = route[r]
                m += 1
        if m == base:
            continue
        close = False
        for a in range(nkept):
            for x in range(base, m):
                for y in range(outer_start[a], outer_start[a + 1]):
                    dx = outer[x, 0] - outer[y, 0]
                    dy = outer[x, 1] - outer[y, 1]
                    dz = outer[x, 2] - outer[y, 2]
                    if math.sqrt(dx * dx + dy * dy + dz * dz) <= route_separation:
                        close = True
                        break
                if close:
                    break
            if close:
                break
        if close:
            continue
        if nkept >= max_legs:
            truncated = True
            break
        kept[nkept] = node
        exits[nkept] = p
        nkept += 1
        outer_start[nkept] = m
    return kept[:nkept], truncated


@njit(cache=True)
def envelope_classes(
    low, spacing, dims, low_threshold, high_threshold, xyz, radii, origin, cell, cell_dims, starts, order
):
    """Which of two thresholds the clearance of each bulk-envelope grid point reaches (0, 1 or 2).

    Grid point ``(i, j, k)`` is ``low + (i, j, k) * spacing`` (the expression of ``BulkEnvelope``). Class 2:
    clearance ``>= high_threshold``; 1: ``>= low_threshold`` only; 0: neither. These two comparisons are all
    the bulk envelope uses its clearances for. Along each grid row the clearance changes by at most
    ``spacing`` per step (it is 1-Lipschitz), so a point whose class follows from the previous bound (with a
    1e-6 Å margin) is not evaluated; every other point is evaluated exactly (the formula of
    ``SpatialIndex.nearest_surface_values``) up to a cap of ``high_threshold + 4 * spacing``, above which only
    the class matters.

    Parameters
    ----------
    low : numpy.ndarray
        Grid corner, Å (envelope frame).
    spacing : float
        Grid spacing, Å.
    dims : numpy.ndarray
        Grid shape.
    low_threshold, high_threshold : float
        The two clearance thresholds, Å.
    xyz, radii, origin, cell, cell_dims, starts, order
        Atom cells of :func:`atom_cells` built from the atoms in the envelope frame.

    Returns
    -------
    numpy.ndarray of int8, shape ``dims``
    """
    out = np.empty((dims[0], dims[1], dims[2]), dtype=np.int8)
    rmax = 0.0
    for n in range(radii.shape[0]):
        rmax = max(rmax, radii[n])
    cap = high_threshold + 4.0 * spacing
    reach = cap + rmax + _REACH_MARGIN
    step = spacing + _REACH_MARGIN
    for i in range(dims[0]):
        px = low[0] + i * spacing
        for j in range(dims[1]):
            py = low[1] + j * spacing
            known = False
            lo = 0.0
            hi = 0.0
            for k in range(dims[2]):
                if known:
                    lo -= step
                    hi += step
                    if lo >= high_threshold:
                        out[i, j, k] = 2
                        continue
                    if hi < low_threshold:
                        out[i, j, k] = 0
                        continue
                    if lo >= low_threshold and hi < high_threshold:
                        out[i, j, k] = 1
                        continue
                pz = low[2] + k * spacing
                value = point_clearance(px, py, pz, reach, xyz, radii, origin, cell, cell_dims, starts, order)
                if value < cap:
                    lo = value
                    hi = value
                else:
                    lo = cap - _REACH_MARGIN
                    hi = math.inf
                known = True
                out[i, j, k] = 2 if value >= high_threshold else (1 if value >= low_threshold else 0)
    return out

"""Nearest-surface and segment-clearance queries against atom spheres.

Every clearance-based analysis in CREVICE asks the same question: how far is
this point (or this segment) from the nearest atom *surface*? Atoms are
modelled as hard spheres whose radii come from
:func:`crevice.radii.atom_vdw_radius` (element-based Bondi-type radii).

Definitions
    *Surface distance* of a point ``p`` to atom ``a`` is
    ``|p - a.coord| - r_a`` in Å. It is negative inside the sphere.
    *Clearance* of ``p`` is the minimum surface distance over all atoms: the
    radius of the largest empty sphere centred at ``p``.

Main entry point: :class:`SpatialIndex`.

Notes
-----
Queries are exact for the sphere model, with or without SciPy. The sphere
model itself is an approximation: it ignores bonding overlap, atom
flexibility and hydrogen atoms that were filtered out before indexing.
"""

from __future__ import annotations

from itertools import chain
import math

import numpy as np

from .geometry import distance, point_segment_distance
from .models import Atom, Coord
from .radii import atom_vdw_radius


class SpatialIndex:
    """Exact atom-surface distance queries for a fixed set of atoms.

    When SciPy is available the atom centres are stored in a
    :class:`scipy.spatial.cKDTree`; otherwise every query is a deterministic
    brute-force loop. Both paths return the same values.

    The KD-tree finds the nearest atom *centre* first. Because atom radii
    differ, a farther centre can have a nearer surface, so the search is then
    widened by the largest atomic radius, which bounds how much any other atom
    can improve on the candidate.

    Parameters
    ----------
    atoms : tuple of Atom
        Atoms to index. Radii are taken from each atom's element.

    Raises
    ------
    ValueError
        If ``atoms`` is empty.

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.spatial import SpatialIndex
    >>> index = SpatialIndex((Atom(1, "C", "UNK", "A", 1, 0.0, 0.0, 0.0, "C"),))
    >>> clearance, atom = index.nearest_surface((3.0, 0.0, 0.0))
    >>> round(clearance, 2), atom.serial      # 3.0 Å minus the 1.70 Å carbon radius
    (1.3, 1)
    >>> round(index.segment_clearance((-2.0, 2.0, 0.0), (2.0, 2.0, 0.0)), 2)
    0.3
    """

    def __init__(self, atoms: tuple[Atom, ...]):
        if not atoms:
            raise ValueError("SpatialIndex requires at least one atom")
        self.atoms = atoms
        self._radii = {atom: atom_vdw_radius(atom) for atom in atoms}
        self._max_radius = max(self._radii.values())
        self._tree = None
        self._xyz = self._radius_array = None
        try:
            from scipy.spatial import cKDTree  # type: ignore[import-not-found]
        except ImportError:
            self._coords = None
        else:
            import numpy as np
            self._coords = [atom.coord for atom in atoms]
            self._tree = cKDTree(self._coords)
            self._xyz = np.array(self._coords, dtype=float).reshape(-1, 3)
            self._radius_array = np.array([self._radii[atom] for atom in atoms], dtype=float)

    def nearest_surface(self, point: Coord) -> tuple[float, Atom]:
        """Clearance of a point and the atom that limits it.

        Parameters
        ----------
        point : Coord
            Query position in Å.

        Returns
        -------
        clearance : float
            Minimum signed surface distance in Å (negative inside an atom).
        atom : Atom
            Atom whose surface is nearest. On an exact tie the first atom found is
            kept; which one that is may differ between the KD-tree and brute-force
            paths.
        """
        if self._tree is not None:
            return self._nearest_surface_kdtree(point)
        return self._nearest_surface_bruteforce(point)

    def atom_surface_distance(self, point: Coord, atom: Atom) -> float:
        """Signed distance from a point to one atom's van der Waals surface.

        Parameters
        ----------
        point : Coord
            Position in Å.
        atom : Atom
            Must be one of the indexed atoms.

        Returns
        -------
        float
            ``|point - atom| - r_atom`` in Å (negative inside the atom).
        """
        return distance(point, atom.coord) - self._radii[atom]

    def segment_clearance(self, start: Coord, end: Coord) -> float:
        """Minimum signed atom-surface clearance over an entire closed segment.

        The exact minimum is found by computing, for each candidate atom, the
        distance from its centre to the segment (see
        :func:`crevice.geometry.point_segment_distance`) minus its radius. This is
        analytic, not sampling along the segment. With a KD-tree, candidates are
        limited to atoms whose centres lie within ``L/2 + best + r_max`` of the
        segment midpoint (``L`` the segment length, ``best`` the clearance found for
        the atom nearest the midpoint, ``r_max`` the largest atomic radius). That
        ball contains every atom that could improve the result.

        Parameters
        ----------
        start, end : Coord
            Segment endpoints in Å. Identical endpoints reduce to
            :meth:`nearest_surface`.

        Returns
        -------
        float
            Segment clearance in Å. No probe radius is subtracted; compare the
            result with the probe radius you require.

        Raises
        ------
        ValueError
            If an endpoint is not a finite 3D coordinate.
        """
        if any(len(p) != 3 or not all(math.isfinite(v) for v in p) for p in (start, end)):
            raise ValueError("Segment endpoints must be finite 3D coordinates")
        if all(a == b for a, b in zip(start, end)):
            return self.nearest_surface(start)[0]
        if self._tree is None:
            return min(point_segment_distance(a.coord, start, end) - self._radii[a] for a in self.atoms)
        midpoint = tuple(a + 0.5 * (b - a) for a, b in zip(start, end))
        _, nearest = self.nearest_surface(midpoint)
        best = point_segment_distance(nearest.coord, start, end) - self._radii[nearest]
        # Any improving atom center lies within half the segment length plus
        # best clearance plus the maximum atom radius of the midpoint.
        radius = 0.5 * distance(start, end) + best + self._max_radius
        radius += max(1e-9, abs(radius) * 1e-12)
        indices = self._tree.query_ball_point(midpoint, radius, eps=0.0)
        if indices:
            # point_segment_distance for every candidate at once, with the same
            # operations in the same order (bit-identical values).
            values = _segment_surface_distances(self._xyz[indices], self._radius_array[indices], start, end)
            best = min(best, float(values.min()))
        return best

    def segment_clearances(self, starts, ends):
        """Clearances of many closed segments at once (values of :meth:`segment_clearance`).

        Bit-identical to ``[segment_clearance(a, b) for a, b in zip(starts,
        ends)]``. Each segment is tested against every atom centre within
        ``L/2 + c_mid + r_max`` (plus a small margin) of its midpoint, where
        ``c_mid`` is the midpoint clearance; this ball contains the one
        :meth:`segment_clearance` searches, and atoms outside that one cannot
        lower the minimum, so the minimum is the same. Distances use the
        operations of :func:`crevice.geometry.point_segment_distance`.

        Parameters
        ----------
        starts, ends : array_like, shape (m, 3)
            Segment endpoints in Å.

        Returns
        -------
        numpy.ndarray, shape (m,)
            Segment clearances in Å.

        Raises
        ------
        ValueError
            If an endpoint is not a finite 3D coordinate.
        """
        a = np.asarray(starts, dtype=float).reshape(-1, 3)
        b = np.asarray(ends, dtype=float).reshape(-1, 3)
        if a.shape != b.shape:
            raise ValueError("starts and ends must have the same shape")
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            raise ValueError("Segment endpoints must be finite 3D coordinates")
        out = np.empty(len(a))
        if len(a) == 0:
            return out
        if self._tree is None:
            return np.array([self.segment_clearance(tuple(map(float, s)), tuple(map(float, e))) for s, e in zip(a, b)])
        same = (a == b).all(axis=1)
        if same.any():
            out[same] = self.nearest_surface_values(a[same])
        rows = np.nonzero(~same)[0]
        if len(rows) == 0:
            return out
        s, e = a[rows], b[rows]
        midpoint = s + 0.5 * (e - s)
        g = s - e
        length = np.sqrt(g[:, 0] * g[:, 0] + g[:, 1] * g[:, 1] + g[:, 2] * g[:, 2])
        radius = 0.5 * length + self.nearest_surface_values(midpoint) + self._max_radius
        radius = np.maximum(radius, 0.0) + 1e-6 + np.abs(radius) * 1e-9
        balls = self._tree.query_ball_point(midpoint, radius, return_sorted=False)
        counts = np.fromiter(map(len, balls), dtype=np.intp, count=len(balls))
        owner = np.repeat(np.arange(len(rows)), counts)
        index = np.fromiter(chain.from_iterable(balls), dtype=np.intp, count=int(counts.sum()))
        values = _segment_surface_distances(self._xyz[index], self._radius_array[index], s[owner], e[owner])
        starts_at = np.cumsum(counts) - counts
        out[rows] = np.minimum.reduceat(values, starts_at)
        return out

    def nearest_surface_values(self, points):
        """Clearances of many points at once (values of :meth:`nearest_surface`).

        Bit-identical to ``[nearest_surface(p)[0] for p in points]``: the same
        KD-tree nearest centre and ball search per point, with every distance
        evaluated by the same floating-point operations in the same order,
        only batched in NumPy.

        Parameters
        ----------
        points : array_like, shape (m, 3)
            Query positions in Å.

        Returns
        -------
        numpy.ndarray, shape (m,)
            Minimum signed surface distance of each point, Å.
        """
        import numpy as np
        pts = np.asarray(points, dtype=float).reshape(-1, 3)
        if len(pts) == 0:
            return np.zeros(0)
        if self._tree is None:
            return np.array([self._nearest_surface_bruteforce(tuple(map(float, p)))[0] for p in pts])
        nearest_distance, nearest_index = self._tree.query(pts, k=1)
        best = nearest_distance - self._radius_array[nearest_index]
        search = np.maximum(np.maximum(nearest_distance, best + self._max_radius), 0.0) + 1e-9
        balls = self._tree.query_ball_point(pts, search, return_sorted=False)
        counts = np.fromiter(map(len, balls), dtype=np.intp, count=len(balls))
        if counts.sum():
            owner = np.repeat(np.arange(len(pts)), counts)
            index = np.fromiter(chain.from_iterable(balls), dtype=np.intp, count=int(counts.sum()))
            values = _point_surface_distances(pts[owner], self._xyz[index], self._radius_array[index])
            nonempty = counts > 0
            starts = (np.cumsum(counts) - counts)[nonempty]
            best[nonempty] = np.minimum(best[nonempty], np.minimum.reduceat(values, starts))
        return best

    def atoms_within_surface_distance(self, point: Coord, cutoff: float) -> tuple[tuple[Atom, float], ...]:
        """Atoms whose surfaces lie within ``cutoff`` of a point.

        Parameters
        ----------
        point : Coord
            Query position in Å.
        cutoff : float
            Maximum signed surface distance in Å (inclusive). Negative values select
            only atoms that enclose the point by at least that depth.

        Returns
        -------
        tuple of (Atom, float)
            ``(atom, surface_distance)`` pairs sorted by increasing distance.
        """
        if self._tree is None:
            return self._atoms_within_surface_distance_bruteforce(point, cutoff)
        search_radius = max(0.0, cutoff + self._max_radius)
        indices = self._tree.query_ball_point(point, search_radius)
        hits = []
        for index in indices:
            atom = self.atoms[index]
            clearance = self.atom_surface_distance(point, atom)
            if clearance <= cutoff:
                hits.append((atom, clearance))
        return tuple(sorted(hits, key=lambda item: item[1]))

    def atoms_within_distance(self, point: Coord, cutoff: float) -> tuple[Atom, ...]:
        """Return atoms whose centres are within ``cutoff`` of a point.

        Parameters
        ----------
        point : Coord
            Query position in Å.
        cutoff : float
            Maximum centre-to-point distance in Å (inclusive). A negative cutoff
            returns an empty tuple.

        Returns
        -------
        tuple of Atom
            Matching atoms. The order is the KD-tree's result order (or input order
            without SciPy) and is not sorted by distance.
        """

        if cutoff < 0:
            return ()
        if self._tree is None:
            return tuple(atom for atom in self.atoms if distance(point, atom.coord) <= cutoff)
        indices = self._tree.query_ball_point(point, cutoff)
        return tuple(self.atoms[int(index)] for index in indices)

    def _nearest_surface_kdtree(self, point: Coord) -> tuple[float, Atom]:
        """KD-tree nearest surface: nearest centre, then a ball search of radius
        ``max(d_centre, best_clearance + max_radius)`` to catch larger atoms.
        """
        nearest_distance, nearest_index = self._tree.query(point, k=1)
        nearest_atom = self.atoms[int(nearest_index)]
        best_clearance = float(nearest_distance) - self._radii[nearest_atom]
        search_radius = max(float(nearest_distance), best_clearance + self._max_radius, 0.0) + 1e-9
        indices = self._tree.query_ball_point(point, search_radius)
        if not indices:
            return best_clearance, nearest_atom
        # Same values as atom_surface_distance, batched; the first atom (in
        # ball order) attaining a strictly smaller clearance wins, exactly as
        # the sequential scan did.
        xyz = self._xyz[indices]
        dx = point[0] - xyz[:, 0]
        dy = point[1] - xyz[:, 1]
        dz = point[2] - xyz[:, 2]
        values = np.sqrt(dx * dx + dy * dy + dz * dz) - self._radius_array[indices]
        k = int(values.argmin())
        if values[k] < best_clearance:
            return float(values[k]), self.atoms[indices[k]]
        return best_clearance, nearest_atom

    def _nearest_surface_bruteforce(self, point: Coord) -> tuple[float, Atom]:
        """Nearest surface by scanning every atom; ties keep the earliest atom."""
        best_atom = self.atoms[0]
        best_clearance = distance(point, best_atom.coord) - self._radii[best_atom]
        for atom in self.atoms[1:]:
            clearance = distance(point, atom.coord) - self._radii[atom]
            if clearance < best_clearance:
                best_atom = atom
                best_clearance = clearance
        return best_clearance, best_atom

    def _atoms_within_surface_distance_bruteforce(self, point: Coord, cutoff: float) -> tuple[tuple[Atom, float], ...]:
        """Brute-force version of :meth:`atoms_within_surface_distance`."""
        hits: list[tuple[Atom, float]] = []
        for atom in self.atoms:
            clearance = self.atom_surface_distance(point, atom)
            if clearance <= cutoff:
                hits.append((atom, clearance))
        return tuple(sorted(hits, key=lambda item: item[1]))


def _point_surface_distances(points, centres, radii):
    """``distance(point, centre) - radius`` row by row, as :func:`crevice.geometry.distance` computes it."""
    import numpy as np
    dx = points[:, 0] - centres[:, 0]
    dy = points[:, 1] - centres[:, 1]
    dz = points[:, 2] - centres[:, 2]
    return np.sqrt(dx * dx + dy * dy + dz * dz) - radii


def _segment_surface_distances(centres, radii, start, end):
    """``point_segment_distance(centre, start, end) - radius`` for many centres.

    The operations mirror :func:`crevice.geometry.point_segment_distance`
    element by element, so the values are bit-identical. ``start`` and
    ``end`` are one segment (3-sequences) or one segment per centre
    (``(m, 3)`` arrays).
    """
    import numpy as np
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)
    d = end - start
    d0, d1, d2 = d[..., 0], d[..., 1], d[..., 2]
    length_squared = d0 * d0 + d1 * d1 + d2 * d2
    w0 = centres[:, 0] - start[..., 0]
    w1 = centres[:, 1] - start[..., 1]
    w2 = centres[:, 2] - start[..., 2]
    zero = length_squared == 0
    safe = np.where(zero, 1.0, length_squared)
    fraction = np.maximum(0.0, np.minimum(1.0, (w0 * d0 + w1 * d1 + w2 * d2) / safe))
    fraction = np.where(zero, 0.0, fraction)
    e0 = w0 - d0 * fraction
    e1 = w1 - d1 * fraction
    e2 = w2 - d2 * fraction
    return np.sqrt(e0 * e0 + e1 * e1 + e2 * e2) - radii

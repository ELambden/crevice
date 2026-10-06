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

import math

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
        try:
            from scipy.spatial import cKDTree  # type: ignore[import-not-found]
        except ImportError:
            self._coords = None
        else:
            self._coords = [atom.coord for atom in atoms]
            self._tree = cKDTree(self._coords)

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
        for index in self._tree.query_ball_point(midpoint, radius, eps=0.0):
            atom = self.atoms[int(index)]
            best = min(best, point_segment_distance(atom.coord, start, end) - self._radii[atom])
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
        best_atom = nearest_atom
        for index in indices:
            atom = self.atoms[index]
            clearance = self.atom_surface_distance(point, atom)
            if clearance < best_clearance:
                best_atom = atom
                best_clearance = clearance
        return best_clearance, best_atom

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

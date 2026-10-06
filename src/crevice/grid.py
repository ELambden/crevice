"""Continuous atom-clearance test for edges between neighbouring grid points.

Grid-based searches (:func:`crevice.voids.void_cast`,
:func:`crevice.tunnels.find_tunnels`) treat two adjacent grid points as
connected only if a probe of radius ``min_radius`` can move along the whole
straight segment between them without entering an atom sphere. Checking only
the two endpoints would let paths jump through thin atom walls; this module
checks the segment analytically.

Definitions
    *Segment clearance* is the minimum, over every point of the closed
    segment, of the signed distance to the nearest atom van der Waals
    surface (Å), computed by
    :meth:`crevice.spatial.SpatialIndex.segment_clearance`. An edge is
    allowed when both endpoint clearances and the segment clearance are at
    least ``min_radius``.

Main entry point: :class:`GridConnectivity`.
"""

from __future__ import annotations

import math

from .models import Coord
from .spatial import SpatialIndex

GridIndex = tuple[int, int, int]


class GridConnectivity:
    """Segment-clearance tests for grid edges, with a bounded cache.

    The graph traversal itself lives in the calling module; this class only
    answers "may a probe of radius ``min_radius`` pass along this edge?".

    Parameters
    ----------
    spatial : SpatialIndex
        Atom-surface query structure for the analysed atoms.
    origin : Coord
        Position of grid index ``(0, 0, 0)`` in Å.
    spacing : float
        Grid spacing in Å.
    min_radius : float
        Required clearance in Å (the probe radius for the search).
    cache_limit : int, default 32768
        Maximum number of exact edge clearances kept in memory. Once full,
        further results are computed but not stored, so memory stays bounded
        and results are unchanged.
    basis : sequence of three Coord, optional
        Rows of an orthonormal grid basis. When given, index ``(i, j, k)`` maps
        to ``origin + spacing * (i*b0 + j*b1 + k*b2)``; otherwise the grid is
        aligned with the x, y and z axes.

    Attributes
    ----------
    segment_queries : int
        Number of exact segment-clearance calculations performed.
    endpoint_certificates : int
        Number of edges accepted by the endpoint bound without an exact query.
    """

    def __init__(self, spatial: SpatialIndex, origin: Coord, spacing: float,
                 min_radius: float, *, cache_limit: int = 32768, basis=None):
        self.basis = basis
        self.spatial = spatial
        self.origin = origin
        self.spacing = spacing
        self.min_radius = min_radius
        self.cache_limit = cache_limit
        self._cache: dict[tuple[GridIndex, GridIndex], float] = {}
        self.segment_queries = 0
        self.endpoint_certificates = 0

    def position(self, index: GridIndex) -> Coord:
        """Cartesian position (Å) of a grid index.

        Parameters
        ----------
        index : tuple of int
            Grid index ``(i, j, k)``.

        Returns
        -------
        Coord
        """
        if self.basis is not None:
            return tuple(self.origin[j] + self.spacing * sum(index[k]*self.basis[k][j] for k in range(3))
                         for j in range(3))
        return tuple(o + self.spacing * i for o, i in zip(self.origin, index))

    def clearance(self, left: GridIndex, right: GridIndex) -> float:
        """Exact minimum atom-surface clearance along the segment between two nodes.

        The pair is ordered before lookup, so ``clearance(a, b) == clearance(b, a)``
        and both share one cache entry.

        Parameters
        ----------
        left, right : tuple of int
            Grid indices of the segment endpoints.

        Returns
        -------
        float
            Signed segment clearance in Å (no probe radius subtracted).
        """
        key = (left, right) if left <= right else (right, left)
        if key not in self._cache:
            value = self.spatial.segment_clearance(self.position(key[0]), self.position(key[1]))
            self.segment_queries += 1
            if len(self._cache) < self.cache_limit:
                self._cache[key] = value
            return value
        return self._cache[key]

    def allows(self, left: GridIndex, right: GridIndex, left_clearance: float,
               right_clearance: float) -> bool:
        """Decide whether a probe of radius ``min_radius`` can traverse an edge.

        Parameters
        ----------
        left, right : tuple of int
            Grid indices of the edge endpoints (normally 26-neighbours).
        left_clearance, right_clearance : float
            Atom-surface clearances (Å) already known at the two endpoints.

        Returns
        -------
        bool
            ``False`` if either endpoint clearance is below ``min_radius``;
            otherwise ``True`` exactly when the segment clearance is at least
            ``min_radius``.

        Notes
        -----
        Distance to a union of spheres is 1-Lipschitz, so every point on a segment
        of length ``L`` has clearance at least ``(c_left + c_right - L) / 2``. When
        that bound already exceeds ``min_radius`` (by more than 1e-12 Å) the edge is
        accepted without an exact query and counted in
        ``endpoint_certificates``. The shortcut never changes the answer.
        """
        if min(left_clearance, right_clearance) < self.min_radius:
            return False
        length = self.spacing * math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
        # Surface clearance is 1-Lipschitz. Combining both endpoint bounds
        # certifies some edges without calculating their exact capacity.
        if (left_clearance + right_clearance - length) / 2 > self.min_radius + 1e-12:
            self.endpoint_certificates += 1
            return True
        return self.clearance(left, right) >= self.min_radius

    def to_dict(self) -> dict:
        """Describe the edge rule and query statistics for result metadata.

        Returns
        -------
        dict
            ``neighbors`` (always 26), ``edge_rule``, ``min_segment_clearance_A``,
            ``segment_queries``, ``endpoint_certificates``, ``cached_edges`` and
            ``cache_limit``.
        """
        return {"neighbors": 26, "edge_rule": "continuous_atom_sphere_clearance",
                "min_segment_clearance_A": self.min_radius,
                "segment_queries": self.segment_queries,
                "endpoint_certificates": self.endpoint_certificates,
                "cached_edges": len(self._cache), "cache_limit": self.cache_limit}

"""Dependency-free vector geometry on 3-tuples.

These helpers operate on :data:`~crevice.models.Coord` tuples ``(x, y, z)``
and use only the standard library, so the parsers and core analyses work
without NumPy. Distances are in the units of the input coordinates (Å
throughout CREVICE).

Main entry points: vector arithmetic (:func:`add`, :func:`sub`,
:func:`dot`, :func:`cross`, :func:`normalize`), distances
(:func:`distance`, :func:`point_segment_distance`), axis handling
(:func:`axis_from_name`, :func:`principal_axes`, :func:`basis_from_direction`,
:func:`molecular_basis`) and :func:`parse_coord` for command-line input.

Examples
--------
>>> from crevice.geometry import add, cross, distance, normalize
>>> add((1.0, 2.0, 3.0), (1.0, 1.0, 1.0))
(2.0, 3.0, 4.0)
>>> cross((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
(0.0, 0.0, 1.0)
>>> distance((0.0, 0.0, 0.0), (3.0, 4.0, 0.0))
5.0
>>> normalize((0.0, 0.0, 2.0))
(0.0, 0.0, 1.0)
"""

from __future__ import annotations

import math
from typing import Iterable

from .models import Atom, Coord


def add(a: Coord, b: Coord) -> Coord:
    """Return the component-wise sum ``a + b``.

    Parameters
    ----------
    a, b : Coord
        3-vectors (or points) in consistent units.

    Returns
    -------
    Coord
    """
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Coord, b: Coord) -> Coord:
    """Return the component-wise difference ``a - b``.

    Parameters
    ----------
    a, b : Coord
        3-vectors (or points) in consistent units.

    Returns
    -------
    Coord
    """
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(a: Coord, factor: float) -> Coord:
    """Return ``a`` multiplied by the scalar ``factor``.

    Parameters
    ----------
    a : Coord
    factor : float

    Returns
    -------
    Coord
    """
    return (a[0] * factor, a[1] * factor, a[2] * factor)


def dot(a: Coord, b: Coord) -> float:
    """Return the dot product of ``a`` and ``b``.

    Parameters
    ----------
    a, b : Coord
        3-vectors (or points) in consistent units.

    Returns
    -------
    float
    """
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Coord, b: Coord) -> Coord:
    """Return the right-handed cross product ``a x b``.

    Parameters
    ----------
    a, b : Coord
        3-vectors (or points) in consistent units.

    Returns
    -------
    Coord
    """
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a: Coord) -> float:
    """Return the Euclidean length of ``a``.

    Parameters
    ----------
    a : Coord

    Returns
    -------
    float
    """
    return math.sqrt(dot(a, a))


def normalize(a: Coord) -> Coord:
    """Return ``a`` scaled to unit length.

    Raises
    ------
    ValueError
        If ``a`` has exactly zero length.
    """
    length = norm(a)
    if length == 0:
        raise ValueError("Cannot normalize a zero-length vector")
    return (a[0] / length, a[1] / length, a[2] / length)


def distance(a: Coord, b: Coord) -> float:
    """Return the Euclidean distance between points ``a`` and ``b``.

    Parameters
    ----------
    a, b : Coord
        3-vectors (or points) in consistent units.

    Returns
    -------
    float
        Same units as the inputs (Å for atom coordinates).
    """
    return norm(sub(a, b))


def point_segment_distance(point: Coord, start: Coord, end: Coord) -> float:
    """Euclidean distance from a point to a closed line segment.

    The point is projected onto the segment and the projection parameter is
    clamped to ``[0, 1]``, so the nearest point may be an endpoint. A
    zero-length segment is treated as the single point ``start``.

    Parameters
    ----------
    point, start, end : Coord
        Query point and segment endpoints.

    Returns
    -------
    float
        Distance in the input units.

    Examples
    --------
    >>> from crevice.geometry import point_segment_distance
    >>> point_segment_distance((0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    1.0
    >>> point_segment_distance((3.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    2.0
    """
    direction = sub(end, start)
    length_squared = dot(direction, direction)
    if length_squared == 0:
        return distance(point, start)
    fraction = max(0.0, min(1.0, dot(sub(point, start), direction) / length_squared))
    return distance(sub(point, start), scale(direction, fraction))


def centroid(points: Iterable[Coord]) -> Coord:
    """Unweighted mean of a set of points.

    Parameters
    ----------
    points : iterable of Coord

    Returns
    -------
    Coord

    Raises
    ------
    ValueError
        If ``points`` is empty.
    """
    coords = tuple(points)
    if not coords:
        raise ValueError("Cannot compute centroid of an empty point set")
    n = float(len(coords))
    return (
        sum(point[0] for point in coords) / n,
        sum(point[1] for point in coords) / n,
        sum(point[2] for point in coords) / n,
    )


def project_onto_axis(point: Coord, origin: Coord, direction: Coord) -> float:
    """Signed coordinate of ``point`` along an axis.

    Returns ``dot(point - origin, direction)``. The result is a length only if
    ``direction`` is a unit vector; it is not normalised here.
    """
    return dot(sub(point, origin), direction)


def point_on_axis(origin: Coord, direction: Coord, t: float) -> Coord:
    """Return ``origin + t * direction``.

    Parameters
    ----------
    origin : Coord
    direction : Coord
        Axis direction; ``t`` is a length only if this is a unit vector.
    t : float
        Coordinate along the axis.

    Returns
    -------
    Coord
    """
    return add(origin, scale(direction, t))


def basis_from_direction(direction: Coord) -> tuple[Coord, Coord]:
    """Two unit vectors perpendicular to ``direction`` and to each other.

    The first vector is ``normalize(direction x ref)`` with ``ref = +z``, or
    ``+y`` when the direction is within about 25.8 degrees of the z axis
    (``|cos| > 0.9``). The second is ``direction x u``. Together with the unit
    direction they form a right-handed orthonormal frame
    ``(u, v, direction)``. The choice depends on the laboratory axes; see
    :func:`molecular_basis` for a molecule-attached alternative.

    Parameters
    ----------
    direction : Coord
        Non-zero axis direction (need not be normalised).

    Returns
    -------
    tuple of Coord
        ``(u, v)``.

    Raises
    ------
    ValueError
        If ``direction`` has zero length.
    """
    unit = normalize(direction)
    reference = (0.0, 0.0, 1.0)
    if abs(dot(unit, reference)) > 0.9:
        reference = (0.0, 1.0, 0.0)
    u = normalize(cross(unit, reference))
    v = normalize(cross(unit, u))
    return u, v


def axis_from_name(axis: str | Coord, atoms: Iterable[Atom]) -> tuple[Coord, Coord, str]:
    """Resolve an axis specification into an origin and a unit direction.

    The origin is always the unweighted centroid of ``atoms``.

    Parameters
    ----------
    axis : {"x", "y", "z", "auto"} or Coord
        ``"x"``, ``"y"`` or ``"z"`` (case-insensitive) select a laboratory
        axis. ``"auto"`` uses the principal axis of largest positional variance
        (:func:`principal_axes`). A 3-vector is normalised and used as given.
    atoms : iterable of Atom
        Atoms that define the centroid (and the shape for ``"auto"``).

    Returns
    -------
    origin : Coord
        Centroid of the atoms, in Å.
    direction : Coord
        Unit axis direction.
    label : str
        ``"x"``, ``"y"``, ``"z"``, ``"vector"`` or ``"auto:principal-shape"``.

    Raises
    ------
    ValueError
        If ``atoms`` is empty, the name is not recognised, or a vector has zero
        length.

    Notes
    -----
    ``"auto"`` is only a shape proposal. :func:`crevice.channels.pore_profile`
    tests all three principal directions when it resolves a connected channel.
    """
    selected = tuple(atoms)
    if not selected:
        raise ValueError("Cannot infer an axis from an empty atom set")
    origin = (
        sum(atom.x for atom in selected) / len(selected),
        sum(atom.y for atom in selected) / len(selected),
        sum(atom.z for atom in selected) / len(selected),
    )
    if not isinstance(axis, str):
        return origin, normalize(axis), "vector"

    name = axis.lower()
    if name == "x":
        return origin, (1.0, 0.0, 0.0), "x"
    if name == "y":
        return origin, (0.0, 1.0, 0.0), "y"
    if name == "z":
        return origin, (0.0, 0.0, 1.0), "z"
    if name != "auto":
        raise ValueError(f"Unknown axis {axis!r}; use x, y, z, auto, or a 3-vector")

    # Shape proposal only; pore profiling checks all principal directions.
    return origin, principal_axes(selected)[0], "auto:principal-shape"


def principal_axes(atoms: Iterable[Atom]) -> tuple[Coord, Coord, Coord]:
    """Principal axes of the atom positions, largest variance first.

    The 3x3 positional covariance matrix (about the centroid, unweighted) is
    diagonalised with cyclic Jacobi rotations (at most 50 sweeps), without
    NumPy. The eigenvectors are returned in order of decreasing eigenvalue.

    Parameters
    ----------
    atoms : iterable of Atom

    Returns
    -------
    tuple of Coord
        Three orthonormal unit vectors.

    Notes
    -----
    Each vector's sign is fixed from atom order: the first atom whose
    projection exceeds half the maximum absolute projection is made positive.
    This makes the axes rotate with the molecule under a rigid motion.
    Degenerate eigenvalues (for example a symmetric channel) do not define a
    unique direction, so the returned vectors within a degenerate subspace are
    arbitrary.

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.geometry import principal_axes
    >>> rod = [Atom(i, "C", "UNK", "A", 1, 0.0, 0.0, float(4 - i), "C") for i in range(5)]
    >>> principal_axes(rod)[0]
    (0.0, 0.0, 1.0)
    """
    atoms = tuple(atoms)
    origin = centroid(a.coord for a in atoms)
    coords = [sub(a.coord, origin) for a in atoms]
    matrix = [[sum(p[i] * p[j] for p in coords) / len(coords)
               for j in range(3)] for i in range(3)]
    vectors = [[float(i == j) for j in range(3)] for i in range(3)]
    for _ in range(50):
        p, q = max(((0, 1), (0, 2), (1, 2)), key=lambda ij: abs(matrix[ij[0]][ij[1]]))
        if abs(matrix[p][q]) <= 1e-13 * max(1.0, sum(abs(matrix[i][i]) for i in range(3))):
            break
        angle = 0.5 * math.atan2(2 * matrix[p][q], matrix[q][q] - matrix[p][p])
        c, s = math.cos(angle), math.sin(angle)
        app, aqq, apq = matrix[p][p], matrix[q][q], matrix[p][q]
        matrix[p][p] = c*c*app - 2*s*c*apq + s*s*aqq
        matrix[q][q] = s*s*app + 2*s*c*apq + c*c*aqq
        matrix[p][q] = matrix[q][p] = 0.0
        for k in range(3):
            if k not in (p, q):
                akp, akq = matrix[k][p], matrix[k][q]
                matrix[k][p] = matrix[p][k] = c*akp - s*akq
                matrix[k][q] = matrix[q][k] = s*akp + c*akq
            vkp, vkq = vectors[k][p], vectors[k][q]
            vectors[k][p], vectors[k][q] = c*vkp - s*vkq, s*vkp + c*vkq
    result = []
    for j in sorted(range(3), key=lambda j: -matrix[j][j]):
        direction = normalize(tuple(vectors[i][j] for i in range(3)))
        projections = [dot(p, direction) for p in coords]
        extent = max(abs(t) for t in projections)
        reference = next((t for t in projections if abs(t) > extent * 0.5), 1.0)
        result.append(scale(direction, -1.0 if reference < 0 else 1.0))
    return tuple(result)


def molecular_basis(direction: Coord, atoms: Iterable[Atom]) -> tuple[Coord, Coord]:
    """Transverse axes attached to the molecule rather than to the laboratory.

    The first atom (in input order) whose offset from the centroid has a
    component perpendicular to ``direction`` longer than 1e-6 Å defines
    ``u`` (that perpendicular component, normalised); ``v = direction x u``.
    Falls back to :func:`basis_from_direction` if every atom lies on the axis.

    Parameters
    ----------
    direction : Coord
        Unit axis direction.
    atoms : iterable of Atom

    Returns
    -------
    tuple of Coord
        ``(u, v)``, both unit length and perpendicular to ``direction``.
    """
    atoms = tuple(atoms)
    origin = centroid(a.coord for a in atoms)
    for atom in atoms:
        offset = sub(atom.coord, origin)
        transverse = sub(offset, scale(direction, dot(offset, direction)))
        if norm(transverse) > 1e-6:
            u = normalize(transverse)
            return u, normalize(cross(direction, u))
    return basis_from_direction(direction)


def parse_coord(text: str) -> Coord:
    """Parse ``"x,y,z"`` text into a coordinate tuple.

    Parameters
    ----------
    text : str
        Three comma-separated numbers; whitespace around them is ignored.

    Returns
    -------
    Coord

    Raises
    ------
    ValueError
        If there are not exactly three fields or a field is not a number.

    Examples
    --------
    >>> from crevice.geometry import parse_coord
    >>> parse_coord("1.5, -2, 0")
    (1.5, -2.0, 0.0)
    """
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 3:
        raise ValueError("Expected a coordinate formatted as x,y,z")
    return (float(parts[0]), float(parts[1]), float(parts[2]))


def padded_box(min_corner: Coord, max_corner: Coord, padding: float) -> tuple[Coord, Coord]:
    """Enlarge an axis-aligned box by ``padding`` on every side.

    Returns
    -------
    tuple of Coord
        ``(min_corner - padding, max_corner + padding)`` component-wise.
    """
    return (
        (min_corner[0] - padding, min_corner[1] - padding, min_corner[2] - padding),
        (max_corner[0] + padding, max_corner[1] + padding, max_corner[2] + padding),
    )

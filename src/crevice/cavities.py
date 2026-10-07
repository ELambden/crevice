"""Detection of buried cavities by flood-filling a clearance grid.

A *cavity* here is a connected set of grid points where a probe of radius
``min_radius`` fits between atoms (atom-surface clearance at least
``min_radius``) and from which no clearance-valid grid path leads to the edge
of the padded box. Enclosure is decided by flood-filling the free grid from
the box boundary; whatever the flood cannot reach is buried.

Main entry points: :func:`detect_cavities` (command ``crevice cavities``),
:func:`classify_voids`, :func:`cavity_surface_points`, :func:`cavity_volume`.

Assumptions and limitations
    Volumes are grid-cell counts times ``spacing**3`` and change with the
    spacing (2 Å by default) and ``min_radius`` (1.4 Å by default).
    "Enclosed" means enclosed at this resolution and probe size: a cavity
    connected to bulk by a channel narrower than the probe, or than the grid
    can resolve, is reported as buried. Open pockets and channels are
    deliberately not cavities here; use :func:`crevice.rolling.rolling_probe_cast`
    (command ``crevice cast``) for pockets and open regions. A cavity is a
    geometric finding, not evidence of function or of water occupancy.
"""

from __future__ import annotations

import math
from typing import Sequence

from .analysis import detect_cavities as _detect_cavities
from .models import Cavity, Coord, StructureFrame
from .radii import RadiusSet, radii_option


@radii_option
def detect_cavities(frame: StructureFrame, *, radii: RadiusSet | str | None = None, **kwargs) -> tuple[Cavity, ...]:
    """Detect enclosed grid voids (optionally also open pocket candidates).

    Thin wrapper around :func:`crevice.analysis.detect_cavities`, which calls
    :func:`crevice.voids.void_cast` with ``mode="cavity"`` (or ``"all"``) and
    converts each resulting component to a :class:`~crevice.models.Cavity`.

    Parameters
    ----------
    frame : StructureFrame
    **kwargs
        Passed to :func:`crevice.analysis.detect_cavities`:

        spacing : float, default 2.0
            Requested grid spacing in Å (coarsened automatically if the grid
            would exceed ``max_grid_points``).
        min_radius : float, default 1.4
            Required atom-surface clearance (probe radius) in Å for every grid
            point and every grid edge.
        max_cavities : int, default 10
            Maximum number of cavities returned, largest first.
        shell_radius : float, default 8.0
            Radius in Å of the neighbourhood used for the octant-burial test and
            for listing nearby residues.
        min_octants : int, default 4
            With ``enclosed_only=False``, a grid point is kept only if atom
            centres within ``shell_radius`` occupy at least this many of the
            eight octants around it. Ignored when ``enclosed_only=True``.
        max_grid_points : int, default 75000
            Grid-size limit that triggers automatic coarsening.
        include_hydrogen : bool, default False
            Treat hydrogen atoms as obstacles.
        include_hetero : bool, default False
            Treat HETATM atoms (ligands, waters, ions) as obstacles.
        enclosed_only : bool, default True
            ``True``: only components the exterior flood fill cannot reach.
            ``False``: every octant-buried free component, including open
            pockets and channels.
        focus_points : sequence of Coord, optional
            Keep only grid points within ``focus_radius`` of at least one of
            these positions (Å).
        focus_radius : float, default 8.0
            Radius (Å) around each focus point.
        min_component_volume : float, default 0.0
            Discard components smaller than this (Å³, inclusive bound).
        max_component_volume : float, optional
            Discard components larger than this (Å³, inclusive bound).
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
    tuple of Cavity
        Largest volume first (ties broken by the smallest grid coordinate).
        ``radius`` is the largest clearance in the cavity, ``volume`` is
        ``grid_points * spacing**3`` in Å³, ``center`` is the
        clearance-weighted centroid and ``score`` equals ``volume`` (a ranking
        value only).

    Raises
    ------
    ValueError
        For invalid parameters (for example ``max_cavities < 1``).

    See Also
    --------
    crevice.voids.void_cast : The underlying grid cast, with full metadata.
    crevice.rolling.rolling_probe_cast : Rolling-probe casts, including pockets.

    Examples
    --------
    A closed spherical shell of 250 carbon atoms (radius 7 Å) has one cavity:

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.cavities import detect_cavities
    >>> def shell_atom(i, n=250, golden=math.pi * (3 - math.sqrt(5))):
    ...     y = 1 - (2 * i + 1) / n
    ...     r = math.sqrt(1 - y * y)
    ...     return Atom(i + 1, "C", "ALA", "A", i + 1, 7 * math.cos(golden * i) * r,
    ...                 7 * y, 7 * math.sin(golden * i) * r, "C")
    >>> shell = StructureFrame(tuple(shell_atom(i) for i in range(250)))
    >>> [(c.id, c.volume, c.grid_points) for c in detect_cavities(shell)]
    [(1, 256.0, 32)]
    """

    return _detect_cavities(frame, **kwargs)


def cavity_volume(cavity: Cavity) -> float:
    """Return the cavity's grid volume in Å³.

    This is ``Cavity.volume``, i.e. ``grid_points * spacing**3``; it depends on
    the grid spacing and probe radius used for detection.
    """

    return cavity.volume


def cavity_surface_points(cavity: Cavity, *, samples: int = 64) -> tuple[Coord, ...]:
    """Evenly spread points on the cavity's representative sphere, for display.

    The sphere has centre ``Cavity.center`` and radius ``Cavity.radius`` (the
    largest clearance in the cavity). Points follow a golden-angle (Fibonacci)
    spiral. The sphere is only a marker; it is not the cavity's shape.

    Parameters
    ----------
    cavity : Cavity
    samples : int, default 64
        Number of points.

    Returns
    -------
    tuple of Coord
        ``samples`` positions in Å.

    Raises
    ------
    ValueError
        If ``samples < 1``.
    """

    if samples < 1:
        raise ValueError("samples must be at least 1")
    points: list[Coord] = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(samples):
        y = 1.0 - (2.0 * index + 1.0) / samples
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden_angle * index
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        points.append(
            (
                cavity.center[0] + cavity.radius * x,
                cavity.center[1] + cavity.radius * y,
                cavity.center[2] + cavity.radius * z,
            )
        )
    return tuple(points)


def classify_voids(
    frame: StructureFrame,
    cavities: Sequence[Cavity],
    *,
    surface_accessibility: float = 6.0,
) -> dict[int, str]:
    """Label cavities by position relative to the structure's bounding box.

    For each cavity, the distance from its centre to the nearest face of the
    axis-aligned bounding box of *all* atom centres in ``frame`` is compared
    with ``surface_accessibility``:

    * distance ``<= surface_accessibility``: ``"surface_proximal"``;
    * otherwise, if ``cavity.radius >= surface_accessibility``:
      ``"large_connected_candidate"``;
    * otherwise: ``"buried_cavity"``.

    Parameters
    ----------
    frame : StructureFrame
    cavities : sequence of Cavity
    surface_accessibility : float, default 6.0
        Threshold in Å, used both as the distance-to-box-face limit and as the
        radius limit.

    Returns
    -------
    dict of int to str
        Label for each cavity id.

    Notes
    -----
    This is a coarse positional heuristic: the bounding box is not the
    molecular surface, and none of the labels tests connectivity to bulk
    (``"large_connected_candidate"`` only means a large radius). Only the
    cavity centre and radius are used. No CREVICE command calls this function.
    """

    min_corner, max_corner = frame.bounding_box()
    classes: dict[int, str] = {}
    for cavity in cavities:
        distance_to_box = min(
            cavity.center[0] - min_corner[0],
            max_corner[0] - cavity.center[0],
            cavity.center[1] - min_corner[1],
            max_corner[1] - cavity.center[1],
            cavity.center[2] - min_corner[2],
            max_corner[2] - cavity.center[2],
        )
        if distance_to_box <= surface_accessibility:
            label = "surface_proximal"
        elif cavity.radius >= surface_accessibility:
            label = "large_connected_candidate"
        else:
            label = "buried_cavity"
        classes[cavity.id] = label
    return classes

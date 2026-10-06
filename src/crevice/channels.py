"""HOLE-like channel radius profiles.

A channel profile answers: along a pore, how large a sphere fits between the
atoms at each position? CREVICE samples points along the channel and records,
for each, the *clearance* (distance from the point to the nearest atom van der
Waals surface, Å) and the *radius* ``max(0, clearance - probe_radius)``.

Method in brief
    :func:`pore_profile` (default mode) proposes axis directions, slices the
    structure into transverse planes, and keeps only free space that is
    laterally enclosed by atoms (an enclosure probe cannot escape to the edge
    of the search disc). Enclosed sections in neighbouring planes are joined
    only by straight segments that clear the probe. The longest joined path
    with two open ends is the channel. Each profile point is the position of
    largest clearance within its section (the best grid point, refined by a
    local search), and the profile stops at the refined mouths. If no such channel is
    found, :class:`crevice.sections.ChannelResolutionError` is raised rather
    than guessing.

Main entry points: :func:`pore_profile`, :func:`profile_along_path`,
:func:`detect_bottlenecks`, :func:`trace_centerline`, :func:`channel_volume`.
The command-line equivalent is ``crevice profile``.

Assumptions and limitations
    Atoms are hard spheres with element radii (:mod:`crevice.radii`).
    Results depend on the section grid spacing, the enclosure probe and the
    axis proposals, and describe one monotone, unbranched passage along the
    chosen axis. A profile is geometry only: a narrow radius is not proof of
    a closed gate, and a wide one is not proof of permeation. No benchmark
    structure in CREVICE has a curated axis, so profiles of real proteins
    are unvalidated measurements.
"""

from __future__ import annotations

from typing import Sequence
import math

from .analysis import pore_profile as _pore_profile
from .geometry import distance, normalize, sub
from .models import ChannelPoint, Coord, PoreProfile, StructureFrame
from .spatial import SpatialIndex
from .radii import RadiusSet, radii_option


@radii_option
def pore_profile(frame: StructureFrame, *, radii: RadiusSet | str | None = None, **kwargs) -> PoreProfile:
    """Resolve a channel and measure its inscribed radius along the axis.

    Thin wrapper around :func:`crevice.analysis.pore_profile`; all keyword
    arguments are passed through unchanged.

    Two modes are selected by ``search_radius``:

    * ``search_radius > 0`` (default 6.0): **connected-channel** mode
      (``method="axial-connected"``). For each candidate direction (the three
      principal axes when ``axis="auto"``, otherwise the given axis), the
      structure is cut into planes about every 0.75 Å. In each plane, grid
      points with spacing ``section_spacing`` inside a square of half-width
      ``search_radius`` are kept where a sphere of radius
      ``max(probe_radius, enclosure_radius)`` fits; connected groups of such
      points that do not reach the square's edge are *enclosed sections*.
      Sections in adjacent planes are joined when the straight segment between
      them clears the same sphere. The joined path with the longest axial span
      that has open space beyond both ends is the channel. It is resampled
      between refined mouths, and where a plane contains the lumen plus a small
      wall pocket the widest full-length path is followed. By default the
      enclosure probe is chosen automatically (see ``enclosure_radius``).
    * ``search_radius == 0``: **fixed-axis scan** (``method="axial-refined"``,
      ``metadata["status"] = "unvalidated_axis_scan"``). ``samples`` points are
      placed evenly on the straight axis through ``origin`` between the extreme
      atom projections (extended by ``padding``), and the clearance is read at
      each point without moving it. No channel connectivity is asserted.

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    **kwargs
        Passed to :func:`crevice.analysis.pore_profile`:

        axis : {"auto", "x", "y", "z"} or Coord, default "auto"
            Channel direction. ``"auto"`` tries the three principal axes of the
            selected atoms and fails if two give comparable spans.
        origin : Coord, optional
            A point inside the intended channel (Å). Defaults to the atom
            centroid. In connected mode a supplied origin must have clearance
            greater than the enclosure and measurement probes.
        samples : int, default 81
            Number of profile points. In connected mode this is a minimum:
            more points are used if needed to keep steps at or below 0.75 Å.
            Must be at least 2.
        padding : float, default 0.0
            Å added beyond both atom extremes; only allowed with
            ``search_radius=0``.
        search_radius : float, default 6.0
            Half-width (Å) of the square transverse search region around the
            axis in connected mode; 0 selects the fixed-axis scan.
        refinement_steps : int, default 4
            Refinement iterations for section centres and mouth positions.
        probe_radius : float, default 0.0
            Probe radius (Å) subtracted from each clearance to give ``radius``.
        include_hydrogen : bool, default False
            Treat hydrogen atoms as obstacles.
        include_hetero : bool, default True
            Treat HETATM atoms (ligands, waters, ions) as obstacles.
        section_spacing : float, default 0.5
            Transverse grid spacing (Å) used to find sections. Smaller values
            resolve narrower passages but cost more.
        enclosure_radius : float or None, default None
            Radius (Å) of the sphere used for the lateral-enclosure and
            connectivity tests (the larger of this and ``probe_radius`` is
            used). ``None`` chooses it automatically: the channel is resolved
            at every probe of a 0.5-3.0 Å ladder, results are grouped into
            channels, and the smallest probe of the channel that resolves at
            the most probes (at least two) is used; the choice, every probe
            tried and the reason are recorded in
            ``metadata["enclosure_probe_selection"]``. A number runs exactly
            that probe (0.8 Å was the fixed default before).
        lateral_exits : bool, default False
            Accept a capped end whose straight axial exit is blocked if a
            lateral exit to bulk exists; every distinct exit leg is reported
            in ``metadata["exits"]`` (see :mod:`crevice.channel_exits`).
        exit_bulk_radius : float, default 6.0
            Rolling bulk-probe radius (Å) for lateral exits; exits whose exit
            points are more than twice this apart are distinct.
        exit_spacing : float, default 0.5
            Lattice spacing (Å) of the lateral-exit search.
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
    PoreProfile
        Points ordered along the axis. ``t`` is the signed axial position (Å)
        from ``axis_origin``; ``radius`` is ``max(0, raw_clearance -
        probe_radius)``; ``nearest_residue`` names the residue whose atom limits
        each point. ``metadata`` records the settings and, in connected mode,
        ``status="resolved"``, mouth positions, section areas (Å²), arc lengths,
        the continuous segment bottleneck and the candidate-axis diagnostics.

    Raises
    ------
    crevice.sections.ChannelResolutionError
        If no laterally enclosed channel with two open ends is found, or the
        axis or region is ambiguous (connected mode only). With the automatic
        probe the error lists every probe tried and carries the selection
        record as ``probe_selection``.
    ValueError
        For invalid parameter values (for example ``samples < 2``, negative
        radii, non-finite origin, or ``padding`` with connected mode).

    See Also
    --------
    profile_along_path : Clearance along a path you supply.
    detect_bottlenecks : Local minima of a profile.
    crevice.trajectory.analyze_trajectory : Profiles for every frame.

    Notes
    -----
    In the returned profile, the smallest ``radius`` among the samples
    (:attr:`~crevice.models.PoreProfile.min_radius`) can be larger than the true constriction
    between samples. Connected mode also records
    ``metadata["continuous_bottleneck_radius_A"]``, the minimum clearance along
    the straight segments joining consecutive samples minus ``probe_radius``,
    which accounts for that.

    Examples
    --------
    A synthetic cylinder of carbon atoms (radius 6 Å, 21 Å long) along z:

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.channels import pore_profile
    >>> atoms = [Atom(12 * k + j + 1, "C", "ALA", "A", k + 1,
    ...               6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6),
    ...               -10.5 + 1.5 * k, "C")
    ...          for k in range(15) for j in range(12)]
    >>> profile = pore_profile(StructureFrame(tuple(atoms)), axis="z")
    >>> profile.method, profile.metadata["status"]
    ('axial-connected', 'resolved')
    >>> round(profile.min_radius, 2)          # 6.0 Å minus the 1.70 Å carbon radius
    4.3
    """

    return _pore_profile(frame, **kwargs)


@radii_option
def trace_centerline(frame: StructureFrame, *, radii: RadiusSet | str | None = None, **kwargs) -> tuple[Coord, ...]:
    """Return the channel centre-line positions found by :func:`pore_profile`.

    Parameters
    ----------
    frame : StructureFrame
    **kwargs
        Passed to :func:`pore_profile` (see there for every option).
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
    tuple of Coord
        The ``position`` of every profile point, in profile order (Å).

    Raises
    ------
    crevice.sections.ChannelResolutionError, ValueError
        As for :func:`pore_profile`.
    """

    return tuple(point.position for point in pore_profile(frame, **kwargs).points)


@radii_option
def profile_along_path(
    frame: StructureFrame,
    path: Sequence[Coord],
    *,
    probe_radius: float = 0.0,
    include_hydrogen: bool = False,
    include_hetero: bool = True,
    radii: RadiusSet | str | None = None,
) -> PoreProfile:
    """Measure atom-surface clearance at each point of a path you supply.

    No centring, smoothing or connectivity test is applied: each coordinate is
    evaluated where it is. Use this to profile a path from another tool, a
    tunnel, or a hand-picked line.

    Parameters
    ----------
    frame : StructureFrame
        Structure whose atoms are the obstacles.
    path : sequence of Coord
        Ordered sample positions in Å; at least one.
    probe_radius : float, default 0.0
        Å subtracted from each clearance to give ``radius`` (clipped at zero).
    include_hydrogen : bool, default False
        Treat hydrogen atoms as obstacles.
    include_hetero : bool, default True
        Treat HETATM atoms as obstacles.
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
    PoreProfile
        ``method="explicit-path"``. ``axis_origin`` is the first path point and
        ``axis_direction`` the unit vector from the first to the last point
        (``(0, 0, 1)`` for a single point). Each point's ``t`` is the cumulative
        path length (Å) from the first point.

    Raises
    ------
    ValueError
        If ``path`` is empty, ``probe_radius`` is negative or non-finite, the
        first and last points coincide in a multi-point path, or the atom
        selection is empty.

    Examples
    --------
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.channels import profile_along_path
    >>> frame = StructureFrame((Atom(1, "C", "ALA", "A", 1, 5.0, 0.0, 0.0, "C"),))
    >>> profile = profile_along_path(frame, [(0.0, 0.0, z) for z in (-3.0, 0.0, 3.0)],
    ...                              probe_radius=1.4)
    >>> [round(p.radius, 2) for p in profile.points]   # 5.0 - 1.70 - 1.4 at z = 0
    [2.73, 1.9, 2.73]
    >>> [p.t for p in profile.points]
    [0.0, 3.0, 6.0]
    """

    if len(path) < 1:
        raise ValueError("path must contain at least one coordinate")
    if not math.isfinite(probe_radius) or probe_radius < 0:
        raise ValueError("probe_radius must be finite and non-negative")
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    spatial = SpatialIndex(atoms)
    points: list[ChannelPoint] = []
    t = 0.0
    previous: Coord | None = None
    for index, coord in enumerate(path):
        if previous is not None:
            t += distance(previous, coord)
        raw_clearance, nearest = spatial.nearest_surface(coord)
        points.append(
            ChannelPoint(
                index=index,
                position=coord,
                radius=max(0.0, raw_clearance - probe_radius),
                raw_clearance=raw_clearance,
                t=t,
                nearest_atom_serial=nearest.serial,
                nearest_residue=nearest.residue_key.label,
            )
        )
        previous = coord
    if len(path) > 1:
        direction = normalize(sub(path[-1], path[0]))
    else:
        direction = (0.0, 0.0, 1.0)
    return PoreProfile(
        axis_origin=path[0],
        axis_direction=direction,
        points=tuple(points),
        probe_radius=probe_radius,
        method="explicit-path",
        metadata={
            "path_points": len(path),
            "include_hydrogen": include_hydrogen,
            "include_hetero": include_hetero,
        },
    )


def detect_bottlenecks(
    profile: PoreProfile,
    *,
    threshold: float | None = None,
    max_points: int | None = None,
) -> tuple[ChannelPoint, ...]:
    """Find local minima of the radius along a profile.

    A point is a local minimum when its radius is less than or equal to both
    neighbours (the profile ends count as having an infinitely wide neighbour
    outside), so plateaus report every tied point.

    Parameters
    ----------
    profile : PoreProfile
    threshold : float, optional
        Keep only minima with ``radius <= threshold`` (Å).
    max_points : int, optional
        Return at most this many minima.

    Returns
    -------
    tuple of ChannelPoint
        Minima sorted by increasing radius, then by point index. If no minimum
        passes ``threshold``, the single global bottleneck is returned instead.
        Empty only when the profile has no points.
    """

    if not profile.points:
        return ()
    candidates: list[ChannelPoint] = []
    points = profile.points
    for index, point in enumerate(points):
        left = points[index - 1].radius if index > 0 else float("inf")
        right = points[index + 1].radius if index < len(points) - 1 else float("inf")
        if point.radius <= left and point.radius <= right:
            if threshold is None or point.radius <= threshold:
                candidates.append(point)
    if not candidates:
        candidates = [profile.bottleneck]
    candidates.sort(key=lambda point: (point.radius, point.index))
    if max_points is not None:
        candidates = candidates[:max_points]
    return tuple(candidates)


def channel_volume(profile: PoreProfile) -> float:
    """Circular-section volume estimate of a profile, in Å³.

    Returns :attr:`~crevice.models.PoreProfile.volume_estimate`: trapezoidal integration of
    ``pi * radius**2`` along the sampled path. This treats each section as a
    circle of the probe-corrected inscribed radius, so it underestimates
    non-circular lumens. Use a cast (``crevice cast``,
    :func:`crevice.rolling.rolling_probe_cast`) for a sampled volume.
    """

    return profile.volume_estimate

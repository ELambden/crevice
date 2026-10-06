"""Engines of the pore profile, cavity, tunnel and residue-contact analyses.

The documented entry points for the first three are the thin wrappers
:func:`crevice.channels.pore_profile`, :func:`crevice.cavities.detect_cavities`
and :func:`crevice.tunnels.find_tunnels` (all exported from ``crevice``); they
pass their keyword arguments here unchanged, and their docstrings give the
full method description. :func:`annotate_residues` and
:func:`nearest_residues` are documented here and re-exported by
:mod:`crevice.residues`. All lengths are in Å, atoms are van der Waals
spheres from the radius set in effect, and every result is geometry at the
chosen resolution, not a validated biological assignment.
"""

from __future__ import annotations

import heapq
import itertools
import math
from collections import deque
from dataclasses import replace
from typing import Callable, Iterable, Sequence

from .geometry import (
    add,
    axis_from_name,
    basis_from_direction,
    centroid,
    distance,
    norm,
    padded_box,
    point_on_axis,
    project_onto_axis,
    scale,
    sub,
)
from .models import Atom, Cavity, ChannelPoint, Coord, PoreProfile, ResidueContact, StructureFrame, Tunnel
from .grid import GridConnectivity
from .radii import RadiusSet, atom_vdw_radius, radii_option, residue_properties
from .spatial import SpatialIndex


@radii_option
def pore_profile(
    frame: StructureFrame,
    *,
    axis: str | Coord = "auto",
    origin: Coord | None = None,
    samples: int = 81,
    padding: float = 0.0,
    search_radius: float = 6.0,
    refinement_steps: int = 4,
    probe_radius: float = 0.0,
    include_hydrogen: bool = False,
    include_hetero: bool = True,
    section_spacing: float = 0.5,
    enclosure_radius: float | None = None,
    lateral_exits: bool = False,
    exit_bulk_radius: float = 6.0,
    exit_spacing: float = 0.5,
    radii: RadiusSet | str | None = None,
) -> PoreProfile:
    """Resolve a connected channel and measure inscribed atom-surface radii.

    Automatic mode tests all principal shape directions against enclosed free
    sections and two open ends. A zero search radius explicitly requests the
    legacy fixed-axis clearance scan, without asserting a channel assignment.
    Closed cavities and unresolved/ambiguous channels require a separate region
    selection. ``section_spacing`` controls transverse discovery resolution.
    Where a plane holds both the lumen and a small enclosed wall pocket, the
    widest connected path (largest bottleneck) is followed; see
    :func:`crevice.sections.connected_profile` for the full selection rule.

    ``lateral_exits=True`` (off by default) accepts a capped end, whose
    straight axial exit is blocked, when a lateral exit to bulk exists; every
    distinct exit leg is reported separately in ``metadata["exits"]`` and
    ``metadata["path_type"]``. Not available for the ``search_radius=0``
    fixed-axis scan.

    ``enclosure_radius`` (Å) is the probe that decides lateral enclosure and
    connectivity. ``None`` (default) chooses it automatically: the smallest
    probe of the channel that resolves at the most probes of a 0.5-3.0 Å
    ladder, recorded in ``metadata["enclosure_probe_selection"]`` (see
    :func:`crevice.sections.connected_profile`). A number runs exactly that
    probe, as before. The ``search_radius=0`` scan does not use it.

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    axis : {"auto", "x", "y", "z"} or Coord, default "auto"
        Channel direction.
    origin : Coord, optional
        Seed inside the intended channel, Å; default the atom centroid.
    samples : int, default 81
        Number of samples (a minimum in connected mode, where axial steps are
        kept at or below 0.75 Å); at least 2.
    padding : float, default 0.0
        Axis extension beyond the atoms, Å; only with ``search_radius=0``.
    search_radius : float, default 6.0
        Transverse search half-width, Å; 0 = unvalidated fixed-axis scan.
    refinement_steps : int, default 4
        Centre and mouth refinement iterations.
    probe_radius : float, default 0.0
        Measurement probe radius, Å, subtracted from each clearance.
    include_hydrogen : bool, default False
        Treat hydrogen atoms as obstacles.
    include_hetero : bool, default True
        Treat HETATM atoms as obstacles.
    section_spacing : float, default 0.5
        Transverse grid spacing, Å.
    enclosure_radius : float or None, default None
        Enclosure probe, Å; ``None`` = automatic (above).
    lateral_exits : bool, default False
        Accept capped ends with a lateral exit (above).
    exit_bulk_radius : float, default 6.0
        Rolling bulk-probe radius, Å, that defines bulk for lateral exits.
    exit_spacing : float, default 0.5
        Lattice spacing, Å, of the lateral-exit search.
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
        See :func:`crevice.channels.pore_profile`.

    Raises
    ------
    crevice.sections.ChannelResolutionError
        When no unique connected channel is resolved (connected mode).
    ValueError
        For invalid settings.
    """

    if samples < 2:
        raise ValueError("samples must be at least 2")
    if not math.isfinite(probe_radius) or probe_radius < 0:
        raise ValueError("probe_radius must be finite and non-negative")
    if not math.isfinite(search_radius) or search_radius < 0:
        raise ValueError("search_radius must be finite and non-negative")
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    if origin is not None and (len(origin) != 3 or not all(math.isfinite(v) for v in origin)):
        raise ValueError("origin must be a finite 3D coordinate")
    if not math.isfinite(section_spacing) or section_spacing <= 0:
        raise ValueError("section_spacing must be finite and positive")
    if not math.isfinite(padding) or padding < 0:
        raise ValueError("padding must be finite and non-negative")
    if enclosure_radius is not None and (not math.isfinite(enclosure_radius) or enclosure_radius < 0):
        raise ValueError("enclosure_radius must be None (automatic) or finite and non-negative")
    if lateral_exits and search_radius == 0:
        raise ValueError("lateral_exits requires a connected channel search (search_radius > 0)")
    if search_radius > 0:
        if padding:
            raise ValueError("Connected channel profiles end at their mouths; padding is only supported for search_radius=0 scans")
        from .sections import connected_profile
        profile = connected_profile(atoms, axis=axis, origin=origin, samples=samples,
                                    search_radius=search_radius, refinement_steps=refinement_steps,
                                    probe_radius=probe_radius, section_spacing=section_spacing,
                                    enclosure_radius=enclosure_radius, lateral_exits=lateral_exits,
                                    exit_bulk_radius=exit_bulk_radius, exit_spacing=exit_spacing)
        return replace(profile, metadata={**profile.metadata, "include_hydrogen": include_hydrogen,
                                          "include_hetero": include_hetero, "padding": padding})
    inferred_origin, direction, axis_label = axis_from_name(axis, atoms)
    axis_origin = origin if origin is not None else inferred_origin
    projections = [project_onto_axis(atom.coord, axis_origin, direction) for atom in atoms]
    start = min(projections) - padding
    end = max(projections) + padding
    step = (end - start) / (samples - 1)
    spatial = SpatialIndex(atoms)
    points: list[ChannelPoint] = []
    for index in range(samples):
        t = start + index * step
        anchor = point_on_axis(axis_origin, direction, t)
        if search_radius > 0 and refinement_steps > 0:
            center = _refine_center(anchor, direction, spatial, search_radius, refinement_steps)
        else:
            center = anchor
        raw_clearance, nearest = spatial.nearest_surface(center)
        radius = max(0.0, raw_clearance - probe_radius)
        points.append(
            ChannelPoint(
                index=index,
                position=center,
                radius=radius,
                raw_clearance=raw_clearance,
                t=t,
                nearest_atom_serial=nearest.serial,
                nearest_residue=nearest.residue_key.label,
            )
        )
    return PoreProfile(
        axis_origin=axis_origin,
        axis_direction=direction,
        points=tuple(points),
        probe_radius=probe_radius,
        metadata={
            "axis": axis_label,
            "status": "unvalidated_axis_scan",
            "selection": "fixed axis; no channel or inter-sample connectivity assertion",
            "samples": samples,
            "padding": padding,
            "search_radius": search_radius,
            "refinement_steps": refinement_steps,
            "include_hydrogen": include_hydrogen,
            "include_hetero": include_hetero,
        },
    )


@radii_option
def detect_cavities(
    frame: StructureFrame,
    *,
    spacing: float = 2.0,
    min_radius: float = 1.4,
    max_cavities: int = 10,
    shell_radius: float = 8.0,
    min_octants: int = 4,
    max_grid_points: int = 75_000,
    include_hydrogen: bool = False,
    include_hetero: bool = False,
    enclosed_only: bool = True,
    focus_points: Sequence[Coord] | None = None,
    focus_radius: float = 8.0,
    min_component_volume: float = 0.0,
    max_component_volume: float | None = None,
    radii: RadiusSet | str | None = None,
) -> tuple[Cavity, ...]:
    """Detect enclosed grid voids; optionally include local pocket candidates.

    Enclosure is evaluated by exterior flood fill at the specified resolution
    and probe clearance. Open tunnels are intentionally separate from cavities.
    Runs :func:`crevice.voids.void_cast` (``mode="cavity"``, or ``"all"`` with
    ``enclosed_only=False``) and converts each component to a
    :class:`~crevice.models.Cavity`.

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    spacing : float, default 2.0
        Requested grid spacing, Å (coarsened if the grid would exceed
        ``max_grid_points``).
    min_radius : float, default 1.4
        Required atom-surface clearance, Å, of every grid point and edge.
    max_cavities : int, default 10
        Maximum number of cavities, largest first; at least 1.
    shell_radius : float, default 8.0
        Neighbourhood radius, Å, of the octant-burial test and residue listing.
    min_octants : int, default 4
        With ``enclosed_only=False``: octants around a point that must hold atom
        centres within ``shell_radius``.
    max_grid_points : int, default 75000
        Grid-size limit that triggers coarsening.
    include_hydrogen : bool, default False
        Treat hydrogen atoms as obstacles.
    include_hetero : bool, default False
        Treat HETATM atoms as obstacles.
    enclosed_only : bool, default True
        Only components the exterior flood fill cannot reach.
    focus_points : sequence of Coord, optional
        Keep only grid points within ``focus_radius`` of one of these, Å.
    focus_radius : float, default 8.0
        Focus radius, Å.
    min_component_volume : float, default 0.0
        Smallest kept component, Å³ (inclusive).
    max_component_volume : float, optional
        Largest kept component, Å³ (inclusive).
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
        Largest first; see :func:`crevice.cavities.detect_cavities`.

    Raises
    ------
    ValueError
        If ``max_cavities < 1`` or a grid setting is invalid.
    """

    from .voids import void_cast

    if max_cavities < 1:
        raise ValueError("max_cavities must be positive")
    cast = void_cast(frame, mode="cavity" if enclosed_only else "all", spacing=spacing,
                     min_radius=min_radius, max_components=max_cavities,
                     shell_radius=shell_radius, min_octants=0 if enclosed_only else min_octants,
                     max_grid_points=max_grid_points, include_hydrogen=include_hydrogen,
                     include_hetero=include_hetero, focus_points=focus_points, focus_radius=focus_radius,
                     min_component_volume=min_component_volume, max_component_volume=max_component_volume)
    return tuple(Cavity(c.id, c.center, c.max_radius, c.volume, c.point_count,
                        c.nearest_residues, c.score) for c in cast.components)



@radii_option
def find_tunnels(
    frame: StructureFrame,
    *,
    start: Coord | None = None,
    spacing: float = 2.0,
    min_radius: float = 0.8,
    max_tunnels: int = 3,
    padding: float = 3.0,
    max_grid_points: int = 120_000,
    include_hydrogen: bool = False,
    include_hetero: bool = False,
    radii: RadiusSet | str | None = None,
) -> tuple[Tunnel, ...]:
    """Find segment-clearance-weighted grid paths to the grid boundary.

    Attach the actual start to its nearest clearance-reachable grid node within
    one cell diagonal. Absence of a path may reflect insufficient resolution.
    Each point's ``t`` is the cumulative path length in Å from the start
    (:func:`crevice.tunnels.with_path_length`) and ``index`` is its step
    number, the same as :func:`crevice.tunnels.find_tunnels`, which is the
    documented entry point (algorithm there).

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    start : Coord, optional
        Start point, Å; default the centroid of the selected atoms.
    spacing : float, default 2.0
        Requested grid spacing, Å (coarsened if the grid would exceed
        ``max_grid_points``; recorded as ``spacing_adjusted``).
    min_radius : float, default 0.8
        Required clearance, Å, of every node and edge (a threshold; reported radii
        are not reduced by it).
    max_tunnels : int, default 3
        Maximum number of tunnels, widest first.
    padding : float, default 3.0
        Margin around the atom bounding box, Å; the box faces stand in for bulk.
    max_grid_points : int, default 120000
        Grid-size limit that triggers coarsening.
    include_hydrogen : bool, default False
        Treat hydrogen atoms as obstacles.
    include_hetero : bool, default False
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
    tuple of Tunnel
        Widest first; empty when the default start is too tight or no path
        reaches the boundary.

    Raises
    ------
    ValueError
        For an explicit start with too little clearance or invalid settings.
    """

    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be finite and positive")
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    base_min, base_max = frame.bounding_box(atoms)
    min_corner, max_corner = padded_box(base_min, base_max, padding)
    requested_spacing = spacing
    spacing = _scale_spacing_for_box(min_corner, max_corner, spacing, max_grid_points)
    dims = _grid_dimensions(min_corner, max_corner, spacing)
    spatial = SpatialIndex(atoms)
    origin = tuple(float(v) for v in start) if start is not None else frame.centroid(atoms)
    if min_radius < 0 or not math.isfinite(min_radius):
        raise ValueError("min_radius must be finite and non-negative")
    if len(origin) != 3 or not all(math.isfinite(v) for v in origin):
        raise ValueError("Tunnel start must be a finite 3D coordinate")
    origin_clearance, origin_nearest = spatial.nearest_surface(origin)
    if origin_clearance < min_radius:
        if start is not None:
            raise ValueError("Tunnel start does not have the required clearance")
        return ()

    nodes: dict[tuple[int, int, int], float] = {}
    nearest_atom_by_node: dict[tuple[int, int, int], Atom] = {}
    for idx in itertools.product(range(dims[0]), range(dims[1]), range(dims[2])):
        point = _grid_point(min_corner, spacing, idx)
        clearance, nearest = spatial.nearest_surface(point)
        if clearance >= min_radius:
            nodes[idx] = clearance
            nearest_atom_by_node[idx] = nearest
    if not nodes:
        return ()

    connectivity = GridConnectivity(spatial, min_corner, spacing, min_radius)
    nearby = sorted((distance(connectivity.position(idx), origin), idx) for idx in nodes
                    if distance(connectivity.position(idx), origin) <= math.sqrt(3) * spacing)
    start_idx = None
    for attachment_length, idx in nearby:
        attachment_clearance = spatial.segment_clearance(origin, connectivity.position(idx))
        if attachment_clearance >= min_radius:
            start_idx = idx
            break
    if start_idx is None:
        return ()
    bottleneck, lengths, previous = _widest_paths(
        nodes, dims, start_idx, spacing, edge_clearance=connectivity.clearance,
        min_radius=min_radius, initial_capacity=attachment_clearance,
        initial_length=attachment_length)
    boundary_nodes = [
        idx
        for idx in nodes
        if idx in bottleneck and _is_boundary(idx, dims)
    ]
    boundary_nodes.sort(key=lambda idx: (bottleneck[idx], -lengths[idx]), reverse=True)

    tunnels: list[Tunnel] = []
    selected_ends: list[Coord] = []
    for idx in boundary_nodes:
        end = _grid_point(min_corner, spacing, idx)
        if any(distance(end, other) < spacing * 3 for other in selected_ends):
            continue
        path_indices = _reconstruct_path(previous, idx)
        points = _channel_points_from_grid(path_indices, min_corner, spacing, nodes, nearest_atom_by_node)
        if not points:
            continue
        if origin != points[0].position:
            seed = ChannelPoint(0, origin, origin_clearance, origin_clearance, 0.0,
                                origin_nearest.serial, origin_nearest.residue_key.label)
            points = (seed,) + tuple(replace(p, index=i + 1, t=float(i + 1)) for i, p in enumerate(points))
        length = _path_length([point.position for point in points])
        mean_radius = sum(point.radius for point in points) / len(points)
        segment_clearances = [connectivity.clearance(a, b) for a, b in zip(path_indices, path_indices[1:])]
        if origin != connectivity.position(start_idx):
            segment_clearances.insert(0, attachment_clearance)
        bottleneck_radius = min([point.radius for point in points] + segment_clearances)
        residues = _path_nearest_residues(frame, points, limit=10)
        score = bottleneck_radius / (1.0 + 0.01 * length)
        tunnels.append(
            Tunnel(
                id=len(tunnels) + 1,
                start=points[0].position,
                end=end,
                points=points,
                length=length,
                bottleneck_radius=bottleneck_radius,
                mean_radius=mean_radius,
                score=score,
                nearest_residues=residues,
                metadata={"spacing_A": spacing, "requested_spacing_A": requested_spacing,
                          "spacing_adjusted": spacing != requested_spacing, "grid_origin": min_corner,
                          "grid_dimensions": dims, "required_clearance_A": min_radius,
                          "connectivity": connectivity.to_dict(),
                          "seed_grid_position": connectivity.position(start_idx),
                          "seed_attachment_length_A": attachment_length,
                          "seed_attachment_clearance_A": attachment_clearance,
                          "segment_clearances_A": segment_clearances,
                          "bottleneck_estimator": "minimum continuous segment clearance including seed",
                          "mean_radius_estimator": "arithmetic mean at returned nodes"},
            )
        )
        selected_ends.append(end)
        if len(tunnels) >= max_tunnels:
            break
    from .tunnels import with_path_length  # deferred: crevice.tunnels imports this module
    return tuple(with_path_length(tunnel) for tunnel in tunnels)


@radii_option
def annotate_residues(
    frame: StructureFrame,
    points: Sequence[ChannelPoint | Coord],
    *,
    cutoff: float = 4.5,
    bottleneck_index: int | None = None,
    identify_bottleneck: bool = True,
    include_hydrogen: bool = False,
    include_hetero: bool = True,
    radii: RadiusSet | str | None = None,
) -> tuple[ResidueContact, ...]:
    """Annotate atomic-surface distance to a sampled region surface.

    ChannelPoint raw clearance defines a local region sphere. Bare coordinates
    retain point-to-atom-surface semantics. Scores describe geometric contacts.

    For each residue, the gap between its atom surfaces and the region surface
    (the sphere of each sample's raw clearance) is measured over every sample
    within ``cutoff``. The role is ``bottleneck`` (gap <= 1.0 Å and contacting a
    sample within two indices of the bottleneck), ``lining`` (gap <= 1.0 Å),
    ``bottleneck-nearby`` (within ``cutoff`` near the bottleneck) or ``nearby``.
    The dimensionless influence score is ``1/(0.25 + gap) + 0.15*ln(1 + n)``
    (n = samples contacted) plus 1.0 for bottleneck roles or 0.25 for lining.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the region.
    points : sequence of ChannelPoint or Coord
        Region samples (for example profile points); bare coordinates are
        treated as zero-radius points.
    cutoff : float, default 4.5
        Largest surface gap, Å, at which a residue is reported.
    bottleneck_index : int, optional
        Index of the bottleneck sample; default the sample of smallest radius
        when ``identify_bottleneck`` and the points are ChannelPoints.
    identify_bottleneck : bool, default True
        Locate the bottleneck automatically when ``bottleneck_index`` is not
        given.
    include_hydrogen : bool, default False
        Include hydrogen atoms.
    include_hetero : bool, default True
        Include HETATM atoms (ligands, waters, ions) as contacting residues.
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
    tuple of ResidueContact
        Highest influence score first (ties: smaller gap first); ``min_distance``
        and ``mean_distance`` are gaps in Å.

    Raises
    ------
    ValueError
        If ``cutoff`` is negative or not finite.
    """

    if not points:
        return ()
    if not math.isfinite(cutoff) or cutoff < 0:
        raise ValueError("cutoff must be finite and non-negative")
    coords = tuple(point.position if isinstance(point, ChannelPoint) else point for point in points)
    region_radii = tuple(max(0.0, point.raw_clearance) if isinstance(point, ChannelPoint) else 0.0
                         for point in points)
    if identify_bottleneck and bottleneck_index is None and isinstance(points[0], ChannelPoint):
        bottleneck_index = min(range(len(points)), key=lambda i: points[i].radius)
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    residues = frame.residues()
    atom_selection = set(atoms)
    spatial = SpatialIndex(atoms)
    hits_by_residue = {}
    for point_index, (coord, region_radius) in enumerate(zip(coords, region_radii)):
        for atom, clearance in spatial.atoms_within_surface_distance(coord, cutoff + region_radius):
            distances, point_hits = hits_by_residue.setdefault(atom.residue_key, ([], set()))
            distances.append(max(0.0, clearance - region_radius))
            point_hits.add(point_index)
    contacts: list[ResidueContact] = []
    for residue in residues:
        residue_atoms = tuple(atom for atom in residue.atoms if atom in atom_selection)
        if not residue_atoms:
            continue
        distances, point_hits = hits_by_residue.get(residue.key, ([], set()))
        if not distances:
            continue
        min_distance = min(distances)
        mean_distance = sum(distances) / len(distances)
        role = _contact_role(min_distance, point_hits, bottleneck_index)
        properties = residue_properties(residue.key.resname)
        influence = _influence_score(min_distance, point_hits, role)
        contacts.append(
            ResidueContact(
                residue=residue.label,
                min_distance=min_distance,
                mean_distance=mean_distance,
                atom_count=len(residue_atoms),
                point_count=len(point_hits),
                role=role,
                properties=properties,
                influence_score=influence,
            )
        )
    return tuple(sorted(contacts, key=lambda contact: (-contact.influence_score, contact.min_distance)))


@radii_option
def nearest_residues(
    frame: StructureFrame,
    point: Coord,
    *,
    limit: int = 8,
    cutoff: float = 8.0,
    radii: RadiusSet | str | None = None,
) -> tuple[str, ...]:
    """Residues with atom surfaces nearest to a point, closest first.

    Parameters
    ----------
    frame : StructureFrame
    point : Coord
        Query point (Å).
    limit : int, default 8
        Maximum number of residues returned.
    cutoff : float, default 8.0
        Surface-distance cutoff in Å.
    radii : RadiusSet, str or None, optional
        Atomic radius set (see :func:`crevice.radii.use_radii`); ``None`` uses
        the set in effect, by default the standard table.

    Returns
    -------
    tuple of str
        Residue labels.
    """
    contacts = annotate_residues(frame, [point], cutoff=cutoff)
    return tuple(contact.residue for contact in contacts[:limit])


def _refine_center(
    anchor: Coord,
    direction: Coord,
    spatial: SpatialIndex,
    search_radius: float,
    refinement_steps: int,
) -> Coord:
    u, v = basis_from_direction(direction)
    best = anchor
    best_clearance, _ = spatial.nearest_surface(best)
    step = search_radius / 2.0
    for _ in range(refinement_steps):
        improved = best
        improved_score = best_clearance - 0.02 * distance(best, anchor)
        for du in (-step, 0.0, step):
            for dv in (-step, 0.0, step):
                candidate = add(best, add(scale(u, du), scale(v, dv)))
                offset = distance(candidate, anchor)
                if offset > search_radius:
                    continue
                clearance, _ = spatial.nearest_surface(candidate)
                score = clearance - 0.02 * offset
                if score > improved_score:
                    improved = candidate
                    improved_score = score
                    best_clearance = clearance
        best = improved
        step *= 0.5
    return best


def _scale_spacing_for_box(min_corner: Coord, max_corner: Coord, spacing: float, max_grid_points: int) -> float:
    dims = _grid_dimensions(min_corner, max_corner, spacing)
    count = dims[0] * dims[1] * dims[2]
    if count <= max_grid_points:
        return spacing
    factor = (count / max_grid_points) ** (1.0 / 3.0)
    return spacing * factor


def _grid_dimensions(min_corner: Coord, max_corner: Coord, spacing: float) -> tuple[int, int, int]:
    return (
        max(2, int(math.floor((max_corner[0] - min_corner[0]) / spacing)) + 1),
        max(2, int(math.floor((max_corner[1] - min_corner[1]) / spacing)) + 1),
        max(2, int(math.floor((max_corner[2] - min_corner[2]) / spacing)) + 1),
    )


def _grid_point(min_corner: Coord, spacing: float, idx: tuple[int, int, int]) -> Coord:
    return (
        min_corner[0] + idx[0] * spacing,
        min_corner[1] + idx[1] * spacing,
        min_corner[2] + idx[2] * spacing,
    )


def _directional_coverage(point: Coord, atoms: tuple[Atom, ...], shell_radius: float) -> int:
    octants: set[tuple[int, int, int]] = set()
    for atom in atoms:
        vector = sub(atom.coord, point)
        if norm(vector) > shell_radius:
            continue
        octants.add(
            (
                1 if vector[0] >= 0 else -1,
                1 if vector[1] >= 0 else -1,
                1 if vector[2] >= 0 else -1,
            )
        )
    return len(octants)


def _cluster_grid(candidates: dict[tuple[int, int, int], tuple[Coord, float]]) -> list[list[tuple[int, int, int]]]:
    unseen = set(candidates)
    clusters: list[list[tuple[int, int, int]]] = []
    neighbor_offsets = [
        offset for offset in itertools.product((-1, 0, 1), repeat=3) if offset != (0, 0, 0)
    ]
    while unseen:
        seed = unseen.pop()
        cluster = [seed]
        queue: deque[tuple[int, int, int]] = deque([seed])
        while queue:
            idx = queue.popleft()
            for offset in neighbor_offsets:
                neighbor = (idx[0] + offset[0], idx[1] + offset[1], idx[2] + offset[2])
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    cluster.append(neighbor)
                    queue.append(neighbor)
        clusters.append(cluster)
    clusters.sort(key=len, reverse=True)
    return clusters


def _weighted_centroid(points: Sequence[Coord], weights: Sequence[float]) -> Coord:
    total = sum(weights)
    if total == 0:
        return centroid(points)
    return (
        sum(point[0] * weight for point, weight in zip(points, weights)) / total,
        sum(point[1] * weight for point, weight in zip(points, weights)) / total,
        sum(point[2] * weight for point, weight in zip(points, weights)) / total,
    )


def _widest_paths(
    nodes: dict[tuple[int, int, int], float],
    dims: tuple[int, int, int],
    start_idx: tuple[int, int, int],
    spacing: float,
    *,
    edge_clearance: Callable[[tuple[int, int, int], tuple[int, int, int]], float],
    min_radius: float,
    initial_capacity: float | None = None,
    initial_length: float = 0.0,
) -> tuple[dict[tuple[int, int, int], float], dict[tuple[int, int, int], float], dict[tuple[int, int, int], tuple[int, int, int] | None]]:
    neighbor_offsets = [
        offset for offset in itertools.product((-1, 0, 1), repeat=3) if offset != (0, 0, 0)
    ]
    capacity = nodes[start_idx] if initial_capacity is None else min(nodes[start_idx], initial_capacity)
    bottleneck = {start_idx: capacity}
    lengths = {start_idx: initial_length}
    previous: dict[tuple[int, int, int], tuple[int, int, int] | None] = {start_idx: None}
    heap: list[tuple[float, float, tuple[int, int, int]]] = [(-capacity, initial_length, start_idx)]
    while heap:
        neg_width, length, idx = heapq.heappop(heap)
        width = -neg_width
        if width < bottleneck.get(idx, -1.0):
            continue
        if length > lengths.get(idx, float("inf")) and abs(width - bottleneck[idx]) < 1e-9:
            continue
        for offset in neighbor_offsets:
            neighbor = (idx[0] + offset[0], idx[1] + offset[1], idx[2] + offset[2])
            if neighbor not in nodes or not _in_grid(neighbor, dims):
                continue
            move_length = spacing * math.sqrt(offset[0] ** 2 + offset[1] ** 2 + offset[2] ** 2)
            edge_width = edge_clearance(idx, neighbor)
            if edge_width < min_radius:
                continue
            new_width = min(width, nodes[neighbor], edge_width)
            new_length = length + move_length
            old_width = bottleneck.get(neighbor, -1.0)
            old_length = lengths.get(neighbor, float("inf"))
            if new_width > old_width + 1e-9 or (abs(new_width - old_width) < 1e-9 and new_length < old_length):
                bottleneck[neighbor] = new_width
                lengths[neighbor] = new_length
                previous[neighbor] = idx
                heapq.heappush(heap, (-new_width, new_length, neighbor))
    return bottleneck, lengths, previous


def _in_grid(idx: tuple[int, int, int], dims: tuple[int, int, int]) -> bool:
    return 0 <= idx[0] < dims[0] and 0 <= idx[1] < dims[1] and 0 <= idx[2] < dims[2]


def _is_boundary(idx: tuple[int, int, int], dims: tuple[int, int, int]) -> bool:
    return (
        idx[0] in {0, dims[0] - 1}
        or idx[1] in {0, dims[1] - 1}
        or idx[2] in {0, dims[2] - 1}
    )


def _reconstruct_path(
    previous: dict[tuple[int, int, int], tuple[int, int, int] | None],
    end: tuple[int, int, int],
) -> list[tuple[int, int, int]]:
    path = [end]
    current = end
    while previous.get(current) is not None:
        current = previous[current]  # type: ignore[assignment]
        path.append(current)
    path.reverse()
    return path


def _channel_points_from_grid(
    path_indices: Sequence[tuple[int, int, int]],
    min_corner: Coord,
    spacing: float,
    nodes: dict[tuple[int, int, int], float],
    nearest_atom_by_node: dict[tuple[int, int, int], Atom],
) -> tuple[ChannelPoint, ...]:
    points: list[ChannelPoint] = []
    for index, idx in enumerate(path_indices):
        position = _grid_point(min_corner, spacing, idx)
        nearest = nearest_atom_by_node[idx]
        points.append(
            ChannelPoint(
                index=index,
                position=position,
                radius=nodes[idx],
                raw_clearance=nodes[idx],
                t=float(index),
                nearest_atom_serial=nearest.serial,
                nearest_residue=nearest.residue_key.label,
            )
        )
    return tuple(points)


def _path_length(points: Sequence[Coord]) -> float:
    if len(points) < 2:
        return 0.0
    return sum(distance(left, right) for left, right in zip(points, points[1:]))


def _path_nearest_residues(frame: StructureFrame, points: Sequence[ChannelPoint], limit: int = 10) -> tuple[str, ...]:
    contacts = annotate_residues(frame, points, cutoff=5.0)
    return tuple(contact.residue for contact in contacts[:limit])


def _contact_role(min_distance: float, point_hits: set[int], bottleneck_index: int | None) -> str:
    near_bottleneck = bottleneck_index is not None and any(abs(index - bottleneck_index) <= 2 for index in point_hits)
    if min_distance <= 1.0 and near_bottleneck:
        return "bottleneck"
    if min_distance <= 1.0:
        return "lining"
    if near_bottleneck:
        return "bottleneck-nearby"
    return "nearby"


def _influence_score(min_distance: float, point_hits: set[int], role: str) -> float:
    contact_strength = 1.0 / (0.25 + max(min_distance, 0.0))
    span_strength = math.log1p(len(point_hits))
    role_bonus = 1.0 if "bottleneck" in role else 0.25 if role == "lining" else 0.0
    return contact_strength + 0.15 * span_strength + role_bonus

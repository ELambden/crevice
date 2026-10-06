"""Clearance-grid casts: fill cavities, pockets and channels with atom-clear points.

A *cast* fills empty space with grid points so it can be measured (volume)
and drawn (as dummy atoms or a surface). Every exported point satisfies a
hard physical bound: its distance to the nearest atom van der Waals surface
(its *clearance*) is at least ``min_radius`` Å. Two neighbouring points
belong to the same component only if a probe of radius ``min_radius`` can
slide along the whole segment between them (checked analytically by
:class:`crevice.grid.GridConnectivity`).

Method in brief
    1. Build a grid over the atom bounding box (plus padding) and keep nodes
       with clearance ``>= min_radius``.
    2. Cavity mode removes every node reachable from the box boundary by a
       flood fill; the rest is buried. Other modes instead require the
       node to be surrounded by atoms in enough octants (``min_octants``).
    3. Group nodes into connected components (26-neighbour, segment-checked)
       and classify each by which box faces it touches:
       ``buried_cavity`` (none), ``surface_connected`` (some) or
       ``through_channel`` (two opposite faces).
    4. Channel mode keeps only the part of the open components near the
       channel centre line; focus points can crop further.
    5. Re-group the surviving points, filter by volume and keep the largest
       ``max_components``.

Main entry points: :func:`void_cast`, :func:`validate_void_cast`,
:func:`component_profile`. For casts that include open pockets bounded by a
rolling outer probe, see :func:`crevice.rolling.rolling_probe_cast`
(command ``crevice cast``).

Assumptions and limitations
    Volumes are ``point count * spacing**3`` and depend on spacing and
    ``min_radius``. Failing to find a grid path is not proof that no
    continuous path exists, and "buried" means enclosed at this resolution
    and probe size. The component kinds are geometric labels; they do not
    assign biological channels or water occupancy.
"""

from __future__ import annotations

import itertools
import math
import operator
from bisect import bisect_left
from collections import deque
from dataclasses import dataclass, replace
from typing import Any, Callable, Sequence

from .geometry import axis_from_name, centroid, distance, norm, project_onto_axis, sub
from .grid import GridConnectivity
from .models import Atom, ChannelPoint, Coord, PoreProfile, StructureFrame, VoidCast, VoidComponent
from .spatial import SpatialIndex
from .radii import RadiusSet, radii_option

GridIndex = tuple[int, int, int]

_NEIGHBOR_OFFSETS: tuple[GridIndex, ...] = tuple(
    offset for offset in itertools.product((-1, 0, 1), repeat=3) if offset != (0, 0, 0)
)
_FACE_NAMES = (
    (0, 0, "x-"),
    (0, -1, "x+"),
    (1, 0, "y-"),
    (1, -1, "y+"),
    (2, 0, "z-"),
    (2, -1, "z+"),
)


@dataclass(frozen=True)
class _GridNode:
    """Candidate grid node: position (Å), clearance (Å), limiting atom and the
    number of atom-occupied octants within the shell radius.
    """
    position: Coord
    clearance: float
    nearest: Atom
    coverage: int


@radii_option
def void_cast(
    frame: StructureFrame,
    *,
    mode: str = "auto",
    axis: str | Coord = "auto",
    spacing: float = 1.5,
    min_radius: float = 0.4,
    max_components: int = 4,
    padding: float = 0.0,
    shell_radius: float = 9.0,
    min_octants: int | None = None,
    min_component_points: int = 4,
    max_grid_points: int = 120_000,
    max_export_points: int = 9_000,
    profile: PoreProfile | None = None,
    centerline_samples: int | None = None,
    centerline_search_radius: float = 6.0,
    channel_radius_scale: float = 1.0,
    channel_radius_padding: float | None = None,
    channel_max_radius: float | None = None,
    focus_points: Sequence[Coord] | None = None,
    focus_radius: float = 8.0,
    min_component_volume: float = 0.0,
    max_component_volume: float | None = None,
    include_hydrogen: bool | None = None,
    include_hetero: bool | None = None,
    radii: RadiusSet | str | None = None,
) -> VoidCast:
    """Fill connected, atom-clear empty space and return it as a cast.

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    mode : {"auto", "channel", "cavity", "pocket", "all"}, default "auto"
        What to fill (case-insensitive; ``"pocket"`` is an alias of
        ``"cavity"``).

        * ``"cavity"``: only buried components (unreachable from the box edge).
        * ``"channel"``: the open space along a channel centre line. The centre
          line comes from ``profile`` or from
          :func:`crevice.channels.pore_profile`. When that profile is a
          connected-channel profile (the default), the cast is delegated to
          :func:`crevice.sections.channel_cast`, which fills the probe-swept
          channel sections instead of the grid procedure described here. The
          grid centre-line fill is used for other profile methods, for example
          ``centerline_search_radius=0`` or a supplied explicit-path profile.
        * ``"auto"``: like ``"channel"``, but if no channel can be resolved
          (and the failure is not an ambiguity) it falls back to ``"cavity"``
          and records the reason in ``metadata["auto_channel_unresolved"]``.
        * ``"all"``: every free component that passes the octant test.
    axis : {"auto", "x", "y", "z"} or Coord, default "auto"
        Axis for the channel profile and for each point's ``t`` coordinate.
    spacing : float, default 1.5
        Requested grid spacing in Å. Enlarged automatically if the grid would
        exceed ``max_grid_points`` (see ``metadata["spacing_adjusted"]``).
    min_radius : float, default 0.4
        Required atom-surface clearance in Å for every point and edge (the
        probe radius of the cast).
    max_components : int, default 4
        Maximum number of final components returned (largest volume first).
    padding : float, default 0.0
        Å added around the atom bounding box. In cavity mode a value of 0 is
        replaced by ``max(3.0, 2 * spacing)`` so the exterior flood fill can
        surround the molecule.
    shell_radius : float, default 9.0
        Radius in Å of the neighbourhood used for the octant test and for
        nearby-residue lists.
    min_octants : int, optional
        Minimum number of the eight octants around a node that must contain an
        atom centre within ``shell_radius``. Defaults to 0 in cavity mode and 3
        otherwise.
    min_component_points : int, default 4
        Discard components with fewer points (raised to at least 8 for
        centre-line channel fills).
    max_grid_points : int, default 120000
        Grid-size limit that triggers automatic coarsening.
    max_export_points : int, default 9000
        Maximum number of display points in :attr:`~crevice.models.VoidCast.points`. Thinning
        is even along ``t`` and always keeps the narrowest point; it does not
        change volumes or selection.
    profile : PoreProfile, optional
        Pre-computed centre line for channel/auto modes. Its recorded
        ``include_hydrogen``/``include_hetero`` settings become the defaults.
    centerline_samples : int, optional
        ``samples`` for the automatically computed profile (81 if not given).
    centerline_search_radius : float, default 6.0
        ``search_radius`` for the automatically computed profile.
    channel_radius_scale : float, default 1.0
        Grid centre-line fills keep points within
        ``local_radius * channel_radius_scale + channel_radius_padding`` of the
        nearest centre-line point. Must stay 1.0 for connected profiles.
    channel_radius_padding : float, optional
        See ``channel_radius_scale`` (Å, default 0).
    channel_max_radius : float, optional
        Upper limit on ``local_radius`` (Å). Defaults to the 85th percentile of
        the profile radii plus one spacing.
    focus_points : sequence of Coord, optional
        Keep only points within ``focus_radius`` of at least one of these
        positions (Å), before components are formed.
    focus_radius : float, default 8.0
    min_component_volume : float, default 0.0
        Discard final components below this volume (Å³, inclusive bound).
    max_component_volume : float, optional
        Discard final components above this volume (Å³, inclusive bound).
    include_hydrogen : bool, optional
        Treat hydrogens as obstacles. Defaults to the profile's setting, else
        ``False``.
    include_hetero : bool, optional
        Treat HETATM atoms as obstacles. Defaults to the profile's setting,
        else ``False``.
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
    VoidCast
        ``components`` are the retained regions with every measurement point;
        ``points`` are the (possibly thinned) display points, whose ``radius``
        and ``raw_clearance`` are both the point's atom-surface clearance (Å).
        ``metadata`` records the grid, the selection order, every candidate's
        keep/reject decision, the enclosure test and the connectivity rule.

    Raises
    ------
    ValueError
        For invalid parameter values, for sphere-scaling options combined with
        a connected profile, or when a channel is required but ambiguous.
    crevice.sections.ChannelResolutionError
        In ``"channel"`` mode when no channel can be resolved.

    See Also
    --------
    validate_void_cast : Re-check the clearance bound of a cast.
    crevice.cavities.detect_cavities : Cavity summaries built on this function.
    crevice.rolling.rolling_probe_cast : Rolling-probe casts including pockets.

    Notes
    -----
    Selection order: source-mode eligibility, centre-line crop, focus crop,
    final connected components, component point threshold, inclusive volume
    bounds, descending volume (ties: lexicographically smallest point
    coordinate), component limit, then display thinning. Because the limit comes
    last, a rejected large source cannot hide a smaller eligible one.

    Examples
    --------
    Fill the interior of a closed spherical shell (radius 7 Å):

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.voids import void_cast, validate_void_cast
    >>> def shell_atom(i, n=250, golden=math.pi * (3 - math.sqrt(5))):
    ...     y = 1 - (2 * i + 1) / n
    ...     r = math.sqrt(1 - y * y)
    ...     return Atom(i + 1, "C", "ALA", "A", i + 1, 7 * math.cos(golden * i) * r,
    ...                 7 * y, 7 * math.sin(golden * i) * r, "C")
    >>> shell = StructureFrame(tuple(shell_atom(i) for i in range(250)))
    >>> cast = void_cast(shell, mode="cavity", spacing=1.0)
    >>> [(c.kind, c.volume) for c in cast.components]
    [('buried_cavity', 496.0)]
    >>> validate_void_cast(shell, cast)["is_physical"]
    True
    """

    if spacing <= 0:
        raise ValueError("spacing must be positive")
    if min_radius < 0:
        raise ValueError("min_radius must be non-negative")
    try:
        if isinstance(max_components, bool):
            raise TypeError("Boolean component limit")
        max_components = operator.index(max_components)
    except TypeError as exc:
        raise ValueError("max_components must be a positive integer") from exc
    if max_components < 1:
        raise ValueError("max_components must be a positive integer")
    if max_export_points < 1:
        raise ValueError("max_export_points must be at least 1")
    for name, value in (("spacing", spacing), ("shell_radius", shell_radius),
                        ("focus_radius", focus_radius), ("channel_radius_scale", channel_radius_scale)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    for name, value in (("min_radius", min_radius), ("padding", padding),
                        ("min_component_volume", min_component_volume),
                        ("channel_radius_padding", channel_radius_padding or 0.0)):
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and non-negative")
    if max_component_volume is not None and (not math.isfinite(max_component_volume)
                                             or max_component_volume < min_component_volume):
        raise ValueError("max_component_volume must be finite and >= min_component_volume")
    if channel_max_radius is not None and (not math.isfinite(channel_max_radius) or channel_max_radius <= 0):
        raise ValueError("channel_max_radius must be finite and positive")
    if max_grid_points < 8 or min_component_points < 1:
        raise ValueError("max_grid_points must be >= 8 and min_component_points >= 1")
    if min_octants is not None and not 0 <= min_octants <= 8:
        raise ValueError("min_octants must be between 0 and 8")
    if focus_points is not None:
        focus_points = tuple(tuple(p) for p in focus_points)
        if not focus_points or any(len(p) != 3 or not all(math.isfinite(v) for v in p) for p in focus_points):
            raise ValueError("focus_points must contain finite 3D coordinates")
    mode = mode.lower().replace("pocket", "cavity")
    if mode not in {"auto", "channel", "cavity", "all"}:
        raise ValueError("mode must be auto, channel, cavity, or all")

    if include_hydrogen is None:
        include_hydrogen = profile.metadata.get("include_hydrogen", False) if profile else False
    if include_hetero is None:
        include_hetero = profile.metadata.get("include_hetero", False) if profile else False
    if mode in {"auto", "channel"} and profile is None:
        from .analysis import pore_profile
        from .sections import ChannelResolutionError
        try:
            # Fixed legacy enclosure probe: void_cast exposes no enclosure option.
            profile = pore_profile(frame, axis=axis, samples=centerline_samples or 81,
                                   search_radius=centerline_search_radius, enclosure_radius=0.8,
                                   include_hydrogen=include_hydrogen, include_hetero=include_hetero)
        except ChannelResolutionError as exc:
            if mode != "auto" or "Ambiguous" in str(exc):
                raise
            cast = void_cast(frame, mode="cavity", axis=axis, spacing=spacing,
                             min_radius=min_radius, max_components=max_components, padding=padding,
                             shell_radius=shell_radius, min_octants=min_octants,
                             min_component_points=min_component_points, max_grid_points=max_grid_points,
                             max_export_points=max_export_points, focus_points=focus_points, focus_radius=focus_radius,
                             min_component_volume=min_component_volume, max_component_volume=max_component_volume,
                             include_hydrogen=include_hydrogen, include_hetero=include_hetero)
            return replace(cast, metadata={**cast.metadata, "auto_channel_unresolved": str(exc)})
    if mode in {"auto", "channel"} and profile is not None and profile.method == "axial-connected":
        if channel_radius_scale != 1.0 or channel_radius_padding or channel_max_radius is not None or padding:
            raise ValueError("Connected section casts use probe-swept boundaries; sphere radius scaling and padding do not apply")
        from .sections import channel_cast
        return channel_cast(frame, profile, spacing=spacing, min_radius=min_radius,
                            max_grid_points=max_grid_points, max_components=max_components,
                            max_export_points=max_export_points, min_component_points=min_component_points,
                            focus_points=focus_points, focus_radius=focus_radius,
                            min_component_volume=min_component_volume, max_component_volume=max_component_volume,
                            include_hydrogen=include_hydrogen, include_hetero=include_hetero)

    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    effective_padding = padding
    if mode == "cavity" and effective_padding == 0.0:
        effective_padding = max(3.0, spacing * 2.0)
    min_corner, max_corner = _padded_box(*frame.bounding_box(atoms), effective_padding)
    requested_spacing = spacing
    spacing = _scale_spacing_for_box(min_corner, max_corner, spacing, max_grid_points)
    dims = _grid_dimensions(min_corner, max_corner, spacing)
    spatial = SpatialIndex(atoms)
    connectivity = GridConnectivity(spatial, min_corner, spacing, min_radius)
    axis_origin, axis_direction, axis_label = axis_from_name(axis, atoms)
    required_octants = min_octants if min_octants is not None else (0 if mode == "cavity" else 3)

    if mode == "cavity":
        candidates = _collect_buried_void_nodes(
            spatial,
            connectivity=connectivity,
            min_corner=min_corner,
            spacing=spacing,
            dims=dims,
            min_radius=min_radius,
            shell_radius=shell_radius,
            min_octants=required_octants,
        )
    else:
        candidates = _collect_void_nodes(
            spatial,
            min_corner=min_corner,
            spacing=spacing,
            dims=dims,
            min_radius=min_radius,
            shell_radius=shell_radius,
            min_octants=required_octants,
        )
    def edge_allowed(left: GridIndex, right: GridIndex) -> bool:
        return connectivity.allows(left, right, candidates[left].clearance, candidates[right].clearance)

    source_clusters = _cluster_indices(candidates, edge_allowed=edge_allowed)
    source_components = _build_components(
        frame,
        source_clusters,
        candidates,
        min_corner=min_corner,
        spacing=spacing,
        dims=dims,
        min_component_points=min_component_points,
        axis_origin=axis_origin,
        axis_direction=axis_direction,
        shell_radius=shell_radius,
    )
    # Keep source topology before cropping. Limit only the measured final pieces,
    # so a rejected or empty high-ranked source cannot hide an eligible region.
    considered_sources = _eligible_source_components(source_components, mode=mode)
    export_mode = _export_mode(mode, considered_sources)
    centerline = None
    if export_mode == "centerline_fill":
        centerline = _ensure_centerline_profile(
            frame,
            profile=profile,
            axis=axis,
            atoms=atoms,
            spacing=spacing,
            samples=centerline_samples,
            search_radius=centerline_search_radius,
            include_hydrogen=include_hydrogen,
            include_hetero=include_hetero,
        )

    export_points, export_metadata = _export_points(
        considered_sources,
        max_points=max_export_points,
        export_mode=export_mode,
        grid_origin=min_corner,
        spacing=spacing,
        dims=dims,
        centerline=centerline,
        min_radius=min_radius,
        channel_radius_scale=channel_radius_scale,
        channel_radius_padding=channel_radius_padding if channel_radius_padding is not None else 0.0,
        channel_max_radius=channel_max_radius,
    )
    if focus_points is not None:
        export_points = tuple(p for p in export_points
                              if any(distance(p.position, focus) <= focus_radius for focus in focus_points))
    cropped_point_count = len(export_points)
    final_min_component_points = max(min_component_points, 8) if export_mode == "centerline_fill" else min_component_points
    final_candidates = _components_from_export_points(
        frame,
        export_points,
        source_components=considered_sources,
        edge_allowed=edge_allowed,
        export_mode=export_mode,
        grid_origin=min_corner,
        spacing=spacing,
        dims=dims,
        shell_radius=shell_radius,
        min_component_points=final_min_component_points,
    )
    cast_components, selection = _select_final_components(
        final_candidates, max_components=max_components,
        min_volume=min_component_volume, max_volume=max_component_volume)
    selection["discarded_small_component_points"] = cropped_point_count - sum(c.point_count for c in final_candidates)
    rejected = [row["candidate_component_id"] for row in selection["candidates"]
                if row["decision"] in {"below_min_volume", "above_max_volume"}]
    selected_source_ids = {c.metadata["source_component_id"] for c in cast_components}
    selected_sources = tuple(c for c in considered_sources if c.id in selected_source_ids)
    full_points = [point for component in cast_components for point in component.points]
    export_points = _renumber_points(_sample_points_by_axis(full_points, max_export_points))
    export_metadata["point_count"] = len(export_points)
    export_metadata["final_component_count"] = len(cast_components)
    export_metadata["final_min_component_points"] = final_min_component_points

    return VoidCast(
        components=cast_components,
        points=export_points,
        spacing=spacing,
        min_radius=min_radius,
        mode=mode,
        metadata={
            "axis": axis_label,
            "axis_origin": axis_origin,
            "axis_direction": axis_direction,
            "grid_origin": min_corner,
            "grid_dimensions": dims,
            "requested_spacing_A": requested_spacing,
            "effective_spacing_A": spacing,
            "spacing_adjusted": spacing != requested_spacing,
            "grid_point_count": dims[0] * dims[1] * dims[2],
            "candidate_point_count": len(candidates),
            "source_cluster_count": len(source_clusters),
            "source_component_count": len(source_components),
            "considered_source_component_ids": [c.id for c in considered_sources],
            "cropped_point_count": cropped_point_count,
            "selected_source_components": [component.to_dict() for component in selected_sources],
            "max_components": max_components,
            "max_export_points": max_export_points,
            "measurement_point_count": len(full_points),
            "display_subsampled": len(export_points) < len(full_points),
            "focus_points": focus_points,
            "focus_radius": focus_radius if focus_points is not None else None,
            "min_component_volume": min_component_volume,
            "max_component_volume": max_component_volume,
            "rejected_component_ids": rejected,
            "selection": selection,
            "enclosure_test": "exterior_flood_fill" if mode == "cavity" else "not_performed",
            "enclosure_scope": "absence of a segment-clearance-valid grid path to the box boundary",
            "connectivity": connectivity.to_dict(),
            "selection_note": "Geometric selection; biological relevance requires a supplied site or evidence",
            "export_mode": export_mode,
            "export": export_metadata,
            "padding": effective_padding,
            "shell_radius": shell_radius,
            "min_octants": required_octants,
            "min_component_points": min_component_points,
            "final_min_component_points": final_min_component_points,
            "include_hydrogen": include_hydrogen,
            "include_hetero": include_hetero,
        },
    )


def component_profile(component: VoidComponent) -> tuple[ChannelPoint, ...]:
    """Return a component's points sorted by ``t``, then by position.

    Parameters
    ----------
    component : VoidComponent

    Returns
    -------
    tuple of ChannelPoint
        Same points, ordered along the cast's axis coordinate.
    """

    return tuple(sorted(component.points, key=lambda point: (point.t, point.position)))


def validate_void_cast(
    frame: StructureFrame,
    cast: VoidCast,
    *,
    tolerance: float = 1e-6,
    include_hydrogen: bool = False,
    include_hetero: bool = False,
) -> dict[str, float | int | bool]:
    """Check that a cast's display points satisfy its clearance bound.

    The clearance of every point in ``cast.points`` is recomputed from the atoms
    of ``frame`` and compared with ``cast.min_radius``.

    Parameters
    ----------
    frame : StructureFrame
        Structure the cast was computed from.
    cast : VoidCast
    tolerance : float, default 1e-6
        Allowed shortfall in Å before a point counts as a violation.
    include_hydrogen : bool, default False
    include_hetero : bool, default False
        Atom selection for the check. Use the same selection as the cast, or
        the check will not be comparable.

    Returns
    -------
    dict
        ``point_count``, ``min_clearance`` (Å, 0.0 when empty),
        ``required_min_radius``, ``violation_count`` and ``is_physical``
        (``True`` when there are no violations).

    Notes
    -----
    Only the display points are checked. Rolling-probe casts record
    ``min_radius=0``, so for them this only checks that no point is inside an
    atom.
    """

    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    spatial = SpatialIndex(atoms)
    clearances = [spatial.nearest_surface(point.position)[0] for point in cast.points]
    minimum = min(clearances, default=0.0)
    violations = sum(1 for clearance in clearances if clearance + tolerance < cast.min_radius)
    return {
        "point_count": len(cast.points),
        "min_clearance": minimum,
        "required_min_radius": cast.min_radius,
        "violation_count": violations,
        "is_physical": violations == 0,
    }


def _collect_void_nodes(
    spatial: SpatialIndex,
    *,
    min_corner: Coord,
    spacing: float,
    dims: tuple[int, int, int],
    min_radius: float,
    shell_radius: float,
    min_octants: int,
) -> dict[GridIndex, _GridNode]:
    """Grid nodes with clearance >= ``min_radius`` and at least ``min_octants``
    occupied octants within ``shell_radius``.
    """
    candidates: dict[GridIndex, _GridNode] = {}
    for idx in itertools.product(range(dims[0]), range(dims[1]), range(dims[2])):
        point = _grid_point(min_corner, spacing, idx)
        clearance, nearest = spatial.nearest_surface(point)
        if clearance < min_radius:
            continue
        coverage = _directional_coverage(point, spatial, shell_radius)
        if coverage < min_octants:
            continue
        candidates[idx] = _GridNode(point, clearance, nearest, coverage)
    return candidates



def _collect_buried_void_nodes(
    spatial: SpatialIndex,
    *,
    connectivity: GridConnectivity,
    min_corner: Coord,
    spacing: float,
    dims: tuple[int, int, int],
    min_radius: float,
    shell_radius: float,
    min_octants: int,
) -> dict[GridIndex, _GridNode]:
    """Free nodes (clearance >= ``min_radius``) that the exterior flood fill does
    not reach, filtered by ``min_octants``.
    """
    empty: dict[GridIndex, _GridNode] = {}
    for idx in itertools.product(range(dims[0]), range(dims[1]), range(dims[2])):
        point = _grid_point(min_corner, spacing, idx)
        clearance, nearest = spatial.nearest_surface(point)
        if clearance < min_radius:
            continue
        coverage = _directional_coverage(point, spatial, shell_radius)
        empty[idx] = _GridNode(point, clearance, nearest, coverage)

    exterior = _flood_exterior(empty, dims, edge_allowed=lambda a, b:
                              connectivity.allows(a, b, empty[a].clearance, empty[b].clearance))
    return {
        idx: node
        for idx, node in empty.items()
        if idx not in exterior and node.coverage >= min_octants
    }


def _flood_exterior(empty: dict[GridIndex, _GridNode], dims: tuple[int, int, int], *,
                    edge_allowed: Callable[[GridIndex, GridIndex], bool]) -> set[GridIndex]:
    """Breadth-first flood fill from every free node on the box boundary through
    edges accepted by ``edge_allowed``; returns the exterior node set.
    """
    seeds = [idx for idx in empty if any(i == 0 or i == n - 1 for i, n in zip(idx, dims))]
    exterior: set[GridIndex] = set(seeds)
    queue: deque[GridIndex] = deque(seeds)
    while queue:
        idx = queue.popleft()
        for offset in _NEIGHBOR_OFFSETS:
            neighbor = (idx[0] + offset[0], idx[1] + offset[1], idx[2] + offset[2])
            if neighbor in empty and neighbor not in exterior and edge_allowed(idx, neighbor):
                exterior.add(neighbor)
                queue.append(neighbor)
    return exterior

def _build_components(
    frame: StructureFrame,
    clusters: Sequence[Sequence[GridIndex]],
    candidates: dict[GridIndex, _GridNode],
    *,
    min_corner: Coord,
    spacing: float,
    dims: tuple[int, int, int],
    min_component_points: int,
    axis_origin: Coord,
    axis_direction: Coord,
    shell_radius: float,
) -> tuple[VoidComponent, ...]:
    """Turn candidate clusters into source components: clearance-weighted centre,
    boundary faces, kind, score, volume ``n * spacing**3`` and nearby residues.
    Clusters smaller than ``min_component_points`` are dropped.
    """
    components: list[VoidComponent] = []
    for cluster in clusters:
        if len(cluster) < min_component_points:
            continue
        nodes = [candidates[idx] for idx in cluster]
        coords = [node.position for node in nodes]
        clearances = [node.clearance for node in nodes]
        weights = [max(clearance, 0.1) for clearance in clearances]
        center = _weighted_centroid(coords, weights)
        boundary_faces = _boundary_faces(cluster, dims)
        axis_values = [project_onto_axis(coord, axis_origin, axis_direction) for coord in coords]
        spans = _coordinate_spans(coords)
        point_count = len(nodes)
        volume = point_count * spacing**3
        kind = _source_component_kind(boundary_faces)
        points = _component_points(cluster, candidates, axis_values)
        nearest = _nearest_residues(frame, center, limit=10, cutoff=shell_radius)
        score = _component_score(
            kind=kind,
            volume=volume,
            clearances=clearances,
            boundary_faces=boundary_faces,
            axis_span=max(axis_values) - min(axis_values),
        )
        components.append(
            VoidComponent(
                id=len(components) + 1,
                kind=kind,
                center=center,
                volume=volume,
                points=points,
                boundary_faces=boundary_faces,
                nearest_residues=nearest,
                score=score,
                metadata={
                    "axis_span": max(axis_values) - min(axis_values),
                    "span_x": spans[0],
                    "span_y": spans[1],
                    "span_z": spans[2],
                    "mean_local_octants": sum(node.coverage for node in nodes) / point_count,
                    "grid_point_count": point_count,
                    "grid_spacing": spacing,
                    "grid_origin": min_corner,
                },
            )
        )
    return tuple(components)


def _component_points(
    cluster: Sequence[GridIndex],
    candidates: dict[GridIndex, _GridNode],
    axis_values: Sequence[float],
) -> tuple[ChannelPoint, ...]:
    """Cluster nodes as ChannelPoints sorted by axis coordinate, then grid index."""
    rows = sorted(zip(cluster, axis_values), key=lambda item: (item[1], item[0]))
    points: list[ChannelPoint] = []
    for index, (idx, axis_value) in enumerate(rows):
        node = candidates[idx]
        points.append(
            ChannelPoint(
                index=index,
                position=node.position,
                radius=node.clearance,
                raw_clearance=node.clearance,
                t=axis_value,
                nearest_atom_serial=node.nearest.serial,
                nearest_residue=node.nearest.residue_key.label,
            )
        )
    return tuple(points)


def _eligible_source_components(
    components: Sequence[VoidComponent], *, mode: str
) -> tuple[VoidComponent, ...]:
    """Sources a mode may use: cavity mode keeps buried components only; channel
    mode keeps open (through/surface) components when any exist; otherwise all.
    """
    if mode == "cavity":
        return tuple(c for c in components if c.kind == "buried_cavity")
    if mode == "channel":
        open_components = tuple(c for c in components if c.kind in {"through_channel", "surface_connected"})
        if open_components:
            return open_components
    return tuple(components)


def _select_final_components(
    components: Sequence[VoidComponent], *, max_components: int,
    min_volume: float, max_volume: float | None,
) -> tuple[tuple[VoidComponent, ...], dict[str, Any]]:
    """Filter all final candidates before applying the volume-ranked limit.

    Candidates arrive sorted by descending volume. Each is marked
    ``below_min_volume``, ``above_max_volume``, ``component_limit`` or
    ``selected``; selected ones are renumbered from 1.
    """
    selected: list[VoidComponent] = []
    decisions: list[dict[str, Any]] = []
    eligible_count = 0
    for component in components:
        output_id = None
        if component.volume < min_volume:
            decision = "below_min_volume"
        elif max_volume is not None and component.volume > max_volume:
            decision = "above_max_volume"
        else:
            eligible_count += 1
            if len(selected) >= max_components:
                decision = "component_limit"
            else:
                decision = "selected"
                output_id = len(selected) + 1
                selected.append(replace(component, id=output_id, metadata={
                    **component.metadata, "candidate_component_id": component.id}))
        decisions.append({"candidate_component_id": component.id,
                          "source_component_id": component.metadata["source_component_id"],
                          "volume": component.volume, "point_count": component.point_count,
                          "decision": decision, "output_component_id": output_id})
    return tuple(selected), {
        "order": ["source_mode_eligibility", "centerline_crop", "focus_crop",
                  "final_connected_components", "final_component_point_threshold", "inclusive_volume_bounds",
                  "descending_final_volume", "component_limit", "display_sampling"],
        "tie_break": "lexicographic minimum grid coordinate",
        "candidate_component_count": len(components),
        "eligible_component_count": eligible_count,
        "selected_component_count": len(selected),
        "omitted_by_limit_count": eligible_count - len(selected),
        "candidates": decisions,
    }


def _export_mode(mode: str, components: Sequence[VoidComponent]) -> str:
    """``"centerline_fill"`` for channel/auto modes with any open source,
    otherwise ``"component_fill"``.
    """
    if mode in {"channel", "auto"} and any(component.kind != "buried_cavity" for component in components):
        return "centerline_fill"
    return "component_fill"


def _export_points(
    components: Sequence[VoidComponent],
    *,
    max_points: int,
    export_mode: str,
    grid_origin: Coord,
    spacing: float,
    dims: tuple[int, int, int],
    centerline: PoreProfile | None,
    min_radius: float,
    channel_radius_scale: float,
    channel_radius_padding: float,
    channel_max_radius: float | None,
) -> tuple[tuple[ChannelPoint, ...], dict[str, float | int | str | None]]:
    """Points to carry forward: a centre-line crop of the sources, or all source points."""
    if export_mode == "centerline_fill" and centerline is not None:
        selected, metadata = _centerline_fill_points(
            components,
            centerline,
            spacing=spacing,
            min_radius=min_radius,
            radius_scale=channel_radius_scale,
            radius_padding=channel_radius_padding,
            max_radius=channel_max_radius,
        )
        points = _renumber_points(selected)
        metadata["point_count"] = len(points)
        metadata["mode"] = export_mode
        return points, metadata

    all_points = [point for component in components for point in component.points]
    points = _renumber_points(all_points)
    return points, {"mode": "component_fill", "path_count": 0, "point_count": len(points)}


def _centerline_fill_points(
    components: Sequence[VoidComponent],
    centerline: PoreProfile,
    *,
    spacing: float,
    min_radius: float,
    radius_scale: float,
    radius_padding: float,
    max_radius: float | None,
) -> tuple[list[ChannelPoint], dict[str, float | int | str | None]]:
    """Keep source points near the centre line.

    Each point is matched to the profile point nearest in axial projection. It is
    kept when that profile point's raw clearance is at least ``min_radius`` and
    the point lies within ``local_radius * radius_scale + radius_padding`` of it,
    where ``local_radius = min(max(raw_clearance, radius, min_radius), cap)``.
    """
    profile_points = tuple(sorted(centerline.points, key=lambda point: point.t))
    if not profile_points:
        return [], {"centerline_points": 0, "radius_cap": max_radius, "path_count": 0}
    ts = [point.t for point in profile_points]
    radius_cap = max_radius if max_radius is not None else _profile_radius_cap(profile_points, spacing)
    selected: list[ChannelPoint] = []
    for point in (point for component in components for point in component.points):
        t = project_onto_axis(point.position, centerline.axis_origin, centerline.axis_direction)
        nearest = _nearest_profile_point(t, profile_points, ts)
        local_radius = min(max(nearest.raw_clearance, nearest.radius, min_radius), radius_cap)
        if nearest.raw_clearance < min_radius:
            continue
        limit = local_radius * radius_scale + radius_padding
        if distance(point.position, nearest.position) <= limit:
            selected.append(point)
    return selected, {
        "centerline_points": len(profile_points),
        "radius_cap": radius_cap,
        "radius_scale": radius_scale,
        "radius_padding": radius_padding,
        "path_count": 0,
    }


def _components_from_export_points(
    frame: StructureFrame,
    points: Sequence[ChannelPoint],
    *,
    source_components: Sequence[VoidComponent],
    edge_allowed: Callable[[GridIndex, GridIndex], bool],
    export_mode: str,
    grid_origin: Coord,
    spacing: float,
    dims: tuple[int, int, int],
    shell_radius: float,
    min_component_points: int,
) -> tuple[VoidComponent, ...]:
    """Re-cluster the surviving points into final components.

    Cropping can only split sources, never join them. In centre-line mode every
    component is ``"centerline_channel"``; otherwise the kind follows the box
    faces touched, with ``"cropped_open_region"`` when a crop makes part of an
    open source look buried. Score equals volume. Sorted by descending volume,
    then smallest point coordinate.
    """
    if not points:
        return ()
    point_by_idx = {_index_from_point(point.position, grid_origin, spacing): point for point in points}
    source_by_idx = {_index_from_point(p.position, grid_origin, spacing): c
                     for c in source_components for p in c.points}
    clusters = _cluster_point_indices(set(point_by_idx), edge_allowed=edge_allowed)
    components: list[VoidComponent] = []
    for cluster in clusters:
        if len(cluster) < min_component_points:
            continue
        cluster_points = [point_by_idx[idx] for idx in cluster]
        weights = [max(point.radius, 0.1) for point in cluster_points]
        center = _weighted_centroid([point.position for point in cluster_points], weights)
        boundary_faces = _boundary_faces(cluster, dims)
        # Cropping deletes nodes; it cannot join previously separate sources.
        source = source_by_idx[cluster[0]]
        kind = "centerline_channel" if export_mode == "centerline_fill" else _source_component_kind(boundary_faces)
        if kind == "buried_cavity" and source.kind != "buried_cavity":
            kind = "cropped_open_region"
        nearest = _nearest_residues(frame, center, limit=10, cutoff=shell_radius)
        components.append(
            VoidComponent(
                id=len(components) + 1,
                kind=kind,
                center=center,
                volume=len(cluster_points) * spacing**3,
                points=tuple(sorted(cluster_points, key=lambda point: (point.t, point.position))),
                boundary_faces=boundary_faces,
                nearest_residues=nearest,
                score=len(cluster_points) * spacing**3,
                metadata={
                    "grid_point_count": len(cluster_points),
                    "grid_spacing": spacing,
                    "grid_origin": grid_origin,
                    "export_mode": export_mode,
                    "source_component_id": source.id,
                    "source_kind": source.kind,
                    "source_boundary_faces": source.boundary_faces,
                    "selection_cropped": len(cluster_points) < source.point_count,
                },
            )
        )
    components = sorted(components, key=lambda c: (-c.volume, min(p.position for p in c.points)))
    return tuple(_renumber_components(components))



def _renumber_components(components: Sequence[VoidComponent]) -> tuple[VoidComponent, ...]:
    """Copy components with ids 1..n in the given order."""
    return tuple(
        VoidComponent(
            id=index,
            kind=component.kind,
            center=component.center,
            volume=component.volume,
            points=component.points,
            boundary_faces=component.boundary_faces,
            nearest_residues=component.nearest_residues,
            score=component.score,
            metadata=component.metadata,
        )
        for index, component in enumerate(components, start=1)
    )

def _ensure_centerline_profile(
    frame: StructureFrame,
    *,
    profile: PoreProfile | None,
    axis: str | Coord,
    atoms: tuple[Atom, ...],
    spacing: float,
    samples: int | None,
    search_radius: float,
    include_hydrogen: bool,
    include_hetero: bool,
) -> PoreProfile:
    """Return ``profile`` or compute one with :func:`crevice.analysis.pore_profile`;
    the default sample count is ``ceil(span / spacing) + 1`` clamped to 17-161.
    """
    if profile is not None:
        return profile
    axis_origin, axis_direction, _axis_label = axis_from_name(axis, atoms)
    projections = [project_onto_axis(atom.coord, axis_origin, axis_direction) for atom in atoms]
    span = max(projections) - min(projections)
    inferred_samples = max(17, min(161, int(math.ceil(span / max(spacing, 0.1))) + 1))
    from .analysis import pore_profile

    return pore_profile(
        frame,
        axis=axis,
        samples=samples or inferred_samples,
        search_radius=search_radius,
        enclosure_radius=0.8,  # fixed legacy probe; no enclosure option here
        refinement_steps=4,
        include_hydrogen=include_hydrogen,
        include_hetero=include_hetero,
    )


def _profile_radius_cap(points: Sequence[ChannelPoint], spacing: float) -> float:
    """Default centre-line radius cap: the 85th-percentile profile radius plus one
    spacing, and never less than one spacing.
    """
    radii = sorted(max(point.raw_clearance, point.radius, 0.0) for point in points)
    if not radii:
        return spacing
    index = min(len(radii) - 1, max(0, int(round(0.85 * (len(radii) - 1)))))
    return max(spacing, radii[index] + spacing)


def _nearest_profile_point(t: float, points: Sequence[ChannelPoint], ts: Sequence[float]) -> ChannelPoint:
    """Profile point nearest in ``t`` (binary search; ties go to the lower point)."""
    index = bisect_left(ts, t)
    if index <= 0:
        return points[0]
    if index >= len(points):
        return points[-1]
    before = points[index - 1]
    after = points[index]
    return before if abs(t - before.t) <= abs(after.t - t) else after


def _renumber_points(points: Sequence[ChannelPoint]) -> tuple[ChannelPoint, ...]:
    """Copy points with ``index`` set to 0..n-1 in the given order."""
    renumbered: list[ChannelPoint] = []
    for index, point in enumerate(points):
        renumbered.append(
            ChannelPoint(
                index=index,
                position=point.position,
                radius=point.radius,
                raw_clearance=point.raw_clearance,
                t=point.t,
                nearest_atom_serial=point.nearest_atom_serial,
                nearest_residue=point.nearest_residue,
            )
        )
    return tuple(renumbered)



def _index_from_point(point: Coord, grid_origin: Coord, spacing: float) -> GridIndex:
    """Grid index of a position on the axis-aligned grid (rounded to nearest)."""
    return (
        int(round((point[0] - grid_origin[0]) / spacing)),
        int(round((point[1] - grid_origin[1]) / spacing)),
        int(round((point[2] - grid_origin[2]) / spacing)),
    )

def _sample_points_by_axis(points: Sequence[ChannelPoint], max_points: int) -> list[ChannelPoint]:
    """Sort points by ``(t, position)`` and thin to ``max_points`` with :func:`_sample_evenly`."""
    ordered = sorted(points, key=lambda point: (point.t, point.position))
    if len(ordered) <= max_points:
        return list(ordered)
    return _sample_evenly(ordered, max_points)


def _sample_evenly(points: Sequence[ChannelPoint], count: int) -> list[ChannelPoint]:
    """Evenly strided subset of ``count`` points that always keeps the minimum
    radius point. If adding it overflows ``count``, the widest points are
    dropped. ``count <= 1`` returns only the minimum-radius point.
    """
    if count >= len(points):
        return list(points)
    if count <= 1:
        return [min(points, key=lambda point: point.radius)]
    stride = (len(points) - 1) / (count - 1)
    indices = {round(i * stride) for i in range(count)}
    bottleneck_index = min(range(len(points)), key=lambda index: points[index].radius)
    indices.add(bottleneck_index)
    selected = [points[index] for index in sorted(indices)]
    if len(selected) > count:
        selected.sort(key=lambda point: (point.radius, point.t, point.position))
        keep = set(selected[:1])
        for point in selected[1:]:
            if len(keep) >= count:
                break
            keep.add(point)
        selected = [point for point in points if point in keep]
    return selected[:count]


def _cluster_point_indices(indices: set[GridIndex], *,
                           edge_allowed: Callable[[GridIndex, GridIndex], bool]) -> list[list[GridIndex]]:
    """Connected components of grid indices under ``edge_allowed`` (26-neighbour
    breadth-first search from sorted seeds), largest first, ties by smallest
    index.
    """
    unseen = set(indices)
    clusters: list[list[GridIndex]] = []
    for seed in sorted(indices):
        if seed not in unseen:
            continue
        unseen.remove(seed)
        cluster = [seed]
        queue: deque[GridIndex] = deque([seed])
        while queue:
            idx = queue.popleft()
            for offset in _NEIGHBOR_OFFSETS:
                neighbor = (idx[0] + offset[0], idx[1] + offset[1], idx[2] + offset[2])
                if neighbor in unseen and edge_allowed(idx, neighbor):
                    unseen.remove(neighbor)
                    cluster.append(neighbor)
                    queue.append(neighbor)
        clusters.append(cluster)
    clusters.sort(key=lambda cluster: (-len(cluster), min(cluster)))
    return clusters


def _cluster_indices(candidates: dict[GridIndex, _GridNode], *,
                     edge_allowed: Callable[[GridIndex, GridIndex], bool]) -> list[list[GridIndex]]:
    """Connected components of candidate nodes; see :func:`_cluster_point_indices`."""
    return _cluster_point_indices(set(candidates), edge_allowed=edge_allowed)


def _directional_coverage(point: Coord, spatial: SpatialIndex, shell_radius: float) -> int:
    """Number of octants (sign patterns of ``atom - point``, zero counted as
    positive) that contain an atom centre within ``shell_radius``; 0 to 8.
    """
    octants: set[tuple[int, int, int]] = set()
    for atom in spatial.atoms_within_distance(point, shell_radius):
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


def _source_component_kind(boundary_faces: tuple[str, ...]) -> str:
    """``"through_channel"`` if two opposite box faces are touched,
    ``"surface_connected"`` if any face is, otherwise ``"buried_cavity"``.
    """
    if _opposing_face_count(boundary_faces) > 0:
        return "through_channel"
    if boundary_faces:
        return "surface_connected"
    return "buried_cavity"


def _component_score(
    *,
    kind: str,
    volume: float,
    clearances: Sequence[float],
    boundary_faces: tuple[str, ...],
    axis_span: float,
) -> float:
    """Heuristic source ranking score (not used for the final selection order):
    ``topology * faces * volume**(2/3) * (mean + max clearance)/2 * (1 + 0.02 * axis_span)``
    with topology 3.0/1.5/1.0 for through/surface/buried and
    ``faces = 1 + 0.25 * n_faces + 0.5 * n_opposing_pairs``.
    """
    mean_radius = sum(clearances) / len(clearances)
    max_radius = max(clearances)
    topology_bonus = 3.0 if kind == "through_channel" else 1.5 if kind == "surface_connected" else 1.0
    face_bonus = 1.0 + 0.25 * len(boundary_faces) + 0.5 * _opposing_face_count(boundary_faces)
    return topology_bonus * face_bonus * (volume ** (2.0 / 3.0)) * (0.5 * mean_radius + 0.5 * max_radius) * (1.0 + 0.02 * axis_span)


def _boundary_faces(cluster: Sequence[GridIndex], dims: tuple[int, int, int]) -> tuple[str, ...]:
    """Sorted names of the grid-box faces (``x-``, ``x+``, ...) touched by a cluster."""
    faces: set[str] = set()
    edge_values = {0: (0, dims[0] - 1), 1: (0, dims[1] - 1), 2: (0, dims[2] - 1)}
    for idx in cluster:
        for axis, side_index, name in _FACE_NAMES:
            side = edge_values[axis][side_index]
            if idx[axis] == side:
                faces.add(name)
    return tuple(sorted(faces))


def _opposing_face_count(faces: Sequence[str]) -> int:
    """Number of axes (0-3) on which both opposite faces are touched."""
    face_set = set(faces)
    return sum(
        1
        for negative, positive in (("x-", "x+"), ("y-", "y+"), ("z-", "z+"))
        if negative in face_set and positive in face_set
    )


def _nearest_residues(frame: StructureFrame, point: Coord, *, limit: int, cutoff: float) -> tuple[str, ...]:
    """Residue labels near a point via :func:`crevice.analysis.nearest_residues`."""
    from .analysis import nearest_residues

    return nearest_residues(frame, point, limit=limit, cutoff=cutoff)


def _scale_spacing_for_box(min_corner: Coord, max_corner: Coord, spacing: float, max_grid_points: int) -> float:
    """Enlarge ``spacing`` by ``(count / max_grid_points)**(1/3)`` when the grid
    would exceed ``max_grid_points``; otherwise return it unchanged.
    """
    dims = _grid_dimensions(min_corner, max_corner, spacing)
    count = dims[0] * dims[1] * dims[2]
    if count <= max_grid_points:
        return spacing
    factor = (count / max_grid_points) ** (1.0 / 3.0)
    return spacing * factor


def _grid_dimensions(min_corner: Coord, max_corner: Coord, spacing: float) -> tuple[int, int, int]:
    """Points per axis, ``floor(extent / spacing) + 1``, at least 2."""
    return (
        max(2, int(math.floor((max_corner[0] - min_corner[0]) / spacing)) + 1),
        max(2, int(math.floor((max_corner[1] - min_corner[1]) / spacing)) + 1),
        max(2, int(math.floor((max_corner[2] - min_corner[2]) / spacing)) + 1),
    )


def _grid_point(min_corner: Coord, spacing: float, idx: GridIndex) -> Coord:
    """Position of grid index ``idx`` on the axis-aligned grid (Å)."""
    return (
        min_corner[0] + idx[0] * spacing,
        min_corner[1] + idx[1] * spacing,
        min_corner[2] + idx[2] * spacing,
    )


def _padded_box(min_corner: Coord, max_corner: Coord, padding: float) -> tuple[Coord, Coord]:
    """Box enlarged by ``padding`` on all sides (unchanged when padding is 0)."""
    if padding == 0:
        return min_corner, max_corner
    return (
        (min_corner[0] - padding, min_corner[1] - padding, min_corner[2] - padding),
        (max_corner[0] + padding, max_corner[1] + padding, max_corner[2] + padding),
    )


def _weighted_centroid(points: Sequence[Coord], weights: Sequence[float]) -> Coord:
    """Weighted mean position; the plain centroid when all weights are zero."""
    total = sum(weights)
    if total == 0:
        return centroid(points)
    return (
        sum(point[0] * weight for point, weight in zip(points, weights)) / total,
        sum(point[1] * weight for point, weight in zip(points, weights)) / total,
        sum(point[2] * weight for point, weight in zip(points, weights)) / total,
    )


def _coordinate_spans(points: Sequence[Coord]) -> Coord:
    """Extent of the points along x, y and z (Å)."""
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))

"""Engines of the pore profile, cavity and residue-contact analyses.

The documented entry points for the first two are the thin wrappers
:func:`crevice.channels.pore_profile` and :func:`crevice.cavities.detect_cavities`
(both exported from ``crevice``); they pass their keyword arguments here
unchanged, and their docstrings give the full method description. :func:`annotate_residues` and
:func:`nearest_residues` are documented here and re-exported by
:mod:`crevice.residues`. All lengths are in Å, atoms are van der Waals
spheres from the radius set in effect, and every result is geometry at the
chosen resolution, not a validated biological assignment.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Sequence

from .models import Cavity, ChannelPoint, Coord, PoreProfile, ResidueContact, StructureFrame
from .radii import RadiusSet, radii_option, residue_properties
from .spatial import SpatialIndex


@radii_option
def pore_profile(
    frame: StructureFrame,
    *,
    axis: str | Coord = "auto",
    origin: Coord | None = None,
    samples: int = 81,
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

    The channel is followed plane by plane from the seed (:mod:`crevice.adaptive_channel`):
    each plane's centre is climbed from the previous one to the in-plane
    clearance maximum, and the plane belongs to the channel while its lumen is
    laterally enclosed. ``axis="auto"`` tries the principal shape directions
    and keeps a unique channel. Closed cavities and unresolved or ambiguous
    channels need a separate region selection (casts, cavities, paths).

    ``lateral_exits=True`` (off by default) accepts a capped end, whose
    straight axial exit is blocked, when a lateral exit to bulk exists; every
    distinct exit leg is reported separately in ``metadata["exits"]`` and
    ``metadata["path_type"]``.

    ``enclosure_radius`` (Å) is the probe that decides lateral enclosure and
    connectivity. ``None`` (default, also ``"auto"``) chooses it per plane:
    the smallest probe of a 0.30-3.00 Å ladder at which that plane's lumen is
    enclosed (each sample's probe in ``metadata["enclosure_probes_A"]``). A
    number runs exactly that probe in every plane.

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    axis : {"auto", "x", "y", "z"} or Coord, default "auto"
        Channel direction.
    origin : Coord, optional
        Seed inside the intended channel, Å; default the atom centroid.
    samples : int, default 81
        Minimum number of samples (axial steps are kept at or below 0.75 Å);
        at least 2.
    search_radius : float, default 6.0
        Base transverse search half-width, Å (must be positive); a wide lumen
        uses a square twice or four times as wide.
    refinement_steps : int, default 4
        Minimum step halvings of each in-plane centre climb (at least 1); the
        climb converges to the local clearance maximum within 0.0001 Å.
    probe_radius : float, default 0.0
        Measurement probe radius, Å, subtracted from each clearance.
    include_hydrogen : bool, default False
        Treat hydrogen atoms as obstacles.
    include_hetero : bool, default True
        Treat HETATM atoms as obstacles.
    section_spacing : float, default 0.5
        Transverse grid spacing, Å.
    enclosure_radius : float, "auto" or None, default None
        Enclosure probe, Å; ``None`` or ``"auto"`` = per plane (above).
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
    if search_radius == 0:
        raise ValueError(
            "search_radius=0 (the unchecked fixed-axis clearance scan) was removed in CREVICE 0.5; "
            "profiles always follow a connected, laterally enclosed channel"
        )
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    if origin is not None and (len(origin) != 3 or not all(math.isfinite(v) for v in origin)):
        raise ValueError("origin must be a finite 3D coordinate")
    if not math.isfinite(section_spacing) or section_spacing <= 0:
        raise ValueError("section_spacing must be finite and positive")
    if isinstance(enclosure_radius, str):
        if enclosure_radius.lower() not in {"auto", "adaptive"}:
            raise ValueError(
                "enclosure_radius must be None or 'auto' (per-plane probe) or a finite non-negative "
                "radius ('ladder' was removed in CREVICE 0.5)"
            )
        enclosure_radius = None
    elif enclosure_radius is not None and (not math.isfinite(enclosure_radius) or enclosure_radius < 0):
        raise ValueError("enclosure_radius must be None (automatic) or finite and non-negative")
    from .sections import connected_profile

    profile = connected_profile(
        atoms,
        axis=axis,
        origin=origin,
        samples=samples,
        search_radius=search_radius,
        refinement_steps=refinement_steps,
        probe_radius=probe_radius,
        section_spacing=section_spacing,
        enclosure_radius=enclosure_radius,
        lateral_exits=lateral_exits,
        exit_bulk_radius=exit_bulk_radius,
        exit_spacing=exit_spacing,
    )
    return replace(
        profile, metadata={**profile.metadata, "include_hydrogen": include_hydrogen, "include_hetero": include_hetero}
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
    and probe clearance. Open pockets, channels and tunnels are not cavities.
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
    cast = void_cast(
        frame,
        mode="cavity" if enclosed_only else "all",
        spacing=spacing,
        min_radius=min_radius,
        max_components=max_cavities,
        shell_radius=shell_radius,
        min_octants=0 if enclosed_only else min_octants,
        max_grid_points=max_grid_points,
        include_hydrogen=include_hydrogen,
        include_hetero=include_hetero,
        focus_points=focus_points,
        focus_radius=focus_radius,
        min_component_volume=min_component_volume,
        max_component_volume=max_component_volume,
    )
    return tuple(
        Cavity(c.id, c.center, c.max_radius, c.volume, c.point_count, c.nearest_residues, c.score)
        for c in cast.components
    )


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
    region_radii = tuple(max(0.0, point.raw_clearance) if isinstance(point, ChannelPoint) else 0.0 for point in points)
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

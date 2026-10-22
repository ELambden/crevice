"""Static-structure workflow of ``crevice report`` and ``crevice cavities --cast``, split into compute and export.

:func:`compute_static_workflow` and :func:`compute_rolling_cast` measure and
return results without writing any file; :func:`export_static_workflow` and
:func:`export_rolling_cast` write them. The bundle writers
:func:`crevice.figures.write_static_publication_bundle` and
:func:`crevice.figures.write_rolling_probe_bundle` are thin wrappers (compute,
then export), so the three routes give identical numbers and files.

Examples
--------
Measure first, decide what to write afterwards::

    from crevice import load_structure
    from crevice.static_workflow import compute_static_workflow, export_static_workflow
    frame = load_structure("1GRM.cif")
    result = compute_static_workflow(frame, samples=101, section_spacing=0.25, cast_spacing=0.25)
    print(result.status, result.profile.min_radius, result.cast.total_volume)
    export_static_workflow(result, structure_path="1GRM.cif", output_dir="out", prefix="1GRM",
                           outputs="figures,pymol")
"""
from __future__ import annotations

import json
import math
import shutil
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Sequence

from .presentation import annotate_option
from .radii import RadiusSet, radii_fields, radii_option


@dataclass
class RollingCastResult:
    """A measured 3D rolling-probe cast and its display regions (no files).

    Attributes
    ----------
    frame : StructureFrame
        The obstacle atoms actually used (after ``selection``).
    intake : dict
        Atom-selection report (written as ``PREFIX_selection.json``).
    cast : VoidCast
        The measured cast.
    labels : dict of int to str
        Display label of each component (dominant selection).
    regions : dict of str to sequence of ChannelPoint
        Disjoint display regions (empty for a single region).
    region_summary : list of dict
        Region records (written as ``PREFIX_region_annotations.json``).
    selection, selection_mode : str
        The obstacle selection and the component selection mode used.
    """

    frame: Any
    intake: dict
    cast: Any
    labels: dict = field(default_factory=dict)
    regions: dict = field(default_factory=dict)
    region_summary: list = field(default_factory=list)
    selection: str = "protein"
    selection_mode: str = "dominant"


@dataclass
class StaticWorkflowResult:
    """Everything ``crevice report`` measures for one structure (no files).

    Attributes
    ----------
    status : {"resolved", "unresolved", "rolling"}
        ``resolved``: a through-channel profile; ``unresolved``: no unique
        channel, a rolling-probe fallback cast instead (no profile is
        invented); ``rolling``: ``cast_mode="rolling"`` was requested.
    frame : StructureFrame
        The analysed atoms (publication selection).
    profile : PoreProfile or None
        The channel profile (``None`` unless resolved).
    contacts : tuple of ResidueContact
        Residues contacting the profile (or the selected cast).
    cast : VoidCast or None
        The profile-guided cast (resolved only).
    exit_casts, segments : tuple
        Lateral exit-leg casts and ENTRY/EXIT segment casts (resolved channels).
    segment_rows, lining : list of dict
        Segment summary rows and lining-residue rows.
    region : PoreProfile or VoidComponent or None
        Region node of the residue network.
    rolling : RollingCastResult or None
        The rolling-probe cast (``unresolved`` and ``rolling``).
    unresolved_reason : str or None
        Why no channel resolved.
    probe_selection : dict or None
        Automatic enclosure-probe record of an unresolved profile.
    cavities, tunnels, network : object or None
        Optional analyses (``None`` when not requested or not completed).
    connectivity : dict or None
        Lining connectivity report (with ``residue_groups``).
    fingerprint : PoreFingerprint or None
        Pore fingerprint (resolved channels only, when requested).
    extras : dict
        Status records of the optional analyses after an unresolved profile.
    settings : dict
        Settings the export needs (``include_hydrogen``, ``include_hetero``,
        ``cast_probe``, ``entry_end``, ``lining_cutoff``, ``exit_bulk_radius``,
        ``residue_groups``, ``cast_mode``, ``focus``).
    """

    status: str
    frame: Any
    profile: Any = None
    contacts: tuple = ()
    cast: Any = None
    exit_casts: tuple = ()
    segments: Any = None
    segment_rows: list = field(default_factory=list)
    lining: list = field(default_factory=list)
    region: Any = None
    rolling: RollingCastResult | None = None
    unresolved_reason: str | None = None
    probe_selection: dict | None = None
    cavities: Any = None
    tunnels: Any = None
    network: Any = None
    connectivity: dict | None = None
    fingerprint: Any = None
    extras: dict = field(default_factory=dict)
    settings: dict = field(default_factory=dict)


@radii_option
def compute_rolling_cast(frame, *, spacing=.5, probe_radius=.8, outer_radius=6., enclosure_fraction=.9, min_depth=0.,
                         min_component_volume=50., max_components=1, selection_mode="dominant", core_fraction=.08,
                         max_grid_points=16_000_000, selection="protein", seed=None, seed_tolerance=.75,
                         radii: RadiusSet | str | None = None) -> RollingCastResult:
    """Measure a 3D rolling-probe cast and its display regions without writing files.

    Parameters
    ----------
    frame : StructureFrame
        Parsed structure.
    spacing : float, default 0.5
        Grid spacing (Å); never silently coarsened.
    probe_radius : float, default 0.8
        Inner rolling-probe radius (Å).
    outer_radius : float, default 6.0
        Outer envelope probe radius (Å).
    enclosure_fraction : float, default 0.9
        Burial threshold for cavity cores.
    min_depth : float, default 0.0
        Extra depth inside the outer envelope (Å).
    min_component_volume : float, default 50.0
        Minimum region volume (Å³).
    max_components : int, default 1
        Maximum qualifying regions.
    selection_mode : {"dominant", "all"}, default "dominant"
        Substantial interior cores, or every void.
    core_fraction : float, default 0.08
        Minimum core volume relative to the strongest core (dominant mode).
    max_grid_points : int, default 16000000
        Grid allocation guard.
    selection : {"protein", "all"}, default "protein"
        Obstacles: protein heavy atoms, or all atoms including HETATM records.
    seed : tuple of float, optional
        Interior seed (Å).
    seed_tolerance : float, default 0.75
        Maximum seed-to-grid snap distance (Å).
    radii : RadiusSet, str or None, optional
        Atomic radius set (see :func:`crevice.radii.use_radii`); ``None`` uses
        the set in effect.

    Returns
    -------
    RollingCastResult
    """
    from .rolling import chain_enclosed_regions, rolling_probe_cast, select_cast_atoms
    frame, intake = select_cast_atoms(frame, selection)
    cast = rolling_probe_cast(frame, spacing=spacing, probe_radius=probe_radius, outer_radius=outer_radius,
                              enclosure_fraction=enclosure_fraction, min_depth=min_depth,
                              min_component_volume=min_component_volume, max_components=max_components,
                              max_grid_points=max_grid_points, include_hetero=selection == 'all',
                              selection_mode=selection_mode, core_fraction=core_fraction,
                              seed=seed, seed_tolerance=seed_tolerance)
    labels, regions, region_summary = {}, {}, []
    if cast.components:
        if selection_mode == 'dominant':
            rows = {r['region_id']: r for r in cast.metadata['dominance'].get('regions', [])}
            for component in cast.components:
                row = rows[component.metadata['source_component_id']]
                name = 'shared' if row['kind'] == 'shared_interface_core' else f'cavity_{component.id}'
                regions[name] = component.points
                labels[component.id] = name
                region_summary.append({'label': name, 'volume_A3': component.volume, 'selection': row})
            if len(regions) == 1:
                regions = {}
        else:
            regions, region_summary = chain_enclosed_regions(frame, cast)
    return RollingCastResult(frame, intake, cast, labels, regions, region_summary, selection, selection_mode)


def export_rolling_cast(result: RollingCastResult, *, structure_path, output_dir, prefix="crevice",
                        surface_smoothing=None, smooth=None) -> dict:
    """Write a :class:`RollingCastResult` as the ``crevice cavities --cast`` bundle.

    The files are those of :func:`crevice.figures.write_rolling_probe_bundle`,
    filtered by the active output plan (:mod:`crevice.outputs`).

    Parameters
    ----------
    result : RollingCastResult
        From :func:`compute_rolling_cast`.
    structure_path : str or Path
        Source file, copied into the bundle (file kind ``structures``).
    output_dir : str or Path
        Bundle directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    surface_smoothing : float, optional
        Legacy display interpolation width (Å).
    smooth : float, optional
        Display-surface smoothing length (Å).

    Returns
    -------
    dict of str to str
        Output key to path, including ``manifest_json``.
    """
    from .figures import _copy_structure_for_bundle, _skip_json, _viewer_structure_for_bundle, _write_json
    from .io import write_void_cast_csv, write_void_cast_json
    from .outputs import drop_blocked, wants
    from .volume_export import resolve_display_smoothing, write_void_cast_dx, write_volume_viewer_bundle
    if Path(prefix).name != prefix or not prefix or any(c in prefix for c in '\n\r"'):
        raise ValueError("prefix must be a simple filename stem")
    resolve_display_smoothing(smooth, surface_smoothing)
    frame, cast = result.frame, result.cast
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    original = _copy_structure_for_bundle(structure_path, root, prefix)
    viewer = _viewer_structure_for_bundle(frame, root / f"{prefix}_viewer.pdb")
    manifest = {'structure_file': str(viewer), 'original_structure_file': str(original),
                'viewer_identity_json': str(viewer.with_suffix('.identities.json')),
                'void_cast_json': str(root / f'{prefix}_void_cast.json'),
                'void_cast_csv': str(root / f'{prefix}_void_cast.csv'),
                'selection_json': str(root / f'{prefix}_selection.json')}
    write_void_cast_json(cast, manifest['void_cast_json'])
    _write_json(result.intake, manifest['selection_json'])
    if cast.components:
        manifest['region_annotations_json'] = str(root/f'{prefix}_region_annotations.json')
        (_write_json if wants(manifest['region_annotations_json']) else _skip_json)(
            {'regions': result.region_summary,
             'interpretation': 'Disjoint display annotation; does not change connected components or imply a biological channel count.'},
            manifest['region_annotations_json'])
        manifest.update(write_volume_viewer_bundle(cast, structure_path=viewer, output_dir=root, prefix=prefix, frame=frame,
                                                   cast_extension=0, regions=result.regions,
                                                   surface_smoothing=surface_smoothing, smooth=smooth))
        for extension, key in [('pml', 'volume_pml'), ('tcl', 'volume_tcl'), ('vmd', 'volume_vmd'), ('cxc', 'volume_cxc')]:
            alias = root / f'{prefix}_pore_cast.{extension}'
            if wants(alias):
                shutil.copyfile(manifest[key], alias)
                manifest[f'pore_cast_{extension}'] = str(alias)
        end_views = {extension: root/f'{prefix}_end_view.{extension}' for extension in ('pml', 'tcl', 'vmd', 'cxc')}
        end_views = {extension: path for extension, path in end_views.items() if wants(path)}
        if end_views:
            end_turn = json.loads(Path(manifest['scene_json']).read_text(encoding="utf-8")).get('end_view_turn_degrees', 74)
            turns = {'pml': f'turn x, {end_turn:g}', 'tcl': f'rotate x by {end_turn:g}',
                     'vmd': f'rotate x by {end_turn:g}', 'cxc': f'turn x {end_turn:g}'}
            for extension, end_view in end_views.items():
                end_view.write_text(Path(manifest[f'volume_{extension}']).read_text(encoding="utf-8")+turns[extension]+'\n', encoding="utf-8", newline="\n")
                manifest[f'end_view_{extension}'] = str(end_view)
        # Each retained region can also be inspected independently in a viewer.
        region_dir = root / 'regions'
        region_dir.mkdir(exist_ok=True)
        for component in cast.components:
            write_void_cast_dx(cast, region_dir/f'region_{component.id:03d}.dx', component_id=component.id)
    write_void_cast_csv(cast, manifest['void_cast_csv'], labels=result.labels)
    manifest_path = root / f'{prefix}_manifest.json'
    drop_blocked(manifest)
    _write_json({'files': manifest, 'prefix': prefix, 'atom_selection': result.intake,
                 'void_cast': cast.to_dict(), **radii_fields(),
                 'channel_profile_status': 'not_assigned_by_3D_void_detection',
                 'interpretation': 'Geometric void regions. Component count is not a biological channel count.',
                 'native_viewers': 'Scripts generated; running a viewer is a separate validation step.'}, manifest_path)
    manifest['manifest_json'] = str(manifest_path)
    return manifest


def measured_region_samples(cast):
    """Centres of the measured cast samples exactly as read back from ``PREFIX_volume.dx``.

    The DX text stores the origin with 9 and the step vectors with 12
    significant digits; the same rounding is applied here so in-memory results
    equal those of the file-based 0.2 workflow.

    Parameters
    ----------
    cast : VoidCast
        Measured cast.

    Returns
    -------
    numpy.ndarray or None
        ``(n, 3)`` sample centres (Å), or ``None`` for a cast without samples.
    grid volume : float
        Volume of one grid cell (Å³) (``0.0`` without samples).
    """
    import numpy as np
    from .volume_export import cast_grid
    if not cast.components:
        return None, 0.0
    grid, origin = cast_grid(cast)
    origin = np.asarray([float(f"{v:.9g}") for v in origin])
    basis = cast.metadata.get("grid_basis", ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    deltas = np.asarray([[float(f"{cast.spacing*v:.12g}") for v in direction] for direction in basis])
    xyz = np.argwhere(grid) @ deltas + origin
    if not len(xyz):
        return None, 0.0
    return xyz, float(abs(np.linalg.det(deltas)))


def _unresolved_extras(frame, cast, *, include_cavities, include_tunnels, include_network, include_hydrogen,
                       include_hetero, residue_groups=None):
    """Compute the optional analyses after an unresolved profile (no files).

    Returns ``(records, cavities, tunnels, network, connectivity)``; the
    records match those written to the manifest by the 0.2 workflow, except
    that ``network.region_source`` is filled in by the export.
    """
    import numpy as np
    from .analysis import annotate_residues, detect_cavities
    from .tunnels import find_tunnels
    obstacles = {"include_hydrogen": include_hydrogen, "include_hetero": include_hetero,
                 "definition": "Publication atom selection; the same obstacle atoms as the fallback rolling-probe cast"}
    records: dict = {}
    cavities = tunnels = network = connectivity = None
    region_xyz, cell = measured_region_samples(cast)
    no_region = "The fallback cast has no measured region samples"
    if not include_cavities:
        records["cavities"] = {"status": "skipped_by_request"}
    else:
        cavities = detect_cavities(frame, include_hydrogen=include_hydrogen, include_hetero=include_hetero)
        records["cavities"] = {"status": "completed", "cavity_count": len(cavities), "obstacles": obstacles,
                               "definition": "Separate enclosed-cavity search (detect_cavities defaults); independent of the rolling-probe cast"}
    if not include_tunnels:
        records["tunnels"] = {"status": "skipped_by_request"}
    else:
        start = None
        start_definition = "Structure centroid (find_tunnels default); " + no_region.lower()
        if region_xyz is not None:
            from .radii import atom_vdw_radius
            from .rolling import _SphereQueries
            query = _SphereQueries(np.asarray([a.coord for a in frame.atoms]),
                                   np.asarray([atom_vdw_radius(a) for a in frame.atoms]))
            clearance, _ = query.points(region_xyz)
            start = tuple(float(v) for v in region_xyz[int(np.argmax(clearance))])
            start_definition = "Maximum atom-surface-clearance sample of the measured fallback cast"
        record = {"start": list(start) if start is not None else None, "start_definition": start_definition,
                  "obstacles": obstacles}
        try:
            tunnels = find_tunnels(frame, start=start, include_hydrogen=include_hydrogen, include_hetero=include_hetero)
        except ValueError as exc:
            record.update(status="not_completed", reason=str(exc))
        else:
            record.update(status="completed", tunnel_count=len(tunnels))
        records["tunnels"] = record
    if not include_network:
        records["network"] = {"status": "skipped_by_request"}
    elif region_xyz is None:
        records["network"] = {"status": "not_completed", "reason": no_region}
    else:
        from .models import ChannelPoint, VoidComponent
        from .networks import build_cavity_network
        volume = float(len(region_xyz) * cell)
        points = tuple(ChannelPoint(i, tuple(map(float, p)), 0.0, 0.0, float(p[2])) for i, p in enumerate(region_xyz))
        region = VoidComponent(1, "measured_void", tuple(map(float, region_xyz.mean(0))), volume, points)
        coords = tuple(p.position for p in points)
        contacts = annotate_residues(frame, coords, identify_bottleneck=False)
        network = build_cavity_network(frame, region, contacts=contacts, residue_groups=residue_groups)
        if residue_groups is not None:
            from .connectivity import lining_connectivity
            connectivity = lining_connectivity(network)
        records["network"] = {"status": "completed", "node_count": len(network.nodes), "edge_count": len(network.edges),
                              "region_source": None, "region_sample_count": len(points),
                              "definition": "Standard publish network; the region node is every measured fallback-cast sample centre (point-to-atom-surface contacts), not a profile"}
    return records, cavities, tunnels, network, connectivity


@radii_option
def compute_static_workflow(
    frame,
    *,
    axis="auto",
    origin=None,
    section_spacing: float = 0.5,
    enclosure_radius: float | None = None,
    include_hydrogen: bool = False,
    include_hetero: bool = True,
    samples: int = 101,
    lateral_exits: bool = False,
    exit_bulk_radius: float = 6.0,
    exit_spacing: float = 0.5,
    search_radius: float = 6.0,
    refinement_steps: int = 4,
    probe_radius: float = 0.0,
    include_cavities: bool = True,
    include_tunnels: bool = True,
    include_network: bool = True,
    cast_mode: str = "auto",
    cast_outer_radius: float = 6.0,
    cast_enclosure_fraction: float = .9,
    cast_extension: float = 2.0,
    cast_spacing: float = 0.5,
    cast_min_radius: float = 0.0,
    cast_max_components: int = 1,
    cast_selection: str = "dominant",
    cast_core_fraction: float = .08,
    cast_max_grid_points: int = 16_000_000,
    cast_focus_points: Sequence | None = None,
    cast_focus_radius: float = 8.0,
    cast_min_component_volume: float = 0.0,
    cast_max_component_volume: float | None = None,
    residue_groups: dict | None = None,
    entry_end: str = "auto",
    lining_cutoff: float = 3.3,
    fingerprint: bool = False,
    radii: RadiusSet | str | None = None,
) -> StaticWorkflowResult:
    """Measure everything ``crevice report`` reports, without writing any file.

    Parameters are those of :func:`crevice.figures.write_static_publication_bundle`
    that affect measurements (the display-only ones belong to
    :func:`export_static_workflow`, except ``cast_extension``, which adds the
    display context recorded in the profile).

    Parameters
    ----------
    frame : StructureFrame
        Parsed structure.
    axis, origin, section_spacing, enclosure_radius, include_hydrogen, include_hetero : object
        Channel profile and atom selection (see :func:`crevice.analysis.pore_profile`).
    samples, lateral_exits, exit_bulk_radius, exit_spacing, search_radius, refinement_steps, probe_radius : object
        Channel profile settings.
    include_cavities, include_tunnels, include_network : bool, default True
        Run the optional cavity search, tunnel search and residue network.
    cast_mode, cast_outer_radius, cast_enclosure_fraction, cast_extension, cast_spacing : object
        Cast settings (``cast_extension`` is the display-only mouth continuation, Å).
    cast_min_radius, cast_max_components, cast_selection, cast_core_fraction, cast_max_grid_points : object
        Cast settings.
    cast_focus_points, cast_focus_radius, cast_min_component_volume, cast_max_component_volume : object
        Cast crops.
    residue_groups : dict of str to str, optional
        Residue groups for the network connectivity report.
    entry_end : {"auto", "start", "end"}, default "auto"
        Which channel end is labelled ENTRY.
    lining_cutoff : float, default 3.3
        Lining-residue distance (Å).
    fingerprint : bool, default False
        Also compute the pore fingerprint of a resolved channel.
    radii : RadiusSet, str or None, optional
        Atomic radius set; ``None`` uses the set in effect.

    Returns
    -------
    StaticWorkflowResult

    Raises
    ------
    ChannelResolutionError
        When no channel resolves and a crop or non-auto mode forbids the fallback.
    ValueError
        For invalid settings.
    """
    from .analysis import annotate_residues, detect_cavities, pore_profile
    from .networks import build_cavity_network
    from .sections import LEGACY_ENCLOSURE_RADIUS, ChannelResolutionError
    from .tunnels import find_tunnels
    from .voids import void_cast
    if entry_end not in (None, "auto", "start", "end"):
        raise ValueError("entry_end must be 'auto', 'start' or 'end'")
    if not math.isfinite(lining_cutoff) or lining_cutoff <= 0:
        raise ValueError("lining_cutoff must be finite and positive")
    # Rolling-probe casts need a number; the automatic choice applies to the channel profile only.
    cast_probe = LEGACY_ENCLOSURE_RADIUS if enclosure_radius is None else enclosure_radius
    settings = {"include_hydrogen": include_hydrogen, "include_hetero": include_hetero, "cast_probe": cast_probe,
                "entry_end": entry_end, "lining_cutoff": lining_cutoff, "exit_bulk_radius": exit_bulk_radius,
                "residue_groups": residue_groups, "cast_mode": cast_mode, "focus": cast_focus_points,
                "fingerprint_requested": bool(fingerprint)}
    rolling_options = dict(spacing=cast_spacing, probe_radius=cast_probe, outer_radius=cast_outer_radius,
                           enclosure_fraction=cast_enclosure_fraction, min_component_volume=cast_min_component_volume,
                           max_components=cast_max_components, max_grid_points=cast_max_grid_points,
                           selection_mode=cast_selection, core_fraction=cast_core_fraction)
    if cast_mode == 'rolling':
        if cast_focus_points is not None or cast_max_component_volume is not None or cast_min_radius != 0:
            raise ValueError('Rolling publication mode does not support focus, maximum-volume or minimum-swept-clearance crops; use the explicit rolling probe and component settings')
        rolling = compute_rolling_cast(frame, selection='all' if include_hetero else 'protein', **rolling_options)
        return StaticWorkflowResult("rolling", frame, rolling=rolling, settings=settings)
    if not math.isfinite(cast_extension) or cast_extension < 0:
        raise ValueError("cast_extension must be finite and non-negative")
    frame = replace(frame, atoms=frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero))
    try:
        profile = pore_profile(frame, axis=axis, origin=origin, section_spacing=section_spacing,
                               enclosure_radius=enclosure_radius, include_hydrogen=include_hydrogen,
                               include_hetero=include_hetero, samples=samples, search_radius=search_radius,
                               refinement_steps=refinement_steps, probe_radius=probe_radius,
                               lateral_exits=lateral_exits, exit_bulk_radius=exit_bulk_radius, exit_spacing=exit_spacing)
    except ChannelResolutionError as exc:
        if (cast_mode != 'auto' or cast_focus_points is not None or cast_max_component_volume is not None
                or cast_min_radius != 0 or include_hydrogen):
            raise
        # General voids remain useful when no unique through-channel exists. Preserve the
        # explicit selected obstacle set and never fabricate a radius profile.
        rolling = compute_rolling_cast(frame, selection='all', **rolling_options)
        records, cavities, tunnels, network, connectivity = _unresolved_extras(
            frame, rolling.cast, include_cavities=include_cavities, include_tunnels=include_tunnels,
            include_network=include_network, include_hydrogen=include_hydrogen, include_hetero=include_hetero,
            residue_groups=residue_groups)
        return StaticWorkflowResult("unresolved", frame, rolling=rolling, unresolved_reason=str(exc),
                                    probe_selection=getattr(exc, 'probe_selection', None), cavities=cavities,
                                    tunnels=tunnels, network=network, connectivity=connectivity, extras=records,
                                    settings=settings)
    if profile.method == "axial-connected" and cast_extension:
        from .geometry import add, scale
        from .spatial import SpatialIndex
        spatial = SpatialIndex(frame.atoms)
        context = []
        for point, sign in ((profile.points[0], -1), (profile.points[-1], 1)):
            segment = []
            count = max(2, int(math.ceil(cast_extension / 0.25)))
            for k in range(count+1):
                dt = sign * cast_extension * k/count
                position = add(point.position, scale(profile.axis_direction, dt))
                raw = spatial.nearest_surface(position)[0]
                segment.append({"t": point.t+dt, "radius": max(0.0, raw-probe_radius)})
            context.append(segment)
        profile = replace(profile, metadata={**profile.metadata, "display_context": context,
                          "display_context_definition": "clearance along axial exit; outside channel statistics"})
    contacts = annotate_residues(frame, profile.points)
    cast = void_cast(frame, mode=cast_mode, axis=axis, spacing=cast_spacing, min_radius=cast_min_radius,
                     max_components=cast_max_components, max_grid_points=cast_max_grid_points,
                     focus_points=cast_focus_points, focus_radius=cast_focus_radius,
                     min_component_volume=cast_min_component_volume, max_component_volume=cast_max_component_volume,
                     profile=profile)
    exit_casts = ()
    if (profile.metadata.get("exits") and cast.mode == "channel"
            and any(e.get("type") == "lateral" for e in profile.metadata["exits"].values())):
        # Lateral exit legs as probe-swept casts, reported separately; the axial cast is unchanged.
        from .channel_exits import EXIT_CAST_DEFINITION, exit_cast_summary, lateral_exit_casts
        exit_casts = lateral_exit_casts(frame, profile, cast)
        if exit_casts:
            cast = replace(cast, metadata={**cast.metadata, "lateral_exit_casts": exit_cast_summary(exit_casts),
                                           "lateral_exit_cast_definition": EXIT_CAST_DEFINITION,
                                           "lateral_exit_casts_note": "separate volumes; not part of total_volume"})
    segments, lining, segment_rows = None, [], []
    if cast.mode == "channel" and cast.metadata.get("export_mode") == "connected_section_fill" and not cast_focus_points:
        # ENTRY / LUMEN / EXIT: the lumen is the cast above, unchanged; the vestibules beyond each
        # mouth are separate segments.
        from .cast_segments import (LABEL_RULE_DEFINITION, LINING_DEFINITION, SEGMENT_DEFINITION, lining_residues,
                                    segment_casts, segment_roles, segment_summary)
        segment_roles(profile, entry_end)
        segments = segment_casts(frame, profile, cast, entry_end=entry_end, bulk_radius=exit_bulk_radius,
                                 lateral_casts=exit_casts)
        segment_rows = segment_summary(cast, segments)
        lining = lining_residues(frame, cast, segments, cutoff=lining_cutoff)
        cast = replace(cast, metadata={**cast.metadata, "cast_segments": segment_rows,
                                       "cast_segment_definition": SEGMENT_DEFINITION,
                                       "segment_label_rule": {"entry_end": entry_end or "auto",
                                                              "rule": segment_roles(profile, entry_end)[1],
                                                              "definition": LABEL_RULE_DEFINITION},
                                       "segment_bulk_radius_A": exit_bulk_radius,
                                       "segments_total_volume_A3": sum(r["volume_A3"] for r in segment_rows),
                                       "segments_total_note": "lumen + every entry and exit segment; "
                                                              "total_volume stays the lumen only",
                                       "lining_cutoff_A": lining_cutoff, "lining_definition": LINING_DEFINITION})
    if cast_focus_points is not None or cast_mode in {"cavity", "all"}:
        from .models import VoidComponent
        cast_points = tuple(p for c in cast.components for p in c.points)
        region = VoidComponent(0, "selected_void", (0.0, 0.0, 0.0), cast.total_volume, cast_points)
        contacts = annotate_residues(frame, cast_points, identify_bottleneck=False)
    else:
        region = profile
    fp = None
    if fingerprint:
        from .fingerprint import pore_fingerprint
        fp = pore_fingerprint(frame, profile, cast=cast if segments is not None else None, entry_end=entry_end,
                              bulk_radius=exit_bulk_radius, include_hetero=include_hetero)
    cavities = detect_cavities(frame) if include_cavities else None
    tunnels = find_tunnels(frame) if include_tunnels else None
    network = connectivity = None
    if include_network:
        network = build_cavity_network(frame, region, contacts=contacts, residue_groups=residue_groups)
        if residue_groups is not None:
            from .connectivity import lining_connectivity
            connectivity = lining_connectivity(network)
    return StaticWorkflowResult("resolved", frame, profile=profile, contacts=contacts, cast=cast, exit_casts=exit_casts,
                                segments=segments, segment_rows=segment_rows, lining=lining, region=region,
                                cavities=cavities, tunnels=tunnels, network=network, connectivity=connectivity,
                                fingerprint=fp, settings=settings)


def _write_cavities(cavities, root, prefix, manifest, *, dpi):
    from .figures import plot_cavity_summary
    from .io import write_cavities_csv, write_cavities_json, write_cavities_pdb
    from .outputs import wants
    manifest["cavities_json"] = str(root / f"{prefix}_cavities.json")
    manifest["cavities_csv"] = str(root / f"{prefix}_cavities.csv")
    manifest["cavities_pdb"] = str(root / f"{prefix}_cavities.pdb")
    manifest["cavity_summary_png"] = str(root / f"{prefix}_cavity_summary.png")
    write_cavities_json(cavities, manifest["cavities_json"])
    write_cavities_csv(cavities, manifest["cavities_csv"])
    if wants(manifest["cavities_pdb"]):
        write_cavities_pdb(cavities, manifest["cavities_pdb"])
    plot_cavity_summary(cavities, manifest["cavity_summary_png"], dpi=dpi)


def _write_tunnels(tunnels, root, prefix, manifest):
    from .io import write_tunnel_points_csv, write_tunnels_csv, write_tunnels_json, write_tunnels_pdb
    from .outputs import wants
    manifest["tunnels_json"] = str(root / f"{prefix}_tunnels.json")
    manifest["tunnels_csv"] = str(root / f"{prefix}_tunnels.csv")
    manifest["tunnel_points_csv"] = str(root / f"{prefix}_tunnel_points.csv")
    manifest["tunnels_pdb"] = str(root / f"{prefix}_tunnels.pdb")
    write_tunnels_json(tunnels, manifest["tunnels_json"])
    write_tunnels_csv(tunnels, manifest["tunnels_csv"])
    write_tunnel_points_csv(tunnels, manifest["tunnel_points_csv"])
    if wants(manifest["tunnels_pdb"]):
        write_tunnels_pdb(tunnels, manifest["tunnels_pdb"])


def _write_network(network, connectivity, root, prefix, manifest, *, dpi):
    from .figures import _write_json, plot_network_chord, plot_network_summary
    from .io import write_connectivity_csv, write_network_csv
    from .networks import network_metrics
    if connectivity is not None:
        from .connectivity import plot_lining_connectivity
        manifest["connectivity_json"] = str(root / f"{prefix}_connectivity.json")
        manifest["connectivity_csv"] = str(root / f"{prefix}_connectivity.csv")
        manifest["connectivity_groups_csv"] = str(root / f"{prefix}_connectivity_groups.csv")
        manifest["connectivity_png"] = str(root / f"{prefix}_connectivity.png")
        _write_json(connectivity, manifest["connectivity_json"])
        write_connectivity_csv(connectivity, manifest["connectivity_csv"], manifest["connectivity_groups_csv"])
        plot_lining_connectivity(network, manifest["connectivity_png"], dpi=dpi)
    manifest["network_json"] = str(root / f"{prefix}_network.json")
    manifest["network_nodes_csv"] = str(root / f"{prefix}_network_nodes.csv")
    manifest["network_edges_csv"] = str(root / f"{prefix}_network_edges.csv")
    manifest["network_summary_png"] = str(root / f"{prefix}_network_summary.png")
    manifest["network_chord_png"] = str(root / f"{prefix}_network_chord.png")
    _write_json({"network": network.to_dict(), "metrics": network_metrics(network)}, manifest["network_json"])
    write_network_csv(network, manifest["network_nodes_csv"], manifest["network_edges_csv"])
    plot_network_summary(network, manifest["network_summary_png"], dpi=dpi)
    plot_network_chord(network, manifest["network_chord_png"], dpi=dpi)


@annotate_option
def export_static_workflow(result: StaticWorkflowResult, *, structure_path, output_dir, prefix: str = "crevice",
                           dpi: int = 200, residue_label_every: int = 5, cast_extension: float = 2.0,
                           surface_smoothing: float | None = None, smooth: float | None = None, outputs=None,
                           annotate: bool | None = None) -> dict:
    """Write a :class:`StaticWorkflowResult` as the ``crevice report`` bundle.

    The files and manifest are those of
    :func:`crevice.figures.write_static_publication_bundle`, filtered by the
    output plan (``outputs`` or the active plan; everything without one).

    Parameters
    ----------
    result : StaticWorkflowResult
        From :func:`compute_static_workflow`.
    structure_path : str or Path
        Source file, copied into the bundle (file kind ``structures``).
    output_dir : str or Path
        Bundle directory (created).
    prefix : str, default "crevice"
        File-name prefix.
    dpi : int, default 200
        Figure resolution.
    residue_label_every : int, default 5
        Landmark spacing of the annotated radius figure.
    cast_extension : float, default 2.0
        Display-only cast continuation beyond each mouth (Å); use the value
        given to :func:`compute_static_workflow`.
    surface_smoothing, smooth : float, optional
        Display-surface smoothing (Å); see :func:`crevice.volume_export.resolve_display_smoothing`.
    outputs : None, str, sequence of str or OutputPlan, optional
        Output selection for this call (policy of ``crevice report``, see
        :mod:`crevice.outputs`); ``None`` keeps the active plan.
    annotate : bool, optional
        Draw explanatory text in the figures; ``None`` inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        Output key to path, including ``manifest_json``.
    """
    from .outputs import outputs_for
    with outputs_for("publish", outputs):
        return _export_static(result, structure_path=structure_path, output_dir=output_dir, prefix=prefix, dpi=dpi,
                              residue_label_every=residue_label_every, cast_extension=cast_extension,
                              surface_smoothing=surface_smoothing, smooth=smooth)


def _export_static(result, *, structure_path, output_dir, prefix, dpi, residue_label_every, cast_extension,
                   surface_smoothing, smooth):
    from .figures import (_copy_structure_for_bundle, _pdb_resi_for_min_radius, _viewer_structure_for_bundle,
                          _write_json, plot_profile_radius, plot_profile_radius_with_residues, plot_residue_contacts,
                          write_publication_pymol_script)
    from .io import (write_profile_csv, write_profile_json, write_residue_contacts_csv, write_void_cast_csv,
                     write_void_cast_json, write_void_cast_pdb)
    from .outputs import drop_blocked, note, wants
    from .volume_export import resolve_display_smoothing
    resolve_display_smoothing(smooth, surface_smoothing)
    settings = result.settings
    if result.status == "rolling":
        return export_rolling_cast(result.rolling, structure_path=structure_path, output_dir=output_dir, prefix=prefix,
                                   surface_smoothing=surface_smoothing, smooth=smooth)
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    frame = result.frame
    if result.status == "unresolved":
        # The rolling export writes the input copy and viewer PDB of its own atom selection.
        manifest = export_rolling_cast(result.rolling, structure_path=structure_path, output_dir=root, prefix=prefix,
                                       surface_smoothing=surface_smoothing, smooth=smooth)
        saved = json.loads(Path(manifest['manifest_json']).read_text(encoding="utf-8"))
        saved['single_channel_unresolved_reason'] = result.unresolved_reason
        saved['auto_fallback'] = '3D rolling-probe cast; no single-channel profile'
        saved['publication_input_selection'] = {'include_hydrogen': settings['include_hydrogen'],
                                                'include_hetero': settings['include_hetero']}
        fallback = json.loads(json.dumps(result.extras))
        if result.cavities is not None:
            _write_cavities(result.cavities, root, prefix, manifest, dpi=dpi)
        if result.tunnels is not None:
            _write_tunnels(result.tunnels, root, prefix, manifest)
        if result.network is not None:
            _write_network(result.network, result.connectivity, root, prefix, manifest, dpi=dpi)
            fallback["network"]["region_source"] = str(manifest.get("volume_dx"))
        status = {'status': 'unresolved', 'reason': result.unresolved_reason,
                  'interpretation': 'No unique through-channel resolved; no pore-radius profile or profile-derived features assigned.',
                  'fallback_analyses': fallback}
        selection = result.probe_selection
        if selection is not None:
            status['enclosure_probe_selection'] = selection
            saved['enclosure_probe'] = {'mode': 'auto', 'chosen_A': None, 'reason': selection['reason'],
                                        'fallback_cast_probe_A': settings['cast_probe']}
        if settings.get("fingerprint_requested"):
            note('fingerprint', 'not written: no through-channel profile resolved (a fingerprint needs one)')
        manifest['profile_status_json'] = str(root / f'{prefix}_profile_status.json')
        _write_json(status, manifest['profile_status_json'])
        saved['profile_status'] = {k: status[k] for k in ('status', 'reason', 'interpretation')}
        if selection is not None:
            saved['profile_status']['enclosure_probe'] = saved.pop('enclosure_probe')
        saved['fallback_analyses'] = fallback
        saved.update(radii_fields())
        drop_blocked(manifest)
        saved['files'].update({k: v for k, v in manifest.items() if k != 'manifest_json'})
        _write_json(saved, manifest['manifest_json'])
        return manifest

    profile, cast, contacts = result.profile, result.cast, result.contacts
    exit_casts, segments, segment_rows, lining = result.exit_casts, result.segments, result.segment_rows, result.lining
    manifest: dict = {}
    original_copy = _copy_structure_for_bundle(structure_path, root, prefix)
    structure_copy = _viewer_structure_for_bundle(frame, root / f"{prefix}_viewer.pdb")
    manifest["structure_file"] = str(structure_copy)
    manifest["original_structure_file"] = str(original_copy)
    manifest["viewer_identity_json"] = str(structure_copy.with_suffix(".identities.json"))
    manifest["profile_json"] = str(root / f"{prefix}_profile.json")
    manifest["profile_csv"] = str(root / f"{prefix}_profile.csv")
    manifest["pore_cast_pdb"] = str(root / f"{prefix}_pore_cast.pdb")
    manifest["pore_cast_pml"] = str(root / f"{prefix}_pore_cast.pml")
    manifest["void_cast_json"] = str(root / f"{prefix}_void_cast.json")
    manifest["pymol_render_png"] = str(root / f"{prefix}_pore_cast.png")
    manifest["profile_radius_png"] = str(root / f"{prefix}_profile_radius.png")
    manifest["profile_radius_annotated_png"] = str(root / f"{prefix}_profile_radius_annotated.png")
    manifest["void_cast_csv"] = str(root / f"{prefix}_void_cast.csv")
    manifest["residue_contacts_csv"] = str(root / f"{prefix}_residue_contacts.csv")
    manifest["residue_contacts_png"] = str(root / f"{prefix}_residue_contacts.png")
    write_profile_json(profile, manifest["profile_json"])
    write_profile_csv(profile, manifest["profile_csv"])
    write_void_cast_json(cast, manifest["void_cast_json"])
    write_void_cast_csv(cast, manifest["void_cast_csv"], extra_components=exit_casts)
    if wants(manifest["pore_cast_pdb"]):
        write_void_cast_pdb(cast, manifest["pore_cast_pdb"])
    write_residue_contacts_csv(contacts, manifest["residue_contacts_csv"])
    plot_profile_radius(profile, manifest["profile_radius_png"], dpi=dpi)
    plot_profile_radius_with_residues(profile, manifest["profile_radius_annotated_png"], dpi=dpi,
                                      label_every=residue_label_every)
    plot_residue_contacts(contacts, manifest["residue_contacts_png"], dpi=dpi)
    if result.fingerprint is not None:
        from .fingerprint import write_fingerprint_bundle
        manifest.update(write_fingerprint_bundle(result.fingerprint, root, prefix, dpi=dpi))
    sphere_scene = wants(manifest["pore_cast_pml"])
    if not sphere_scene:
        del manifest["pymol_render_png"]   # named only inside the sphere scene
    if sphere_scene:
        write_publication_pymol_script(
            structure_path=structure_copy, pore_cast_path=manifest["pore_cast_pdb"],
            script_path=manifest["pore_cast_pml"], profile=profile, output_png=manifest["pymol_render_png"], dpi=dpi,
            structure_local_name=Path(structure_copy).name, pore_cast_local_name=Path(manifest["pore_cast_pdb"]).name,
            output_png_local_name=Path(manifest["pymol_render_png"]).name,
            bottleneck_resi=_pdb_resi_for_min_radius(cast.points))
    from .volume_export import write_volume_viewer_bundle
    if cast.components:
        manifest.update(write_volume_viewer_bundle(cast, structure_path=structure_copy, output_dir=root, prefix=prefix,
                                                   frame=frame, profile=profile, cast_extension=cast_extension,
                                                   surface_smoothing=surface_smoothing, smooth=smooth,
                                                   exit_casts=exit_casts, segments=segments, lining=lining or None))
        if profile.method == "axial-connected" and sphere_scene:
            # The primary pore-cast scene must show the measured full grid; the dummy-atom
            # representation is kept separately for compatibility.
            legacy = root / f"{prefix}_pore_spheres.pml"
            if wants(legacy):
                shutil.copyfile(manifest["pore_cast_pml"], legacy)
                manifest["pore_spheres_pml"] = str(legacy)
            shutil.copyfile(manifest["volume_pml"], manifest["pore_cast_pml"])
            manifest["pymol_render_png"] = str(root / f"{prefix}_pymol.png")
        for extension in ['vmd', 'tcl', 'cxc']:
            target = root/f'{prefix}_pore_cast.{extension}'
            if wants(target):
                shutil.copyfile(manifest[f'volume_{extension}'], target)
                manifest[f'pore_cast_{extension}'] = str(target)
    if result.cavities is not None:
        _write_cavities(result.cavities, root, prefix, manifest, dpi=dpi)
    if result.tunnels is not None:
        _write_tunnels(result.tunnels, root, prefix, manifest)
    if result.network is not None:
        _write_network(result.network, result.connectivity, root, prefix, manifest, dpi=dpi)
    if segments is not None:
        from .cast_segments import add_lining_columns, write_lining_csv, write_segments_csv
        manifest["cast_segments_csv"] = str(root / f"{prefix}_cast_segments.csv")
        manifest["lining_residues_csv"] = str(root / f"{prefix}_lining_residues.csv")
        write_segments_csv(segment_rows, manifest["cast_segments_csv"])
        write_lining_csv(lining, manifest["lining_residues_csv"])
        for key in ("residue_contacts_csv", "network_nodes_csv"):
            if key in manifest and Path(manifest[key]).is_file():
                add_lining_columns(manifest[key], lining)
    manifest_path = root / f"{prefix}_manifest.json"
    drop_blocked(manifest)
    _write_json({"files": manifest, "dpi": dpi, "prefix": prefix, "void_cast": cast.to_dict(), **radii_fields(),
                 "profile_status": _resolved_status(result, manifest)}, manifest_path)
    manifest["manifest_json"] = str(manifest_path)
    return manifest


def _resolved_status(result, manifest) -> dict:
    profile, cast = result.profile, result.cast
    exit_casts, segments, segment_rows, lining = result.exit_casts, result.segments, result.segment_rows, result.lining
    status = {"status": "resolved", "method": profile.method, "profile_json": manifest["profile_json"]}
    if "path_type" in profile.metadata:
        status.update({"path_type": profile.metadata["path_type"],
                       "exit_types": {k: e["type"] for k, e in profile.metadata["exits"].items()},
                       "exit_leg_counts": {k: e.get("leg_count", 0) for k, e in profile.metadata["exits"].items()}})
    if exit_casts:
        status.update({"lateral_exit_cast_volumes_A3": {f"{c.metadata['end']}_{c.metadata['leg']}": c.volume
                                                        for c in exit_casts},
                       "lateral_exit_cast_note": "separate from void_cast.total_volume (axial channel only)"})
    if segments is not None:
        status.update({"cast_segment_volumes_A3": {r["segment"]: r["volume_A3"] for r in segment_rows},
                       "segment_label_rule": cast.metadata["segment_label_rule"]["rule"],
                       "lining_residue_counts": {r["segment"]: sum(1 for x in lining if x["segment"] == r["segment"])
                                                 for r in segment_rows},
                       "cast_segment_note": "lumen = void_cast.total_volume; entry/exit volumes separate"})
    if "enclosure_probe_selection" in profile.metadata:
        selection = profile.metadata["enclosure_probe_selection"]
        status["enclosure_probe"] = {"mode": "auto", "chosen_A": selection["chosen_A"], "reason": selection["reason"],
                                     "candidates_tried": [r["probe_A"] for r in selection["candidates"]
                                                          if r["status"] != "not_tried"]}
    return status

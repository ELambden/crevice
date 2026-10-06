"""Standard-analysis hydration integration; solvent stays outside pore geometry."""
from pathlib import Path
import json

from .presentation import annotate_option
from .radii import RadiusSet, radii_option


def hydration_dependencies():
    """Return required hydration libraries that cannot be imported (empty on a correct install)."""
    from importlib.util import find_spec
    return [name for name in ['numpy','scipy','matplotlib'] if find_spec(name) is None]


@radii_option
@annotate_option
def static_hydration_bundle(frame,output_dir,*,prefix='crevice',source_path=None,volume_path=None,residue_ids=None,roles=None,
                            cutoff=3.5,sasa_probe=1.4,sasa_samples=256,dpi=240,analysis_options=None,annotate=None,
                            radii: RadiusSet | str | None = None):
    """Measure and write hydration for one structure (``crevice hydration`` without a trajectory).

    The residue focus is ``residue_ids``, else the boundary residues and
    nonlocal partners of ``volume_path`` (see
    :func:`crevice.residue_evidence.boundary_residue_evidence`), else every
    selected protein residue. Deposited waters with positive occupancy are
    counted once; a structure without explicit water reports hydration as
    unobserved, never as dry. For GRO input the full system and its box are
    reloaded so solvent is available. Optional density, typed interactions and
    the coordinated scenes, evidence table and overview figures follow
    ``analysis_options`` (see :func:`analysis_options_from_args`).

    Parameters
    ----------
    frame : StructureFrame
        Parsed structure.
    output_dir : str or Path
        Output directory.
    prefix : str, default "crevice"
        File-name prefix.
    source_path : str or Path, optional
        Source file; a ``.gro`` path triggers reloading the full solvated system.
    volume_path : str or Path, optional
        Measured binary cast DX defining the residue focus and scene cavity.
    residue_ids : sequence of str, optional
        Explicit focus residue IDs; an empty list writes an ``unavailable`` report.
    roles : dict of str to str, optional
        Residue role labels.
    cutoff : float, default 3.5
        Water oxygen to protein heavy-atom contact distance (Å).
    sasa_probe : float, default 1.4
        SASA probe radius (Å).
    sasa_samples : int, default 256
        Quadrature points per atom.
    dpi : int, default 240
        Figure resolution.
    analysis_options : dict, optional
        Density spacing/smoothing/level and skip flags (see
        :func:`analysis_options_from_args`).
    annotate : bool, optional
        Draw titles and scene text (see
        :func:`crevice.hydration_export.write_hydration_bundle`). ``None``
        (default) inherits the surrounding setting.
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
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_hydration.* and optional density, interaction and analysis files : bundle
        See :func:`crevice.hydration_export.write_hydration_bundle`.

    Colours and representation
    --------------------------
    As :func:`crevice.hydration_export.write_hydration_bundle` and
    :func:`crevice.analysis_scene_templates.write_scene_scripts`.
    """
    from .hydration import static_hydration
    from .hydration_export import write_hydration_bundle
    from .residue_evidence import read_binary_dx,boundary_residue_evidence
    box=None; evidence=None
    if volume_path:
        grid,origin,deltas=read_binary_dx(volume_path);evidence=boundary_residue_evidence(frame,grid,origin,deltas)
        if residue_ids is None:residue_ids=[r['residue'] for r in evidence['residues'] if r['role']!='other']
        if roles is None:roles={r['residue']:r['role'] for r in evidence['residues']}
    if source_path and Path(source_path).suffix.lower()=='.gro':
        from .parser import load_structure_report
        from .trajectory import _open_universe
        full,_=load_structure_report(source_path,model_index=frame.model_index,md_selection='all')
        universe=_open_universe(source_path);universe.trajectory[frame.model_index];box=universe.dimensions
        frame=full
    if residue_ids is not None and not residue_ids:
        root=Path(output_dir);root.mkdir(parents=True,exist_ok=True);p=root/(prefix+'_hydration.json')
        p.write_text(json.dumps({'status':'unavailable_no_identified_residues','residues':[],'interpretation':'No pore/cavity residue focus was available; not evidence of absent hydration'},indent=2)+'\n')
        return {'hydration_json':str(p)}
    result=static_hydration(frame,residue_ids=residue_ids,box=box,cutoff=cutoff,sasa_probe=sasa_probe,sasa_samples=sasa_samples)
    if evidence is not None:result['boundary_area_by_residue_A2']={r['residue']:r['boundary_area_A2'] for r in evidence['residues']}
    from .hydration_extensions import attach_static_extensions
    attach_static_extensions(result,frame,box=box,volume_path=volume_path,options=analysis_options)
    return write_hydration_bundle(result,output_dir,prefix=prefix,roles=roles,dpi=dpi)


def append_static_hydration(args,frame,source_path,files,*,residue_ids=None,roles=None):
    """Add the standard hydration bundle to a CLI command's outputs.

    Used by ``cast``, ``publish``, ``analyze``, ``residue-evidence`` and
    ``trajectory``: unless ``--skip-hydration`` is set, runs
    :func:`static_hydration_bundle` with the command's hydration options and
    merges its files into ``files`` and the command manifest. Figure and scene
    text follows the current annotation setting (CLI ``--annotate``).

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    frame : StructureFrame
        The analysed structure.
    source_path : str or Path
        Structure file.
    files : dict of str to str
        Command outputs; updated in place.
    residue_ids : sequence of str, optional
        Residue focus.
    roles : dict of str to str, optional
        Residue roles.

    Returns
    -------
    dict of str to str
        The hydration outputs.
    """
    if getattr(args,'skip_hydration',False):return {}
    prefix=getattr(args,'prefix',None) or Path(source_path).stem
    missing=hydration_dependencies()
    if missing:
        root=Path(args.out_dir);root.mkdir(parents=True,exist_ok=True);p=root/(prefix+'_hydration.json')
        p.write_text(json.dumps({'status':'unavailable_missing_required_dependencies','missing':missing,'install':'reinstall CREVICE (pip install --force-reinstall crevice)','other_analyses_retained':True},indent=2)+'\n')
        extra={'hydration_json':str(p)}
    else:
        extra=static_hydration_bundle(frame,args.out_dir,prefix=prefix,source_path=source_path,volume_path=files.get('volume_dx'),
          residue_ids=residue_ids,roles=roles,cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,dpi=args.dpi,analysis_options=analysis_options_from_args(args))
    files.update(extra)
    manifest=files.get('manifest_json')
    if manifest:
        p=Path(manifest);data=json.loads(p.read_text());data['files'].update(extra);p.write_text(json.dumps(data,indent=2)+'\n')
    return extra


def analysis_options_from_args(args):
    """Collect the hydration-analysis display options from parsed CLI arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments of ``hydration``, ``publish``, ``cavity-trajectory`` or
        ``region-trajectory``; missing attributes take their defaults.

    Returns
    -------
    dict
        ``density_spacing`` (``--water-density-spacing``, default 0.5 Å),
        ``density_smoothing`` (``--water-density-smoothing``, default 0.5 Å, the
        Gaussian standard deviation of the display density map),
        ``density_level`` (``--water-density-level``, default 0.05 oxygens per
        Å³, the contour level of the density surface), ``skip_density``
        (``--skip-water-density``), ``skip_interactions``
        (``--skip-interactions``) and ``skip_views`` (``--skip-analysis-views``).
    """
    return {'density_spacing':getattr(args,'water_density_spacing',.5),
            'density_smoothing':getattr(args,'water_density_smoothing',.5),
            'density_level':getattr(args,'water_density_level',.05),
            'skip_density':getattr(args,'skip_water_density',False),
            'skip_interactions':getattr(args,'skip_interactions',False),
            'skip_views':getattr(args,'skip_analysis_views',False)}

"""Standard-analysis hydration integration; solvent stays outside pore geometry."""
from pathlib import Path
from ._portable_paths import dumps_record
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
                            cutoff=3.5,sasa_probe=1.4,sasa_samples=256,dpi=200,analysis_options=None,annotate=None,
                            radii: RadiusSet | str | None = None):
    """Measure and write hydration for one structure (``crevice trajectory --water`` without a trajectory).

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
    dpi : int, default 200
        Figure resolution. The default suits screens; use 400 or more for print.
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
        p.write_text(dumps_record({'status':'unavailable_no_identified_residues','residues':[],'interpretation':'No pore/cavity residue focus was available; not evidence of absent hydration'},p,indent=2)+'\n', encoding="utf-8", newline="\n")
        return {'hydration_json':str(p)}
    result=static_hydration(frame,residue_ids=residue_ids,box=box,cutoff=cutoff,sasa_probe=sasa_probe,sasa_samples=sasa_samples)
    if evidence is not None:result['boundary_area_by_residue_A2']={r['residue']:r['boundary_area_A2'] for r in evidence['residues']}
    from .hydration_extensions import attach_static_extensions
    attach_static_extensions(result,frame,box=box,volume_path=volume_path,options=analysis_options)
    return write_hydration_bundle(result,output_dir,prefix=prefix,roles=roles,dpi=dpi)


def append_static_hydration(args,frame,source_path,files,*,residue_ids=None,roles=None):
    """Add the standard hydration bundle to a CLI command's outputs.

    Used by ``cavities --cast``, ``report``, ``report --analyze`` and ``residues --evidence``: when
    the active output plan includes ``hydration`` (``--outputs hydration`` or
    ``--all``; without a plan, unless ``--skip-hydration``), runs
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
    from .outputs import active_plan
    plan=active_plan()
    if plan is not None and not plan.runs('hydration'):
        return append_residue_interactions(args,frame,source_path,files,residue_ids=residue_ids) if plan.runs('interactions') else {}
    if plan is None and getattr(args,'skip_hydration',False):return {}
    prefix=getattr(args,'prefix',None) or Path(source_path).stem
    missing=hydration_dependencies()
    if missing:
        root=Path(args.out_dir);root.mkdir(parents=True,exist_ok=True);p=root/(prefix+'_hydration.json')
        p.write_text(dumps_record({'status':'unavailable_missing_required_dependencies','missing':missing,'install':'reinstall CREVICE (pip install --force-reinstall crevice)','other_analyses_retained':True},p,indent=2)+'\n', encoding="utf-8", newline="\n")
        extra={'hydration_json':str(p)}
    else:
        extra=static_hydration_bundle(frame,args.out_dir,prefix=prefix,source_path=source_path,volume_path=files.get('volume_dx'),
          residue_ids=residue_ids,roles=roles,cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,dpi=args.dpi,analysis_options=analysis_options_from_args(args))
    files.update(extra)
    manifest=files.get('manifest_json')
    if manifest:
        p=Path(manifest);data=json.loads(p.read_text(encoding="utf-8"));data['files'].update(extra);p.write_text(dumps_record(data,p,indent=2)+'\n', encoding="utf-8", newline="\n")
    return extra


def append_residue_interactions(args,frame,source_path,files,*,residue_ids=None):
    """Add the typed residue–residue interactions (no water analyses) to a CLI command's outputs.

    Runs :func:`crevice.residue_interactions.static_residue_interactions` with
    the residue focus the hydration analysis would use (``residue_ids``, else
    the boundary residues of the command's ``volume_dx``) and writes
    ``PREFIX_interactions.json`` and the interaction/conformation tables. The
    residue-residue numbers equal those of the hydration workflow; water
    bridges and residue-water hydrogen bonds are not reported.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (``out_dir``, ``prefix``).
    frame : StructureFrame
        The analysed structure.
    source_path : str or Path
        Structure file.
    files : dict of str to str
        Command outputs; updated in place (and the command manifest).
    residue_ids : sequence of str, optional
        Residue focus.

    Returns
    -------
    dict of str to str
        The interaction outputs.
    """
    from .residue_interactions import static_residue_interactions,write_residue_interactions
    prefix=getattr(args,'prefix',None) or Path(source_path).stem
    from .outputs import note
    try:
        if residue_ids is not None and not residue_ids:
            result=None
        else:
            result=static_residue_interactions(frame,source_path=source_path,volume_path=None if residue_ids is not None else files.get('volume_dx'),
                                               residue_ids=residue_ids)
    except ValueError as exc:
        # A default-on analysis: record why it could not run instead of failing the command.
        note('interactions',f'not completed: {exc}')
        return {}
    extra=write_residue_interactions(result,args.out_dir,prefix)
    files.update(extra)
    manifest=files.get('manifest_json')
    if manifest:
        p=Path(manifest);data=json.loads(p.read_text(encoding="utf-8"));data['files'].update(extra);p.write_text(dumps_record(data,p,indent=2)+'\n', encoding="utf-8", newline="\n")
    return extra


def analysis_options_from_args(args):
    """Collect the hydration-analysis display options from parsed CLI arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments of ``hydration``, ``report``, ``trajectory --cavities`` or
        ``trajectory --region``; missing attributes take their defaults.

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

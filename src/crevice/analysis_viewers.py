"""Portable coordinated scenes for measured hydration and residue chemistry."""
from pathlib import Path
from itertools import combinations
import json
import math
import shutil

from .presentation import annotate_option


def hydration_rgb(occupancy):
    """Absolute linear red-to-blue scale; missing observations are neutral grey.

    Used for residue colouring by water-contact occupancy in the unified analysis
    scenes: 0 (never hydrated) is red, 1 (hydrated in every frame) is blue, and
    values between are mixed linearly.

    Parameters
    ----------
    occupancy : float or None
        Water-contact occupancy in [0, 1]; ``None`` = not observed (for example
        no explicit water).

    Returns
    -------
    list of float
        RGB in [0, 1]: ``[1 - occupancy, 0, occupancy]``, or
        ``[0.55, 0.55, 0.55]`` for ``None``.

    Raises
    ------
    ValueError
        For a value outside [0, 1].

    Examples
    --------
    >>> from crevice.analysis_viewers import hydration_rgb
    >>> hydration_rgb(0.25), hydration_rgb(None)
    ([0.75, 0.0, 0.25], [0.55, 0.55, 0.55])
    """
    if occupancy is None:
        return [.55,.55,.55]
    if not math.isfinite(float(occupancy)) or not 0<=occupancy<=1:
        raise ValueError('Hydration colour requires occupancy in [0,1] or None')
    return [1-float(occupancy),0.,float(occupancy)]


def _camera(frame, points, basis):
    import numpy as np
    from .models import VoidCast, VoidComponent, ChannelPoint
    from .presentation import scene_geometry
    points = np.asarray(points)
    sample = points[::max(1,int(len(points)//1000))]
    cloud = tuple(ChannelPoint(i,tuple(map(float,p)),0.,0.,float(p@basis[2])) for i,p in enumerate(sample))
    component = VoidComponent(1,'reference',tuple(map(float,points.mean(0))),float(len(points)),cloud)
    cast = VoidCast((component,),cloud,1.,0.,'rolling',{'grid_basis':basis.tolist()})
    geometry = scene_geometry(cast,frame)
    geometry['vertical_span'] += 10.0
    geometry['radius'] = max(geometry['radius'],geometry['vertical_span']*.6)
    return geometry


def write_analysis_assets(analysis, report, output_dir, prefix):
    """Write the shared scene definition and assets of the coordinated hydration views.

    One scene definition drives the PyMOL, VMD and ChimeraX scripts written by
    :func:`write_analysis_views`, so every viewer shows the same colours,
    residues, surfaces and camera. Residues are coloured by their all-frame
    water-contact occupancy on an absolute scale (:func:`hydration_rgb`);
    direct boundary contributors are shown as sticks; the reference cavity and
    the mean water density are display surfaces; up to 12 typed contacts
    observed in the displayed reference and at 50% occupancy or more are drawn
    as dashes; trajectories supply first, middle and last aligned snapshots.
    With water membership, each snapshot's member-water oxygens are written as
    their own group. This function writes no text into any scene.

    Parameters
    ----------
    analysis : dict
        Hydration analysis with ``reference_frame``, reference atom indices,
        snapshots, optional density, interactions and ``member_waters``.
    report : dict
        The hydration report (per-residue statistics and water status).
    output_dir : str or Path
        Output directory.
    prefix : str
        File-name prefix.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_analysis_scene.json : JSON
        Residues (roles, occupancy, colours, atom serials), surfaces, contacts,
        camera, snapshots, ``member_waters`` (``null`` without membership),
        display definitions and limitations. :func:`write_analysis_views` adds
        an ``annotations`` entry recording the scene text.
    PREFIX_analysis_protein.pdb, PREFIX_analysis_protein.identities.json : structure
        Renumbered reference protein and its identity map.
    PREFIX_analysis_frame_NNNNNN.pdb (+ ``.identities.json``) : snapshots
        Aligned protein of each displayed source frame.
    PREFIX_analysis_member_waters_NNNNNN.pdb : member waters
        Member water oxygens of each displayed snapshot (membership runs only;
        none for a snapshot without members).
    PREFIX_analysis_geometry.npz : arrays
        Reference coordinates and pre-computed surface meshes.
    PREFIX_analysis_cavity.dx, PREFIX_analysis_cavity_display.dx, PREFIX_analysis_cavity_local.dx : grids
        Reference cavity (measured and display) for the viewers.

    Colours and representation
    --------------------------
    See :func:`crevice.analysis_scene_templates.write_scene_scripts`: residues
    red (0%) to blue (100%) water-contact occupancy, grey unobserved; pale
    reference cavity (18% opaque); teal water density (42% opaque); gold, teal,
    purple and blue typed-contact dashes; orange 0.7 Å member-water spheres.
    """
    import numpy as np
    from .presentation import write_viewer_structure
    from .hydration import minimum_vectors
    from .volume_export import presentation_grid
    from .water_density import write_scalar_dx
    root = Path(output_dir).resolve()
    frame = analysis['reference_frame']
    files = {}
    def path(key,suffix):
        p=root/(prefix+suffix); files[key]=str(p); return p
    structure = write_viewer_structure(frame,path('analysis_protein_pdb','_analysis_protein.pdb'))
    files['analysis_protein_identities_json'] = str(structure.with_suffix('.identities.json'))
    records = {r['residue']:r for r in report['residues']}
    identities = json.loads(structure.with_suffix('.identities.json').read_text())
    residue_mapping = {r['original']:r for r in identities['residue_mapping']}
    protein_xyz = np.asarray([a.coord for a in frame.atoms])
    atom_indices = np.asarray(analysis['reference_atom_indices'])
    by_residue = {}
    for key,row in records.items():
        serials = [i+1 for i,a in enumerate(frame.atoms) if a.residue_key.label==key]
        occupancy = row['water_contact_occupancy']['mean'] if report['water_population'] else None
        by_residue[key] = {'residue':key,'role':row['role'],'hydration_occupancy':occupancy,
            'color_rgb':hydration_rgb(occupancy),'atom_serials':serials,
            'mean_water_count':row['statistics']['water_count']['mean'],
            'mean_SASA_A2':row['statistics']['sasa_A2']['mean'],
            'boundary_area_A2':row.get('mean_boundary_area_A2'),
            'centroid_A':protein_xyz[np.asarray(serials)-1].mean(0).tolist(),
            **{k:residue_mapping[key][k] for k in ['viewer_chain','viewer_resid']}}
    direct = [key for key,r in by_residue.items() if r['role'] in {'boundary_lining','ever_boundary'}]
    displayed = direct or list(by_residue)
    partners = [key for key in by_residue if key not in displayed]
    field_meshes = {}
    surfaces = []
    basis = np.eye(3)
    site_points = np.asarray([by_residue[k]['centroid_A'] for k in displayed])
    try:
        from skimage.measure import marching_cubes
    except ImportError:
        marching_cubes = None
    def mesh(field,origin,spacing,field_basis,key,level):
        if marching_cubes is None or not float(field.min()) < level < float(field.max()):
            return False
        vertices,faces,normals,_=marching_cubes(field,level=level,spacing=(spacing,)*3)
        field_meshes.update({key+'_vertices':vertices@field_basis+origin,
                            key+'_faces':faces,
                            key+'_normals':normals@field_basis})
        return True
    cavity_path=analysis.get('reference_volume_path')
    if cavity_path:
        from .residue_evidence import read_binary_dx
        grid,origin,deltas = read_binary_dx(cavity_path)
        spacing=float(np.linalg.norm(deltas[0])); basis=deltas/spacing
        if not np.allclose(basis@basis.T,np.eye(3),atol=1e-7):
            raise ValueError('Unified cavity scenes require orthogonal, equally spaced map samples')
        site_points=np.argwhere(grid)@deltas+origin
        measured=path('analysis_cavity_dx','_analysis_cavity.dx')
        if Path(cavity_path).resolve()!=measured.resolve():shutil.copy2(cavity_path,measured)
        field=presentation_grid(grid,spacing,.6)
        local=path('analysis_cavity_local_dx','_analysis_cavity_local.dx')
        write_scalar_dx(field,origin@basis.T,spacing,local)
        display=path('analysis_cavity_display_dx','_analysis_cavity_display.dx')
        write_scalar_dx(field,origin,spacing,display,basis=basis)
        available=mesh(field,origin,spacing,basis,'cavity',.5)
        surfaces.append({'key':'cavity','name':'crevice_cavity','dx':display.name,'local_dx':local.name,
                         'transform':basis.T.tolist(),'level':.5,'color_rgb':[.69,.73,.73],
                         'opacity':.18,'mesh_available':available,'meaning':'reference cavity display; measured occupancy retained separately'})
    density=analysis.get('water_density')
    if density is not None and density['metadata']['status']=='observed':
        level=analysis['analysis_options']['density_level'];field=density['display_density_A3']
        available=mesh(field,density['origin'],density['spacing_A'],np.eye(3),'water',level)
        if float(field.min()) < level < float(field.max()):
            surfaces.append({'key':'water','name':'crevice_water_density','dx':prefix+'_water_density_display.dx',
                             'local_dx':prefix+'_water_density_display.dx','transform':np.eye(3).tolist(),
                             'level':level,'color_rgb':[.08,.66,.59],'opacity':.42,
                             'mesh_available':available,'meaning':'Gaussian display of mean water-oxygen density; units oxygen / A^3'})
    geometry=_camera(frame,site_points,basis)
    contacts=[]
    interaction_path=root/(prefix+'_interactions.json')
    if analysis.get('interactions') is not None and interaction_path.is_file():
        interactions=json.loads(interaction_path.read_text())
        colors={'hydrogen_bond':[.83,.61,.12],'water_bridge_hbond':[.05,.66,.59],
                'salt_bridge_candidate':[.65,.30,.65],
                'aromatic_parallel_candidate':[.30,.45,.70],'aromatic_edge_face_candidate':[.30,.45,.70]}
        eligible=[e for e in interactions['edges'] if e['reference_observation'] and e['observed_occupancy']>=.5
                  and (e['source'] in displayed or e['target'] in displayed)]
        eligible.sort(key=lambda e:(-e['observed_occupancy'],e['kind'],e['source'],e['target']))
        full=analysis['reference_full_coordinates'];rotation=analysis['reference_rotation'];box=analysis['reference_box']
        for edge in eligible[:12]:
            observation=edge['reference_observation'];indices=observation['atom_indices']
            if edge['kind']=='hydrogen_bond':indices=[indices[0],indices[-1]]
            coordinates=[np.asarray(full[indices[0]])]
            for index in indices[1:]:
                vector=(full[index]-coordinates[-1])@rotation.T
                coordinates.append(coordinates[-1]+minimum_vectors(vector[None,:],box)[0]@rotation)
            contacts.append({'kind':edge['kind'],'source':edge['source'],'target':edge['target'],
                'occupancy':edge['observed_occupancy'],'points_A':[p.tolist() for p in coordinates],
                'color_rgb':colors[edge['kind']]})
    snapshots=[]
    for snapshot,time in zip(analysis.get('snapshot_frames',[frame]),analysis.get('snapshot_times_ps',[None])):
        name=prefix+f'_analysis_frame_{snapshot.frame_index:06d}.pdb';snapshot_path=root/name
        write_viewer_structure(snapshot,snapshot_path)
        snapshots.append({'pdb':name,'source_frame_index':snapshot.frame_index,'time_ps':time})
        files[f'analysis_frame_{snapshot.frame_index}_pdb']=str(snapshot_path)
        files[f'analysis_frame_{snapshot.frame_index}_identities_json']=str(snapshot_path.with_suffix('.identities.json'))
    member_waters=None
    if analysis.get('member_waters') is not None:
        # Own named group, visible by default; hidden independently of view modes.
        member_waters={'object':'crevice_member_waters','color_rgb':[1.,.5,0.],'sphere_radius_A':.7,'visible_by_default':True,
            'definition':"aligned water oxygens whose centre lies in an occupied voxel of that snapshot's measured region; geometric membership, not binding or hydrogen bonding",
            'snapshots':[]}
        for snapshot in snapshots:
            record=analysis['member_waters'].get(int(snapshot['source_frame_index']))
            entry={'source_frame_index':snapshot['source_frame_index'],'count':None if record is None else len(record['water_ids']),
                   'water_ids':None if record is None else record['water_ids'],'pdb':None,
                   'status':'unresolved_geometry' if record is None else 'measured'}
            if record and record['water_ids']:
                name=prefix+f"_analysis_member_waters_{snapshot['source_frame_index']:06d}.pdb"
                lines=[f"HETATM{i:5d}  O   HOH W{i:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00           O"
                       for i,p in enumerate(record['positions_A'],1)]
                (root/name).write_text('REMARK   CREVICE member water oxygens; residue number is the list position, IDs are in the scene JSON\n'+'\n'.join(lines)+'\nEND\n')
                entry['pdb']=name;files[f"analysis_member_waters_{snapshot['source_frame_index']}_pdb"]=str(root/name)
            member_waters['snapshots'].append(entry)
    np.savez_compressed(path('analysis_geometry_npz','_analysis_geometry.npz'),protein_coords_A=protein_xyz,
                        source_atom_indices=atom_indices,**field_meshes)
    scene={'prefix':prefix,'protein_pdb':structure.name,'geometry_npz':prefix+'_analysis_geometry.npz',
           'region_definition':report['settings'].get('region_definition'),
           'frame_count':report['frame_count'],'source_reference_frame_index':frame.frame_index,
           'hydration_definition':'fraction of analysed frames with at least one water oxygen within the configured residue contact cutoff',
           'static_definition':'single observed snapshot contacts; 0/100% colours are not equilibrium hydration probabilities' if report['mode'].startswith('static') else None,
           'water_status':report['water_status'],'hydration_cutoff_A':report['settings']['contact_cutoff_A'],
           'color_scale':{'minimum':0,'maximum':1,'red_at_zero':[1,0,0],'blue_at_one':[0,0,1],
                          'missing_rgb':hydration_rgb(None),'interpolation':'linear RGB; absolute contact occupancy, not maximum-normalized water counts'},
           'residues':list(by_residue.values()),'display_residues':displayed,'partner_residues':partners,
           'display_definition':'all identified direct boundary contributors as sticks; all selected residues if no boundary roles; partner sticks can be enabled separately',
           'surfaces':surfaces,'contacts':contacts,
           'contact_display_definition':'up to 12 edges observed in the displayed reference and at >=50% saved-frame occupancy, sorted by occupancy; present only in network view',
           'camera':geometry,'snapshots':snapshots,'member_waters':member_waters,
           'density_domain':density['metadata']['domain'] if density is not None else None,
           'pymol_surface_status':'available' if marching_cubes is not None else 'unavailable_optional_scikit_image',
           'limitations':['Identified geometric contributors are not proven crucial functional residues.',
                          'Density is a mean spatial observation; a displayed snapshot can differ from that mean.',
                          'Static missing waters are unobserved; chemical classification has separate coverage limits.']}
    path('analysis_scene_json','_analysis_scene.json').write_text(json.dumps(scene,indent=2,allow_nan=False)+'\n')
    return files


@annotate_option
def write_analysis_views(output_dir,prefix,files,*,dpi=240,annotate=None):
    """Write the native viewer scripts, evidence table and overview figures for a scene definition.

    Calls :func:`crevice.analysis_scene_templates.write_scene_scripts` and
    :func:`crevice.analysis_report.write_analysis_report` on the assets written
    by :func:`write_analysis_assets`.

    Parameters
    ----------
    output_dir : str or Path
        Directory containing ``PREFIX_analysis_scene.json`` and its assets.
    prefix : str
        File-name prefix.
    files : dict of str to str
        Files written so far; used to find optional interaction and density records.
    dpi : int, default 240
        Resolution of the overview figures.
    annotate : bool, optional
        Draw the scene legend, captions and banners and the figure titles.
        ``None`` (default) inherits the surrounding setting, which is off unless
        enabled with :func:`crevice.presentation.figure_annotations` or ``--annotate``.

    Returns
    -------
    dict of str to str
        Output keys and paths of both writers.

    Outputs
    -------
    PREFIX_analysis*.{pml,vmd,tcl,cxc}, PREFIX_analysis*_chimerax.py : scenes
        See :func:`crevice.analysis_scene_templates.write_scene_scripts`.
    PREFIX_analysis_evidence.csv, PREFIX_analysis_overview.png, PREFIX_hydration_legend.png, ... : table and figures
        See :func:`crevice.analysis_report.write_analysis_report`. No HTML is written.

    Colours and representation
    --------------------------
    As documented by the two writers.
    """
    from .analysis_scene_templates import write_scene_scripts
    from .analysis_report import write_analysis_report
    extra=write_scene_scripts(output_dir,prefix)
    extra.update(write_analysis_report(output_dir,prefix,{**files,**extra},dpi=dpi))
    return extra

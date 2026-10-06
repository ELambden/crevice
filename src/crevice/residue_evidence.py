"""Explain a measured void boundary and its nonlocal residue contacts.

This is the engine of ``crevice residue-evidence`` (and of the boundary
evidence inside ``cavity-trajectory``). The 0.5 surface of a measured binary
volume map (a ``*_volume.dx`` from ``cast`` or ``publish``) is triangulated;
each triangle's area is assigned to the residue(s) whose van der Waals sphere
is nearest (:func:`boundary_residue_evidence`). Residues with assigned area
*line* the boundary; residues that do not line it but make a nonlocal contact
with a lining residue are *partners*. :func:`write_residue_evidence_bundle`
writes the tables, figures and a stick overlay on the cast scenes;
:func:`select_stick_residues` decides which residues are drawn as sticks
(``--stick-residues``). :func:`profile_constriction_evidence` lists the
atoms that limit a channel profile's continuous path.

Boundary area is geometric evidence. Contact partners and correlated motion
support follow-up hypotheses; neither establishes functional control.
"""
from __future__ import annotations
import csv,json,math
from pathlib import Path

from .presentation import FigureText, annotate_option
from .radii import RadiusSet, radii_option


def read_binary_dx(path):
    """Read a regular binary occupancy DX; refuse smoothed display/density maps.

    Parameters
    ----------
    path : str or Path
        OpenDX scalar file as written by
        :func:`crevice.volume_export.write_void_cast_dx` without smoothing
        (``*_volume.dx``).

    Returns
    -------
    grid : numpy.ndarray
        Values (0 or 1) with the file's ``counts`` shape (z fastest).
    origin : numpy.ndarray, shape (3,)
        Position of grid point (0, 0, 0), Å.
    deltas : numpy.ndarray, shape (3, 3)
        Grid step vectors as rows, Å.

    Raises
    ------
    ValueError
        For inconsistent counts, a singular or non-finite basis, or any value
        other than 0 and 1 (display fields and densities are refused).
    """
    import numpy as np
    lines=Path(path).read_text().splitlines()
    counts=next(line for line in lines if 'class gridpositions counts' in line)
    shape=tuple(int(v) for v in counts.split()[-3:])
    origin=np.asarray([float(v) for v in next(line for line in lines if line.startswith('origin ')).split()[1:]])
    deltas=np.asarray([[float(v) for v in line.split()[1:]] for line in lines if line.startswith('delta ')])
    start=next(i for i,line in enumerate(lines) if 'data follows' in line)
    chunks=[]
    for line in lines[start+1:]:
        if line.startswith(('attribute','object','component')):break
        chunks.append(line)
    values=np.fromstring(' '.join(chunks),sep=' ')
    if len(shape)!=3 or any(v<2 for v in shape) or values.size!=math.prod(shape):
        raise ValueError('DX counts do not match a nonempty 3D scalar grid')
    if deltas.shape!=(3,3) or not np.isfinite(deltas).all() or abs(np.linalg.det(deltas))<1e-12:
        raise ValueError('DX deltas must define a finite nonsingular basis')
    if not np.isfinite(origin).all() or not np.isin(values,[0.,1.]).all():
        raise ValueError('Use the measured binary volume DX, not a smoothed display map or density field')
    return values.reshape(shape),origin,deltas


@radii_option
def boundary_residue_evidence(frame,grid,origin,deltas,*,lining_distance=1.5,
                              contact_cutoff=4.5,residue_groups=None,boundary_exclusion=None,obstacle_atoms=None,
                              radii: RadiusSet | str | None = None):
    """Area-weighted nearest-sphere attribution of a binary map's 0.5 surface.

    Each triangle centroid contributes its area. Ties across residues share that
    area; faces farther than lining_distance from protein are left unassigned.
    Artificial mouth/crop faces can remain: assignment is not molecular SASA.

    Algorithm: marching cubes at 0.5 on the binary grid; for each triangle, the
    surface gap to every nearby atom sphere (centre distance minus radius); the
    smallest gap within 1e-7 Å decides the owner residue(s), provided it is at
    most ``lining_distance`` and the face is not excluded. A residue contact
    network (atom-centre cutoff ``contact_cutoff``) then gives, for every
    residue, its nonlocal contacts with lining residues (nonlocal = different
    chain or more than two positions apart in residue order).

    Parameters
    ----------
    frame : StructureFrame
        Structure in the coordinates of the map.
    grid : array_like
        Binary 3D occupancy grid with an empty one-voxel border.
    origin : array_like, shape (3,)
        Grid origin, Å.
    deltas : array_like, shape (3, 3)
        Grid step vectors, Å.
    lining_distance : float, default 1.5
        Largest surface gap (Å) at which a boundary face is assigned to a
        residue.
    contact_cutoff : float, default 4.5
        Residue contact cutoff (Å, minimum heavy-atom centre distance).
    residue_groups : dict of str to str, optional
        Residue ID to group label (for example helix names), copied into each
        record's ``group``.
    boundary_exclusion : callable, optional
        Called with the triangle centroids (n x 3); returns one boolean per
        triangle, True for faces to leave unassigned (for example crop faces of
        a reference region).
    obstacle_atoms : sequence of Atom, optional
        Atoms that own boundary faces; default the frame's non-HETATM heavy
        atoms. HETATM obstacles (lipids, ions, ligands) are reported separately.
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
    dict
        ``method``, ``units``, settings, ``boundary_area_A2``,
        ``excluded_boundary_area_A2``, ``assigned_boundary_area_A2`` and
        fraction, protein and other-obstacle assigned areas, ``triangle_count``,
        ``boundary_overlap_area_fraction``, ``grid_voxel_volume_A3``,
        ``measurement_volume_A3``, ``residues`` (one record per residue:
        ``group``, ``role`` = ``boundary_lining``, ``nonlocal_lining_partner`` or
        ``other``, ``boundary_area_A2``, ``boundary_fraction``,
        ``minimum_boundary_gap_A``, ``nonlocal_lining_contact_count``,
        ``contacted_lining_boundary_fraction``, ``nonlocal_lining_partners``;
        sorted by area, then contacted fraction), ``network``,
        ``lining_residues``, ``ranking`` and ``limitations``.

    Raises
    ------
    ValueError
        For a negative ``lining_distance``, a grid that is not binary, empty or
        full, occupancy touching the grid border, or an exclusion of the wrong
        length.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    from skimage.measure import marching_cubes
    from .networks import build_residue_network
    from .radii import atom_vdw_radius
    if not math.isfinite(lining_distance) or lining_distance<0:
        raise ValueError('lining_distance must be finite and nonnegative')
    grid=np.asarray(grid);origin=np.asarray(origin);deltas=np.asarray(deltas)
    if grid.ndim!=3 or not np.isin(grid,[0,1]).all() or not (grid==1).any() or not (grid==0).any():
        raise ValueError('Boundary evidence requires a nonempty binary occupancy grid with exterior zeros')
    if any(np.any(np.take(grid,index,axis=axis)) for axis in range(3) for index in (0,-1)):
        raise ValueError('Binary occupancy needs an empty exterior border; boundary truncation would be ambiguous')
    vertices,faces,_,_=marching_cubes(grid.astype('float32'),level=.5)
    triangles=(vertices@deltas+origin)[faces]
    areas=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)/2
    centers=triangles.mean(axis=1)
    excluded=np.zeros(len(centers),dtype=bool) if boundary_exclusion is None else np.asarray(boundary_exclusion(centers),dtype=bool)
    if excluded.shape!=(len(centers),):raise ValueError('Boundary exclusion must return one boolean per triangle')
    atoms=frame.selected_atoms(include_hydrogen=False,include_hetero=False) if obstacle_atoms is None else tuple(obstacle_atoms)
    xyz=np.asarray([a.coord for a in atoms]);radii=np.asarray([atom_vdw_radius(a) for a in atoms])
    tree=cKDTree(xyz);nearest,index=tree.query(centers)
    candidates=tree.query_ball_point(centers,nearest-radii[index]+radii.max()+1e-7)
    rows=np.repeat(np.arange(len(centers)),[len(c) for c in candidates]);cols=np.concatenate(candidates).astype(int)
    gaps=np.linalg.norm(centers[rows]-xyz[cols],axis=1)-radii[cols]
    minimum=np.full(len(centers),np.inf);np.minimum.at(minimum,rows,gaps)
    owners={}
    for row,col,gap in zip(rows,cols,gaps):
        if not excluded[row] and gap<=minimum[row]+1e-7 and gap<=lining_distance:
            owners[(int(row),atoms[int(col)].residue_key.label)]=float(gap)
    owner_count=np.bincount([row for row,_ in owners],minlength=len(centers))
    area_by_residue={};gap_by_residue={}
    for (row,residue),gap in owners.items():
        area_by_residue[residue]=area_by_residue.get(residue,0.)+float(areas[row]/owner_count[row])
        gap_by_residue[residue]=min(gap_by_residue.get(residue,float('inf')),gap)
    obstacle_labels={a.residue_key.label for a in atoms if a.hetero}
    obstacle_areas={k:v for k,v in area_by_residue.items() if k in obstacle_labels}
    protein_assigned=sum(v for k,v in area_by_residue.items() if k not in obstacle_labels)
    total=float(areas.sum());assigned=float(areas[owner_count>0].sum())
    network=build_residue_network(frame,cutoff=contact_cutoff,distance_metric='center',residue_groups=residue_groups)
    partners={node.id:[] for node in network.nodes}
    for edge in network.edges:
        if not edge.metadata.get('nonlocal_contact'):continue
        for a,b in [(edge.source,edge.target),(edge.target,edge.source)]:
            if b in area_by_residue:
                partners[a].append({'residue':b,'center_distance_A':edge.distance,
                                    'boundary_fraction':area_by_residue[b]/total,
                                    'closest_pair':edge.metadata['closest_pair'],
                                    'interaction_candidate':edge.interaction})
    records=[]
    for node in network.nodes:
        area=area_by_residue.get(node.id,0.);support=partners[node.id]
        records.append({'residue':node.id,'group':node.metadata.get('group','unassigned'),
                        'role':'boundary_lining' if area else 'nonlocal_lining_partner' if support else 'other',
                        'boundary_area_A2':area,'boundary_fraction':area/total,
                        'minimum_boundary_gap_A':gap_by_residue.get(node.id),
                        'nonlocal_lining_contact_count':len(support),
                        'contacted_lining_boundary_fraction':sum(p['boundary_fraction'] for p in support),
                        'nonlocal_lining_partners':sorted(support,key=lambda p:(-p['boundary_fraction'],p['residue']))})
    records.sort(key=lambda r:(-r['boundary_area_A2'],-r['contacted_lining_boundary_fraction'],r['residue']))
    return {'method':'binary_boundary_triangle_centroid_nearest_vdw_sphere',
            'units':{'length':'angstrom','area':'angstrom squared'},'lining_distance_A':lining_distance,
            'contact_cutoff_A':contact_cutoff,'contact_metric':'minimum heavy-atom center distance',
            'boundary_area_A2':total,'excluded_boundary_area_A2':float(areas[excluded].sum()),'assigned_boundary_area_A2':assigned,'assigned_boundary_fraction':assigned/total,
            'protein_assigned_boundary_area_A2':float(protein_assigned),
            'other_obstacle_boundary_area_A2':float(sum(obstacle_areas.values())),
            'other_obstacle_boundary_areas_A2':obstacle_areas,
            'triangle_count':len(areas),'boundary_overlap_area_fraction':float(areas[minimum<0].sum()/total),
            'grid_voxel_volume_A3':float(abs(np.linalg.det(deltas))),
            'measurement_volume_A3':float(grid.sum()*abs(np.linalg.det(deltas))),
            'residues':records,'network':network.to_dict(),
            'lining_residues':[r['residue'] for r in records if r['boundary_area_A2']>0],
            'ranking':'Boundary contributors by assigned area; supporting partners by contacted lining fraction. No combined causal importance score.',
            'limitations':['Grid-boundary quadrature, not molecular solvent-accessible surface area.',
                          'Artificial mouth/crop faces are not distinguished anatomically.',
                          'Nonlocal means different chain or >2 positions in observed residue order.',
                          'Static geometry cannot establish functional control, energetic interactions or mutational effects.']}


def restore_viewer_identities(frame,structure_path):
    """Restore analysis IDs when reading a CREVICE viewer PDB and its sidecar.

    Viewer PDBs (:func:`crevice.presentation.write_viewer_structure`) renumber
    chains, residues and atoms; their ``.identities.json`` sidecar records the
    originals. When the sidecar exists, every atom gets back its original chain,
    residue number, insertion code and serial.

    Parameters
    ----------
    frame : StructureFrame
        Structure read from ``structure_path``.
    structure_path : str or Path
        Path of the viewer PDB; the sidecar is looked for beside it.

    Returns
    -------
    frame : StructureFrame
        With original identities, or unchanged when there is no sidecar.
    record : dict
        ``definition`` and ``sidecar`` (resolved path or ``None``).

    Raises
    ------
    ValueError
        If the sidecar does not match the structure.
    """
    from dataclasses import replace
    import re
    path=Path(structure_path).with_suffix('.identities.json')
    if not path.is_file():return frame,{'definition':'input structure identities','sidecar':None}
    data=json.loads(path.read_text());serials=data.get('original_serials_in_viewer_order',[])
    if len(serials)!=len(frame.atoms) or data.get('atom_count')!=len(frame.atoms):
        raise ValueError('Viewer identity sidecar atom count does not match the structure')
    residues={(r['viewer_chain'],r['viewer_resid']):r['original'] for r in data['residue_mapping']}
    atoms=[]
    for atom in frame.atoms:
        label=residues.get((atom.chain_id,atom.resid))
        if label is None or not 1<=atom.serial<=len(serials):
            raise ValueError('Viewer atom/residue identity is absent from its sidecar')
        chain,tail=label.rsplit(':',1)
        if not tail.startswith(atom.resname):raise ValueError('Viewer residue name differs from its identity sidecar')
        match=re.fullmatch(r'(-?\d+)([^\d]*)',tail[len(atom.resname):])
        if match is None:raise ValueError('Cannot decode the original residue identifier')
        atoms.append(replace(atom,chain_id='' if chain=='_' else chain,resid=int(match[1]),icode=match[2],serial=serials[atom.serial-1]))
    return replace(frame,atoms=tuple(atoms)),{'definition':'original identities restored from CREVICE viewer sidecar','sidecar':str(path.resolve())}


def _viewer_serial_mapping(frame,scene):
    """Validate direct PDB IDs or the exporter's reversible source-ID sidecar."""
    import ast
    import numpy as np
    source=scene.read_text()
    body=source.split('python\n',1)[1].rsplit('python end',1)[0]
    tree=ast.parse(body)
    loader=next((node for node in ast.walk(tree) if isinstance(node,ast.Call)
                 and isinstance(node.func,ast.Attribute) and node.func.attr=='load'
                 and isinstance(node.func.value,ast.Name) and node.func.value.id=='cmd'
                 and node.args and isinstance(node.args[0],ast.Call)
                 and isinstance(node.args[0].func,ast.Name) and node.args[0].func.id=='crevice_input'),None)
    if loader is None:raise ValueError('Cannot identify the CREVICE scene protein input')
    local,fallback=(ast.literal_eval(v) for v in loader.args[0].args[:2])
    protein=scene.parent/local
    if not protein.is_file():protein=Path(fallback)
    coordinates={}
    for line in protein.read_text().splitlines():
        if line[:6].strip() in {'ATOM','HETATM'}:
            serial=int(line[6:11])
            if serial in coordinates:raise ValueError('Viewer PDB atom serials must be unique')
            coordinates[serial]=[float(line[i:i+8]) for i in [30,38,46]]
    ids=[a.serial for a in frame.atoms]
    if len(set(ids))!=len(ids):raise ValueError('Analysis atom serials must be unique for a residue overlay')
    candidates=[('direct PDB serial identities',ids)]
    sidecar=protein.with_suffix('.identities.json')
    if sidecar.is_file():
        originals=json.loads(sidecar.read_text()).get('original_serials_in_viewer_order',[])
        remap={serial:i for i,serial in enumerate(originals,1)}
        if len(remap)==len(originals) and set(ids)==set(remap):
            candidates.append(('reversible viewer identity sidecar',[remap[i] for i in ids]))
    expected=np.asarray([a.coord for a in frame.atoms]);valid=[]
    for definition,mapped in candidates:
        if len(mapped)==len(coordinates) and set(mapped)==set(coordinates):
            actual=np.asarray([coordinates[i] for i in mapped])
            if np.allclose(actual,expected,atol=.001,rtol=0):valid.append((definition,mapped))
    if not valid:raise ValueError('Scene protein coordinates and atom identities do not match the analysis reference')
    if any(v[1]!=valid[0][1] for v in valid):raise ValueError('Ambiguous scene atom identity mapping')
    return valid[0][1],protein,valid[0][0]


STICK_RESIDUE_DEFAULT='top:12'


def select_stick_residues(report,stick_residues=STICK_RESIDUE_DEFAULT):
    """Choose the boundary/partner rows shown as viewer sticks.

    ``'top:N'`` keeps the N highest-ranked boundary contributors and nonlocal
    partners (default N=12), ``'all'`` keeps every boundary/partner residue and
    an explicit selection is residue IDs (a comma-separated string, a sequence
    or ``@file`` with one or more IDs per line; ``#`` starts a comment line).
    Selected residues keep their reported role and colour; an ID outside those
    two roles is an error.

    Parameters
    ----------
    report : dict
        :func:`boundary_residue_evidence` result.
    stick_residues : str or sequence of str, default "top:12"
        The request (CLI ``--stick-residues``); ``None`` means the default.

    Returns
    -------
    lining_rows : list of dict
        Selected boundary residues, in report order (largest area first).
    partner_rows : list of dict
        Selected nonlocal partners, by contacted lining fraction.
    selection_record : dict
        ``request``, ``mode`` (``top``, ``all`` or ``selected``),
        ``definition``, ``lining_count`` and ``partner_count``.

    Raises
    ------
    ValueError
        For a malformed ``top:N``, an empty selection, or IDs absent from the
        report or neither lining nor partner.

    Examples
    --------
    >>> from crevice.residue_evidence import select_stick_residues
    >>> report = {"residues": [
    ...     {"residue": "A:LEU5", "role": "boundary_lining", "contacted_lining_boundary_fraction": 0.0},
    ...     {"residue": "A:ASP9", "role": "nonlocal_lining_partner", "contacted_lining_boundary_fraction": 0.2},
    ...     {"residue": "A:GLY1", "role": "other", "contacted_lining_boundary_fraction": 0.0}]}
    >>> lining, partners, record = select_stick_residues(report, "all")
    >>> [r["residue"] for r in lining], [r["residue"] for r in partners], record["mode"]
    (['A:LEU5'], ['A:ASP9'], 'all')
    """
    rows=report['residues']
    lining=[r for r in rows if r['role']=='boundary_lining']
    support=sorted((r for r in rows if r['role']=='nonlocal_lining_partner'),
                   key=lambda r:(-r['contacted_lining_boundary_fraction'],r['residue']))
    request=STICK_RESIDUE_DEFAULT if stick_residues is None else stick_residues
    if isinstance(request,str) and request.strip().lower()=='all':
        mode='all';definition='All boundary residues and all nonlocal nonlining partners'
    elif isinstance(request,str) and request.strip().lower().startswith('top:'):
        text=request.strip()[4:]
        if not text.isdigit():raise ValueError('stick_residues top:N needs a nonnegative integer N')
        count=int(text);lining,support=lining[:count],support[:count]
        mode='top';definition=f'Top {count} direct boundary contributors and top {count} nonlocal nonlining partners'
    else:
        if isinstance(request,str):
            text=request.strip()
            if text.startswith('@'):text=Path(text[1:]).read_text()
            requested=[v.strip() for line in text.splitlines() for v in line.split(',') if v.strip() and not v.strip().startswith('#')]
        else:
            requested=[str(v).strip() for v in request]
        if not requested:raise ValueError('stick_residues selection is empty')
        roles={r['residue']:r['role'] for r in rows}
        unknown=[v for v in requested if v not in roles]
        if unknown:raise ValueError('Selected stick residues are absent from the evidence report: '+', '.join(unknown))
        other=[v for v in requested if roles[v] not in {'boundary_lining','nonlocal_lining_partner'}]
        if other:raise ValueError('Selected stick residues are neither boundary residues nor nonlocal partners: '+', '.join(other))
        chosen=set(requested);lining=[r for r in lining if r['residue'] in chosen];support=[r for r in support if r['residue'] in chosen]
        mode='selected';definition=f'{len(chosen)} user-selected boundary/partner residues'
    return lining,support,{'request':request if isinstance(request,str) else list(request),'mode':mode,'definition':definition,
                           'lining_count':len(lining),'partner_count':len(support)}


@annotate_option
def write_residue_evidence_bundle(report,frame,output_dir,*,prefix='crevice',scene_path=None,volume_path=None,dpi=300,
                                  stick_residues=STICK_RESIDUE_DEFAULT,annotate=None):
    """Write the tables, figure and viewer overlays of ``crevice residue-evidence``.

    ``report`` (from :func:`boundary_residue_evidence`) assigns every boundary
    face of a measured binary cast to the nearest atom within
    ``lining_distance``: residues owning boundary area are direct boundary
    contributors, ranked by their fraction of assigned boundary area.
    Non-lining residues in contact (``contact_cutoff``) with lining residues
    are nonlocal structural partners, ranked by the boundary fraction of the
    lining residues they contact. The figure ranks the top 12 residues per
    role. With a scene (``scene_path``, or one built from ``volume_path``),
    PyMOL, VMD and ChimeraX overlays show the selected residues as opaque
    sticks on the cartoon with the cast. These are geometric candidates for
    follow-up, not functional assignments.

    Parameters
    ----------
    report : dict
        Evidence report from :func:`boundary_residue_evidence`.
    frame : StructureFrame
        The structure the report was computed on.
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    scene_path : str or Path, optional
        Existing CREVICE volume scene (``*_volume.pml`` etc.) to overlay.
    volume_path : str or Path, optional
        Binary DX used to build a scene when ``scene_path`` is not given.
    dpi : int, default 300
        Figure resolution.
    stick_residues : str, default "top:12"
        Overlay membership: ``top:N`` per role, ``all``, comma-separated IDs or
        ``@FILE`` (see :func:`select_stick_residues`). The bar chart and tables
        are unchanged by this option.
    annotate : bool, optional
        Draw the panel titles ("Direct boundary contributors", "Nonlocal
        structural partners") and the figure title. The overlays contain no
        text. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_residue_evidence.json : JSON
        Full report, display membership, colours and overlay definition; kept
        because ``trajectory --lining-evidence`` and the viewer overlays read it.
    PREFIX_residue_evidence.csv : table
        Per residue: group, role, boundary area (Å²) and fraction, minimum gap
        (Å), nonlocal lining contacts, contacted-lining boundary fraction.
    PREFIX_residue_evidence_partners.csv : table
        One row per nonlocal contact between a residue and a boundary-lining
        residue: ``residue``, ``lining_residue``, ``interaction_candidate``,
        ``center_distance_A`` (edge distance under the network metric),
        ``surface_distance_A``, ``closest_atoms`` and
        ``lining_boundary_fraction`` (the lining residue's assigned boundary
        fraction).
    PREFIX_residue_evidence.png : figure
        Two bar charts: "Assigned boundary area (%)" of direct contributors and
        "Boundary fraction of contacted lining residues (%)" of partners
        (residue IDs as tick labels).
    PREFIX_residue_context.{pml,vmd,tcl,cxc} : scenes
        Overlays on the volume scene, with ``PREFIX_residue_reference.npz``
        and copied scene assets (only with a scene).

    Colours and representation
    --------------------------
    * Figure: orange-gold bars (``#d98b32``) for boundary contributors, magenta
      (``#a34477``) for nonlocal partners.
    * Overlays: gold boundary residues (``[0.85, 0.55, 0.20]``) and magenta
      nonlocal partners (``[0.64, 0.27, 0.47]``) as opaque sticks (radius
      0.19 Å); the rest of the protein as cartoon; the cast surface 55% opaque
      (45% transparent). No labels are drawn.
    """
    import numpy as np
    from .figures import _pyplot,_save_figure
    root=Path(output_dir).resolve();root.mkdir(parents=True,exist_ok=True)
    if Path(prefix).name!=prefix or not prefix:raise ValueError('prefix must be a filename stem')
    paths={'residue_evidence_json':str(root/f'{prefix}_residue_evidence.json'),
           'residue_evidence_csv':str(root/f'{prefix}_residue_evidence.csv'),
           'residue_evidence_png':str(root/f'{prefix}_residue_evidence.png')}
    fields=['residue','group','role','boundary_area_A2','boundary_fraction','minimum_boundary_gap_A',
            'nonlocal_lining_contact_count','contacted_lining_boundary_fraction']
    with Path(paths['residue_evidence_csv']).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows({key:row[key] for key in fields} for row in report['residues'])
    paths['residue_evidence_partners_csv']=str(root/f'{prefix}_residue_evidence_partners.csv')
    partner_fields=['residue','lining_residue','interaction_candidate','center_distance_A','surface_distance_A','closest_atoms','lining_boundary_fraction']
    with Path(paths['residue_evidence_partners_csv']).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=partner_fields);writer.writeheader()
        for row in report['residues']:
            for partner in row.get('nonlocal_lining_partners') or []:
                pair=partner.get('closest_pair') or {}
                writer.writerow({'residue':row['residue'],'lining_residue':partner['residue'],
                                 'interaction_candidate':partner.get('interaction_candidate'),
                                 'center_distance_A':partner.get('center_distance_A'),
                                 'surface_distance_A':pair.get('surface_distance'),
                                 'closest_atoms':'-'.join(pair.get('atoms',[])),
                                 'lining_boundary_fraction':partner.get('boundary_fraction')})
    lining,support,_=select_stick_residues(report,STICK_RESIDUE_DEFAULT)
    stick_lining,stick_support,stick_selection=select_stick_residues(report,stick_residues)
    text=FigureText(f"Top {len(lining)} direct boundary contributors by assigned boundary area and top {len(support)} nonlocal "
                    "structural partners by boundary fraction of the lining residues they contact")
    plt=_pyplot();fig,axes=plt.subplots(1,2,figsize=(11,5.5),layout='constrained')
    for ax,rows,key,color,title,label in [(axes[0],lining,'boundary_fraction','#d98b32','Direct boundary contributors','Assigned boundary area (%)'),
              (axes[1],support,'contacted_lining_boundary_fraction','#a34477','Nonlocal structural partners','Boundary fraction of contacted lining residues (%)')]:
        ax.barh([r['residue'] for r in rows],[100*r[key] for r in rows],color=color)
        ax.invert_yaxis();ax.set_xlabel(label);text.title(ax,title);ax.spines[['top','right']].set_visible(False)
        if not rows:text.notice(ax,.5,.5,'No qualifying residues',transform=ax.transAxes,ha='center')
    text.suptitle(fig,'Geometric residue evidence — candidates for follow-up',fontsize=13)
    _save_figure(fig,Path(paths['residue_evidence_png']),dpi=dpi,text=text)
    report={**report,'display_lining_residues':[r['residue'] for r in stick_lining],
            'display_partner_residues':[r['residue'] for r in stick_support],
            'display_colors':{'boundary_lining':'gold','nonlocal_lining_partner':'magenta'},
            'stick_residue_selection':stick_selection}
    if scene_path is None and (volume_path or report.get('source_volume_dx')):
        from .residue_viewers import scene_from_binary_volume
        base_files=scene_from_binary_volume(frame,volume_path or report['source_volume_dx'],root,prefix)
        scene_path=base_files['volume_pml']
        paths.update(base_files)
    if scene_path is not None:
        scene=Path(scene_path).resolve(strict=True)
        viewer_ids,protein_path,mapping_definition=_viewer_serial_mapping(frame,scene)
        report['overlay_atom_mapping']=mapping_definition
        if scene.parent!=root:
            import shutil
            stem=scene.stem
            for ending in ('_volume','_pore_cast','_end_view'):
                if stem.endswith(ending):stem=stem[:-len(ending)];break
            for source in scene.parent.iterdir():
                if source.is_file() and source.name.startswith(stem+'_') and source.suffix in {'.pdb','.npz','.json','.dx','.bild'}:
                    target=root/source.name
                    if target.exists() and target.read_bytes()!=source.read_bytes():
                        raise FileExistsError(f'Refusing to replace different scene asset: {target}')
                    if not target.exists():shutil.copy2(source,target)
        import shutil
        for source in [protein_path,protein_path.with_suffix('.identities.json')]:
            if source.is_file() and source.parent!=root:
                target=root/source.name
                if target.exists() and target.read_bytes()!=source.read_bytes():
                    raise FileExistsError(f'Refusing to replace different protein asset: {target}')
                if not target.exists():shutil.copy2(source,target)
        membership={role:[serial for a,serial in zip(frame.atoms,viewer_ids) if a.residue_key.label in {r['residue'] for r in rows}]
                    for role,rows in [('lining',stick_lining),('partners',stick_support)]}
        coords_path=root/f'{prefix}_residue_reference.npz'
        np.savez_compressed(coords_path,protein_coords=np.asarray([a.coord for a in frame.atoms]),atom_serials=np.asarray(viewer_ids))
        evidence_path=Path(paths['residue_evidence_json'])
        extra=f"""\npython
# Background protein stays cartoon; identified residues also show opaque sticks.
crevice_reference = np.load(crevice_input({coords_path.name!r}, {str(coords_path)!r}), allow_pickle=False)
crevice_expected = crevice_reference["protein_coords"]
crevice_ids = crevice_reference["atom_serials"].tolist()
crevice_model = cmd.get_model("crevice_protein", state=1)
crevice_actual_by_id = {{a.id: a.coord for a in crevice_model.atom}}
if len(crevice_actual_by_id)!=len(crevice_model.atom) or set(crevice_actual_by_id)!=set(crevice_ids):
    raise ValueError("Residue overlay requires matching unique PDB atom serials")
crevice_actual = np.asarray([crevice_actual_by_id[i] for i in crevice_ids])
if not np.allclose(crevice_actual, crevice_expected, atol=.001, rtol=0):
    raise ValueError("Residue overlay requires matching protein coordinates by atom serial")
crevice_reference.close()
cmd.set("cartoon_color", -1, "crevice_protein")
cmd.set_color("crevice_boundary_evidence", [0.85, 0.55, 0.20])
cmd.set_color("crevice_partner_evidence", [0.64, 0.27, 0.47])
crevice_membership = {membership!r}
for role, color in [("partners", "crevice_partner_evidence"), ("lining", "crevice_boundary_evidence")]:
    atom_ids = crevice_membership[role]
    if atom_ids:
        crevice_selected = "crevice_protein and id " + "+".join(str(i) for i in atom_ids)
        cmd.color(color, crevice_selected)
        cmd.show("sticks", crevice_selected)
        cmd.set("stick_radius", 0.19, crevice_selected)
        cmd.set("stick_transparency", 0.0, crevice_selected)
cmd.set("stick_quality", 24, "crevice_protein")
for crevice_object in cmd.get_names("objects"):
    if crevice_object.startswith("crevice_volume"):
        cmd.set("cgo_transparency", 0.45, crevice_object)
python end
"""
        overlay=root/f'{prefix}_residue_context.pml'
        body=extra.split('python\n',1)[1].rsplit('python end',1)[0]
        base=scene.read_text()
        if base.count('python end')!=1:
            raise ValueError('The evidence overlay requires one CREVICE volume-scene Python block')
        overlay.write_text(base.replace('python end',body+'\npython end',1))
        report['display_atom_serials']=membership
        paths['residue_context_pml']=str(overlay);paths['residue_reference_npz']=str(coords_path)
        if stick_selection['mode']=='top' and stick_selection['request']==STICK_RESIDUE_DEFAULT:
            report['overlay_definition']='Top 12 direct contributors gold, top 12 nonlining nonlocal partners magenta, shown as opaque sticks; background cartoon and 55%-opaque cast.'
        else:
            report['overlay_definition']=(f"{stick_selection['definition']}: {len(stick_lining)} boundary residues gold, "
                f"{len(stick_support)} nonlining nonlocal partners magenta, shown as opaque sticks; background cartoon and 55%-opaque cast.")
        report['overlay_representation']='cartoon_with_evidence_sticks'
        report['overlay_cast_opacity']=.55
        from .residue_viewers import write_companion_contexts
        paths.update(write_companion_contexts(scene,root,prefix,protein_path,membership,viewer_ids,frame))
    Path(paths['residue_evidence_json']).write_text(json.dumps(report,indent=2)+'\n')
    return paths


@radii_option
def profile_constriction_evidence(frame,profile,*,tolerance=.25, radii: RadiusSet | str | None = None):
    """Atoms/residues limiting the continuous piecewise-linear profile path.

    Minimize distance to each selected atom sphere along every entire segment.
    These are geometric constraints on this path, not functional gate assignments.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile (same atom selection as the profile).
    profile : PoreProfile
        Resolved profile with at least two samples.
    tolerance : float, default 0.25
        Atoms whose smallest clearance is within this many Å of the global
        minimum are listed as limiting (one row per residue, its closest atom).
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
    dict
        ``status`` (``measured_path_constraints``), ``bottleneck_radius_A``
        (minimum clearance minus the profile's probe radius), ``raw_clearance_A``,
        ``tolerance_A``, ``probe_radius_A``, ``residues`` (``residue``, ``atom``,
        ``atom_serial``, ``raw_clearance_A``, ``path_position_A``,
        ``segment_index``; narrowest first), ``definition`` and
        ``interpretation``. For an unvalidated fixed-axis scan:
        ``{"status": "unavailable_unvalidated_axis_scan", "residues": []}``.

    Raises
    ------
    ValueError
        For a negative ``tolerance`` or a profile with fewer than two samples.
    """
    import numpy as np
    from .radii import atom_vdw_radius
    if not math.isfinite(tolerance) or tolerance<0:raise ValueError('tolerance must be finite and nonnegative')
    if profile.metadata.get('status')=='unvalidated_axis_scan':
        return {'status':'unavailable_unvalidated_axis_scan','residues':[]}
    atoms=frame.selected_atoms(include_hydrogen=profile.metadata.get('include_hydrogen',False),
                               include_hetero=profile.metadata.get('include_hetero',True))
    xyz=np.asarray([a.coord for a in atoms]);radii=np.asarray([atom_vdw_radius(a) for a in atoms])
    path=np.asarray([p.position for p in profile.points])
    if len(path)<2:raise ValueError('A continuous profile needs at least two points')
    minimum=np.full(len(atoms),np.inf);locations=np.zeros_like(xyz);segments=np.zeros(len(atoms),dtype=int)
    for j,(a,b) in enumerate(zip(path,path[1:])):
        v=b-a;den=float(v@v)
        fraction=np.clip((xyz-a)@v/den,0,1) if den>0 else np.zeros(len(atoms))
        closest=a+fraction[:,None]*v;gap=np.linalg.norm(xyz-closest,axis=1)-radii
        better=gap<minimum;minimum[better]=gap[better];locations[better]=closest[better];segments[better]=j
    raw=float(minimum.min());owners={}
    for i in np.flatnonzero(minimum<=raw+tolerance):
        atom=atoms[int(i)];row={'residue':atom.residue_key.label,'atom':atom.name,'atom_serial':atom.serial,
           'raw_clearance_A':float(minimum[i]),'path_position_A':locations[i].tolist(),'segment_index':int(segments[i])}
        if row['residue'] not in owners or row['raw_clearance_A']<owners[row['residue']]['raw_clearance_A']:
            owners[row['residue']]=row
    return {'status':'measured_path_constraints','bottleneck_radius_A':max(0.,raw-profile.probe_radius),
            'raw_clearance_A':raw,'tolerance_A':tolerance,'probe_radius_A':profile.probe_radius,
            'residues':sorted(owners.values(),key=lambda r:(r['raw_clearance_A'],r['residue'])),
            'definition':'Global minimum atom-sphere clearance over every continuous profile segment; tolerance includes near-limiting residues',
            'interpretation':'Geometric path constraints, not evidence of causal gating or functional control'}

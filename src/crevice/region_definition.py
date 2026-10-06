"""Versioned, identity-anchored regional observables and explicit review provenance.

A region definition (JSON) names a regional observable: the topology and
trajectory it belongs to (with SHA-256 hashes), the protein and obstacle atom
selections, the alignment scaffold, the grid geometry, the reference region
(an existing DX map, or a sphere around residue landmarks to be proposed),
the landmarks and their evidence, and a review status. ``crevice
region-init`` writes a candidate definition, ``crevice region-prepare``
(:func:`prepare_region`) fixes its reference region in the declared frame and
writes the review scenes, and ``crevice region-trajectory``
(:func:`trajectory_arguments`) measures the fixed region in every frame.

Schema versions
---------------
* ``schema_version: 1``: no radius set is stored. A version-1 definition is
  read as using :data:`crevice.radii.DEFAULT_RADII` (the standard element
  table plus CHARMM36 radii for recognised ion residues), the default since
  2 Oct 2026. Definitions prepared before then were prepared with the default
  of their day (before 30 Sep 2026 a CHARMM ``SOD`` was read as sulfur,
  1.80 Å; 30 Sep-1 Oct as neutral sodium, 2.27 Å); their fixed reference map
  is not recomputed, only the per-frame measurements use the new radii.
* ``schema_version: 2``: as version 1 plus a required top-level ``radii``
  object, the provenance of the radius set
  (:meth:`crevice.radii.RadiusSet.provenance`). ``region-prepare`` always
  writes version 2. Later steps rebuild that set
  (:func:`crevice.radii.radius_set_from_provenance`) and use it, or refuse a
  different explicitly chosen set unless ``allow_radii_mismatch`` /
  ``--allow-radii-mismatch`` is given (:func:`resolve_region_radii`).

A candidate is never promoted to an anatomical or functional assignment by this
module. Preparing a reference fixes its domain; it is not re-centred per frame.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import math
import re

from .presentation import annotate_option
from .radii import RadiusSet, radii_option

#: Region-definition schema versions this module reads. Version 2 adds the
#: required ``radii`` record; version 1 is read as the default radius set.
SUPPORTED_SCHEMA_VERSIONS = (1, 2)
#: The version :func:`prepare_region` writes.
CURRENT_SCHEMA_VERSION = 2


def write_json(path, value):
    """Write ``value`` as indented JSON (``allow_nan=False``) with a trailing newline.

    Parameters
    ----------
    path : str or Path
        Output file; overwritten.
    value : object
        JSON-serialisable value. NaN or infinity raises ``ValueError``.
    """
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def sha256(path):
    """SHA-256 hex digest of a file, read in 1 MiB blocks.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    str
        64 lower-case hex characters.
    """
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def file_record(path, base=None):
    """A definition asset record ``{"path", "sha256"}`` for an existing file.

    Parameters
    ----------
    path : str or Path
        Existing file.
    base : str or Path, optional
        Directory the recorded path is made relative to (normally the
        definition's directory); absolute path when omitted.

    Returns
    -------
    dict
        ``path`` (relative to ``base`` or absolute) and ``sha256``.
    """
    import os
    p = Path(path).resolve(strict=True)
    return {'path': os.path.relpath(p, base) if base else str(p), 'sha256': sha256(p)}


def _keys(value, required, optional=(), label='definition'):
    if not isinstance(value, dict):
        raise ValueError(f'{label} must be an object')
    missing = set(required)-value.keys()
    unknown = value.keys()-set(required)-set(optional)
    if missing or unknown:
        raise ValueError(f'{label}: missing fields {sorted(missing)}; unknown fields {sorted(unknown)}')


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be nonempty text')


def _positive(value, label, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int,float)) or not math.isfinite(value) or (value<0 if zero else value<=0):
        raise ValueError(f'{label} must be finite and {"nonnegative" if zero else "positive"}')


def _vector(value, label):
    if not isinstance(value, list) or len(value)!=3:
        raise ValueError(f'{label} requires three finite numbers')
    for v in value:
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
            raise ValueError(f'{label} requires three finite numbers')


def validate_definition(data):
    """Check a region definition against its schema; raise on the first problem.

    Parameters
    ----------
    data : dict
        Parsed definition JSON. Required keys: ``schema_version`` (1 or 2),
        ``region_id`` (portable file stem), ``label``, ``review`` (status
        ``candidate`` or ``reviewed``, rationale, sources; a reviewed region
        also needs reviewer, date and scope), ``system`` (topology and
        trajectory records, reference frame, selections, ``pbc`` ``check`` or
        ``none``, assembly, preparation and membrane notes), ``alignment``,
        ``geometry`` (grid spacing, probe radius, margin, axis, grid phase,
        profile padding, grid-point limit), ``reference`` (``existing_dx``
        with an asset record, or ``landmark_sphere``) and ``landmarks``.
        Optional: ``coordinate_landmarks`` and ``prepared``. Version 2 also
        requires ``radii`` (a radius-set provenance record with ``name`` and
        a 64-hex ``table_sha256``); version 1 must not have it.

    Returns
    -------
    dict
        ``data`` unchanged.

    Raises
    ------
    ValueError
        Naming the first missing, unknown or invalid field.
    """
    version=data.get('schema_version') if isinstance(data,dict) else None
    if type(version) is not int or version not in SUPPORTED_SCHEMA_VERSIONS:
        raise ValueError('Unsupported region schema_version; expected 1 or 2')
    required=['schema_version','region_id','label','review','system','alignment','geometry','reference','landmarks']
    _keys(data, required+(['radii'] if version>=2 else []), ['coordinate_landmarks','prepared'])
    if version>=2:
        radii=data['radii']
        if not isinstance(radii,dict) or not isinstance(radii.get('name'),str) or not re.fullmatch('[0-9a-f]{64}',str(radii.get('table_sha256'))):
            raise ValueError('radii must be a radius-set provenance record with a name and a table_sha256')
    if not isinstance(data['region_id'],str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}',data['region_id']):
        raise ValueError('region_id requires a portable alphanumeric filename stem')
    _text(data['label'],'label')
    review=data['review']
    _keys(review,['status','rationale','sources'],['reviewer','reviewed_at','scope'],'review')
    if review['status'] not in {'candidate','reviewed'}:raise ValueError('review status must be candidate or reviewed')
    _text(review['rationale'],'review rationale')
    if not isinstance(review['sources'],list):raise ValueError('review sources must be a list')
    for s in review['sources']:_text(s,'review source')
    if review['status']=='reviewed':
        for key in ['reviewer','reviewed_at','scope']:_text(review.get(key),f'reviewed regions require {key}')
        if not review['sources']:raise ValueError('Reviewed regions require supporting sources')
    system=data['system']
    _keys(system,['topology','trajectory','reference_frame','protein_selection','obstacle_selection','pbc','assembly','preparation','membrane'],label='system')
    if type(system['reference_frame']) is not int or system['reference_frame']<0:raise ValueError('reference_frame must be a nonnegative integer')
    if system['pbc'] not in {'check','none'}:raise ValueError('Region definitions require preprocessed whole structures; pbc is check or none')
    for k in ['protein_selection','obstacle_selection','assembly','preparation','membrane']:_text(system[k],k)
    alignment=data['alignment'];_keys(alignment,['residues','rationale'],label='alignment')
    _text(alignment['rationale'],'alignment rationale')
    ids=alignment['residues']
    if ids is not None:
        if not isinstance(ids,list) or not ids or any(not isinstance(x,str) or not x for x in ids) or len(set(ids))!=len(ids):
            raise ValueError('alignment residues must be null or unique full residue IDs')
    g=data['geometry']
    _keys(g,['spacing_A','probe_radius_A','region_margin_A','axis','grid_phase','profile_padding_A','max_grid_points'],label='geometry')
    for k in ['spacing_A','probe_radius_A']:_positive(g[k],k)
    for k in ['region_margin_A','profile_padding_A']:_positive(g[k],k,zero=True)
    if type(g['max_grid_points']) is not int or g['max_grid_points']<1:raise ValueError('max_grid_points must be a positive integer')
    _vector(g['axis'],'axis');_vector(g['grid_phase'],'grid_phase')
    if sum(x*x for x in g['axis'])==0:raise ValueError('axis must be nonzero')
    if any(not 0<=x<1 for x in g['grid_phase']):raise ValueError('grid_phase fractions must be in [0,1)')
    ref=data['reference']
    if not isinstance(ref,dict):raise ValueError('reference must be an object')
    if ref.get('kind')=='existing_dx':_keys(ref,['kind','asset'],label='reference')
    elif ref.get('kind')=='landmark_sphere':
        _keys(ref,['kind','radius_A','max_seed_distance_A','ambiguity_distance_A'],label='reference')
        for k in ['radius_A','max_seed_distance_A']:_positive(ref[k],k)
        _positive(ref['ambiguity_distance_A'],'ambiguity_distance_A',zero=True)
        if ref['max_seed_distance_A']>ref['radius_A']:raise ValueError('Seed search distance must not exceed the declared region radius')
        if review['status']!='candidate':raise ValueError('Automatic landmark proposals must remain candidate; review a prepared definition explicitly')
    else:raise ValueError('reference kind must be existing_dx or landmark_sphere')
    if not isinstance(data['landmarks'],list):raise ValueError('landmarks must be a list')
    seen=set()
    for row in data['landmarks']:
        _keys(row,['residue','atoms','evidence'],label='landmark')
        _text(row['residue'],'landmark residue');_text(row['evidence'],'landmark evidence')
        names=row['atoms']
        if not isinstance(names,list) or not names or any(not isinstance(x,str) or not x for x in names) or len(names)!=len(set(names)):
            raise ValueError('Each landmark requires unique explicit atom names')
        if row['residue'] in seen:raise ValueError('Duplicate landmark residue')
        seen.add(row['residue'])
    if ref['kind']=='landmark_sphere' and not seen:raise ValueError('Landmark proposal requires residue landmarks')
    for row in data.get('coordinate_landmarks',[]):
        _keys(row,['label','points_A','evidence','transfer_rmsd_A'],label='coordinate landmark')
        _text(row['label'],'coordinate landmark label');_text(row['evidence'],'coordinate landmark evidence')
        if not isinstance(row['points_A'],list) or not row['points_A']:raise ValueError('Coordinate landmark points must be nonempty')
        for p in row['points_A']:_vector(p,'coordinate landmark point')
        if row['transfer_rmsd_A'] is not None:_positive(row['transfer_rmsd_A'],'transfer RMSD',zero=True)
    return data


def load_definition(path):
    """Read and validate a region definition and verify its input files.

    Duplicate JSON keys are refused. The topology, trajectory (unless
    ``null``) and, for an ``existing_dx`` reference, the reference map are
    resolved relative to the definition's directory and their SHA-256 must
    equal the recorded digests.

    Parameters
    ----------
    path : str or Path
        Definition JSON.

    Returns
    -------
    data : dict
        The validated definition.
    assets : dict of str to Path or None
        Resolved ``topology``, ``trajectory`` (``None`` for a static input)
        and, if present, ``reference`` paths.

    Raises
    ------
    ValueError
        For an invalid definition or a file whose hash differs.
    """
    path=Path(path).resolve(strict=True)
    def unique(pairs):
        out={}
        for key,value in pairs:
            if key in out:raise ValueError(f'Duplicate JSON key {key}')
            out[key]=value
        return out
    data=json.loads(path.read_text(),object_pairs_hook=unique)
    validate_definition(data)
    assets={}
    records={'topology':data['system']['topology'],'trajectory':data['system']['trajectory']}
    if data['reference']['kind']=='existing_dx':records['reference']=data['reference']['asset']
    for key,record in records.items():
        if key=='trajectory' and record is None:assets[key]=None;continue
        _keys(record,['path','sha256'],label=key+' asset')
        _text(record['path'],key+' path')
        if not isinstance(record['sha256'],str) or not re.fullmatch('[0-9a-f]{64}',record['sha256']):raise ValueError(key+' asset needs a SHA256 digest')
        p=(path.parent/record['path']).resolve(strict=True)
        if sha256(p)!=record['sha256']:raise ValueError(key+' asset SHA256 differs from the region definition')
        assets[key]=p
    return data,assets


def declared_radii(data):
    """The radius-set record a definition declares, and how it was obtained.

    Parameters
    ----------
    data : dict
        A validated definition.

    Returns
    -------
    record : dict or None
        The ``radii`` provenance of a version-2 definition; for a prepared
        version-1 definition, the provenance of
        :data:`crevice.radii.DEFAULT_RADII` (the documented meaning of an
        older definition); ``None`` for an unprepared version-1 request, which
        has not chosen a set yet.
    origin : str
        ``"definition"``, ``"legacy_v1_default"`` or ``"not_declared"``.
    """
    from .radii import DEFAULT_RADII
    if data['schema_version']>=2:return data['radii'],'definition'
    if 'prepared' in data:return DEFAULT_RADII.provenance(),'legacy_v1_default'
    return None,'not_declared'


def resolve_region_radii(data, requested=None, *, allow_radii_mismatch=False, base_dir=None):
    """Choose the radius set for work on a region, checking it against the definition.

    Rules:

    * The definition declares a set (version 2, or a prepared version 1,
      which means :data:`crevice.radii.DEFAULT_RADII`) and nothing else was
      chosen: that set is used (rebuilt with
      :func:`crevice.radii.radius_set_from_provenance`).
    * A set was chosen explicitly (``requested``: ``--radii`` or ``radii=``)
      and has the same table hash as the declared one: it is used.
    * The two differ: ``ValueError`` naming both, unless
      ``allow_radii_mismatch`` is true, in which case the requested set is
      used and the mismatch is recorded.
    * The definition declares nothing (an unprepared version-1 request): the
      requested set, or the default set, is used.

    Parameters
    ----------
    data : dict
        A validated definition.
    requested : RadiusSet, str, os.PathLike or None, optional
        The explicitly chosen set; ``None`` means none was chosen (for the
        CLI this is :func:`crevice.radii.active_radii`).
    allow_radii_mismatch : bool, default False
        Use ``requested`` even if it differs from the declared set.
    base_dir : str or Path, optional
        Directory for a relative radius-file path in the declared record.

    Returns
    -------
    radii : RadiusSet
        The set to use.
    check : dict
        ``status`` (``"definition_set_used"``, ``"requested_matches_definition"``,
        ``"mismatch_overridden"``, ``"no_declared_set"``), ``declared_origin``,
        ``declared`` and ``used`` (name and table hash of each).

    Raises
    ------
    ValueError
        On a mismatch without ``allow_radii_mismatch``, or if the declared set
        cannot be rebuilt (pass the same set explicitly then).
    """
    from .radii import effective_radii,radius_set_from_provenance,radius_set_label,resolve_radii,same_radius_set
    declared,origin=declared_radii(data)
    chosen=resolve_radii(requested)
    def short(value):
        if value is None:return None
        if isinstance(value,RadiusSet):return {'name':value.name,'table_sha256':value.table_sha256}
        return {'name':value.get('name'),'table_sha256':value.get('table_sha256')}
    if declared is None:
        used=effective_radii(chosen);status='no_declared_set'
    elif chosen is None:
        used=radius_set_from_provenance(declared,base_dir=base_dir);status='definition_set_used'
    elif same_radius_set(chosen,declared):
        used=chosen;status='requested_matches_definition'
    elif allow_radii_mismatch:
        used=chosen;status='mismatch_overridden'
    else:
        raise ValueError(f'Radius set mismatch: the region definition declares {radius_set_label(declared)} '
                         f'({origin}) but {radius_set_label(chosen)} was requested. Omit --radii to use the '
                         'definition\'s set, prepare a new definition with the other set, or pass '
                         '--allow-radii-mismatch to measure with the requested set anyway (recorded).')
    return used,{'status':status,'declared_origin':origin,'declared':short(declared),'used':short(used)}


def frame_contract(frame):
    """Fingerprint of a reference frame's atoms and coordinates.

    Parameters
    ----------
    frame : StructureFrame
        The reference frame (protein selection).

    Returns
    -------
    dict
        ``source_frame_index``, ``atom_count``, ``identities_sha256`` (residue
        label, atom name, altloc, element and hetero flag of every atom, in
        order), ``coordinates_float32_sha256`` (little-endian float32
        coordinates) and ``coordinate_convention``. ``region-trajectory``
        requires the current inputs to reproduce it exactly.
    """
    import numpy as np
    identities=[[a.residue_key.label,a.name,a.altloc,a.element,a.hetero] for a in frame.atoms]
    xyz=np.asarray([a.coord for a in frame.atoms],dtype='<f4')
    return {'source_frame_index':frame.frame_index,'atom_count':len(identities),
            'identities_sha256':hashlib.sha256(json.dumps(identities,separators=(',',':')).encode()).hexdigest(),
            'coordinates_float32_sha256':hashlib.sha256(xyz.tobytes()).hexdigest(),
            'coordinate_convention':'Angstrom, input reference frame; no coordinate transfer applied'}


def resolve_landmarks(frame, rows):
    """Resolve residue landmarks to atoms by exact identity.

    Each landmark's residue label (chain, residue name, number, insertion
    code) and each declared atom name must match exactly one atom of
    ``frame``; hydrogens are refused.

    Parameters
    ----------
    frame : StructureFrame
        Reference frame.
    rows : list of dict
        Definition landmarks (``residue``, ``atoms``, ``evidence``).

    Returns
    -------
    list of dict
        Each row plus ``points_A`` (declared atom coordinates),
        ``centroid_A``, ``anchor_atom_serials`` and ``residue_atom_serials``
        (1-based positions in ``frame``).

    Raises
    ------
    ValueError
        If an atom is missing or ambiguous, or is a hydrogen.
    """
    import numpy as np
    result=[]
    for row in rows:
        selected=[];serials=[]
        for name in row['atoms']:
            matches=[(i,a) for i,a in enumerate(frame.atoms,1) if a.residue_key.label==row['residue'] and a.name==name]
            if len(matches)!=1:raise ValueError(f'Landmark {row["residue"]}/{name} resolves to {len(matches)} atoms; explicit identity required')
            i,a=matches[0]
            if a.element.upper() in {'H','D'}:raise ValueError('Region landmarks require heavy atoms')
            selected.append(a.coord);serials.append(i)
        all_serials=[i for i,a in enumerate(frame.atoms,1) if a.residue_key.label==row['residue']]
        result.append({**row,'points_A':selected,'centroid_A':np.mean(selected,axis=0).tolist(),
                       'anchor_atom_serials':serials,'residue_atom_serials':all_serials})
    return result


def propose_landmark_region(frame, landmarks, settings, geometry, obstacle_source):
    """One locally seeded probe component inside an explicitly bounded sphere.

    The anchor is the equal-weight mean of the landmark centroids (each
    residue counts once, whatever its atom count). On a grid of the
    definition's spacing, axis and phase around the anchor, probe centres
    (clearance at least the probe radius to the obstacle atoms, with radii from
    the radius set in effect) inside the declared sphere are connected along
    atom-clear segments; the component nearest the anchor is kept, swept by the
    probe and clipped to the sphere. Seeds are candidate empty-space samples,
    not ligand poses or access claims.

    Parameters
    ----------
    frame : StructureFrame
        Reference frame (protein).
    landmarks : list of dict
        Resolved landmarks (:func:`resolve_landmarks`).
    settings : dict
        The ``landmark_sphere`` reference: ``radius_A`` (sphere radius),
        ``max_seed_distance_A`` (largest allowed anchor-to-component
        distance) and ``ambiguity_distance_A`` (two components closer to each
        other's distance than this are ambiguous).
    geometry : dict
        The definition's ``geometry`` (spacing, probe radius, axis, grid
        phase, grid-point limit).
    obstacle_source : object
        Obstacle atoms with ``atoms``, ``radii`` and ``snapshot(index)``.

    Returns
    -------
    cast : VoidCast
        The proposed region (one component, ``landmark_sphere_candidate``).
    diagnostic : dict
        Anchor, candidates with distances, the selected one, grid size, crop
        boundary samples, minimum clearance, obstacle report, definition and
        limitations.

    Raises
    ------
    ValueError
        If the grid exceeds ``max_grid_points``, no admissible component lies
        within ``max_seed_distance_A``, the nearest two are ambiguous, or the
        swept region is empty.
    """
    import numpy as np
    from scipy import ndimage
    from .geometry import basis_from_direction
    from .rolling import _SphereQueries,_connected
    from .radii import atom_vdw_radius
    from .cavity_obstacles import aligned_obstacle_images
    from .models import ChannelPoint,VoidComponent,VoidCast
    anchor=np.asarray([r['centroid_A'] for r in landmarks]).mean(0)
    h=geometry['spacing_A'];probe=geometry['probe_radius_A'];radius=settings['radius_A']
    axis=np.asarray(geometry['axis'],float);axis/=np.linalg.norm(axis)
    u,v=basis_from_direction(tuple(axis));basis=np.array([u,v,axis])
    half=math.ceil((radius+probe+2*h)/h);shape=(2*half+1,)*3
    if math.prod(shape)>geometry['max_grid_points']:raise ValueError('Landmark proposal exceeds max_grid_points')
    index=np.indices(shape).reshape(3,-1).T-half
    phase=np.asarray(geometry['grid_phase']);xyz=(index+phase)*h@basis+anchor
    dist=np.linalg.norm(xyz-anchor,axis=1).reshape(shape)
    # Every obstacle that could alter sample clearance or an eligible grid edge.
    pq=_SphereQueries(np.asarray([a.coord for a in frame.atoms]),np.asarray([atom_vdw_radius(a) for a in frame.atoms]))
    pad=max(float(pq.points(xyz)[0].max()),probe+h)+float(obstacle_source.radii.max())+1e-6
    coords,box=obstacle_source.snapshot(frame.frame_index)
    obstacles,obstacle_report=aligned_obstacle_images(obstacle_source.atoms,coords,xyz.min(0),xyz.max(0),rotation=np.eye(3),translation=np.zeros(3),padding=pad,box=box)
    queries=_SphereQueries(np.asarray([a.coord for a in obstacles])@basis.T,np.asarray([atom_vdw_radius(a) for a in obstacles]))
    clearance,nearest=queries.points(xyz@basis.T);clearance=clearance.reshape(shape)
    origin=xyz[0]@basis.T
    centers=(dist<=radius)&(clearance>=probe)
    labels=_connected(centers,clearance,origin,h,probe,queries)
    candidates=[]
    for label in np.unique(labels[centers]):
        mask=labels==label;flat=np.flatnonzero(mask);j=int(flat[np.argmin(dist.ravel()[flat])])
        candidates.append({'component':int(label),'nearest_distance_A':float(dist.ravel()[j]),'seed_A':xyz[j].tolist(),'seed_clearance_A':float(clearance.ravel()[j]),'probe_centers':int(mask.sum())})
    candidates.sort(key=lambda r:(r['nearest_distance_A'],r['component']))
    if not candidates or candidates[0]['nearest_distance_A']>settings['max_seed_distance_A']:
        raise ValueError('No admissible probe center within the declared landmark seed distance; no reference created')
    if len(candidates)>1 and candidates[1]['nearest_distance_A']-candidates[0]['nearest_distance_A']<=settings['ambiguity_distance_A']:
        raise ValueError('Ambiguous landmark seed: competing disconnected components are similarly close; no reference created')
    selected=labels==candidates[0]['component']
    swept=ndimage.distance_transform_edt(~selected,sampling=h)<=probe+1e-9
    filled=swept&(dist<=radius+1e-9)&(clearance>=0)
    flat=np.flatnonzero(filled);points=tuple(ChannelPoint(i,tuple(map(float,xyz[j])),float(clearance.ravel()[j]),float(clearance.ravel()[j]),float((xyz[j]-anchor)@axis),obstacles[int(nearest[j])].serial,obstacles[int(nearest[j])].residue_key.label) for i,j in enumerate(flat))
    if not points:raise ValueError('Landmark proposal has no atom-clear swept samples')
    component=VoidComponent(1,'landmark_sphere_candidate',tuple(np.mean([p.position for p in points],axis=0)),len(points)*h**3,points)
    diagnostic={'anchor_A':anchor.tolist(),'anchor_definition':'equal-weight mean of local residue landmark atom centroids',
      'radius_A':radius,'candidates':candidates,'selected':candidates[0],'grid_points':math.prod(shape),
      'crop_boundary_samples':int(np.sum(filled&(dist>=radius-.9*h))),
      'minimum_sample_clearance_A':float(clearance[filled].min()),'obstacles':obstacle_report,
      'definition':'Single nearest admissible probe-center component in a landmark-centred sphere; probe sweeps clipped to that sphere and atom-clear sample centers',
      'limitations':['The sphere is a declared crop, not an anatomical cavity extent or mouth.',
                     'A nearby unoccupied probe center is not a ligand pose or evidence of substrate access.',
                     'No absent sampled connection proves continuum disconnection.']}
    return VoidCast((component,),points,h,0.,'rolling',{'grid_basis':basis.tolist(),'algorithm':'landmark_sphere_candidate'}),diagnostic


class StaticMmcifRegionSource:
    """Use canonical mmCIF identities, with MDAnalysis only for atom selections.

    Deposited crystal symmetry and assembly operators are never applied silently.
    This adapter deliberately does not treat the crystal cell as an MD box.

    Parameters
    ----------
    topology : str or Path
        mmCIF file (read with :func:`crevice.parser.load_structure_report`).
    system : dict
        The definition's ``system`` block; ``reference_frame`` must be 0, and
        ``protein_selection``/``obstacle_selection`` are MDAnalysis selections
        evaluated on an in-memory universe built from the parsed atoms.

    Attributes
    ----------
    frame : StructureFrame
        The selected protein atoms.
    atoms : tuple of Atom
        Obstacle heavy atoms (protein atoms first keep their identity; other
        atoms get an ``/obstacle/`` chain label and ``hetero=True``).
    radii : numpy.ndarray
        Obstacle radii from the radius set in effect when the source is built.
    coords : numpy.ndarray
        Obstacle coordinates (Å).
    metadata : dict
        Obstacle selection, atom counts and periodicity.
    reader : dict
        Trajectory-reader report.

    Raises
    ------
    ValueError
        For a reference frame other than 0, an empty protein selection, or an
        obstacle selection that omits analysed protein heavy atoms.
    """
    def __init__(self, topology, system):
        import numpy as np
        import MDAnalysis as mda
        from dataclasses import replace
        from collections import Counter
        from .parser import load_structure_report
        from .radii import atom_vdw_radius
        if system['reference_frame']!=0:raise ValueError('Static mmCIF region input currently requires model/reference frame 0')
        original,reader=load_structure_report(topology)
        keys=list(dict.fromkeys(a.residue_key for a in original.atoms));lookup={k:i for i,k in enumerate(keys)}
        chains=list(dict.fromkeys(k.chain_id for k in keys))
        u=mda.Universe.empty(len(original.atoms),n_residues=len(keys),n_segments=len(chains),
            atom_resindex=[lookup[a.residue_key] for a in original.atoms],
            residue_segindex=[chains.index(k.chain_id) for k in keys],trajectory=True)
        u.add_TopologyAttr('names',[a.name for a in original.atoms]);u.add_TopologyAttr('elements',[a.element for a in original.atoms])
        u.add_TopologyAttr('resnames',[k.resname for k in keys]);u.add_TopologyAttr('resids',[k.resid for k in keys])
        u.add_TopologyAttr('icodes',[k.icode for k in keys]);u.add_TopologyAttr('segids',chains)
        u.add_TopologyAttr('chainIDs',[a.chain_id for a in original.atoms]);u.atoms.positions=[a.coord for a in original.atoms]
        indices=list(map(int,u.select_atoms(system['protein_selection']).indices))
        if not indices:raise ValueError('Region protein selection is empty')
        protein={i:replace(original.atoms[i],serial=j+1) for j,i in enumerate(indices)}
        self.frame=replace(original,atoms=tuple(protein.values()),frame_index=0)
        chosen=[int(i) for i in u.select_atoms(system['obstacle_selection']).indices if original.atoms[int(i)].element.upper() not in {'H','D'}]
        if not set(indices).issubset(chosen):raise ValueError('Obstacle selection must retain every analysed protein heavy atom')
        self.atoms=tuple(protein[i] if i in protein else replace(original.atoms[i],serial=i+1,
            chain_id=f'{original.atoms[i].chain_id}/obstacle/{lookup[original.atoms[i].residue_key]}',hetero=True) for i in chosen)
        self.radii=np.asarray([atom_vdw_radius(a) for a in self.atoms]);self.coords=np.asarray([a.coord for a in self.atoms])
        self.metadata={'selection':system['obstacle_selection'],'selected_heavy_atom_count':len(chosen),
            'analysed_protein_heavy_atom_count':len(indices),'other_heavy_atom_count':len(chosen)-len(indices),
            'residue_atom_counts':dict(Counter(a.resname for a in self.atoms)),
            'periodicity':'Static deposited coordinates; crystal symmetry and assembly operators not applied'}
        self.reader={**reader,'selection':system['protein_selection'],'selected_atom_count':len(indices),
            'selection_engine':'MDAnalysis over canonical parser atoms; original residue identities retained',
            'periodicity':self.metadata['periodicity']}

    def snapshot(self,index):
        """Obstacle coordinates of frame ``index`` (only 0 exists) and no box."""
        if index!=0:raise ValueError('Static obstacle source has only reference frame 0')
        return self.coords.copy(),None

    def prepare(self,*args,**kwargs):
        """Obstacle atoms around the reference, as :meth:`crevice.cavity_obstacles.ObstacleSource.prepare`.

        Parameters
        ----------
        *args, **kwargs
            Passed unchanged to
            :meth:`crevice.cavity_obstacles.ObstacleSource.prepare` (reference
            and aligned frames, fit, reference region, margin, probe radius,
            grid-point limit).

        Returns
        -------
        tuple
            Obstacle atoms and their report, as that method returns them.
        """
        from .cavity_obstacles import ObstacleSource
        return ObstacleSource.prepare(self,*args,**kwargs)


@radii_option
@annotate_option
def prepare_region(path, output_dir, *, annotate=None, radii: RadiusSet | str | None = None,
                   allow_radii_mismatch=False):
    """Validate a region definition, fix its reference and write review scenes.

    The definition's topology/trajectory hashes are checked, the declared
    reference frame is read, landmarks are resolved to atoms and the reference
    region is either copied (``existing_dx``) or proposed around the landmarks
    (``landmark_sphere``). Diagnostics (volume, obstacle overlap, landmark
    distances, warnings) and a prepared, versioned definition are written, and
    :func:`crevice.region_review.write_region_review` writes the native review
    scenes and page (``crevice region-prepare``).

    The radius set is chosen by :func:`resolve_region_radii`: an unprepared
    request (from ``region-init``) takes ``radii`` or the set in effect; a
    definition that already declares a set uses it, and a different
    explicitly chosen set is refused unless ``allow_radii_mismatch``. The
    prepared definition is written as schema version 2 with that set's
    provenance in ``radii``, so every later step measures with the same radii.

    Parameters
    ----------
    path : str or Path
        Region definition JSON (for example from ``crevice region-init``).
    output_dir : str or Path
        New or empty output directory.
    annotate : bool, optional
        Draw residue-landmark labels, the status caption and the legend line in
        the review scenes. ``None`` (default) inherits the surrounding setting,
        which is off unless enabled with
        :func:`crevice.presentation.figure_annotations` or ``--annotate``.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`. The set used is stored in the prepared
        definition (``radii``, schema version 2).
    allow_radii_mismatch : bool, default False
        If the input definition already declares a radius set (schema 2, or a
        prepared schema 1, which means the default set), measure with ``radii``
        even though it differs; the mismatch is recorded as ``radii_check`` in
        the review JSON and manifest.

    Returns
    -------
    dict of str to str
        Output keys and paths, including ``review_scene_json`` and ``review_landmarks_csv``.

    Raises
    ------
    ValueError
        For an invalid definition, a non-empty output directory, unresolvable
        landmarks, an ambiguous or missing landmark seed, or a radius-set
        mismatch without ``allow_radii_mismatch``.

    Outputs
    -------
    request.json, reader.json : JSON
        The submitted definition and trajectory-reader provenance.
    REGION_region.json, REGION_review.json, REGION_reference.dx : definition
        Prepared definition, review diagnostics and fixed reference grid.
    REGION_review.{pml,vmd,tcl,cxc}, REGION_review_landmarks.csv : review
        See :func:`crevice.region_review.write_region_review`.
    manifest.json : JSON
        Provenance: region ID, review status, every file, ``radii`` (the radius set used)
        and ``radii_check`` (how it was chosen).

    Colours and representation
    --------------------------
    As :func:`crevice.region_review.write_region_review`: gold landmark sticks,
    purple transferred points, green nearby obstacles, teal reference surface,
    grey-blue cartoon; no text unless annotating.
    """
    from .radii import active_radii,use_radii
    data,_=load_definition(path)
    used,check=resolve_region_radii(data,active_radii(),allow_radii_mismatch=allow_radii_mismatch,
                                    base_dir=Path(path).resolve().parent)
    with use_radii(used):
        return _prepare_region(path,output_dir,used,check)


def _prepare_region(path, output_dir, used, check):
    """Body of :func:`prepare_region`, run with the chosen radius set ``used`` in effect."""
    import os
    import numpy as np
    from scipy.spatial import cKDTree
    from .trajectory import load_trajectory_report,align_frame_report
    from .cavity_trajectory import CavityReference
    from .cavity_obstacles import ObstacleSource
    from .region_review import write_region_review
    data,assets=load_definition(path);root=Path(output_dir).resolve()
    if root.exists() and any(root.iterdir()):raise ValueError('Region preparation requires an empty output directory; existing evidence is preserved')
    system=data['system'];index=system['reference_frame']
    if assets['trajectory'] is None and assets['topology'].suffix.lower() in {'.cif','.mmcif'}:
        source=StaticMmcifRegionSource(assets['topology'],system);frame=source.frame;reader=source.reader
    else:
        trajectory,reader=load_trajectory_report(assets['topology'],assets['trajectory'],selection=system['protein_selection'],start=index,stop=index+1,pbc=system['pbc'])
        frame=trajectory.frames[0]
        source=ObstacleSource(assets['topology'],assets['trajectory'],protein_selection=system['protein_selection'],obstacle_selection=system['obstacle_selection'],reference_frame=frame,frame_times={index:trajectory.times()[0]},pbc=system['pbc'])
    if any(a.hetero or a.element.upper() in {'H','D'} for a in frame.atoms):raise ValueError('Region protein_selection must contain protein heavy atoms only')
    landmarks=resolve_landmarks(frame,data['landmarks'])
    _,fit=align_frame_report(frame,frame,residues=data['alignment']['residues'])
    root.mkdir(parents=True,exist_ok=True)
    write_json(root/'request.json',data)
    write_json(root/'reader.json',reader)
    prefix=data['region_id'];g=data['geometry']
    if data['reference']['kind']=='landmark_sphere':
        from .volume_export import write_void_cast_dx
        cast,proposal=propose_landmark_region(frame,landmarks,data['reference'],g,source)
        reference_path=write_void_cast_dx(cast,root/(prefix+'_reference.dx'))
    else:
        import shutil
        proposal=None;reference_path=root/(prefix+'_reference.dx');shutil.copy2(assets['reference'],reference_path)
    reference=CavityReference.from_dx(reference_path,spacing=g['spacing_A'],axis=g['axis'],grid_phase=g['grid_phase'],profile_padding=g['profile_padding_A'])
    obstacle_atoms,obstacle_report=source.prepare(frame,frame,fit,reference,margin=g['region_margin_A'],probe_radius=g['probe_radius_A'],max_grid_points=g['max_grid_points'])
    from .rolling import _SphereQueries
    from .radii import atom_vdw_radius
    queries=_SphereQueries(np.asarray([a.coord for a in obstacle_atoms]),np.asarray([atom_vdw_radius(a) for a in obstacle_atoms]))
    clearance,_=queries.points(reference.points)
    tree=cKDTree(reference.points)
    for row in landmarks:row['minimum_distance_to_reference_samples_A']=float(tree.query(row['points_A'])[0].min())
    coordinates=[]
    for row in data.get('coordinate_landmarks',[]):
        distances=tree.query(row['points_A'])[0]
        coordinates.append({**row,'minimum_distance_to_reference_samples_A':float(distances.min()),
                            'points_inside_reference_neighbourhood':int(np.sum(distances<=g['region_margin_A']))})
    diagnostics={'schema_version':1,'region_id':prefix,'label':data['label'],'review':data['review'],
        'reference_frame':frame_contract(frame),'alignment':fit,'landmarks':landmarks,'coordinate_landmarks':coordinates,
        'reference_volume_A3':float(np.sum(reference.grid)*abs(np.linalg.det(reference.deltas))),
        'reference_samples_inside_obstacle_spheres':int(np.sum(clearance < -1e-6)),
        'reference_minimum_clearance_A':float(clearance.min()),'proposal':proposal,'obstacles':obstacle_report,
        'assembly':system['assembly'],'preparation':system['preparation'],'membrane':system['membrane'],
        'warnings':['Anatomical identity, membrane sides and ligand access require independent review.',
                    'Reference coordinates stay fixed after protein alignment; landmarks are not used to chase a void each frame.']}
    if diagnostics['reference_samples_inside_obstacle_spheres']:
        diagnostics['warnings'].append('Reference samples overlap the selected obstacles; they are an anchor domain, not a current free-volume measurement. Each trajectory frame is recalculated.')
    if any(r['transfer_rmsd_A'] is not None and r['transfer_rmsd_A']>2 for r in coordinates):
        diagnostics['warnings'].append('Coordinate transfer RMSD exceeds the 2 A display-review flag; this threshold is a review prompt, not a biological acceptance criterion.')
    diagnostics['radii']=used.provenance();diagnostics['radii_check']=check
    prepared=deepcopy(data)
    prepared['schema_version']=CURRENT_SCHEMA_VERSION;prepared['radii']=used.provenance()
    for key in ['topology','trajectory']:
        if assets[key] is not None:prepared['system'][key]={**system[key],'path':os.path.relpath(assets[key],root)}
    prepared['reference']={'kind':'existing_dx','asset':file_record(reference_path,root)}
    prepared['prepared']={'request_sha256':sha256(path),'frame_contract':frame_contract(frame),'proposal':proposal,
                           'definition':'Reference fixed in the declared source frame; preserve separate region IDs for alternative assignments'}
    definition=root/(prefix+'_region.json');write_json(definition,prepared)
    write_json(root/(prefix+'_review.json'),diagnostics)
    files=write_region_review(frame,reference,diagnostics,obstacle_atoms,root,prefix)
    files.update(definition_json=str(definition),review_json=str(root/(prefix+'_review.json')),reference_dx=str(reference_path))
    from .radii import radii_fields
    write_json(root/'manifest.json',{'region_id':prefix,'review_status':data['review']['status'],'files':files,
                                     **radii_fields(),'radii_check':check})
    return files


def trajectory_arguments(path, output_dir, *, radii=None, allow_radii_mismatch=False):
    """Resolve a prepared region into the existing cavity trajectory CLI.

    The definition must be prepared (``existing_dx`` reference, ``prepared``
    record), declare a trajectory, use reference frame 0, and its frame
    contract (atom identities and float32 coordinates of the reference frame)
    must match the current input files exactly.

    Parameters
    ----------
    path : str or Path
        Prepared definition (``REGION_region.json`` from ``region-prepare``).
    output_dir : str or Path
        Output directory passed on as ``--out-dir``.
    radii : RadiusSet, str, os.PathLike or None, optional
        An explicitly chosen radius set (``--radii``); ``None`` uses the set the
        definition declares (version 2) or, for version 1, the default set.
    allow_radii_mismatch : bool, default False
        Measure with ``radii`` even if it differs from the declared set.

    Returns
    -------
    argv : list of str
        ``crevice cavity-trajectory`` arguments: topology, trajectory,
        reference map, selections, PBC, ``--geometry-mode reference-region``
        and the definition's grid geometry.
    metadata : dict
        Region provenance (schema version, region ID, label, review status,
        definition SHA-256 and content, ``meaning``), plus ``radii`` (the set
        to measure with, provenance), ``radii_check`` (see
        :func:`resolve_region_radii`) and ``radius_set`` (the
        :class:`~crevice.radii.RadiusSet` object itself; remove it before
        writing JSON).

    Raises
    ------
    ValueError
        If the definition is not prepared, has no trajectory, uses another
        reference frame, its frame contract differs, or the radius sets differ
        without ``allow_radii_mismatch``.
    """
    data,assets=load_definition(path)
    if 'prepared' not in data or data['reference']['kind']!='existing_dx':raise ValueError('Run region-prepare before region-trajectory')
    if assets['trajectory'] is None:raise ValueError('region-trajectory requires a declared trajectory')
    if data['system']['reference_frame']!=0:raise ValueError('region-trajectory currently requires source reference frame 0')
    from .trajectory import load_trajectory_report
    t,_=load_trajectory_report(assets['topology'],assets['trajectory'],selection=data['system']['protein_selection'],stop=1,pbc=data['system']['pbc'])
    resolve_landmarks(t.frames[0],data['landmarks'])
    if frame_contract(t.frames[0])!=data['prepared']['frame_contract']:raise ValueError('Prepared reference atom identities or coordinate frame differ')
    g=data['geometry'];s=data['system']
    argv=['cavity-trajectory',str(assets['topology']),str(assets['trajectory']),'--reference-volume-dx',str(assets['reference']),
          '--out-dir',str(output_dir),'--prefix',data['region_id'],'--selection',s['protein_selection'],
          '--obstacle-selection',s['obstacle_selection'],'--pbc',s['pbc'],'--geometry-mode','reference-region']
    for flag,key in [('spacing','spacing_A'),('probe-radius','probe_radius_A'),('region-margin','region_margin_A'),('profile-padding','profile_padding_A'),('max-grid-points','max_grid_points')]:argv.extend(['--'+flag,str(g[key])])
    argv.extend(['--axis='+','.join(map(str,g['axis'])),'--grid-phase='+','.join(map(str,g['grid_phase']))])
    used,check=resolve_region_radii(data,radii,allow_radii_mismatch=allow_radii_mismatch,base_dir=Path(path).resolve().parent)
    metadata={'schema_version':1,'region_id':data['region_id'],'label':data['label'],'review':data['review'],
              'definition_sha256':sha256(path),'definition':data,
              'meaning':'Named fixed regional observable; review status is declared provenance, not automatic biological validation',
              'radii':used.provenance(),'radii_check':check,'radius_set':used}
    return argv,metadata

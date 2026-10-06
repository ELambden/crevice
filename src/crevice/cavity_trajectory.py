"""Aligned, reference-matched cavity geometry and changing residue evidence.

This is the engine of ``crevice cavity-trajectory`` and
``crevice region-trajectory``. A fixed reference cavity
(:class:`CavityReference`, read from a binary DX map in the coordinates of
the first trajectory frame) defines one lattice, axis and anchor for every
frame. Each frame is rigidly aligned to the first frame, its cavity is
measured on that lattice (:func:`analyze_cavity_frame`) and either matched to
the reference by sampled overlap (``geometry_mode="rolling"``) or measured
inside the fixed reference neighbourhood (``"reference-region"``). The
matched cavity's volume, planar sections along the reference axis
(:func:`cavity_section_profile`), boundary residues and residue contacts are
recorded; :func:`summarize_cavity_trajectory` and
:func:`write_cavity_trajectory_bundle` turn them into statistics, tables and
figures.

Profiles describe measured planar sections, including one-sided/closed cavities.
They are not a through-path, permeation radius, or functional gate assignment.
Unresolved frames are missing values, never zeros.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
from pathlib import Path
import math

from .presentation import NON_CHANNEL_CAST_HEX, FigureText, annotate_option
from .radii import RadiusSet, active_radii, radii_fields, radii_option


@dataclass
class CavityReference:
    """Fixed reference cavity: occupied samples, lattice and axis shared by every frame.

    Build it with :meth:`from_dx`. All coordinates are in Å, in the frame of the
    first trajectory frame (the alignment reference).

    Attributes
    ----------
    points : numpy.ndarray, shape (n, 3)
        Centres of the occupied reference voxels.
    basis : numpy.ndarray, shape (3, 3)
        Right-handed orthonormal lattice axes as rows; the third row is the
        reference (profile) axis.
    anchor : numpy.ndarray, shape (3,)
        Lattice origin: the mean of ``points`` shifted by ``grid_phase``.
    spacing : float
        Lattice spacing, Å.
    axial_indices : numpy.ndarray of int
        Lattice planes along the axis covered by the profile (the reference extent
        plus ``profile_padding`` on each side); position = index x spacing.
    origin, grid, deltas
        The DX map as read: origin (Å), boolean occupancy grid and the three grid
        step vectors (Å).
    source_path : str or None
        Resolved path of the DX file.
    source_sha256 : str or None
        SHA-256 of the DX file, recorded in every output.
    grid_phase : tuple of float
        Lattice phase as fractions of a spacing in [0, 1) along each basis axis.
    """
    points: object
    basis: object
    anchor: object
    spacing: float
    axial_indices: object
    origin: object
    grid: object
    deltas: object
    source_path: str | None = None
    source_sha256: str | None = None
    grid_phase: tuple = (0.,0.,0.)

    @classmethod
    def from_dx(cls,path,*,spacing=None,axis=None,profile_padding=4.,grid_phase=(0.,0.,0.)):
        """Read a reference cavity from a binary DX map.

        Parameters
        ----------
        path : str or Path
            DX file whose nonzero values mark occupied voxels (for example a
            ``*_volume.dx`` written by ``crevice cast``). It must already be in the
            coordinates of the first trajectory frame.
        spacing : float, optional
            Lattice spacing, Å; default the length of the map's first grid step.
        axis : sequence of float, optional
            Reference axis direction (normalised here). Default: the map's own grid
            axes, with the third grid axis as the profile axis.
        profile_padding : float, default 4.0
            Extra axial range (Å) on each side of the reference cavity covered by
            the section profile, so a growing cavity is still measured.
        grid_phase : sequence of float, default (0, 0, 0)
            Lattice phase as fractions of a spacing in [0, 1).

        Returns
        -------
        CavityReference

        Raises
        ------
        ValueError
            For an empty map, a nonpositive spacing, negative padding, an invalid
            axis or phase, or a non-orthogonal grid.
        """
        import numpy as np
        from .residue_evidence import read_binary_dx
        from .geometry import basis_from_direction
        grid,origin,deltas=read_binary_dx(path)
        points=np.argwhere(grid)@deltas+origin
        if not len(points):raise ValueError('Reference volume is empty')
        h=float(np.linalg.norm(deltas[0])) if spacing is None else float(spacing)
        if not math.isfinite(h) or h<=0:raise ValueError('spacing must be positive')
        if not math.isfinite(profile_padding) or profile_padding<0:raise ValueError('profile_padding must be nonnegative')
        if axis is None:
            basis=deltas/np.linalg.norm(deltas,axis=1)[:,None]
        else:
            direction=np.asarray(axis,dtype=float)
            if direction.shape!=(3,) or not np.isfinite(direction).all() or np.linalg.norm(direction)==0:raise ValueError('axis must be a finite nonzero direction')
            direction/=np.linalg.norm(direction);u,v=basis_from_direction(tuple(direction));basis=np.array([u,v,direction])
        if not np.allclose(basis@basis.T,np.eye(3),atol=1e-8):raise ValueError('Reference grid must be orthogonal')
        if np.linalg.det(basis)<0:basis[0]*=-1
        phase=np.asarray(grid_phase,dtype=float)
        if phase.shape!=(3,) or not np.isfinite(phase).all() or np.any((phase<0)|(phase>=1)):
            raise ValueError("grid_phase must contain three finite fractions in [0,1)")
        anchor=points.mean(0)+(phase*h)@basis
        t=(points-anchor)@basis[2]
        axial=np.arange(math.floor((t.min()-profile_padding)/h),math.ceil((t.max()+profile_padding)/h)+1)
        import hashlib
        return cls(points,basis,anchor,h,axial,origin,grid,deltas,str(Path(path).resolve()),hashlib.sha256(Path(path).read_bytes()).hexdigest(),tuple(map(float,phase)))


def match_reference_component(cast,reference,*,minimum_overlap=.1,ambiguity_ratio=.8):
    """Symmetric sampled overlap, with explicit rejection of ambiguous identity.

    For each cast component, *precision* is the fraction of its samples within
    ``0.9 x max(reference spacing, cast spacing)`` of a reference sample and
    *recall* the fraction of reference samples within that distance of the
    component; the overlap score is their harmonic mean.

    Parameters
    ----------
    cast : VoidCast
        Cast of one aligned frame on the reference lattice.
    reference : CavityReference
        Fixed reference cavity.
    minimum_overlap : float, default 0.1
        Smallest acceptable best score, in (0, 1].
    ambiguity_ratio : float, default 0.8
        The match is refused as ambiguous when the second-best score is at least
        this fraction of the best, in (0, 1].

    Returns
    -------
    component : VoidComponent or None
        The matched component, or ``None`` when unresolved.
    report : dict
        Every candidate's ``overlap_score``, ``reference_recall``,
        ``candidate_precision`` and ``volume_A3``; the thresholds and distance
        tolerance; ``status`` (``matched``, ``unresolved_no_reference_match`` or
        ``unresolved_ambiguous_reference_match``) and, when matched, ``selected``.

    Raises
    ------
    ValueError
        If a threshold is outside (0, 1].
    """
    import numpy as np
    from scipy.spatial import cKDTree
    if not 0<minimum_overlap<=1 or not 0<ambiguity_ratio<=1:raise ValueError('Matching thresholds must lie in (0,1]')
    target=np.asarray(reference.points);tree=cKDTree(target);rows=[]
    tolerance=.9*max(reference.spacing,cast.spacing)
    for component in cast.components:
        points=np.asarray([p.position for p in component.points])
        precision=float(np.mean(tree.query(points)[0]<=tolerance))
        recall=float(np.mean(cKDTree(points).query(target)[0]<=tolerance)) if precision else 0.
        score=2*precision*recall/(precision+recall) if precision+recall else 0.
        rows.append({'component_id':component.id,'overlap_score':score,'reference_recall':recall,
                     'candidate_precision':precision,'volume_A3':component.volume})
    rows.sort(key=lambda r:(-r['overlap_score'],r['component_id']))
    report={'candidates':rows,'minimum_overlap':minimum_overlap,'ambiguity_ratio':ambiguity_ratio,
            'distance_tolerance_A':tolerance,'definition':'harmonic mean of sampled nearest-neighbor overlap fractions; matching the fixed aligned reference, not a biological identity proof'}
    if not rows or rows[0]['overlap_score']<minimum_overlap:
        return None,{**report,'status':'unresolved_no_reference_match'}
    if len(rows)>1 and rows[1]['overlap_score']>=ambiguity_ratio*rows[0]['overlap_score']:
        return None,{**report,'status':'unresolved_ambiguous_reference_match'}
    winner=next(c for c in cast.components if c.id==rows[0]['component_id'])
    return winner,{**report,'status':'matched','selected':rows[0]}


def cavity_section_profile(component,reference):
    """Measured lattice planes: area-equivalent diameters and local sphere size.

    For every reference axial plane, the component's samples in that plane give
    the total section area (samples x spacing²), the largest face-connected
    section area, their area-equivalent diameters ``2*sqrt(A/pi)``, the number
    of connected pieces and the diameter of the largest atom-clear sphere
    centred in the plane (twice the largest sample clearance). Empty planes of a
    successfully matched cavity have zero area. An unresolved cavity is handled
    separately as missing, never as a zero-width observation.

    Parameters
    ----------
    component : VoidComponent
        Matched component whose samples lie on the reference lattice.
    reference : CavityReference
        Fixed reference cavity (axis, anchor, spacing, axial planes).

    Returns
    -------
    dict
        Arrays over ``reference.axial_indices``: ``total_area_A2``,
        ``largest_area_A2``, ``equivalent_diameter_A``,
        ``largest_component_diameter_A``,
        ``maximum_atom_clear_sphere_diameter_A``, ``section_component_count``;
        and ``outside_profile_volume_A3`` (volume outside the profiled planes)
        and ``extent_A`` (first and last occupied plane position).

    Raises
    ------
    ValueError
        If the samples are not on the reference lattice.
    """
    import numpy as np
    from scipy import ndimage
    xyz=np.asarray([p.position for p in component.points]);h=reference.spacing
    local=(xyz-reference.anchor)@reference.basis.T/h
    index=np.rint(local).astype(int)
    if not np.allclose(local,index,atol=1e-5,rtol=0):raise ValueError('Profile requires the fixed reference lattice')
    clear=np.asarray([p.raw_clearance for p in component.points]);n=len(reference.axial_indices)
    total=np.zeros(n);largest=np.zeros(n);sphere=np.zeros(n);pieces=np.zeros(n,dtype=int)
    for j,k in enumerate(reference.axial_indices):
        selected=index[:,2]==k
        if not selected.any():continue
        xy=index[selected,:2];xy=xy-xy.min(0)
        plane=np.zeros(tuple(xy.max(0)+3),bool);plane[tuple((xy+1).T)]=True
        labels,count=ndimage.label(plane);sizes=np.bincount(labels.ravel())[1:]
        total[j]=selected.sum()*h*h;largest[j]=sizes.max()*h*h;pieces[j]=count
        sphere[j]=2*max(0.,float(clear[selected].max()))
    return {'total_area_A2':total,'largest_area_A2':largest,
            'equivalent_diameter_A':2*np.sqrt(total/np.pi),
            'largest_component_diameter_A':2*np.sqrt(largest/np.pi),
            'maximum_atom_clear_sphere_diameter_A':sphere,'section_component_count':pieces,
            'outside_profile_volume_A3':float(np.sum(~np.isin(index[:,2],reference.axial_indices))*h**3),
            'extent_A':[float(index[:,2].min()*h),float(index[:,2].max()*h)]}


def analyze_cavity_frame(frame,reference_frame,reference,*,geometry=None,alignment_residues=None,
                         minimum_overlap=.1,ambiguity_ratio=.8,lining_distance=1.5,contact_cutoff=4.5,
                         geometry_mode="rolling",region_margin=2.,obstacles=None,waters=None):
    """Align one frame and measure its cavity, boundary residues and contacts.

    1. Rigidly align ``frame`` onto ``reference_frame``
       (:func:`crevice.trajectory.align_frame_report`).
    2. Build the residue contact network of the aligned frame (atom-centre
       distance <= ``contact_cutoff``); contacts are defined even when the
       cavity is unresolved.
    3. Measure the cavity on the reference lattice: ``"rolling"`` runs
       :func:`crevice.rolling.rolling_probe_cast` and matches a component with
       :func:`match_reference_component`; ``"reference-region"`` measures
       probe-swept free volume inside the fixed reference neighbourhood
       (:func:`crevice.regional_volume.regional_probe_cast`), excluding crop
       faces from residue attribution.
    4. For a measured component: volume, :func:`cavity_section_profile` and
       :func:`crevice.residue_evidence.boundary_residue_evidence`. A measured
       zero volume gives zero areas, not missing values.
    5. Optionally, water membership (:mod:`crevice.water_membership`).

    Parameters
    ----------
    frame : StructureFrame
        Frame to analyse (original coordinates).
    reference_frame : StructureFrame
        Alignment reference (the first trajectory frame).
    reference : CavityReference
        Fixed reference cavity.
    geometry : dict, optional
        Keyword arguments for the cavity measurement, for example
        ``probe_radius`` (default 1.4 Å), ``outer_radius``, ``selection_mode``
        (default ``"dominant"``), ``enclosure_fraction``,
        ``min_component_volume`` (default 1 Å³) or ``max_grid_points``. The
        lattice keys ``grid_basis``, ``grid_anchor``, ``spacing`` and
        ``max_components`` are controlled by the reference and refused.
    alignment_residues : sequence of str, optional
        Residue IDs whose atoms define the fit; default all matching protein
        heavy atoms.
    minimum_overlap, ambiguity_ratio : float
        Matching thresholds (``"rolling"`` mode); see
        :func:`match_reference_component`.
    lining_distance : float, default 1.5
        Atom-surface distance (Å) from the cavity boundary within which a residue
        lines it.
    contact_cutoff : float, default 4.5
        Residue contact cutoff (Å, atom centres).
    geometry_mode : {"rolling", "reference-region"}, default "rolling"
        Cavity measurement mode (see above).
    region_margin : float, default 2.0
        Margin (Å) of the fixed reference neighbourhood (``"reference-region"``).
    obstacles : ObstacleSource, optional
        Extra obstacle atoms such as lipids or ligands
        (:class:`crevice.cavity_obstacles.ObstacleSource`); requires
        ``"reference-region"``.
    waters : WaterSource, optional
        When given, water membership is measured with this frame's grid and fit.

    Returns
    -------
    dict
        ``frame_index``, ``alignment`` (fit record), ``matching``, ``network``,
        ``obstacles``, ``status``, ``volume_A3``, ``profile``, ``boundary`` and,
        with ``waters``, ``water_membership``. Unresolved geometry leaves
        ``volume_A3``, ``profile`` and ``boundary`` as ``None``.

    Raises
    ------
    ValueError
        For an unknown ``geometry_mode``, a reference-controlled ``geometry`` key
        or obstacles outside ``"reference-region"`` mode.
    """
    import numpy as np
    from .trajectory import align_frame_report
    from .rolling import rolling_probe_cast
    from .volume_export import cast_grid
    from .residue_evidence import boundary_residue_evidence
    from .networks import build_residue_network
    aligned,fit=align_frame_report(frame,reference_frame,residues=alignment_residues)
    # Contacts remain defined even if cavity correspondence is unresolved.
    network=build_residue_network(aligned,cutoff=contact_cutoff,distance_metric='center')
    settings=dict(geometry or {})
    for key in ['grid_basis','grid_anchor','spacing','max_components']:
        if key in settings:raise ValueError(f'{key} is controlled by the trajectory reference')
    settings.setdefault('selection_mode','dominant');settings.setdefault('min_component_volume',1.)
    obstacle_atoms=None;obstacle_report=None
    if obstacles is not None:
        if geometry_mode!="reference-region":raise ValueError("Explicit obstacles currently require reference-region geometry")
        obstacle_atoms,obstacle_report=obstacles.prepare(frame,aligned,fit,reference,margin=region_margin,probe_radius=settings.get("probe_radius",1.4),max_grid_points=settings.get("max_grid_points",2_000_000))
    exclusion=None
    if geometry_mode=='reference-region':
        from .regional_volume import regional_probe_cast
        from scipy.spatial import cKDTree
        cast,component,regional=regional_probe_cast(aligned,reference,margin=region_margin,probe_radius=settings.get('probe_radius',1.4),max_grid_points=settings.get('max_grid_points',2_000_000),obstacle_atoms=obstacle_atoms)
        tree=cKDTree(reference.points)
        exclusion=lambda centers:tree.query(centers)[0]>=max(0.,region_margin-.9*reference.spacing)
        matching={'status':'regional_measured','regional_definition':regional,'selected':{}}
    elif geometry_mode=='rolling':
        cast=rolling_probe_cast(aligned,spacing=reference.spacing,grid_basis=reference.basis,
                                grid_anchor=reference.anchor,max_components=None,**settings)
        component,matching=match_reference_component(cast,reference,minimum_overlap=minimum_overlap,ambiguity_ratio=ambiguity_ratio)
    else:raise ValueError('geometry_mode must be rolling or reference-region')
    result={'frame_index':frame.frame_index,'alignment':fit,'matching':matching,'network':network,
            'obstacles':obstacle_report,'status':matching['status'],'volume_A3':None,'profile':None,'boundary':None}
    if component is not None and not component.points:
        ids=sorted({a.residue_key.label for a in aligned.atoms if not a.hetero and a.element!='H'})
        blank={'lining_residues':[],'boundary_area_A2':0.,'excluded_boundary_area_A2':0.,'assigned_boundary_area_A2':0.,'measurement_volume_A3':0.,
               'residues':[{'residue':key,'role':'other','boundary_area_A2':0.,'boundary_fraction':0.,'contacted_lining_boundary_fraction':0.,'nonlocal_lining_partners':[]} for key in ids]}
        empty={key:np.zeros(len(reference.axial_indices)) for key in ['total_area_A2','largest_area_A2','equivalent_diameter_A','largest_component_diameter_A','maximum_atom_clear_sphere_diameter_A','section_component_count']}
        empty.update(outside_profile_volume_A3=0.,extent_A=[None,None])
        result.update(volume_A3=0.,profile=empty,boundary=blank)
    elif component is not None:
        measured=replace(cast,components=(component,),points=component.points)
        grid,origin=cast_grid(measured);delta=reference.basis*reference.spacing
        boundary=boundary_residue_evidence(aligned,grid,origin,delta,lining_distance=lining_distance,contact_cutoff=contact_cutoff,boundary_exclusion=exclusion,obstacle_atoms=obstacle_atoms)
        result.update(volume_A3=component.volume,profile=cavity_section_profile(component,reference),boundary=boundary)
    if waters is not None:
        # Paired with this frame's measured grid and fit; the grid itself is not retained.
        from .water_membership import frame_water_membership
        result['water_membership']=frame_water_membership(frame,fit,reference,component,waters,geometry_mode=geometry_mode,
            region_margin=region_margin,probe_radius=settings.get('probe_radius',1.4),max_grid_points=settings.get('max_grid_points',2_000_000))
    return result


_WORKER=None

def _initialize_worker(frames,reference_frame,reference,settings,radii=None):
    # The radius set in effect in the parent is passed explicitly, so forked and
    # spawned workers measure with the same atomic radii (see crevice.radii.use_radii).
    global _WORKER
    _WORKER=frames,reference_frame,reference,settings
    if radii is not None:
        from .radii import _ACTIVE_RADII
        _ACTIVE_RADII.set(radii)


def _analyze_worker(index):
    frames,reference_frame,reference,settings=_WORKER
    return index,analyze_cavity_frame(frames[index],reference_frame,reference,**settings)


@radii_option
def analyze_cavity_trajectory(trajectory,reference,*,geometry=None,alignment_residues=None,
                               minimum_overlap=.1,ambiguity_ratio=.8,lining_distance=1.5,
                               contact_cutoff=4.5,workers=1,progress_callback=None,geometry_mode="rolling",region_margin=2.,obstacles=None,waters=None,
                               radii: RadiusSet | str | None = None):
    """Analyze every input frame; numerical/programming errors propagate.

    Each frame is processed by :func:`analyze_cavity_frame` against the first
    frame and the fixed reference. Residue contacts are accumulated in time
    order by :class:`crevice.residue_dynamics.ResidueDynamics`, seeded with the
    residues lining the reference cavity in the first frame. Parallel
    completion order does not alter frame ordering or statistics. Missing
    anatomical correspondence remains an explicit per-frame result.

    Parameters
    ----------
    trajectory : Trajectory
        Frames with finite, strictly increasing timestamps (ps); see
        :func:`crevice.trajectory.load_trajectory`.
    reference : CavityReference
        Fixed reference cavity in the first frame's coordinates.
    geometry : dict, optional
        Cavity measurement settings (see :func:`analyze_cavity_frame`).
    alignment_residues : sequence of str, optional
        Residues defining the rigid fit; default all protein heavy atoms.
    minimum_overlap : float, default 0.1
        Smallest acceptable reference overlap score (``"rolling"`` mode).
    ambiguity_ratio : float, default 0.8
        Second-best/best score ratio that makes a match ambiguous.
    lining_distance : float, default 1.5
        Boundary lining distance, Å.
    contact_cutoff : float, default 4.5
        Residue contact cutoff, Å.
    workers : int, default 1
        Worker processes (fork where available, else spawn). Results are
        identical for any number of workers.
    progress_callback : callable, optional
        Called as ``progress_callback(index, frame_record)`` as each frame
        finishes (completion order).
    geometry_mode : {"rolling", "reference-region"}, default "rolling"
        Cavity measurement mode.
    region_margin : float, default 2.0
        Reference neighbourhood margin, Å (``"reference-region"``).
    obstacles : ObstacleSource, optional
        Extra obstacle atoms (``"reference-region"`` only).
    waters : WaterSource, optional
        Measure water membership in every frame.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call (including every worker process) comes from that set,
        and results with a ``metadata`` dict record it as
        ``metadata["radii"]``; see :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict
        ``frames`` (one :func:`analyze_cavity_frame` record per frame, with
        ``time_ps``), ``reference``, ``reference_boundary`` (first-frame boundary
        evidence of the reference map), ``residue_dynamics``,
        ``residue_centroids`` (aligned residue centroids per frame, Å),
        ``settings``, ``frame_count``, ``aligned`` and ``time_ps``.

    Raises
    ------
    ValueError
        For an invalid ``workers``, missing or non-increasing timestamps, or
        obstacles outside ``"reference-region"`` mode.
    """
    from .residue_dynamics import ResidueDynamics
    from .residue_evidence import boundary_residue_evidence
    from .trajectory import align_frame_report
    import numpy as np
    if isinstance(workers,bool) or int(workers)!=workers or workers<1:raise ValueError('workers must be a positive integer')
    frames=trajectory.frames;first=frames[0];n=len(frames);times=trajectory.times()
    # Validate timestamps up front for all subsequent temporal statistics.
    if any(t is None or not math.isfinite(t) for t in times):raise ValueError('Cavity trajectory requires actual finite timestamps')
    if n>1 and np.any(np.diff(times)<=0):raise ValueError('Frame times must increase strictly')
    initial_obstacles=None
    if obstacles is not None:
        if geometry_mode!="reference-region":raise ValueError("Explicit obstacles currently require reference-region geometry")
        aligned,fit=align_frame_report(first,first,residues=alignment_residues)
        initial_obstacles,_=obstacles.prepare(first,aligned,fit,reference,margin=region_margin,probe_radius=(geometry or {}).get("probe_radius",1.4),max_grid_points=(geometry or {}).get("max_grid_points",2_000_000))
    initial=boundary_residue_evidence(first,reference.grid,reference.origin,reference.deltas,
                                      lining_distance=lining_distance,contact_cutoff=contact_cutoff,obstacle_atoms=initial_obstacles)
    dynamics=ResidueDynamics(lining_residues=initial['lining_residues'],min_occupancy=0)
    settings=dict(geometry=geometry,alignment_residues=alignment_residues,minimum_overlap=minimum_overlap,
                  ambiguity_ratio=ambiguity_ratio,lining_distance=lining_distance,contact_cutoff=contact_cutoff,geometry_mode=geometry_mode,region_margin=region_margin,obstacles=obstacles)
    if waters is not None:settings['waters']=waters
    records=[None]*n
    if workers==1:
        results=((i,analyze_cavity_frame(frames[i],first,reference,**settings)) for i in range(n))
        pool=None
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor,as_completed
        # Fork shares the read-only compact coordinate trajectory on POSIX.
        context=mp.get_context('fork' if 'fork' in mp.get_all_start_methods() else 'spawn')
        pool=ProcessPoolExecutor(max_workers=int(workers),mp_context=context,
                                initializer=_initialize_worker,initargs=(frames,first,reference,settings,active_radii()))
        futures=[pool.submit(_analyze_worker,i) for i in range(n)]
        results=(future.result() for future in as_completed(futures))
    try:
        for i,row in results:
            row['time_ps']=float(times[i]);records[i]=row
            if progress_callback:progress_callback(i,row)
    finally:
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
    # Sequential accumulation preserves physical time order for split-half motion.
    for i,row in enumerate(records):
        aligned,_=align_frame_report(frames[i],first,residues=alignment_residues)
        dynamics.add(aligned,row.pop('network'),frame_index=frames[i].frame_index,time=times[i])
    return {'frames':records,'reference':reference,'reference_boundary':initial,
            'residue_dynamics':dynamics.finish(aligned=True,include_contact_frames=True),
            'residue_centroids':np.asarray(dynamics.positions),
            'settings':{**settings,'obstacles':obstacles.metadata if obstacles is not None else None,**({'waters':waters.metadata} if waters is not None else {}),'workers':int(workers),'spacing_A':reference.spacing,**radii_fields()},
            'frame_count':n,'aligned':True,'time_ps':np.asarray(times,dtype=float)}


def describe_series(values,*,confidence=.95,block_length=None,replicates=2000,seed=20260913,regular=True):
    """Descriptive coverage plus autocorrelation-aware uncertainty for the mean.

    NaN values are missing frames. The record always reports the frame counts
    and coverage; with observations it adds the mean, sample SD, median,
    2.5/97.5% frame quantiles (a description of fluctuations, not an interval),
    extremes and first/second-half means. A ``confidence_interval`` for the mean
    (:func:`crevice.uncertainty.block_mean_confidence`) is attempted only for a
    complete, regularly sampled, non-constant series; otherwise its ``status``
    says why it was withheld (``no_observations``, ``unavailable_missing_frames``,
    ``unavailable_irregular_timestamps``, ``unavailable_observed_constant``).

    Parameters
    ----------
    values : array_like, shape (n,)
        One value per frame; NaN = missing.
    confidence : float, default 0.95
        Nominal level of the mean interval.
    block_length : int, optional
        Bootstrap block length in frames; default from the autocorrelation.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260913
        Bootstrap random seed.
    regular : bool, default True
        Whether the frames are evenly spaced in time; ``False`` withholds the
        interval.

    Returns
    -------
    dict
        ``total_frames``, ``observed_frames``, ``missing_frames``,
        ``coverage_fraction``, ``mean``, ``sample_sd``, ``median``,
        ``quantile_025``, ``quantile_975``, ``minimum``, ``maximum``,
        ``first_half_mean``, ``second_half_mean``, ``statistical_inefficiency``,
        ``effective_frames`` and ``confidence_interval``. Flattened to CSV by
        :func:`crevice.io.summary_statistics_rows`.

    Examples
    --------
    >>> from crevice.cavity_trajectory import describe_series
    >>> stats = describe_series([1.0, 2.0, float("nan"), 3.0])
    >>> stats["observed_frames"], stats["mean"], stats["confidence_interval"]["status"]
    (3, 2.0, 'unavailable_missing_frames')
    """
    import numpy as np
    from .uncertainty import block_mean_confidence,statistical_inefficiency
    values=np.asarray(values,dtype=float);observed=values[np.isfinite(values)];n=len(values)
    result={'total_frames':n,'observed_frames':len(observed),'missing_frames':n-len(observed),
            'coverage_fraction':len(observed)/n,'mean':None,'sample_sd':None,'median':None,
            'quantile_025':None,'quantile_975':None,'minimum':None,'maximum':None,
            'confidence_interval':{'status':'no_observations'},'statistical_inefficiency':None,'effective_frames':None}
    if len(observed):
        q=np.quantile(observed,[.025,.5,.975]);result.update(mean=float(observed.mean()),sample_sd=float(observed.std(ddof=1)) if len(observed)>1 else None,
           median=float(q[1]),quantile_025=float(q[0]),quantile_975=float(q[2]),minimum=float(observed.min()),maximum=float(observed.max()))
        result['first_half_mean']=float(np.nanmean(values[:n//2])) if np.isfinite(values[:n//2]).any() else None
        result['second_half_mean']=float(np.nanmean(values[n//2:])) if np.isfinite(values[n//2:]).any() else None
        if len(observed)==n and n>=2:
            g=statistical_inefficiency(values);result.update(statistical_inefficiency=g,effective_frames=n/g)
        if len(observed)==n and np.ptp(observed)==0:
            result['confidence_interval']={'status':'unavailable_observed_constant','reason':'No variation observed; a zero-width sampling interval is not justified for unobserved transitions'}
        elif len(observed)!=n:result['confidence_interval']={'status':'unavailable_missing_frames'}
        elif not regular:result['confidence_interval']={'status':'unavailable_irregular_timestamps'}
        else:result['confidence_interval']=block_mean_confidence(values[:,None],confidence=confidence,block_length=block_length,replicates=replicates,seed=seed)
    return result


def summarize_cavity_trajectory(analysis,*,confidence=.95,block_length=None,replicates=2000,seed=20260913):
    """All-frame arrays and statistics of a cavity-trajectory analysis.

    Volumes, section widths, residue boundary areas and fractions, lining
    occupancy and contact occupancy are gathered into frame x quantity arrays
    (NaN for unresolved frames) and described with :func:`describe_series`. The
    width profile uses one joint block-bootstrap family over all axial planes
    (pointwise and simultaneous bounds, clipped at zero; constant planes get no
    interval). Per residue, lining frequency is reported over observed frames and
    as all-frame bounds (unresolved frames counted as never/always lining).
    Pearson correlations of volume with each residue's boundary area and
    centroid distance are descriptive only.

    Parameters
    ----------
    analysis : dict
        Result of :func:`analyze_cavity_trajectory`.
    confidence : float, default 0.95
        Nominal level of the mean intervals.
    block_length : int, optional
        Bootstrap block length in frames; default from the autocorrelation.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260913
        Bootstrap random seed.

    Returns
    -------
    report : dict
        Frame counts, status counts, time range, settings, reference provenance,
        ``volume_A3``, ``profile`` (rows per axial position with observed
        2.5-97.5% frame range and mean bounds), ``residues``,
        ``residue_dynamics`` and ``limitations`` (written as
        ``PREFIX_cavity_statistics.json``).
    arrays : dict of str to numpy.ndarray
        ``time_ps``, ``volume_A3``, ``residue_ids``, ``boundary_area_A2``,
        ``boundary_fraction``, ``nonlocal_contacted_boundary_fraction``,
        ``nonlining_partner`` and the six section-profile arrays (written as
        ``PREFIX_cavity_frames.npz``).
    """
    import numpy as np
    from collections import Counter
    from copy import deepcopy
    from .uncertainty import block_mean_confidence
    frames=analysis['frames'];n=len(frames);reference=analysis['reference'];times=analysis['time_ps']
    dynamics=deepcopy(analysis['residue_dynamics'])
    regular=n<3 or np.allclose(np.diff(times),np.diff(times)[0],rtol=1e-7,atol=1e-7)
    options=dict(confidence=confidence,block_length=block_length,replicates=replicates,seed=seed,regular=regular)
    volume=np.asarray([r['volume_A3'] if r['volume_A3'] is not None else np.nan for r in frames])
    ids=sorted(r['residue'] for r in analysis['reference_boundary']['residues']);lookup={k:i for i,k in enumerate(ids)}
    area=np.full((n,len(ids)),np.nan);support=area.copy();fraction=area.copy();partner=area.copy();pairs={}
    keys=['total_area_A2','largest_area_A2','equivalent_diameter_A','largest_component_diameter_A','maximum_atom_clear_sphere_diameter_A','section_component_count']
    profiles={key:np.full((n,len(reference.axial_indices)),np.nan) for key in keys}
    for i,row in enumerate(frames):
        if row['boundary'] is None:continue
        for key in keys:profiles[key][i]=row['profile'][key]
        for residue in row['boundary']['residues']:
            j=lookup[residue['residue']];area[i,j]=residue['boundary_area_A2'];fraction[i,j]=residue['boundary_fraction']
            support[i,j]=residue['contacted_lining_boundary_fraction'];partner[i,j]=float(residue['role']=='nonlocal_lining_partner')
            if residue['role']=='nonlocal_lining_partner':
                for contact in residue['nonlocal_lining_partners']:
                    key=tuple(sorted((residue['residue'],contact['residue'])))
                    pairs.setdefault(key,set()).add(i)
    valid=np.isfinite(volume);count=int(valid.sum())
    # One joint family for widths, keeping missing frames/planes intact.
    width=profiles['largest_component_diameter_A']
    ci=block_mean_confidence(width,confidence=confidence,block_length=block_length,replicates=replicates,seed=seed) if regular else {'status':'unavailable_irregular_timestamps'}
    constant=np.isfinite(width).all(axis=0) & (np.nan_to_num(width).max(axis=0)==np.nan_to_num(width).min(axis=0))
    ci['constant_columns_without_sampling_interval']=np.flatnonzero(constant).tolist()
    for key in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
        if key in ci:
            ci[key]=[None if constant[j] else max(0.,v) if v is not None else None for j,v in enumerate(ci[key])]
    profile_rows=[]
    for j,k in enumerate(reference.axial_indices):
        x=width[:,j];y=x[np.isfinite(x)]
        row={'position_A':float(k*reference.spacing),'observed_frames':len(y),'coverage_fraction':len(y)/n,
             'occupied_fraction_of_observed':float(np.mean(y>0)) if len(y) else None,
             'mean_diameter_A':float(y.mean()) if len(y) else None,'median_diameter_A':float(np.median(y)) if len(y) else None,
             'fluctuation_lower_A':float(np.quantile(y,.025)) if len(y) else None,'fluctuation_upper_A':float(np.quantile(y,.975)) if len(y) else None}
        for out,key in [('mean_lower_A','pointwise_lower'),('mean_upper_A','pointwise_upper'),('simultaneous_mean_lower_A','simultaneous_lower'),('simultaneous_mean_upper_A','simultaneous_upper')]:
            row[out]=ci.get(key,[None]*len(reference.axial_indices))[j]
        profile_rows.append(row)
    residues=[]
    for j,key in enumerate(ids):
        lining=np.where(np.isfinite(area[:,j]),(area[:,j]>0).astype(float),np.nan)
        direct_count=int(np.nansum(lining));partner_count=int(np.nansum(partner[:,j]))
        residues.append({'residue':key,'observed_geometry_frames':count,'total_frames':n,
            'lining_frames':direct_count,'lining_frequency_observed':direct_count/count if count else None,
            'lining_frequency_all_frames_lower':direct_count/n,'lining_frequency_all_frames_upper':(direct_count+n-count)/n,
            'nonlining_partner_frames':partner_count,'partner_frequency_observed':partner_count/count if count else None,
            'boundary_area_A2':describe_series(area[:,j],**options),'boundary_fraction':describe_series(fraction[:,j],**options),
            'lining_occupancy':describe_series(lining,**options),'nonlocal_contacted_boundary_fraction':describe_series(support[:,j],**options)})
    # Associate motion with geometry using the same aligned frames; descriptive
    # correlation only, conditional on matched frames, no p-values/causal claim.
    for row in dynamics['contacts']:
        key=tuple(sorted((row['source'],row['target'])));dynamic=pairs.get(key,set())
        row['dynamic_nonlining_boundary_partner_frames']=len(dynamic)
        row['dynamic_boundary_eligible_frames']=count
        row['dynamic_partner_frequency_observed']=len(dynamic)/count if count else None
        occurrence=np.zeros(n);occurrence[row.pop('contact_frame_rows')]=1.
        row['occupancy_statistics']=describe_series(occurrence,**options)
        for key in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
            interval=row['occupancy_statistics']['confidence_interval']
            if key in interval:interval[key]=[float(np.clip(v,0,1)) if v is not None else None for v in interval[key]]
    centroids=analysis.get('residue_centroids')
    def correlation(a,b):
        keep=np.isfinite(a)&np.isfinite(b)
        if keep.sum()<3 or np.std(a[keep])<1e-10 or np.std(b[keep])<1e-10:return None
        return float(np.corrcoef(a[keep],b[keep])[0,1])
    for j,row in enumerate(residues):
        row['volume_boundary_area_correlation']=correlation(volume,area[:,j])
        row['volume_centroid_distance_correlation']=correlation(volume,np.linalg.norm(centroids[:,j]-reference.anchor,axis=1)) if centroids is not None else None
        row['correlation_interpretation']='Descriptive Pearson association on matched frames; no uncertainty or causal interpretation; area and volume are geometrically related'
    report={'frame_count':n,'matched_frames':count,'unresolved_frames':n-count,'status_counts':dict(Counter(r['status'] for r in frames)),
            'time_range_ps':[float(times[0]),float(times[-1])],'regular_timestamps':bool(regular),'alignment_applied':True,
            'settings':analysis['settings'],'reference':{'source_path':reference.source_path,'source_sha256':reference.source_sha256,
              'grid_phase_fraction':list(reference.grid_phase),'grid_basis':reference.basis.tolist(),'grid_anchor_A':reference.anchor.tolist(),'spacing_A':reference.spacing,
              'coordinate_contract':'Reference map must already be in the first trajectory frame coordinates; geometric matching is not anatomical validation'},
            'volume_A3':describe_series(volume,**options),
            'profile':{'rows':profile_rows,'confidence_interval':ci,'axis_direction':reference.basis[2].tolist(),'axis_origin':reference.anchor.tolist(),
              'definition':'2*sqrt(A/pi), where A is the largest face-connected measured planar section area; zero means no section of a successfully matched cavity',
              'fluctuation_band':'2.5–97.5% empirical frame quantiles, conditional on matched geometry; not a confidence interval',
              'mean_band':'autocorrelation-aware joint batch bootstrap; may be unavailable with missing frames or insufficient independent batches',
              'limitation':'area-equivalent sectional width, not an inscribed permeation diameter or a connected through-path'},
            'residues':residues,'residue_dynamics':dynamics,
            'limitations':['Reference matching is geometric and cannot establish biological cavity identity.','Statistical intervals do not include probe/grid, force-field or anatomical errors.','All selected frames are attempted; unresolved geometry is missing, never zero.','Single-trajectory stationarity and convergence are not established.','Residue boundary/contact/motion evidence describes associations, not functional or causal control.']}
    arrays={'time_ps':times,'volume_A3':volume,'residue_ids':np.asarray(ids),'boundary_area_A2':area,'boundary_fraction':fraction,
            'nonlocal_contacted_boundary_fraction':support,'nonlining_partner':partner,**profiles}
    return report,arrays


@annotate_option
def write_cavity_trajectory_bundle(analysis,output_dir,*,prefix='crevice',confidence=.95,
                                   block_length=None,replicates=2000,seed=20260913,dpi=240,annotate=None):
    """Write the statistics, tables and figures of ``crevice cavity-trajectory``.

    ``analysis`` (from :func:`analyze_cavity_trajectory`) holds, for every
    aligned frame, the measured cavity matched to the fixed reference, its
    boundary residues and its sectional widths along the fixed reference axis.
    :func:`summarize_cavity_trajectory` turns these into all-frame arrays and
    statistics: observed frame ranges (2.5-97.5% quantiles, a description of
    fluctuations) and, separately, pointwise confidence intervals for means from
    a block bootstrap over frames, withheld when sampling safeguards fail.
    Unresolved frames are missing values, never zeros. Every frame is kept.

    Parameters
    ----------
    analysis : dict
        Result of :func:`analyze_cavity_trajectory`.
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    confidence : float, default 0.95
        Nominal level of the mean intervals; 0 gives descriptive statistics only.
    block_length : int, optional
        Bootstrap block length in frames; default from the autocorrelation.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260913
        Bootstrap random seed.
    dpi : int, default 240
        Resolution of every figure.
    annotate : bool, optional
        Draw panel titles and figure titles (frame counts, mean-interval status
        and interpretation reminders) in the figures. ``None`` (default) inherits
        the surrounding setting, which is off unless enabled with
        :func:`crevice.presentation.figure_annotations` or ``--annotate``.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_cavity_statistics.json : JSON
        Settings, volume, width-profile, residue and contact statistics, limitations.
    PREFIX_cavity_statistics_summary.csv : table
        The JSON's interval statistics flattened, one row per quantity (volume,
        each residue's boundary area/fraction/occupancy, each contact's
        occupancy) and per correlation; columns as in
        :func:`crevice.io.summary_statistics_rows`.
    PREFIX_cavity_frames.csv : table
        One row per frame: ``frame_index``, ``time_ps``, ``status``,
        ``volume_A3``, ``alignment_rmsd_A``, matching scores, axial extent
        (``axial_start_A``/``axial_end_A``) and boundary counts.
    PREFIX_cavity_section_frames.csv : table
        One row per frame with measured geometry and reference-axis position:
        ``frame_index``, ``time_ps``, ``position_A``, ``total_area_A2``,
        ``largest_area_A2``, ``equivalent_diameter_A``,
        ``largest_component_diameter_A``,
        ``maximum_atom_clear_sphere_diameter_A``, ``section_component_count``
        (empty = unresolved section).
    PREFIX_cavity_residue_frames.csv : table
        One row per frame and residue on that frame's measured boundary:
        ``frame_index``, ``time_ps``, ``residue``, ``boundary_area_A2``,
        ``boundary_fraction``, ``nonlocal_contacted_boundary_fraction`` and
        ``nonlining_partner`` (1 = nonlocal lining partner).
    PREFIX_cavity_frames.npz : arrays
        The same all-frame arrays (the array format read by ``region-compare``).
    PREFIX_cavity_frame_diagnostics.json : JSON
        Full per-frame alignment (rotation, translation) and matching
        diagnostics, kept as provenance.
    PREFIX_cavity_width_profile.csv, PREFIX_cavity_residues.csv, PREFIX_cavity_contacts.csv : tables
        Width statistics per axial position; residue boundary statistics; contacts.
    PREFIX_cavity_volume_timeseries.png : figure
        Three panels against time (ns): "Cavity volume (Å³)" with the observed
        mean and mean interval, "Boundary residues" and "CA fit RMSD (Å)".
    PREFIX_cavity_diameter_profile.png, PREFIX_cavity_radius_profile.png : figures
        Area-equivalent diameter/radius (Å) of the largest connected section
        against "Position along reference axis (Å)", with a coverage panel ("Fraction").
    PREFIX_cavity_residues.png : figure
        Top 20 residues: "Boundary frequency | matched geometry" and "Mean
        assigned boundary area (Å²)".
    PREFIX_cavity_width_heatmap.png : figure
        Area-equivalent diameter by reference-axis position (Å) and time (ns).
    PREFIX_cavity_section_statistics.csv, PREFIX_cavity_free_sphere_statistics.csv, PREFIX_cavity_profile_statistics.json : tables
        From :func:`crevice.cavity_reports.write_additional_cavity_reports`.
    PREFIX_cavity_section_radius_statistics.png, PREFIX_cavity_free_sphere_radius_statistics.png, PREFIX_cavity_dynamic_partners.png, PREFIX_cavity_residue_heatmap.png : figures
        See :func:`crevice.cavity_reports.write_additional_cavity_reports`.
    PREFIX_cavity_partner_statistics.csv, PREFIX_cavity_contact_statistics.csv : tables
        Dynamic partner and contact statistics.
    PREFIX_cavity_manifest.json : JSON
        File index, frame and matched-frame counts, interpretation.

    Every PNG embeds a ``Description`` (read it with
    :func:`crevice.presentation.figure_description`) that records the frame
    counts, interval status and every omitted title.

    Colours and representation
    --------------------------
    * Volume timeseries: teal (``#178c91``) volume line, brown (``#a36a22``)
      observed-mean line, orange (``#d98c33``) mean-interval band at 25%;
      orange boundary-residue line; blue-grey (``#547c9b``) RMSD line.
    * Width profiles (a non-channel cast): violet (``#a855f7``, the
      non-channel cast colour of :data:`crevice.presentation.NON_CHANNEL_CAST_HEX`)
      band (20%) = observed 2.5-97.5% frame range, violet line = mean over
      matched frames, orange band (40%) = pointwise mean
      interval; coverage panel blue-grey = matched coverage, magenta
      (``#a34477``) = section present given matched. Legends name these series.
    * Residue bars: orange boundary frequency, blue-grey mean boundary area.
    * Width heatmap: ``viridis`` colour scale with a colour bar
      "Area-equivalent diameter (Å)"; light grey = unresolved geometry.
    """
    import csv,json
    import numpy as np
    from .figures import _pyplot,_save_figure
    root=Path(output_dir);root.mkdir(parents=True,exist_ok=True)
    if not prefix or Path(prefix).name!=prefix:raise ValueError('prefix must be a simple filename stem')
    report,arrays=summarize_cavity_trajectory(analysis,confidence=confidence,block_length=block_length,replicates=replicates,seed=seed)
    files={}
    def path(key,suffix):
        target=root/(prefix+suffix);files[key]=str(target);return target
    path('statistics_json','_cavity_statistics.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    from .io import write_summary_statistics_csv
    write_summary_statistics_csv(report,path('statistics_summary_csv','_cavity_statistics_summary.csv'))
    np.savez_compressed(path('all_frame_arrays_npz','_cavity_frames.npz'),**arrays)
    frame_rows=[]
    for row in analysis['frames']:
        selected=row['matching'].get('selected',{})
        frame_rows.append({'frame_index':row['frame_index'],'time_ps':row['time_ps'],'status':row['status'],
          'volume_A3':row['volume_A3'],'alignment_rmsd_A':row['alignment']['rmsd_after_A'],
          'matching_score':selected.get('overlap_score'),'reference_recall':selected.get('reference_recall'),
          'candidate_precision':selected.get('candidate_precision'),
          'profile_outside_volume_A3':row['profile']['outside_profile_volume_A3'] if row['profile'] else None,
          'axial_start_A':row['profile']['extent_A'][0] if row['profile'] else None,
          'axial_end_A':row['profile']['extent_A'][1] if row['profile'] else None,
          'lining_residues':len(row['boundary']['lining_residues']) if row['boundary'] else None,
          'excluded_crop_boundary_A2':row['boundary'].get('excluded_boundary_area_A2',0.) if row['boundary'] else None,
          'region_volume_components':row['matching'].get('regional_definition',{}).get('sampled_volume_components')})
    def table(target,rows):
        if not rows:target.write_text('');return
        with target.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    table(path('frames_csv','_cavity_frames.csv'),frame_rows)
    def cell(v):
        v=float(v);return v if np.isfinite(v) else ''
    frame_ids=[row['frame_index'] for row in analysis['frames']];frame_times=[row['time_ps'] for row in analysis['frames']]
    positions=[r['position_A'] for r in report['profile']['rows']]
    section_keys=[k for k in ['total_area_A2','largest_area_A2','equivalent_diameter_A','largest_component_diameter_A',
                              'maximum_atom_clear_sphere_diameter_A','section_component_count'] if k in arrays]
    section_rows=[{'frame_index':frame_ids[i],'time_ps':frame_times[i],'position_A':positions[j],**{k:cell(arrays[k][i,j]) for k in section_keys}}
                  for i in range(len(frame_ids)) if any(np.isfinite(arrays[k][i]).any() for k in section_keys)
                  for j in range(len(positions))]
    table(path('section_frames_csv','_cavity_section_frames.csv'),section_rows)
    residue_keys=[k for k in ['boundary_area_A2','boundary_fraction','nonlocal_contacted_boundary_fraction','nonlining_partner'] if k in arrays]
    residue_ids=[str(v) for v in arrays['residue_ids']]
    residue_frame_rows=[{'frame_index':frame_ids[i],'time_ps':frame_times[i],'residue':residue_ids[j],**{k:cell(arrays[k][i,j]) for k in residue_keys}}
                        for i in range(len(frame_ids)) for j in range(len(residue_ids)) if np.isfinite(arrays['boundary_area_A2'][i,j])]
    table(path('residue_frames_csv','_cavity_residue_frames.csv'),residue_frame_rows)
    path('frame_diagnostics_json','_cavity_frame_diagnostics.json').write_text(json.dumps([
        {k:v for k,v in row.items() if k not in {'boundary','profile','water_membership'}} for row in analysis['frames']],indent=2,allow_nan=False)+'\n')
    table(path('width_profile_csv','_cavity_width_profile.csv'),report['profile']['rows'])
    residue_rows=[]
    for row in report['residues']:
        area=row['boundary_area_A2'];ci=area['confidence_interval'];freq=row['lining_occupancy']
        residue_rows.append({'residue':row['residue'],'lining_frames':row['lining_frames'],'observed_geometry_frames':row['observed_geometry_frames'],
          'total_frames':row['total_frames'],'lining_frequency_observed':row['lining_frequency_observed'],
          'lining_frequency_all_frames_lower':row['lining_frequency_all_frames_lower'],'lining_frequency_all_frames_upper':row['lining_frequency_all_frames_upper'],
          'partner_frequency_observed':row['partner_frequency_observed'],'mean_boundary_area_A2':area['mean'],
          'sd_boundary_area_A2':area['sample_sd'],'boundary_mean_lower_A2':ci.get('pointwise_lower',[None])[0],
          'boundary_mean_upper_A2':ci.get('pointwise_upper',[None])[0],'mean_CI_status':ci['status'],
          'lining_occupancy_CI_status':freq['confidence_interval']['status']})
    table(path('residue_statistics_csv','_cavity_residues.csv'),residue_rows)
    table(path('contacts_csv','_cavity_contacts.csv'),[
        {k:v for k,v in row.items() if not isinstance(v,(dict,list))} for row in report['residue_dynamics']['contacts']])
    plt=_pyplot();times=np.asarray(arrays['time_ps'])/1000;valid=np.isfinite(arrays['volume_A3'])
    text=FigureText(f"Cavity volume, boundary-residue count and CA fit RMSD against time for all {len(times)} attempted frames "
                    f"({int(valid.sum())} matched); unresolved frames are gaps")
    fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True,layout='constrained')
    axes[0].plot(times,arrays['volume_A3'],color='#178c91',lw=.9);axes[0].set_ylabel('Cavity volume (Å³)')
    stat=report['volume_A3'];bounds=stat['confidence_interval']
    if stat['mean'] is not None:axes[0].axhline(stat['mean'],color='#a36a22',lw=1.4,label='Observed mean')
    if bounds.get('pointwise_lower',[None])[0] is not None:
        axes[0].axhspan(bounds['pointwise_lower'][0],bounds['pointwise_upper'][0],color='#d98c33',alpha=.25,label=f'{confidence:.0%} mean interval')
    text.title(axes[0],f"All {len(times)} frames attempted · {int(valid.sum())} matched · mean interval: {bounds['status']}",fontsize=10)
    axes[1].plot(times,[r['lining_residues'] if r['lining_residues'] is not None else np.nan for r in frame_rows],color='#d98c33',lw=.9);axes[1].set_ylabel('Boundary residues')
    axes[2].plot(times,[r['alignment_rmsd_A'] for r in frame_rows],color='#547c9b',lw=.9);axes[2].set_ylabel('CA fit RMSD (Å)');axes[2].set_xlabel('Time (ns)')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    _save_figure(fig,path('volume_timeseries_png','_cavity_volume_timeseries.png'),dpi=dpi,text=text)
    rows=report['profile']['rows'];x=np.array([r['position_A'] for r in rows])
    def values(key):return np.array([r[key] if r[key] is not None else np.nan for r in rows])
    for quantity,scale in [('diameter',1.),('radius',.5)]:
        text=FigureText(f"Area-equivalent {quantity} of the largest connected planar cavity section along the fixed reference axis: "
                        "observed frame range, mean and pointwise mean interval, with frame coverage below")
        fig,axes=plt.subplots(2,1,figsize=(10,6.5),sharex=True,height_ratios=[4,1],layout='constrained')
        ax=axes[0];ax.fill_between(x,scale*values('fluctuation_lower_A'),scale*values('fluctuation_upper_A'),color=NON_CHANNEL_CAST_HEX,alpha=.20,label='Observed 2.5–97.5% frame range')
        ax.plot(x,scale*values('mean_diameter_A'),color=NON_CHANNEL_CAST_HEX,lw=1.7,label='Mean over matched frames')
        lower=values('mean_lower_A');upper=values('mean_upper_A')
        if np.isfinite(lower).any():ax.fill_between(x,scale*lower,scale*upper,color='#d98c33',alpha=.4,label=f'{confidence:.0%} pointwise mean interval')
        ax.set_ylabel(f'Area-equivalent {quantity} (Å)');ax.set_ylim(bottom=0);ax.legend(frameon=False,fontsize=9)
        text.title(ax,'Largest connected planar cavity section · aligned fixed coordinate\nNot an inscribed permeation radius or through-channel assignment',fontsize=11)
        axes[1].plot(x,values('coverage_fraction'),color='#547c9b',label='Matched coverage')
        axes[1].plot(x,values('occupied_fraction_of_observed'),color='#a34477',label='Section present | matched')
        axes[1].set_ylim(-.02,1.02);axes[1].set_ylabel('Fraction');axes[1].set_xlabel('Position along reference axis (Å)');axes[1].legend(frameon=False,fontsize=8,ncol=2)
        for axis in axes:axis.spines[['top','right']].set_visible(False)
        text.suptitle(fig,'Mean-band status: '+report['profile']['confidence_interval']['status'],fontsize=10)
        _save_figure(fig,path(quantity+'_profile_png','_cavity_'+quantity+'_profile.png'),dpi=dpi,text=text)
    ranked=sorted(residue_rows,key=lambda r:(-(r['lining_frequency_observed'] or 0),-(r['mean_boundary_area_A2'] or 0),r['residue']))[:20]
    text=FigureText(f"Boundary frequency and mean assigned boundary area of the {len(ranked)} residues most often on the cavity boundary")
    fig,axes=plt.subplots(1,2,figsize=(12,7),sharey=True,layout='constrained');labels=[r['residue'] for r in ranked]
    axes[0].barh(labels,[r['lining_frequency_observed'] or 0 for r in ranked],color='#d98c33');axes[0].invert_yaxis();axes[0].set_xlabel('Boundary frequency | matched geometry');axes[0].set_xlim(0,1.02)
    axes[1].barh(labels,[r['mean_boundary_area_A2'] or 0 for r in ranked],color='#547c9b');axes[1].set_xlabel('Mean assigned boundary area (Å²)')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    text.suptitle(fig,'Changing boundary involvement · geometric evidence, not functional control')
    _save_figure(fig,path('residue_statistics_png','_cavity_residues.png'),dpi=dpi,text=text)
    # Time-position map distinguishes fluctuating shape from missing geometry.
    text=FigureText("Area-equivalent diameter of the largest connected section at each reference-axis position in every frame")
    fig,ax=plt.subplots(figsize=(11,5),layout='constrained')
    im=ax.pcolormesh(x,times,arrays['largest_component_diameter_A'],shading='nearest',cmap='viridis');ax.set_facecolor('#dddddd')
    fig.colorbar(im,ax=ax,label='Area-equivalent diameter (Å)');ax.set_xlabel('Reference-axis position (Å)');ax.set_ylabel('Time (ns)');text.title(ax,'All-frame cavity width · grey = unresolved geometry')
    _save_figure(fig,path('width_heatmap_png','_cavity_width_heatmap.png'),dpi=dpi,text=text)
    from .cavity_reports import write_additional_cavity_reports
    files.update(write_additional_cavity_reports(report,arrays,root,prefix=prefix,confidence=confidence,
                 block_length=block_length,replicates=replicates,seed=seed,dpi=dpi))
    path('manifest_json','_cavity_manifest.json').write_text(json.dumps({'files':files,'frame_count':len(times),'matched_frames':int(valid.sum()),'interpretation':report['limitations'],**radii_fields()},indent=2)+'\n')
    return files

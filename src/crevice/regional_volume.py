"""Probe-swept volume in an explicit, fixed reference neighbourhood.

This bounded geometric observable is defined even when cavities split/merge.
It does not assert ligand accessibility or replace global cavity detection.

Trajectory cavity analysis can either track "the" cavity from frame to frame
(which becomes ill-defined when a cavity splits or merges) or measure the
free volume inside a fixed region defined once from a reference cavity map.
This module implements the second option (``geometry_mode="reference-region"``
in :func:`crevice.cavity_trajectory.analyze_cavity_frame`).

Definition
    The *domain* is every grid point within ``margin`` Å of a reference
    cavity sample. Probe centres are grid points within ``margin +
    probe_radius`` of the reference whose atom-surface clearance is at least
    ``probe_radius``; they are grouped into connected components with
    segment-checked face edges. Components that contain a centre within
    ``0.9 * spacing`` of a reference sample are retained. The measured
    volume is the part of the domain within ``probe_radius`` of a retained
    centre and outside every atom, i.e. the probe-swept space clipped to
    the fixed neighbourhood.

Main entry points: :func:`regional_probe_cast` and :func:`regional_grid`.
"""
import math
from .radii import RadiusSet, radii_option


def regional_grid(reference,margin=2.,probe_radius=1.4,max_grid_points=2_000_000):
    """Grid covering a reference neighbourhood, cached on the reference.

    The grid uses the reference's basis, spacing and anchor, and spans the
    reference samples padded by ``margin + 2 * probe_radius + 2 * spacing``.
    The result is cached in ``reference._regional_cache`` for each parameter
    combination, so every frame reuses one grid.

    Parameters
    ----------
    reference : CavityReference
        Reference region (``points``, ``basis``, ``anchor``, ``spacing``).
    margin : float, default 2.0
        Neighbourhood margin in Å (finite, non-negative).
    probe_radius : float, default 1.4
        Probe radius in Å (finite, positive).
    max_grid_points : int, default 2000000

    Returns
    -------
    tuple
        ``(index, xyz, shape, distance, tree, low)``: integer grid indices,
        world coordinates of every node (Å), grid shape, distance (Å) from
        each node to the nearest reference sample, a KD-tree of the reference
        samples and the lowest grid index.

    Raises
    ------
    ValueError
        For an invalid margin or probe radius, or a grid over
        ``max_grid_points``.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    if not math.isfinite(margin) or margin<0:raise ValueError('region margin must be finite and nonnegative')
    if not math.isfinite(probe_radius) or probe_radius<=0:raise ValueError('region probe radius must be finite and positive')
    h=reference.spacing;basis=reference.basis;anchor=reference.anchor
    cache=getattr(reference,'_regional_cache',{})
    key=(margin,probe_radius,max_grid_points,h,tuple(anchor),tuple(basis.ravel()))
    if key not in cache:
        local=(reference.points-anchor)@basis.T
        padding=margin+2*probe_radius+2*h
        low=np.floor((local.min(0)-padding)/h).astype(int);high=np.ceil((local.max(0)+padding)/h).astype(int)
        shape=tuple(high-low+1)
        if math.prod(shape)>max_grid_points:raise ValueError('Reference neighborhood exceeds max_grid_points')
        index=np.indices(shape).reshape(3,-1).T+low
        xyz=index*h@basis+anchor;tree=cKDTree(reference.points)
        distance=tree.query(xyz)[0].reshape(shape)
        cache[key]=(index,xyz,shape,distance,tree,low);reference._regional_cache=cache
    return cache[key]


@radii_option
def regional_probe_cast(frame,reference,*,margin=2.,probe_radius=1.4,max_grid_points=2_000_000,obstacle_atoms=None, radii: RadiusSet | str | None = None):
    """Measure probe-swept free volume inside a fixed reference neighbourhood.

    Parameters
    ----------
    frame : StructureFrame
        Frame already aligned to the reference.
    reference : CavityReference
        Fixed reference region.
    margin : float, default 2.0
        Neighbourhood margin in Å around the reference samples.
    probe_radius : float, default 1.4
        Probe radius in Å.
    max_grid_points : int, default 2000000
    obstacle_atoms : sequence of Atom, optional
        Atoms that bound the space. Defaults to the frame's non-HETATM heavy
        atoms; see :class:`crevice.cavity_obstacles.ObstacleSource` for adding
        lipids, ligands or water.
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
    cast : VoidCast
        ``mode="rolling"``, ``min_radius=0``, with the grid basis and origin and
        the report in ``metadata``.
    component : VoidComponent
        Kind ``"fixed_reference_neighborhood"``; volume ``n * spacing**3``
        (Å³). ``t`` of each point is its coordinate along the reference's
        third basis vector from the anchor. Empty (volume 0) when nothing is
        retained.
    report : dict
        Definition, margin, probe radius, component counts, grid size, volume
        and limitations. A zero volume means no retained probe-swept samples in
        this region, not evidence of a closed transport pathway.

    Raises
    ------
    ValueError
        If there are no obstacle atoms, or as for :func:`regional_grid`.
    """
    import numpy as np
    from scipy import ndimage
    from .rolling import _SphereQueries,_connected
    from .radii import atom_vdw_radius
    from .models import ChannelPoint,VoidComponent,VoidCast
    h=reference.spacing;basis=reference.basis;anchor=reference.anchor
    index,xyz,shape,distance,tree,low=regional_grid(reference,margin,probe_radius,max_grid_points)
    atoms=frame.selected_atoms(include_hydrogen=False,include_hetero=False) if obstacle_atoms is None else tuple(obstacle_atoms)
    if not atoms:raise ValueError("Regional geometry needs physical obstacles")
    coords=np.asarray([a.coord for a in atoms])@basis.T
    queries=_SphereQueries(coords,np.asarray([atom_vdw_radius(a) for a in atoms]))
    clearance,nearest=queries.points(xyz@basis.T);clearance=clearance.reshape(shape)
    domain=distance<=margin+1e-9
    centers=(distance<=margin+probe_radius+1e-9)&(clearance>=probe_radius)
    origin=anchor@basis.T+low*h
    labels=_connected(centers,clearance,origin,h,probe_radius,queries)
    seeds=centers&(distance<=.9*h)
    retained=np.unique(labels[seeds]);retained=retained[retained>0]
    selected=np.isin(labels,retained)&(labels>0)
    if selected.any():
        dist=ndimage.distance_transform_edt(~selected,sampling=h)
        filled=domain&(clearance>=0)&(dist<=probe_radius+1e-9)
    else:filled=np.zeros(shape,bool)
    flat=filled.ravel();raw=clearance.ravel()[flat];near=nearest[flat];positions=xyz[flat]
    points=tuple(ChannelPoint(i,tuple(map(float,p)),float(r),float(r),float((p-anchor)@basis[2]),atoms[int(j)].serial,atoms[int(j)].residue_key.label)
                 for i,(p,r,j) in enumerate(zip(positions,raw,near)))
    pieces,n=ndimage.label(filled)
    component=VoidComponent(1,'fixed_reference_neighborhood',tuple(map(float,positions.mean(0) if len(positions) else anchor)),len(points)*h**3,points)
    report={'definition':'Union of probe sweeps from all atom-clear center components intersecting the fixed reference samples, clipped to the reference-neighborhood domain',
            'reference_margin_A':margin,'probe_radius_A':probe_radius,'reference_seed_distance_A':.9*h,
            'probe_center_components':len(retained),'sampled_volume_components':int(n),'grid_points':int(math.prod(shape)),
            'volume_A3':component.volume,'zero_definition':'No retained probe-swept samples in this specified region, not evidence of a closed transport pathway',
            'crop_boundary_exclusion_band_A':.9*h,
            'obstacle_definition':'protein heavy atoms' if obstacle_atoms is None else 'explicit physical obstacle atoms',
            'obstacle_atom_images':len(atoms),
            'limitations':'A fixed local measurement domain; chosen obstacles and spherical probes do not establish substrate or membrane accessibility'}
    cast=VoidCast((component,) if points else (),points,h,0.,'rolling',{'grid_basis':basis.tolist(),'grid_origin':(origin@basis).tolist(),'algorithm':'fixed_reference_neighborhood_probe_sweep',**report})
    return cast,component,report

"""Residue water proximity and geometric probe-accessible surface area.

Explicit-water observations and geometric exposure are separate measurements.
Distances alone are not hydrogen bonds, hydration energies or solvent access.

Two per-residue measurements are made for a set of protein heavy atoms:

Water contacts
    The number of distinct water molecules whose oxygen lies within
    ``cutoff`` (default 3.5 Å, inclusive) of any heavy atom of the residue,
    also split into backbone atoms (``N CA C O OXT OT1 OT2``), side-chain
    atoms, and nitrogen/oxygen atoms. Only waters present in the input are
    counted; a structure without deposited waters reports missing (NaN)
    hydration, not zero.
Probe-accessible surface area (SASA)
    Shrake-Rupley numerical quadrature: each atom's sphere of radius
    ``r_vdw + probe`` (probe 1.4 Å by default) is covered with evenly spread
    points (256 by default); the fraction not inside any neighbouring
    expanded sphere, times the sphere area, is the atom's SASA (Å²). All
    non-water heavy atoms (including ligands, lipids and ions) occlude.
    Waters do not.

Periodic systems are handled with MDAnalysis minimum-image distances.

Main entry points: :func:`static_hydration` (command ``crevice hydration``
on one structure), :class:`HydrationTopology` (per-frame measurement, used by
:func:`crevice.hydration_trajectory.analyze_hydration_trajectory`) and
:func:`sampled_sasa`.

Examples
--------
>>> from crevice.models import Atom, StructureFrame
>>> from crevice.hydration import static_hydration
>>> frame = StructureFrame((
...     Atom(1, "N", "SER", "A", 1, 0.0, 0.0, 0.0, "N"),
...     Atom(2, "OG", "SER", "A", 1, 1.5, 0.0, 0.0, "O"),
...     Atom(3, "O", "HOH", "W", 10, 4.0, 0.0, 0.0, "O", hetero=True),
... ))
>>> result = static_hydration(frame)
>>> result["residue_ids"], result["metrics"]["water_count"].tolist()
(('A:SER1',), [1.0])
>>> result["metrics"]["sidechain_water_count"].tolist(), result["water_status"]
([1.0], 'observed_explicit_water')
"""
from dataclasses import dataclass
import math
from .radii import RadiusSet, radii_option

#: Residue names treated as water (compared in upper case).
WATER_RESNAMES = frozenset({'HOH','WAT','SOL','TIP3','TIP3P','TIP4','TIP4P','TIP5','TIP5P','SPC','SPCE','H2O','DOD'})
#: MDAnalysis selection built from :data:`WATER_RESNAMES`.
DEFAULT_WATER_SELECTION = 'resname ' + ' '.join(sorted(WATER_RESNAMES))
#: Atom names counted as backbone for the backbone/side-chain split.
BACKBONE_NAMES = frozenset({'N','CA','C','O','OXT','OT1','OT2'})
#: Names of the per-residue hydration metrics, i.e. distinct-water counts (all,
#: backbone, side-chain and N/O atoms) and SASA in Å² (the same four atom sets).
METRICS = ('water_count','backbone_water_count','sidechain_water_count','NO_water_count',
           'sasa_A2','backbone_sasa_A2','sidechain_sasa_A2','NO_sasa_A2')


def checked_box(box, cutoff=0.):
    """Validate an MDAnalysis unit cell for periodic distance calculations.

    Parameters
    ----------
    box : array_like of 6 floats or None
        ``[a, b, c, alpha, beta, gamma]`` in Å and degrees.
    cutoff : float, default 0.0
        Largest distance (Å) that will be searched; it must be smaller than half
        the minimum cell height so minimum-image distances are unambiguous.

    Returns
    -------
    ndarray or None
        The box as a float array, or ``None`` when ``box`` is ``None``.

    Raises
    ------
    ValueError
        For a malformed or singular cell, or a cutoff that is too large.
    """
    import numpy as np
    if box is None:return None
    box=np.asarray(box,dtype=float)
    if box.shape!=(6,) or not np.isfinite(box).all() or np.any(box[:3]<=0) or np.any((box[3:]<=0)|(box[3:]>=180)):
        raise ValueError('Periodic hydration requires a finite, positive, valid six-value cell')
    from MDAnalysis.lib.mdamath import triclinic_vectors
    vectors=triclinic_vectors(box).astype(float);volume=abs(np.linalg.det(vectors))
    areas=np.linalg.norm(np.cross(vectors[[1,2,0]],vectors[[2,0,1]]),axis=1)
    if volume<=0 or np.any(areas<=0):raise ValueError('Periodic cell is singular')
    if cutoff>=.5*min(volume/areas):raise ValueError('Hydration cutoff/surface neighbourhood must be smaller than half the minimum cell height')
    return box


def minimum_vectors(vectors,box):
    """Apply the minimum-image convention to displacement vectors.

    Returns the vectors unchanged when ``box`` is ``None``; otherwise uses
    :func:`MDAnalysis.lib.distances.minimize_vectors`.
    """
    import numpy as np
    vectors=np.asarray(vectors,dtype=float)
    if box is None:return vectors
    from MDAnalysis.lib.distances import minimize_vectors
    return minimize_vectors(vectors,box)


def distance_pairs(left,right,cutoff,box=None):
    """All pairs of points closer than a cutoff.

    Parameters
    ----------
    left, right : array_like, shape (n, 3) and (m, 3)
        Coordinates in Å.
    cutoff : float
        Maximum distance in Å (inclusive).
    box : array_like, optional
        Unit cell for minimum-image distances (MDAnalysis ``capped_distance``);
        a KD-tree search is used without a box.

    Returns
    -------
    pairs : ndarray of int, shape (k, 2)
        Index pairs ``(i, j)`` into ``left`` and ``right``.
    distances : ndarray, shape (k,)
        Distances in Å.
    """
    import numpy as np
    left=np.asarray(left,dtype=float).reshape((-1,3));right=np.asarray(right,dtype=float).reshape((-1,3))
    if not len(left) or not len(right):return np.empty((0,2),int),np.empty(0)
    if box is not None:
        from MDAnalysis.lib.distances import capped_distance
        return capped_distance(left,right,max_cutoff=cutoff,box=box,return_distances=True)
    from scipy.spatial import cKDTree
    neighbours=cKDTree(left).query_ball_tree(cKDTree(right),cutoff)
    pairs=np.asarray([(i,j) for i,ids in enumerate(neighbours) for j in ids],dtype=int).reshape((-1,2))
    return pairs,np.linalg.norm(left[pairs[:,0]]-right[pairs[:,1]],axis=1)


def sphere_directions(count):
    """Evenly spread unit vectors on a sphere (golden-angle spiral).

    Parameters
    ----------
    count : int
        Number of directions, at least 32.

    Returns
    -------
    ndarray, shape (count, 3)

    Raises
    ------
    ValueError
        If ``count`` is not an integer of at least 32.
    """
    import numpy as np
    if isinstance(count,bool) or int(count)!=count or count<32:raise ValueError('SASA requires at least 32 sphere samples')
    k=np.arange(int(count));z=1-2*(k+.5)/count;phi=k*math.pi*(3-math.sqrt(5));r=np.sqrt(1-z*z)
    return np.column_stack([r*np.cos(phi),r*np.sin(phi),z])


def sampled_sasa(target_xyz,environment_xyz,environment_radii,self_indices,*,probe=1.4,samples=256,box=None,orientation=None):
    """Shrake-Rupley quadrature with nonwater heavy-atom obstacles.

    Periodic neighbour vectors are minimized before optional reference-frame
    rotation; rotating the dot grid with the fit avoids rigid-motion noise.

    Parameters
    ----------
    target_xyz : array_like, shape (n, 3)
        Atoms whose SASA is wanted (Å).
    environment_xyz : array_like, shape (m, 3)
        All occluding atoms (Å); includes the targets.
    environment_radii : array_like, shape (m,)
        Van der Waals radii of the environment atoms (Å).
    self_indices : array_like of int, shape (n,)
        Row of each target atom in the environment arrays (so an atom does not
        occlude itself).
    probe : float, default 1.4
        Probe radius in Å.
    samples : int, default 256
        Points per atom sphere (at least 32). More points reduce quadrature
        noise (roughly ``1/samples`` of the sphere area per point).
    box : array_like, optional
        Unit cell for periodic neighbours.
    orientation : array_like, shape (3, 3), optional
        Orthonormal rotation applied to neighbour vectors, so that the sample
        points rotate with an alignment fit.

    Returns
    -------
    ndarray, shape (n,)
        SASA per target atom in Å².

    Raises
    ------
    ValueError
        For a non-positive probe, invalid radii, a non-orthonormal orientation,
        or an invalid box.
    """
    import numpy as np
    if not math.isfinite(probe) or probe<=0:raise ValueError('SASA probe must be finite and positive')
    directions=sphere_directions(samples);xyz=np.asarray(target_xyz,dtype=float);env=np.asarray(environment_xyz,dtype=float)
    expanded=np.asarray(environment_radii,dtype=float)+probe;self_indices=np.asarray(self_indices,dtype=int)
    if np.any(expanded<=0) or not np.isfinite(expanded).all():raise ValueError('Invalid SASA radii')
    box=checked_box(box,2*expanded.max())
    rotation=np.eye(3) if orientation is None else np.asarray(orientation,dtype=float)
    if rotation.shape!=(3,3) or not np.allclose(rotation@rotation.T,np.eye(3),atol=1e-8):raise ValueError('SASA orientation must be orthonormal')
    pairs,_=distance_pairs(xyz,env,2*expanded.max(),box)
    neighbours=[[] for _ in xyz]
    for i,j in pairs:
        if j!=self_indices[i]:neighbours[i].append(j)
    area=np.zeros(len(xyz))
    for i,ids in enumerate(neighbours):
        radius=expanded[self_indices[i]]
        if not ids:area[i]=4*math.pi*radius**2;continue
        vectors=minimum_vectors(env[ids]-xyz[i],box)@rotation
        thresholds=(radius**2+np.sum(vectors*vectors,axis=1)-expanded[ids]**2)/(2*radius)
        free=~np.any(directions@vectors.T>thresholds[None,:]+1e-10,axis=1)
        area[i]=free.mean()*4*math.pi*radius**2
    return area


@dataclass
class HydrationTopology:
    """Immutable atom metadata; coordinates are supplied one frame at a time.

    Attributes
    ----------
    target_atoms : tuple of Atom
        Protein heavy atoms whose residues are measured.
    environment_radii : array_like
        Radii (Å) of all occluding (non-water heavy) atoms.
    self_indices : array_like of int
        Row of each target atom in the environment arrays.
    water_ids : tuple
        One stable identity per water oxygen; must be unique.
    residue_ids : tuple of str
        Sorted residue labels of the target atoms (derived).

    Raises
    ------
    ValueError
        If there are no target atoms or water identities repeat.
    """
    target_atoms: tuple
    environment_radii: object
    self_indices: object
    water_ids: tuple

    def __post_init__(self):
        import numpy as np
        if not self.target_atoms:raise ValueError('Hydration requires at least one selected protein heavy atom')
        self.residue_ids=tuple(sorted({a.residue_key.label for a in self.target_atoms}))
        lookup={key:i for i,key in enumerate(self.residue_ids)}
        self.atom_residue=np.asarray([lookup[a.residue_key.label] for a in self.target_atoms])
        self.backbone=np.asarray([a.name.upper() in BACKBONE_NAMES for a in self.target_atoms])
        self.NO=np.asarray([a.element.upper() in {'N','O'} for a in self.target_atoms])
        if len(set(self.water_ids))!=len(self.water_ids):raise ValueError('Each water identity must have exactly one selected oxygen')

    def measure(self,xyz,environment_xyz,water_xyz,*,cutoff=3.5,sasa_probe=1.4,sasa_samples=256,box=None,orientation=None):
        """Measure water contacts and SASA for one set of coordinates.

        Parameters
        ----------
        xyz : array_like, shape (n_targets, 3)
            Target atom coordinates (Å).
        environment_xyz : array_like, shape (n_environment, 3)
            Occluding atom coordinates (Å).
        water_xyz : array_like, shape (n_waters, 3)
            Water oxygen coordinates (Å).
        cutoff : float, default 3.5
            Water-oxygen contact distance in Å (inclusive).
        sasa_probe : float, default 1.4
            SASA probe radius in Å.
        sasa_samples : int, default 256
            Points per atom sphere.
        box : array_like, optional
            Unit cell for periodic distances.
        orientation : array_like, shape (3, 3), optional
            Rotation passed to :func:`sampled_sasa`.

        Returns
        -------
        dict
            ``metrics``: arrays over ``residue_ids`` for every name in
            :data:`METRICS` (counts of distinct waters, NaN when no waters are
            present; SASA in Å²). ``memberships``: set of water IDs per residue.
            ``bridges``: residue index pairs that share a water, mapped to the
            shared water IDs (a shared water within the cutoff of both, not a
            hydrogen-bond bridge). ``nearest_contact_oxygen_distance_A``,
            ``water_population`` and ``water_status``.

        Raises
        ------
        ValueError
            For array sizes that do not match the topology, non-finite
            coordinates, or an invalid cutoff or box.
        """
        import numpy as np
        if not math.isfinite(cutoff) or cutoff<=0:raise ValueError('Water contact cutoff must be finite and positive')
        xyz=np.asarray(xyz,dtype=float).reshape((-1,3));env=np.asarray(environment_xyz,dtype=float).reshape((-1,3));waters=np.asarray(water_xyz,dtype=float).reshape((-1,3))
        if len(xyz)!=len(self.target_atoms) or len(waters)!=len(self.water_ids) or len(env)!=len(self.environment_radii):raise ValueError('Hydration topology and coordinates differ')
        if not all(np.isfinite(a).all() for a in [xyz,env,waters]):raise ValueError('Hydration coordinates must be finite')
        box=checked_box(box,cutoff);pairs,dist=distance_pairs(xyz,waters,cutoff,box)
        groups=[[set() for _ in self.residue_ids] for _ in range(4)]
        closest=np.full(len(self.residue_ids),np.nan)
        for (i,j),distance in zip(pairs,dist):
            residue=self.atom_residue[i];groups[0][residue].add(int(j));groups[1 if self.backbone[i] else 2][residue].add(int(j))
            if self.NO[i]:groups[3][residue].add(int(j))
            if not np.isfinite(closest[residue]) or distance<closest[residue]:closest[residue]=distance
        metrics={key:np.asarray([len(s) for s in groups[k]],dtype=float) if len(waters) else np.full(len(self.residue_ids),np.nan)
                 for k,key in enumerate(METRICS[:4])}
        areas=sampled_sasa(xyz,env,self.environment_radii,self.self_indices,probe=sasa_probe,samples=sasa_samples,box=box,orientation=orientation)
        for key,mask in [('sasa_A2',np.ones(len(xyz),bool)),('backbone_sasa_A2',self.backbone),('sidechain_sasa_A2',~self.backbone),('NO_sasa_A2',self.NO)]:
            metrics[key]=np.bincount(self.atom_residue[mask],weights=areas[mask],minlength=len(self.residue_ids))
        memberships=[{self.water_ids[j] for j in s} for s in groups[0]]
        by_water={}
        for i,ids in enumerate(memberships):
            for water in ids:by_water.setdefault(water,[]).append(i)
        bridges={}
        from itertools import combinations
        for water,ids in by_water.items():
            for pair in combinations(sorted(ids),2):bridges.setdefault(pair,[]).append(water)
        return {'metrics':metrics,'memberships':memberships,'bridges':bridges,'nearest_contact_oxygen_distance_A':closest,
                'water_population':len(waters),'water_status':'observed_explicit_water' if len(waters) else 'unavailable_no_explicit_water'}


@radii_option
def static_hydration(frame,*,residue_ids=None,box=None,cutoff=3.5,sasa_probe=1.4,sasa_samples=256,water_resnames=WATER_RESNAMES, radii: RadiusSet | str | None = None):
    """Use observed waters; absence of deposited waters is missing hydration.

    Waters are the oxygen atoms of residues named in ``water_resnames`` with
    occupancy absent or positive. Occluding atoms are every non-hydrogen,
    non-water atom (ligands, lipids and ions included). Measured atoms are the
    non-HETATM heavy atoms of the focus residues.

    Parameters
    ----------
    frame : StructureFrame
    residue_ids : iterable of str, optional
        Residue labels to measure; all protein residues by default.
    box : array_like, optional
        Unit cell for periodic distances (for example from an MD snapshot).
    cutoff : float, default 3.5
        Water-oxygen contact distance in Å.
    sasa_probe : float, default 1.4
        SASA probe radius in Å.
    sasa_samples : int, default 256
        Points per atom sphere.
    water_resnames : iterable of str, default WATER_RESNAMES
        Residue names treated as water (compared in upper case).
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call (SASA environment radii) comes from that set, and
        results with a ``metadata`` dict record it as ``metadata["radii"]``;
        see :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict
        The :meth:`HydrationTopology.measure` result plus ``residue_ids`` (the
        measured residues), ``untyped_focus_residues``, ``source``,
        ``partial_occupancy_waters`` and ``settings`` (which repeats
        ``untyped_focus_residues`` so it reaches the written JSON).

    Raises
    ------
    ValueError
        If ``residue_ids`` contains labels absent from the frame, if none of
        them is a measurable protein residue, or as for
        :meth:`HydrationTopology.measure`.

    Notes
    -----
    A requested residue that is present but has no non-HETATM heavy atom is
    not measured and is listed in ``untyped_focus_residues`` with a reason:
    ``"water"``, ``"hetero_record"`` (for example modified or D-amino-acid
    residues deposited as HETATM in a PDB file, such as DLE/DVA in 1GRM,
    ligands and ions) or ``"no_heavy_atoms"``. Its heavy atoms still occlude
    SASA unless it is water. Load mmCIF with the Gemmi reader to treat
    modified polymer residues as protein (see :func:`crevice.mmcif.load_mmcif`).
    """
    import numpy as np
    from .radii import atom_vdw_radius
    water_resnames=set(water_resnames)
    waters=[a for a in frame.atoms if a.resname.upper() in water_resnames and a.element.upper()=='O' and (a.occupancy is None or a.occupancy>0)]
    environment=[a for a in frame.atoms if a.element.upper() not in {'H','D'} and a.resname.upper() not in water_resnames]
    protein=[(i,a) for i,a in enumerate(environment) if not a.hetero]
    known={a.residue_key.label for _,a in protein};untyped=[]
    if residue_ids is None:focus=known
    else:
        requested=set(residue_ids);present={a.residue_key.label for a in frame.atoms}
        if requested-present:raise ValueError('Hydration focus residue identities are not present in the structure: '+str(sorted(requested-present)))
        # Present but not measurable as protein: reported explicitly, never silently dropped.
        for label in sorted(requested-known):
            members=[a for a in frame.atoms if a.residue_key.label==label]
            reason=('water' if any(a.resname.upper() in water_resnames for a in members) else
                    'hetero_record' if any(a.hetero and a.element.upper() not in {'H','D'} for a in members) else 'no_heavy_atoms')
            untyped.append({'residue':label,'reason':reason})
        focus=requested&known
        if not focus:raise ValueError('Hydration focus has no protein residue to measure; untyped residues: '+str([u['residue'] for u in untyped]))
    targets=[(i,a) for i,a in protein if a.residue_key.label in focus]
    model=HydrationTopology(tuple(a for i,a in targets),np.asarray([atom_vdw_radius(a) for a in environment]),
                            np.asarray([i for i,a in targets]),tuple(a.residue_key.label for a in waters))
    result=model.measure([a.coord for _,a in targets],[a.coord for a in environment],[a.coord for a in waters],
                         cutoff=cutoff,sasa_probe=sasa_probe,sasa_samples=sasa_samples,box=box)
    result.update(residue_ids=model.residue_ids,source=frame.source,untyped_focus_residues=untyped,
                  partial_occupancy_waters=sum(a.occupancy is not None and a.occupancy<1 for a in waters),
                  settings={'contact_cutoff_A':cutoff,'sasa_probe_A':sasa_probe,'sasa_sphere_points':sasa_samples,
                            'periodic':box is not None,'source_model_index':frame.model_index,'source_frame_index':frame.frame_index,'environment_heavy_atoms':len(environment),
                            'environment_policy':'all retained nonwater heavy atoms, including ligands/lipids/ions',
                            'water_resnames':sorted(water_resnames),'focus':'all protein' if residue_ids is None else 'provided residue identities',
                            'untyped_focus_residues':untyped,
                            'untyped_focus_policy':'requested residues without non-HETATM heavy atoms (HETATM-marked modified or non-standard residues, ligands, ions, waters) are not measured; they stay in the occluding environment unless water'})
    return result

"""Three-dimensional two-probe voids, including solvent-facing pockets.

The outer probe defines a finite molecular envelope. The inner probe selects
empty space inside it; enclosed pockets are retained even without bulk access.
All coordinates are in Å. This is sampled atom-sphere geometry, not a
claim about water occupancy, permeation, or the number of biological channels.

Method in brief
    1. Lay a grid (default 0.5 Å) in a molecule-attached, right-handed basis
       (the principal axes of the atoms unless ``grid_basis`` is given) and
       compute every node's atom-surface clearance.
    2. *Outer probe* (radius ``outer_radius``, default 6 Å): nodes where it
       fits and that connect to the grid boundary are bulk solvent. The
       space it sweeps from there is outside the molecule; everything
       further in is the *envelope*. The outer probe cannot enter narrow
       clefts, so pockets and channels stay inside the envelope, and the
       edge of its sweep defines each pocket's mouth.
    3. *Inner probe* (radius ``probe_radius``, default 0.8 Å): its centres
       must have clearance ``>= probe_radius`` and lie deep enough in the
       envelope that the whole sphere fits inside it.
    4. Inner-probe centres are grouped into connected components
       (face-neighbour edges, each checked analytically along the segment),
       and the volume the probe sweeps from them is measured.
    5. ``selection_mode="dominant"`` instead qualifies substantial buried
       cores and fills them with variable-radius atom-clear spheres
       (:func:`crevice.cavity_selection.dominant_sphere_fill`).

Main entry points: :func:`rolling_probe_cast` (command ``crevice cast``),
:func:`select_cast_atoms`, :func:`chain_enclosed_regions`.

Assumptions and limitations
    Results depend on the grid spacing and both probe radii. Boundary voxels
    approximate a surface: their centres are atom-clear, but a voxel can
    still overlap an atom. A missing sampled path does not prove the
    continuum is disconnected. Component counts and volumes are geometric
    measurements; they do not identify biological channels.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import replace
import math
import operator

from .models import ChannelPoint, VoidCast, VoidComponent
from .radii import RadiusSet, atom_vdw_radius, radii_option


class _SphereQueries:
    """Exact batched point/axis-segment distances to unequal atom spheres.

    Atoms are grouped by radius, with one KD-tree per group, so nearest-surface
    queries are exact without a radius-bound search.

    Parameters
    ----------
    coords : ndarray, shape (n, 3)
        Atom centres (Å), in the grid's coordinate frame.
    radii : ndarray, shape (n,)
        Atomic radii (Å).
    """
    def __init__(self, coords, radii):
        import numpy as np
        from scipy.spatial import cKDTree
        self.groups = []
        for radius in np.unique(radii):
            ids = np.flatnonzero(radii == radius)
            self.groups.append((float(radius), ids, cKDTree(coords[ids])))

    def points(self, points):
        """Nearest-surface clearance and nearest-atom index for many points.

        Returns
        -------
        clearance : ndarray
            Signed distance to the nearest sphere surface (Å).
        nearest : ndarray of int32
            Index of that atom in the input arrays.
        """
        import numpy as np
        best = np.full(len(points), np.inf)
        nearest = np.zeros(len(points), dtype=np.int32)
        for radius, ids, tree in self.groups:
            distance, index = tree.query(points)
            take = distance-radius < best
            best[take] = distance[take]-radius
            nearest[take] = ids[index[take]]
        return best, nearest

    def segments(self, midpoints, axis, half_length, required):
        """Whether axis-aligned segments clear every sphere inflated by ``required``.

        Parameters
        ----------
        midpoints : ndarray, shape (m, 3)
            Segment midpoints.
        axis : int
            Axis (0, 1 or 2) the segments are parallel to.
        half_length : float
            Half the segment length (Å).
        required : float
            Probe radius added to every atom radius.

        Returns
        -------
        ndarray of bool
            ``True`` where no inflated sphere intersects the segment.
        """
        import numpy as np
        # Any intersecting sphere center lies within R + half_length of the
        # midpoint. Query all candidates, then minimize along the full segment.
        allowed = np.ones(len(midpoints), dtype=bool)
        for radius, _, tree in self.groups:
            active = np.flatnonzero(allowed)
            for start in range(0, len(active), 8192):
                rows = active[start:start+8192]
                hits = tree.query_ball_point(midpoints[rows], radius+required+half_length+1e-10)
                counts = np.fromiter(map(len, hits), dtype=np.int64, count=len(hits))
                if not counts.sum():
                    continue
                source = np.repeat(np.arange(len(rows)), counts)
                targets = np.concatenate([h for h in hits if h]).astype(int)
                delta = tree.data[targets]-midpoints[rows[source]]
                delta[:, axis] = np.maximum(np.abs(delta[:, axis])-half_length, 0)
                overlap = np.einsum('ij,ij->i', delta, delta) < (radius+required)**2-1e-12
                allowed[rows[source[overlap]]] = False
        return allowed


def _connected(mask, clearance, origin, spacing, probe, queries):
    """Face-connected labels with analytic segment clearance on every edge.

    Safely distant nodes are collapsed with ndimage; only nodes close to an
    atom boundary need explicit edges. The 1-Lipschitz clearance bound proves
    edges between endpoints with clearance >= probe + spacing/2 are valid.
    """
    import numpy as np
    from scipy import ndimage, sparse
    from scipy.sparse.csgraph import connected_components
    safe = mask & (clearance >= probe+spacing/2+1e-10)
    nodes, count = ndimage.label(safe)
    uncertain = mask & ~safe
    n = int(uncertain.sum())
    if not n:
        return nodes
    nodes[uncertain] = np.arange(count+1, count+n+1, dtype=np.int32)
    lefts, rights = [], []
    for axis in range(3):
        lo, hi = [slice(None)]*3, [slice(None)]*3
        lo[axis], hi[axis] = slice(None, -1), slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        connected = mask[lo] & mask[hi] & (uncertain[lo] | uncertain[hi])
        index = np.argwhere(connected)
        if not len(index):
            continue
        mids = origin+index*spacing
        mids[:, axis] += spacing/2
        valid = queries.segments(mids, axis, spacing/2, probe)
        index = index[valid]
        other = index.copy()
        other[:, axis] += 1
        lefts.append(nodes[tuple(index.T)])
        rights.append(nodes[tuple(other.T)])
    if lefts:
        a, b = np.concatenate(lefts), np.concatenate(rights)
        graph = sparse.coo_matrix((np.ones(len(a), dtype=np.int8), (a, b)),
                                  shape=(count+n+1, count+n+1)).tocsr()
        _, labels = connected_components(graph, directed=False)
        # Graph node zero is isolated, so its component remains background zero.
        return labels[nodes].astype(np.int32)
    return nodes


def _ray_enclosure(blocked):
    """Count atom-blocked sampled rays in the 26 lattice directions.

    For each voxel, the count is the number of the 26 lattice directions along
    which a ray through the grid meets a blocked voxel (a blocked voxel counts
    for itself). This optional burial crop separates interior space from
    surface grooves. It is a discrete enclosure descriptor, not an
    accessibility/path test.
    """
    import itertools
    import numpy as np
    count = np.zeros(blocked.shape, dtype=np.uint8)
    for direction in itertools.product((-1,0,1), repeat=3):
        if direction == (0,0,0):
            continue
        axis = next(i for i,d in enumerate(direction) if d)
        sign = direction[axis]
        hit = blocked.copy()
        indices = range(1,blocked.shape[axis]) if sign == 1 else range(blocked.shape[axis]-2,-1,-1)
        for i in indices:
            source, target = [slice(None)]*3, [slice(None)]*3
            source[axis], target[axis] = i-sign, i
            for j,d in enumerate(direction):
                if j != axis and d:
                    source[j] = slice(None,-1) if d == 1 else slice(1,None)
                    target[j] = slice(1,None) if d == 1 else slice(None,-1)
            hit[tuple(target)] |= hit[tuple(source)]
        count += hit
    return count


def select_cast_atoms(frame, selection='protein'):
    """Explicit input policy; preserve nonstandard polymer ATOM records.

    ``"protein"`` means non-HETATM heavy atoms for PDB/mmCIF. GRO input has
    already applied its MDAnalysis selection. ``"all"`` retains every heavy
    atom. This does not infer that a bound ligand is absent in the experimental
    conformation.

    Parameters
    ----------
    frame : StructureFrame
    selection : {"protein", "all"}, default "protein"
        Which atoms act as obstacles for a cast.

    Returns
    -------
    frame : StructureFrame
        Copy of ``frame`` containing only the selected atoms.
    report : dict
        ``selection``, ``input_atoms``, ``selected_atoms``,
        ``excluded_residue_atom_counts`` (by residue name) and ``definition``.

    Raises
    ------
    ValueError
        For an unknown selection or an empty result.
    """
    if selection not in {'protein', 'all'}:
        raise ValueError('selection must be protein or all')
    atoms = frame.selected_atoms(include_hydrogen=False, include_hetero=selection == 'all')
    keep = set(atoms)
    return replace(frame, atoms=atoms), {
        'selection': selection, 'input_atoms': len(frame.atoms), 'selected_atoms': len(atoms),
        'excluded_residue_atom_counts': dict(sorted(Counter(a.resname for a in frame.atoms if a not in keep).items())),
        'definition': 'non-HETATM heavy atoms' if selection == 'protein' else 'all heavy atoms',
    }


@radii_option
def rolling_probe_cast(frame, *, spacing=0.5, probe_radius=0.8, outer_radius=6.0,
                       min_component_volume=20.0, max_components=None, min_depth=0.0,
                       enclosure_fraction=0.0, selection_mode="all", core_fraction=.08,
                       max_grid_points=16_000_000, max_export_points=9000,
                       include_hydrogen=False, include_hetero=False, seed=None, seed_tolerance=.75,
                       grid_basis=None, grid_anchor=None,
                       radii: RadiusSet | str | None = None):
    """Fill voids inside a bulk-large-probe envelope.

    ``selection_mode='all'`` retains the legacy fixed-small-probe sweep.
    ``'dominant'`` qualifies substantial buried cores and fills variable-radius
    atom-clear spheres. The core and interface rules are reported in metadata;
    ``max_components`` is an upper limit applied after independent
    qualification.

    A larger rolling sphere can sweep exterior space but cannot enter narrow
    clefts. Its exterior-connected sweep defines the mouth boundary. Small
    spheres must fit wholly within the remaining envelope, including closed
    cavities. Their swept sample volume is clipped to that envelope and to
    non-negative atom clearance. Components refer to the small probe-centre
    graph; neighbouring swept regions are assigned without overlap.

    An optional seed retains its atom-clear probe-centre component within the
    envelope. No two-open-end requirement is imposed. No axis, circular cross
    section, channel count, or silent coarsening is imposed. Six-neighbour paths
    have exact segment checks; a missing sampled path does not prove continuum
    disconnection. Boundary voxels approximate a surface and can intersect atoms
    even though their centres are atom-clear.

    Parameters
    ----------
    frame : StructureFrame
        Structure; see :func:`select_cast_atoms` for choosing obstacles.
    spacing : float, default 0.5
        Grid spacing in Å. Never coarsened: a grid larger than
        ``max_grid_points`` raises instead.
    probe_radius : float, default 0.8
        Inner probe radius in Å (the space being measured).
    outer_radius : float, default 6.0
        Outer probe radius in Å that defines bulk solvent and pocket mouths.
        Must exceed ``probe_radius``. Larger values keep wider, shallower
        pockets inside the envelope.
    min_component_volume : float, default 20.0
        Components with a smaller swept volume (Å³) are discarded.
    max_components : int, optional
        Keep at most this many components (largest first). No limit by
        default. ``crevice cast`` uses 1.
    min_depth : float, default 0.0
        Extra distance (Å) that inner-probe spheres must keep from the outer
        envelope boundary; larger values trim pocket mouths.
    enclosure_fraction : float, default 0.0
        Fraction (0-1) of the 26 lattice ray directions that must hit an atom
        (inflated by ``probe_radius``) for a probe centre to count as buried.
        In ``"all"`` mode this is a hard crop of the probe centres; in
        ``"dominant"`` mode it qualifies cores only, and the sphere fill can
        extend beyond them. ``crevice cast`` uses 0.9.
    selection_mode : {"all", "dominant"}, default "all"
        ``"all"``: every inner-probe component. ``"dominant"``: substantial
        buried cores filled with spheres of radius ``min(clearance, depth)``
        (see :func:`crevice.cavity_selection.dominant_sphere_fill`).
        ``crevice cast`` uses ``"dominant"``.
    core_fraction : float, default 0.08
        Dominant mode: minimum core size relative to the largest core.
    max_grid_points : int, default 16000000
        Maximum number of grid nodes.
    max_export_points : int, default 9000
        Maximum number of display points (evenly strided subset).
    include_hydrogen : bool, default False
    include_hetero : bool, default False
        Atoms that act as obstacles.
    seed : Coord, optional
        Point (Å) inside the cavity of interest. Only the probe-centre component
        reached from it is kept. The seed must itself fit the probe, and must be
        joined by an atom-clear segment to a probe centre within
        ``seed_tolerance``.
    seed_tolerance : float, default 0.75
        Maximum seed-to-grid-centre distance in Å.
    grid_basis : array_like, shape (3, 3), optional
        Orthonormal grid axes as rows (made right-handed if necessary). Used
        to keep one grid across trajectory frames.
    grid_anchor : Coord, optional
        A point (Å) that is made a grid node, fixing the grid phase.
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
    VoidCast
        ``mode="rolling"`` and ``min_radius=0``. Components have kind
        ``"rolling_probe_void"`` or ``"dominant_cavity_region"`` and are ordered
        by volume (ties by source label). Each point's ``radius`` and
        ``raw_clearance`` are the sample's atom-surface clearance (Å) and ``t``
        its coordinate along the third grid axis from the atom centroid.
        ``metadata`` records the grid, probes, every candidate's
        ``retained``/``below_min_volume``/``component_limit`` decision, the
        seed report and definitions of mouth, volume and connectivity.

    Raises
    ------
    ValueError
        For invalid parameters, a grid over ``max_grid_points``, an invalid
        basis or anchor, or a seed that does not fit or cannot be connected.

    See Also
    --------
    select_cast_atoms : Obstacle-atom policy used by ``crevice cast``.
    crevice.voids.void_cast : Clearance-grid casts without an outer probe.

    Examples
    --------
    The enclosed interior of a closed spherical shell (radius 7 Å), sampled at
    0.5 Å:

    >>> import math
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.rolling import rolling_probe_cast
    >>> def shell_atom(i, n=250, golden=math.pi * (3 - math.sqrt(5))):
    ...     y = 1 - (2 * i + 1) / n
    ...     r = math.sqrt(1 - y * y)
    ...     return Atom(i + 1, "C", "ALA", "A", i + 1, 7 * math.cos(golden * i) * r,
    ...                 7 * y, 7 * math.sin(golden * i) * r, "C")
    >>> shell = StructureFrame(tuple(shell_atom(i) for i in range(250)))
    >>> cast = rolling_probe_cast(shell, spacing=0.5)
    >>> [(c.kind, c.volume) for c in cast.components]
    [('rolling_probe_void', 586.875)]
    """
    try:
        import numpy as np
        from scipy import ndimage
    except ImportError as exc:
        raise ValueError("Rolling-probe analysis requires NumPy and SciPy, which are required CREVICE dependencies but could not be imported; reinstall CREVICE (pip install --force-reinstall crevice)") from exc
    from .geometry import principal_axes
    for name, value in [('spacing', spacing), ('probe_radius', probe_radius), ('outer_radius', outer_radius)]:
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be finite and positive')
    if outer_radius <= probe_radius:
        raise ValueError('outer_radius must exceed probe_radius')
    if not math.isfinite(enclosure_fraction) or not 0 <= enclosure_fraction <= 1:
        raise ValueError('enclosure_fraction must be between 0 and 1')
    if not math.isfinite(min_depth) or min_depth < 0:
        raise ValueError('min_depth must be finite and non-negative')
    if not math.isfinite(min_component_volume) or min_component_volume < 0:
        raise ValueError('min_component_volume must be finite and non-negative')
    for name, value in [('max_grid_points', max_grid_points), ('max_export_points', max_export_points),
                         ('max_components', max_components)]:
        if value is not None and (isinstance(value, bool) or operator.index(value) < 1):
            raise ValueError(f'{name} must be a positive integer')
    if selection_mode not in {'all', 'dominant'}:
        raise ValueError('selection_mode must be all or dominant')
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    xyz = np.array([a.coord for a in atoms])
    basis = np.array(principal_axes(atoms) if grid_basis is None else grid_basis,dtype=float)
    if basis.shape!=(3,3) or not np.isfinite(basis).all() or not np.allclose(basis@basis.T,np.eye(3),atol=1e-8):
        raise ValueError('grid_basis must contain three orthonormal directions')
    # A right-handed, molecule-attached grid preserves rigid transformations
    # for nondegenerate principal directions, without rotating the input file.
    if np.linalg.det(basis) < 0:
        basis[2] *= -1
    coords = xyz @ basis.T
    radii = np.array([atom_vdw_radius(a) for a in atoms])
    padding = radii.max()+outer_radius+2*spacing
    origin = coords.min(axis=0)-padding
    if grid_anchor is not None:
        anchor=np.asarray(grid_anchor,dtype=float)
        if anchor.shape!=(3,) or not np.isfinite(anchor).all():
            raise ValueError('grid_anchor must be a finite three-dimensional coordinate')
        anchor=anchor@basis.T
        origin=anchor+np.floor((origin-anchor)/spacing)*spacing
    dims = tuple(int(n) for n in np.ceil((coords.max(axis=0)+padding-origin)/spacing).astype(int)+1)
    total = math.prod(dims)
    if total > max_grid_points:
        raise ValueError(f'Rolling-probe grid needs {total:,} points at {spacing:g} A; '
                         f'max_grid_points={max_grid_points:,}. Increase the limit or explicitly choose a coarser spacing.')
    queries = _SphereQueries(coords, radii)
    clearance = np.empty(total, dtype=np.float64)
    nearest = np.empty(total, dtype=np.int32)
    for start in range(0, total, 131072):
        stop = min(start+131072, total)
        ijk = np.array(np.unravel_index(np.arange(start, stop), dims)).T
        clearance[start:stop], nearest[start:stop] = queries.points(origin+ijk*spacing)
    clearance, nearest = clearance.reshape(dims), nearest.reshape(dims)
    large_labels = _connected(clearance >= outer_radius, clearance, origin, spacing, outer_radius, queries)
    boundary = np.unique(np.concatenate([large_labels[0].ravel(), large_labels[-1].ravel(),
                             large_labels[:,0].ravel(), large_labels[:,-1].ravel(),
                             large_labels[:,:,0].ravel(), large_labels[:,:,-1].ravel()]))
    bulk = np.isin(large_labels, boundary[boundary > 0])
    del large_labels
    distance_bulk = ndimage.distance_transform_edt(~bulk, sampling=spacing)
    envelope = distance_bulk > outer_radius+1e-9
    del bulk
    # The entire small sphere must fit inside the outer envelope; center-only
    # clipping creates thin surface artifacts where the probes overlap.
    centers = (distance_bulk > outer_radius+probe_radius+min_depth+1e-9) & (clearance >= probe_radius)
    seed_report=None
    if seed is not None:
        from .spatial import SpatialIndex
        seed_world=np.asarray(seed,dtype=float)
        if seed_world.shape!=(3,) or not np.isfinite(seed_world).all():
            raise ValueError('seed must be a finite three-dimensional coordinate')
        if not math.isfinite(seed_tolerance) or seed_tolerance<=0:
            raise ValueError('seed_tolerance must be finite and positive')
        seed_local=seed_world@basis.T
        raw,_=queries.points(seed_local[None,:])
        if raw[0]<probe_radius-1e-9:
            raise ValueError('Seed does not accommodate the requested probe radius')
        seed_mask=centers.copy()
        if selection_mode=='all' and enclosure_fraction:
            seed_mask &= _ray_enclosure(clearance < probe_radius)>=math.ceil(26*enclosure_fraction)
        center_labels=_connected(seed_mask,clearance,origin,spacing,probe_radius,queries)
        index=np.argwhere(seed_mask)
        distances=np.linalg.norm(origin+index*spacing-seed_local,axis=1)
        candidates=np.flatnonzero(distances<=seed_tolerance)
        if not len(candidates):raise ValueError('No interior probe center is close enough to the cavity seed')
        spatial=SpatialIndex(atoms);selected=None
        for candidate in candidates[np.argsort(distances[candidates],kind='stable')]:
            point=(origin+index[candidate]*spacing)@basis
            if spatial.segment_clearance(tuple(seed_world),tuple(point))>=probe_radius-1e-9:
                selected=int(candidate);break
        if selected is None:raise ValueError('No atom-clear segment connects the seed to the sampled interior')
        chosen=int(center_labels[tuple(index[selected])]);centers &= center_labels==chosen
        seed_report={'coordinate_A':seed_world.tolist(),'matched_center_A':((origin+index[selected]*spacing)@basis).tolist(),
                     'snap_distance_A':float(distances[selected]),'maximum_snap_distance_A':seed_tolerance,
                     'seed_clearance_A':float(raw[0]),'retained_probe_centers':int(centers.sum()),
                     'definition':'One interior probe-center component connected to the seed by atom-clear segments; no two-open-end requirement',
                     'accessibility_interpretation':'Connectivity within the chosen geometric envelope, not proof of substrate passage or solvent access through a membrane'}
        del center_labels,seed_mask,index,distances
    dominance = None
    if selection_mode == 'dominant':
        from .cavity_selection import dominant_sphere_fill
        swept_labels, center_counts, dominance = dominant_sphere_fill(
            atoms=atoms, coords=coords, clearance=clearance, centers=centers,
            depth=distance_bulk-outer_radius-min_depth, origin=origin, spacing=spacing,
            probe_radius=probe_radius, enclosure_fraction=enclosure_fraction, queries=queries,
            min_component_volume=min_component_volume, core_fraction=core_fraction)
        del centers
    else:
        if enclosure_fraction:
            centers &= _ray_enclosure(clearance < probe_radius) >= math.ceil(26*enclosure_fraction)
        labels = _connected(centers, clearance, origin, spacing, probe_radius, queries)
        center_counts = np.bincount(labels.ravel())
        del centers
        if labels.max():
            dist, nearest_center = ndimage.distance_transform_edt(labels == 0, sampling=spacing, return_indices=True)
            swept = envelope & (clearance >= 0) & (dist <= probe_radius+1e-9)
            swept_labels = np.where(swept, labels[tuple(nearest_center)], 0)
            del dist, nearest_center, swept
        else:
            swept_labels = labels.copy()
    del distance_bulk
    counts = np.bincount(swept_labels.ravel())
    order = sorted((i for i in range(1,len(counts)) if counts[i]), key=lambda i: (-counts[i], i))
    kept = [i for i in order if counts[i]*spacing**3 >= min_component_volume]
    if max_components is not None:
        kept = kept[:max_components]
    components, decisions = [], []
    atom_center = xyz.mean(0)
    objects = ndimage.find_objects(swept_labels)
    for label in order:
        decisions.append({'source_component_id': int(label), 'volume_A3': int(counts[label])*spacing**3,
                          'center_points': int(center_counts[label]),
                          'decision': 'retained' if label in kept else
                          'below_min_volume' if counts[label]*spacing**3 < min_component_volume else 'component_limit'})
    for new_id, label in enumerate(kept, 1):
        box = objects[label-1]
        index = np.argwhere(swept_labels[box] == label)+np.array([s.start for s in box])
        world = (origin+index*spacing) @ basis
        rr = clearance[tuple(index.T)]
        nn = nearest[tuple(index.T)]
        points = tuple(ChannelPoint(i, tuple(map(float,p)), float(r), float(r), float((p-atom_center)@basis[2]),
                                    atoms[int(n)].serial, atoms[int(n)].residue_key.label)
                       for i,(p,r,n) in enumerate(zip(world,rr,nn)))
        near_chains = Counter(atoms[int(n)].chain_id for n in nn)
        components.append(VoidComponent(new_id, 'dominant_cavity_region' if dominance is not None else 'rolling_probe_void', tuple(map(float,world.mean(0))),
                          len(points)*spacing**3, points,
                          nearest_residues=tuple(dict.fromkeys(p.nearest_residue for p in points))[:20],
                          metadata={'source_component_id': int(label), 'probe_center_points': int(center_counts[label]),
                                    'nearest_chain_point_counts': dict(near_chains),
                                    'span_world_A': np.ptp(world,axis=0).tolist(),
                                    'interpretation': 'selected geometric cavity region, not an assigned biological channel'}))
    full = tuple(p for c in components for p in c.points)
    stride = max(1, math.ceil(len(full)/max_export_points))
    display = tuple(replace(p,index=i) for i,p in enumerate(full[::stride]))
    return VoidCast(tuple(components), display, spacing, 0.0, 'rolling', {
        'algorithm': 'dominant_core_variable_sphere_fill' if dominance is not None else 'two_probe_exterior_envelope',
        'export_mode': 'rolling_probe_swept_fill', 'selection_mode': selection_mode, 'dominance': dominance, 'seed_selection':seed_report,
        'probe_radius_A': probe_radius, 'outer_radius_A': outer_radius,
        'requested_spacing_A': spacing, 'effective_spacing_A': spacing, 'spacing_adjusted': False,
        'grid_basis': basis.tolist(), 'grid_origin': (origin@basis).tolist(),
        'grid_anchor':list(map(float,grid_anchor)) if grid_anchor is not None else None,
        'fixed_grid_basis':grid_basis is not None, 'grid_dimensions': [int(n) for n in dims], 'grid_point_count': total,
        'axis_origin': xyz.mean(0).tolist(), 'axis_direction': basis[2].tolist(),
        'axis_definition': 'grid direction only; no through-channel assignment',
        'measurement_point_count': len(full), 'display_subsampled': stride > 1,
        'include_hydrogen': include_hydrogen, 'include_hetero': include_hetero,
        'connectivity': ('analytic 6-neighbor core/support paths and unions of atom-clear balls; region segmentation is distinct from whole-union connectivity' if dominance is not None else
                         '6-neighbor probe-center paths with analytic atom-sphere segment clearance'),
        'mouth_definition': 'boundary of exterior-connected outer-probe swept volume',
        'volume_definition': (dominance['fill_definition'] if dominance is not None else
                              'sampled small-probe swept void clipped to outer envelope; atom-clear sample centers'),
        'uncertainty': 'grid and probe dependent; voxel surfaces are approximate; no convergence or functional assignment implied',
        'selection': decisions, 'min_component_volume_A3': min_component_volume,
        'max_components': max_components, 'min_depth_A': min_depth,
        'enclosure_fraction': enclosure_fraction, 'enclosure_directions': 26,
        'enclosure_required_rays': math.ceil(26*enclosure_fraction),
        'enclosure_definition': _enclosure_definition(selection_mode, enclosure_fraction), 'source_component_count': len(order),
    })


def _enclosure_definition(selection_mode, enclosure_fraction):
    """Metadata text stating exactly how ``enclosure_fraction`` was used."""
    rays = math.ceil(26*enclosure_fraction)
    test = (f'a probe centre is buried when at least {rays} of the 26 lattice ray directions '
            '(ceil(26 * enclosure_fraction)), traced through the grid to its edge, meet a node whose '
            'atom clearance is below the probe radius (probe-inflated atoms sampled on the grid)')
    if not rays:
        return ('enclosure_fraction is 0, so the burial ray test (26 lattice ray directions meeting '
                'probe-inflated atoms sampled on the grid) is not applied: every probe centre inside the envelope is used')
    if selection_mode == 'dominant':
        return (test + '; in dominant mode it qualifies buried cores only (core and support selection); '
                'the variable-radius sphere fill can extend beyond those centres, so it is not a crop of the filled volume')
    return test + '; in all mode it is a hard crop of the probe centres before connectivity and the swept volume are formed'


def chain_enclosed_regions(frame, cast, *, fraction=.75):
    """Disjoint display annotation by single-chain ray enclosure.

    The connected components and measured union are untouched. Colours identify
    space enclosed by each chain, with a separate shared/unassigned region.
    These labels are geometric annotations, not a biological channel
    classifier.

    For each non-HETATM chain, the chain's atom spheres are rasterised on a grid
    in the cast's basis, and each cast sample receives the number of the 26 ray
    directions blocked by that chain alone. A sample belongs to the chain with
    the highest count, provided it reaches ``ceil(26 * fraction)``. Other
    samples are ``"shared"``; with three or more chains, samples that fewer than
    three chains each block in at least 6 directions become ``"other_void"``.

    Parameters
    ----------
    frame : StructureFrame
        The structure the cast was computed from.
    cast : VoidCast
        A cast whose metadata contains ``grid_basis`` (identity otherwise).
    fraction : float, default 0.75
        Required fraction of blocked rays for single-chain assignment.

    Returns
    -------
    regions : dict of str to tuple of ChannelPoint
        Keys ``"chain_<id>"``, ``"shared"`` and ``"other_void"`` (only those
        with points).
    summary : list of dict
        Name, point count, volume (Å³), definition and centre of each region.
        Both are empty for fewer than two chains or an empty cast.

    Raises
    ------
    ValueError
        If the annotation grid would exceed 16,000,000 points.
    """
    import numpy as np
    chains = sorted(set(a.chain_id for a in frame.atoms if not a.hetero))
    if len(chains) < 2 or not cast.components:
        return {}, []
    basis = np.asarray(cast.metadata.get('grid_basis', np.eye(3)))
    spacing = cast.spacing
    xyz = np.asarray([a.coord for a in frame.atoms]) @ basis.T
    origin = xyz.min(0)-4
    dims = tuple(int(n) for n in np.ceil((xyz.max(0)+4-origin)/spacing)+1)
    if math.prod(dims) > 16_000_000:
        raise ValueError('Chain annotation grid exceeds 16,000,000 points')
    points = tuple(p for c in cast.components for p in c.points)
    ijk = np.rint((np.asarray([p.position for p in points])@basis.T-origin)/spacing).astype(int)
    inside = np.all((ijk>=0) & (ijk<np.asarray(dims)),axis=1)
    best = np.zeros(len(points),dtype=np.uint8)
    assignment = np.full(len(points), -1, dtype=np.int16)
    contributors = np.zeros(len(points), dtype=np.uint8)
    for ci,chain in enumerate(chains):
        blocked = np.zeros(dims,dtype=bool)
        for atom,coord in zip(frame.atoms,xyz):
            if atom.chain_id != chain:
                continue
            radius = atom_vdw_radius(atom)
            low = np.maximum(np.ceil((coord-radius-origin)/spacing).astype(int),0)
            high = np.minimum(np.floor((coord+radius-origin)/spacing).astype(int)+1,dims)
            if np.any(high<=low):
                continue
            x,y,z = [origin[i]+np.arange(low[i],high[i])*spacing-coord[i] for i in range(3)]
            blocked[tuple(slice(a,b) for a,b in zip(low,high))] |= x[:,None,None]**2+y[None,:,None]**2+z[None,None,:]**2 < radius**2
        count = _ray_enclosure(blocked)
        scores = np.zeros(len(points),dtype=np.uint8)
        scores[inside] = count[tuple(ijk[inside].T)]
        contributors += scores >= 6
        take = (scores>=math.ceil(26*fraction)) & (scores>best)
        assignment[take], best[take] = ci, scores[take]
    if len(chains) >= 3:
        assignment[(assignment < 0) & (contributors < 3)] = -2
    regions = {}
    summary = []
    for ci,name in [*enumerate(chains),(-1,'shared'),(-2,'other_void')]:
        selected = tuple(p for i,p in enumerate(points) if assignment[i]==ci)
        if not selected:
            continue
        label = name if ci < 0 else f'chain_{name}'
        regions[label] = selected
        summary.append({'name':label,'points':len(selected),'volume_A3':len(selected)*spacing**3,
                        'definition': ('unassigned space with fewer than three contributing chains' if ci==-2 else
                                       'not enclosed by one chain; for trimers or larger, at least three chains each block >=6 rays') if ci<0 else
                        f'at least {math.ceil(26*fraction)} / 26 rays blocked by this chain alone',
                        'center':np.asarray([p.position for p in selected]).mean(0).tolist()})
    return regions, summary

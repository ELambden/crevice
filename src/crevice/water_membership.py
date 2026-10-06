"""Water-oxygen membership in each frame's measured cavity or region.

Two separate observables are reported for the same aligned oxygens:

* fixed reference neighbourhood: oxygen within the configured margin of the
  fixed occupied reference samples (the existing ``region_water_count``
  definition, computed by the same shared function);
* instantaneous measured region: oxygen centre inside an occupied voxel of the
  grid measured in that same frame.

This module is used by ``crevice cavity-trajectory --water-membership`` and
``crevice region-trajectory --water-membership``. Per frame,
:func:`frame_water_membership` (inside each geometry worker) reads the
oxygens with :class:`WaterSource` and calls :func:`water_membership`;
:func:`summarize_water_membership` builds all-frame arrays and statistics;
:func:`write_water_membership_bundle` writes the JSON, CSV tables and figure;
:func:`viewer_member_waters` supplies the orange member-water spheres of the
viewer scenes (group ``crevice_member_waters``).

Periodic images use the minimum image about the aligned reference anchor;
regions that reach half the minimum cell height are refused, and
``--pbc unwrap`` is refused because it moves protein independently of
solvent. Membership in a probe-swept geometric region is not binding,
hydrogen bonding, hydration free energy or permeation. Protein motion through
a fixed crop is not a hydration transition.
"""
from pathlib import Path
import json
import os

from .presentation import FigureText, annotate_option

MEMBERSHIP_TOLERANCES_A = (0.5, 1.0)
_KEY_OFFSET = 1 << 19
_KEY_BASE = 1 << 20

LIMITATIONS = [
    'Membership is an oxygen-centre-in-occupied-voxel test against the probe-swept measured region; it is not binding, hydrogen bonding, hydration free energy or permeation.',
    'The measured region depends on probe radius, grid spacing, crop margin, obstacles and alignment; oxygens in sites too narrow for the probe are not members.',
    'Fixed-reference and instantaneous counts are different observables; protein motion through a fixed crop is not a hydration transition.',
    'Unresolved geometry is missing, never zero; measured zero volume gives zero members and an undefined density.',
    '100 ps saved-frame spacing does not resolve water exchange; membership persistence is a saved-frame population, not a residence lifetime.',
    'Mean intervals are approximate block-bootstrap intervals conditional on representative stationary sampling; strata by geometry are descriptive only.',
]


def image_waters_around_anchor(waters, box, rotation, translation, anchor):
    """Aligned oxygen vectors from the reference anchor, periodic minimum image.

    Exactly the operation used by the existing fixed-region water count: the
    aligned anchor is mapped into the original frame by the inverse fit, oxygen
    vectors are minimum-imaged there and rotated into the reference frame.

    Parameters
    ----------
    waters : array_like, shape (n, 3)
        Oxygen coordinates in the original (unaligned) frame, Å.
    box : array_like or None
        MDAnalysis cell ``[a, b, c, alpha, beta, gamma]`` (Å, degrees); ``None``
        disables periodic imaging.
    rotation : array_like, shape (3, 3)
        Rotation of the frame's rigid fit (row-vector convention:
        ``aligned = original @ rotation + translation``).
    translation : array_like, shape (3,)
        Translation of the same fit, Å.
    anchor : array_like, shape (3,)
        Reference anchor in the aligned (reference) frame, Å.

    Returns
    -------
    numpy.ndarray, shape (n, 3)
        Oxygen vectors relative to ``anchor`` in the aligned frame, Å (add
        ``anchor`` for aligned positions).
    """
    import numpy as np
    from .hydration import minimum_vectors
    rotation = np.asarray(rotation, dtype=float); translation = np.asarray(translation, dtype=float)
    anchor = np.asarray(anchor, dtype=float)
    original_anchor = (anchor-translation)@rotation.T
    return minimum_vectors(np.asarray(waters, dtype=float)-original_anchor, box)@rotation


def fixed_reference_water_mask(relative, anchor, tree, margin, bound):
    """Existing fixed-neighbourhood rule: aligned oxygen within margin of a reference sample.

    Parameters
    ----------
    relative : array_like, shape (n, 3)
        Oxygen vectors from the anchor (:func:`image_waters_around_anchor`), Å.
    anchor : array_like, shape (3,)
        Reference anchor, Å.
    tree : scipy.spatial.cKDTree
        Tree of the fixed occupied reference sample centres (aligned frame).
    margin : float
        Region margin, Å: an oxygen counts when its nearest reference sample is
        at most this far (plus 1e-9 Å).
    bound : float
        Pre-filter radius, Å: oxygens farther than this from the anchor are not
        tested (it must cover every reference sample plus ``margin``).

    Returns
    -------
    numpy.ndarray of bool, shape (n,)
        True for oxygens in the fixed reference neighbourhood.
    """
    import numpy as np
    relative = np.asarray(relative, dtype=float).reshape((-1, 3))
    mask = np.zeros(len(relative), bool)
    near = np.flatnonzero(np.linalg.norm(relative, axis=1) <= bound)
    if len(near):
        mask[near] = tree.query(relative[near]+anchor)[0] <= margin+1e-9
    return mask


def lattice_indices(points, anchor, basis, spacing):
    """Half-up voxel index on the fixed reference lattice: voxel k owns [k-1/2, k+1/2).

    Parameters
    ----------
    points : array_like, shape (n, 3)
        Aligned positions, Å.
    anchor : array_like, shape (3,)
        Lattice origin (the reference anchor), Å.
    basis : array_like, shape (3, 3)
        Orthonormal lattice axes as rows (the third row is the reference axis).
    spacing : float
        Lattice spacing h, Å.

    Returns
    -------
    numpy.ndarray of int64, shape (n, 3)
        ``floor((points - anchor) @ basis.T / h + 0.5)``: exact ties go up.

    Examples
    --------
    >>> import numpy as np
    >>> from crevice.water_membership import lattice_indices
    >>> lattice_indices([[0.25, -0.25, 0.74]], np.zeros(3), np.eye(3), 0.5).tolist()
    [[1, 0, 1]]
    """
    import numpy as np
    local = (np.asarray(points, dtype=float).reshape((-1, 3))-np.asarray(anchor))@np.asarray(basis).T/float(spacing)
    return np.floor(local+.5).astype(np.int64)


def voxel_keys(index):
    """Encode integer lattice indices as one sortable int64 key per voxel.

    Each index component is offset by 2**19 and packed in base 2**20, so keys are
    unique for components in (-2**19, 2**19); they are used for fast membership
    tests with :func:`numpy.searchsorted`.

    Parameters
    ----------
    index : array_like of int, shape (n, 3)
        Lattice indices, for example from :func:`lattice_indices`.

    Returns
    -------
    numpy.ndarray of int64, shape (n,)
        Keys; the same voxel always gives the same key.

    Raises
    ------
    ValueError
        If an index component is at least 2**19 in magnitude.

    Examples
    --------
    >>> from crevice.water_membership import voxel_keys
    >>> keys = voxel_keys([[0, 0, 0], [0, 0, 1], [0, 0, 0]])
    >>> bool(keys[0] == keys[2]), bool(keys[1] > keys[0])
    (True, True)
    """
    import numpy as np
    index = np.asarray(index, dtype=np.int64).reshape((-1, 3))
    if len(index) and (np.abs(index).max() >= _KEY_OFFSET):
        raise ValueError('Lattice index exceeds the supported membership key range')
    shifted = index+_KEY_OFFSET
    return (shifted[:, 0]*_KEY_BASE+shifted[:, 1])*_KEY_BASE+shifted[:, 2]


def _contains(sorted_keys, keys):
    import numpy as np
    if not len(sorted_keys) or not len(keys):
        return np.zeros(len(keys), bool)
    position = np.clip(np.searchsorted(sorted_keys, keys), 0, len(sorted_keys)-1)
    return sorted_keys[position] == keys


def _minimum_cell_height(box):
    import numpy as np
    from MDAnalysis.lib.mdamath import triclinic_vectors
    from .hydration import checked_box
    checked_box(box, 0.)
    vectors = triclinic_vectors(box).astype(float); volume = abs(np.linalg.det(vectors))
    areas = np.linalg.norm(np.cross(vectors[[1, 2, 0]], vectors[[2, 0, 1]]), axis=1)
    return float(min(volume/areas))


def _domain_keys(reference, margin, probe_radius, max_grid_points):
    """Sorted voxel keys of the fixed discretised reference-region measurement domain."""
    import numpy as np
    from .regional_volume import regional_grid
    cache = getattr(reference, '_membership_domain_cache', {})
    key = (margin, probe_radius, max_grid_points)
    if key not in cache:
        index, _, _, distance, _, _ = regional_grid(reference, margin, probe_radius, max_grid_points)
        cache[key] = np.sort(voxel_keys(index[distance.ravel() <= margin+1e-9]))
        reference._membership_domain_cache = cache
    return cache[key]


def _reference_tree(reference, margin):
    import numpy as np
    from scipy.spatial import cKDTree
    cache = getattr(reference, '_membership_tree_cache', None)
    if cache is None or cache[0] != margin:
        points = np.asarray(reference.points, dtype=float)
        bound = float(np.linalg.norm(points-reference.anchor, axis=1).max()+margin)
        cache = (margin, cKDTree(points), bound); reference._membership_tree_cache = cache
    return cache[1], cache[2]


def water_membership(waters, water_ids, box, fit, reference, component, *, region_margin=2.,
                     geometry_mode='reference-region', probe_radius=1.4, max_grid_points=2_000_000,
                     tolerances=MEMBERSHIP_TOLERANCES_A):
    """Paired fixed-reference and instantaneous-region membership for one frame.

    1. Oxygens are minimum-imaged about the anchor and aligned with the frame's
       own rigid fit (:func:`image_waters_around_anchor`).
    2. Fixed reference count: :func:`fixed_reference_water_mask`.
    3. Instantaneous members: each oxygen centre is assigned to its half-up
       voxel on the fixed reference lattice (:func:`lattice_indices`); it is a
       member when that voxel is one of the occupied samples of ``component``.
    4. Cross-tabulation (both / fixed only / cavity only), sensitivity counts
       (voxel rule OR within each tolerance of an occupied sample centre), the
       count per axial lattice plane and, in ``reference-region`` mode, the count
       in the fixed discretised measurement domain.

    Parameters
    ----------
    waters : array_like, shape (n, 3)
        Original-frame oxygen coordinates, Å.
    water_ids : array_like of int, shape (n,)
        Stable water identities (MD topology residue indices).
    box : array_like or None
        Periodic cell; ``None`` disables imaging.
    fit : dict
        The rigid transform already used for this frame's geometry:
        ``rotation_rows`` (3 x 3) and ``translation_A`` (3,).
    reference : CavityReference
        Fixed reference cavity (anchor, basis, spacing, occupied samples, axial
        plane indices).
    component : VoidComponent or None
        The measured region of this frame; ``None`` when the geometry is
        unresolved, no points when the measured volume is zero. Its samples must
        lie on the reference lattice.
    region_margin : float, default 2.0
        Margin of the fixed reference neighbourhood, Å.
    geometry_mode : str, default "reference-region"
        ``"reference-region"`` also counts oxygens in the fixed measurement
        domain; any other mode reports ``domain_voxel_count`` as ``None``.
    probe_radius : float, default 1.4
        Probe radius of the regional grid, Å (selects the cached domain).
    max_grid_points : int, default 2000000
        Grid-size limit of the regional domain.
    tolerances : sequence of float, default (0.5, 1.0)
        Distance tolerances, Å, of the sensitivity counts (sensitivity only).

    Returns
    -------
    dict
        ``status`` (``measured`` or ``unresolved_geometry``),
        ``water_population``, ``fixed_reference_count``, ``domain_voxel_count``,
        ``cavity_member_count``, sorted ``member_water_ids`` and their aligned
        ``member_positions_A``, ``both_count``, ``fixed_only_count``,
        ``cavity_only_count``, ``domain_noncavity_count``, ``tolerance_counts``
        (keyed by tolerance as a string), ``axial_plane_counts``,
        ``outside_profile_count`` and ``volume_A3``. Unresolved geometry gives
        ``None`` for every instantaneous field, never zero.

    Raises
    ------
    ValueError
        For mismatched or non-finite inputs, samples off the reference lattice
        or duplicated, a region reaching half the minimum periodic cell height,
        or members outside the fixed measurement domain.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    waters = np.asarray(waters, dtype=float).reshape((-1, 3)); water_ids = np.asarray(water_ids, dtype=np.int64)
    if len(waters) != len(water_ids):
        raise ValueError('Water coordinates and identities differ in length')
    if not np.isfinite(waters).all():
        raise ValueError('Water coordinates must be finite')
    h = float(reference.spacing); anchor = np.asarray(reference.anchor, dtype=float); basis = np.asarray(reference.basis, dtype=float)
    tree, bound = _reference_tree(reference, region_margin)
    positions = np.asarray([p.position for p in component.points], dtype=float).reshape((-1, 3)) if component is not None else np.empty((0, 3))
    # Every sample, domain voxel and tolerance shell tested lies within this reach of the anchor.
    reach = max(bound, float(np.linalg.norm(positions-anchor, axis=1).max()) if len(positions) else 0.)+h+max(map(float, tolerances), default=0.)
    if box is not None and reach >= .5*_minimum_cell_height(box):
        raise ValueError('Measured region extends beyond half the minimum periodic cell height; periodic water membership would be ambiguous')
    rotation = np.asarray(fit['rotation_rows'], dtype=float); translation = np.asarray(fit['translation_A'], dtype=float)
    relative = image_waters_around_anchor(waters, box, rotation, translation, anchor) if len(waters) else np.empty((0, 3))
    fixed = fixed_reference_water_mask(relative, anchor, tree, region_margin, bound)
    near = np.flatnonzero(np.linalg.norm(relative, axis=1) <= reach)
    result = {'fixed_reference_count': int(fixed.sum()), 'water_population': int(len(waters))}
    near_keys = voxel_keys(lattice_indices(relative[near]+anchor, anchor, basis, h))
    if geometry_mode == 'reference-region':
        domain = _contains(_domain_keys(reference, region_margin, probe_radius, max_grid_points), near_keys)
        result['domain_voxel_count'] = int(domain.sum())
    else:
        domain = None; result['domain_voxel_count'] = None
    if component is None:
        result.update(status='unresolved_geometry', cavity_member_count=None, member_water_ids=None, member_positions_A=None,
                      both_count=None, fixed_only_count=None, cavity_only_count=None,
                      domain_noncavity_count=None, tolerance_counts=None, axial_plane_counts=None,
                      outside_profile_count=None, volume_A3=None)
        return result
    local = (positions-anchor)@basis.T/h; index = np.rint(local).astype(np.int64)
    if len(index) and not np.allclose(local, index, atol=1e-5, rtol=0):
        raise ValueError('Water membership requires the measured region on the fixed reference lattice')
    occupied = np.unique(voxel_keys(index))
    if len(occupied) != len(index):
        raise ValueError('Measured region contains duplicate lattice samples')
    member_near = _contains(occupied, near_keys)
    member = np.zeros(len(waters), bool); member[near[member_near]] = True
    counts = {}
    if len(positions):
        distance = cKDTree(positions).query(relative[near]+anchor)[0]
    for tolerance in tolerances:
        counts[str(float(tolerance))] = int(np.sum(member_near | (distance <= float(tolerance)+1e-9))) if len(positions) else 0
    planes = lattice_indices(relative[member]+anchor, anchor, basis, h)[:, 2]
    axial = np.asarray(reference.axial_indices, dtype=np.int64)
    lookup = {int(k): j for j, k in enumerate(axial)}; plane_counts = np.zeros(len(axial), dtype=np.int32); outside = 0
    for k in planes:
        j = lookup.get(int(k))
        if j is None:
            outside += 1
        else:
            plane_counts[j] += 1
    order = np.argsort(water_ids[member], kind='stable')
    result.update(status='measured', cavity_member_count=int(member.sum()), member_water_ids=water_ids[member][order],
                  member_positions_A=(relative[member]+anchor)[order],
                  both_count=int(np.sum(member & fixed)), fixed_only_count=int(np.sum(fixed & ~member)),
                  cavity_only_count=int(np.sum(member & ~fixed)),
                  domain_noncavity_count=int(result['domain_voxel_count']-member.sum()) if domain is not None else None,
                  tolerance_counts=counts, axial_plane_counts=plane_counts, outside_profile_count=int(outside),
                  volume_A3=float(len(index)*h**3))
    if domain is not None and np.any(member_near & ~domain):
        raise ValueError('Measured region extends outside its fixed measurement domain')
    return result


class WaterSource:
    """Stream one water-oxygen snapshot at a time; reopened after fork.

    The MDAnalysis universe is opened lazily and dropped when pickled, so each
    ``cavity-trajectory`` worker process reopens its own reader. One oxygen per
    water residue is selected; its residue index is the stable water identity.

    Parameters
    ----------
    topology, trajectory : str or Path
        MD topology and trajectory read by MDAnalysis.
    protein_selection : str
        MDAnalysis selection of the protein atoms used for the geometry; its atom
        count must equal that of ``reference_frame``.
    reference_frame : StructureFrame
        Reference protein frame of the analysis.
    frame_times : mapping of int to float
        Source frame index to time (ps) of every analysed frame; snapshots check
        the trajectory timestamps against it.
    water_selection : str, optional
        MDAnalysis selection of the water residues; default
        :data:`crevice.hydration.DEFAULT_WATER_SELECTION` (common water residue
        names such as HOH, SOL, TIP3, WAT).
    pbc : {"check", "none"}, default "check"
        ``"check"`` requires cell dimensions and applies minimum imaging;
        ``"none"`` disables imaging explicitly. ``"unwrap"`` is refused.

    Attributes
    ----------
    oxygen_indices : numpy.ndarray
        Topology indices of the selected oxygens.
    water_ids : numpy.ndarray of int64
        Residue index of each oxygen.
    metadata : dict
        Water selection, population, identity, imaging, alignment and storage
        descriptions recorded in the outputs.

    Raises
    ------
    ValueError
        For an unsupported ``pbc``, a protein atom-count mismatch or a water
        selection that does not give exactly one oxygen per residue.
    """

    def __init__(self, topology, trajectory, *, protein_selection, reference_frame, frame_times,
                 water_selection=None, pbc='check'):
        import numpy as np
        from .trajectory import _open_universe, _element_names
        from .hydration import DEFAULT_WATER_SELECTION
        if pbc not in {'check', 'none'}:
            raise ValueError('Water membership requires pbc check or none; unwrap moves protein independently of solvent')
        self.topology, self.trajectory, self.pbc = str(topology), str(trajectory), pbc
        self.frame_times = dict(frame_times)
        self._universe = _open_universe(topology, trajectory); self._pid = os.getpid()
        u = self._universe; atoms = u.atoms; elements = np.asarray(_element_names(atoms))
        protein = u.select_atoms(protein_selection)
        if len(protein) != len(reference_frame.atoms):
            raise ValueError('Water source and protein reference atom counts differ')
        self.protein_indices = protein.indices.copy()
        self.water_selection = water_selection or DEFAULT_WATER_SELECTION
        group = u.select_atoms(self.water_selection)
        oxygen = group.indices[elements[group.indices] == 'O'] if len(group) else np.empty(0, int)
        resindices = np.asarray(atoms.resindices)[oxygen]
        if len(set(resindices)) != len(resindices) or len(set(resindices)) != len(set(np.asarray(atoms.resindices)[group.indices])):
            raise ValueError('Water selection must identify exactly one oxygen per water residue')
        self.oxygen_indices = oxygen.copy(); self.water_ids = resindices.astype(np.int64)
        self.metadata = {'water_selection': self.water_selection, 'water_population': int(len(oxygen)),
            'water_identity': 'stable MD topology residue index (water-resindex); one oxygen per residue',
            'periodic_imaging': 'minimum image about the aligned reference anchor (shared with the fixed-region count); refused beyond half the minimum cell height' if pbc != 'none' else 'disabled explicitly',
            'alignment': 'same rigid fit as the frame geometry',
            'storage': 'one water snapshot per process; same source indices/times as protein'}

    def __getstate__(self):
        state = self.__dict__.copy(); state['_universe'] = None; state['_pid'] = None
        return state

    def snapshot(self, frame):
        """Oxygen coordinates and periodic cell of one analysed frame.

        Parameters
        ----------
        frame : StructureFrame
            The protein frame being analysed; its ``frame_index`` is the source
            trajectory index.

        Returns
        -------
        coords : numpy.ndarray, shape (n, 3)
            Original-frame oxygen coordinates, Å.
        box : numpy.ndarray or None
            Cell dimensions, or ``None`` with ``pbc="none"``.

        Raises
        ------
        ValueError
            If the frame was not declared, the timestamps differ, the protein
            coordinates of the two readers differ (inconsistent imaging) or the cell
            is missing with ``pbc="check"``.
        """
        import numpy as np
        from .trajectory import _open_universe
        if self._pid != os.getpid() or self._universe is None:
            self._universe = _open_universe(self.topology, self.trajectory); self._pid = os.getpid()
        index = frame.frame_index
        if index not in self.frame_times:
            raise ValueError('Water frame is absent from the declared protein analysis')
        ts = self._universe.trajectory[index]
        if not np.isclose(float(ts.time), self.frame_times[index], rtol=0, atol=1e-5):
            raise ValueError('Water and protein timestamps differ')
        coords = self._universe.atoms.positions
        if not np.allclose(coords[self.protein_indices], np.asarray([a.coord for a in frame.atoms]), atol=1e-5, rtol=0):
            raise ValueError('Solvent source and selected protein coordinates differ; use a consistently imaged input')
        box = None
        if self.pbc != 'none':
            if ts.dimensions is None:
                raise ValueError('Periodic water membership needs cell dimensions; explicitly use pbc none')
            box = ts.dimensions.copy()
        return coords[self.oxygen_indices].astype(float), box


def frame_water_membership(frame, fit, reference, component, waters, *, geometry_mode, region_margin,
                           probe_radius, max_grid_points):
    """Membership of one frame, reading its oxygens from a :class:`WaterSource`.

    Parameters
    ----------
    frame : StructureFrame
        The analysed protein frame.
    fit : dict
        The frame's rigid fit (see :func:`water_membership`).
    reference : CavityReference
        Fixed reference cavity.
    component : VoidComponent or None
        The frame's measured region (``None`` = unresolved).
    waters : WaterSource
        Oxygen source.
    geometry_mode, region_margin, probe_radius, max_grid_points
        As in :func:`water_membership`.

    Returns
    -------
    dict
        The :func:`water_membership` record, or
        ``{"status": "unavailable_no_explicit_water", "water_population": 0}``
        when the selection has no water.
    """
    if not len(waters.oxygen_indices):
        return {'status': 'unavailable_no_explicit_water', 'water_population': 0}
    coords, box = waters.snapshot(frame)
    return water_membership(coords, waters.water_ids, box, fit, reference, component, region_margin=region_margin,
                            geometry_mode=geometry_mode, probe_radius=probe_radius, max_grid_points=max_grid_points)


def _correlation(a, b):
    import numpy as np
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 3 or np.std(a[keep]) < 1e-10 or np.std(b[keep]) < 1e-10:
        return None
    return float(np.corrcoef(a[keep], b[keep])[0, 1])


def summarize_water_membership(analysis, *, confidence=.95, block_length=None, replicates=2000,
                               seed=20260925, axial_bin_A=1., hydration_region_counts=None):
    """All-frame arrays and statistics for paired fixed/instantaneous membership.

    Every per-frame series (fixed-reference count, domain count, member count,
    cross-tabulation, tolerance counts, number density = members / measured
    volume, fixed minus member count) is described with
    :func:`crevice.cavity_trajectory.describe_series`: frame quantiles and,
    when sampling allows, an approximate block-bootstrap interval for the mean
    (bounds clipped at zero except for the difference). Unresolved frames are
    NaN and withhold intervals. Also computed: volume-tertile strata
    (descriptive, no intervals), axial bins along the fixed axis (mean members,
    wet fraction, volume, density and width), Pearson correlations with
    volume, saved-frame membership persistence
    (:func:`crevice.hydration_trajectory.sampled_survival`) and exact
    bookkeeping identities, which are recorded rather than assumed.

    Parameters
    ----------
    analysis : dict
        :func:`crevice.cavity_trajectory.analyze_cavity_trajectory` result
        computed with water membership (every frame has ``water_membership``).
    confidence : float, default 0.95
        Nominal level of the mean intervals; 0 gives descriptive statistics only.
    block_length : int, optional
        Bootstrap block length in frames; default chosen from the
        autocorrelation.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260925
        Bootstrap random seed.
    axial_bin_A : float, default 1.0
        Width (Å) of the axial bins.
    hydration_region_counts : array_like, optional
        Per-frame ``region_water_count`` from the hydration pass; when given it is
        compared with the fixed-reference count (``fixed_reference_crosscheck``).

    Returns
    -------
    report : dict
        ``status``, frame count, settings, definitions, units,
        ``bookkeeping_checks``, ``statistics``, ``geometry_conditional_strata``,
        ``correlations``, ``axial_rows``, ``membership_persistence``,
        ``fixed_reference_crosscheck`` and ``limitations`` (written as
        ``PREFIX_water_membership.json``).
    arrays : dict of str to numpy.ndarray
        Per-frame arrays, member identities (``member_water_offsets`` +
        ``member_water_resindex``) and aligned member positions (written as
        ``PREFIX_water_membership_frames.npz``).

    Raises
    ------
    ValueError
        If a frame lacks membership or ``axial_bin_A`` is not positive.
    """
    import numpy as np
    from .cavity_trajectory import describe_series
    from .hydration_trajectory import sampled_survival
    frames = analysis['frames']; n = len(frames); reference = analysis['reference']; times = np.asarray(analysis['time_ps'], dtype=float)
    rows = [r.get('water_membership') for r in frames]
    if any(r is None for r in rows):
        raise ValueError('Water membership was not measured for every frame')
    if not (float(axial_bin_A) > 0 and np.isfinite(axial_bin_A)):
        raise ValueError('axial_bin_A must be positive')
    settings = analysis['settings']; mode = settings.get('geometry_mode')
    available = all(r['status'] != 'unavailable_no_explicit_water' for r in rows)
    def series(key):
        return np.asarray([r.get(key) if r.get(key) is not None else np.nan for r in rows], dtype=float)
    arrays = {'time_ps': times, 'frame_indices': np.asarray([r['frame_index'] for r in frames]),
              'volume_A3': np.asarray([r['volume_A3'] if r['volume_A3'] is not None else np.nan for r in frames], dtype=float)}
    for key in ['fixed_reference_count', 'domain_voxel_count', 'cavity_member_count', 'both_count', 'fixed_only_count',
                'cavity_only_count', 'domain_noncavity_count', 'outside_profile_count']:
        arrays[key] = series(key)
    measured = np.isfinite(arrays['cavity_member_count'])
    for tolerance in MEMBERSHIP_TOLERANCES_A:
        arrays[f'cavity_member_count_tolerance_{tolerance:g}A'] = np.asarray(
            [r['tolerance_counts'][str(float(tolerance))] if r.get('tolerance_counts') else np.nan for r in rows], dtype=float)
    volume = arrays['volume_A3']
    with np.errstate(divide='ignore', invalid='ignore'):
        arrays['cavity_number_density_per_A3'] = np.where(measured & (volume > 0), arrays['cavity_member_count']/volume, np.nan)
    arrays['fixed_minus_cavity_count'] = arrays['fixed_reference_count']-arrays['cavity_member_count']
    axial = np.asarray(reference.axial_indices); planes = np.full((n, len(axial)), np.nan)
    for i, r in enumerate(rows):
        if r.get('axial_plane_counts') is not None:
            planes[i] = r['axial_plane_counts']
    arrays['axial_plane_positions_A'] = axial*reference.spacing; arrays['axial_plane_member_counts'] = planes
    ids = [r.get('member_water_ids') for r in rows]
    arrays['member_water_offsets'] = np.cumsum([0]+[len(x) if x is not None else 0 for x in ids]).astype(np.int64)
    arrays['member_water_resindex'] = np.concatenate([np.asarray(x, dtype=np.int64) for x in ids if x is not None] or [np.empty(0, np.int64)])
    positions = [r.get('member_positions_A') for r in rows]
    arrays['member_water_positions_A'] = np.concatenate([np.asarray(x, float).reshape((-1, 3)) for x in positions if x is not None] or [np.empty((0, 3))])
    arrays['membership_measured'] = measured
    # Exact bookkeeping identities, recorded rather than assumed.
    m = measured
    checks = {'every_frame_attempted': True, 'measured_frames': int(m.sum()), 'unresolved_frames': int(np.sum([r['status'] == 'unresolved_geometry' for r in rows])),
              'zero_volume_frames': int(np.sum(m & (volume == 0))),
              'fixed_equals_both_plus_fixed_only': bool(np.all(arrays['fixed_reference_count'][m] == arrays['both_count'][m]+arrays['fixed_only_count'][m])),
              'cavity_equals_both_plus_cavity_only': bool(np.all(arrays['cavity_member_count'][m] == arrays['both_count'][m]+arrays['cavity_only_count'][m])),
              'axial_planes_plus_outside_equal_cavity': bool(np.all(np.nansum(planes[m], axis=1)+arrays['outside_profile_count'][m] == arrays['cavity_member_count'][m])),
              'member_ids_equal_cavity_count': bool(np.all(np.diff(arrays['member_water_offsets'])[m] == arrays['cavity_member_count'][m])),
              'zero_volume_has_zero_members': bool(np.all(arrays['cavity_member_count'][m & (volume == 0)] == 0))}
    if mode == 'reference-region':
        checks['domain_equals_cavity_plus_domain_noncavity'] = bool(np.all(arrays['domain_voxel_count'][m] == arrays['cavity_member_count'][m]+arrays['domain_noncavity_count'][m]))
        checks['cavity_within_domain'] = bool(np.all(arrays['domain_noncavity_count'][m] >= 0))
    regular = n < 3 or np.allclose(np.diff(times), np.diff(times)[0], rtol=1e-7, atol=1e-7)
    options = dict(confidence=confidence or .95, block_length=block_length, replicates=replicates, seed=seed, regular=regular and confidence != 0)
    def describe(x, lower=0.):
        result = describe_series(x, **options)
        if confidence == 0:
            result['confidence_interval'] = {'status': 'not_requested'}
        interval = result['confidence_interval']
        for field in ['pointwise_lower', 'pointwise_upper', 'simultaneous_lower', 'simultaneous_upper']:
            if field in interval and lower is not None:
                interval[field] = [max(lower, v) if v is not None else None for v in interval[field]]
        return result
    statistics = {}
    if available:
        for key in ['fixed_reference_count', 'domain_voxel_count', 'cavity_member_count', 'domain_noncavity_count',
                    'cavity_number_density_per_A3', 'both_count', 'fixed_only_count', 'cavity_only_count']+[
                    f'cavity_member_count_tolerance_{t:g}A' for t in MEMBERSHIP_TOLERANCES_A]:
            if np.isfinite(arrays[key]).any():
                statistics[key] = describe(arrays[key])
        statistics['fixed_minus_cavity_count'] = describe(arrays['fixed_minus_cavity_count'], lower=None)
    # Geometry-conditional strata: descriptive only, frames are not contiguous blocks.
    strata = []
    observed = np.flatnonzero(m)
    if available and len(observed) >= 3:
        edges = np.quantile(volume[observed], [0, 1/3, 2/3, 1])
        labels = np.clip(np.searchsorted(edges[1:-1], volume[observed], side='right'), 0, 2)
        for s, name in enumerate(['lower_volume_tertile', 'middle_volume_tertile', 'upper_volume_tertile']):
            chosen = observed[labels == s]
            if not len(chosen):
                continue
            density = arrays['cavity_number_density_per_A3'][chosen]
            strata.append({'stratum': name, 'frames': int(len(chosen)), 'volume_min_A3': float(volume[chosen].min()),
                'volume_max_A3': float(volume[chosen].max()), 'mean_volume_A3': float(volume[chosen].mean()),
                'mean_cavity_member_count': float(arrays['cavity_member_count'][chosen].mean()),
                'mean_fixed_reference_count': float(arrays['fixed_reference_count'][chosen].mean()),
                'mean_number_density_per_A3': float(np.nanmean(density)) if np.isfinite(density).any() else None,
                'pooled_number_density_per_A3': float(arrays['cavity_member_count'][chosen].sum()/volume[chosen].sum()) if volume[chosen].sum() > 0 else None,
                'interval': 'not estimated: strata are not contiguous time blocks'})
    # Linked axial water/width summaries on declared bins along the fixed axis.
    positions = axial*reference.spacing; bins = np.floor(positions/float(axial_bin_A)+1e-9).astype(int)
    area = np.full((n, len(axial)), np.nan); width = area.copy()
    for i, r in enumerate(frames):
        if r.get('profile') is not None and r['volume_A3'] is not None:
            area[i] = r['profile']['total_area_A2']; width[i] = r['profile']['largest_component_diameter_A']
    arrays['axial_plane_total_area_A2'] = area; arrays['axial_plane_largest_diameter_A'] = width
    axial_rows = []
    if available:
        for b in np.unique(bins):
            cols = bins == b; ok = m & np.isfinite(area[:, cols]).all(axis=1)
            if not ok.any():
                continue
            water = np.nansum(planes[:, cols], axis=1)[ok]; vol = (area[:, cols].sum(axis=1)*reference.spacing)[ok]
            wid = width[:, cols].mean(axis=1)[ok]; wet = water > 0
            if vol.sum() == 0 and water.sum() == 0:
                continue
            axial_rows.append({'bin_start_A': float(b*axial_bin_A), 'bin_end_A': float((b+1)*axial_bin_A), 'planes': int(cols.sum()),
                'measured_frames': int(ok.sum()), 'mean_member_count': float(water.mean()), 'wet_frames': int(wet.sum()),
                'wet_fraction': float(wet.mean()), 'mean_volume_A3': float(vol.mean()),
                'pooled_number_density_per_A3': float(water.sum()/vol.sum()) if vol.sum() > 0 else None,
                'mean_width_A': float(wid.mean()), 'mean_width_wet_A': float(wid[wet].mean()) if wet.any() else None,
                'mean_width_dry_A': float(wid[~wet].mean()) if (~wet).any() else None,
                'water_width_correlation': _correlation(water.astype(float), wid)})
    histories = [set(map(int, x)) if x is not None else None for x in ids]
    if not available:
        persistence = {'status': 'unavailable_no_explicit_water', 'rows': []}
    elif any(x is None for x in histories):
        persistence = {'status': 'unavailable_missing_frames', 'rows': []}
    else:
        persistence = sampled_survival(histories, times)
    crosscheck = {'status': 'not_available'}
    if hydration_region_counts is not None:
        other = np.asarray(hydration_region_counts, dtype=float)
        if other.shape != arrays['fixed_reference_count'].shape:
            crosscheck = {'status': 'frame_count_differs'}
        else:
            both = np.isfinite(other) & np.isfinite(arrays['fixed_reference_count'])
            difference = np.abs(other[both]-arrays['fixed_reference_count'][both])
            crosscheck = {'status': 'identical' if both.any() and not difference.any() else 'differs',
                          'frames_compared': int(both.sum()), 'max_abs_difference': float(difference.max()) if both.any() else None,
                          'definition': 'Per-frame comparison with the unchanged hydration-pass region_water_count'}
    report = {'status': 'measured' if available else 'unavailable_no_explicit_water', 'frame_count': n,
        'geometry_mode': mode, 'time_range_ps': [float(times[0]), float(times[-1])], 'regular_timestamps': bool(regular),
        'settings': {'region_margin_A': settings.get('region_margin'), 'spacing_A': reference.spacing,
            'probe_radius_A': (settings.get('geometry') or {}).get('probe_radius'), 'axial_bin_A': float(axial_bin_A),
            'tolerance_sensitivity_A': list(MEMBERSHIP_TOLERANCES_A), 'waters': settings.get('waters'),
            'alignment_residues': settings.get('alignment_residues'),
            'region_definition': (settings.get('region_definition') or {}).get('region_id')},
        'definitions': {
            'cavity_member_count': 'Aligned, minimum-imaged water oxygens whose centre lies in an occupied voxel of the region measured in that frame; voxel k owns [k-1/2,k+1/2) of each fixed lattice coordinate (spacing h)',
            'fixed_reference_count': 'Existing definition: aligned oxygens within the region margin of any fixed occupied reference sample (identical function to hydration region_water_count)',
            'domain_voxel_count': 'Oxygens whose voxel lies in the fixed discretised reference-region measurement domain (reference-region mode only)',
            'domain_noncavity_count': 'domain_voxel_count minus cavity_member_count: oxygens in the fixed domain but outside the measured region',
            'both/fixed_only/cavity_only': 'Cross-tabulation of fixed-reference and instantaneous membership for the same oxygens',
            'tolerance_counts': 'Sensitivity only: voxel rule OR within the stated distance of an occupied sample centre',
            'cavity_number_density_per_A3': 'cavity_member_count / measured volume; undefined when volume is zero or geometry unresolved',
            'unresolved': 'Missing (NaN), never zero; mean intervals of series with missing frames are withheld',
            'axial': 'Members binned by fixed-axis lattice plane (position relative to the reference anchor); widths are largest planar-component area-equivalent diameters'},
        'units': {'counts': 'water oxygens', 'volume': 'Å³', 'density': 'oxygens per Å³', 'positions': 'Å', 'widths': 'Å'},
        'bookkeeping_checks': checks, 'statistics': statistics, 'geometry_conditional_strata': strata,
        'correlations': {'volume_vs_cavity_member_count': _correlation(volume, arrays['cavity_member_count']),
                         'volume_vs_fixed_reference_count': _correlation(volume, arrays['fixed_reference_count']),
                         'interpretation': 'Descriptive Pearson association on measured frames; no p-value or causal interpretation; count and volume share geometry'},
        'axial_rows': axial_rows, 'membership_persistence': persistence,
        'fixed_reference_crosscheck': crosscheck, 'limitations': LIMITATIONS}
    return report, arrays


@annotate_option
def write_water_membership_bundle(analysis, output_dir, *, prefix='crevice', confidence=.95, block_length=None,
                                  replicates=2000, seed=20260925, dpi=240, axial_bin_A=1., hydration_region_counts=None,
                                  annotate=None):
    """Write instantaneous-cavity water-membership statistics, tables and a figure.

    A water oxygen is a member in a frame if its centre lies in an occupied
    voxel of that frame's measured region on the fixed reference lattice (see
    :func:`water_membership`). Members are paired per frame with the
    fixed-reference neighbourhood count; both series, their overlap, number
    density, axial distribution and volume-stratified statistics are summarised
    by :func:`summarize_water_membership`. Unresolved frames are missing, not zero.

    Parameters
    ----------
    analysis : dict
        Cavity-trajectory analysis computed with water membership.
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    confidence : float, default 0.95
        Nominal level of mean intervals.
    block_length : int, optional
        Bootstrap block length in frames.
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260925
        Bootstrap random seed.
    dpi : int, default 240
        Figure resolution.
    axial_bin_A : float, default 1.0
        Axial bin width (Å) of the axial table and figure panel.
    hydration_region_counts : array_like, optional
        Fixed-region counts from the hydration pass, cross-checked per frame.
    annotate : bool, optional
        Draw the panel title and figure title (interpretation reminders).
        ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_water_membership.json : JSON
        Definitions, statistics, strata, correlations, persistence, cross-check.
    PREFIX_water_membership_summary.csv : table
        The JSON's interval statistics (one row per count or density series)
        and correlations, flattened as in :func:`crevice.io.summary_statistics_rows`.
    PREFIX_water_membership_strata.csv : table
        Volume-tertile strata: ``stratum``, ``frames``, volume range and mean
        (Å³), mean member and fixed-reference counts, mean and pooled number
        density (Å⁻³). An empty file when there are no strata.
    PREFIX_water_membership_frames.csv, PREFIX_water_membership_frames.npz : per frame
        Counts, density, member IDs and aligned member-oxygen positions (Å).
    PREFIX_water_membership_axial.csv : table
        Per axial bin: mean members, wet fraction, volume, density, width.
    PREFIX_water_membership.png : figure
        Three panels: "Water oxygens" against "Time (ns)" for both series;
        "Members | measured" against "Measured volume (Å³)"; mean members per bin
        ("Water oxygens") and mean "Area-equivalent width (Å)" against "Fixed
        reference-axis position (Å)". Measured runs only.

    Colours and representation
    --------------------------
    * Time panel: blue-grey (``#547c9b``) fixed-reference count, teal
      (``#178c91``) member count, with a two-entry legend.
    * Scatter: translucent teal points.
    * Axial panel: translucent teal bars (members) and an orange (``#d98c33``)
      width line on the right axis.
    * In the viewer scenes, members of each displayed snapshot are orange
      0.7 Å spheres in their own group ``crevice_member_waters``.
    """
    import csv
    import numpy as np
    from .figures import _pyplot, _save_figure
    root = Path(output_dir); root.mkdir(parents=True, exist_ok=True)
    if not prefix or Path(prefix).name != prefix:
        raise ValueError('prefix must be a simple filename stem')
    report, arrays = summarize_water_membership(analysis, confidence=confidence, block_length=block_length, replicates=replicates,
                                                seed=seed, axial_bin_A=axial_bin_A, hydration_region_counts=hydration_region_counts)
    files = {}
    def path(key, suffix):
        target = root/(prefix+suffix); files[key] = str(target); return target
    def table(target, rows):
        if not rows:
            target.write_text(''); return
        with target.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    path('water_membership_json', '_water_membership.json').write_text(json.dumps(report, indent=2, allow_nan=False, ensure_ascii=False)+'\n', encoding='utf-8')
    from .io import write_summary_statistics_csv
    write_summary_statistics_csv(report, path('water_membership_summary_csv', '_water_membership_summary.csv'))
    table(path('water_membership_strata_csv', '_water_membership_strata.csv'), report.get('geometry_conditional_strata') or [])
    np.savez_compressed(path('water_membership_npz', '_water_membership_frames.npz'), **arrays)
    def value(x):
        return None if not np.isfinite(x) else (int(x) if float(x).is_integer() else float(x))
    keys = ['volume_A3', 'fixed_reference_count', 'domain_voxel_count', 'cavity_member_count', 'domain_noncavity_count',
            'both_count', 'fixed_only_count', 'cavity_only_count', 'outside_profile_count', 'cavity_number_density_per_A3']+[
            f'cavity_member_count_tolerance_{t:g}A' for t in MEMBERSHIP_TOLERANCES_A]
    rows = [{'frame_index': int(arrays['frame_indices'][i]), 'time_ps': float(arrays['time_ps'][i]),
             'status': (analysis['frames'][i]['water_membership'] or {}).get('status'), **{k: value(arrays[k][i]) for k in keys}}
            for i in range(len(arrays['time_ps']))]
    table(path('water_membership_frames_csv', '_water_membership_frames.csv'), rows)
    table(path('water_membership_axial_csv', '_water_membership_axial.csv'), report['axial_rows'])
    if report['status'] == 'measured':
        plt = _pyplot(); t = arrays['time_ps']/1000
        text = FigureText("Fixed-reference and instantaneous-region water counts against time, members against measured volume, "
                          "and mean members and width along the fixed reference axis")
        fig, axes = plt.subplots(3, 1, figsize=(10, 8.5), sharex=False, layout='constrained')
        axes[0].plot(t, arrays['fixed_reference_count'], color='#547c9b', lw=.8, label='Fixed reference neighbourhood')
        axes[0].plot(t, arrays['cavity_member_count'], color='#178c91', lw=.8, label="In this frame's measured region")
        axes[0].set_ylabel('Water oxygens'); axes[0].set_xlabel('Time (ns)'); axes[0].legend(frameon=False, fontsize=8)
        text.title(axes[0], 'Two separate observables · unresolved frames are gaps, not zero', fontsize=10)
        ok = np.isfinite(arrays['cavity_member_count'])
        axes[1].scatter(arrays['volume_A3'][ok], arrays['cavity_member_count'][ok], s=5, color='#178c91', alpha=.5)
        axes[1].set_xlabel('Measured volume (Å³)'); axes[1].set_ylabel('Members | measured')
        if report['axial_rows']:
            x = [(r['bin_start_A']+r['bin_end_A'])/2 for r in report['axial_rows']]
            axes[2].bar(x, [r['mean_member_count'] for r in report['axial_rows']], width=.9*axial_bin_A, color='#178c91', alpha=.7, label='Mean members per bin')
            twin = axes[2].twinx(); twin.plot(x, [r['mean_width_A'] for r in report['axial_rows']], color='#d98c33', lw=1.4, label='Mean width')
            twin.set_ylabel('Area-equivalent width (Å)'); twin.set_ylim(bottom=0)
        axes[2].set_xlabel('Fixed reference-axis position (Å)'); axes[2].set_ylabel('Water oxygens')
        for ax in axes:
            ax.spines[['top', 'right']].set_visible(False)
        text.suptitle(fig, 'Geometric membership, not binding, permeation or functional state', fontsize=10)
        _save_figure(fig, path('water_membership_png', '_water_membership.png'), dpi=dpi, text=text)
    return files


def viewer_member_waters(analysis):
    """Aligned member-oxygen coordinates per source frame, for molecular viewer groups.

    Parameters
    ----------
    analysis : dict
        Cavity-trajectory analysis computed with water membership.

    Returns
    -------
    dict of int to dict
        For each measured frame index: ``water_ids`` (sorted) and
        ``positions_A`` (aligned oxygen centres, Å). The viewer scenes draw
        them as orange 0.7 Å spheres in the group ``crevice_member_waters``.
        Frames with unresolved geometry or no water are omitted.
    """
    frames = {}
    for row in analysis['frames']:
        record = row.get('water_membership') or {}
        if record.get('member_water_ids') is None:
            continue
        frames[int(row['frame_index'])] = {'water_ids': [int(x) for x in record['member_water_ids']],
                                          'positions_A': [list(map(float, p)) for p in record['member_positions_A']]}
    return frames

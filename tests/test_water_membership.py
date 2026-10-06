"""Instantaneous-region water membership: analytic, conservation, rigid/PBC and
moving-boundary controls. Software/synthetic validation only; no biological claim."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pytest
from crevice.models import ChannelPoint, VoidComponent
from crevice.cavity_trajectory import CavityReference, cavity_section_profile, analyze_cavity_trajectory
from crevice.water_membership import (water_membership, summarize_water_membership, write_water_membership_bundle,
                                    WaterSource, lattice_indices)

IDENTITY = {'rotation_rows': np.eye(3).tolist(), 'translation_A': [0., 0., 0.]}


def block(lo, hi, spacing=1.):
    r = np.arange(lo, hi+1)*spacing
    return np.array([(x, y, z) for x in r for y in r for z in r], dtype=float)


def reference_for(points, spacing=1.):
    points = np.asarray(points, float)
    axial = np.arange(int(np.floor(points[:, 2].min()/spacing))-2, int(np.ceil(points[:, 2].max()/spacing))+3)
    return CavityReference(points, np.eye(3), np.zeros(3), spacing, axial, None, None, None)


def component(xyz):
    xyz = np.asarray(xyz, float).reshape((-1, 3))
    points = tuple(ChannelPoint(i, tuple(map(float, p)), 1., 1., float(p[2])) for i, p in enumerate(xyz))
    return VoidComponent(1, 'test', (0, 0, 0), float(len(points)), points)


def test_voxel_rule_ties_crosstab_and_tolerance_sensitivity_are_analytic():
    ref = reference_for(block(-1, 1))
    waters = np.array([[0, 0, 0], [.49, -.49, .2], [1.49, 0, 0], [1.5, 0, 0], [-1.5, 0, 0], [-1.51, 0, 0], [0, 0, 5]], float)
    got = water_membership(waters, np.arange(10, 17), None, IDENTITY, ref, component(block(-1, 1)),
                           region_margin=.5, geometry_mode='rolling')
    # Half-up ties: +1.5 belongs to voxel 2 (outside), -1.5 to voxel -1 (inside).
    assert got['member_water_ids'].tolist() == [10, 11, 12, 14]
    assert got['cavity_member_count'] == 4 and got['volume_A3'] == 27
    # Fixed rule: distance <= 0.5 A to a reference sample -> waters 0, 2, 3, 4.
    assert got['fixed_reference_count'] == 4
    assert (got['both_count'], got['fixed_only_count'], got['cavity_only_count']) == (3, 1, 1)
    assert got['tolerance_counts'] == {'0.5': 5, '1.0': 6}
    plane = list(ref.axial_indices).index(0)
    assert got['axial_plane_counts'][plane] == 4 and got['axial_plane_counts'].sum() == 4
    assert got['domain_voxel_count'] is None
    assert lattice_indices([[.5, -.5, 0]], np.zeros(3), np.eye(3), 1.).tolist() == [[1, 0, 0]]


def test_sphere_counts_match_independent_voxel_enumeration_and_volume():
    h = .25; radius = 3.
    grid = block(-14, 14, h); occupied = grid[np.linalg.norm(grid, axis=1) <= radius]
    ref = reference_for(occupied, h)
    rng = np.random.default_rng(25); waters = rng.uniform(-4, 4, size=(40000, 3))
    got = water_membership(waters, np.arange(len(waters)), None, IDENTITY, ref, component(occupied),
                           region_margin=0., geometry_mode='rolling', tolerances=())
    index = np.floor(waters/h+.5)
    expected = np.flatnonzero(np.linalg.norm(index*h, axis=1) <= radius+1e-12)
    assert got['member_water_ids'].tolist() == expected.tolist()
    # Uniform samples: member fraction estimates the reported voxel volume.
    fraction = len(expected)/len(waters); volume_fraction = len(occupied)*h**3/8**3
    assert abs(fraction-volume_fraction) < 4*np.sqrt(volume_fraction*(1-volume_fraction)/len(waters))
    assert got['volume_A3'] == pytest.approx(len(occupied)*h**3)


@pytest.mark.parametrize('box', [[40, 42, 44, 90, 90, 90], [40, 42, 44, 75, 80, 65]])
def test_rigid_transform_and_periodic_image_shifts_leave_membership_unchanged(box):
    from scipy.spatial.transform import Rotation
    from MDAnalysis.lib.mdamath import triclinic_vectors
    points = block(-2, 2, .5); ref = reference_for(points, .5)
    from scipy.spatial import cKDTree
    # Voxel-interior positions away from every distance threshold, so the
    # control tests transforms rather than floating-point face ties.
    rng = np.random.default_rng(7)
    aligned = (rng.integers(-8, 9, size=(1200, 3))+rng.uniform(-.4, .4, size=(1200, 3)))*.5
    distance = cKDTree(points).query(aligned)[0]
    aligned = aligned[np.min(np.abs(distance[:, None]-np.array([.5, 1., 2.])), axis=1) > 1e-3][:400]
    ids = np.arange(len(aligned))
    kwargs = dict(region_margin=2., geometry_mode='reference-region', probe_radius=1.4)
    base = water_membership(aligned, ids, np.asarray(box, float), IDENTITY, ref, component(points), **kwargs)
    assert base['cavity_member_count'] > 0 and base['fixed_only_count'] > 0
    R = Rotation.from_rotvec([.4, -.9, .3]).as_matrix(); t = np.array([31.2, -7.5, 12.])
    cell = triclinic_vectors(box).astype(float)
    original = (aligned-t)@R.T+rng.integers(-2, 3, size=(len(aligned), 3))@cell
    moved = water_membership(original, ids, np.asarray(box, float), {'rotation_rows': R.tolist(), 'translation_A': t.tolist()},
                             ref, component(points), **kwargs)
    keys = ['fixed_reference_count', 'domain_voxel_count', 'cavity_member_count', 'both_count', 'fixed_only_count',
            'cavity_only_count', 'domain_noncavity_count', 'tolerance_counts']
    assert {k: moved[k] for k in keys} == {k: base[k] for k in keys}
    assert moved['member_water_ids'].tolist() == base['member_water_ids'].tolist()
    assert moved['axial_plane_counts'].tolist() == base['axial_plane_counts'].tolist()
    # Nonperiodic input without image shifts agrees too.
    plain = water_membership(aligned, ids, None, IDENTITY, ref, component(points), **kwargs)
    assert plain['member_water_ids'].tolist() == base['member_water_ids'].tolist()
    # A region reaching half the minimum cell height is refused, not guessed.
    with pytest.raises(ValueError, match='half the minimum periodic cell'):
        water_membership(aligned, ids, np.array([9, 9, 9, 90, 90, 90.]), IDENTITY, ref, component(points), **kwargs)


def test_moving_boundary_and_moving_water_separate_the_two_observables():
    points = block(-2, 2, .5); ref = reference_for(points, .5)
    kwargs = dict(region_margin=2., geometry_mode='reference-region', probe_radius=1.4, tolerances=())
    waters = np.array([[0, 0, -.8], [0, 0, .8], [.3, -.4, .2], [0, 0, 2.], [5, 5, 5]], float); ids = np.arange(5)
    full = water_membership(waters, ids, None, IDENTITY, ref, component(points), **kwargs)
    upper = points[points[:, 2] >= 0]
    shrunk = water_membership(waters, ids, None, IDENTITY, ref, component(upper), **kwargs)
    # Cavity changes, water does not: fixed count identical, instantaneous count follows the cavity.
    assert shrunk['fixed_reference_count'] == full['fixed_reference_count'] == 4
    assert shrunk['domain_voxel_count'] == full['domain_voxel_count'] == 4
    assert full['member_water_ids'].tolist() == [0, 1, 2]
    assert shrunk['member_water_ids'].tolist() == [1, 2]
    assert shrunk['domain_noncavity_count'] == full['domain_noncavity_count']+1
    # Water leaves both definitions: both counts drop together.
    away = waters.copy(); away[1] = [30, 30, 30]
    left = water_membership(away, ids, None, IDENTITY, ref, component(points), **kwargs)
    assert (left['fixed_reference_count'], left['cavity_member_count']) == (3, 2)
    # Water moves from the measured region to the fixed domain outside it:
    # fixed-reference count unchanged, instantaneous count decreases.
    inside = waters.copy(); inside[0] = [0, 0, 1.9]
    moved = water_membership(inside, ids, None, IDENTITY, ref, component(points), **kwargs)
    assert moved['fixed_reference_count'] == full['fixed_reference_count']
    assert moved['cavity_member_count'] == full['cavity_member_count']-1
    for row in [full, shrunk, left, moved]:
        assert row['fixed_reference_count'] == row['both_count']+row['fixed_only_count']
        assert row['cavity_member_count'] == row['both_count']+row['cavity_only_count']
        assert row['domain_voxel_count'] == row['cavity_member_count']+row['domain_noncavity_count']
    with pytest.raises(ValueError, match='outside its fixed measurement domain'):
        water_membership(np.array([[0, 0, 4.]]), [0], None, IDENTITY, ref, component(np.array([[0, 0, 4.]])), **kwargs)


def synthetic_analysis(rows, ref, volumes, profiles, times):
    frames = [{'frame_index': i, 'volume_A3': v, 'profile': p, 'water_membership': r, 'status': r['status']}
              for i, (r, v, p) in enumerate(zip(rows, volumes, profiles))]
    return {'frames': frames, 'reference': ref, 'time_ps': np.asarray(times, float),
            'settings': {'geometry_mode': 'rolling', 'region_margin': .5, 'geometry': {'probe_radius': 1.4}}}


def test_zero_volume_unresolved_and_absent_water_frames_are_explicit(tmp_path):
    ref = reference_for(block(-1, 1)); waters = np.array([[0, 0, 0], [0, 0, 1.2], [9, 9, 9]], float); ids = np.arange(3)
    measured = component(block(-1, 1)); empty = component(np.empty((0, 3)))
    kw = dict(region_margin=.5, geometry_mode='rolling')
    rows = [water_membership(waters, ids, None, IDENTITY, ref, measured, **kw),
            water_membership(waters, ids, None, IDENTITY, ref, empty, **kw),
            water_membership(waters, ids, None, IDENTITY, ref, None, **kw)]
    zeros = {k: np.zeros(len(ref.axial_indices)) for k in ['total_area_A2', 'largest_component_diameter_A']}
    analysis = synthetic_analysis(rows, ref, [27., 0., None], [cavity_section_profile(measured, ref), zeros, None], [0, 100, 200])
    report, arrays = summarize_water_membership(analysis, replicates=200)
    assert arrays['cavity_member_count'][:2].tolist() == [2, 0] and np.isnan(arrays['cavity_member_count'][2])
    assert arrays['fixed_reference_count'].tolist() == [2, 2, 2]
    assert arrays['cavity_number_density_per_A3'][0] == pytest.approx(2/27)
    assert np.isnan(arrays['cavity_number_density_per_A3'][1:]).all()
    stat = report['statistics']['cavity_member_count']
    assert stat['observed_frames'] == 2 and stat['confidence_interval']['status'] == 'unavailable_missing_frames'
    assert report['statistics']['fixed_reference_count']['observed_frames'] == 3
    assert report['membership_persistence']['status'] == 'unavailable_missing_frames'
    checks = report['bookkeeping_checks']
    assert checks['unresolved_frames'] == 1 and checks['zero_volume_frames'] == 1
    assert all(v for k, v in checks.items() if isinstance(v, bool))
    assert arrays['member_water_offsets'].tolist() == [0, 2, 2, 2]
    json.dumps(report, allow_nan=False)
    files = write_water_membership_bundle(analysis, tmp_path, replicates=200, dpi=50)
    assert all(Path(p).is_file() for p in files.values())
    summary = Path(files['water_membership_summary_csv']).read_text().splitlines()
    assert summary[0].startswith('quantity,item,value,') and len(summary) - 1 >= len(report['statistics'])
    assert Path(files['water_membership_strata_csv']).is_file()
    absent = [{'status': 'unavailable_no_explicit_water', 'water_population': 0}]*3
    report, _ = summarize_water_membership(synthetic_analysis(absent, ref, [27., 0., None], [None]*3, [0, 100, 200]), replicates=200)
    assert report['status'] == 'unavailable_no_explicit_water' and report['statistics'] == {}


class StreamedWaters:
    """Duck-typed WaterSource: per-frame original coordinates and a periodic cell."""
    def __init__(self, coords, box=None):
        self.coords = coords; self.box = box; n = len(coords[0])
        self.oxygen_indices = np.arange(n); self.water_ids = np.arange(100, 100+n)
        self.metadata = {'water_selection': 'synthetic', 'water_population': n}

    def snapshot(self, frame):
        return np.asarray(self.coords[frame.frame_index], float), self.box


def test_end_to_end_moving_boundary_rigid_motion_and_periodic_shift(tmp_path):
    from crevice.models import Atom
    from crevice.trajectory import Trajectory
    from crevice.volume_export import write_void_cast_dx
    from crevice.rolling import rolling_probe_cast
    from test_channel_coordinate import rectangular_channel
    channel = rectangular_channel(half_x=4, half_y=4, half_length=4, capped=True)
    cast = rolling_probe_cast(channel, spacing=.5, probe_radius=.8, selection_mode='all', enclosure_fraction=0, min_component_volume=1,
                              grid_basis=np.eye(3), grid_anchor=(0, 0, 0))
    main = max(cast.components, key=lambda c: c.volume)
    ref = CavityReference.from_dx(write_void_cast_dx(replace(cast, components=(main,), points=main.points), tmp_path/'ref.dx'), axis=(0, 0, 1))
    scaffold = sorted({a.residue_key.label for a in channel.atoms})
    plug = Atom(9000, 'C', 'LEU', 'A', 9000, 0., 0., 40., 'C')
    def frame(index, atoms):
        return replace(channel, atoms=tuple(atoms), frame_index=index)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.]]); shift = np.array([11.3, -4.2, 7.9])
    def move(atom, p):
        return replace(atom, x=float(p[0]), y=float(p[1]), z=float(p[2]))
    base_atoms = channel.atoms+(plug,)
    moved_atoms = tuple(move(a, np.asarray(a.coord)@rotation.T+shift) for a in base_atoms)
    plugged = channel.atoms+(replace(plug, z=0.),)
    frames = [frame(0, base_atoms), frame(1, moved_atoms), frame(2, plugged), frame(3, base_atoms)]
    waters = np.array([[0, 0, 0.1], [.4, .3, -1.1], [-.3, .2, 1.3], [20, 20, 20]], float)
    box = np.array([60, 60, 60, 90, 90, 90.])
    departed = waters.copy(); departed[0] = [0, 25, 0]
    coords = [waters, waters@rotation.T+shift+np.array([60, 0, -60]), waters, departed]
    analysis = analyze_cavity_trajectory(Trajectory(frames, time_step=100.), ref, geometry_mode='reference-region',
                                         geometry={'probe_radius': 1.4}, alignment_residues=scaffold,
                                         waters=StreamedWaters(coords, box))
    rows = [f['water_membership'] for f in analysis['frames']]
    assert rows[0]['cavity_member_count'] >= 2 and 100 in rows[0]['member_water_ids']
    # Rigid motion plus a periodic image shift: identical observations.
    assert rows[1]['member_water_ids'].tolist() == rows[0]['member_water_ids'].tolist()
    assert rows[1]['fixed_reference_count'] == rows[0]['fixed_reference_count']
    assert analysis['frames'][1]['volume_A3'] == analysis['frames'][0]['volume_A3']
    # Boundary moves into the fixed region, waters unchanged.
    assert analysis['frames'][2]['volume_A3'] < analysis['frames'][0]['volume_A3']
    assert rows[2]['fixed_reference_count'] == rows[0]['fixed_reference_count']
    assert 100 not in rows[2]['member_water_ids'] and rows[2]['cavity_member_count'] < rows[0]['cavity_member_count']
    # Water leaves, boundary unchanged: both observables lose it.
    assert rows[3]['fixed_reference_count'] == rows[0]['fixed_reference_count']-1
    assert rows[3]['cavity_member_count'] == rows[0]['cavity_member_count']-1
    report, arrays = summarize_water_membership(analysis, replicates=200)
    assert all(v for v in report['bookkeeping_checks'].values() if isinstance(v, bool))
    assert report['settings']['waters']['water_population'] == 4
    json.dumps(report, allow_nan=False)


def test_cli_membership_is_additive_and_matches_unchanged_fixed_count(tmp_path):
    from crevice.cli import main
    from test_hydration import md_fixture, reference_dx
    source, trajectory = md_fixture(tmp_path); dx = reference_dx(tmp_path/'ref.dx')
    common = ['cavity-trajectory', str(source), str(trajectory), '--reference-volume-dx', str(dx), '--geometry-mode', 'reference-region',
              '--max-frames', '3', '--dpi', '50', '--bootstrap-replicates', '200', '--hydration-sasa-points', '32',
              '--skip-water-density', '--skip-interactions']
    assert main(common+['--out-dir', str(tmp_path/'plain')]) == 0
    assert main(common+['--out-dir', str(tmp_path/'member'), '--water-membership', '--workers', '2']) == 0
    assert not list((tmp_path/'plain').glob('*water_membership*'))
    report = json.loads((tmp_path/'member/crevice_water_membership.json').read_text())
    assert report['fixed_reference_crosscheck']['status'] == 'identical'
    with np.load(tmp_path/'member/crevice_water_membership_frames.npz') as a:
        assert a['fixed_reference_count'].tolist() == [1, 0, 1]
        assert a['cavity_member_count'][0] == a['cavity_member_count'][2] and a['cavity_member_count'][1] == 0
    with np.load(tmp_path/'plain/crevice_cavity_frames.npz') as a, np.load(tmp_path/'member/crevice_cavity_frames.npz') as b:
        assert all(np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == 'f') for k in a.files)
    with np.load(tmp_path/'plain/crevice_hydration_frames.npz') as a, np.load(tmp_path/'member/crevice_hydration_frames.npz') as b:
        assert np.array_equal(a['region_water_count'], b['region_water_count'])
    # Viewer group: own object, on by default, one entry per snapshot, coordinates from the membership pass.
    plain_scene = json.loads((tmp_path/'plain/crevice_analysis_scene.json').read_text())
    scene = json.loads((tmp_path/'member/crevice_analysis_scene.json').read_text())
    assert plain_scene['member_waters'] is None
    group = scene['member_waters']
    assert group['object'] == 'crevice_member_waters' and group['visible_by_default'] is True
    assert [e['source_frame_index'] for e in group['snapshots']] == [s['source_frame_index'] for s in scene['snapshots']]
    with np.load(tmp_path/'member/crevice_water_membership_frames.npz') as a:
        offsets = a['member_water_offsets']; positions = a['member_water_positions_A']; frames = a['frame_indices'].tolist()
    for entry in group['snapshots']:
        row = frames.index(entry['source_frame_index']); expected = positions[offsets[row]:offsets[row+1]]
        assert entry['count'] == len(expected)
        if entry['pdb']:
            lines = [l for l in (tmp_path/'member'/entry['pdb']).read_text().splitlines() if l.startswith('HETATM')]
            got = np.array([[float(l[k:k+8]) for k in (30, 38, 46)] for l in lines])
            assert np.allclose(got, expected, atol=6e-4)
        else:
            assert entry['count'] == 0
    for script in ['crevice_analysis.pml', 'crevice_analysis.vmd', 'crevice_analysis_chimerax.py']:
        assert 'crevice_member_waters' in (tmp_path/'member'/script).read_text()
    manifest = json.loads((tmp_path/'member/crevice_cavity_manifest.json').read_text())
    assert 'water_membership_json' in manifest['files']
    with pytest.raises(ValueError, match='unwrap'):
        WaterSource(source, trajectory, protein_selection='protein', reference_frame=None, frame_times={}, pbc='unwrap')

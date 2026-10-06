"""``--smooth X`` display surfaces: constraints, option handling and measurement invariance."""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
import json
import math
from pathlib import Path

import numpy as np
import pytest

from crevice import pore_profile, void_cast
from crevice.presentation import extend_channel_cast, write_viewer_structure
from crevice.volume_export import (DEFAULT_SMOOTH_A, SMOOTH_METHOD, DisplaySmoothing, cast_grid,
                                 constrained_display_field, resolve_display_smoothing,
                                 validate_display_width, write_volume_viewer_bundle)
from test_channel_coordinate import rectangular_channel


def two_rooms_with_neck():
    grid = np.zeros((19, 19, 19), dtype=np.float32)
    grid[2:8, 5:14, 5:14] = 1
    grid[11:17, 5:14, 5:14] = 1
    grid[7:12, 9, 9] = 1   # one-sample neck
    grid[4, 8, 8] = 0      # one-sample empty pocket
    return grid


def diagonal_node_bounds(field_shape, grid, k):
    """Independent (loop) reference for the cell/face/edge agreement clamp."""
    lo = np.empty(field_shape, dtype=np.float32)
    hi = np.empty(field_shape, dtype=np.float32)
    for index in np.ndindex(*field_shape):
        corners = [(i//k, i//k + (i % k != 0)) for i in index]
        values = grid[corners[0][0]:corners[0][1]+1, corners[1][0]:corners[1][1]+1, corners[2][0]:corners[2][1]+1]
        lo[index], hi[index] = values.min(), values.max()
    return lo, hi


@pytest.mark.parametrize('width', [.2, .3, .6, 1.2])
@pytest.mark.parametrize('k', [1, 2, 3])
def test_constrained_field_keeps_samples_necks_and_pockets(width, k):
    from scipy.ndimage import label
    grid = two_rooms_with_neck()
    field = constrained_display_field(grid, .5, width, k)
    assert field.shape == tuple(k*(n-1)+1 for n in grid.shape)
    assert field.dtype == np.float32
    measured = field[::k, ::k, ::k]
    assert np.array_equal(measured > .5, grid > .5)
    assert label(field > .5)[1] == 1  # neck remains connected at display resolution
    assert measured[4, 8, 8] < .5 and measured[9, 9, 9] > .5
    # The pocket stays a closed empty region inside the room.
    empty, count = label(field < .5)
    assert count == 2 and empty[4*k, 8*k, 8*k] != empty[0, 0, 0]
    assert np.any((field > 0) & (field < 1))


def test_constrained_field_clamps_every_agreeing_cell_node():
    rng = np.random.default_rng(7)
    grid = (rng.random((6, 5, 7)) > .55).astype(np.float32)
    grid[0], grid[-1], grid[:, 0], grid[:, -1], grid[:, :, 0], grid[:, :, -1] = 0, 0, 0, 0, 0, 0
    k = 3
    field = constrained_display_field(grid, .25, .6, k)
    lo, hi = diagonal_node_bounds(field.shape, grid, k)
    assert np.all(field[lo > .5] >= .501)
    assert np.all(field[hi < .5] <= .499)


def test_display_surface_stays_within_boundary_cells():
    from scipy.spatial import cKDTree
    from skimage.measure import marching_cubes
    grid = two_rooms_with_neck()
    h, k = .5, 2
    field = constrained_display_field(grid, h, 1.2, k)
    vertices = marching_cubes(field, .5, spacing=(h/k,)*3)[0]
    occupied = np.argwhere(grid > .5)*h
    empty = np.argwhere(grid < .5)*h
    bound = h*math.sqrt(3) + 1e-6
    assert cKDTree(occupied).query(vertices)[0].max() <= bound
    assert cKDTree(empty).query(vertices)[0].max() <= bound


def test_zero_width_and_invalid_values():
    grid = two_rooms_with_neck()
    assert np.array_equal(constrained_display_field(grid, .5, 0, 2), grid)
    raw = resolve_display_smoothing(smooth=0)
    field, spacing, factor = raw.field(grid, .5)
    assert np.array_equal(field, grid) and spacing == .5 and factor == 1
    assert raw.metadata(.5, grid.shape)['method'] == 'none'
    for bad in (-1, -1e-9, float('nan'), float('inf'), '-0.5', 'wide', None, True):
        with pytest.raises(ValueError, match='smooth'):
            if bad is None:
                validate_display_width(bad)
            else:
                resolve_display_smoothing(smooth=bad)
        if bad is not None:
            with pytest.raises(ValueError, match='surface_smoothing'):
                resolve_display_smoothing(surface_smoothing=bad)
    with pytest.raises(ValueError, match='binary'):
        constrained_display_field(grid*.7, .5, .3)
    with pytest.raises(ValueError, match='smooth_supersample'):
        resolve_display_smoothing(smooth=.3, smooth_supersample=0)


def test_option_resolution_and_conflict():
    default = resolve_display_smoothing()
    assert (default.method, default.width_A, default.option) == (SMOOTH_METHOD, DEFAULT_SMOOTH_A, 'default')
    chosen = resolve_display_smoothing(smooth=.9)
    assert (chosen.method, chosen.width_A, chosen.option) == (SMOOTH_METHOD, .9, 'smooth')
    legacy = resolve_display_smoothing(surface_smoothing=.6)
    assert (legacy.method, legacy.width_A) == ('legacy_clamped_gaussian', .6)
    assert resolve_display_smoothing(surface_smoothing=0).method == 'none'
    with pytest.raises(ValueError, match='not both'):
        resolve_display_smoothing(smooth=.3, surface_smoothing=.6)
    with pytest.raises(ValueError, match='not both'):
        resolve_display_smoothing(smooth=0, surface_smoothing=0)
    # Legacy output is exactly the previous sample-clamped presentation field.
    from crevice.volume_export import presentation_grid
    grid = two_rooms_with_neck()
    assert np.array_equal(legacy.field(grid, .5)[0], presentation_grid(grid, .5, .6))
    # Automatic supersampling respects the display-point budget.
    small = DisplaySmoothing(SMOOTH_METHOD, .3, 'smooth', max_display_points=40**3)
    assert small.factor((30, 30, 30), .5) == 1
    assert resolve_display_smoothing(smooth=.3).factor((30, 30, 30), .25) == 2


def test_cli_smooth_options():
    from crevice.cli import build_parser
    parser = build_parser()
    for command in (['cast', 'input.pdb', '--out-dir', 'out'], ['publish', 'input.pdb', '--out-dir', 'out']):
        args = parser.parse_args(command)
        assert args.smooth is None and args.surface_smoothing is None
        assert parser.parse_args(command+['--smooth', '0.9']).smooth == .9
        assert parser.parse_args(command+['--smooth', '0']).smooth == 0
        assert parser.parse_args(command+['--surface-smoothing', '.6']).surface_smoothing == .6
        for bad in (['--smooth', '-0.1'], ['--smooth', 'nan'], ['--smooth', 'inf'], ['--surface-smoothing', '-1'],
                    ['--smooth', '.3', '--surface-smoothing', '.6']):
            with pytest.raises(SystemExit):
                parser.parse_args(command+bad)


@pytest.fixture(scope='module')
def channel_display():
    frame = rectangular_channel()
    profile = pore_profile(frame, enclosure_radius=0.8, search_radius=8, samples=41)
    cast = void_cast(frame, profile=profile, spacing=.5, min_radius=0, max_grid_points=200000)
    return frame, profile, cast


def test_bundle_metadata_mesh_and_measured_dx(channel_display, tmp_path):
    frame, profile, cast = channel_display
    structure = write_viewer_structure(frame, tmp_path/'viewer.pdb')
    before = cast.to_dict(include_points=True)
    runs = {}
    for name, options in {'raw': {'smooth': 0}, 'default': {}, 'wide': {'smooth': 1.2},
                          'legacy': {'surface_smoothing': .6}}.items():
        runs[name] = write_volume_viewer_bundle(cast, structure_path=structure, output_dir=tmp_path/name,
                                                frame=frame, profile=profile, **options)
    measured = {Path(files['volume_dx']).read_bytes() for files in runs.values()}
    assert len(measured) == 1
    assert cast.to_dict(include_points=True) == before
    display_shape = cast_grid(extend_channel_cast(cast, frame, 2))[0].shape
    for name, files in runs.items():
        meta = json.loads(Path(files['volume_metadata_json']).read_text())
        smoothing = meta['display_smoothing']
        text = Path(files['display_dx']).read_text().splitlines()
        counts = tuple(int(v) for v in text[1].split()[-3:])
        delta = float(np.linalg.norm([float(v) for v in text[3].split()[1:]]))
        k = smoothing['supersample_factor']
        assert counts == tuple(k*(n-1)+1 for n in display_shape)
        assert delta == pytest.approx(cast.spacing/k)
        assert meta['unsmoothed_reference_dx'] == files['volume_dx']
        assert meta['measured_volume_A3'] == cast.total_volume
        if name == 'raw':
            assert smoothing['method'] == 'none' and k == 1 and meta['smoothing'] is False
        elif name == 'legacy':
            assert smoothing['method'] == 'legacy_clamped_gaussian' and meta['surface_smoothing_A'] == .6
        else:
            assert smoothing['method'] == SMOOTH_METHOD and k == 2
            assert smoothing['width_A'] == (1.2 if name == 'wide' else DEFAULT_SMOOTH_A)
            assert meta['smooth_A'] == smoothing['width_A']
            assert 'standard deviation' in smoothing['definition'] and meta['smoothing_definition'] == smoothing['definition']
        with np.load(files['pymol_mesh_npz']) as mesh:
            vertices = mesh['volume_vertices']
        assert len(vertices) > 0
    with pytest.raises(ValueError, match='not both'):
        write_volume_viewer_bundle(cast, structure_path=structure, output_dir=tmp_path/'both', frame=frame,
                                   profile=profile, smooth=.3, surface_smoothing=.6)
    assert not (tmp_path/'both').exists()


MEASUREMENT_FILES = ('_profile.json', '_profile.csv', '_void_cast.json', '_pore_cast.pdb', '_volume.dx',
                     '_residue_contacts.csv', '_void_cast.csv', '_network.json', '_network_nodes.csv',
                     '_network_edges.csv', '_scene.json')


def test_publication_measurements_identical_with_and_without_smoothing(tmp_path):
    from crevice.figures import write_static_publication_bundle
    frame = rectangular_channel()
    source = write_viewer_structure(frame, tmp_path/'channel.pdb')
    outputs = {}
    for name, options in {'raw': {'smooth': 0}, 'default': {}, 'wide': {'smooth': 1.2},
                          'legacy': {'surface_smoothing': .6}}.items():
        files = write_static_publication_bundle(frame, structure_path=source, output_dir=tmp_path/name, prefix='c',
                                                samples=41, search_radius=8, cast_spacing=.5, dpi=60,
                                                include_cavities=False, include_tunnels=False, **options)
        outputs[name] = files
    for suffix in MEASUREMENT_FILES:
        blobs = {(tmp_path/name/('c'+suffix)).read_bytes() for name in outputs}
        assert len(blobs) == 1, suffix
    displays = {(tmp_path/name/'c_display.dx').read_bytes() for name in outputs}
    assert len(displays) == 4
    with pytest.raises(ValueError, match='not both'):
        write_static_publication_bundle(frame, enclosure_radius=0.8, structure_path=source, output_dir=tmp_path/'both',
                                        smooth=.3, surface_smoothing=.6)


def test_display_field_has_no_exact_isovalue_nodes():
    from skimage.measure import marching_cubes
    grid = two_rooms_with_neck()
    for width in (.2, .4, .8):
        field = constrained_display_field(grid, .25, width, 2)
        assert not np.any(field == .5)
        vertices, faces, _, _ = marching_cubes(field, .5)
        edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
        _, uses = np.unique(edges, axis=0, return_counts=True)
        assert np.all(uses == 2)  # closed two-manifold contour

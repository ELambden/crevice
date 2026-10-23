"""Streaming of cavity/region trajectory frames (memory contract of ``trajectory --cavities``).

Software tests on synthetic geometry only; no biological result is implied.
They pin that frames reach the residue dynamics in time order whatever order
the workers finish in, that at most ``4 * workers`` frames are in flight, that
the compact per-frame records return exactly the rows of
:func:`crevice.cavity_trajectory.analyze_cavity_frame`, and that the written
bundle equals one written from the full, unstreamed records.
"""
import copy
import math
import pickle
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import crevice.cavity_trajectory as ct
from crevice.cavity_trajectory import (CavityReference, _BoundaryResidues, _uninvolved_row,
                                       analyze_cavity_frame, analyze_cavity_trajectory,
                                       write_cavity_trajectory_bundle)
from crevice.residue_dynamics import ResidueDynamics
from crevice.rolling import rolling_probe_cast
from crevice.trajectory import Trajectory
from crevice.volume_export import write_void_cast_dx
from test_channel_coordinate import rectangular_channel

OPTIONS = dict(geometry_mode='reference-region', geometry={'probe_radius': 1.4})


def _frames(count):
    frame = rectangular_channel(half_x=3, half_y=4, half_length=5, capped=True)
    # A separate cluster of residues far from the cavity: their boundary rows are 'other'.
    far = tuple(replace(frame.atoms[0], serial=len(frame.atoms) + k + 1, resid=len(frame.atoms) + k + 1,
                        x=30. + 1.6 * (k % 4), y=1.6 * (k // 4), z=0.) for k in range(16))
    frame = replace(frame, atoms=frame.atoms + far)
    rng = np.random.default_rng(7)
    frames = []
    for index in range(count):
        jitter = rng.normal(scale=.15, size=(len(frame.atoms), 3)) if index else np.zeros((len(frame.atoms), 3))
        atoms = tuple(replace(a, x=a.x + float(d[0]), y=a.y + float(d[1]), z=a.z + float(d[2]))
                      for a, d in zip(frame.atoms, jitter))
        frames.append(replace(frame, atoms=atoms, frame_index=index))
    return frames


@pytest.fixture(scope='module')
def system(tmp_path_factory):
    frames = _frames(12)
    cast = rolling_probe_cast(frames[0], spacing=.75, selection_mode='all', enclosure_fraction=0, min_component_volume=1)
    main = max(cast.components, key=lambda c: c.volume)
    cast = replace(cast, components=(main,), points=main.points)
    ref = CavityReference.from_dx(write_void_cast_dx(cast, tmp_path_factory.mktemp('ref') / 'ref.dx'), axis=(0, 0, 1))
    return Trajectory(frames, time_step=10.), ref


@pytest.fixture(scope='module')
def serial(system):
    trajectory, ref = system
    return analyze_cavity_trajectory(trajectory, ref, workers=1, **OPTIONS)


@pytest.fixture(scope='module')
def unstreamed(system):
    return _unstreamed(*system)


def _unstreamed(trajectory, reference):
    """Records as 0.4.0 kept them: full analyze_cavity_frame rows, network removed."""
    rows = []
    for frame, t in zip(trajectory.frames, trajectory.times()):
        row = analyze_cavity_frame(frame, trajectory.frames[0], reference, **OPTIONS)
        row.pop('network'); row['time_ps'] = float(t); rows.append(row)
    return rows


def test_records_are_compact_and_return_the_full_rows(system, serial, unstreamed):
    trajectory, ref = system
    full = copy.deepcopy(unstreamed)
    for row, expected in zip(serial['frames'], full):
        assert 'network' not in row and 'network' not in row['boundary']
        residues = row['boundary']['residues']
        assert isinstance(residues, _BoundaryResidues)
        rows = list(residues)
        assert rows == expected['boundary']['residues'] and residues == expected['boundary']['residues']
        # Same key order and value types, not only equal values.
        for got, want in zip(rows, expected['boundary']['residues']):
            assert list(got.items()) == list(want.items())
            assert [type(v) for v in got.values()] == [type(v) for v in want.values()]
        assert sum(isinstance(item, tuple) for item in residues._items) > 0   # uninvolved rows are compact
        expected['boundary'].pop('network')
        assert {k: v for k, v in row['boundary'].items() if k != 'residues'} == \
               {k: v for k, v in expected['boundary'].items() if k != 'residues'}
    short = Trajectory(trajectory.frames[:2], time_step=10.)
    kept = analyze_cavity_trajectory(short, ref, workers=1, retain_networks=True, **OPTIONS)
    assert all('network' in row['boundary'] and 'network' not in row for row in kept['frames'])


def test_uninvolved_rows_are_shared_between_frames_and_fresh_on_access(serial):
    first, second = (row['boundary']['residues'] for row in serial['frames'][:2])
    shared = {item for item in first._items if isinstance(item, tuple)}
    assert shared and any(any(item is key for key in shared) for item in second._items)
    index = next(i for i, item in enumerate(first._items) if isinstance(item, tuple))
    edited = first[index]; edited['boundary_area_A2'] = 99.
    assert first[index]['boundary_area_A2'] == 0.
    assert pickle.loads(pickle.dumps(first)) == first and copy.deepcopy(first) == list(first)


def test_only_exact_uninvolved_rows_are_compacted():
    template = _uninvolved_row('A:ALA:1', 'unassigned')
    near = [dict(template, boundary_area_A2=-0.), dict(template, contacted_lining_boundary_fraction=0.),
            dict(template, minimum_boundary_gap_A=1.2), dict(template, group=None), {**template, 'extra': 1},
            dict(template, role='boundary_lining'), dict(template, nonlocal_lining_partners=[{'residue': 'A:GLY:2'}])]
    blank = _uninvolved_row('A:ALA:2')
    table = _BoundaryResidues([template, *near, blank])
    assert isinstance(table._items[0], tuple) and isinstance(table._items[-1], tuple)
    assert all(isinstance(item, dict) for item in table._items[1:-1])
    assert table == [template, *near, blank] and table[1:3] == near[:2]
    assert math.copysign(1., table[1]['boundary_area_A2']) < 0


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='patches the worker function, which needs fork')
def test_out_of_order_completion_feeds_dynamics_in_time_order(system, serial, monkeypatch):
    trajectory, ref = system
    monkeypatch.setenv('CREVICE_MP_START_METHOD', 'fork')
    original = ct._trajectory_frame

    def slow_first(frame, *args, **kwargs):
        if frame.frame_index == 0:
            time.sleep(1.5)
        return original(frame, *args, **kwargs)

    monkeypatch.setattr(ct, '_trajectory_frame', slow_first)
    added = []
    add = ResidueDynamics.add

    def recording(self, frame, network, *, frame_index, time):
        added.append(frame_index)
        return add(self, frame, network, frame_index=frame_index, time=time)

    monkeypatch.setattr(ResidueDynamics, 'add', recording)
    completed = []
    parallel = analyze_cavity_trajectory(trajectory, ref, workers=2, progress_callback=lambda i, row: completed.append(i), **OPTIONS)
    assert completed[0] != 0 and sorted(completed) == list(range(12))   # frame 0 finished late
    # A frame is submitted only within 4 x workers of the oldest unfinished frame.
    finished = set()
    for index in completed:
        oldest = min(set(range(12)) - finished)
        assert index < oldest + 8
        finished.add(index)
    assert added == list(range(12))
    assert parallel['residue_dynamics'] == serial['residue_dynamics']
    assert np.array_equal(parallel['residue_centroids'], serial['residue_centroids'])
    assert [r['boundary']['residues'] for r in parallel['frames']] == [r['boundary']['residues'] for r in serial['frames']]


def test_bundle_equals_the_bundle_of_unstreamed_records(system, serial, unstreamed, tmp_path):
    trajectory, ref = system
    parallel = analyze_cavity_trajectory(trajectory, ref, workers=3, **OPTIONS)
    options = dict(replicates=200, dpi=40)
    a = write_cavity_trajectory_bundle(serial, tmp_path / 'a', **options)
    b = write_cavity_trajectory_bundle(dict(serial, frames=unstreamed), tmp_path / 'b', **options)
    parallel['settings']['workers'] = 1
    c = write_cavity_trajectory_bundle(parallel, tmp_path / 'c', **options)
    assert set(a) == set(b) == set(c)
    for key, path in a.items():
        if path.endswith(('.csv', '.json')) and not path.endswith('_manifest.json'):
            text = Path(path).read_bytes()
            assert text == Path(b[key]).read_bytes() == Path(c[key]).read_bytes(), key
        elif path.endswith('.npz'):
            with np.load(path) as x, np.load(b[key]) as y, np.load(c[key]) as z:
                for name in x.files:
                    assert np.array_equal(x[name], y[name], equal_nan=x[name].dtype.kind == 'f')
                    assert np.array_equal(x[name], z[name], equal_nan=x[name].dtype.kind == 'f')
    assert Path(a['residue_frames_csv']).stat().st_size > 0

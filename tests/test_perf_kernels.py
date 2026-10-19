"""The batched geometry kernels return exactly (bit for bit) what the scalar code returns.

Software checks on synthetic atoms: the vectorised section-plane graph, the
connected-component labelling, the batched point and segment clearances and
the shared probe-ladder frame cache must not change any value.
"""
import math
import random

import numpy as np
import pytest

from crevice.geometry import point_segment_distance
from crevice.models import Atom
from crevice.sections import SectionGeometry, _FrameCache, _STEPS, connected_profile
from crevice.spatial import SpatialIndex


def _random_atoms(seed, count=400, size=14.0):
    rng = random.Random(seed)
    elements = ("C", "N", "O", "S", "H")
    return tuple(Atom(i + 1, "X", "UNK", "A", i + 1, rng.uniform(0, size), rng.uniform(0, size),
                      rng.uniform(0, size), rng.choice(elements)) for i in range(count))


def _ring_channel(radius=4.2, layers=14, per_ring=12):
    atoms, serial = [], 1
    for k in range(layers):
        for m in range(per_ring):
            angle = 2 * math.pi * (m + 0.5 * (k % 2)) / per_ring
            atoms.append(Atom(serial, "C", "UNK", "A", serial, radius * math.cos(angle),
                              radius * math.sin(angle), 1.5 * k - 10.0, "C"))
            serial += 1
    return tuple(atoms)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_batched_point_and_segment_clearances_are_bit_identical(seed):
    atoms = _random_atoms(seed)
    index = SpatialIndex(atoms)
    rng = random.Random(seed + 100)
    points = [(rng.uniform(-2, 16), rng.uniform(-2, 16), rng.uniform(-2, 16)) for _ in range(600)]
    assert index.nearest_surface_values(points).tolist() == [index.nearest_surface(p)[0] for p in points]
    starts, ends = points[:300], points[300:]
    ends[0] = starts[0]  # zero-length segment
    reference = [min(point_segment_distance(a.coord, s, e) - index._radii[a] for a in atoms)
                 if s != e else index.nearest_surface(s)[0] for s, e in zip(starts, ends)]
    assert [index.segment_clearance(s, e) for s, e in zip(starts, ends)] == reference
    assert index.segment_clearances(starts, ends).tolist() == reference


@pytest.mark.parametrize("probe", [0.0, 0.5, 0.8, 1.4, 3.0])
@pytest.mark.parametrize("atoms", [_random_atoms(7, 250), _ring_channel()], ids=["random", "ring"])
def test_plane_graph_matches_scalar_plane_masks(atoms, probe):
    centre = tuple(sum(a.coord[k] for a in atoms) / len(atoms) for k in range(3))
    geometry = SectionGeometry(atoms, centre, (0.0, 0.0, 1.0), spacing=0.5, extent=8.0, probe=probe)
    h = geometry.half
    n = 2 * h + 1
    for t in np.linspace(-9.0, 9.0, 13):
        free_ref, blocked_ref = geometry._plane_masks(float(t))
        free, blocked = geometry._plane_graph(float(t))
        assert np.array_equal(free, free_ref)
        for k, (di, dj) in enumerate(_STEPS):
            i, j = np.nonzero(free)
            ti, tj = i + di, j + dj
            ok = (ti >= 0) & (ti < n) & (tj >= 0) & (tj < n)
            i, j, ti, tj = i[ok], j[ok], ti[ok], tj[ok]
            both = free[ti, tj]
            assert np.array_equal(blocked[k, i[both], j[both]], blocked_ref[k, i[both], j[both]])


def test_shared_frame_cache_gives_the_independent_result():
    atoms = _ring_channel()
    kwargs = dict(axis=(0.0, 0.0, 1.0), origin=(0.0, 0.0, 0.0), samples=21, search_radius=7.0,
                  refinement_steps=4, probe_radius=0.0, section_spacing=0.5)
    alone = connected_profile(atoms, enclosure_radius=0.8, **kwargs)
    from crevice import sections
    token = sections._ACTIVE_FRAMES.set(_FrameCache())
    try:
        connected_profile(atoms, enclosure_radius=0.6, **kwargs)
        shared = connected_profile(atoms, enclosure_radius=0.8, **kwargs)
    finally:
        sections._ACTIVE_FRAMES.reset(token)
    assert shared == alone


def test_dx_value_rows_equal_formatting_every_value():
    from crevice.volume_export import dx_value_rows

    def reference(values, spec):
        return [" ".join(format(v, spec) for v in values[i:i + 3]) for i in range(0, len(values), 3)]

    rng = np.random.default_rng(7)
    for dtype in (np.float32, np.float64, np.int64):
        for n in (0, 1, 2, 3, 10, 1001):
            values = (rng.random(n) * 3).round(3).astype(dtype)
            if n > 3 and dtype is not np.int64:
                values[:3] = (-0.0, np.nan, np.inf)
            for spec in (".6g", ".12g"):
                assert dx_value_rows(values, spec) == reference(values, spec)


def test_pruned_maximum_clearance_equals_the_full_maximum():
    from crevice.cavity_obstacles import _max_clearance
    from crevice.rolling import _SphereQueries
    rng = np.random.default_rng(3)
    coords = rng.uniform(0, 12, (300, 3))
    queries = _SphereQueries(coords, rng.choice([1.5, 1.7, 1.8], 300))
    for shape in ((17, 23, 9), (5, 4, 3), (1, 1, 1)):
        xyz = np.indices(shape).reshape(3, -1).T * 0.37 - 1.0
        assert _max_clearance(queries, xyz, shape) == float(queries.points(xyz)[0].max())


def test_near_only_regional_distances_are_the_full_ones_within_the_bound():
    from types import SimpleNamespace
    from crevice.regional_volume import regional_grid
    rng = np.random.default_rng(5)
    points = np.unique(np.rint(rng.normal(0, 1.5, (400, 3)) / 0.4), axis=0) * 0.4
    make = lambda: SimpleNamespace(points=points, basis=np.eye(3), anchor=np.zeros(3), spacing=0.4)
    full = regional_grid(make(), 2.0, 1.4, 10_000_000)[3]
    near = regional_grid(make(), 2.0, 1.4, 10_000_000, near_only=True)[3]
    within = full < 2.0 + 1.4 + 0.4
    assert np.array_equal(full[within], near[within]) and np.isinf(near[~within]).all()


def test_lateral_exit_legs_do_not_depend_on_the_prefetch_block(monkeypatch):
    # Segment clearances are prepared per lattice block; the block size (and
    # so what is prepared ahead) must not change any leg.
    import sys
    from pathlib import Path
    import crevice.channel_exits as exits
    from crevice.analysis import pore_profile
    sys.path.insert(0, str(Path(__file__).parent))
    from test_capped_channels import OPTIONS, capped_channel
    frame = capped_channel(second_portal=True)
    results = []
    for block in (1, 2, 4):
        monkeypatch.setattr(exits, "_BLOCK", block)
        results.append(pore_profile(frame, lateral_exits=True, **OPTIONS).metadata["exits"])
    assert results[0] == results[1] == results[2]


def test_route_overlap_test_matches_the_plain_definition():
    from crevice.paths import _min_pair_distance, _routes_overlap

    def plain(pos_a, pos_b, ids_a, ids_b, separation, route_separation):
        if np.linalg.norm(pos_a[-1] - pos_b[-1]) <= separation:
            return True
        shared = 0
        while shared < min(len(ids_a), len(ids_b)) and ids_a[shared] == ids_b[shared]:
            shared += 1
        junction = pos_a[shared - 1] if shared else None
        one, two = pos_a[shared:], pos_b[shared:]
        if junction is not None:
            one = one[np.linalg.norm(one - junction, axis=1) > route_separation]
            two = two[np.linalg.norm(two - junction, axis=1) > route_separation]
        return not len(one) or not len(two) or _min_pair_distance(one, two) <= route_separation

    rng = np.random.default_rng(11)
    seen = set()
    for _ in range(300):
        prefix = rng.integers(0, 6)
        steps_a = rng.choice([-0.5, 0.0, 0.5], (rng.integers(1, 40), 3))
        steps_b = rng.choice([-0.5, 0.0, 0.5], (rng.integers(1, 40), 3))
        prefix = min(prefix, len(steps_a), len(steps_b))
        steps_b[:prefix] = steps_a[:prefix]
        pos_a, pos_b = np.cumsum(steps_a, axis=0), np.cumsum(steps_b, axis=0)
        ids_a = list(range(len(pos_a)))
        ids_b = list(range(prefix)) + list(range(1000, 1000 + len(pos_b) - prefix))
        for separation, route_separation in ((3.0, 1.0), (12.0, 6.0), (1.0, 0.5)):
            result = _routes_overlap(pos_a, pos_b, ids_a, ids_b, separation, route_separation)
            assert result == plain(pos_a, pos_b, ids_a, ids_b, separation, route_separation)
            seen.add(result)
    assert seen == {True, False}

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

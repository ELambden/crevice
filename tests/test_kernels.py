"""The compiled kernels (crevice._kernels, numba) give exactly the NumPy results.

Software tests on synthetic structures. Skipped where numba is not installed;
the NumPy path itself is covered by every other profile test.
"""
import math

import numpy as np
import pytest

from crevice import Atom, StructureFrame, pore_profile
from crevice import _kernels, adaptive_channel
from crevice.sections import SectionGeometry, _STEPS

numba_only = pytest.mark.skipif(not _kernels.HAVE_NUMBA, reason="numba not installed")


def _lattice(*, waist=None, tube_z=None):
    """2 Å carbon lattice with a 3 Å lumen along z, an optional waist and side tube (as test_adaptive_enclosure)."""
    def segment_distance(p, a, b):
        ab = [b[i] - a[i] for i in range(3)]
        ap = [p[i] - a[i] for i in range(3)]
        t = max(0.0, min(1.0, sum(x * y for x, y in zip(ab, ap)) / sum(x * x for x in ab)))
        return math.dist(p, [a[i] + t * ab[i] for i in range(3)])
    coords = []
    for i in range(-4, 5):
        for j in range(-4, 5):
            for k in range(-4, 5):
                p = (2.0 * i, 2.0 * j, 2.0 * k)
                if math.hypot(p[0], p[1]) < 3.0:
                    continue
                if tube_z is not None and segment_distance(p, (0.0, 0.0, tube_z), (12.0, 0.0, tube_z)) < 2.5:
                    continue
                coords.append(p)
    if waist is not None:
        coords += [(waist, 0.0, 4.0), (-waist, 0.0, 4.0), (0.0, waist, 4.0), (0.0, -waist, 4.0)]
    coords.sort(key=lambda p: (p != (8.0, 0.0, -8.0), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


@numba_only
@pytest.mark.parametrize("probe", [0.3, 0.7, 1.2])
@pytest.mark.parametrize("t", [-4.0, -1.3, 0.0, 4.0, 5.1])
def test_plane_graph_kernel_equals_numpy(probe, t):
    frame = _lattice(waist=2.3, tube_z=-4.0)
    geometry = SectionGeometry(frame.atoms, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), spacing=0.25, extent=7.0, probe=probe)
    state = geometry._state
    free, blocked = _kernels.plane_graph(state.px, state.py, state.pz, state.pr, t, geometry.half,
                                         geometry.spacing, probe, _kernels.STEPS)
    ref_free, ref_blocked = geometry._plane_graph(t)
    assert np.array_equal(free, ref_free)
    assert np.array_equal(blocked, ref_blocked)


def test_kernel_steps_match_the_sections_order():
    assert [tuple(row) for row in _kernels.STEPS.tolist()] == list(_STEPS)
    assert all(tuple(_kernels.STEPS[7 - k]) == (-di, -dj) for k, (di, dj) in enumerate(_STEPS))


@numba_only
@pytest.mark.parametrize("probe", [0.3, 0.7, 1.2])
def test_component_kernel_equals_numpy_labels(probe, monkeypatch):
    frame = _lattice(waist=2.3, tube_z=-4.0)
    tracer = adaptive_channel._Tracer(frame.atoms, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), spacing=0.25, extent=7.0,
                                      ladder=adaptive_channel.ADAPTIVE_ENCLOSURE_LADDER, refinement_steps=4,
                                      shared=None)
    for t in (-4.0, 0.0, 4.0):
        geometry = tracer.geometry(probe, 1)
        compiled = tracer._component(geometry, t, 0, 0)
        monkeypatch.setattr(_kernels, "HAVE_NUMBA", False)
        reference = tracer._component(geometry, t, 0, 0)
        monkeypatch.setattr(_kernels, "HAVE_NUMBA", True)
        if reference is None:
            assert compiled is None
        else:
            assert np.array_equal(compiled[0], reference[0]) and np.array_equal(compiled[1], reference[1])


@numba_only
def test_profile_is_identical_with_and_without_numba(monkeypatch):
    frame = _lattice(waist=2.3, tube_z=-4.0)
    options = dict(axis="z", origin=(0.0, 0.0, -1.0), samples=41, search_radius=7.0, section_spacing=0.25)
    compiled = pore_profile(frame, **options)
    monkeypatch.setattr(_kernels, "HAVE_NUMBA", False)
    reference = pore_profile(frame, **options)
    assert compiled.points == reference.points
    assert compiled.metadata == reference.metadata


# -- lateral-exit search ------------------------------------------------------------------------------------


def test_round10_is_python_round():
    rng = np.random.default_rng(7)
    values = np.concatenate(
        [
            rng.uniform(-10.0, 10.0, 20000),
            rng.uniform(0.0, 1e-6, 2000),
            10.0 ** rng.uniform(-12, 4.9, 5000),
            # exact decimal ties (odd multiples of 2**-11 are ties at ten decimals) and their neighbours
            np.arange(1, 4001, 2) / 2048.0,
            np.nextafter(np.arange(1, 4001, 2) / 2048.0, 1.0),
            np.nextafter(np.arange(1, 4001, 2) / 2048.0, -1.0),
            np.round(rng.uniform(0, 5, 2000), 10) + 5e-11,
        ]
    )
    for x in values.tolist() + [0.0, -0.0, math.inf, 1e-300, 99999.99999999995]:
        assert _kernels.round10(x) == round(x, 10), x
    assert _kernels.round10(-math.inf) == -math.inf


def _mixed_atoms():
    from test_capped_channels import capped_channel

    frame = capped_channel(second_portal=True)
    elements = ("C", "N", "O", "S")
    return tuple(
        Atom(a.serial, elements[a.serial % 4], "ALA", "A", a.serial, *a.coord, elements[a.serial % 4])
        for a in frame.atoms
    )


def test_compiled_clearances_equal_the_spatial_index():
    from crevice.spatial import SpatialIndex

    atoms = _mixed_atoms()
    spatial = SpatialIndex(atoms)
    cells = _kernels.atom_cells(spatial._xyz, spatial._radius_array)
    xyz, radii, origin, cell, dims, starts, order = cells
    rmax = float(radii.max())
    rng = np.random.default_rng(3)
    points = rng.uniform(-12.0, 12.0, (3000, 3))
    points[:500] = np.round(points[:500] * 2) / 2  # lattice-like points, exact ties between atoms
    reference = spatial.nearest_surface_values(points)
    for p, ref in zip(points, reference):
        value = _kernels.point_clearance(p[0], p[1], p[2], max(ref, 0.0) + 0.9 + rmax, *cells)
        assert value == ref
    ends = points + rng.uniform(-0.6, 0.6, points.shape)
    segments = spatial.segment_clearances(points, ends)
    probe = 0.3
    for a, b, ref in zip(points, ends, segments):
        for limit in (probe + 1e-6, 0.8, 2.5, 50.0):
            value = _kernels.segment_clearance_below(*a, *b, limit, probe, rmax, *cells)
            if ref <= probe:
                assert value <= probe
            else:
                assert value == min(limit, ref)


def test_bulk_state_defers_to_numpy_at_rounding_ties():
    distance = np.zeros(8)
    distance[1::2] = 10.0  # cells alternate bulk / not bulk along the last axis
    common = (np.zeros(3), np.eye(3), np.zeros(3), 1.0, np.array([2, 2, 2], dtype=np.int64), distance, 6.0)
    assert _kernels._bulk_state(np.array([0.2, 0.2, 0.5]), *common) == 2
    assert _kernels._bulk_state(np.array([0.2, 0.2, 0.2]), *common) == 1
    assert _kernels._bulk_state(np.array([0.2, 0.2, 0.8]), *common) == 0
    assert _kernels._bulk_state(np.array([0.2, 0.2, -0.7]), *common) == 1  # off the grid


@numba_only
@pytest.mark.parametrize(
    "variant, probe, spacing",
    [
        ({}, 0.8, 0.5),
        ({}, None, 0.5),
        ({"lumen_top": 6.0}, 0.8, 0.5),
        ({"second_portal": True}, 0.8, 0.5),
        ({"second_portal": True}, 0.8, 0.35),
        ({"second_portal": True}, None, 0.5),
        ({"portal": False}, 0.8, 0.5),
    ],
)
def test_lateral_exit_kernel_equals_numpy(variant, probe, spacing, monkeypatch):
    from crevice.sections import ChannelResolutionError
    from test_capped_channels import capped_channel

    frame = capped_channel(**variant)
    options = dict(axis="z", origin=(0.0, 0.0, 0.0), samples=25, search_radius=7.0, section_spacing=0.5)
    options.update(enclosure_radius=probe, lateral_exits=True, exit_spacing=spacing)

    def run():
        try:
            return pore_profile(frame, **options).to_dict()
        except ChannelResolutionError as exc:
            return str(exc)

    compiled = run()
    monkeypatch.setattr(_kernels, "HAVE_NUMBA", False)
    reference = run()
    assert compiled == reference
    if isinstance(reference, dict):
        assert reference["metadata"]["path_type"] == "capped_channel_lateral_exit"


@numba_only
def test_lateral_exit_legs_kernel_equals_numpy_for_many_legs(monkeypatch):
    from crevice.channel_exits import lateral_exit_legs
    from crevice.sections import SectionGeometry, Section
    from test_capped_channels import capped_channel

    frame = capped_channel(second_portal=True)
    geometry = SectionGeometry(frame.atoms, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), spacing=0.5, extent=7.0, probe=0.5)
    start = Section(2.0, (), (0.0, 0.0, 2.0), geometry.spatial.nearest_surface((0.0, 0.0, 2.0))[0])
    for kwargs in ({}, {"max_legs": 1}, {"separation": 1.0, "route_separation": 1.0}, {"max_nodes": 400}):
        compiled = lateral_exit_legs(geometry, start, 1, end="upper", **kwargs)
        monkeypatch.setattr(_kernels, "HAVE_NUMBA", False)
        reference = lateral_exit_legs(geometry, start, 1, end="upper", **kwargs)
        monkeypatch.setattr(_kernels, "HAVE_NUMBA", True)
        assert [leg.to_dict() for leg in compiled[0]] == [leg.to_dict() for leg in reference[0]]
        assert compiled[1] == reference[1]
    with pytest.raises(ValueError, match="exceeded 5 lattice nodes"):
        lateral_exit_legs(geometry, start, 1, end="upper", max_nodes=5)


@numba_only
@pytest.mark.parametrize("spacing, bulk_radius, sign", [(1.0, 6.0, 1), (0.5, 3.0, -1), (0.7, 4.5, 1)])
def test_bulk_envelope_kernel_equals_numpy(spacing, bulk_radius, sign, monkeypatch):
    from crevice.channel_exits import BulkEnvelope

    atoms = _mixed_atoms()
    args = (atoms, (0.3, -0.2, 0.1), ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), 2.0, sign)
    compiled = BulkEnvelope(*args, spacing=spacing, bulk_radius=bulk_radius)
    monkeypatch.setattr(_kernels, "HAVE_NUMBA", False)
    reference = BulkEnvelope(*args, spacing=spacing, bulk_radius=bulk_radius)
    assert compiled.dims == reference.dims
    assert np.array_equal(compiled.distance, reference.distance)
    assert (reference.distance <= bulk_radius).any() and (reference.distance > bulk_radius).any()

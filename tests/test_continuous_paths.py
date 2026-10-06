"""Atom-sphere segment clearance and grid-path regression contracts."""

import math
from dataclasses import replace

import numpy as np
import pytest

from crevice import StructureFrame, find_tunnels, void_cast
from crevice import analysis, voids
from crevice.spatial import SpatialIndex
from test_ensemble_cast_network import atom, cage, ring


def diagonal_fixture():
    frame = StructureFrame((atom(1, (1.5, 1.5, 1.5)),
                            atom(2, (0, 0, 6)), atom(3, (6, 6, 0))))
    spatial = SpatialIndex(frame.atoms)
    nodes = {}
    for idx, position in (((0, 0, 0), (0, 0, 0)), ((1, 1, 1), (3, 3, 3))):
        clearance, nearest = spatial.nearest_surface(position)
        assert clearance > 0.2
        nodes[idx] = voids._GridNode(position, clearance, nearest, 8)
    return frame, spatial, nodes


def test_clear_diagonal_endpoints_do_not_merge_across_atom(monkeypatch):
    frame, _spatial, nodes = diagonal_fixture()
    monkeypatch.setattr(voids, "_collect_void_nodes", lambda *args, **kwargs: nodes)
    cast = void_cast(frame, mode="all", axis="z", spacing=3,
                     min_radius=0.2, min_component_points=1)
    assert cast.component_count == 2
    assert [c.point_count for c in cast.components] == [1, 1]
    assert len({c.metadata["source_component_id"] for c in cast.components}) == 2


def test_tunnel_retains_actual_off_grid_start():
    frame = ring()
    start = (0.2, 0.3, 0.4)
    tunnels = find_tunnels(frame, start=start, spacing=2, min_radius=0.5, max_tunnels=1)
    assert tunnels
    assert tunnels[0].start == start
    assert tunnels[0].points[0].position == start
    spatial = SpatialIndex(frame.atoms)
    tunnel = tunnels[0]
    clearances = [spatial.segment_clearance(a.position, b.position)
                  for a, b in zip(tunnel.points, tunnel.points[1:])]
    assert min(clearances) >= 0.5
    assert tunnel.bottleneck_radius == pytest.approx(min(clearances))
    assert tunnel.metadata["segment_clearances_A"] == pytest.approx(clearances)
    assert tunnel.length == pytest.approx(sum(math.dist(a.position, b.position)
                                            for a, b in zip(tunnel.points, tunnel.points[1:])))


@pytest.mark.parametrize("start,end,expected", [
    ((-3, 0, 0), (3, 0, 0), -1.7),
    ((-3, -3, 0), (3, 3, 0), -1.7),
    ((-3, -3, -3), (3, 3, 3), -1.7),
    ((-3, 2, 0), (3, 2, 0), 0.3),
    ((2, 0, 0), (4, 0, 0), 0.3),
    ((0, 2, 0), (0, 2, 0), 0.3),
])
def test_analytic_segment_clearance(start, end, expected):
    spatial = SpatialIndex((atom(1, (0, 0, 0)),))
    assert spatial.segment_clearance(start, end) == pytest.approx(expected, abs=1e-12)
    assert spatial.segment_clearance(end, start) == pytest.approx(expected, abs=1e-12)
    spatial._tree = None
    assert spatial.segment_clearance(start, end) == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("endpoint", [(math.nan, 0, 0), (0, math.inf, 0), (0, 0), (0, 0, 0, 0)])
def test_nonfinite_or_malformed_segment_rejected(endpoint):
    spatial = SpatialIndex((atom(1, (0, 0, 0)),))
    with pytest.raises(ValueError, match="finite 3D"):
        spatial.segment_clearance((0, 0, 0), endpoint)


def test_tree_and_fallback_match_independent_scalar_minimization():
    from scipy.optimize import minimize_scalar
    from crevice.radii import atom_vdw_radius

    rng = np.random.default_rng(913)
    atoms = tuple(atom(i, xyz, element=("C", "O", "H", "K")[i % 4])
                  for i, xyz in enumerate(rng.uniform(-10, 10, (24, 3)), 1))
    spatial = SpatialIndex(atoms)
    tree = spatial._tree
    for start, end in rng.uniform(-15, 15, (12, 2, 3)):
        # Minimize each atom's convex distance along the segment independently.
        best = []
        for a in atoms:
            def objective(t):
                return np.linalg.norm(start + t * (end - start) - a.coord) - atom_vdw_radius(a)
            optimum = minimize_scalar(objective, bounds=(0, 1), method="bounded",
                                       options={"xatol": 1e-13})
            best.append(min(objective(0), objective(1), optimum.fun))
        spatial._tree = tree
        accelerated = spatial.segment_clearance(tuple(start), tuple(end))
        spatial._tree = None
        assert accelerated == pytest.approx(min(best), abs=1e-7)
        assert spatial.segment_clearance(tuple(start), tuple(end)) == pytest.approx(accelerated, abs=1e-12)


def test_segment_rigid_transform_and_large_radius_endpoint_candidate():
    atoms = (atom(1, (0, 2, 0)), atom(2, (9, 2.6, 0), element="K"))
    start, end = (-10, 0, 0), (10, 0, 0)
    expected = SpatialIndex(atoms).segment_clearance(start, end)
    assert expected == pytest.approx(-0.15)
    def transform(point):
        x, y, z = point
        return (-y + 12, z - 7, -x + 3)
    transformed = tuple(replace(a, x=transform(a.coord)[0], y=transform(a.coord)[1],
                                z=transform(a.coord)[2]) for a in atoms)
    assert SpatialIndex(transformed).segment_clearance(transform(start), transform(end)) == pytest.approx(expected)


def test_flood_fill_does_not_escape_through_blocked_diagonal():
    from crevice.grid import GridConnectivity
    _frame, spatial, nodes = diagonal_fixture()
    graph = GridConnectivity(spatial, (0, 0, 0), 3, 0.2)
    exterior = voids._flood_exterior(nodes, (3, 3, 3), edge_allowed=lambda a, b:
                                    graph.allows(a, b, nodes[a].clearance, nodes[b].clearance))
    assert exterior == {(0, 0, 0)}


@pytest.mark.parametrize("required, reachable", [(0.2, True), (0.4, False)])
def test_widest_path_uses_between_node_capacity(required, reachable):
    from crevice.grid import GridConnectivity
    spatial = SpatialIndex((atom(1, (0, 0, 0)),))
    graph = GridConnectivity(spatial, (-3, 2, 0), 6, required)
    left, right = (0, 0, 0), (1, 0, 0)
    nodes = {p: spatial.nearest_surface(graph.position(p))[0] for p in (left, right)}
    assert min(nodes.values()) > 1.9
    widths, _, _ = analysis._widest_paths(nodes, (2, 1, 1), left, 6,
                                         edge_clearance=graph.clearance, min_radius=required)
    assert (right in widths) == reachable
    if reachable:
        assert widths[right] == pytest.approx(0.3)


def test_valid_seed_does_not_teleport_out_of_unsampled_cage():
    frame = StructureFrame(cage())
    assert SpatialIndex(frame.atoms).nearest_surface((0, 0, 0))[0] > 3
    assert find_tunnels(frame, start=(0, 0, 0), spacing=20, min_radius=0.5) == ()


def test_edge_certificates_and_bounded_cache_match_exact_clearance():
    from crevice.grid import GridConnectivity
    spatial = SpatialIndex((atom(1, (0, 0, 0)),))
    graph = GridConnectivity(spatial, (-3, -3, -3), 1, 0.2, cache_limit=2)
    for left, right in [((x, 3, 3), (x + 1, 3, 3)) for x in range(6)]:
        a, b = graph.position(left), graph.position(right)
        ca, cb = spatial.nearest_surface(a)[0], spatial.nearest_surface(b)[0]
        assert graph.allows(left, right, ca, cb) == (spatial.segment_clearance(a, b) >= 0.2)
        assert graph.clearance(right, left) == pytest.approx(spatial.segment_clearance(a, b))
    assert graph.to_dict()["cached_edges"] <= 2
    left, right = (20, 20, 20), (21, 20, 20)
    before = graph.segment_queries
    assert graph.allows(left, right, 20, 20)
    assert graph.segment_queries == before


def test_numpy_segment_coordinates_are_supported():
    spatial = SpatialIndex((atom(1, (0, 0, 0)),))
    assert spatial.segment_clearance(np.array([-3, 2, 0]), np.array([3, 2, 0])) == pytest.approx(0.3)


@pytest.mark.parametrize("wall,spacing,phase,expected", [
    ("closed", 1, 0.5, False), ("wide", 1, 0.5, True),
    ("narrow", 1, 0, True), ("narrow", 1, 0.5, False), ("narrow", 0.125, 0.5, True),
])
def test_wall_resolution_and_grid_phase(wall, spacing, phase, expected):
    from pathlib import Path
    import runpy
    demo = runpy.run_path(str(Path(__file__).resolve().parents[1] / "examples/connectivity_validation.py"))
    result = demo["measure_gate"](wall, spacing, phase)
    assert result["connected"] is expected

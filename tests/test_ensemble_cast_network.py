"""Scientific regression tests for ensemble, focused-void and network behavior."""

import csv
import itertools
import json
import math
from dataclasses import replace

import numpy as np
import pytest

from crevice import (Atom, StructureFrame, pore_profile, void_cast, detect_cavities,
                   analyze_trajectory, build_residue_network, build_cavity_network,
                   profile_distribution, write_profile_distribution_csv,
                   plot_trajectory_profiles, persistent_network, lining_connectivity,
                   write_volume_viewer_bundle, write_void_cast_dx)
from crevice.models import (ChannelPoint, PoreProfile, FrameAnalysis, NetworkNode, NetworkEdge,
                          ResidueInteractionNetwork, VoidComponent, VoidCast)
from crevice.volume_export import cast_grid
from crevice.trajectory import align_frame


def atom(i, xyz, **kw):
    return Atom(i, kw.pop("name", "CA"), kw.pop("resname", "ALA"), kw.pop("chain_id", "A"),
                kw.pop("resid", i), *xyz, element=kw.pop("element", "C"), **kw)


def ring(radius=5):
    return StructureFrame(tuple(atom(1 + i * 32 + j,
        (radius * math.cos(j * math.tau / 32), radius * math.sin(j * math.tau / 32), z))
        for i, z in enumerate(range(-6, 7, 2)) for j in range(32)))


def cage(shift=0):
    coords = set()
    for axis in range(3):
        for side in (-5, 5):
            for a, b in itertools.product(range(-5, 6, 2), repeat=2):
                p = [a, b]
                p.insert(axis, side)
                coords.add((p[0] + shift, p[1], p[2]))
    return tuple(atom(i, p) for i, p in enumerate(sorted(coords), 1))


def frame_profile(radius, origin=(0, 0, 0), ts=(-2, 0, 2)):
    points = tuple(ChannelPoint(i, (origin[0], origin[1], origin[2] + t), radius, radius, t)
                   for i, t in enumerate(ts))
    return FrameAnalysis(0, profile=PoreProfile(origin, (0, 0, 1), points))


def graph(weights=(0.8, 0.8, 0.1)):
    nodes = tuple(NetworkNode(n, "residue", n, metadata={"group": g})
                  for n, g in (("a", "TM1"), ("b", "TM2"), ("c", "TM3"), ("isolated", "loop")))
    edges = (NetworkEdge("a", "b", "contact", 3, weights[0]),
             NetworkEdge("b", "c", "contact", 3, weights[1]),
             NetworkEdge("a", "c", "contact", 10, weights[2]))
    return ResidueInteractionNetwork(nodes, edges)


def test_exact_frame_quantiles_and_csv(tmp_path):
    distribution = profile_distribution([frame_profile(1), frame_profile(3)], samples=3)
    for row in distribution["rows"]:
        assert row["mean_radius"] == 2
        assert row["lower_radius"] == pytest.approx(1.2)
        assert row["upper_radius"] == pytest.approx(2.8)
        assert row["n_frames"] == 2
    write_profile_distribution_csv(distribution, tmp_path / "band.csv")
    with (tmp_path / "band.csv").open() as handle:
        assert len(list(csv.DictReader(handle))) == 3


def test_common_coordinate_grid_has_no_extrapolation():
    frames = [frame_profile(1), frame_profile(3, origin=(0, 0, 1))]
    rows = profile_distribution(frames, samples=3)["rows"]
    assert rows[0]["coordinate"] == -1
    assert rows[-1]["coordinate"] == 2
    rows = profile_distribution(frames, coordinate_grid=(-5, -2, 0, 3, 5))["rows"]
    assert [r["n_frames"] for r in rows] == [0, 1, 2, 1, 0]
    assert rows[0]["mean_radius"] is None


@pytest.mark.parametrize("quantiles", [(0.9, 0.1), (-0.1, 0.9), (0, math.nan)])
def test_invalid_quantiles_rejected(quantiles):
    with pytest.raises(ValueError):
        profile_distribution([frame_profile(1)], quantiles=quantiles)


def test_axis_and_probe_mismatch_rejected():
    frame = frame_profile(1)
    for profile in (replace(frame.profile, probe_radius=1),
                    replace(frame.profile, axis_direction=(1, 0, 0))):
        with pytest.raises(ValueError):
            profile_distribution([frame, replace(frame, profile=profile)])


def test_rigid_alignment_and_fixed_reference():
    reference = ring()
    mobile = replace(reference, atoms=tuple(replace(a, x=-a.z + 8, y=a.y - 5, z=a.x + 3)
                                             for a in reference.atoms))
    aligned = align_frame(mobile, reference)
    assert np.asarray([a.coord for a in aligned.atoms]) == pytest.approx(
        np.asarray([a.coord for a in reference.atoms]), abs=1e-6)
    result = analyze_trajectory([reference, mobile], align=True,
                               profile_kwargs={"axis": "auto", "search_radius": 0, "samples": 9})
    rows = profile_distribution(result)["rows"]
    assert all(r["upper_radius"] - r["lower_radius"] < 1e-6 for r in rows)


def test_distribution_and_volume_figures_are_nonblank(tmp_path):
    from PIL import Image
    figure = tmp_path / "band.png"
    frames = [frame_profile(1), frame_profile(3)]
    plot_trajectory_profiles(frames, figure, dpi=90)
    image = np.asarray(Image.open(figure).convert("RGB"))
    teal = (image[:, :, 1] > image[:, :, 0] + 10) & (image[:, :, 2] > image[:, :, 0] + 10)
    assert teal.sum() > 1000
    plot_trajectory_profiles(frames, tmp_path / "band.svg")
    plot_trajectory_profiles(frames, tmp_path / "timeseries.png", view="timeseries", dpi=80)
    plot_trajectory_profiles(frames, tmp_path / "missing.png", coordinate_grid=(-5, -2, 0, 2, 5), dpi=80)


def test_enclosed_cavities_and_sparse_negative():
    frame = StructureFrame(cage())
    cavities = detect_cavities(frame, spacing=1, min_radius=0.5)
    assert len(cavities) == 1
    assert cavities[0].volume > 50
    sparse = StructureFrame(tuple(atom(i, p) for i, p in enumerate(itertools.product((-4, 4), repeat=3))))
    assert detect_cavities(sparse, spacing=1, min_radius=0.5) == ()
    assert detect_cavities(sparse, spacing=1, min_radius=0.5, enclosed_only=False)


def test_two_cavities_focus_and_volume_filters():
    frame = StructureFrame(cage() + tuple(replace(a, serial=a.serial + 1000, resid=a.resid + 1000) for a in cage(16)))
    full = void_cast(frame, mode="cavity", spacing=1, min_radius=0.5)
    assert len(full.components) == 2
    focused = void_cast(frame, mode="cavity", spacing=1, min_radius=0.5,
                        focus_points=[(0, 0, 0)], focus_radius=7)
    assert len(focused.components) == 1
    assert focused.total_volume * 2 == full.total_volume
    assert not void_cast(frame, mode="cavity", spacing=1, min_radius=0.5,
                         max_component_volume=1).components


def test_focus_miss_never_returns_unfocused_cast():
    frame = ring()
    profile = pore_profile(frame, search_radius=0, samples=13)
    far_profile = replace(profile, points=tuple(replace(p, position=(100, 100, p.t)) for p in profile.points))
    assert not void_cast(frame, mode="channel", spacing=1, profile=far_profile).points
    assert not void_cast(frame, mode="all", spacing=1, focus_points=[(100, 100, 100)]).points


def test_sampling_does_not_change_volume_topology_or_maps():
    frame = ring()
    profile = pore_profile(frame, search_radius=0, samples=13)
    full = void_cast(frame, mode="channel", spacing=1, profile=profile)
    thinned = void_cast(frame, mode="channel", spacing=1, profile=profile, max_export_points=5)
    assert len(thinned.points) == 5
    assert thinned.total_volume == full.total_volume > 0
    assert len(thinned.components) == len(full.components)
    assert np.array_equal(cast_grid(thinned)[0], cast_grid(full)[0])
    assert all(math.hypot(*p.position[:2]) <= 3.3 + 1e-6 for p in full.points)


def test_map_preserves_branches_and_components(tmp_path):
    points = tuple(ChannelPoint(i, p, 1, 1, p[2]) for i, p in enumerate(
        [(0, 0, z) for z in range(4)] + [(x, 0, 2) for x in (1, 2, 3)] + [(8, 0, 0)]))
    component = VoidComponent(1, "test", (0, 0, 0), 8, points)
    cast = VoidCast((component,), points[:1], 1, 0.4, "all")
    grid, origin = cast_grid(cast)
    assert grid.sum() == 8
    from scipy.ndimage import label
    assert label(grid)[1] == 2
    path = write_void_cast_dx(cast, tmp_path / "cast.dx")
    text = path.read_text()
    values = text.split("data follows\n")[1].split("attribute")[0].split()
    assert np.array([float(v) for v in values]).reshape(grid.shape).tolist() == grid.tolist()
    assert origin.tolist() == [-1, -1, -1]


def test_all_viewers_share_map(tmp_path):
    frame = ring()
    cast = void_cast(frame, mode="channel", axis="z", spacing=1,
                     profile=pore_profile(frame, search_radius=0, samples=13))
    manifest = write_volume_viewer_bundle(cast, structure_path=tmp_path / "protein.pdb",
                                         output_dir=tmp_path / "with spaces")
    metadata = json.loads(__import__("pathlib").Path(manifest["volume_metadata_json"]).read_text())
    assert metadata["boundary_validation"] == "not_continuously_validated"
    assert metadata["connectivity"]["edge_rule"] == "continuous_atom_sphere_clearance"
    for path in manifest.values():
        assert __import__("pathlib").Path(path).exists()
    pml = __import__("pathlib").Path(manifest["volume_pml"]).read_text()
    compile(pml.split("python\n", 1)[1].split("python end")[0], "scene", "exec")
    assert 'cmd.load_cgo' in pml
    assert metadata['pymol_representation']=='standalone_triangle_mesh'
    assert "Isosurface 0.5" in __import__("pathlib").Path(manifest["volume_tcl"]).read_text()
    assert "surfaceSmoothing false" in __import__("pathlib").Path(manifest["volume_cxc"]).read_text()


def test_sidechain_evidence_and_full_identity():
    frame = StructureFrame((atom(1, (0, 0, 0), name="NZ", resname="LYS", element="N"),
                            atom(2, (3, 0, 0), name="OD1", resname="ASP", element="O")))
    network = build_residue_network(frame, residue_groups={"A:LYS1": "TM1", "A:ASP2": "TM2"})
    assert network.edges[0].interaction == "salt_bridge_candidate"
    assert network.edges[0].metadata["charged_atom_pairs"]
    with pytest.raises(ValueError):
        build_residue_network(frame, residue_groups={"B:LYS1": "TM1"})


def test_wide_pore_network_keeps_lining():
    frame = ring(radius=10)
    profile = pore_profile(frame, search_radius=0, samples=7)
    network = build_cavity_network(frame, profile)
    assert any(e.interaction in {"region_lining", "region_bottleneck"} for e in network.edges)


def test_weighted_paths_ignore_geometric_distance_and_region_shortcuts():
    network = graph()
    region = NetworkNode("region:channel", "region", "pore")
    network = replace(network, nodes=network.nodes + (region,), edges=network.edges +
                      (NetworkEdge(region.id, "a", "region_lining", 0),
                       NetworkEdge(region.id, "c", "region_nearby", 1, 100)))
    report = lining_connectivity(network)
    rows = {r["residue"]: r for r in report["residues"]}
    assert rows["c"]["path"] == ["a", "b", "c"]
    assert rows["c"]["cost_to_lining"] == pytest.approx(2.5)
    assert not rows["isolated"]["reachable"]
    assert len(report["group_contacts"]) == 3
    removed = lining_connectivity(network, removed_residues=["b"])
    c = next(r for r in removed["residues"] if r["residue"] == "c")
    assert c["cost_to_lining"] == 10


def test_articulation_removal_disconnects_distal_node():
    network = graph()
    network = replace(network, edges=network.edges[:2])
    before = lining_connectivity(network, lining_residues=["a"])
    assert next(r for r in before["residues"] if r["residue"] == "b")["articulation_point"]
    after = lining_connectivity(network, lining_residues=["a"], removed_residues=["b"])
    assert not next(r for r in after["residues"] if r["residue"] == "c")["reachable"]


def test_contact_occupancy_missing_nodes_and_threshold():
    a = graph()
    b = replace(a, edges=a.edges[1:])
    c = replace(a, nodes=tuple(n for n in a.nodes if n.id != "a"), edges=a.edges[1:2])
    pooled = persistent_network([a, b, c], min_occupancy=0.5)
    ab = next(e for e in pooled.edges if {e.source, e.target} == {"a", "b"})
    assert ab.weight == 0.5
    assert ab.metadata["eligible_frames"] == 2
    assert not any({e.source, e.target} == {"a", "b"} for e in persistent_network([a, b, c], min_occupancy=0.6).edges)
    with pytest.raises(ValueError):
        persistent_network([replace(a, metadata={"cutoff": 3}), replace(b, metadata={"cutoff": 4})])


def test_lining_persistence_combines_lining_and_bottleneck_roles():
    network = graph()
    region = NetworkNode("region:channel", "region", "pore")
    frames = [replace(network, nodes=network.nodes + (region,), edges=network.edges +
                      (NetworkEdge(region.id, "a", role, 0),))
              for role in ("region_lining", "region_bottleneck")]
    pooled = persistent_network(frames, min_occupancy=0.75)
    assert lining_connectivity(pooled)["lining_residues"] == ["a"]

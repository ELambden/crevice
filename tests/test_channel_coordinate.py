"""Analytic/synthetic contracts for channel assignment, separate from real data."""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
from dataclasses import replace
import math

import numpy as np
import pytest

from crevice import Atom, StructureFrame, pore_profile, void_cast
from crevice.geometry import dot, distance
from crevice.spatial import SpatialIndex
from crevice.volume_export import cast_grid, write_void_cast_dx


def rectangular_channel(*, half_x=4, half_y=6, half_length=10, curved=False, capped=False):
    coords = set()
    for z in range(-half_length, half_length+1):
        shift = 0.9*math.cos(z/half_length*math.pi/2) if curved else 0.0
        for x in range(-half_x, half_x+1):
            coords.update(((x+shift, -half_y, z), (x+shift, half_y, z)))
        for y in range(-half_y, half_y+1):
            coords.update(((-half_x+shift, y, z), (half_x+shift, y, z)))
    if capped:
        coords.update((x, y, half_length) for x in range(-half_x, half_x+1)
                      for y in range(-half_y, half_y+1))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C")
                                for i, p in enumerate(sorted(coords), 1)))


def transform(frame):
    from scipy.spatial.transform import Rotation
    rotation = Rotation.from_rotvec([0.73, -0.51, 0.28]).as_matrix()
    translation = np.array([13.1, -29.4, 8.6])
    coords = np.array([a.coord for a in frame.atoms]) @ rotation.T + translation
    moved = replace(frame, atoms=tuple(replace(a, x=float(p[0]), y=float(p[1]), z=float(p[2]))
                                       for a, p in zip(frame.atoms, coords)))
    return moved, rotation, translation


@pytest.fixture(scope="module")
def channel():
    frame = rectangular_channel()
    return frame, pore_profile(frame, enclosure_radius=0.8, samples=41, search_radius=8)


def test_automatic_axis_is_rigid_transform_equivariant(channel):
    frame, profile = channel
    moved, rotation, translation = transform(frame)
    other = pore_profile(moved, enclosure_radius=0.8, samples=41, search_radius=8)
    assert np.allclose(other.axis_direction, rotation @ profile.axis_direction, atol=1e-9)
    assert np.allclose([p.radius for p in other.points], [p.radius for p in profile.points], atol=1e-8)
    assert np.allclose([p.position for p in other.points],
                       np.array([p.position for p in profile.points]) @ rotation.T + translation, atol=1e-8)
    assert profile.metadata["axis"] == "auto:connected-sections"


def test_nonlongest_principal_direction_can_be_the_channel():
    frame = rectangular_channel(half_y=10, half_length=5)
    profile = pore_profile(frame, enclosure_radius=0.8, samples=21, search_radius=12)
    assert abs(dot(profile.axis_direction, (0, 0, 1))) > 0.999
    assert profile.metadata["axis_candidates"][0]["through_span_A"] == 0


def test_radius_is_inscribed_clearance_not_a_remote_width(channel):
    frame, profile = channel
    # Carbon walls at +/-4 give a 2.3 A inscribed sphere; the other dimension
    # is wider. Measure actual 3D atom surfaces, without averaging wall widths.
    middle = min(profile.points, key=lambda p: abs(p.t))
    assert middle.radius == pytest.approx(4-1.7, abs=0.04)
    assert max(p.radius for p in profile.points if abs(p.position[2]) <= 10) < 2.5
    assert min(profile.metadata["segment_clearances_A"]) >= 0.8
    assert profile.metadata["continuous_bottleneck_radius_A"] <= profile.min_radius
    assert profile.metadata["arc_lengths_A"][-1] == pytest.approx(profile.length)


def test_cast_fills_non_circular_section_and_rotates_with_structure(channel, tmp_path):
    frame, profile = channel
    options = dict(spacing=.5, min_radius=0, max_grid_points=200_000)
    cast = void_cast(frame, profile=profile, **options)
    assert cast.component_count == 1
    assert cast.metadata["status"] == "resolved"
    assert cast.metadata["export_mode"] == "connected_section_fill"
    # Near-wall corners lie outside every old centerline inscribed sphere.
    assert any(abs(p.position[0]) > 1.0 and abs(p.position[1]) > 3.2
               and abs(p.position[2]) < 1 for p in cast.components[0].points)
    assert all(abs(p.position[0]) < 2.5 and abs(p.position[1]) < 4.5 for p in cast.components[0].points if abs(p.position[2]) <= 10)
    spatial = SpatialIndex(frame.atoms)
    assert min(spatial.nearest_surface(p.position)[0] for c in cast.components for p in c.points) >= 0
    moved, rotation, translation = transform(frame)
    other_profile = pore_profile(moved, enclosure_radius=0.8, samples=41, search_radius=8)
    other = void_cast(moved, profile=other_profile, **options)
    assert other.total_volume == cast.total_volume
    assert other.component_count == cast.component_count
    assert np.array_equal(cast_grid(cast)[0], cast_grid(other)[0])
    dx = write_void_cast_dx(other, tmp_path / "rotated.dx").read_text()
    deltas = np.array([[float(v) for v in line.split()[1:]] for line in dx.splitlines() if line.startswith("delta ")])
    assert np.allclose(deltas, np.array(other.metadata["grid_basis"])*other.spacing)
    assert not np.allclose(deltas, np.eye(3)*other.spacing)
    thinned = void_cast(frame, profile=profile, max_export_points=3, **options)
    assert thinned.total_volume == cast.total_volume
    assert cast_grid(thinned)[0].sum() == sum(c.point_count for c in cast.components)


def test_curved_channel_continuity_and_solvent_termination():
    frame = rectangular_channel(curved=True)
    profile = pore_profile(frame, enclosure_radius=0.8, search_radius=8, samples=41)
    spatial = SpatialIndex(frame.atoms)
    assert profile.method == "axial-connected"
    assert all(spatial.segment_clearance(a.position, b.position) > 0.8
               for a, b in zip(profile.points, profile.points[1:]))
    assert max(p.position[0] for p in profile.points)-min(p.position[0] for p in profile.points) > .4
    # Mouths lie past the end atom centers, where inflated end spheres cease to enclose.
    assert all(abs(p.position[2]) < 12.5 for p in profile.points)
    assert max(abs(p.position[2]) for p in profile.points) > 10.0
    assert max(p.radius for p in profile.points if abs(p.position[2]) <= 10) < 2.5


def test_wall_and_open_space_cannot_be_reported_as_a_through_channel():
    with pytest.raises(ValueError, match="No connected channel"):
        pore_profile(rectangular_channel(capped=True), enclosure_radius=0.8, axis="z", search_radius=8)
    frame = StructureFrame((Atom(1, "C", "ALA", "A", 1, 0, 0, 0, "C"),))
    with pytest.raises(ValueError, match="No connected channel"):
        pore_profile(frame, enclosure_radius=0.8)


def test_explicit_seed_must_be_feasible_and_selects_one_parallel_channel():
    one = rectangular_channel(half_x=3, half_y=3)
    both = StructureFrame(tuple(replace(a, x=a.x-5) for a in one.atoms) +
                          tuple(replace(a, x=a.x+5, serial=a.serial+10000, resid=a.resid+10000)
                                for a in one.atoms))
    with pytest.raises(ValueError, match="Ambiguous"):
        pore_profile(both, enclosure_radius=0.8, axis="z", search_radius=10)
    profile = pore_profile(both, enclosure_radius=0.8, axis="z", origin=(-5, 0, 0), search_radius=6)
    assert all(abs(p.position[0]+5) < .5 for p in profile.points)
    with pytest.raises(ValueError, match="origin must have clearance"):
        pore_profile(one, enclosure_radius=0.8, axis="z", origin=one.atoms[0].coord)


def test_no_silent_coarsening_and_atom_selection_mismatch(channel):
    frame, profile = channel
    with pytest.raises(ValueError, match="No automatic coarsening"):
        void_cast(frame, profile=profile, spacing=.25, max_grid_points=100)
    with pytest.raises(ValueError, match="same atom selection"):
        void_cast(frame, profile=profile, include_hetero=False)


@pytest.mark.parametrize("options", [{"section_spacing": 0}, {"section_spacing": float("nan")},
                                     {"enclosure_radius": -1}, {"origin": (float("inf"), 0, 0)}])
def test_invalid_channel_controls(channel, options):
    with pytest.raises(ValueError):
        pore_profile(channel[0], **options)


def test_fixed_axis_scan_is_explicitly_unvalidated(channel):
    p = pore_profile(channel[0], enclosure_radius=0.8, axis="z", search_radius=0, samples=7)
    assert len(p.points) == 7
    assert p.metadata["status"] == "unvalidated_axis_scan"


def test_cli_uses_auto_axis_and_map_scene(tmp_path):
    from crevice.cli import main, build_parser
    import json
    frame = rectangular_channel(half_length=7)
    path = tmp_path / "wall.pdb"
    path.write_text("\n".join(
        f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
        for a in frame.atoms) + "\nEND\n")
    args = build_parser().parse_args(["publish", str(path), "--out-dir", str(tmp_path / "out")])
    assert args.axis == "auto"
    assert main(["publish", str(path), "--enclosure-radius", "0.8", "--out-dir", str(tmp_path / "out"), "--search-radius", "8",
                 "--samples", "21", "--section-spacing", "0.5", "--enclosure-radius", "0.8",
                 "--skip-cavities", "--skip-tunnels", "--skip-network", "--dpi", "60"]) == 0
    root = tmp_path / "out"
    scene = (root / "wall_pore_cast.pml").read_text()
    assert "wall_volume.dx" in scene and "cmd.load_cgo" in scene
    assert (root / "wall_pymol_mesh.npz").is_file()
    assert "alter " not in scene and "show surface" not in scene
    profile = json.loads((root / "wall_profile.json").read_text())
    assert profile["method"] == "axial-connected"


def test_discontinuous_enclosure_cannot_select_just_one_channel_fragment():
    frame = rectangular_channel()
    cuffs = replace(frame, atoms=tuple(a for a in frame.atoms if abs(a.z) >= 4))
    with pytest.raises(ValueError, match="No connected channel"):
        pore_profile(cuffs, enclosure_radius=0.8, axis="z", search_radius=8)


def test_rotated_connectivity_measures_the_actual_world_segment():
    from crevice.grid import GridConnectivity
    spatial = SpatialIndex((Atom(1, 'C', 'ALA', 'A', 1, 0, 0, 0, 'C'),))
    basis = ((0, 1, 0), (-1, 0, 0), (0, 0, 1))
    grid = GridConnectivity(spatial, (0, -3, 0), 6, .2, basis=basis)
    assert grid.position((1, 0, 0)) == (0, 3, 0)
    assert not grid.allows((0, 0, 0), (1, 0, 0), 1.3, 1.3)


def pocket_channel(*, slit=2.4999):
    """Straight lumen in a solid carbon lattice, plus four small wall pockets.

    Lattice sites every 2 A (a solid wall for a 0.8 A enclosure probe) are
    removed within 3 A of the z axis (lumen) and within 2.6 A of four short
    tubes in the z = 0 plane running from the axis to (+/-6, 0, 0) and
    (0, +/-6, 0). A pair of atoms at 4 A from the axis, 2*slit apart, closes
    each tube in the z = 0 plane only: their inflated 2.5 A disks overlap there
    (slit < 2.5) but not a few hundredths of an angstrom off that plane. So at
    z = 0 each pocket is a separate enclosed section of <= 4 grid points
    (0.5 A grid; clearance ~1.18 A, lumen ~3.0 A) that straight segments from
    the lumen centres just above and below still reach with more than 0.8 A
    clearance. The pocket at -x has the smallest grid index, so section
    discovery lists it before the lumen. The layout has four-fold symmetry and
    the first atom fixes the transverse basis to (+x, +y).
    """
    tubes = ((-1, 0), (1, 0), (0, -1), (0, 1))

    def tube_distance(p, d):
        s = max(0.0, min(6.0, p[0]*d[0] + p[1]*d[1]))
        return math.dist(p, (s*d[0], s*d[1], 0.0))

    coords = [(2.0*i, 2.0*j, 2.0*k) for i in range(-4, 5) for j in range(-4, 5) for k in range(-4, 5)]
    coords = [p for p in coords if math.hypot(p[0], p[1]) >= 3.0
              and all(tube_distance(p, d) >= 2.6 for d in tubes)]
    coords += [(4*dx - dy*sign*slit, 4*dy + dx*sign*slit, 0.0) for dx, dy in tubes for sign in (-1, 1)]
    coords.sort(key=lambda p: (p != (8.0, 0.0, -8.0), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


POCKET_OPTIONS = dict(axis="z", search_radius=7.0, section_spacing=0.5, enclosure_radius=0.8)


def test_pocket_fixture_offers_a_tied_wall_pocket_before_the_lumen():
    # Guard: the fixture must present exactly the tie that routed profiles
    # through one-grid-point wall interstices in 1GRM before the widest-path rule.
    from crevice.sections import SectionGeometry
    frame = pocket_channel()
    geometry = SectionGeometry(frame.atoms, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), spacing=0.5, extent=7.0, probe=0.8)
    assert geometry.u == (1.0, 0.0, 0.0)
    sections = [geometry.refine(s, 4) for s in geometry.sections(0.0)]
    first, lumen = sections[0], max(sections, key=lambda s: s.clearance)
    assert len(sections) == 5 and len(first.indices) <= 4 and first.clearance < 1.3 < 2.9 < lumen.clearance
    for dz in (-0.2, 0.2):
        neighbour = max((geometry.refine(s, 4) for s in geometry.sections(dz)), key=lambda s: s.clearance)
        assert geometry.spatial.segment_clearance(neighbour.center, first.center) > 0.8


@pytest.mark.parametrize("origin", [None, (0.0, 0.0, 0.0)])
def test_profile_follows_the_lumen_not_a_tied_wall_pocket(origin):
    # Before the fix: without a seed, 25/81/101 samples reported the 1.18 A
    # pocket as the minimum; with the seed the discovery path went through the
    # pocket, failed the seed-segment test and the channel was refused.
    frame = pocket_channel()
    profile = pore_profile(frame, origin=origin, samples=25, **POCKET_OPTIONS)
    middle = profile.points[len(profile.points)//2]
    assert abs(middle.t) < 1e-6  # a sample lies in the z = 0 pocket plane
    assert math.hypot(*middle.position[:2]) < 0.1
    assert profile.min_radius == pytest.approx(2.3, abs=0.01)
    assert min(profile.metadata["section_areas_A2"]) > 4 * 0.5**2
    assert profile.metadata["continuous_bottleneck_radius_A"] == pytest.approx(2.3, abs=1e-6)


def test_minimum_radius_is_stable_across_sample_counts():
    # Odd counts put a sample on the pocket plane, even counts do not. The
    # tolerance is the final refinement step (section_spacing / 2**4): a
    # sampled minimum may land anywhere within one refinement move of the
    # exact in-plane maximum, and axial sampling of this wall changes the
    # sampled minimum by less than that.
    frame = pocket_channel()
    tolerance = 0.5 / 2**4
    minima, bottlenecks = [], []
    for samples in (81, 101, 118, 150):
        profile = pore_profile(frame, samples=samples, **POCKET_OPTIONS)
        assert min(profile.metadata["section_areas_A2"]) > 4 * 0.5**2
        minima.append(profile.min_radius)
        bottlenecks.append(profile.metadata["continuous_bottleneck_radius_A"])
    assert max(minima) - min(minima) <= tolerance
    assert max(bottlenecks) - min(bottlenecks) <= tolerance

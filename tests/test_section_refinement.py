"""Synthetic checks of section-centre refinement (``SectionGeometry.refine``).

Before 0.4.0 the centre took one hill-climbing move per halving step (four
steps) and could stop short of the in-plane clearance maximum: on 1GRM by up
to 0.03 Å within 0.1 Å of the centre, which made the minimum radius depend
slightly on the sample count. The refined centre is now an exact local maximum: the
three-atom clearance vertex at the top of the hill that holds the grid's best
node. Synthetic geometry only; not a real-system validation.
"""
import math

import numpy as np
import pytest

from crevice.models import Atom
from crevice.radii import atom_vdw_radius
from crevice.sections import REFINEMENT_TOLERANCE_A, SectionGeometry, connected_profile

CENTRE = (0.137, -0.081)
# An irregular ring: distances from CENTRE vary, so the maximum is a three-atom vertex off the grid.
RING = ((0, 3.7), (37, 4.3), (71, 3.9), (110, 4.4), (146, 3.8), (183, 4.2), (219, 4.0),
        (255, 4.5), (291, 3.75), (326, 4.25))


def _ring(z=0.0, first=1):
    return tuple(Atom(first + k, "C", "ALA", "A", 1, CENTRE[0] + d*math.cos(math.radians(a)),
                      CENTRE[1] + d*math.sin(math.radians(a)), z, "C") for k, (a, d) in enumerate(RING))


def _brute_max(atoms, z=0.0):
    """In-plane maximum clearance by brute force (0.005 Å grid, then 2e-5 and 2e-7 Å)."""
    xyz = np.array([a.coord for a in atoms])
    radii = np.array([atom_vdw_radius(a) for a in atoms])
    best = np.array([CENTRE[0], CENTRE[1], z])
    for half, h in ((1.5, 0.005), (0.006, 0.00002), (0.00004, 2e-7)):
        g = np.arange(-half, half + h/2, h)
        a, b = np.meshgrid(g, g, indexing="ij")
        pts = best + np.stack([a.ravel(), b.ravel(), np.zeros(a.size)], axis=1)
        vals = (np.linalg.norm(pts[:, None, :] - xyz[None], axis=2) - radii).min(axis=1)
        best = pts[vals.argmax()]
    return float(vals.max())


def test_refined_centre_is_the_exact_in_plane_maximum():
    atoms = _ring()
    geometry = SectionGeometry(atoms, (0, 0, 0), (0, 0, 1), spacing=0.25, extent=6.0, probe=0.5)
    (section,) = geometry.sections(0.0)
    refined = geometry.refine(section, 4)
    expected = _brute_max(atoms)
    assert refined.clearance > section.clearance + 0.005          # the grid node is clearly off the maximum
    assert refined.clearance == pytest.approx(expected, abs=1e-6)
    assert refined.clearance >= expected - 1e-9
    # At the maximum three atom surfaces are equidistant (a clearance vertex).
    surface = sorted(math.dist(refined.center, a.coord) - atom_vdw_radius(a) for a in atoms)
    assert surface[2] - surface[0] < 1e-9
    assert refined.t == section.t and refined.indices == section.indices
    # Deterministic and memoised.
    assert geometry.refine(section, 4) == refined


def test_zero_refinement_steps_keep_the_grid_node():
    geometry = SectionGeometry(_ring(), (0, 0, 0), (0, 0, 1), spacing=0.25, extent=6.0, probe=0.5)
    (section,) = geometry.sections(0.0)
    assert geometry.refine(section, 0) == section


def test_profile_records_the_refinement_and_reaches_the_in_plane_maximum():
    atoms = tuple(a for k, z in enumerate((-6.0, -3.0, 0.0, 3.0, 6.0)) for a in _ring(z, first=1 + 10*k))
    profiles = [connected_profile(atoms, axis=(0, 0, 1), origin=(0, 0, 0), samples=n, search_radius=6.0,
                                  refinement_steps=4, probe_radius=0.0, section_spacing=0.25, enclosure_radius=0.5)
                for n in (21, 29, 37)]
    record = profiles[0].metadata["section_refinement"]
    assert record["position_tolerance_A"] == REFINEMENT_TOLERANCE_A
    assert record["minimum_halvings"] == 4
    # In the ring planes the radius is the in-plane maximum (other rings are too far to limit it).
    in_plane, checked = _brute_max(_ring()), 0
    for profile in profiles:
        ring_points = [p for p in profile.points if min(abs(p.t - z) for z in (-3.0, 0.0, 3.0)) < 1e-9]
        for point in ring_points:
            assert point.radius == pytest.approx(in_plane, abs=1e-6)
            checked += 1
    assert checked


def _peanut(dx=0.0, dy=0.0):
    """Two overlapping circles of carbon atoms (lobes) joined by a neck, shifted by ``(dx, dy)``.

    Every atom of a lobe is on its circle, so the lobe's clearance maximum is
    its centre, exactly ``R - 1.7`` Å: 2.70 Å for lobe A at the origin, 2.73 Å
    for lobe B at (6.2, 0.05). Atoms of one circle inside the other are left out.
    """
    lobes = ((0.0, 0.0, 4.4), (6.2, 0.05, 4.43))
    points = []
    for cx, cy, radius in lobes:
        count = int(2*math.pi*radius)
        for k in range(count):
            angle = 2*math.pi*k/count + 0.3
            x, y = cx + radius*math.cos(angle), cy + radius*math.sin(angle)
            if all(math.hypot(x - ox, y - oy) > other - 0.01 for ox, oy, other in lobes if (ox, oy) != (cx, cy)):
                points.append((x, y))
    return tuple(Atom(k + 1, "C", "ALA", "A", 1, x + dx, y + dy, 0.0, "C") for k, (x, y) in enumerate(points))


def _peanut_section(dx, dy):
    geometry = SectionGeometry(_peanut(dx, dy), (0, 0, 0), (0, 0, 1), spacing=0.25, extent=12.0, probe=0.5)
    (section,) = geometry.sections(0.0)
    return geometry, section


OFFSETS = [(dx, dy) for dx in (0.0, 0.05, 0.1, 0.15, 0.2) for dy in (0.0, 0.05, 0.1, 0.15, 0.2)]


def test_refinement_reaches_the_maximum_of_the_grid_centres_lobe_at_every_grid_offset():
    """The refined centre is the exact maximum of the lobe that holds the grid's best node.

    At some offsets that node lies on the lower lobe (2.70 A); refinement does
    not jump to the higher lobe (2.73 A), by definition (see
    ``SectionGeometry.refine``). Both cases occur here.
    """
    assert atom_vdw_radius(_peanut()[0]) == 1.7
    lobes = {"A": ((0.0, 0.0), 4.4 - 1.7), "B": ((6.2, 0.05), 4.43 - 1.7)}
    seen = set()
    for dx, dy in OFFSETS:
        geometry, section = _peanut_section(dx, dy)
        lobe = "A" if section.center[0] - dx < 3.1 else "B"
        (cx, cy), expected = lobes[lobe]
        refined = geometry.refine(section, 4)
        assert refined.clearance == pytest.approx(expected, abs=1e-9)
        assert math.dist(refined.center, (cx + dx, cy + dy, 0.0)) < 1e-6
        seen.add(lobe)
    assert seen == {"A", "B"}


def _is_clearance_vertex(xyz, radii, centre, u, v, clearance):
    """True when three atoms limit ``centre`` and their in-plane projections enclose it (a strict local maximum)."""
    import itertools
    surface = np.linalg.norm(xyz - centre, axis=1) - radii
    active = np.flatnonzero(surface <= clearance + 1e-6)
    projected = np.stack([(xyz[active] - centre) @ u, (xyz[active] - centre) @ v], axis=1)
    for a, b, c in itertools.combinations(range(len(active)), 3):
        corners = projected[[a, b, c]]
        matrix = (corners[1:] - corners[0]).T
        if abs(np.linalg.det(matrix)) > 1e-12:
            weights = np.linalg.solve(matrix, -corners[0])
            if weights.min() >= -1e-9 and weights.sum() <= 1 + 1e-9:
                return True
    return False


def _1grm_profile():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from workspace_paths import pdb_cache_dir
    from crevice.channels import pore_profile
    from crevice.parser import load_structure
    source = pdb_cache_dir() / "1GRM.cif"
    if not source.is_file():
        pytest.skip(f"{source} not cached")
    frame = load_structure(str(source), parser="gemmi")
    atoms = frame.selected_atoms(include_hydrogen=False, include_hetero=True)
    xyz = np.array([a.coord for a in atoms])
    radii = np.array([atom_vdw_radius(a) for a in atoms])
    profile = pore_profile(frame, axis="auto", samples=101, section_spacing=0.25, enclosure_radius=0.8)
    u, v = (np.array(w) for w in profile.metadata["section_basis"])
    stalled = [p.index for p in profile.points
               if not _is_clearance_vertex(xyz, radii, np.array(p.position), u, v, p.raw_clearance)]
    return profile, stalled


def test_every_1grm_section_centre_is_an_exact_clearance_vertex():
    """Real structure (optional; skipped without the cached mmCIF), checked for the software property only."""
    profile, stalled = _1grm_profile()
    assert len(profile.points) == 101 and stalled == []


def test_without_the_linearised_move_a_1grm_centre_stalls_on_a_ridge(monkeypatch):
    """Negative control: eight fixed directions alone stop short of a vertex on a two-atom ridge."""
    monkeypatch.setattr(SectionGeometry, "_linear_move", lambda self, *args: None)
    _, stalled = _1grm_profile()
    assert stalled


@pytest.mark.parametrize("samples", [153, 159, 189])
def test_4hfi_profile_does_not_detour_through_a_wall_pocket(samples):
    """Optional (curated 4HFI benchmark input and HOLE's simple.rad; skipped when absent).

    At these sample counts one plane of the 4HFI pore grazes the search square
    and a wall pocket 5 A off the axis in that plane can be joined to its
    neighbours. Without window bridging the only continuous path detoured
    through the pocket (minimum radius 1.84 A instead of 2.59 A); a path that
    much narrower than the discovery path now competes with the bridged one.
    """
    import json
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from workspace_paths import workspace_path
    from crevice.channels import pore_profile
    from crevice.parser import load_structure
    from crevice.radii import use_radii
    entry = workspace_path("benchmarks/channels/4HFI")
    radii = workspace_path("tools/hole2/share/hole2/rad/simple.rad")
    if not (entry / "curated.pdb").is_file() or not radii.is_file():
        pytest.skip("curated 4HFI benchmark input or HOLE simple.rad not available")
    axis = json.loads((entry / "axis.json").read_text())
    with use_radii(str(radii)):
        frame = load_structure(str(entry / "curated.pdb"), radii=str(radii))
        profile = pore_profile(frame, axis=tuple(axis["axis_unit"]), origin=tuple(axis["origin"]), samples=samples,
                               section_spacing=0.25, enclosure_radius=1.8, radii=str(radii))
    direction = np.array(axis["axis_unit"])
    steps = np.diff(np.array([p.position for p in profile.points]), axis=0)
    lateral = np.linalg.norm(steps - np.outer(steps @ direction, direction), axis=1)
    assert lateral.max() < 0.5
    assert profile.min_radius > 2.5
    assert profile.metadata["window_bridged_samples"]

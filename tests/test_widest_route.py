"""Synthetic contracts for the widest-route choice of the adaptive tracer's centres.

Software/synthetic tests only; they say nothing about any real channel. The fixture is a solid 2 Å carbon
lattice with a wide tube carved along z; for |z| <= 2 a fin of atoms at x = -1 Å splits the tube into a narrow
lobe (x < -1) and a wide branch (x > -1), joined in the plane through a gap at the fin's end and above and below
the fin. The seed lies in the narrow lobe. The straight in-plane segment between the two lobes' tops crosses the
fin, so a per-plane rule cannot tell that the wide branch belongs to the same channel; the route across planes can.
"""
import math

import pytest

from crevice import Atom, StructureFrame, pore_profile
from crevice.adaptive_channel import _widest_route
from crevice.spatial import SpatialIndex


def _fin_tube():
    coords = []
    for i in range(-6, 7):
        for j in range(-6, 7):
            for k in range(-5, 6):
                p = (2.0 * i, 2.0 * j, 2.0 * k)
                if math.hypot(p[0], p[1]) >= 6.7:
                    coords.append(p)
    for z in (-2.0, -1.0, 0.0, 1.0, 2.0):
        coords += [(-1.0, -7.0 + y, z) for y in range(9)]  # the fin, y = -7 ... 1 Å
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


def _best_clearance(frame, z, side):
    """Largest clearance on a 0.1 Å grid of the plane at ``z``, on one side of the fin."""
    index = SpatialIndex(frame.atoms)
    grid = [0.1 * i for i in range(-70, 71)]
    return max(index.nearest_surface((x, y, z))[0] for x in grid if side * (x + 1.0) > 0 for y in grid)


@pytest.fixture(scope="module")
def split():
    frame = _fin_tube()
    profile = pore_profile(frame, axis="z", origin=(-3.85, 0.0, 0.0), search_radius=9.0, section_spacing=0.25,
                           samples=41)
    return frame, profile


def test_a_split_lumen_is_reported_along_its_wider_branch(split):
    frame, profile = split
    narrow, wide = _best_clearance(frame, 0.0, -1), _best_clearance(frame, 0.0, 1)
    assert wide > narrow + 0.9                                    # the fixture: 2.80 Å against 1.84 Å
    # The seed lies in the narrow lobe, yet every sample beside the fin is centred in the wide branch.
    beside = [p for p in profile.points if abs(p.t) <= 2.0]
    assert beside and all(p.position[0] > -1.0 for p in beside)
    # The profile's minimum is the wide branch's clearance (within the grid), not the narrow lobe's.
    assert profile.min_radius == pytest.approx(wide, abs=0.05)
    assert profile.metadata["unjoined_samples"] == [] and profile.metadata["unenclosed_samples"] == []


def test_the_split_lumen_profile_does_not_depend_on_the_sample_count(split):
    frame, profile = split
    other = pore_profile(frame, axis="z", origin=(-3.85, 0.0, 0.0), search_radius=9.0, section_spacing=0.25,
                         samples=57)
    assert (other.points[0].t, other.points[-1].t) == pytest.approx((profile.points[0].t, profile.points[-1].t))
    assert other.min_radius == pytest.approx(profile.min_radius, abs=0.02)


class _C:
    """A stand-in candidate: only ``clearance`` is read by the route choice."""

    def __init__(self, name, clearance=1.0):
        self.name, self.clearance = name, clearance


def _route(layers, weights):
    return [c.name for c in _widest_route(layers, weights)]


def test_the_route_takes_the_branch_whose_narrowest_step_is_widest():
    # Two parallel branches a and b between single-candidate planes s and e; a is wider at every step but one.
    layers = [[_C("s")], [_C("a1"), _C("b1")], [_C("a2"), _C("b2")], [_C("e")]]
    weights = [[[2.0, 1.5]], [[2.5, -math.inf], [-math.inf, 1.6]], [[0.9], [1.4]]]
    assert _route(layers, weights) == ["s", "b1", "b2", "e"]       # a's narrowest step is 0.9, b's 1.4


def test_the_route_is_widest_on_each_side_of_its_bottleneck():
    # Every route passes the step m1 -> a2 (0.5), its bottleneck, so a0 and b0 tie on the narrowest step alone;
    # the wider first step (b0 -> m1, 1.2) decides. x1 is reached from neither start.
    layers = [[_C("a0"), _C("b0")], [_C("m1"), _C("x1")], [_C("a2"), _C("b2")]]
    weights = [[[0.8, -math.inf], [1.2, -math.inf]], [[0.5, 0.3], [-math.inf, -math.inf]]]
    assert _route(layers, weights) == ["b0", "m1", "a2"]
    # Two branches between single-candidate planes: the problem splits there, and each part is solved alone.
    layers = [[_C("a0"), _C("b0")], [_C("m")], [_C("a2"), _C("b2")], [_C("e")]]
    weights = [[[0.7], [0.9]], [[3.0, 2.0]], [[0.4], [2.5]]]
    assert _route(layers, weights) == ["b0", "m", "b2", "e"]


def test_a_single_plane_takes_its_highest_candidate():
    assert _route([[_C("x", 1.0), _C("y", 2.0)]], []) == ["y"]

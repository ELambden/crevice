"""Synthetic contracts for the lateral-branch rule and for ends that already lie in bulk.

Software/synthetic tests only; they say nothing about any real channel. The channel is a 2 Å carbon lattice
with a wide lumen along z, a narrow filter (four atoms, 0.6 Å clearance) below the seed, a side portal leaving
the lumen at z = 2 towards +x, and a hole into the cap above it that ends at z = 8, where a side tube leads
out towards -x (so the upper end is capped). With ``lateral_exits`` the capped upper end moves to the portal
when the portal is clearly wider than the hole and its tube, and stays where it is otherwise. An open end
(the hole through the whole cap) never moves.
"""
import math

import pytest

from crevice import Atom, StructureFrame, adaptive_channel, pore_profile
from crevice.sections import ChannelResolutionError


def _segment_distance(p, a, b):
    ab = [b[i] - a[i] for i in range(3)]
    ap = [p[i] - a[i] for i in range(3)]
    t = max(0.0, min(1.0, sum(x * y for x, y in zip(ab, ap)) / sum(x * x for x in ab)))
    return math.dist(p, [a[i] + t * ab[i] for i in range(3)])


def branched_channel(*, portal_radius, hole_radius, hole_top=8.0, n=5):
    """Lattice sites within 4.5 Å of the axis removed below z = 4 (lumen), within ``hole_radius`` of it up to
    ``hole_top`` (cap hole; 2.1 Å leaves a 1.13 Å clearance, 3.0 Å about 2.3 Å), within ``hole_radius`` of the
    segment (0, 0, hole_top)-(-14, 0, hole_top) (the cap's side tube, when the hole stops inside the block) and
    within ``portal_radius`` of the segment (0, 0, 2)-(14, 0, 2) (side portal; same clearances); filter atoms at
    2.3 Å from the axis at z = -4."""
    coords = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                p = (2.0 * i, 2.0 * j, 2.0 * k)
                r = math.hypot(p[0], p[1])
                if (p[2] <= 4.0 and r < 4.5) or (hole_radius and r < hole_radius and p[2] <= hole_top):
                    continue
                tube = (0.0, 0.0, hole_top), (-2.0 * n - 4, 0.0, hole_top)
                if hole_radius and hole_top < 2 * n and _segment_distance(p, *tube) < hole_radius:
                    continue
                if portal_radius and _segment_distance(p, (0.0, 0.0, 2.0), (2.0 * n + 4, 0.0, 2.0)) < portal_radius:
                    continue
                coords.append(p)
    coords += [(2.3, 0.0, -4.0), (-2.3, 0.0, -4.0), (0.0, 2.3, -4.0), (0.0, -2.3, -4.0)]
    coords.sort(key=lambda p: (p != (2.0 * n, 0.0, -2.0 * n), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


OPTIONS = dict(axis="z", origin=(0.0, 0.0, 0.0), samples=41, search_radius=8.0)


def test_a_wider_side_portal_ends_the_lumen_and_is_reported(monkeypatch):
    frame = branched_channel(portal_radius=3.0, hole_radius=2.1)
    with pytest.raises(ChannelResolutionError):
        pore_profile(frame, **OPTIONS)  # the upper end is capped
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    with monkeypatch.context() as patch:  # the trace alone follows the 1.13 A hole into the cap
        patch.setattr(adaptive_channel, "_lateral_branches", lambda trace, **kwargs: trace)
        traced = pore_profile(frame, lateral_exits=True, **OPTIONS)
    assert traced.points[-1].t > 6.0
    assert traced.metadata["exits"]["upper"]["bottleneck_raw_clearance_A"] < 1.2
    meta = profile.metadata
    upper = meta["exits"]["upper"]
    assert meta["channel_mouths"]["upper"]["reason"] == "lateral_branch"
    assert upper["type"] == "lateral" and upper["leg_count"] >= 1
    assert -1.0 < profile.points[-1].t < 3.0  # the lumen ends where the portal leaves (z = 2)
    branch = upper["lateral_branch"]
    assert branch["lateral_route_width_A"] > 2.0 > 1.2 > branch["continuation_width_A"]
    assert upper["bottleneck_raw_clearance_A"] > 2.0
    assert upper["legs"][0]["exit_point"][0] > 8.0  # through the +x portal
    assert meta["exits"]["lower"]["type"] == "axial"
    assert len(profile.points) == 41 and min(p.raw_clearance for p in profile.points) == pytest.approx(0.6, abs=0.02)


def test_a_narrower_branch_leaves_the_capped_end_in_place(monkeypatch):
    frame = branched_channel(portal_radius=2.5, hole_radius=3.0)
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    monkeypatch.setattr(adaptive_channel, "_lateral_branches", lambda trace, **kwargs: trace)
    unchanged = pore_profile(frame, lateral_exits=True, **OPTIONS)
    assert profile.to_dict() == unchanged.to_dict()
    assert profile.metadata["exits"]["upper"]["type"] == "lateral"
    assert profile.metadata["channel_mouths"]["upper"]["reason"] != "lateral_branch"
    assert profile.points[-1].t > 6.0


@pytest.mark.parametrize("portal_radius", [3.0, 2.5, 0.0])
def test_an_open_end_never_moves(portal_radius):
    frame = branched_channel(portal_radius=portal_radius, hole_radius=2.1, hole_top=99.0)
    plain = pore_profile(frame, **OPTIONS)
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    assert profile.points == plain.points
    assert profile.metadata["channel_mouths"] == plain.metadata["channel_mouths"]
    assert {e["type"] for e in profile.metadata["exits"].values()} == {"axial"}
    assert profile.metadata["path_type"] == "through_channel"


def test_a_capped_end_already_in_bulk_is_open(monkeypatch):
    from test_capped_channels import OPTIONS as CAPPED, capped_channel

    frame = capped_channel(portal=False)
    with pytest.raises(ChannelResolutionError):
        pore_profile(frame, **CAPPED)
    with pytest.raises(ChannelResolutionError):
        pore_profile(frame, lateral_exits=True, **CAPPED)
    # Software check of the plumbing: the closed cap's end sample is declared to lie in bulk.
    monkeypatch.setattr(adaptive_channel, "_end_in_bulk", lambda envelope, lumen: True)
    plain = pore_profile(frame, **CAPPED)
    assert plain.metadata["ends_in_bulk"]["ends"] == ["upper"]
    opted = pore_profile(frame, lateral_exits=True, **CAPPED)
    assert opted.metadata["exits"]["upper"]["type"] == "bulk"
    assert opted.metadata["path_type"] == "through_channel"
    assert opted.points == plain.points


def test_end_in_bulk_uses_the_bulk_envelope_depth():
    class Envelope:
        def __init__(self, depth):
            self.value = depth

        def depth(self, point):
            return self.value

    lumen = adaptive_channel.Lumen(0.0, (0.0, 0.0, 0.0), 1.0, 0.5, None)
    assert adaptive_channel._end_in_bulk(Envelope(0.0), lumen)
    assert not adaptive_channel._end_in_bulk(Envelope(0.3), lumen)

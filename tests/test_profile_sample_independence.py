"""Synthetic regression: a connected profile must not depend on the sample count.

A tube of stacked atom rings has a lumen that, in one thin slab between two
rings (about 0.1 Å thick), just reaches the edge of the 6 Å transverse search
square while protein still encloses it. Before 0.4.0, resampling planes that
fell in that slab broke the channel ("Channel continuity failed"), so the
same channel resolved at 21 samples and failed at 29. Synthetic geometry
only; not a real-system validation.
"""
import math

import pytest

from crevice.models import Atom
from crevice.sections import ChannelResolutionError, connected_profile

RINGS = (-8.0, -5.4, -2.8, 0.1, 2.7, 5.3, 7.9)  # one wider gap (-2.8 to 0.1): the grazing slab near -1.35


def _tube(radius=8.3, rings=RINGS, per_ring=48, gap=None):
    atoms = []
    for ring, z in enumerate(rings):
        for k in range(per_ring):
            if gap is not None and ring in gap[0] and k in gap[1]:
                continue
            a = 2 * math.pi * k / per_ring
            atoms.append(Atom(len(atoms) + 1, "C", "ALA", "A", ring + 1,
                              radius * math.cos(a), radius * math.sin(a), z, "C"))
    return tuple(atoms)


def _profile(atoms, samples):
    return connected_profile(atoms, axis=(0, 0, 1), origin=(0, 0, 0), samples=samples, search_radius=6.0,
                             refinement_steps=4, probe_radius=0.0, section_spacing=0.25, enclosure_radius=1.0)


def test_grazing_slab_is_bridged_only_where_the_lumen_is_enclosed():
    from crevice.adaptive_channel import _Tracer
    atoms = _tube()
    tracer = _Tracer(atoms, (0, 0, 0), (0, 0, 1), spacing=0.25, extent=6.0, ladder=(1.0,), refinement_steps=4,
                     shared=None)
    center = (0.0, 0.0, -1.35)
    found = tracer.enclosed_at(-1.35, center, 1.0)
    assert found is not None and found[1]                 # enclosed only in the wider square
    found = tracer.enclosed_at(-3.0, (0.0, 0.0, -3.0), 1.0)
    assert found is not None and not found[1]             # enclosed in the plane's own square
    # A real wall opening (eight atoms removed from the rings at z = -2.8, 0.1 and 2.7) stays open in
    # the wider square too.
    tracer = _Tracer(_tube(gap=((2, 3, 4), range(0, 8))), (0, 0, 0), (0, 0, 1), spacing=0.25, extent=6.0,
                     ladder=(1.0,), refinement_steps=4, shared=None)
    assert tracer.enclosed_at(0.1, (0.0, 0.0, 0.1), 1.0) is None


def test_profile_resolves_at_every_sample_count_and_agrees():
    atoms = _tube()
    plain, crossing = _profile(atoms, 21), _profile(atoms, 29)   # 29 samples place a plane in the slab
    # Planes in the grazing slab are enclosed in the wider square (bridged); none breaks the channel.
    for profile in (plain, crossing):
        assert all(abs(profile.points[i].t + 1.35) < 1.0 for i in profile.metadata["window_bridged_samples"])
    assert (plain.points[0].t, plain.points[-1].t) == pytest.approx((crossing.points[0].t, crossing.points[-1].t), abs=1e-9)
    assert crossing.min_radius == pytest.approx(plain.min_radius, abs=0.005)
    assert plain.metadata["continuous_bottleneck_radius_A"] == pytest.approx(
        crossing.metadata["continuous_bottleneck_radius_A"], abs=1e-9)


def test_a_real_lateral_opening_still_breaks_the_channel():
    open_wall = _tube(gap=((2, 3, 4), range(0, 8)))
    with pytest.raises(ChannelResolutionError):
        _profile(open_wall, 29)

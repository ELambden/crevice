"""Synthetic contracts for the per-plane enclosure probe (the automatic profile mode, 0.5).

Software/synthetic tests only; they say nothing about any real channel. The
fixtures are solid 2 Å carbon lattices with a carved lumen along z (as in
``test_auto_enclosure.py``):

* ``both``: a narrow constriction (four atoms 2.3 Å from the axis at z = +4,
  0.6 Å clearance) and, at z = -4, a side tube through the wall (clearance
  ~1.13 Å). One probe for the whole channel would have to be below 0.6 Å (to
  pass the constriction) and above 1.13 Å (to close the side tube), so no
  fixed probe resolves it; the per-plane probe can.
* ``waist(d)``: four atoms at distance ``d`` from the axis in the z = 0 plane,
  an exactly known constriction of clearance ``d - 1.7`` Å, used for the
  smallest resolvable radius against the section grid.
"""
import json
import math

import pytest

from crevice import Atom, StructureFrame, pore_profile
from crevice.adaptive_channel import ADAPTIVE_ENCLOSURE_LADDER, _smallest_true
from crevice.sections import ChannelResolutionError


def _segment_distance(p, a, b):
    ab = [b[i]-a[i] for i in range(3)]
    ap = [p[i]-a[i] for i in range(3)]
    t = max(0.0, min(1.0, sum(x*y for x, y in zip(ab, ap))/sum(x*x for x in ab)))
    return math.dist(p, [a[i]+t*ab[i] for i in range(3)])


def lattice(*, lumen=3.0, waist=None, waist_z=0.0, tube_z=None, sealed=False):
    coords = []
    for i in range(-4, 5):
        for j in range(-4, 5):
            for k in range(-4, 5):
                p = (2.0*i, 2.0*j, 2.0*k)
                if math.hypot(p[0], p[1]) < lumen and not (sealed and k == 0):
                    continue
                if tube_z is not None and _segment_distance(p, (0.0, 0.0, tube_z), (12.0, 0.0, tube_z)) < 2.5:
                    continue
                coords.append(p)
    if waist is not None:
        coords += [(waist, 0.0, waist_z), (-waist, 0.0, waist_z), (0.0, waist, waist_z), (0.0, -waist, waist_z)]
    coords.sort(key=lambda p: (p != (8.0, 0.0, -8.0), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


OPTIONS = dict(axis="z", origin=(0.0, 0.0, -1.0), samples=41, search_radius=7.0, section_spacing=0.25)


@pytest.fixture(scope="module")
def both():
    frame = lattice(waist=2.3, waist_z=4.0, tube_z=-4.0)
    return frame, pore_profile(frame, **OPTIONS)


def test_a_constriction_and_a_leaky_wall_together_need_a_per_plane_probe(both):
    frame, profile = both
    # No fixed probe covers both: 0.5 A leaks at the side tube, 1.2 A cannot pass the waist.
    assert profile.points[0].t < -6.5 and profile.points[-1].t > 6.0
    for fixed in (0.5, 1.2):
        try:
            other = pore_profile(frame, **{**OPTIONS, "enclosure_radius": fixed})
        except ChannelResolutionError:
            continue
        assert not (other.points[0].t < -6.5 and other.points[-1].t > 6.0)
    md = profile.metadata
    assert md["status"] == "resolved" and md["enclosure_mode"] == "adaptive"
    probes = list(zip((p.t for p in profile.points), md["enclosure_probes_A"]))
    # Planes at the constriction need probes below its 0.6 A clearance ...
    assert max(q for t, q in probes if abs(t-4.0) <= 1.0) < 0.6
    # ... and planes cutting the side tube (z = -4 +- 2.5) probes that close it (about 1 A).
    assert max(q for t, q in probes if abs(t+4.0) <= 1.5) >= 1.0
    # The measured radius does not depend on the probes: the exact waist clearance.
    assert md["continuous_bottleneck_radius_A"] == pytest.approx(0.6, abs=0.01)
    # Every measured sample is enclosed at its own probe and joins its neighbours.
    assert all(c > p for c, p in zip((pt.raw_clearance for pt in profile.points), md["enclosure_probes_A"]))
    assert all(s > min(a, b) for s, a, b in zip(md["segment_clearances_A"], md["enclosure_probes_A"],
                                              md["enclosure_probes_A"][1:]))


def test_per_plane_metadata_and_selection_record(both):
    _, profile = both
    md = profile.metadata
    assert len(md["enclosure_probes_A"]) == len(profile.points) == len(md["section_search_radii_A"])
    assert set(md["enclosure_probes_A"]) <= set(ADAPTIVE_ENCLOSURE_LADDER)
    assert md["enclosure_radius_A"] == max(md["enclosure_probes_A"])
    assert md["enclosure_probe_range_A"] == [min(md["enclosure_probes_A"]), max(md["enclosure_probes_A"])]
    selection = md["enclosure_probe_selection"]
    assert selection["mode"] == "adaptive" and selection["chosen_A"] is None
    assert selection["ladder_A"] == list(ADAPTIVE_ENCLOSURE_LADDER)


def test_the_result_is_deterministic(both):
    frame, profile = both
    again = pore_profile(frame, **OPTIONS)
    assert again.points == profile.points and again.metadata == profile.metadata


@pytest.mark.parametrize("spacing", [0.25, 0.5])
@pytest.mark.parametrize("samples", [41, 57])
def test_smallest_resolvable_constriction_against_the_section_grid(spacing, samples):
    """A waist resolves when its clearance exceeds the smallest probe (planes narrower than the grid refine it).

    The limit is the smallest ladder probe plus half a diagonal of the finest local grid
    (``section_spacing / MAX_GRID_REFINEMENT``), whatever the sample count, so whether a waist
    resolves does not depend on where the samples fall.
    """
    from crevice.adaptive_channel import MAX_GRID_REFINEMENT
    limit = ADAPTIVE_ENCLOSURE_LADDER[0] + spacing/MAX_GRID_REFINEMENT*math.sqrt(2.0)/2
    options = {**OPTIONS, "origin": (0.0, 0.0, -4.0), "section_spacing": spacing, "samples": samples}
    for clearance in (limit + 0.02, 0.5, 1.0):
        profile = pore_profile(lattice(waist=clearance + 1.7), **options)
        assert profile.metadata["continuous_bottleneck_radius_A"] == pytest.approx(clearance, abs=2e-3)
    with pytest.raises(ChannelResolutionError):
        pore_profile(lattice(waist=ADAPTIVE_ENCLOSURE_LADDER[0] - 0.02 + 1.7), **options)


def test_a_sealed_channel_is_refused_with_an_adaptive_record():
    with pytest.raises(ChannelResolutionError) as info:
        pore_profile(lattice(sealed=True), **{**OPTIONS, "origin": (0.0, 0.0, -4.0)})
    assert info.value.probe_selection["mode"] == "adaptive" and info.value.probe_selection["chosen_A"] is None


def test_profile_ends_inside_the_van_der_waals_envelope():
    from test_channel_coordinate import rectangular_channel
    profile = pore_profile(rectangular_channel(), samples=41, search_radius=6)
    assert -11.7 - 1e-9 <= profile.points[0].t and profile.points[-1].t <= 11.7 + 1e-9


def test_monotone_search_returns_the_smallest_true_for_any_hint():
    values = [0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4, 2.0, 3.0]
    for n in range(len(values) + 1):
        v = values[:n]
        for cut in range(n + 1):
            expected = (v[cut], True) if cut < n else None
            for hint in (None, 0.1, 0.3, 0.45, 0.8, 5.0):
                assert _smallest_true(v, lambda x: cut < n and x >= v[cut], hint) == expected


def test_channel_cast_uses_each_plane_probe_and_fills_every_plane(both):
    import warnings
    from crevice.sections import channel_cast
    frame, profile = both
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        cast = channel_cast(frame, profile, spacing=0.5, min_radius=0.0, max_grid_points=2_000_000,
                            max_components=1, max_export_points=5000, min_component_points=1, focus_points=None,
                            focus_radius=0.0, min_component_volume=0.0, max_component_volume=None,
                            include_hydrogen=False, include_hetero=True)
    md = cast.metadata
    assert md["status"] == "resolved" and md["unresolved_section_coordinates_A"] == []
    assert len(md["plane_enclosure_radii_A"]) == md["grid_dimensions"][2]
    assert min(md["plane_enclosure_radii_A"]) < 0.6 < 1.0 <= max(md["plane_enclosure_radii_A"])


def test_trajectory_finds_the_per_plane_probe_in_every_frame(both):
    from crevice.trajectory import analyze_trajectory
    frame, profile = both
    result = analyze_trajectory([frame, frame], profile_kwargs=dict(OPTIONS), align=False)
    probe = result.metadata["enclosure_probe"]
    assert probe["mode"] == "adaptive" and probe["chosen_A"] is None
    assert probe["frames"] == {"resolved_frames": 2, "unresolved_frames": 0}
    assert [f.metadata["enclosure_probe"]["source"] for f in result.frames] == ["per_plane", "per_plane"]
    assert all(f.profile.metadata["enclosure_probes_A"] == profile.metadata["enclosure_probes_A"]
               for f in result.frames)


def _write_pdb(frame, path):
    path.write_text("\n".join(
        f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
        for a in frame.atoms) + "\nEND\n", encoding="utf-8", newline="\n")
    return path


def test_cli_default_is_per_plane_and_says_what_it_found(tmp_path, capsys, monkeypatch):
    from crevice.cli import main
    from crevice import progress
    path = _write_pdb(lattice(waist=2.3, waist_z=4.0, tube_z=-4.0), tmp_path / "both.pdb")
    common = ["--axis", "z", "--origin", "0,0,-1", "--search-radius", "7", "--samples", "41",
              "--section-spacing", "0.25"]
    monkeypatch.delenv(progress.PROGRESS_ENV, raising=False)
    assert main(["profile", str(path), "-o", str(tmp_path / "p.json"), *common]) == 0
    captured = capsys.readouterr()
    assert "enclosure_radius=per-plane" in captured.out
    lines = captured.err.splitlines()
    assert any(line.startswith("crevice: initialising the enclosure probe for both.pdb (per-plane probe") for line in lines)
    found = next(line for line in lines if line.startswith("crevice: pore radius found in "))
    assert "minimum" in found and "per-plane enclosure probes" in found
    assert lines[-1].startswith("crevice: profile: finished in ")
    data = json.loads((tmp_path / "p.json").read_text(encoding="utf-8"))
    assert data["metadata"]["enclosure_probe_selection"]["mode"] == "adaptive"
    # --quiet silences every status line; stdout is unchanged.
    assert main(["profile", str(path), "-o", str(tmp_path / "q.json"), "--quiet", *common]) == 0
    quiet = capsys.readouterr()
    assert quiet.err == "" and quiet.out == captured.out.replace("p.json", "q.json")
    # CREVICE_PROGRESS=0 silences them too.
    monkeypatch.setenv(progress.PROGRESS_ENV, "0")
    assert main(["profile", str(path), "-o", str(tmp_path / "r.json"), *common]) == 0
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("samples", [41, 97, 203])
def test_channel_extent_does_not_depend_on_the_sample_count(both, samples):
    """The mouths come from the discovery trace; resampling only adds points between them."""
    frame, reference = both
    profile = pore_profile(frame, **{**OPTIONS, "samples": samples})
    assert (profile.points[0].t, profile.points[-1].t) == (reference.points[0].t, reference.points[-1].t)
    assert profile.metadata["channel_mouths"] == reference.metadata["channel_mouths"]
    assert profile.metadata["continuous_bottleneck_radius_A"] == pytest.approx(0.6, abs=0.01)
    flagged = set(profile.metadata["unenclosed_samples"]) | set(profile.metadata["unjoined_samples"])
    assert flagged <= set(range(len(profile.points)))

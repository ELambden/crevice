"""Synthetic contracts for the automatic enclosure probe (``enclosure_radius=None``).

Software/synthetic tests only; they say nothing about any real channel. Each
fixture is a solid 2 Å carbon lattice (x, y, z in [-8, 8]; the lattice itself
is sealed for any probe above ~0.04 Å) with a carved lumen along z:

* narrow: lumen radius 3 Å with four extra atoms at 2.3 Å from the axis in the
  z = 0 plane, a constriction of 0.6 Å clearance, narrower than the legacy
  0.8 Å probe;
* leaky: a wide lumen (5 Å, on-axis clearance ~3.96 Å) with a side tube at
  z = 0 through the wall to the +x face (clearance ~1.13 Å), wider than 0.8 Å,
  so at 0.8 Å the lumen is not laterally enclosed there;
* closed: the narrow lumen sealed by a full lattice layer, which no probe
  resolves.
"""
import json
import math

import pytest

from crevice import Atom, StructureFrame, pore_profile
from crevice.sections import AUTO_ENCLOSURE_LADDER, ChannelResolutionError


def _segment_distance(p, a, b):
    ab = [b[i]-a[i] for i in range(3)]
    ap = [p[i]-a[i] for i in range(3)]
    t = max(0.0, min(1.0, sum(x*y for x, y in zip(ab, ap))/sum(x*x for x in ab)))
    return math.dist(p, [a[i]+t*ab[i] for i in range(3)])


def lattice_channel(*, lumen=3.0, constriction=None, side_tube=False, sealed=False):
    coords = []
    for i in range(-4, 5):
        for j in range(-4, 5):
            for k in range(-4, 5):
                p = (2.0*i, 2.0*j, 2.0*k)
                if math.hypot(p[0], p[1]) < lumen and not (sealed and k == 0):
                    continue
                if side_tube and _segment_distance(p, (0.0, 0.0, 0.0), (12.0, 0.0, 0.0)) < 2.5:
                    continue
                coords.append(p)
    if constriction is not None:
        coords += [(constriction, 0.0, 0.0), (-constriction, 0.0, 0.0),
                   (0.0, constriction, 0.0), (0.0, -constriction, 0.0)]
    # A fixed first atom keeps the transverse basis deterministic.
    coords.sort(key=lambda p: (p != (8.0, 0.0, -8.0), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


OPTIONS = dict(axis="z", origin=(0.0, 0.0, -4.0), samples=25, search_radius=7.0, section_spacing=0.5)


@pytest.fixture(scope="module")
def narrow():
    frame = lattice_channel(constriction=2.3)
    return frame, pore_profile(frame, **OPTIONS)


@pytest.fixture(scope="module")
def leaky():
    frame = lattice_channel(lumen=5.0, side_tube=True)
    return frame, pore_profile(frame, **OPTIONS)


def test_narrow_channel_resolves_below_the_legacy_probe(narrow):
    frame, profile = narrow
    with pytest.raises(ChannelResolutionError):
        pore_profile(frame, enclosure_radius=0.8, **OPTIONS)
    selection = profile.metadata["enclosure_probe_selection"]
    assert selection["mode"] == "auto" and selection["stability"] == "plateau"
    chosen = selection["chosen_A"]
    # The probe must stay below the constriction (0.6 A clearance).
    assert chosen < profile.metadata["continuous_bottleneck_radius_A"] < 0.61
    assert profile.metadata["continuous_bottleneck_radius_A"] == pytest.approx(0.6, abs=0.01)
    assert 0.6 <= profile.min_radius < 0.7   # samples need not land exactly on z = 0
    assert profile.metadata["enclosure_radius_A"] == chosen
    # Smallest probe of the winning group, which has at least two probes.
    group = next(g for g in selection["groups"] if g["group"] == selection["chosen_group"])
    assert group["support"] >= 2 and chosen == min(group["probes_A"])
    assert [row["probe_A"] for row in selection["candidates"] if row["role"] == "ladder"] == list(AUTO_ENCLOSURE_LADDER)


def test_auto_profile_equals_the_explicit_run_at_the_chosen_probe(narrow):
    frame, profile = narrow
    chosen = profile.metadata["enclosure_probe_selection"]["chosen_A"]
    explicit = pore_profile(frame, enclosure_radius=chosen, **OPTIONS)
    assert "enclosure_probe_selection" not in explicit.metadata
    assert explicit.points == profile.points
    assert {k: v for k, v in profile.metadata.items() if k != "enclosure_probe_selection"} == explicit.metadata


def test_wide_channel_with_a_leaky_wall_closes_the_gap_but_not_the_lumen(leaky):
    frame, profile = leaky
    # At the legacy probe the side tube (clearance ~1.13 A) breaks enclosure.
    with pytest.raises(ChannelResolutionError):
        pore_profile(frame, enclosure_radius=0.8, **OPTIONS)
    selection = profile.metadata["enclosure_probe_selection"]
    chosen = selection["chosen_A"]
    assert 2*math.sqrt(2)-1.7 < chosen < profile.metadata["continuous_bottleneck_radius_A"]
    # Radii are atom-surface distances, independent of the probe: the lumen's
    # on-axis clearance is the (4, 4) lattice site at 5.657 A minus 1.7 A.
    assert profile.min_radius == pytest.approx(math.hypot(4, 4)-1.7, abs=0.02)
    smaller = [r for r in selection["candidates"] if r["probe_A"] < 1.1]
    assert smaller and all(r["status"] == "unresolved" for r in smaller)


def test_wide_channel_uses_the_smallest_probe_of_its_plateau(leaky):
    _, profile = leaky
    selection = profile.metadata["enclosure_probe_selection"]
    resolved = [r["probe_A"] for r in selection["candidates"] if r["status"] == "resolved"]
    assert selection["chosen_A"] == min(resolved)
    assert len(selection["groups"]) == 1 and selection["groups"][0]["support"] == len(resolved) >= 2


def test_closed_channel_is_refused_with_the_record_attached():
    frame = lattice_channel(constriction=2.3, sealed=True)
    with pytest.raises(ChannelResolutionError, match="No enclosure probe resolved a stable channel") as info:
        pore_profile(frame, **OPTIONS)
    record = info.value.probe_selection
    assert record["chosen_A"] is None and record["groups"] == []
    assert {row["status"] for row in record["candidates"]} <= {"unresolved", "not_tried"}


def test_explicit_origin_bounds_the_candidates(narrow):
    _, profile = narrow
    rows = profile.metadata["enclosure_probe_selection"]["candidates"]
    skipped = [r for r in rows if r["status"] == "not_tried"]
    # The seed at z = -4 has 1.3 A clearance (lumen radius 3 A): larger probes cannot seed.
    assert skipped and min(r["probe_A"] for r in skipped) >= 1.3-1e-9
    assert all("origin" in r["reason"] for r in skipped)


def test_invalid_enclosure_radius_is_rejected():
    frame = lattice_channel()
    with pytest.raises(ValueError, match="enclosure_radius"):
        pore_profile(frame, enclosure_radius=-1.0, **OPTIONS)


def test_cli_auto_is_the_default_and_is_reported(tmp_path, capsys):
    from crevice.cli import main
    frame = lattice_channel(constriction=2.3)
    path = tmp_path / "narrow.pdb"
    path.write_text("\n".join(
        f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
        for a in frame.atoms) + "\nEND\n")
    common = ["--axis", "z", "--origin", "0,0,-4", "--search-radius", "7", "--samples", "25",
              "--section-spacing", "0.5"]
    assert main(["profile", str(path), "-o", str(tmp_path / "auto.json"), *common]) == 0
    out = capsys.readouterr().out
    assert "enclosure_radius=" in out and "(auto:" in out
    data = json.loads((tmp_path / "auto.json").read_text())
    assert data["metadata"]["enclosure_probe_selection"]["mode"] == "auto"
    assert main(["profile", str(path), "-o", str(tmp_path / "auto2.json"), "--enclosure-radius", "auto", *common]) == 0
    assert json.loads((tmp_path / "auto2.json").read_text()) == data
    with pytest.raises(SystemExit):
        main(["profile", str(path), "--enclosure-radius", "wide", *common])
    capsys.readouterr()
    # An explicit probe below the constriction resolves without any selection record.
    assert main(["profile", str(path), "-o", str(tmp_path / "fixed.json"), "--enclosure-radius", "0.5", *common]) == 0
    assert "(auto:" not in capsys.readouterr().out
    assert "enclosure_probe_selection" not in json.loads((tmp_path / "fixed.json").read_text())["metadata"]


def test_publish_records_the_choice_or_the_refusal(tmp_path):
    from crevice.cli import main
    for name, frame, resolved in (("narrow", lattice_channel(constriction=2.3), True),
                                  ("closed", lattice_channel(constriction=2.3, sealed=True), False)):
        path = tmp_path / f"{name}.pdb"
        path.write_text("\n".join(
            f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
            for a in frame.atoms) + "\nEND\n")
        out = tmp_path / name
        with pytest.warns(RuntimeWarning) if not resolved else _no_warning_check():
            assert main(["publish", str(path), "--out-dir", str(out), "--axis", "z", "--origin", "0,0,-4",
                         "--search-radius", "7", "--samples", "25", "--section-spacing", "0.5",
                         "--skip-cavities", "--skip-tunnels", "--skip-network", "--dpi", "40"]) == 0
        manifest = json.loads((out / f"{name}_manifest.json").read_text())
        probe = manifest["profile_status"]["enclosure_probe"]
        assert probe["mode"] == "auto"
        if resolved:
            assert probe["chosen_A"] is not None and probe["candidates_tried"]
        else:
            assert probe["chosen_A"] is None and probe["fallback_cast_probe_A"] == 0.8
            status = json.loads((out / f"{name}_profile_status.json").read_text())
            assert status["enclosure_probe_selection"]["chosen_A"] is None


class _no_warning_check:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


# --- Rollout to static-suite, benchmark and trajectory (1 Oct 2026) -----------
# A short ladder keeps these tests fast; the selection rule is unchanged.
SHORT_LADDER = (0.50, 0.55, 0.60, 0.80, 1.00)


@pytest.fixture
def short_ladder(monkeypatch):
    import crevice.sections as sections
    monkeypatch.setattr(sections, "AUTO_ENCLOSURE_LADDER", SHORT_LADDER)
    calls = []
    original = sections._auto_enclosure_profile

    def counted(*args, **kwargs):
        calls.append(kwargs.get("origin"))
        return original(*args, **kwargs)
    monkeypatch.setattr(sections, "_auto_enclosure_profile", counted)
    return calls


def _perturbed(frame, dz=0.05):
    """The same lattice with the constriction atoms moved slightly along z."""
    from dataclasses import replace
    return replace(frame, atoms=tuple(replace(a, z=a.z+dz) if abs(a.z) < 1e-9 and math.hypot(a.x, a.y) < 3
                                      else a for a in frame.atoms))


def test_trajectory_chooses_the_probe_once_on_the_reference_and_pins_it(short_ladder):
    from crevice.trajectory import analyze_trajectory
    frame = lattice_channel(constriction=2.3)
    reference = pore_profile(frame, **OPTIONS)
    short_ladder.clear()
    chosen = reference.metadata["enclosure_probe_selection"]["chosen_A"]
    result = analyze_trajectory([frame, _perturbed(frame), frame], profile_kwargs=dict(OPTIONS))
    # One automatic choice (the reference frame); every frame is an explicit run.
    assert len(short_ladder) == 1
    probe = result.metadata["enclosure_probe"]
    assert probe["mode"] == "auto" and probe["chosen_A"] == chosen
    assert probe["selection"]["chosen_A"] == chosen and probe["reason"] == probe["selection"]["reason"]
    assert probe["applied"] == "fixed for every frame"
    assert all(f.profile is not None for f in result.frames)
    assert {f.profile.metadata["enclosure_radius_A"] for f in result.frames} == {chosen}
    assert all("enclosure_probe_selection" not in f.profile.metadata for f in result.frames)


def test_trajectory_explicit_probe_overrides_the_automatic_choice(short_ladder):
    from crevice.trajectory import analyze_trajectory
    frame = lattice_channel(constriction=2.3)
    result = analyze_trajectory([frame, frame], profile_kwargs={**OPTIONS, "enclosure_radius": 0.5})
    assert short_ladder == []
    probe = dict(result.metadata["enclosure_probe"])
    fallback = probe.pop("fallback")
    assert probe == {"mode": "explicit", "chosen_A": 0.5, "reason": "enclosure_radius given explicitly",
                     "applied": "fixed for every frame"}
    assert fallback["enabled"] is False and fallback["pinned_frames"] == 2 and fallback["fallback_frames"] == 0
    assert {f.profile.metadata["enclosure_radius_A"] for f in result.frames} == {0.5}


def test_trajectory_records_a_failed_automatic_choice(short_ladder):
    from crevice.trajectory import analyze_trajectory
    frame = lattice_channel(constriction=2.3, sealed=True)
    result = analyze_trajectory([frame, frame], profile_kwargs=dict(OPTIONS))
    probe = result.metadata["enclosure_probe"]
    assert probe["mode"] == "auto" and probe["chosen_A"] is None and probe["selection"]["chosen_A"] is None
    assert len(short_ladder) == 1
    assert {f.metadata["profile_status"] for f in result.frames} == {"not_attempted_reference_unresolved"}
    with pytest.raises(ChannelResolutionError):
        analyze_trajectory([frame], profile_kwargs=dict(OPTIONS), on_profile_error="raise")


def test_trajectory_fixed_axis_scan_records_no_probe():
    from crevice.trajectory import analyze_trajectory
    frame = lattice_channel()
    result = analyze_trajectory([frame], profile_kwargs={"axis": "z", "search_radius": 0, "samples": 5})
    assert result.metadata["enclosure_probe"]["mode"] == "not_used"
    assert analyze_trajectory([frame], analyses=("cavities",),
                              cavity_kwargs={"max_grid_points": 20_000}).metadata["enclosure_probe"] is None


def _write_pdb(frame, path):
    path.write_text("\n".join(
        f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
        for a in frame.atoms) + "\nEND\n")
    return path


def test_trajectory_cli_defaults_to_auto_and_records_the_choice(tmp_path, capsys, short_ladder):
    import csv
    from crevice.cli import main
    frame = lattice_channel(constriction=2.3)
    paths = [str(_write_pdb(f, tmp_path / f"f{i}.pdb")) for i, f in enumerate((frame, _perturbed(frame), frame))]
    common = ["--axis", "z", "--origin", "0,0,-4", "--search-radius", "7", "--samples", "25",
              "--section-spacing", "0.5", "--no-align", "--confidence", "0", "--skip-hydration"]
    assert main(["trajectory", *paths, "-o", str(tmp_path / "auto.json"), *common]) == 0
    out = capsys.readouterr().out
    assert "auto, chosen once on the reference frame and fixed for every frame" in out
    assert len(short_ladder) == 1
    data = json.loads((tmp_path / "auto.json").read_text())
    probe = data["metadata"]["enclosure_probe"]
    assert probe["mode"] == "auto" and probe["chosen_A"] in SHORT_LADDER
    with (tmp_path / "auto_frames.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3 and {float(r["enclosure_radius_A"]) for r in rows} == {probe["chosen_A"]}
    short_ladder.clear()
    assert main(["trajectory", *paths, "-o", str(tmp_path / "fixed.json"), "--enclosure-radius", "0.5", *common]) == 0
    assert "chosen once" not in capsys.readouterr().out and short_ladder == []
    fixed = json.loads((tmp_path / "fixed.json").read_text())["metadata"]["enclosure_probe"]
    assert fixed["mode"] == "explicit" and fixed["chosen_A"] == 0.5


def test_static_suite_and_benchmark_default_to_auto(tmp_path, short_ladder):
    import csv
    from crevice.cli import main
    from crevice.presentation import write_viewer_structure
    from crevice.static_suite import run_static_benchmark_suite
    from test_channel_coordinate import rectangular_channel
    cache = tmp_path / "cache"
    cache.mkdir()
    write_viewer_structure(rectangular_channel(half_length=6), cache / "1GRM.pdb")
    rows = {}
    for label, enclosure in (("auto", None), ("fixed", 0.8)):
        result = run_static_benchmark_suite(
            cache_dir=cache, output_dir=tmp_path / label, pdb_ids=("1GRM",), extension="pdb",
            include_network=False, include_tunnels=False, include_cavities=False, cast_spacing=1,
            profile_samples=11, search_radius=8, probe_radii=(0.0, 3.0), dpi=40, enclosure_radius=enclosure)
        system = result["systems"][0]
        # The 3 A measurement probe exceeds the ~2.3 A clearance: an unresolved
        # sensitivity row, not a failed system.
        assert system["status"] == "ok" and system["profile_sensitivity_unresolved_probe_radii"] == [3.0]
        with open(tmp_path / label / "1GRM" / "1GRM_profile_sensitivity.csv") as handle:
            sensitivity = list(csv.DictReader(handle))
        assert [r["status"] for r in sensitivity] == ["resolved", "unresolved"]
        assert {float(r["enclosure_radius_A"]) for r in sensitivity} == {system["enclosure_radius"]}
        assert result["parameters"]["enclosure_radius"] == ("auto" if enclosure is None else 0.8)
        with open(result["summary_csv"]) as handle:
            rows[label] = next(csv.DictReader(handle))
    assert rows["auto"]["enclosure_probe_mode"] == "auto" and float(rows["auto"]["enclosure_radius_A"]) in SHORT_LADDER
    assert rows["fixed"]["enclosure_probe_mode"] == "explicit" and float(rows["fixed"]["enclosure_radius_A"]) == 0.8
    assert len(short_ladder) == 1   # sensitivity rows reuse the system's probe
    for label, extra in (("auto", []), ("fixed", ["--enclosure-radius", "0.8"])):
        output = tmp_path / f"benchmark_{label}.json"
        assert main(["benchmark", "1GRM", "--cache-dir", str(cache), "--format", "pdb", "--samples", "11",
                     "--search-radius", "8", "--output", str(output), *extra]) == 0
        row = json.loads(output.read_text())["benchmarks"][0]
        assert row["status"] == "ok" and row["enclosure_probe_mode"] == ("auto" if not extra else "explicit")
        assert row["enclosure_radius"] == (float(rows["auto"]["enclosure_radius_A"]) if not extra else 0.8)


# --- Per-frame fallback (2 Oct 2026) ------------------------------------------
# Reference: the narrow lattice (automatic probe ~0.5 Å). Second frame: the
# leaky lattice, whose side tube (clearance ~1.13 Å) is open to that small
# probe, so the pinned probe cannot resolve it but a larger per-frame probe can.
# align=False: the two lattices have different atoms.

def test_pinned_auto_probe_falls_back_per_frame_and_records_it(narrow, leaky, tmp_path):
    from crevice.io import write_trajectory_csv
    from crevice.trajectory import analyze_trajectory
    import csv
    pinned = narrow[1].metadata["enclosure_probe_selection"]["chosen_A"]
    rescued = leaky[1].metadata["enclosure_probe_selection"]["chosen_A"]
    assert rescued > pinned
    result = analyze_trajectory([narrow[0], leaky[0], narrow[0]], profile_kwargs=dict(OPTIONS), align=False)
    sources = [f.metadata["enclosure_probe"]["source"] for f in result.frames]
    assert sources == ["pinned", "fallback", "pinned"]
    middle = result.frames[1]
    assert middle.profile is not None and middle.metadata["profile_status"] == "resolved"
    assert middle.metadata["enclosure_probe"]["radius_A"] == middle.profile.metadata["enclosure_radius_A"] > pinned
    assert "pinned_failure" in middle.metadata["enclosure_probe"]
    fallback = result.metadata["enclosure_probe"]["fallback"]
    assert fallback["enabled"] and fallback["fallback_frames"] == 1 and fallback["pinned_frames"] == 2
    assert fallback["fallback_frame_indices"] == [middle.frame_index]
    write_trajectory_csv(result, tmp_path / "f.csv", tmp_path / "p.csv")
    rows = list(csv.DictReader(open(tmp_path / "f.csv")))
    assert [r["enclosure_probe_source"] for r in rows] == sources
    assert float(rows[1]["enclosure_radius_A"]) == middle.profile.metadata["enclosure_radius_A"]


def test_explicit_probe_stays_strict_unless_fallback_is_requested(narrow, leaky):
    from crevice.trajectory import analyze_trajectory
    pinned = narrow[1].metadata["enclosure_probe_selection"]["chosen_A"]
    frames = [narrow[0], leaky[0]]
    strict = analyze_trajectory(frames, profile_kwargs={**OPTIONS, "enclosure_radius": pinned}, align=False)
    assert [f.metadata["enclosure_probe"]["source"] for f in strict.frames] == ["pinned", "none"]
    assert strict.frames[1].profile is None and strict.frames[1].metadata["enclosure_probe"]["fallback"] == "disabled"
    assert strict.metadata["enclosure_probe"]["fallback"]["enabled"] is False
    assert strict.metadata["enclosure_probe"]["fallback"]["unresolved_frames"] == 1
    with pytest.raises(ChannelResolutionError):
        analyze_trajectory(frames, profile_kwargs={**OPTIONS, "enclosure_radius": pinned}, align=False,
                           on_profile_error="raise")
    opted = analyze_trajectory(frames, profile_kwargs={**OPTIONS, "enclosure_radius": pinned}, align=False,
                               probe_fallback=True)
    assert [f.metadata["enclosure_probe"]["source"] for f in opted.frames] == ["pinned", "fallback"]
    off = analyze_trajectory(frames, profile_kwargs=dict(OPTIONS), align=False, probe_fallback=False)
    assert [f.metadata["enclosure_probe"]["source"] for f in off.frames] == ["pinned", "none"]


def test_fallback_frames_pool_into_the_distribution_and_are_counted(narrow, leaky):
    from crevice.ensemble import profile_distribution
    from crevice.trajectory import analyze_trajectory
    result = analyze_trajectory([narrow[0], leaky[0]], profile_kwargs=dict(OPTIONS), align=False)
    distribution = profile_distribution(result, assume_aligned=True)
    assert distribution["frame_count"] == 2 and distribution["enclosure_probe_fallback_frame_count"] == 1
    assert len(distribution["enclosure_radii_A"]) == 2

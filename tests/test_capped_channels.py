"""Synthetic contracts for capped channel ends with lateral exits (``lateral_exits``).

Software/synthetic tests only; they say nothing about any real channel.
"""
import json
import math

import pytest

from crevice import Atom, StructureFrame, pore_profile
from crevice.geometry import dot, sub
from crevice.sections import ChannelResolutionError
from crevice.spatial import SpatialIndex

from test_channel_coordinate import rectangular_channel


def _segment_distance(p, a, b):
    ab = [b[i]-a[i] for i in range(3)]
    ap = [p[i]-a[i] for i in range(3)]
    t = max(0.0, min(1.0, sum(x*y for x, y in zip(ab, ap))/sum(x*x for x in ab)))
    return math.dist(p, [a[i]+t*ab[i] for i in range(3)])


def capped_channel(*, portal=True, lumen_top=2.0, second_portal=False):
    """Solid 2 A carbon lattice (x, y, z in [-8, 8]) with a lumen and a cap.

    The lumen (lattice sites within 3 A of the z axis removed) runs from the
    open bottom face up to z = 2; the lattice layers above it are a cap, so
    the straight upward exit is blocked. With ``portal`` a side tube (sites
    within 2.5 A of the segment (0, 0, 2)-(12, 0, 2) removed) leads from the
    top of the lumen to the +x face: a lateral exit with ~1.13 A clearance.
    The first atom fixes the transverse basis to (+x, +y). ``lumen_top=6``
    extends the lumen into a chamber above the portal, closed by the top layer:
    enclosed sections beyond the mouth, as under the TWIK-1 cap.
    ``second_portal`` adds the mirror-image tube to the -x face, so the cap
    has two symmetric side openings that tie exactly (like the two
    C2-related portals of TWIK-1).
    """
    coords = []
    for i in range(-4, 5):
        for j in range(-4, 5):
            for k in range(-4, 5):
                p = (2.0*i, 2.0*j, 2.0*k)
                if math.hypot(p[0], p[1]) < 3.0 and p[2] <= lumen_top:
                    continue
                if portal and _segment_distance(p, (0.0, 0.0, 2.0), (12.0, 0.0, 2.0)) < 2.5:
                    continue
                if second_portal and _segment_distance(p, (0.0, 0.0, 2.0), (-12.0, 0.0, 2.0)) < 2.5:
                    continue
                coords.append(p)
    coords.sort(key=lambda p: (p != (8.0, 0.0, -8.0), p))
    return StructureFrame(tuple(Atom(i, "C", "ALA", "A", i, *p, "C") for i, p in enumerate(coords, 1)))


OPTIONS = dict(axis="z", origin=(0.0, 0.0, 0.0), samples=25, search_radius=7.0, section_spacing=0.5,
               enclosure_radius=0.8)


@pytest.fixture(scope="module")
def capped():
    frame = capped_channel()
    return frame, pore_profile(frame, lateral_exits=True, **OPTIONS)


def test_capped_end_is_refused_without_the_option():
    with pytest.raises(ChannelResolutionError, match="No connected channel"):
        pore_profile(capped_channel(), **OPTIONS)


def test_capped_end_leaves_through_the_lateral_portal(capped):
    frame, profile = capped
    md = profile.metadata
    assert md["path_type"] == "capped_channel_lateral_exit"
    assert md["exits"]["lower"]["type"] == "axial"
    end = md["exits"]["upper"]
    assert end["type"] == "lateral" and end["straight_exit_clearance_A"] <= 0.8
    # One portal: exactly one distinct leg, and the end reports its bottleneck.
    assert end["leg_count"] == 1 and len(end["legs"]) == 1
    upper = end["legs"][0]
    assert upper["leg"] == 1 and end["bottleneck_radius_A"] == upper["bottleneck_radius_A"]
    # The leg starts at the capped end of the axial profile, never re-enters
    # the channel below that plane, and ends outside the +x face in bulk.
    start = profile.points[-1]
    assert upper["points"][0] == pytest.approx(list(start.position))
    assert min(dot(sub(p, start.position), profile.axis_direction) for p in upper["points"]) >= -1e-9
    assert upper["exit_point"][0] > 8.0
    spatial = SpatialIndex(frame.atoms)
    vertices = upper["vertices"]
    assert all(spatial.segment_clearance(a, b) > 0.8 for a, b in zip(vertices, vertices[1:]))
    # Bottleneck of the leg is the portal (lattice sites at 2.83 A from its axis).
    assert upper["bottleneck_radius_A"] == pytest.approx(2*math.sqrt(2)-1.7, abs=0.06)
    assert md["path_bottleneck_radius_A"] == pytest.approx(
        min(md["continuous_bottleneck_radius_A"], upper["bottleneck_radius_A"]))
    # Axial samples, and hence min_radius, stay the planar section profile.
    assert profile.min_radius > upper["bottleneck_radius_A"]
    assert md["exit_policy"]["capped_ends"] == ["upper"]


def test_capped_end_without_a_connected_exit_is_refused():
    with pytest.raises(ChannelResolutionError, match="upper end is capped"):
        pore_profile(capped_channel(portal=False), lateral_exits=True, **OPTIONS)


def test_through_channel_is_unchanged_by_the_option():
    frame = rectangular_channel()
    plain = pore_profile(frame, samples=41, search_radius=8)
    opted = pore_profile(frame, samples=41, search_radius=8, lateral_exits=True)
    assert [p.position for p in opted.points] == [p.position for p in plain.points]
    assert [p.radius for p in opted.points] == [p.radius for p in plain.points]
    assert opted.metadata["path_type"] == "through_channel"
    assert {e["type"] for e in opted.metadata["exits"].values()} == {"axial"}
    assert "exits" not in plain.metadata and "path_type" not in plain.metadata
    assert {k: v for k, v in opted.metadata.items() if k in plain.metadata} == plain.metadata


def test_lateral_exit_controls_are_validated():
    frame = capped_channel()
    with pytest.raises(ValueError, match="search_radius"):
        pore_profile(frame, axis="z", search_radius=0, lateral_exits=True)
    with pytest.raises(ValueError, match="exit_bulk_radius"):
        pore_profile(frame, lateral_exits=True, exit_bulk_radius=0.5, **OPTIONS)


def _publish_capped(tmp_path, *extra):
    from crevice.cli import main
    frame = capped_channel()
    path = tmp_path / "capped.pdb"
    path.write_text("\n".join(
        f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
        for a in frame.atoms) + "\nEND\n")
    out = tmp_path / "out"
    assert main(["publish", str(path), "--out-dir", str(out), "--axis", "z", "--origin", "0,0,0",
                 "--search-radius", "7", "--samples", "25", "--section-spacing", "0.5",
                 "--enclosure-radius", "0.8", "--lateral-exits", "--skip-cavities", "--skip-tunnels", "--skip-network",
                 "--skip-hydration", "--dpi", "60", *extra]) == 0
    return out


def test_publish_reports_path_type_and_draws_the_exit_leg_as_a_cast(tmp_path):
    import csv
    out = _publish_capped(tmp_path)
    manifest = json.loads((out / "capped_manifest.json").read_text())
    assert manifest["profile_status"]["path_type"] == "capped_channel_lateral_exit"
    assert manifest["profile_status"]["exit_types"] == {"lower": "axial", "upper": "lateral"}
    profile = json.loads((out / "capped_profile.json").read_text())
    assert profile["metadata"]["exits"]["upper"]["type"] == "lateral"
    assert manifest["profile_status"]["exit_leg_counts"] == {"lower": 0, "upper": 1}
    # The leg is a separate cast component; the axial volume excludes it.
    cast = json.loads((out / "capped_void_cast.json").read_text())
    legs = cast["metadata"]["lateral_exit_casts"]
    assert [(leg["end"], leg["leg"], leg["kind"]) for leg in legs] == [("upper", 1, "lateral_exit")]
    assert legs[0]["volume_A3"] > 0
    assert cast["total_volume"] == sum(c["volume"] for c in cast["components"])
    assert all(c["kind"] != "lateral_exit" for c in cast["components"])
    assert manifest["profile_status"]["lateral_exit_cast_volumes_A3"] == {"upper_1": legs[0]["volume_A3"]}
    rows = list(csv.DictReader((out / "capped_void_cast.csv").open()))
    assert [r["label"] for r in rows if r["kind"] == "lateral_exit"] == ["lateral_exit_upper_1"]
    assert (out / "capped_exit_casts.dx").exists()
    metadata = json.loads((out / "capped_volume_metadata.json").read_text())
    assert metadata["measured_volume_A3"] == cast["total_volume"]
    assert metadata["display_points"] > metadata["measurement_points"]
    # Lines and rings are opt-in; their coordinates stay in the scene JSON.
    geometry = json.loads((out / "capped_scene.json").read_text())
    assert [(leg["end"], leg["leg"]) for leg in geometry["exit_paths"]] == [("upper", 1)]
    assert geometry["mouth_rings"] and not geometry["mouth_rings_drawn"] and not geometry["exit_paths_drawn"]
    scene = (out / "capped_volume.pml").read_text()
    assert "crevice_exits" not in scene and 'crevice_scene.get("mouth_rings_drawn")' in scene
    assert "0.16 0.47 0.84" not in (out / "capped_mouths.bild").read_text()
    assert "color 23" not in (out / "capped_volume.tcl").read_text()
    assert "capped_mouths.bild" not in (out / "capped_volume.cxc").read_text()
    assert "color #2 #1494a1" in (out / "capped_volume.cxc").read_text()


def test_publish_guides_are_restored_by_the_options(tmp_path):
    out = _publish_capped(tmp_path, "--exit-centre-lines", "--mouth-guides")
    assert "crevice_exits" in (out / "capped_volume.pml").read_text()
    bild = (out / "capped_mouths.bild").read_text()
    assert "0.16 0.47 0.84" in bild and bild.count(".cylinder") > 64
    tcl = (out / "capped_volume.tcl").read_text()
    assert "color 23" in tcl and "radius 0.065" in tcl
    assert "capped_mouths.bild" in (out / "capped_volume.cxc").read_text()
    geometry = json.loads((out / "capped_scene.json").read_text())
    assert geometry["mouth_rings_drawn"] and geometry["exit_paths_drawn"]


def test_enclosed_chamber_beyond_a_capped_mouth_does_not_block_the_exit():
    # Without lateral exits the chamber above the portal is "another reachable
    # enclosed section beyond a gap", and the run is refused as a fragment.
    frame = capped_channel(lumen_top=6.0)
    with pytest.raises(ChannelResolutionError):
        pore_profile(frame, **OPTIONS)
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    assert profile.metadata["exits"]["upper"]["type"] == "lateral"
    assert profile.points[-1].t < 3.0
    assert profile.metadata["exits"]["upper"]["legs"][0]["exit_point"][0] > 8.0


def test_symmetric_portals_are_both_reported():
    # Two mirror-image side tubes tie exactly; both are distinct exits (exit
    # points far more than twice the bulk radius apart), reported widest
    # first with the same bottleneck, instead of one winning a tie-break.
    frame = capped_channel(second_portal=True)
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    end = profile.metadata["exits"]["upper"]
    assert end["leg_count"] == 2 and [leg["leg"] for leg in end["legs"]] == [1, 2]
    first, second = end["legs"]
    assert {math.copysign(1, first["exit_point"][0]), math.copysign(1, second["exit_point"][0])} == {-1.0, 1.0}
    assert math.dist(first["exit_point"], second["exit_point"]) > end["exit_separation_A"] == 12.0
    assert first["bottleneck_radius_A"] == pytest.approx(second["bottleneck_radius_A"], abs=1e-9)
    assert end["bottleneck_radius_A"] == first["bottleneck_radius_A"]
    spatial = SpatialIndex(frame.atoms)
    for leg in end["legs"]:
        assert all(spatial.segment_clearance(a, b) > 0.8 for a, b in zip(leg["vertices"], leg["vertices"][1:]))
    assert profile.metadata["path_bottleneck_radius_A"] == pytest.approx(
        min(profile.metadata["continuous_bottleneck_radius_A"], first["bottleneck_radius_A"]))


def test_first_leg_is_the_single_exit_search(capped):
    # lateral_exit (single leg) is the first leg of lateral_exit_legs.
    from crevice.channel_exits import lateral_exit, lateral_exit_legs
    from crevice.sections import SectionGeometry, Section
    frame, profile = capped
    md = profile.metadata
    atoms = frame.selected_atoms(include_hydrogen=False, include_hetero=True)
    geometry = SectionGeometry(atoms, profile.axis_origin, profile.axis_direction, spacing=0.5,
                               extent=7.0, probe=0.8, basis=md["section_basis"])
    top = profile.points[-1]
    start = Section(top.t, (), top.position, top.raw_clearance)
    single = lateral_exit(geometry, start, 1, end="upper")
    legs, info = lateral_exit_legs(geometry, start, 1, end="upper")
    assert single.points == legs[0].points and single.bottleneck_clearance == legs[0].bottleneck_clearance
    assert len(legs) == 1 and not info["legs_truncated"] and not info["search_truncated"]


def _captured_axes(monkeypatch, plot, *args, **kwargs):
    """Run a plot function and return the axes of the figure it saves."""
    from crevice import figures
    captured = {}
    original = figures._save_figure

    def keep(fig, path, **options):
        captured["axes"] = list(fig.axes)
        captured["lines"] = [(line.get_label(), line.get_marker(), line.get_linestyle(), line.get_color(),
                              list(line.get_xdata()), list(line.get_ydata())) for line in fig.axes[0].lines]
        captured["collections"] = len(fig.axes[0].collections)
        original(fig, path, **options)

    monkeypatch.setattr(figures, "_save_figure", keep)
    plot(*args, **kwargs)
    return captured


@pytest.mark.parametrize("name", ["plot_profile_radius", "plot_profile_radius_with_residues"])
def test_profile_plots_draw_every_exit_leg_and_no_sample_markers(capped, monkeypatch, tmp_path, name):
    from crevice import figures
    from crevice.presentation import figure_description
    frame = capped_channel(second_portal=True)
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    out = tmp_path / f"{name}.png"
    captured = _captured_axes(monkeypatch, getattr(figures, name), profile, out, dpi=50)
    lines = {label: rest for label, *rest in captured["lines"]}
    # The radius curve is a plain line: no per-sample markers.
    marker, style, color, xs, ys = lines["Pore radius"]
    assert marker in (None, "None", "", " ") and style == "-"
    assert xs == [p.t for p in profile.points]
    # Each distinct leg is a bluish-green dash-dot continuation beyond the upper mouth.
    legs = profile.metadata["exits"]["upper"]["legs"]
    exit_lines = [row for row in captured["lines"] if row[3] == figures.CREVICE_EXIT]
    assert len(exit_lines) == len(legs) == 2
    assert exit_lines[0][0] == "Lateral exit (path length beyond mouth)" and exit_lines[1][0].startswith("_")
    top = profile.points[-1].t
    for (_, marker, style, _, xs, ys), leg in zip(exit_lines, legs):
        assert style == "-." and marker in (None, "None", "", " ")
        assert xs[0] == pytest.approx(top) and xs[-1] == pytest.approx(top + leg["arc_lengths_A"][-1])
        assert ys == leg["radii_A"]
    # The bottleneck dot is kept (a scatter collection besides the fill).
    assert captured["collections"] >= 2
    described = figure_description(out)
    assert "Lateral exit leg 2 at the upper end" in described["Description"]


def test_profile_plot_without_exits_has_no_exit_line(monkeypatch, tmp_path):
    from crevice import figures
    from test_channel_coordinate import rectangular_channel
    profile = pore_profile(rectangular_channel(), samples=21, search_radius=8, enclosure_radius=0.8)
    captured = _captured_axes(monkeypatch, figures.plot_profile_radius, profile, tmp_path / "p.png", dpi=50)
    assert not [row for row in captured["lines"] if row[3] == figures.CREVICE_EXIT]
    assert [row[1] for row in captured["lines"] if row[0] == "Pore radius"][0] in (None, "None", "", " ")


def test_exit_leg_casts_are_disjoint_from_the_channel_and_stop_at_bulk():
    from crevice import void_cast
    from crevice.channel_exits import BulkEnvelope, lateral_exit_casts
    frame = capped_channel()
    profile = pore_profile(frame, lateral_exits=True, **OPTIONS)
    cast = void_cast(frame, profile=profile, spacing=0.5, min_radius=0.0, max_grid_points=400000)
    legs = lateral_exit_casts(frame, profile, cast)
    assert [(c.kind, c.metadata["end"], c.metadata["leg"]) for c in legs] == [("lateral_exit", "upper", 1)]
    channel = {tuple(round(v, 6) for v in p.position) for c in cast.components for p in c.points}
    leg = legs[0]
    assert channel.isdisjoint(tuple(round(v, 6) for v in p.position) for p in leg.points)
    assert all(p.t > cast.metadata["mouth_planes_A"][1] for p in leg.points)
    assert leg.volume == pytest.approx(leg.point_count*0.5**3)
    record = profile.metadata["exits"]["upper"]["legs"][0]
    envelope = BulkEnvelope(frame.atoms, profile.axis_origin, cast.metadata["grid_basis"], profile.points[-1].t, 1,
                            spacing=record["bulk_grid_spacing_A"], bulk_radius=record["bulk_radius_A"])
    assert all(envelope.depth(p.position) > 0 for p in leg.points)
    # No lateral exit, no leg casts: through channels are untouched.
    through = rectangular_channel()
    plain = pore_profile(through, samples=41, search_radius=8, lateral_exits=True)
    assert lateral_exit_casts(through, plain, void_cast(through, profile=plain, spacing=0.5, min_radius=0.0)) == ()

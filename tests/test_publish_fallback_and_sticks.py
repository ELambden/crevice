"""Software/synthetic contracts for publish fallback coverage and stick selection.

No biological channel, cavity or residue role is validated by these fixtures.
"""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
import json
import re
from pathlib import Path

import numpy as np
import pytest

from crevice.presentation import write_viewer_structure
from test_channel_coordinate import rectangular_channel


EXTRA_FILES = ("cavities_json", "cavities_csv", "cavities_pdb", "cavity_summary_png", "tunnels_json", "tunnels_csv",
               "tunnel_points_csv", "tunnels_pdb", "network_json", "network_nodes_csv", "network_edges_csv",
               "network_summary_png", "network_chord_png")


@pytest.fixture(scope="module")
def unresolved_bundle(tmp_path_factory):
    from crevice.figures import write_static_publication_bundle
    root = tmp_path_factory.mktemp("unresolved")
    frame = rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True)
    source = write_viewer_structure(frame, root / "pocket.pdb")
    with pytest.warns(RuntimeWarning, match="No unique through-channel"):
        files = write_static_publication_bundle(frame, enclosure_radius=0.8, structure_path=source, output_dir=root / "published",
            prefix="pocket", cast_spacing=.5, cast_max_grid_points=500000, cast_min_component_volume=30,
            include_hetero=False, dpi=60)
    return frame, files


def test_unresolved_profile_completes_requested_extras_without_a_profile(unresolved_bundle):
    from crevice.residue_evidence import read_binary_dx
    from crevice.rolling import _SphereQueries
    from crevice.radii import atom_vdw_radius
    frame, files = unresolved_bundle
    manifest = json.loads(Path(files["manifest_json"]).read_text())
    for key in EXTRA_FILES:
        assert Path(files[key]).is_file(), key
        assert manifest["files"][key] == files[key]
    # No profile or profile-derived output is invented.
    assert not any(k.startswith(("profile_json", "profile_csv", "profile_radius")) for k in manifest["files"])
    assert manifest["profile_status"]["status"] == "unresolved"
    assert manifest["profile_status"]["reason"] == manifest["single_channel_unresolved_reason"]
    status = json.loads(Path(files["profile_status_json"]).read_text())
    assert status["status"] == "unresolved" and "no pore-radius profile" in status["interpretation"]
    fallback = manifest["fallback_analyses"]
    assert status["fallback_analyses"] == fallback
    assert {name: row["status"] for name, row in fallback.items()} == {
        "cavities": "completed", "tunnels": "completed", "network": "completed"}
    # Tunnel seed is the maximum-clearance sample of the complete measured map.
    grid, origin, deltas = read_binary_dx(files["volume_dx"])
    samples = np.argwhere(grid) @ deltas + origin
    query = _SphereQueries(np.asarray([a.coord for a in frame.atoms]),
                           np.asarray([atom_vdw_radius(a) for a in frame.atoms]))
    clearance, _ = query.points(samples)
    assert fallback["tunnels"]["start"] == pytest.approx(samples[int(np.argmax(clearance))].tolist(), abs=0)
    tunnels = json.loads(Path(files["tunnels_json"]).read_text())["tunnels"]
    assert tunnels and len(tunnels) == fallback["tunnels"]["tunnel_count"]
    assert all(t["start"] == pytest.approx(fallback["tunnels"]["start"]) for t in tunnels)
    assert json.loads(Path(files["cavities_json"]).read_text())["cavities"] == []
    network = json.loads(Path(files["network_json"]).read_text())["network"]
    region = [n for n in network["nodes"] if n["kind"] == "region"]
    assert len(region) == 1 and region[0]["id"] == "region:void:1"
    assert fallback["network"]["region_sample_count"] == len(samples)
    assert any(e["source"] == region[0]["id"] for e in network["edges"])


def test_unresolved_profile_records_skipped_requests(tmp_path):
    from crevice.figures import write_static_publication_bundle
    frame = rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True)
    source = write_viewer_structure(frame, tmp_path / "pocket.pdb")
    with pytest.warns(RuntimeWarning):
        files = write_static_publication_bundle(frame, enclosure_radius=0.8, structure_path=source, output_dir=tmp_path / "out",
            cast_spacing=.5, cast_max_grid_points=500000, cast_min_component_volume=30, include_hetero=False,
            include_cavities=False, include_tunnels=False, include_network=False, dpi=60)
    manifest = json.loads(Path(files["manifest_json"]).read_text())
    assert not any(key in files for key in EXTRA_FILES)
    assert {k: v["status"] for k, v in manifest["fallback_analyses"].items()} == {
        "cavities": "skipped_by_request", "tunnels": "skipped_by_request", "network": "skipped_by_request"}


def test_resolved_profile_keeps_extras_and_reports_resolved(tmp_path):
    from crevice import find_tunnels, pore_profile
    from crevice.figures import write_static_publication_bundle
    frame = rectangular_channel()
    source = write_viewer_structure(frame, tmp_path / "channel.pdb")
    files = write_static_publication_bundle(frame, enclosure_radius=0.8, structure_path=source, output_dir=tmp_path / "out",
        prefix="channel", samples=41, search_radius=8, cast_spacing=.5, dpi=60)
    manifest = json.loads(Path(files["manifest_json"]).read_text())
    for key in EXTRA_FILES + ("profile_json", "profile_csv"):
        assert Path(files[key]).is_file(), key
    assert manifest["profile_status"]["status"] == "resolved"
    assert "profile_status_json" not in files and "fallback_analyses" not in manifest
    profile = json.loads(Path(files["profile_json"]).read_text())
    assert profile["points"] and profile["min_radius"] == pytest.approx(
        pore_profile(frame, enclosure_radius=0.8, samples=41, search_radius=8).min_radius)
    # The resolved path's tunnel call is unchanged: default start and obstacles.
    expected = find_tunnels(frame)
    assert len(json.loads(Path(files["tunnels_json"]).read_text())["tunnels"]) == len(expected)


def test_publish_cli_prints_unresolved_profile(tmp_path, capsys):
    from crevice.cli import main
    frame = rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True)
    source = write_viewer_structure(frame, tmp_path / "pocket.pdb")
    with pytest.warns(RuntimeWarning):
        assert main(["publish", str(source), "--enclosure-radius", "0.8", "--out-dir", str(tmp_path / "out"), "--exclude-hetero",
                     "--cast-spacing", "0.5", "--min-component-volume", "30", "--skip-hydration",
                     "--skip-network", "--dpi", "60"]) == 0
    assert "profile=unresolved" in capsys.readouterr().out
    manifest = json.loads((tmp_path / "out/pocket_manifest.json").read_text())
    assert manifest["fallback_analyses"]["network"]["status"] == "skipped_by_request"
    assert manifest["fallback_analyses"]["tunnels"]["status"] == "completed"
    assert Path(manifest["files"]["input_report_json"]).is_file()


# --- full-list / selected-list residue sticks --------------------------------

@pytest.fixture(scope="module")
def evidence(unresolved_bundle):
    from crevice.residue_evidence import boundary_residue_evidence, read_binary_dx
    frame, files = unresolved_bundle
    grid, origin, deltas = read_binary_dx(files["volume_dx"])
    return frame, files, boundary_residue_evidence(frame, grid, origin, deltas)


def _write(evidence, root, **options):
    from crevice.residue_evidence import write_residue_evidence_bundle
    frame, files, report = evidence
    paths = write_residue_evidence_bundle(report, frame, root, prefix="ctx", scene_path=files["volume_pml"],
                                          dpi=50, **options)
    return paths, json.loads(Path(paths["residue_evidence_json"]).read_text())


def _serials(frame, residues):
    return sorted(a.serial for a in frame.atoms if a.residue_key.label in set(residues))


def _check_scenes(paths, frame, lining, partners):
    membership = {"lining": _serials(frame, lining), "partners": _serials(frame, partners)}
    pml = Path(paths["residue_context_pml"]).read_text()
    assert "crevice_membership = " + repr(membership) in pml
    assert 'cmd.set("stick_transparency", 0.0, crevice_selected)' in pml
    assert "cmd.set_color(\"crevice_boundary_evidence\", [0.85, 0.55, 0.20])" in pml
    tcl = Path(paths["residue_context_tcl"]).read_text()
    for role, ids in membership.items():
        assert f"set crevice_evidence_serials({role}) {{{' '.join(map(str, ids))}}}" in tcl
    assert "material change opacity CREVICEResidues 1.0" in tcl
    cxc = Path(paths["residue_context_cxc"]).read_text()
    for role, color in [("partners", "#a34478"), ("lining", "#d98c33")]:
        residues = lining if role == "lining" else partners
        match = re.search(r"^color (\S+) " + color + " target a$", cxc, flags=re.M)
        if not residues:
            assert match is None
            continue
        viewer_ids = {":" + r.rsplit(":", 1)[1][3:] for r in residues}  # ALA<resid>
        assert {"".join(s.split("/A")[-1:]) for s in match[1].split("|")} == viewer_ids
    assert cxc.count(" stickRadius 0.19") == sum(bool(x) for x in (lining, partners))
    assert "transparency" in cxc and " 0 target a" in cxc
    return membership


def test_default_sticks_remain_top_twelve(evidence, tmp_path):
    frame, _, report = evidence
    paths, saved = _write(evidence, tmp_path)
    lining = [r["residue"] for r in report["residues"] if r["role"] == "boundary_lining"][:12]
    partners = [r["residue"] for r in sorted((r for r in report["residues"] if r["role"] == "nonlocal_lining_partner"),
                key=lambda r: (-r["contacted_lining_boundary_fraction"], r["residue"]))][:12]
    assert saved["display_lining_residues"] == lining and saved["display_partner_residues"] == partners
    assert saved["overlay_definition"].startswith("Top 12 direct contributors gold")
    assert saved["stick_residue_selection"]["mode"] == "top"
    assert saved["display_atom_serials"] == _check_scenes(paths, frame, lining, partners)


def test_all_sticks_export_every_boundary_and_partner_residue(evidence, tmp_path):
    frame, _, report = evidence
    paths, saved = _write(evidence, tmp_path, stick_residues="all")
    lining = {r["residue"] for r in report["residues"] if r["role"] == "boundary_lining"}
    partners = {r["residue"] for r in report["residues"] if r["role"] == "nonlocal_lining_partner"}
    assert len(lining) > 12 and len(partners) > 12
    assert set(saved["display_lining_residues"]) == lining and set(saved["display_partner_residues"]) == partners
    assert saved["stick_residue_selection"] == {"request": "all", "mode": "all", "lining_count": len(lining),
        "partner_count": len(partners), "definition": "All boundary residues and all nonlocal nonlining partners"}
    assert f"{len(lining)} boundary residues gold" in saved["overlay_definition"]
    membership = _check_scenes(paths, frame, sorted(lining), sorted(partners))
    assert saved["display_atom_serials"] == membership
    # The ranked chart/table stay complete and unchanged by the stick choice.
    rows = Path(paths["residue_evidence_csv"]).read_text().splitlines()
    assert len(rows) == len(report["residues"]) + 1


def test_top_n_and_selected_sticks(evidence, tmp_path):
    frame, _, report = evidence
    _, top = _write(evidence, tmp_path / "top", stick_residues="top:3")
    assert len(top["display_lining_residues"]) == 3 and len(top["display_partner_residues"]) == 3
    boundary = [r["residue"] for r in report["residues"] if r["role"] == "boundary_lining"]
    partner = [r["residue"] for r in report["residues"] if r["role"] == "nonlocal_lining_partner"]
    chosen = [boundary[-1], partner[-1], boundary[20]]
    paths, saved = _write(evidence, tmp_path / "sel", stick_residues=",".join(chosen))
    assert saved["display_lining_residues"] == [r for r in boundary if r in chosen]
    assert saved["display_partner_residues"] == [partner[-1]]
    assert saved["stick_residue_selection"]["mode"] == "selected"
    _check_scenes(paths, frame, saved["display_lining_residues"], saved["display_partner_residues"])
    listing = tmp_path / "ids.txt"
    listing.write_text("# chosen residues\n" + "\n".join(chosen) + "\n")
    _, from_file = _write(evidence, tmp_path / "file", stick_residues="@" + str(listing))
    assert from_file["display_atom_serials"] == saved["display_atom_serials"]
    _, from_list = _write(evidence, tmp_path / "list", stick_residues=chosen)
    assert from_list["display_atom_serials"] == saved["display_atom_serials"]


@pytest.mark.parametrize("request_,match", [("top:x", "nonnegative integer"), ("A:ALA999999", "absent"), ("", "empty")])
def test_invalid_stick_selection_fails(evidence, request_, match):
    from crevice.residue_evidence import select_stick_residues
    with pytest.raises(ValueError, match=match):
        select_stick_residues(evidence[2], request_)


def test_selected_residue_outside_both_roles_is_rejected():
    from crevice.residue_evidence import select_stick_residues
    report = {"residues": [
        {"residue": "A:ALA1", "role": "boundary_lining", "contacted_lining_boundary_fraction": 0.},
        {"residue": "A:ALA9", "role": "nonlocal_lining_partner", "contacted_lining_boundary_fraction": .2},
        {"residue": "A:ALA5", "role": "other", "contacted_lining_boundary_fraction": 0.}]}
    with pytest.raises(ValueError, match="neither boundary"):
        select_stick_residues(report, ["A:ALA1", "A:ALA5"])
    lining, partners, record = select_stick_residues(report, "A:ALA9, A:ALA1")
    assert [r["residue"] for r in lining] == ["A:ALA1"] and [r["residue"] for r in partners] == ["A:ALA9"]
    assert record["mode"] == "selected" and select_stick_residues(report, "top:0")[:2] == ([], [])


def test_residue_evidence_cli_stick_option(evidence, tmp_path):
    from crevice.cli import build_parser, main
    assert build_parser().parse_args(["residue-evidence", "x.pdb", "--volume-dx", "v.dx", "--out-dir", "o"]).stick_residues == "top:12"
    frame, files, report = evidence
    assert main(["residue-evidence", files["structure_file"], "--volume-dx", files["volume_dx"],
                 "--out-dir", str(tmp_path), "--prefix", "full", "--scene", files["volume_pml"],
                 "--stick-residues", "all", "--skip-hydration", "--dpi", "50"]) == 0
    saved = json.loads((tmp_path / "full_residue_evidence.json").read_text())
    assert saved["stick_residue_selection"]["mode"] == "all"
    assert len(saved["display_lining_residues"]) == sum(r["role"] == "boundary_lining" for r in report["residues"])
    for extension in ("pml", "tcl", "vmd", "cxc"):
        assert (tmp_path / f"full_residue_context.{extension}").is_file()

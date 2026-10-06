from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import crevice
from crevice import (
    Atom,
    StructureFrame,
    analyze_trajectory,
    annotate_residues,
    build_cavity_network,
    custom_radii_set,
    network_metrics,
    plot_network_chord,
    plot_profile_radius,
    plot_profile_radius_with_residues,
    pore_profile,
    radii_preset,
    result_envelope,
    write_publication_pymol_script,
    write_static_publication_bundle,
    trajectory_from_frames,
    validate_void_cast,
    void_cast,
    write_void_cast_json,
    write_void_cast_pdb,
)
from crevice.cli import main


class CreviceApiTests(unittest.TestCase):
    def test_cavity_network_includes_region_node(self) -> None:
        frame = _channel_frame()
        profile = pore_profile(frame, axis="z", samples=5, search_radius=0.0)
        contacts = annotate_residues(frame, profile.points, cutoff=4.0)

        network = build_cavity_network(frame, profile, contacts=contacts, cutoff=4.0)
        metrics = network_metrics(network)

        self.assertEqual(network.nodes[0].kind, "region")
        self.assertGreater(metrics["node_count"], 1)
        self.assertTrue(any(edge.source == "region:channel" for edge in network.edges))

    def test_static_trajectory_analysis_extracts_features(self) -> None:
        trajectory = trajectory_from_frames([_channel_frame(), _channel_frame()], time_step=2.0)
        result = analyze_trajectory(
            trajectory,
            analyses=("profile", "residues", "network", "features"),
            profile_kwargs={"axis": "z", "samples": 5, "search_radius": 0.0},
            align=False,  # Legacy toy fixture repeats atom identities within residues.
        )

        self.assertEqual(len(result.frames), 2)
        self.assertEqual(result.frames[1].time, 2.0)
        self.assertIn("min_radius", result.frames[0].features)
        self.assertIn("network_density", result.frames[0].features)

    def test_radii_presets_and_schema_envelope(self) -> None:
        preset = radii_preset("hole")
        custom = custom_radii_set(element_radii={"C": 1.9})
        atom = _channel_frame().atoms[0]
        envelope = result_envelope("profile", {"ok": True})

        self.assertEqual(preset.name, "hole")
        self.assertEqual(custom.radius_for_atom(atom), 1.9)
        self.assertEqual(envelope["schema"], "crevice.profile.v0.1")

    def test_cli_network_and_static_trajectory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "mini.pdb"
            network_output = root / "network.json"
            trajectory_output = root / "trajectory.json"
            structure.write_text(_pdb_fixture())

            network_status = main([
                "network",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "--search-radius",
                "0",
                "-o",
                str(network_output),
            ])
            trajectory_status = main([
                "trajectory",
                "--no-align", "--confidence", "0",  # Explicit legacy scan of a nonphysical toy fixture.
                str(structure),
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "--search-radius",
                "0",
                "--analysis",
                "profile",
                "--analysis",
                "features",
                "-o",
                str(trajectory_output),
            ])

            self.assertEqual(network_status, 0)
            self.assertEqual(trajectory_status, 0)
            self.assertIn("network", json.loads(network_output.read_text()))
            self.assertEqual(json.loads(trajectory_output.read_text())["frame_count"], 2)

    def test_void_cast_exports_only_physical_empty_space(self) -> None:
        # Isolated atoms test clearance, but do not establish channel enclosure.
        from test_channel_coordinate import rectangular_channel
        frame = rectangular_channel()
        cast = void_cast(
            frame,
            mode="channel",
            axis="z",
            spacing=1.0,
            min_radius=1.0,
            max_components=2,
            centerline_search_radius=8.0,
            max_grid_points=20_000,
        )
        validation = validate_void_cast(frame, cast)

        self.assertGreater(cast.point_count, 0)
        self.assertEqual(cast.metadata["export_mode"], "connected_section_fill")
        self.assertTrue(validation["is_physical"])
        self.assertEqual(validation["violation_count"], 0)
        self.assertGreaterEqual(validation["min_clearance"], cast.min_radius - 1e-6)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            json_path = root / "void_cast.json"
            pdb_path = root / "void_cast.pdb"
            write_void_cast_json(cast, json_path)
            write_void_cast_pdb(cast, pdb_path)

            data = json.loads(json_path.read_text())
            self.assertEqual(data["method"], "clearance-grid-fill")
            self.assertIn("components", data)
            self.assertIn("HETATM", pdb_path.read_text())

    def test_publication_profile_png_and_pymol_script(self) -> None:
        from PIL import Image

        frame = _channel_frame()
        profile = pore_profile(frame, axis="z", samples=5, search_radius=0.0)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "profile.png"
            pml = root / "pore_cast.pml"
            pdb = root / "pore_cast.pdb"
            structure = root / "mini.pdb"
            structure.write_text(_pdb_fixture())
            crevice.write_profile_pdb(profile, pdb)

            plot_profile_radius(profile, png, dpi=400)
            plot_profile_radius_with_residues(profile, root / "profile_annotated.png", dpi=400)
            write_publication_pymol_script(
                structure_path=structure,
                pore_cast_path=pdb,
                script_path=pml,
                profile=profile,
                output_png=root / "pore_cast.png",
                dpi=400,
            )

            self.assertTrue(png.exists())
            with Image.open(png) as image:
                self.assertEqual(image.format, "PNG")
                self.assertGreaterEqual(round(image.info.get("dpi", (0, 0))[0]), 399)
            script = pml.read_text()
            self.assertIn("cmd.load(_crevice_input", script)
            self.assertIn("cmd.ray(2400, 1800)", script)
            self.assertIn("cmd.png(_png_target, dpi=400)", script)
            self.assertIn("crevice_cartoon", script)
            self.assertIn("color crevice_cyan, crevice_cartoon", script)
            self.assertIn("crevice_pore_fill", script)
            self.assertIn("show surface, crevice_pore_fill", script)
            self.assertIn("color crevice_pink, crevice_pore_fill", script)
            self.assertIn("crevice_bottleneck", script)
            self.assertIn("spectrum b, blue_white_red", script)
            self.assertTrue((root / "profile_annotated.png").exists())


    def test_publication_summary_figures_handle_empty_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cavity_png = root / "cavities.png"
            chord_png = root / "network_chord.png"
            frame = _channel_frame()
            profile = pore_profile(frame, axis="z", samples=5, search_radius=0.0)
            contacts = annotate_residues(frame, profile.points, cutoff=4.0)
            network = build_cavity_network(frame, profile, contacts=contacts, cutoff=4.0)

            crevice.plot_cavity_summary((), cavity_png, dpi=400)
            plot_network_chord(network, chord_png, dpi=400)

            self.assertTrue(cavity_png.exists())
            self.assertTrue(chord_png.exists())

    def test_cli_profile_publication_outputs_and_publish_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "mini.pdb"
            output = root / "profile.json"
            csv_output = root / "profile.csv"
            pore_cast = root / "pore_cast.pdb"
            profile_png = root / "profile.png"
            annotated_png = root / "profile_annotated.png"
            pml = root / "pore_cast.pml"
            network_output = root / "network.json"
            chord_png = root / "network_chord.png"
            bundle_dir = root / "bundle"
            structure.write_text(_pdb_fixture())

            profile_status = main([
                "profile",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "--search-radius",
                "0",
                "-o",
                str(output),
                "--csv",
                str(csv_output),
                "--pdb",
                str(pore_cast),
                "--png",
                str(profile_png),
                "--annotated-png",
                str(annotated_png),
                "--pymol",
                str(pml),
                "--pymol-png",
                str(root / "render.png"),
            ])
            network_status = main([
                "network",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "--search-radius",
                "0",
                "-o",
                str(network_output),
                "--chord-png",
                str(chord_png),
            ])
            publish_status = main([
                "publish",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "--search-radius",
                "0",
                "--out-dir",
                str(bundle_dir),
                "--prefix",
                "mini",
                "--skip-cavities",
                "--skip-tunnels",
            ])

            self.assertEqual(profile_status, 0)
            self.assertEqual(network_status, 0)
            self.assertEqual(publish_status, 0)
            self.assertTrue(profile_png.exists())
            self.assertTrue(annotated_png.exists())
            self.assertTrue(chord_png.exists())
            self.assertIn("cmd.ray(2400, 1800)", pml.read_text())
            manifest = json.loads((bundle_dir / "mini_manifest.json").read_text())
            self.assertIn("structure_file", manifest["files"])
            self.assertIn("profile_radius_png", manifest["files"])
            self.assertIn("profile_radius_annotated_png", manifest["files"])
            self.assertIn("network_chord_png", manifest["files"])
            self.assertTrue((bundle_dir / "mini_structure.pdb").exists())
            self.assertTrue((bundle_dir / "mini_profile_radius.png").exists())
            self.assertTrue((bundle_dir / "mini_profile_radius_annotated.png").exists())
            self.assertTrue((bundle_dir / "mini_network_chord.png").exists())
            bundle_pml = (bundle_dir / "mini_pore_cast.pml").read_text()
            self.assertIn("mini_viewer.pdb", bundle_pml)
            self.assertIn("mini_pore_cast.pdb", bundle_pml)
            self.assertIn("crevice_cartoon", bundle_pml)
            self.assertIn("crevice_pore_fill", bundle_pml)
            self.assertIn("color crevice_cyan, crevice_cartoon", bundle_pml)
            self.assertIn("color crevice_pink, crevice_pore_fill", bundle_pml)
            self.assertNotIn("results/", bundle_pml)

    def test_cli_static_suite_generates_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache"
            out_dir = root / "suite"
            cache.mkdir()
            (cache / "1GRM.pdb").write_text(_pdb_fixture())

            status = main([
                "static-suite",
                "1GRM",
                "--cache-dir",
                str(cache),
                "--out-dir",
                str(out_dir),
                # The cache holds PDB-format files; mmCIF is the default.
                "--format",
                "pdb",
                "--samples",
                "5",
                "--search-radius",
                "0",
                "--probe-radius",
                "0",
                "--skip-cavities",
                "--skip-tunnels",
            ])

            self.assertEqual(status, 0)
            summary = json.loads((out_dir / "static_benchmark_summary.json").read_text())
            self.assertEqual(summary["systems"][0]["status"], "ok")
            self.assertTrue((out_dir / "1GRM" / "1GRM_manifest.json").exists())
            self.assertTrue((out_dir / "1GRM" / "1GRM_network_chord.png").exists())


def _channel_frame() -> StructureFrame:
    atoms: list[Atom] = []
    serial = 1
    residues = [
        ("LEU", "CD1", (5.0, 0.0)),
        ("ASP", "OD1", (-5.0, 0.0)),
        ("LYS", "NZ", (0.0, 5.0)),
        ("PHE", "CZ", (0.0, -5.0)),
    ]
    for z in (-4.0, 0.0, 4.0):
        for resid, (resname, name, xy) in enumerate(residues, start=1):
            atoms.append(
                Atom(
                    serial=serial,
                    name=name,
                    resname=resname,
                    chain_id="A",
                    resid=resid,
                    x=xy[0],
                    y=xy[1],
                    z=z,
                    element=name[0],
                )
            )
            serial += 1
    return StructureFrame(tuple(atoms), source="channel")


def _pdb_fixture() -> str:
    lines = []
    serial = 1
    for z in (-2.0, 2.0):
        for resid, (resname, name, x, y) in enumerate(
            [
                ("LEU", "CD1", 5.0, 0.0),
                ("ASP", "OD1", -5.0, 0.0),
                ("LYS", "NZ", 0.0, 5.0),
                ("PHE", "CZ", 0.0, -5.0),
            ],
            start=1,
        ):
            element = "C" if name.startswith("C") else name[0]
            lines.append(
                f"ATOM  {serial:5d} {name:>4s} {resname:>3s} A{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{10.00:6.2f}          {element:>2s}"
            )
            serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    unittest.main()


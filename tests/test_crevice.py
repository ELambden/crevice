from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import crevice
from channel_fixtures import channel_frame, channel_pdb_text

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
        profile = pore_profile(frame, axis="z", samples=5)
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
            profile_kwargs={"axis": "z", "samples": 5},
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
            structure.write_text(_pdb_fixture(), encoding="utf-8", newline="\n")

            network_status = main([
                "residues", "--network",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
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
                "--analysis",
                "profile",
                "--analysis",
                "features",
                "-o",
                str(trajectory_output),
            ])

            self.assertEqual(network_status, 0)
            self.assertEqual(trajectory_status, 0)
            self.assertIn("network", json.loads(network_output.read_text(encoding="utf-8")))
            self.assertEqual(json.loads(trajectory_output.read_text(encoding="utf-8"))["frame_count"], 2)

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

            data = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertEqual(data["method"], "clearance-grid-fill")
            self.assertIn("components", data)
            self.assertIn("HETATM", pdb_path.read_text(encoding="utf-8"))

    def test_publication_profile_png_and_pymol_script(self) -> None:
        from PIL import Image

        frame = _channel_frame()
        profile = pore_profile(frame, axis="z", samples=5)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "profile.png"
            pml = root / "pore_cast.pml"
            pdb = root / "pore_cast.pdb"
            structure = root / "mini.pdb"
            structure.write_text(_pdb_fixture(), encoding="utf-8", newline="\n")
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
            script = pml.read_text(encoding="utf-8")
            self.assertIn("cmd.load(crevice_input('mini.pdb')", script)
            self.assertNotIn(str(root), script)
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
            profile = pore_profile(frame, axis="z", samples=5)
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
            structure.write_text(_pdb_fixture(), encoding="utf-8", newline="\n")

            profile_status = main([
                "profile",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
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
                "residues", "--network",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "-o",
                str(network_output),
                "--chord-png",
                str(chord_png),
            ])
            publish_status = main([
                "report",
                str(structure),
                "--axis",
                "z",
                "--samples",
                "5",
                "--out-dir",
                str(bundle_dir),
                "--prefix",
                "mini",
                "--skip-hydration",
                "--all",   # 0.3: the complete bundle is opt-in
            ])

            self.assertEqual(profile_status, 0)
            self.assertEqual(network_status, 0)
            self.assertEqual(publish_status, 0)
            self.assertTrue(profile_png.exists())
            self.assertTrue(annotated_png.exists())
            self.assertTrue(chord_png.exists())
            self.assertIn("cmd.ray(2400, 1800)", pml.read_text(encoding="utf-8"))
            manifest = json.loads((bundle_dir / "mini_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("structure_file", manifest["files"])
            self.assertIn("profile_radius_png", manifest["files"])
            self.assertIn("profile_radius_annotated_png", manifest["files"])
            self.assertIn("network_chord_png", manifest["files"])
            self.assertTrue((bundle_dir / "mini_structure.pdb").exists())
            self.assertTrue((bundle_dir / "mini_profile_radius.png").exists())
            self.assertTrue((bundle_dir / "mini_profile_radius_annotated.png").exists())
            self.assertTrue((bundle_dir / "mini_network_chord.png").exists())
            bundle_pml = (bundle_dir / "mini_pore_cast.pml").read_text(encoding="utf-8")
            # A resolved channel: the pore-cast script draws the measured lumen mesh.
            self.assertIn("mini_viewer.pdb", bundle_pml)
            self.assertIn("mini_pymol_mesh.npz", bundle_pml)
            self.assertIn("crevice_surface('crevice_lumen', 'lumen')", bundle_pml)
            self.assertNotIn("results/", bundle_pml)

    def test_cli_static_suite_generates_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache"
            out_dir = root / "suite"
            cache.mkdir()
            (cache / "1GRM.pdb").write_text(_pdb_fixture(), encoding="utf-8", newline="\n")

            status = main([
                "report", "--suite",
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
                "--probe-radius",
                "0",
                "--skip-cavities",
            ])

            self.assertEqual(status, 0)
            summary = json.loads((out_dir / "static_benchmark_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["systems"][0]["status"], "ok")
            self.assertTrue((out_dir / "1GRM" / "1GRM_manifest.json").exists())
            self.assertTrue((out_dir / "1GRM" / "1GRM_network_chord.png").exists())


def _channel_frame() -> StructureFrame:
    return channel_frame()


def _pdb_fixture() -> str:
    return channel_pdb_text()


if __name__ == "__main__":
    unittest.main()

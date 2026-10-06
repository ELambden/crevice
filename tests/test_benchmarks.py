from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from crevice.benchmarks import (
    DEFAULT_BENCHMARK_IDS,
    benchmark_ids_with_status,
    benchmark_path,
    benchmark_systems,
    curated_benchmark_ids,
    fetch_benchmark_structure,
    get_benchmark_system,
)
from crevice.parser import load_structure
from crevice.static_suite import _cast_mode_for_expected


class BenchmarkMetadataTests(unittest.TestCase):
    def test_expected_benchmark_ids_are_registered(self) -> None:
        self.assertEqual(DEFAULT_BENCHMARK_IDS, ("1GRM", "2CHB", "6MVY", "1AF6", "4PYP"))
        systems = benchmark_systems()
        self.assertEqual(len(systems), 5)
        self.assertEqual(get_benchmark_system("2chb").pdb_id, "2CHB")

    def test_benchmark_path_normalizes_extension(self) -> None:
        path = benchmark_path("cache", "1grm", extension="mmcif")
        self.assertEqual(path, Path("cache") / "1GRM.cif")

    def test_misidentified_entries_no_longer_carry_their_old_labels(self) -> None:
        """4PYP is human GLUT1; 2CHB is a cholera toxin B-pentamer complex."""
        glut1 = get_benchmark_system("4pyp")
        self.assertNotIn("photoactive", glut1.name.lower())
        self.assertIn("glut1", glut1.name.lower())
        toxin = get_benchmark_system("2chb")
        self.assertNotIn("channel benchmark system", toxin.name)

    def test_withdrawn_geometry_labels_fall_back_to_automatic_mode(self) -> None:
        """A label derived from a misidentification is withdrawn, not reassigned."""
        for pdb_id in ("4PYP", "2CHB"):
            with self.subTest(pdb_id=pdb_id):
                system = get_benchmark_system(pdb_id)
                self.assertEqual(system.expected_geometry, "unassigned")
                self.assertEqual(_cast_mode_for_expected(system.expected_geometry), "auto")

    def test_no_benchmark_system_is_curated_yet(self) -> None:
        """Curation gates biological interpretation; none has been done."""
        self.assertEqual(curated_benchmark_ids(), ())
        self.assertEqual(benchmark_ids_with_status("set_aside"), ("2CHB",))
        self.assertEqual(
            set(benchmark_ids_with_status("uncurated")),
            {"1GRM", "6MVY", "1AF6", "4PYP"},
        )

    def test_unknown_curation_status_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "status must be one of"):
            benchmark_ids_with_status("validated")


@unittest.skipUnless(
    os.environ.get("CREVICE_RUN_PDB_INTEGRATION") == "1",
    "set CREVICE_RUN_PDB_INTEGRATION=1 to run real PDB integration checks",
)
class PDBIntegrationTests(unittest.TestCase):
    def test_requested_pdb_systems_load(self) -> None:
        cache_dir = os.environ.get("CREVICE_PDB_DIR")
        fetch = os.environ.get("CREVICE_FETCH_PDB") == "1"
        if cache_dir is None:
            if not fetch:
                self.skipTest("set CREVICE_PDB_DIR or CREVICE_FETCH_PDB=1")
            temp = tempfile.TemporaryDirectory()
            self.addCleanup(temp.cleanup)
            cache_dir = temp.name

        for pdb_id in DEFAULT_BENCHMARK_IDS:
            path = benchmark_path(cache_dir, pdb_id)
            if not path.exists():
                if not fetch:
                    self.skipTest(f"{path} is missing; set CREVICE_FETCH_PDB=1 to download")
                path = fetch_benchmark_structure(pdb_id, cache_dir)
            with self.subTest(pdb_id=pdb_id):
                frame = load_structure(path)
                self.assertGreater(len(frame.atoms), 0)
                self.assertGreater(len(frame.residues()), 0)


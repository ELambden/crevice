"""Independent scientific contract tests.

Each test states a geometric or intake property with a known answer (analytic
clearances, model and chain selection, tunnel and network definitions) and
checks CREVICE against it without relying on CREVICE's own reference outputs.
"""

import math
import random
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from crevice import Atom, StructureFrame, load_structure, pore_profile
from crevice.analysis import annotate_residues, find_tunnels
from crevice.channels import profile_along_path
from crevice.models import ChannelPoint, PoreProfile
from crevice.networks import build_residue_network
from crevice.radii import atom_vdw_radius
from crevice.spatial import SpatialIndex


def atom(serial, xyz, element="C", **kwargs):
    return Atom(serial=serial, name=kwargs.pop("name", element),
                resname=kwargs.pop("resname", "ALA"),
                chain_id=kwargs.pop("chain_id", "A"),
                resid=kwargs.pop("resid", serial),
                x=xyz[0], y=xyz[1], z=xyz[2], element=element, **kwargs)


def ring(radius=5.0):
    # At ring planes, exact axial clearance is radius - carbon VDW radius.
    return StructureFrame(tuple(
        atom(1 + iz * 32 + j,
             (radius * math.cos(j * math.tau / 32),
              radius * math.sin(j * math.tau / 32), z))
        for iz, z in enumerate((-4.0, 0.0, 4.0)) for j in range(32)))


def pdb_line(x=5.0):
    return (f"ATOM  {1:5d} {'CA':>4s} {'ALA':>3s} A{1:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}          {'C':>2s}\n")


class ScientificContracts(unittest.TestCase):
    def test_atomistic_ring_clearance(self):
        profile = pore_profile(ring(), axis="z", samples=3, search_radius=0)
        for point in profile.points:
            self.assertAlmostEqual(point.radius, 3.3, places=10)

    def test_translation_invariance(self):
        frame = ring()
        moved = StructureFrame(tuple(replace(a, x=a.x + 31, y=a.y - 17, z=a.z + 11)
                                     for a in frame.atoms))
        left = pore_profile(frame, axis="z", samples=9, search_radius=0)
        right = pore_profile(moved, axis="z", samples=9, search_radius=0)
        for a, b in zip(left.points, right.points):
            self.assertAlmostEqual(a.radius, b.radius, places=10)

    def test_explicit_axis_rotation_invariance(self):
        frame = ring()
        moved = StructureFrame(tuple(replace(a, x=a.z, z=-a.x) for a in frame.atoms))
        left = pore_profile(frame, axis="z", samples=9, search_radius=0)
        right = pore_profile(moved, axis="x", samples=9, search_radius=0)
        for a, b in zip(left.points, right.points):
            self.assertAlmostEqual(a.radius, b.radius, places=10)

    def test_surface_query_matches_independent_bruteforce(self):
        rng = random.Random(42)
        atoms = tuple(atom(i, tuple(rng.uniform(-10, 10) for _ in range(3)),
                           element=rng.choice(["C", "O", "S", "K"]))
                      for i in range(70))
        index = SpatialIndex(atoms)
        for _ in range(30):
            point = tuple(rng.uniform(-10, 10) for _ in range(3))
            expected = min(math.dist(point, a.coord) - atom_vdw_radius(a) for a in atoms)
            self.assertAlmostEqual(index.nearest_surface(point)[0], expected, places=10)
        index._tree = None
        self.assertAlmostEqual(index.nearest_surface(point)[0], expected, places=10)

    def test_repeated_serials_do_not_change_atomic_radii(self):
        atoms = (atom(1, (0, 0, 0), "C"), atom(1, (20, 0, 0), "O"))
        self.assertAlmostEqual(SpatialIndex(atoms).nearest_surface((3, 0, 0))[0], 1.3)

    def test_closed_slice_is_retained(self):
        frame = ring()
        closed = StructureFrame(frame.atoms + (atom(1000, (0, 0, 0)),))
        profile = pore_profile(closed, axis="z", samples=3, search_radius=0)
        self.assertEqual(profile.min_radius, 0.0)

    def test_probe_subtraction_contract(self):
        profile = pore_profile(ring(), axis="z", samples=3, search_radius=0, probe_radius=1.4)
        self.assertAlmostEqual(profile.min_radius, 1.9)
        self.assertAlmostEqual(profile.bottleneck.raw_clearance, 3.3)

    def test_negative_probe_is_rejected(self):
        with self.assertRaises(ValueError):
            pore_profile(ring(), probe_radius=-1.0)

    def test_constant_radius_volume(self):
        points = tuple(ChannelPoint(i, (0, 0, float(i)), 2.0, 2.0, float(i), 1, "A:ALA1")
                       for i in range(3))
        profile = PoreProfile((0, 0, 0), (0, 0, 1), points)
        self.assertAlmostEqual(profile.volume_estimate, math.pi * 4 * 2)

    def test_pdb_model_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "models.pdb"
            path.write_text("MODEL        1\n" + pdb_line(5) + "ENDMDL\nMODEL        2\n" +
                            pdb_line(50) + "ENDMDL\nEND\n")
            first = load_structure(path, model_index=0)
            second = load_structure(path, model_index=1)
            self.assertEqual([a.x for a in first.atoms], [5.0])
            self.assertEqual([a.x for a in second.atoms], [50.0])

    def test_mmcif_model_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "models.cif"
            tags = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id",
                    "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
                    "pdbx_PDB_model_num"]
            path.write_text("data_test\nloop_\n" +
                            "".join("_atom_site." + tag + "\n" for tag in tags) +
                            "ATOM 1 C CA ALA A 1 5 0 0 1\nATOM 2 C CA ALA A 1 50 0 0 2\n#\n")
            self.assertEqual(len(load_structure(path, model_index=0).atoms), 1)

    def test_malformed_coordinates_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "invalid.pdb"
            line = pdb_line()
            path.write_text(line[:30] + " invalid" + line[38:])
            with self.assertRaises(ValueError):
                load_structure(path)

    def test_chain_and_insertion_identity(self):
        frame = StructureFrame((atom(1, (0, 0, 0), resid=1),
                                atom(2, (3, 0, 0), resid=1, chain_id="B"),
                                atom(3, (0, 3, 0), resid=1, icode="A")))
        graph = build_residue_network(frame)
        self.assertEqual(len({n.id for n in graph.nodes}), 3)

    def test_network_cutoff_is_surface_distance(self):
        frame = StructureFrame((atom(1, (0, 0, 0)), atom(2, (7.8, 0, 0))))
        graph = build_residue_network(frame, cutoff=4.5)
        self.assertEqual(len(graph.edges), 1)
        self.assertAlmostEqual(graph.edges[0].distance, 4.4)

    def test_backbone_only_pair_is_not_a_salt_bridge(self):
        frame = StructureFrame((atom(1, (0, 0, 0), resname="LYS", name="CA"),
                                atom(2, (4, 0, 0), resname="ASP", name="CA")))
        graph = build_residue_network(frame)
        self.assertFalse(any("salt_bridge" in edge.interaction for edge in graph.edges))

    def test_wide_pore_lining_is_not_lost(self):
        frame = ring(radius=10.0)
        profile = pore_profile(frame, axis="z", samples=3, search_radius=0)
        contacts = annotate_residues(frame, profile.points)
        self.assertGreater(len(contacts), 0, "All ring atoms bound the pore but are beyond a centerline cutoff")

    def test_tunnel_start_inside_atom_is_rejected(self):
        frame = StructureFrame((atom(1, (0, 0, 0)),))
        try:
            tunnels = find_tunnels(frame, start=(0, 0, 0), spacing=1, padding=4,
                                   min_radius=0.8, max_tunnels=1)
        except ValueError:
            return
        self.assertEqual(tunnels, (), "An inaccessible start must not teleport to an exterior grid point")

    def test_explicit_path_exact_clearance(self):
        frame = StructureFrame((atom(1, (0, 0, 0)),))
        profile = profile_along_path(frame, [(3, 0, 0), (4, 0, 0)])
        self.assertAlmostEqual(profile.points[0].radius, 1.3)
        self.assertAlmostEqual(profile.points[1].radius, 2.3)


if __name__ == "__main__":
    unittest.main(verbosity=2)

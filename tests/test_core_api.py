from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from crevice import (
    Atom,
    StructureFrame,
    annotate_residues,
    cluster_profiles,
    detect_cavities,
    find_tunnels,
    load_structure,
    pore_profile,
    profile_features,
    representative_profile,
    write_profile_csv,
    write_profile_json,
    write_profile_pdb,
    write_residue_contacts_csv,
)
from crevice.cli import main


class CoreApiTests(unittest.TestCase):
    def test_loads_pdb_and_groups_residues(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mini.pdb"
            path.write_text(_pdb_fixture())

            frame = load_structure(path)

            self.assertEqual(len(frame.atoms), 8)
            self.assertEqual(len(frame.residues()), 4)
            self.assertEqual(frame.atoms[0].element, "C")

    def test_pore_profile_and_residue_contacts(self) -> None:
        frame = _channel_frame()
        profile = pore_profile(frame, axis="z", samples=5, search_radius=0.0)
        contacts = annotate_residues(frame, profile.points, cutoff=4.0)

        self.assertEqual(len(profile.points), 5)
        self.assertGreater(profile.min_radius, 3.0)
        self.assertTrue(any(contact.role in {"lining", "bottleneck", "bottleneck-nearby"} for contact in contacts))
        self.assertTrue(any("hydrophobic" in contact.properties for contact in contacts))

    def test_sparse_box_is_an_open_candidate_not_an_enclosed_void(self) -> None:
        frame = _box_frame()
        self.assertEqual(detect_cavities(frame, spacing=2.0, min_radius=1.0), ())
        cavities = detect_cavities(frame, spacing=2.0, min_radius=1.0, max_cavities=3, enclosed_only=False)

        self.assertGreaterEqual(len(cavities), 1)
        self.assertGreater(cavities[0].volume, 0.0)

    def test_tunnel_search_reaches_boundary(self) -> None:
        frame = _channel_frame()
        tunnels = find_tunnels(frame, start=(0.0, 0.0, 0.0), spacing=2.5, min_radius=1.0, max_tunnels=1)

        self.assertEqual(len(tunnels), 1)
        self.assertGreater(tunnels[0].length, 0.0)
        self.assertGreater(tunnels[0].bottleneck_radius, 1.0)

    def test_io_and_ml_helpers(self) -> None:
        frame = _channel_frame()
        profile = pore_profile(frame, axis="z", samples=5, search_radius=0.0)
        contacts = annotate_residues(frame, profile.points, cutoff=4.0)
        features = profile_features(profile, contacts)
        labels = cluster_profiles([profile, profile], k=2)
        index, representative = representative_profile([profile, profile])

        self.assertIn("min_radius", features)
        self.assertEqual(labels, [0, 0])
        self.assertEqual(index, 0)
        self.assertIs(representative, profile)

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            write_profile_json(profile, out / "profile.json")
            write_profile_csv(profile, out / "profile.csv")
            write_profile_pdb(profile, out / "profile.pdb")
            write_residue_contacts_csv(contacts, out / "contacts.csv")

            data = json.loads((out / "profile.json").read_text())
            self.assertEqual(data["point_count"], 5)
            self.assertIn("HETATM", (out / "profile.pdb").read_text())
            self.assertIn("residue", (out / "contacts.csv").read_text().splitlines()[0])

    def test_cli_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "mini.pdb"
            output = root / "profile.json"
            csv_output = root / "profile.csv"
            pdb_output = root / "profile.pdb"
            structure.write_text(_pdb_fixture())

            status = main([
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
                str(pdb_output),
            ])

            self.assertEqual(status, 0)
            self.assertTrue(output.exists())
            self.assertTrue(csv_output.exists())
            self.assertTrue(pdb_output.exists())

    def test_cli_profile_writes_its_csv_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "mini.pdb"
            structure.write_text(_pdb_fixture())
            status = main(["profile", str(structure), "--axis", "z", "--samples", "5", "--search-radius", "0",
                           "-o", str(root / "out" / "run.json")])
            self.assertEqual(status, 0)
            header = (root / "out" / "run.csv").read_text().splitlines()[0]
            self.assertTrue(header.startswith("index,x_A,y_A,z_A,radius_A"))
            self.assertEqual(len((root / "out" / "run.csv").read_text().splitlines()), 6)


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


def _box_frame() -> StructureFrame:
    atoms: list[Atom] = []
    serial = 1
    for x in (-4.0, 4.0):
        for y in (-4.0, 4.0):
            for z in (-4.0, 4.0):
                atoms.append(
                    Atom(
                        serial=serial,
                        name="C",
                        resname="ALA",
                        chain_id="A",
                        resid=serial,
                        x=x,
                        y=y,
                        z=z,
                        element="C",
                    )
                )
                serial += 1
    return StructureFrame(tuple(atoms), source="box")


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
            lines.append(_pdb_atom(serial, name, resname, "A", resid, x, y, z, element))
            serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def _pdb_atom(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"ATOM  {serial:5d} {name:>4s} {resname:>3s} {chain}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{10.00:6.2f}          {element:>2s}"
    )


if __name__ == "__main__":
    unittest.main()


"""Atoms without a radius in the set in force, and AlphaFold DB identifiers (offline).

Software/synthetic checks: hand-written PDB records (a designed residue with a
gadolinium atom, a modified residue keyed by RESNAME:ATOMNAME) and synthetic
AlphaFold payloads served through a monkeypatched downloader.
"""
import json

import pytest

from crevice import pdb as pdbmod
from crevice.channels import profile_along_path
from crevice.models import Atom, StructureFrame
from crevice.parser import load_structure_report, predicted_model_report
from crevice.radii import (UnrecognisedElementWarning, atom_vdw_radius, load_radius_file, radii_preset,
                           use_radii, write_radius_template)


def line(serial, name, resname, element, x, record="HETATM", resid=1):
    return (f"{record:<6}{serial:5d} {name:<4} {resname:>3} A{resid:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          {element:>2}\n")


@pytest.fixture
def designed_pdb(tmp_path):
    """Alanine plus a designed residue XGD carrying a Gd atom (no tabulated radius)."""
    path = tmp_path / "designed.pdb"
    path.write_text(line(1, "N", "ALA", "N", 0.0, "ATOM") + line(2, "CA", "ALA", "C", 1.458, "ATOM")
                    + line(3, "C1", "XGD", "C", 4.0, resid=2) + line(4, "GD1", "XGD", "GD", 6.0, resid=2) + "END\n")
    return path


def test_designed_residue_without_radius_is_excluded_and_reported(designed_pdb):
    with pytest.warns(UnrecognisedElementWarning) as caught:
        frame, report = load_structure_report(designed_pdb)
    message = str(caught[0].message)
    for hint in ("GD (1 atom)", "XGD:GD1", "--radii FILE", "radii=", "write_radius_template", "RESNAME:ATOMNAME"):
        assert hint in message
    assert [a.name for a in frame.atoms] == ["N", "CA", "C1"]
    assert report["excluded_no_radius_atoms"] == {"GD": 1}
    assert report["excluded_no_radius_names"] == {"XGD:GD1": 1}
    assert report["unsupported_radius_elements"] == ["GD"]
    assert report["radius_set"] == "default" and report["atom_count"] == 3


def test_user_radius_file_includes_the_element_and_is_recorded(designed_pdb, tmp_path):
    radii = write_radius_template(tmp_path / "gd.csv", elements=["GD"])
    radii.write_text(radii.read_text().replace("# element,GD,<radius in angstrom>", "element,GD,2.10"))
    frame, report = load_structure_report(designed_pdb, radii=radii)
    assert [a.name for a in frame.atoms] == ["N", "CA", "C1", "GD1"]
    assert report["excluded_no_radius_atoms"] == {} and report["radius_set"] == "gd"
    profile = profile_along_path(frame, [(6.0, 3.0, 0.0)], radii=radii)
    recorded = profile.metadata["radii"]
    assert recorded["element_overrides_A"] == {"GD": 2.1}
    assert recorded["source_sha256"] and recorded["source"].endswith("gd.csv")
    assert profile.points[0].raw_clearance == pytest.approx(3.0 - 2.1)
    # Same exclusion decision when the set is chosen for a whole block (as --radii does).
    with use_radii(radii):
        assert len(load_structure_report(designed_pdb)[0].atoms) == 4


def test_modified_residue_keyed_by_resname_atomname(designed_pdb, tmp_path):
    path = tmp_path / "modified.json"
    path.write_text(json.dumps({"base": "default", "atom_radii": {"XGD:GD1": 2.30}}))
    frame, report = load_structure_report(designed_pdb, radii=path)
    gd = frame.atoms[-1]
    assert gd.name == "GD1" and atom_vdw_radius(gd, path) == 2.30
    assert load_radius_file(path).provenance()["atom_radii_A"] == {"XGD:GD1": 2.3}
    # The key is residue-specific: a GD1 atom in another residue still has no radius.
    assert not load_radius_file(path).has_radius("GD1", "OTH", "GD")


def test_has_radius_follows_the_set_precedence():
    from crevice.radii import RadiusSet
    carbon_only = RadiusSet("carbon-only", hole_records=(("C???", "???", 1.85),))
    assert carbon_only.has_radius("CA", "ALA", "C") and not carbon_only.has_radius("GD1", "XGD", "GD")
    assert radii_preset("default").has_radius("SOD", "SOD", "")  # recognised ion residue
    assert not radii_preset("default").has_radius("GD1", "XGD", "GD")


def test_deuterium_is_treated_as_hydrogen():
    atom = Atom(1, "D1", "DOD", "W", 1, 0.0, 0.0, 0.0, "D")
    assert atom_vdw_radius(atom) == atom_vdw_radius(Atom(2, "H1", "HOH", "W", 1, 0.0, 0.0, 0.0, "H")) == 1.2
    frame = StructureFrame((atom, Atom(3, "O", "DOD", "W", 1, 1.0, 0.0, 0.0, "O")))
    assert [a.name for a in frame.selected_atoms()] == ["O"]
    assert len(frame.selected_atoms(include_hydrogen=True)) == 2


def test_md_selection_with_radius_less_atom_is_refused(designed_pdb):
    pytest.importorskip("MDAnalysis")
    from crevice.trajectory import load_trajectory_report
    with pytest.raises(ValueError, match=r"GD \(1 atom\).*\(all\) and not resname XGD"):
        load_trajectory_report(designed_pdb, selection="all")
    _, report = load_trajectory_report(designed_pdb, selection="(all) and not resname XGD")
    assert report["selected_atom_count"] == 2


def test_radius_template_requires_csv(tmp_path):
    with pytest.raises(ValueError, match="csv"):
        write_radius_template(tmp_path / "radii.json")


# --- AlphaFold DB ------------------------------------------------------------

AF_PDB = (
    "HEADER                                            01-AUG-25                     \n"
    "TITLE     ALPHAFOLD MONOMER V2.0 PREDICTION FOR LACTOSE PERMEASE (P02920)       \n"
    "DBREF  XXXX A    1   417  UNP    P02920   LACY_ECOLI       1    417             \n"
    "ATOM      1  N   MET A   1      -9.000   1.000   2.000  1.00 57.62           N  \n"
    "END\n"
)
AF_API = json.dumps([{"modelEntityId": "AF-P02920-F1", "latestVersion": 6}])


@pytest.fixture
def fake_alphafold(monkeypatch):
    calls = []

    def download(url, accession, *, timeout, max_bytes):
        calls.append(url)
        if "/api/prediction/" in url:
            return AF_API.encode()
        if url.endswith("AF-P02920-F1-model_v6.pdb"):
            return AF_PDB.encode()
        raise RuntimeError(f"unexpected {url}")

    monkeypatch.setattr(pdbmod, "_download", download)
    return calls


@pytest.mark.parametrize("text, expected", [
    ("AF-P02920-F1", ("AF-P02920-F1", None, None)),
    ("af-p02920-f1-model_v6", ("AF-P02920-F1", 6, None)),
    ("AF-P02920-F1-model_v6.pdb", ("AF-P02920-F1", 6, "pdb")),
    ("/some/dir/AF-A0A024RBG1-F2-model_v4.cif", ("AF-A0A024RBG1-F2", 4, "cif")),
])
def test_alphafold_identifiers_are_recognised(text, expected):
    model = pdbmod.parse_alphafold_id(text)
    assert (model.entry, model.version, model.fmt) == expected
    assert pdbmod.is_alphafold_id(text)


@pytest.mark.parametrize("text", ["1GRM", "pdb_00001grm", "AF-XXXXXX-F1", "AF-P02920", "AF-P02920-F0"])
def test_non_alphafold_identifiers_are_rejected(text):
    assert pdbmod.parse_alphafold_id(text) is None


def test_alphafold_fetch_resolve_and_cache(tmp_path, fake_alphafold):
    record = pdbmod.fetch_alphafold_model("AF-P02920-F1", tmp_path, fmt="pdb")
    assert record.path.name == "AF-P02920-F1-model_v6.pdb"
    assert record.identity_verified and "P02920" in record.identity_note
    assert pdbmod.read_provenance(record.path).sha256 == record.sha256
    calls = len(fake_alphafold)
    source = pdbmod.resolve_structure("AF-P02920-F1", cache_dir=tmp_path, fmt="pdb", offline=True)
    assert source.kind == "cached_accession" and source.path == record.path
    assert pdbmod.fetch_structure("af-p02920-f1", tmp_path, fmt="pdb").path == record.path
    assert len(fake_alphafold) == calls  # served from the cache
    frame, report = load_structure_report(source.path)
    assert report["predictor"] == "AlphaFold" and "pLDDT" in report["bfactor_meaning"]
    assert frame.atoms[0].bfactor == 57.62


def test_alphafold_mismatch_and_assembly_are_errors(tmp_path, fake_alphafold, monkeypatch):
    with pytest.raises(ValueError, match="assemblies"):
        pdbmod.resolve_structure("AF-P02920-F1", cache_dir=tmp_path, assembly=1)
    with pytest.raises(RuntimeError, match="not in the cache"):
        pdbmod.resolve_structure("AF-P02920-F1", cache_dir=tmp_path, offline=True)
    monkeypatch.setattr(pdbmod, "_download", lambda url, *a, **k: AF_PDB.replace("P02920", "P00001").encode())
    with pytest.raises(RuntimeError, match="P00001"):
        pdbmod.fetch_alphafold_model("AF-P02920-F1-model_v6", tmp_path, fmt="pdb")
    assert not (tmp_path / "AF-P02920-F1-model_v6.pdb").exists()


def test_predicted_model_report_ignores_experimental_headers():
    assert predicted_model_report("HEADER    MEMBRANE PROTEIN\nTITLE     GRAMICIDIN A\n") == {}
    assert predicted_model_report("data_AF-P02920-F1\n_entry.id AF-P02920-F1\n")["model_type"] == "predicted"


def test_extended_accession_matches_the_classic_id_stated_in_the_file():
    verified, note = pdbmod.check_identity("data_1GRM\n_entry.id   1GRM\n", "pdb_00001grm", "cif",
                                           "test", strict=True)
    assert verified and "pdb_00001grm" in note
    with pytest.raises(RuntimeError, match="not the requested"):
        pdbmod.check_identity("_entry.id 2CHB\n", "pdb_00001grm", "cif", "test", strict=True)

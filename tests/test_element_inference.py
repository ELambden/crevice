"""Element inference and radius provenance contracts (synthetic records only).

An explicit element that is a real element symbol must be kept as written, so a
mercury atom (``HG``) is never truncated to hydrogen and then dropped by the
hydrogen filter. Name-based inference is only a fallback for files without an
element column, and must not turn carbon atoms into calcium or sodium.
"""

from __future__ import annotations

import warnings

import pytest

from crevice.models import Atom
from crevice.parser import load_structure_report
from crevice.radii import (
    ELEMENT_SYMBOLS,
    VDW_RADII,
    default_radius_elements,
    infer_element,
    is_element_symbol,
    radius_source,
    vdw_radius,
)


def pdb_line(serial, name, resname, element, *, record="HETATM", x=0.0):
    return (f"{record:<6}{serial:5d} {name:<4} {resname:>3} A{serial:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          {element:>2}")


# --- explicit element column ------------------------------------------------

@pytest.mark.parametrize("symbol", ["HG", "Hg", "hg", " HG", "CD", "PT", "GD", "YB", "AU"])
def test_explicit_two_letter_elements_are_kept(symbol):
    assert infer_element("X1", symbol) == symbol.strip().upper()


def test_mercury_is_not_hydrogen_and_survives_the_hydrogen_filter(tmp_path):
    path = tmp_path / "hg.pdb"
    path.write_text("\n".join([
        pdb_line(1, "CA", "ALA", "C", record="ATOM"),
        pdb_line(2, "HG", "HG", "HG", x=5.0),
        "END",
    ]) + "\n")
    frame, _ = load_structure_report(path)
    mercury = [a for a in frame.atoms if a.resname == "HG"]
    assert [a.element for a in mercury] == ["HG"]
    kept = frame.selected_atoms(include_hydrogen=False, include_hetero=True)
    assert "HG" in [a.element for a in kept]
    assert vdw_radius("HG") == 1.55


def test_explicit_element_wins_over_the_atom_name():
    # A calcium ion named CA in residue CA, and a C-alpha whose column says C.
    assert infer_element("CA", "CA", "CA") == "CA"
    assert infer_element("CA", "C", "ALA") == "C"
    assert infer_element("CA", "CA", "ALA") == "CA"  # the file says calcium


def test_invalid_explicit_symbol_falls_back_to_the_name():
    assert not is_element_symbol("XX")
    assert infer_element("OG1", "XX", "SER") == "O"
    assert infer_element("OG1", "", "SER") == "O"
    assert infer_element("OG1", ".", "SER") == "O"


def test_charge_suffix_in_element_field_is_ignored():
    assert infer_element("FE", "FE2+", "HEM") == "FE"
    assert infer_element("O1", "O1-", "SO4") == "O"


# --- name-based fallback ----------------------------------------------------

def test_c_alpha_is_never_calcium_without_an_element_column():
    for resname in ("ALA", "GLY", "DLE", "MSE", "UNK", "HEM"):
        assert infer_element("CA", "", resname) == "C"
    assert infer_element("CA", "", "CA") == "CA"


def test_ligand_carbons_and_nitrogens_are_not_metal_ions():
    assert infer_element("CAB", "", "HEM") == "C"
    assert infer_element("CA1", "", "LIG") == "C"
    assert infer_element("NA", "", "HEM") == "N"
    assert infer_element("NA", "", "NA") == "NA"
    assert infer_element("CO", "", "CO") == "CO"
    assert infer_element("CO1", "", "LIG") == "C"


def test_unambiguous_two_letter_names_still_resolve():
    assert infer_element("CL1", "", "LIG") == "CL"
    assert infer_element("FE", "", "HEM") == "FE"
    assert infer_element("ZN", "", "ZN") == "ZN"
    assert infer_element("SE", "", "MSE") == "SE"
    assert infer_element("MG", "", "MG") == "MG"
    assert infer_element("CLA", "", "CLA") == "CL"


# --- symbol list and radii --------------------------------------------------

def test_periodic_table_is_complete_and_matches_gemmi_when_available():
    assert len(ELEMENT_SYMBOLS) == 118
    assert {"H", "HE", "HG", "OG", "TS"} <= ELEMENT_SYMBOLS
    gemmi = pytest.importorskip("gemmi")
    assert ELEMENT_SYMBOLS == {gemmi.Element(n).name.upper() for n in range(1, 119)}


def test_every_tabulated_radius_is_an_element_and_has_a_source():
    for element, radius in VDW_RADII.items():
        assert element in ELEMENT_SYMBOLS
        assert 1.0 <= radius <= 3.5
        assert radius_source(element) != "default"


@pytest.mark.parametrize("element, radius", [
    ("HG", 1.55), ("CD", 1.58), ("AG", 1.72), ("AU", 1.66), ("PT", 1.72), ("PD", 1.63),
    ("PB", 2.02), ("TL", 1.96), ("SI", 2.10), ("AS", 1.85), ("RB", 3.03), ("CS", 3.43),
    ("SR", 2.49), ("BA", 2.68), ("AL", 1.84),
])
def test_heavy_and_metal_ion_radii(element, radius):
    assert vdw_radius(element) == radius


def test_existing_radii_are_unchanged():
    expected = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
                "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90, "B": 1.92, "FE": 1.94,
                "MG": 1.73, "MN": 1.97, "ZN": 1.39, "CA": 2.31, "NA": 2.27, "K": 2.75,
                "LI": 1.82, "CU": 1.40, "NI": 1.63, "CO": 1.67}
    assert {key: VDW_RADII[key] for key in expected} == expected


def test_element_without_radius_is_excluded_and_reported(tmp_path):
    # The element table has no Gd; loaded structures now leave such atoms out
    # (with a warning) instead of giving them the 1.70 A default.
    from crevice.radii import UnrecognisedElementWarning
    assert vdw_radius("GD") == 1.70
    assert radius_source("GD") == "default"
    atoms = [Atom(1, "GD", "GD", "A", 1, 0.0, 0.0, 0.0, "GD"),
             Atom(2, "CA", "ALA", "A", 2, 5.0, 0.0, 0.0, "C")]
    assert default_radius_elements(atoms) == ["GD"]
    path = tmp_path / "gd.pdb"
    path.write_text(pdb_line(1, "GD", "GD", "GD") + "\n"
                    + pdb_line(2, "CA", "ALA", "C", record="ATOM", x=5.0) + "\nEND\n")
    with pytest.warns(UnrecognisedElementWarning, match=r"GD \(1 atom\).*--radii"):
        frame, report = load_structure_report(path)
    assert [a.element for a in frame.atoms] == ["C"]
    assert report["unsupported_radius_elements"] == ["GD"]
    assert report["excluded_no_radius_atoms"] == {"GD": 1}
    assert report["excluded_no_radius_names"] == {"GD:GD": 1}
    assert report["atom_count"] == 1


# --- force-field ion names --------------------------------------------------

@pytest.mark.parametrize("resname, name, element", [
    # CHARMM36 toppar_water_ions.str
    ("SOD", "SOD", "NA"), ("POT", "POT", "K"), ("CLA", "CLA", "CL"), ("CAL", "CAL", "CA"),
    ("MG", "MG", "MG"), ("LIT", "LIT", "LI"), ("CES", "CES", "CS"), ("ZN2", "ZN", "ZN"),
    ("BAR", "BAR", "BA"), ("RUB", "RUB", "RB"), ("CD2", "CD", "CD"),
    # GROMACS ions.itp / AMBER
    ("NA", "NA", "NA"), ("K", "K", "K"), ("CL", "CL", "CL"), ("CA", "CA", "CA"),
    ("Na+", "Na+", "NA"), ("K+", "K+", "K"), ("Cl-", "Cl-", "CL"), ("RB", "RB", "RB"),
    ("CS", "CS", "CS"), ("CU1", "CU", "CU"), ("Li+", "Li+", "LI"),
])
def test_force_field_ion_names_map_to_their_element(resname, name, element):
    assert infer_element(name, "", resname) == element


def test_ion_rule_needs_both_names_and_never_touches_protein_atoms():
    # C-alpha, proline CD, heme NA and a ligand atom called SOD stay as before.
    assert infer_element("CA", "", "ALA") == "C"
    assert infer_element("CD", "", "PRO") == "C"
    assert infer_element("NA", "", "HEM") == "N"
    assert infer_element("SOD", "", "LIG") == "S"
    assert infer_element("POT", "", "LIG") == "P"


def test_explicit_element_column_beats_the_ion_table():
    assert infer_element("SOD", "S", "SOD") == "S"
    assert infer_element("CA", "C", "CA") == "C"


def test_every_ion_entry_is_a_real_element_with_a_cited_name_source():
    from crevice.radii import ION_NAME_ELEMENTS, ion_name_source
    for (resname, name), element in ION_NAME_ELEMENTS.items():
        assert element in ELEMENT_SYMBOLS
        assert ion_name_source(resname, name)
    assert ion_name_source("SOD", "SOD") == "CHARMM36 toppar_water_ions.str"
    assert ion_name_source("ALA", "CA") is None

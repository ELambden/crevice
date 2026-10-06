"""Synthetic mmCIF intake contracts for the optional Gemmi-backed loader.

These are parser and policy contracts on constructed text. They do not validate
biological assembly choice, cavity assignment, or any real deposited structure.
"""

from __future__ import annotations

import builtins

import pytest

pytest.importorskip("gemmi", reason="mmCIF intake contracts require the structures extra")

from crevice import mmcif


ATOM_COLUMNS = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id",
    "label_asym_id", "label_entity_id", "label_seq_id", "auth_seq_id",
    "auth_comp_id", "auth_atom_id", "auth_asym_id", "label_alt_id",
    "pdbx_PDB_ins_code", "occupancy", "B_iso_or_equiv",
    "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num",
)

DEFAULT_ROW = {
    "group_PDB": "ATOM", "type_symbol": "C", "label_atom_id": "CA",
    "label_comp_id": "ALA", "label_asym_id": "A", "label_entity_id": "1",
    "label_seq_id": "1", "auth_seq_id": "1", "auth_comp_id": ".",
    "auth_atom_id": ".", "auth_asym_id": "A", "label_alt_id": ".",
    "pdbx_PDB_ins_code": ".", "occupancy": "1.00", "B_iso_or_equiv": "10.0",
    "Cartn_x": "0.0", "Cartn_y": "0.0", "Cartn_z": "0.0",
    "pdbx_PDB_model_num": "1",
}

DEFAULT_ENTITIES = (("1", "polymer"), ("2", "non-polymer"))


def row(serial, **overrides):
    """One atom-site row; unset columns take a valid default."""
    data = dict(DEFAULT_ROW, id=str(serial))
    data.update({key: str(value) for key, value in overrides.items()})
    return data


def mmcif_text(rows, *, entities=DEFAULT_ENTITIES, entry_id="TEST",
               columns=ATOM_COLUMNS, extra="", blocks=1):
    lines = [f"data_{entry_id}", f"_entry.id {entry_id}"]
    if entities is not None:
        lines += ["loop_", "_entity.id", "_entity.type"]
        lines += [f"{eid} {etype}" for eid, etype in entities]
        lines.append("#")
    lines.append("loop_")
    lines += [f"_atom_site.{column}" for column in columns]
    for item in rows:
        lines.append(" ".join(_token(item[column]) for column in columns))
    lines.append("#")
    if extra:
        lines.append(extra)
    text = "\n".join(lines) + "\n"
    if blocks > 1:
        text += text.replace(f"data_{entry_id}", f"data_{entry_id}_2", 1)
    return text


def _token(value: str) -> str:
    return value if value and " " not in value else (value or ".")


def write(tmp_path, rows, **kwargs):
    path = tmp_path / "structure.cif"
    path.write_text(mmcif_text(rows, **kwargs))
    return path


def protein_rows():
    """Two polymer residues plus a non-polymer water."""
    return [
        row(1, label_atom_id="N", type_symbol="N", label_seq_id=1, auth_seq_id=1),
        row(2, label_atom_id="CA", type_symbol="C", label_seq_id=1, auth_seq_id=1),
        row(3, label_atom_id="CB", type_symbol="C", label_seq_id=2, auth_seq_id=2,
            label_comp_id="SER", Cartn_x=3.0),
        row(4, group_PDB="HETATM", label_atom_id="O", type_symbol="O",
            label_comp_id="HOH", label_asym_id="B", label_entity_id=2,
            label_seq_id=".", auth_seq_id=50, Cartn_x=20.0),
    ]


# --- selection policy -------------------------------------------------------

def test_polymer_heavy_keeps_polymer_and_drops_nonpolymer(tmp_path):
    frame, report = mmcif.load_mmcif(write(tmp_path, protein_rows()))
    assert [atom.name for atom in frame.atoms] == ["N", "CA", "CB"]
    assert report["excluded_atoms"]["nonpolymer"] == 1
    assert report["polymer_chain_count"] == 1


def test_polymer_hetatm_is_retained_by_entity_membership(tmp_path):
    """A modified polymer residue deposited as HETATM stays in the polymer."""
    rows = protein_rows() + [
        row(5, group_PDB="HETATM", label_atom_id="CD", type_symbol="C",
            label_comp_id="MSE", label_entity_id=1, label_seq_id=3,
            auth_seq_id=3, Cartn_x=6.0),
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows))
    assert "CD" in [atom.name for atom in frame.atoms]
    assert report["polymer_hetatm_retained"] == 1


def test_all_heavy_keeps_nonpolymer_but_not_hydrogen(tmp_path):
    rows = protein_rows() + [
        row(5, label_atom_id="HA", type_symbol="H", label_seq_id=2, auth_seq_id=2),
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows), selection="all-heavy")
    assert len(frame.atoms) == 4
    assert report["excluded_atoms"]["hydrogen_or_deuterium"] == 1
    assert "nonpolymer" not in report["excluded_atoms"]


def test_all_selection_retains_hydrogen(tmp_path):
    rows = protein_rows() + [
        row(5, label_atom_id="HA", type_symbol="H", label_seq_id=2, auth_seq_id=2),
    ]
    frame, _ = mmcif.load_mmcif(write(tmp_path, rows), selection="all")
    assert len(frame.atoms) == 5


def test_deuterium_is_excluded_with_hydrogen(tmp_path):
    rows = protein_rows() + [
        row(5, label_atom_id="DA", type_symbol="D", label_seq_id=2, auth_seq_id=2),
    ]
    _, report = mmcif.load_mmcif(write(tmp_path, rows), selection="all-heavy")
    assert report["excluded_atoms"]["hydrogen_or_deuterium"] == 1


def test_invalid_selection_and_namespace_are_rejected(tmp_path):
    path = write(tmp_path, protein_rows())
    with pytest.raises(ValueError, match="selection must be"):
        mmcif.load_mmcif(path, selection="backbone")
    with pytest.raises(ValueError, match="chain_ids must be"):
        mmcif.load_mmcif(path, chain_ids="entity")


def test_empty_selection_reports_the_selection_not_a_frame_error(tmp_path):
    """Excluding every atom is a selection outcome, not a malformed frame."""
    rows = [row(1, group_PDB="HETATM", label_comp_id="HOH", label_entity_id=2,
                label_asym_id="B", type_symbol="O", label_seq_id=".", auth_seq_id=1)]
    with pytest.raises(ValueError, match="selection retained no atoms"):
        mmcif.load_mmcif(write(tmp_path, rows))


# --- model choice -----------------------------------------------------------

def test_model_index_selects_the_requested_model(tmp_path):
    rows = [
        row(1, Cartn_x=0.0, pdbx_PDB_model_num=1),
        row(2, Cartn_x=1.0, pdbx_PDB_model_num=1, label_atom_id="CB"),
        row(3, Cartn_x=7.0, pdbx_PDB_model_num=2),
        row(4, Cartn_x=8.0, pdbx_PDB_model_num=2, label_atom_id="CB"),
    ]
    path = write(tmp_path, rows)
    first, first_report = mmcif.load_mmcif(path, model_index=0)
    second, second_report = mmcif.load_mmcif(path, model_index=1)
    assert [atom.x for atom in first.atoms] == [0.0, 1.0]
    assert [atom.x for atom in second.atoms] == [7.0, 8.0]
    assert first_report["model_ids"] == ["1", "2"]
    assert second_report["selected_model_id"] == "2"
    assert second.model_index == 1


def test_model_index_out_of_range_is_rejected(tmp_path):
    path = write(tmp_path, protein_rows())
    with pytest.raises(ValueError, match="out of range"):
        mmcif.load_mmcif(path, model_index=3)


@pytest.mark.parametrize("bad", [-1, True, 1.5, "0"])
def test_model_index_rejects_non_index_values(tmp_path, bad):
    path = write(tmp_path, protein_rows())
    with pytest.raises(ValueError, match="model_index must be"):
        mmcif.load_mmcif(path, model_index=bad)


def test_absent_model_column_is_a_single_model(tmp_path):
    columns = tuple(c for c in ATOM_COLUMNS if c != "pdbx_PDB_model_num")
    path = write(tmp_path, protein_rows(), columns=columns)
    _, report = mmcif.load_mmcif(path)
    assert report["model_ids"] == ["1"]


# --- chain namespaces and copy identity -------------------------------------

def test_label_chains_separate_copies_sharing_one_auth_chain(tmp_path):
    """Assembly copies share an auth chain; label IDs must keep them distinct."""
    rows = [
        row(1, label_asym_id="A", auth_asym_id="A"),
        row(2, label_asym_id="B", auth_asym_id="A", Cartn_x=30.0),
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows), chain_ids="label")
    assert sorted(report["chains"]) == ["A", "B"]
    assert len(frame.atoms) == 2
    assert report["chains"]["B"]["auth_asym_ids"] == ["A"]


def test_auth_namespace_rejects_colliding_copies(tmp_path):
    rows = [
        row(1, label_asym_id="A", auth_asym_id="A"),
        row(2, label_asym_id="B", auth_asym_id="A", Cartn_x=30.0),
    ]
    with pytest.raises(ValueError, match="Duplicate selected atom identity"):
        mmcif.load_mmcif(write(tmp_path, rows), chain_ids="auth")


def test_missing_chain_id_in_requested_namespace_is_rejected(tmp_path):
    rows = [row(1, auth_asym_id=".")]
    with pytest.raises(ValueError, match="Missing auth chain ID"):
        mmcif.load_mmcif(write(tmp_path, rows), chain_ids="auth")


def test_insertion_codes_keep_residues_distinct(tmp_path):
    rows = [
        row(1, auth_seq_id=52, pdbx_PDB_ins_code="."),
        row(2, auth_seq_id=52, pdbx_PDB_ins_code="A", Cartn_x=4.0),
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows))
    assert len(frame.residues()) == 2
    assert report["residue_count"] == 2


def test_auth_seq_id_is_preferred_with_label_fallback(tmp_path):
    rows = [
        row(1, label_seq_id=1, auth_seq_id=101),
        row(2, label_seq_id=2, auth_seq_id=".", label_atom_id="CB", Cartn_x=4.0),
    ]
    frame, _ = mmcif.load_mmcif(write(tmp_path, rows))
    assert [atom.resid for atom in frame.atoms] == [101, 2]


# --- alternate locations ----------------------------------------------------

def test_first_lexical_altloc_is_kept_with_blank_atoms(tmp_path):
    rows = [
        row(1, label_atom_id="N", type_symbol="N", label_alt_id="."),
        row(2, label_atom_id="CA", label_alt_id="B", Cartn_x=2.0),
        row(3, label_atom_id="CA", label_alt_id="A", Cartn_x=1.0),
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows))
    assert [(atom.name, atom.altloc) for atom in frame.atoms] == [("N", ""), ("CA", "A")]
    assert report["altloc_residue_count"] == 1
    assert report["excluded_atoms"]["alternate_location"] == 1


def test_blank_and_labelled_altloc_on_one_atom_is_rejected(tmp_path):
    """A file mixing conventions for one atom has no deterministic conformer."""
    rows = [
        row(1, label_atom_id="CA", label_alt_id="."),
        row(2, label_atom_id="CA", label_alt_id="B", Cartn_x=2.0),
    ]
    with pytest.raises(ValueError, match="blank and labelled alternate"):
        mmcif.load_mmcif(write(tmp_path, rows))


def test_altloc_choice_is_per_residue(tmp_path):
    rows = [
        row(1, auth_seq_id=1, label_alt_id="A", Cartn_x=1.0),
        row(2, auth_seq_id=1, label_alt_id="B", Cartn_x=2.0),
        row(3, auth_seq_id=2, label_alt_id="B", Cartn_x=3.0),
        row(4, auth_seq_id=2, label_alt_id="C", Cartn_x=4.0),
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows))
    assert [atom.x for atom in frame.atoms] == [1.0, 3.0]
    assert report["altloc_residue_count"] == 2


# --- occupancy --------------------------------------------------------------

def test_zero_occupancy_atoms_are_excluded(tmp_path):
    rows = protein_rows() + [
        row(5, label_atom_id="CG", label_seq_id=2, auth_seq_id=2, occupancy="0.00")
    ]
    frame, report = mmcif.load_mmcif(write(tmp_path, rows))
    assert "CG" not in [atom.name for atom in frame.atoms]
    assert report["excluded_atoms"]["zero_occupancy"] == 1


def test_negative_occupancy_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="Negative atom occupancy"):
        mmcif.load_mmcif(write(tmp_path, [row(1, occupancy="-0.5")]))


def test_absent_occupancy_column_is_accepted(tmp_path):
    columns = tuple(c for c in ATOM_COLUMNS if c != "occupancy")
    frame, _ = mmcif.load_mmcif(write(tmp_path, protein_rows(), columns=columns))
    assert all(atom.occupancy is None for atom in frame.atoms)


# --- malformed input --------------------------------------------------------

def test_missing_required_columns_are_rejected(tmp_path):
    columns = tuple(c for c in ATOM_COLUMNS if c != "Cartn_z")
    with pytest.raises(ValueError, match="complete atom-site coordinate table"):
        mmcif.load_mmcif(write(tmp_path, protein_rows(), columns=columns))


def test_absent_atom_site_category_is_rejected(tmp_path):
    path = tmp_path / "empty.cif"
    path.write_text("data_TEST\n_entry.id TEST\n")
    with pytest.raises(ValueError, match="complete atom-site coordinate table"):
        mmcif.load_mmcif(path)


def test_multi_block_input_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="single-block"):
        mmcif.load_mmcif(write(tmp_path, protein_rows(), blocks=2))


def test_unparsable_text_is_rejected(tmp_path):
    path = tmp_path / "broken.cif"
    path.write_text("data_TEST\nloop_\n_atom_site.id\n'unterminated\n")
    with pytest.raises(ValueError):
        mmcif.load_mmcif(path)


def test_duplicate_atom_site_ids_are_rejected(tmp_path):
    rows = [row(1), row(1, label_atom_id="CB", Cartn_x=2.0)]
    with pytest.raises(ValueError, match="duplicate mmCIF atom-site ID"):
        mmcif.load_mmcif(write(tmp_path, rows))


def test_missing_atom_site_id_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="duplicate mmCIF atom-site ID"):
        mmcif.load_mmcif(write(tmp_path, [row(".")]))


def test_noninteger_atom_site_id_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="atom-site ID"):
        mmcif.load_mmcif(write(tmp_path, [row("1a")]))


def test_unsupported_group_pdb_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="group_PDB"):
        mmcif.load_mmcif(write(tmp_path, [row(1, group_PDB="ANISOU")]))


@pytest.mark.parametrize("missing", [".", "?"])
def test_missing_coordinate_is_rejected_not_defaulted(tmp_path, missing):
    """An absent coordinate must fail, never silently become the origin."""
    rows = [row(1), row(2, label_atom_id="CB", Cartn_y=missing)]
    with pytest.raises(ValueError, match="[Cc]oordinate"):
        mmcif.load_mmcif(write(tmp_path, rows))


def test_nonnumeric_coordinate_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="[Cc]oordinate"):
        mmcif.load_mmcif(write(tmp_path, [row(1, Cartn_x="abc")]))


def test_nonfinite_coordinate_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="finite"):
        mmcif.load_mmcif(write(tmp_path, [row(1, Cartn_x="NaN")]))


def test_missing_residue_number_is_rejected(tmp_path):
    rows = [row(1, label_seq_id=".", auth_seq_id=".")]
    with pytest.raises(ValueError, match="residue number"):
        mmcif.load_mmcif(write(tmp_path, rows))


def test_noninteger_residue_number_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="residue number"):
        mmcif.load_mmcif(write(tmp_path, [row(1, auth_seq_id="12A")]))


def test_missing_element_symbol_is_rejected(tmp_path):
    """Element drives the van der Waals radius; it cannot default to blank."""
    with pytest.raises(ValueError, match="element"):
        mmcif.load_mmcif(write(tmp_path, [row(1, type_symbol=".")]))


@pytest.mark.parametrize("column", ["label_atom_id", "label_comp_id"])
def test_blank_required_label_value_is_rejected(tmp_path, column):
    with pytest.raises(ValueError, match=f"Missing required mmCIF atom-site value {column}"):
        mmcif.load_mmcif(write(tmp_path, [row(1, **{column: "."})]))


def test_nonnumeric_occupancy_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="occupancy"):
        mmcif.load_mmcif(write(tmp_path, [row(1, occupancy="high")]))


# --- entity information -----------------------------------------------------

def test_polymer_selection_requires_entity_types(tmp_path):
    with pytest.raises(ValueError, match="entity type"):
        mmcif.load_mmcif(write(tmp_path, protein_rows(), entities=None))


def test_absent_entity_table_does_not_claim_zero_polymers(tmp_path):
    """Without _entity, polymer counts are unknown rather than zero."""
    _, report = mmcif.load_mmcif(
        write(tmp_path, protein_rows(), entities=None), selection="all-heavy"
    )
    assert report["entity_types_available"] is False
    assert report["polymer_chain_count"] is None
    assert report["polymer_hetatm_retained"] is None


def test_entity_table_present_reports_polymer_counts(tmp_path):
    _, report = mmcif.load_mmcif(write(tmp_path, protein_rows()), selection="all-heavy")
    assert report["entity_types_available"] is True
    assert report["polymer_chain_count"] == 1


# --- intake report ----------------------------------------------------------

def test_report_records_parser_provenance_and_policy(tmp_path):
    import gemmi

    _, report = mmcif.load_mmcif(write(tmp_path, protein_rows()))
    assert report["parser"] == "gemmi"
    assert report["parser_version"] == gemmi.__version__
    assert report["entry_id"] == "TEST"
    assert report["selection"] == "polymer-heavy"
    assert report["chain_id_namespace"] == "label"
    assert report["coordinate_units"] == "Angstrom"
    assert report["completeness_validated"] is False
    assert "first lexical" in report["altloc_policy"]


def test_report_atom_accounting_is_exact(tmp_path):
    rows = protein_rows() + [
        row(5, label_atom_id="HA", type_symbol="H", label_seq_id=2, auth_seq_id=2),
        row(6, label_atom_id="CG", label_seq_id=2, auth_seq_id=2, occupancy="0.00"),
        row(7, label_atom_id="OG", type_symbol="O", label_seq_id=2, auth_seq_id=2,
            label_comp_id="SER", label_alt_id="A", Cartn_y=1.0),
        row(8, label_atom_id="OG", type_symbol="O", label_seq_id=2, auth_seq_id=2,
            label_comp_id="SER", label_alt_id="B", Cartn_y=2.0),
    ]
    _, report = mmcif.load_mmcif(write(tmp_path, rows))
    assert report["input_model_atom_count"] == 8
    assert report["atom_count"] == 4
    assert report["excluded_atoms"] == {
        "alternate_location": 1, "zero_occupancy": 1,
        "nonpolymer": 1, "hydrogen_or_deuterium": 1,
    }
    assert report["atom_count"] + sum(report["excluded_atoms"].values()) == 8


def test_report_counts_only_the_selected_model(tmp_path):
    rows = [row(1), row(2, pdbx_PDB_model_num=2), row(3, pdbx_PDB_model_num=2,
                                                     label_atom_id="CB", Cartn_x=2.0)]
    _, report = mmcif.load_mmcif(write(tmp_path, rows), model_index=1)
    assert report["input_model_atom_count"] == 2
    assert report["atom_count"] == 2


def test_report_flags_elements_without_radii(tmp_path):
    rows = [row(1), row(2, label_atom_id="XX", type_symbol="XX", Cartn_x=4.0)]
    _, report = mmcif.load_mmcif(write(tmp_path, rows), selection="all")
    assert report["unsupported_radius_elements"] == ["XX"]


def test_report_records_bounds_and_assembly_definitions(tmp_path):
    extra = "\n".join([
        "loop_", "_pdbx_struct_assembly.id", "_pdbx_struct_assembly.details",
        "1 author_defined_assembly", "#",
    ])
    rows = [row(1, Cartn_x=-2.0), row(2, label_atom_id="CB", Cartn_x=5.0, Cartn_z=1.0)]
    _, report = mmcif.load_mmcif(write(tmp_path, rows, extra=extra))
    assert report["bounds_A"] == ((-2.0, 0.0, 0.0), (5.0, 0.0, 1.0))
    assert [item["id"] for item in report["assembly_definitions"]] == ["1"]


def test_report_counts_unobserved_records(tmp_path):
    extra = "\n".join([
        "loop_", "_pdbx_unobs_or_zero_occ_residues.id",
        "_pdbx_unobs_or_zero_occ_residues.auth_seq_id", "1 41", "2 42", "#",
    ])
    _, report = mmcif.load_mmcif(write(tmp_path, protein_rows(), extra=extra))
    assert report["unobserved_residue_records"] == 2
    assert report["unobserved_atom_records"] == 0


def test_report_is_json_serializable(tmp_path):
    import json

    _, report = mmcif.load_mmcif(write(tmp_path, protein_rows()))
    assert json.loads(json.dumps(report))["atom_count"] == 3


def test_entry_id_placeholder_is_not_treated_as_an_accession(tmp_path):
    """RCSB-expanded assemblies carry a placeholder entry ID."""
    path = write(tmp_path, protein_rows(), entry_id="XXXX")
    _, report = mmcif.load_mmcif(path)
    assert report["entry_id"] == "XXXX"
    assert report["entry_id_is_placeholder"] is True


# --- optional dependency ----------------------------------------------------

def test_absent_gemmi_reports_the_install_extra(tmp_path, monkeypatch):
    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "gemmi":
            raise ImportError("no gemmi")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(RuntimeError, match=r"crevice\[structures\]"):
        mmcif.load_mmcif(write(tmp_path, protein_rows()))

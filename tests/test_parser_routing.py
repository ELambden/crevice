"""Which parser handles which input, and what changes when it does.

These are routing and equivalence contracts on constructed files. They do not
validate either parser against a real deposited structure.
"""

from __future__ import annotations

import json

import pytest

from crevice import mmcif
from crevice.parser import load_structure, load_structure_report

from test_mmcif_intake import protein_rows, row, write


PDB_TEXT = """\
HEADER    TRANSPORT PROTEIN                       01-JAN-00   1GRM
ATOM      1  N   VAL A   1      -1.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  VAL A   1       0.000   0.000   0.000  1.00 10.00           C
HETATM    3  O   HOH A 100       9.000   9.000   9.000  1.00 10.00           O
END
"""


def pdb_file(tmp_path):
    path = tmp_path / "structure.pdb"
    path.write_text(PDB_TEXT)
    return path


# --- routing ----------------------------------------------------------------

def test_mmcif_uses_the_validated_parser_when_available(tmp_path):
    pytest.importorskip("gemmi")
    _, report = load_structure_report(write(tmp_path, protein_rows()))
    assert report["parser"] == "gemmi"
    assert report["input_format"] == "mmcif"


def test_pdb_always_uses_the_builtin_parser(tmp_path):
    _, report = load_structure_report(pdb_file(tmp_path))
    assert report["parser"] == "builtin"
    assert report["input_format"] == "pdb"


def test_builtin_parser_can_be_requested_for_mmcif(tmp_path):
    frame, report = load_structure_report(write(tmp_path, protein_rows()), parser="builtin")
    assert report["parser"] == "builtin"
    assert report["input_format"] == "mmcif"
    assert len(frame.atoms) > 0


def test_auto_falls_back_to_builtin_without_the_structures_extra(tmp_path, monkeypatch):
    monkeypatch.setattr(mmcif, "gemmi_available", lambda: False)
    _, report = load_structure_report(write(tmp_path, protein_rows()))
    assert report["parser"] == "builtin"


def test_requesting_gemmi_without_the_extra_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setattr(mmcif, "gemmi_available", lambda: False)
    with pytest.raises(RuntimeError, match=r"crevice\[structures\]"):
        load_structure_report(write(tmp_path, protein_rows()), parser="gemmi")


def test_requesting_gemmi_for_a_pdb_file_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="reads mmCIF"):
        load_structure_report(pdb_file(tmp_path), parser="gemmi")


def test_unknown_parser_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="parser must be one of"):
        load_structure_report(pdb_file(tmp_path), parser="biopython")


def test_negative_model_index_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="model_index must be non-negative"):
        load_structure_report(pdb_file(tmp_path), model_index=-1)


# --- what routing preserves -------------------------------------------------

def test_both_parsers_load_the_same_atoms_for_a_plain_mmcif(tmp_path):
    """Routing must not change coordinates or atom count for ordinary input."""
    pytest.importorskip("gemmi")
    path = write(tmp_path, protein_rows())
    validated, _ = load_structure_report(path)
    builtin, _ = load_structure_report(path, parser="builtin")
    assert len(validated.atoms) == len(builtin.atoms)
    assert [a.coord for a in validated.atoms] == [a.coord for a in builtin.atoms]
    assert [a.name for a in validated.atoms] == [a.name for a in builtin.atoms]


def test_load_structure_retains_hydrogen_and_hetero_for_downstream_filtering(tmp_path):
    """Filtering stays with selected_atoms, so nothing is dropped at load time."""
    pytest.importorskip("gemmi")
    rows = protein_rows() + [
        row(5, label_atom_id="HA", type_symbol="H", label_seq_id=2, auth_seq_id=2),
    ]
    frame = load_structure(write(tmp_path, rows))
    assert len(frame.atoms) == 5
    assert any(atom.element == "H" for atom in frame.atoms)
    assert any(atom.hetero for atom in frame.atoms)


def test_modified_polymer_residue_survives_hetero_filtering(tmp_path):
    """An MSE deposited as HETATM is polymer, so --no-hetero must keep it."""
    pytest.importorskip("gemmi")
    rows = protein_rows() + [
        row(5, group_PDB="HETATM", label_atom_id="SE", type_symbol="SE",
            label_comp_id="MSE", label_entity_id=1, label_seq_id=3,
            auth_seq_id=3, Cartn_x=6.0),
    ]
    path = write(tmp_path, rows)
    frame, report = load_structure_report(path)
    assert report["hetero_definition"] == "nonpolymer"
    kept = frame.selected_atoms(include_hydrogen=False, include_hetero=False)
    assert "SE" in [atom.name for atom in kept]
    assert "O" not in [atom.name for atom in kept]

    legacy, legacy_report = load_structure_report(path, parser="builtin")
    assert legacy_report["hetero_definition"] == "group_pdb"
    legacy_kept = legacy.selected_atoms(include_hydrogen=False, include_hetero=False)
    assert "SE" not in [atom.name for atom in legacy_kept]


def test_entity_free_mmcif_falls_back_to_the_record_type(tmp_path):
    pytest.importorskip("gemmi")
    path = write(tmp_path, [row(1), row(2, group_PDB="HETATM", label_atom_id="O",
                                 type_symbol="O", label_comp_id="HOH", Cartn_x=9.0)],
                 entities=None)
    _, report = load_structure_report(path)
    assert report["hetero_definition"] == "group_pdb"


def test_model_index_is_honoured_through_routing(tmp_path):
    pytest.importorskip("gemmi")
    rows = [row(1, Cartn_x=0.0, pdbx_PDB_model_num=1),
            row(2, Cartn_x=7.0, pdbx_PDB_model_num=2)]
    path = write(tmp_path, rows)
    assert load_structure(path, model_index=1).atoms[0].x == 7.0
    assert load_structure(path, model_index=1, parser="builtin").atoms[0].x == 7.0


def test_chain_id_namespace_is_selectable(tmp_path):
    pytest.importorskip("gemmi")
    rows = [row(1, label_asym_id="A", auth_asym_id="P"),
            row(2, label_asym_id="B", auth_asym_id="Q", Cartn_x=30.0)]
    path = write(tmp_path, rows)
    assert {a.chain_id for a in load_structure(path).atoms} == {"A", "B"}
    assert {a.chain_id for a in load_structure(path, chain_ids="auth").atoms} == {"P", "Q"}


def test_reports_are_json_serializable(tmp_path):
    pytest.importorskip("gemmi")
    for path in (write(tmp_path, protein_rows()), pdb_file(tmp_path)):
        _, report = load_structure_report(path)
        assert json.loads(json.dumps(report))["completeness_validated"] is False

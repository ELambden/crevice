"""The builtin and Gemmi mmCIF readers agree where a difference was accidental.

Synthetic files only (software tests). Chain namespace and alternate-location
choice now follow one rule in both readers; the deliberate differences
(heteroatom meaning, zero-occupancy atoms) are stated in every report.
"""
from __future__ import annotations

import pytest

from crevice.parser import load_structure_report

from test_mmcif_intake import row, write


def chain_rows():
    return [row(1, label_asym_id="A", auth_asym_id="P"),
            row(2, label_asym_id="B", auth_asym_id="Q", label_seq_id=2, auth_seq_id=2, Cartn_x=30.0),
            row(3, label_asym_id="C", auth_asym_id=".", label_seq_id=3, auth_seq_id=3, Cartn_x=60.0)]


def test_builtin_mmcif_honours_the_chain_namespace():
    import tempfile, pathlib
    path = write(pathlib.Path(tempfile.mkdtemp()), chain_rows())
    frame, report = load_structure_report(path, parser="builtin")
    assert [a.chain_id for a in frame.atoms] == ["A", "B", "C"]
    assert report["chain_id_namespace"] == "label"
    frame, report = load_structure_report(path, parser="builtin", chain_ids="auth")
    # A null auth_asym_id falls back to the label ID for that row.
    assert [a.chain_id for a in frame.atoms] == ["P", "Q", "C"]
    assert report["chain_id_namespace"] == "auth"


def test_both_readers_give_the_same_chains(tmp_path):
    pytest.importorskip("gemmi")
    path = write(tmp_path, chain_rows()[:2])
    for namespace in ("label", "auth"):
        builtin, _ = load_structure_report(path, parser="builtin", chain_ids=namespace)
        gemmi, _ = load_structure_report(path, parser="gemmi", chain_ids=namespace)
        assert [a.chain_id for a in builtin.atoms] == [a.chain_id for a in gemmi.atoms]


def altloc_rows():
    return [
        # residue 1: A and B -> A
        row(1, label_atom_id="CA", label_alt_id="A", occupancy="0.6"),
        row(2, label_atom_id="CA", label_alt_id="B", occupancy="0.4", Cartn_x=0.5),
        # residue 2: only B and C deposited -> B (was dropped entirely by the builtin reader)
        row(3, label_atom_id="CA", label_alt_id="B", occupancy="0.5", label_seq_id=2, auth_seq_id=2, Cartn_x=4.0),
        row(4, label_atom_id="CA", label_alt_id="C", occupancy="0.5", label_seq_id=2, auth_seq_id=2, Cartn_x=4.5),
        # residue 3: blank atom plus numbered locations -> blank and 1
        row(5, label_atom_id="N", label_seq_id=3, auth_seq_id=3, Cartn_x=8.0),
        row(6, label_atom_id="CA", label_alt_id="1", occupancy="0.7", label_seq_id=3, auth_seq_id=3, Cartn_x=8.5),
        row(7, label_atom_id="CA", label_alt_id="2", occupancy="0.3", label_seq_id=3, auth_seq_id=3, Cartn_x=9.0),
    ]


EXPECTED_ALTLOC = [("CA", 1, "A"), ("CA", 2, "B"), ("N", 3, ""), ("CA", 3, "1")]


def test_builtin_mmcif_keeps_the_first_alternate_location_per_residue(tmp_path):
    frame, report = load_structure_report(write(tmp_path, altloc_rows()), parser="builtin")
    assert [(a.name, a.resid, a.altloc) for a in frame.atoms] == EXPECTED_ALTLOC
    assert report["altloc_policy"] == "blank plus first lexical nonblank per residue"


def test_builtin_pdb_keeps_the_first_alternate_location_per_residue(tmp_path):
    def line(serial, name, alt, resid, x):
        return (f"ATOM  {serial:5d} {name:<4}{alt:1}ALA A{resid:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}           C")
    path = tmp_path / "alt.pdb"
    path.write_text("\n".join([line(1, "CA", "A", 1, 0), line(2, "CA", "B", 1, .5),
                               line(3, "CA", "B", 2, 4), line(4, "CA", "C", 2, 4.5),
                               line(5, "N", " ", 3, 8), line(6, "CA", "1", 3, 8.5),
                               line(7, "CA", "2", 3, 9), "END"]) + "\n")
    frame, report = load_structure_report(path)
    assert [(a.name, a.resid, a.altloc) for a in frame.atoms] == EXPECTED_ALTLOC
    assert report["altloc_policy"] == "blank plus first lexical nonblank per residue"


def test_both_mmcif_readers_choose_the_same_alternate_locations(tmp_path):
    pytest.importorskip("gemmi")
    path = write(tmp_path, altloc_rows())
    builtin, _ = load_structure_report(path, parser="builtin")
    gemmi, _ = load_structure_report(path, parser="gemmi")
    assert [(a.name, a.resid, a.altloc) for a in builtin.atoms] == [(a.name, a.resid, a.altloc) for a in gemmi.atoms]


def test_deliberate_differences_are_stated_in_the_report(tmp_path):
    rows = [row(1), row(2, label_atom_id="CB", occupancy="0.0", Cartn_x=1.5)]
    frame, report = load_structure_report(write(tmp_path, rows), parser="builtin")
    assert len(frame.atoms) == 2
    assert report["zero_occupancy_policy"] == "retained"
    assert report["hetero_definition"] == "group_pdb"
    gemmi = pytest.importorskip("gemmi")
    frame, report = load_structure_report(write(tmp_path, rows), parser="gemmi")
    assert len(frame.atoms) == 1
    assert report["zero_occupancy_policy"] == "excluded"

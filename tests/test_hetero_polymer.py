"""Polymer residues deposited as HETATM count as protein (``hetero_polymer``)."""

import pytest

from crevice.parser import load_structure_report, polymer_linked_residues
from crevice.rolling import select_cast_atoms


def _line(record, serial, name, resname, chain, resid, x, y, z, element):
    return (f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain}{resid:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}")


def _residue(record, chain, resid, resname, x0, serial, y0=0.0):
    """Backbone N, CA, C, O; consecutive residues 3.33 A apart give a 1.33 A C-N bond."""
    atoms = (("N", 0.0, 0.0), ("CA", 1.0, 1.0), ("C", 2.0, 0.0), ("O", 2.5, 1.1))
    return [_line(record, serial + k, name, resname, chain, resid, x0 + dx, y0 + dy, 0.0, name[0])
            for k, (name, dx, dy) in enumerate(atoms)]


def _write(tmp_path, lines, name="model.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def _chain_with_ligands():
    """ALA1 - DLE2 (HETATM) - ALA3 in chain A, a ligand with N and C atoms, a water."""
    return (_residue("ATOM", "A", 1, "ALA", 0.0, 1) + _residue("HETATM", "A", 2, "DLE", 3.33, 5)
            + _residue("ATOM", "A", 3, "ALA", 6.66, 9)
            # Ligand in the same chain: its N sits 1.9 A from ALA3's C (no bond), its C is far away.
            + [_line("HETATM", 13, "N", "LIG", "A", 101, 10.56, 0.0, 0.0, "N"),
               _line("HETATM", 14, "C", "LIG", "A", 101, 11.9, 0.5, 0.0, "C"),
               _line("HETATM", 15, "O", "HOH", "A", 201, 20.0, 0.0, 0.0, "O")])


def _hetero_names(frame):
    return sorted({a.resname for a in frame.atoms if a.hetero})


def test_peptide_bonded_hetatm_residue_is_polymer_and_ligands_stay_hetero(tmp_path):
    path = _write(tmp_path, _chain_with_ligands())
    frame, report = load_structure_report(path)
    assert _hetero_names(frame) == ["HOH", "LIG"]
    assert report["hetero_definition"] == "group_pdb_except_polymer"
    assert report["hetero_polymer"] == "include"
    assert report["polymer_hetatm_residues"] == {"DLE": 1}
    assert report["polymer_hetatm_retained"] == 4
    cast_frame, selection = select_cast_atoms(frame, "protein")
    assert {a.resname for a in cast_frame.atoms} == {"ALA", "DLE"}
    assert selection["polymer_hetero_residues_excluded"] == {}


def test_opt_out_restores_the_record_type(tmp_path):
    path = _write(tmp_path, _chain_with_ligands())
    frame, report = load_structure_report(path, hetero_polymer="exclude")
    assert _hetero_names(frame) == ["DLE", "HOH", "LIG"]
    assert report["hetero_definition"] == "group_pdb"
    assert report["polymer_hetatm_residues"] == {}
    _, selection = select_cast_atoms(frame, "protein")
    assert selection["polymer_hetero_residues_excluded"] == {"DLE": 1}
    with pytest.raises(ValueError, match="hetero_polymer"):
        load_structure_report(path, hetero_polymer="maybe")


def test_hetatm_runs_join_through_each_other_and_at_chain_ends(tmp_path):
    # FVA1 (HETATM) - DLE2 (HETATM) - ALA3 - ETA-like cap (HETATM, only an N): all polymer.
    lines = (_residue("HETATM", "A", 1, "FVA", 0.0, 1) + _residue("HETATM", "A", 2, "DLE", 3.33, 5)
             + _residue("ATOM", "A", 3, "ALA", 6.66, 9)
             + [_line("HETATM", 13, "N", "ETA", "A", 4, 9.99, 0.0, 0.0, "N"),
                _line("HETATM", 14, "CA", "ETA", "A", 4, 10.9, 1.0, 0.0, "C")])
    frame, report = load_structure_report(_write(tmp_path, lines))
    assert _hetero_names(frame) == []
    assert report["polymer_hetatm_residues"] == {"DLE": 1, "ETA": 1, "FVA": 1}


def test_bonds_to_another_chain_or_through_side_chains_do_not_count(tmp_path):
    # A HETATM residue whose N is 1.33 A from a C of chain B, and a retinal-like ligand
    # 1.4 A from a lysine NZ (Schiff base; not a peptide bond): both stay hetero.
    lines = (_residue("ATOM", "B", 1, "ALA", 0.0, 1)
             + _residue("HETATM", "A", 2, "DLE", 3.33, 5)
             + [_line("ATOM", 9, "NZ", "LYS", "A", 5, 0.0, 10.0, 0.0, "N"),
                _line("HETATM", 10, "C15", "RET", "A", 301, 1.4, 10.0, 0.0, "C"),
                _line("HETATM", 11, "C", "RET", "A", 301, 2.9, 10.0, 0.0, "C")])
    frame, report = load_structure_report(_write(tmp_path, lines))
    assert _hetero_names(frame) == ["DLE", "RET"]
    assert report["polymer_hetatm_residues"] == {}


def test_modres_and_seqres_mark_polymer_without_an_atom_record_neighbour(tmp_path):
    # MSE listed in MODRES but separated from the chain (gap); a two-residue HETATM peptide
    # whose names are in SEQRES; the same peptide without SEQRES stays hetero (a peptide ligand).
    header = ["SEQRES   1 C    2  DLE DVA",
              "MODRES 9XYZ MSE A    5  MET  SELENOMETHIONINE"]
    body = (_residue("ATOM", "A", 1, "ALA", 0.0, 1) + _residue("HETATM", "A", 5, "MSE", 30.0, 5)
            + _residue("HETATM", "C", 1, "DLE", 0.0, 9, y0=20.0)
            + _residue("HETATM", "C", 2, "DVA", 3.33, 13, y0=20.0))
    frame, report = load_structure_report(_write(tmp_path, header + body))
    assert _hetero_names(frame) == []
    assert report["polymer_hetatm_residues"] == {"DLE": 1, "DVA": 1, "MSE": 1}
    frame, report = load_structure_report(_write(tmp_path, body, "plain.pdb"))
    assert _hetero_names(frame) == ["DLE", "DVA", "MSE"]


def test_builtin_mmcif_reader_applies_the_same_rule(tmp_path):
    rows = []
    for line in _chain_with_ligands():
        record, serial, name = line[:6].strip(), int(line[6:11]), line[12:16].strip()
        resname, chain, resid = line[17:20].strip(), line[21], int(line[22:26])
        x, y, z = (float(line[c:c + 8]) for c in (30, 38, 46))
        rows.append(f"{record} {serial} {line[76:78].strip()} {name} {resname} {chain} {resid} {x} {y} {z} 1")
    text = "\n".join(["data_test", "loop_"] + [f"_atom_site.{tag}" for tag in (
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id",
        "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num")] + rows) + "\n#\n"
    path = tmp_path / "model.cif"
    path.write_text(text)
    frame, report = load_structure_report(path, parser="builtin")
    assert _hetero_names(frame) == ["HOH", "LIG"]
    assert report["polymer_hetatm_residues"] == {"DLE": 1}


def test_mdanalysis_reader_adds_peptide_bonded_residues_to_protein(tmp_path):
    pytest.importorskip("MDAnalysis")
    from crevice.trajectory import load_trajectory_report
    path = _write(tmp_path, _chain_with_ligands())
    traj, report = load_trajectory_report(path)
    frame = traj.frames[0]
    assert sorted({a.resname for a in frame.atoms}) == ["ALA", "DLE"]
    assert not any(a.hetero for a in frame.atoms)
    assert report["polymer_hetatm_residues"] == {"DLE": 1}
    traj, report = load_trajectory_report(path, hetero_polymer="exclude")
    assert sorted({a.resname for a in traj.frames[0].atoms}) == ["ALA"]
    assert report["polymer_hetatm_residues"] == {}


def test_cli_option_and_note(tmp_path, capsys):
    from crevice.cli import _add_input_args, _load_input
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("structure")
    _add_input_args(parser)
    path = _write(tmp_path, _chain_with_ligands())
    frame, _, report = _load_input(parser.parse_args([str(path)]))
    assert "DLE x1" in capsys.readouterr().err
    assert _hetero_names(frame) == ["HOH", "LIG"]
    frame, _, report = _load_input(parser.parse_args([str(path), "--hetero-polymer", "exclude"]))
    assert capsys.readouterr().err == ""
    assert _hetero_names(frame) == ["DLE", "HOH", "LIG"]


def test_polymer_linked_residues_needs_a_polymer_anchor():
    keys = [1, 1, 2, 2]
    xyz = [(0, 0, 0), (2, 0, 0), (3.33, 0, 0), (5, 0, 0)]
    names = ["N", "C", "N", "C"]
    assert polymer_linked_residues(keys, "AAAA", names, xyz, polymer=set(), candidates={1, 2}) == set()
    assert polymer_linked_residues(keys, "AAAA", names, xyz, polymer=set(), candidates={1, 2},
                                   resnames={1: "DLE", 2: "DVA"}, sequences={"A": {"DLE", "DVA"}}) == {1, 2}

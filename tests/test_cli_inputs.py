"""CLI contracts for structure input, analysis selection and trajectory flags.

Downloads are mocked. These check argument plumbing and what each command
records, not the numerical correctness of any analysis.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from crevice import pdb
from crevice.cli import main

ANALYSES = ("profile", "residues", "cavities", "tunnels", "network", "features")

from test_download_integrity import FakeResponse  # noqa: F401


CIF_TEXT = """\
data_1GRM
_entry.id 1GRM
loop_
_entity.id
_entity.type
1 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
%s
#
"""


def _cif_rows():
    rows = []
    serial = 1
    # A ring of atoms around the z axis, wide enough for a profile to run.
    for index, (x, y, z) in enumerate([
        (6, 0, -3), (-6, 0, -3), (0, 6, -3), (0, -6, -3),
        (6, 0, 3), (-6, 0, 3), (0, 6, 3), (0, -6, 3),
    ]):
        rows.append(
            f"ATOM {serial} C CA ALA A 1 {index + 1} {index + 1} A . . 1.00 10.0 "
            f"{float(x)} {float(y)} {float(z)} 1"
        )
        serial += 1
    return "\n".join(rows)


CIF_PAYLOAD = (CIF_TEXT % _cif_rows()).encode()


@pytest.fixture
def network(monkeypatch):
    calls: list[str] = []

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        return FakeResponse(CIF_PAYLOAD)

    monkeypatch.setattr(pdb, "urlopen", fake_urlopen)
    return calls


@pytest.fixture
def local_cif(tmp_path):
    path = tmp_path / "local.cif"
    path.write_bytes(CIF_PAYLOAD)
    return path


def manifest_for(out_dir, prefix):
    return json.loads((Path(out_dir) / f"{prefix}_manifest.json").read_text())


# --- structure input --------------------------------------------------------

def test_analyze_accepts_a_local_file(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--samples", "5", "--search-radius", "0"]) == 0
    data = manifest_for(out, "local")
    assert data["source"]["kind"] == "local_file"
    assert data["source"]["download"] is None
    assert network == []


def test_analyze_accepts_an_accession_and_fetches_it(tmp_path, network):
    out = tmp_path / "out"
    assert main(["analyze", "1GRM", "--cache-dir", str(tmp_path / "cache"),
                 "--out-dir", str(out), "--samples", "5", "--search-radius", "0"]) == 0
    data = manifest_for(out, "1GRM")
    assert data["source"]["kind"] == "downloaded"
    assert data["source"]["accession"] == "1GRM"
    assert network == ["https://files.rcsb.org/download/1GRM.cif"]


def test_a_cached_accession_is_not_refetched(tmp_path, network):
    cache = tmp_path / "cache"
    for index in range(2):
        assert main(["analyze", "1GRM", "--cache-dir", str(cache),
                     "--out-dir", str(tmp_path / f"out{index}"),
                     "--samples", "5", "--search-radius", "0"]) == 0
    assert len(network) == 1


def test_offline_uncached_accession_fails_with_advice(tmp_path, network, capsys):
    assert main(["analyze", "1GRM", "--cache-dir", str(tmp_path / "cache"),
                 "--out-dir", str(tmp_path / "out"), "--offline"]) == 2
    assert "--offline" in capsys.readouterr().err
    assert network == []


def test_an_unresolvable_argument_is_reported(tmp_path, network, capsys):
    assert main(["analyze", "not-a-thing", "--out-dir", str(tmp_path / "out")]) == 2
    assert "neither an existing file nor a PDB accession" in capsys.readouterr().err


def test_assembly_flag_reaches_the_download_url(tmp_path, network):
    assert main(["analyze", "1GRM", "--cache-dir", str(tmp_path / "cache"),
                 "--assembly", "1", "--out-dir", str(tmp_path / "out"),
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert network == ["https://files.rcsb.org/download/1GRM-assembly1.cif"]


def test_parser_choice_is_recorded_in_the_manifest(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out), "--parser", "builtin",
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert manifest_for(out, "local")["parser_report"]["parser"] == "builtin"


def test_default_parser_is_the_validated_one(tmp_path, local_cif, network):
    pytest.importorskip("gemmi")
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert manifest_for(out, "local")["parser_report"]["parser"] == "gemmi"


# --- analysis selection -----------------------------------------------------

def test_default_analyses_are_profile_and_residues(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--samples", "5", "--search-radius", "0"]) == 0
    data = manifest_for(out, "local")
    assert data["analyses"] == ["profile", "residues"]
    assert {"profile.json", "profile.csv", "residues.csv"} <= set(data["files"])
    assert "residues.json" not in data["files"]  # the residue table is the CSV
    assert {"hydration_json", "hydration_residues_csv", "hydration_frames_npz"} <= set(data["files"])


def test_requested_analyses_select_the_outputs(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--analysis", "cavities", "--analysis", "features",
                 "--samples", "5", "--search-radius", "0",
                 "--cavity-spacing", "3"]) == 0
    data = manifest_for(out, "local")
    assert data["analyses"] == ["cavities", "features"]
    assert "cavities.json" in data["files"]
    assert "profile.json" not in data["files"]
    for path in data["files"].values():
        assert Path(path).is_file()


def test_repeated_analyses_are_deduplicated(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--analysis", "profile", "--analysis", "profile",
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert manifest_for(out, "local")["analyses"] == ["profile"]


@pytest.mark.parametrize("analysis", ANALYSES)
def test_every_analysis_choice_runs_and_writes_its_outputs(tmp_path, local_cif, network, analysis):
    """Each choice must actually execute; a subset check hid a missing import."""
    out = tmp_path / analysis
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--analysis", analysis, "--samples", "5", "--search-radius", "0",
                 "--cavity-spacing", "3", "--tunnel-spacing", "3"]) == 0
    data = manifest_for(out, "local")
    assert data["analyses"] == [analysis]
    assert data["files"], f"{analysis} wrote no files"
    for path in data["files"].values():
        assert Path(path).is_file(), f"{analysis} listed a file it did not write"


def test_all_analyses_together_run_in_one_invocation(tmp_path, local_cif, network):
    out = tmp_path / "all"
    command = ["analyze", str(local_cif), "--out-dir", str(out),
               "--samples", "5", "--search-radius", "0",
               "--cavity-spacing", "3", "--tunnel-spacing", "3"]
    for analysis in ANALYSES:
        command += ["--analysis", analysis]
    assert main(command) == 0
    data = manifest_for(out, "local")
    assert data["analyses"] == list(ANALYSES)
    for path in data["files"].values():
        assert Path(path).is_file()


def test_unknown_analysis_is_rejected_by_the_parser(tmp_path, local_cif):
    with pytest.raises(SystemExit):
        main(["analyze", str(local_cif), "--out-dir", str(tmp_path), "--analysis", "magic"])


def test_manifest_states_that_nothing_is_curated(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out),
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert "no curated axis" in manifest_for(out, "local")["interpretation"]


def test_prefix_overrides_the_output_names(tmp_path, local_cif, network):
    out = tmp_path / "out"
    assert main(["analyze", str(local_cif), "--out-dir", str(out), "--prefix", "run7",
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert (out / "run7_manifest.json").is_file()


# --- fetch ------------------------------------------------------------------

def test_fetch_accepts_an_arbitrary_accession(tmp_path, network, capsys):
    assert main(["fetch", "7XYZ", "--cache-dir", str(tmp_path)]) == 2
    # 7XYZ is not the payload's accession, so identity checking must reject it.
    assert "reports accession 1GRM" in capsys.readouterr().err


def test_fetch_writes_provenance_records(tmp_path, network):
    report = tmp_path / "downloads.json"
    assert main(["fetch", "1GRM", "--cache-dir", str(tmp_path / "cache"),
                 "--json", str(report)]) == 0
    downloads = json.loads(report.read_text())["downloads"]
    assert downloads[0]["pdb_id"] == "1GRM"
    assert downloads[0]["identity_verified"] is True


def test_fetch_list_still_describes_the_benchmarks(network, capsys):
    assert main(["fetch", "--list"]) == 0
    out = capsys.readouterr().out
    assert "set_aside" in out and "2CHB" in out
    assert network == []


def test_fetch_rejects_a_malformed_accession(tmp_path, network, capsys):
    assert main(["fetch", "nope", "--cache-dir", str(tmp_path)]) == 2
    assert "not a PDB accession" in capsys.readouterr().err
    assert network == []


# --- trajectory -------------------------------------------------------------

@pytest.fixture
def gromacs_pair(tmp_path):
    from test_trajectory_input import build_universe
    import MDAnalysis as mda
    import numpy as np

    universe = build_universe()
    count = len(universe.atoms)
    topology = tmp_path / "system.gro"
    trajectory = tmp_path / "run.xtc"
    universe.atoms.write(str(topology))
    base = np.arange(count * 3, dtype=np.float32).reshape(count, 3)
    with mda.Writer(str(trajectory), n_atoms=count) as writer:
        for index in range(5):
            universe.atoms.positions = base + index * 10.0
            writer.write(universe.atoms)
    return topology, trajectory


def test_trajectory_inspect_reports_without_analyzing(tmp_path, gromacs_pair, capsys):
    topology, trajectory = gromacs_pair
    report = tmp_path / "info.json"
    assert main(["trajectory", str(trajectory), "--topology", str(topology),
                 "--inspect", "-o", str(tmp_path / "unused.json"),
                 "--report-json", str(report)]) == 0
    assert json.loads(report.read_text())["frame_count"] == 5
    assert not (tmp_path / "unused.json").exists()
    assert "frame_count" in capsys.readouterr().out


def test_trajectory_inspect_requires_a_topology(tmp_path, gromacs_pair, capsys):
    _, trajectory = gromacs_pair
    assert main(["trajectory", str(trajectory), "--inspect",
                 "-o", str(tmp_path / "x.json")]) == 2
    assert "--topology" in capsys.readouterr().err


def test_trajectory_inspect_needs_no_output_path(gromacs_pair, capsys):
    """The README example: --inspect prints to stdout, so -o is not required."""
    topology, trajectory = gromacs_pair
    assert main(["trajectory", str(trajectory), "--topology", str(topology), "--inspect"]) == 0
    out = capsys.readouterr().out
    assert "frame_count\t5" in out


def test_trajectory_analysis_still_requires_an_output_path(gromacs_pair, capsys):
    topology, trajectory = gromacs_pair
    assert main(["trajectory", str(trajectory), "--topology", str(topology)]) == 2
    assert "-o/--output" in capsys.readouterr().err


def test_trajectory_flags_control_the_frames_read(tmp_path, gromacs_pair):
    topology, trajectory = gromacs_pair
    report = tmp_path / "reader.json"
    assert main(["trajectory", str(trajectory), "--topology", str(topology),
                 "--selection", "name CA", "--start", "1", "--stop", "5", "--stride", "2",
                 "-o", str(tmp_path / "traj.json"), "--report-json", str(report),
                 "--samples", "5", "--search-radius", "0"]) == 0
    data = json.loads(report.read_text())
    assert data["loaded_frame_count"] == 2
    assert data["selection"] == "name CA"
    assert data["selected_atom_count"] == 4
    assert data["available_frame_count"] == 5


def test_trajectory_max_frames_guard_is_reachable_from_the_cli(tmp_path, gromacs_pair, capsys):
    topology, trajectory = gromacs_pair
    assert main(["trajectory", str(trajectory), "--topology", str(topology),
                 "--max-frames", "2", "-o", str(tmp_path / "traj.json")]) == 2
    assert "raise stride" in capsys.readouterr().err


def test_trajectory_stride_is_not_applied_twice(tmp_path, gromacs_pair):
    """The reader already subsamples, so the analysis must not stride again."""
    topology, trajectory = gromacs_pair
    output = tmp_path / "traj.json"
    assert main(["trajectory", str(trajectory), "--topology", str(topology),
                 "--stride", "2", "-o", str(output),
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert len(json.loads(output.read_text())["frames"]) == 3


def test_static_frame_sequence_still_works(tmp_path, local_cif, network):
    output = tmp_path / "traj.json"
    assert main(["trajectory", str(local_cif), str(local_cif), "-o", str(output),
                 "--samples", "5", "--search-radius", "0"]) == 0
    assert len(json.loads(output.read_text())["frames"]) == 2

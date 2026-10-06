"""Accession grammar, general fetching and structure-input resolution.

Every network call is mocked. These contracts cover which entry is requested,
where it lands and how it is described; they do not contact RCSB.
"""

from __future__ import annotations

import pytest

from crevice import benchmarks, pdb

from test_download_integrity import CIF_PAYLOAD, PDB_PAYLOAD, FakeResponse, serve, transport  # noqa: F401


# --- accession grammar ------------------------------------------------------

@pytest.mark.parametrize("text,expected", [
    ("1GRM", "1GRM"), ("1grm", "1GRM"), ("4pyp", "4PYP"), ("101M", "101M"),
    ("7A1B", "7A1B"), (" 1grm ", "1GRM"),
    ("pdb_00001grm", "pdb_00001grm"), ("PDB_00001GRM", "pdb_00001grm"),
])
def test_valid_accessions_normalize(text, expected):
    assert pdb.is_accession(text)
    assert pdb.normalize_accession(text) == expected


@pytest.mark.parametrize("text", [
    "GRM1", "1GR", "1GRMX", "", "   ", "structure.pdb", "1GR M", "protein",
    "./1GRM", "pdb_0001grm", "pdb_00001grmx",
])
def test_invalid_accessions_are_rejected(text):
    assert not pdb.is_accession(text)
    with pytest.raises(ValueError, match="not a PDB accession"):
        pdb.normalize_accession(text)


# --- locations --------------------------------------------------------------

@pytest.mark.parametrize("fmt,assembly,expected", [
    ("cif", None, "https://files.rcsb.org/download/1GRM.cif"),
    ("mmcif", None, "https://files.rcsb.org/download/1GRM.cif"),
    ("pdb", None, "https://files.rcsb.org/download/1GRM.pdb"),
    ("cif", 1, "https://files.rcsb.org/download/1GRM-assembly1.cif"),
    ("cif", 2, "https://files.rcsb.org/download/1GRM-assembly2.cif"),
    ("pdb", 1, "https://files.rcsb.org/download/1GRM.pdb1"),
])
def test_download_urls_follow_rcsb_naming(fmt, assembly, expected):
    assert pdb.structure_url("1grm", fmt=fmt, assembly=assembly) == expected


@pytest.mark.parametrize("fmt,assembly,expected", [
    ("cif", None, "1GRM.cif"), ("pdb", None, "1GRM.pdb"),
    ("cif", 1, "1GRM-assembly1.cif"), ("pdb", 2, "1GRM-assembly2.pdb"),
])
def test_cache_filenames_always_carry_the_format_suffix(tmp_path, fmt, assembly, expected):
    """The legacy .pdb<n> assembly name is normalised locally."""
    path = pdb.cached_structure_path(tmp_path, "1grm", fmt=fmt, assembly=assembly)
    assert path == tmp_path / expected


@pytest.mark.parametrize("bad", [0, -1, True, 1.5, "1"])
def test_assembly_must_be_a_positive_integer(bad):
    with pytest.raises(ValueError, match="assembly must be a positive integer"):
        pdb.structure_url("1GRM", assembly=bad)


def test_unknown_format_is_rejected():
    with pytest.raises(ValueError, match="format must be pdb, cif, or mmcif"):
        pdb.structure_url("1GRM", fmt="xyz")


# --- general fetching -------------------------------------------------------

def test_any_accession_can_be_fetched_not_only_benchmarks(tmp_path, transport):
    """The five named benchmarks are metadata, not a download allow-list."""
    record = pdb.fetch_structure("7XYZ", tmp_path, fmt="pdb")
    assert record.path == tmp_path / "7XYZ.pdb"
    assert record.identity_verified is True
    assert transport.calls[0]["url"] == "https://files.rcsb.org/download/7XYZ.pdb"


def test_invalid_accession_is_rejected_before_any_request(tmp_path, transport):
    with pytest.raises(ValueError, match="not a PDB accession"):
        pdb.fetch_structure("not-an-id", tmp_path)
    assert transport.calls == []


def test_assembly_download_uses_the_assembly_url_and_records_it(tmp_path, transport):
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id XXXX"))
    record = pdb.fetch_structure("1GRM", tmp_path, fmt="cif", assembly=1)
    assert transport.calls[0]["url"].endswith("1GRM-assembly1.cif")
    assert record.path == tmp_path / "1GRM-assembly1.cif"
    assert record.assembly == 1


def test_assembly_placeholder_identity_is_recorded_not_assumed(tmp_path, transport):
    """An expanded assembly states no accession, so it cannot be verified."""
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id XXXX"))
    record = pdb.fetch_structure("1GRM", tmp_path, fmt="cif", assembly=1)
    assert record.identity_verified is False
    assert "placeholder" in record.identity_note
    assert record.path.exists()


def test_assembly_accession_mismatch_does_not_abort_the_fetch(tmp_path, transport):
    """Assembly files are not required to state the requested accession."""
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id 9ZZZ"))
    record = pdb.fetch_structure("1GRM", tmp_path, fmt="cif", assembly=1)
    assert record.identity_verified is False
    assert "9ZZZ" in record.identity_note


def test_asymmetric_unit_accession_mismatch_still_aborts(tmp_path, transport):
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id 9ZZZ"))
    with pytest.raises(RuntimeError, match="reports accession 9ZZZ"):
        pdb.fetch_structure("1GRM", tmp_path, fmt="cif")


def test_asymmetric_and_assembly_files_do_not_collide_in_the_cache(tmp_path, transport):
    serve(transport, CIF_PAYLOAD)
    unit = pdb.fetch_structure("1GRM", tmp_path, fmt="cif")
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id XXXX"))
    assembly = pdb.fetch_structure("1GRM", tmp_path, fmt="cif", assembly=1)
    assert unit.path != assembly.path
    assert unit.path.exists() and assembly.path.exists()


# --- input resolution -------------------------------------------------------

def test_existing_local_file_resolves_without_network(tmp_path, transport):
    structure = tmp_path / "local.pdb"
    structure.write_bytes(PDB_PAYLOAD)
    source = pdb.resolve_structure(structure, cache_dir=tmp_path / "cache")
    assert source.kind == "local_file"
    assert source.path == structure
    assert source.accession is None
    assert transport.calls == []


def test_a_local_file_named_like_an_accession_wins(tmp_path, transport, monkeypatch):
    """A path that exists is never shadowed by a same-named PDB entry."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "1GRM").write_bytes(PDB_PAYLOAD)
    source = pdb.resolve_structure("1GRM", cache_dir=tmp_path / "cache")
    assert source.kind == "local_file"
    assert transport.calls == []


def test_uncached_accession_is_downloaded(tmp_path, transport):
    source = pdb.resolve_structure("1GRM", cache_dir=tmp_path, fmt="pdb")
    assert source.kind == "downloaded"
    assert source.accession == "1GRM"
    assert source.record.sha256
    assert len(transport.calls) == 1


def test_cached_accession_is_reused_without_network(tmp_path, transport):
    pdb.fetch_structure("1GRM", tmp_path, fmt="pdb")
    source = pdb.resolve_structure("1GRM", cache_dir=tmp_path, fmt="pdb")
    assert source.kind == "cached_accession"
    assert len(transport.calls) == 1


def test_offline_uncached_accession_names_the_cache_directory(tmp_path, transport):
    with pytest.raises(RuntimeError, match=str(tmp_path)):
        pdb.resolve_structure("1GRM", cache_dir=tmp_path, fmt="pdb", offline=True)
    assert transport.calls == []


def test_offline_cached_accession_still_resolves(tmp_path, transport):
    pdb.fetch_structure("1GRM", tmp_path, fmt="pdb")
    source = pdb.resolve_structure("1GRM", cache_dir=tmp_path, fmt="pdb", offline=True)
    assert source.kind == "cached_accession"


def test_unusable_spec_names_both_interpretations(tmp_path, transport):
    with pytest.raises(ValueError, match="neither an existing file nor a PDB accession"):
        pdb.resolve_structure("missing.pdb", cache_dir=tmp_path)
    assert transport.calls == []


def test_directory_spec_is_rejected(tmp_path, transport):
    with pytest.raises(ValueError, match="is a directory"):
        pdb.resolve_structure(tmp_path, cache_dir=tmp_path)


def test_source_is_json_serializable(tmp_path, transport):
    import json

    source = pdb.resolve_structure("1GRM", cache_dir=tmp_path, fmt="pdb")
    data = json.loads(json.dumps(source.to_dict()))
    assert data["kind"] == "downloaded"
    assert data["download"]["pdb_id"] == "1GRM"


def test_local_source_reports_no_download(tmp_path, transport):
    structure = tmp_path / "local.pdb"
    structure.write_bytes(PDB_PAYLOAD)
    source = pdb.resolve_structure(structure, cache_dir=tmp_path)
    assert source.to_dict()["download"] is None


# --- benchmark layer --------------------------------------------------------

def test_benchmark_fetch_still_restricts_to_registered_systems(tmp_path, transport):
    with pytest.raises(ValueError, match="Unknown benchmark system"):
        benchmarks.fetch_benchmark_record("7XYZ", tmp_path, extension="pdb")
    assert transport.calls == []


def test_unknown_benchmark_error_points_at_the_general_fetcher(tmp_path):
    with pytest.raises(ValueError, match=r"crevice\.pdb\.fetch_structure"):
        benchmarks.get_benchmark_system("7XYZ")


def test_benchmark_fetch_delegates_to_the_general_transport(tmp_path, transport):
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    assert isinstance(record, pdb.DownloadRecord)
    assert record.path == tmp_path / "1GRM.pdb"


# --- extended accessions and the stated entry ID ------------------------------
# Offline, with synthetic mmCIF text: RCSB currently states classic IDs in
# `_entry.id` (1GRM) and will state extended ones (pdb_00001grm) for new
# entries, so a request in either form must match either stated form.

def _cif(entry_id: str) -> str:
    return f"data_x\n_entry.id {entry_id}\nloop_\n_atom_site.id\n1\n#\n"


@pytest.mark.parametrize("requested", ["pdb_00001grm", "PDB_00001GRM", "1GRM", "1grm"])
@pytest.mark.parametrize("stated", ["1GRM", "pdb_00001grm", "PDB_00001GRM"])
def test_identity_matches_classic_and_extended_forms(requested, stated):
    from crevice.pdb import check_identity

    verified, note = check_identity(_cif(stated), requested, "cif", "synthetic", strict=True)
    assert verified, note


def test_extended_identity_mismatch_is_still_rejected():
    from crevice.pdb import check_identity

    with pytest.raises(RuntimeError, match="reports accession"):
        check_identity(_cif("9ZZZ"), "pdb_00001grm", "cif", "synthetic", strict=True)
    with pytest.raises(RuntimeError, match="reports accession"):
        check_identity(_cif("pdb_00009zzz"), "1GRM", "cif", "synthetic", strict=True)
    # A genuinely extended ID (non-zero prefix) has no classic equivalent.
    with pytest.raises(RuntimeError, match="reports accession"):
        check_identity(_cif("1GRM"), "pdb_10001grm", "cif", "synthetic", strict=True)


def test_strict_fetch_of_an_extended_accession_verifies(tmp_path, monkeypatch):
    from crevice import pdb as pdb_module

    class Response:
        def __init__(self, payload):
            self.payload, self.done = payload, False

        def read(self, size=-1):
            if self.done:
                return b""
            self.done = True
            return self.payload

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    payload = (_cif("1GRM").replace("loop_\n_atom_site.id\n1\n",
                                    "loop_\n_atom_site.group_PDB\n_atom_site.id\nATOM 1\n")).encode()
    monkeypatch.setattr(pdb_module, "urlopen", lambda url, timeout=None: Response(payload))
    record = pdb_module.fetch_structure("pdb_00001grm", tmp_path, fmt="cif")
    assert record.identity_verified is True
    assert record.pdb_id == "pdb_00001grm"
    # The cached-entry path must normalise the same way.
    again = pdb_module.cached_record("pdb_00001grm", record.path, "cif")
    assert again.pdb_id == "pdb_00001grm"
    record.path.with_name(record.path.name + ".provenance.json").unlink()
    fresh = pdb_module.cached_record("pdb_00001grm", record.path, "cif")
    assert fresh.pdb_id == "pdb_00001grm" and fresh.identity_verified is True


# --- default formats ------------------------------------------------------------

def test_benchmark_helpers_default_to_the_same_format_as_fetch_and_the_cli():
    import inspect
    from crevice.cli import build_parser

    expected = inspect.signature(pdb.fetch_structure).parameters["fmt"].default
    assert expected == "cif"
    for function in (benchmarks.benchmark_path, benchmarks.fetch_benchmark_structure,
                     benchmarks.fetch_benchmark_record, benchmarks.fetch_all_benchmark_structures):
        assert inspect.signature(function).parameters["extension"].default == expected, function.__name__
    parser = build_parser()
    for command in (["fetch", "1GRM"], ["benchmark"], ["static-suite"]):
        assert parser.parse_args(command).format == expected

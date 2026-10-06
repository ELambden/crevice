"""Transport, cache and provenance contracts for benchmark structure downloads.

Every network call here is mocked. These contracts cover payload integrity and
recorded provenance; they do not contact RCSB and do not validate that a real
deposited entry is a correct benchmark for any biological question.
"""

from __future__ import annotations

import hashlib
import json
from urllib.error import HTTPError, URLError

import pytest

from crevice import benchmarks, pdb


PDB_PAYLOAD = b"""HEADER    TRANSPORT PROTEIN                       01-JAN-00   1GRM
ATOM      1  N   VAL A   1      -1.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  VAL A   1       0.000   0.000   0.000  1.00 10.00           C
END
"""

CIF_PAYLOAD = b"""data_1GRM
_entry.id 1GRM
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.Cartn_x
ATOM 1 0.0
#
"""

HTML_PAYLOAD = b"<html><head><title>404 Not Found</title></head><body>No entry</body></html>"


class FakeResponse:
    """Minimal urlopen stand-in supporting chunked reads and mid-stream failure."""

    def __init__(self, payload: bytes, *, fail_after: int | None = None):
        self.payload = payload
        self.fail_after = fail_after
        self.offset = 0
        self.reads = 0

    def read(self, size: int = -1) -> bytes:
        self.reads += 1
        if self.fail_after is not None and self.reads > self.fail_after:
            raise URLError("connection reset")
        if size is None or size < 0:
            size = len(self.payload) - self.offset
        chunk = self.payload[self.offset : self.offset + size]
        self.offset += len(chunk)
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


@pytest.fixture
def transport(monkeypatch):
    """Record every request and serve a scripted response."""

    calls: list[dict] = []
    state = {"response": None}

    def fake_urlopen(url, timeout=None):
        calls.append({"url": url, "timeout": timeout})
        factory = state["response"]
        result = factory() if factory is not None else FakeResponse(_payload_for(url))
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(pdb, "urlopen", fake_urlopen)
    return type("Transport", (), {"calls": calls, "state": state})()


def _payload_for(url: str) -> bytes:
    """Serve each accession its own payload, as RCSB would."""
    accession = url.rsplit("/", 1)[-1].split(".")[0].encode()
    return PDB_PAYLOAD.replace(b"1GRM", accession)


def serve(transport, payload=None, **kwargs):
    transport.state["response"] = lambda: FakeResponse(payload, **kwargs)


def cached(tmp_path, pdb_id="1GRM", extension="pdb"):
    return benchmarks.benchmark_path(tmp_path, pdb_id, extension=extension)


# --- payload integrity ------------------------------------------------------

def test_valid_pdb_payload_is_cached(tmp_path, transport):
    path = benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert path == cached(tmp_path)
    assert path.read_bytes() == PDB_PAYLOAD
    assert transport.calls[0]["url"] == "https://files.rcsb.org/download/1GRM.pdb"


def test_html_error_page_is_not_cached_as_a_structure(tmp_path, transport):
    serve(transport, HTML_PAYLOAD)
    with pytest.raises(RuntimeError, match="markup"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert not cached(tmp_path).exists()


def test_pdb_payload_without_coordinates_is_rejected(tmp_path, transport):
    serve(transport, b"HEADER    NOTHING HERE\nEND\n")
    with pytest.raises(RuntimeError, match="no ATOM or HETATM"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert not cached(tmp_path).exists()


def test_cif_payload_without_atom_site_is_rejected(tmp_path, transport):
    serve(transport, b"data_1GRM\n_entry.id 1GRM\n")
    with pytest.raises(RuntimeError, match="atom_site"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="cif")
    assert not cached(tmp_path, extension="cif").exists()


def test_valid_cif_payload_is_cached(tmp_path, transport):
    serve(transport, CIF_PAYLOAD)
    path = benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="cif")
    assert path.read_bytes() == CIF_PAYLOAD
    assert transport.calls[0]["url"].endswith("1GRM.cif")


def test_empty_payload_is_rejected(tmp_path, transport):
    serve(transport, b"   \n")
    with pytest.raises(RuntimeError, match="empty payload"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert not cached(tmp_path).exists()


def test_oversized_payload_is_rejected_before_writing(tmp_path, transport):
    serve(transport, b"ATOM  " + b"x" * 4096)
    with pytest.raises(RuntimeError, match="exceeds the 1024 byte limit"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, max_bytes=1024, extension="pdb")
    assert not cached(tmp_path).exists()


def test_max_bytes_must_be_positive(tmp_path, transport):
    with pytest.raises(ValueError, match="max_bytes must be a positive integer"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, max_bytes=0, extension="pdb")


# --- failure handling -------------------------------------------------------

def test_interrupted_download_leaves_no_partial_file(tmp_path, transport):
    serve(transport, PDB_PAYLOAD, fail_after=1)
    with pytest.raises(RuntimeError, match="Could not fetch 1GRM"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, max_bytes=64, extension="pdb")
    assert list(tmp_path.iterdir()) == []


def test_http_error_names_the_url(tmp_path, transport):
    transport.state["response"] = lambda: HTTPError(
        "https://files.rcsb.org/download/1GRM.pdb", 404, "Not Found", {}, None
    )
    with pytest.raises(RuntimeError, match="https://files.rcsb.org/download/1GRM.pdb"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert not cached(tmp_path).exists()


def test_failed_download_preserves_an_existing_cache_entry(tmp_path, transport):
    benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    transport.state["response"] = lambda: URLError("offline")
    with pytest.raises(RuntimeError):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, overwrite=True, extension="pdb")
    assert cached(tmp_path).read_bytes() == PDB_PAYLOAD


def test_non_https_download_template_is_rejected(tmp_path, transport, monkeypatch):
    monkeypatch.setattr(pdb, "RCSB_DOWNLOAD", "ftp://files.rcsb.org/{filename}")
    with pytest.raises(ValueError, match="https"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert transport.calls == []


# --- caching ----------------------------------------------------------------

def test_existing_cache_entry_is_not_refetched(tmp_path, transport):
    benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert len(transport.calls) == 1


def test_overwrite_refetches_and_replaces(tmp_path, transport):
    benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    replacement = PDB_PAYLOAD.replace(b"VAL", b"ALA")
    serve(transport, replacement)
    benchmarks.fetch_benchmark_structure("1GRM", tmp_path, overwrite=True, extension="pdb")
    assert cached(tmp_path).read_bytes() == replacement
    assert len(transport.calls) == 2


def test_timeout_is_passed_to_the_transport(tmp_path, transport):
    benchmarks.fetch_benchmark_structure("1GRM", tmp_path, timeout=7.5, extension="pdb")
    assert transport.calls[0]["timeout"] == 7.5


def test_fetch_all_returns_one_path_per_system(tmp_path, transport):
    paths = benchmarks.fetch_all_benchmark_structures(tmp_path, extension="pdb")
    assert len(paths) == len(benchmarks.DEFAULT_BENCHMARK_IDS)
    assert {path.stem for path in paths} == set(benchmarks.DEFAULT_BENCHMARK_IDS)


# --- provenance -------------------------------------------------------------

def test_download_records_hash_size_and_source(tmp_path, transport):
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    assert record.pdb_id == "1GRM"
    assert record.sha256 == hashlib.sha256(PDB_PAYLOAD).hexdigest()
    assert record.byte_count == len(PDB_PAYLOAD)
    assert record.url == "https://files.rcsb.org/download/1GRM.pdb"
    assert record.content_format == "pdb"
    assert record.source == "download"
    assert record.retrieved_utc.endswith("Z")


def test_provenance_sidecar_is_written_next_to_the_structure(tmp_path, transport):
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    sidecar = tmp_path / "1GRM.pdb.provenance.json"
    assert json.loads(sidecar.read_text())["sha256"] == record.sha256
    assert benchmarks.read_provenance(record.path).sha256 == record.sha256


def test_verify_cached_structure_accepts_an_untouched_file(tmp_path, transport):
    path = benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert benchmarks.verify_cached_structure(path).sha256 == hashlib.sha256(PDB_PAYLOAD).hexdigest()


def test_verify_cached_structure_detects_modified_bytes(tmp_path, transport):
    path = benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    path.write_bytes(PDB_PAYLOAD.replace(b"VAL", b"ALA"))
    with pytest.raises(RuntimeError, match="does not match its recorded SHA-256"):
        benchmarks.verify_cached_structure(path)


def test_verify_cached_structure_requires_a_sidecar(tmp_path):
    path = cached(tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(PDB_PAYLOAD)
    with pytest.raises(RuntimeError, match="no recorded provenance"):
        benchmarks.verify_cached_structure(path)


def test_record_for_a_preexisting_cache_entry_is_marked_as_such(tmp_path):
    """Files cached before provenance recording are reported, not re-labelled."""
    path = cached(tmp_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(PDB_PAYLOAD)
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    assert record.source == "existing_cache"
    assert record.url == ""
    assert record.sha256 == hashlib.sha256(PDB_PAYLOAD).hexdigest()


def test_record_is_json_serializable(tmp_path, transport):
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    assert json.loads(json.dumps(record.to_dict()))["pdb_id"] == "1GRM"


# --- accession identity -----------------------------------------------------

def test_pdb_header_identity_mismatch_is_rejected(tmp_path, transport):
    serve(transport, PDB_PAYLOAD.replace(b"1GRM", b"9ZZZ"))
    with pytest.raises(RuntimeError, match="reports accession 9ZZZ"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="pdb")
    assert not cached(tmp_path).exists()


def test_matching_pdb_header_verifies_identity(tmp_path, transport):
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    assert record.identity_verified is True
    assert "HEADER" in record.identity_note


def test_matching_cif_entry_id_verifies_identity(tmp_path, transport):
    serve(transport, CIF_PAYLOAD)
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="cif")
    assert record.identity_verified is True
    assert "_entry.id" in record.identity_note


def test_cif_entry_id_mismatch_is_rejected(tmp_path, transport):
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id 9ZZZ"))
    with pytest.raises(RuntimeError, match="reports accession 9ZZZ"):
        benchmarks.fetch_benchmark_structure("1GRM", tmp_path, extension="cif")


def test_placeholder_entry_id_is_cached_but_not_identity_verified(tmp_path, transport):
    """An expanded assembly carries no accession; that is recorded, not assumed."""
    serve(transport, CIF_PAYLOAD.replace(b"_entry.id 1GRM", b"_entry.id XXXX"))
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="cif")
    assert record.path.exists()
    assert record.identity_verified is False
    assert "placeholder" in record.identity_note


def test_absent_identity_field_is_recorded_as_unverified(tmp_path, transport):
    serve(transport, b"ATOM      1  N   VAL A   1       0.000   0.000   0.000\nEND\n")
    record = benchmarks.fetch_benchmark_record("1GRM", tmp_path, extension="pdb")
    assert record.identity_verified is False
    assert "no accession" in record.identity_note

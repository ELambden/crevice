"""Offline tests of the example-data fetcher (crevice.examples) and `crevice fetch-example`.

No network access: a small fake dataset is registered and served from a local
directory or a file:// URL.
"""
import gzip
import hashlib

import pytest

from crevice import examples
from crevice.cli import main


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def fake_dataset(tmp_path, monkeypatch):
    source = tmp_path / "assets"
    source.mkdir()
    plain = b"frame data\n" * 10
    unpacked = b"topology line\n" * 20
    packed = gzip.compress(unpacked, mtime=0)
    (source / "t.xtc").write_bytes(plain)
    (source / "t.gro.gz").write_bytes(packed)
    dataset = examples.ExampleDataset(
        name="fake",
        description="test data",
        licence="CC0",
        files=(
            examples.ExampleFile("t.gro.gz", _sha(packed), len(packed), "packed",
                                 unpacked_name="t.gro", unpacked_sha256=_sha(unpacked)),
            examples.ExampleFile("t.xtc", _sha(plain), len(plain), "plain"),
        ),
    )
    monkeypatch.setitem(examples.EXAMPLES, "fake", dataset)
    for name in (examples.BASE_URL_ENV, examples.SOURCE_DIR_ENV, examples.CACHE_DIR_ENV):
        monkeypatch.delenv(name, raising=False)
    return source, plain, unpacked


def test_fetch_from_source_dir_verifies_and_unpacks(fake_dataset, tmp_path):
    source, plain, unpacked = fake_dataset
    fetched = examples.fetch_example("fake", cache_dir=tmp_path / "cache", source_dir=source)
    assert fetched.directory == (tmp_path / "cache" / "fake").resolve()
    assert fetched["t.xtc"].read_bytes() == plain
    assert fetched["t.gro"].read_bytes() == unpacked
    assert sorted(fetched.downloaded) == ["t.gro", "t.gro.gz", "t.xtc"]
    assert not list(fetched.directory.glob("*.part"))
    again = examples.fetch_example("fake", cache_dir=tmp_path / "cache", source_dir=source)
    assert again.downloaded == []
    forced = examples.fetch_example("fake", cache_dir=tmp_path / "cache", source_dir=source, overwrite=True)
    assert sorted(forced.downloaded) == ["t.gro", "t.gro.gz", "t.xtc"]


def test_fetch_from_file_url_and_environment(fake_dataset, tmp_path, monkeypatch):
    source, plain, _ = fake_dataset
    fetched = examples.fetch_example("fake", cache_dir=tmp_path / "c1", base_url=source.as_uri())
    assert fetched["t.xtc"].read_bytes() == plain
    monkeypatch.setenv(examples.BASE_URL_ENV, source.as_uri() + "/")
    monkeypatch.setenv(examples.CACHE_DIR_ENV, str(tmp_path / "c2"))
    fetched = examples.fetch_example("fake")
    assert fetched.directory == (tmp_path / "c2" / "fake").resolve()
    assert fetched["t.gro"].is_file()


def test_hash_mismatch_is_refused_and_removed(fake_dataset, tmp_path):
    source, _, _ = fake_dataset
    (source / "t.xtc").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="SHA-256 mismatch for t.xtc"):
        examples.fetch_example("fake", cache_dir=tmp_path / "cache", source_dir=source)
    assert not (tmp_path / "cache" / "fake" / "t.xtc").exists()


def test_corrupt_cache_is_refetched(fake_dataset, tmp_path):
    source, plain, _ = fake_dataset
    cache = tmp_path / "cache"
    examples.fetch_example("fake", cache_dir=cache, source_dir=source)
    (cache / "fake" / "t.xtc").write_bytes(b"corrupt")
    fetched = examples.fetch_example("fake", cache_dir=cache, source_dir=source)
    assert fetched.downloaded == ["t.xtc"]
    assert fetched["t.xtc"].read_bytes() == plain


def test_unknown_name_and_registry():
    with pytest.raises(KeyError, match="glut1-excerpt"):
        examples.fetch_example("no-such-example")
    glut1 = examples.EXAMPLES["glut1-excerpt"]
    names = [item.name for item in glut1.files]
    assert {"glut1_excerpt.xtc", "glut1_excerpt.gro.gz", "glut1_site_candidate_region.json",
            "glut1_site_candidate_reference.dx"} <= set(names)
    assert all(len(item.sha256) == 64 for item in glut1.files)
    assert glut1.licence == "CC-BY-4.0"
    assert sum(item.size for item in glut1.files) < 15_000_000
    assert examples.DEFAULT_BASE_URL.startswith("https://github.com/") and examples.DEFAULT_BASE_URL.endswith("/")


def test_default_cache_dir_uses_xdg(monkeypatch, tmp_path):
    monkeypatch.delenv(examples.CACHE_DIR_ENV, raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    assert examples.default_cache_dir() == tmp_path / "crevice" / "examples"


def test_cli_fetch_example(fake_dataset, tmp_path, capsys):
    source, _, _ = fake_dataset
    assert main(["fetch-example", "--list"]) == 0
    assert "glut1-excerpt" in capsys.readouterr().out
    assert main(["fetch-example", "fake", "--cache-dir", str(tmp_path / "cache"), "--source-dir", str(source)]) == 0
    out = capsys.readouterr().out
    assert "t.gro\tfetched" in out and "t.xtc\tfetched\tsha256 verified" in out
    assert f"example directory={(tmp_path / 'cache' / 'fake').resolve()}" in out

"""Download example datasets by name, with SHA-256 verification and a local cache.

The example data are too large for the source repository, so they are
distributed as release assets and fetched on demand. Each dataset is a fixed
list of files with known sizes and SHA-256 hashes; a file is used only when
its hash matches, and a ``.gz`` asset is decompressed after download and its
decompressed hash is checked too.

The download location is ``base_url + file name``. It defaults to the GitHub
Release asset of the CREVICE release that introduced the data and can be
changed with the ``base_url`` argument or the ``CREVICE_EXAMPLE_BASE_URL``
environment variable. For offline use, ``source_dir`` (or
``CREVICE_EXAMPLE_SOURCE_DIR``) names a local directory holding the same
files, which are then copied and verified instead of downloaded.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import tempfile
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_BASE_URL = "https://github.com/ELambden/crevice/releases/download/v0.1.1/"
"""Default location of the example assets (the v0.1.1 GitHub Release)."""

BASE_URL_ENV = "CREVICE_EXAMPLE_BASE_URL"
SOURCE_DIR_ENV = "CREVICE_EXAMPLE_SOURCE_DIR"
CACHE_DIR_ENV = "CREVICE_EXAMPLE_CACHE_DIR"


@dataclass(frozen=True)
class ExampleFile:
    """One downloadable file of an example dataset.

    Attributes
    ----------
    name : str
        Asset file name, appended to the base URL and used in the cache.
    sha256 : str
        Expected SHA-256 of the downloaded asset.
    size : int
        Expected size of the downloaded asset in bytes.
    description : str
        Short description of the file's contents.
    unpacked_name : str or None
        For a gzip asset, the name of the decompressed file written next to it.
    unpacked_sha256 : str or None
        Expected SHA-256 of the decompressed file.
    """

    name: str
    sha256: str
    size: int
    description: str
    unpacked_name: str | None = None
    unpacked_sha256: str | None = None


@dataclass(frozen=True)
class ExampleDataset:
    """A named example dataset: its files and a description.

    Attributes
    ----------
    name : str
        Dataset name accepted by :func:`fetch_example`.
    description : str
        One-line description of the data.
    files : tuple of ExampleFile
        The files of the dataset.
    licence : str
        Licence of the data.
    """

    name: str
    description: str
    files: tuple[ExampleFile, ...]
    licence: str


@dataclass
class FetchedExample:
    """Local copies of a fetched example dataset.

    Attributes
    ----------
    name : str
        Dataset name.
    directory : pathlib.Path
        Cache directory holding the verified files.
    paths : dict of str to pathlib.Path
        Local path of every file by name, including decompressed files.
    downloaded : list of str
        Names of the files copied, downloaded or decompressed by this call
        (the others were already cached with the right hash).
    """

    name: str
    directory: Path
    paths: dict[str, Path] = field(default_factory=dict)
    downloaded: list[str] = field(default_factory=list)

    def __getitem__(self, name: str) -> Path:
        """Return the local path of the file called ``name``.

        Parameters
        ----------
        name : str
            File name, e.g. ``"glut1_excerpt.xtc"``.

        Returns
        -------
        pathlib.Path
            Local path of the verified file.
        """
        return self.paths[name]


EXAMPLES: dict[str, ExampleDataset] = {
    "glut1-excerpt": ExampleDataset(
        name="glut1-excerpt",
        description=(
            "21 frames (every 5 ns over 100 ns) of an all-atom membrane MD simulation of a "
            "human GLUT1 mutant (N45Q/E329Q/N411Q): protein, lipids, water and ions (the "
            "simulation's Tris molecules omitted), CHARMM naming, with a prepared region "
            "definition for the candidate sugar site"
        ),
        licence="CC-BY-4.0",
        files=(
            ExampleFile(
                "glut1_excerpt.gro.gz",
                "e427bd5def55f553092795e68478e41cfac4be68b07cd767defc10081cf4a3e1",
                1594631,
                "topology and first-frame coordinates (GRO, gzip)",
                unpacked_name="glut1_excerpt.gro",
                unpacked_sha256="246647c455e83ee8365e8e3788b07e5a3ea733f0d1767eb80af22e02d7656f22",
            ),
            ExampleFile(
                "glut1_excerpt.xtc",
                "1de05da86b0ccfbab03120b8dd16695110fa2f100c4c93c2ddc3055a37888058",
                9718520,
                "21-frame trajectory (XTC)",
            ),
            ExampleFile(
                "glut1_site_candidate_region.json",
                "065ea77977a46587b5969f43595a6b73ff45cb35b4b50b2cb2c23974929a64b1",
                6719,
                "prepared region definition (schema 2) bound to the files above",
            ),
            ExampleFile(
                "glut1_site_candidate_reference.dx",
                "7b231d4cf615a464ecc6588b4b789e1c956a6448705cc5167a020b795fa08e6d",
                376223,
                "fixed reference map of the region",
            ),
            ExampleFile(
                "glut1_excerpt_provenance.json",
                "a5cae6dc712bc1a97e9b9094700e35d14c82e13e4a39dd2339158352a589247f",
                3802,
                "source hashes, frame indices, times and selection",
            ),
            ExampleFile(
                "glut1_excerpt_README.md",
                "8f20c74bb1e718b486870624fa972cfb0def3e0e2a4c6325f3bbfd38cd6f6f24",
                3960,
                "data card: system, files and licence",
            ),
        ),
    ),
}
"""Registry of the example datasets, by name."""


def list_examples() -> list[ExampleDataset]:
    """Return the registered example datasets.

    Returns
    -------
    list of ExampleDataset
        Every dataset that :func:`fetch_example` can fetch, sorted by name.
    """
    return [EXAMPLES[name] for name in sorted(EXAMPLES)]


def default_cache_dir() -> Path:
    """Return the cache root used when no ``cache_dir`` is given.

    Returns
    -------
    pathlib.Path
        ``$CREVICE_EXAMPLE_CACHE_DIR`` if set, otherwise
        ``$XDG_CACHE_HOME/crevice/examples`` (``~/.cache/crevice/examples``).
    """
    configured = os.environ.get(CACHE_DIR_ENV)
    if configured:
        return Path(configured).expanduser()
    root = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(root).expanduser() / "crevice" / "examples"


def sha256_file(path: str | os.PathLike) -> str:
    """Return the SHA-256 hex digest of a file.

    Parameters
    ----------
    path : str or os.PathLike
        File to hash.

    Returns
    -------
    str
        Lower-case hexadecimal digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _matches(path: Path, sha256: str) -> bool:
    return path.is_file() and sha256_file(path) == sha256


def _retrieve(source: str, target: Path, timeout: float) -> None:
    """Copy ``source`` (local path or URL) to ``target`` via a temporary file."""
    target.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=target.parent, prefix=target.name + ".", suffix=".part")
    os.close(handle)
    try:
        if "://" in source:
            request = urllib.request.Request(source, headers={"User-Agent": "crevice-fetch-example"})
            with urllib.request.urlopen(request, timeout=timeout) as response, open(temporary, "wb") as out:
                shutil.copyfileobj(response, out)
        else:
            shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def fetch_example(
    name: str,
    cache_dir: str | os.PathLike | None = None,
    base_url: str | None = None,
    source_dir: str | os.PathLike | None = None,
    overwrite: bool = False,
    timeout: float = 60.0,
) -> FetchedExample:
    """Fetch an example dataset into the local cache and verify every file.

    Files already in the cache with the expected SHA-256 are reused. Others
    are copied from ``source_dir`` or downloaded from ``base_url``, and the
    dataset is returned only if every hash matches.

    Parameters
    ----------
    name : str
        Dataset name, e.g. ``"glut1-excerpt"`` (see :func:`list_examples`).
    cache_dir : str or os.PathLike, optional
        Cache root; the files go in ``cache_dir/name``. Defaults to
        :func:`default_cache_dir`.
    base_url : str, optional
        URL prefix of the assets (a trailing ``/`` is added if missing);
        ``file://`` URLs work. Defaults to ``$CREVICE_EXAMPLE_BASE_URL`` or
        :data:`DEFAULT_BASE_URL`.
    source_dir : str or os.PathLike, optional
        Local directory holding the asset files, used instead of any URL
        (offline use). Defaults to ``$CREVICE_EXAMPLE_SOURCE_DIR`` if set.
    overwrite : bool, default False
        Fetch every file again even if a verified copy is cached.
    timeout : float, default 60.0
        Network timeout in seconds per file.

    Returns
    -------
    FetchedExample
        Cache directory and the local path of every file.

    Raises
    ------
    KeyError
        If ``name`` is not a registered dataset.
    ValueError
        If a fetched or decompressed file does not have the expected hash; the
        bad file is removed.
    """
    if name not in EXAMPLES:
        raise KeyError(f"Unknown example dataset {name!r}; available: {', '.join(sorted(EXAMPLES))}")
    dataset = EXAMPLES[name]
    directory = Path(cache_dir).expanduser() if cache_dir is not None else default_cache_dir()
    directory = (directory / name).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    local = source_dir if source_dir is not None else os.environ.get(SOURCE_DIR_ENV)
    prefix = base_url or os.environ.get(BASE_URL_ENV) or DEFAULT_BASE_URL
    if not prefix.endswith("/"):
        prefix += "/"
    result = FetchedExample(name=name, directory=directory)
    for item in dataset.files:
        target = directory / item.name
        if overwrite or not _matches(target, item.sha256):
            source = str(Path(local).expanduser() / item.name) if local else prefix + urllib.parse.quote(item.name)
            _retrieve(source, target, timeout)
            result.downloaded.append(item.name)
            actual = sha256_file(target)
            if actual != item.sha256:
                target.unlink()
                raise ValueError(f"SHA-256 mismatch for {item.name} from {source}: expected {item.sha256}, got {actual}")
        result.paths[item.name] = target
        if item.unpacked_name:
            unpacked = directory / item.unpacked_name
            if overwrite or not _matches(unpacked, item.unpacked_sha256 or ""):
                handle, temporary = tempfile.mkstemp(dir=directory, prefix=item.unpacked_name + ".", suffix=".part")
                os.close(handle)
                try:
                    with gzip.open(target, "rb") as src, open(temporary, "wb") as out:
                        shutil.copyfileobj(src, out)
                    actual = sha256_file(temporary)
                    if actual != item.unpacked_sha256:
                        raise ValueError(f"SHA-256 mismatch for decompressed {item.unpacked_name}: "
                                         f"expected {item.unpacked_sha256}, got {actual}")
                    os.replace(temporary, unpacked)
                    result.downloaded.append(item.unpacked_name)
                finally:
                    if os.path.exists(temporary):
                        os.unlink(temporary)
            result.paths[item.unpacked_name] = unpacked
    return result

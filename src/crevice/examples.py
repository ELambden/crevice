"""Download example datasets by name, with SHA-256 verification and a local cache.

The example data are too large for the source repository, so they are
distributed as release assets and fetched on demand. Each dataset is a fixed
list of files with known sizes and SHA-256 hashes; a file is used only when
its hash matches, and a ``.gz`` asset is decompressed after download and its
decompressed hash is checked too.

The download location is ``base_url + file name``. It defaults to the GitHub
Release assets of the CREVICE release that introduced the dataset
(:attr:`ExampleDataset.base_url`) and can be changed with the ``base_url``
argument or the ``CREVICE_EXAMPLE_BASE_URL`` environment variable. For offline use, ``source_dir`` (or
``CREVICE_EXAMPLE_SOURCE_DIR``) names a local directory holding the same
files, which are then copied and verified instead of downloaded.
"""

from __future__ import annotations

import gzip
import hashlib
import re
import os
import shutil
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_BASE_URL = "https://github.com/ELambden/crevice/releases/download/v0.1.1/"
"""Default location of the example assets (the v0.1.1 GitHub Release)."""

BASE_URL_ENV = "CREVICE_EXAMPLE_BASE_URL"
RELEASE_052 = "https://github.com/ELambden/crevice/releases/download/v0.5.2/"
"""Release assets of the datasets added in 0.5.2 (AQP4 and OmpF excerpts)."""
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
    base_url : str or None
        URL prefix of the dataset's release assets; ``None`` means
        :data:`DEFAULT_BASE_URL` (the first datasets' release).
    """

    name: str
    description: str
    files: tuple[ExampleFile, ...]
    licence: str
    base_url: str | None = None


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
    "aqp4-excerpt": ExampleDataset(
        name="aqp4-excerpt",
        description=(
            "200 consecutive frames (30-50 ns, every 100 ps) of an all-atom membrane MD simulation of the "
            "human aquaporin-4 monomer: protein, lipid phosphorus atoms, ions and the 3,093 waters that come "
            "near the pore, CHARMM naming; the water permeation example"
        ),
        licence="CC-BY-4.0",
        files=(
            ExampleFile(
                "aqp4_excerpt.gro.gz",
                "87ab2188f2c50ebb15f0a300239a800e6c1ea0d2b046aef689ed08ca71bb66a5",
                167665,
                "topology and first-frame coordinates (GRO, gzip)",
                unpacked_name="aqp4_excerpt.gro",
                unpacked_sha256="49f7d1aa77fe4d5df8b23717c9d83683a96d5054c4557c46bddbe8a9610f32b1",
            ),
            ExampleFile(
                "aqp4_excerpt.xtc",
                "fb4bf2780d3bd7f014c3a432abf726c3c955038ae7bc438b4f7d3a061fbadfd8",
                9777740,
                "200-frame trajectory, 30-50 ns every 100 ps (XTC)",
            ),
            ExampleFile(
                "aqp4_excerpt_provenance.json",
                "a9f2eedb84a9b9f44049af7d4c7b5b71990454a2d3d69a5fb25b7ef70c998066",
                1431,
                "source hashes, frame indices, times and selections",
            ),
            ExampleFile(
                "aqp4_excerpt_README.md",
                "40f45cb799136e822cecfefcb76c9ab42c3e94f11b5808073a2214cdf27b4abc",
                3964,
                "data card: system, what was kept, protocol and licence",
            ),
        ),
        base_url=RELEASE_052,
    ),
    "ompf-excerpt": ExampleDataset(
        name="ompf-excerpt",
        description=(
            "100 consecutive frames (50-60 ns, every 100 ps) of an all-atom membrane MD simulation of the "
            "OmpF porin trimer in 1 M KCl under an applied field: protein, lipid phosphorus atoms and every "
            "ion, CHARMM naming, with the run parameters; the ion permeation and conductance example"
        ),
        licence="CC-BY-4.0",
        files=(
            ExampleFile(
                "ompf_excerpt.gro.gz",
                "2c44ea59d9111efc8077a0c1cb6d631552bbe950e6764805aadffadc6c9d988a",
                226958,
                "topology and first-frame coordinates (GRO, gzip)",
                unpacked_name="ompf_excerpt.gro",
                unpacked_sha256="64a0836442998e9e1e664605f47a521be550c8851cb1dc085530da5b4dc11c87",
            ),
            ExampleFile(
                "ompf_excerpt.xtc",
                "af65496fa3b513dc805c34ff7a5fcf900f67be17b5bfcbeba59bc637ba0674c2",
                7028820,
                "100-frame trajectory, 50-60 ns every 100 ps (XTC)",
            ),
            ExampleFile(
                "ompf_excerpt_production.mdp",
                "1b48a5df4a66a1235a87d19ae972977b3b150e1163df390130c3ef6e35f48093",
                1350,
                "production run parameters, with the applied field (GROMACS .mdp)",
            ),
            ExampleFile(
                "ompf_excerpt_provenance.json",
                "896c06df24a3b577d8b3e2b638642109ac774bd9b45437edac400b6912691765",
                1379,
                "source hashes, frame indices, times and selections",
            ),
            ExampleFile(
                "ompf_excerpt_README.md",
                "751546a7611f6eeb3df0dd2b032cf0003b65727474d6f917a5ab2e404b8570e7",
                3624,
                "data card: system, what was kept, protocol and licence",
            ),
        ),
        base_url=RELEASE_052,
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


def _release_tag(url: str) -> str | None:
    """Return the release tag of a GitHub ``releases/download/TAG/`` URL, if it is one."""
    match = re.search(r"/releases/download/([^/]+)/", url)
    return urllib.parse.unquote(match.group(1)) if match else None


def _retrieve_error_message(name: str, source: str, exc: OSError, local: bool) -> str:
    """One-line explanation of a failed example-asset copy or download."""
    if local:
        if isinstance(exc, FileNotFoundError):
            return (
                f"Example file {name} is not in the source directory ({source} does not exist); "
                "put the release assets there or drop --source-dir to download them"
            )
        return f"Could not copy example file {name} from {source}: {exc}"
    if isinstance(exc, urllib.error.HTTPError):
        if exc.code == 404:
            tag = _release_tag(source)
            where = f"the GitHub release {tag}" if tag else "the base URL"
            return (
                f"Example file {name} was not found at {source} (HTTP 404). The asset is probably not "
                f"attached to {where}; check --base-url (or ${BASE_URL_ENV}), or download the assets "
                "yourself and pass --source-dir DIR"
            )
        return (
            f"Could not download example file {name} from {source}: HTTP {exc.code} {exc.reason}; "
            "try again later or pass --source-dir DIR with a local copy"
        )
    reason = exc.reason if isinstance(exc, urllib.error.URLError) else exc
    return (
        f"Could not download example file {name} from {source}: {reason}. Check the network "
        "connection; when offline, pass --source-dir DIR with a local copy of the assets"
    )


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
        ``file://`` URLs work. Defaults to ``$CREVICE_EXAMPLE_BASE_URL``, then
        the dataset's own release (:attr:`ExampleDataset.base_url`), then
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
    RuntimeError
        If a file cannot be copied or downloaded. The message names the file,
        the source and the likely cause (for an HTTP 404, an asset missing from
        the GitHub release; for a network error, being offline and the
        ``source_dir`` alternative).
    """
    if name not in EXAMPLES:
        raise KeyError(f"Unknown example dataset {name!r}; available: {', '.join(sorted(EXAMPLES))}")
    dataset = EXAMPLES[name]
    directory = Path(cache_dir).expanduser() if cache_dir is not None else default_cache_dir()
    directory = (directory / name).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    local = source_dir if source_dir is not None else os.environ.get(SOURCE_DIR_ENV)
    prefix = base_url or os.environ.get(BASE_URL_ENV) or dataset.base_url or DEFAULT_BASE_URL
    if not prefix.endswith("/"):
        prefix += "/"
    result = FetchedExample(name=name, directory=directory)
    for item in dataset.files:
        target = directory / item.name
        if overwrite or not _matches(target, item.sha256):
            source = str(Path(local).expanduser() / item.name) if local else prefix + urllib.parse.quote(item.name)
            try:
                _retrieve(source, target, timeout)
            except OSError as exc:  # URLError and HTTPError are OSError subclasses
                raise RuntimeError(_retrieve_error_message(item.name, source, exc, local is not None)) from exc
            result.downloaded.append(item.name)
            actual = sha256_file(target)
            if actual != item.sha256:
                target.unlink()
                raise ValueError(
                    f"SHA-256 mismatch for {item.name} from {source}: expected {item.sha256}, got {actual}"
                )
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
                        raise ValueError(
                            f"SHA-256 mismatch for decompressed {item.unpacked_name}: "
                            f"expected {item.unpacked_sha256}, got {actual}"
                        )
                    os.replace(temporary, unpacked)
                    result.downloaded.append(item.unpacked_name)
                finally:
                    if os.path.exists(temporary):
                        os.unlink(temporary)
            result.paths[item.unpacked_name] = unpacked
    return result

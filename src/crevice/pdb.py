"""Fetching and resolving structures from RCSB or the local filesystem.

This module owns accession grammar, download transport, payload validation and
cache provenance for any PDB entry. :mod:`crevice.benchmarks` layers metadata
for a named subset on top of it. Nothing here validates the biological content
of an entry: a successful fetch means the bytes arrived intact and stated the
expected accession, not that the entry is a correct input for any question.

Workflow
    :func:`resolve_structure` turns whatever the user typed (a path or an
    accession such as ``1GRM``) into a local file. Accessions are served from
    a cache directory (default ``.crevice/pdb``) or downloaded over HTTPS
    from ``files.rcsb.org`` by :func:`fetch_structure`. Each download is
    size-limited, checked to be the requested format, compared with the
    accession the file states about itself, written atomically, and
    described by a ``<file>.provenance.json`` sidecar holding its SHA-256
    and retrieval time (UTC).

AlphaFold DB models
    Identifiers such as ``AF-P02920-F1`` (optionally ``-model_v6``, or the full
    download name ``AF-P02920-F1-model_v6.pdb``) are also recognised
    (:func:`parse_alphafold_id`). They are fetched from
    ``alphafold.ebi.ac.uk`` by :func:`fetch_alphafold_model` (the latest
    version, looked up through the AlphaFold DB API, unless one is given),
    cached under their AlphaFold DB file name with the same provenance
    sidecar, and checked against ``_entry.id`` (mmCIF) or the ``DBREF``
    UniProt accession (PDB). AlphaFold files are predictions: their B-factor
    column holds pLDDT confidence, and
    :func:`crevice.parser.load_structure_report` flags them as
    ``model_type="predicted"``.

The command-line equivalent is ``crevice fetch``; every command that takes a
structure argument also accepts a PDB accession or an AlphaFold DB identifier.

Examples
--------
>>> from crevice.pdb import is_accession, normalize_accession, structure_url
>>> is_accession("1grm"), is_accession("pdb_00001grm"), is_accession("GRM1")
(True, True, False)
>>> normalize_accession(" 1grm ")
'1GRM'
>>> structure_url("1GRM", fmt="pdb")
'https://files.rcsb.org/download/1GRM.pdb'
>>> structure_url("1af6", assembly=1)
'https://files.rcsb.org/download/1AF6-assembly1.cif'
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlsplit
from urllib.request import urlopen

RCSB_DOWNLOAD = "https://files.rcsb.org/download/{filename}"

#: Classic four-character PDB accession, e.g. 1GRM.
ACCESSION_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
#: Extended PDB accession, e.g. pdb_00001grm.
EXTENDED_ACCESSION_PATTERN = re.compile(r"^pdb_[0-9A-Za-z]{8}$")

#: Cap on a single download, large enough for deposited assemblies.
MAX_DOWNLOAD_BYTES = 128 * 1024 * 1024
DOWNLOAD_CHUNK_BYTES = 65536
PROVENANCE_SUFFIX = ".provenance.json"

#: Accession values that carry no entry identity, e.g. expanded assemblies.
PLACEHOLDER_ACCESSIONS = {"", ".", "?", "XXXX"}

STRUCTURE_FORMATS = ("pdb", "cif")

#: Grammar of an AlphaFold DB model identifier, ``AF-<UniProt>-F<fragment>``, optionally
#: followed by ``-model_v<version>`` and a ``.pdb``/``.cif`` suffix (the
#: AlphaFold DB download file name), case-insensitive.
ALPHAFOLD_PATTERN = re.compile(r"^AF-([0-9A-Z]{6,10})-F([1-9][0-9]*)(?:-MODEL_V([1-9][0-9]*))?(?:\.(PDB|CIF))?$")
#: UniProt accession grammar (as published by UniProt).
UNIPROT_PATTERN = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
ALPHAFOLD_DOWNLOAD = "https://alphafold.ebi.ac.uk/files/{filename}"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"


@dataclass(frozen=True)
class DownloadRecord:
    """What was retrieved, from where, and what was checked about it.

    Written to, and read back from, the ``.provenance.json`` sidecar next to a
    cached structure.

    Attributes
    ----------
    pdb_id : str
        Normalised accession (upper-case classic ID or lower-case extended ID).
    path : pathlib.Path
        Local file.
    url : str
        Source URL, or ``""`` for a cache entry this code did not download.
    sha256 : str
        Hex SHA-256 of the file bytes.
    byte_count : int
        File size in bytes.
    retrieved_utc : str
        ISO-8601 UTC timestamp (``YYYY-MM-DDTHH:MM:SSZ``) of the download, or
        the file modification time for a cache entry without provenance.
    content_format : str
        ``"pdb"`` or ``"cif"``.
    identity_verified : bool
        ``True`` only if the file states the requested accession. Assembly files
        carry a placeholder entry ID, so they are never verified.
    identity_note : str
        Human-readable reason for ``identity_verified``.
    source : str
        ``"download"`` or ``"existing_cache"``.
    assembly : int or None, default None
        Biological-assembly number, or ``None`` for the deposited unit.
    """

    pdb_id: str
    path: Path
    url: str
    sha256: str
    byte_count: int
    retrieved_utc: str
    content_format: str
    identity_verified: bool
    identity_note: str
    source: str
    assembly: int | None = None

    def to_dict(self) -> dict:
        """Return the fields as a JSON-friendly dictionary (``path`` as a string)."""
        data = asdict(self)
        data["path"] = str(self.path)
        return data


@dataclass(frozen=True)
class StructureSource:
    """Where a resolved structure came from, whatever the caller typed.

    Attributes
    ----------
    path : pathlib.Path
        Local file to read.
    spec : str
        The original argument (path or accession text).
    kind : str
        ``"local_file"``, ``"cached_accession"`` or ``"downloaded"``.
    accession : str or None, default None
        Normalised accession when ``spec`` was one.
    record : DownloadRecord or None, default None
        Provenance for cached or downloaded accessions.
    """

    path: Path
    spec: str
    kind: str
    accession: str | None = None
    record: DownloadRecord | None = None

    def to_dict(self) -> dict:
        """Return a JSON-friendly dictionary, with the download record nested under ``download``."""
        return {
            "path": str(self.path), "spec": self.spec, "kind": self.kind,
            "accession": self.accession,
            "download": self.record.to_dict() if self.record is not None else None,
        }


# --- accessions -------------------------------------------------------------

def is_accession(text: str) -> bool:
    """Return ``True`` when the text is a syntactically valid PDB accession.

    Accepted forms are the classic four-character ID (a digit followed by three
    letters or digits, any case) and the extended ``pdb_`` plus eight
    characters. Surrounding whitespace is ignored. Existence at RCSB is not
    checked.
    """

    candidate = str(text).strip()
    return bool(ACCESSION_PATTERN.match(candidate)
                or EXTENDED_ACCESSION_PATTERN.match(candidate.lower()))


def normalize_accession(text: str) -> str:
    """Return the canonical form of an accession, or raise for invalid input.

    Parameters
    ----------
    text : str

    Returns
    -------
    str
        Classic IDs in upper case (``"1GRM"``); extended IDs in lower case
        (``"pdb_00001grm"``).

    Raises
    ------
    ValueError
        If the text is not a valid accession.
    """

    candidate = str(text).strip()
    if EXTENDED_ACCESSION_PATTERN.match(candidate.lower()):
        return candidate.lower()
    if ACCESSION_PATTERN.match(candidate):
        return candidate.upper()
    raise ValueError(
        f"{text!r} is not a PDB accession; expected four characters starting "
        "with a digit, such as 1GRM, or an extended ID such as pdb_00001grm"
    )


def accession_identity(text: str) -> str:
    """Comparison key under which classic and extended forms of one entry agree.

    An extended ID whose first four characters after ``pdb_`` are ``0000`` is
    the extended form of a classic ID (``pdb_00001grm`` and ``1GRM`` are the
    same entry) and maps to the upper-case classic ID. Other extended IDs map
    to their lower-case form, classic IDs to upper case, and anything else to
    its stripped upper-case text.

    Examples
    --------
    >>> from crevice.pdb import accession_identity
    >>> accession_identity("pdb_00001grm"), accession_identity("PDB_00001GRM"), accession_identity("1grm")
    ('1GRM', '1GRM', '1GRM')
    >>> accession_identity("pdb_10001grm")
    'pdb_10001grm'
    """
    candidate = str(text).strip()
    lowered = candidate.lower()
    if EXTENDED_ACCESSION_PATTERN.match(lowered):
        body = lowered[4:]
        if body.startswith("0000") and ACCESSION_PATTERN.match(body[4:]):
            return body[4:].upper()
        return lowered
    return candidate.upper()


def normalize_format(extension: str) -> str:
    """Normalise a format name to ``"pdb"`` or ``"cif"``.

    Leading dots are removed and ``"mmcif"`` maps to ``"cif"``
    (case-insensitive).

    Parameters
    ----------
    extension : str
        Format name or file extension, for example ``".CIF"``, ``"mmcif"`` or
        ``"pdb"``.

    Returns
    -------
    str
        ``"pdb"`` or ``"cif"``.

    Raises
    ------
    ValueError
        For any other format.

    Examples
    --------
    >>> from crevice.pdb import normalize_format
    >>> normalize_format(".mmCIF"), normalize_format("PDB")
    ('cif', 'pdb')
    """
    fmt = str(extension).lower().lstrip(".")
    if fmt == "mmcif":
        fmt = "cif"
    if fmt not in STRUCTURE_FORMATS:
        raise ValueError("format must be pdb, cif, or mmcif")
    return fmt


def _normalize_assembly(assembly) -> int | None:
    """Validate an assembly number: ``None`` or a positive ``int`` (``bool`` rejected)."""
    if assembly is None:
        return None
    if isinstance(assembly, bool) or not isinstance(assembly, int) or assembly < 1:
        raise ValueError("assembly must be a positive integer, such as 1")
    return assembly


# --- locations --------------------------------------------------------------

def remote_filename(accession: str, fmt: str, assembly: int | None = None) -> str:
    """The filename RCSB serves for this entry, unit and format.

    The deposited unit is ``<id>.<fmt>``. Biological assemblies are
    ``<id>-assembly<n>.cif`` for mmCIF and ``<id>.pdb<n>`` for legacy PDB
    format. Inputs are not normalised here.

    Parameters
    ----------
    accession : str
        Normalised accession, for example ``"1GRM"``.
    fmt : {"pdb", "cif"}
        Normalised format.
    assembly : int, optional
        Biological assembly number; ``None`` = the deposited unit.

    Returns
    -------
    str
        Remote file name.

    Examples
    --------
    >>> from crevice.pdb import remote_filename
    >>> remote_filename("3UKM", "cif"), remote_filename("3UKM", "cif", 1), remote_filename("3UKM", "pdb", 1)
    ('3UKM.cif', '3UKM-assembly1.cif', '3UKM.pdb1')
    """

    if assembly is None:
        return f"{accession}.{fmt}"
    # RCSB serves mmCIF assemblies as <id>-assembly<n>.cif and legacy PDB
    # assemblies as <id>.pdb<n>.
    return f"{accession}-assembly{assembly}.cif" if fmt == "cif" else f"{accession}.pdb{assembly}"


def local_filename(accession: str, fmt: str, assembly: int | None = None) -> str:
    """The cache filename, kept uniform so the format is always the suffix.

    Parameters
    ----------
    accession : str
    fmt : str
    assembly : int, optional

    Returns
    -------
    str
        ``<id>.<fmt>`` or ``<id>-assembly<n>.<fmt>`` (inputs are not
        normalised here).
    """

    stem = accession if assembly is None else f"{accession}-assembly{assembly}"
    return f"{stem}.{fmt}"


def structure_url(accession: str, *, fmt: str = "cif", assembly: int | None = None) -> str:
    """Download URL for an entry at ``files.rcsb.org``.

    Parameters
    ----------
    accession : str
        Any valid accession; it is normalised.
    fmt : {"cif", "pdb", "mmcif"}, default "cif"
    assembly : int, optional
        Biological-assembly number (positive).

    Returns
    -------
    str

    Raises
    ------
    ValueError
        For an invalid accession, format or assembly number.
    """
    accession = normalize_accession(accession)
    fmt = normalize_format(fmt)
    return RCSB_DOWNLOAD.format(filename=remote_filename(accession, fmt, _normalize_assembly(assembly)))


def cached_structure_path(
    cache_dir: str | Path, accession: str, *, fmt: str = "cif", assembly: int | None = None
) -> Path:
    """Where a structure is (or would be) cached.

    Parameters
    ----------
    cache_dir : str or pathlib.Path
    accession : str
    fmt : {"cif", "pdb", "mmcif"}, default "cif"
    assembly : int, optional

    Returns
    -------
    pathlib.Path
        ``cache_dir / local_filename(...)``. The file need not exist.

    Raises
    ------
    ValueError
        For an invalid accession, format or assembly number.
    """
    accession = normalize_accession(accession)
    fmt = normalize_format(fmt)
    return Path(cache_dir) / local_filename(accession, fmt, _normalize_assembly(assembly))


def provenance_path(path: str | Path) -> Path:
    """Path of the provenance sidecar for a structure file.

    Parameters
    ----------
    path : str or pathlib.Path
        Structure file.

    Returns
    -------
    pathlib.Path
        ``<path>.provenance.json`` in the same directory.
    """
    structure = Path(path)
    return structure.with_name(structure.name + PROVENANCE_SUFFIX)


# --- fetching ---------------------------------------------------------------

def fetch_structure(
    accession: str,
    cache_dir: str | Path,
    *,
    fmt: str = "cif",
    assembly: int | None = None,
    overwrite: bool = False,
    timeout: float = 30.0,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
) -> DownloadRecord:
    """Fetch any RCSB entry into a local cache and record what was checked.

    If the destination already exists and ``overwrite`` is false, nothing is
    downloaded and :func:`cached_record` describes the existing file.
    Otherwise the file is downloaded over HTTPS, validated, written atomically
    (a partial download never replaces a good file) and a provenance sidecar is
    written.

    Parameters
    ----------
    accession : str
        Classic or extended PDB accession. An AlphaFold DB identifier is passed
        to :func:`fetch_alphafold_model` instead (``assembly`` must then be
        ``None``).
    cache_dir : str or pathlib.Path
        Directory for cached files (created if needed).
    fmt : {"cif", "pdb", "mmcif"}, default "cif"
        File format to request.
    assembly : int, optional
        Biological-assembly number; ``None`` fetches the deposited unit.
    overwrite : bool, default False
        Download even if a cached file exists.
    timeout : float, default 30.0
        Network timeout in seconds.
    max_bytes : int, default 128 MiB
        Maximum download size; larger payloads are rejected.

    Returns
    -------
    DownloadRecord

    Raises
    ------
    ValueError
        For an invalid accession, format, assembly or ``max_bytes``.
    RuntimeError
        If the download fails, exceeds ``max_bytes``, is empty, is markup
        rather than a coordinate file, lacks coordinate records, or (for the
        deposited unit) states a different accession.

    Notes
    -----
    The payload must parse as the requested format and, where the file states
    an accession, must match the requested one. Assembly files carry a
    placeholder entry ID, so they are cached with ``identity_verified`` false
    rather than assumed correct.

    Examples
    --------
    >>> from crevice.pdb import fetch_structure
    >>> record = fetch_structure("1GRM", ".crevice/pdb", fmt="pdb")  # doctest: +SKIP
    >>> record.identity_verified, record.path.name                   # doctest: +SKIP
    (True, '1GRM.pdb')
    """

    if is_alphafold_id(accession):
        if assembly is not None:
            raise ValueError("AlphaFold DB models have no biological assemblies; omit assembly")
        return fetch_alphafold_model(accession, cache_dir, fmt=fmt, overwrite=overwrite,
                                     timeout=timeout, max_bytes=max_bytes)
    accession = normalize_accession(accession)
    fmt = normalize_format(fmt)
    assembly = _normalize_assembly(assembly)
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes <= 0:
        raise ValueError("max_bytes must be a positive integer")
    destination = cached_structure_path(cache_dir, accession, fmt=fmt, assembly=assembly)
    if destination.exists() and not overwrite:
        return cached_record(accession, destination, fmt, assembly=assembly)

    url = structure_url(accession, fmt=fmt, assembly=assembly)
    require_https(url)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = _download(url, accession, timeout=timeout, max_bytes=max_bytes)
    text = validate_payload(payload, fmt, accession, url)
    verified, note = check_identity(text, accession, fmt, url, strict=assembly is None)
    _atomic_write(destination, payload)
    record = DownloadRecord(
        pdb_id=accession, path=destination, url=url,
        sha256=hashlib.sha256(payload).hexdigest(), byte_count=len(payload),
        retrieved_utc=_utc_now(), content_format=fmt,
        identity_verified=verified, identity_note=note, source="download",
        assembly=assembly,
    )
    write_provenance(record)
    return record


def resolve_structure(
    spec: str | Path,
    *,
    cache_dir: str | Path = ".crevice/pdb",
    fmt: str = "cif",
    assembly: int | None = None,
    offline: bool = False,
    overwrite: bool = False,
    timeout: float = 30.0,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
) -> StructureSource:
    """Resolve a user-supplied structure argument to a local file.

    An existing path always wins, so a local file is never shadowed by an
    accession that happens to share its name. Otherwise a syntactically valid
    PDB accession or AlphaFold DB identifier (:func:`parse_alphafold_id`) is
    served from the cache, or downloaded when not ``offline``.

    Parameters
    ----------
    spec : str or pathlib.Path
        File path, PDB accession or AlphaFold DB identifier.
    cache_dir : str or pathlib.Path, default ".crevice/pdb"
        Cache directory for accessions (relative to the working directory).
    fmt : {"cif", "pdb", "mmcif"}, default "cif"
        Format used for accessions.
    assembly : int, optional
        Biological-assembly number for PDB accessions (an error for AlphaFold DB
        identifiers).
    offline : bool, default False
        Never download; a missing cache entry is an error.
    overwrite : bool, default False
        Re-download a cached accession.
    timeout : float, default 30.0
        Network timeout in seconds.
    max_bytes : int, default 128 MiB
        Maximum download size.

    Returns
    -------
    StructureSource

    Raises
    ------
    ValueError
        If ``spec`` is a directory, or is neither an existing file nor a valid
        accession.
    RuntimeError
        If ``offline`` is set and the accession is not cached, or a download
        fails (see :func:`fetch_structure`).
    """

    text = str(spec)
    candidate = Path(text)
    if candidate.exists():
        if candidate.is_dir():
            raise ValueError(f"{candidate} is a directory, not a structure file")
        return StructureSource(path=candidate, spec=text, kind="local_file")
    model = parse_alphafold_id(text)
    if model is not None:
        if assembly is not None:
            raise ValueError("AlphaFold DB models have no biological assemblies; omit --assembly")
        af_fmt = normalize_format(model.fmt or fmt)
        cached = (Path(cache_dir) / model.filename(model.version, af_fmt) if model.version
                  else _cached_alphafold(cache_dir, model, af_fmt))
        if cached is not None and cached.exists() and not overwrite:
            return StructureSource(path=cached, spec=text, kind="cached_accession",
                                   accession=cached.stem, record=_cached_alphafold_record(model, cached, af_fmt))
        if offline:
            raise RuntimeError(
                f"{model.entry} is not in the cache at {Path(cache_dir)} and downloads are "
                "disabled; drop --offline or fetch it first"
            )
        record = fetch_alphafold_model(text, cache_dir, fmt=af_fmt, overwrite=overwrite,
                                       timeout=timeout, max_bytes=max_bytes)
        return StructureSource(path=record.path, spec=text, kind="downloaded",
                               accession=record.pdb_id, record=record)
    if not is_accession(text):
        raise ValueError(
            f"{text!r} is neither an existing file nor a PDB accession (such as 1GRM) or "
            "AlphaFold DB identifier (such as AF-P02920-F1); pass a path to a local structure or one of those IDs"
        )
    accession = normalize_accession(text)
    cached = cached_structure_path(cache_dir, accession, fmt=fmt, assembly=assembly)
    if cached.exists() and not overwrite:
        return StructureSource(
            path=cached, spec=text, kind="cached_accession", accession=accession,
            record=cached_record(accession, cached, normalize_format(fmt),
                                 assembly=_normalize_assembly(assembly)),
        )
    if offline:
        raise RuntimeError(
            f"{accession} is not in the cache at {Path(cache_dir)} and downloads are "
            "disabled; drop --offline or fetch it first"
        )
    record = fetch_structure(accession, cache_dir, fmt=fmt, assembly=assembly,
                             overwrite=overwrite, timeout=timeout, max_bytes=max_bytes)
    return StructureSource(path=record.path, spec=text, kind="downloaded",
                           accession=accession, record=record)


# --- AlphaFold DB -----------------------------------------------------------

@dataclass(frozen=True)
class AlphaFoldId:
    """A parsed AlphaFold DB model identifier.

    Attributes
    ----------
    uniprot : str
        UniProt accession, for example ``"P02920"``.
    fragment : int
        Fragment number (``F1`` for most proteins; long proteins are split).
    version : int or None
        Model version, or ``None`` for "the latest".
    fmt : str or None
        ``"pdb"`` or ``"cif"`` when the text carried a file suffix.
    """

    uniprot: str
    fragment: int
    version: int | None = None
    fmt: str | None = None

    @property
    def entry(self) -> str:
        """Entry ID without version, for example ``"AF-P02920-F1"``."""
        return f"AF-{self.uniprot}-F{self.fragment}"

    def filename(self, version: int, fmt: str) -> str:
        """AlphaFold DB file name, for example ``"AF-P02920-F1-model_v6.pdb"``.

        Parameters
        ----------
        version : int
            Model version.
        fmt : {"pdb", "cif"}
            File format (the suffix).

        Returns
        -------
        str
        """
        return f"{self.entry}-model_v{version}.{fmt}"


def parse_alphafold_id(text: str) -> AlphaFoldId | None:
    """Parse an AlphaFold DB identifier or download file name.

    Accepts ``AF-P02920-F1``, ``AF-P02920-F1-model_v6`` and
    ``AF-P02920-F1-model_v6.pdb`` (any case; a leading directory is ignored).

    Parameters
    ----------
    text : str

    Returns
    -------
    AlphaFoldId or None
        ``None`` when the text is not an AlphaFold identifier or the UniProt
        part is not a valid UniProt accession.

    Examples
    --------
    >>> from crevice.pdb import parse_alphafold_id
    >>> model = parse_alphafold_id("AF-P02920-F1-model_v6.pdb")
    >>> model.entry, model.version, model.fmt
    ('AF-P02920-F1', 6, 'pdb')
    >>> parse_alphafold_id("1GRM") is None
    True
    """
    match = ALPHAFOLD_PATTERN.match(Path(str(text).strip()).name.upper())
    if not match or not UNIPROT_PATTERN.match(match.group(1)):
        return None
    uniprot, fragment, version, suffix = match.groups()
    return AlphaFoldId(uniprot, int(fragment), int(version) if version else None,
                       suffix.lower() if suffix else None)


def is_alphafold_id(text: str) -> bool:
    """True when the text is an AlphaFold DB identifier (see :func:`parse_alphafold_id`)."""
    return parse_alphafold_id(text) is not None


def alphafold_latest_version(model: AlphaFoldId, *, timeout: float = 30.0) -> int:
    """Ask the AlphaFold DB API for the latest version of a model.

    Parameters
    ----------
    model : AlphaFoldId
    timeout : float, default 30.0
        Network timeout in seconds.

    Returns
    -------
    int

    Raises
    ------
    RuntimeError
        If the API cannot be reached, returns invalid data, or has no entry
        for this UniProt accession and fragment.
    """
    url = ALPHAFOLD_API.format(uniprot=model.uniprot)
    require_https(url)
    payload = _download(url, model.entry, timeout=timeout, max_bytes=4 * 1024 * 1024)
    try:
        entries = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"AlphaFold DB returned invalid JSON for {model.uniprot}") from exc
    for item in entries if isinstance(entries, list) else []:
        if isinstance(item, dict) and str(item.get("modelEntityId", "")).upper() == model.entry.upper():
            version = item.get("latestVersion")
            if isinstance(version, int) and version > 0:
                return version
    raise RuntimeError(f"AlphaFold DB has no model {model.entry} (UniProt {model.uniprot})")


def check_alphafold_identity(text: str, model: AlphaFoldId, fmt: str) -> tuple[bool, str]:
    """Check that a file states the requested AlphaFold model.

    mmCIF files must have ``_entry.id`` equal to the entry ID; PDB files must
    carry a ``DBREF`` record naming the UniProt accession.

    Parameters
    ----------
    text : str
        Decoded file content.
    model : AlphaFoldId
        The requested model.
    fmt : {"pdb", "cif"}
        File format of ``text``.

    Returns
    -------
    verified : bool
    note : str

    Raises
    ------
    RuntimeError
        If the file names a different entry or accession.
    """
    if fmt == "cif":
        stated, field = payload_accession(text, fmt)
        if stated is None or stated in PLACEHOLDER_ACCESSIONS:
            return False, f"no entry ID recorded in {field}"
        if stated != model.entry.upper():
            raise RuntimeError(f"file states entry {stated}, not the requested {model.entry}")
        return True, f"matched {model.entry} from _entry.id"
    for line in text.splitlines():
        if line.startswith("DBREF"):
            accession = line[33:41].strip().upper()
            if accession and accession != model.uniprot:
                raise RuntimeError(f"DBREF names UniProt {accession}, not the requested {model.uniprot}")
            if accession:
                return True, f"matched UniProt {model.uniprot} from DBREF"
    return False, "no DBREF record naming a UniProt accession"


def _cached_alphafold(cache_dir: str | Path, model: AlphaFoldId, fmt: str) -> Path | None:
    """Newest cached version of a model in ``fmt``, or ``None``."""
    versions = []
    for path in Path(cache_dir).glob(f"{model.entry}-model_v*.{fmt}"):
        parsed = parse_alphafold_id(path.name)
        if parsed is not None and parsed.entry == model.entry and parsed.version:
            versions.append((parsed.version, path))
    return max(versions)[1] if versions else None


def fetch_alphafold_model(
    identifier: str,
    cache_dir: str | Path,
    *,
    fmt: str = "cif",
    overwrite: bool = False,
    timeout: float = 30.0,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
) -> DownloadRecord:
    """Fetch an AlphaFold DB model into the local cache and record what was checked.

    The file is cached under its AlphaFold DB name
    (``AF-<UniProt>-F<n>-model_v<k>.<fmt>``). Without a version in
    ``identifier``, a cached copy of any version is reused (the newest), and
    otherwise the latest version is looked up through the AlphaFold DB API.
    Downloads are validated like RCSB files (format, HTTPS, size limit, atomic
    write, provenance sidecar) and checked against the requested entry.

    Parameters
    ----------
    identifier : str
        ``AF-P02920-F1``, ``AF-P02920-F1-model_v6`` or the download file name.
    cache_dir : str or pathlib.Path
    fmt : {"cif", "pdb", "mmcif"}, default "cif"
        Format to fetch; a ``.pdb``/``.cif`` suffix in ``identifier`` wins.
    overwrite : bool, default False
        Download even if a cached copy exists.
    timeout : float, default 30.0
        Network timeout in seconds.
    max_bytes : int, default 128 MiB

    Returns
    -------
    DownloadRecord
        ``pdb_id`` holds the versioned model ID, for example
        ``"AF-P02920-F1-model_v6"``.

    Raises
    ------
    ValueError
        For an invalid identifier, format or ``max_bytes``.
    RuntimeError
        If the API or download fails or the file names a different model.

    Notes
    -----
    AlphaFold models are predictions: the B-factor column holds pLDDT
    confidence, and low-confidence regions should not be read as structure.

    Examples
    --------
    >>> from crevice.pdb import fetch_alphafold_model
    >>> record = fetch_alphafold_model("AF-P02920-F1", ".crevice/pdb", fmt="pdb")  # doctest: +SKIP
    >>> record.path.name, record.identity_verified                                 # doctest: +SKIP
    ('AF-P02920-F1-model_v6.pdb', True)
    """
    model = parse_alphafold_id(identifier)
    if model is None:
        raise ValueError(f"{identifier!r} is not an AlphaFold DB identifier such as AF-P02920-F1")
    fmt = normalize_format(model.fmt or fmt)
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes <= 0:
        raise ValueError("max_bytes must be a positive integer")
    if model.version is None and not overwrite:
        cached = _cached_alphafold(cache_dir, model, fmt)
        if cached is not None:
            return _cached_alphafold_record(model, cached, fmt)
    version = model.version or alphafold_latest_version(model, timeout=timeout)
    destination = Path(cache_dir) / model.filename(version, fmt)
    if destination.exists() and not overwrite:
        return _cached_alphafold_record(model, destination, fmt)
    url = ALPHAFOLD_DOWNLOAD.format(filename=model.filename(version, fmt))
    require_https(url)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = _download(url, model.entry, timeout=timeout, max_bytes=max_bytes)
    text = validate_payload(payload, fmt, model.entry, url)
    verified, note = check_alphafold_identity(text, model, fmt)
    _atomic_write(destination, payload)
    record = DownloadRecord(
        pdb_id=destination.stem, path=destination, url=url,
        sha256=hashlib.sha256(payload).hexdigest(), byte_count=len(payload),
        retrieved_utc=_utc_now(), content_format=fmt,
        identity_verified=verified, identity_note=note, source="download",
    )
    write_provenance(record)
    return record


def _cached_alphafold_record(model: AlphaFoldId, destination: Path, fmt: str) -> DownloadRecord:
    """Provenance of a cached AlphaFold file (sidecar, or re-hashed and re-checked)."""
    recorded = read_provenance(destination)
    if recorded is not None:
        return recorded
    payload = destination.read_bytes()
    try:
        verified, note = check_alphafold_identity(payload.decode("utf-8", errors="replace"), model, fmt)
    except RuntimeError as exc:
        verified, note = False, str(exc)
    mtime = datetime.fromtimestamp(destination.stat().st_mtime, timezone.utc)
    return DownloadRecord(
        pdb_id=destination.stem, path=destination, url="",
        sha256=hashlib.sha256(payload).hexdigest(), byte_count=len(payload),
        retrieved_utc=mtime.strftime("%Y-%m-%dT%H:%M:%SZ"), content_format=fmt,
        identity_verified=verified, identity_note=note, source="existing_cache",
    )


# --- transport --------------------------------------------------------------

def _download(url: str, accession: str, *, timeout: float, max_bytes: int) -> bytes:
    """Stream a URL in chunks, enforcing ``max_bytes``; wrap network errors in ``RuntimeError``."""
    chunks: list[bytes] = []
    total = 0
    try:
        with urlopen(url, timeout=timeout) as response:
            while True:
                chunk = response.read(min(DOWNLOAD_CHUNK_BYTES, max_bytes))
                if not chunk:
                    break
                total += len(chunk)
                if total > max_bytes:
                    raise RuntimeError(
                        f"{url} exceeds the {max_bytes} byte limit; "
                        "raise max_bytes to accept it"
                    )
                chunks.append(chunk)
    except (URLError, OSError) as exc:
        raise RuntimeError(f"Could not fetch {accession} from {url}: {exc}") from exc
    return b"".join(chunks)


def validate_payload(payload: bytes, fmt: str, accession: str, url: str) -> str:
    """Reject anything that is not the requested coordinate format.

    A PDB payload must contain an ``ATOM``/``HETATM`` record; an mmCIF payload
    must contain a ``data_`` block and an ``_atom_site.`` table. Empty payloads
    and markup (HTML/XML error pages) are rejected.

    Parameters
    ----------
    payload : bytes
        Downloaded bytes.
    fmt : {"pdb", "cif"}
        Requested format.
    accession : str
        Requested accession (for the error message).
    url : str
        Source URL (for the error message).

    Returns
    -------
    str
        The payload decoded as UTF-8 (undecodable bytes replaced).

    Raises
    ------
    RuntimeError
        If the payload fails a check.
    """

    if not payload.strip():
        raise RuntimeError(f"RCSB returned an empty payload for {accession}")
    text = payload.decode("utf-8", errors="replace")
    if text.lstrip()[:1] == "<":
        raise RuntimeError(f"{url} returned markup, not a {fmt} coordinate file")
    lines = text.splitlines()
    if fmt == "pdb":
        if not any(line.startswith(("ATOM  ", "HETATM")) for line in lines):
            raise RuntimeError(f"{url} contains no ATOM or HETATM coordinate records")
    else:
        if not any(line.startswith("data_") for line in lines):
            raise RuntimeError(f"{url} contains no mmCIF data block")
        if "_atom_site." not in text:
            raise RuntimeError(f"{url} contains no _atom_site coordinate table")
    return text


def check_identity(text: str, accession: str, fmt: str, source: str, *, strict: bool) -> tuple[bool, str]:
    """Compare the requested accession with the one the file states.

    Parameters
    ----------
    text : str
        Decoded file content.
    accession : str
        Requested accession.
    fmt : {"pdb", "cif"}
    source : str
        URL or path, used in messages.
    strict : bool
        Raise on a mismatch instead of returning ``False``.

    Returns
    -------
    verified : bool
        ``True`` only when the stated accession names the requested entry.
        Case is ignored, and a classic ID matches its extended form
        (``1GRM`` and ``pdb_00001grm``; see :func:`accession_identity`). A
        missing or placeholder accession gives ``False``.
    note : str
        Explanation.

    Raises
    ------
    RuntimeError
        On a mismatch when ``strict`` is true.
    """

    stated, field = payload_accession(text, fmt)
    if stated is None:
        return False, f"no accession recorded: {field}"
    if stated in PLACEHOLDER_ACCESSIONS:
        return False, f"no accession recorded: {field} holds the placeholder {stated}"
    if accession_identity(stated) != accession_identity(accession):
        message = f"{source} reports accession {stated}, not the requested {accession}"
        if strict:
            raise RuntimeError(message)
        return False, message
    return True, f"matched {accession} from {field}"


def payload_accession(text: str, fmt: str) -> tuple[str | None, str]:
    r"""Read the accession the file states about itself, if it states one.

    PDB files use columns 63-66 of the first ``HEADER`` record; mmCIF files use
    the ``_entry.id`` item.

    Parameters
    ----------
    text : str
        Decoded file contents.
    fmt : {"pdb", "cif"}
        Format of ``text``.

    Returns
    -------
    accession : str or None
        Upper-case accession, ``""`` for an empty field, or ``None`` when the
        record/item is absent.
    field : str
        Description of where the value was (or was not) found.

    Examples
    --------
    >>> from crevice.pdb import payload_accession
    >>> payload_accession("data_1GRM\n_entry.id 1grm\n", "cif")
    ('1GRM', '_entry.id')
    """

    if fmt == "pdb":
        for line in text.splitlines():
            if line.startswith("HEADER"):
                return line[62:66].strip().upper() or "", "the HEADER record, columns 63-66"
        return None, "the file has no HEADER record"
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("_entry.id"):
            parts = stripped.split(None, 1)
            value = parts[1].strip().strip("'\"").upper() if len(parts) > 1 else ""
            return value, "_entry.id"
    return None, "the file has no _entry.id value"


# --- cache provenance -------------------------------------------------------

def cached_record(
    accession: str, destination: Path, fmt: str, *, assembly: int | None = None
) -> DownloadRecord:
    """Describe a cache entry that this call did not download.

    Returns the recorded provenance if a sidecar exists. Otherwise the file is
    hashed, its stated accession checked non-strictly, and the modification time
    used as ``retrieved_utc`` with ``source="existing_cache"``.

    Parameters
    ----------
    accession : str
        Requested accession.
    destination : Path
        Cached structure file.
    fmt : {"pdb", "cif"}
        Its format.
    assembly : int, optional
        Biological assembly number, recorded in the result.

    Returns
    -------
    DownloadRecord
        Provenance of the cached file (SHA-256, size, retrieval time, identity
        check result and note).
    """

    recorded = read_provenance(destination)
    if recorded is not None:
        return recorded
    payload = destination.read_bytes()
    text = payload.decode("utf-8", errors="replace")
    verified, note = check_identity(text, accession, fmt, str(destination), strict=False)
    mtime = datetime.fromtimestamp(destination.stat().st_mtime, timezone.utc)
    return DownloadRecord(
        pdb_id=normalize_accession(accession), path=destination, url="",
        sha256=hashlib.sha256(payload).hexdigest(), byte_count=len(payload),
        retrieved_utc=mtime.strftime("%Y-%m-%dT%H:%M:%SZ"), content_format=fmt,
        identity_verified=verified, identity_note=note, source="existing_cache",
        assembly=assembly,
    )


def read_provenance(path: str | Path) -> DownloadRecord | None:
    """Return the recorded provenance for a cached structure, if any.

    Parameters
    ----------
    path : str or pathlib.Path
        Structure file (not the sidecar).

    Returns
    -------
    DownloadRecord or None
        ``None`` when there is no sidecar.
    """

    sidecar = provenance_path(path)
    if not sidecar.exists():
        return None
    data = json.loads(sidecar.read_text())
    return DownloadRecord(**{**data, "path": Path(data["path"])})


def write_provenance(record: DownloadRecord) -> None:
    """Write a record's JSON sidecar next to its structure file.

    Parameters
    ----------
    record : DownloadRecord
        Written to :func:`provenance_path` of ``record.path`` (sorted keys,
        indented, overwriting any existing sidecar).
    """
    sidecar = provenance_path(record.path)
    sidecar.write_text(json.dumps(record.to_dict(), indent=2, sort_keys=True) + "\n")


def verify_cached_structure(path: str | Path) -> DownloadRecord:
    """Re-hash a cached structure and compare it with its recorded provenance.

    Parameters
    ----------
    path : str or Path
        Cached structure with a ``.provenance.json`` sidecar.

    Returns
    -------
    DownloadRecord
        The recorded provenance, when the SHA-256 matches.

    Raises
    ------
    RuntimeError
        If there is no provenance sidecar or the hash differs.
    """

    structure = Path(path)
    record = read_provenance(structure)
    if record is None:
        raise RuntimeError(f"{structure} has no recorded provenance to verify against")
    digest = hashlib.sha256(structure.read_bytes()).hexdigest()
    if digest != record.sha256:
        raise RuntimeError(
            f"{structure} does not match its recorded SHA-256 "
            f"({digest}, recorded {record.sha256})"
        )
    return record


# --- helpers ----------------------------------------------------------------

def _atomic_write(destination: Path, payload: bytes) -> None:
    """Replace the destination only once the whole payload is on disk."""

    handle = tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=destination.name + ".", suffix=".part", delete=False
    )
    temporary = Path(handle.name)
    try:
        with handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def require_https(url: str) -> None:
    """Refuse any download URL that is not HTTPS.

    Parameters
    ----------
    url : str

    Raises
    ------
    ValueError
        If the URL scheme is not ``https``.
    """
    scheme = urlsplit(url).scheme
    if scheme != "https":
        raise ValueError(f"Structure downloads require https, not {scheme!r}: {url}")


def _utc_now() -> str:
    """Current UTC time as ``YYYY-MM-DDTHH:MM:SSZ``."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

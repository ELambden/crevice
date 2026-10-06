"""Named benchmark structures and their curation status.

This module registers a small set of RCSB entries used for end-to-end and
scalability checks (``crevice benchmark``, ``crevice static-suite``) and
restricts benchmark downloads to that set. Fetching is delegated to
:mod:`crevice.pdb`.

Curation status
    Every entry carries a ``curation_status``. An entry is ``"curated"``
    only once its assembly, atom selection, channel axis and expected
    regions have been fixed and checked. **No entry is curated yet**
    (:func:`curated_benchmark_ids` returns an empty tuple), so any cavity
    count, volume or radius computed on these structures is a provisional
    geometric measurement, not a validated result.

Examples
--------
>>> from crevice.benchmarks import benchmark_systems, curated_benchmark_ids
>>> [system.pdb_id for system in benchmark_systems()]
['1GRM', '2CHB', '6MVY', '1AF6', '4PYP']
>>> curated_benchmark_ids()
()
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .pdb import (
    MAX_DOWNLOAD_BYTES,
    DownloadRecord,
    cached_structure_path,
    fetch_structure,
    provenance_path,
    read_provenance,
    verify_cached_structure,
)


@dataclass(frozen=True)
class BenchmarkSystem:
    """A named entry in the intended static benchmark set.

    ``expected_geometry`` selects a workflow cast mode; it is not a validated
    geometric result for the entry. Values other than ``channel`` and ``cavity``
    fall back to automatic mode selection, and ``unassigned`` marks a label that
    was withdrawn rather than replaced. ``curation_status`` records whether the
    assembly, atom selection, axis and expected regions have been curated:
    nothing in this set is curated yet, so every result from it is provisional.

    Attributes
    ----------
    pdb_id : str
        RCSB accession.
    name : str
        Short descriptive name of the entry.
    expected_geometry : str
        ``"channel"``, ``"cavity"`` or ``"unassigned"`` (see above).
    notes : str
        Caveats about identity, assembly and scope.
    curation_status : str, default "uncurated"
        One of :data:`CURATION_STATUSES`: ``"curated"``, ``"uncurated"`` or
        ``"set_aside"``.
    identity_source : str, default ""
        URL of the RCSB entry page used to check the entry's identity.
    """

    pdb_id: str
    name: str
    expected_geometry: str
    notes: str
    curation_status: str = "uncurated"
    identity_source: str = ""

    @property
    def normalized_id(self) -> str:
        """The accession in upper case."""
        return self.pdb_id.upper()


RCSB_ENTRY = "https://www.rcsb.org/structure/"

#: Registered benchmark systems keyed by upper-case accession.
BENCHMARK_SYSTEMS: dict[str, BenchmarkSystem] = {
    "1GRM": BenchmarkSystem(
        pdb_id="1GRM",
        name="gramicidin channel system",
        expected_geometry="channel",
        notes="Small membrane-channel benchmark for fast end-to-end checks. Its "
              "modified polymer residues are deposited as HETATM and must be "
              "retained by entity membership, not dropped as ligands.",
        identity_source=RCSB_ENTRY + "1GRM",
    ),
    "2CHB": BenchmarkSystem(
        pdb_id="2CHB",
        name="cholera toxin B-pentamer complex",
        expected_geometry="unassigned",
        notes="Set aside pending a scope decision. The earlier 'channel benchmark "
              "system' label misidentified this entry; no pore geometry is claimed "
              "for it here and it is not part of any curated benchmark set.",
        curation_status="set_aside",
        identity_source=RCSB_ENTRY + "2CHB",
    ),
    "1AF6": BenchmarkSystem(
        pdb_id="1AF6",
        name="Escherichia coli maltoporin–sucrose trimer",
        expected_geometry="channel",
        notes="Identity and deposited trimeric assembly checked against RCSB on "
              "14 September 2026. Per-protomer channels and the central interfacial "
              "space still require separate geometric/accessibility validation.",
        identity_source=RCSB_ENTRY + "1AF6",
    ),
    "4PYP": BenchmarkSystem(
        pdb_id="4PYP",
        name="human glucose transporter GLUT1",
        expected_geometry="unassigned",
        notes="Corrected from an earlier 'photoactive yellow protein' label, which "
              "was a misidentification. Its inward-open annotation is not "
              "equivalent to a buried, disconnected cavity, so the previous "
              "cavity mode is withdrawn rather than reassigned.",
        identity_source=RCSB_ENTRY + "4PYP",
    ),
    "6MVY": BenchmarkSystem(
        pdb_id="6MVY",
        name="large membrane-channel benchmark system",
        expected_geometry="channel",
        notes="Large, complex benchmark for scalability and robust exports. The "
              "deposited two-chain coordinates are not the biological tetramer; an "
              "assembly must be chosen explicitly before any geometric claim.",
        identity_source=RCSB_ENTRY + "6MVY",
    ),
}

#: Allowed values of :attr:`BenchmarkSystem.curation_status`.
CURATION_STATUSES = ("curated", "uncurated", "set_aside")

#: Declared order of the benchmark set, used by every benchmark command.
DEFAULT_BENCHMARK_IDS = ("1GRM", "2CHB", "6MVY", "1AF6", "4PYP")


def benchmark_systems() -> tuple[BenchmarkSystem, ...]:
    """All registered benchmark systems in the declared order.

    Returns
    -------
    tuple of BenchmarkSystem
        Ordered as :data:`DEFAULT_BENCHMARK_IDS`.
    """
    return tuple(BENCHMARK_SYSTEMS[pdb_id] for pdb_id in DEFAULT_BENCHMARK_IDS)


def benchmark_ids_with_status(status: str) -> tuple[str, ...]:
    """IDs whose curation status matches, in the declared benchmark order.

    Parameters
    ----------
    status : {"curated", "uncurated", "set_aside"}

    Returns
    -------
    tuple of str

    Raises
    ------
    ValueError
        For an unknown status.
    """

    if status not in CURATION_STATUSES:
        raise ValueError(f"status must be one of {', '.join(CURATION_STATUSES)}")
    return tuple(
        pdb_id for pdb_id in DEFAULT_BENCHMARK_IDS
        if BENCHMARK_SYSTEMS[pdb_id].curation_status == status
    )


def curated_benchmark_ids() -> tuple[str, ...]:
    """IDs whose assembly, selection, axis and expected regions are curated.

    This is currently empty: no entry in the set has been curated, so no
    benchmark result from it should be read as a validated geometric outcome.
    """

    return benchmark_ids_with_status("curated")


def get_benchmark_system(pdb_id: str) -> BenchmarkSystem:
    """Look up a registered benchmark system by accession (case-insensitive).

    Parameters
    ----------
    pdb_id : str

    Returns
    -------
    BenchmarkSystem

    Raises
    ------
    ValueError
        If the accession is not registered. Use
        :func:`crevice.pdb.fetch_structure` for any other entry.
    """
    key = str(pdb_id).upper()
    try:
        return BENCHMARK_SYSTEMS[key]
    except KeyError as exc:
        valid = ", ".join(DEFAULT_BENCHMARK_IDS)
        raise ValueError(
            f"Unknown benchmark system {pdb_id!r}; registered IDs are {valid}. "
            "Any other accession can be fetched with crevice.pdb.fetch_structure."
        ) from exc


def benchmark_path(
    cache_dir: str | Path, pdb_id: str, *, extension: str = "cif", assembly: int | None = None
) -> Path:
    """Cache path of a registered benchmark structure.

    Parameters
    ----------
    cache_dir : str or pathlib.Path
    pdb_id : str
        Registered accession.
    extension : {"cif", "pdb", "mmcif"}, default "cif"
    assembly : int, optional

    Returns
    -------
    pathlib.Path
        The file need not exist.

    Raises
    ------
    ValueError
        For an unregistered accession or invalid format/assembly.
    """
    system = get_benchmark_system(pdb_id)
    return cached_structure_path(cache_dir, system.normalized_id, fmt=extension, assembly=assembly)


def fetch_benchmark_structure(
    pdb_id: str,
    cache_dir: str | Path,
    *,
    extension: str = "cif",
    overwrite: bool = False,
    timeout: float = 30.0,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
    assembly: int | None = None,
) -> Path:
    """Fetch a named benchmark structure from RCSB into a local cache.

    Same as :func:`fetch_benchmark_record` but returns only the local path.

    Parameters
    ----------
    pdb_id : str
        Registered benchmark accession (see :func:`benchmark_systems`).
    cache_dir : str or pathlib.Path
        Cache directory (created); a valid cached file is reused.
    extension : {"cif", "pdb", "mmcif"}, default "cif"
        File format.
    overwrite : bool, default False
        Download again even when a cached file exists.
    timeout : float, default 30.0
        Network timeout, seconds.
    max_bytes : int, default 128 MiB
        Largest accepted download.
    assembly : int, optional
        Biological assembly number; default the deposited unit.

    Returns
    -------
    pathlib.Path
        The cached structure file.

    Raises
    ------
    ValueError, RuntimeError
        As for :func:`fetch_benchmark_record`.
    """

    return fetch_benchmark_record(
        pdb_id, cache_dir, extension=extension, overwrite=overwrite,
        timeout=timeout, max_bytes=max_bytes, assembly=assembly,
    ).path


def fetch_benchmark_record(
    pdb_id: str,
    cache_dir: str | Path,
    *,
    extension: str = "cif",
    overwrite: bool = False,
    timeout: float = 30.0,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
    assembly: int | None = None,
) -> DownloadRecord:
    """Fetch a named benchmark structure, restricted to the registered systems.

    Use :func:`crevice.pdb.fetch_structure` for any other accession; this
    wrapper exists so benchmark commands cannot silently analyse an unregistered
    entry.

    Parameters
    ----------
    pdb_id : str
        Registered accession.
    cache_dir : str or pathlib.Path
        Cache directory (created); a valid cached file is reused.
    extension : {"cif", "pdb", "mmcif"}, default "cif"
        File format; the same default as :func:`crevice.pdb.fetch_structure`
        and the ``--format`` option of ``crevice fetch``/``benchmark``.
    overwrite : bool, default False
        Download again even when a cached file exists.
    timeout : float, default 30.0
        Network timeout, seconds.
    max_bytes : int, default 128 MiB
        Largest accepted download.
    assembly : int, optional
        Biological assembly number; default the deposited unit.

    Returns
    -------
    DownloadRecord
        Path, URL, SHA-256, size, retrieval time and identity check.

    Raises
    ------
    ValueError, RuntimeError
        As for :func:`get_benchmark_system` and
        :func:`crevice.pdb.fetch_structure`.
    """

    system = get_benchmark_system(pdb_id)
    return fetch_structure(
        system.normalized_id, cache_dir, fmt=extension, assembly=assembly,
        overwrite=overwrite, timeout=timeout, max_bytes=max_bytes,
    )


def fetch_all_benchmark_structures(
    cache_dir: str | Path,
    *,
    extension: str = "cif",
    overwrite: bool = False,
    timeout: float = 30.0,
    max_bytes: int = MAX_DOWNLOAD_BYTES,
) -> tuple[Path, ...]:
    """Fetch every registered benchmark structure (deposited units only).

    Parameters
    ----------
    cache_dir : str or pathlib.Path
    extension : {"cif", "pdb", "mmcif"}, default "cif"
    overwrite : bool, default False
    timeout : float, default 30.0
    max_bytes : int, default 128 MiB

    Returns
    -------
    tuple of pathlib.Path
        Local paths in :data:`DEFAULT_BENCHMARK_IDS` order.
    """
    return tuple(
        fetch_benchmark_structure(
            pdb_id,
            cache_dir,
            extension=extension,
            overwrite=overwrite,
            timeout=timeout,
            max_bytes=max_bytes,
        )
        for pdb_id in DEFAULT_BENCHMARK_IDS
    )


__all__ = [
    "BENCHMARK_SYSTEMS", "CURATION_STATUSES", "DEFAULT_BENCHMARK_IDS",
    "BenchmarkSystem", "DownloadRecord", "benchmark_ids_with_status",
    "benchmark_path", "benchmark_systems", "curated_benchmark_ids",
    "fetch_all_benchmark_structures", "fetch_benchmark_record",
    "fetch_benchmark_structure", "get_benchmark_system", "provenance_path",
    "read_provenance", "verify_cached_structure",
]

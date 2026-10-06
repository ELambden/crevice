"""Explicit, entity-aware mmCIF intake using the optional Gemmi parser.

Install with ``pip install 'crevice[structures]'``. Compared with the
built-in reader in :mod:`crevice.parser`, this loader:

* validates the whole ``_atom_site`` table (unique integer IDs, supported
  record types, finite coordinates in *every* model);
* uses ``_entity`` types so that modified polymer residues deposited as
  ``HETATM`` are kept as polymer, not dropped as ligands;
* drops zero-occupancy atoms (the built-in reader keeps them);
* rejects an atom with both blank and labelled alternate locations;
* returns a report of every policy applied and every excluded atom class.

Both readers use label chain IDs by default (which keep assembly copies and
ligands distinct) and the same deterministic alternate-location rule: blank
plus the first lexical non-blank location per residue.

Main entry point: :func:`load_mmcif` (used automatically by
:func:`crevice.parser.load_structure` for ``.cif`` files when Gemmi is
installed).

Notes
-----
Values are read with ``block.get_mmcif_category``, which returns ``None``
for ``?`` and ``False`` for ``.``; both are treated as missing. The loader
checks file integrity, not biological completeness: unobserved residues are
counted in the report but not modelled.
"""

from __future__ import annotations

from collections import Counter
import math
import operator
from pathlib import Path

from .models import Atom, StructureFrame
from .radii import VDW_RADII

#: Values RCSB writes when an entry has no accession, e.g. expanded assemblies.
PLACEHOLDER_ENTRY_IDS = {"", ".", "?", "XXXX"}

REQUIRED_ATOM_SITE_COLUMNS = frozenset({
    "id", "group_PDB", "type_symbol", "label_atom_id", "label_comp_id",
    "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
})


def gemmi_available() -> bool:
    """True when the optional structures extra (Gemmi) is installed."""

    try:
        import gemmi  # noqa: F401
    except ImportError:
        return False
    return True


def read_mmcif(text: str):
    """Parse mmCIF text and return its single data block.

    Parameters
    ----------
    text : str
        Complete mmCIF file content.

    Returns
    -------
    gemmi.cif.Block

    Raises
    ------
    RuntimeError
        If Gemmi is not installed.
    ValueError
        If the text does not parse or does not contain exactly one data block.
    """
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("Validated mmCIF intake requires pip install 'crevice[structures]'") from exc
    try:
        return gemmi.cif.read_string(text).sole_block()
    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"Invalid single-block mmCIF: {exc}") from exc


def category_rows(block, category: str) -> list[dict]:
    """Rows of an mmCIF category as dictionaries.

    Parameters
    ----------
    block : gemmi.cif.Block
    category : str
        Category name including the trailing dot, for example ``"_entity."``.

    Returns
    -------
    list of dict
        One dictionary per row, keyed by item name without the category prefix
        (``"id"``, ``"type"``, ...). Empty if the category is absent.
    """
    columns = block.get_mmcif_category(category)
    return [dict(zip(columns, row)) for row in zip(*columns.values())] if columns else []


def atom_site_rows(block) -> list[dict]:
    """Return atom-site rows after checking table-wide structural integrity.

    Coordinates and identifiers are validated across every model, so a defect
    outside the requested model is reported rather than silently skipped.

    Parameters
    ----------
    block : gemmi.cif.Block

    Returns
    -------
    list of dict
        Every ``_atom_site`` row.

    Raises
    ------
    ValueError
        If required columns (``id``, ``group_PDB``, ``type_symbol``,
        ``label_atom_id``, ``label_comp_id``, ``label_asym_id``, ``Cartn_x/y/z``)
        are missing, an ID is missing, non-integral or duplicated, a record type
        is not ``ATOM``/``HETATM``, or a coordinate is missing, non-numeric or
        non-finite.
    """
    rows = category_rows(block, "_atom_site.")
    if not rows or not REQUIRED_ATOM_SITE_COLUMNS.issubset(rows[0]):
        raise ValueError("mmCIF requires a complete atom-site coordinate table")
    identifiers = set()
    for row in rows:
        identifier = atom_site_id(row)
        if identifier in identifiers:
            raise ValueError(f"Missing or duplicate mmCIF atom-site ID: {identifier}")
        identifiers.add(identifier)
        if _value(row, "group_PDB") not in {"ATOM", "HETATM"}:
            raise ValueError(f"Unsupported mmCIF atom-site group_PDB: {_value(row, 'group_PDB')!r}")
        for key in ("Cartn_x", "Cartn_y", "Cartn_z"):
            coordinate(row, key)
    return rows


def load_mmcif(
    path: str | Path, *, model_index: int = 0,
    selection: str = "polymer-heavy", chain_ids: str = "label",
    hetero_policy: str = "group_pdb",
) -> tuple[StructureFrame, dict]:
    """Load one model, retaining polymer HETATM via entity membership.

    Label chain IDs preserve assembly copies and separate ligands from polymers.
    Residue numbers remain author numbers, with label numbers only as a fallback.
    Blank alternate locations plus the first lexical non-blank location per
    residue are retained. This is a deterministic policy, not conformer
    inference.

    Parameters
    ----------
    path : str or pathlib.Path
        mmCIF file (single data block, uncompressed).
    model_index : int, default 0
        Zero-based index into the distinct ``pdbx_PDB_model_num`` values in
        order of appearance (a missing model number counts as ``"1"``).
    selection : {"polymer-heavy", "all-heavy", "all"}, default "polymer-heavy"
        ``"polymer-heavy"``: atoms of polymer entities without H/D (requires
        entity types for every atom). ``"all-heavy"``: all entities without
        H/D. ``"all"``: every atom, including hydrogens.
    chain_ids : {"label", "auth"}, default "label"
        Which chain-ID column becomes ``Atom.chain_id``.
    hetero_policy : {"group_pdb", "nonpolymer", "auto"}, default "group_pdb"
        What ``Atom.hetero`` means. ``"group_pdb"`` mirrors the deposited record
        type. ``"nonpolymer"`` uses entity membership instead, so a modified
        polymer residue deposited as HETATM is not treated as a ligand; it
        requires entity types. ``"auto"`` uses ``"nonpolymer"`` when entity
        types are available and falls back to ``"group_pdb"``. The report
        records which applied.

    Returns
    -------
    frame : StructureFrame
    report : dict
        Parser and version, entry and model IDs (with a flag for placeholder
        entry IDs in assembly files), the selection and policies, counts of
        excluded atoms by reason (``alternate_location``, ``zero_occupancy``,
        ``nonpolymer``, ``hydrogen_or_deuterium``), per-chain summaries,
        polymer HETATM retained, elements without a tabulated radius, the
        coordinate bounds (Å), assembly definitions and counts of unobserved
        residue/atom records. ``completeness_validated`` is always ``False``.

    Raises
    ------
    RuntimeError
        If Gemmi is not installed.
    ValueError
        For invalid options, an out-of-range model, integrity failures (see
        :func:`atom_site_rows`), an atom with both blank and labelled alternate
        locations, negative occupancy, a missing element or chain ID, duplicate
        ``(residue, atom name)`` identities (try ``chain_ids="label"``), or an
        empty selection.

    Examples
    --------
    >>> from crevice.mmcif import load_mmcif
    >>> frame, report = load_mmcif("1grm.cif", hetero_policy="auto")  # doctest: +SKIP
    >>> report["hetero_definition"], report["polymer_hetatm_retained"] > 0  # doctest: +SKIP
    ('nonpolymer', True)
    """
    model_index = _model_index(model_index)
    if selection not in {"polymer-heavy", "all-heavy", "all"}:
        raise ValueError("selection must be polymer-heavy, all-heavy, or all")
    if chain_ids not in {"label", "auth"}:
        raise ValueError("chain_ids must be label or auth")
    if hetero_policy not in {"group_pdb", "nonpolymer", "auto"}:
        raise ValueError("hetero_policy must be group_pdb, nonpolymer, or auto")
    source = Path(path)
    block = read_mmcif(source.read_text())
    rows = atom_site_rows(block)
    model_ids = list(dict.fromkeys(_value(row, "pdbx_PDB_model_num") or "1" for row in rows))
    if model_index >= len(model_ids):
        raise ValueError(f"model_index {model_index} out of range for {len(model_ids)} models")
    rows = [row for row in rows if (_value(row, "pdbx_PDB_model_num") or "1") == model_ids[model_index]]
    entity_types = {_value(row, "id"): _value(row, "type") for row in category_rows(block, "_entity.")}
    entity_types_available = bool(entity_types)
    if hetero_policy == "auto":
        hetero_definition = "nonpolymer" if entity_types_available else "group_pdb"
    elif hetero_policy == "nonpolymer" and not entity_types_available:
        raise ValueError("hetero_policy 'nonpolymer' requires _entity types for every atom")
    else:
        hetero_definition = hetero_policy
    if selection == "polymer-heavy" and any(_value(row, "label_entity_id") not in entity_types for row in rows):
        raise ValueError("polymer-heavy selection requires entity type for every atom")

    alternatives = _residue_alternatives(rows)
    chosen = {key: min(values) for key, values in alternatives.items()}
    atoms = []
    seen = set()
    chains: dict[str, dict] = {}
    excluded = Counter()
    polymer_hetero = 0
    for row in rows:
        alt = _value(row, "label_alt_id")
        if alt and alt != chosen[residue_id(row)]:
            excluded["alternate_location"] += 1
            continue
        occupancy = _number(row, "occupancy")
        if occupancy is not None and occupancy < 0:
            raise ValueError("Negative atom occupancy")
        if occupancy == 0:
            excluded["zero_occupancy"] += 1
            continue
        element = _value(row, "type_symbol").upper()
        if not element:
            raise ValueError(f"Missing type_symbol element for atom-site ID {atom_site_id(row)}")
        is_polymer = entity_types.get(_value(row, "label_entity_id")) == "polymer"
        if selection == "polymer-heavy" and not is_polymer:
            excluded["nonpolymer"] += 1
            continue
        if selection != "all" and element in {"H", "D"}:
            excluded["hydrogen_or_deuterium"] += 1
            continue
        chain = _value(row, f"{chain_ids}_asym_id")
        if not chain:
            raise ValueError(f"Missing {chain_ids} chain ID")
        atom = Atom(
            serial=atom_site_id(row), name=atom_name(row),
            resname=_value(row, "auth_comp_id") or _required(row, "label_comp_id"),
            chain_id=chain, resid=residue_number(row), icode=residue_id(row)[2],
            x=coordinate(row, "Cartn_x"), y=coordinate(row, "Cartn_y"),
            z=coordinate(row, "Cartn_z"),
            element=element, altloc=alt, occupancy=occupancy,
            bfactor=_number(row, "B_iso_or_equiv"),
            hetero=(not is_polymer) if hetero_definition == "nonpolymer"
                   else _value(row, "group_PDB") == "HETATM",
        )
        identity = (atom.residue_key, atom.name)
        if identity in seen:
            raise ValueError(f"Duplicate selected atom identity: {identity}; use label chain IDs")
        seen.add(identity)
        atoms.append(atom)
        polymer_hetero += int(is_polymer and _value(row, "group_PDB") == "HETATM")
        info = chains.setdefault(chain, {"atom_count": 0, "auth_asym_ids": set(),
                                         "label_asym_ids": set(), "entity_ids": set(),
                                         "polymer": is_polymer})
        info["atom_count"] += 1
        info["polymer"] = info["polymer"] or is_polymer
        for field in ("auth_asym", "label_asym", "entity"):
            column = "label_entity_id" if field == "entity" else field + "_id"
            info[field + "_ids"].add(_value(row, column))
    if not atoms:
        raise ValueError(f"The {selection} selection retained no atoms from {source}")
    frame = StructureFrame(tuple(atoms), source=str(source), model_index=model_index)
    for info in chains.values():
        for field, value in info.items():
            if isinstance(value, set):
                info[field] = sorted(value)
    entry_id = block.find_value("_entry.id")
    import gemmi
    report = {
        "parser": "gemmi", "parser_version": gemmi.__version__,
        "entry_id": entry_id, "model_ids": model_ids,
        "entry_id_is_placeholder": _is_placeholder(entry_id),
        "model_index": model_index, "selected_model_id": model_ids[model_index],
        "selection": selection, "chain_id_namespace": chain_ids,
        "residue_numbering": "auth_seq_id, label_seq_id fallback",
        "altloc_policy": "blank plus first lexical nonblank per residue",
        "zero_occupancy_policy": "excluded",
        "hetero_definition": hetero_definition,
        "altloc_residue_count": len(alternatives), "excluded_atoms": dict(excluded),
        "input_model_atom_count": len(rows), "atom_count": len(atoms),
        "residue_count": len(frame.residues()), "chains": chains,
        "entity_types_available": entity_types_available,
        "polymer_chain_count": (sum(info["polymer"] for info in chains.values())
                                if entity_types_available else None),
        "polymer_hetatm_retained": polymer_hetero if entity_types_available else None,
        "unsupported_radius_elements": sorted({a.element for a in atoms} - VDW_RADII.keys()),
        "bounds_A": frame.bounding_box(), "coordinate_units": "Angstrom",
        "assembly_definitions": category_rows(block, "_pdbx_struct_assembly."),
        "unobserved_residue_records": len(category_rows(block, "_pdbx_unobs_or_zero_occ_residues.")),
        "unobserved_atom_records": len(category_rows(block, "_pdbx_unobs_or_zero_occ_atoms.")),
        "completeness_validated": False,
    }
    return frame, report


def residue_id(row: dict) -> tuple[str, str, str]:
    """Group atoms by residue without requiring a parsable residue number.

    Parameters
    ----------
    row : dict of str to str
        One ``_atom_site`` row keyed by item name without the category prefix
        (for example ``"label_asym_id"``); absent or null (non-string) items
        count as empty.

    Returns
    -------
    tuple of str
        ``(label_asym_id, auth_seq_id or label_seq_id, pdbx_PDB_ins_code)``.

    Examples
    --------
    >>> from crevice.mmcif import residue_id
    >>> residue_id({"label_asym_id": "A", "label_seq_id": "7", "auth_seq_id": "107"})
    ('A', '107', '')
    """
    return (_value(row, "label_asym_id"),
            _value(row, "auth_seq_id") or _value(row, "label_seq_id"),
            _value(row, "pdbx_PDB_ins_code"))


def atom_name(row: dict) -> str:
    """Atom name from ``auth_atom_id``, falling back to ``label_atom_id``.

    Parameters
    ----------
    row : dict of str to str
        One ``_atom_site`` row (see :func:`residue_id`).

    Returns
    -------
    str
        The atom name.

    Raises
    ------
    ValueError
        If neither is present.
    """
    return _value(row, "auth_atom_id") or _required(row, "label_atom_id")


def atom_site_id(row: dict) -> int:
    """The integer ``_atom_site.id`` of a row.

    Raises
    ------
    ValueError
        If the ID is missing or not an integer.
    """
    text = _value(row, "id")
    if not text:
        raise ValueError("Missing or duplicate mmCIF atom-site ID")
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"Nonintegral mmCIF atom-site ID: {text!r}") from exc


def residue_number(row: dict) -> int:
    """Integer residue number: ``auth_seq_id``, falling back to ``label_seq_id``.

    Parameters
    ----------
    row : dict of str to str
        One ``_atom_site`` row (see :func:`residue_id`).

    Returns
    -------
    int
        The residue number.

    Raises
    ------
    ValueError
        If neither is present or the value is not an integer.
    """
    text = residue_id(row)[1]
    if not text:
        raise ValueError(f"Missing residue number for atom-site ID {atom_site_id(row)}")
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"Nonintegral residue number: {text!r}") from exc


def coordinate(row: dict, key: str) -> float:
    """One Cartesian coordinate (Å) of a row.

    Parameters
    ----------
    row : dict
    key : {"Cartn_x", "Cartn_y", "Cartn_z"}

    Returns
    -------
    float

    Raises
    ------
    ValueError
        If the value is missing, non-numeric or not finite.
    """
    text = _value(row, key)
    if not text:
        raise ValueError(f"Missing {key} coordinate in mmCIF atom-site table")
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(f"Nonnumeric {key} coordinate: {text!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Coordinates must be finite in every model; {key} is {text!r}")
    return value


def _residue_alternatives(rows: list[dict]) -> dict[tuple, set[str]]:
    """Collect nonblank alternate locations, rejecting mixed conventions.

    An atom carrying both a blank and a labelled alternate location has no
    deterministic conformer, so the file is rejected instead of guessed at.
    """
    alternatives: dict[tuple, set[str]] = {}
    per_atom: dict[tuple, set[str]] = {}
    for row in rows:
        alt = _value(row, "label_alt_id")
        if alt:
            alternatives.setdefault(residue_id(row), set()).add(alt)
        per_atom.setdefault((residue_id(row), atom_name(row)), set()).add(alt)
    for (residue, name), alts in per_atom.items():
        if "" in alts and len(alts) > 1:
            raise ValueError(
                f"Atom {name} of residue {residue} has both blank and labelled "
                "alternate locations"
            )
    return alternatives


def _model_index(model_index) -> int:
    """Validate a model index as a non-negative integer (``bool`` rejected)."""
    if isinstance(model_index, bool):
        raise ValueError("model_index must be a non-negative integer")
    try:
        model_index = operator.index(model_index)
    except TypeError as exc:
        raise ValueError("model_index must be a non-negative integer") from exc
    if model_index < 0:
        raise ValueError("model_index must be a non-negative integer")
    return model_index


def _is_placeholder(entry_id) -> bool:
    """True for a missing, non-string or placeholder (``.``, ``?``, ``XXXX``) entry ID."""
    if not isinstance(entry_id, str):
        return True
    return entry_id.strip().strip("'\"").upper() in PLACEHOLDER_ENTRY_IDS


def _value(row: dict, key: str) -> str:
    """String value of a row item, or ``""`` for absent or null items."""
    value = row.get(key)
    return value if isinstance(value, str) else ""


def _required(row: dict, key: str) -> str:
    """String value of a required item; raise ``ValueError`` if absent."""
    value = _value(row, key)
    if not value:
        raise ValueError(f"Missing required mmCIF atom-site value {key}")
    return value


def _number(row: dict, key: str) -> float | None:
    """Finite float value of an item, ``None`` if absent; raise ``ValueError`` otherwise."""
    value = _value(row, key)
    if not value:
        return None
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"Nonnumeric {key}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"Nonfinite {key}")
    return number

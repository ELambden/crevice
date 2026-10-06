"""Load PDB, mmCIF and GRO files into a :class:`~crevice.models.StructureFrame`.

Dispatch by file suffix
    * ``.gro``: read through MDAnalysis (:mod:`crevice.trajectory`) with an
      explicit atom selection (``"protein"`` by default); coordinates are
      converted from nm to Å.
    * ``.cif`` / ``.mmcif`` (uncompressed): the entity-aware Gemmi loader
      (:func:`crevice.mmcif.load_mmcif`) when Gemmi is installed, otherwise
      the built-in reader.
    * anything else: the built-in fixed-column PDB reader.

What the readers keep
    All readers keep blank alternate locations plus the first lexical
    non-blank location of each residue (:data:`ALTLOC_POLICY`), and both mmCIF
    readers honour ``chain_ids`` (label by default). Hydrogens, waters, ions
    and ligands are kept; filtering happens later through
    :meth:`~crevice.models.StructureFrame.selected_atoms`.

Deliberate differences between the mmCIF readers (recorded in each report)
    * ``Atom.hetero``: the built-in readers use the ``HETATM`` record type
      (``hetero_definition="group_pdb"``), so a modified polymer residue
      deposited as HETATM counts as a heteroatom. The Gemmi path uses entity
      membership when ``_entity`` types are present
      (``hetero_definition="nonpolymer"``), so it stays in the polymer.
    * Zero-occupancy atoms: kept by the built-in readers
      (``zero_occupancy_policy="retained"``), dropped by Gemmi
      (``"excluded"``).
    * Integrity: Gemmi validates the whole atom-site table and rejects an atom
      with both blank and labelled alternate locations; the built-in readers
      keep both.

The built-in readers infer missing elements from atom names
(:func:`crevice.radii.infer_element`). No reader checks completeness (missing
residues or atoms are not reported).

Main entry points: :func:`load_structure` and :func:`load_structure_report`.
"""

from __future__ import annotations

import re
import math
import shlex
import warnings
from dataclasses import replace
from pathlib import Path

from .models import Atom, StructureFrame
from .radii import infer_element, is_recognised_ion, partition_atoms_with_radii


PARSERS = ("auto", "gemmi", "builtin")

#: Alternate-location rule shared by every structure reader.
ALTLOC_POLICY = "blank plus first lexical nonblank per residue"


def is_mmcif_path(path: str | Path) -> bool:
    """Return ``True`` when the file name ends in ``.cif`` or ``.mmcif``.

    All suffixes are joined before the test, so ``model.cif`` qualifies but
    ``model.cif.gz`` does not.

    Parameters
    ----------
    path : str or Path
        File name or path (the file need not exist).

    Returns
    -------
    bool
    """
    suffixes = "".join(Path(path).suffixes).lower()
    return suffixes.endswith(".cif") or suffixes.endswith(".mmcif")


def load_structure(
    path: str | Path,
    *,
    model_index: int = 0,
    parser: str = "auto",
    chain_ids: str = "label",
    md_selection: str = "protein",
    radii=None,
) -> StructureFrame:
    r"""Load a structure file and return only the frame.

    Equivalent to ``load_structure_report(...)[0]``; see
    :func:`load_structure_report` for the full description of each option.

    Parameters
    ----------
    path : str or Path
        PDB, mmCIF or GRO file.
    model_index : int, default 0
        Zero-based model (GRO: frame) to read.
    parser : {"auto", "gemmi", "builtin"}, default "auto"
        mmCIF reader choice; GRO requires ``"auto"``.
    chain_ids : {"label", "auth"}, default "label"
        mmCIF chain-ID namespace.
    md_selection : str, default "protein"
        MDAnalysis selection for GRO input.
    radii : RadiusSet, str, os.PathLike or None, default None
        Radius set used to decide which atoms have a radius (atoms without
        one are excluded with a warning); ``None`` uses the set in effect.

    Returns
    -------
    StructureFrame

    Examples
    --------
    >>> import pathlib, tempfile
    >>> from crevice.parser import load_structure
    >>> path = pathlib.Path(tempfile.mkdtemp()) / "tiny.pdb"
    >>> _ = path.write_text(
    ...     "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n")
    >>> frame = load_structure(path)
    >>> frame.atoms[0].residue_key.label, frame.atoms[0].coord
    ('A:ALA1', (1.0, 2.0, 3.0))
    """

    return load_structure_report(
        path, model_index=model_index, parser=parser, chain_ids=chain_ids, md_selection=md_selection,
        radii=radii,
    )[0]


def load_structure_report(
    path: str | Path,
    *,
    model_index: int = 0,
    parser: str = "auto",
    chain_ids: str = "label",
    md_selection: str = "protein",
    radii=None,
) -> tuple[StructureFrame, dict]:
    r"""Load a structure and report which parser ran and what it did.

    For mmCIF, ``"auto"`` prefers the validated entity-aware Gemmi loader and
    falls back to the builtin parser when the ``structures`` extra is absent;
    the report always names the parser that actually ran. Hydrogen and
    heteroatom filtering is left to :meth:`~crevice.models.StructureFrame.selected_atoms`.
    GRO uses the MDAnalysis reader and defaults to ``md_selection="protein"``;
    its report records excluded atoms and the conversion from nm to Å.

    Parameters
    ----------
    path : str or pathlib.Path
        Structure file.
    model_index : int, default 0
        Zero-based model to read (for GRO, the frame index).
    parser : {"auto", "gemmi", "builtin"}, default "auto"
        ``"gemmi"`` requires an mmCIF file and the Gemmi package; ``"builtin"``
        forces the dependency-free reader; GRO input requires ``"auto"``.
    chain_ids : {"label", "auth"}, default "label"
        Chain-ID namespace for mmCIF input, used by both mmCIF readers
        (``label_asym_id`` or ``auth_asym_id``; the built-in reader falls back
        to the other column where the requested one is null). Ignored for PDB
        and GRO input.
    md_selection : str, default "protein"
        MDAnalysis selection for GRO input.
    radii : RadiusSet, str, os.PathLike or None, default None
        Radius set in force for this load (:func:`crevice.radii.resolve_radii`;
        ``None`` uses the set in effect, normally the CLI ``--radii`` choice or
        :data:`~crevice.radii.DEFAULT_RADII`). Atoms that get no radius from
        it (:meth:`~crevice.radii.RadiusSet.has_radius`) are excluded with an
        :class:`~crevice.radii.UnrecognisedElementWarning`. Pass the same set
        you will analyse with.

    Returns
    -------
    frame : StructureFrame
    report : dict
        ``parser``, ``input_format``, ``model_index``, atom and residue counts,
        the policies applied (``chain_id_namespace``, ``altloc_policy``,
        ``zero_occupancy_policy``, ``hetero_definition``, element source),
        ``unsupported_radius_elements`` (elements of the excluded atoms),
        ``excluded_no_radius_atoms`` (count per element),
        ``excluded_no_radius_names`` (count per ``"RESNAME:ATOMNAME"``),
        ``radius_set`` (name of the set checked), ``radius_note``,
        AlphaFold flags from :func:`predicted_model_report` and
        ``recognised_ion_atoms`` (count per ``"RESNAME:ATOMNAME"`` of the
        recognised monatomic ion residues, which get ion radii in the default
        radius set; see :func:`crevice.radii.is_recognised_ion`).
        ``completeness_validated`` is always ``False``.

    Raises
    ------
    ValueError
        For an unknown parser, negative ``model_index``, a non-mmCIF file with
        ``parser="gemmi"``, a non-``"auto"`` parser for GRO, or a file with no
        readable atoms.
    RuntimeError
        If ``parser="gemmi"`` is requested but Gemmi is not installed.

    Notes
    -----
    Not every atom in the file is returned. Every reader drops the
    alternate locations that :data:`ALTLOC_POLICY` does not select and atoms
    of other models. The Gemmi path (``selection="all"`` in
    :func:`crevice.mmcif.load_mmcif`) also drops zero-occupancy atoms, and its
    ``Atom.hetero`` follows entity membership where entity types are
    available. Hydrogens, waters, ions and ligands are kept by every reader.
    The two mmCIF paths can therefore give different atom sets and hetero
    flags for the same file; check ``report["parser"]``,
    ``report["hetero_definition"]`` and ``report["zero_occupancy_policy"]``.

    Examples
    --------
    >>> import pathlib, tempfile
    >>> from crevice.parser import load_structure_report
    >>> path = pathlib.Path(tempfile.mkdtemp()) / "tiny.pdb"
    >>> _ = path.write_text(
    ...     "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n"
    ...     "HETATM    2  O   HOH W   5       6.000   0.000   0.000  1.00  0.00           O\n")
    >>> frame, report = load_structure_report(path)
    >>> report["parser"], report["atom_count"], [a.hetero for a in frame.atoms]
    ('builtin', 2, [False, True])
    """

    from .radii import use_radii

    with use_radii(radii):
        return _load_structure_report(Path(path), model_index=model_index, parser=parser,
                                      chain_ids=chain_ids, md_selection=md_selection)


def _load_structure_report(source: Path, *, model_index: int, parser: str, chain_ids: str,
                           md_selection: str) -> tuple[StructureFrame, dict]:
    """Body of :func:`load_structure_report`, run inside the chosen radius set."""
    if parser not in PARSERS:
        raise ValueError(f"parser must be one of {', '.join(PARSERS)}")
    if model_index < 0:
        raise ValueError("model_index must be non-negative")
    if source.suffix.lower() == ".gro":
        if parser != "auto":
            raise ValueError("GRO input uses MDAnalysis; use parser=auto")
        from .trajectory import load_trajectory_report
        trajectory, report = load_trajectory_report(source, selection=md_selection,
                                                    start=model_index, stop=model_index+1)
        frame = trajectory.frames[0]
        return _with_radius_note(frame, {**report, "parser": "mdanalysis", "input_format": "gro",
                                                "coordinate_units": "Angstrom (converted from GRO nm)"})
    mmcif_input = is_mmcif_path(source)
    if parser == "gemmi" and not mmcif_input:
        raise ValueError(
            f"The gemmi parser reads mmCIF; {source.name} is not a .cif or .mmcif file"
        )
    if mmcif_input and parser != "builtin":
        from . import mmcif

        available = mmcif.gemmi_available()
        if parser == "gemmi" and not available:
            raise RuntimeError(
                "Validated mmCIF intake requires pip install 'crevice[structures]'"
            )
        if available:
            frame, report = mmcif.load_mmcif(
                source, model_index=model_index, selection="all",
                chain_ids=chain_ids, hetero_policy="auto",
            )
            with source.open(errors="replace") as handle:
                header = handle.read(65536)
            return _with_radius_note(frame, {**report, "input_format": "mmcif",
                                             **predicted_model_report(header)})

    text = source.read_text()
    if chain_ids not in {"label", "auth"}:
        raise ValueError("chain_ids must be label or auth")
    if mmcif_input:
        atoms = _parse_mmcif(text, model_index=model_index, chain_ids=chain_ids)
    else:
        atoms = _parse_pdb(text, model_index=model_index)
    if not atoms:
        raise ValueError(f"No atoms could be loaded from {source}")
    frame = StructureFrame(tuple(atoms), source=str(source), model_index=model_index)
    return _with_radius_note(frame, {
        **predicted_model_report(text),
        "parser": "builtin",
        "input_format": "mmcif" if mmcif_input else "pdb",
        "model_index": model_index,
        "atom_count": len(frame.atoms),
        "residue_count": len(frame.residues()),
        "chain_id_namespace": chain_ids if mmcif_input else "PDB chain column",
        "altloc_policy": ALTLOC_POLICY,
        "zero_occupancy_policy": "retained",
        "hetero_definition": "group_pdb",
        "element_source": "type_symbol or name inference",
        "completeness_validated": False,
    })


def _with_radius_note(frame: StructureFrame, report: dict) -> tuple[StructureFrame, dict]:
    """Exclude atoms without a radius, warn, and add the radius entries to a report.

    Parameters
    ----------
    frame : StructureFrame
        The loaded structure.
    report : dict
        The reader report; modified in place.

    Returns
    -------
    frame : StructureFrame
        ``frame`` without the atoms that get no radius from the set in effect
        (:func:`~crevice.radii.partition_atoms_with_radii`).
    report : dict
        ``report`` with ``unsupported_radius_elements`` (elements of the
        excluded atoms), ``excluded_no_radius_atoms`` (count per element),
        ``excluded_no_radius_names`` (count per ``"RESNAME:ATOMNAME"``),
        ``radius_set``, ``radius_note``, an updated ``atom_count`` and
        ``recognised_ion_atoms`` (``{"RESNAME:ATOMNAME": count}`` of recognised
        monatomic ion residues, :func:`~crevice.radii.is_recognised_ion`).

    Raises
    ------
    ValueError
        If no atom has a radius.
    """
    from collections import Counter
    from .radii import UnrecognisedElementWarning, effective_radii, unrecognised_element_message
    kept, by_element, by_name = partition_atoms_with_radii(frame.atoms)
    report["radius_set"] = effective_radii().name
    report["unsupported_radius_elements"] = sorted(by_element)
    report["excluded_no_radius_atoms"] = by_element
    report["excluded_no_radius_names"] = by_name
    if by_element:
        note = unrecognised_element_message(by_element, by_name)
        warnings.warn(f"{frame.source}: {note}", UnrecognisedElementWarning, stacklevel=4)
        if not kept:
            raise ValueError(f"No atom of {frame.source} has a radius in the radius set in force. {note}")
        frame = replace(frame, atoms=kept)
        if "atom_count" in report:
            report["atom_count"] = len(kept)
    else:
        note = f"every atom has a radius in radius set {effective_radii().name!r}"
    report["radius_note"] = note
    ions = Counter(f"{atom.resname.strip().upper()}:{atom.name.strip().upper()}" for atom in frame.atoms
                   if is_recognised_ion(atom.name, atom.resname, atom.element))
    report["recognised_ion_atoms"] = dict(sorted(ions.items()))
    return frame, report


def predicted_model_report(text: str) -> dict:
    """Flag a computationally predicted model (currently AlphaFold) from its header.

    AlphaFold files store per-residue confidence (pLDDT, 0-100) in the
    B-factor column, and their coordinates are a prediction rather than an
    experimental structure. The first 64 kB are searched for ``ALPHAFOLD``
    (PDB ``TITLE``, mmCIF software/title records) or an ``AF-`` entry ID.

    Parameters
    ----------
    text : str
        File content (only the start is inspected).

    Returns
    -------
    dict
        Empty for other files; otherwise ``model_type="predicted"``,
        ``predictor="AlphaFold"`` and ``bfactor_meaning``.

    Examples
    --------
    >>> from crevice.parser import predicted_model_report
    >>> predicted_model_report("TITLE     ALPHAFOLD MONOMER V2.0 PREDICTION FOR X")["predictor"]
    'AlphaFold'
    >>> predicted_model_report("HEADER    MEMBRANE PROTEIN")
    {}
    """
    head = text[:65536]
    if "ALPHAFOLD" in head.upper() or re.search(r"(?m)^_entry\.id\s+['\"]?AF-", head):
        return {"model_type": "predicted", "predictor": "AlphaFold",
                "bfactor_meaning": "pLDDT per-residue confidence (0-100), not a crystallographic B-factor; "
                                   "treat low-confidence regions with caution"}
    return {}


def _parse_pdb(text: str, *, model_index: int = 0) -> list[Atom]:
    """Fixed-column PDB reader for one model.

    Without ``MODEL`` records every ``ATOM``/``HETATM`` belongs to model 0.
    Alternate locations follow :data:`ALTLOC_POLICY` (see :func:`_select_altlocs`).
    """
    atoms: list[Atom] = []
    has_models = any(line[:6].strip() == "MODEL" for line in text.splitlines())
    current_model = -1 if has_models else 0
    active = not has_models
    for line in text.splitlines():
        record = line[0:6].strip()
        if record == "MODEL":
            current_model += 1
            active = True
            continue
        if record == "ENDMDL":
            active = False
        if not active or current_model != model_index:
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        altloc = _slice(line, 16, 17).strip()
        serial = _safe_int(_slice(line, 6, 11), -1)
        name = _slice(line, 12, 16).strip()
        resname = _slice(line, 17, 20).strip() or "UNK"
        chain_id = _slice(line, 21, 22).strip()
        resid = _safe_int(_slice(line, 22, 26), 0)
        icode = _slice(line, 26, 27).strip()
        x = _coordinate(_slice(line, 30, 38))
        y = _coordinate(_slice(line, 38, 46))
        z = _coordinate(_slice(line, 46, 54))
        occupancy = _optional_float(_slice(line, 54, 60))
        bfactor = _optional_float(_slice(line, 60, 66))
        element = infer_element(name, _slice(line, 76, 78), resname)
        atoms.append(
            Atom(
                serial=serial,
                name=name,
                resname=resname,
                chain_id=chain_id,
                resid=resid,
                icode=icode,
                x=x,
                y=y,
                z=z,
                element=element,
                altloc=altloc,
                occupancy=occupancy,
                bfactor=bfactor,
                hetero=record == "HETATM",
            )
        )
    return _select_altlocs(atoms)


def _select_altlocs(atoms: list[Atom]) -> list[Atom]:
    """Keep blank alternate locations plus the first lexical non-blank one per residue.

    Residues are grouped by ``(chain_id, resid, icode)``. This is the rule of
    the Gemmi loader (:func:`crevice.mmcif.load_mmcif`), a deterministic choice
    rather than conformer inference; a residue deposited only as ``B``/``C``
    keeps ``B``. Atoms parsed without a serial (marked ``-1``) are numbered by
    their position among the kept atoms.
    """
    chosen: dict[tuple, str] = {}
    for atom in atoms:
        if atom.altloc:
            key = (atom.chain_id, atom.resid, atom.icode)
            chosen[key] = min(chosen.get(key, atom.altloc), atom.altloc)
    kept: list[Atom] = []
    for atom in atoms:
        if atom.altloc and atom.altloc != chosen[(atom.chain_id, atom.resid, atom.icode)]:
            continue
        if atom.serial == -1:
            atom = replace(atom, serial=len(kept) + 1)
        kept.append(atom)
    return kept


def _parse_mmcif(text: str, *, model_index: int = 0, chain_ids: str = "label") -> list[Atom]:
    """Minimal ``_atom_site`` loop reader for one model.

    Chain IDs come from ``label_asym_id`` or ``auth_asym_id`` as ``chain_ids``
    requests; a null value falls back to the other column for that row. For
    atom name, residue name and residue number the author columns are
    preferred, with the same per-row fallback to the label columns (as in the
    Gemmi loader). Models are numbered in order of first appearance of
    ``pdbx_PDB_model_num``. Alternate locations follow :data:`ALTLOC_POLICY`.
    Quoted multi-line values are not supported.
    """
    lines = text.splitlines()
    atoms: list[Atom] = []
    model_ids: list[str] = []
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped != "loop_":
            i += 1
            continue
        i += 1
        tags: list[str] = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            tags.append(lines[i].strip())
            i += 1
        if not tags or not all(tag.startswith("_atom_site.") for tag in tags):
            while i < len(lines) and lines[i].strip() not in {"loop_", "#"}:
                i += 1
            continue

        values: list[str] = []
        while i < len(lines):
            line = lines[i].strip()
            if not line or line == "#":
                i += 1
                break
            if line == "loop_" or line.startswith("_") or line.startswith("data_"):
                break
            values.extend(_tokenize_cif_line(line))
            i += 1

        width = len(tags)
        rows = [values[start : start + width] for start in range(0, len(values), width)]
        tag_index = {tag: idx for idx, tag in enumerate(tags)}

        def cols(*names: str) -> tuple[int, ...]:
            """Every present column for a field, in preference order."""
            return tuple(
                tag_index[f"_atom_site.{name}"]
                for name in names
                if f"_atom_site.{name}" in tag_index
            )

        columns = {
            "record": cols("group_PDB"),
            "serial": cols("id"),
            "element": cols("type_symbol"),
            "name": cols("auth_atom_id", "label_atom_id"),
            "altloc": cols("label_alt_id"),
            "resname": cols("auth_comp_id", "label_comp_id"),
            "chain": (cols("label_asym_id", "auth_asym_id") if chain_ids == "label"
                      else cols("auth_asym_id", "label_asym_id")),
            "resid": cols("auth_seq_id", "label_seq_id"),
            "icode": cols("pdbx_PDB_ins_code"),
            "x": cols("Cartn_x"),
            "y": cols("Cartn_y"),
            "z": cols("Cartn_z"),
            "occupancy": cols("occupancy"),
            "bfactor": cols("B_iso_or_equiv"),
            "model": cols("pdbx_PDB_model_num"),
        }
        required = {"name", "resname", "chain", "resid", "x", "y", "z"}
        if any(not columns[name] for name in required):
            continue

        def value(row: list[str], key: str, default: str = "") -> str:
            """First non-null value for a field, falling back per row.

            An auth column present but null (`.` or `?`) must not shadow the
            label column; the fallback is per row, not per table.
            """
            for index in columns[key]:
                cleaned = _clean(row[index])
                if cleaned:
                    return cleaned
            return default

        for row in rows:
            if len(row) != width:
                raise ValueError("Incomplete mmCIF atom-site row")
            model_id = value(row, "model", "1")
            if model_id not in model_ids:
                model_ids.append(model_id)
            if model_ids.index(model_id) != model_index:
                continue
            record = value(row, "record", "ATOM")
            if record not in {"ATOM", "HETATM"}:
                continue
            altloc = value(row, "altloc")
            name = value(row, "name")
            resname = value(row, "resname", "UNK")
            atoms.append(
                Atom(
                    serial=_safe_int(value(row, "serial"), -1),
                    name=name,
                    resname=resname,
                    chain_id=value(row, "chain"),
                    resid=_safe_int(value(row, "resid"), 0),
                    icode=value(row, "icode"),
                    x=_coordinate(value(row, "x")),
                    y=_coordinate(value(row, "y")),
                    z=_coordinate(value(row, "z")),
                    element=infer_element(name, value(row, "element"), resname),
                    altloc=altloc,
                    occupancy=_optional_float(value(row, "occupancy")),
                    bfactor=_optional_float(value(row, "bfactor")),
                    hetero=record == "HETATM",
                )
            )
    return _select_altlocs(atoms)


def _coordinate(text: str) -> float:
    """Parse a coordinate; raise ``ValueError`` if it is not a finite number."""
    value = float(_clean(text))
    if not math.isfinite(value):
        raise ValueError("Coordinates must be finite")
    return value


def _slice(line: str, start: int, end: int) -> str:
    """Fixed-column slice that tolerates short lines."""
    return line[start:end] if len(line) >= start else ""


def _safe_int(text: str, default: int) -> int:
    """Parse an integer, else the first signed digit run, else ``default``."""
    cleaned = _clean(text)
    try:
        return int(cleaned)
    except ValueError:
        match = re.search(r"-?\d+", cleaned)
        return int(match.group(0)) if match else default


def _safe_float(text: str, default: float) -> float:
    """Parse a float, else return ``default``."""
    cleaned = _clean(text)
    try:
        return float(cleaned)
    except ValueError:
        return default


def _optional_float(text: str) -> float | None:
    """Parse a float, or ``None`` for blank, null or non-numeric text."""
    cleaned = _clean(text)
    if not cleaned:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def _clean(text: str) -> str:
    """Strip whitespace and quotes; mmCIF null markers ``.`` and ``?`` become ``""``."""
    value = text.strip().strip("'\"")
    return "" if value in {".", "?"} else value


def _tokenize_cif_line(line: str) -> list[str]:
    """Split one mmCIF data line into cleaned tokens, honouring quotes."""
    lexer = shlex.shlex(line, posix=False)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return [_clean(token) for token in lexer]

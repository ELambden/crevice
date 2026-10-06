"""Van der Waals radii, element inference and coarse residue classes.

Every clearance in CREVICE is measured to atom *surfaces*, so the atomic
radius directly sets every reported radius, volume and contact distance.
This module supplies:

* :data:`VDW_RADII`: element radii in Å from Bondi (1964), with Mantina et
  al. (2009) for main-group elements Bondi did not tabulate (for example B,
  Ca, Rb, Cs, Sr, Ba) and three legacy values (Fe, Mn, Co);
  :func:`radius_source` names the source of each. Any element not listed
  (most transition metals and lanthanides) keeps its symbol and gets
  :data:`DEFAULT_RADIUS` (1.70 Å, the carbon value); the structure readers
  report such elements (:func:`default_radius_elements`).
* :data:`ION_RADII`: CHARMM36 Lennard-Jones Rmin/2 radii of monatomic ions
  (``toppar_water_ions.str``), used by default for recognised ion residues
  only (:func:`is_recognised_ion`).
* :data:`ELEMENT_SYMBOLS` and :func:`infer_element`: element symbols from an
  explicit element column (any real symbol is kept) or, failing that, from
  PDB/mmCIF/MD atom names.
* :data:`ION_NAME_ELEMENTS`: residue/atom names of monatomic MD ions
  (CHARMM36, GROMACS, AMBER) mapped to their elements.
* :func:`atom_vdw_radius`: the radius that the analyses actually use.
* :func:`residue_properties`: residue-name classes (hydrophobic, polar,
  charged, aromatic, small) used to label contacts and network edges.
* :class:`RadiusSet`, :data:`DEFAULT_RADII`, :func:`radii_preset`,
  :func:`custom_radii_set`, :func:`load_radius_file`: radius tables (presets
  ``default``, ``bondi``, ``hole``, ``charmm_like``; JSON, CSV or HOLE
  ``.rad`` files).
* :func:`use_radii`, :func:`effective_radii`, :func:`resolve_radii`,
  :func:`radii_provenance`, :func:`radius_set_from_provenance`,
  :func:`same_radius_set`: choose the set the analyses use, record it, and
  rebuild or compare recorded sets.

Notes
-----
Every analysis obtains atomic radii from :func:`atom_vdw_radius`. With no
radius set chosen it uses :data:`DEFAULT_RADII`: the standard
:data:`VDW_RADII` element table for every atom, except that a recognised
monatomic ion residue (for example CHARMM ``SOD``/``SOD``, GROMACS ``NA``/``NA``
or a crystallographic ``ZN``/``ZN``) gets the CHARMM36 radius of
:data:`ION_RADII`. Before 2 Oct 2026 such ions got their neutral-atom
:data:`VDW_RADII` value (preset ``bondi``). A set is chosen with an analysis's
``radii=`` argument, a :func:`use_radii` block or the CLI ``--radii`` option;
all three set one :mod:`contextvars` variable that :func:`atom_vdw_radius`
reads, so nested analyses and trajectory worker processes use the same set.
Outputs always record the set they used (``radii`` in result metadata and
manifests: name, source, file SHA-256, overrides, ion radii and a table
hash), including default runs.

Residue classes are based on residue names only; they do not account for
protonation state, local environment or modified residues.

Examples
--------
>>> from crevice.radii import infer_element, vdw_radius, residue_properties
>>> infer_element("CA", resname="ALA"), infer_element("CA", resname="CA")
('C', 'CA')
>>> vdw_radius("O")
1.52
>>> residue_properties("HIS")
('polar', 'positive', 'aromatic')
"""

from __future__ import annotations

import contextlib
import contextvars
import functools
import os
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Mapping

from .models import Atom

#: Every element symbol (Z = 1-118), upper case. Used to decide whether an
#: explicit element column holds a real symbol; the list is embedded so the
#: result does not depend on optional packages.
ELEMENT_SYMBOLS: frozenset[str] = frozenset("""
H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI CU ZN
GA GE AS SE BR KR RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I XE CS BA LA
CE PR ND PM SM EU GD TB DY HO ER TM YB LU HF TA W RE OS IR PT AU HG TL PB BI PO
AT RN FR RA AC TH PA U NP PU AM CM BK CF ES FM MD NO LR RF DB SG BH HS MT DS RG
CN NH FL MC LV TS OG
""".split())

# Sources of the tabulated radii (Å):
#   Bondi 1964: A. Bondi, J. Phys. Chem. 68, 441-451 (1964), Table I.
#   Mantina 2009: M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer and
#     D. G. Truhlar, J. Phys. Chem. A 113, 5806-5812 (2009), Table 12
#     (main-group elements that Bondi did not tabulate).
#   legacy: values carried over from the first CREVICE table whose source
#     was not recorded. They are kept unchanged so results stay reproducible.
_RADIUS_SOURCES: dict[str, tuple[str, float]] = {
    # Bondi 1964
    "H": ("Bondi 1964", 1.20), "C": ("Bondi 1964", 1.70), "N": ("Bondi 1964", 1.55),
    "O": ("Bondi 1964", 1.52), "F": ("Bondi 1964", 1.47), "P": ("Bondi 1964", 1.80),
    "S": ("Bondi 1964", 1.80), "CL": ("Bondi 1964", 1.75), "BR": ("Bondi 1964", 1.85),
    "I": ("Bondi 1964", 1.98), "SE": ("Bondi 1964", 1.90), "MG": ("Bondi 1964", 1.73),
    "ZN": ("Bondi 1964", 1.39), "NA": ("Bondi 1964", 2.27), "K": ("Bondi 1964", 2.75),
    "LI": ("Bondi 1964", 1.82), "CU": ("Bondi 1964", 1.40), "NI": ("Bondi 1964", 1.63),
    "HE": ("Bondi 1964", 1.40), "NE": ("Bondi 1964", 1.54), "AR": ("Bondi 1964", 1.88),
    "KR": ("Bondi 1964", 2.02), "XE": ("Bondi 1964", 2.16), "SI": ("Bondi 1964", 2.10),
    "AS": ("Bondi 1964", 1.85), "TE": ("Bondi 1964", 2.06), "GA": ("Bondi 1964", 1.87),
    "IN": ("Bondi 1964", 1.93), "SN": ("Bondi 1964", 2.17), "TL": ("Bondi 1964", 1.96),
    "PB": ("Bondi 1964", 2.02), "PD": ("Bondi 1964", 1.63), "AG": ("Bondi 1964", 1.72),
    "CD": ("Bondi 1964", 1.58), "PT": ("Bondi 1964", 1.72), "AU": ("Bondi 1964", 1.66),
    "HG": ("Bondi 1964", 1.55), "U": ("Bondi 1964", 1.86),
    # Mantina 2009
    "B": ("Mantina 2009", 1.92), "CA": ("Mantina 2009", 2.31), "BE": ("Mantina 2009", 1.53),
    "AL": ("Mantina 2009", 1.84), "GE": ("Mantina 2009", 2.11), "SB": ("Mantina 2009", 2.06),
    "BI": ("Mantina 2009", 2.07), "PO": ("Mantina 2009", 1.97), "AT": ("Mantina 2009", 2.02),
    "RN": ("Mantina 2009", 2.20), "RB": ("Mantina 2009", 3.03), "CS": ("Mantina 2009", 3.43),
    "SR": ("Mantina 2009", 2.49), "BA": ("Mantina 2009", 2.68), "FR": ("Mantina 2009", 3.48),
    "RA": ("Mantina 2009", 2.83),
    # legacy, source not recorded
    "FE": ("legacy", 1.94), "MN": ("legacy", 1.97), "CO": ("legacy", 1.67),
}

#: Radius in Å used for any element not in :data:`VDW_RADII` (the carbon value).
DEFAULT_RADIUS = 1.70

#: Van der Waals radii in Å by upper-case element symbol. Values are from
#: Bondi (1964), Mantina et al. (2009) for main-group elements Bondi lacked,
#: and three legacy values (Fe, Mn, Co) whose source was not recorded; see
#: :func:`radius_source`. Elements not listed (most transition metals and
#: lanthanides) get :data:`DEFAULT_RADIUS`.
VDW_RADII: dict[str, float] = {key: value for key, (_, value) in _RADIUS_SOURCES.items()}

PROTEIN_CARBON_NAMES = {
    "C",
    "CA",
    "CB",
    "CD",
    "CD1",
    "CD2",
    "CE",
    "CE1",
    "CE2",
    "CE3",
    "CG",
    "CG1",
    "CG2",
    "CH2",
    "CZ",
    "CZ2",
    "CZ3",
}

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "HIS", "TRP"}
POSITIVE = {"ARG", "LYS", "HIS"}
NEGATIVE = {"ASP", "GLU"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
SMALL = {"GLY", "ALA", "SER", "CYS"}

#: Two-letter element prefixes of atom names that are inferred from the name
#: alone (no element column). None of them starts a common carbon or nitrogen
#: atom name.
UNAMBIGUOUS_NAME_ELEMENTS = frozenset({"CL", "BR", "FE", "MG", "MN", "ZN", "LI", "CU", "SE"})
#: Two-letter element symbols that are also common carbon or nitrogen atom
#: names (C-alpha ``CA``, heme ``NA``, ligand ``CO1``/``NI``...). From the name
#: alone they are read as the metal only for a single-atom ion residue whose
#: residue name and atom name are both that symbol.
AMBIGUOUS_NAME_ELEMENTS = frozenset({"CA", "NA", "CO", "NI"})

# Monatomic ions of common MD force fields, keyed by (residue name, atom name),
# both upper case. A rule applies only when BOTH names match, so C-alpha "CA"
# of ALA, carbon "CD" of PRO or heme nitrogen "NA" are never read as ions.
# Sources (checked against the files where a local copy was available):
#   CHARMM36 toppar_water_ions.str (RESI/ATOM records; the copy checked is
#     from a CHARMM-GUI toppar set and carries no version string): LIT, SOD,
#     MG, POT, CAL, RUB, CES, BAR, ZN2 (atom ZN), CD2 (atom CD), CLA.
#   GROMACS 2022.4 share/top/*.ff/ions.itp (amber03/94/96/99/99sb/99sb-ildn/GS,
#     charmm27, gromos43a1-54a7, oplsaa) and share/top/residuetypes.dat ("Ion"
#     entries): NA, K, CL, CA, MG, ZN, LI, RB, CS, CU, CU1 (atom CU), F, BR, I,
#     NA+ and CL- (older GROMACS names).
#   AMBER legacy leap names (Joung and Cheatham, J. Phys. Chem. B 112, 9020,
#     2008; AmberTools atomic_ions.lib before the names were changed to the
#     element symbols): Li+, Na+, K+, Rb+, Cs+, F-, Cl-, Br-, I-. No AMBER file
#     was available offline, so these names were not checked against one.
# Pseudo-atoms such as GROMACS "IB+" (big positive ion) have no element and are
# deliberately absent.
_CHARMM36_IONS = {
    ("LIT", "LIT"): "LI", ("SOD", "SOD"): "NA", ("MG", "MG"): "MG", ("POT", "POT"): "K",
    ("CAL", "CAL"): "CA", ("RUB", "RUB"): "RB", ("CES", "CES"): "CS", ("BAR", "BAR"): "BA",
    ("ZN2", "ZN"): "ZN", ("CD2", "CD"): "CD", ("CLA", "CLA"): "CL",
}
_GROMACS_IONS = {
    ("NA", "NA"): "NA", ("K", "K"): "K", ("CL", "CL"): "CL", ("CA", "CA"): "CA",
    ("MG", "MG"): "MG", ("ZN", "ZN"): "ZN", ("LI", "LI"): "LI", ("RB", "RB"): "RB",
    ("CS", "CS"): "CS", ("CU", "CU"): "CU", ("CU1", "CU"): "CU", ("F", "F"): "F",
    ("BR", "BR"): "BR", ("I", "I"): "I", ("NA+", "NA+"): "NA", ("CL-", "CL-"): "CL",
}
_AMBER_LEGACY_IONS = {
    ("LI+", "LI+"): "LI", ("NA+", "NA+"): "NA", ("K+", "K+"): "K", ("RB+", "RB+"): "RB",
    ("CS+", "CS+"): "CS", ("F-", "F-"): "F", ("CL-", "CL-"): "CL", ("BR-", "BR-"): "BR",
    ("I-", "I-"): "I",
}
#: Element of each recognised monatomic MD ion, keyed by upper-case
#: ``(residue name, atom name)``. Used by :func:`infer_element` only when a
#: file has no usable element column. Sources: CHARMM36
#: ``toppar_water_ions.str``, GROMACS 2022.4 ``ions.itp``/``residuetypes.dat``
#: and the legacy AMBER (Joung-Cheatham) names; see :func:`ion_name_source`.
ION_NAME_ELEMENTS: dict[tuple[str, str], str] = {**_AMBER_LEGACY_IONS, **_GROMACS_IONS, **_CHARMM36_IONS}


# CHARMM36 Lennard-Jones Rmin/2 (Å) of the monatomic ions, from the NONBONDED
# section of toppar_water_ions.str (columns: atom type, ignored, epsilon,
# Rmin/2). The copy read is the one CHARMM-GUI v3.7 wrote into a membrane
# builder job on 25 Aug 2026 (toppar/ directory; the file itself carries no
# version string; SHA-256 41c6474470fa66e927198582a314f9e390cd6bf6fbd1449475d35e3346d72bdc).
# The GROMACS SOD.itp/CLA.itp of the same job give sigma = 0.251367 and
# 0.404468 nm, i.e. Rmin/2 = sigma * 2**(1/6) / 2 = 1.41075 and 2.27 Å.
# Keyed by element; the value is (CHARMM atom type, Rmin/2).
_CHARMM36_ION_RMIN_HALF: dict[str, tuple[str, float]] = {
    "LI": ("LIT", 1.2975), "NA": ("SOD", 1.41075), "MG": ("MG", 1.185), "K": ("POT", 1.76375),
    "CA": ("CAL", 1.367), "RB": ("RUB", 1.90), "CS": ("CES", 2.100), "BA": ("BAR", 1.890),
    "ZN": ("ZN", 1.09), "CD": ("CAD", 1.357), "CL": ("CLA", 2.27),
}

#: Label of the source of the default ion radii (:data:`ION_RADII`).
ION_RADII_SOURCE = "CHARMM36 toppar_water_ions.str"

#: Full citation of :data:`ION_RADII`, recorded in radius-set provenance.
ION_RADII_CITATION = (
    "CHARMM36 toppar_water_ions.str, NONBONDED Rmin/2 of LIT, SOD, MG, POT, CAL, RUB, CES, BAR, ZN, "
    "CAD, CLA (copy distributed by CHARMM-GUI v3.7, 25 Aug 2026; the file has no version string; "
    "SHA-256 41c6474470fa66e927198582a314f9e390cd6bf6fbd1449475d35e3346d72bdc). Ion parameters: "
    "Beglov and Roux, J. Chem. Phys. 100, 9050 (1994) and later Roux-group values cited in the file; "
    "Zn: Stote and Karplus, Proteins 23, 12 (1995)")

#: Default radii (Å) of recognised monatomic ions, keyed by element. Each value is
#: the CHARMM36 Lennard-Jones Rmin/2 of that ion type (``toppar_water_ions.str``). Applied
#: only to atoms that :func:`is_recognised_ion` accepts (an ion *residue* of
#: :data:`ION_NAME_ELEMENTS`); every other atom of the element keeps its
#: :data:`VDW_RADII` value. Ions the CHARMM36 file does not parameterise (Cu, F,
#: Br, I) keep :data:`VDW_RADII` as well.
ION_RADII: dict[str, float] = {element: value for element, (_, value) in _CHARMM36_ION_RMIN_HALF.items()}


def is_recognised_ion(atom_name: str, resname: str, element: str = "") -> bool:
    """Whether an atom is a recognised monatomic ion (an ion residue).

    An atom is a recognised ion when its ``(residue name, atom name)`` pair is
    in :data:`ION_NAME_ELEMENTS` *and* its element (explicit column if it holds
    a real symbol, otherwise inferred by :func:`infer_element`) is the element
    that table gives. Both names must match, so C-alpha ``CA`` of any amino
    acid, carbon ``CD`` of PRO or heme nitrogen ``NA`` are never ions, and a
    metal atom inside a larger residue (heme FE, a zinc in a ligand residue
    such as ``ZNH``) is never an ion either. A single-atom PDB residue whose
    residue and atom names are both the symbol (``ZN``/``ZN``, ``NA``/``NA``,
    ``CA``/``CA``, ``MG``/``MG``, ``K``/``K``, ``CL``/``CL``) *is* a
    recognised ion, whether it comes from an MD topology or a crystal
    structure.

    Parameters
    ----------
    atom_name : str
        Atom name as written in the file (case-insensitive).
    resname : str
        Residue name (case-insensitive).
    element : str, default ""
        Explicit element column, if any. An explicit element that differs from
        the table's element (for example element ``C`` on atom ``CA`` of
        residue ``CA``) makes the atom not an ion.

    Returns
    -------
    bool

    Examples
    --------
    >>> from crevice.radii import is_recognised_ion
    >>> is_recognised_ion("SOD", "SOD"), is_recognised_ion("CA", "ALA"), is_recognised_ion("ZN", "ZN", "ZN")
    (True, False, True)
    >>> is_recognised_ion("FE", "HEM", "FE")
    False
    """
    ion = ION_NAME_ELEMENTS.get((resname.strip().upper(), atom_name.strip().upper()))
    return ion is not None and infer_element(atom_name, element, resname) == ion


def ion_name_source(resname: str, atom_name: str) -> str | None:
    """Which force-field naming convention an ion name pair comes from.

    Parameters
    ----------
    resname, atom_name : str
        Residue and atom names (case-insensitive).

    Returns
    -------
    str or None
        ``"CHARMM36 toppar_water_ions.str"``, ``"GROMACS 2022.4 ions.itp"``
        and/or ``"AMBER legacy (Joung-Cheatham)"`` joined by ``"; "``, or
        ``None`` if the pair is not a recognised ion.

    Examples
    --------
    >>> from crevice.radii import ion_name_source
    >>> ion_name_source("SOD", "SOD")
    'CHARMM36 toppar_water_ions.str'
    >>> ion_name_source("ALA", "CA") is None
    True
    """
    key = (resname.strip().upper(), atom_name.strip().upper())
    sources = [label for label, table in (("CHARMM36 toppar_water_ions.str", _CHARMM36_IONS),
                                          ("GROMACS 2022.4 ions.itp", _GROMACS_IONS),
                                          ("AMBER legacy (Joung-Cheatham)", _AMBER_LEGACY_IONS))
               if key in table]
    return "; ".join(sources) or None


def is_element_symbol(text: str) -> bool:
    """Return ``True`` if the text is a periodic-table symbol (case-insensitive).

    Examples
    --------
    >>> from crevice.radii import is_element_symbol
    >>> is_element_symbol("Hg"), is_element_symbol("XX")
    (True, False)
    """
    return text.strip().upper() in ELEMENT_SYMBOLS


def infer_element(atom_name: str, explicit_element: str = "", resname: str = "") -> str:
    """Infer an element symbol from an atom record.

    Rules, applied in order:

    1. If ``explicit_element`` (with any charge digits or signs removed, for
       example ``"FE2+"``) is a periodic-table symbol, return it in upper case,
       whatever the atom name. An explicit ``HG`` is mercury, never hydrogen,
       and an explicit ``CA`` is calcium even on a residue called ALA.
       Blank, ``"."``, ``"?"`` or a value that is not an element symbol falls
       through to the name rules below.
    2. A recognised monatomic MD ion, matched on the residue name *and* the
       atom name together (:data:`ION_NAME_ELEMENTS`): CHARMM36 ``SOD``/``SOD``
       is Na, ``POT`` K, ``CLA`` Cl, ``CAL`` Ca, ``MG`` Mg, ``LIT`` Li,
       ``CES`` Cs, ``ZN2``/``ZN`` Zn, ``BAR`` Ba, ``RUB`` Rb, ``CD2``/``CD``
       Cd; GROMACS/AMBER ``NA``, ``K``, ``CL``, ``CA``, ``MG``, ``ZN``,
       ``LI``, ``RB``, ``CS``, ``CU``, ``F``, ``BR``, ``I`` (residue and atom
       name equal), ``NA+``, ``K+``, ``CL-`` and the other legacy AMBER names.
       ``CA`` is calcium only in a residue named ``CA``; C-alpha is never
       affected.
    3. Strip the atom name, upper-case it and drop leading digits (so
       ``"1HB"`` becomes ``"HB"``). An empty name gives ``"C"``.
    4. Standard protein carbon names (``CA``, ``CB``, ``CD1``, ...) give
       ``"C"``, except when the residue name is ``"CA"`` (a calcium ion).
    5. Names starting with ``CL``, ``BR``, ``FE``, ``MG``, ``MN``, ``ZN``,
       ``LI``, ``CU`` or ``SE`` give that two-letter symbol.
    6. The ambiguous symbols ``CA``, ``NA``, ``CO`` and ``NI`` are returned
       only when the atom name *and* the residue name both equal the symbol
       (a single-atom ion such as ``CA`` in residue ``CA``). Otherwise, for
       example ligand carbon ``CAB`` or heme nitrogen ``NA``, rule 7 applies.
    7. Otherwise the first letter of the name.

    Parameters
    ----------
    atom_name : str
        Atom name as written in the file.
    explicit_element : str, default ""
        Element column from the file, if any.
    resname : str, default ""
        Residue name; used to recognise monatomic ions (rule 2) and single-atom
        ions named like carbon or nitrogen atoms (rules 4 and 6).

    Returns
    -------
    str
        Upper-case element symbol.

    Notes
    -----
    Rules 2-7 apply only to files without a usable element column (for
    example GRO files); an explicit element always wins. Rule 2 knows only the
    ion names listed in :data:`ION_NAME_ELEMENTS` (sources in
    :func:`ion_name_source`); any other ion residue falls through to the name
    heuristics (rules 3-7), which can mistake it for a light element. Before
    rule 2 was added, CHARMM ``SOD`` and ``POT`` were read as S and P.
    Supply element columns where possible.
    Elements without a tabulated radius keep their symbol and use
    :data:`DEFAULT_RADIUS` (see :func:`default_radius_elements`).

    Examples
    --------
    >>> from crevice.radii import infer_element
    >>> infer_element("1HB2"), infer_element("OG1"), infer_element("CL1")
    ('H', 'O', 'CL')
    >>> infer_element("X", explicit_element="zn"), infer_element("HG", explicit_element="HG")
    ('ZN', 'HG')
    >>> infer_element("CA", resname="ALA"), infer_element("CAB", resname="HEM")
    ('C', 'C')
    >>> infer_element("SOD", resname="SOD"), infer_element("POT", resname="POT")
    ('NA', 'K')
    >>> infer_element("CD", resname="CD2"), infer_element("CD", resname="PRO")
    ('CD', 'C')
    """

    explicit = "".join(ch for ch in explicit_element.strip().upper() if ch.isalpha())
    if explicit in ELEMENT_SYMBOLS:
        return explicit

    name = atom_name.strip().upper()
    residue = resname.strip().upper()
    ion = ION_NAME_ELEMENTS.get((residue, name))
    if ion is not None:
        return ion
    if not name:
        return "C"
    while name and name[0].isdigit():
        name = name[1:]
    if name in PROTEIN_CARBON_NAMES and residue != "CA":
        return "C"
    prefix = name[:2]
    if len(name) >= 2 and prefix in UNAMBIGUOUS_NAME_ELEMENTS:
        return prefix
    if name in AMBIGUOUS_NAME_ELEMENTS and residue == name:
        return name
    return name[0] if name else "C"


def radius_source(element: str, resname: str = "", atom_name: str = "") -> str:
    """Where the default radius of an element (or of one atom) comes from.

    With only ``element`` this describes the :data:`VDW_RADII` entry. With
    ``resname`` and ``atom_name`` as well, it describes the radius the default
    set (:data:`DEFAULT_RADII`) gives that atom: a recognised ion residue
    (:func:`is_recognised_ion`) of an element in :data:`ION_RADII` reports the
    CHARMM36 source; anything else reports the element's table source.

    Parameters
    ----------
    element : str
        Element symbol (case-insensitive).
    resname, atom_name : str, default ""
        Residue and atom names of a specific atom (optional).

    Returns
    -------
    str
        ``"CHARMM36 toppar_water_ions.str (Rmin/2 of SOD)"`` (the CHARMM atom
        type in brackets) for a recognised ion, otherwise ``"Bondi 1964"``,
        ``"Mantina 2009"``, ``"legacy"`` (an original CREVICE value whose
        source was not recorded) or ``"default"`` (not tabulated;
        :data:`DEFAULT_RADIUS` is used).

    Examples
    --------
    >>> from crevice.radii import radius_source
    >>> radius_source("HG"), radius_source("GD")
    ('Bondi 1964', 'default')
    >>> radius_source("NA", "SOD", "SOD"), radius_source("C", "ALA", "CA")
    ('CHARMM36 toppar_water_ions.str (Rmin/2 of SOD)', 'Bondi 1964')
    """
    symbol = element.strip().upper()
    if atom_name and symbol in _CHARMM36_ION_RMIN_HALF and is_recognised_ion(atom_name, resname, symbol):
        return f"{ION_RADII_SOURCE} (Rmin/2 of {_CHARMM36_ION_RMIN_HALF[symbol][0]})"
    entry = _RADIUS_SOURCES.get(symbol)
    return entry[0] if entry else "default"


def default_radius_elements(atoms) -> list[str]:
    """Sorted elements of ``atoms`` that have no tabulated radius.

    These atoms are kept, with their element symbol, and use
    :data:`DEFAULT_RADIUS` (1.70 Å). The structure readers list them in their
    report as ``unsupported_radius_elements`` and warn once per load.

    Parameters
    ----------
    atoms : iterable of Atom

    Returns
    -------
    list of str
    """
    return sorted({atom.element.strip().upper() for atom in atoms} - VDW_RADII.keys())


#: Hydrogen isotopes, which get the ``H`` radius (unless a set lists them) and the
#: hydrogen filter of :meth:`crevice.models.StructureFrame.selected_atoms`
#: removes them like hydrogen.
HYDROGEN_ISOTOPES = frozenset({"D"})


class UnrecognisedElementWarning(UserWarning):
    """Atoms were excluded because the radius set in force gives them no radius."""


def partition_atoms_with_radii(atoms, radii: "RadiusSet | str | os.PathLike | None" = None):
    """Split atoms into those with a radius in a set and a tally of the rest.

    Parameters
    ----------
    atoms : iterable of Atom
    radii : RadiusSet, str, os.PathLike or None, default None
        The set to check (:func:`resolve_radii`); ``None`` uses the set in
        effect (:func:`effective_radii`).

    Returns
    -------
    kept : tuple of Atom
        Atoms with a radius (:meth:`RadiusSet.has_radius`), in input order.
    by_element : dict of str to int
        Excluded atoms per element symbol (``"?"`` when the element is
        unknown), sorted.
    by_name : dict of str to int
        Excluded atoms per ``"RESNAME:ATOMNAME"``, sorted; these are the keys a
        radius file can use for them.

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.radii import partition_atoms_with_radii
    >>> atoms = [Atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
    ...          Atom(2, "GD1", "XGD", "A", 2, 4.0, 0.0, 0.0, "GD")]
    >>> kept, by_element, by_name = partition_atoms_with_radii(atoms)
    >>> [a.name for a in kept], by_element, by_name
    (['CA'], {'GD': 1}, {'XGD:GD1': 1})
    """
    from collections import Counter

    radius_set = effective_radii(radii)
    kept, elements, names = [], Counter(), Counter()
    for atom in atoms:
        if radius_set.has_radius(atom.name, atom.resname, atom.element):
            kept.append(atom)
            continue
        elements[infer_element(atom.name, atom.element, atom.resname) or "?"] += 1
        names[f"{atom.resname.strip().upper()}:{atom.name.strip().upper()}"] += 1
    return tuple(kept), dict(sorted(elements.items())), dict(sorted(names.items()))


def unrecognised_element_message(by_element: Mapping[str, int], by_name: Mapping[str, int] | None = None,
                                 source: str = "", *, radii: "RadiusSet | str | os.PathLike | None" = None,
                                 excluded: bool = True) -> str:
    """Explain which atoms have no radius and how to supply one.

    Parameters
    ----------
    by_element : mapping of str to int
        Atoms without a radius per element.
    by_name : mapping of str to int, optional
        The same atoms per ``"RESNAME:ATOMNAME"`` (up to five are listed).
    source : str, default ""
        File name to prefix.
    radii : RadiusSet, str, os.PathLike or None, default None
        The set that was checked (``None``: the set in effect).
    excluded : bool, default True
        ``True`` for "were excluded" (structure files), ``False`` for "cannot
        be analysed" (MD selections, which are refused instead).

    Returns
    -------
    str
    """
    listing = ", ".join(f"{element} ({count} atom{'s' if count != 1 else ''})"
                        for element, count in sorted(by_element.items()))
    names = sorted(by_name or {})
    keys = ", ".join(names[:5]) + (", ..." if len(names) > 5 else "")
    radius_set = effective_radii(radii)
    outcome = "were excluded from the analysis" if excluded else "cannot be analysed"
    return ((f"{source}: " if source else "")
            + f"atoms with no van der Waals radius in radius set {radius_set.name!r} {outcome}: {listing}"
            + (f" [{keys}]" if keys else "") + ". "
            "To include them, give their radii in a radius file and pass it with --radii FILE "
            "(Python: radii=FILE, also when loading the structure). Start from "
            "crevice.radii.write_radius_template('radii.csv', elements=[...]); entries can be keyed by "
            "element, atom name or RESNAME:ATOMNAME, so designed or modified residues can be given "
            "their own radii.")


def write_radius_template(path, *, elements=(), base: str = "default") -> Path:
    """Write an editable CSV radius file that extends a built-in set.

    The file starts from ``base`` (``base,<preset>,``, so every built-in
    radius is kept) and lists commented example rows for ``elements`` and for
    atom-name keys. Remove the ``#`` and fill in radii (Å) for the atoms that
    need them, then pass the file with ``--radii FILE`` or ``radii=FILE``.

    Parameters
    ----------
    path : str or os.PathLike
        Output ``.csv`` file (overwritten).
    elements : iterable of str, default ()
        Element symbols to list as rows to fill in.
    base : str, default "default"
        Preset to start from (:func:`radii_preset`).

    Returns
    -------
    pathlib.Path

    Raises
    ------
    ValueError
        For an unknown preset or a path that does not end in ``.csv``.

    Examples
    --------
    >>> import pathlib, tempfile
    >>> from crevice.radii import load_radius_file, write_radius_template
    >>> path = write_radius_template(pathlib.Path(tempfile.mkdtemp()) / "radii.csv", elements=["GD"])
    >>> load_radius_file(path).has_radius("C", "ALA", "C")
    True
    """
    target = Path(path)
    if target.suffix.lower() != ".csv":
        raise ValueError("write_radius_template writes a .csv radius file")
    preset = radii_preset(base)
    lines = ["# CREVICE radius file (see the Atomic radii page of the documentation).",
             f"# Starts from the built-in {preset.name!r} set; rows below add or override radii in angstrom.",
             "# Keys: element (GD), atom name (GD1) or RESNAME:ATOMNAME (XGD:GD1, most specific wins).",
             "kind,key,radius", f"base,{preset.name},"]
    lines += [f"# element,{str(e).strip().upper()},<radius in angstrom>" for e in elements]
    lines += ["# atom,RESNAME:ATOMNAME,<radius in angstrom>"]
    target.write_text("\n".join(lines) + "\n")
    return target


def vdw_radius(element: str) -> float:
    """Van der Waals radius of an element from the standard table, in Å.

    This is the element table only. It does not know about ions: the radius
    the analyses use for an atom, including the CHARMM36 radius of a
    recognised ion residue, comes from :func:`atom_vdw_radius`.

    Parameters
    ----------
    element : str
        Element symbol (case-insensitive, surrounding whitespace ignored).

    Returns
    -------
    float
        The :data:`VDW_RADII` value, or :data:`DEFAULT_RADIUS` (1.70 Å) for an
        element not in the table.
    """

    return VDW_RADII.get(element.strip().upper(), DEFAULT_RADIUS)


def atom_vdw_radius(atom: Atom, radii: "RadiusSet | str | os.PathLike | None" = None) -> float:
    """Radius used for an atom in every clearance calculation, in Å.

    Every analysis obtains atomic radii from this function. It returns
    :meth:`RadiusSet.radius_for_atom` of ``radii`` or, if that is ``None``, of
    the set in effect (:func:`use_radii`, an analysis's ``radii=`` argument or
    the CLI ``--radii`` option). When no set was chosen anywhere, the set in
    effect is :data:`DEFAULT_RADII`: a recognised monatomic ion residue
    (:func:`is_recognised_ion`, for example CHARMM ``SOD``) gets its CHARMM36
    radius from :data:`ION_RADII` and every other atom gets
    ``vdw_radius(infer_element(atom.name, atom.element, atom.resname))``.

    Parameters
    ----------
    atom : Atom
        The atom; its ``name``, ``resname`` and ``element`` are used.
    radii : RadiusSet, str, os.PathLike or None, default None
        Radius set for this lookup (see :func:`resolve_radii`); ``None`` uses
        the set in effect.

    Returns
    -------
    float
        Radius in Å.

    Raises
    ------
    ValueError
        If the set has HOLE records and none matches the atom.

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.radii import atom_vdw_radius
    >>> sodium = Atom(1, "SOD", "SOD", "", 1, 0.0, 0.0, 0.0, "")
    >>> calpha = Atom(2, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "")
    >>> atom_vdw_radius(sodium), atom_vdw_radius(calpha), atom_vdw_radius(sodium, "bondi")
    (1.41075, 1.7, 2.27)
    """
    return effective_radii(radii).radius_for_names(atom.name, atom.resname, atom.element)


def radius_for_names(name: str, resname: str = "", element: str = "",
                     radii: "RadiusSet | str | os.PathLike | None" = None) -> float:
    """:func:`atom_vdw_radius` for an atom given by its names (no :class:`Atom` needed).

    Parameters
    ----------
    name, resname : str
        Atom and residue names.
    element : str, default ""
        Explicit element column or an already inferred symbol.
    radii : RadiusSet, str, os.PathLike or None, default None
        As for :func:`atom_vdw_radius`.

    Returns
    -------
    float
        Radius in Å.
    """
    return effective_radii(radii).radius_for_names(name, resname, element)


def residue_properties(resname: str) -> tuple[str, ...]:
    """Coarse chemical classes of a residue, from its name.

    A residue can belong to several classes. Membership sets:

    * ``hydrophobic``: ALA VAL LEU ILE MET PHE TRP PRO TYR
    * ``polar``: SER THR ASN GLN CYS TYR HIS TRP
    * ``positive``: ARG LYS HIS
    * ``negative``: ASP GLU
    * ``aromatic``: PHE TYR TRP HIS
    * ``small``: GLY ALA SER CYS

    Parameters
    ----------
    resname : str
        Residue name (case-insensitive). Only the standard three-letter names
        above are recognised; protonation variants such as ``HSD`` or ``HIE``
        are not mapped to their parent residue here.

    Returns
    -------
    tuple of str
        Classes in the order listed above, or ``("other",)`` when none apply
        (for example GLY is only ``small``; water and ligands are ``other``).

    Examples
    --------
    >>> from crevice.radii import residue_properties
    >>> residue_properties("LYS"), residue_properties("HOH")
    (('positive',), ('other',))
    """
    residue = resname.strip().upper()
    props: list[str] = []
    if residue in HYDROPHOBIC:
        props.append("hydrophobic")
    if residue in POLAR:
        props.append("polar")
    if residue in POSITIVE:
        props.append("positive")
    if residue in NEGATIVE:
        props.append("negative")
    if residue in AROMATIC:
        props.append("aromatic")
    if residue in SMALL:
        props.append("small")
    if not props:
        props.append("other")
    return tuple(props)

# ---------------------------------------------------------------------------
# Radius sets
# ---------------------------------------------------------------------------

#: Lookup order used by :meth:`RadiusSet.radius_for_atom`; also recorded in
#: provenance.
RADIUS_PRECEDENCE = ("RESNAME:ATOMNAME", "ATOMNAME", "recognised ion residue (ion_radii by element)",
                     "HOLE VDWR records (first match)", "element", "default")


def _hole_match(pattern: str, text: str) -> bool:
    """HOLE's ``LMATCH``: fixed width, case-insensitive, ``?`` in either string matches anything."""
    return all(p == "?" or t == "?" or p == t for p, t in zip(pattern, text))


def _canonical_sha256(value) -> str:
    import hashlib
    import json

    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True)
class RadiusSet:
    """A named table of atomic radii that the analyses can use.

    Every clearance, cast, section, contact distance and SASA in CREVICE is
    measured to atom surfaces, so the radius set changes the numbers. When no
    set is given anywhere, the analyses use :data:`DEFAULT_RADII` (the
    standard :data:`VDW_RADII` element table, plus CHARMM36 radii for
    recognised ion residues); pass a set as ``radii=`` to an analysis, enter
    :func:`use_radii`, or use ``--radii`` on the command line to use another
    one. Every output records the set it used (:meth:`provenance`).

    Lookup precedence for an atom (first hit wins):

    1. ``"RESNAME:ATOMNAME"`` in :attr:`atom_radii`;
    2. ``"ATOMNAME"`` in :attr:`atom_radii`;
    3. a recognised monatomic ion residue (:func:`is_recognised_ion`) whose
       element is in :attr:`ion_radii`;
    4. :attr:`hole_records`, in file order (HOLE ``VDWR`` semantics; see
       :func:`load_radius_file`). If the set has records and none matches, a
       ``ValueError`` is raised, as HOLE stops on an unmatched atom;
    5. the element (:func:`infer_element`) in :attr:`element_radii`;
    6. :attr:`default_radius`.

    Lookups are cached per ``(atom name, residue name, element)``, so a set
    costs one dictionary lookup per atom after the first occurrence of a name.

    Attributes
    ----------
    name : str
        Identifier of the set, recorded in provenance.
    element_radii : mapping of str to float
        Element symbol (upper case) to radius in Å.
    atom_radii : mapping of str to float
        Atom-name or ``"RESNAME:ATOMNAME"`` keys (upper case) to radius in Å.
    default_radius : float, default 1.70
        Radius in Å for anything not matched.
    hole_records : tuple of (str, str, float)
        Ordered HOLE ``VDWR`` records ``(atom pattern, residue pattern,
        radius)``: a 4-character atom-name pattern and a 3-character residue
        pattern in which ``?`` matches any character.
    ion_radii : mapping of str to float
        Element symbol (upper case) to radius in Å, applied only to recognised
        monatomic ion residues (:func:`is_recognised_ion`). Empty means ions
        are treated like any other atom of their element.
    source : str
        Where the values come from (preset description or absolute file path).
    source_sha256 : str or None
        SHA-256 of the radius file the set was read from, if any.
    citation : str
        Literature source of the values, if known.

    See Also
    --------
    radii_preset, custom_radii_set, load_radius_file, use_radii

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.radii import radii_preset
    >>> atom = Atom(1, "CB", "ALA", "A", 1, 0.0, 0.0, 0.0, "C")
    >>> radii_preset("bondi").radius_for_atom(atom), radii_preset("hole").radius_for_atom(atom)
    (1.7, 1.85)
    >>> sodium = Atom(2, "SOD", "SOD", "", 1, 0.0, 0.0, 0.0, "")
    >>> radii_preset("default").radius_for_atom(sodium), radii_preset("bondi").radius_for_atom(sodium)
    (1.41075, 2.27)
    """

    name: str
    element_radii: Mapping[str, float] = field(default_factory=dict)
    atom_radii: Mapping[str, float] = field(default_factory=dict)
    default_radius: float = DEFAULT_RADIUS
    hole_records: tuple = ()
    source: str = ""
    source_sha256: str | None = None
    citation: str = ""
    ion_radii: Mapping[str, float] = field(default_factory=dict)
    _cache: dict = field(default_factory=dict, init=False, repr=False, compare=False)

    def radius_for_element(self, element: str) -> float:
        """Radius in Å for an element symbol.

        Parameters
        ----------
        element : str
            Case-insensitive symbol.

        Returns
        -------
        float
            The table value, or :attr:`default_radius` if absent. Deuterium
            (``D``) uses the ``H`` entry unless the table lists ``D`` itself.
        """
        symbol = element.strip().upper()
        if symbol not in self.element_radii and symbol in HYDROGEN_ISOTOPES and "H" in self.element_radii:
            symbol = "H"
        return self.element_radii.get(symbol, self.default_radius)

    def radius_for_names(self, name: str, resname: str = "", element: str = "") -> float:
        """Radius in Å for an atom given by its names.

        Parameters
        ----------
        name : str
            Atom name.
        resname : str, default ""
            Residue name.
        element : str, default ""
            Explicit element column, if any (otherwise inferred from the names
            by :func:`infer_element`).

        Returns
        -------
        float

        Raises
        ------
        ValueError
            If the set has HOLE records and none matches the atom.
        """
        key = (name, resname, element)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        atom_key = name.strip().upper()
        residue = resname.strip().upper()
        exact_key = f"{residue}:{atom_key}"
        if exact_key in self.atom_radii:
            value = self.atom_radii[exact_key]
        elif atom_key in self.atom_radii:
            value = self.atom_radii[atom_key]
        elif self.ion_radii and self._ion_element(name, resname, element) in self.ion_radii:
            value = self.ion_radii[self._ion_element(name, resname, element)]
        elif self.hole_records:
            atom_field, residue_field = atom_key.ljust(4)[:4], residue.ljust(3)[:3]
            for atom_pattern, residue_pattern, radius in self.hole_records:
                if _hole_match(atom_pattern, atom_field) and _hole_match(residue_pattern, residue_field):
                    value = radius
                    break
            else:
                raise ValueError(
                    f"Radius set {self.name!r}: no VDWR record matches atom {atom_key!r} in residue "
                    f"{residue!r} (HOLE stops on the same atom). Add a specific record, or a final "
                    f"catch-all such as 'VDWR ???? ??? 1.70'.")
        else:
            value = self.radius_for_element(infer_element(name, element, resname))
        self._cache[key] = float(value)
        return float(value)

    @staticmethod
    def _ion_element(name: str, resname: str, element: str) -> str | None:
        """Element of a recognised ion residue, or ``None`` for any other atom."""
        ion = ION_NAME_ELEMENTS.get((resname.strip().upper(), name.strip().upper()))
        if ion is None or infer_element(name, element, resname) != ion:
            return None
        return ion

    def radius_for_atom(self, atom: Atom) -> float:
        """Radius in Å for an atom, using the precedence described on the class.

        Parameters
        ----------
        atom : Atom

        Returns
        -------
        float
        """
        return self.radius_for_names(atom.name, atom.resname, atom.element)

    def has_radius(self, name: str, resname: str = "", element: str = "") -> bool:
        """Whether this set gives the atom a radius of its own.

        True when one of the first five precedence rules on the class matches:
        an ``RESNAME:ATOMNAME`` or ``ATOMNAME`` entry, a recognised ion residue
        in :attr:`ion_radii`, a matching HOLE record, or the atom's element
        (deuterium counts as hydrogen) in :attr:`element_radii`. The catch-all
        :attr:`default_radius` does not count: the structure readers exclude
        atoms without a radius (:func:`partition_atoms_with_radii`).

        Parameters
        ----------
        name, resname : str
            Atom and residue names.
        element : str, default ""
            Explicit element column, if any.

        Returns
        -------
        bool

        Examples
        --------
        >>> from crevice.radii import radii_preset, custom_radii_set
        >>> radii_preset("default").has_radius("GD1", "XGD", "GD")
        False
        >>> custom_radii_set(base="default", atom_radii={"XGD:GD1": 2.0}).has_radius("GD1", "XGD", "GD")
        True
        """
        key = ("__has__", name, resname, element)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        atom_key = name.strip().upper()
        residue = resname.strip().upper()
        if f"{residue}:{atom_key}" in self.atom_radii or atom_key in self.atom_radii:
            found = True
        elif self.ion_radii and self._ion_element(name, resname, element) in self.ion_radii:
            found = True
        elif self.hole_records:
            atom_field, residue_field = atom_key.ljust(4)[:4], residue.ljust(3)[:3]
            found = any(_hole_match(a, atom_field) and _hole_match(r, residue_field)
                        for a, r, _ in self.hole_records)
        else:
            symbol = infer_element(name, element, resname)
            found = symbol in self.element_radii or (symbol in HYDROGEN_ISOTOPES and "H" in self.element_radii)
        self._cache[key] = found
        return found

    @property
    def table_sha256(self) -> str:
        """SHA-256 of the effective table (all radii, records and the default)."""
        table = {
            "element_radii": {k: float(v) for k, v in self.element_radii.items()},
            "atom_radii": {k: float(v) for k, v in self.atom_radii.items()},
            "hole_records": [list(r) for r in self.hole_records],
            "default_radius": float(self.default_radius),
        }
        if self.ion_radii:
            # Only present when used, so sets without ion radii keep the hash
            # they had before ion radii existed.
            table["ion_radii"] = {k: float(v) for k, v in self.ion_radii.items()}
        return _canonical_sha256(table)

    def provenance(self) -> dict:
        """What to record about this set in an output manifest.

        Every CREVICE output records this dictionary (as ``radii``), including
        runs with the default set. It identifies the set completely:
        :func:`radius_set_from_provenance` rebuilds the set from it and
        :attr:`table_sha256` verifies the rebuild.

        Returns
        -------
        dict
            ``name``, ``source``, ``source_sha256`` (radius file, or ``None``),
            ``citation``, ``table_sha256``, ``default_radius_A``,
            ``element_overrides_A`` (element radii that differ from, or are
            not in, the standard table), ``standard_elements_missing`` (standard
            elements this set does not list, which get the default radius;
            empty for a set of HOLE records), ``atom_radii_A``,
            ``ion_radii_A`` (radii of recognised ion residues by element;
            empty when ions are treated as ordinary atoms), ``ion_scope`` (the
            rule that decides which atoms are ions), ``hole_records`` and
            ``precedence``.
        """
        cached = self._cache.get("__provenance__")
        if cached is None:
            elements = {k: float(v) for k, v in sorted(self.element_radii.items())}
            cached = {
                "name": self.name,
                "source": self.source or "custom",
                "source_sha256": self.source_sha256,
                "citation": self.citation or None,
                "table_sha256": self.table_sha256,
                "default_radius_A": float(self.default_radius),
                "element_overrides_A": {k: v for k, v in elements.items() if VDW_RADII.get(k) != v},
                "standard_elements_missing": sorted(VDW_RADII.keys() - elements.keys()) if not self.hole_records else [],
                "atom_radii_A": {k: float(v) for k, v in sorted(self.atom_radii.items())},
                "ion_radii_A": {k: float(v) for k, v in sorted(self.ion_radii.items())},
                "ion_scope": ION_SCOPE if self.ion_radii else None,
                "hole_records": [[a, r, float(v)] for a, r, v in self.hole_records],
                "precedence": list(RADIUS_PRECEDENCE),
            }
            self._cache["__provenance__"] = cached
        import copy

        return copy.deepcopy(cached)


# HOLE 2.3.1 rad/simple.rad, VDWR records in file order (remark and BOND lines
# omitted). Its header reads "van der Waals radii: AMBER united atom, from
# Weiner et al. (1984), JACS, vol 106 pp765-768". Checked against the copy in
# the HOLE 2.3.1 source tree (github.com/osmart/hole2, tag v2.3.1).
_HOLE_SIMPLE_RAD = (
    ("C???", "???", 1.85),
    ("O???", "???", 1.65),
    ("S???", "???", 2.00),
    ("N???", "???", 1.75),
    ("H???", "???", 1.00),
    ("P???", "???", 2.10),
    ("E2? ", "GLN", 1.00),
    ("D2? ", "ASN", 1.00),
    ("LP??", "???", 0.00),
)

_BONDI_CITATION = ("Bondi, J. Phys. Chem. 68, 441 (1964); Mantina et al., J. Phys. Chem. A 113, 5806 "
                   "(2009) for main-group elements Bondi lacked; Fe, Mn, Co legacy (see radius_source)")

#: The rule that decides which atoms get :attr:`RadiusSet.ion_radii`; recorded
#: in provenance.
ION_SCOPE = ("recognised monatomic ion residues only: (residue name, atom name) in "
             "crevice.radii.ION_NAME_ELEMENTS and the atom's element equal to that entry's element; "
             "metals inside other residues, C-alpha CA, PRO CD and heme NA are not ions")

#: The radius set used when none is chosen. It combines the standard
#: :data:`VDW_RADII` element table with :data:`ION_RADII` (CHARMM36 Rmin/2) for recognised ion
#: residues. Also available as the preset ``"default"``.
DEFAULT_RADII = RadiusSet(
    "default", VDW_RADII, ion_radii=ION_RADII,
    source="crevice.radii.VDW_RADII (the standard table) with crevice.radii.ION_RADII for recognised ion residues",
    citation=f"{_BONDI_CITATION}. Ion residues: {ION_RADII_CITATION}")

#: Built-in radius tables; see :func:`radii_preset`.
RADIUS_PRESETS: dict[str, RadiusSet] = {
    "default": DEFAULT_RADII,
    "bondi": RadiusSet("bondi", VDW_RADII, source="crevice.radii.VDW_RADII (the standard table) for every atom, "
                                                  "ions included",
                       citation=_BONDI_CITATION),
    "hole": RadiusSet("hole", hole_records=_HOLE_SIMPLE_RAD, source="HOLE 2.3.1 rad/simple.rad",
                      citation="HOLE: Smart et al., J. Mol. Graph. 14, 354 (1996); simple.rad values: AMBER "
                               "united-atom radii, Weiner et al., J. Am. Chem. Soc. 106, 765 (1984)"),
    "charmm_like": RadiusSet("charmm_like", {**VDW_RADII, "C": 2.00, "N": 1.85, "O": 1.70, "S": 2.00},
                             ion_radii=ION_RADII,
                             source="CREVICE approximation: VDW_RADII with larger C, N, O, S (not read "
                                    "from a CHARMM parameter file); ION_RADII for recognised ion residues",
                             citation=f"Ion residues: {ION_RADII_CITATION}"),
}


def radii_preset(name: str = "default") -> RadiusSet:
    """Return one of the built-in radius tables.

    Parameters
    ----------
    name : {"default", "bondi", "hole", "charmm_like"}, default "default"
        Case-insensitive preset name.

        * ``"default"``: :data:`DEFAULT_RADII`, the set used when none is
          chosen. Every atom gets its :data:`VDW_RADII` element radius
          (Bondi 1964; Mantina et al. 2009 where Bondi has none), except
          recognised monatomic ion residues (:func:`is_recognised_ion`), which
          get the CHARMM36 Lennard-Jones Rmin/2 of :data:`ION_RADII`: Li
          1.2975, Na 1.41075, Mg 1.185, K 1.76375, Ca 1.367, Rb 1.90, Cs 2.10,
          Ba 1.89, Zn 1.09, Cd 1.357, Cl 2.27 Å. Choosing it explicitly gives
          the same numbers and the same provenance as choosing nothing.
        * ``"bondi"``: :data:`VDW_RADII` for every atom, ions included (a
          CHARMM ``SOD`` is a neutral-atom sodium, 2.27 Å). This was the
          default between the ion-name table (30 Sep 2026) and the CHARMM36
          ion radii; it differs from ``"default"`` only for ion residues.
        * ``"hole"``: HOLE 2.3.1's default radius file ``simple.rad`` (AMBER
          united-atom radii of Weiner et al. 1984): any atom whose *name*
          starts with C 1.85 Å, O 1.65, S 2.00, N 1.75, H 1.00, P 2.10 Å;
          GLN ``E2?`` and ASN ``D2?`` 1.00 Å; ``LP??`` 0.00 Å. Matching is by
          atom name as in HOLE, not by element (a CHARMM ``SOD`` ion is
          2.00 Å, ``CLA`` 1.85 Å), and an atom no record matches (for example
          ``ZN``) raises ``ValueError`` as HOLE does.
        * ``"charmm_like"``: :data:`VDW_RADII` with larger C, N, O and S radii
          (2.00, 1.85, 1.70 and 2.00 Å) and the :data:`ION_RADII` of the
          default set for ion residues. The C, N, O, S values are an
          approximation, not read from a force-field parameter file.

    Returns
    -------
    RadiusSet
        The shared preset object (immutable).

    Raises
    ------
    ValueError
        For an unknown preset name.

    Examples
    --------
    >>> from crevice.radii import radii_preset
    >>> radii_preset("hole").radius_for_names("OG", "SER")
    1.65
    >>> radii_preset().radius_for_names("CLA", "CLA"), radii_preset("bondi").radius_for_names("CLA", "CLA")
    (2.27, 1.75)
    """

    key = name.strip().lower()
    try:
        return RADIUS_PRESETS[key]
    except KeyError as exc:
        valid = ", ".join(sorted(RADIUS_PRESETS))
        raise ValueError(f"Unknown radii preset {name!r}; valid presets are {valid}") from exc


def _radius_value(value, label: str) -> float:
    import math

    try:
        radius = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{label}: radius {value!r} is not a number") from None
    if not math.isfinite(radius) or radius < 0:
        raise ValueError(f"{label}: radius must be finite and non-negative, got {value!r}")
    return radius


def _element_key(key: str, label: str) -> str:
    symbol = str(key).strip().upper()
    if symbol not in ELEMENT_SYMBOLS:
        raise ValueError(f"{label}: {key!r} is not an element symbol")
    return symbol


def _atom_key(key: str, label: str) -> str:
    text = str(key).strip().upper()
    parts = text.split(":")
    if not text or len(parts) > 2 or not all(part.strip() for part in parts):
        raise ValueError(f"{label}: atom key {key!r} must be ATOMNAME or RESNAME:ATOMNAME")
    return ":".join(part.strip() for part in parts)


def custom_radii_set(
    *,
    name: str = "user",
    element_radii: Mapping[str, float] | None = None,
    atom_radii: Mapping[str, float] | None = None,
    default_radius: float = DEFAULT_RADIUS,
    base: str | RadiusSet | None = None,
    source: str = "",
    citation: str = "",
    ion_radii: Mapping[str, float] | None = None,
) -> RadiusSet:
    """Build a :class:`RadiusSet` from user-supplied radii.

    Keys are stripped and upper-cased; values must be finite and non-negative.

    Parameters
    ----------
    name : str, default "user"
        Identifier recorded in provenance.
    element_radii : mapping of str to float, optional
        Element symbol to radius in Å.
    atom_radii : mapping of str to float, optional
        ``"ATOMNAME"`` or ``"RESNAME:ATOMNAME"`` to radius in Å. These win
        over everything else, including ion radii.
    default_radius : float, default 1.70
        Radius in Å for anything not matched.
    base : str or RadiusSet, optional
        Preset name or set to start from; its element radii, atom radii, ion
        radii and HOLE records are kept and the entries given here override
        them. Without a base, only the entries given here (and
        ``default_radius``) exist, so every element not listed gets
        ``default_radius`` and ions are ordinary atoms unless ``ion_radii`` is
        given. ``base="default"`` keeps the CHARMM36 ion radii;
        ``base="bondi"`` treats ions as neutral atoms.
    source, citation : str, optional
        Recorded in provenance.
    ion_radii : mapping of str to float, optional
        Element symbol to radius in Å for recognised ion residues
        (:func:`is_recognised_ion`), for example ``{"NA": 1.02}`` for an
        ionic (Shannon) sodium radius.

    Returns
    -------
    RadiusSet

    Raises
    ------
    ValueError
        For an invalid key or radius.

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.radii import custom_radii_set
    >>> radii = custom_radii_set(element_radii={"c": 1.9}, atom_radii={"ser:og": 1.6})
    >>> radii.radius_for_atom(Atom(1, "OG", "SER", "A", 1, 0.0, 0.0, 0.0, "O"))
    1.6
    >>> radii.radius_for_element("N")
    1.7
    >>> custom_radii_set(base="bondi", atom_radii={"SOD:SOD": 1.02}).radius_for_names("SOD", "SOD")
    1.02
    >>> custom_radii_set(base="default", ion_radii={"na": 1.02}).radius_for_names("NA", "NA")
    1.02
    """

    start = radii_preset(base) if isinstance(base, str) else base
    elements = dict(start.element_radii) if start is not None else {}
    atoms = dict(start.atom_radii) if start is not None else {}
    elements.update({_element_key(k, name): _radius_value(v, f"{name} element {k}")
                     for k, v in dict(element_radii or {}).items()})
    atoms.update({_atom_key(k, name): _radius_value(v, f"{name} atom {k}")
                  for k, v in dict(atom_radii or {}).items()})
    ions = dict(start.ion_radii) if start is not None else {}
    ions.update({_element_key(k, name): _radius_value(v, f"{name} ion {k}")
                 for k, v in dict(ion_radii or {}).items()})
    return RadiusSet(
        name=name,
        element_radii=elements,
        atom_radii=atoms,
        default_radius=_radius_value(default_radius, f"{name} default"),
        hole_records=start.hole_records if start is not None else (),
        source=source or (f"custom, based on {start.name}" if start is not None else "custom"),
        citation=citation or (start.citation if start is not None else ""),
        ion_radii=ions,
    )


def _file_sha256(path) -> str:
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_hole_radius_file(path) -> tuple:
    """Read the ``VDWR`` records of a HOLE radius file (``.rad``).

    The file is read exactly as HOLE 2.3.1 reads it (``tsradr.f``): each line
    is upper-cased; a line starting ``VDWR`` in columns 1-4 gives the atom-name
    pattern in columns 6-9, the residue-name pattern in columns 11-13 and the
    radius in columns 14-23. ``BOND`` records (bond radii for HOLE's other
    options), ``remark`` lines and anything else are ignored. Records are kept
    in file order because the first matching record wins.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    tuple of (str, str, float)
        ``(atom pattern, residue pattern, radius)`` per ``VDWR`` record.

    Raises
    ------
    ValueError
        If the file has no ``VDWR`` record, a record contains a tab, or a
        radius field is blank or has no decimal point (HOLE's ``F10.3`` read
        would silently give zero or divide it by 1000).
    """
    records = []
    for number, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.upper()
        if not line.startswith("VDWR"):
            continue
        if "\t" in line:
            raise ValueError(f"{path}:{number}: tab in a VDWR record; HOLE reads fixed columns")
        line = line.ljust(23)
        value = line[13:23].strip()
        if not value or "." not in value:
            raise ValueError(f"{path}:{number}: radius field (columns 14-23) {value!r} needs a decimal "
                             "point; HOLE's F10.3 read would change its value")
        records.append((line[5:9], line[10:13], _radius_value(value, f"{path}:{number}")))
    if not records:
        raise ValueError(f"{path}: no VDWR records found")
    return tuple(records)


def load_radius_file(path) -> RadiusSet:
    """Read a radius set from a JSON, CSV or HOLE ``.rad`` file.

    The format is chosen by the extension (case-insensitive):

    ``.json``
        An object with optional keys ``name`` (default: file stem),
        ``base`` (a preset name to start from; without it only the listed
        radii exist), ``default_radius`` (Å, default 1.70),
        ``element_radii`` (``{"C": 1.7, ...}``), ``atom_radii``
        (``{"OW": 1.5, "SOD:SOD": 1.02, ...}``), ``ion_radii``
        (``{"NA": 1.02, ...}``, applied to recognised ion residues) and
        ``citation``. Other keys are an error.
    ``.csv``
        A header row containing ``kind,key,radius`` and one row per entry:
        ``element,C,1.70``; ``atom,SOD:SOD,1.02`` (``ATOMNAME`` or
        ``RESNAME:ATOMNAME``); ``ion,NA,1.02`` (element of a recognised ion
        residue); ``default,,1.70``; ``base,default,`` (radius empty).
        Blank lines and lines starting ``#`` are ignored.
    ``.rad``
        A HOLE radius file; only ``VDWR`` records are used, with HOLE's
        column layout and first-match wildcard rules
        (:func:`read_hole_radius_file`).

    Atom names are compared after stripping blanks and upper-casing. For
    ``.rad`` records CREVICE left-justifies the atom name in the 4-character
    field; HOLE instead takes PDB columns 14-17 (prefixing column 13 only when
    it holds ``H``). The two agree for every standard name of up to three
    characters and for 4-character hydrogen names; they differ for other
    4-character names, names with a leading digit (``1HB2``) and two-letter
    element names written from column 13 (HOLE reads PDB ``ZN`` as ``N``).
    Residue names are truncated to three characters as in HOLE.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    RadiusSet
        With ``source`` set to the absolute path and ``source_sha256`` to the
        file's SHA-256, both recorded in output provenance.

    Raises
    ------
    ValueError
        For an unsupported extension or invalid content.
    """
    import csv
    import json

    path = Path(path).expanduser().resolve(strict=True)
    suffix = path.suffix.lower()
    digest = _file_sha256(path)
    if suffix == ".rad":
        return RadiusSet(path.stem, hole_records=read_hole_radius_file(path), source=str(path),
                         source_sha256=digest, citation="")
    if suffix == ".json":
        data = json.loads(path.read_text())
        if not isinstance(data, dict):
            raise ValueError(f"{path}: expected a JSON object")
        unknown = set(data) - {"name", "base", "default_radius", "element_radii", "atom_radii", "ion_radii",
                               "citation"}
        if unknown:
            raise ValueError(f"{path}: unknown keys {sorted(unknown)}")
        for key in ("element_radii", "atom_radii", "ion_radii"):
            if not isinstance(data.get(key, {}), dict):
                raise ValueError(f"{path}: {key} must be an object")
        spec = dict(name=str(data.get("name") or path.stem), element_radii=data.get("element_radii"),
                    atom_radii=data.get("atom_radii"), default_radius=data.get("default_radius", DEFAULT_RADIUS),
                    base=data.get("base"), citation=str(data.get("citation") or ""), ion_radii=data.get("ion_radii"))
    elif suffix == ".csv":
        lines = [line for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
        rows = list(csv.DictReader(lines))
        if not lines or not {"kind", "key", "radius"} <= {h.strip().lower() for h in (csv.DictReader(lines).fieldnames or [])}:
            raise ValueError(f"{path}: CSV needs a header with kind,key,radius")
        elements, atoms, ions, base, default = {}, {}, {}, None, DEFAULT_RADIUS
        for number, row in enumerate(rows, 2):
            row = {str(k).strip().lower(): (v or "").strip() for k, v in row.items() if k is not None}
            kind = row["kind"].lower()
            if kind == "element":
                elements[row["key"]] = row["radius"]
            elif kind == "atom":
                atoms[row["key"]] = row["radius"]
            elif kind == "ion":
                ions[row["key"]] = row["radius"]
            elif kind == "default":
                default = row["radius"]
            elif kind == "base":
                base = row["key"]
            else:
                raise ValueError(f"{path}: row {number}: kind must be element, atom, ion, default or base")
        spec = dict(name=path.stem, element_radii=elements, atom_radii=atoms, default_radius=default, base=base,
                    citation="", ion_radii=ions)
    else:
        raise ValueError(f"{path}: radius files must be .json, .csv or .rad")
    radii = custom_radii_set(source=str(path), **spec)
    return replace(radii, source=str(path), source_sha256=digest)


def resolve_radii(radii: "RadiusSet | str | os.PathLike | None") -> RadiusSet | None:
    """Turn a ``radii=`` argument into a :class:`RadiusSet`.

    Parameters
    ----------
    radii : RadiusSet, str, os.PathLike or None
        ``None`` stays ``None`` (no set chosen: the enclosing set, by default
        :data:`DEFAULT_RADII`, applies); a :class:`RadiusSet` is returned
        unchanged; a preset name (``default``, ``bondi``, ``hole``,
        ``charmm_like``, case-insensitive) gives :func:`radii_preset`; any other
        string or path must be a radius file (:func:`load_radius_file`). A
        preset name wins over a file of the same name; write ``./hole`` for
        the file.

    Returns
    -------
    RadiusSet or None
        ``None`` only for ``radii=None``.

    Raises
    ------
    ValueError
        For an unknown name that is not an existing file.
    """
    if radii is None or isinstance(radii, RadiusSet):
        return radii
    if isinstance(radii, str) and radii.strip().lower() in RADIUS_PRESETS:
        return radii_preset(radii)
    path = Path(radii).expanduser()
    if path.is_file():
        return load_radius_file(path)
    valid = ", ".join(sorted(RADIUS_PRESETS))
    raise ValueError(f"Unknown radius set {str(radii)!r}: use a preset ({valid}) or the path of a "
                     ".json, .csv or HOLE .rad radius file")


# ---------------------------------------------------------------------------
# The radius set in effect
# ---------------------------------------------------------------------------

_ACTIVE_RADII: contextvars.ContextVar = contextvars.ContextVar("crevice_radii", default=None)


def active_radii() -> RadiusSet | None:
    """The radius set explicitly chosen here, or ``None`` if none was chosen.

    Returns
    -------
    RadiusSet or None
        The set of the innermost :func:`use_radii` block (or ``radii=``
        argument, or ``--radii`` option) of this thread or process. ``None``
        means no set was chosen, so :data:`DEFAULT_RADII` applies; use
        :func:`effective_radii` for the set that is actually used.
    """
    return _ACTIVE_RADII.get()


def effective_radii(radii: "RadiusSet | str | os.PathLike | None" = None) -> RadiusSet:
    """The radius set the analyses use here.

    Parameters
    ----------
    radii : RadiusSet, str, os.PathLike or None, default None
        A set to resolve (:func:`resolve_radii`); ``None`` means the set in
        effect.

    Returns
    -------
    RadiusSet
        ``radii`` if given, otherwise the set chosen by the innermost
        :func:`use_radii` block, otherwise :data:`DEFAULT_RADII`. Never
        ``None``.

    Examples
    --------
    >>> from crevice.radii import effective_radii, use_radii
    >>> effective_radii().name
    'default'
    >>> with use_radii("hole"):
    ...     effective_radii().name
    'hole'
    """
    chosen = resolve_radii(radii) if radii is not None else _ACTIVE_RADII.get()
    return DEFAULT_RADII if chosen is None else chosen


@contextlib.contextmanager
def use_radii(radii: "RadiusSet | str | os.PathLike | None"):
    """Use a radius set for every analysis run inside a ``with`` block.

    This is the one mechanism behind every ``radii=`` argument and the CLI
    ``--radii`` option: :func:`atom_vdw_radius`, which all analyses call,
    reads the set in effect. The setting is a :mod:`contextvars` variable, so
    it is private to the current thread (and asyncio task) and is restored when
    the block ends, even on error.

    Worker processes: ``crevice cavity-trajectory``/``region-trajectory``
    worker pools receive the set explicitly (as an initialiser argument) and
    enter it in each worker, so forked and spawned workers use the same set.
    Threads you start yourself do not inherit it; pass ``radii=`` to the
    analysis inside the thread, or run it with ``contextvars.copy_context()``.

    Parameters
    ----------
    radii : RadiusSet, str, os.PathLike or None
        See :func:`resolve_radii`. ``None`` changes nothing: the enclosing
        setting (by default :data:`DEFAULT_RADII`) stays in effect.

    Yields
    ------
    RadiusSet or None
        The set chosen inside the block (``None`` if none was chosen, meaning
        :data:`DEFAULT_RADII`).

    Examples
    --------
    >>> from crevice.models import Atom
    >>> from crevice.radii import atom_vdw_radius, use_radii
    >>> atom = Atom(1, "OG", "SER", "A", 1, 0.0, 0.0, 0.0, "O")
    >>> with use_radii("hole"):
    ...     atom_vdw_radius(atom)
    1.65
    >>> atom_vdw_radius(atom)
    1.52
    """
    resolved = resolve_radii(radii)
    if resolved is None:
        yield _ACTIVE_RADII.get()
        return
    token = _ACTIVE_RADII.set(resolved)
    try:
        yield resolved
    finally:
        _ACTIVE_RADII.reset(token)


def radii_provenance(radii: "RadiusSet | str | os.PathLike | None" = None) -> dict:
    """Provenance of the given set, or of the set in effect.

    Parameters
    ----------
    radii : RadiusSet, str, os.PathLike or None, default None
        A set to describe; ``None`` describes :func:`effective_radii`.

    Returns
    -------
    dict
        :meth:`RadiusSet.provenance`. Since 2 Oct 2026 this is never ``None``:
        runs with the default set record :data:`DEFAULT_RADII` (name
        ``"default"``). Outputs written before then recorded nothing for the
        default set; for those, absence meant the standard table *at the time*
        (see :doc:`/methods/atomic-radii`).
    """
    return effective_radii(radii).provenance()


def radii_fields() -> dict:
    """``{"radii": provenance}`` of the set in effect, for a manifest or JSON output.

    Returns
    -------
    dict
        One key, ``"radii"``, holding :func:`radii_provenance`. Always
        present, including for default runs.
    """
    return {"radii": radii_provenance()}


def _record_radii(result, provenance: dict) -> None:
    items = result if isinstance(result, (list, tuple)) else (result,)
    for item in items:
        metadata = getattr(item, "metadata", None)
        if isinstance(metadata, dict):
            metadata.setdefault("radii", provenance)


def radii_option(func):
    """Decorator applying an analysis's keyword-only ``radii`` argument.

    ``radii=None`` inherits the set in effect (:data:`DEFAULT_RADII` unless a
    caller chose one); anything else is resolved with :func:`resolve_radii`
    and used for the whole call, including nested analyses. Results with a
    ``metadata`` dict (profiles, casts, tunnels, networks, trajectory
    analyses) get ``metadata["radii"]`` provenance of the set used, default
    set included; an existing entry is kept.

    Parameters
    ----------
    func : callable
        An analysis with a keyword-only ``radii`` parameter.

    Returns
    -------
    callable
        The wrapped analysis.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with use_radii(kwargs.get("radii")) as active:
            if "radii" in kwargs:
                kwargs["radii"] = active
            result = func(*args, **kwargs)
            _record_radii(result, (DEFAULT_RADII if active is None else active).provenance())
        return result
    return wrapper


# ---------------------------------------------------------------------------
# Rebuilding and comparing recorded sets
# ---------------------------------------------------------------------------

def radius_set_from_provenance(provenance: Mapping, *, base_dir=None) -> RadiusSet:
    """Rebuild the :class:`RadiusSet` that a ``radii`` provenance record describes.

    Used where a later step must measure with the same radii as an earlier
    one (``region-trajectory`` with a prepared region definition,
    ``hydration --cavity-results``). Candidates are tried in order and the
    first whose :attr:`RadiusSet.table_sha256` equals the recorded
    ``table_sha256`` is returned:

    1. the preset of the recorded ``name``;
    2. the recorded ``source`` file, if it exists and its SHA-256 equals the
       recorded ``source_sha256`` (a relative path is taken relative to
       ``base_dir``);
    3. a set rebuilt from the recorded tables: element radii = the standard
       table without ``standard_elements_missing``, updated with
       ``element_overrides_A`` (or no element radii for a set of HOLE
       records), plus ``atom_radii_A``, ``ion_radii_A``, ``hole_records`` and
       ``default_radius_A``.

    Parameters
    ----------
    provenance : mapping
        A ``radii`` record (:meth:`RadiusSet.provenance`).
    base_dir : str or Path, optional
        Directory against which a relative ``source`` path is resolved.

    Returns
    -------
    RadiusSet

    Raises
    ------
    ValueError
        If no candidate reproduces the recorded table hash (for example, a
        record from an older CREVICE whose standard table differed).

    Examples
    --------
    >>> from crevice.radii import radii_preset, radius_set_from_provenance
    >>> radius_set_from_provenance(radii_preset("hole").provenance()).name
    'hole'
    """
    if not isinstance(provenance, Mapping) or "table_sha256" not in provenance:
        raise ValueError("Radius provenance must be a recorded 'radii' object with a table_sha256")
    target = provenance["table_sha256"]
    candidates = []
    name = str(provenance.get("name") or "")
    if name.lower() in RADIUS_PRESETS:
        candidates.append(RADIUS_PRESETS[name.lower()])
    source = provenance.get("source")
    if source and provenance.get("source_sha256"):
        path = Path(source)
        if not path.is_absolute() and base_dir is not None:
            path = Path(base_dir) / path
        if path.is_file() and _file_sha256(path) == provenance["source_sha256"]:
            try:
                candidates.append(load_radius_file(path))
            except ValueError:
                pass
    missing = set(provenance.get("standard_elements_missing") or ())
    overrides = {k: float(v) for k, v in (provenance.get("element_overrides_A") or {}).items()}
    common = dict(atom_radii={k: float(v) for k, v in (provenance.get("atom_radii_A") or {}).items()},
                  default_radius=float(provenance.get("default_radius_A", DEFAULT_RADIUS)),
                  hole_records=tuple((a, r, float(v)) for a, r, v in (provenance.get("hole_records") or ())),
                  ion_radii={k: float(v) for k, v in (provenance.get("ion_radii_A") or {}).items()},
                  source=str(source or ""), source_sha256=provenance.get("source_sha256"),
                  citation=str(provenance.get("citation") or ""))
    elements = {k: v for k, v in VDW_RADII.items() if k not in missing}
    elements.update(overrides)
    candidates.append(RadiusSet(name or "recorded", elements, **common))
    if common["hole_records"]:
        candidates.append(RadiusSet(name or "recorded", overrides, **common))
    for candidate in candidates:
        if candidate.table_sha256 == target:
            return candidate
    raise ValueError(f"Cannot rebuild the recorded radius set {name!r} (table {target[:12]}...): no preset, "
                     "source file or recorded table reproduces its hash. Pass the same radius set explicitly "
                     "with --radii / radii=.")


def same_radius_set(first: Mapping | RadiusSet | None, second: Mapping | RadiusSet | None) -> bool:
    """Whether two radius sets (or their provenance records) have the same table.

    Two sets are the same when their :attr:`RadiusSet.table_sha256` values are
    equal: every element, atom, ion and HOLE radius and the default radius
    agree. Names, sources and citations are not compared.

    Parameters
    ----------
    first, second : RadiusSet, mapping or None
        Sets or ``radii`` provenance records. ``None`` never matches.

    Returns
    -------
    bool

    Examples
    --------
    >>> from crevice.radii import radii_preset, same_radius_set
    >>> same_radius_set(radii_preset("default"), radii_preset("default").provenance())
    True
    >>> same_radius_set(radii_preset("default"), radii_preset("bondi"))
    False
    """
    def digest(value):
        if value is None:
            return None
        if isinstance(value, RadiusSet):
            return value.table_sha256
        return value.get("table_sha256") if isinstance(value, Mapping) else None
    a, b = digest(first), digest(second)
    return a is not None and a == b


def radius_set_label(value: Mapping | RadiusSet | None) -> str:
    """Short human-readable name of a set for error messages: ``name (table 1a2b3c4d5e6f)``.

    Parameters
    ----------
    value : RadiusSet, mapping or None
        A set or a ``radii`` provenance record.

    Returns
    -------
    str
        Name and the first 12 hex digits of the table hash, or
        ``"no recorded radius set"`` for ``None``.
    """
    if value is None:
        return "no recorded radius set"
    if isinstance(value, RadiusSet):
        value = {"name": value.name, "table_sha256": value.table_sha256}
    return f"{value.get('name')!s} (table {str(value.get('table_sha256'))[:12]})"


def reconcile_recorded_radii(recorded: Mapping | None, requested: "RadiusSet | str | os.PathLike | None" = None, *,
                             allow_mismatch: bool = False, what: str = "the earlier results",
                             base_dir=None) -> tuple[RadiusSet, dict]:
    """Choose the radius set for a step that reuses recorded results.

    Rules:

    * ``recorded`` is a ``radii`` record and nothing was chosen explicitly
      (``requested`` is ``None``): the recorded set is rebuilt
      (:func:`radius_set_from_provenance`) and used.
    * ``requested`` has the same table hash as ``recorded``: it is used.
    * They differ: ``ValueError``, unless ``allow_mismatch``, in which case
      ``requested`` is used and the mismatch is recorded.
    * Nothing is recorded (results written before 2 Oct 2026 with the default
      set, which recorded no ``radii``): the requested set, or the default
      set, is used and the check says it could not be verified.

    Parameters
    ----------
    recorded : mapping or None
        The ``radii`` record of the earlier results.
    requested : RadiusSet, str, os.PathLike or None, optional
        The explicitly chosen set (``--radii``/``radii=``), or ``None``.
    allow_mismatch : bool, default False
        Use ``requested`` even if it differs.
    what : str, default "the earlier results"
        Wording for the error message.
    base_dir : str or Path, optional
        Directory for a relative radius-file path in ``recorded``.

    Returns
    -------
    radii : RadiusSet
        The set to use.
    check : dict
        ``status`` (``"recorded_set_used"``, ``"requested_matches_recorded"``,
        ``"mismatch_overridden"`` or ``"unverified_no_recorded_set"``),
        ``recorded`` and ``used`` (name and table hash).

    Raises
    ------
    ValueError
        On a mismatch without ``allow_mismatch``.
    """
    chosen = resolve_radii(requested)

    def short(value):
        if value is None:
            return None
        if isinstance(value, RadiusSet):
            return {"name": value.name, "table_sha256": value.table_sha256}
        return {"name": value.get("name"), "table_sha256": value.get("table_sha256")}

    if not recorded:
        used, status = effective_radii(chosen), "unverified_no_recorded_set"
    elif chosen is None:
        used, status = radius_set_from_provenance(recorded, base_dir=base_dir), "recorded_set_used"
    elif same_radius_set(chosen, recorded):
        used, status = chosen, "requested_matches_recorded"
    elif allow_mismatch:
        used, status = chosen, "mismatch_overridden"
    else:
        raise ValueError(f"Radius set mismatch: {what} used {radius_set_label(recorded)} but "
                         f"{radius_set_label(chosen)} was requested. Omit --radii to use the recorded set, or pass "
                         "--allow-radii-mismatch to proceed with the requested set (recorded).")
    return used, {"status": status, "recorded": short(recorded), "used": short(used)}

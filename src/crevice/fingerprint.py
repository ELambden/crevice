"""Pore fingerprint: radius, lining hydrophobicity, chemistry, water and conductance along a channel.

One table along a single path coordinate ``s`` (Å) of a resolved channel
combines what CREVICE measures separately: the radius profile, the
hydrophobicity and chemistry of the residues lining the path, the ENTRY /
LUMEN / EXIT segments (:mod:`crevice.cast_segments`), and, for MD
trajectories, the radius distribution and the water-oxygen count along the
path. A geometric conductance estimate summarises the radius profile.

Definitions (all distances in Å; atoms are van der Waals spheres of the
radius set in force, :func:`crevice.radii.effective_radii`):

* **Path coordinate s.** LUMEN samples are the channel profile samples, at
  ``s = t - t_first`` (the profile's axial coordinate ``t`` for axial
  profiles, the cumulative path length for explicit and tunnel paths), so the
  lumen spans ``0 <= s <= L``. Beyond each mouth the path continues along one
  centre line: at an open (axial) end the straight vestibule centre line of
  :func:`crevice.cast_segments.axial_vestibule_leg` (needs the channel cast;
  every 0.25 Å out to bulk), at a capped end the widest lateral exit leg
  (rank 1) of :func:`crevice.channel_exits.lateral_exit_legs`. Their samples
  sit at ``s = -arc`` (lower end) or ``s = L + arc`` (upper end), ``arc``
  being the path length along that centre line. Other lateral legs are
  reported (rank, length, bottleneck) but not drawn on ``s``.
* **Segments.** ENTRY/EXIT labels follow
  :func:`crevice.cast_segments.segment_roles` (geometric labels only, no
  transport direction). The vestibule or leg samples of an end take that
  end's label (``entry_1``, ``exit_1``); profile samples are ``lumen``.
* **Radius** ``radius_A``: the profile radius (free-sphere clearance minus the
  profile probe radius) for LUMEN samples; the centre-line clearance minus the
  same probe radius for vestibule and leg samples (a free-sphere radius along
  a 3D path, not a planar section radius).
* **Lining atoms and residues.** An atom lines a sample when its van der Waals
  surface lies within ``lining_shell`` (default 1.4 Å, one water-probe radius)
  of the sample's free sphere: ``|x - c| - r_atom - R <= lining_shell``, with
  ``c`` the sample centre and ``R`` its raw clearance. Hydrogens, water and
  recognised monatomic ions are excluded. A residue lines a sample when any of
  its atoms does; its **anchor** ``s_j`` is the sample at which its surface
  gap is smallest.
* **Hydrophobicity** ``hydrophobicity_kcal_mol``: Nadaraya-Watson Gaussian
  kernel average along ``s`` of the residue values ``h_j`` of all lining
  residues, ``H(s) = sum_j h_j K(s - s_j) / sum_j K(s - s_j)`` with
  ``K(d) = exp(-d**2 / (2 sigma**2))``, one weight per residue, and
  ``sigma = hydrophobicity_bandwidth`` (default 3.0 Å, about one helical turn
  of rise). Where ``sum_j K < 0.01`` (no residue within about 3 sigma) the
  value is empty. The same kernel form is used by CHAP (Klesse et al. 2019).
  ``h_j`` is the **Wimley-White whole-residue octanol scale** (Wimley,
  Creamer & White 1996, *Biochemistry* 35:5109; White & Wimley 1999, *Annu.
  Rev. Biophys. Biomol. Struct.* 28:319) with its sign reversed,
  ``h = -dG(water -> octanol)``, so positive values are hydrophobic. The
  octanol scale was chosen because it is experimentally determined for whole
  residues (side chain plus peptide unit), includes charged and neutral forms,
  and describes partitioning into a low-dielectric environment, the relevant
  contrast for a pore wall. Asp, Glu, Lys and Arg are taken as charged and His
  as neutral (pH 7); protonation-variant names (HSD, HSE, HIE, HID, HSP, HIP,
  ASH, GLH, LYN, CYX ...) and common D-amino-acid and modified residues
  (DLE, DVA, MSE ...) map to their parent and charge state
  (:data:`RESIDUE_PARENTS`); other residues (ligands, lipids, ethanolamine
  caps) have no value and are reported as ``other``. Residue-level values
  cannot see that a residue lines the pore with its backbone only (for
  example gramicidin A, whose lumen is lined by carbonyl oxygens of
  hydrophobic residues).
* **Chemistry** classes come from :func:`crevice.radii.residue_properties`
  of the parent residue: ``hydrophobic``, ``polar``, ``charged`` (positive or
  negative), ``aromatic``; residues may belong to several (Tyr is
  hydrophobic, polar and aromatic), so counts per class do not add up to the
  residue count. Residues with no class are ``other``.
* **Water (trajectories)** along the fixed reference path: in every aligned
  frame (same rigid fit as the protein, periodic minimum image about the path
  centroid, :func:`crevice.water_membership.image_waters_around_anchor`) a
  water oxygen belongs to sample ``i`` when ``i`` is the nearest path sample
  and the oxygen centre lies inside that frame's free sphere at the sample
  (distance to the centre at most the frame's raw clearance there). The bin
  width ``ds_i`` is half the distance to the neighbouring samples. Reported:
  mean count per frame, linear density (mean count / ``ds_i``, molecules per
  Å), number density (mean count over the mean of ``pi R**2 ds_i``, Å⁻³,
  bulk water is about 0.0334 Å⁻³) and occupancy (fraction of frames with at
  least one oxygen in the bin). The free sphere underestimates a wide,
  non-circular lumen, so densities are lower bounds of what fills the cross
  section.
* **Conductance** (static and per frame; LUMEN only): the HOLE geometric
  estimate (Smart et al. 1993, *Biophys. J.* 65:2455; Smart, Breed, Smith &
  Sansom 1997, *Biophys. J.* 72:1109), ``G_macro = 1 / (rho * F)`` with
  ``F = integral ds / (pi r(s)**2)`` over the lumen (trapezoid rule), and
  ``rho`` the bulk resistivity of the electrolyte (default 150 mM KCl,
  0.53 Ohm m at 25 °C, :data:`ELECTROLYTES`). Also reported: ``G_access``,
  which adds the access resistance ``rho / (4 a)`` of each mouth of radius
  ``a`` (Hall 1975, *J. Gen. Physiol.* 66:531; Hille, *Ion Channels of
  Excitable Membranes*, 3rd ed., 2001, ch. 11); and ``G_hole_corrected =
  G_macro / 5.594``, HOLE 2.3.1's second-generation constant empirical
  divisor (source ``hcapgr.f``, calibrated on eight channels in 1 M KCl with
  the capsule option, Smart et al. 1997), given for orientation only. All are
  **order-of-magnitude geometric estimates, not measurements**: they ignore
  ion-wall interactions, selectivity, dehydration and surface charge, and
  macroscopic continuum resistivity is not valid in a single-file pore. A
  zero radius anywhere in the lumen gives ``G = 0`` (occluded).
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from .presentation import FigureText, annotate_option

#: Wimley-White whole-residue octanol scale, dG(water -> octanol) in kcal/mol
#: (Wimley, Creamer & White 1996; White & Wimley 1999). Keys: parent name plus
#: charge-state suffix where the scale has two forms.
WIMLEY_WHITE_OCTANOL = {
    "ALA": 0.50, "ARG+": 1.81, "ASN": 0.85, "ASP-": 3.64, "ASP0": 0.43, "CYS": -0.02, "GLN": 0.77,
    "GLU-": 3.63, "GLU0": 0.11, "GLY": 1.15, "HIS+": 2.33, "HIS0": 0.11, "ILE": -1.12, "LEU": -1.25,
    "LYS+": 2.80, "MET": -0.67, "PHE": -1.71, "PRO": 0.14, "SER": 0.46, "THR": 0.25, "TRP": -2.09,
    "TYR": -0.71, "VAL": -0.46,
}
#: Default charge state of each parent residue at pH 7 (key into the scale).
DEFAULT_STATE = {"ARG": "ARG+", "ASP": "ASP-", "GLU": "GLU-", "HIS": "HIS0", "LYS": "LYS+"}
#: Residue-name variants mapped to (parent residue, scale key).
RESIDUE_PARENTS = {
    "HSD": ("HIS", "HIS0"), "HSE": ("HIS", "HIS0"), "HIE": ("HIS", "HIS0"), "HID": ("HIS", "HIS0"),
    "HSP": ("HIS", "HIS+"), "HIP": ("HIS", "HIS+"), "HIS1": ("HIS", "HIS0"), "HISA": ("HIS", "HIS0"),
    "HISB": ("HIS", "HIS0"), "HISH": ("HIS", "HIS+"), "ASH": ("ASP", "ASP0"), "ASPP": ("ASP", "ASP0"),
    "ASPH": ("ASP", "ASP0"), "GLH": ("GLU", "GLU0"), "GLUP": ("GLU", "GLU0"), "GLUH": ("GLU", "GLU0"),
    "LYN": ("LYS", "LYS+"), "LSN": ("LYS", "LYS+"), "CYX": ("CYS", "CYS"), "CYM": ("CYS", "CYS"),
    "MSE": ("MET", "MET"), "DAL": ("ALA", "ALA"), "DLE": ("LEU", "LEU"), "DVA": ("VAL", "VAL"),
    "DIL": ("ILE", "ILE"), "DPN": ("PHE", "PHE"), "DTR": ("TRP", "TRP"), "DTY": ("TYR", "TYR"),
    "DSN": ("SER", "SER"), "DTH": ("THR", "THR"), "DGN": ("GLN", "GLN"), "DSG": ("ASN", "ASN"),
    "DAS": ("ASP", "ASP-"), "DGL": ("GLU", "GLU-"), "DLY": ("LYS", "LYS+"), "DAR": ("ARG", "ARG+"),
    "DHI": ("HIS", "HIS0"), "DPR": ("PRO", "PRO"), "DCY": ("CYS", "CYS"), "MED": ("MET", "MET"),
    "FVA": ("VAL", "VAL"),
}
HYDROPHOBICITY_SCALE = {
    "name": "wimley_white_octanol",
    "values_kcal_mol": "-dG(water->octanol) of the whole residue; positive = hydrophobic",
    "citation": "Wimley WC, Creamer TP, White SH (1996) Biochemistry 35:5109-5124; "
                "White SH, Wimley WC (1999) Annu Rev Biophys Biomol Struct 28:319-365",
    "charge_states": "Asp, Glu, Lys, Arg charged; His neutral (pH 7); protonation-variant names map to their state",
}
#: Bulk resistivities (Ohm m) of common electrolytes at 25 °C, from the
#: electrical conductivities tabulated in the CRC Handbook of Chemistry and
#: Physics ("Electrical Conductivity of Aqueous Solutions"; 150 mM values
#: interpolated): 150 mM KCl 1.88 S/m, 150 mM NaCl 1.56 S/m, 1 M KCl 11.2 S/m,
#: 1 M NaCl 8.6 S/m. HOLE uses 1/12 Ohm m for 1 M KCl.
ELECTROLYTES = {"kcl-150mm": 0.532, "nacl-150mm": 0.641, "kcl-1m": 0.0894, "nacl-1m": 0.116}
#: HOLE 2.3.1 second-generation constant empirical divisor (hcapgr.f; Smart et al. 1997).
HOLE_EMPIRICAL_DIVISOR = 5.594
#: Bulk liquid-water number density, molecules per Å³ (about 0.997 g/cm³ at 25 °C).
BULK_WATER_DENSITY_A3 = 0.0334
#: Lining chemistry classes, in display order.
CHEMISTRY_CLASSES = ("hydrophobic", "polar", "charged", "aromatic")
CONDUCTANCE_DEFINITION = (
    "G_macro = 1/(rho*F), F = integral over the lumen of ds/(pi r(s)^2) (trapezoid); G_access adds rho/(4a) at "
    "each mouth (Hall 1975; Hille 2001); G_hole_corrected = G_macro/5.594 (HOLE 2.3.1 constant divisor, 1 M KCl "
    "calibration, Smart et al. 1997). Order-of-magnitude geometric estimates, not measurements")
CONDUCTANCE_CITATIONS = (
    "Smart OS, Goodfellow JM, Wallace BA (1993) Biophys J 65:2455-2460",
    "Smart OS, Breed J, Smith GR, Sansom MSP (1997) Biophys J 72:1109-1126",
    "Hall JE (1975) J Gen Physiol 66:531-532",
    "Hille B (2001) Ion Channels of Excitable Membranes, 3rd ed., Sinauer, ch. 11")
#: Columns of the per-sample fingerprint CSV.
SAMPLE_FIELDS = (
    "s_A", "segment", "segment_type", "source", "t_A", "x_A", "y_A", "z_A", "radius_A", "raw_clearance_A",
    "hydrophobicity_kcal_mol", "hydrophobicity_support", "lining_residue_count", "n_hydrophobic", "n_polar",
    "n_charged", "n_aromatic", "nearest_lining_residue")
#: Extra per-sample columns of a trajectory fingerprint.
ENSEMBLE_FIELDS = (
    "radius_mean_A", "radius_sd_A", "radius_q_low_A", "radius_q_high_A", "water_count_mean",
    "water_linear_density_per_A", "water_number_density_per_A3", "water_occupancy")


def residue_hydrophobicity(resname: str) -> tuple[str | None, float | None]:
    """Parent residue and Wimley-White hydrophobicity of a residue name.

    Parameters
    ----------
    resname : str
        Residue name (case-insensitive), including protonation variants and the
        D-amino-acid and modified names of :data:`RESIDUE_PARENTS`.

    Returns
    -------
    parent : str or None
        Standard parent residue name, or ``None`` when not an amino acid of the
        scale.
    value : float or None
        ``-dG(water -> octanol)`` in kcal/mol (positive = hydrophobic), or
        ``None``.

    Examples
    --------
    >>> residue_hydrophobicity("LEU"), residue_hydrophobicity("asp")
    (('LEU', 1.25), ('ASP', -3.64))
    >>> residue_hydrophobicity("HSD"), residue_hydrophobicity("POPC")
    (('HIS', -0.11), (None, None))
    """
    name = resname.strip().upper()
    if name in RESIDUE_PARENTS:
        parent, key = RESIDUE_PARENTS[name]
    else:
        parent, key = name, DEFAULT_STATE.get(name, name)
    if key not in WIMLEY_WHITE_OCTANOL:
        return None, None
    return parent, -WIMLEY_WHITE_OCTANOL[key]


def residue_classes(resname: str) -> tuple[str, ...]:
    """Chemistry classes (hydrophobic, polar, charged, aromatic) of a residue.

    Parameters
    ----------
    resname : str
        Residue name; variants map to their parent first
        (:func:`residue_hydrophobicity`).

    Returns
    -------
    tuple of str
        Subset of :data:`CHEMISTRY_CLASSES` from
        :func:`crevice.radii.residue_properties` (``positive`` and
        ``negative`` become ``charged``), or ``("other",)``.

    Examples
    --------
    >>> residue_classes("TYR"), residue_classes("HSP"), residue_classes("GLY")
    (('hydrophobic', 'polar', 'aromatic'), ('polar', 'charged', 'aromatic'), ('other',))
    """
    from .radii import residue_properties
    parent = residue_hydrophobicity(resname)[0] or resname
    props = set(residue_properties(parent))
    if props & {"positive", "negative"}:
        props.add("charged")
    out = tuple(c for c in CHEMISTRY_CLASSES if c in props)
    return out or ("other",)


def resistance_integral(s, radius) -> float:
    """Geometric factor ``F = integral ds / (pi r(s)**2)`` by the trapezoid rule.

    Parameters
    ----------
    s : sequence of float
        Path coordinate, Å, increasing.
    radius : sequence of float
        Radius at each ``s``, Å.

    Returns
    -------
    float
        ``F`` in Å⁻¹; ``inf`` if any radius is zero or negative; ``0.0`` for
        fewer than two samples.

    Examples
    --------
    >>> round(resistance_integral([0.0, 10.0, 20.0], [2.0, 2.0, 2.0]), 6)   # 20 / (pi * 4)
    1.591549
    """
    if any(r <= 0 for r in radius):
        return math.inf
    inv = [1.0/(math.pi*r*r) for r in radius]
    return sum(0.5*(inv[i]+inv[i+1])*(s[i+1]-s[i]) for i in range(len(s)-1))


def conductance_estimate(s, radius, *, resistivity_ohm_m: float = ELECTROLYTES["kcl-150mm"],
                         electrolyte: str | None = "kcl-150mm") -> dict:
    """HOLE-style geometric conductance estimate of a pore profile.

    Parameters
    ----------
    s : sequence of float
        Path coordinate of the lumen samples, Å, increasing.
    radius : sequence of float
        Radius of each sample, Å.
    resistivity_ohm_m : float, default 0.532
        Bulk resistivity of the electrolyte, Ohm m (150 mM KCl at 25 °C).
    electrolyte : str or None, default "kcl-150mm"
        Name recorded with the value (``None`` for a custom resistivity).

    Returns
    -------
    dict
        ``geometric_factor_per_A`` (``F``), ``length_A``, ``min_radius_A``,
        ``g_macro_pS`` (``1/(rho F)``), ``g_access_pS`` (with the Hall access
        resistance ``rho/(4a)`` of both end radii), ``g_hole_corrected_pS``
        (``g_macro / 5.594``), ``occluded``, ``resistivity_ohm_m``,
        ``electrolyte``, ``definition``, ``citations`` and ``caveat``.

    Raises
    ------
    ValueError
        If the resistivity is not finite and positive, or the arrays differ in
        length or have fewer than two samples.

    Examples
    --------
    >>> g = conductance_estimate([0.0, 25.0], [1.5, 1.5], resistivity_ohm_m=1/12)
    >>> round(g["g_macro_pS"], 1)   # 12 S/m * pi (1.5 Å)^2 / 25 Å
    339.3
    """
    if not math.isfinite(resistivity_ohm_m) or resistivity_ohm_m <= 0:
        raise ValueError("resistivity must be finite and positive")
    if len(s) != len(radius) or len(s) < 2:
        raise ValueError("conductance needs at least two samples with one radius each")
    factor = resistance_integral(s, radius)
    occluded = not math.isfinite(factor)
    rho = resistivity_ohm_m
    # R [Ohm] = rho [Ohm m] * F [1/Å] * 1e10; G [pS] = 1e12 / R
    g_macro = 0.0 if occluded else 100.0/(rho*factor)
    if occluded:
        g_access = 0.0
    else:
        access = 1.0/(4*radius[0]) + 1.0/(4*radius[-1])
        g_access = 100.0/(rho*(factor+access))
    return {"geometric_factor_per_A": factor if not occluded else None, "length_A": s[-1]-s[0],
            "min_radius_A": min(radius), "g_macro_pS": g_macro, "g_access_pS": g_access,
            "g_hole_corrected_pS": g_macro/HOLE_EMPIRICAL_DIVISOR, "occluded": occluded,
            "resistivity_ohm_m": rho, "electrolyte": electrolyte, "definition": CONDUCTANCE_DEFINITION,
            "citations": list(CONDUCTANCE_CITATIONS),
            "caveat": "order-of-magnitude geometric estimate, not a measurement; continuum resistivity ignores "
                      "ion-wall interactions, selectivity, dehydration and surface charge"}


def resolve_resistivity(electrolyte: str | None = "kcl-150mm", resistivity: float | None = None) -> tuple[float, str | None]:
    """Resistivity (Ohm m) and its label from a preset name or an explicit value.

    Parameters
    ----------
    electrolyte : str or None, default "kcl-150mm"
        Key of :data:`ELECTROLYTES`.
    resistivity : float, optional
        Explicit resistivity, Ohm m; overrides ``electrolyte``.

    Returns
    -------
    tuple of (float, str or None)
        Resistivity and the label (``"custom"`` for an explicit value).

    Raises
    ------
    ValueError
        For an unknown electrolyte name.
    """
    if resistivity is not None:
        return float(resistivity), "custom"
    if electrolyte not in ELECTROLYTES:
        raise ValueError(f"unknown electrolyte {electrolyte!r}; choose from {', '.join(ELECTROLYTES)}")
    return ELECTROLYTES[electrolyte], electrolyte


@dataclass
class PoreFingerprint:
    """Fingerprint of one channel along its path coordinate ``s``.

    Attributes
    ----------
    samples : list of dict
        One row per path sample in increasing ``s``; keys
        :data:`SAMPLE_FIELDS`, plus :data:`ENSEMBLE_FIELDS` for trajectories.
    lining : list of dict
        One row per lining residue: ``residue``, ``resname``, ``parent``,
        ``hydrophobicity_kcal_mol``, ``classes``, ``anchor_s_A``,
        ``anchor_segment``, ``min_gap_A``, ``samples_lined``, ``segments``.
    segments : list of dict
        One row per segment in ENTRY, LUMEN, EXIT order: ``segment``,
        ``segment_type``, ``source``, ``s_start_A``, ``s_end_A``,
        ``length_A``, ``min_radius_A``, ``residue_count``, per-class counts
        ``n_<class>``, ``n_other``, ``mean_hydrophobicity_kcal_mol`` and
        ``residues``.
    exit_legs : list of dict
        Every lateral exit leg: ``end``, ``segment``, ``leg``, ``length_A``,
        ``bottleneck_radius_A``, ``on_path`` (rank 1 legs are on ``s``).
    conductance : dict
        :func:`conductance_estimate` of the lumen; for trajectories also
        ``per_frame_g_macro_pS`` and summary statistics.
    parameters : dict
        Definitions, scale, bandwidth, lining shell and source settings.
    water : dict or None
        Water selection, frame count and definitions (trajectories only).
    frames : list of dict
        Per-frame rows (trajectories only): ``frame_index``, ``time_ps``,
        ``fit_rmsd_A``, ``lumen_min_radius_A``, ``g_macro_pS``,
        ``g_access_pS`` and ``pore_water_count``.
    """

    samples: list
    lining: list
    segments: list
    exit_legs: list
    conductance: dict
    parameters: dict
    water: dict | None = None
    frames: list = field(default_factory=list)

    def to_dict(self) -> dict:
        """JSON summary (everything except the per-sample rows and per-frame rows).

        Returns
        -------
        dict
            ``conductance``, ``segments``, ``lining_residues``, ``exit_legs``,
            ``parameters``, ``water``, ``frames`` and ``sample_count``.
        """
        return {"conductance": self.conductance, "segments": self.segments, "lining_residues": self.lining,
                "exit_legs": self.exit_legs, "parameters": self.parameters, "water": self.water,
                "frames": self.frames, "sample_count": len(self.samples)}


def _excluded(atom) -> bool:
    from .hydration import WATER_RESNAMES
    from .radii import is_recognised_ion
    return (atom.element.upper() in {"H", "D"} or atom.resname.strip().upper() in WATER_RESNAMES
            or is_recognised_ion(atom.name, atom.resname, atom.element))


def _path_track(frame, profile, cast, entry_end, bulk_radius):
    """Path samples (dicts) along s and the lateral exit-leg report."""
    from .cast_segments import axial_vestibule_leg, segment_roles
    roles, rule = segment_roles(profile, entry_end)
    t0 = profile.points[0].t
    length = profile.points[-1].t - t0
    rows = [{"s_A": p.t - t0, "segment": "lumen", "segment_type": "lumen", "source": "channel_profile",
             "t_A": p.t, "position": tuple(p.position), "radius_A": p.radius, "raw_clearance_A": p.raw_clearance}
            for p in profile.points]
    exits = profile.metadata.get("exits") or {}
    legs_report = []
    probe = profile.probe_radius
    for end in ("lower", "upper"):
        kind = roles[end]
        record = exits.get(end) or {}
        sign = -1.0 if end == "lower" else 1.0
        base = 0.0 if end == "lower" else length
        centre = None
        if record.get("type") == "lateral":
            legs = sorted(record.get("legs", ()), key=lambda leg: leg.get("leg", 1))
            for leg in legs:
                rank = leg.get("leg", 1)
                legs_report.append({"end": end, "segment": f"{kind}_{rank}", "leg": rank,
                                    "length_A": leg.get("length_A"),
                                    "bottleneck_radius_A": leg.get("bottleneck_radius_A"), "on_path": rank == legs[0].get("leg", 1)})
            if legs:
                leg = legs[0]
                points = [tuple(p) for p in leg["points"]]
                raw = list(leg["raw_clearances_A"])
                arc = list(leg.get("arc_lengths_A") or [])
                if len(arc) != len(points):
                    arc = [0.0]
                    for a, b in zip(points, points[1:]):
                        arc.append(arc[-1]+math.dist(a, b))
                centre = (points, raw, arc, "lateral_exit_leg", f"{kind}_{legs[0].get('leg', 1)}")
        elif cast is not None and profile.method == "axial-connected" and cast.metadata.get("grid_basis") is not None:
            leg = axial_vestibule_leg(frame, profile, cast, end, bulk_radius=bulk_radius)
            if leg is not None:
                points = [tuple(p) for p in leg["points"]]
                arc = [0.25*(k+1) for k in range(len(points))]
                centre = (points, list(leg["raw_clearances_A"]), arc, "axial_vestibule", f"{kind}_1")
        if centre is None:
            continue
        points, raw, arc, source, name = centre
        for p, c, a in zip(points, raw, arc):
            if a <= 1e-9:
                continue  # the mouth sample itself is already the lumen end
            rows.append({"s_A": base+sign*a, "segment": name, "segment_type": kind, "source": source, "t_A": None,
                         "position": p, "radius_A": max(0.0, c-probe), "raw_clearance_A": c})
    rows.sort(key=lambda r: r["s_A"])
    return rows, legs_report, rule


def _lining(frame, rows, shell, include_hetero=True):
    """Lining atoms per sample: (per-sample residue sets, residue records)."""
    import numpy as np
    from scipy.spatial import cKDTree
    from .radii import effective_radii
    radii = effective_radii()
    atoms = [a for a in frame.atoms if not _excluded(a) and (include_hetero or not a.hetero)]
    if not atoms:
        return [set() for _ in rows], {}
    xyz = np.asarray([a.coord for a in atoms], dtype=float)
    rad = np.asarray([radii.radius_for_atom(a) for a in atoms], dtype=float)
    tree = cKDTree(xyz)
    rmax = float(rad.max())
    per_sample, records = [], {}
    for i, row in enumerate(rows):
        centre = np.asarray(row["position"], dtype=float)
        reach = max(row["raw_clearance_A"], 0.0) + shell + rmax
        found = tree.query_ball_point(centre, reach)
        keys = set()
        for k in found:
            gap = float(np.linalg.norm(xyz[k]-centre) - rad[k] - max(row["raw_clearance_A"], 0.0))
            if gap > shell:
                continue
            key = atoms[k].residue_key
            keys.add(key)
            rec = records.setdefault(key, {"gap": math.inf, "index": i, "samples": set()})
            rec["samples"].add(i)
            if gap < rec["gap"]:
                rec["gap"], rec["index"] = gap, i
        per_sample.append(keys)
    return per_sample, records


def pore_fingerprint(frame, profile, *, cast=None, entry_end: str | None = "auto", bulk_radius: float = 6.0,
                     hydrophobicity_bandwidth: float = 3.0, lining_shell: float = 1.4,
                     electrolyte: str | None = "kcl-150mm", resistivity: float | None = None,
                     include_hetero: bool = True) -> PoreFingerprint:
    """Fingerprint of a resolved channel in one structure.

    Definitions: module docstring.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile.
    profile : PoreProfile
        Resolved channel profile (at least two samples).
    cast : VoidCast, optional
        Connected channel cast of the profile; needed for the straight
        vestibule centre lines beyond open (axial) mouths. Without it only
        lateral exit legs extend the path beyond the lumen.
    entry_end : {"auto", "start", "end"} or None, default "auto"
        ENTRY/EXIT labelling rule (:func:`crevice.cast_segments.segment_roles`).
    bulk_radius : float, default 6.0
        Bulk-probe radius (Å) that ends the vestibule centre lines.
    hydrophobicity_bandwidth : float, default 3.0
        Gaussian kernel standard deviation ``sigma`` along ``s``, Å.
    lining_shell : float, default 1.4
        Surface-gap shell around each sample's free sphere, Å.
    electrolyte : str or None, default "kcl-150mm"
        Resistivity preset (:data:`ELECTROLYTES`).
    resistivity : float, optional
        Explicit bulk resistivity, Ohm m; overrides ``electrolyte``.
    include_hetero : bool, default True
        Count HETATM residues as lining (use the profile's setting; CLI
        ``--exclude-hetero`` turns it off).

    Returns
    -------
    PoreFingerprint

    Raises
    ------
    ValueError
        For fewer than two profile samples, a non-positive bandwidth or a
        negative shell.
    """
    if len(profile.points) < 2:
        raise ValueError("a fingerprint needs a resolved profile with at least two samples")
    if not math.isfinite(hydrophobicity_bandwidth) or hydrophobicity_bandwidth <= 0:
        raise ValueError("hydrophobicity bandwidth must be finite and positive")
    if not math.isfinite(lining_shell) or lining_shell < 0:
        raise ValueError("lining shell must be finite and non-negative")
    rho, label = resolve_resistivity(electrolyte, resistivity)
    rows, legs, rule = _path_track(frame, profile, cast, entry_end, bulk_radius)
    per_sample, records = _lining(frame, rows, lining_shell, include_hetero)
    s = [r["s_A"] for r in rows]
    lining = []
    by_key = {}
    for key in sorted(records, key=lambda k: (records[k]["index"], k)):
        rec = records[key]
        parent, value = residue_hydrophobicity(key.resname)
        segs = sorted({rows[i]["segment"] for i in rec["samples"]}, key=lambda n: (n.split("_")[0] != "entry", n != "lumen", n))
        entry = {"residue": key.label, "resname": key.resname, "parent": parent, "hydrophobicity_kcal_mol": value,
                 "classes": list(residue_classes(key.resname)), "anchor_s_A": s[rec["index"]],
                 "anchor_segment": rows[rec["index"]]["segment"], "min_gap_A": round(rec["gap"], 4),
                 "samples_lined": len(rec["samples"]), "segments": segs}
        lining.append(entry)
        by_key[key] = entry
    sigma = hydrophobicity_bandwidth
    valued = [(e["anchor_s_A"], e["hydrophobicity_kcal_mol"]) for e in lining if e["hydrophobicity_kcal_mol"] is not None]
    for row, keys in zip(rows, per_sample):
        num = den = 0.0
        for sj, h in valued:
            w = math.exp(-0.5*((row["s_A"]-sj)/sigma)**2)
            num += w*h
            den += w
        row["hydrophobicity_kcal_mol"] = num/den if den >= 0.01 else None
        row["hydrophobicity_support"] = den
        row["lining_residue_count"] = len(keys)
        for c in CHEMISTRY_CLASSES:
            row[f"n_{c}"] = sum(c in by_key[k]["classes"] for k in keys)
        nearest = min(keys, key=lambda k: records[k]["gap"], default=None)
        row["nearest_lining_residue"] = nearest.label if nearest else ""
    segments = []
    names = list(dict.fromkeys(r["segment"] for r in rows))
    names.sort(key=lambda n: (n.split("_")[0] != "entry", n != "lumen", n))
    for name in names:
        idx = [i for i, r in enumerate(rows) if r["segment"] == name]
        keys = set().union(*(per_sample[i] for i in idx))
        entries = [by_key[k] for k in sorted(keys, key=lambda k: (records[k]["index"], k))]
        values = [e["hydrophobicity_kcal_mol"] for e in entries if e["hydrophobicity_kcal_mol"] is not None]
        segment = {"segment": name, "segment_type": rows[idx[0]]["segment_type"], "source": rows[idx[0]]["source"],
                   "s_start_A": s[idx[0]], "s_end_A": s[idx[-1]], "length_A": s[idx[-1]]-s[idx[0]],
                   "min_radius_A": min(rows[i]["radius_A"] for i in idx), "residue_count": len(entries)}
        for c in CHEMISTRY_CLASSES:
            segment[f"n_{c}"] = sum(c in e["classes"] for e in entries)
        segment["n_other"] = sum(e["classes"] == ["other"] for e in entries)
        segment["mean_hydrophobicity_kcal_mol"] = sum(values)/len(values) if values else None
        segment["residues"] = [e["residue"] for e in entries]
        segments.append(segment)
    lumen = [r for r in rows if r["segment"] == "lumen"]
    conductance = conductance_estimate([r["s_A"] for r in lumen], [r["radius_A"] for r in lumen],
                                       resistivity_ohm_m=rho, electrolyte=label)
    for row in rows:
        x, y, z = row.pop("position")
        row.update({"x_A": x, "y_A": y, "z_A": z})
    parameters = {
        "path_coordinate": "s = t - t_first on the channel profile (lumen 0..L); vestibule/leg samples at -arc "
                           "(lower end) or L + arc (upper end) along their centre line",
        "profile_method": profile.method, "probe_radius_A": profile.probe_radius, "label_rule": rule,
        "entry_end": entry_end or "auto", "bulk_radius_A": bulk_radius, "vestibules": cast is not None,
        "lining_shell_A": lining_shell, "include_hetero": include_hetero,
        "lining_definition": "atom (no H, water or ions) whose vdW surface is within lining_shell_A of the sample's "
                             "free sphere; residue anchored at the sample of smallest surface gap",
        "hydrophobicity_scale": HYDROPHOBICITY_SCALE, "hydrophobicity_bandwidth_A": sigma,
        "hydrophobicity_definition": "Nadaraya-Watson Gaussian kernel average over lining-residue anchors along s, "
                                     "one weight per residue; empty where the kernel sum < 0.01",
        "chemistry_classes": "crevice.radii.residue_properties of the parent residue; positive/negative -> charged; "
                             "a residue may be in several classes",
        "radii": _radii_record(),
    }
    return PoreFingerprint(rows, lining, segments, legs, conductance, parameters)


def _radii_record():
    from .radii import radii_provenance
    try:
        return radii_provenance()
    except Exception:  # pragma: no cover - provenance is informative only
        return None


def fingerprint_structure(frame, *, profile=None, cast_spacing: float = 0.5, cast_max_grid_points: int = 16_000_000,
                          entry_end: str | None = "auto", bulk_radius: float = 6.0, hydrophobicity_bandwidth: float = 3.0,
                          lining_shell: float = 1.4, electrolyte: str | None = "kcl-150mm",
                          resistivity: float | None = None, **profile_kwargs) -> tuple:
    """Resolve the channel, its cast and its fingerprint in one structure.

    Parameters
    ----------
    frame : StructureFrame
        Structure to analyse.
    profile : PoreProfile, optional
        Already resolved profile; otherwise :func:`crevice.channels.pore_profile`
        runs with ``profile_kwargs``.
    cast_spacing : float, default 0.5
        Channel-cast grid spacing, Å (as ``publish``).
    cast_max_grid_points : int, default 16000000
        Grid-size guard of the cast.
    entry_end, bulk_radius, hydrophobicity_bandwidth, lining_shell, electrolyte, resistivity
        As :func:`pore_fingerprint`.
    **profile_kwargs
        Passed to :func:`crevice.channels.pore_profile` (``axis``, ``origin``,
        ``lateral_exits`` ...).

    Returns
    -------
    profile : PoreProfile
    cast : VoidCast or None
        The channel cast, or ``None`` when it is not a connected channel cast.
    fingerprint : PoreFingerprint

    Raises
    ------
    ValueError
        If the profile is unresolved.
    """
    from .analysis import pore_profile
    from .voids import void_cast
    if profile is None:
        profile = pore_profile(frame, **profile_kwargs)
    if profile.metadata.get("status", "resolved") != "resolved" or len(profile.points) < 2:
        raise ValueError("no resolved channel profile; a fingerprint needs one (see crevice profile)")
    cast = None
    if profile.method == "axial-connected":
        cast = void_cast(frame, mode="channel", spacing=cast_spacing, min_radius=0.0, max_components=1,
                         max_grid_points=cast_max_grid_points, profile=profile,
                         include_hydrogen=profile_kwargs.get("include_hydrogen", False),
                         include_hetero=profile_kwargs.get("include_hetero", True))
        if cast.metadata.get("export_mode") != "connected_section_fill":
            cast = None
    fp = pore_fingerprint(frame, profile, cast=cast, entry_end=entry_end, bulk_radius=bulk_radius,
                          hydrophobicity_bandwidth=hydrophobicity_bandwidth, lining_shell=lining_shell,
                          electrolyte=electrolyte, resistivity=resistivity,
                          include_hetero=profile_kwargs.get("include_hetero", True))
    return profile, cast, fp


def _bin_widths(s):
    n = len(s)
    if n == 1:
        return [1.0]
    out = []
    for i in range(n):
        lo = s[i] - (s[i]-s[i-1])/2 if i > 0 else s[i] - (s[1]-s[0])/2
        hi = s[i] + (s[i+1]-s[i])/2 if i < n-1 else s[i] + (s[-1]-s[-2])/2
        out.append(max(hi-lo, 1e-6))
    return out


def count_path_waters(centres, raw_clearances, oxygens, *, tree=None):
    """Count water oxygens per path sample (nearest sample, inside its free sphere).

    Parameters
    ----------
    centres : array_like, shape (n, 3)
        Path sample centres, Å.
    raw_clearances : array_like, shape (n,)
        Free-sphere radius (raw clearance) at each sample in this frame, Å.
    oxygens : array_like, shape (m, 3)
        Aligned water-oxygen positions, Å.
    tree : scipy.spatial.cKDTree, optional
        Prebuilt tree of ``centres``.

    Returns
    -------
    numpy.ndarray, shape (n,)
        Oxygens assigned to each sample: nearest sample ``i`` with centre
        distance at most ``raw_clearances[i]``.

    Examples
    --------
    >>> count_path_waters([(0, 0, 0), (0, 0, 2)], [1.0, 1.0], [(0, 0, 0.5), (0, 0, 2.2), (3, 0, 0)]).tolist()
    [1.0, 1.0]
    """
    import numpy as np
    from scipy.spatial import cKDTree
    centres = np.asarray(centres, dtype=float)
    raw = np.asarray(raw_clearances, dtype=float)
    counts = np.zeros(len(centres))
    oxygens = np.asarray(oxygens, dtype=float).reshape(-1, 3)
    if not len(oxygens) or not len(centres):
        return counts
    tree = tree if tree is not None else cKDTree(centres)
    dist, idx = tree.query(oxygens, distance_upper_bound=float(max(raw.max(), 0.0)) + 1e-6)
    ok = np.isfinite(dist)
    idx, d = idx[ok], dist[ok]
    inside = d <= raw[idx]
    np.add.at(counts, idx[inside], 1)
    return counts


def fingerprint_trajectory(topology, trajectory, *, selection: str = "protein", start: int = 0, stop: int | None = None,
                           stride: int = 1, max_frames: int | None = None, pbc: str = "check",
                           water_selection: str | None = None, include_water: bool = True,
                           quantiles=(0.1, 0.9), progress=None, **fingerprint_kwargs) -> PoreFingerprint:
    """Fingerprint along a fixed reference path through every frame of an MD trajectory.

    The channel, its cast, segments and lining residues are resolved once on
    the first loaded frame (the reference). Every frame is rigidly fitted to
    it (:func:`crevice.trajectory.align_frame_report`); the clearance is then
    measured at the fixed reference path samples
    (:func:`crevice.channels.profile_along_path`) and water oxygens are counted
    along the path (module docstring). The path is not re-resolved per frame,
    so a frame in which the channel moves away from the reference centre line
    reports the clearance at the reference positions.

    Parameters
    ----------
    topology, trajectory : str or Path
        MD topology and trajectory read by MDAnalysis.
    selection : str, default "protein"
        MDAnalysis selection of the protein atoms (the geometry).
    start, stop, stride, max_frames : int
        Frame range and guard (:func:`crevice.trajectory.load_trajectory_report`).
    pbc : {"check", "none"}, default "check"
        Periodic handling; ``"unwrap"`` is refused when water is counted.
    water_selection : str, optional
        MDAnalysis water selection; default
        :data:`crevice.hydration.DEFAULT_WATER_SELECTION`.
    include_water : bool, default True
        Count water oxygens along the path.
    quantiles : tuple of float, default (0.1, 0.9)
        Lower and upper frame quantiles of the radius band.
    progress : callable, optional
        Called with the number of frames done after each frame.
    **fingerprint_kwargs
        Passed to :func:`fingerprint_structure` (profile settings, entry end,
        bandwidth, shell, electrolyte, resistivity).

    Returns
    -------
    PoreFingerprint
        Reference-frame fingerprint whose samples also carry
        :data:`ENSEMBLE_FIELDS`, ``conductance`` gains per-frame statistics,
        ``frames`` lists every frame's time, fit RMSD, lumen minimum radius,
        ``G_macro`` and pore water count, and ``water`` records the water
        definitions.

    Raises
    ------
    ValueError
        If the reference channel is unresolved or ``pbc="unwrap"`` with water.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    from .channels import profile_along_path
    from .trajectory import align_frame_report, load_trajectory_report
    traj, reader = load_trajectory_report(topology, trajectory, selection=selection, start=start, stop=stop,
                                          stride=stride, max_frames=max_frames, pbc=pbc)
    reference = traj.frames[0]
    profile, _cast, fp = fingerprint_structure(reference, **fingerprint_kwargs)
    rows = fp.samples
    centres = [(r["x_A"], r["y_A"], r["z_A"]) for r in rows]
    s = [r["s_A"] for r in rows]
    widths = np.asarray(_bin_widths(s))
    lumen_idx = [i for i, r in enumerate(rows) if r["segment"] == "lumen"]
    waters = None
    water_meta = None
    if include_water:
        if pbc == "unwrap":
            raise ValueError("water counting needs pbc check or none; unwrap moves protein independently of solvent")
        from .water_membership import WaterSource, image_waters_around_anchor
        waters = WaterSource(topology, trajectory, protein_selection=selection, reference_frame=reference,
                             frame_times=dict(zip(range(*reader["frame_indices"]), traj.times())),
                             water_selection=water_selection, pbc=pbc)
        if not len(waters.oxygen_indices):
            waters = None
        water_meta = {**(waters.metadata if waters else {"water_population": 0}),
                      "definition": "oxygen assigned to its nearest path sample and counted when inside that frame's "
                                    "free sphere there (centre distance <= raw clearance); ds = half the distance to "
                                    "neighbouring samples; number density = mean count / mean(pi R^2 ds)",
                      "bulk_water_density_per_A3": BULK_WATER_DENSITY_A3}
    tree = cKDTree(np.asarray(centres))
    anchor = np.asarray(centres).mean(0)
    radius_frames, raw_frames, count_frames, frames = [], [], [], []
    rho = fp.conductance["resistivity_ohm_m"]
    include_h = fingerprint_kwargs.get("include_hydrogen", False)
    include_het = fingerprint_kwargs.get("include_hetero", True)
    for n, (frame, time) in enumerate(zip(traj.frames, traj.times())):
        aligned, report = align_frame_report(frame, reference)
        measured = profile_along_path(aligned, centres, probe_radius=profile.probe_radius,
                                      include_hydrogen=include_h, include_hetero=include_het)
        radius = np.asarray([p.radius for p in measured.points])
        raw = np.asarray([p.raw_clearance for p in measured.points])
        counts = np.zeros(len(rows))
        if waters is not None:
            o, box = waters.snapshot(frame)
            rel = image_waters_around_anchor(o, box, np.asarray(report["rotation_rows"]),
                                             np.asarray(report["translation_A"]), anchor) + anchor
            counts = count_path_waters(centres, raw, rel, tree=tree)
        g = conductance_estimate([s[i] for i in lumen_idx], [float(radius[i]) for i in lumen_idx],
                                 resistivity_ohm_m=rho, electrolyte=fp.conductance["electrolyte"])
        radius_frames.append(radius)
        raw_frames.append(raw)
        count_frames.append(counts)
        frames.append({"frame_index": frame.frame_index, "time_ps": time, "fit_rmsd_A": report["rmsd_after_A"],
                       "lumen_min_radius_A": float(radius[lumen_idx].min()), "g_macro_pS": g["g_macro_pS"],
                       "g_access_pS": g["g_access_pS"],
                       "pore_water_count": int(counts[lumen_idx].sum()) if waters is not None else None})
        if progress is not None:
            progress(n+1)
    R = np.asarray(radius_frames)
    raw_all = np.clip(np.asarray(raw_frames), 0.0, None)
    C = np.asarray(count_frames)
    lo, hi = np.quantile(R, quantiles[0], axis=0), np.quantile(R, quantiles[1], axis=0)
    volume = (math.pi*raw_all**2*widths).mean(0)
    for i, row in enumerate(rows):
        row.update({"radius_mean_A": float(R[:, i].mean()), "radius_sd_A": float(R[:, i].std()),
                    "radius_q_low_A": float(lo[i]), "radius_q_high_A": float(hi[i])})
        if waters is not None:
            mean = float(C[:, i].mean())
            row.update({"water_count_mean": mean, "water_linear_density_per_A": mean/widths[i],
                        "water_number_density_per_A3": mean/volume[i] if volume[i] > 0 else None,
                        "water_occupancy": float((C[:, i] > 0).mean())})
        else:
            row.update({"water_count_mean": None, "water_linear_density_per_A": None,
                        "water_number_density_per_A3": None, "water_occupancy": None})
    gs = np.asarray([f["g_macro_pS"] for f in frames])
    mean_profile = conductance_estimate([s[i] for i in lumen_idx], [float(R[:, i].mean()) for i in lumen_idx],
                                        resistivity_ohm_m=rho, electrolyte=fp.conductance["electrolyte"])
    fp.conductance = {**fp.conductance, "reference_frame": "first loaded frame",
                      "per_frame_g_macro_pS": {"mean": float(gs.mean()), "median": float(np.median(gs)),
                                               "q_low": float(np.quantile(gs, quantiles[0])),
                                               "q_high": float(np.quantile(gs, quantiles[1])),
                                               "min": float(gs.min()), "max": float(gs.max()),
                                               "occluded_frames": int((gs == 0).sum())},
                      "g_macro_of_mean_radius_pS": mean_profile["g_macro_pS"]}
    fp.frames = frames
    fp.water = water_meta if waters is not None else ({**water_meta, "note": "no water oxygens selected"} if water_meta else None)
    fp.parameters = {**fp.parameters, "trajectory": {
        "topology": str(topology), "trajectory": str(trajectory), "selection": selection, "frames": len(frames),
        "quantiles": list(quantiles), "reader": reader,
        "path": "fixed reference path of the first loaded frame; every frame fitted to it (Kabsch, CA atoms)"}}
    return fp


def write_fingerprint_csv(fp: PoreFingerprint, path) -> Path:
    """Write the per-sample fingerprint table.

    Parameters
    ----------
    fp : PoreFingerprint
        Result of :func:`pore_fingerprint` or :func:`fingerprint_trajectory`.
    path : str or Path
        Output CSV; columns :data:`SAMPLE_FIELDS` (plus :data:`ENSEMBLE_FIELDS`
        when present), unit-suffixed; empty cells are missing values.

    Returns
    -------
    Path
        The written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(SAMPLE_FIELDS) + [f for f in ENSEMBLE_FIELDS if fp.samples and f in fp.samples[0]]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in fp.samples:
            writer.writerow({k: ("" if row.get(k) is None else round(row[k], 6) if isinstance(row.get(k), float) else row.get(k))
                             for k in fields})
    return path


def write_fingerprint_bundle(fp: PoreFingerprint, output_dir, prefix: str = "crevice", *, dpi: int = 200) -> dict:
    """Write ``PREFIX_fingerprint.csv``, ``.json`` and ``.png``.

    Parameters
    ----------
    fp : PoreFingerprint
        The fingerprint.
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix.
    dpi : int, default 200
        Figure resolution.

    Returns
    -------
    dict
        ``fingerprint_csv``, ``fingerprint_json``, ``fingerprint_png`` (and
        ``fingerprint_frames_csv`` for trajectories) as path strings.
    """
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    out = {"fingerprint_csv": str(write_fingerprint_csv(fp, root/f"{prefix}_fingerprint.csv"))}
    path = root/f"{prefix}_fingerprint.json"
    path.write_text(json.dumps(fp.to_dict(), indent=2, default=str)+"\n")
    out["fingerprint_json"] = str(path)
    if fp.frames:
        frames = root/f"{prefix}_fingerprint_frames.csv"
        with frames.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fp.frames[0]))
            writer.writeheader()
            writer.writerows(fp.frames)
        out["fingerprint_frames_csv"] = str(frames)
    out["fingerprint_png"] = str(plot_pore_fingerprint(fp, root/f"{prefix}_fingerprint.png", dpi=dpi))
    return out


#: Fingerprint panel colours. Radius line blue, hydrophobic and hydrophilic
#: fills orange and blue (a diverging pair around zero), water blue, strip
#: ticks dark ink.
FINGERPRINT_COLOURS = {"radius": "#2a78d6", "hydrophobic": "#eb6834", "hydrophilic": "#2a78d6",
                       "water": "#2a78d6", "ink": "#3d3d3a", "zero": "#8d8d88"}


@annotate_option
def plot_pore_fingerprint(fp: PoreFingerprint, path, *, dpi: int = 200, annotate: bool | None = None) -> Path:
    """Stacked panels along ``s``: radius, hydrophobicity, water and chemistry.

    Parameters
    ----------
    fp : PoreFingerprint
        The fingerprint.
    path : str or Path
        Output image (PNG recommended).
    dpi : int, default 200
        Raster resolution.
    annotate : bool, optional
        Draw the title and the conductance call-out; ``None`` inherits the
        surrounding setting (off unless ``--annotate``).

    Returns
    -------
    Path
        The written image.

    Outputs
    -------
    path : image
        7.0 inch wide figure; panels share the x axis "Path coordinate s (Å)":
        radius (Å), lining hydrophobicity (kcal/mol), water (molecules per Å,
        trajectories only) and a four-row chemistry strip.

    Colours and representation
    --------------------------
    * Segment shading across every panel: ENTRY yellow ``#f0d43a``, LUMEN teal
      ``#1494a1``, EXIT rust ``#c2410c`` (the cast-segment colours), light
      alpha; a legend names the segments present.
    * Radius: blue line; for trajectories the frame mean with a light band
      between the lower and upper frame quantiles, and the reference frame as
      a thin grey line.
    * Hydrophobicity: dark line; fill above zero orange (hydrophobic), below
      zero blue (hydrophilic); grey zero line.
    * Water: blue line of the mean linear density.
    * Chemistry strip: one row per class (hydrophobic, polar, charged,
      aromatic) with a dark tick at each lining residue's anchor ``s``; the
      row label carries identity, not colour.
    """
    from .cast_segments import SEGMENT_HEX
    from .figures import _ensure_parent, _pyplot, _save_figure, _style_axes
    plt = _pyplot()
    from matplotlib.patches import Patch
    output = _ensure_parent(path)
    rows = fp.samples
    s = [r["s_A"] for r in rows]
    g = fp.conductance
    text = FigureText(
        "Pore fingerprint along the path coordinate: radius, lining hydrophobicity (Wimley-White octanol, "
        f"Gaussian sigma {fp.parameters['hydrophobicity_bandwidth_A']} A)"
        + (", water linear density" if fp.water and fp.water.get("water_population") else "")
        + f" and lining chemistry; geometric G_macro {g['g_macro_pS']:.3g} pS in {g['electrolyte']} "
          f"(order-of-magnitude estimate)")
    has_water = bool(rows) and rows[0].get("water_linear_density_per_A") is not None
    ratios = [2.2, 1.6] + ([1.4] if has_water else []) + [1.0]
    fig, axes = plt.subplots(len(ratios), 1, figsize=(7.0, 1.55*sum(ratios)+0.6), sharex=True,
                             gridspec_kw={"height_ratios": ratios}, constrained_layout=True)
    c = FINGERPRINT_COLOURS
    present = []
    lumen_span = next(((seg["s_start_A"], seg["s_end_A"]) for seg in fp.segments if seg["segment"] == "lumen"), None)
    for seg in fp.segments:
        kind = seg["segment_type"]
        lo, hi = seg["s_start_A"], seg["s_end_A"]
        # Shade contiguous spans: an end segment meets the lumen at its mouth.
        if lumen_span and kind == "entry":
            hi = lumen_span[0] if hi < lumen_span[0] else hi
            lo = lumen_span[1] if lo > lumen_span[1] else lo
        elif lumen_span and kind == "exit":
            lo = lumen_span[1] if lo > lumen_span[1] else lo
            hi = lumen_span[0] if hi < lumen_span[0] else hi
        if kind not in present:
            present.append(kind)
        for ax in axes:
            ax.axvspan(lo, hi, color=SEGMENT_HEX[kind], alpha=0.16 if kind != "lumen" else 0.10, linewidth=0, zorder=0)
    ax = axes[0]
    if "radius_mean_A" in rows[0]:
        ax.fill_between(s, [r["radius_q_low_A"] for r in rows], [r["radius_q_high_A"] for r in rows],
                        color=c["radius"], alpha=0.22, linewidth=0, label="Frame quantile band")
        ax.plot(s, [r["radius_A"] for r in rows], color=c["zero"], linewidth=1.0, label="Reference frame")
        ax.plot(s, [r["radius_mean_A"] for r in rows], color=c["radius"], linewidth=2.0, label="Frame mean")
    else:
        ax.plot(s, [r["radius_A"] for r in rows], color=c["radius"], linewidth=2.0, label="Radius")
    ax.set_ylabel("Radius (Å)")
    ax.set_ylim(bottom=0)
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) < 2:
        handles, labels = [], []
    handles += [Patch(color=SEGMENT_HEX[k], alpha=0.35 if k != "lumen" else 0.25) for k in present]
    labels += [k.upper() for k in present]
    ax.legend(handles, labels, frameon=False, fontsize=8, ncol=min(len(labels), 3), loc="upper left")
    text.title(ax, "CREVICE pore fingerprint")
    text.text(ax, 0.99, 0.95, f"G_macro ≈ {g['g_macro_pS']:.2g} pS ({g['electrolyte']}); geometric estimate",
              transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax = axes[1]
    h = [r["hydrophobicity_kcal_mol"] for r in rows]
    hs = [v if v is not None else float("nan") for v in h]
    ax.axhline(0, color=c["zero"], linewidth=0.8)
    import numpy as np
    arr = np.asarray(hs, dtype=float)
    ax.fill_between(s, arr, 0, where=arr >= 0, color=c["hydrophobic"], alpha=0.35, linewidth=0, interpolate=True)
    ax.fill_between(s, arr, 0, where=arr < 0, color=c["hydrophilic"], alpha=0.30, linewidth=0, interpolate=True)
    ax.plot(s, arr, color=c["ink"], linewidth=1.4)
    ax.set_ylabel("Hydrophobicity\n(kcal/mol)")
    finite = arr[np.isfinite(arr)]
    if finite.size:
        bound = max(abs(finite).max(), 0.5)*1.15
        ax.set_ylim(-bound, bound)
    else:
        text.notice(ax, 0.5, 0.5, "no lining residues with a scale value", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8)
    k = 2
    if has_water:
        ax = axes[k]
        ax.plot(s, [r["water_linear_density_per_A"] for r in rows], color=c["water"], linewidth=1.8)
        ax.set_ylabel("Water\n(molecules/Å)")
        ax.set_ylim(bottom=0)
        k += 1
    ax = axes[k]
    for row_index, cls in enumerate(CHEMISTRY_CLASSES):
        xs = [e["anchor_s_A"] for e in fp.lining if cls in e["classes"]]
        y = len(CHEMISTRY_CLASSES) - 1 - row_index
        ax.vlines(xs, y-0.32, y+0.32, color=c["ink"], linewidth=1.3)
    ax.set_yticks(range(len(CHEMISTRY_CLASSES)))
    ax.set_yticklabels([cls.capitalize() for cls in reversed(CHEMISTRY_CLASSES)], fontsize=8)
    ax.set_ylim(-0.6, len(CHEMISTRY_CLASSES)-0.4)
    ax.set_xlabel("Path coordinate s (Å)")
    ax.set_xlim(min(s), max(s))
    for ax in axes:
        _style_axes(ax)
    for ax in axes[:-1]:
        ax.grid(True, axis="y", color="#e2e2df", linewidth=0.6)
    text.note(f"G_macro {g['g_macro_pS']:.4g} pS, G_access {g['g_access_pS']:.4g} pS, HOLE-corrected "
              f"{g['g_hole_corrected_pS']:.4g} pS ({g['electrolyte']}, rho {g['resistivity_ohm_m']} Ohm m)")
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


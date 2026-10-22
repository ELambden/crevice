"""Permeation events, water permeabilities and ion current through a channel in an MD trajectory.

The analysis follows one representative atom per molecule (the water oxygen
or Martini water bead, the ion) of every requested species and assigns it, in
every frame, to one of four states relative to each channel:

* ``below`` (axial coordinate under the channel extent),
* ``inside`` the lumen (within the extent and within the lumen's lateral
  limit),
* ``above`` (axial coordinate over the extent),
* ``out``: within the axial extent but outside the lumen laterally (for
  example in the protein or the bilayer next to the pore).

A *permeation event* is a complete crossing: the molecule leaves one bulk
side, stays inside the lumen without interruption and leaves to the opposite
side (below -> inside -> above is ``up``, the reverse is ``down``), the
definition of Zhu, Tajkhorshid and Schulten (Biophys. J. 86, 50-57, 2004).
A molecule that enters and leaves on the same side, or that leaves the lumen
laterally (``out``), does not permeate; ``out`` also clears the side the
molecule came from. Periodic boundaries are handled per molecule: positions
are taken as the minimum image about the channel centre, and the frame to
frame displacement is unwrapped, so a molecule that reaches the other side
through the periodic bulk is not a permeation, and one that crosses the whole
extent between two saved frames is reported as an unresolved crossing (not
counted as an event).

Channel definitions (all in the coordinates of the reference frame, the first
analysed frame after least-squares fitting of the protein alignment atoms):

* :class:`CylinderChannel`: an explicit axial extent along an axis, with an
  optional lateral radius around an axis line (several of them describe
  several pores, for example the four pores of an aquaporin tetramer);
* :class:`ProfileChannel`: the per-frame channel profiles of
  :func:`crevice.trajectory.analyze_trajectory`; the extent is the profile's
  first and last sample along its axis and the lumen radius is the profile
  radius (plus a tolerance) at the molecule's axial position;
* :class:`RegionChannel`: a fixed region map (OpenDX), for example the
  reference map of a prepared named region.

Derived quantities (:func:`analyze_permeation`):

* counts per direction, net flux, unidirectional rates with exact Poisson
  (Garwood) intervals, occupancy (mean with a block-bootstrap interval) and
  dwell/transit times;
* water: the collective coordinate n(t) (dn = sum of axial displacements
  inside the channel divided by the channel length; only the part of a
  displacement inside the extent counts), its diffusion coefficient D_n from
  the mean-square displacement, the osmotic permeability ``pf = v_w D_n`` and
  the diffusive permeability ``pd = v_w q0`` with ``q0`` the unidirectional
  event rate (Zhu, Tajkhorshid and Schulten, Biophys. J. 2004; Phys. Rev.
  Lett. 93, 224501, 2004);
* ions: the current ``I = e * sum(q * (up - down)) / T`` with Poisson and
  block-bootstrap intervals, and the conductance ``G = I / V`` only when a
  voltage is known (given, or an applied field times the box length).

Martini input: ``W`` (and ``WF``, ``PW``) water beads stand for four water
molecules; counts, ``q0`` and ``pd`` are multiplied by the mapping and the
collective coordinate is expressed in water molecules (``D_n`` multiplied by
the mapping squared) before ``pf = v_w D_n``. Martini 2 ``ION`` residues take
their charge from the atom name (``NA+``, ``CL-``). Coarse-grained input is
flagged in the summary. Coarse-grained time scales are not converted.

Outputs (:func:`write_permeation_outputs`): ``PREFIX_permeation.json`` (the
summary with definitions, units and uncertainties), ``PREFIX_permeation_events.csv``,
``PREFIX_permeation_visits.csv``, ``PREFIX_permeation_frames.csv`` (per-frame
occupancy, cumulative net crossings and n(t)) and the figure
``PREFIX_permeation.png`` (:func:`plot_permeation`).
"""
from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np

from .presentation import FigureText, annotate_option

#: Molar volume of liquid water used by Zhu et al. (2004), cm^3/mol.
WATER_MOLAR_VOLUME_CM3 = 18.0
#: Avogadro constant, 1/mol.
AVOGADRO = 6.02214076e23
#: Volume of one water molecule ``v_w``, cm^3.
WATER_MOLECULE_VOLUME_CM3 = WATER_MOLAR_VOLUME_CM3 / AVOGADRO
#: Elementary charge, C.
ELEMENTARY_CHARGE_C = 1.602176634e-19

#: Residue names read as water. The representative atom is the oxygen (or the bead).
WATER_RESNAMES = frozenset({"TIP3", "TIP3P", "TP3", "T3P", "TIP4", "TIP4P", "T4P", "TIP5", "TIP5P", "T5P", "SOL",
                            "WAT", "HOH", "H2O", "SPC", "SPCE", "OPC", "OPC3", "SWM4", "W", "WF", "PW"})
#: Martini water residues and the number of water molecules per bead.
MARTINI_WATER_MAPPING = {"W": 4, "WF": 4, "PW": 4}
#: Residue names of Martini 2 ions, whose charge is read from the atom name.
MARTINI_ION_RESNAMES = frozenset({"ION"})
#: Atom names preferred as the representative atom of a multi-atom molecule.
REPRESENTATIVE_NAMES = ("OH2", "OW", "O", "OW1", "OH", "W")

#: Default charges (e) by residue name (or, for Martini ``ION``, atom name).
#: CHARMM (SOD, POT, CLA, CAL, ...), AMBER (Na+, K+, Cl-, ...), GROMACS/OPLS and
#: Martini names. Upper-case keys; lookup is case-insensitive.
DEFAULT_CHARGES = {
    "NA": 1.0, "SOD": 1.0, "NA+": 1.0, "NA1+": 1.0,
    "K": 1.0, "POT": 1.0, "K+": 1.0,
    "LI": 1.0, "LIT": 1.0, "LI+": 1.0,
    "RB": 1.0, "RB+": 1.0, "CS": 1.0, "CES": 1.0, "CS+": 1.0,
    "CL": -1.0, "CLA": -1.0, "CL-": -1.0, "CL1-": -1.0,
    "BR": -1.0, "BR-": -1.0, "IOD": -1.0, "I-": -1.0,
    "MG": 2.0, "MG2+": 2.0, "MG+2": 2.0,
    "CA": 2.0, "CAL": 2.0, "CA2+": 2.0, "CA+2": 2.0,
    "ZN": 2.0, "ZN2": 2.0, "ZN2+": 2.0,
    "BA": 2.0, "BAR": 2.0, "BA2+": 2.0,
}

_STATE_BELOW, _STATE_INSIDE, _STATE_ABOVE, _STATE_OUT = -1, 0, 1, 2

LIMITATIONS = [
    "An event needs the molecule to be seen inside the lumen in at least one saved frame; a crossing completed "
    "between two saved frames is reported as unresolved and not counted.",
    "Entry, exit and transit times are resolved only to the saved-frame interval.",
    "The channel extent and lumen are geometric definitions; counts depend on them (and on the lateral tolerance).",
    "pf from the collective-diffusion model and pd from the event rate assume a stationary equilibrium trajectory; "
    "their intervals are sampling intervals only (no force-field or model error).",
    "The ion current counts complete crossings of the channel extent; conductance is reported only for a known "
    "voltage and is not a validated single-channel conductance.",
    "Coarse-grained (Martini) time scales are not converted; rates and permeabilities are in simulation time.",
]


# ---------------------------------------------------------------------------
# Channel definitions
# ---------------------------------------------------------------------------

def _unit(vector) -> np.ndarray:
    v = np.asarray(vector, dtype=float).reshape(3)
    norm = float(np.linalg.norm(v))
    if not math.isfinite(norm) or norm == 0:
        raise ValueError("channel axis must be a finite non-zero vector")
    return v / norm


def _lateral(vectors: np.ndarray, axis: np.ndarray) -> np.ndarray:
    return np.linalg.norm(vectors - np.outer(vectors @ axis, axis), axis=1)


@dataclass
class CylinderChannel:
    """A channel with an explicit axial extent and an optional lateral radius.

    Attributes
    ----------
    name : str
        Label used in the outputs (for example ``"pore1"``).
    axis : tuple of float
        Channel axis (normalised); ``above`` is the ``+axis`` side.
    zmin, zmax : float
        Axial extent, Å, as projections ``x . axis`` of reference-frame
        coordinates.
    point : tuple of float or None, default None
        A point on the axis line (reference frame, Å). ``None`` uses the
        origin passed by the caller (:func:`analyze_permeation` uses the
        centroid of the alignment atoms).
    radius : float or None, default None
        Lateral radius of the lumen around the axis line, Å. ``None`` means no
        lateral limit: every molecule within the extent is inside.
    """

    name: str
    axis: tuple
    zmin: float
    zmax: float
    point: tuple | None = None
    radius: float | None = None

    def __post_init__(self):
        self.axis = tuple(_unit(self.axis))
        if not (math.isfinite(self.zmin) and math.isfinite(self.zmax)) or self.zmax <= self.zmin:
            raise ValueError("channel extent needs finite zmin < zmax")
        if self.radius is not None and not (math.isfinite(self.radius) and self.radius > 0):
            raise ValueError("channel radius must be positive")

    kind = "cylinder"

    def bind(self, default_point) -> None:
        """Fix the axis line through ``default_point`` when no point was given.

        Parameters
        ----------
        default_point : array_like, shape (3,)
            Point used when :attr:`point` is ``None``.
        """
        if self.point is None:
            self.point = tuple(float(v) for v in default_point)

    def centre(self) -> np.ndarray:
        """Point of the axis line at the middle of the extent (reference frame, Å).

        Returns
        -------
        numpy.ndarray, shape (3,)
        """
        a = np.asarray(self.axis)
        p = np.asarray(self.point, dtype=float)
        return p + ((self.zmin + self.zmax) / 2 - p @ a) * a

    def extent(self, frame: int) -> tuple:
        """Axial extent ``(zmin, zmax)`` in Å for analysed frame ``frame``.

        Parameters
        ----------
        frame : int
            Analysed-frame index (ignored: the extent is fixed).
        """
        return self.zmin, self.zmax

    def fixed_extent(self) -> tuple:
        """Extent used for the collective coordinate n(t) (fixed), Å.

        Returns
        -------
        tuple of float
        """
        return self.zmin, self.zmax

    def lumen(self, positions: np.ndarray, frame: int) -> np.ndarray:
        """Whether reference-frame ``positions`` inside the extent are in the lumen.

        Parameters
        ----------
        positions : numpy.ndarray, shape (n, 3)
            Positions already known to lie within the axial extent, Å.
        frame : int
            Analysed-frame index (ignored).

        Returns
        -------
        numpy.ndarray of bool
        """
        if self.radius is None:
            return np.ones(len(positions), dtype=bool)
        return _lateral(positions - np.asarray(self.point), np.asarray(self.axis)) <= self.radius

    def describe(self) -> dict:
        """JSON record of the definition.

        Returns
        -------
        dict
        """
        return {"kind": self.kind, "name": self.name, "axis": list(self.axis), "extent_A": [self.zmin, self.zmax],
                "axis_point_A": list(self.point) if self.point is not None else None,
                "lateral_radius_A": self.radius,
                "inside_rule": ("zmin <= x.axis <= zmax and distance from the axis line <= lateral_radius_A"
                                if self.radius is not None else "zmin <= x.axis <= zmax (no lateral limit)")}


class ProfileChannel:
    """A channel defined by per-frame radius profiles.

    The extent of a frame is the axial range of its profile samples; a
    molecule within it is inside the lumen when its distance from the profile
    centre line (interpolated at its axial position) is at most the profile
    radius there plus ``tolerance``. Analysed frames without a resolved
    profile use the most recent earlier profile (or the first one), and those
    frames are counted (``held_frames``).

    Parameters
    ----------
    name : str
        Label used in the outputs.
    profiles : sequence of PoreProfile or None
        Profiles in the reference frame, in trajectory order. ``None`` marks a
        frame without a profile.
    frame_numbers : sequence of int, optional
        Trajectory frame number of each profile; analysed frames are mapped to
        the latest profile at or before them. Default ``0, 1, 2, ...``.
    tolerance : float, default 1.0
        Added to the profile radius, Å.
    axis : array_like, optional
        Channel axis; default the first profile's ``axis_direction``.

    Raises
    ------
    ValueError
        If no profile is resolved.
    """

    kind = "profile"

    def __init__(self, name, profiles, frame_numbers=None, *, tolerance=1.0, axis=None):
        self.name = name
        frames = list(range(len(profiles))) if frame_numbers is None else [int(f) for f in frame_numbers]
        if len(frames) != len(profiles):
            raise ValueError("frame_numbers must match profiles")
        resolved = [(f, p) for f, p in zip(frames, profiles) if p is not None and len(p.points) >= 2]
        if not resolved:
            raise ValueError("no resolved channel profile: give an explicit extent (--channel-extent) instead")
        self.axis = _unit(axis if axis is not None else resolved[0][1].axis_direction)
        self.tolerance = float(tolerance)
        self._numbers = np.asarray([f for f, _ in resolved])
        self._tables = []
        for _, profile in resolved:
            pos = np.asarray([p.position for p in profile.points], dtype=float)
            s = pos @ self.axis
            order = np.argsort(s, kind="stable")
            s, pos = s[order], pos[order]
            radius = np.asarray([p.radius for p in profile.points], dtype=float)[order]
            keep = np.concatenate([[True], np.diff(s) > 1e-9])
            self._tables.append((s[keep], pos[keep], radius[keep]))
        self.total_profiles = len(profiles)
        self.resolved_profiles = len(resolved)
        self._frame_numbers = None
        self.held_frames = 0
        first = self._tables[0]
        self._centre = first[1][len(first[1]) // 2]

    def bind_frames(self, frame_numbers: Sequence[int]) -> None:
        """Map each analysed frame to a profile (latest at or before it).

        Parameters
        ----------
        frame_numbers : sequence of int
            Trajectory frame number of each analysed frame.
        """
        numbers = np.asarray(frame_numbers, dtype=int)
        index = np.searchsorted(self._numbers, numbers, side="right") - 1
        index[index < 0] = 0
        self._map = index
        self.held_frames = int(np.sum(self._numbers[index] != numbers))

    def bind(self, default_point) -> None:
        """No-op (the profile fixes its own centre line).

        Parameters
        ----------
        default_point : array_like
            Ignored.
        """

    def _table(self, frame):
        return self._tables[int(self._map[frame]) if hasattr(self, "_map") else 0]

    def centre(self) -> np.ndarray:
        """Middle sample of the first resolved profile (reference frame, Å).

        Returns
        -------
        numpy.ndarray, shape (3,)
        """
        return np.asarray(self._centre)

    def extent(self, frame: int) -> tuple:
        """Axial range ``(zmin, zmax)`` of the frame's profile, Å.

        Parameters
        ----------
        frame : int
            Analysed-frame index.
        """
        s = self._table(frame)[0]
        return float(s[0]), float(s[-1])

    def fixed_extent(self) -> tuple:
        """Median extent over the resolved profiles (for n(t)), Å.

        Returns
        -------
        tuple of float
        """
        lows = [t[0][0] for t in self._tables]
        highs = [t[0][-1] for t in self._tables]
        return float(np.median(lows)), float(np.median(highs))

    def lumen(self, positions: np.ndarray, frame: int) -> np.ndarray:
        """Lateral test against the frame's profile radius plus the tolerance.

        Parameters
        ----------
        positions : numpy.ndarray, shape (n, 3)
            Reference-frame positions within the extent, Å.
        frame : int
            Analysed-frame index.

        Returns
        -------
        numpy.ndarray of bool
        """
        s_table, centres, radii = self._table(frame)
        s = positions @ self.axis
        centre = np.column_stack([np.interp(s, s_table, centres[:, k]) for k in range(3)])
        radius = np.interp(s, s_table, radii)
        return _lateral(positions - centre, self.axis) <= radius + self.tolerance

    def describe(self) -> dict:
        """JSON record of the definition.

        Returns
        -------
        dict
        """
        lo, hi = self.fixed_extent()
        return {"kind": self.kind, "name": self.name, "axis": self.axis.tolist(),
                "median_extent_A": [lo, hi], "lumen_tolerance_A": self.tolerance,
                "resolved_profiles": self.resolved_profiles, "total_profiles": self.total_profiles,
                "held_frames": self.held_frames,
                "inside_rule": "within the frame's profile extent and within profile radius + tolerance of the "
                               "profile centre line at the same axial position",
                "missing_profile_rule": "latest earlier resolved profile (or the first one)"}


class RegionChannel:
    """A channel defined by a fixed region map (OpenDX grid) and an axis.

    Grid points with a value above ``level`` form the region; it is dilated by
    ``tolerance`` (Å). A molecule is inside when its reference-frame position
    falls in a region voxel; ``below``/``above`` are the axial ranges outside
    the region's projected extent.

    Parameters
    ----------
    name : str
        Label used in the outputs.
    path : str or Path
        OpenDX map in the coordinates of the reference frame.
    axis : array_like
        Channel axis.
    tolerance : float, default 0.0
        Dilation of the region, Å.
    level : float, default 0.5
        Threshold of the map.
    """

    kind = "region"

    def __init__(self, name, path, axis, *, tolerance=0.0, level=0.5):
        from .viewer import read_dx_field
        values, origin, deltas = read_dx_field(path)
        if not np.allclose(deltas, np.diag(np.diag(deltas))):
            raise ValueError("region maps must use an axis-aligned grid")
        self.name, self.path, self.axis = name, str(path), _unit(axis)
        self.tolerance, self.level = float(tolerance), float(level)
        self.origin, self.spacing = np.asarray(origin, float), np.diag(deltas).astype(float)
        mask = values > level
        if not mask.any():
            raise ValueError(f"{path}: the region map has no point above {level}")
        if tolerance > 0:
            from scipy.ndimage import binary_dilation
            steps = np.ceil(tolerance / self.spacing).astype(int)
            grid = np.ogrid[tuple(slice(-k, k + 1) for k in steps)]
            ball = sum((g * h) ** 2 for g, h in zip(grid, self.spacing)) <= tolerance ** 2
            mask = binary_dilation(np.pad(mask, [(k, k) for k in steps]), structure=ball)
            self.origin = self.origin - steps * self.spacing
        self.mask = mask
        points = np.argwhere(mask) * self.spacing + self.origin
        s = points @ self.axis
        self._extent = (float(s.min() - 0.5 * np.abs(self.spacing @ np.abs(self.axis))),
                        float(s.max() + 0.5 * np.abs(self.spacing @ np.abs(self.axis))))
        self._centre = points.mean(axis=0)
        self.voxels = int(mask.sum())

    def bind(self, default_point) -> None:
        """No-op (the region fixes its own centre).

        Parameters
        ----------
        default_point : array_like
            Ignored.
        """

    def centre(self) -> np.ndarray:
        """Centroid of the region points (reference frame, Å).

        Returns
        -------
        numpy.ndarray, shape (3,)
        """
        return self._centre

    def extent(self, frame: int) -> tuple:
        """Projected axial extent of the region, Å.

        Parameters
        ----------
        frame : int
            Analysed-frame index (ignored).
        """
        return self._extent

    def fixed_extent(self) -> tuple:
        """Extent used for n(t), Å.

        Returns
        -------
        tuple of float
        """
        return self._extent

    def lumen(self, positions: np.ndarray, frame: int) -> np.ndarray:
        """Whether positions fall in a region voxel.

        Parameters
        ----------
        positions : numpy.ndarray, shape (n, 3)
            Reference-frame positions, Å.
        frame : int
            Analysed-frame index (ignored).

        Returns
        -------
        numpy.ndarray of bool
        """
        index = np.rint((positions - self.origin) / self.spacing).astype(int)
        ok = np.all((index >= 0) & (index < np.asarray(self.mask.shape)), axis=1)
        out = np.zeros(len(positions), dtype=bool)
        out[ok] = self.mask[tuple(index[ok].T)]
        return out

    def describe(self) -> dict:
        """JSON record of the definition.

        Returns
        -------
        dict
        """
        return {"kind": self.kind, "name": self.name, "map": self.path, "axis": self.axis.tolist(),
                "extent_A": list(self._extent), "dilation_A": self.tolerance, "level": self.level,
                "region_voxels": self.voxels,
                "inside_rule": "reference-frame position in a (dilated) region voxel"}


def region_channel_from_definition(path, *, name=None, tolerance=0.0) -> RegionChannel:
    """Build a :class:`RegionChannel` from a prepared named-region definition.

    Parameters
    ----------
    path : str or Path
        A prepared region definition (``REGION_region.json`` written by
        ``crevice region-prepare``); its reference map and axis are used.
    name : str, optional
        Channel label; default the region ID.
    tolerance : float, default 0.0
        Dilation of the region, Å.

    Returns
    -------
    RegionChannel

    Raises
    ------
    ValueError
        If the definition has no prepared reference map.
    """
    from .region_definition import load_definition
    data, assets = load_definition(path)
    if data["reference"]["kind"] != "existing_dx" or "reference" not in assets:
        raise ValueError("the region definition has no prepared reference map; run region-prepare first")
    channel = RegionChannel(name or data["region_id"], assets["reference"], data["geometry"]["axis"],
                            tolerance=tolerance)
    channel.region_definition = {"region_id": data["region_id"], "label": data["label"],
                                 "review_status": data["review"]["status"],
                                 "reference_frame": data["system"]["reference_frame"]}
    return channel


# ---------------------------------------------------------------------------
# Species
# ---------------------------------------------------------------------------

@dataclass
class Species:
    """One permeating species: representative atoms, charges and kind.

    Attributes
    ----------
    label : str
        Name in the outputs (the token, or ``RESNAME:ATOMNAME`` for Martini ions).
    token : str
        The ``--permeation`` token or selection it came from.
    kind : {"water", "ion", "neutral"}
        Water molecules get pf/pd; charged species the current.
    atom_indices : numpy.ndarray
        Universe atom index of the representative atom of each molecule.
    charge : float
        Charge per molecule, e.
    charge_source : str
        ``override``, ``table``, ``topology`` or ``none``.
    water_mapping : int
        Water molecules per representative (4 for Martini water beads, else 1).
    coarse_grained : bool
        Martini water bead or ion.
    representative : str
        How the representative atom was chosen.
    molecules : list of dict
        Residue identity of each molecule (``resid``, ``resname``, ``segid``).
    """

    label: str
    token: str
    kind: str
    atom_indices: np.ndarray
    charge: float
    charge_source: str
    water_mapping: int = 1
    coarse_grained: bool = False
    representative: str = ""
    molecules: list = field(default_factory=list)

    def describe(self) -> dict:
        """JSON record of the species (without per-molecule data).

        Returns
        -------
        dict
        """
        return {"label": self.label, "token": self.token, "kind": self.kind, "molecules": len(self.atom_indices),
                "charge_e": self.charge, "charge_source": self.charge_source, "water_mapping": self.water_mapping,
                "coarse_grained": self.coarse_grained, "representative_atom": self.representative}


def parse_permeation_tokens(values) -> list:
    """Split ``--permeation`` values into species tokens.

    A value containing white space is one MDAnalysis selection (optionally
    ``LABEL=selection``); otherwise it is a comma-separated list of names.

    Parameters
    ----------
    values : str or sequence of str
        ``--permeation`` values.

    Returns
    -------
    list of str

    Examples
    --------
    >>> parse_permeation_tokens(["W,NA,CL", "water=resname TIP3 and name OH2"])
    ['W', 'NA', 'CL', 'water=resname TIP3 and name OH2']
    """
    if isinstance(values, str):
        values = [values]
    tokens = []
    for value in values or ():
        value = str(value).strip()
        if not value:
            continue
        if any(c.isspace() for c in value):
            tokens.append(value)
        else:
            tokens += [t.strip() for t in value.split(",") if t.strip()]
    return list(dict.fromkeys(tokens))


def parse_charge_overrides(values) -> dict:
    """Parse ``NAME=q`` charge overrides.

    Parameters
    ----------
    values : sequence of str
        Items such as ``"SOD=1"`` or ``"NA+=1,CL-=-1"``.

    Returns
    -------
    dict of str to float
        Upper-case names to charges (e).

    Raises
    ------
    ValueError
        For an item without ``=`` or a non-numeric charge.
    """
    out = {}
    for value in values or ():
        for item in str(value).split(","):
            item = item.strip()
            if not item:
                continue
            name, sep, charge = item.rpartition("=")
            if not sep or not name.strip():
                raise ValueError(f"--charge expects NAME=q, got {item!r}")
            try:
                out[name.strip().upper()] = float(charge)
            except ValueError as exc:
                raise ValueError(f"--charge expects a number after '=', got {item!r}") from exc
    return out


def _topology_charges(atoms):
    try:
        return np.asarray(atoms.charges, dtype=float)
    except (AttributeError, Exception):
        return None


def resolve_species(universe, tokens, *, charges=None, water_mapping=None) -> tuple:
    """Resolve species tokens into :class:`Species` records.

    A plain token matches residue names first and, if none match, atom names
    of non-protein atoms. Each residue is one molecule; its representative
    atom is the only selected atom, else the first named ``OH2``, ``OW``,
    ``O``, ``OW1``, ``OH`` or ``W``, else the first selected atom. A molecule
    already claimed by an earlier token is not counted twice. Martini ``ION``
    residues are split by atom name (``ION:NA+``, ``ION:CL-``).

    Charges: ``charges`` overrides (by label, residue name or atom name), then
    :data:`DEFAULT_CHARGES` (by residue name, or atom name for Martini ions),
    then the topology charge of the molecule (rounded to 0.001 e), else 0.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The system.
    tokens : sequence of str
        From :func:`parse_permeation_tokens`.
    charges : dict, optional
        Overrides from :func:`parse_charge_overrides`.
    water_mapping : int, optional
        Water molecules per Martini water bead (default 4).

    Returns
    -------
    species : list of Species
        Species with at least one molecule.
    absent : list of str
        Tokens that matched nothing.

    Raises
    ------
    ValueError
        If no token matches any atom, or a selection is invalid.
    """
    charges = {k.upper(): v for k, v in (charges or {}).items()}
    atoms = universe.atoms
    names = np.asarray([str(n) for n in atoms.names])
    resnames = np.asarray([str(r) for r in atoms.resnames])
    upper_res = np.char.upper(resnames.astype(str))
    upper_names = np.char.upper(names.astype(str))
    try:
        protein_mask = np.zeros(len(atoms), dtype=bool)
        protein_mask[universe.select_atoms("protein").indices] = True
    except Exception:
        protein_mask = np.zeros(len(atoms), dtype=bool)
    claimed = np.zeros(len(universe.residues), dtype=bool)
    species, absent = [], []
    for token in tokens:
        label, selection = token, None
        if any(c.isspace() for c in token):
            head, sep, tail = token.partition("=")
            if sep and head.strip() and not any(c.isspace() for c in head.strip()):
                label, selection = head.strip(), tail.strip()
            else:
                selection = token
            try:
                group = universe.select_atoms(selection)
            except Exception as exc:
                raise ValueError(f"invalid --permeation selection {selection!r}: {exc}") from exc
            indices = group.indices
            how = "selection"
        else:
            key = token.upper()
            indices = np.flatnonzero(upper_res == key)
            how = "resname"
            if not len(indices):
                indices = np.flatnonzero((upper_names == key) & ~protein_mask)
                how = "atom name"
        if not len(indices):
            absent.append(token)
            continue
        group = atoms[indices]
        resindex = group.resindices
        # Split Martini ions by atom name; otherwise one species per token.
        if how != "selection" and set(np.char.upper(group.resnames.astype(str))) <= MARTINI_ION_RESNAMES:
            parts = [(f"{token}:{n}", np.flatnonzero(np.asarray(group.names) == n)) for n in dict.fromkeys(group.names)]
        else:
            parts = [(label, np.arange(len(group)))]
        for part_label, members in parts:
            sub = group[members]
            reps, mols, rule = [], [], set()
            order = np.argsort(sub.resindices, kind="stable")
            sub = sub[order]
            boundaries = np.flatnonzero(np.diff(sub.resindices)) + 1
            for chunk in np.split(np.arange(len(sub)), boundaries):
                res = int(sub.resindices[chunk[0]])
                if claimed[res]:
                    continue
                if len(chunk) == 1:
                    pick, why = chunk[0], "single atom"
                else:
                    named = [i for i in chunk if str(sub.names[i]).upper() in REPRESENTATIVE_NAMES]
                    pick, why = (named[0], "named oxygen/bead") if named else (chunk[0], "first selected atom")
                claimed[res] = True
                rule.add(why)
                reps.append(int(sub.indices[pick]))
                r = universe.residues[res]
                mols.append({"resid": int(r.resid), "resname": str(r.resname),
                             "segid": str(getattr(r, "segid", "") or "")})
            if not reps:
                continue
            rep_atoms = atoms[np.asarray(reps)]
            resname = str(rep_atoms.resnames[0]).upper()
            atomname = str(rep_atoms.names[0]).upper()
            martini_ion = resname in MARTINI_ION_RESNAMES
            charge, source = 0.0, "none"
            for key in (part_label.upper(), token.upper(), resname, atomname):
                if key in charges:
                    charge, source = charges[key], "override"
                    break
            else:
                key = atomname if martini_ion else resname
                if key in DEFAULT_CHARGES and len(set(np.char.upper(rep_atoms.resnames.astype(str)))) == 1:
                    charge, source = DEFAULT_CHARGES[key], "table"
                else:
                    topo = _topology_charges(universe.residues[rep_atoms.resindices[0]].atoms)
                    if topo is not None and abs(topo.sum()) > 1e-3:
                        charge, source = round(float(topo.sum()), 3), "topology"
            is_water = resname in WATER_RESNAMES and charge == 0
            mapping = 1
            cg = martini_ion
            if is_water and resname in MARTINI_WATER_MAPPING:
                mapping = int(water_mapping or MARTINI_WATER_MAPPING[resname])
                cg = True
            kind = "water" if is_water else ("ion" if charge != 0 else "neutral")
            species.append(Species(part_label, token, kind, np.asarray(reps), float(charge), source,
                                   mapping, cg, "; ".join(sorted(rule)) + f" (matched by {how})", mols))
    if not species:
        raise ValueError("no --permeation species matched any atom: " + ", ".join(tokens))
    return species, absent


# ---------------------------------------------------------------------------
# Applied voltage
# ---------------------------------------------------------------------------

def read_applied_field(path) -> dict:
    """Read a static applied electric field from a GROMACS ``.mdp`` or ``.tpr``.

    ``.mdp``: ``electric-field-x/y/z = E0 omega t0 sigma`` (GROMACS 2018 and
    later) or the older ``E-x = n E0 phase``. ``.tpr``: the
    ``electric-field`` block of ``gmx dump -s`` (needs ``gmx`` or ``gmx_mpi``
    on the PATH).

    Parameters
    ----------
    path : str or Path
        Parameter or run-input file.

    Returns
    -------
    dict
        ``E_V_per_nm`` (x, y, z amplitudes), ``static`` (no oscillation or
        pulse), ``source`` and ``reader``.

    Raises
    ------
    ValueError
        If the file cannot be read or ``gmx`` is not available for a ``.tpr``.
    """
    path = Path(path)
    field_v = [0.0, 0.0, 0.0]
    static = True
    if path.suffix.lower() == ".tpr":
        gmx = shutil.which("gmx") or shutil.which("gmx_mpi") or shutil.which("gmx_d")
        if gmx is None:
            raise ValueError("reading the field from a .tpr needs GROMACS (gmx) on the PATH; give --efield V/nm "
                             "or the .mdp file instead")
        run = subprocess.run([gmx, "-quiet", "dump", "-s", str(path)], capture_output=True, text=True, timeout=300)
        if run.returncode != 0:
            raise ValueError(f"gmx dump failed for {path}: {run.stderr.strip()[-300:]}")
        text = run.stdout
        block = text.partition("electric-field:")[2]
        if not block:
            raise ValueError(f"{path}: no electric-field block in gmx dump output")
        for k, axis in enumerate("xyz"):
            part = re.search(rf"^\s*{axis}:\s*\n((?:\s+\w+\s*=\s*\S+\s*\n)+)", block, flags=re.M)
            if not part:
                continue
            values = dict(re.findall(r"(\w+)\s*=\s*(\S+)", part.group(1)))
            field_v[k] = float(values.get("E0", 0))
            if float(values.get("omega", 0)) != 0 or float(values.get("sigma", 0)) != 0:
                static = False
        reader = "gmx dump"
    else:
        text = path.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            line = line.split(";", 1)[0].strip()
            m = re.match(r"(?i)(electric[-_]field[-_]|E[-_])([xyz])\s*=\s*(.*)", line)
            if not m:
                continue
            k = "xyz".index(m.group(2).lower())
            words = [float(w) for w in m.group(3).split()]
            if m.group(1).lower().startswith("electric"):
                if words:
                    field_v[k] = words[0]
                    if len(words) > 1 and words[1] != 0 or len(words) > 3 and words[3] != 0:
                        static = False
            elif len(words) >= 2 and words[0] >= 1:
                field_v[k] = words[1]
                if words[0] > 1:
                    static = False
        reader = "mdp"
    return {"E_V_per_nm": field_v, "static": static, "source": str(path), "reader": reader}


def voltage_record(*, voltage_mV=None, efield=None, box_lengths_A=None, axis=None, source=None) -> dict:
    """The transmembrane voltage along the channel axis, when one is known.

    Sign convention: ``V`` is positive when it drives cations along ``+axis``
    (the potential on the ``below`` side minus the potential on the
    ``above`` side). A field ``E`` along a box axis gives ``V = E * L`` with
    ``L`` the mean box length along it, signed by the axis direction.

    Parameters
    ----------
    voltage_mV : float, optional
        Given voltage, mV (takes precedence).
    efield : sequence of 3 float, optional
        Applied field (V/nm) along x, y, z; at most one non-zero component.
    box_lengths_A : sequence of 3 float, optional
        Mean box lengths along x, y, z, Å.
    axis : array_like, optional
        Channel axis.
    source : str, optional
        Where the field came from.

    Returns
    -------
    dict
        ``status`` (``given``, ``from_field``, ``undefined``), ``voltage_mV``
        and how it was obtained.
    """
    if voltage_mV is not None:
        return {"status": "given", "voltage_mV": float(voltage_mV),
                "sign_convention": "positive drives cations along +axis (below side minus above side)"}
    if efield is None or not any(efield):
        return {"status": "undefined", "voltage_mV": None,
                "reason": "no voltage given (--voltage) and no applied field (--efield/--efield-from)"
                          + (f"; the field read from {source} is zero" if source else "")}
    e = np.asarray(efield, dtype=float)
    nonzero = np.flatnonzero(e)
    if len(nonzero) != 1:
        return {"status": "undefined", "voltage_mV": None, "E_V_per_nm": e.tolist(),
                "reason": "only a field along one box axis is supported"}
    k = int(nonzero[0])
    if box_lengths_A is None:
        return {"status": "undefined", "voltage_mV": None, "reason": "box length unknown"}
    a = _unit(axis) if axis is not None else np.eye(3)[k]
    cosine = float(a[k])
    length_nm = float(box_lengths_A[k]) / 10
    record = {"E_V_per_nm": e.tolist(), "field_axis": "xyz"[k], "mean_box_length_nm": length_nm,
              "axis_field_cosine": cosine, "source": source,
              "sign_convention": "positive drives cations along +axis (below side minus above side)",
              "formula": "V = E * L_box * sign(axis . field direction)"}
    if abs(cosine) < math.cos(math.radians(45)):
        record.update(status="undefined", voltage_mV=None,
                      reason="the channel axis is more than 45 degrees from the field direction")
        return record
    record.update(status="from_field", voltage_mV=1000 * e[k] * length_nm * math.copysign(1, cosine))
    return record


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def poisson_interval(count, confidence=0.95) -> tuple:
    """Exact (Garwood) two-sided interval for a Poisson count.

    Parameters
    ----------
    count : int
        Observed count.
    confidence : float, default 0.95
        Two-sided level.

    Returns
    -------
    tuple of float
        ``(lower, upper)`` bounds on the expected count.

    Examples
    --------
    >>> [round(v, 3) for v in poisson_interval(0)]
    [0.0, 3.689]
    """
    from scipy.stats import chi2
    alpha = 1 - confidence
    lower = 0.0 if count == 0 else 0.5 * chi2.ppf(alpha / 2, 2 * count)
    upper = 0.5 * chi2.ppf(1 - alpha / 2, 2 * count + 2)
    return float(lower), float(upper)


def collective_diffusion(n, dt_ps, *, blocks=5, max_lag=None, min_lag=1, confidence=0.95) -> dict:
    """Diffusion coefficient of a collective coordinate from its mean-square displacement.

    ``MSD(tau) = <(n(t + tau) - n(t))^2>`` over all time origins is fitted by
    least squares with ``a + 2 D tau`` for lags ``min_lag..max_lag`` frames.
    The interval comes from the same fit on ``blocks`` contiguous segments:
    estimate +/- t(blocks - 1) * SD(segments) / sqrt(blocks).

    Parameters
    ----------
    n : array_like
        Collective coordinate per frame (dimensionless).
    dt_ps : float
        Frame interval, ps.
    blocks : int, default 5
        Number of segments for the interval.
    max_lag : int, optional
        Largest fitted lag in frames; default a tenth of a segment
        (at least 2).
    min_lag : int, default 1
        Smallest fitted lag in frames.
    confidence : float, default 0.95
        Two-sided level.

    Returns
    -------
    dict
        ``status``, ``D_per_ps``, ``D_per_s``, ``interval_per_s``, ``intercept``,
        per-segment estimates, the fitted lags and the MSD curve. The lower
        bound of ``interval_per_s`` is clipped at zero
        (``interval_unclipped_per_s`` keeps the raw value).

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(1)
    >>> n = np.concatenate([[0], np.cumsum(rng.normal(0, np.sqrt(2 * 0.01 * 1.0), 20000))])
    >>> round(collective_diffusion(n, 1.0)["D_per_ps"], 2)
    0.01
    """
    from scipy.stats import t as student
    n = np.asarray(n, dtype=float)
    frames = len(n)
    segment = frames // max(1, int(blocks))
    lag_hi = int(max_lag) if max_lag is not None else max(2, segment // 10)
    lag_lo = int(min_lag)
    record = {"method": "least-squares fit MSD(tau) = a + 2 D tau over all time origins",
              "frames": frames, "frame_interval_ps": float(dt_ps), "fit_lags_frames": [lag_lo, lag_hi],
              "fit_lags_ps": [lag_lo * dt_ps, lag_hi * dt_ps], "blocks": int(blocks),
              "block_frames": segment, "confidence": confidence}
    if frames < 2 * lag_hi + 2 or lag_hi <= lag_lo:
        record.update(status="insufficient_frames", D_per_ps=None, D_per_s=None, interval_per_s=None)
        return record

    def fit(series):
        lags = np.arange(lag_lo, lag_hi + 1)
        msd = np.asarray([np.mean((series[lag:] - series[:-lag]) ** 2) for lag in lags])
        tau = lags * dt_ps
        slope, intercept = np.polyfit(tau, msd, 1)
        return slope / 2, intercept, tau, msd

    d, intercept, tau, msd = fit(n)
    record.update(status="estimated", D_per_ps=float(d), D_per_s=float(d * 1e12), intercept=float(intercept),
                  msd_curve={"lag_ps": tau.tolist(), "msd": msd.tolist()})
    estimates = []
    if blocks >= 2 and segment >= 2 * lag_hi + 2:
        for b in range(int(blocks)):
            estimates.append(fit(n[b * segment:(b + 1) * segment])[0])
    record["block_D_per_ps"] = [float(v) for v in estimates]
    if len(estimates) >= 2:
        se = float(np.std(estimates, ddof=1) / math.sqrt(len(estimates)))
        crit = float(student.ppf(0.5 + confidence / 2, len(estimates) - 1))
        record["standard_error_per_s"] = se * 1e12
        raw = [float((d - crit * se) * 1e12), float((d + crit * se) * 1e12)]
        record["interval_unclipped_per_s"] = raw
        # D_n cannot be negative: the reported lower bound is clipped at zero.
        record["interval_per_s"] = [max(0.0, raw[0]), raw[1]]
        record["interval_clipped_at_zero"] = bool(raw[0] < 0)
        record["interval_method"] = f"estimate +/- t({len(estimates) - 1}) x SD of {len(estimates)} segment estimates / sqrt({len(estimates)})"
    else:
        record["interval_per_s"] = None
        record["interval_method"] = "unavailable: segments shorter than twice the largest fitted lag"
    return record


def _block_bootstrap(per_frame, *, block_frames, replicates=2000, seed=20261020, confidence=0.95):
    """Percentile interval of the per-frame mean from resampled contiguous blocks."""
    x = np.asarray(per_frame, dtype=float)
    k = len(x) // block_frames
    if k < 2:
        return {"status": "insufficient_blocks", "blocks": k}
    sums = x[:k * block_frames].reshape(k, block_frames).mean(axis=1)
    if np.ptp(sums) == 0:
        return {"status": "no_variation_between_blocks", "blocks": k, "block_frames": block_frames,
                "reason": "every block has the same value (for example no events); a bootstrap interval is not informative"}
    rng = np.random.default_rng(seed)
    boot = sums[rng.integers(0, k, size=(int(replicates), k))].mean(axis=1)
    alpha = (1 - confidence) / 2
    return {"status": "estimated", "blocks": k, "block_frames": block_frames, "replicates": int(replicates),
            "seed": seed, "lower_per_frame": float(np.quantile(boot, alpha)),
            "upper_per_frame": float(np.quantile(boot, 1 - alpha)),
            "method": "percentile interval of the mean over resampled contiguous blocks (frames beyond the last whole block dropped)"}


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

@dataclass
class PermeationResult:
    """Result of :func:`analyze_permeation`.

    Attributes
    ----------
    summary : dict
        Definitions, settings, per-channel and per-species statistics (the
        JSON summary).
    events : list of dict
        One row per permeation event or unresolved crossing.
    visits : list of dict
        One row per lumen visit (contiguous inside frames).
    frames : dict
        Per-frame arrays: ``frame``, ``time_ps`` and, per ``"channel/species"``
        key, ``occupancy``, ``net_crossings`` (cumulative) and ``n`` (the
        collective coordinate).
    """

    summary: dict
    events: list
    visits: list
    frames: dict


def _open(topology, trajectory):
    if hasattr(topology, "trajectory") and hasattr(topology, "atoms"):
        return topology
    from .trajectory import _open_universe
    return _open_universe(topology, trajectory)


def _kabsch(mobile, reference):
    xc, yc = mobile.mean(axis=0), reference.mean(axis=0)
    h = (mobile - xc).T @ (reference - yc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt)) or 1.0
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    rmsd = float(np.sqrt(np.mean(np.sum(((mobile - xc) @ rotation + yc - reference) ** 2, axis=1))))
    return rotation, xc, yc, rmsd


def analyze_permeation(topology, trajectory=None, *, species=("W", "NA", "CL"), channels=None, start=0, stop=None,
                       stride=1, align=True, alignment_selection="protein and name CA", charges=None,
                       water_mapping=None, voltage_mV=None, efield=None, efield_source=None, msd_blocks=5,
                       msd_max_lag=None, bootstrap_replicates=2000, bootstrap_seed=20261020, confidence=0.95,
                       progress_callback=None) -> PermeationResult:
    """Count permeation events and derive water permeabilities and ion current.

    Parameters
    ----------
    topology : str, Path or MDAnalysis.Universe
        Topology file, or an open Universe (then ``trajectory`` is ignored).
    trajectory : str or Path, optional
        Trajectory file.
    species : sequence of str, default ("W", "NA", "CL")
        Species tokens (residue names, atom names or MDAnalysis selections;
        :func:`parse_permeation_tokens`). Tokens that match nothing are listed
        as absent.
    channels : sequence
        :class:`CylinderChannel`, :class:`ProfileChannel` or
        :class:`RegionChannel` objects (reference-frame coordinates). Several
        channels are analysed independently (for example one per pore).
    start, stop, stride : int
        Frames analysed: ``range(start, stop, stride)``. The first analysed
        frame is the reference frame.
    align : bool, default True
        Fit the alignment atoms of every frame onto the reference frame;
        ``False`` uses the raw coordinates.
    alignment_selection : str, default "protein and name CA"
        MDAnalysis selection of the fit atoms (and of the default axis point).
    charges : dict, optional
        Charge overrides (:func:`parse_charge_overrides`).
    water_mapping : int, optional
        Water molecules per Martini water bead (default 4).
    voltage_mV : float, optional
        Known transmembrane voltage (see :func:`voltage_record` for the sign).
    efield : sequence of 3 float, optional
        Applied static field, V/nm; with the mean box length it gives the voltage.
    efield_source : str, optional
        Where ``efield`` came from (recorded).
    msd_blocks : int, default 5
        Segments for the D_n interval (:func:`collective_diffusion`).
    msd_max_lag : int, optional
        Largest MSD lag in frames.
    bootstrap_replicates : int, default 2000
        Block-bootstrap replicates of the current and occupancy intervals.
    bootstrap_seed : int, default 20261020
        Seed of the block bootstrap.
    confidence : float, default 0.95
        Two-sided level of every interval.
    progress_callback : callable, optional
        Called as ``(index, total)`` after each frame.

    Returns
    -------
    PermeationResult

    Raises
    ------
    ValueError
        Without channels, for an empty frame range, or when no species matches.
    """
    if not channels:
        raise ValueError("analyze_permeation needs at least one channel definition")
    names = [c.name for c in channels]
    if len(set(names)) != len(names):
        raise ValueError("channel names must be unique")
    universe = _open(topology, trajectory)
    tokens = parse_permeation_tokens(species)
    found, absent = resolve_species(universe, tokens, charges=charges, water_mapping=water_mapping)
    total_frames = len(universe.trajectory)
    from .trajectory import _frame_range
    frame_numbers = list(_frame_range(total_frames, start, stop, stride))
    if len(frame_numbers) < 2:
        raise ValueError("permeation analysis needs at least two frames")
    fit_atoms = universe.select_atoms(alignment_selection) if alignment_selection else universe.atoms[:0]
    if align and len(fit_atoms) < 3:
        raise ValueError(f"alignment selection {alignment_selection!r} has fewer than 3 atoms; use --no-align or "
                         "another selection (Martini backbone beads: 'name BB')")
    coarse = any(s.coarse_grained for s in found) or bool(len(universe.select_atoms("name BB")))
    all_indices = np.concatenate([s.atom_indices for s in found])
    offsets = np.cumsum([0] + [len(s.atom_indices) for s in found])
    rep = universe.atoms[all_indices]
    n_mol = len(rep)

    universe.trajectory[frame_numbers[0]]
    reference_fit = fit_atoms.positions.astype(float).copy() if len(fit_atoms) else None
    default_point = reference_fit.mean(axis=0) if reference_fit is not None else np.zeros(3)
    for channel in channels:
        channel.bind(default_point)
        if hasattr(channel, "bind_frames"):
            channel.bind_frames(frame_numbers)

    F, C, S = len(frame_numbers), len(channels), len(found)
    times = np.zeros(F)
    occupancy = np.zeros((C, S, F), dtype=np.int32)
    n_coord = np.zeros((C, S, F))
    net_cum = np.zeros((C, S, F), dtype=np.int64)
    box_lengths = []
    fit_rmsd = []
    hist_bins = []
    hist = []
    for channel in channels:
        lo, hi = channel.fixed_extent()
        edges = np.arange(math.floor(lo), math.ceil(hi) + 0.5, 0.5)
        hist_bins.append(edges)
        hist.append(np.zeros((S, len(edges) - 1)))
    species_of = np.repeat(np.arange(S), np.diff(offsets))
    # Per-channel online state.
    st = [{"prev_state": None, "prev_d": None, "prev_s": None,
           "side": np.zeros(n_mol, dtype=np.int8), "inside": np.zeros(n_mol, dtype=bool),
           "visit_start": np.full(n_mol, -1, dtype=np.int64), "entry_side": np.zeros(n_mol, dtype=np.int8),
           "net": np.zeros(S, dtype=np.int64), "n": np.zeros(S)} for _ in channels]
    events, visits = [], []

    def close_visits(c, mask, frame_index, exit_to):
        s_state = st[c]
        for i in np.flatnonzero(mask):
            begin = int(s_state["visit_start"][i])
            entered = int(s_state["entry_side"][i])
            visits.append({"channel": c, "molecule": int(i), "start": begin, "end": frame_index - 1,
                           "entered_from": entered, "exited_to": exit_to if isinstance(exit_to, str) else int(exit_to[i])})

    for k, frame in enumerate(frame_numbers):
        ts = universe.trajectory[frame]
        times[k] = float(ts.time)
        box = ts.dimensions
        has_box = box is not None and np.all(np.asarray(box[:3]) > 0)
        if has_box:
            box_lengths.append(np.asarray(box[:3], dtype=float))
        x = rep.positions.astype(float)
        if align:
            rotation, xc, yc, rmsd = _kabsch(fit_atoms.positions.astype(float), reference_fit)
            fit_rmsd.append(rmsd)
        else:
            rotation, xc, yc = np.eye(3), np.zeros(3), np.zeros(3)
        for c, channel in enumerate(channels):
            s_state = st[c]
            axis = np.asarray(channel.axis, dtype=float)
            centre = channel.centre()
            centre_lab = (centre - yc) @ rotation.T + xc
            d_lab = x - centre_lab
            if has_box:
                from MDAnalysis.lib.distances import minimize_vectors
                d_lab = minimize_vectors(d_lab, np.asarray(box, dtype=np.float64))
            d_ref = d_lab @ rotation
            pos = centre + d_ref
            s = pos @ axis
            zmin, zmax = channel.extent(k)
            state = np.where(s < zmin, _STATE_BELOW, np.where(s > zmax, _STATE_ABOVE, _STATE_OUT)).astype(np.int8)
            slab = np.flatnonzero(state == _STATE_OUT)
            if len(slab):
                inside = channel.lumen(pos[slab], k)
                state[slab[inside]] = _STATE_INSIDE
            is_in = state == _STATE_INSIDE
            if s_state["prev_state"] is not None:
                prev = s_state["prev_state"]
                delta = d_lab - s_state["prev_d"]
                if has_box:
                    delta = minimize_vectors(delta, np.asarray(box, dtype=np.float64))
                ds = (delta @ rotation) @ axis
                s_prev = s_state["prev_s"]
                s_new = s_prev + ds
                # n(t): only the part of each displacement inside the fixed extent counts.
                lo, hi = channel.fixed_extent()
                through = ((prev == _STATE_BELOW) & (s_new > zmax)) | ((prev == _STATE_ABOVE) & (s_new < zmin))
                moving = (prev == _STATE_INSIDE) | is_in | through
                if moving.any():
                    contrib = np.clip(s_new[moving], lo, hi) - np.clip(s_prev[moving], lo, hi)
                    s_state["n"] += np.bincount(species_of[moving], weights=contrib, minlength=S) / (hi - lo)
                # Unresolved crossings: the whole extent crossed between two saved frames.
                for i in np.flatnonzero(through & ~(prev == _STATE_INSIDE)):
                    events.append({"channel": c, "molecule": int(i), "direction": "up" if prev[i] == _STATE_BELOW else "down",
                                   "entry_frame": k - 1, "exit_frame": k, "status": "unresolved_single_interval"})
                # Leaving the lumen.
                leaving = s_state["inside"] & ~is_in
                to_side = leaving & (state != _STATE_OUT)
                if to_side.any():
                    close_visits(c, to_side, k, state)
                    for i in np.flatnonzero(to_side):
                        entered = int(s_state["entry_side"][i])
                        if entered != 0 and entered == -int(state[i]):
                            direction = "up" if entered == _STATE_BELOW else "down"
                            events.append({"channel": c, "molecule": int(i), "direction": direction,
                                           "entry_frame": int(s_state["visit_start"][i]), "exit_frame": k,
                                           "status": "complete"})
                            s_state["net"][species_of[i]] += 1 if direction == "up" else -1
                lateral = leaving & (state == _STATE_OUT)
                if lateral.any():
                    close_visits(c, lateral, k, "out")
                # Entering the lumen.
                entering = is_in & ~s_state["inside"]
                s_state["visit_start"][entering] = k
                s_state["entry_side"][entering] = s_state["side"][entering]
            else:
                s_state["visit_start"][is_in] = 0
                s_state["entry_side"][is_in] = 0
            s_state["inside"] = is_in
            bulk = (state == _STATE_BELOW) | (state == _STATE_ABOVE)
            s_state["side"][bulk] = state[bulk]
            s_state["side"][state == _STATE_OUT] = 0
            s_state["prev_state"], s_state["prev_d"], s_state["prev_s"] = state, d_lab, s
            occupancy[c, :, k] = np.bincount(species_of[is_in], minlength=S)
            n_coord[c, :, k] = s_state["n"]
            net_cum[c, :, k] = s_state["net"]
            if is_in.any():
                for j in range(S):
                    sel = is_in & (species_of == j)
                    if sel.any():
                        hist[c][j] += np.histogram(s[sel], bins=hist_bins[c])[0]
        if progress_callback is not None:
            progress_callback(k, F)
    for c in range(C):
        open_mask = st[c]["inside"]
        for i in np.flatnonzero(open_mask):
            visits.append({"channel": c, "molecule": int(i), "start": int(st[c]["visit_start"][i]), "end": F - 1,
                           "entered_from": int(st[c]["entry_side"][i]), "exited_to": "open_at_end"})

    return _assemble(universe, found, absent, channels, frame_numbers, times, occupancy, n_coord, net_cum, events,
                     visits, offsets, species_of, box_lengths, fit_rmsd, hist_bins, hist, coarse,
                     dict(align=align, alignment_selection=alignment_selection, start=start, stop=stop, stride=stride,
                          voltage_mV=voltage_mV, efield=efield, efield_source=efield_source, msd_blocks=msd_blocks,
                          msd_max_lag=msd_max_lag, replicates=bootstrap_replicates, seed=bootstrap_seed,
                          confidence=confidence, total_frames=total_frames))


def _side_name(value):
    return {-1: "below", 1: "above", 0: "unknown"}.get(value, value)


def _assemble(universe, found, absent, channels, frame_numbers, times, occupancy, n_coord, net_cum, events, visits,
              offsets, species_of, box_lengths, fit_rmsd, hist_bins, hist, coarse, opts) -> PermeationResult:
    from .cavity_trajectory import describe_series
    F = len(frame_numbers)
    diffs = np.diff(times)
    dt = float(np.median(diffs)) if len(diffs) else 0.0
    regular = bool(len(diffs) and dt > 0 and np.allclose(diffs, dt, rtol=1e-3, atol=1e-3))
    duration_ps = float(times[-1] - times[0]) if F > 1 else 0.0
    if not duration_ps > 0:
        raise ValueError("frame times must increase; the trajectory has no usable time axis")
    duration_s = duration_ps * 1e-12
    conf = opts["confidence"]
    mean_box = np.mean(box_lengths, axis=0).tolist() if box_lengths else None
    block_frames = max(1, F // max(2, opts["msd_blocks"] * 2))

    def mol_row(i):
        j = int(species_of[i])
        local = i - offsets[j]
        info = found[j].molecules[local]
        return j, {"species": found[j].label, "atom_index": int(found[j].atom_indices[local]), **info}

    event_rows = []
    for e in events:
        j, info = mol_row(e["molecule"])
        entry, exit_ = e["entry_frame"], e["exit_frame"]
        event_rows.append({"channel": channels[e["channel"]].name, **info, "direction": e["direction"],
                           "status": e["status"], "entry_frame": frame_numbers[entry], "exit_frame": frame_numbers[exit_],
                           "entry_time_ps": float(times[entry]), "exit_time_ps": float(times[exit_]),
                           "transit_time_ps": float(times[exit_] - times[entry]), "_species_index": j,
                           "_channel_index": e["channel"], "_exit_index": exit_})
    visit_rows = []
    for v in visits:
        j, info = mol_row(v["molecule"])
        b, end = v["start"], v["end"]
        visit_rows.append({"channel": channels[v["channel"]].name, **info, "start_frame": frame_numbers[b],
                           "end_frame": frame_numbers[end], "frames_inside": end - b + 1,
                           "dwell_time_ps": float((end - b + 1) * dt), "entered_from": _side_name(v["entered_from"]),
                           "exited_to": _side_name(v["exited_to"]),
                           "censored": v["entered_from"] == 0 and b == 0 or v["exited_to"] == "open_at_end",
                           "_species_index": j, "_channel_index": v["channel"]})

    voltage = voltage_record(voltage_mV=opts["voltage_mV"], efield=opts["efield"], box_lengths_A=mean_box,
                             axis=channels[0].axis, source=opts["efield_source"])

    def stats_for(c, j, sp):
        complete = [e for e in event_rows if e["_channel_index"] == c and e["_species_index"] == j and e["status"] == "complete"]
        unresolved = sum(1 for e in event_rows if e["_channel_index"] == c and e["_species_index"] == j and e["status"] != "complete")
        up = sum(e["direction"] == "up" for e in complete)
        down = len(complete) - up
        lo_up, hi_up = poisson_interval(up, conf)
        lo_dn, hi_dn = poisson_interval(down, conf)
        transit = np.asarray([e["transit_time_ps"] for e in complete], float)
        mine = [v for v in visit_rows if v["_channel_index"] == c and v["_species_index"] == j]
        dwell = np.asarray([v["dwell_time_ps"] for v in mine if not v["censored"]], float)
        occ = describe_series(occupancy[c, j].astype(float), confidence=conf, replicates=opts["replicates"],
                              seed=opts["seed"], regular=regular)
        per_ns = 1e3 / duration_ps
        record = {"species": sp.label, "kind": sp.kind,
                  "events_up": up, "events_down": down, "events_total": up + down, "net_up_minus_down": up - down,
                  "unresolved_crossings": unresolved,
                  "rate_up_per_ns": up * per_ns, "rate_up_interval_per_ns": [lo_up * per_ns, hi_up * per_ns],
                  "rate_down_per_ns": down * per_ns, "rate_down_interval_per_ns": [lo_dn * per_ns, hi_dn * per_ns],
                  "net_flux_per_ns": (up - down) * per_ns,
                  "transit_time_ps": ({"mean": float(transit.mean()), "median": float(np.median(transit)),
                                       "min": float(transit.min()), "max": float(transit.max())} if len(transit) else None),
                  "lumen_visits": len(mine), "uncensored_visits": int(len(dwell)),
                  "dwell_time_ps": ({"mean": float(dwell.mean()), "median": float(np.median(dwell)),
                                     "max": float(dwell.max())} if len(dwell) else None),
                  "occupancy": {k: occ.get(k) for k in ("mean", "sample_sd", "median", "minimum", "maximum",
                                                        "statistical_inefficiency", "effective_frames")}
                               | {"confidence_interval": _compact_interval(occ.get("confidence_interval"))},
                  "axial_density": {"bin_edges_A": hist_bins[c].tolist(),
                                    "mean_molecules_per_bin": (hist[c][j] / F).tolist()}}
        m = sp.water_mapping
        if sp.kind == "water":
            q0_count = (up + down) / 2
            lo_t, hi_t = poisson_interval(up + down, conf)
            q0 = m * q0_count / duration_s
            pd = WATER_MOLECULE_VOLUME_CM3 * q0
            dn = collective_diffusion(m * n_coord[c, j], dt, blocks=opts["msd_blocks"], max_lag=opts["msd_max_lag"],
                                      confidence=conf)
            pf = WATER_MOLECULE_VOLUME_CM3 * dn["D_per_s"] if dn.get("D_per_s") is not None else None
            pf_int = ([WATER_MOLECULE_VOLUME_CM3 * v for v in dn["interval_per_s"]]
                      if dn.get("interval_per_s") else None)
            record["water"] = {
                "water_molecules_per_representative": m,
                "q0_per_ns": q0 * 1e-9, "q0_interval_per_ns": [m * lo_t / 2 / duration_s * 1e-9, m * hi_t / 2 / duration_s * 1e-9],
                "pd_cm3_per_s": pd,
                "pd_interval_cm3_per_s": [WATER_MOLECULE_VOLUME_CM3 * m * lo_t / 2 / duration_s,
                                          WATER_MOLECULE_VOLUME_CM3 * m * hi_t / 2 / duration_s],
                "pd_interval_method": "exact Poisson (Garwood) interval on the total event count, halved",
                "collective_diffusion": dn,
                "D_n_per_s": dn.get("D_per_s"),
                "pf_cm3_per_s": pf, "pf_interval_cm3_per_s": pf_int,
                "pf_over_pd": (pf / pd) if pf is not None and pd > 0 else None,
                "v_w_cm3": WATER_MOLECULE_VOLUME_CM3,
                "definitions": {"q0": "unidirectional permeation rate = (events up + events down) / (2 T)",
                                "pd": "v_w q0 (Zhu, Tajkhorshid & Schulten, Biophys. J. 2004, Eq. 3)",
                                "n": "collective coordinate: dn = sum of axial displacements inside the extent / L",
                                "D_n": "MSD(n)/(2 tau) at equilibrium",
                                "pf": "v_w D_n (collective-diffusion model, Zhu, Tajkhorshid & Schulten 2004)",
                                "v_w": "18 cm^3/mol / N_A"},
                "coarse_grained_mapping": (f"each bead is {m} water molecules: event counts, q0 and pd are multiplied by "
                                           f"{m}; n is expressed in water molecules ({m} x bead n), so D_n is "
                                           f"multiplied by {m * m}" if m != 1 else None)}
        if sp.charge != 0:
            record["ion"] = _current(c, j, sp, up, down, event_rows, F, block_frames, duration_s, opts)
        return record

    channel_records = []
    for c, channel in enumerate(channels):
        rows = [stats_for(c, j, sp) for j, sp in enumerate(found)]
        lo, hi = channel.fixed_extent()
        total_current = None
        charged = [r for r in rows if "ion" in r]
        if charged:
            total_current = _total_current(c, found, event_rows, F, block_frames, duration_s, opts, voltage)
        channel_records.append({"definition": channel.describe(), "collective_extent_A": [lo, hi],
                                "species": rows, "ion_current_total": total_current})
    totals = []
    for j, sp in enumerate(found):
        rows = [cr["species"][j] for cr in channel_records]
        entry = {"species": sp.label, "events_up": sum(r["events_up"] for r in rows),
                 "events_down": sum(r["events_down"] for r in rows),
                 "unresolved_crossings": sum(r["unresolved_crossings"] for r in rows)}
        entry["net_up_minus_down"] = entry["events_up"] - entry["events_down"]
        if sp.kind == "water":
            pfs = [r["water"]["pf_cm3_per_s"] for r in rows]
            entry["pd_cm3_per_s_sum_over_channels"] = sum(r["water"]["pd_cm3_per_s"] for r in rows)
            entry["pf_cm3_per_s_sum_over_channels"] = sum(pfs) if all(v is not None for v in pfs) else None
        totals.append(entry)
    for record in channel_records:
        if record["ion_current_total"] is not None:
            record["ion_current_total"]["voltage"] = voltage
    rmsd = np.asarray(fit_rmsd) if fit_rmsd else None
    warnings = []
    if absent:
        warnings.append("species tokens that matched no atom: " + ", ".join(absent))
    if not regular:
        warnings.append("frame times are irregular: occupancy intervals withheld; D_n uses the median interval")
    if rmsd is not None and rmsd.max() > 3.0:
        warnings.append(f"alignment RMSD up to {rmsd.max():.2f} A: check that the protein is whole and the selection is stable")
    unresolved_total = sum(t["unresolved_crossings"] for t in totals)
    if unresolved_total:
        warnings.append(f"{unresolved_total} crossings completed between two saved frames (not counted); save frames more often")
    for channel in channels:
        if getattr(channel, "kind", "") == "cylinder" and channel.radius is None:
            warnings.append(f"channel {channel.name}: no lateral radius; every molecule within the axial extent counts as inside")
        if getattr(channel, "kind", "") == "profile" and channel.held_frames:
            warnings.append(f"channel {channel.name}: {channel.held_frames} analysed frames used an earlier profile")
    summary = {
        "schema": "crevice.permeation/1",
        "coarse_grained": bool(coarse),
        "frames": {"analysed": F, "available": opts["total_frames"], "first": frame_numbers[0], "last": frame_numbers[-1],
                   "stride": opts["stride"], "frame_interval_ps": dt, "regular_interval": regular,
                   "duration_ps": duration_ps, "first_time_ps": float(times[0])},
        "reference_frame": {"frame": frame_numbers[0], "aligned": bool(opts["align"]),
                            "alignment_selection": opts["alignment_selection"] if opts["align"] else None,
                            "fit_rmsd_A": ({"mean": float(rmsd.mean()), "max": float(rmsd.max())} if rmsd is not None else None)},
        "mean_box_A": mean_box,
        "species": [s.describe() for s in found],
        "absent_tokens": absent,
        "channels": channel_records,
        "totals_over_channels": totals,
        "voltage": voltage,
        "conductance_defined": voltage["status"] != "undefined",
        "definitions": {
            "states": "below: x.axis < zmin; above: x.axis > zmax; inside: within the extent and the lumen; "
                      "out: within the extent but outside the lumen (clears the remembered side)",
            "event": "complete crossing below -> inside (uninterrupted) -> above ('up') or the reverse ('down'); "
                     "Zhu, Tajkhorshid & Schulten, Biophys. J. 86, 50-57 (2004)",
            "unresolved_crossing": "the whole extent crossed between two consecutive analysed frames (not counted)",
            "periodic_boundaries": "minimum image about the channel centre; frame-to-frame displacement unwrapped per molecule",
            "entry_frame": "first frame inside the lumen", "exit_frame": "first frame on the far side",
            "transit_time": "exit time minus entry time (resolution: one frame interval)",
            "dwell_time": "frames inside x frame interval for one lumen visit; visits cut by the start or end of the trajectory are censored",
            "positive_direction": "+axis ('up', from below to above)",
            "current": "I = e sum_species q (up - down) / T, positive for net positive charge moving along +axis",
            "conductance": "G = I / V, only when the voltage is known",
        },
        "units": {"positions": "A", "time": "ps", "rates": "1/ns", "permeability": "cm^3/s", "D_n": "1/s",
                  "current": "pA", "conductance": "pS", "voltage": "mV", "charge": "e"},
        "confidence": conf,
        "warnings": warnings,
        "limitations": LIMITATIONS,
    }
    frames = {"frame": list(frame_numbers), "time_ps": times.tolist(), "series": {}}
    for c, channel in enumerate(channels):
        for j, sp in enumerate(found):
            frames["series"][f"{channel.name}/{sp.label}"] = {
                "occupancy": occupancy[c, j].tolist(), "net_crossings": net_cum[c, j].tolist(),
                "n": (sp.water_mapping * n_coord[c, j]).tolist() if sp.kind == "water" else n_coord[c, j].tolist()}
    for row in event_rows + visit_rows:
        for key in [k for k in row if k.startswith("_")]:
            del row[key]
    return PermeationResult(summary, event_rows, visit_rows, frames)


def _compact_interval(ci):
    if not isinstance(ci, dict):
        return ci
    if ci.get("status") != "estimated":
        return {"status": ci.get("status"), "reason": ci.get("reason")}
    lower, upper = ci.get("pointwise_lower") or [None], ci.get("pointwise_upper") or [None]
    return {"status": "estimated", "lower": lower[0], "upper": upper[0], "method": ci.get("method"),
            "block_length_frames": ci.get("block_length_frames")}


def _charge_series(c, j_list, charges, event_rows, F):
    series = np.zeros(F)
    for e in event_rows:
        if e["_channel_index"] == c and e["status"] == "complete" and e["_species_index"] in j_list:
            q = charges[e["_species_index"]]
            series[e["_exit_index"]] += q if e["direction"] == "up" else -q
    return series


def _current_record(series, var_counts, max_q, total_events, F, block_frames, duration_s, opts):
    conf = opts["confidence"]
    from scipy.stats import norm
    to_pA = ELEMENTARY_CHARGE_C / duration_s * 1e12
    net_charge = float(series.sum())
    current = net_charge * to_pA
    if total_events:
        half = float(norm.ppf(0.5 + conf / 2) * math.sqrt(var_counts)) * to_pA
        poisson = {"interval_pA": [current - half, current + half],
                   "method": "normal approximation to the difference of Poisson counts (variance sum q^2 (up + down))"}
    else:
        upper = poisson_interval(0, conf)[1] * max_q * to_pA
        poisson = {"interval_pA": [-upper, upper],
                   "method": "no events: exact Poisson upper bound on the total crossing count times the largest |q|"}
    boot = _block_bootstrap(series, block_frames=block_frames, replicates=opts["replicates"], seed=opts["seed"],
                            confidence=conf)
    if boot["status"] == "estimated":
        # current = (mean net charge per frame) x frames x e / T
        boot["interval_pA"] = [boot["lower_per_frame"] * F * to_pA, boot["upper_per_frame"] * F * to_pA]
    return {"net_charge_transferred_e": net_charge, "current_pA": current, "poisson": poisson, "block_bootstrap": boot}


def _current(c, j, sp, up, down, event_rows, F, block_frames, duration_s, opts):
    charges = {j: sp.charge}
    series = _charge_series(c, [j], charges, event_rows, F)
    record = _current_record(series, sp.charge ** 2 * (up + down), abs(sp.charge), up + down, F, block_frames,
                             duration_s, opts)
    record["charge_e"] = sp.charge
    record["charge_source"] = sp.charge_source
    return record


def _total_current(c, found, event_rows, F, block_frames, duration_s, opts, voltage):
    charges = {j: s.charge for j, s in enumerate(found) if s.charge != 0}
    series = _charge_series(c, list(charges), charges, event_rows, F)
    complete = [e for e in event_rows if e["_channel_index"] == c and e["status"] == "complete" and e["_species_index"] in charges]
    var = sum(charges[e["_species_index"]] ** 2 for e in complete)
    record = _current_record(series, var, max(abs(q) for q in charges.values()), len(complete), F, block_frames,
                             duration_s, opts)
    v = voltage.get("voltage_mV")
    if voltage["status"] != "undefined" and v:
        g = record["current_pA"] / (v / 1000)
        record["conductance_pS"] = g
        record["conductance_poisson_interval_pS"] = sorted(x / (v / 1000) for x in record["poisson"]["interval_pA"])
        if "interval_pA" in record["block_bootstrap"]:
            record["conductance_bootstrap_interval_pS"] = sorted(x / (v / 1000) for x in record["block_bootstrap"]["interval_pA"])
        record["conductance_status"] = "defined"
    else:
        record["conductance_pS"] = None
        record["conductance_status"] = "undefined: " + voltage.get("reason", "no voltage")
    return record


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

EVENT_COLUMNS = ("channel", "species", "resid", "resname", "segid", "atom_index", "direction", "status", "entry_frame",
                 "exit_frame", "entry_time_ps", "exit_time_ps", "transit_time_ps")
VISIT_COLUMNS = ("channel", "species", "resid", "resname", "segid", "atom_index", "start_frame", "end_frame",
                 "frames_inside", "dwell_time_ps", "entered_from", "exited_to", "censored")


def _write_csv(path, columns, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_permeation_outputs(result: PermeationResult, output_dir, prefix="crevice", *, dpi=200, annotate=None) -> dict:
    """Write the permeation JSON summary, CSV tables and figure.

    Parameters
    ----------
    result : PermeationResult
        From :func:`analyze_permeation`.
    output_dir : str or Path
        Directory (created if needed).
    prefix : str, default "crevice"
        File-name prefix.
    dpi : int, default 200
        Figure resolution.
    annotate : bool or None, default None
        Draw explanatory text in the figure (``None`` inherits ``--annotate``).

    Returns
    -------
    dict
        Paths of the files written (``summary_json``, ``events_csv``,
        ``visits_csv``, ``frames_csv`` and, if written, ``figure_png``).
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    files = {"summary_json": out / f"{prefix}_permeation.json", "events_csv": out / f"{prefix}_permeation_events.csv",
             "visits_csv": out / f"{prefix}_permeation_visits.csv", "frames_csv": out / f"{prefix}_permeation_frames.csv"}
    _write_csv(files["events_csv"], EVENT_COLUMNS, result.events)
    _write_csv(files["visits_csv"], VISIT_COLUMNS, result.visits)
    keys = list(result.frames["series"])
    columns = ["frame", "time_ps"] + [f"{k}:{q}" for k in keys for q in ("occupancy", "net_crossings", "n")]
    with files["frames_csv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for i, (frame, time) in enumerate(zip(result.frames["frame"], result.frames["time_ps"])):
            row = [frame, f"{time:.6g}"]
            for k in keys:
                s = result.frames["series"][k]
                row += [s["occupancy"][i], s["net_crossings"][i], f"{s['n'][i]:.6g}"]
            writer.writerow(row)
    figure = out / f"{prefix}_permeation.png"
    from .outputs import wants
    if wants(figure):
        plot_permeation(result, figure, dpi=dpi, annotate=annotate)
        files["figure_png"] = figure
    summary = dict(result.summary)
    summary["files"] = {k: v.name for k, v in files.items()}
    from ._portable_paths import dumps_record  # 0.4: recorded paths relative to this file
    files["summary_json"].write_text(dumps_record(summary, files["summary_json"], indent=2, sort_keys=True, allow_nan=False) + "\n",
                                     encoding="utf-8", newline="\n")
    return {k: str(v) for k, v in files.items()}


@annotate_option
def plot_permeation(result: PermeationResult, path, *, dpi=200, annotate=None):
    """Figure: lumen occupancy, cumulative net crossings and the water collective coordinate.

    Three panels share the time axis (ns): the number of molecules inside the
    lumen per species, the cumulative net crossings (up minus down) per
    species, and n(t) of the water species (in water molecules). Series are
    summed over channels. Clean by default; ``annotate=True`` (or
    ``--annotate``) adds a title and the event counts. The description is
    always embedded in the PNG metadata.

    Parameters
    ----------
    result : PermeationResult
        From :func:`analyze_permeation`.
    path : str or Path
        Output PNG.
    dpi : int, default 200
        Resolution.
    annotate : bool or None, default None
        Draw explanatory text (``None`` inherits the surrounding setting).

    Returns
    -------
    matplotlib.figure.Figure
        The figure (already saved and closed).
    """
    from .figures import CREVICE_BLUE, CREVICE_GRAY, CREVICE_GREEN, CREVICE_ORANGE, CREVICE_PURPLE, CREVICE_RED, \
        CREVICE_TEAL, _pyplot, _save_figure, _style_axes
    plt = _pyplot()
    colours = [CREVICE_TEAL, CREVICE_ORANGE, CREVICE_PURPLE, CREVICE_BLUE, CREVICE_RED, CREVICE_GREEN, CREVICE_GRAY]
    time_ns = np.asarray(result.frames["time_ps"]) / 1000
    labels = [s["label"] for s in result.summary["species"]]
    kinds = {s["label"]: s["kind"] for s in result.summary["species"]}
    channels = [c["definition"]["name"] for c in result.summary["channels"]]

    def summed(label, quantity):
        return np.sum([result.frames["series"][f"{c}/{label}"][quantity] for c in channels], axis=0)

    water = [lab for lab in labels if kinds[lab] == "water"]
    rows = 3 if water else 2
    fig, axes = plt.subplots(rows, 1, figsize=(7.4, 2.1 * rows + 0.6), sharex=True, layout="constrained")
    text = FigureText(f"Permeation through {len(channels)} channel(s) ({', '.join(channels)}): molecules inside the "
                      "lumen, cumulative net crossings (up minus down) and the water collective coordinate n(t) "
                      "against time, summed over channels")
    for k, lab in enumerate(labels):
        colour = colours[k % len(colours)]
        axes[0].plot(time_ns, summed(lab, "occupancy"), color=colour, lw=1.2, label=lab)
        axes[1].step(time_ns, summed(lab, "net_crossings"), where="post", color=colour, lw=1.6, label=lab)
    for k, lab in enumerate(water):
        axes[2].plot(time_ns, summed(lab, "n"), color=colours[labels.index(lab) % len(colours)], lw=1.2, label=lab)
    text.axis_label(axes[0].set_ylabel, "Inside lumen")
    text.axis_label(axes[1].set_ylabel, "Net crossings")
    if water:
        text.axis_label(axes[2].set_ylabel, "n(t) (waters)")
    text.axis_label(axes[-1].set_xlabel, "Time (ns)")
    for ax in axes:
        ax.grid(True, color="#d8d8d8", lw=.6, alpha=.75)
        _style_axes(ax)
    axes[1].axhline(0, color=CREVICE_GRAY, lw=.8)
    if len(labels) > 1:
        axes[0].legend(frameon=False, fontsize=8, ncol=min(4, len(labels)))
    totals = result.summary["totals_over_channels"]
    text.title(axes[0], "Permeation events · " + ", ".join(f"{t['species']} {t['events_up']} up / {t['events_down']} down"
                                                          for t in totals), fontsize=10)
    _save_figure(fig, Path(path), dpi=dpi, text=text)
    return fig

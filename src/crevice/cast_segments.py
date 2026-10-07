"""Channel cast segments (ENTRY, LUMEN, EXIT) and the residues lining them.

A resolved channel cast (:func:`crevice.voids.void_cast` of an axial
profile) is the **lumen**: the probe-swept space between the two channel
mouths. Its volume, ``total_volume`` and every profile value are unchanged
here. Beyond each mouth this module adds the space-filling vestibule out to
bulk solvent, so the whole route is reported as separate segments.

Segments (all distances in Å, atoms as van der Waals spheres of the radius set
in force):

* **LUMEN** (``lumen``): the axial channel cast between the mouth planes.
* **Open axial mouth**: one vestibule segment. Its centre line runs straight
  from the mouth sample's centre along the channel axis, away from the lumen,
  every 0.25 Å until the first sample in bulk (at most 60 Å, and only while
  the sample is outside every atom).
* **Capped mouth** (``lateral_exits``): one segment per lateral exit leg; the
  leg's own centre line (:func:`crevice.channel_exits.lateral_exit_legs`).
* Both kinds are filled by the same probe-swept rule as the lateral exit
  casts (:func:`crevice.channel_exits.lateral_exit_casts`): on the lumen's
  lattice, only planes beyond the mouth plane; core nodes inside the free
  sphere of a centre-line sample with clearance at least the enclosure probe,
  swept by (enclosure probe - min_radius), connected to the centre line.
* **Outer boundary.** Nodes in bulk are excluded. Bulk is the space swept by a
  rolling probe of ``bulk_radius`` (default 6 Å, CLI ``--exit-bulk-radius``)
  that is connected to the outside of the protein
  (:class:`crevice.channel_exits.BulkEnvelope`). This radius sets *where the
  vestibule ends* (the solvent boundary), not how wide the segment may be:
  inside that boundary the segment fills all probe-accessible space attached
  to its centre line. A larger bulk probe moves the boundary outward.

Segments are disjoint from the lumen and from each other (different
half-spaces; nodes shared by legs at one end go to the nearest centre line)
and are reported separately; their volumes are never added to
``total_volume``. ``segments_total_volume_A3`` (lumen plus every segment) is
a separate, clearly named sum.

Labelling rule (geometric labels only; **"entry" and "exit" do not claim a
transport direction**):

* ``entry_end="auto"`` (default): if exactly one end is capped and leaves
  through lateral exit legs, those legs are EXITs and the open axial end is
  the ENTRY; otherwise the axis start (lower end, first profile sample) is
  the ENTRY and the axis end is the EXIT.
* ``entry_end="start"`` or ``"end"`` forces the ENTRY to that end of the
  axis (CLI ``--entry-end``).
* Several segments of one type are numbered by leg rank: ``exit_1``,
  ``exit_2``, ...

Lining residues: a residue lines a segment when any of its atoms has its
centre within ``cutoff`` (default 3.3 Å, a hydrogen-bond donor-acceptor
distance; CLI ``--lining-cutoff``) of a measured grid sample of that segment
(the cast nodes, which are probe-centre positions clear of every atomic
van der Waals sphere of the radius set in force). The distance is atom centre
to nearest grid node, so on a grid of spacing ``h`` it can overestimate the
distance to the continuous cast by at most ``h * sqrt(3) / 2``.

Viewer objects and colours (see
:func:`crevice.volume_export.write_volume_viewer_bundle`): PyMOL objects
``crevice_entry_N`` (yellow ``#f0d43a``), ``crevice_lumen`` (teal ``#1494a1``)
and ``crevice_exit_N`` (rust ``#c2410c``); lining residues as stick objects
``crevice_<segment>_lining`` (disabled by default) and selections
``crevice_<segment>_lining_sel``.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

from .presentation import CHANNEL_CAST_RGB

#: Segment types in display order.
SEGMENT_TYPES = ("entry", "lumen", "exit")
#: Default segment surface colours as RGB in 0-1. Entry yellow ``#f0d43a``,
#: lumen teal ``#1494a1`` (the channel colour), exit rust ``#c2410c``.
SEGMENT_RGB = {"entry": (0.941, 0.831, 0.227), "lumen": CHANNEL_CAST_RGB, "exit": (0.761, 0.255, 0.047)}
#: The same colours as hex strings.
SEGMENT_HEX = {"entry": "#f0d43a", "lumen": "#1494a1", "exit": "#c2410c"}
#: Default lining cutoff, Å (atom centre to nearest segment grid node).
DEFAULT_LINING_CUTOFF_A = 3.3
#: Longest straight vestibule centre line beyond an open axial mouth, Å.
MAX_VESTIBULE_LENGTH_A = 60.0

VESTIBULE_DEFINITION = (
    "probe-swept vestibule beyond an open axial mouth: centre line straight along the channel axis from the "
    "mouth sample's centre, away from the lumen, every 0.25 A until the first sample in bulk; filled with the "
    "lateral-exit cast rule (lumen lattice, beyond the mouth plane, outside bulk, core nodes in centre-line free "
    "spheres with clearance >= enclosure probe, swept by enclosure probe - min_radius, connected to the centre "
    "line). Reported separately; never part of the axial channel volume")
SEGMENT_DEFINITION = (
    "LUMEN = axial channel cast between the mouth planes (total_volume); ENTRY/EXIT = space-filling probe-swept "
    "casts beyond each mouth out to the bulk-solvent boundary (rolling bulk probe of bulk_radius_A connected to "
    "the outside; an outer boundary, not a width limit): the axial vestibule at an open mouth, one segment per "
    "lateral exit leg at a capped mouth. Disjoint; volumes reported separately")
LABEL_RULE_DEFINITION = (
    "auto: if exactly one end is capped with lateral exit legs, those legs are EXITs and the open axial end is the "
    "ENTRY; otherwise the axis start (lower end) is the ENTRY and the axis end the EXIT. entry_end='start'|'end' "
    "forces the ENTRY end. Geometric labels only, not a claim about transport direction")
LINING_DEFINITION = (
    "residue with any atom centre within cutoff_A of a measured grid node of the segment (probe-centre positions "
    "clear of every atomic van der Waals sphere of the radius set in force); distance = atom centre to nearest "
    "node, overestimating the distance to the continuous cast by at most spacing*sqrt(3)/2")


def segment_roles(profile, entry_end: str | None = "auto") -> tuple[dict[str, str], str]:
    """Assign ENTRY and EXIT to the two ends of a channel profile.

    Parameters
    ----------
    profile : PoreProfile
        Resolved channel profile; ``metadata["exits"]`` (when present) tells
        which ends are capped (``"lateral"``).
    entry_end : {"auto", "start", "end"} or None, default "auto"
        ``None`` means ``"auto"``. See the module docstring for the rule.

    Returns
    -------
    roles : dict of str to str
        ``{"lower": "entry" | "exit", "upper": "entry" | "exit"}``; ``lower``
        is the axis start (first profile sample).
    rule : str
        ``"capped_end_exits"``, ``"axis_start_entry"``, ``"override_start"``
        or ``"override_end"``.

    Raises
    ------
    ValueError
        For any other ``entry_end``.

    Examples
    --------
    >>> from crevice.models import PoreProfile
    >>> segment_roles(PoreProfile((0, 0, 0), (0, 0, 1), ()))
    ({'lower': 'entry', 'upper': 'exit'}, 'axis_start_entry')
    >>> segment_roles(PoreProfile((0, 0, 0), (0, 0, 1), (), metadata={"exits": {
    ...     "lower": {"type": "lateral"}, "upper": {"type": "axial"}}}))
    ({'lower': 'exit', 'upper': 'entry'}, 'capped_end_exits')
    """
    choice = "auto" if entry_end is None else entry_end
    exits = profile.metadata.get("exits") or {}
    if choice == "auto":
        capped = [end for end in ("lower", "upper") if (exits.get(end) or {}).get("type") == "lateral"]
        if len(capped) == 1:
            entry, rule = ("upper" if capped[0] == "lower" else "lower"), "capped_end_exits"
        else:
            entry, rule = "lower", "axis_start_entry"
    elif choice == "start":
        entry, rule = "lower", "override_start"
    elif choice == "end":
        entry, rule = "upper", "override_end"
    else:
        raise ValueError("entry_end must be 'auto', 'start' or 'end'")
    return {end: ("entry" if end == entry else "exit") for end in ("lower", "upper")}, rule


def axial_vestibule_leg(frame, profile, cast, end: str, *, bulk_radius: float = 6.0,
                        bulk_spacing: float = 1.0, step: float = 0.25) -> dict | None:
    """Straight centre line beyond an open axial mouth, out to bulk.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile and cast (the cast's atom selection is used).
    profile : PoreProfile
        Resolved channel profile.
    cast : VoidCast
        Its connected channel cast.
    end : {"lower", "upper"}
        Mouth to continue from.
    bulk_radius : float, default 6.0
        Rolling bulk-probe radius, Å.
    bulk_spacing : float, default 1.0
        Bulk-envelope grid spacing, Å.
    step : float, default 0.25
        Sample spacing along the centre line, Å.

    Returns
    -------
    dict or None
        Leg record for :func:`crevice.channel_exits.probe_swept_leg_casts`
        (``points``, ``raw_clearances_A``, ``leg`` = 1, ``bulk_radius_A``,
        ``bulk_grid_spacing_A``, ``bottleneck_radius_A`` = smallest raw
        clearance, ``length_A``, ``reaches_bulk``), or ``None`` when no sample
        beyond the mouth is clear of the atoms.
    """
    from .channel_exits import BulkEnvelope
    from .geometry import add, scale
    from .spatial import SpatialIndex
    meta = cast.metadata
    atoms = frame.selected_atoms(include_hydrogen=meta.get("include_hydrogen", False),
                                 include_hetero=meta.get("include_hetero", True))
    spatial = SpatialIndex(atoms)
    sign = -1 if end == "lower" else 1
    mouth = profile.points[0] if end == "lower" else profile.points[-1]
    basis = meta["grid_basis"]
    bulk = BulkEnvelope(spatial.atoms, profile.axis_origin, basis, mouth.t, sign,
                        spacing=bulk_spacing, bulk_radius=bulk_radius)
    direction = scale(profile.axis_direction, sign)
    points, raw = [], []
    reached = False
    for k in range(1, int(MAX_VESTIBULE_LENGTH_A/step)+1):
        position = add(mouth.position, scale(direction, k*step))
        clearance = spatial.nearest_surface(position)[0]
        if clearance <= 0:
            break
        points.append(tuple(position))
        raw.append(float(clearance))
        if bulk.depth(position) <= 0.0:
            reached = True
            break
    if not points:
        return None
    return {"points": [list(p) for p in points], "raw_clearances_A": raw, "leg": 1,
            "bulk_radius_A": bulk_radius, "bulk_grid_spacing_A": bulk_spacing,
            "bottleneck_radius_A": min(raw), "length_A": len(points)*step, "reaches_bulk": reached}


def segment_casts(frame, profile, cast, *, entry_end: str | None = "auto", bulk_radius: float | None = None,
                  lateral_casts=None) -> tuple:
    """ENTRY and EXIT segment casts beyond the two mouths of a channel cast.

    The lumen is ``cast`` itself and is not modified. Rules: module docstring.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the profile and cast.
    profile : PoreProfile
        Resolved channel profile (``method == "axial-connected"``).
    cast : VoidCast
        Its connected channel cast (the LUMEN).
    entry_end : {"auto", "start", "end"} or None, default "auto"
        Labelling rule (:func:`segment_roles`).
    bulk_radius : float, optional
        Rolling bulk-probe radius (Å) for open axial ends; default the radius
        recorded by the lateral exit search, else 6.0.
    lateral_casts : sequence of VoidComponent, optional
        Already computed :func:`crevice.channel_exits.lateral_exit_casts`; computed
        here when omitted.

    Returns
    -------
    tuple of VoidComponent
        Entries first, then exits, each numbered by leg rank. Lateral legs keep
        kind ``"lateral_exit"``; axial vestibules have kind
        ``"axial_vestibule"``. Each component's metadata gains ``segment``
        (e.g. ``"entry_1"``), ``segment_type``, ``segment_index``, ``source``,
        ``label_rule`` and ``segment_definition``. Empty when the cast is not
        a connected channel cast or is focus-cropped.

    Examples
    --------
    >>> from crevice.models import PoreProfile
    >>> segment_casts(None, PoreProfile((0, 0, 0), (0, 0, 1), ()), None)
    ()
    """
    from dataclasses import replace
    from .channel_exits import lateral_exit_casts, probe_swept_leg_casts
    if (cast is None or cast.metadata.get("export_mode") != "connected_section_fill"
            or cast.metadata.get("focus_points") or not profile.points):
        return ()
    roles, rule = segment_roles(profile, entry_end)
    exits = profile.metadata.get("exits") or {}
    lateral = tuple(lateral_casts) if lateral_casts is not None else lateral_exit_casts(frame, profile, cast)
    if bulk_radius is None:
        bulk_radius = next((leg.get("bulk_radius_A") for e in exits.values() for leg in (e or {}).get("legs", ())
                            if leg.get("bulk_radius_A")), 6.0)
    next_id = max((c.id for c in (*cast.components, *lateral)), default=0)+1
    found = {"lower": [], "upper": []}
    for c in lateral:
        found[c.metadata["end"]].append(c)
    for end in ("lower", "upper"):
        if (exits.get(end) or {}).get("type") == "lateral":
            continue
        leg = axial_vestibule_leg(frame, profile, cast, end, bulk_radius=bulk_radius)
        if leg is None:
            continue
        made = probe_swept_leg_casts(frame, profile, cast, [(end, leg)], kind="axial_vestibule",
                                     definition=VESTIBULE_DEFINITION, first_id=next_id)
        for c in made:
            c.metadata.update({"centre_line_reaches_bulk": leg["reaches_bulk"], "bulk_radius_A": bulk_radius})
        found[end].extend(made)
        next_id += len(made)
    result = []
    for kind in ("entry", "exit"):
        for end in ("lower", "upper"):
            if roles[end] != kind:
                continue
            for index, c in enumerate(sorted(found[end], key=lambda c: c.metadata["leg"]), 1):
                result.append(replace(c, metadata={**c.metadata, "segment": f"{kind}_{index}",
                                                   "segment_type": kind, "segment_index": index,
                                                   "source": c.kind, "label_rule": rule,
                                                   "segment_definition": SEGMENT_DEFINITION}))
    return tuple(result)


def segment_summary(cast, segments) -> list[dict]:
    """One JSON/CSV row per segment, LUMEN included.

    Parameters
    ----------
    cast : VoidCast
        The channel cast (LUMEN).
    segments : sequence of VoidComponent
        From :func:`segment_casts`.

    Returns
    -------
    list of dict
        Rows in ENTRY, LUMEN, EXIT order: ``segment``, ``segment_type``,
        ``segment_index``, ``end`` (``lower``/``upper``; ``both`` for the
        lumen), ``source`` (``axial_lumen``, ``axial_vestibule`` or
        ``lateral_exit``), ``component_ids``, ``volume_A3``, ``point_count``,
        ``min_clearance_A``, ``max_clearance_A`` and ``label_rule``.

    Examples
    --------
    >>> from crevice.models import VoidCast
    >>> rows = segment_summary(VoidCast((), (), 0.5, 0.0, "channel"), ())
    >>> rows[0]["segment"], rows[0]["volume_A3"]
    ('lumen', 0)
    """
    points = [p for c in cast.components for p in c.points]
    lumen = {"segment": "lumen", "segment_type": "lumen", "segment_index": 1, "end": "both",
             "source": "axial_lumen", "component_ids": [c.id for c in cast.components],
             "volume_A3": cast.total_volume, "point_count": len(points),
             "min_clearance_A": min((p.raw_clearance for p in points), default=None),
             "max_clearance_A": max((p.raw_clearance for p in points), default=None),
             "label_rule": segments[0].metadata["label_rule"] if segments else None}
    rows = [{"segment": c.metadata["segment"], "segment_type": c.metadata["segment_type"],
             "segment_index": c.metadata["segment_index"], "end": c.metadata["end"],
             "source": c.metadata["source"], "component_ids": [c.id], "volume_A3": c.volume,
             "point_count": c.point_count, "min_clearance_A": c.min_radius, "max_clearance_A": c.max_radius,
             "label_rule": c.metadata["label_rule"]} for c in segments]
    return ([r for r in rows if r["segment_type"] == "entry"] + [lumen]
            + [r for r in rows if r["segment_type"] == "exit"])


def segment_points(cast, segments) -> dict:
    """Measured grid samples of each segment, LUMEN first.

    Parameters
    ----------
    cast : VoidCast
        The channel cast (LUMEN).
    segments : sequence of VoidComponent
        From :func:`segment_casts`.

    Returns
    -------
    dict of str to tuple of ChannelPoint
        ``{"lumen": ..., "entry_1": ..., ...}``.
    """
    out = {"lumen": tuple(p for c in cast.components for p in c.points)}
    for c in segments:
        out[c.metadata["segment"]] = tuple(c.points)
    return out


def viewer_residue_ids(frame) -> dict:
    """Viewer chain and residue number of every residue of ``frame``.

    The same numbering as :func:`crevice.presentation.write_viewer_structure`
    (residues renumbered 1, 2, ... per chain; multi-character chains mapped
    to unused single characters).

    Parameters
    ----------
    frame : StructureFrame
        The structure written to the viewer PDB.

    Returns
    -------
    dict of ResidueKey to (str, int)
        Viewer chain and residue number.
    """
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    chains = list(dict.fromkeys(a.chain_id for a in frame.atoms))
    chain_map = {c: (c if len(c) == 1 and c in alphabet else None) for c in chains}
    available = [c for c in alphabet if c not in chain_map.values()]
    for c in chains:
        if chain_map[c] is None:
            chain_map[c] = available.pop(0)
    out, counts = {}, {}
    for a in frame.atoms:
        if a.residue_key not in out:
            counts[a.chain_id] = counts.get(a.chain_id, 0)+1
            out[a.residue_key] = (chain_map[a.chain_id], counts[a.chain_id])
    return out


def lining_residues(frame, cast, segments=(), *, cutoff: float = DEFAULT_LINING_CUTOFF_A) -> list[dict]:
    """Residues lining the LUMEN and each ENTRY/EXIT segment.

    Parameters
    ----------
    frame : StructureFrame
        Structure of the cast; its atoms (filtered by the cast's
        ``include_hydrogen``/``include_hetero`` settings) are tested, and its
        full atom list defines the viewer numbering (:func:`viewer_residue_ids`).
    cast : VoidCast
        The channel cast (LUMEN).
    segments : sequence of VoidComponent, optional
        From :func:`segment_casts`.
    cutoff : float, default 3.3
        Atom centre to nearest segment grid node, Å (``LINING_DEFINITION``).

    Returns
    -------
    list of dict
        One row per (segment, residue), segments in ENTRY, LUMEN, EXIT order
        and residues in structure order: ``segment``, ``segment_type``,
        ``role`` (``"<segment>_lining"``), ``residue`` (CREVICE label),
        ``chain_id``, ``resname``, ``resid``, ``icode``, ``hetero``,
        ``min_distance_A``, ``closest_atom``, ``atoms_within``, ``cutoff_A``,
        ``viewer_chain``, ``viewer_resid`` and ``viewer_atom_serials``
        (1-based positions in ``frame.atoms``, as in the viewer PDB).

    Raises
    ------
    ValueError
        If ``cutoff`` is not finite and positive.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("lining cutoff must be finite and positive")
    meta = cast.metadata
    selected = set(map(id, frame.selected_atoms(include_hydrogen=meta.get("include_hydrogen", False),
                                                include_hetero=meta.get("include_hetero", True))))
    atoms = [(i, a) for i, a in enumerate(frame.atoms, 1) if id(a) in selected]
    if not atoms:
        return []
    xyz = np.asarray([a.coord for _, a in atoms], dtype=float)
    viewer = viewer_residue_ids(frame)
    serials = {}
    for i, a in enumerate(frame.atoms, 1):
        serials.setdefault(a.residue_key, []).append(i)
    order = {key: n for n, key in enumerate(viewer)}
    groups = segment_points(cast, segments)
    names = ([n for n in groups if n.startswith("entry")] + ["lumen"] + [n for n in groups if n.startswith("exit")])
    rows = []
    for name in names:
        points = groups[name]
        if not points:
            continue
        distance, _ = cKDTree(np.asarray([p.position for p in points])).query(xyz, distance_upper_bound=cutoff)
        best = {}
        for (_, atom), d in zip(atoms, distance):
            if d <= cutoff:
                key = atom.residue_key
                record = best.setdefault(key, [math.inf, "", 0, atom])
                record[2] += 1
                if d < record[0]:
                    record[0], record[1] = float(d), atom.name
        for key in sorted(best, key=order.get):
            d, closest, count, atom = best[key]
            rows.append({"segment": name, "segment_type": name.split("_")[0], "role": f"{name}_lining",
                         "residue": key.label, "chain_id": key.chain_id, "resname": key.resname,
                         "resid": key.resid, "icode": key.icode.strip(), "hetero": bool(atom.hetero),
                         "min_distance_A": round(d, 4), "closest_atom": closest, "atoms_within": count,
                         "cutoff_A": cutoff, "viewer_chain": viewer[key][0], "viewer_resid": viewer[key][1],
                         "viewer_atom_serials": serials[key]})
    return rows


LINING_FIELDS = ("segment", "segment_type", "role", "residue", "chain_id", "resname", "resid", "icode", "hetero",
                 "min_distance_A", "closest_atom", "atoms_within", "cutoff_A", "viewer_chain", "viewer_resid")
SEGMENT_FIELDS = ("segment", "segment_type", "segment_index", "end", "source", "component_ids", "volume_A3",
                  "point_count", "min_clearance_A", "max_clearance_A", "label_rule")


def write_segments_csv(rows, path) -> None:
    """Write :func:`segment_summary` rows as CSV.

    Columns: ``SEGMENT_FIELDS`` (``component_ids`` semicolon-separated;
    volumes in Å³, clearances in Å).

    Parameters
    ----------
    rows : sequence of dict
        From :func:`segment_summary`.
    path : str or Path
        Output CSV.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SEGMENT_FIELDS), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "component_ids": ";".join(map(str, row["component_ids"]))})


def write_lining_csv(rows, path) -> None:
    """Write :func:`lining_residues` rows as CSV.

    Columns: ``LINING_FIELDS`` (distances in Å; ``role`` is
    ``"<segment>_lining"``; ``viewer_chain``/``viewer_resid`` address the
    bundle's viewer PDB).

    Parameters
    ----------
    rows : sequence of dict
        From :func:`lining_residues`.
    path : str or Path
        Output CSV.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(LINING_FIELDS), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def add_lining_columns(path, rows, *, residue_field: str = "residue") -> bool:
    """Append ``lumen_lining`` and ``lining_segments`` columns to a residue CSV.

    ``lumen_lining`` is ``True``/``False`` for rows with a residue label and
    empty otherwise (for example region nodes of a network table);
    ``lining_segments`` lists every segment the residue lines, semicolon-separated.

    Parameters
    ----------
    path : str or Path
        Existing CSV with a residue-label column.
    rows : sequence of dict
        From :func:`lining_residues`.
    residue_field : str, default "residue"
        Column holding CREVICE residue labels.

    Returns
    -------
    bool
        ``False`` (file untouched) when the CSV lacks ``residue_field`` or
        already has the columns.
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        table = list(reader)
    if residue_field not in fields or "lumen_lining" in fields:
        return False
    lines = {}
    for row in rows:
        lines.setdefault(row["residue"], []).append(row["segment"])
    for row in table:
        label = row.get(residue_field) or ""
        found = lines.get(label, [])
        row["lumen_lining"] = ("lumen" in found) if label else ""
        row["lining_segments"] = ";".join(found)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields+["lumen_lining", "lining_segments"])
        writer.writeheader()
        writer.writerows(table)
    return True

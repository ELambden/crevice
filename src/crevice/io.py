"""Writers for CREVICE results: JSON records, CSV tables, PDB dummy atoms and a PyMOL script.

Every analysis result object (:mod:`crevice.models`) has a JSON writer here
that stores its full ``to_dict()`` record, and every tabular dataset has a CSV
writer whose column names carry their unit (``_A`` = Å, ``_A2`` = Å²,
``_A3`` = Å³, ``_ps``, ``_degrees``). The PDB writers store sample points as
dummy ``HETATM`` records with the radius (Å) in the B-factor column, so a
molecular viewer can draw them as spheres of that radius.

Writers create missing parent directories and overwrite existing files. JSON
files are indented with sorted keys. Nothing is plotted here; figures are in
:mod:`crevice.figures`.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from .models import Cavity, ChannelPoint, PoreProfile, ResidueContact, ResidueInteractionNetwork, TrajectoryAnalysis, Tunnel, VoidCast


def write_profile_json(profile: PoreProfile, path: str | Path) -> None:
    """Write a pore profile as JSON.

    The file holds :meth:`crevice.models.PoreProfile.to_dict`: the method, axis
    origin and direction, probe radius, minimum/maximum/mean radius (Å), length
    (Å), volume estimate (Å³), the bottleneck sample, every sample (centre ``x``,
    ``y``, ``z``, ``radius``, ``raw_clearance``, ``axis_position`` in Å, nearest
    atom serial and residue) and the profile metadata (enclosure-probe choice,
    exits, settings, radius set).

    Parameters
    ----------
    profile : PoreProfile
        Resolved profile, for example from :func:`crevice.channels.pore_profile`.
    path : str or Path
        Output JSON path.
    """
    _write_json(profile.to_dict(), path)


# CSV column names carry their unit (ASCII: _A = Å, _A2 = Å², _A3 = Å³, _ps, _degrees).
POINT_COLUMNS = {"x": "x_A", "y": "y_A", "z": "z_A", "radius": "radius_A",
                 "raw_clearance": "raw_clearance_A", "axis_position": "axis_position_A"}


def _point_row(point: ChannelPoint) -> dict[str, Any]:
    return {POINT_COLUMNS.get(key, key): value for key, value in point.to_dict().items()}


def write_profile_csv(profile: PoreProfile, path: str | Path) -> None:
    """Write one row per profile sample.

    Columns: ``index``, ``x_A``, ``y_A``, ``z_A`` (sample centre, Å), ``radius_A``
    (reported radius, Å), ``raw_clearance_A`` (atom-surface clearance before
    any probe subtraction, Å), ``axis_position_A`` (signed axial coordinate
    ``t``, Å), ``nearest_atom_serial`` and ``nearest_residue``.

    Parameters
    ----------
    profile : PoreProfile
        Resolved profile.
    path : str or Path
        Output CSV (parent directories are created).
    """
    _write_rows(path, [_point_row(point) for point in profile.points])


def write_profile_pdb(profile: PoreProfile, path: str | Path) -> None:
    """Write the profile samples as dummy atoms for a molecular viewer.

    One ``HETATM`` record per sample, residue name ``PRF``, atom name ``X``,
    chain ``A``, residue number = serial, coordinates = sample centre (Å) and
    B-factor = reported radius (Å, clipped to 0-99.99). Load it with the protein
    and show the atoms as spheres scaled by B-factor to see the pore lumen (see
    :func:`export_pymol_script`).

    Parameters
    ----------
    profile : PoreProfile
        Resolved profile.
    path : str or Path
        Output PDB path.
    """
    write_points_pdb(profile.points, path, resname="PRF")


def write_cavities_json(cavities: Sequence[Cavity], path: str | Path) -> None:
    """Write detected cavities as JSON.

    The file holds ``{"cavities": [...]}`` with each
    :meth:`crevice.models.Cavity.to_dict`: ``id``, ``center`` (Å), ``radius``
    (largest atom-surface clearance, Å), ``volume`` (Å³), ``grid_points``,
    ``nearest_residues``, ``score`` and metadata.

    Parameters
    ----------
    cavities : sequence of Cavity
        From :func:`crevice.analysis.detect_cavities`.
    path : str or Path
        Output JSON path.
    """
    _write_json({"cavities": [cavity.to_dict() for cavity in cavities]}, path)


def write_cavities_pdb(cavities: Sequence[Cavity], path: str | Path) -> None:
    """Write one dummy atom per cavity at its centre.

    Residue name ``CAV``, residue number = serial, B-factor = the cavity's largest
    clearance radius (Å). The cavity ``id`` is stored as the point index; the
    cavity shape itself is not represented (use a void cast for that).

    Parameters
    ----------
    cavities : sequence of Cavity
        From :func:`crevice.analysis.detect_cavities`.
    path : str or Path
        Output PDB path.
    """
    points = [
        ChannelPoint(
            index=cavity.id,
            position=cavity.center,
            radius=cavity.radius,
            raw_clearance=cavity.radius,
            t=float(cavity.id),
        )
        for cavity in cavities
    ]
    write_points_pdb(points, path, resname="CAV")


def write_tunnels_json(tunnels: Sequence[Tunnel], path: str | Path) -> None:
    """Write tunnels as JSON.

    The file holds ``{"tunnels": [...]}`` with each
    :meth:`crevice.models.Tunnel.to_dict`: ``id``, ``start``/``end`` (Å),
    ``length`` (Å), ``bottleneck_radius`` and ``mean_radius`` (Å), ``score``,
    ``nearest_residues``, every path node (``axis_position`` = cumulative path
    length from the start, Å) and metadata (grid spacing, whether it was
    enlarged).

    Parameters
    ----------
    tunnels : sequence of Tunnel
        From :func:`crevice.tunnels.find_tunnels`.
    path : str or Path
        Output JSON path.
    """
    _write_json({"tunnels": [tunnel.to_dict() for tunnel in tunnels]}, path)


def write_tunnels_pdb(tunnels: Sequence[Tunnel], path: str | Path) -> None:
    """Write every tunnel path node as a dummy atom.

    Residue name ``TUN``; atoms of all tunnels are numbered consecutively, in
    tunnel order then path order, with B-factor = node clearance radius (Å). The
    tunnel ID is not stored; use the CSV or JSON to separate tunnels.

    Parameters
    ----------
    tunnels : sequence of Tunnel
        From :func:`crevice.tunnels.find_tunnels`.
    path : str or Path
        Output PDB path.
    """
    points: list[ChannelPoint] = []
    serial = 1
    for tunnel in tunnels:
        for point in tunnel.points:
            points.append(
                ChannelPoint(
                    index=serial,
                    position=point.position,
                    radius=point.radius,
                    raw_clearance=point.raw_clearance,
                    t=point.t,
                    nearest_atom_serial=point.nearest_atom_serial,
                    nearest_residue=point.nearest_residue,
                )
            )
            serial += 1
    write_points_pdb(points, path, resname="TUN")


def write_void_cast_json(cast: VoidCast, path: str | Path, *, include_points: bool = False) -> None:
    """Write a void cast as JSON.

    The file holds :meth:`crevice.models.VoidCast.to_dict`: method, mode, grid
    spacing and minimum clearance (Å), component count, total volume (Å³), each
    component's summary and the cast metadata.

    Parameters
    ----------
    cast : VoidCast
        From :func:`crevice.voids.void_cast` or :func:`crevice.rolling.rolling_probe_cast`.
    path : str or Path
        Output JSON path.
    include_points : bool, default False
        Also store every retained grid sample. This can make the file large; the
        measured volume itself is better exported as a DX map
        (:func:`crevice.volume_export.write_void_cast_dx`).
    """
    _write_json(cast.to_dict(include_points=include_points), path)


def write_void_cast_pdb(cast: VoidCast, path: str | Path) -> None:
    """Write every retained cast sample as a dummy atom.

    Residue name ``VOI``, B-factor = the sample's atom-surface clearance (Å). Large
    casts give large files; the DX maps of :mod:`crevice.volume_export` are the
    preferred volume representation.

    Parameters
    ----------
    cast : VoidCast
        Void cast whose ``points`` are written.
    path : str or Path
        Output PDB path.
    """
    write_points_pdb(cast.points, path, resname="VOI")


def write_network_json(network: ResidueInteractionNetwork, path: str | Path) -> None:
    """Write a residue interaction network as JSON.

    The file holds :meth:`crevice.models.ResidueInteractionNetwork.to_dict`:
    ``nodes`` (id, kind ``residue`` or ``region``, label, residue, properties,
    metadata), ``edges`` (source, target, interaction class, distance in Å,
    weight, metadata) and the construction settings.

    Parameters
    ----------
    network : ResidueInteractionNetwork
        From :func:`crevice.networks.build_residue_network` or
        :func:`crevice.networks.build_cavity_network`.
    path : str or Path
        Output JSON path.
    """
    _write_json(network.to_dict(), path)


def write_trajectory_json(analysis: TrajectoryAnalysis, path: str | Path) -> None:
    """Write a trajectory analysis as JSON.

    The file holds :meth:`crevice.models.TrajectoryAnalysis.to_dict`: per frame
    the profile (or its unresolved status and reason), cavities, tunnels, network
    counts, alignment record and features, plus the analysis metadata. The same
    per-frame values are tabulated by :func:`write_trajectory_csv`.

    Parameters
    ----------
    analysis : TrajectoryAnalysis
        From :func:`crevice.trajectory.analyze_trajectory`.
    path : str or Path
        Output JSON path.
    """
    _write_json(analysis.to_dict(), path)


def write_residue_contacts_csv(contacts: Sequence[ResidueContact], path: str | Path) -> None:
    """Write one row per contacting residue.

    Columns: ``residue``, ``min_distance_A`` and ``mean_distance_A`` (smallest
    and mean atom-surface gap to the region surface, Å), ``atom_count``,
    ``point_count`` (region samples contacted), ``role`` (``bottleneck``,
    ``bottleneck-nearby``, ``lining`` or ``nearby``), ``properties``
    (semicolon-separated residue classes) and ``influence_score`` (dimensionless).

    Parameters
    ----------
    contacts : sequence of ResidueContact
        From :func:`crevice.analysis.annotate_residues`.
    path : str or Path
        Output CSV (parent directories are created).
    """
    fields = ["residue", "min_distance_A", "mean_distance_A", "atom_count", "point_count",
              "role", "properties", "influence_score"]
    rows = []
    for contact in contacts:
        row = contact.to_dict()
        rows.append({"residue": row["residue"], "min_distance_A": row["min_distance"],
                     "mean_distance_A": row["mean_distance"], "atom_count": row["atom_count"],
                     "point_count": row["point_count"], "role": row["role"],
                     "properties": ";".join(row["properties"]), "influence_score": row["influence_score"]})
    _write_rows(path, rows, fields)


def write_residue_contacts_json(contacts: Sequence[ResidueContact], path: str | Path) -> None:
    """Write residue contacts as JSON.

    The file holds ``{"residue_contacts": [...]}`` with each
    :meth:`crevice.models.ResidueContact.to_dict`. The same fields are written as
    CSV by :func:`write_residue_contacts_csv`, which the command-line tools use.

    Parameters
    ----------
    contacts : sequence of ResidueContact
        From :func:`crevice.analysis.annotate_residues`.
    path : str or Path
        Output JSON path.
    """
    _write_json({"residue_contacts": [contact.to_dict() for contact in contacts]}, path)


def write_void_cast_csv(cast: VoidCast, path: str | Path, *, labels: dict[int, str] | None = None,
                        extra_components=()) -> None:
    """Write one row per retained cast region.

    Columns: ``region_id``, ``label`` (display-region label, empty when the
    cast has a single region), ``kind``, ``volume_A3``, ``point_count``,
    ``center_x_A``/``center_y_A``/``center_z_A``, ``min_radius_A``,
    ``mean_radius_A``, ``max_radius_A`` (sample atom-surface clearances, Å),
    ``spacing_A``, ``score``, ``boundary_faces`` and ``nearest_residues``
    (semicolon-separated). The header is written even when there are no regions.

    Parameters
    ----------
    cast : VoidCast
        Cast whose components are written.
    path : str or Path
        Output CSV.
    labels : dict of int to str, optional
        Display-region label per component ID.
    extra_components : sequence of VoidComponent, optional
        Components reported beside the cast but not part of it, such as
        lateral exit leg casts (:func:`crevice.channel_exits.lateral_exit_casts`).
        One row each after the cast's rows, with ``kind`` ``lateral_exit`` and
        ``label`` ``lateral_exit_<end>_<leg>``; their volumes are not part of
        ``cast.total_volume``.
    """
    fields = ["region_id", "label", "kind", "volume_A3", "point_count", "center_x_A", "center_y_A",
              "center_z_A", "min_radius_A", "mean_radius_A", "max_radius_A", "spacing_A", "score",
              "boundary_faces", "nearest_residues"]
    labels = labels or {}
    rows = [{"region_id": c.id, "label": labels.get(c.id, ""), "kind": c.kind, "volume_A3": c.volume,
             "point_count": c.point_count, "center_x_A": c.center[0], "center_y_A": c.center[1],
             "center_z_A": c.center[2], "min_radius_A": c.min_radius, "mean_radius_A": c.mean_radius,
             "max_radius_A": c.max_radius, "spacing_A": cast.spacing, "score": c.score,
             "boundary_faces": ";".join(c.boundary_faces), "nearest_residues": ";".join(c.nearest_residues)}
            for c in cast.components]
    rows += [{"region_id": c.id, "label": f"lateral_exit_{c.metadata.get('end')}_{c.metadata.get('leg')}",
              "kind": c.kind, "volume_A3": c.volume, "point_count": c.point_count, "center_x_A": c.center[0],
              "center_y_A": c.center[1], "center_z_A": c.center[2], "min_radius_A": c.min_radius,
              "mean_radius_A": c.mean_radius, "max_radius_A": c.max_radius, "spacing_A": cast.spacing,
              "score": c.score, "boundary_faces": ";".join(c.boundary_faces),
              "nearest_residues": ";".join(c.nearest_residues)} for c in extra_components]
    _write_rows(path, rows, fields)


def write_cavities_csv(cavities: Sequence[Cavity], path: str | Path) -> None:
    """Write one row per detected cavity.

    Columns: ``cavity_id``, ``center_x_A``/``center_y_A``/``center_z_A``,
    ``max_radius_A`` (largest atom-surface clearance, Å), ``volume_A3``,
    ``grid_points``, ``score`` and ``nearest_residues`` (semicolon-separated).
    Together with the manifest this holds every field of the cavities JSON.

    Parameters
    ----------
    cavities : sequence of Cavity
        From :func:`crevice.analysis.detect_cavities`.
    path : str or Path
        Output CSV.
    """
    fields = ["cavity_id", "center_x_A", "center_y_A", "center_z_A", "max_radius_A", "volume_A3",
              "grid_points", "score", "nearest_residues"]
    rows = [{"cavity_id": c.id, "center_x_A": c.center[0], "center_y_A": c.center[1], "center_z_A": c.center[2],
             "max_radius_A": c.radius, "volume_A3": c.volume, "grid_points": c.grid_points, "score": c.score,
             "nearest_residues": ";".join(c.nearest_residues)} for c in cavities]
    _write_rows(path, rows, fields)


def write_tunnels_csv(tunnels: Sequence[Tunnel], path: str | Path) -> None:
    """Write one row per tunnel.

    Columns: ``tunnel_id``, ``start_{x,y,z}_A``, ``end_{x,y,z}_A``,
    ``length_A`` (polyline length), ``bottleneck_radius_A`` (minimum clearance
    over the connecting segments), ``mean_radius_A``, ``score``,
    ``point_count``, ``grid_spacing_A`` and ``nearest_residues``
    (semicolon-separated). Per-node samples are in :func:`write_tunnel_points_csv`.

    Parameters
    ----------
    tunnels : sequence of Tunnel
        From :func:`crevice.tunnels.find_tunnels`.
    path : str or Path
        Output CSV.
    """
    fields = ["tunnel_id", "start_x_A", "start_y_A", "start_z_A", "end_x_A", "end_y_A", "end_z_A",
              "length_A", "bottleneck_radius_A", "mean_radius_A", "score", "point_count", "grid_spacing_A",
              "nearest_residues"]
    rows = [{"tunnel_id": t.id, "start_x_A": t.start[0], "start_y_A": t.start[1], "start_z_A": t.start[2],
             "end_x_A": t.end[0], "end_y_A": t.end[1], "end_z_A": t.end[2], "length_A": t.length,
             "bottleneck_radius_A": t.bottleneck_radius, "mean_radius_A": t.mean_radius, "score": t.score,
             "point_count": len(t.points), "grid_spacing_A": t.metadata.get("spacing_A", ""),
             "nearest_residues": ";".join(t.nearest_residues)} for t in tunnels]
    _write_rows(path, rows, fields)


def write_tunnel_points_csv(tunnels: Sequence[Tunnel], path: str | Path) -> None:
    """Write one row per tunnel path node.

    Columns: ``tunnel_id`` followed by the profile point columns of
    :func:`write_profile_csv`; ``index`` is the step number and
    ``axis_position_A`` the cumulative path length from the start (Å).

    Parameters
    ----------
    tunnels : sequence of Tunnel
        From :func:`crevice.tunnels.find_tunnels`.
    path : str or Path
        Output CSV.
    """
    fields = ["tunnel_id", "index", "x_A", "y_A", "z_A", "radius_A", "raw_clearance_A", "axis_position_A",
              "nearest_atom_serial", "nearest_residue"]
    rows = [{"tunnel_id": t.id, **_point_row(point)} for t in tunnels for point in t.points]
    _write_rows(path, rows, fields)


def write_network_csv(network: ResidueInteractionNetwork, nodes_path: str | Path, edges_path: str | Path) -> None:
    """Write the network as a node table and an edge table.

    ``nodes_path`` columns: ``node_id``, ``kind`` (``residue`` or ``region``),
    ``label``, ``residue``, ``resname``, ``chain_id``, ``resid``, ``icode``,
    ``group``, ``properties`` (semicolon-separated), ``x_A``/``y_A``/``z_A``
    (residue centroid, Å), ``degree`` (all incident edges, including a region
    edge) and ``residue_degree`` (residue-residue edges only).

    ``edges_path`` columns: ``source``, ``target``, ``interaction``,
    ``distance_A`` (under the network's metric; see ``distance_metric``),
    ``distance_metric``, ``cutoff_A``, ``weight``, ``closest_atoms``,
    ``center_distance_A``, ``surface_distance_A``, ``charged_atom_pairs``,
    ``sequence_separation_observed``, ``different_chains``,
    ``nonlocal_contact`` and ``point_count`` (region edges). Empty cells mean
    the value does not apply to that edge.

    Parameters
    ----------
    network : ResidueInteractionNetwork
        Residue (and optional region) network.
    nodes_path, edges_path : str or Path
        Output CSVs.
    """
    degree = network.degree()
    regions = {node.id for node in network.nodes if node.kind == "region"}
    residue_degree = {node.id: 0 for node in network.nodes}
    for edge in network.edges:
        if edge.source in regions or edge.target in regions:
            continue
        residue_degree[edge.source] = residue_degree.get(edge.source, 0) + 1
        residue_degree[edge.target] = residue_degree.get(edge.target, 0) + 1
    node_fields = ["node_id", "kind", "label", "residue", "resname", "chain_id", "resid", "icode", "group",
                   "properties", "x_A", "y_A", "z_A", "degree", "residue_degree"]
    node_rows = []
    for node in network.nodes:
        meta = node.metadata or {}
        position = meta.get("position") or (None, None, None)
        node_rows.append({"node_id": node.id, "kind": node.kind, "label": node.label, "residue": node.residue or "",
                          "resname": meta.get("resname", ""), "chain_id": meta.get("chain_id", ""),
                          "resid": meta.get("resid", ""), "icode": meta.get("icode", ""),
                          "group": meta.get("group") or "", "properties": ";".join(node.properties),
                          "x_A": position[0], "y_A": position[1], "z_A": position[2],
                          "degree": degree.get(node.id, 0), "residue_degree": residue_degree.get(node.id, 0)})
    _write_rows(nodes_path, node_rows, node_fields)
    edge_fields = ["source", "target", "interaction", "distance_A", "distance_metric", "cutoff_A", "weight",
                   "closest_atoms", "center_distance_A", "surface_distance_A", "charged_atom_pairs",
                   "sequence_separation_observed", "different_chains", "nonlocal_contact", "point_count"]
    metric = network.metadata.get("distance_metric", "")
    edge_rows = []
    for edge in network.edges:
        meta = edge.metadata or {}
        pair = meta.get("closest_pair") or {}
        region = edge.source in regions or edge.target in regions
        edge_rows.append({
            "source": edge.source, "target": edge.target, "interaction": edge.interaction,
            "distance_A": edge.distance,
            "distance_metric": "region_surface_gap" if region else metric,
            "cutoff_A": network.metadata.get("cutoff", ""), "weight": edge.weight,
            "closest_atoms": "-".join(pair.get("atoms", [])),
            "center_distance_A": pair.get("center_distance", ""), "surface_distance_A": pair.get("surface_distance", ""),
            "charged_atom_pairs": ";".join("-".join(str(v) for v in item[:2]) + f"@{item[2]:.3f}"
                                           for item in meta.get("charged_atom_pairs", [])),
            "sequence_separation_observed": meta.get("sequence_separation_observed", ""),
            "different_chains": meta.get("different_chains", ""), "nonlocal_contact": meta.get("nonlocal_contact", ""),
            "point_count": meta.get("point_count", "")})
    _write_rows(edges_path, edge_rows, edge_fields)


def write_connectivity_csv(report: dict[str, Any], residues_path: str | Path, groups_path: str | Path) -> None:
    """Write a :func:`crevice.connectivity.lining_connectivity` report as two tables.

    ``residues_path``: one row per residue with ``residue``, ``group``,
    ``is_lining``, ``removed``, ``reachable``, ``cost_to_lining``
    (summed ``1/weight`` along the cheapest path; dimensionless),
    ``hops_to_lining``, ``path`` (residue IDs joined by ``>``), ``path_hops``,
    ``lining_occupancy``, ``bridge_endpoint`` and ``articulation_point``.
    ``groups_path``: one row per pair of residue groups in contact with
    ``source``, ``target``, ``contact_count`` and ``total_strength``.

    Parameters
    ----------
    report : dict
        Result of :func:`crevice.connectivity.lining_connectivity`.
    residues_path, groups_path : str or Path
        Output CSVs.
    """
    fields = ["residue", "group", "is_lining", "removed", "reachable", "cost_to_lining", "hops_to_lining",
              "path", "path_hops", "lining_occupancy", "bridge_endpoint", "articulation_point"]
    rows = [{**row, "path": ">".join(row.get("path") or [])} for row in report.get("residues", [])]
    _write_rows(residues_path, rows, fields)
    _write_rows(groups_path, report.get("group_contacts", []), ["source", "target", "contact_count", "total_strength"])


FEATURE_UNITS = {"min_radius": "Å", "max_radius": "Å", "mean_radius": "Å", "length": "Å",
                 "volume_estimate": "Å³", "bottleneck_position": "Å", "radius_variance": "Å²"}


def write_features_csv(features: dict[str, Any], path: str | Path) -> None:
    """Write numeric features as ``feature,value,unit`` rows.

    ``unit`` is ``Å``, ``Å²`` or ``Å³`` for the geometric features of
    :func:`crevice.ml.profile_features` and empty for counts, fractions and the
    dimensionless mean influence score.

    Parameters
    ----------
    features : dict of str to float
        For example :func:`crevice.ml.profile_features`.
    path : str or Path
        Output CSV.
    """
    rows = [{"feature": key, "value": value, "unit": FEATURE_UNITS.get(key, "")} for key, value in features.items()]
    _write_rows(path, rows, ["feature", "value", "unit"])


#: Column order of :func:`summary_statistics_rows` and
#: :func:`write_summary_statistics_csv`.
SUMMARY_STATISTICS_FIELDS = (
    "quantity", "item", "value", "total_frames", "observed_frames", "missing_frames", "coverage_fraction",
    "mean", "sample_sd", "median", "quantile_025", "quantile_975", "minimum", "maximum",
    "first_half_mean", "second_half_mean", "statistical_inefficiency", "effective_frames",
    "interval_status", "interval_method", "interval_confidence", "mean_lower", "mean_upper",
    "simultaneous_lower", "simultaneous_upper", "block_length_frames", "replicates",
)
_SERIES_FIELDS = SUMMARY_STATISTICS_FIELDS[3:18]


def _is_series_statistic(value: Any) -> bool:
    return isinstance(value, dict) and "mean" in value and "confidence_interval" in value


def _item_label(entry: Any, index: int) -> str:
    if isinstance(entry, dict):
        if entry.get("residue") is not None:
            return str(entry["residue"])
        if entry.get("source") is not None and entry.get("target") is not None:
            return f"{entry['source']}|{entry['target']}"
        if entry.get("stratum") is not None:
            return str(entry["stratum"])
    return str(index)


def _single(values: Any) -> Any:
    """The value of a one-element interval list; empty for scalars absent or vectors."""
    if isinstance(values, list) and len(values) == 1:
        return values[0]
    return None


def summary_statistics_rows(report: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten the nested interval statistics of a summary report into table rows.

    The trajectory summaries (``*_cavity_statistics.json``,
    ``*_water_membership.json``, ``*_hydration.json`` and ``comparison.json``)
    describe each quantity with the same nested record: frame counts, mean,
    SD, median, 2.5/97.5% frame quantiles, extremes, half-trajectory means,
    statistical inefficiency and a nested ``confidence_interval`` for the mean.
    This function walks ``report`` depth-first in key order and returns one
    flat row per such record, plus one row per numeric value whose key (or an
    enclosing key) contains ``correlation``.

    Parameters
    ----------
    report : dict
        A summary report as written to one of the JSON files above.

    Returns
    -------
    list of dict
        Rows with the keys of :data:`SUMMARY_STATISTICS_FIELDS`:

        ``quantity``
            Dotted JSON path of the record with list positions dropped, for
            example ``volume_A3``, ``residues.boundary_area_A2`` or
            ``statistics.cavity_member_count``. The unit is the key's suffix
            (``_A3`` = Å³, ``_A2`` = Å², ``_A`` = Å; counts and fractions
            have none).
        ``item``
            For a record inside a list, the entry's residue (``residue``), pair
            (``source|target``), stratum, or else its position; empty otherwise.
        ``value``
            Only for correlation rows: the Pearson coefficient (empty
            when it was not estimable). All other columns are then empty.
        ``total_frames`` ... ``effective_frames``
            The record's own fields, unchanged (empty when absent or null).
        ``interval_status``, ``interval_method``, ``interval_confidence``
            ``confidence_interval`` status (for example ``complete``,
            ``insufficient_independent_blocks`` or
            ``unavailable_observed_constant``), method and nominal level.
        ``mean_lower``, ``mean_upper``, ``simultaneous_lower``, ``simultaneous_upper``
            Pointwise and simultaneous bounds on the mean; empty when the
            interval was withheld. Bounds already clipped in the JSON (for
            example occupancies to [0, 1]) are copied as clipped.
        ``block_length_frames``, ``replicates``
            Bootstrap settings, when an interval was attempted.

    Examples
    --------
    >>> from crevice.io import summary_statistics_rows
    >>> report = {"volume_A3": {"mean": 10.0, "total_frames": 3,
    ...           "confidence_interval": {"status": "unavailable_missing_frames"}},
    ...           "residues": [{"residue": "A:ALA1",
    ...                         "area_A2": {"mean": 2.0, "confidence_interval": {"status": "complete",
    ...                                     "pointwise_lower": [1.5], "pointwise_upper": [2.5]}},
    ...                         "volume_correlation": 0.5}]}
    >>> [(r["quantity"], r["item"], r["mean"], r["value"], r["mean_lower"], r["interval_status"])
    ...  for r in summary_statistics_rows(report)]  # doctest: +NORMALIZE_WHITESPACE
    [('volume_A3', '', 10.0, '', '', 'unavailable_missing_frames'),
     ('residues.area_A2', 'A:ALA1', 2.0, '', 1.5, 'complete'),
     ('residues.volume_correlation', 'A:ALA1', '', 0.5, '', '')]
    """
    rows: list[dict[str, Any]] = []

    def blank(value: Any) -> Any:
        return "" if value is None else value

    def visit(node: Any, path: tuple[str, ...], item: str) -> None:
        if _is_series_statistic(node):
            interval = node.get("confidence_interval") or {}
            row = {"quantity": ".".join(path), "item": item, "value": ""}
            row.update({key: blank(node.get(key)) for key in _SERIES_FIELDS})
            row.update(interval_status=blank(interval.get("status")), interval_method=blank(interval.get("method")),
                       interval_confidence=blank(interval.get("confidence")),
                       mean_lower=blank(_single(interval.get("pointwise_lower"))),
                       mean_upper=blank(_single(interval.get("pointwise_upper"))),
                       simultaneous_lower=blank(_single(interval.get("simultaneous_lower"))),
                       simultaneous_upper=blank(_single(interval.get("simultaneous_upper"))),
                       block_length_frames=blank(interval.get("block_length_frames")),
                       replicates=blank(interval.get("replicates")))
            rows.append(row)
        elif isinstance(node, dict):
            for key, value in node.items():
                named = any("correlation" in part for part in path + (str(key),))
                if named and (value is None or isinstance(value, (int, float))) and not isinstance(value, bool):
                    rows.append({**{field: "" for field in SUMMARY_STATISTICS_FIELDS},
                                 "quantity": ".".join(path + (str(key),)), "item": item, "value": blank(value)})
                else:
                    visit(value, path + (str(key),), item)
        elif isinstance(node, list):
            for index, entry in enumerate(node):
                visit(entry, path, _item_label(entry, index) if not item else f"{item}/{_item_label(entry, index)}")

    visit(report, (), "")
    return rows


def write_summary_statistics_csv(report: dict[str, Any], path: str | Path) -> list[dict[str, Any]]:
    """Write the nested interval statistics of a summary report as one CSV table.

    The JSON report is unchanged and remains the complete record (settings,
    definitions, limitations, full per-coordinate interval vectors); this table
    repeats its per-quantity statistics in a flat form for spreadsheets and
    data frames. See :func:`summary_statistics_rows` for the columns.

    Parameters
    ----------
    report : dict
        A summary report (``*_cavity_statistics.json``,
        ``*_water_membership.json``, ``*_hydration.json`` or
        ``comparison.json`` content).
    path : str or Path
        Output CSV. A header-only file is written when the report holds no
        interval statistics.

    Returns
    -------
    list of dict
        The rows written.
    """
    rows = summary_statistics_rows(report)
    _write_rows(path, rows, SUMMARY_STATISTICS_FIELDS)
    return rows


def write_trajectory_csv(analysis: TrajectoryAnalysis, frames_path: str | Path, profiles_path: str | Path) -> None:
    """Write per-frame summaries and all per-frame profile samples of a trajectory analysis.

    ``frames_path``: one row per analysed frame with ``frame_index``,
    ``source_index``, ``time_ps``, ``profile_status``, ``profile_reason``,
    ``min_radius_A``, ``mean_radius_A``, ``max_radius_A``,
    ``bottleneck_axis_position_A``, ``bottleneck_residue``,
    ``profile_point_count``, ``enclosure_radius_A`` (the effective enclosure
    probe of that frame's profile: the pinned probe, or the frame's own
    automatic choice for a fallback frame), ``enclosure_probe_source``
    (``pinned``, ``fallback``, ``none`` when no probe resolved the frame,
    ``not_attempted``, or empty when no profile was requested; see
    :func:`crevice.trajectory.analyze_trajectory`),
    ``residue_contact_count``, ``cavity_count``,
    ``tunnel_count``, ``network_node_count``, ``network_edge_count``,
    ``alignment_status``, ``alignment_rmsd_before_A``,
    ``alignment_rmsd_after_A`` and one ``feature_<name>`` column per numeric
    feature (unit appended for the geometric features: ``_A``, ``_A2``, ``_A3``).

    ``profiles_path``: one row per profile sample of every frame with a
    profile: ``frame_index``, ``time_ps`` and the point columns of
    :func:`write_profile_csv`. The header is written even when no frame has a profile.

    Parameters
    ----------
    analysis : TrajectoryAnalysis
        From :func:`crevice.trajectory.analyze_trajectory`.
    frames_path, profiles_path : str or Path
        Output CSVs.
    """
    suffix = {"Å": "_A", "Å²": "_A2", "Å³": "_A3"}
    frame_rows, profile_rows = [], []
    for frame in analysis.frames:
        meta = frame.metadata or {}
        alignment = meta.get("alignment") or {}
        profile = frame.profile
        row = {"frame_index": frame.frame_index, "source_index": meta.get("source_index", ""),
               "time_ps": frame.time if frame.time is not None else "",
               "profile_status": meta.get("profile_status", "resolved" if profile else ""),
               "profile_reason": meta.get("profile_reason", ""),
               "min_radius_A": profile.min_radius if profile else "",
               "mean_radius_A": profile.mean_radius if profile else "",
               "max_radius_A": profile.max_radius if profile else "",
               "bottleneck_axis_position_A": profile.bottleneck.t if profile else "",
               "bottleneck_residue": (profile.bottleneck.nearest_residue or "") if profile else "",
               "profile_point_count": len(profile.points) if profile else 0,
               "enclosure_radius_A": profile.metadata.get("enclosure_radius_A", "") if profile else "",
               "enclosure_probe_source": (meta.get("enclosure_probe") or {}).get("source", ""),
               "residue_contact_count": len(frame.residue_contacts), "cavity_count": len(frame.cavities),
               "tunnel_count": len(frame.tunnels),
               "network_node_count": len(frame.network.nodes) if frame.network else meta.get("network_node_count", ""),
               "network_edge_count": len(frame.network.edges) if frame.network else meta.get("network_edge_count", ""),
               "alignment_status": alignment.get("status", ""),
               "alignment_rmsd_before_A": alignment.get("rmsd_before_A", ""),
               "alignment_rmsd_after_A": alignment.get("rmsd_after_A", "")}
        for key, value in (frame.features or {}).items():
            row[f"feature_{key}{suffix.get(FEATURE_UNITS.get(key, ''), '')}"] = value
        frame_rows.append(row)
        if profile:
            profile_rows.extend({"frame_index": frame.frame_index, "time_ps": row["time_ps"], **_point_row(point)}
                                for point in profile.points)
    _write_rows(frames_path, frame_rows)
    _write_rows(profiles_path, profile_rows, ["frame_index", "time_ps", "index", "x_A", "y_A", "z_A", "radius_A",
                                              "raw_clearance_A", "axis_position_A", "nearest_atom_serial",
                                              "nearest_residue"])


def _write_rows(path: str | Path, rows: Sequence[dict[str, Any]], fields: Sequence[str] | None = None) -> None:
    """Write dictionaries as CSV; the field order is ``fields`` or first-seen key order."""
    _ensure_parent(path)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        writer.writerows(rows)


def write_points_pdb(points: Iterable[ChannelPoint], path: str | Path, *, resname: str = "DUM") -> None:
    """Write points as dummy ``HETATM`` atoms with the radius in the B-factor column.

    Each point becomes atom ``X`` of residue ``resname``, chain ``A``, with atom
    serial and residue number both equal to its position in ``points`` (from 1),
    occupancy 1.00, B-factor = ``point.radius`` clipped to 0-99.99 Å and element
    ``X``. The file ends with ``END``.

    Parameters
    ----------
    points : iterable of ChannelPoint
        Points to write; only ``position`` and ``radius`` are used.
    path : str or Path
        Output PDB path.
    resname : str, default "DUM"
        Three-character residue name identifying the point set (``PRF``
        profile, ``CAV`` cavities, ``TUN`` tunnels, ``VOI`` cast samples).

    Examples
    --------
    >>> import tempfile, pathlib
    >>> from crevice.models import ChannelPoint
    >>> from crevice.io import write_points_pdb
    >>> path = pathlib.Path(tempfile.mkdtemp()) / "points.pdb"
    >>> write_points_pdb([ChannelPoint(index=0, position=(1.0, 2.0, 3.0), radius=1.5,
    ...                                raw_clearance=1.5, t=0.0)], path, resname="PRF")
    >>> print(path.read_text().splitlines()[0].rstrip())
    HETATM    1  X   PRF A   1       1.000   2.000   3.000  1.00  1.50          X
    """
    _ensure_parent(path)
    lines: list[str] = []
    for serial, point in enumerate(points, start=1):
        x, y, z = point.position
        radius = max(0.0, min(point.radius, 99.99))
        lines.append(
            f"HETATM{serial:5d}  X   {resname:>3s} A{serial:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{radius:6.2f}          X "
        )
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def export_pymol_script(
    *,
    structure_path: str | Path,
    dummy_path: str | Path,
    script_path: str | Path,
    object_name: str = "crevice_pore_cast",
    output_png: str | Path | None = None,
    width: int = 2400,
    height: int = 1800,
    dpi: int = 400,
) -> None:
    """Write a PyMOL script showing profile dummy atoms inside the protein.

    Used by ``crevice profile --pymol``. The dummy-atom PDB (one atom per
    profile sample, B-factor = radius in Å) is shown as translucent spheres of
    that radius coloured by radius, with a surface of the same spheres, inside
    the protein cartoon. Files are read from the working directory when
    present there, else from their absolute paths. No text is drawn.

    Parameters
    ----------
    structure_path : str or Path
        Protein structure.
    dummy_path : str or Path
        Profile dummy-atom PDB (``crevice profile --pdb``).
    script_path : str or Path
        Output ``.pml`` path.
    object_name : str, default "crevice_pore_cast"
        PyMOL object name of the dummy atoms.
    output_png : str or Path, optional
        When given, the script ray-traces this PNG.
    width, height : int, default 2400, 1800
        Rendered image size in pixels.
    dpi : int, default 400
        DPI recorded in the rendered PNG.

    Outputs
    -------
    script_path : PyMOL script
        Run with ``pymol script.pml``.
    output_png : PNG, optional
        Written by PyMOL when the script runs.

    Colours and representation
    --------------------------
    Cyan cartoon (``[0.00, 0.82, 0.92]``, 12% transparent); dummy-atom spheres
    55% transparent on PyMOL's ``blue_white_red`` spectrum by radius (narrow
    blue, wide red); pink (``[1.00, 0.05, 0.55]``) surface of the spheres,
    18% transparent; white background.
    """

    _ensure_parent(script_path)
    structure_abs = Path(structure_path).resolve()
    dummy_abs = Path(dummy_path).resolve()
    structure_local = structure_abs.name
    dummy_local = dummy_abs.name
    png_block = ""
    if output_png is not None:
        png_abs = Path(output_png).resolve()
        png_local = png_abs.name
        png_block = f"""python
_required = [{structure_local!r}, {dummy_local!r}]
_png_target = _crevice_output({png_local!r}, {str(png_abs)!r}, _required)
cmd.ray({width}, {height})
cmd.png(_png_target, dpi={dpi})
python end
"""
    script = f"""reinitialize
set retain_order, 1
bg_color white
set ray_opaque_background, on
set antialias, 4
set ray_trace_mode, 1
set ray_shadows, off
set ambient, 0.36
set direct, 0.72
set specular, 0.18
set shininess, 24
set depth_cue, 0
set cartoon_fancy_helices, on
set cartoon_sampling, 14
set_color crevice_cyan, [0.00, 0.82, 0.92]
set_color crevice_pink, [1.00, 0.05, 0.55]
set_color crevice_deep_blue, [0.08, 0.18, 0.55]
python
from pathlib import Path
from pymol import cmd

def _crevice_input(local_name, fallback):
    local = Path.cwd() / local_name
    return str(local if local.exists() else Path(fallback))

def _crevice_output(local_name, fallback, required_local_inputs):
    if all((Path.cwd() / name).exists() for name in required_local_inputs):
        return str(Path.cwd() / local_name)
    return str(Path(fallback))

cmd.load(_crevice_input({structure_local!r}, {str(structure_abs)!r}), 'crevice_cartoon')
cmd.load(_crevice_input({dummy_local!r}, {str(dummy_abs)!r}), {object_name!r})
python end
hide everything, crevice_cartoon
show cartoon, crevice_cartoon
color crevice_cyan, crevice_cartoon
set cartoon_transparency, 0.12, crevice_cartoon
hide everything, {object_name}
alter {object_name}, vdw=max(0.25, b)
rebuild
create crevice_pore_fill, {object_name}
hide everything, crevice_pore_fill
show surface, crevice_pore_fill
color crevice_pink, crevice_pore_fill
set surface_color, crevice_pink, crevice_pore_fill
set surface_quality, 2, crevice_pore_fill
set transparency, 0.18, crevice_pore_fill
show spheres, {object_name}
set sphere_quality, 4, {object_name}
set sphere_transparency, 0.55, {object_name}
spectrum b, blue_white_red, {object_name}
orient crevice_cartoon or {object_name} or crevice_pore_fill
zoom crevice_cartoon or {object_name} or crevice_pore_fill, 6
{png_block}"""
    Path(script_path).write_text(script)

def _write_json(data: object, path: str | Path) -> None:
    _ensure_parent(path)
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def _ensure_parent(path: str | Path) -> None:
    parent = Path(path).parent
    if parent != Path("."):
        parent.mkdir(parents=True, exist_ok=True)

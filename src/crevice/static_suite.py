"""Static benchmark suite for CREVICE validation systems."""

from __future__ import annotations

import csv
import json
import shutil
import time
from pathlib import Path
from typing import Sequence

from .analysis import annotate_residues, detect_cavities, pore_profile
from .tunnels import find_tunnels
from .voids import validate_void_cast, void_cast
from .benchmarks import DEFAULT_BENCHMARK_IDS, benchmark_path, fetch_benchmark_structure, get_benchmark_system
from .figures import (
    plot_cavity_summary,
    plot_network_chord,
    plot_network_summary,
    plot_profile_radius,
    plot_profile_radius_with_residues,
    plot_residue_contacts,
    write_publication_pymol_script,
)
from .io import (
    write_cavities_json,
    write_cavities_pdb,
    write_profile_csv,
    write_profile_json,
    write_profile_pdb,
    write_residue_contacts_csv,
    write_tunnels_json,
    write_tunnels_pdb,
    write_void_cast_json,
    write_void_cast_pdb,
    write_cavities_csv,
    write_network_csv,
    write_tunnel_points_csv,
    write_tunnels_csv,
    write_void_cast_csv,
)
from .networks import build_cavity_network, network_metrics
from .parser import load_structure
from .presentation import annotate_option
from .radii import RadiusSet, radii_fields, radii_option


@radii_option
@annotate_option
def run_static_benchmark_suite(
    *,
    cache_dir: str | Path = ".crevice/pdb",
    output_dir: str | Path = "results/static-benchmark-suite",
    pdb_ids: Sequence[str] = DEFAULT_BENCHMARK_IDS,
    extension: str = "pdb",
    fetch: bool = False,
    profile_samples: int = 61,
    search_radius: float = 3.0,
    refinement_steps: int = 4,
    probe_radii: Sequence[float] = (0.0, 1.0, 1.4),
    residue_label_every: int = 5,
    cavity_spacing: float = 3.0,
    tunnel_spacing: float = 4.0,
    cast_spacing: float = 1.5,
    cast_min_radius: float = 0.4,
    cast_max_components: int = 4,
    cast_max_grid_points: int = 120_000,
    max_grid_points: int = 60_000,
    include_cavities: bool = True,
    include_tunnels: bool = True,
    include_network: bool = True,
    dpi: int = 400,
    annotate: bool | None = None,
    enclosure_radius: float | None = None,
    radii: RadiusSet | str | None = None,
) -> dict[str, object]:
    """Run the static analyses of ``crevice static-suite`` across benchmark systems.

    For each configured benchmark entry, the structure is loaded (optionally
    fetched), profiled at each probe radius, and its residue contacts, void
    cast, PyMOL scene, cavities, tunnels and network are written with figures.
    Per-system failures are recorded rather than raised. The benchmark axes are
    assumed, not curated: outputs describe geometry at the stated parameters
    and are not validated biological results.

    Parameters
    ----------
    cache_dir : str or Path, default ".crevice/pdb"
        Structure cache.
    output_dir : str or Path, default "results/static-benchmark-suite"
        Suite directory; each system gets a subdirectory.
    pdb_ids : sequence of str
        Benchmark IDs (default :data:`crevice.benchmarks.DEFAULT_BENCHMARK_IDS`).
    extension : {"pdb", "cif", "mmcif"}, default "pdb"
        Structure format.
    fetch : bool, default False
        Download missing structures.
    profile_samples : int, default 61
        Profile samples.
    search_radius : float, default 3.0
        Centre-refinement radius (Å).
    refinement_steps : int, default 4
        Centre-refinement iterations.
    probe_radii : sequence of float, default (0.0, 1.0, 1.4)
        Measurement probe radii (Å) for the sensitivity rows. Each row keeps
        the system's enclosure probe fixed; a row whose probe does not pass
        the channel is recorded as ``unresolved`` with its reason instead of
        failing the system.
    residue_label_every : int, default 5
        Landmark spacing in the residue-landmark profile figure.
    cavity_spacing, tunnel_spacing, cast_spacing : float
        Grid spacings (Å): 3.0, 4.0 and 1.5 by default.
    cast_min_radius : float, default 0.4
        Minimum clearance of cast samples (Å).
    cast_max_components : int, default 4
        Maximum cast components.
    cast_max_grid_points, max_grid_points : int
        Grid allocation guards for the cast and for cavity/tunnel searches.
    include_cavities, include_tunnels, include_network : bool, default True
        Optional analyses.
    dpi : int, default 400
        Figure resolution.
    annotate : bool, optional
        Draw titles, value call-outs and residue landmark names in every figure.
        ``None`` (default) inherits the surrounding setting, which is off unless
        enabled with :func:`crevice.presentation.figure_annotations` or ``--annotate``.
    enclosure_radius : float, optional
        Enclosure probe (Å) of each system's channel profile (see
        :func:`crevice.channels.pore_profile`). ``None`` (default; CLI
        ``--enclosure-radius auto``) chooses it per system with the automatic
        rule (smallest probe of the most persistent channel, 0.5-3.0 Å) and
        uses the chosen value for that system's sensitivity rows. A number
        runs that probe everywhere (``0.8`` reproduces the suite outputs
        recorded before the automatic default). The choice is recorded in
        each system summary (``enclosure_probe``) and the summary CSV
        (``enclosure_radius_A``, ``enclosure_probe_mode``).
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict
        Suite summary with per-system rows, ``summary_json`` and ``summary_csv``.

    Outputs
    -------
    static_benchmark_summary.csv : table
        One row per system. Columns carry units (``profile_min_radius_A``,
        ``void_cast_total_volume_A3``, ``top_tunnel_length_A``, ...); ``error``
        holds the failure message of a system that did not complete.
    static_benchmark_summary.json : JSON
        Kept as the machine-readable suite record: parameters, per-system
        summaries with nested timings, output validation and top residues, and
        each system's manifest path.
    <ID>/<ID>_manifest.json, <ID>/<ID>_system_summary.json : JSON
        File index with the cast record; the system's summary record.
    <ID>/<ID>_profile.csv, <ID>/<ID>_profile.json : profile
        Per-sample table (see :func:`crevice.io.write_profile_csv`) and the full
        profile record (metadata, mouths, settings).
    <ID>/<ID>_profile_sensitivity.csv : table
        One row per probe radius: ``probe_radius_A``, ``enclosure_radius_A``
        (the system's enclosure probe, fixed across rows), ``status``
        (``resolved`` or ``unresolved``), ``reason`` (why an unresolved row
        failed), ``min_radius_A``, ``mean_radius_A``, ``max_radius_A``,
        ``bottleneck_index``, ``bottleneck_t_A``, ``bottleneck_residue``,
        ``length_A``, ``volume_estimate_A3``, ``seconds``.
    <ID>/<ID>_residue_contacts.csv : table
        :func:`crevice.io.write_residue_contacts_csv`.
    <ID>/<ID>_void_cast.csv, <ID>/<ID>_void_cast.json, <ID>/<ID>_pore_cast.pdb : cast
        Per-region table, full cast record, dummy atoms.
    <ID>/<ID>_cavities.{csv,json,pdb}, <ID>/<ID>_tunnels.{csv,json,pdb}, <ID>/<ID>_tunnel_points.csv : tables
        Cavity and tunnel tables (see :mod:`crevice.io`), JSON records and PDB dummy atoms.
    <ID>/<ID>_network_nodes.csv, <ID>/<ID>_network_edges.csv, <ID>/<ID>_network.json : network
        Node and edge tables and the network record with metrics.
    <ID>/... : figures and scenes
        :func:`crevice.figures.plot_profile_radius`,
        :func:`~crevice.figures.plot_profile_radius_with_residues`,
        :func:`~crevice.figures.plot_residue_contacts`,
        :func:`~crevice.figures.plot_cavity_summary`,
        :func:`~crevice.figures.plot_network_summary` and
        :func:`~crevice.figures.plot_network_chord`, the volume viewer scenes and
        a PyMOL scene (:func:`~crevice.figures.write_publication_pymol_script`).
        No tunnel summary figure is written; tunnel data are in the CSV/PDB files.

    Colours and representation
    --------------------------
    As documented for each figure function and the PyMOL scene.
    """

    out_root = Path(output_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    manifests: dict[str, str] = {}
    for pdb_id in pdb_ids:
        system = get_benchmark_system(pdb_id)
        structure_path = benchmark_path(cache_dir, system.pdb_id, extension=extension)
        if not structure_path.exists():
            if not fetch:
                summaries.append(_error_row(
                    system.pdb_id, system.expected_geometry, structure_path,
                    f"missing structure: no {extension} file at {structure_path}; "
                    f"rerun with --fetch, or --format to match the cache",
                    system.curation_status))
                continue
            structure_path = fetch_benchmark_structure(system.pdb_id, cache_dir, extension=extension)
        try:
            summary, manifest_path = _run_system(
                pdb_id=system.pdb_id,
                expected_geometry=system.expected_geometry,
                curation_status=system.curation_status,
                structure_path=structure_path,
                output_dir=out_root / system.pdb_id,
                profile_samples=profile_samples,
                search_radius=search_radius,
                refinement_steps=refinement_steps,
                probe_radii=probe_radii,
                residue_label_every=residue_label_every,
                cavity_spacing=cavity_spacing,
                tunnel_spacing=tunnel_spacing,
                cast_spacing=cast_spacing,
                cast_min_radius=cast_min_radius,
                cast_max_components=cast_max_components,
                cast_max_grid_points=cast_max_grid_points,
                max_grid_points=max_grid_points,
                include_cavities=include_cavities,
                include_tunnels=include_tunnels,
                include_network=include_network,
                dpi=dpi,
                enclosure_radius=enclosure_radius,
            )
        except (OSError, ValueError, RuntimeError) as exc:
            row = _error_row(system.pdb_id, system.expected_geometry, structure_path,
                             f"{type(exc).__name__}: {exc}", system.curation_status)
            selection = getattr(exc, "probe_selection", None)
            if selection is not None:
                row["enclosure_probe"] = {"mode": "auto", "chosen_A": None, "reason": selection["reason"],
                                          "selection": selection}
                row["enclosure_probe_mode"] = "auto"
            summaries.append(row)
            continue
        summaries.append(summary)
        manifests[system.pdb_id] = str(manifest_path)

    summary_json = out_root / "static_benchmark_summary.json"
    summary_csv = out_root / "static_benchmark_summary.csv"
    result = {
        "systems": summaries,
        "manifests": manifests,
        "parameters": {
            "profile_samples": profile_samples,
            "search_radius": search_radius,
            "refinement_steps": refinement_steps,
            "probe_radii": list(probe_radii),
            "cavity_spacing": cavity_spacing,
            "tunnel_spacing": tunnel_spacing,
            "cast_spacing": cast_spacing,
            "cast_min_radius": cast_min_radius,
            "cast_max_components": cast_max_components,
            "cast_max_grid_points": cast_max_grid_points,
            "max_grid_points": max_grid_points,
            "include_cavities": include_cavities,
            "include_tunnels": include_tunnels,
            "include_network": include_network,
            "dpi": dpi,
            "enclosure_radius": "auto" if enclosure_radius is None else enclosure_radius,
        },
        **radii_fields(),
    }
    _write_json(result, summary_json)
    _write_summary_csv(summaries, summary_csv)
    return {**result, "summary_json": str(summary_json), "summary_csv": str(summary_csv)}


def _run_system(
    *,
    pdb_id: str,
    expected_geometry: str,
    curation_status: str,
    structure_path: Path,
    output_dir: Path,
    profile_samples: int,
    search_radius: float,
    refinement_steps: int,
    probe_radii: Sequence[float],
    residue_label_every: int,
    cavity_spacing: float,
    tunnel_spacing: float,
    cast_spacing: float,
    cast_min_radius: float,
    cast_max_components: int,
    cast_max_grid_points: int,
    max_grid_points: int,
    include_cavities: bool,
    include_tunnels: bool,
    include_network: bool,
    dpi: int,
    enclosure_radius: float | None = None,
) -> tuple[dict[str, object], Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, str] = {}
    timings: dict[str, float] = {}
    prefix = pdb_id
    structure_copy = _copy_structure(structure_path, output_dir, prefix)
    manifest["structure_file"] = str(structure_copy)

    started = time.perf_counter()
    frame = load_structure(structure_path)
    timings["load_seconds"] = time.perf_counter() - started

    started = time.perf_counter()
    profile = pore_profile(
        frame,
        axis="auto",
        samples=profile_samples,
        search_radius=search_radius,
        refinement_steps=refinement_steps,
        enclosure_radius=enclosure_radius,
    )
    timings["profile_seconds"] = time.perf_counter() - started
    enclosure_probe = _enclosure_probe_record(profile, enclosure_radius)

    started = time.perf_counter()
    contacts = annotate_residues(frame, profile.points)
    timings["residue_annotation_seconds"] = time.perf_counter() - started

    started = time.perf_counter()
    cast = void_cast(
        frame,
        mode=_cast_mode_for_expected(expected_geometry),
        axis="auto",
        spacing=cast_spacing,
        min_radius=cast_min_radius,
        max_components=cast_max_components,
        max_grid_points=cast_max_grid_points,
        profile=profile,
    )
    timings["void_cast_seconds"] = time.perf_counter() - started
    cast_validation = validate_void_cast(frame, cast)

    manifest.update(
        {
            "profile_json": str(output_dir / f"{prefix}_profile.json"),
            "profile_csv": str(output_dir / f"{prefix}_profile.csv"),
            "pore_cast_pdb": str(output_dir / f"{prefix}_pore_cast.pdb"),
            "pore_cast_pml": str(output_dir / f"{prefix}_pore_cast.pml"),
            "void_cast_json": str(output_dir / f"{prefix}_void_cast.json"),
            "void_cast_csv": str(output_dir / f"{prefix}_void_cast.csv"),
            "pymol_render_png": str(output_dir / f"{prefix}_pore_cast.png"),
            "profile_radius_png": str(output_dir / f"{prefix}_profile_radius.png"),
            "profile_radius_annotated_png": str(output_dir / f"{prefix}_profile_radius_annotated.png"),
            "residue_contacts_csv": str(output_dir / f"{prefix}_residue_contacts.csv"),
            "residue_contacts_png": str(output_dir / f"{prefix}_residue_contacts.png"),
        }
    )
    write_profile_json(profile, manifest["profile_json"])
    write_profile_csv(profile, manifest["profile_csv"])
    write_void_cast_json(cast, manifest["void_cast_json"])
    write_void_cast_csv(cast, manifest["void_cast_csv"])
    write_void_cast_pdb(cast, manifest["pore_cast_pdb"])
    write_residue_contacts_csv(contacts, manifest["residue_contacts_csv"])
    plot_profile_radius(profile, manifest["profile_radius_png"], dpi=dpi)
    plot_profile_radius_with_residues(
        profile,
        manifest["profile_radius_annotated_png"],
        dpi=dpi,
        label_every=residue_label_every,
    )
    plot_residue_contacts(contacts, manifest["residue_contacts_png"], dpi=dpi)
    write_publication_pymol_script(
        structure_path=structure_copy,
        pore_cast_path=manifest["pore_cast_pdb"],
        script_path=manifest["pore_cast_pml"],
        profile=profile,
        output_png=manifest["pymol_render_png"],
        dpi=dpi,
        structure_local_name=Path(manifest["structure_file"]).name,
        pore_cast_local_name=Path(manifest["pore_cast_pdb"]).name,
        output_png_local_name=Path(manifest["pymol_render_png"]).name,
        bottleneck_resi=_pdb_resi_for_min_radius(cast.points),
    )
    if cast.components:
        from .volume_export import write_volume_viewer_bundle
        manifest.update(write_volume_viewer_bundle(cast, structure_path=structure_copy,
                        output_dir=output_dir, prefix=prefix))
        if profile.method == "axial-connected":
            legacy = output_dir / f"{prefix}_pore_spheres.pml"
            shutil.copyfile(manifest["pore_cast_pml"], legacy)
            manifest["pore_spheres_pml"] = str(legacy)
            shutil.copyfile(manifest["volume_pml"], manifest["pore_cast_pml"])


    sensitivity_rows = _profile_sensitivity(
        frame,
        probe_radii=probe_radii,
        samples=profile_samples,
        search_radius=search_radius,
        refinement_steps=refinement_steps,
        enclosure_radius=enclosure_probe["chosen_A"],
    )
    manifest["sensitivity_csv"] = str(output_dir / f"{prefix}_profile_sensitivity.csv")
    _write_rows_csv(sensitivity_rows, manifest["sensitivity_csv"])

    network = None
    network_data: dict[str, object] = {}
    if include_network:
        started = time.perf_counter()
        network = build_cavity_network(frame, profile, contacts=contacts)
        timings["network_seconds"] = time.perf_counter() - started
        network_data = {"network": network.to_dict(), "metrics": network_metrics(network)}
        manifest["network_json"] = str(output_dir / f"{prefix}_network.json")
        manifest["network_nodes_csv"] = str(output_dir / f"{prefix}_network_nodes.csv")
        manifest["network_edges_csv"] = str(output_dir / f"{prefix}_network_edges.csv")
        manifest["network_summary_png"] = str(output_dir / f"{prefix}_network_summary.png")
        manifest["network_chord_png"] = str(output_dir / f"{prefix}_network_chord.png")
        _write_json(network_data, manifest["network_json"])
        write_network_csv(network, manifest["network_nodes_csv"], manifest["network_edges_csv"])
        plot_network_summary(network, manifest["network_summary_png"], dpi=dpi)
        plot_network_chord(network, manifest["network_chord_png"], dpi=dpi)

    cavities = ()
    if include_cavities:
        started = time.perf_counter()
        cavities = detect_cavities(
            frame,
            spacing=cavity_spacing,
            max_cavities=8,
            max_grid_points=max_grid_points,
        )
        timings["cavity_seconds"] = time.perf_counter() - started
        manifest["cavities_json"] = str(output_dir / f"{prefix}_cavities.json")
        manifest["cavities_csv"] = str(output_dir / f"{prefix}_cavities.csv")
        manifest["cavities_pdb"] = str(output_dir / f"{prefix}_cavities.pdb")
        manifest["cavity_summary_png"] = str(output_dir / f"{prefix}_cavity_summary.png")
        write_cavities_json(cavities, manifest["cavities_json"])
        write_cavities_csv(cavities, manifest["cavities_csv"])
        write_cavities_pdb(cavities, manifest["cavities_pdb"])
        plot_cavity_summary(cavities, manifest["cavity_summary_png"], dpi=dpi)

    tunnels = ()
    if include_tunnels:
        started = time.perf_counter()
        tunnels = find_tunnels(
            frame,
            spacing=tunnel_spacing,
            min_radius=1.0,
            max_tunnels=3,
            max_grid_points=max_grid_points,
        )
        timings["tunnel_seconds"] = time.perf_counter() - started
        manifest["tunnels_json"] = str(output_dir / f"{prefix}_tunnels.json")
        manifest["tunnels_csv"] = str(output_dir / f"{prefix}_tunnels.csv")
        manifest["tunnel_points_csv"] = str(output_dir / f"{prefix}_tunnel_points.csv")
        manifest["tunnels_pdb"] = str(output_dir / f"{prefix}_tunnels.pdb")
        write_tunnels_json(tunnels, manifest["tunnels_json"])
        write_tunnels_csv(tunnels, manifest["tunnels_csv"])
        write_tunnel_points_csv(tunnels, manifest["tunnel_points_csv"])
        write_tunnels_pdb(tunnels, manifest["tunnels_pdb"])

    validation = _validate_manifest(manifest)
    summary = {
        "pdb_id": pdb_id,
        "expected_geometry": expected_geometry,
        "curation_status": curation_status,
        "status": "ok" if validation["missing_required"] == [] else "missing_outputs",
        "atom_count": len(frame.atoms),
        "residue_count": len(frame.residues()),
        "profile_min_radius": profile.min_radius,
        "profile_mean_radius": profile.mean_radius,
        "profile_max_radius": profile.max_radius,
        "profile_bottleneck_index": profile.bottleneck.index,
        "profile_bottleneck_t": profile.bottleneck.t,
        "profile_bottleneck_residue": profile.bottleneck.nearest_residue,
        "enclosure_radius": enclosure_probe["chosen_A"],
        "enclosure_probe_mode": enclosure_probe["mode"],
        "enclosure_probe": enclosure_probe,
        "profile_sensitivity_unresolved_probe_radii": [row["probe_radius_A"] for row in sensitivity_rows
                                                       if row["status"] != "resolved"],
        "residue_contact_count": len(contacts),
        "void_cast_mode": cast.mode,
        "void_cast_component_count": cast.component_count,
        "void_cast_point_count": cast.point_count,
        "void_cast_total_volume": cast.total_volume,
        "void_cast_min_clearance": cast.min_clearance,
        "void_cast_mean_clearance": cast.mean_clearance,
        "void_cast_max_clearance": cast.max_clearance,
        "void_cast_component_kinds": [component.kind for component in cast.components],
        "void_cast_validation": cast_validation,
        "top_residues": [contact.to_dict() for contact in contacts[:10]],
        "network_node_count": len(network.nodes) if network else 0,
        "network_edge_count": len(network.edges) if network else 0,
        "network_density": network_metrics(network)["density"] if network else 0.0,
        "cavity_count": len(cavities),
        "top_cavity_volume": cavities[0].volume if cavities else 0.0,
        "top_cavity_radius": cavities[0].radius if cavities else 0.0,
        "tunnel_count": len(tunnels),
        "top_tunnel_bottleneck_radius": tunnels[0].bottleneck_radius if tunnels else 0.0,
        "top_tunnel_length": tunnels[0].length if tunnels else 0.0,
        "timings": timings,
        "output_validation": validation,
    }
    manifest_path = output_dir / f"{prefix}_manifest.json"
    _write_json({"files": manifest, "summary": summary, "void_cast": cast.to_dict(), **radii_fields()}, manifest_path)
    manifest["manifest_json"] = str(manifest_path)
    _write_json(summary, output_dir / f"{prefix}_system_summary.json")
    return summary, manifest_path



def _cast_mode_for_expected(expected_geometry: str) -> str:
    normalized = expected_geometry.lower()
    if normalized in {"channel", "cavity"}:
        return normalized
    return "auto"


def _pdb_resi_for_min_radius(points: Sequence) -> int | None:
    if not points:
        return None
    return min(enumerate(points, start=1), key=lambda item: item[1].radius)[0]

def _enclosure_probe_record(profile, enclosure_radius: float | None) -> dict[str, object]:
    """Which enclosure probe a suite profile used and why (``mode``, ``chosen_A``, ``reason``)."""
    if profile.method != "axial-connected":
        return {"mode": "not_used", "chosen_A": None,
                "reason": "search_radius=0 fixed-axis scan does not use an enclosure probe"}
    if enclosure_radius is not None:
        return {"mode": "explicit", "chosen_A": enclosure_radius, "reason": "enclosure_radius given explicitly"}
    selection = profile.metadata["enclosure_probe_selection"]
    return {"mode": "auto", "chosen_A": selection["chosen_A"], "reason": selection["reason"],
            "selection": selection}


def _profile_sensitivity(
    frame,
    *,
    probe_radii: Sequence[float],
    samples: int,
    search_radius: float,
    refinement_steps: int,
    enclosure_radius: float | None = 0.8,
) -> list[dict[str, object]]:
    """One row per measurement probe, all at the same (explicit) enclosure probe.

    ``enclosure_radius`` is ``None`` only for ``search_radius=0`` scans, which
    do not use it.

    A probe that does not pass the channel gives an ``unresolved`` row with
    the resolver's reason instead of an exception.
    """
    from .sections import ChannelResolutionError
    rows = []
    for probe_radius in probe_radii:
        started = time.perf_counter()
        try:
            profile = pore_profile(
                frame,
                axis="auto",
                samples=samples,
                search_radius=search_radius,
                refinement_steps=refinement_steps,
                probe_radius=probe_radius,
                enclosure_radius=enclosure_radius,
            )
        except ChannelResolutionError as exc:
            rows.append({"probe_radius_A": probe_radius, "enclosure_radius_A": enclosure_radius,
                         "status": "unresolved", "reason": str(exc).split(". ")[0],
                         "seconds": time.perf_counter() - started})
            continue
        rows.append(
            {
                "probe_radius_A": probe_radius,
                "enclosure_radius_A": enclosure_radius,
                "status": "resolved",
                "reason": "",
                "min_radius_A": profile.min_radius,
                "mean_radius_A": profile.mean_radius,
                "max_radius_A": profile.max_radius,
                "bottleneck_index": profile.bottleneck.index,
                "bottleneck_t_A": profile.bottleneck.t,
                "bottleneck_residue": profile.bottleneck.nearest_residue,
                "length_A": profile.length,
                "volume_estimate_A3": profile.volume_estimate,
                "seconds": time.perf_counter() - started,
            }
        )
    return rows


def _copy_structure(source: Path, output_dir: Path, prefix: str) -> Path:
    suffix = "".join(source.suffixes) or ".pdb"
    destination = output_dir / f"{prefix}_structure{suffix}"
    if source.resolve() != destination.resolve():
        shutil.copyfile(source, destination)
    return destination


def _validate_manifest(manifest: dict[str, str]) -> dict[str, object]:
    optional_missing_allowed = {"pymol_render_png"}
    existing = []
    missing_required = []
    missing_optional = []
    for key, value in manifest.items():
        exists = Path(value).exists()
        if exists:
            existing.append(key)
        elif key in optional_missing_allowed:
            missing_optional.append(key)
        else:
            missing_required.append(key)
    return {
        "existing": sorted(existing),
        "missing_required": sorted(missing_required),
        "missing_optional": sorted(missing_optional),
    }


def _error_row(pdb_id: str, expected_geometry: str, path: Path, error: str,
               curation_status: str = "uncurated") -> dict[str, object]:
    return {
        "pdb_id": pdb_id,
        "expected_geometry": expected_geometry,
        "curation_status": curation_status,
        "status": "error",
        "path": str(path),
        "error": error,
    }


# Summary-row key -> CSV column (the unit is part of the column name: _A = Å, _A3 = Å³).
SUMMARY_COLUMNS = {
    "pdb_id": "pdb_id", "expected_geometry": "expected_geometry", "curation_status": "curation_status",
    "status": "status", "atom_count": "atom_count", "residue_count": "residue_count",
    "profile_min_radius": "profile_min_radius_A", "profile_mean_radius": "profile_mean_radius_A",
    "profile_max_radius": "profile_max_radius_A", "profile_bottleneck_index": "profile_bottleneck_index",
    "profile_bottleneck_t": "profile_bottleneck_t_A", "profile_bottleneck_residue": "profile_bottleneck_residue",
    "enclosure_radius": "enclosure_radius_A", "enclosure_probe_mode": "enclosure_probe_mode",
    "residue_contact_count": "residue_contact_count", "void_cast_mode": "void_cast_mode",
    "void_cast_component_count": "void_cast_component_count", "void_cast_point_count": "void_cast_point_count",
    "void_cast_total_volume": "void_cast_total_volume_A3", "void_cast_min_clearance": "void_cast_min_clearance_A",
    "void_cast_mean_clearance": "void_cast_mean_clearance_A", "void_cast_max_clearance": "void_cast_max_clearance_A",
    "network_node_count": "network_node_count", "network_edge_count": "network_edge_count",
    "network_density": "network_density", "cavity_count": "cavity_count",
    "top_cavity_volume": "top_cavity_volume_A3", "top_cavity_radius": "top_cavity_radius_A",
    "tunnel_count": "tunnel_count", "top_tunnel_bottleneck_radius": "top_tunnel_bottleneck_radius_A",
    "top_tunnel_length": "top_tunnel_length_A", "total_seconds": "total_seconds", "error": "error",
}


def _write_summary_csv(rows: Sequence[dict[str, object]], path: str | Path) -> None:
    """One row per system; column names carry units (see ``SUMMARY_COLUMNS``); ``error`` holds a failure message."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SUMMARY_COLUMNS.values()))
        writer.writeheader()
        for row in rows:
            csv_row = {column: row.get(key, "") for key, column in SUMMARY_COLUMNS.items()}
            timings = row.get("timings")
            if isinstance(timings, dict):
                csv_row["total_seconds"] = sum(float(value) for value in timings.values())
            writer.writerow(csv_row)


def _write_rows_csv(rows: Sequence[dict[str, object]], path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_json(data: object, path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

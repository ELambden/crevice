"""CREVICE public API.

CREVICE stands for Cavity and Residue Environment Visualisation, Interaction and
Connectivity Evaluation.

CREVICE measures empty space in protein structures and MD trajectories
(channel radius profiles, cavities, pockets, tunnels and cast volumes) and
relates it to the residues and waters around it (lining residues, contact
networks, hydration and typed interactions). All lengths are in ångström,
areas in Å², volumes in Å³ and MD times in ps.

The names imported here are the supported entry points; each is documented
in its defining module. The main groups are:

* loading: :func:`~crevice.parser.load_structure`, :func:`~crevice.pdb.resolve_structure`,
  :func:`~crevice.trajectory.load_trajectory`;
* geometry: :func:`~crevice.channels.pore_profile`, :func:`~crevice.rolling.rolling_probe_cast`,
  :func:`~crevice.voids.void_cast`, :func:`~crevice.cavities.detect_cavities`, :func:`~crevice.tunnels.find_tunnels`;
* residues and networks: :func:`~crevice.analysis.annotate_residues`,
  :func:`~crevice.networks.build_residue_network`, :func:`~crevice.networks.build_cavity_network`;
* trajectories and statistics: :func:`~crevice.trajectory.analyze_trajectory`,
  :func:`~crevice.ensemble.profile_distribution`, :func:`~crevice.uncertainty.block_mean_confidence`;
* outputs: the ``write_*`` exporters and ``plot_*`` figure functions.

The same analyses are available from the ``crevice`` command line.

Results are geometric measurements of the supplied coordinates. They depend
on atomic radii, probe sizes and grid spacing, and are not by themselves
evidence of biological function, permeation or solvent access.

Examples
--------
>>> from crevice import Atom, StructureFrame, profile_along_path
>>> frame = StructureFrame((Atom(1, "CA", "ALA", "A", 1, 4.0, 0.0, 0.0, "C"),))
>>> round(profile_along_path(frame, [(0.0, 0.0, 0.0)]).min_radius, 2)
2.3
"""

from .benchmarks import (
    BENCHMARK_SYSTEMS,
    DEFAULT_BENCHMARK_IDS,
    BenchmarkSystem,
    DownloadRecord,
    benchmark_path,
    benchmark_systems,
    fetch_all_benchmark_structures,
    fetch_benchmark_record,
    fetch_benchmark_structure,
    get_benchmark_system,
    provenance_path,
    read_provenance,
    verify_cached_structure,
)
from .pdb import (
    DownloadRecord,
    StructureSource,
    cached_structure_path,
    fetch_structure,
    is_accession,
    normalize_accession,
    resolve_structure,
    structure_url,
)
from .cavities import cavity_surface_points, cavity_volume, classify_voids, detect_cavities
from .channels import channel_volume, detect_bottlenecks, pore_profile, profile_along_path, trace_centerline
from .figures import (
    plot_cavity_summary,
    plot_network_chord,
    plot_network_summary,
    plot_profile_radius,
    plot_profile_radius_with_residues,
    plot_residue_contacts,
    plot_trajectory_profiles,
    write_publication_pymol_script,
    write_static_publication_bundle,
)
from .io import (
    export_pymol_script,
    write_cavities_csv,
    write_connectivity_csv,
    write_features_csv,
    write_network_csv,
    write_trajectory_csv,
    write_tunnel_points_csv,
    write_tunnels_csv,
    write_void_cast_csv,
    write_cavities_json,
    write_cavities_pdb,
    write_profile_csv,
    write_network_json,
    write_profile_json,
    write_profile_pdb,
    write_residue_contacts_csv,
    write_residue_contacts_json,
    write_tunnels_json,
    write_trajectory_json,
    write_tunnels_pdb,
    write_void_cast_json,
    write_void_cast_pdb,
)
from .ml import (
    cluster_conformations,
    cluster_profiles,
    extract_features,
    profile_features,
    reduce_dimensions,
    representative_frames,
    representative_profile,
    residue_importance,
    transition_summary,
)
from .models import (
    Atom,
    Cavity,
    ChannelPoint,
    FrameAnalysis,
    NetworkEdge,
    NetworkNode,
    PoreProfile,
    ResidueContact,
    ResidueInteractionNetwork,
    StructureFrame,
    TrajectoryAnalysis,
    Tunnel,
    VoidCast,
    VoidComponent,
)
from .networks import build_cavity_network, build_residue_network, compare_networks, network_metrics
from .parser import load_structure, load_structure_report
from .volume_export import write_void_cast_dx, write_volume_viewer_bundle
from .ensemble import profile_distribution, write_profile_distribution_csv
from .connectivity import persistent_network, lining_connectivity, plot_lining_connectivity
from .radii import RadiusSet, custom_radii_set, load_radius_file, radii_preset, use_radii
from .residues import annotate_residues, bottleneck_residues, lining_residues, nearest_residues, residue_contact_table
from .schemas import SCHEMA_VERSION, result_envelope
from .static_suite import run_static_benchmark_suite
from .trajectory import (
    Trajectory,
    analyze_trajectory,
    cavity_timeseries,
    inspect_trajectory,
    load_static_frames,
    load_trajectory,
    load_trajectory_report,
    network_timeseries,
    profile_timeseries,
    trajectory_from_frames,
)
from .tunnels import cluster_tunnels, find_tunnels, rank_tunnels, tunnel_profile
from .voids import component_profile, validate_void_cast, void_cast

__all__ = [
    "write_void_cast_dx", "write_volume_viewer_bundle",
    "profile_distribution", "write_profile_distribution_csv", "persistent_network",
    "lining_connectivity", "plot_lining_connectivity",
    "Atom",
    "BENCHMARK_SYSTEMS",
    "BenchmarkSystem",
    "Cavity",
    "ChannelPoint",
    "DEFAULT_BENCHMARK_IDS",
    "FrameAnalysis",
    "NetworkEdge",
    "NetworkNode",
    "PoreProfile",
    "RadiusSet",
    "ResidueContact",
    "ResidueInteractionNetwork",
    "SCHEMA_VERSION",
    "StructureFrame",
    "Trajectory",
    "TrajectoryAnalysis",
    "Tunnel",
    "VoidComponent",
    "VoidCast",
    "analyze_trajectory",
    "annotate_residues",
    "benchmark_path",
    "benchmark_systems",
    "bottleneck_residues",
    "build_cavity_network",
    "build_residue_network",
    "cavity_surface_points",
    "cavity_timeseries",
    "cavity_volume",
    "channel_volume",
    "classify_voids",
    "cluster_conformations",
    "cluster_profiles",
    "cluster_tunnels",
    "void_cast",
    "validate_void_cast",
    "component_profile",
    "compare_networks",
    "custom_radii_set",
    "detect_bottlenecks",
    "detect_cavities",
    "export_pymol_script",
    "extract_features",
    "fetch_all_benchmark_structures",
    "fetch_benchmark_record",
    "fetch_benchmark_structure",
    "provenance_path",
    "read_provenance",
    "verify_cached_structure",
    "StructureSource",
    "cached_structure_path",
    "fetch_structure",
    "is_accession",
    "normalize_accession",
    "resolve_structure",
    "structure_url",
    "DownloadRecord",
    "find_tunnels",
    "get_benchmark_system",
    "lining_residues",
    "load_static_frames",
    "load_structure_report",
    "load_structure",
    "load_trajectory_report",
    "inspect_trajectory",
    "load_trajectory",
    "nearest_residues",
    "network_metrics",
    "network_timeseries",
    "plot_cavity_summary",
    "plot_network_chord",
    "plot_profile_radius_with_residues",
    "plot_network_summary",
    "plot_profile_radius",
    "plot_residue_contacts",
    "plot_trajectory_profiles",
    "write_publication_pymol_script",
    "write_static_publication_bundle",
    "pore_profile",
    "profile_along_path",
    "profile_features",
    "profile_timeseries",
    "radii_preset",
    "load_radius_file",
    "use_radii",
    "rank_tunnels",
    "reduce_dimensions",
    "representative_frames",
    "representative_profile",
    "residue_contact_table",
    "residue_importance",
    "result_envelope",
    "run_static_benchmark_suite",
    "trace_centerline",
    "trajectory_from_frames",
    "transition_summary",
    "tunnel_profile",
    "write_cavities_json",
    "write_cavities_pdb",
    "write_network_json",
    "write_trajectory_json",
    "write_profile_csv",
    "write_profile_json",
    "write_profile_pdb",
    "write_residue_contacts_csv",
    "write_residue_contacts_json",
    "write_tunnels_json",
    "write_tunnels_pdb",
    "write_void_cast_pdb",
    "write_void_cast_json",
    "write_cavities_csv",
    "write_connectivity_csv",
    "write_features_csv",
    "write_network_csv",
    "write_trajectory_csv",
    "write_tunnel_points_csv",
    "write_tunnels_csv",
    "write_void_cast_csv",
]

from .rolling import rolling_probe_cast
__all__.append("rolling_probe_cast")

__version__ = "0.1.0"

from .residue_evidence import boundary_residue_evidence,read_binary_dx,write_residue_evidence_bundle
from .residue_dynamics import ResidueDynamics,write_dynamics_summary
from .uncertainty import block_mean_confidence
from .trajectory import align_frame_report
__all__.extend(["boundary_residue_evidence","read_binary_dx","write_residue_evidence_bundle",
                "ResidueDynamics","write_dynamics_summary","block_mean_confidence","align_frame_report"])

from .residue_evidence import profile_constriction_evidence
__all__.append("profile_constriction_evidence")

from .residue_evidence import restore_viewer_identities
__all__.append("restore_viewer_identities")

from .residue_evidence import select_stick_residues
__all__.append("select_stick_residues")

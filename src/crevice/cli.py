"""Command-line interface for CREVICE (the ``crevice`` command).

``crevice <command> --help`` lists every option with its unit and default;
the same reference is generated in the documentation (CLI reference page).
Commands:

* static structures: ``profile``, ``residues``, ``cavities``, ``tunnels``,
  ``network``, ``features``, ``cast``, ``residue-evidence``, ``hydration``
  (without ``--trajectory``), ``analyze`` (several analyses at once) and
  ``publish`` (the complete bundle of tables, figures and viewer scenes);
* trajectories: ``trajectory``, ``cavity-trajectory``, ``hydration
  --trajectory`` and the named-region workflow ``region-init``,
  ``region-prepare``, ``region-trajectory``, ``region-compare``;
* utilities: ``fetch`` (download from RCSB), ``benchmark`` and
  ``static-suite``.

A structure argument is a local PDB/mmCIF/GRO file, an RCSB accession or an
AlphaFold DB ID such as ``AF-P02920-F1`` (downloaded and cached; an existing
path wins). Options shared by every
command that writes figures or scenes: ``--annotate`` draws titles, captions,
labels and viewer legends (clean output by default; the text is always
recorded in PNG metadata and scene JSON). Options shared by every command
that measures geometry: ``--radii`` selects the atomic radius set for the
whole command, including worker processes.

Each command function converts its arguments and calls the same public
Python functions documented in the API reference; a single-analysis command
writes its JSON to ``-o`` and its tables as CSV beside it (``<stem>.csv``).
"""

from __future__ import annotations

import argparse
import math
import json
import sys
from pathlib import Path

from .analysis import annotate_residues, detect_cavities, pore_profile
from .tunnels import find_tunnels
from .benchmarks import (
    DEFAULT_BENCHMARK_IDS,
    benchmark_path,
    benchmark_systems,
    fetch_benchmark_structure,
    get_benchmark_system,
)
from .figures import (
    plot_cavity_summary,
    plot_network_chord,
    plot_network_summary,
    plot_profile_radius,
    plot_profile_radius_with_residues,
    plot_residue_contacts,
    plot_trajectory_profiles,
    write_static_publication_bundle,
)
from .geometry import parse_coord
from .io import (
    export_pymol_script,
    write_cavities_csv,
    write_connectivity_csv,
    write_features_csv,
    write_network_csv,
    write_tunnel_points_csv,
    write_tunnels_csv,
    write_cavities_json,
    write_cavities_pdb,
    write_network_json,
    write_profile_csv,
    write_profile_json,
    write_profile_pdb,
    write_trajectory_json,
    write_residue_contacts_csv,
    write_residue_contacts_json,
    write_tunnels_json,
    write_tunnels_pdb,
)
from .ml import profile_features
from .networks import build_cavity_network, build_residue_network, network_metrics
from .parser import PARSERS, load_structure_report
from .pdb import fetch_structure, resolve_structure
from .static_suite import run_static_benchmark_suite
from .trajectory import analyze_trajectory, load_static_frames
from .trajectory import DEFAULT_MAX_FRAMES, inspect_trajectory, load_trajectory_report
from .ensemble import profile_distribution, write_profile_distribution_csv
from .connectivity import persistent_network, lining_connectivity, plot_lining_connectivity


def build_parser() -> argparse.ArgumentParser:
    """Build the ``crevice`` argument parser with every subcommand.

    Each subparser stores its handler as ``func``; :func:`main` calls it. The
    parser is also read by the documentation build to generate the CLI
    reference.

    Returns
    -------
    argparse.ArgumentParser
        The top-level parser (program name ``crevice``).

    Examples
    --------
    >>> from crevice.cli import build_parser
    >>> args = build_parser().parse_args(["profile", "1GRM", "-o", "out.json"])
    >>> args.output, args.enclosure_radius, args.csv   # None = automatic probe
    ('out.json', None, None)
    """
    parser = argparse.ArgumentParser(
        prog="crevice",
        description=(
            "CREVICE: Cavity and Residue Environment Visualisation, Interaction and "
            "Connectivity Evaluation, for static structures and trajectories."
        ),
        epilog=(
            "Figures and molecular viewer scenes are clean by default: no titles, captions, "
            "value call-outs, residue labels or viewer legends. Add --annotate to any command "
            "that writes figures or scenes to draw them; the text is always recorded in each "
            "PNG's metadata and in the scene JSON files."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    region_init = subparsers.add_parser("region-init", help="Create a candidate region definition with hashed inputs and explicit landmarks")
    region_init.add_argument("topology", help="MD topology (or a static structure) of the system, hashed into the definition")
    region_init.add_argument("--trajectory", help="MD trajectory of the system, hashed into the definition (required later by region-trajectory)")
    region_init.add_argument("--region-id", required=True, help="Short region identifier (a simple file stem); used as the output prefix of region-trajectory")
    region_init.add_argument("--landmarks-json", help="List of full residue IDs, anchor atom names and evidence")
    region_init.add_argument("--reference-volume-dx", help="Existing reference in source frame 0; otherwise propose around local landmarks")
    region_init.add_argument("--output", required=True, help="New definition JSON path; an existing file is refused (use a new versioned name)")
    region_init.add_argument("--selection", default="protein and not name H*", help="MDAnalysis selection of the analysed protein atoms")
    region_init.add_argument("--obstacle-selection", default="protein and not name H*", help="MDAnalysis selection of every atom that bounds the region (must include the analysed protein heavy atoms)")
    region_init.add_argument("--spacing", type=float, default=.5, help="Region grid spacing in Å")
    region_init.add_argument("--probe-radius", type=float, default=1.4, help="Probe radius in Å that defines the probe-swept region volume")
    region_init.add_argument("--radius", type=float, default=6., help="Without --reference-volume-dx: radius in Å of the landmark sphere in which the reference region is proposed")
    region_init.set_defaults(func=_cmd_region_init)
    region_prepare = subparsers.add_parser("region-prepare", help="Validate a named region, bind its reference frame and write native review overlays")
    region_prepare.add_argument("definition", help="Region definition JSON from region-init (edited and reviewed)")
    region_prepare.add_argument("--out-dir", required=True, help="Output directory for the prepared definition, reference map and review overlays")
    region_prepare.set_defaults(func=_cmd_region_prepare)
    region_compare = subparsers.add_parser("region-compare", help="Compare separate named regions on identical trajectory frames")
    region_compare.add_argument("left_statistics", help="The <region-id>_cavity_statistics.json of the first region-trajectory run")
    region_compare.add_argument("right_statistics", help="The <region-id>_cavity_statistics.json of the second region-trajectory run (same trajectory and frames)")
    region_compare.add_argument("--left-definition", help="Prepared definition required for legacy observations lacking region provenance")
    region_compare.add_argument("--right-definition", help="As --left-definition, for the right run")
    region_compare.add_argument("--left-input-provenance", help="Historical input SHA256 record required when binding legacy observations")
    region_compare.add_argument("--right-input-provenance", help="As --left-input-provenance, for the right run")
    region_compare.add_argument("--out-dir", required=True, help="New or empty output directory")
    region_compare.set_defaults(func=_cmd_region_compare)
    region_md = subparsers.add_parser("region-trajectory", help="Run geometry and hydration using a prepared, versioned region definition")
    region_md.add_argument("definition", help="Prepared region definition JSON from region-prepare")
    region_md.add_argument("--out-dir", required=True, help="Output directory")
    region_md.add_argument("--workers", type=int, default=1, help="Worker processes for the per-frame geometry; results do not depend on it")
    region_md.add_argument("--max-frames", type=int, default=DEFAULT_MAX_FRAMES, help="Refuse to load more frames than this")
    region_md.add_argument("--stop", type=int, help="Stop before this trajectory frame (default: all frames)")
    region_md.add_argument("--confidence", type=float, default=.95, help="Nominal level of the approximate mean intervals; 0 gives descriptive statistics only")
    region_md.add_argument("--block-length", type=int, help="Bootstrap block length in frames; default chosen from the autocorrelation")
    region_md.add_argument("--bootstrap-replicates", type=int, default=2000, help="Block-bootstrap replicates of the mean intervals")
    region_md.add_argument("--bootstrap-seed", type=int, default=20260913, help="Random seed of the block bootstrap")
    region_md.add_argument("--dpi", type=int, default=240, help="Figure resolution in dots per inch")
    region_md.set_defaults(func=_cmd_region_trajectory)

    hydration = subparsers.add_parser("hydration", help="Observed residue water contacts, geometric SASA and sampled trajectory hydration")
    _add_structure_arg(hydration)
    hydration.add_argument("--trajectory", help="Optional MD trajectory; structure becomes its topology")
    hydration.add_argument("--cavity-results", help="Existing cavity trajectory statistics JSON, to reuse its residue focus and aligned geometry")
    hydration.add_argument("--volume-dx", help="Static measured cavity DX to select boundary residues and nonlocal partners")
    hydration.add_argument("--residues-json", help="Optional list of full residue IDs")
    hydration.add_argument("--out-dir", required=True, help="Output directory")
    hydration.add_argument("--prefix", default="crevice", help="Output file-name prefix (a simple file stem)")
    hydration.add_argument("--selection", default="protein and not name H*", help="MDAnalysis selection of the analysed protein atoms (trajectory input)")
    hydration.add_argument("--water-selection", help="MDAnalysis selection of water residues (one oxygen per residue)")
    hydration.add_argument("--start", type=int, default=0, help="First trajectory frame to read")
    hydration.add_argument("--stop", type=int, help="Stop before this trajectory frame (default: all frames)")
    hydration.add_argument("--stride", type=int, default=1, help="Read every Nth trajectory frame")
    hydration.add_argument("--max-frames", type=int, default=DEFAULT_MAX_FRAMES, help="Refuse to load more frames than this")
    hydration.add_argument("--pbc", choices=["check","none"], default="check", help="check: minimum-image water distances in the periodic cell (cell dimensions required); none: no periodic imaging")
    hydration.add_argument("--confidence", type=float, default=.95, help="Mean confidence level; zero requests descriptive statistics only")
    hydration.add_argument("--block-length", type=int, help="Bootstrap block length in frames; default chosen from the autocorrelation")
    hydration.add_argument("--bootstrap-replicates", type=int, default=2000, help="Block-bootstrap replicates of the mean intervals (at least 200)")
    hydration.add_argument("--bootstrap-seed", type=int, default=20260914, help="Random seed of the block bootstrap")
    hydration.add_argument("--dpi", type=int, default=240, help="Figure resolution in dots per inch")
    hydration.set_defaults(func=_cmd_hydration)

    cavity_md = subparsers.add_parser("cavity-trajectory", help="Track measured cavities in every aligned trajectory frame; volumes, changing residue evidence and sectional width statistics")
    cavity_md.add_argument("topology", help="MD topology (.gro, .pdb, .psf, .tpr)")
    cavity_md.add_argument("trajectory", help="MD trajectory (.xtc, .dcd, ...)")
    cavity_md.add_argument("--reference-volume-dx", required=True, help="Measured binary cavity map in the first trajectory frame's coordinates")
    cavity_md.add_argument("--out-dir", required=True, help="Output directory")
    cavity_md.add_argument("--prefix", default="crevice", help="Output file-name prefix (a simple file stem)")
    cavity_md.add_argument("--selection", default="protein and not name H*", help="MDAnalysis selection of the analysed protein atoms")
    cavity_md.add_argument("--obstacle-selection", help="Independent MDAnalysis physical-obstacle selection; must retain analysed protein heavy atoms. Reference-region mode only; default protein-only geometry")
    cavity_md.add_argument("--grid-phase", type=parse_coord, default=(0.,0.,0.), help="Lattice phase fractions x,y,z in [0,1) for numerical sensitivity")
    cavity_md.add_argument("--pbc", choices=['check','none','unwrap'], default='check', help="check: screen for split molecules and use cell dimensions; none: no periodic handling; unwrap: make the selection whole (refused with --water-membership)")
    cavity_md.add_argument("--max-frames", type=int, default=DEFAULT_MAX_FRAMES, help="Refuse to load more frames than this")
    cavity_md.add_argument("--stop", type=int, help="Stop before this trajectory frame (default: all frames)")
    cavity_md.add_argument("--workers", type=int, default=1, help="Worker processes for the per-frame geometry; results do not depend on it")
    cavity_md.add_argument("--max-grid-points", type=int, default=2_000_000, help="Explicit per-frame grid allocation guard; no automatic coarsening")
    cavity_md.add_argument("--spacing", type=float, help="Fixed lattice spacing in Å; default the reference map's own spacing")
    cavity_md.add_argument("--axis", type=parse_coord, help="Fixed analysis direction x,y,z in the reference frame; default measured-map third direction")
    cavity_md.add_argument("--profile-padding", type=float, default=4., help="Extra axial range in Å beyond the reference cavity covered by the section profile")
    cavity_md.add_argument("--geometry-mode",choices=["rolling","reference-region"],default="rolling",help="Global component matching or explicitly bounded reference-neighborhood volume")
    cavity_md.add_argument("--region-margin",type=float,default=2.,help="Reference-neighborhood growth distance in Å; reference-region mode only")
    cavity_md.add_argument("--probe-radius", type=float, default=.8, help="Probe radius in Å of the per-frame cavity measurement")
    cavity_md.add_argument("--outer-radius", type=float, default=6., help="Rolling mode: outer envelope probe radius in Å")
    cavity_md.add_argument("--cavity-selection", choices=['dominant','all'], default='dominant', help="Rolling mode: dominant interior cores, or all fixed-probe voids")
    cavity_md.add_argument("--enclosure-fraction", type=float, default=.9, help="Rolling mode: burial threshold (fraction of 26 lattice rays)")
    cavity_md.add_argument("--min-component-volume", type=float, default=1., help="Rolling mode: smallest component kept, in Å³")
    cavity_md.add_argument("--minimum-overlap", type=float, default=.1, help="Rolling mode: smallest accepted overlap score with the reference, in (0, 1]")
    cavity_md.add_argument("--ambiguity-ratio", type=float, default=.8, help="Rolling mode: a match is ambiguous (unresolved) when the second-best score is at least this fraction of the best")
    cavity_md.add_argument("--lining-distance", type=float, default=1.5, help="Largest atom-surface gap in Å at which a boundary face is assigned to a residue")
    cavity_md.add_argument("--contact-cutoff", type=float, default=4.5, help="Residue contact cutoff in Å (heavy-atom centres)")
    cavity_md.add_argument("--alignment-residues-json", help="JSON list of full residue IDs defining the rigid fit; default all matching protein atoms")
    cavity_md.add_argument("--confidence", type=float, default=.95, help="Nominal level of the approximate mean intervals; 0 gives descriptive statistics only")
    cavity_md.add_argument("--block-length", type=int, help="Bootstrap block length in frames; default chosen from the autocorrelation")
    cavity_md.add_argument("--bootstrap-replicates", type=int, default=2000, help="Block-bootstrap replicates of the mean intervals")
    cavity_md.add_argument("--bootstrap-seed", type=int, default=20260913, help="Random seed of the block bootstrap")
    cavity_md.add_argument("--dpi", type=int, default=240, help="Figure resolution in dots per inch")
    cavity_md.set_defaults(func=_cmd_cavity_trajectory)

    cast = subparsers.add_parser("cast", help="Fill protein-interior voids with a 3D rolling probe, including open pockets and multiple regions")
    _add_structure_arg(cast)
    cast.add_argument("--out-dir", required=True, help="Output directory for the cast bundle")
    cast.add_argument("--prefix", help="Output file-name prefix; defaults to the structure file stem")
    cast.add_argument("--spacing", type=float, default=.5, help="Requested grid spacing in Å; never silently coarsened")
    cast.add_argument("--probe-radius", type=float, default=.8,
                      help="Inner rolling probe radius in Å (fixed default 0.8, no 'auto'). Unlike the profile's "
                           "--enclosure-radius, which only decides which wall gaps count as closed while the radii "
                           "come from atom-surface clearance, this probe defines the cast itself: its centres and "
                           "swept spheres set the reported volume and connectivity, so it is a measurement setting, "
                           "not chosen automatically")
    cast.add_argument("--seed", type=parse_coord, help="Interior cavity seed x,y,z in Å; retain only its probe-connected interior")
    cast.add_argument("--seed-tolerance", type=float, default=.75, help="Maximum atom-clear seed-to-grid snap distance in Å")
    cast.add_argument("--outer-radius", type=float, default=6., help="Outer envelope probe radius in Å")
    cast.add_argument("--enclosure-fraction", type=float, default=.9, help="Burial threshold (fraction of 26 lattice rays). Dominant selection: qualifies cavity cores and the sphere fill can extend beyond them; --cavity-selection all: hard crop of probe centres")
    cast.add_argument("--min-depth", type=float, default=0., help="Additional depth inside the outer envelope in Å")
    cast.add_argument("--min-component-volume", type=float, default=50., help="Minimum sampled region volume in Å³")
    cast.add_argument("--max-cavities", "--max-components", dest="max_components", type=int, default=1,
                      help="Maximum qualifying cavity regions (default 1); increasing it never forces extra results")
    cast.add_argument("--cavity-selection", choices=['dominant','all'], default='dominant',
                      help="Dominant: substantial interior cores. All: legacy fixed-probe void enumeration")
    cast.add_argument("--core-fraction", type=float, default=.08, help="Minimum core volume relative to the strongest core (dominant mode)")
    _add_smoothing_args(cast)
    cast.add_argument("--max-grid-points", type=int, default=16_000_000, help="Largest allowed grid; a larger request is an error (no automatic coarsening)")
    cast.add_argument("--selection", choices=['protein','all'], default='protein', help="Protein: non-HETATM heavy atoms. All: retain bound species/solvent as obstacles")
    cast.add_argument("--dpi", type=int, default=250, help="Resolution of the hydration/analysis figures written with the cast (the cast itself has no Matplotlib figure; use the viewer scenes)")
    cast.set_defaults(func=_cmd_cast)

    profile = subparsers.add_parser("profile", help="Measure the radius profile of a channel that runs through the protein")
    _add_structure_arg(profile)
    _add_profile_args(profile)
    profile.add_argument("-o", "--output", required=True,
                         help="Output JSON path; the per-sample table is also written as <stem>.csv beside it")
    profile.add_argument("--csv", help="Per-sample CSV path instead of the default <stem>.csv beside the JSON")
    profile.add_argument("--pdb", help="Optional profile-trace dummy-atom PDB output")
    profile.add_argument("--png", help="Optional 400 dpi pore-radius profile PNG output")
    profile.add_argument("--annotated-png", help="Optional 400 dpi pore-radius PNG with a lane of nearest-residue landmarks (residue names are drawn with --annotate)")
    profile.add_argument("--residue-label-every", type=int, default=5, help="Landmark spacing in profile points for the residue-landmark radius figure")
    profile.add_argument("--dpi", type=int, default=400, help="Figure DPI for PNG outputs")
    profile.add_argument("--pymol", help="Optional publication-ready PyMOL script path; requires --pdb")
    profile.add_argument("--pymol-png", help="Optional PyMOL-rendered PNG path written into the .pml script")
    profile.set_defaults(func=_cmd_profile)

    residues = subparsers.add_parser("residues", help="Annotate channel-lining residues")
    _add_structure_arg(residues)
    _add_profile_args(residues)
    residues.add_argument("-o", "--output", required=True, help="Output CSV path")
    residues.add_argument("--json", help="Optional residue JSON output")
    residues.add_argument("--png", help="Optional 400 dpi residue influence PNG output")
    residues.add_argument("--dpi", type=int, default=400, help="Figure DPI for PNG outputs")
    residues.add_argument("--cutoff", type=float, default=4.5, help="Residue surface-distance cutoff")
    residues.set_defaults(func=_cmd_residues)

    cavities = subparsers.add_parser("cavities", help="Find buried cavities on a clearance grid")
    _add_structure_arg(cavities)
    cavities.add_argument("-o", "--output", required=True,
                          help="Output JSON path; the cavity table is also written as <stem>.csv beside it")
    cavities.add_argument("--pdb", help="Optional dummy-atom PDB output")
    cavities.add_argument("--png", help="Optional 400 dpi cavity summary PNG output")
    cavities.add_argument("--dpi", type=int, default=400, help="Figure DPI for PNG outputs")
    cavities.add_argument("--spacing", type=float, default=2.0, help="Grid spacing in Å")
    cavities.add_argument("--min-radius", type=float, default=1.4, help="Minimum void clearance radius in Å")
    cavities.add_argument("--max-cavities", type=int, default=10, help="Maximum cavities to report")
    cavities.add_argument("--include-hydrogen", action="store_true", help="Include hydrogen atoms")
    cavities.add_argument("--include-hetero", action="store_true", help="Include HETATM records")
    cavities.set_defaults(func=_cmd_cavities)

    tunnels = subparsers.add_parser("tunnels", help="Find the widest access tunnels from a start point to the outside")
    _add_structure_arg(tunnels)
    tunnels.add_argument("-o", "--output", required=True,
                         help="Output JSON path; tunnel and path-node tables are also written as <stem>.csv and <stem>_points.csv")
    tunnels.add_argument("--pdb", help="Optional dummy-atom PDB output")
    tunnels.add_argument("--start", help="Start coordinate formatted as x,y,z in Å; defaults to structure centroid")
    tunnels.add_argument("--spacing", type=float, default=2.0, help="Grid spacing in Å")
    tunnels.add_argument("--min-radius", type=float, default=0.8, help="Minimum traversable clearance radius in Å")
    tunnels.add_argument("--max-tunnels", type=int, default=3, help="Maximum tunnels to report")
    tunnels.add_argument("--include-hydrogen", action="store_true", help="Include hydrogen atoms")
    tunnels.add_argument("--include-hetero", action="store_true", help="Include HETATM records")
    tunnels.set_defaults(func=_cmd_tunnels)

    network = subparsers.add_parser("network", help="Build residue and channel/cavity interaction networks")
    _add_structure_arg(network)
    _add_profile_args(network)
    network.add_argument("-o", "--output", required=True,
                         help="Output network JSON path; node and edge tables are also written as <stem>_nodes.csv and <stem>_edges.csv")
    network.add_argument("--cutoff", type=float, default=4.5, help="Residue/residue and region/residue cutoff in Å")
    network.add_argument("--contact-metric",choices=["center","surface"],default="center",help="Atom-center distance by default; surface retains the legacy van der Waals gap. The networks inside publish/analyze use the surface gap with heteroatoms excluded")
    network.add_argument("--region", choices=["none", "profile"], default="profile", help="Add a region pseudo-node")
    network.add_argument("--png", help="Optional 400 dpi network summary PNG output")
    network.add_argument("--chord-png", help="Optional 400 dpi chord-style network PNG output")
    network.add_argument("--dpi", type=int, default=400, help="Figure DPI for PNG outputs")
    network.set_defaults(func=_cmd_network)
    network.add_argument("--residue-groups", help="JSON mapping full residue IDs to helix/strand names")
    network.add_argument("--connectivity-json", help="Paths and group contact record; the tables are also written as <stem>.csv and <stem>_groups.csv")
    network.add_argument("--connectivity-png", help="Group-colored graph-distance plot")
    network.add_argument("--remove-residue", action="append", default=[], help="Residue ID for graph removal sensitivity")

    evidence=subparsers.add_parser("residue-evidence",help="Explain a measured cast boundary and nonlocal residue partners")
    _add_structure_arg(evidence)
    evidence.add_argument("--volume-dx",required=True,help="Measured binary volume DX in the protein coordinate frame")
    evidence.add_argument("--out-dir",required=True, help="Output directory")
    evidence.add_argument("--prefix",default="crevice", help="Output file-name prefix (a simple file stem)")
    evidence.add_argument("--lining-distance",type=float,default=1.5, help="Largest atom-surface gap in Å at which a boundary face is assigned to a residue")
    evidence.add_argument("--contact-cutoff",type=float,default=4.5, help="Residue contact cutoff in Å (heavy-atom centres) for nonlocal partners")
    evidence.add_argument("--residue-groups",help="JSON mapping full residue IDs to helix/strand groups")
    evidence.add_argument("--scene",help="Optional matching volume PML for a cartoon-color evidence overlay")
    evidence.add_argument("--stick-residues",default="top:12",metavar="all|top:N|IDS|@FILE",
                          help="Residues shown as opaque sticks in the PyMOL/VMD/ChimeraX overlays: 'top:N' per role "
                               "(default top:12), 'all' boundary residues and nonlocal partners, or comma-separated "
                               "residue IDs / @file of IDs (each must be a boundary residue or nonlocal partner). "
                               "The ranked bar chart and tables are unchanged")
    evidence.add_argument("--dpi",type=int,default=300, help="Figure resolution in dots per inch")
    evidence.set_defaults(func=_cmd_residue_evidence)

    features = subparsers.add_parser("features", help="Export profile and residue feature JSON")
    _add_structure_arg(features)
    _add_profile_args(features)
    features.add_argument("-o", "--output", required=True,
                          help="Output JSON path; the features are also written as <stem>.csv (feature, value, unit)")
    features.add_argument("--cutoff", type=float, default=4.5, help="Residue surface-distance cutoff in Å")
    features.set_defaults(func=_cmd_features)

    trajectory = subparsers.add_parser("trajectory", help="Analyze a sequence of static structure frames")
    trajectory.add_argument("structures", nargs="+", help="Input PDB/mmCIF files, one per frame")
    trajectory.add_argument("-o", "--output", help="Output trajectory-analysis JSON path (required unless --inspect); "
                            "per-frame summaries and profiles are also written as <stem>_frames.csv and <stem>_profiles.csv")
    trajectory.add_argument(
        "--analysis",
        action="append",
        choices=["profile", "residues", "cavities", "tunnels", "network", "features"],
        help="Analysis to run; may be repeated. Defaults to profile and features.",
    )
    trajectory.add_argument("--stride", type=int, default=1, help="Analyze every Nth input frame")
    trajectory.add_argument("--axis", default="auto", help="Channel axis: auto, x, y, z, or vector")
    trajectory.add_argument("--origin", type=parse_coord, help="Reference channel seed x,y,z")
    trajectory.add_argument("--section-spacing", type=float, default=0.5, help="Transverse channel grid spacing in Å")
    trajectory.add_argument("--enclosure-radius", type=_enclosure_radius_arg, default=None, metavar="{auto,Å}",
                            help="Probe radius (Å) for channel enclosure, fixed for every frame. 'auto' (default) "
                                 "chooses it once on the first (reference) frame with the profile rule (smallest "
                                 "probe of the most persistent channel, 0.5-3.0 Å) and then uses that value for every "
                                 "frame; the choice and reason are recorded in the output JSON "
                                 "(metadata.enclosure_probe) and printed. A number runs that probe for every frame")
    trajectory.add_argument("--probe-fallback", action=argparse.BooleanOptionalAction, default=None,
                            help="Per-frame fallback: if the pinned enclosure probe does not resolve a frame, retry "
                                 "that frame with the automatic choice made on it (same axis and origin). Default: "
                                 "on with --enclosure-radius auto, off with an explicit probe (kept strict); "
                                 "--probe-fallback enables it for an explicit probe, --no-probe-fallback disables "
                                 "it. Each frame's source (pinned/fallback/none) is in <stem>_frames.csv "
                                 "(enclosure_probe_source) and the JSON; fallbacks are counted and printed")
    trajectory.add_argument("--samples", type=int, default=81, help="Minimum number of profile samples; connected profiles add points to keep axial steps ≤ 0.75 Å")
    trajectory.add_argument("--search-radius", type=float, default=6.0, help="Transverse search half-width (Å); 0 requests an unvalidated fixed-axis scan")
    trajectory.add_argument("--png", help="Axial mean, frame fluctuations and approximate confidence band")
    trajectory.add_argument("--topology", help="MD topology (.gro, .pdb, .psf, .tpr); pass one trajectory file as structures")
    trajectory.add_argument("--selection", default="protein",
                            help="MDAnalysis atom selection applied when reading a trajectory")
    trajectory.add_argument("--start", type=int, default=0, help="First trajectory frame to read")
    trajectory.add_argument("--stop", type=int, help="Stop before this trajectory frame")
    trajectory.add_argument("--max-frames", type=int, default=DEFAULT_MAX_FRAMES,
                            help="Refuse to load more frames than this; 0 disables the guard")
    trajectory.add_argument("--inspect", action="store_true",
                            help="Report frame and selection counts without loading or analyzing")
    trajectory.add_argument("--report-json", help="Optional path for the reader provenance report")
    trajectory.add_argument("--align",action=argparse.BooleanOptionalAction,default=True,help="Fit matching CA atoms to the first frame (default); --no-align disables fitting")
    trajectory.add_argument("--alignment-residue",action="append",help="Full residue ID in a stable alignment scaffold; repeat to select multiple residues")
    trajectory.add_argument("--alignment-residues-json",help="JSON list of full residue IDs for a preselected alignment scaffold")
    trajectory.add_argument("--pbc",choices=["check","none","unwrap"],default="check",help="Screen split bonds before fitting; unwrap requires one reliably bonded fragment")
    trajectory.add_argument("--profile-errors",choices=["record","raise"],default="record",help="Retain unresolved profile status by default")
    trajectory.add_argument("--confidence",type=float,default=.95,help="Approximate mean confidence level; 0 disables confidence estimation")
    trajectory.add_argument("--block-length",type=int,help="Target contiguous batch length in analyzed frames; default uses autocorrelation")
    trajectory.add_argument("--bootstrap-replicates",type=int,default=2000, help="Block-bootstrap replicates of the mean intervals")
    trajectory.add_argument("--bootstrap-seed",type=int,default=20260911, help="Random seed of the block bootstrap")
    trajectory.add_argument("--profile-quantity",choices=["radius","diameter"],default="radius", help="Plot radius or diameter (Å) in the --png figure")
    trajectory.add_argument("--contact-metric",choices=["center","surface"],default="center",help="Residue-contact distance: heavy-atom centres (default) or van der Waals surface gap; the Python analyze_trajectory default is the surface gap")
    trajectory.add_argument("--contact-cutoff",type=float,default=4.5, help="Residue contact cutoff in Å for --network-json")
    trajectory.add_argument("--lining-evidence",help="Reference residue-evidence JSON; residue identities must match the trajectory")
    trajectory.add_argument("--view", choices=["distribution", "timeseries"], default="distribution", help="--png figure: distribution (mean, frame quantile band and mean interval against position) or timeseries (each frame's minimum, mean and maximum against time)")
    trajectory.add_argument("--quantiles", type=float, nargs=2, default=(0.1, 0.9), help="Lower and upper frame quantiles of the fluctuation band")
    trajectory.add_argument("--distribution-csv", help="Coordinate-wise radius statistics (Å) and frame counts; the interval record is also written as <stem>.json")
    trajectory.add_argument("--network-json", help="Persistent physical contacts and descriptive aligned residue motion (also writes CSV/PNG)")
    trajectory.add_argument("--occupancy", type=float, default=0.75, help="Minimum contact occupancy")
    trajectory.add_argument("--residue-groups", help="JSON mapping residue IDs to helix/strand names")
    trajectory.add_argument("--dpi", type=int, default=400, help="Figure DPI for PNG outputs")
    trajectory.set_defaults(func=_cmd_trajectory)

    publish = subparsers.add_parser("publish", help="Generate a publication-ready static-structure output bundle")
    _add_structure_arg(publish)
    _add_profile_args(publish)
    publish.add_argument("--out-dir", required=True, help="Directory for the publication-ready output bundle")
    publish.add_argument("--prefix", help="Output filename prefix; defaults to the structure stem")
    publish.add_argument("--dpi", type=int, default=400, help="Figure DPI for publication PNG outputs")
    publish.add_argument("--residue-label-every", type=int, default=5, help="Landmark spacing in profile points for the residue-landmark radius figure")
    publish.add_argument("--cast-mode", choices=["auto", "channel", "cavity", "all", "rolling"], default="auto", help="Void-cast component selection mode for PyMOL pore/cavity fill")
    publish.add_argument("--cast-outer-radius", type=float, default=6., help="Outer envelope probe for rolling cast mode")
    publish.add_argument("--cast-enclosure-fraction", type=float, default=.9, help="Rolling cast: burial threshold (fraction of 26 lattice rays)")
    publish.add_argument("--cast-extension", type=float, default=2.0,
                         help="Display-only cast continuation beyond each mouth (Å); excluded from measured volume")
    publish.add_argument("--cast-spacing", type=float, default=0.5, help="Void-cast grid spacing in Å")
    publish.add_argument("--cast-min-radius", type=float, default=0.0, help="Minimum atom-surface clearance in Å for void-cast points")
    publish.add_argument("--cast-max-cavities", "--cast-max-components", dest="cast_max_components", type=int, default=1,
                         help="Maximum qualifying cast regions (default 1)")
    publish.add_argument("--cast-selection", choices=['dominant','all'], default='dominant', help="Selection for rolling casts and automatic rolling fallback")
    publish.add_argument("--cast-core-fraction", type=float, default=.08, help="Rolling cast (dominant selection): minimum core volume relative to the strongest core")
    _add_smoothing_args(publish)
    publish.add_argument("--cast-max-grid-points", type=int, default=16000000, help="Maximum grid points for void-cast analysis")
    publish.add_argument("--focus", type=float, nargs=3, action="append", help="Focus coordinate x y z in Å; repeat for branches")
    publish.add_argument("--focus-radius", type=float, default=8.0, help="Radius in Å around each --focus point")
    publish.add_argument("--min-component-volume", type=float, default=0.0, help="Inclusive minimum final component volume in Å³ after cropping")
    publish.add_argument("--max-component-volume", type=float, help="Inclusive maximum final component volume in Å³ after cropping")
    publish.add_argument("--residue-groups", help="JSON mapping residue IDs to helix/strand names")
    publish.add_argument("--skip-cavities", action="store_true", help="Do not run cavity detection for the bundle")
    publish.add_argument("--skip-tunnels", action="store_true", help="Do not run tunnel detection for the bundle")
    publish.add_argument("--skip-network", action="store_true", help="Do not build the residue/channel network for the bundle")
    publish.add_argument("--entry-end", choices=["auto", "start", "end"], default="auto",
                         help="Which channel end is labelled ENTRY for the cast segments (ENTRY, LUMEN, EXIT). auto: "
                              "if exactly one end is capped with lateral exits, those legs are EXITs and the open end "
                              "is the ENTRY; otherwise the axis start. Geometric labels, not a transport direction")
    publish.add_argument("--lining-cutoff", type=float, default=3.3,
                         help="Lining residues: atom centre within this distance (Å) of a segment's measured grid "
                              "nodes (default 3.3, a hydrogen-bond distance); written per segment to "
                              "PREFIX_lining_residues.csv and as hidden viewer selections")
    publish.set_defaults(func=_cmd_publish)

    fetch = subparsers.add_parser("fetch", help="Fetch structures from RCSB by accession, or AlphaFold DB models by AF- identifier")
    fetch.add_argument(
        "pdb_ids",
        nargs="*",
        help=f"Accessions to fetch, e.g. 1GRM 4PYP; defaults to {', '.join(DEFAULT_BENCHMARK_IDS)}",
    )
    fetch.add_argument("--cache-dir", default=".crevice/pdb", help="Local structure cache directory")
    fetch.add_argument("--format", choices=["pdb", "cif", "mmcif"], default="cif", help="Downloaded structure format")
    fetch.add_argument("--assembly", type=int, help="Fetch this biological assembly instead of the deposited unit")
    fetch.add_argument("--overwrite", action="store_true", help="Overwrite existing cached files")
    fetch.add_argument("--timeout", type=float, default=30.0, help="Network timeout in seconds")
    fetch.add_argument("--list", action="store_true", help="List configured benchmark systems without downloading")
    fetch.add_argument("--json", help="Optional path for the download provenance records")
    fetch.set_defaults(func=_cmd_fetch)

    fetch_example = subparsers.add_parser("fetch-example", help="Download an example dataset by name into a local cache, verifying SHA-256 hashes")
    fetch_example.add_argument("name", nargs="?", help="Dataset name, e.g. glut1-excerpt (see --list)")
    fetch_example.add_argument("--cache-dir", help="Cache root; files go in CACHE_DIR/NAME (default: $CREVICE_EXAMPLE_CACHE_DIR or ~/.cache/crevice/examples)")
    fetch_example.add_argument("--base-url", help="URL prefix of the assets (default: $CREVICE_EXAMPLE_BASE_URL or the v0.1.1 GitHub Release)")
    fetch_example.add_argument("--source-dir", help="Local directory holding the asset files, for offline use (default: $CREVICE_EXAMPLE_SOURCE_DIR)")
    fetch_example.add_argument("--overwrite", action="store_true", help="Fetch every file again even if a verified copy is cached")
    fetch_example.add_argument("--timeout", type=float, default=60.0, help="Network timeout in seconds per file")
    fetch_example.add_argument("--list", action="store_true", help="List the available datasets and their files without downloading")
    fetch_example.set_defaults(func=_cmd_fetch_example)

    analyze = subparsers.add_parser("analyze", help="Run selected analyses on one structure")
    _add_structure_arg(analyze)
    _add_profile_args(analyze)
    analyze.add_argument(
        "--analysis", action="append",
        choices=["profile", "residues", "cavities", "tunnels", "network", "features"],
        help="Analysis to run; may be repeated. Defaults to profile and residues.",
    )
    analyze.add_argument("--out-dir", required=True, help="Directory for analysis outputs")
    analyze.add_argument("--prefix", help="Output filename prefix; defaults to the structure stem")
    analyze.add_argument("--cutoff", type=float, default=4.5, help="Residue surface-distance cutoff in Å")
    analyze.add_argument("--cavity-spacing", type=float, default=2.0, help="Cavity grid spacing in Å")
    analyze.add_argument("--cavity-min-radius", type=float, default=1.4, help="Minimum cavity clearance radius in Å")
    analyze.add_argument("--max-cavities", type=int, default=10, help="Maximum cavities to report")
    analyze.add_argument("--cavity-max-grid-points", type=int, default=75_000,
                         help="Grid-point cap for cavity detection; the requested spacing is "
                              "coarsened to stay under it, so raise this to honour a fine spacing")
    analyze.add_argument("--tunnel-spacing", type=float, default=2.0, help="Tunnel grid spacing in Å")
    analyze.add_argument("--tunnel-min-radius", type=float, default=0.8, help="Minimum traversable clearance in Å")
    analyze.add_argument("--max-tunnels", type=int, default=3, help="Maximum tunnels to report")
    analyze.add_argument("--tunnel-max-grid-points", type=int, default=120_000,
                         help="Grid-point cap for tunnel search; see --cavity-max-grid-points")
    analyze.add_argument("--start", help="Tunnel start coordinate as x,y,z; defaults to the centroid")
    analyze.add_argument("--png", action="store_true", help="Also write figures for the profile, residue, cavity and network analyses (tunnel data are tables only)")
    analyze.add_argument("--dpi", type=int, default=400, help="Figure DPI for PNG outputs")
    analyze.set_defaults(func=_cmd_analyze)

    benchmark = subparsers.add_parser("benchmark", help="Run static checks on configured CREVICE PDB benchmarks")
    benchmark.add_argument("pdb_ids", nargs="*", help=f"PDB IDs to check; defaults to {', '.join(DEFAULT_BENCHMARK_IDS)}")
    benchmark.add_argument("--cache-dir", default=".crevice/pdb", help="Local structure cache directory")
    benchmark.add_argument("--format", choices=["pdb", "cif", "mmcif"], default="cif", help="Structure format")
    benchmark.add_argument("--fetch", action="store_true", help="Download missing structures before checking")
    benchmark.add_argument("--output", help="Optional JSON summary path; the rows are also written as <stem>.csv")
    benchmark.add_argument("--samples", type=int, default=41, help="Profile samples for benchmark summaries")
    benchmark.add_argument("--search-radius", type=float, default=3.0, help="Transverse search half-width (Å) of each benchmark profile")
    benchmark.add_argument("--enclosure-radius", type=_enclosure_radius_arg, default=None, metavar="{auto,Å}",
                           help="Enclosure probe (Å) of each profile. 'auto' (default) chooses it per system with the "
                                "profile rule (smallest probe of the most persistent channel, 0.5-3.0 Å) and records "
                                "it (enclosure_radius_A, enclosure_probe_mode, enclosure_probe_reason); a number "
                                "runs that probe (0.8 reproduces the summaries recorded before auto became the default)")
    benchmark.set_defaults(func=_cmd_benchmark)

    static_suite = subparsers.add_parser("static-suite", help="Run publication-output static analyses across benchmark systems")
    static_suite.add_argument("pdb_ids", nargs="*", help=f"PDB IDs to analyze; defaults to {', '.join(DEFAULT_BENCHMARK_IDS)}")
    static_suite.add_argument("--cache-dir", default=".crevice/pdb", help="Local structure cache directory")
    static_suite.add_argument("--out-dir", default="results/static-benchmark-suite", help="Output directory for suite results")
    static_suite.add_argument("--format", choices=["pdb", "cif", "mmcif"], default="cif", help="Structure format")
    static_suite.add_argument("--fetch", action="store_true", help="Download missing structures before analyzing")
    static_suite.add_argument("--samples", type=int, default=61, help="Profile samples for static analyses")
    static_suite.add_argument("--search-radius", type=float, default=3.0, help="Transverse search half-width (Å) of each channel profile")
    static_suite.add_argument("--refinement-steps", type=int, default=4, help="Profile center-refinement iterations")
    static_suite.add_argument("--probe-radius", action="append", type=float,
                              help="Measurement probe radius (Å) for the sensitivity rows; may be repeated (default "
                                   "0, 1.0 and 1.4). Rows keep the system's enclosure probe; a probe that does not "
                                   "pass the channel is recorded as an unresolved row")
    static_suite.add_argument("--enclosure-radius", type=_enclosure_radius_arg, default=None, metavar="{auto,Å}",
                              help="Enclosure probe (Å) of each system's channel profile. 'auto' (default) chooses it "
                                   "per system with the profile rule (smallest probe of the most persistent channel, "
                                   "0.5-3.0 Å), uses it for that system's sensitivity rows and records it in the "
                                   "summaries (enclosure_radius_A, enclosure_probe_mode); a number runs that probe "
                                   "(0.8 reproduces the suite outputs recorded before auto became the default)")
    static_suite.add_argument("--cavity-spacing", type=float, default=3.0, help="Cavity grid spacing in Å")
    static_suite.add_argument("--tunnel-spacing", type=float, default=4.0, help="Tunnel grid spacing in Å")
    static_suite.add_argument("--cast-spacing", type=float, default=1.5, help="Void-cast grid spacing in Å")
    static_suite.add_argument("--cast-min-radius", type=float, default=0.4, help="Minimum atom-surface clearance in Å for void-cast points")
    static_suite.add_argument("--cast-max-components", type=int, default=4, help="Maximum final connected components exported in each pore cast")
    static_suite.add_argument("--cast-max-grid-points", type=int, default=120000, help="Maximum grid points for void-cast analysis")
    static_suite.add_argument("--max-grid-points", type=int, default=60000, help="Maximum grid points for cavity/tunnel analyses")
    static_suite.add_argument("--residue-label-every", type=int, default=5, help="Landmark spacing in profile points for the residue-landmark radius figure")
    static_suite.add_argument("--dpi", type=int, default=400, help="Figure DPI")
    static_suite.add_argument("--skip-cavities", action="store_true", help="Skip cavity detection")
    static_suite.add_argument("--skip-tunnels", action="store_true", help="Skip tunnel detection")
    static_suite.add_argument("--skip-network", action="store_true", help="Skip residue-network analysis")
    static_suite.set_defaults(func=_cmd_static_suite)
    for command in [hydration, cavity_md, region_md, evidence, cast, publish, analyze, trajectory]:
        if command is not hydration:command.add_argument("--skip-hydration", action="store_true", help="Skip standard water-contact/SASA analysis")
        if command in (cavity_md, region_md):command.add_argument("--water-membership", action="store_true", help="Also count aligned water oxygens inside each frame's measured region, paired with the fixed-reference count")
        command.add_argument("--hydration-cutoff", type=float, default=3.5, help="Water oxygen–protein heavy-atom contact distance in Å")
        command.add_argument("--hydration-probe", type=float, default=1.4, help="Geometric SASA probe radius in Å")
        command.add_argument("--hydration-sasa-points", type=int, default=256, help="Equal-area quadrature points per atom for SASA")
        command.add_argument('--water-density-spacing', type=float, default=.5, help='Aligned water-density bin size in Å')
        command.add_argument('--water-density-smoothing', type=float, default=.5, help='Display-only Gaussian sigma in Å; zero keeps raw density')
        command.add_argument('--water-density-level', type=float, default=.05, help='Fixed display isovalue in water oxygens per Å³')
        command.add_argument('--skip-water-density', action='store_true', help="Skip the aligned water-oxygen density maps")
        command.add_argument('--skip-interactions', action='store_true', help='Skip typed chemical interactions and residue conformation measurements')
        command.add_argument('--skip-analysis-views', action='store_true', help='Keep numerical analyses but omit the unified molecular scenes, evidence table and overview figures')
    for command in [region_prepare, region_compare, region_md, hydration, cavity_md, cast, profile, residues, cavities,
                    network, evidence, trajectory, publish, analyze, static_suite]:
        _add_annotate_arg(command)
    for command in [cast, profile, trajectory, publish, analyze, static_suite]:
        _add_display_guide_args(command)
    for command in [region_prepare, region_md, hydration, cavity_md, cast, profile, residues, cavities, tunnels,
                    network, evidence, features, trajectory, publish, analyze, benchmark, static_suite]:
        _add_radii_arg(command)
    for command in [region_prepare, region_md, hydration, region_compare]:
        command.add_argument(
            "--allow-radii-mismatch", action="store_true",
            help="Proceed even though the radius set differs from the one recorded in the region definition "
                 "or earlier results (the mismatch is recorded). Without it a mismatch is an error")
    return parser


def _add_radii_arg(parser: argparse.ArgumentParser) -> None:
    """The atomic radius set used by every geometric measurement of the command."""
    parser.add_argument(
        "--radii", metavar="{default,bondi,hole,charmm_like,PATH}",
        help="Atomic radius set for every clearance, cast, contact and SASA radius. default (when omitted): "
             "the standard element table (crevice.radii.VDW_RADII: Bondi 1964, Mantina 2009) with CHARMM36 "
             "Rmin/2 radii for recognised monatomic ion residues (SOD, CLA, POT, NA, ZN...); bondi: the "
             "element table for every atom, ions included; hole: HOLE 2.3.1 simple.rad (AMBER united-atom "
             "radii, matched by atom name); charmm_like: larger C/N/O/S approximation plus the ion radii. "
             "PATH: a .json, .csv or HOLE .rad radius file. The set used is always recorded (name, source, "
             "file SHA-256, overrides, ion radii, table hash) in the output manifests/metadata as 'radii'. "
             "region-trajectory and hydration --cavity-results use the set recorded in the region "
             "definition or earlier results when this is omitted")


def _add_display_guide_args(parser: argparse.ArgumentParser) -> None:
    """Opt-in display guides: channel-mouth guides and exit-leg centre lines."""
    parser.add_argument(
        "--mouth-guides", action="store_true",
        help="Draw channel-mouth guides: orange rings in the PyMOL/VMD/ChimeraX scenes and orange mouth "
             "brackets/lines on profile plots. Off by default; mouth positions are always recorded in the "
             "profile JSON, scene JSON and manifests")
    parser.add_argument(
        "--exit-centre-lines", action="store_true",
        help="With --lateral-exits: also draw each lateral exit leg as a thin blue centre-line tube in the "
             "viewer scenes. Off by default; the legs are drawn as cast surfaces")


def _add_annotate_arg(parser: argparse.ArgumentParser) -> None:
    """The single opt-in for explanatory text in figures and viewer scenes."""
    parser.add_argument(
        "--annotate", action="store_true",
        help="Draw explanatory text in figures and viewer scenes: titles, captions, value call-outs, "
             "residue/landmark labels, viewer legends and 2D labels. Off by default (clean figures); "
             "the text is always recorded in PNG metadata and scene JSON")


def main(argv: list[str] | None = None) -> int:
    """Run the ``crevice`` command line.

    Parses ``argv``, validates the shared hydration settings, then runs the
    selected command inside :func:`crevice.presentation.figure_annotations`
    (``--annotate``), :func:`crevice.presentation.display_guides`
    (``--mouth-guides``, ``--exit-centre-lines``) and
    :func:`crevice.radii.use_radii` (``--radii``), so one setting applies to
    every figure, scene and atomic radius of the command.

    Parameters
    ----------
    argv : list of str, optional
        Arguments without the program name; default ``sys.argv[1:]``.

    Returns
    -------
    int
        0 on success, 2 when the command raised an error (the message is
        printed to standard error as ``crevice: error: ...``). Argument-parsing
        errors exit through :mod:`argparse` (status 2) and ``--help`` exits with 0.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if hasattr(args,'hydration_cutoff') and not getattr(args,'skip_hydration',False) and not getattr(args,'inspect',False):
            if not math.isfinite(args.hydration_cutoff) or args.hydration_cutoff<=0:raise ValueError('hydration cutoff must be finite and positive')
            if not math.isfinite(args.hydration_probe) or args.hydration_probe<=0:raise ValueError('hydration SASA probe must be finite and positive')
            if args.hydration_sasa_points<32:raise ValueError('hydration SASA requires at least 32 points per atom')
            if not all(math.isfinite(v) and v>0 for v in [args.water_density_spacing,args.water_density_level]):raise ValueError('Water density spacing and level must be finite and positive')
            if not math.isfinite(args.water_density_smoothing) or args.water_density_smoothing<0:raise ValueError('Water density smoothing must be finite and nonnegative')
        from .presentation import display_guides, figure_annotations
        from .radii import use_radii
        # One radius set for the whole command, including trajectory worker processes.
        with figure_annotations(getattr(args, "annotate", False)), use_radii(getattr(args, "radii", None)), \
                display_guides(mouth_guides=getattr(args, "mouth_guides", False),
                               exit_centre_lines=getattr(args, "exit_centre_lines", False)):
            args.func(args)
    except Exception as exc:
        print(f"crevice: error: {exc}", file=sys.stderr)
        return 2
    return 0


def _display_width(text: str) -> float:
    from .volume_export import validate_display_width
    try:
        return validate_display_width(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from None


def _add_smoothing_args(parser: argparse.ArgumentParser) -> None:
    """Display-surface smoothing options; measurements are never smoothed."""
    from .volume_export import DEFAULT_SMOOTH_A
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--smooth", type=_display_width, metavar="X",
                       help=f"Display-surface smoothing length X in Å (default {DEFAULT_SMOOTH_A:g}; 0 shows the raw voxel "
                            "boundary). Gaussian twicing G*(2B-G*B), sigma X, on a 2x supersampled display grid, clamped so every measured sample "
                            "keeps its inside/outside state and the surface moves only within boundary grid cells. "
                            "Measured grids, volumes and profiles are unchanged. See docs/SURFACE_SMOOTHING.md")
    group.add_argument("--surface-smoothing", type=_display_width, metavar="W",
                       help="Legacy display interpolation width in Å on the measured grid (sample-clamped Gaussian; "
                            "0 gives the raw voxel boundary). Cannot be combined with --smooth")


def _add_structure_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "structure",
        help="Local PDB/mmCIF file, an RCSB accession such as 1GRM, or an AlphaFold DB ID such as AF-P02920-F1. "
             "An existing path always wins over an accession of the same name.",
    )
    _add_input_args(parser)


def _add_input_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("structure input")
    group.add_argument("--cache-dir", default=".crevice/pdb",
                       help="Where fetched structures are cached")
    group.add_argument("--format", dest="input_format", choices=["pdb", "cif", "mmcif"],
                       default="cif", help="Format to fetch when given an accession")
    group.add_argument("--assembly", type=int,
                       help="Fetch this biological assembly instead of the deposited unit")
    group.add_argument("--offline", action="store_true",
                       help="Never download; use the cache or fail")
    group.add_argument("--refetch", action="store_true",
                       help="Re-download even when the structure is already cached")
    group.add_argument("--timeout", type=float, default=30.0,
                       help="Network timeout in seconds")
    group.add_argument("--model", type=int, default=0,
                       help="Model index for multi-model structures")
    group.add_argument("--parser", choices=list(PARSERS), default="auto",
                       help="mmCIF parser; auto prefers the validated gemmi loader. builtin marks HETATM "
                            "records (including modified polymer residues) as heteroatoms and keeps "
                            "zero-occupancy atoms; gemmi uses entity types and drops them")
    group.add_argument("--md-selection", default="protein",
                       help="MDAnalysis selection for local GRO files (default protein removes solvent, ions and lipids)")
    group.add_argument("--chain-ids", choices=["label", "auth"], default="label",
                       help="mmCIF chain-ID namespace (label_asym_id or auth_asym_id), for either parser")


def _resolve_input(args: argparse.Namespace, spec: str | None = None):
    """Resolve a structure argument to a local file, fetching if needed."""

    return resolve_structure(
        spec if spec is not None else args.structure,
        cache_dir=args.cache_dir, fmt=args.input_format, assembly=args.assembly,
        offline=args.offline, overwrite=args.refetch, timeout=args.timeout,
    )


def _load_input(args: argparse.Namespace, spec: str | None = None):
    """Resolve and load a structure, returning the frame and both reports."""

    source = _resolve_input(args, spec)
    frame, report = load_structure_report(
        source.path, model_index=args.model, parser=args.parser, chain_ids=args.chain_ids,
        md_selection=args.md_selection
    )
    if source.kind == "downloaded":
        print(f"fetched {source.accession} -> {source.path}")
    return frame, source, report


def _enclosure_radius_arg(value: str) -> float | None:
    """``--enclosure-radius``: ``auto`` (None, automatic choice) or a radius in Å."""
    if value.strip().lower() == "auto":
        return None
    try:
        radius = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("expected 'auto' or a radius in Å") from None
    if not math.isfinite(radius) or radius < 0:
        raise argparse.ArgumentTypeError("the enclosure radius must be finite and non-negative")
    return radius


def _add_profile_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--axis", default="auto", help="Channel direction: auto, x, y, z, or x,y,z vector")
    parser.add_argument("--origin", type=parse_coord, help="Seed inside the intended channel, x,y,z (Å)")
    parser.add_argument("--section-spacing", type=float, default=0.5, help="Transverse channel grid spacing (Å)")
    parser.add_argument("--enclosure-radius", type=_enclosure_radius_arg, default=None, metavar="{auto,Å}",
                        help="Probe radius (Å) that decides lateral enclosure and connectivity, separate from the "
                             "measured radius. 'auto' (default) tries 0.5-3.0 Å and uses the smallest probe of "
                             "the channel that resolves at the most probes (recorded in the metadata); a number "
                             "runs that probe only")
    parser.add_argument("--lateral-exits", action="store_true",
                        help="Accept a capped end (blocked straight axial exit) if a lateral exit to bulk exists; "
                             "every distinct exit leg is reported separately (see docs: profile, capped channels)")
    parser.add_argument("--exit-bulk-radius", type=float, default=6.0,
                        help="Rolling probe radius (Å) that defines bulk solvent: the outer boundary of the lateral "
                             "exit legs (--lateral-exits) and, in publish, of every ENTRY/EXIT cast segment. Not a width "
                             "limit: inside the boundary a segment fills all probe-accessible space. Exit legs whose exit "
                             "points are more than twice this apart are reported as distinct")
    parser.add_argument("--exit-spacing", type=float, default=0.5,
                        help="With --lateral-exits: lattice spacing (Å) of the exit-path search")
    parser.add_argument("--samples", type=int, default=81,
                        help="Minimum number of profile samples; connected profiles add points to keep axial "
                             "steps ≤ 0.75 Å")
    parser.add_argument("--padding", type=float, default=0.0, help="Axis padding (Å); only with --search-radius 0")
    parser.add_argument("--search-radius", type=float, default=6.0, help="Transverse search half-width (Å); 0 requests an unvalidated fixed-axis scan")
    parser.add_argument("--refinement-steps", type=int, default=4, help="Center refinement iterations")
    parser.add_argument("--probe-radius", type=float, default=0.0, help="Probe radius (Å) subtracted from each measured radius")
    parser.add_argument("--include-hydrogen", action="store_true", help="Include hydrogen atoms")
    parser.add_argument("--exclude-hetero", action="store_true", help="Exclude HETATM records")


def _profile_axis(value):
    return parse_coord(value) if isinstance(value, str) and "," in value else value


def _profile_from_args(args: argparse.Namespace):
    frame, _source, _report = _load_input(args)
    return frame, pore_profile(
        frame,
        axis=_profile_axis(args.axis),
        origin=args.origin,
        section_spacing=args.section_spacing,
        enclosure_radius=args.enclosure_radius,
        lateral_exits=args.lateral_exits,
        exit_bulk_radius=args.exit_bulk_radius,
        exit_spacing=args.exit_spacing,
        samples=args.samples,
        padding=args.padding,
        search_radius=args.search_radius,
        refinement_steps=args.refinement_steps,
        probe_radius=args.probe_radius,
        include_hydrogen=args.include_hydrogen,
        include_hetero=not args.exclude_hetero,
    )


def _cmd_profile(args: argparse.Namespace) -> None:
    frame, profile = _profile_from_args(args)
    write_profile_json(profile, args.output)
    write_profile_csv(profile, args.csv or _sibling(args.output, ".csv"))
    if args.pdb:
        write_profile_pdb(profile, args.pdb)
    if args.png:
        plot_profile_radius(profile, args.png, dpi=args.dpi)
    if args.annotated_png:
        plot_profile_radius_with_residues(
            profile,
            args.annotated_png,
            dpi=args.dpi,
            label_every=args.residue_label_every,
        )
    if args.pymol:
        if not args.pdb:
            raise ValueError("--pymol requires --pdb so the dummy atoms can be loaded")
        export_pymol_script(
            structure_path=frame.source,
            dummy_path=args.pdb,
            script_path=args.pymol,
            output_png=args.pymol_png,
            dpi=args.dpi,
        )
    print(
        "profile "
        f"points={len(profile.points)} min_radius={profile.min_radius:.3f} "
        f"mean_radius={profile.mean_radius:.3f} bottleneck_index={profile.bottleneck.index}"
    )
    selection = profile.metadata.get("enclosure_probe_selection")
    if selection:
        print(f"enclosure_radius={selection['chosen_A']:.3f} (auto: {selection['reason']})")


def _cmd_residues(args: argparse.Namespace) -> None:
    frame, profile = _profile_from_args(args)
    contacts = annotate_residues(
        frame,
        profile.points,
        cutoff=args.cutoff,
        include_hydrogen=args.include_hydrogen,
        include_hetero=not args.exclude_hetero,
    )
    write_residue_contacts_csv(contacts, args.output)
    if args.json:
        write_residue_contacts_json(contacts, args.json)
        _add_radii_provenance(args.json)
    else:
        _write_csv_provenance(args.output, "residues")
    if args.png:
        plot_residue_contacts(contacts, args.png, dpi=args.dpi)
    print(f"residues contacts={len(contacts)} output={args.output}")


def _cmd_cavities(args: argparse.Namespace) -> None:
    frame, _source, _report = _load_input(args)
    cavities = detect_cavities(
        frame,
        spacing=args.spacing,
        min_radius=args.min_radius,
        max_cavities=args.max_cavities,
        include_hydrogen=args.include_hydrogen,
        include_hetero=args.include_hetero,
    )
    write_cavities_json(cavities, args.output)
    _add_radii_provenance(args.output)
    write_cavities_csv(cavities, _sibling(args.output, ".csv"))
    if args.pdb:
        write_cavities_pdb(cavities, args.pdb)
    if args.png:
        plot_cavity_summary(cavities, args.png, dpi=args.dpi)
    print(f"cavities count={len(cavities)} output={args.output}")


def _cmd_tunnels(args: argparse.Namespace) -> None:
    frame, _source, _report = _load_input(args)
    start = parse_coord(args.start) if args.start else None
    tunnels = find_tunnels(
        frame,
        start=start,
        spacing=args.spacing,
        min_radius=args.min_radius,
        max_tunnels=args.max_tunnels,
        include_hydrogen=args.include_hydrogen,
        include_hetero=args.include_hetero,
    )
    write_tunnels_json(tunnels, args.output)
    _add_radii_provenance(args.output)
    write_tunnels_csv(tunnels, _sibling(args.output, ".csv"))
    write_tunnel_points_csv(tunnels, _sibling(args.output, "_points.csv"))
    if args.pdb:
        write_tunnels_pdb(tunnels, args.pdb)
    print(f"tunnels count={len(tunnels)} output={args.output}")


def _cmd_network(args: argparse.Namespace) -> None:
    if args.region=="profile":
        frame,profile=_profile_from_args(args)
    else:
        frame,_,_=_load_input(args)
        profile=None
    groups = json.loads(Path(args.residue_groups).read_text()) if args.residue_groups else None
    contacts = annotate_residues(
        frame,
        profile.points if profile else (),
        cutoff=args.cutoff,
        include_hydrogen=args.include_hydrogen,
        include_hetero=not args.exclude_hetero,
    )
    if args.region == "profile":
        network = build_cavity_network(
            frame,
            profile,
            contacts=contacts,
            residue_groups=groups,
            distance_metric=args.contact_metric,
            cutoff=args.cutoff,
            include_hydrogen=args.include_hydrogen,
            include_hetero=not args.exclude_hetero,
        )
    else:
        network = build_residue_network(
            frame,
            residue_groups=groups,
            distance_metric=args.contact_metric,
            cutoff=args.cutoff,
            include_hydrogen=args.include_hydrogen,
            include_hetero=not args.exclude_hetero,
        )
    from .radii import radii_fields
    data = {"network": network.to_dict(), "metrics": network_metrics(network), **radii_fields()}
    _write_json(data, args.output)
    write_network_csv(network, _sibling(args.output, "_nodes.csv"), _sibling(args.output, "_edges.csv"))
    if args.connectivity_json:
        report = lining_connectivity(network)
        write_connectivity_csv(report, _sibling(args.connectivity_json, ".csv"),
                               _sibling(args.connectivity_json, "_groups.csv"))
        if args.remove_residue:
            report["after_removal"] = lining_connectivity(network, removed_residues=args.remove_residue)
            write_connectivity_csv(report["after_removal"], _sibling(args.connectivity_json, "_after_removal.csv"),
                                   _sibling(args.connectivity_json, "_after_removal_groups.csv"))
        from .radii import radii_fields
        _write_json({**report, **radii_fields()}, args.connectivity_json)
    if args.connectivity_png:
        plot_lining_connectivity(network, args.connectivity_png, dpi=args.dpi)
    if args.png:
        plot_network_summary(network, args.png, dpi=args.dpi)
    if args.chord_png:
        plot_network_chord(network, args.chord_png, dpi=args.dpi)
    print(f"network nodes={len(network.nodes)} edges={len(network.edges)} output={args.output}")


def _cmd_features(args: argparse.Namespace) -> None:
    frame, profile = _profile_from_args(args)
    contacts = annotate_residues(
        frame,
        profile.points,
        cutoff=args.cutoff,
        include_hydrogen=args.include_hydrogen,
        include_hetero=not args.exclude_hetero,
    )
    data = {
        "profile": profile.to_dict(),
        "features": profile_features(profile, contacts),
        "top_residues": [contact.to_dict() for contact in contacts[:25]],
    }
    from .radii import radii_fields
    data.update(radii_fields())
    _write_json(data, args.output)
    write_features_csv(data["features"], _sibling(args.output, ".csv"))
    print(f"features output={args.output}")


def _cmd_trajectory(args: argparse.Namespace) -> None:
    if not 0<=args.confidence<1:
        raise ValueError("--confidence must be in [0,1)")
    if not args.align and args.confidence and (args.png or args.distribution_csv):
        raise ValueError("Confidence bands require --align; use --confidence 0 for explicitly unaligned descriptive outputs")
    if args.topology and len(args.structures) != 1:
        raise ValueError("--topology requires exactly one trajectory file")
    if args.inspect:
        if not args.topology:
            raise ValueError("--inspect describes an MD trajectory; pass --topology")
        info = inspect_trajectory(args.topology, args.structures[0], selection=args.selection)
        for key, value in info.items():
            print(f"{key}\t{value}")
        if args.report_json:
            _write_json(info, args.report_json)
        return
    if not args.output:
        raise ValueError("trajectory requires -o/--output (only --inspect prints without one)")
    reader_report = None
    if args.topology:
        trajectory, reader_report = load_trajectory_report(
            args.topology, args.structures[0], selection=args.selection,
            start=args.start, stop=args.stop, stride=args.stride,
            max_frames=args.max_frames or None,pbc=args.pbc,
        )
        print(f"read {reader_report['loaded_frame_count']} of "
              f"{reader_report['available_frame_count']} frames, "
              f"{reader_report['selected_atom_count']} atoms per frame")
    else:
        trajectory = load_static_frames(args.structures)
    if args.report_json:
        _write_json(reader_report or {"reader": "static-file-sequence",
                                      "frames": [str(s) for s in args.structures]},
                    args.report_json)
    analyses = tuple(args.analysis) if args.analysis else ("profile", "features")
    if args.network_json:
        analyses = tuple(dict.fromkeys((*analyses, "network")))
    groups = json.loads(Path(args.residue_groups).read_text()) if args.residue_groups else None
    if args.alignment_residue and args.alignment_residues_json:
        raise ValueError("Choose --alignment-residue or --alignment-residues-json")
    alignment_residues=json.loads(Path(args.alignment_residues_json).read_text()) if args.alignment_residues_json else args.alignment_residue
    if alignment_residues is not None and (not isinstance(alignment_residues,list) or not all(isinstance(r,str) for r in alignment_residues)):
        raise ValueError("Alignment scaffold must be a JSON list of full residue IDs")
    lining=json.loads(Path(args.lining_evidence).read_text())["lining_residues"] if args.lining_evidence else ()
    result = analyze_trajectory(
        trajectory,
        analyses=analyses,
        stride=1 if args.topology else args.stride,
        align=args.align,alignment_residues=alignment_residues,
        on_profile_error=args.profile_errors,retain_networks=not bool(args.network_json),
        lining_residues=lining,network_occupancy=args.occupancy,
        network_kwargs={"residue_groups":groups,"distance_metric":args.contact_metric,"cutoff":args.contact_cutoff},
        profile_kwargs={"axis": _profile_axis(args.axis), "origin": args.origin, "samples": args.samples,
                        "search_radius": args.search_radius, "section_spacing": args.section_spacing,
                        "enclosure_radius": args.enclosure_radius},
        probe_fallback=args.probe_fallback,
    )
    probe = result.metadata.get("enclosure_probe")
    if probe and probe.get("mode") == "auto":
        if probe["chosen_A"] is None:
            print(f"enclosure_radius=none (auto: no probe could be chosen on the reference frame, so no frame was "
                  f"profiled: {probe['reason']})")
        else:
            print(f"enclosure_radius={probe['chosen_A']:.3f} (auto, chosen once on the reference frame and fixed "
                  f"for every frame: {probe['reason']})")
    fallback = (probe or {}).get("fallback")
    if fallback:
        print(f"enclosure probe per frame: pinned={fallback['pinned_frames']} fallback={fallback['fallback_frames']} "
              f"unresolved={fallback['unresolved_frames']} (fallback {'on' if fallback['enabled'] else 'off'})"
              + (f"; fallback probes by source frame: {fallback['fallback_probes_A']}" if fallback['fallback_frames'] else ""))
    write_trajectory_json(result, args.output)
    from .io import write_trajectory_csv
    write_trajectory_csv(result, _sibling(args.output, "_frames.csv"), _sibling(args.output, "_profiles.csv"))
    distribution=None
    if args.png or args.distribution_csv:
        if any(f.profile is not None for f in result.frames):
            distribution=profile_distribution(result,quantiles=args.quantiles,confidence=args.confidence or None,
                     block_length=args.block_length,bootstrap_replicates=args.bootstrap_replicates,seed=args.bootstrap_seed)
        else:
            distribution={"status":"no_profiles","rows":[],"total_frame_count":len(result.frames),
                          "confidence_interval":{"status":"unavailable_no_reference_profile"},
                          "reason":result.metadata.get("reference_profile_reason")}
    if args.png:
        plot_trajectory_profiles(result,args.png,dpi=args.dpi,view=args.view,quantiles=args.quantiles,
                                 distribution=distribution,quantity=args.profile_quantity)
    if args.distribution_csv:
        write_profile_distribution_csv(distribution,args.distribution_csv)
        from .radii import radii_fields
        _write_json({**distribution,**radii_fields()},str(Path(args.distribution_csv).with_suffix(".json")))
    if args.network_json:
        from .radii import radii_fields
        _write_json({**result.metadata['residue_dynamics'],**radii_fields()},args.network_json)
        from .residue_dynamics import write_dynamics_summary
        write_dynamics_summary(result.metadata['residue_dynamics'],Path(args.network_json).with_suffix(''),dpi=args.dpi)
    if not args.skip_hydration:
        from .hydration_workflow import static_hydration_bundle,hydration_dependencies,analysis_options_from_args
        missing_hydration=hydration_dependencies()
        focus=sorted(set(lining)|{c.residue for f in result.frames for c in f.residue_contacts})
        output=Path(args.output);hydration_files={}
        hydration_times=[f.time for f in result.frames]
        invalid_hydration_times=any(t is None or not math.isfinite(t) for t in hydration_times) or any(b<=a for a,b in zip(hydration_times,hydration_times[1:]) if a is not None and b is not None)
        if missing_hydration:
            hydration_files={"status":"unavailable_missing_required_dependencies","missing":missing_hydration,"install":"reinstall CREVICE (pip install --force-reinstall crevice)"}
        elif args.topology and focus and invalid_hydration_times:
            hydration_files={"status":"unavailable_invalid_timestamps","reason":"Hydration dynamics requires actual finite increasing frame times; other analyses are retained"}
        elif args.topology and focus:
            from .hydration_trajectory import analyze_hydration_trajectory
            from .hydration_export import write_hydration_bundle
            hydrated=analyze_hydration_trajectory(args.topology,args.structures[0],residue_ids=focus,selection=args.selection,
                start=args.start,stop=args.stop,stride=args.stride,max_frames=args.max_frames or None,pbc=args.pbc,alignment_residues=alignment_residues,align=args.align,
                cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,analysis_options=analysis_options_from_args(args))
            hydration_files=write_hydration_bundle(hydrated,output.parent,prefix=output.stem,dpi=args.dpi,confidence=args.confidence,block_length=args.block_length,replicates=args.bootstrap_replicates,seed=args.bootstrap_seed)
        elif not args.topology:
            for i,frame in enumerate(trajectory.frames[::args.stride]):
                hydration_files[str(i)]=static_hydration_bundle(frame,output.parent/(output.stem+'_hydration_frames'),prefix=f'frame_{i:04d}',
                    residue_ids=focus,cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,dpi=args.dpi,analysis_options=analysis_options_from_args(args))
        else:
            hydration_files={'status':'unavailable_no_identified_residues','reason':'Supply lining evidence or use the cavity-trajectory workflow'}
        result.metadata['hydration_files']=hydration_files;write_trajectory_json(result,args.output)
    print(f"trajectory frames={len(result.frames)} analyses={','.join(result.analyses)} output={args.output}")


def _cmd_cast(args: argparse.Namespace) -> None:
    from .figures import write_rolling_probe_bundle
    frame, source, report = _load_input(args)
    prefix = args.prefix or source.path.stem
    manifest = write_rolling_probe_bundle(frame, structure_path=source.path, output_dir=args.out_dir,
        prefix=prefix, spacing=args.spacing, probe_radius=args.probe_radius, outer_radius=args.outer_radius,
        enclosure_fraction=args.enclosure_fraction, min_depth=args.min_depth,
        min_component_volume=args.min_component_volume, max_components=args.max_components,
        max_grid_points=args.max_grid_points, selection=args.selection, dpi=args.dpi,
        selection_mode=args.cavity_selection, core_fraction=args.core_fraction, surface_smoothing=args.surface_smoothing, smooth=args.smooth,
        seed=args.seed, seed_tolerance=args.seed_tolerance)
    path = Path(args.out_dir)/f'{prefix}_input_report.json'
    _write_json({'source': source.to_dict(), 'intake': report}, path)
    saved = json.loads(Path(manifest['manifest_json']).read_text())
    saved['files']['input_report_json'] = str(path)
    _write_json(saved, manifest['manifest_json'])
    from .hydration_workflow import append_static_hydration
    append_static_hydration(args,frame,source.path,manifest)
    cast = saved['void_cast']
    print(f"cast regions={cast['component_count']} volume_A3={cast['total_volume']:.3f} spacing_A={cast['spacing']} out_dir={args.out_dir}")


def _cmd_publish(args: argparse.Namespace) -> None:
    frame, source, _report = _load_input(args)
    prefix = args.prefix or source.path.stem
    manifest = write_static_publication_bundle(
        frame,
        structure_path=source.path,
        output_dir=args.out_dir,
        prefix=prefix,
        axis=_profile_axis(args.axis),
        origin=args.origin,
        section_spacing=args.section_spacing,
        enclosure_radius=args.enclosure_radius,
        lateral_exits=args.lateral_exits,
        exit_bulk_radius=args.exit_bulk_radius,
        exit_spacing=args.exit_spacing,
        samples=args.samples,
        search_radius=args.search_radius,
        refinement_steps=args.refinement_steps,
        probe_radius=args.probe_radius,
        include_hydrogen=args.include_hydrogen,
        include_hetero=not args.exclude_hetero,
        dpi=args.dpi,
        include_cavities=not args.skip_cavities,
        include_tunnels=not args.skip_tunnels,
        include_network=not args.skip_network,
        residue_label_every=args.residue_label_every,
        cast_mode=args.cast_mode,
        cast_outer_radius=args.cast_outer_radius,
        cast_enclosure_fraction=args.cast_enclosure_fraction,
        cast_extension=args.cast_extension,
        cast_spacing=args.cast_spacing,
        cast_min_radius=args.cast_min_radius,
        cast_max_components=args.cast_max_components,
        cast_selection=args.cast_selection, cast_core_fraction=args.cast_core_fraction,
        surface_smoothing=args.surface_smoothing,
        smooth=args.smooth,
        cast_max_grid_points=args.cast_max_grid_points,
        cast_focus_points=args.focus,
        cast_focus_radius=args.focus_radius,
        cast_min_component_volume=args.min_component_volume,
        cast_max_component_volume=args.max_component_volume,
        residue_groups=json.loads(Path(args.residue_groups).read_text()) if args.residue_groups else None,
        entry_end=args.entry_end,
        lining_cutoff=args.lining_cutoff,
    )
    input_report_path = Path(args.out_dir) / f"{prefix}_input_report.json"
    _write_json({"source": source.to_dict(), "intake": _report}, input_report_path)
    saved_manifest = json.loads(Path(manifest["manifest_json"]).read_text())
    saved_manifest["files"]["input_report_json"] = str(input_report_path)
    _write_json(saved_manifest, manifest["manifest_json"])
    from .hydration_workflow import append_static_hydration
    append_static_hydration(args,frame,source.path,manifest)
    profile_status = "unresolved" if "profile_status_json" in manifest else "resolved" if "profile_json" in manifest else "not_assigned"
    print(f"publish files={len(manifest)} out_dir={args.out_dir} profile={profile_status}")


def _cmd_fetch(args: argparse.Namespace) -> None:
    if args.list:
        for system in benchmark_systems():
            print(
                f"{system.pdb_id}\t{system.expected_geometry}\t{system.curation_status}\t"
                f"{system.name}\t{system.notes}"
            )
        return

    pdb_ids = tuple(args.pdb_ids) if args.pdb_ids else DEFAULT_BENCHMARK_IDS
    records = [
        fetch_structure(
            pdb_id,
            args.cache_dir,
            fmt=args.format,
            assembly=args.assembly,
            overwrite=args.overwrite,
            timeout=args.timeout,
        )
        for pdb_id in pdb_ids
    ]
    for record in records:
        identity = "verified" if record.identity_verified else f"unverified ({record.identity_note})"
        print(f"{record.path}\t{record.source}\t{record.byte_count} bytes\t{identity}")
    if args.json:
        _write_json({"downloads": [record.to_dict() for record in records]}, args.json)


def _cmd_fetch_example(args: argparse.Namespace) -> None:
    from .examples import fetch_example, list_examples
    if args.list or not args.name:
        for dataset in list_examples():
            print(f"{dataset.name}\t{dataset.description}\t{dataset.licence}")
            for item in dataset.files:
                print(f"  {item.name}\t{item.size} bytes\t{item.description}")
        return
    fetched = fetch_example(args.name, cache_dir=args.cache_dir, base_url=args.base_url, source_dir=args.source_dir,
                            overwrite=args.overwrite, timeout=args.timeout)
    for name, path in fetched.paths.items():
        state = "fetched" if name in fetched.downloaded else "cached"
        print(f"{path}\t{state}\tsha256 verified")
    print(f"example directory={fetched.directory}")


def _cmd_analyze(args: argparse.Namespace) -> None:
    frame, source, parse_report = _load_input(args)
    analyses = tuple(dict.fromkeys(args.analysis or ("profile", "residues")))
    out_dir = Path(args.out_dir)
    prefix = args.prefix or source.path.stem
    files: dict[str, str] = {}

    def target(name: str) -> str:
        path = str(out_dir / f"{prefix}_{name}")
        files[name] = path
        return path

    profile = None
    contacts: tuple = ()
    if {"profile", "residues", "network", "features"} & set(analyses):
        profile = pore_profile(
            frame, axis=_profile_axis(args.axis), origin=args.origin, samples=args.samples, padding=args.padding,
            section_spacing=args.section_spacing, enclosure_radius=args.enclosure_radius,
            lateral_exits=args.lateral_exits, exit_bulk_radius=args.exit_bulk_radius,
            exit_spacing=args.exit_spacing,
            search_radius=args.search_radius, refinement_steps=args.refinement_steps,
            probe_radius=args.probe_radius, include_hydrogen=args.include_hydrogen,
            include_hetero=not args.exclude_hetero,
        )
    if {"residues", "network", "features"} & set(analyses):
        contacts = annotate_residues(frame, profile.points, cutoff=args.cutoff)

    if "profile" in analyses:
        write_profile_json(profile, target("profile.json"))
        write_profile_csv(profile, target("profile.csv"))
        if args.png:
            plot_profile_radius(profile, target("profile.png"), dpi=args.dpi)
    if "residues" in analyses:
        write_residue_contacts_csv(contacts, target("residues.csv"))
        if args.png:
            plot_residue_contacts(contacts, target("residues.png"), dpi=args.dpi)
    if "cavities" in analyses:
        cavities = detect_cavities(
            frame, spacing=args.cavity_spacing, min_radius=args.cavity_min_radius,
            max_cavities=args.max_cavities, max_grid_points=args.cavity_max_grid_points,
            include_hydrogen=args.include_hydrogen,
            include_hetero=not args.exclude_hetero,
        )
        write_cavities_json(cavities, target("cavities.json"))
        write_cavities_csv(cavities, target("cavities.csv"))
        if args.png:
            plot_cavity_summary(cavities, target("cavities.png"), dpi=args.dpi)
    if "tunnels" in analyses:
        tunnels = find_tunnels(
            frame, start=parse_coord(args.start) if args.start else None,
            spacing=args.tunnel_spacing, min_radius=args.tunnel_min_radius,
            max_tunnels=args.max_tunnels, max_grid_points=args.tunnel_max_grid_points,
            include_hydrogen=args.include_hydrogen,
            include_hetero=not args.exclude_hetero,
        )
        write_tunnels_json(tunnels, target("tunnels.json"))
        write_tunnels_csv(tunnels, target("tunnels.csv"))
        write_tunnel_points_csv(tunnels, target("tunnel_points.csv"))
    if "network" in analyses:
        network = build_cavity_network(frame, profile, contacts=contacts, cutoff=args.cutoff)
        write_network_json(network, target("network.json"))
        write_network_csv(network, target("network_nodes.csv"), target("network_edges.csv"))
        if args.png:
            plot_network_summary(network, target("network.png"), dpi=args.dpi)
    if "features" in analyses:
        write_features_csv(profile_features(profile, contacts), target("features.csv"))

    from .hydration_workflow import append_static_hydration
    append_static_hydration(args,frame,source.path,files,residue_ids=[c.residue for c in contacts])
    manifest = {
        "analyses": list(analyses),
        "files": files,
        "source": source.to_dict(),
        "parser_report": parse_report,
        "interpretation": "Geometric outputs at the stated parameters; no curated "
                          "axis, assembly or biological validation is implied.",
    }
    from .radii import radii_fields
    manifest.update(radii_fields())
    _write_json(manifest, str(out_dir / f"{prefix}_manifest.json"))
    print(f"analyze analyses={','.join(analyses)} files={len(files)} out_dir={args.out_dir}")


def _cmd_static_suite(args: argparse.Namespace) -> None:
    result = run_static_benchmark_suite(
        cache_dir=args.cache_dir,
        output_dir=args.out_dir,
        pdb_ids=tuple(args.pdb_ids) if args.pdb_ids else DEFAULT_BENCHMARK_IDS,
        extension=args.format,
        fetch=args.fetch,
        profile_samples=args.samples,
        search_radius=args.search_radius,
        refinement_steps=args.refinement_steps,
        probe_radii=tuple(args.probe_radius) if args.probe_radius else (0.0, 1.0, 1.4),
        residue_label_every=args.residue_label_every,
        cavity_spacing=args.cavity_spacing,
        tunnel_spacing=args.tunnel_spacing,
        cast_spacing=args.cast_spacing,
        cast_min_radius=args.cast_min_radius,
        cast_max_components=args.cast_max_components,
        cast_max_grid_points=args.cast_max_grid_points,
        max_grid_points=args.max_grid_points,
        include_cavities=not args.skip_cavities,
        include_tunnels=not args.skip_tunnels,
        include_network=not args.skip_network,
        dpi=args.dpi,
        enclosure_radius=args.enclosure_radius,
    )
    ok = sum(system.get("status") == "ok" for system in result["systems"])
    print(f"static-suite ok={ok}/{len(result['systems'])} summary={result['summary_json']}")


def _cmd_benchmark(args: argparse.Namespace) -> None:
    pdb_ids = tuple(pdb_id.upper() for pdb_id in args.pdb_ids) if args.pdb_ids else DEFAULT_BENCHMARK_IDS
    rows = []
    for pdb_id in pdb_ids:
        system = get_benchmark_system(pdb_id)
        path = benchmark_path(args.cache_dir, system.pdb_id, extension=args.format)
        row = {
            "pdb_id": system.pdb_id,
            "expected_geometry": system.expected_geometry,
            "curation_status": system.curation_status,
            "path": str(path),
            "status": "pending",
        }
        try:
            if not path.exists():
                if not args.fetch:
                    raise FileNotFoundError(f"{path} is missing; rerun with --fetch to download")
                path = fetch_benchmark_structure(system.pdb_id, args.cache_dir, extension=args.format)
                row["path"] = str(path)
            frame = load_structure_report(path)[0]
            profile = pore_profile(frame, axis="auto", samples=args.samples, search_radius=args.search_radius,
                                   enclosure_radius=args.enclosure_radius)
            selection = profile.metadata.get("enclosure_probe_selection")
            if profile.method != "axial-connected":
                row.update({"enclosure_radius": None, "enclosure_probe_mode": "not_used",
                            "enclosure_probe_reason": "search_radius=0 fixed-axis scan does not use an enclosure probe"})
            else:
                row.update({"enclosure_radius": selection["chosen_A"] if selection else args.enclosure_radius,
                            "enclosure_probe_mode": "auto" if selection else "explicit",
                            "enclosure_probe_reason": selection["reason"] if selection else
                            "enclosure_radius given explicitly"})
            contacts = annotate_residues(frame, profile.points)
            row.update(
                {
                    "status": "ok",
                    "atom_count": len(frame.atoms),
                    "residue_count": len(frame.residues()),
                    "profile_min_radius": profile.min_radius,
                    "profile_mean_radius": profile.mean_radius,
                    "profile_bottleneck_index": profile.bottleneck.index,
                    "residue_contact_count": len(contacts),
                }
            )
            if system.expected_geometry == "cavity":
                cavities = detect_cavities(frame, spacing=2.5, max_cavities=5)
                row["cavity_count"] = len(cavities)
                row["top_cavity_volume"] = cavities[0].volume if cavities else 0.0
        except Exception as exc:  # benchmark output should capture per-system failures
            row.update({"status": "error", "error": str(exc)})
            selection = getattr(exc, "probe_selection", None)
            if selection is not None:
                row.update({"enclosure_radius": None, "enclosure_probe_mode": "auto",
                            "enclosure_probe_reason": selection["reason"]})
        rows.append(row)
        print(f"{row['pdb_id']}\t{row['status']}\t{row.get('atom_count', '-')} atoms\t{row.get('error', '')}")
    from .radii import radii_fields
    summary = {"benchmarks": rows, "ok": sum(row["status"] == "ok" for row in rows), "total": len(rows),
               **radii_fields()}
    if args.output:
        _write_json(summary, args.output)
        _write_benchmark_csv(rows, _sibling(args.output, ".csv"))


def _add_radii_provenance(path: str | Path) -> None:
    """Add ``radii`` provenance (the radius set in effect, default included) to a JSON output.

    The file is re-written with the same ``indent=2, sort_keys=True`` layout
    as the writers use, so only the ``radii`` key is added.
    """
    from .radii import radii_fields
    target = Path(path)
    data = json.loads(target.read_text())
    data.update(radii_fields())
    target.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def _write_csv_provenance(csv_path: str | Path, command: str) -> None:
    """Write ``<csv stem>_provenance.json`` beside a CSV that has no JSON companion.

    It records the command, the CSV file name and ``radii`` (the radius set in
    effect), so a CSV-only output still says which atomic radii produced it.
    """
    from .radii import radii_fields
    target = Path(csv_path)
    _write_json({"command": command, "csv": target.name, **radii_fields()},
                target.with_name(target.stem + "_provenance.json"))


def _sibling(path: str | Path, suffix: str) -> Path:
    """``out.json`` -> ``out<suffix>``: the CSV companion of a named JSON output."""
    target = Path(path)
    return target.with_name(target.stem + suffix) if target.suffix else target.with_name(target.name + suffix)


BENCHMARK_COLUMNS = {"pdb_id": "pdb_id", "expected_geometry": "expected_geometry", "curation_status": "curation_status",
                     "path": "path", "status": "status", "atom_count": "atom_count", "residue_count": "residue_count",
                     "profile_min_radius": "profile_min_radius_A", "profile_mean_radius": "profile_mean_radius_A",
                     "profile_bottleneck_index": "profile_bottleneck_index",
                     "enclosure_radius": "enclosure_radius_A", "enclosure_probe_mode": "enclosure_probe_mode",
                     "enclosure_probe_reason": "enclosure_probe_reason",
                     "residue_contact_count": "residue_contact_count", "cavity_count": "cavity_count",
                     "top_cavity_volume": "top_cavity_volume_A3", "error": "error"}


def _write_benchmark_csv(rows: list[dict], path: Path) -> None:
    """One row per benchmark system; radii in Å (``_A``), volume in Å³ (``_A3``)."""
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(BENCHMARK_COLUMNS.values()))
        writer.writeheader()
        writer.writerows({column: row.get(key, "") for key, column in BENCHMARK_COLUMNS.items()} for row in rows)


def _write_json(data: object, path: str | Path) -> None:
    target = Path(path)
    if target.parent != Path("."):
        target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def _cmd_residue_evidence(args):
    from .residue_evidence import read_binary_dx,boundary_residue_evidence,write_residue_evidence_bundle,restore_viewer_identities
    frame,source,_=_load_input(args)
    frame,identity_provenance=restore_viewer_identities(frame,args.structure)
    grid,origin,deltas=read_binary_dx(args.volume_dx)
    groups=json.loads(Path(args.residue_groups).read_text()) if args.residue_groups else None
    report=boundary_residue_evidence(frame,grid,origin,deltas,lining_distance=args.lining_distance,
                                    contact_cutoff=args.contact_cutoff,residue_groups=groups)
    report['source_volume_dx']=str(Path(args.volume_dx).resolve())
    report['residue_identity_provenance']=identity_provenance
    files=write_residue_evidence_bundle(report,frame,args.out_dir,prefix=args.prefix,scene_path=args.scene,dpi=args.dpi,
                                        stick_residues=args.stick_residues)
    from .hydration_workflow import append_static_hydration
    append_static_hydration(args,frame,source.path,files,residue_ids=[r['residue'] for r in report['residues'] if r['role']!='other'],roles={r['residue']:r['role'] for r in report['residues']})
    from .radii import radii_fields
    _write_json({'files':files,**radii_fields()},Path(args.out_dir)/f'{args.prefix}_evidence_manifest.json')
    print(f"residue-evidence lining={len(report['lining_residues'])} out_dir={args.out_dir}")


def _cmd_cavity_trajectory(args):
    from .cavity_trajectory import CavityReference,analyze_cavity_trajectory,write_cavity_trajectory_bundle
    if args.obstacle_selection and args.geometry_mode!="reference-region":
        raise ValueError("--obstacle-selection requires --geometry-mode reference-region")
    trajectory,reader=load_trajectory_report(args.topology,args.trajectory,selection=args.selection,
                                            stop=args.stop,max_frames=args.max_frames,pbc=args.pbc)
    root=Path(args.out_dir);root.mkdir(parents=True,exist_ok=True)
    _write_json(reader,root/(args.prefix+'_reader.json'))
    reference=CavityReference.from_dx(args.reference_volume_dx,spacing=args.spacing,axis=args.axis,profile_padding=args.profile_padding,grid_phase=args.grid_phase)
    scaffold=json.loads(Path(args.alignment_residues_json).read_text()) if args.alignment_residues_json else None
    geometry=dict(probe_radius=args.probe_radius,outer_radius=args.outer_radius,selection_mode=args.cavity_selection,
                  enclosure_fraction=args.enclosure_fraction,min_component_volume=args.min_component_volume,max_grid_points=args.max_grid_points)
    obstacles=None
    if args.obstacle_selection:
        from .cavity_obstacles import ObstacleSource
        obstacles=ObstacleSource(args.topology,args.trajectory,protein_selection=args.selection,obstacle_selection=args.obstacle_selection,
            reference_frame=trajectory.frames[0],frame_times=dict(zip(range(*reader["frame_indices"]),trajectory.times())),pbc=args.pbc)
    waters=None
    if getattr(args,'water_membership',False):
        from .water_membership import WaterSource
        waters=WaterSource(args.topology,args.trajectory,protein_selection=args.selection,reference_frame=trajectory.frames[0],
            frame_times=dict(zip(range(*reader["frame_indices"]),trajectory.times())),pbc=args.pbc)
    completed=0
    def progress(index,row):
        nonlocal completed
        completed+=1
        _write_json({'completed_frames':completed,'total_frames':len(trajectory.frames),
                     'last_source_index':index,'last_status':row['status']},root/(args.prefix+'_progress.json'))
        if completed%25==0:print(f"cavity frames completed={completed}/{len(trajectory.frames)}",flush=True)
    analysis=analyze_cavity_trajectory(trajectory,reference,geometry=geometry,alignment_residues=scaffold,
              minimum_overlap=args.minimum_overlap,ambiguity_ratio=args.ambiguity_ratio,lining_distance=args.lining_distance,
              contact_cutoff=args.contact_cutoff,workers=args.workers,progress_callback=progress,geometry_mode=args.geometry_mode,region_margin=args.region_margin,obstacles=obstacles,waters=waters)
    if getattr(args,'region_definition_metadata',None):
        analysis['settings']['region_definition']=args.region_definition_metadata
    files=write_cavity_trajectory_bundle(analysis,root,prefix=args.prefix,confidence=args.confidence,
             block_length=args.block_length,replicates=args.bootstrap_replicates,seed=args.bootstrap_seed,dpi=args.dpi)
    if not args.skip_hydration:
        from .hydration_trajectory import analyze_hydration_trajectory
        from .hydration_workflow import analysis_options_from_args
        from .hydration_export import write_hydration_bundle
        hydration=analyze_hydration_trajectory(args.topology,args.trajectory,context_json=files['statistics_json'],max_frames=args.max_frames,pbc=args.pbc,
           cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,analysis_options=analysis_options_from_args(args),
           progress_callback=lambda i,n: print(f"hydration frames completed={i+1}/{n}",flush=True) if (i+1)%25==0 else None)
        if waters is not None:
            from .water_membership import viewer_member_waters
            hydration['member_waters']=viewer_member_waters(analysis)
        extra=write_hydration_bundle(hydration,root,prefix=args.prefix,confidence=args.confidence,block_length=args.block_length,replicates=args.bootstrap_replicates,seed=args.bootstrap_seed,dpi=args.dpi)
        manifest=Path(files['manifest_json']);data=json.loads(manifest.read_text());data['files'].update(extra);_write_json(data,manifest)
    if waters is not None:
        from .water_membership import write_water_membership_bundle
        import numpy as np
        counts=None
        if not args.skip_hydration and extra.get('hydration_frames_npz'):
            with np.load(extra['hydration_frames_npz']) as saved:counts=saved['region_water_count'].copy() if 'region_water_count' in saved else None
        extra=write_water_membership_bundle(analysis,root,prefix=args.prefix,confidence=args.confidence,block_length=args.block_length,
            replicates=args.bootstrap_replicates,seed=args.bootstrap_seed,dpi=args.dpi,hydration_region_counts=counts)
        manifest=Path(files['manifest_json']);data=json.loads(manifest.read_text());data['files'].update(extra);_write_json(data,manifest)
    print(f"cavity trajectory frames={len(trajectory.frames)} statistics={files['statistics_json']}")


def _cmd_hydration(args):
    from .hydration_workflow import static_hydration_bundle,hydration_dependencies,analysis_options_from_args
    if hydration_dependencies():raise RuntimeError("Hydration needs NumPy, SciPy and Matplotlib, which are required CREVICE dependencies but could not be imported: reinstall CREVICE (pip install --force-reinstall crevice)")
    from .hydration_trajectory import analyze_hydration_trajectory
    from .hydration_export import write_hydration_bundle
    residues=json.loads(Path(args.residues_json).read_text()) if args.residues_json else None
    if args.trajectory:
        result=analyze_hydration_trajectory(args.structure,args.trajectory,context_json=args.cavity_results,residue_ids=residues,
            allow_radii_mismatch=args.allow_radii_mismatch,
            selection=args.selection,water_selection=args.water_selection,start=args.start,stop=args.stop,stride=args.stride,
            max_frames=args.max_frames,pbc=args.pbc,cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,analysis_options=analysis_options_from_args(args),
            progress_callback=lambda i,n: print(f"hydration frames completed={i+1}/{n}",flush=True) if (i+1)%25==0 else None)
        files=write_hydration_bundle(result,args.out_dir,prefix=args.prefix,dpi=args.dpi,confidence=args.confidence,block_length=args.block_length,replicates=args.bootstrap_replicates,seed=args.bootstrap_seed)
    else:
        if args.cavity_results:raise ValueError('--cavity-results requires --trajectory')
        frame,source,_=_load_input(args)
        from .residue_evidence import restore_viewer_identities
        frame,_=restore_viewer_identities(frame,source.path)
        files=static_hydration_bundle(frame,args.out_dir,prefix=args.prefix,source_path=source.path,volume_path=args.volume_dx,residue_ids=residues,
            cutoff=args.hydration_cutoff,sasa_probe=args.hydration_probe,sasa_samples=args.hydration_sasa_points,dpi=args.dpi,analysis_options=analysis_options_from_args(args))
    print('hydration statistics='+files['hydration_json'])


def _cmd_region_init(args):
    from .region_definition import file_record,validate_definition,write_json
    target=Path(args.output).resolve()
    if target.exists():raise ValueError('Definition exists; use a new versioned filename')
    landmarks=json.loads(Path(args.landmarks_json).read_text()) if args.landmarks_json else []
    reference={'kind':'existing_dx','asset':file_record(args.reference_volume_dx,target.parent)} if args.reference_volume_dx else {
        'kind':'landmark_sphere','radius_A':args.radius,'max_seed_distance_A':min(3.,args.radius),'ambiguity_distance_A':.5}
    data={'schema_version':1,'region_id':args.region_id,'label':args.region_id,
      'review':{'status':'candidate','rationale':'Candidate regional observable; edit source-backed identity, assembly and landmark evidence before anatomical review.','sources':[]},
      'system':{'topology':file_record(args.topology,target.parent),'trajectory':file_record(args.trajectory,target.parent) if args.trajectory else None,
        'reference_frame':0,'protein_selection':args.selection,'obstacle_selection':args.obstacle_selection,'pbc':'check',
        'assembly':'Unreviewed input assembly; no operators applied.','preparation':'Starting-model and preparation provenance unprovided.',
        'membrane':'Membrane orientation and anatomical sides unassigned; profile axis is a coordinate convention.'},
      'alignment':{'residues':None,'rationale':'All matching protein CA atoms (heavy-atom fallback); stable scaffold not independently established.'},
      'geometry':{'spacing_A':args.spacing,'probe_radius_A':args.probe_radius,'region_margin_A':2.,'axis':[0.,0.,1.],
                  'grid_phase':[0.,0.,0.],'profile_padding_A':4.,'max_grid_points':8000000},
      'reference':reference,'landmarks':landmarks}
    validate_definition(data);target.parent.mkdir(parents=True,exist_ok=True);write_json(target,data)
    print(f'candidate region definition={target}')


def _cmd_region_prepare(args):
    from .region_definition import prepare_region
    files=prepare_region(args.definition,args.out_dir,allow_radii_mismatch=args.allow_radii_mismatch)
    print('region review='+files['review_scene_json']+' landmarks='+files['review_landmarks_csv'])


def _cmd_region_trajectory(args):
    from .region_definition import trajectory_arguments,write_json
    root=Path(args.out_dir).resolve()
    if root.exists() and any(root.iterdir()):raise ValueError('Region trajectory requires an empty output directory; choose a new --out-dir')
    from .radii import active_radii,use_radii
    argv,metadata=trajectory_arguments(args.definition,root,radii=active_radii(),
                                       allow_radii_mismatch=args.allow_radii_mismatch)
    radius_set=metadata.pop('radius_set')
    resolved=build_parser().parse_args(argv)
    for key,value in vars(args).items():
        if key not in {'definition','out_dir','command','func'}:setattr(resolved,key,value)
    root.mkdir(parents=True,exist_ok=True)
    scaffold=metadata['definition']['alignment']['residues']
    if scaffold is not None:
        path=root/(metadata['region_id']+'_alignment_residues.json');write_json(path,scaffold)
        resolved.alignment_residues_json=str(path)
    resolved.region_definition_metadata=metadata
    write_json(root/(metadata['region_id']+'_region_provenance.json'),metadata)
    if metadata['radii_check']['status']!='requested_matches_definition':
        print(f"radii={metadata['radii']['name']} ({metadata['radii_check']['status']}, "
              f"table {metadata['radii']['table_sha256'][:12]})")
    # The definition's radius set (or the explicitly allowed override) for every frame and worker.
    with use_radii(radius_set):
        _cmd_cavity_trajectory(resolved)


def _cmd_region_compare(args):
    from .region_comparison import compare_regions
    result=compare_regions(args.left_statistics,args.right_statistics,args.out_dir,
        left_definition=args.left_definition,right_definition=args.right_definition,
        left_input_provenance=args.left_input_provenance,right_input_provenance=args.right_input_provenance,
        allow_radii_mismatch=args.allow_radii_mismatch)
    print(f"paired region frames={result['frame_count']} comparison={Path(args.out_dir)/'comparison.json'}")

# Python API

All names below can be imported directly from the top-level package, for example
`from crevice import pore_profile`. Each entry links to the documentation in its
defining module. The [module index](modules.rst) lists every submodule, including
workflow functions that are not exported at the top level (see
[Submodule-only workflows](#submodule-only-workflows)).

:::{admonition} Early interface
:class: note

Before version 1.0, a minor release may change some signatures. Low-level
defaults that differ from the CLI are listed in
[CLI or Python?](../getting-started/cli-vs-python.md#defaults-that-differ-between-the-cli-and-python).
:::

## Loading structures and trajectories

Read PDB/mmCIF/GRO files, fetch RCSB entries with provenance, and open MD trajectories.

```{autosummary}
:nosignatures:

~crevice.parser.load_structure
~crevice.parser.load_structure_report
~crevice.pdb.resolve_structure
~crevice.pdb.fetch_structure
~crevice.pdb.is_accession
~crevice.pdb.normalize_accession
~crevice.pdb.structure_url
~crevice.pdb.cached_structure_path
~crevice.pdb.StructureSource
~crevice.pdb.DownloadRecord
~crevice.pdb.verify_cached_structure
~crevice.pdb.read_provenance
~crevice.pdb.provenance_path
~crevice.trajectory.load_trajectory
~crevice.trajectory.load_trajectory_report
~crevice.trajectory.inspect_trajectory
~crevice.trajectory.load_static_frames
~crevice.trajectory.trajectory_from_frames
~crevice.trajectory.Trajectory
~crevice.trajectory.align_frame_report
```

## Channel profiles

Resolve a two-ended channel and measure radii along it.

```{autosummary}
:nosignatures:

~crevice.channels.pore_profile
~crevice.channels.profile_along_path
~crevice.channels.trace_centerline
~crevice.channels.detect_bottlenecks
~crevice.channels.channel_volume
~crevice.residue_evidence.profile_constriction_evidence
```

## Casts and voids

Fill cavities, pockets and channels with atom-clear samples.

```{autosummary}
:nosignatures:

~crevice.rolling.rolling_probe_cast
~crevice.voids.void_cast
~crevice.voids.validate_void_cast
~crevice.voids.component_profile
~crevice.cavities.classify_voids
~crevice.cavities.cavity_volume
~crevice.cavities.cavity_surface_points
```

## Cavities and tunnels

Grid searches for enclosed voids and widest exit paths.

```{autosummary}
:nosignatures:

~crevice.cavities.detect_cavities
~crevice.tunnels.find_tunnels
~crevice.tunnels.rank_tunnels
~crevice.tunnels.cluster_tunnels
~crevice.tunnels.tunnel_profile
```

## Residues and networks

Residue roles, boundary attribution and contact graphs.

```{autosummary}
:nosignatures:

~crevice.analysis.annotate_residues
~crevice.residues.lining_residues
~crevice.residues.bottleneck_residues
~crevice.analysis.nearest_residues
~crevice.residues.residue_contact_table
~crevice.residue_evidence.boundary_residue_evidence
~crevice.residue_evidence.read_binary_dx
~crevice.residue_evidence.restore_viewer_identities
~crevice.residue_evidence.select_stick_residues
~crevice.networks.build_residue_network
~crevice.networks.build_cavity_network
~crevice.networks.network_metrics
~crevice.networks.compare_networks
~crevice.connectivity.persistent_network
~crevice.connectivity.lining_connectivity
```

## Trajectory analysis and statistics

Frame-by-frame analysis, distributions and approximate mean uncertainty.

```{autosummary}
:nosignatures:

~crevice.trajectory.analyze_trajectory
~crevice.trajectory.profile_timeseries
~crevice.trajectory.cavity_timeseries
~crevice.trajectory.network_timeseries
~crevice.ensemble.profile_distribution
~crevice.uncertainty.block_mean_confidence
~crevice.residue_dynamics.ResidueDynamics
```

## Features and exploratory utilities

Descriptive features, clustering and dimensionality reduction. No trained models are included.

```{autosummary}
:nosignatures:

~crevice.ml.extract_features
~crevice.ml.profile_features
~crevice.ml.cluster_profiles
~crevice.ml.cluster_conformations
~crevice.ml.reduce_dimensions
~crevice.ml.representative_frames
~crevice.ml.representative_profile
~crevice.ml.residue_importance
~crevice.ml.transition_summary
```

## Atomic radii

Van der Waals radius sets used for every clearance (see
[Atomic radii](../methods/atomic-radii.md)).

```{autosummary}
:nosignatures:

~crevice.radii.RadiusSet
~crevice.radii.radii_preset
~crevice.radii.custom_radii_set
~crevice.radii.load_radius_file
~crevice.radii.read_hole_radius_file
~crevice.radii.resolve_radii
~crevice.radii.use_radii
~crevice.radii.active_radii
~crevice.radii.radii_provenance
~crevice.radii.atom_vdw_radius
~crevice.radii.infer_element
~crevice.radii.ion_name_source
```

## Writers and bundles

Write results to disk. Analysis functions never write files themselves.

```{autosummary}
:nosignatures:

~crevice.io.write_profile_json
~crevice.io.write_profile_csv
~crevice.io.write_profile_pdb
~crevice.io.write_residue_contacts_csv
~crevice.io.write_residue_contacts_json
~crevice.io.write_cavities_json
~crevice.io.write_cavities_pdb
~crevice.io.write_tunnels_json
~crevice.io.write_tunnels_pdb
~crevice.io.write_network_json
~crevice.io.write_trajectory_json
~crevice.io.write_void_cast_json
~crevice.io.write_void_cast_pdb
~crevice.io.write_void_cast_csv
~crevice.io.write_cavities_csv
~crevice.io.write_tunnels_csv
~crevice.io.write_tunnel_points_csv
~crevice.io.write_network_csv
~crevice.io.write_connectivity_csv
~crevice.io.write_features_csv
~crevice.io.write_trajectory_csv
~crevice.io.summary_statistics_rows
~crevice.io.write_summary_statistics_csv
~crevice.volume_export.write_void_cast_dx
~crevice.ensemble.write_profile_distribution_csv
~crevice.residue_dynamics.write_dynamics_summary
~crevice.residue_evidence.write_residue_evidence_bundle
~crevice.volume_export.write_volume_viewer_bundle
~crevice.figures.write_static_publication_bundle
~crevice.figures.write_publication_pymol_script
~crevice.io.export_pymol_script
~crevice.schemas.result_envelope
~crevice.schemas.SCHEMA_VERSION
```

## Plotting

Matplotlib figures. Each function writes one image.

```{autosummary}
:nosignatures:

~crevice.figures.plot_profile_radius
~crevice.figures.plot_profile_radius_with_residues
~crevice.figures.plot_residue_contacts
~crevice.figures.plot_cavity_summary
~crevice.figures.plot_network_summary
~crevice.figures.plot_network_chord
~crevice.figures.plot_trajectory_profiles
~crevice.connectivity.plot_lining_connectivity
```

## Data models

Frozen dataclasses returned by the analyses.

```{autosummary}
:nosignatures:

~crevice.models.Atom
~crevice.models.StructureFrame
~crevice.models.ChannelPoint
~crevice.models.PoreProfile
~crevice.models.Cavity
~crevice.models.Tunnel
~crevice.models.VoidComponent
~crevice.models.VoidCast
~crevice.models.ResidueContact
~crevice.models.NetworkNode
~crevice.models.NetworkEdge
~crevice.models.ResidueInteractionNetwork
~crevice.models.FrameAnalysis
~crevice.models.TrajectoryAnalysis
```

## Benchmark registry

Configured (uncurated) benchmark systems and suite runner.

```{autosummary}
:nosignatures:

~crevice.benchmarks.BENCHMARK_SYSTEMS
~crevice.benchmarks.DEFAULT_BENCHMARK_IDS
~crevice.benchmarks.BenchmarkSystem
~crevice.benchmarks.benchmark_systems
~crevice.benchmarks.get_benchmark_system
~crevice.benchmarks.benchmark_path
~crevice.benchmarks.fetch_benchmark_structure
~crevice.benchmarks.fetch_benchmark_record
~crevice.benchmarks.fetch_all_benchmark_structures
~crevice.static_suite.run_static_benchmark_suite
```

(submodule-only-workflows)=
## Submodule-only workflows

These functions drive the `cast`, `cavity-trajectory`, `hydration` and `region-*`
commands but are not re-exported from `crevice`. Import them from their modules.

```{autosummary}
:nosignatures:

~crevice.figures.write_rolling_probe_bundle
~crevice.rolling.select_cast_atoms
~crevice.cavity_trajectory.CavityReference
~crevice.cavity_trajectory.analyze_cavity_trajectory
~crevice.cavity_trajectory.summarize_cavity_trajectory
~crevice.cavity_trajectory.write_cavity_trajectory_bundle
~crevice.hydration.static_hydration
~crevice.hydration_trajectory.analyze_hydration_trajectory
~crevice.hydration_export.write_hydration_bundle
~crevice.hydration_workflow.static_hydration_bundle
~crevice.region_definition.load_definition
~crevice.region_definition.prepare_region
~crevice.region_definition.trajectory_arguments
~crevice.region_comparison.compare_regions
```

# API reference

CREVICE's Python interface mirrors the command line: each command converts its
options and calls the functions below, so anything you can do with `crevice`
you can also do from a script or a notebook.

```python
import crevice

frame = crevice.load_structure("1GRM.cif")
profile = crevice.pore_profile(frame)
print(profile.min_radius)
```

Names shown without a module path can be imported straight from `crevice`
(for example `from crevice import pore_profile`). Names shown with their module
path, such as `crevice.hydration.static_hydration`, are imported from that
module. Every entry links to its full documentation, and the
[module index](modules.rst) lists every submodule.

(api-mdanalysis)=
## MDAnalysis objects as input

Every function that reads a structure or a trajectory from a path also accepts
an MDAnalysis `Universe` or `AtomGroup` in place of the topology path (see
{mod}`crevice.mdinput`). Static loaders read the current frame; trajectory
loaders iterate the trajectory with the usual `start`/`stop`/`stride` and put
the Universe back on the frame it was on. An `AtomGroup` is the selection
itself (its atoms replace `selection=`); waters and obstacles are still taken
from its whole Universe. The results are identical to reading the same files
by path.

```python
import MDAnalysis as mda
import crevice

u = mda.Universe("system.gro", "production.xtc")
u.trajectory[100]
frame = crevice.load_structure(u)                      # frame 100, md_selection="protein"
traj = crevice.load_trajectory(u, stride=10)           # every 10th frame
chain_a = crevice.load_trajectory(u.select_atoms("protein and segid PROA"), stride=10)
result = crevice.analyze_trajectory(traj)

from crevice.hydration_trajectory import analyze_hydration_trajectory
hydration = analyze_hydration_trajectory(u, residue_ids=["A:SER83"])  # residue IDs as CREVICE labels them
```

The same applies to {func}`~crevice.trajectory.inspect_trajectory`,
{func}`~crevice.fingerprint.fingerprint_trajectory`,
{func}`~crevice.tracking.track_cavities`,
{func}`~crevice.residue_interactions.trajectory_residue_interactions` and the
obstacle and water sources of the cavity-trajectory workflow.

## Figures

Every `plot_*` function returns the `matplotlib.figure.Figure`. Pass a path to
also save it (the saved file is unchanged from earlier releases), or omit it to
get the figure only:

```python
fig = crevice.plot_profile_radius(profile)            # nothing written
fig.axes[0].set_xlim(-10, 10)
fig.savefig("profile.svg")
crevice.plot_profile_radius(profile, "profile.png")   # saved and returned
```

The returned figure is closed in pyplot, so loops do not accumulate open
figures; in a notebook, evaluate it to display it. `fig.crevice_description`
holds the description that saved PNG/PDF/SVG files carry in their metadata.

## Uncertainty and sensitivity

{func}`crevice.uncertainty.sampling_uncertainty` gives the effective sample
size, the blocking standard error and a moving-block bootstrap interval of the
mean of one per-frame series, with the block sizes used;
{func}`crevice.uncertainty.trajectory_frame_statistics` applies it to every
per-frame quantity of a trajectory analysis. See
[Uncertainty, sensitivity and grid resolution](../methods/uncertainty.md).

```python
from crevice.uncertainty import sampling_uncertainty, trajectory_frame_statistics

stats = trajectory_frame_statistics(result)
s = stats["quantities"]["min_radius_A"]["sampling_uncertainty"]
print(s["effective_sample_size"], s["block_bootstrap_lower"], s["block_bootstrap_upper"],
      s["block_bootstrap_block_length_frames"])
```

Signatures may still change between minor releases before version 1.0. Where
a Python default differs from the command line, it is listed in
[CLI or Python?](../getting-started/cli-vs-python.md#defaults-that-differ-between-the-cli-and-python).

(api-structures)=
## Structures and inputs

Read PDB, mmCIF and GRO files, fetch entries from the RCSB PDB with provenance, and the atom and frame records everything else works on.

```{eval-rst}
.. autosummary::
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
   ~crevice.models.Atom
   ~crevice.models.StructureFrame
```

(api-profiles)=
## Profiles and casts

Measure a channel's radius along its length, and fill pockets, cavities and channels with atom-clear probe spheres.

```{eval-rst}
.. autosummary::
   :nosignatures:

   ~crevice.channels.pore_profile
   ~crevice.channels.profile_along_path
   ~crevice.channels.trace_centerline
   ~crevice.channels.detect_bottlenecks
   ~crevice.channels.channel_volume
   ~crevice.rolling.rolling_probe_cast
   ~crevice.voids.void_cast
   ~crevice.voids.validate_void_cast
   ~crevice.voids.component_profile
   ~crevice.cavities.classify_voids
   ~crevice.cavities.cavity_volume
   ~crevice.cavities.cavity_surface_points
   crevice.rolling.select_cast_atoms
   ~crevice.models.ChannelPoint
   ~crevice.models.PoreProfile
   ~crevice.models.VoidComponent
   ~crevice.models.VoidCast
   ~crevice.fingerprint.fingerprint_structure
   ~crevice.fingerprint.pore_fingerprint
   ~crevice.fingerprint.fingerprint_trajectory
   ~crevice.fingerprint.conductance_estimate
   crevice.paths.find_paths
   crevice.paths.PathNetwork
   crevice.paths.PathBranch
   crevice.paths.residue_seed
   crevice.paths.region_seed
   crevice.paths.FreeSpace
   crevice.paths.write_paths_outputs
   crevice.paths.plot_paths
   crevice.paths.write_path_scenes
   crevice.paths.paths_view
   ~crevice.fingerprint.write_fingerprint_bundle
   ~crevice.fingerprint.plot_pore_fingerprint
   ~crevice.fingerprint.PoreFingerprint
```

(api-cavities)=
## Cavities and tunnels

Grid searches for buried cavities and for the widest routes out of the protein.

```{eval-rst}
.. autosummary::
   :nosignatures:

   ~crevice.cavities.detect_cavities
   ~crevice.tunnels.find_tunnels
   ~crevice.tunnels.rank_tunnels
   ~crevice.tunnels.cluster_tunnels
   ~crevice.tunnels.tunnel_profile
   ~crevice.models.Cavity
   ~crevice.models.Tunnel
```

(api-residues)=
## Residues and networks

Which residues line a channel or form a cavity boundary, and how residues contact each other.

```{eval-rst}
.. autosummary::
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
   ~crevice.residue_evidence.profile_constriction_evidence
   ~crevice.models.ResidueContact
   ~crevice.models.NetworkNode
   ~crevice.models.NetworkEdge
   ~crevice.models.ResidueInteractionNetwork
```

(api-hydration)=
## Hydration

Water contacts, solvent-accessible area and water density, for one structure or a trajectory.

```{eval-rst}
.. autosummary::
   :nosignatures:

   crevice.hydration.static_hydration
   crevice.hydration_workflow.static_hydration_bundle
   crevice.hydration_trajectory.analyze_hydration_trajectory
   crevice.hydration_export.write_hydration_bundle
```

(api-trajectories)=
## Trajectories and regions

Open MD trajectories, analyse every frame, summarise with approximate uncertainty, and work with named regions.

```{eval-rst}
.. autosummary::
   :nosignatures:

   ~crevice.trajectory.load_trajectory
   ~crevice.trajectory.load_trajectory_report
   ~crevice.trajectory.inspect_trajectory
   ~crevice.trajectory.load_static_frames
   ~crevice.trajectory.trajectory_from_frames
   ~crevice.trajectory.Trajectory
   ~crevice.trajectory.OriginClearanceWarning
   ~crevice.trajectory.align_frame_report
   ~crevice.trajectory.analyze_trajectory
   ~crevice.trajectory.profile_timeseries
   ~crevice.trajectory.cavity_timeseries
   ~crevice.trajectory.network_timeseries
   ~crevice.ensemble.profile_distribution
   ~crevice.uncertainty.block_mean_confidence
   crevice.uncertainty.sampling_uncertainty
   crevice.uncertainty.trajectory_frame_statistics
   crevice.uncertainty.integrated_autocorrelation_time
   crevice.uncertainty.blocking_standard_error
   crevice.uncertainty.moving_block_bootstrap
   crevice.sensitivity.run_sensitivity
   crevice.grid_guard.grid_resolution_check
   ~crevice.residue_dynamics.ResidueDynamics
   crevice.cavity_trajectory.CavityReference
   crevice.cavity_trajectory.analyze_cavity_trajectory
   crevice.cavity_trajectory.summarize_cavity_trajectory
   crevice.cavity_trajectory.write_cavity_trajectory_bundle
   crevice.tracking.track_cavities
   crevice.tracking.CavityTracker
   crevice.tracking.TrackedComponent
   crevice.tracking.measure_frame
   crevice.tracking.voxel_keys
   crevice.tracking.write_tracking_outputs
   crevice.tracking.plot_tracking
   crevice.region_definition.load_definition
   crevice.region_definition.prepare_region
   crevice.region_definition.trajectory_arguments
   crevice.region_comparison.compare_regions
   ~crevice.models.FrameAnalysis
   ~crevice.models.TrajectoryAnalysis
```

(api-radii)=
## Atomic radii

The radius sets behind every clearance (see [Atomic radii](../methods/atomic-radii.md)).

```{eval-rst}
.. autosummary::
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

(api-outputs)=
## Outputs and viewers

Write results, figures and PyMOL/VMD/ChimeraX scenes. Analysis functions never write files themselves, so these are the only functions that touch the disk.

```{eval-rst}
.. autosummary::
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
   crevice.figures.write_rolling_probe_bundle
   ~crevice.figures.plot_profile_radius
   ~crevice.figures.plot_profile_radius_with_residues
   ~crevice.figures.plot_residue_contacts
   ~crevice.figures.plot_cavity_summary
   ~crevice.figures.plot_network_summary
   ~crevice.figures.plot_network_chord
   ~crevice.figures.plot_trajectory_profiles
   ~crevice.connectivity.plot_lining_connectivity
```

(api-view)=
## Interactive views

Show results in 3D in Jupyter (with the `view` extra) or write them as one self-contained HTML file. See [Explore results interactively](../examples/explore-interactively.md).

```{eval-rst}
.. autosummary::
   :nosignatures:

   ~crevice.viewer.view
   ~crevice.viewer.CreviceView
   ~crevice.viewer.Layer
   ~crevice.viewer.field_mesh
   ~crevice.viewer.cast_mesh
   ~crevice.viewer.read_dx_field
```

(api-more)=
## Features and benchmarks

Descriptive features and clustering for exploratory work (no trained models are included), and the registry of benchmark structures, which are not curated references.

```{eval-rst}
.. autosummary::
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

```{toctree}
:hidden:

modules
```

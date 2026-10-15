# CLI or Python?

You can use CREVICE from the command line or from Python; both run the same
calculations. This page helps you choose, shows which function sits behind
each command, and lists the few defaults that differ.

| | Command line (`crevice ...`) | Python (`import crevice`) |
|---|---|---|
| Best for | standard runs, batch scripts, reproducible bundles | custom pipelines, notebooks, new analyses |
| Input | a file path, an RCSB accession or an AlphaFold DB identifier | a {class}`~crevice.models.StructureFrame` or {class}`~crevice.trajectory.Trajectory` you load first |
| Output | files (JSON, CSV, PNG, DX, viewer scripts) plus a manifest | result objects (dataclasses and dictionaries); no files unless you call a writer |
| Provenance | recorded automatically in manifests and input reports | your responsibility; the result objects carry method metadata |
| Selections | flags such as `--md-selection`, `--selection`, `--exclude-hetero` | keyword arguments such as `md_selection=`, `include_hetero=` |
| Atomic radii | `--radii {bondi,hole,charmm_like,PATH}` for the whole command | `radii=` on each analysis, or {func}`crevice.radii.use_radii` for a block |
| Figure text | clean by default; `--annotate` | clean by default; `annotate=True` or {func}`crevice.presentation.figure_annotations` |
| Display guides | off by default; `--mouth-guides`, `--exit-centre-lines` | off by default; {func}`crevice.presentation.display_guides` or `mouth_guides=`/`exit_centre_lines=` on the scene writer |

## Mapping commands to functions

The commands are thin wrappers that load the input, call the analysis, then
write files and figures. The main public functions are:

| Command | Main Python entry point | Writers / plots |
|---|---|---|
| `profile` | {func}`crevice.pore_profile` | {func}`crevice.write_profile_json`, {func}`crevice.write_profile_csv`, {func}`crevice.plot_profile_radius` |
| `residues` | {func}`crevice.annotate_residues` | {func}`crevice.write_residue_contacts_csv`, {func}`crevice.plot_residue_contacts` |
| `cavities` | {func}`crevice.detect_cavities` | {func}`crevice.write_cavities_json`, {func}`crevice.write_cavities_csv`, {func}`crevice.plot_cavity_summary` |
| `tunnels` | {func}`crevice.find_tunnels` | {func}`crevice.write_tunnels_json`, {func}`crevice.write_tunnels_csv`, {func}`crevice.write_tunnel_points_csv` |
| `network` | {func}`crevice.build_cavity_network`, {func}`crevice.build_residue_network` | {func}`crevice.write_network_json`, {func}`crevice.write_network_csv`, {func}`crevice.plot_network_summary` |
| `cast` | {func}`crevice.rolling_probe_cast` | {func}`crevice.write_volume_viewer_bundle`, {func}`crevice.write_void_cast_dx`, {func}`crevice.write_void_cast_csv` |
| `residue-evidence` | {func}`crevice.boundary_residue_evidence`, {func}`crevice.select_stick_residues` | {func}`crevice.write_residue_evidence_bundle` |
| `trajectory` | {func}`crevice.load_trajectory`, {func}`crevice.analyze_trajectory`, {func}`crevice.profile_distribution` | {func}`crevice.plot_trajectory_profiles`, {func}`crevice.write_trajectory_json`, {func}`crevice.write_trajectory_csv` |
| `cavity-trajectory` | {func}`crevice.cavity_trajectory.analyze_cavity_trajectory` | {func}`crevice.cavity_trajectory.write_cavity_trajectory_bundle` |
| `hydration` | {func}`crevice.hydration.static_hydration`, {func}`crevice.hydration_trajectory.analyze_hydration_trajectory` | {func}`crevice.hydration_export.write_hydration_bundle` |
| `region-*` | {func}`crevice.region_definition.prepare_region`, {func}`crevice.region_comparison.compare_regions` | written by the same functions |
| `publish` | {func}`crevice.static_workflow.compute_static_workflow` (no files) | {func}`crevice.static_workflow.export_static_workflow`; {func}`crevice.write_static_publication_bundle` does both |

The first eight rows can be imported straight from `crevice`. The trajectory
cavity, hydration and region workflows are reached through their submodules,
and some of their helpers take CLI argument objects.

## The same run from both sides

{func}`crevice.config.settings` builds the typed settings of any command
(defaults are read from the command-line parser, so they cannot drift) and
{func}`crevice.config.run` runs them through exactly the code the command line
uses, including the output policy and the manifest's run record:

```python
from crevice.config import settings, run

s = settings("publish", "1GRM", out_dir="results/1GRM", prefix="1GRM",
             samples=101, outputs="default,hydration")
print(s.to_dict()["geometry"]["samples"])   # the same JSON as --print-config
record = run(s)                             # same files as the command line
```

To measure first and decide what to write later, split the `publish` workflow:

```python
from crevice import load_structure
from crevice.static_workflow import compute_static_workflow, export_static_workflow

frame = load_structure("1GRM.cif")
result = compute_static_workflow(frame, samples=101)   # writes nothing
print(result.status, result.profile.min_radius, result.cast.total_volume)
export_static_workflow(result, structure_path="1GRM.cif", output_dir="out",
                       prefix="1GRM", outputs="figures,pymol")
```

The bundle writers ({func}`crevice.write_static_publication_bundle`,
{func}`crevice.figures.write_rolling_probe_bundle`) accept the same
`outputs=` values as `--outputs` (plus `"all"` and `"legacy"`); without it they
write the complete bundle, as in 0.2.

## Defaults that differ between the CLI and Python

Some low-level functions have different defaults from the command line. If you
want results that match the command line, pass these values explicitly:

- Residue-network distance: `crevice network` and `crevice trajectory` use the
  heavy-atom **centre** distance (4.5 Å). {func}`crevice.build_residue_network`
  and {func}`crevice.build_cavity_network` default to the van der Waals
  **surface-gap** metric. Pass `distance_metric="center"` to match
  `crevice network`. The networks built *inside* `publish` and `analyze` still
  use those low-level defaults (surface gap, heteroatoms excluded). The same
  named analysis can therefore give different numbers depending on the command
  that ran it. This split is deliberate for now: unifying it is an open
  interface-design decision, and each network records its metric and
  heteroatom choice in `metadata`.
- Trajectories: `crevice trajectory` adds a 0.95 mean confidence band
  (`--confidence`); {func}`crevice.profile_distribution` adds none unless you
  pass `confidence=`. Alignment is on by default in both (`--align`,
  `align=True`). For a plain list of frames, rather than a
  `TrajectoryAnalysis`, a confidence band needs `assume_aligned=True`.
  `analyze_trajectory` builds networks with the surface-gap default; pass
  `network_kwargs={"distance_metric": "center"}` to match the command.
- Profiles: `--samples`/`samples=` (default 81) is a minimum for connected
  profiles in both interfaces; extra points keep axial steps at or below
  0.75 Å (`metadata["requested_samples"]` and `["samples"]`). The enclosure
  probe is chosen automatically in both (`--enclosure-radius auto`,
  `enclosure_radius=None`; see [`profile`](../tools/profile.md)), including
  `static-suite`/{func}`crevice.static_suite.run_static_benchmark_suite` and
  `benchmark`. For trajectories, `crevice trajectory` and
  {func}`crevice.analyze_trajectory` both choose once on the reference frame
  and use that probe for every frame (`metadata["enclosure_probe"]`). A number
  fixes the probe in every interface. `crevice cast` and
  {func}`crevice.rolling_probe_cast` have no automatic probe: their
  `probe_radius` (default 0.8 Å) defines the cast itself (see
  [`cast`](../tools/cast.md#why-the-probe-radius-is-fixed)).
- Benchmarks and fetching: the Python helpers and the commands all default to
  mmCIF (`fmt="cif"`/`extension="cif"`, `--format cif`).
- Cavity casts: `crevice cast` selects the dominant cavity (`--max-cavities 1`).
  {func}`crevice.rolling_probe_cast` defaults to `selection_mode="all"` with no
  count limit. Pass `selection_mode="dominant", max_components=1,
  enclosure_fraction=0.9` to reproduce the CLI.
- Heteroatoms: `profile`, `residues`, `network`, `analyze` and `publish` keep
  HETATM records unless you pass `--exclude-hetero`. The standalone `cavities`
  and `tunnels` commands drop them unless you pass `--include-hetero`, but
  `analyze --analysis cavities/tunnels` keeps them by default. `crevice network`
  keeps heteroatoms, whereas {func}`crevice.build_residue_network` and
  {func}`crevice.build_cavity_network` default to `include_hetero=False`. Check
  `include_hetero` in each signature, and record the command you ran.

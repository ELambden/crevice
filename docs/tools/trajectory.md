# `trajectory`: profiles and contacts across frames

## Scientific question

How does a through-channel profile vary across MD frames (or a set of static
structures)? How precisely is its mean known? Which residue contacts persist,
and how do residues move together? `trajectory` runs the static analyses frame
by frame after alignment and summarises them.

For one-sided pockets and cavities without two open ends, use
[`cavity-trajectory`](cavity-trajectory.md).

## Method

- **Reading.** Either one trajectory file with `--topology`, or several static
  files, one per frame. `--selection` (MDAnalysis syntax) filters atoms while
  reading. `--start/--stop/--stride` choose frames, and `--max-frames` (default
  500) guards against loading too much. `--inspect` prints frame, atom and
  residue counts to the terminal without loading any frames; it needs no `-o`
  (every analysis run does), and `--report-json` also saves the counts.
- **Alignment and periodicity.** Matching CA atoms are fitted to the first
  frame (turn this off with `--no-align`). A scaffold can be predeclared with
  `--alignment-residue` or `--alignment-residues-json`. `--pbc check|unwrap|none`
  controls the split-bond screen. See
  [Inputs and selections](../getting-started/inputs-and-selections.md#trajectories-frames-alignment-and-periodic-boundaries).
- **Profiles.** All frames share the reference axis and seed. `--samples`
  (default 81) is a minimum for connected profiles: more points are used where
  needed to keep axial steps at or below 0.75 Å, and the profile metadata
  records `requested_samples` and `samples`. With `--enclosure-radius auto`
  (the default) the enclosure probe is chosen **once**, on the first
  (reference) frame, by the rule of
  [`profile`](profile.md#choosing-the-enclosure-probe---enclosure-radius-auto),
  and that value is *pinned*: it is passed explicitly for every frame, so all
  frames share one definition of the wall. The command prints the chosen
  probe; the output JSON records it in `metadata.enclosure_probe` (`mode`,
  `chosen_A`, `reason`, and the full reference-frame `selection` record). A
  number (for example `--enclosure-radius 0.8`) runs that probe for every
  frame and records `mode: explicit`. The choice costs one reference-frame
  profile per ladder probe (up to 21), not one per frame.
- **Per-frame fallback.** A pinned probe can fail on a frame whose own
  resolving window does not contain it (the wall moved). With the automatic
  probe, such a frame is then profiled again with the automatic choice made on
  that frame alone (same axis and origin). Each frame records which probe it
  used: `<stem>_frames.csv` has `enclosure_radius_A` (the frame's effective
  probe) and `enclosure_probe_source` (`pinned`, `fallback`, or `none` when
  neither resolved it); the JSON has `frames[].metadata.enclosure_probe`
  (`source`, `radius_A`, `pinned_A`, the pinned failure and the fallback's
  reason or failure), and `metadata.enclosure_probe.fallback` counts pinned,
  fallback and unresolved frames and lists the fallback probe of each
  fallback frame. The command prints the counts. An explicit probe stays
  strict: a frame it cannot resolve is unresolved unless you add
  `--probe-fallback`; `--no-probe-fallback` disables the fallback for the
  automatic probe. A fallback frame's wall is defined by a different probe.
  Inside a channel's resolving window the measured radii barely depend on the
  probe, but its mouths, and so its axial extent, can move; distributions pool
  fallback frames with the others and record how many there were
  (`enclosure_probe_fallback_frame_count`, `enclosure_radii_A`). Each fallback
  costs up to 21 extra profiles of that frame.
- **Unresolved frames.** A frame that does
  not resolve is recorded as missing (`--profile-errors record`). If the
  reference frame does not resolve (including when no probe gives a stable
  channel), later profiles are marked not attempted, `enclosure_probe` records
  `chosen_A: null` with the reason, and other analyses continue.
  `cavity-trajectory` and `region-trajectory` do not build channel profiles and
  have no enclosure probe; their `--probe-radius` is the rolling-probe radius
  of the cast (see [`cast`](cast.md#why-the-probe-radius-is-fixed)).
- **Distributions.** Radii are interpolated onto a shared axial grid without
  extrapolation. A light band shows frame quantiles (`--quantiles`, default
  10th–90th percentile). A dark band shows an approximate confidence interval
  for the **mean**, from a jointly resampled contiguous-batch bootstrap
  (`--confidence`, `--block-length`, `--bootstrap-replicates`). The interval is
  withheld when fewer than eight effective batches exist or coverage is
  incomplete.
- **Residue dynamics** (`--network-json`). Physical contact occupancy (heavy-atom
  centres, `--contact-cutoff` 4.5 Å, persistent threshold `--occupancy` 0.75),
  pooled across chemical labels. Also aligned residue-centroid RMSF and a
  dynamic cross-correlation matrix (DCCM) with split-half diagnostics. With
  `--lining-evidence`, lining membership comes from a reference
  [`residue-evidence`](residue-evidence.md) JSON of the **same** structure.

## Assumptions

- Sampling is stationary and representative, and the batch length captures
  the time correlation. Neither is established automatically.
- Molecules are whole before fitting. CREVICE does not repair periodic images.
- The default all-CA fit is a convenience. It does not claim that those atoms
  form a rigid biological core.

`--radii` selects the atomic radius set for every geometric measurement (default: the standard table); see [Atomic radii](../methods/atomic-radii.md).

## Outputs

- `-o trajectory.json`: per-frame results, reader provenance, fit diagnostics
  and profile status.
- `<stem>_frames.csv` (beside `-o`): one row per frame: `frame_index`,
  `source_index`, `time_ps`, `profile_status`, `profile_reason`,
  `min_radius_A`, `mean_radius_A`, `max_radius_A`,
  `bottleneck_axis_position_A`, `bottleneck_residue`, `profile_point_count`,
  `enclosure_radius_A`, `enclosure_probe_source`, counts, alignment
  RMSD (`_A`) and `feature_*` columns.
- `<stem>_profiles.csv`: every profile sample of every frame (`frame_index`,
  `time_ps` and the [`profile`](profile.md#outputs) CSV columns).
- `--png`: radius distribution (`--view distribution`: light-teal frame
  quantile band, darker teal approximate mean confidence band, teal mean, with
  a legend, and frame coverage below) or time series (`--view timeseries`:
  red minimum, blue mean, grey maximum radius per frame, with a legend). The
  status title (resolved/total frames) and the "mean confidence unavailable"
  note are drawn only with `--annotate`; both are in the PNG metadata.
  `--distribution-csv`: radius statistics and frame counts at each coordinate
  (`coordinate_A`, `mean_radius_A`, `median_radius_A`, quantile and interval
  columns, all in Å even when `--profile-quantity diameter` doubles the
  plotted values); the interval method record is written beside it as
  `<stem>.json`.
- `--network-json`: compact residue-dynamics report, with companion
  `<stem>_contacts.csv`, `<stem>_residues.csv` (`centroid_rmsf_A`, lining
  partners) and partner plot (blue-grey occupancy bars; magenta dots and bars for the
  whole-trajectory and half-trajectory motion correlations).
- `--report-json`: reader provenance.

## Python equivalent

```python
from crevice import analyze_trajectory, load_trajectory, profile_distribution

traj = load_trajectory("system.gro", "run.xtc", selection="protein",
                       stride=10, max_frames=500)
analysis = analyze_trajectory(traj, analyses=("profile",), align=True)
stats = profile_distribution(analysis, confidence=0.95, assume_aligned=True)
```

Python and command-line defaults differ in three places; pass them explicitly
to reproduce a command:

- **Confidence band.** `crevice trajectory` requests a 0.95 mean band
  (`--confidence`, 0 disables it); {func}`crevice.profile_distribution`
  computes none unless `confidence=` is given.
- **Alignment.** Both fit CA atoms to the first frame by default
  (`--align` / `align=True`). {func}`crevice.profile_distribution` trusts the
  `aligned_to_first_frame` record of a `TrajectoryAnalysis`; for a plain list
  of frames a confidence band needs `assume_aligned=True`.
- **Residue contacts.** `crevice trajectory` uses heavy-atom centre distances
  (`--contact-metric center`); `analyze_trajectory` passes `network_kwargs` to
  {func}`crevice.build_cavity_network`, whose default is the surface gap with
  heteroatoms excluded. Pass `network_kwargs={"distance_metric": "center"}` to
  match the command.

## Testing and validation status

- **Software / synthetic:** coverage checks of the correlated-data estimator,
  and fit, PBC-screen and identity contracts.
- **Per-frame fallback (software and one real ensemble):** synthetic tests check
  that a frame the pinned probe cannot resolve is rescued and recorded, that an
  explicit probe stays strict unless `--probe-fallback` is given, and that the
  distribution counts fallback frames. On the five NMR models of 1GRM read as
  a five-frame trajectory, the automatic probe (0.75 Å, chosen on model 1)
  resolves models 1, 3 and 4; the fallback rescues model 5 (0.80 Å, chosen on
  that model); model 2 resolves at none of the probes tried with that model's
  axis and origin. This shows the bookkeeping on real coordinates; it is not a
  validation of any model.
- **Real-input use:** an all-atom membrane-transporter MD trajectory is read and
  aligned, but its reference does **not** give a monotone two-ended channel, so
  its profiles are recorded as unresolved, not as zero. All 64 rigidly
  transformed copies of 1GRM coordinates resolve; that is a positive control on
  real coordinates, not an MD test.
- **Biological / functional:** not established. Contact networks and DCCM are
  associations, not causal pathways.

Method note: [Trajectory profiles and residue evidence](../TRAJECTORY_RESIDUE_METHODS.md).

## Known limitations

- Short or strongly correlated trajectories often give no mean interval, or
  intervals whose nominal coverage is not guaranteed.
- DCCM has no uncertainty estimate.
- The confidence bands leave out force-field error, unsampled states, grid and
  probe error, and wrong channel assignment.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: trajectory
```

# Channel profiles over a trajectory

`crevice trajectory` follows a channel through a molecular-dynamics run, or
through any set of structures such as an NMR ensemble. It aligns the frames,
measures the [profile](profile.md) in each one, and tells you how the radius
varies, how well the mean is known, which residue contacts persist and how
residues move together.

For pockets and cavities that don't run right through the protein, use
[cavities over a trajectory](cavity-trajectory.md).

## Quick start

```bash
crevice trajectory run.xtc --topology system.gro --inspect          # count frames and atoms first
crevice trajectory run.xtc --topology system.gro --stride 10 \
    -o traj.json --png traj_radius.png --distribution-csv traj_distribution.csv
```

```python
from crevice import analyze_trajectory, load_trajectory, profile_distribution

traj = load_trajectory("system.gro", "run.xtc", selection="protein",
                       stride=10, max_frames=500)
analysis = analyze_trajectory(traj, analyses=("profile",), align=True)
stats = profile_distribution(analysis, confidence=0.95, assume_aligned=True)
```

```{figure} ../examples/images/1grm_nmr_profiles.png
:alt: Radius distribution of the gramicidin A channel across its NMR models.
:width: 520px

The five NMR models of 1GRM read as a short trajectory: the frame range, the
interval for the mean and the mean radius along the channel.
```

## How it works

- **Reading frames.** Give one trajectory with `--topology`, or several
  structure files, one per frame. `--selection` (MDAnalysis syntax) chooses the
  atoms, `--start`, `--stop` and `--stride` choose the frames, and
  `--max-frames` (500) stops you loading more than you meant to. `--inspect`
  prints the frame, atom and residue counts without loading anything.
- **Alignment.** Matching CA atoms are fitted to the first frame (turn this
  off with `--no-align`, or choose the fitting residues with
  `--alignment-residue` or `--alignment-residues-json`). `--pbc` controls the
  check for molecules split across the periodic box; see
  [Inputs and selections](../getting-started/inputs-and-selections.md#trajectories-frames-alignment-and-periodic-boundaries).
- **Profiles.** Every frame uses the reference frame's axis and seed. The
  enclosure probe is chosen once, on the reference frame, and then used for
  every frame, so all frames share one definition of the wall
  (`metadata.enclosure_probe`). A number such as `--enclosure-radius 0.8`
  fixes it instead.
- **When a frame doesn't fit.** If the wall moves so that the shared probe
  can't resolve a frame, that frame is profiled again with a probe chosen on
  the frame itself. Each frame records which probe it used
  (`enclosure_radius_A` and `enclosure_probe_source` in the frames table:
  `pinned`, `fallback` or `none`), and the command prints the counts. A
  fixed probe stays strict unless you add `--probe-fallback`;
  `--no-probe-fallback` turns the rescue off. Frames that still don't resolve
  are recorded as missing, never as zero. If the reference frame itself doesn't
  resolve, the later profiles are skipped and the other analyses carry on.
- **Distributions.** Radii are put on a shared axial grid. A light band shows
  the spread across frames (`--quantiles`, 10th to 90th percentile by
  default), and a darker band an approximate confidence interval for the
  **mean**, from a block bootstrap that respects time correlation
  (`--confidence`, `--block-length`, `--bootstrap-replicates`). The interval is
  withheld when there are fewer than eight effective blocks.
- **Residue dynamics** (`--network-json`). How often residue pairs are in
  contact (heavy-atom centres within 4.5 Å; persistent above 75 % of frames),
  how much each residue moves (RMSF) and which residues move together (a
  dynamic cross-correlation matrix with split-half checks). With
  `--lining-evidence`, the lining residues come from a
  [`residue-evidence`](residue-evidence.md) report of the same structure.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--topology` | none | topology for a trajectory file |
| `--selection` | `protein` | MDAnalysis selection of the atoms to read |
| `--start`, `--stop`, `--stride`, `--max-frames` | 0, end, 1, 500 | which frames to analyse |
| `--no-align` | aligned | skip fitting to the first frame |
| `--view` | `distribution` | plot the distribution or a `timeseries` |
| `--confidence` | 0.95 | level of the mean interval; `0` turns it off |
| `--network-json` | none | write the residue-dynamics report |

Channel profile settings, atomic radii and figure text are
[shared options](../reference/cli/common-options.md); everything is listed
under [`crevice trajectory`](../reference/cli/trajectory.md).

## What you get

- **`-o trajectory.json`**: every frame's results, how the frames were read
  and fitted, and the profile status of each.
- **`<stem>_frames.csv`**: one row per frame with its time, profile status,
  minimum, mean and maximum radius, bottleneck position and residue, probe
  used and alignment RMSD.
- **`<stem>_profiles.csv`**: every profile point of every frame.
- **`--png`**: the radius distribution (teal bands and mean, with frame
  coverage below) or, with `--view timeseries`, the minimum (red), mean (blue)
  and maximum (grey) radius per frame.
- **`--distribution-csv`**: the radius statistics at each position along the
  channel, in Å.
- **`--network-json`**: the residue-dynamics report, with contact and
  residue tables and a partner plot.

## Matching the command in Python

Three defaults differ between the command and the Python functions:

- `crevice trajectory` adds a 0.95 mean interval;
  {func}`crevice.profile_distribution` adds none unless you pass
  `confidence=`.
- Both align by default, but for a plain list of frames (rather than a
  `TrajectoryAnalysis`) the interval needs `assume_aligned=True`.
- The command measures residue contacts between heavy-atom centres; pass
  `network_kwargs={"distance_metric": "center"}` to
  {func}`crevice.analyze_trajectory` to match.

:::{admonition} Interpreting results
:class: crevice-interpret

The mean interval assumes your sampling is stationary and representative and
that the blocks are long enough to capture the time correlation; CREVICE can't
check either for you. Short or strongly correlated runs often get no interval,
and the interval never includes force-field error, unsampled states, grid
error or a wrong channel choice. Molecules must be whole before fitting, and
the default CA fit is a convenience, not a claim about a rigid core. Frames
rescued with their own probe pool with the rest, but their mouths can sit
slightly differently. Contact networks and correlated motion are associations,
not pathways. The detailed method is in
[Trajectory profiles and residue evidence](../TRAJECTORY_RESIDUE_METHODS.md).
:::

The [NMR ensemble tutorial](../examples/trajectory-nmr-ensemble.md) runs all
of this on a small public example.

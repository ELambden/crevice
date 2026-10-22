# Uncertainty, sensitivity and grid resolution

Consecutive MD frames are correlated, so a trajectory of N frames holds fewer
than N independent samples. CREVICE reports three complementary,
correlation-aware estimates for the mean of every per-frame summary it writes
for a trajectory, and separately shows how its headline numbers depend on
measurement definitions (`--sensitivity`) and on the grid spacing (the
grid-resolution warning). All of these are statistical or numerical
descriptions of the supplied trajectory; none is a force-field, sampling-state
or biological uncertainty.

## Where the estimates appear

Since 0.4 every statistics record written by `describe_series` carries a
`sampling_uncertainty` block, and the flat summary tables carry the same values
as extra columns (earlier columns and values are unchanged):

| Output | Quantities |
|---|---|
| `PREFIX_cavity_statistics.json` / `_cavity_statistics_summary.csv` (`trajectory --cavities`, `trajectory --region`) | cavity or regional volume, residue boundary areas and fractions, lining occupancy, contact occupancy |
| `*_hydration.json` / `_hydration_summary.csv`, `*_water_membership.json`, interaction summaries | water counts, occupancies, interaction frequencies |
| `comparison.json` (`region-compare`) | paired differences |
| `OUTPUT_frame_statistics.json` / `.csv` (`trajectory`, new) | per-frame minimum, mean and maximum radius, bottleneck position, residue-contact, cavity, tunnel and network counts |
| `--distribution-csv` rows (`trajectory`) | per axial position of the profile distribution |

From Python: {func}`crevice.uncertainty.sampling_uncertainty` (one series),
{func}`crevice.uncertainty.trajectory_frame_statistics` (a
{class}`~crevice.models.TrajectoryAnalysis`), and the lower-level
{func}`~crevice.uncertainty.integrated_autocorrelation_time`,
{func}`~crevice.uncertainty.blocking_standard_error` and
{func}`~crevice.uncertainty.moving_block_bootstrap`.

## Methods

**Effective sample size.** The normalised autocorrelation ρ(t) is computed by
FFT, and the integrated autocorrelation time is
τ_int(M) = 1 + 2 Σ_{t=1..M} ρ(t), with Sokal's automatic window: the
smallest M with M ≥ 5 τ_int(M) (A. D. Sokal, *Monte Carlo methods in
statistical mechanics: foundations and new algorithms*, lecture notes, 1997;
the same estimator as `emcee.autocorr.integrated_time`). The effective sample
size is N / τ_int (`effective_sample_size`; τ_int is reported in frames, and
`autocorrelation_window_converged` is false when no window satisfied the rule).

**Blocking.** Flyvbjerg–Petersen blocking (H. Flyvbjerg and H. G. Petersen,
*J. Chem. Phys.* 91, 461 (1989)) averages adjacent frames pairwise, level by
level; the standard error of the mean is sqrt(s_k / n_k) at the chosen level.
The level is chosen automatically with the lag-one test of M. Jonsson
(*Phys. Rev. E* 98, 043304 (2018)) exactly as in the author's reference
implementation (`github.com/computative/block`, `python/tictoc.py`: the first
level k with M_k < χ²₀.₉₉(k + 1)); we read that code, not the paper's text. Unlike
the reference code, an odd trailing frame is dropped at each level instead of
truncating the series to a power of two. The block size used (2^k frames) is
reported as `blocking_block_size_frames`.

**Moving-block bootstrap.** ceil(N / L) overlapping blocks of L frames are
drawn with replacement (H. R. Künsch, *Ann. Statist.* 17, 1217 (1989)); the
percentile interval of the replicate means is reported as
`block_bootstrap_lower` / `block_bootstrap_upper`, their SD as
`block_bootstrap_standard_error`. The block length is L = ceil(5 τ_int)
unless `--block-length` is given, and is always reported
(`block_bootstrap_block_length_frames`); at least 8 blocks (N / L ≥ 8) are
required, otherwise the status is `insufficient_blocks` and no interval is
given. Replicates and seed follow `--bootstrap-replicates` and
`--bootstrap-seed`, so results are reproducible.

These sit beside, and do not replace, the existing studentized batch-bootstrap
`confidence_interval` (with its Geyer-sequence `statistical_inefficiency`).
Estimates are withheld (with a status) for series with missing frames,
irregular time steps, fewer than 20 frames or no variation.

## How well they work (synthetic check)

On AR(1) series with known autocorrelation (200 replicate series per row;
a synthetic check run during development) the estimates are close to
the exact values when the trajectory is long compared with its correlation
time, and biased low when it is not:

| φ | N | true τ | mean τ_int | true SE | blocking SE | bootstrap SE | 95 % coverage |
|---|---|---|---|---|---|---|---|
| 0.0 | 1000 | 1.0 | 1.06 | 0.032 | 0.031 | 0.032 | 0.955 |
| 0.5 | 1000 | 3.0 | 2.99 | 0.063 | 0.056 | 0.060 | 0.940 |
| 0.9 | 1000 | 19 | 16.9 | 0.316 | 0.231 | 0.275 | 0.853 (177/200 estimated) |
| 0.97 | 1000 | 66 | 43.5 | 1.054 | 0.629 | 0.589 | 0.735 (34/200 estimated) |

Rule of thumb: trust the intervals when N / τ_int is at least about 50; below
that every estimate is optimistic, and a non-converged window or blocking level,
a small effective sample size or an `insufficient_blocks` status says so. The
statistics assume a stationary trajectory; they cannot detect a slow process
the trajectory never sampled. Cost: about 3 ms per 1000-frame series.

## Definition sensitivity (`--sensitivity`)

`--sensitivity` (opt-in) reruns the headline numbers with one setting changed
at a time and writes `PREFIX_sensitivity.csv` (one row per quantity and
variant: value, difference and relative difference from the baseline, status,
seconds) and `PREFIX_sensitivity.json` (adds the spread per quantity and the
total cost). The manifest lists both files.

| Command | Headline quantities | Variants |
|---|---|---|
| `profile` | minimum and mean radius | neighbouring enclosure-probe ladder radii; section grid 0.25 / 0.5 Å |
| `cavities --cast` | cast volume, region count | neighbouring ladder radii around the rolling probe; grid 0.25 / 0.5 Å |
| `report` | minimum radius, cast volume | neighbouring enclosure-probe ladder radii; cast grid 0.25 / 0.5 Å |
| `trajectory --cavities`, `trajectory --region` | mean volume, resolved-frame fraction | margin ±1 Å (reference-region mode); neighbouring probe radii; grid 0.25 / 0.5 Å; on `--sensitivity-frames` (default 10) evenly spaced frames |

The baseline grid value is skipped. Each variant reruns the whole measurement,
so the cost is roughly the number of variants (two to five) times the
baseline run time; on 1GRM the extra time was about 1 s for `profile`, 5 s for
`cavities --cast` and 9 s for `report`. A variant that fails (for example, no channel at
a smaller probe) is recorded with its error and left out of the spread. The
spread is descriptive: the variants are not equally plausible alternatives,
and it is not a confidence interval.

## Grid-resolution warning

Grid-based measurements lose narrow openings and volume as the grid coarsens:
on 1GRM the rolling-probe cast volume was 289.9 / 282.8 / 257.0 / 238.4 Å³ at
0.25 / 0.35 / 0.5 / 0.7 Å, and synthetic walls lost a narrow opening at coarse
grids (CREVICE 0.3 development runs). When the narrowest clearance a
run has to resolve is below two grid spacings, CREVICE prints a warning and
records `grid_resolution` (status, spacing, clearance and its source,
threshold, the recommendation of at most 0.35 Å for narrow channels) in the
manifest; for `paths` it is recorded in the paths JSON settings. Defaults and
results are not changed. Only measured clearances are checked: the profile
minimum radius (`report`) and the narrowest route bottleneck (`paths`). Casts
and per-frame regions (`cavities --cast`, `trajectory --cavities`,
`trajectory --region`) measure no narrowest width, so their check is recorded
as not assessed and no warning is printed; choose a fine grid (at most
0.35 Å) yourself for narrow channels.

### Cast volume and region selection on a grid

The default cast (`cavities --cast`: dominant region, 0.8 Å probe) is defined
by a union of atom-clear spheres centred on grid nodes of a high-clearance
core. Two grid effects follow from that definition.

*Smooth bias.* Sphere centres are a grid sample of the core, so the union
misses a layer of order one spacing at its surface and the volume grows as the
grid is refined, roughly linearly in the spacing. On the 1GRM mmCIF (all
272 heavy atoms; mean of four grid offsets) the volume is 297, 291, 281, 259
and 235 Å³ at 0.2, 0.25, 0.35, 0.5 and 0.7 Å; on 4PYP it is 3377, 3247 and
3116 Å³ at 0.35, 0.5 and 0.7 Å. Grid offsets move a single run by about
±0.5 % at 0.25 Å and up to ±3 % at 0.7 Å. Compare volumes only at the same spacing,
and use 0.25 Å when the volume itself matters.

*Jumps.* Two decisions are made on the grid and can change a volume
discontinuously. (1) With `--max-cavities 1`, one of two similar regions is
kept. (2) Cores whose probe-radius supports touch are filled as one region
spanning both lobes; a neck whose clearance is within about
`spacing * sqrt(3) / 2` of the probe radius can be found at one spacing or
grid offset and missed at another. On T4 lysozyme L99A (1L83) the dominant
region is about 1100 Å³ at 0.2 Å (two lobes joined) and 670–700 Å³ at 0.5 Å
(one lobe), with intermediate spacings switching between the two by grid
offset; even at 0.2 Å a neck to a third chamber lies within the tolerance,
so every 1L83 run is flagged. These are properties of the definition on a grid, not a convergence
failure: at fine enough spacing the decision settles. CREVICE flags both
cases: the cast metadata record `selection_stability` (`near_tie`,
`grid_sensitive_grouping`, a `message`) and the command prints a warning;
nothing is changed. Rerun at a finer spacing (and, for a key number, a few
grid offsets) before reporting a flagged volume.


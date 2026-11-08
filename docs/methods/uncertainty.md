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
| `OUTPUT_frame_statistics.json` / `.csv` (`trajectory`, new) | per-frame minimum, mean and maximum radius, bottleneck position, residue-contact, cavity and network counts |
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
truncating the series to a power of two, and the lag-one autocovariance of each
level is taken about the mean of that level's blocks (for 2^d frames this is the
reference code's result). The block size used (2^k frames) is
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

`benchmarks/nongeom/uncertainty_coverage.py` draws stationary series whose
autocorrelation is known exactly (Gaussian AR(1) with coefficient φ, and a
two-state Markov chain with p = 0.2 as a model of a binary occupancy), 400
replicate series per row, and counts how often each 95 % interval contains the
true mean. Coverage is conditional on an interval being reported; the
fraction of series that got one is in brackets. Selected AR(1) rows (all rows
in `benchmarks/results/nongeom/uncertainty_coverage.csv`):

| φ | N | true τ | N/τ | mean τ_int | true SE | blocking SE | bootstrap SE | batch interval | block bootstrap | mean ± 1.96 blocking SE |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 1000 | 1 | 1000 | 1.06 | 0.032 | 0.032 | 0.032 | 0.98 (1.00) | 0.96 (1.00) | 0.96 |
| 0.5 | 1000 | 3 | 333 | 2.96 | 0.055 | 0.049 | 0.052 | 0.98 (1.00) | 0.94 (1.00) | 0.92 |
| 0.9 | 1000 | 19 | 53 | 17.2 | 0.137 | 0.103 | 0.119 | 0.94 (0.70) | 0.87 (0.88) | 0.85 |
| 0.95 | 1000 | 39 | 26 | 28.8 | 0.196 | 0.133 | 0.141 | 0.88 (0.14) | 0.87 (0.47) | 0.83 |
| 0.98 | 1000 | 99 | 10 | 51.1 | 0.307 | 0.164 | 0.136 | – (0.00) | 0.62 (0.09) | 0.71 |
| 0.9 | 5000 | 19 | 263 | 18.9 | 0.062 | 0.053 | 0.058 | 0.95 (1.00) | 0.91 (1.00) | 0.89 |
| 0.98 | 5000 | 99 | 51 | 88.8 | 0.140 | 0.102 | 0.118 | 0.96 (0.70) | 0.91 (0.85) | 0.88 |

The Markov-chain rows behave the same way at the same N/τ. tau_int and the
Geyer statistical inefficiency agree with `emcee.autocorr.integrated_time`
(to rounding, once emcee's values below 1 are clipped to 1 as here) and with
`pymbar.timeseries.statistical_inefficiency` (means within about 3 % except
where the series is far too short), and the blocking estimate equals a
transcription of Jonsson's reference code on 2^d-frame series
(`benchmarks/results/nongeom/uncertainty_reference.csv`).

Rule of thumb: trust the intervals when N / τ_int is at least about 50; below
that every estimate is optimistic and the batch interval is usually withheld
(a non-converged window or blocking level, a small effective sample size or
an `insufficient_blocks` status says so). Even with long series the
moving-block bootstrap and the blocking standard error cover slightly less
than nominal (about 0.89–0.93 at N/τ of 100–1000); the studentized batch
interval stays at or above about 0.95. The statistics assume a stationary
trajectory; they cannot detect a slow process the trajectory never sampled.
Cost: about 3 ms per 1000-frame series.

## Event counts

Permeation event counts are not frame series, so the estimates above do not
apply to them. Their rates, q0 and pd use a quasi-Poisson interval: the exact
Poisson (Garwood) interval widened by the dispersion index of the counts in
10 contiguous time blocks, never narrower than Garwood. The ion current uses
the block variance of the net charge when it exceeds the Poisson variance.
Coverage on synthetic single-file water and gated currents is in
`benchmarks/nongeom/README.md`; details are in
[Permeation and conductance](../tools/permeation.md#intervals-for-correlated-events).

## Definition sensitivity (`--sensitivity`)

`--sensitivity` (opt-in) reruns the headline numbers with one setting changed
at a time and writes `PREFIX_sensitivity.csv` (one row per quantity and
variant: value, difference and relative difference from the baseline, status,
seconds) and `PREFIX_sensitivity.json` (adds the spread per quantity and the
total cost). The manifest lists both files.

| Command | Headline quantities | Variants |
|---|---|---|
| `profile` | minimum and mean radius | section grid 0.25 / 0.5 Å; with a fixed `--enclosure-radius`, also the neighbouring probe radii |
| `cavities --cast` | cast volume, region count | neighbouring ladder radii around the rolling probe; grid 0.25 / 0.5 Å |
| `report` | minimum radius, cast volume | cast grid 0.25 / 0.5 Å; with a fixed `--enclosure-radius`, also the neighbouring probe radii |
| `trajectory --cavities`, `trajectory --region` | mean volume, resolved-frame fraction | margin ±1 Å (reference-region mode); neighbouring probe radii; grid 0.25 / 0.5 Å; on `--sensitivity-frames` (default 10) evenly spaced frames, with the run's obstacles |

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


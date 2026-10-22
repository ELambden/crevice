# A small public "trajectory": the 1GRM NMR ensemble

The trajectory commands need a topology and a sequence of frames. A solution
NMR entry with several models is a small public stand-in: `trajectory` reads
the five models of [PDB 1GRM](https://www.rcsb.org/structure/1GRM) (legacy PDB
format) as five frames, aligns them, profiles the channel in every frame and
summarises the profiles and residue contacts across frames. It takes under a
minute. NMR models are not a time series: the "time" here is just the model
order (1 ps per model is assumed because the file has no time step), so
correlations along the frame order have no physical meaning. With a real MD run
the same commands take a topology such as `system.gro` and a trajectory such as
`run.xtc`.

## Inspect before analysing

```bash
crevice fetch 1GRM --format pdb --cache-dir .crevice/pdb
crevice trajectory --topology .crevice/pdb/1GRM.pdb --selection all --inspect .crevice/pdb/1GRM.pdb
```

```text
total_atom_count	272
selected_atom_count	272
selected_residue_count	32
frame_count	5
time_step_ps	1.0
box_dimensions	None
```

`--inspect` reports what would be read without analysing anything. `--selection`
is an MDAnalysis selection; `all` keeps every atom. `protein` gives the same
272 atoms here: gramicidin's D-amino acids and modified termini are not in the
MDAnalysis `protein` selection, but CREVICE adds residues peptide-bonded into
the chain to it.

## Profiles in every frame

```bash
crevice trajectory .crevice/pdb/1GRM.pdb --topology .crevice/pdb/1GRM.pdb --selection all \
    -o 1GRM_nmr.json --png 1GRM_nmr_profiles.png \
    --distribution-csv 1GRM_nmr_distribution.csv --network-json 1GRM_nmr_network.json --dpi 150
```

```text
read 5 of 5 frames, 272 atoms per frame
enclosure_radius=0.750 (auto, chosen once on the reference frame and fixed for every frame: smallest probe of channel group 1, which resolved at 5 of 21 probes tried)
enclosure probe per frame: pinned=3 fallback=1 unresolved=1 (fallback on); fallback probes by source frame: {'4': 0.8}
trajectory frames=5 analyses=profile,features,network output=1GRM_nmr.json
```

The water analysis is off by default (since 0.3) and this ensemble has no
explicit water. The typed residue interactions, on by default, need a periodic
box, which an NMR ensemble does not have: the command prints
`crevice: note: interactions not completed (...)`, records the reason under
`metadata.interaction_files` in `1GRM_nmr.json` and completes everything else.
`--skip-interactions` leaves them out.

Every frame is aligned on matching CA atoms to the first frame, whose profile
also fixes the axis and the enclosure probe (0.75 Å, chosen once and pinned
for every frame). A frame where the pinned probe resolves no channel falls back
to the automatic probe choice on that frame. `1GRM_nmr_frames.csv` has one row
per frame:

| Frame | Profile | Probe (Å) | Minimum radius (Å) | Mean radius (Å) |
|---:|---|---:|---:|---:|
| 0 | resolved | 0.75 (pinned) | 1.337 | 1.571 |
| 1 | unresolved | | | |
| 2 | resolved | 0.75 (pinned) | 1.222 | 1.478 |
| 3 | resolved | 0.75 (pinned) | 1.138 | 1.489 |
| 4 | resolved | 0.8 (fallback) | 1.214 | 1.581 |

Model 4 did not resolve with the pinned 0.75 Å probe and was resolved by the
per-frame fallback at 0.8 Å (`enclosure_probe_source` = `fallback`); its radii
were measured with a different probe from the other frames, which is why the
probe is recorded per frame. In model 1 neither the pinned probe nor the
fallback resolved a connected channel with two open ends; the reason is
recorded in the `profile_reason` column. That frame is **missing**, not zero,
and everything below is conditional on the four resolved frames.
`1GRM_nmr_profiles.csv` holds every sample of every resolved frame.

```{figure} images/1grm_nmr_profiles.png
:alt: Mean pore radius and frame range of the 1GRM NMR models along the aligned axis
:width: 85%

`1GRM_nmr_profiles.png` (default `--view distribution`). Upper panel: mean
pore radius (Å, dark teal line) of the resolved frames at each aligned axial
coordinate (Å), and the 10–90% frame range (teal band; `--quantiles`), which
describes fluctuations, not uncertainty. No confidence band for the mean is
drawn: with four frames the block bootstrap has too few independent blocks
(status `insufficient_independent_blocks` in `1GRM_nmr_distribution.json`), so
the interval is withheld rather than reported. Lower panel: fraction of all
frames whose profile covers each coordinate (0.8 = 4 of 5).
```

`1GRM_nmr_distribution.csv` tabulates the same curves per coordinate
(`coordinate_A`, `n_frames`, `coverage_fraction`, `mean_radius_A`,
`median_radius_A`, `lower_radius_A`, `upper_radius_A` and the mean-interval
columns, empty here). `--network-json` writes residue contacts present across
frames (`1GRM_nmr_network_contacts.csv`: occupancy, mean distance in Å, contact
role and the dynamic cross-correlation of residue centroids, which is
descriptive only), a per-residue table and a figure of the persistent contacts
between lining residues and their partners.

## The same in Python

```python
from crevice import load_trajectory, analyze_trajectory, profile_distribution

trajectory = load_trajectory(".crevice/pdb/1GRM.pdb", ".crevice/pdb/1GRM.pdb", selection="all")  # topology, frames
result = analyze_trajectory(trajectory, analyses=("profile",))
print([frame.profile is not None for frame in result.frames])      # [True, False, True, True, True]
distribution = profile_distribution(result, quantiles=(0.1, 0.9))
```

## Real MD trajectories

For an MD run, the cavity- and region-based commands follow the same pattern:
[`trajectory --cavities`](../tools/cavity-trajectory.md) tracks a cavity defined by
a measured reference map (`*_volume.dx` from `cavities --cast` or `report`) through every
frame, and [`trajectory --water`](../tools/hydration.md) and the
[`region-*`](../tools/regions.md) commands add explicit-water analyses. The
command shape (with placeholder file names) is

```bash
crevice trajectory --cavities system.gro run.xtc --reference-volume-dx pocket_volume.dx \
    --geometry-mode reference-region --probe-radius 1.4 --region-margin 2 \
    --out-dir trajectory-results --prefix system
```

These commands need explicit water and periodic-cell information, which the
NMR ensemble does not have, so they are not run on this page.

:::{admonition} Interpreting the results
:class: crevice-interpret

The profile values and the figure come from a CREVICE 0.4.0 run, and the
Python block gave the same pattern of resolved frames. NMR models aren't a time series, so this is a
demonstration of the workflow on public coordinates, not a statement about
gramicidin dynamics.
:::

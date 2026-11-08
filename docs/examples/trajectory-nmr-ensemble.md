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
topology	.crevice/pdb/1GRM.pdb
trajectory	.crevice/pdb/1GRM.pdb
selection	all
total_atom_count	272
selected_atom_count	272
selected_residue_count	32
frame_count	5
time_step_ps	1.0
total_time_ps	4.0
box_dimensions	None
has_elements	True
has_chain_ids	True
segment_ids	['A', 'B']
```

(MDAnalysis also warns on standard error that the reader has no time step and
sets it to 1 ps.)

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
enclosure_radius=per-plane (auto, found again in every frame): resolved=5 unresolved=0
trajectory frames=5 analyses=profile,features,network output=1GRM_nmr.json
```

and on standard error, among the status lines,

```text
crevice: pore radius over 5 frames found in 2.3 s: minimum radius 1.217 A (mean of 5 resolved frame minima, range 1.139-1.337 A; no uncertainty: too few frames)
```

The water analysis is off by default and this ensemble has no explicit water.
The typed residue interactions, on by default, need a periodic box, which an
NMR ensemble does not have: the command prints
`crevice: note: interactions not completed (...)`, records the reason under
`metadata.interaction_files` in `1GRM_nmr.json` and completes everything else.
`--skip-interactions` leaves them out.

Every frame is aligned on matching CA atoms to the first frame, whose profile
fixes the axis and the seed. The enclosure probe is chosen again, slice by
slice, in every frame (`metadata.enclosure_probe` in `1GRM_nmr.json`: mode
`adaptive`, with the per-frame counts). `1GRM_nmr_frames.csv` has one row per
frame; `enclosure_radius_A` is the largest slice probe of that frame and
`enclosure_probe_source` is `per_plane`:

| Frame | Profile | Largest slice probe (Å) | Minimum radius (Å) | Mean radius (Å) |
|---:|---|---:|---:|---:|
| 0 | resolved | 1.4 | 1.337 | 1.582 |
| 1 | resolved | 1.6 | 1.171 | 1.545 |
| 2 | resolved | 1.4 | 1.222 | 1.503 |
| 3 | resolved | 1.4 | 1.139 | 1.506 |
| 4 | resolved | 1.4 | 1.215 | 1.595 |

All five models resolve a channel. A frame that did not would be recorded with
its reason in the `profile_reason` column, and it would be **missing**, not
zero: every summary would then be conditional on the resolved frames.
`1GRM_nmr_profiles.csv` holds every sample of every resolved frame.

```{figure} images/1grm_nmr_profiles.png
:alt: Mean pore radius and frame range of the 1GRM NMR models along the aligned axis
:width: 85%

`1GRM_nmr_profiles.png` (default `--view distribution`). Upper panel: mean
pore radius (Å, dark teal line) of the resolved frames at each aligned axial
coordinate (Å), and the 10–90% frame range (teal band; `--quantiles`), which
describes fluctuations, not uncertainty. No confidence band for the mean is
drawn: with five frames the block bootstrap has too few independent blocks
(status `insufficient_independent_blocks` in `1GRM_nmr_distribution.json`), so
the interval is withheld rather than reported. Lower panel: fraction of all
frames whose profile covers each coordinate (1.0 everywhere here).
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
print([frame.profile is not None for frame in result.frames])      # [True, True, True, True, True]
distribution = profile_distribution(result, quantiles=(0.1, 0.9))
```

## Real MD trajectories

For an MD run, the cavity- and region-based commands follow the same pattern:
[`trajectory --cavities`](../tools/cavity-trajectory.md) tracks a cavity defined by
a measured reference map (`*_volume.dx` from `cavities --cast` or `report`) through every
frame, and [`trajectory --water`](../tools/hydration.md) and the
[named-region workflow](../tools/regions.md) (`crevice region`,
`trajectory --region`) add explicit-water analyses. The
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

The profile values and the figure come from one run of the commands above,
and the Python block gave the same pattern of resolved frames. NMR models aren't a time series, so this is a
demonstration of the workflow on public coordinates, not a statement about
gramicidin dynamics.
:::

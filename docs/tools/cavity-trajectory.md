# `cavity-trajectory`: all-frame cavity statistics

## Scientific question

How do a cavity's volume, cross-sectional widths and boundary residues change
over an aligned MD trajectory? How precisely are their means known? The command
works on one-sided pockets and cavities and does **not** need a through-channel.

## Method

The reference is a measured binary cavity map (`--reference-volume-dx`) in the
**first frame's** coordinates, for example from [`cast`](cast.md) or a prepared
[region](regions.md). Each frame is aligned to the first with matching CA atoms.
Two geometric observables are available:

- `--geometry-mode rolling` (default) recomputes global rolling-probe
  components on a fixed reference lattice in every frame. The component that
  best matches the reference is chosen by an overlap score. A match needs a
  score ≥ `--minimum-overlap` (0.1). The frame is **ambiguous** when a second
  candidate scores ≥ `--ambiguity-ratio` (0.8) times the best. An absent or
  ambiguous match is recorded as missing, never as zero volume.
- `--geometry-mode reference-region` measures a declared local region. It keeps
  every probe-centre component that touches the reference, within
  `--region-margin` (2 Å) of the occupied reference samples, and clips the
  result to that domain. A zero means no free space was measured in that region.
  It does not show that transport is closed.

For each frame, CREVICE records the volume, area-equivalent and local atom-clear
sphere widths along a fixed axis, per-residue boundary area (nearest-VDW
attribution within `--lining-distance` 1.5 Å, excluding crop faces), nonlocal
partners, and all physical contact pairs. Scalar series get mean, SD, empirical
frame range and, when estimable, approximate batch-bootstrap mean intervals.

`--water-membership` also counts, in every frame, the water oxygens inside the
region measured in that frame, paired with the fixed-reference count. It is
computed inside the geometry worker because the per-frame grid is not saved.
See [instantaneous-cavity membership](hydration.md#fixed-reference-counts-and-instantaneous-cavity-membership).

`--obstacle-selection` keeps chosen non-protein atoms (for example lipids) as
walls in `reference-region` mode, with periodic images enumerated. See
[Physical obstacles and numerical sensitivity](../CAVITY_OBSTACLES.md).

## Assumptions

- The reference map is binary, not a smoothed display map, and is in frame-0
  coordinates.
- Geometry uses protein heavy atoms unless `--obstacle-selection` is given.
  Leaving out membrane and ligands can open space that is physically blocked.
- The lattice basis, origin and phase are fixed after alignment. Changing
  `--spacing`, `--axis` or `--grid-phase` resamples the geometry.
- Mean intervals assume representative, stationary sampling.

## Key parameters and units

| Option | Default | Unit |
|---|---:|---|
| `--geometry-mode` | `rolling` | `rolling` or `reference-region` |
| `--spacing` | from the reference map | Å |
| `--probe-radius` | 0.8 | Å |
| `--outer-radius` | 6.0 | Å |
| `--region-margin` | 2.0 | Å (reference-region only) |
| `--axis` | third direction of the map | `x,y,z` |
| `--contact-cutoff` | 4.5 | Å |
| `--max-frames` | 500 | frames (a guard, not a stride) |
| `--workers` | 1 | processes |
| `--water-membership` | off | adds instantaneous-cavity water counts |
| `--radii` | standard table | radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file ([Atomic radii](../methods/atomic-radii.md)) |

## Outputs

CSV tables: `*_cavity_frames.csv` (one row per frame: `time_ps`,
`volume_A3`, `alignment_rmsd_A`, matching scores, axial extent),
`*_cavity_section_frames.csv` (per frame and axis position: section areas
`_A2`, diameters `_A`, component count), `*_cavity_residue_frames.csv` (per
frame and boundary residue: `boundary_area_A2`, fractions, partner flag),
`*_cavity_width_profile.csv`, `*_cavity_section_statistics.csv`,
`*_cavity_free_sphere_statistics.csv`, `*_cavity_residues.csv`,
`*_cavity_contacts.csv`, `*_cavity_partner_statistics.csv`,
`*_cavity_contact_statistics.csv` and `*_cavity_statistics_summary.csv` (every
interval statistic of `*_cavity_statistics.json` on one row: the volume, each
residue's boundary area, fraction and lining occupancy, each contact's
occupancy, and the volume correlations; columns as in
{func}`crevice.io.summary_statistics_rows`). The same arrays are in
`*_cavity_frames.npz`; `*_cavity_statistics.json` and
`*_cavity_profile_statistics.json` are the structured reports read by
`hydration --cavity-results` and `region-compare`, and
`*_cavity_frame_diagnostics.json` keeps the full per-frame alignment and
matching provenance. Figures: shaded width profiles, time–position and residue-area heatmaps, and
partner and motion summaries. In the width profiles a violet band (`#a855f7`, the
non-channel cast colour) is the observed 2.5–97.5% frame range, the violet line
the mean and an orange band the
pointwise mean interval (legends name them); heatmaps use `viridis` (width) and
`cividis` (boundary area) with labelled colour bars, grey = unresolved. Frame
counts, interval status and interpretation reminders are drawn as titles only
with `--annotate`; they are always in the PNG metadata. Standard hydration outputs are appended unless
`--skip-hydration`. With `--water-membership`, the
`PREFIX_water_membership.{json,png}`, `_frames.{csv,npz}`, `_axial.csv`,
`_summary.csv` and `_strata.csv` files
are added; without it no membership file is written and the other outputs are
unchanged. The statistics JSON can later be reused by
[`hydration --cavity-results`](hydration.md) and
[`region-compare`](regions.md).

## Python equivalent

```python
from crevice import load_trajectory
from crevice.cavity_trajectory import (CavityReference, analyze_cavity_trajectory,
                                     write_cavity_trajectory_bundle)

traj = load_trajectory("system.gro", "run.xtc",
                       selection="protein and not name H*", stop=5)
reference = CavityReference.from_dx("reference_volume.dx")
analysis = analyze_cavity_trajectory(traj, reference)
```

Instantaneous-cavity membership is requested with
`analyze_cavity_trajectory(..., waters=WaterSource(...))` and written with
{func}`crevice.water_membership.write_water_membership_bundle`; the CLI flag
assembles the `WaterSource` from the same topology, trajectory and frame range.

This workflow is reached through {mod}`crevice.cavity_trajectory`. Its geometry
settings are passed as a dictionary (`geometry=`), unlike the keyword arguments
used by the static functions. See [CLI or Python?](../getting-started/cli-vs-python.md).

## Testing and validation status

- **Software / synthetic:** matching and ambiguity rules, reference-region
  clipping, crop-face exclusion, periodic obstacle imaging and statistical
  suppression rules are covered by tests. Instantaneous-cavity membership has
  its own synthetic controls (`tests/test_water_membership.py`).
- **Real-input use:** the workflow has been run on an all-atom membrane-protein
  MD trajectory, where the regional volume depended on every geometric setting
  (region margin, probe, grid spacing and phase). Check sensitivity on your own
  system ([how](../CAVITY_OBSTACLES.md#checking-numerical-sensitivity)).
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL, VMD and ChimeraX for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established. A reference pocket is an
  anatomical **candidate** until it is curated, and one simulation is not a
  replica set.

Method note: [Cavity trajectories: measurement and interpretation](../CAVITY_TRAJECTORY_METHODS.md).

## Known limitations

- Split, merge and expansion events can break the match in rolling mode.
- Default volumes are not converged with respect to grid spacing, and domain
  and probe choices change what is measured.
- Independent replicas, anatomical review and finer convergence are still
  needed.
- `--obstacle-selection` works only in `reference-region` mode.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: cavity-trajectory
```

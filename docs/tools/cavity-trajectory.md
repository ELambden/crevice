# Cavities over a trajectory

`crevice cavity-trajectory` follows a cavity through every frame of a
simulation. It tells you how the cavity's volume and widths change, which
residues form its wall and how often, and how precisely the averages are
known. It works on pockets and one-sided cavities, so the cavity doesn't need
to run right through the protein.

## Quick start

You need a reference: a measured cast map of the cavity in the coordinates of
the first frame, for example from [`cast`](cast.md) run on that frame, or from
a prepared [named region](regions.md).

```bash
crevice cavity-trajectory system.gro run.xtc \
    --reference-volume-dx reference_volume.dx --out-dir cavity --workers 4
```

```python
from crevice import load_trajectory
from crevice.cavity_trajectory import (CavityReference, analyze_cavity_trajectory,
                                       write_cavity_trajectory_bundle)

traj = load_trajectory("system.gro", "run.xtc",
                       selection="protein and not name H*", stop=5)
reference = CavityReference.from_dx("reference_volume.dx")
analysis = analyze_cavity_trajectory(traj, reference)
```

In Python the geometry settings go in one dictionary (`geometry=`), unlike the
keyword arguments of the static functions.

## How it works

Each frame is aligned to the first on matching CA atoms, and the cavity is
measured on a lattice fixed to the reference. There are two ways to measure
it:

- **`--geometry-mode rolling`** (the default) casts the whole protein in each
  frame and picks the region that best overlaps the reference. A match needs
  an overlap score of at least `--minimum-overlap` (0.1), and a frame is
  *ambiguous* if a second region scores at least `--ambiguity-ratio` (0.8)
  times the best. Missing or ambiguous matches are recorded as missing, never
  as zero volume.
- **`--geometry-mode reference-region`** measures a fixed local region: every
  bit of free space that touches the reference, within `--region-margin`
  (2 Å) of it. Here a zero means no free space was found in that region.

In every frame CREVICE records the volume, the widths along a fixed axis, each
wall residue's share of the boundary (as in
[boundary residues](residue-evidence.md)), the partner residues and all residue
contacts. Each series gets a mean, spread, frame range and, where it can be
estimated, an approximate interval for the mean.

Two options are worth knowing about:

- **`--water-membership`** also counts the waters inside the cavity as
  measured in each frame, next to the fixed-reference count; see
  [water in the cavity](hydration.md#water-in-the-cavity-two-ways-to-count).
- **`--obstacle-selection`** keeps chosen non-protein atoms, such as lipids,
  as walls (reference-region mode only), including their periodic images; see
  [Physical obstacles and numerical sensitivity](../CAVITY_OBSTACLES.md).

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--reference-volume-dx` | required | the measured reference map, in frame-0 coordinates |
| `--geometry-mode` | `rolling` | `rolling` or `reference-region` |
| `--spacing` | from the reference | grid spacing in Å |
| `--probe-radius`, `--outer-radius` | 0.8, 6.0 Å | the cast probes |
| `--region-margin` | 2.0 Å | how far around the reference to look (reference-region) |
| `--axis` | from the map | the axis for widths, `x,y,z` |
| `--max-frames`, `--workers` | 500, 1 | frame guard and parallel processes |
| `--water-membership` | off | count waters inside the cavity in each frame |

Everything is listed under
[`crevice cavity-trajectory`](../reference/cli/cavity-trajectory.md).

## What you get

- **Per frame**: `*_cavity_frames.csv` (time, volume, alignment RMSD, match
  scores, axial extent), `*_cavity_section_frames.csv` (widths along the
  axis) and `*_cavity_residue_frames.csv` (each wall residue's area).
- **Summaries**: width profiles, section and free-sphere statistics, residue,
  partner and contact tables, and `*_cavity_statistics_summary.csv`, which puts
  every summary statistic on its own row.
- **`*_cavity_statistics.json`**: the structured report, which
  [`hydration --cavity-results`](hydration.md) and
  [`region-compare`](regions.md) read later.
- **Figures**: width profiles (a violet band for the frame range, a violet
  line for the mean and an orange band for the mean interval),
  time-by-position and residue-area heatmaps (unresolved frames in grey), and
  partner and motion summaries.
- Standard hydration outputs unless `--skip-hydration`, and the water
  membership files with `--water-membership`.

:::{admonition} Interpreting results
:class: crevice-interpret

Cavity volumes over a trajectory depend on every geometric choice: in our own
membrane-protein run the regional volume changed with the margin, probe, grid
spacing and grid phase, so check sensitivity on your system
([how](../CAVITY_OBSTACLES.md#checking-numerical-sensitivity)). Leaving out
membrane and ligands can open space that is really blocked. In rolling mode,
cavities that split, merge or grow can break the match. Mean intervals assume
stationary, representative sampling, and one simulation is not a replica set.
A reference pocket is a candidate until it has been reviewed, and a zero
volume in a region doesn't show that transport is closed. The method is
described in [Cavity trajectories](../CAVITY_TRAJECTORY_METHODS.md).
:::

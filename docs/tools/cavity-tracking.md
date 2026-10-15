# Cavity tracking

`crevice track-cavities` follows every cavity near a site through a
trajectory and gives each one an identity that persists from frame to frame.
It records when a cavity appears (birth) or disappears (death), splits into
pieces or merges with another, and reports each identity's volume, position
and lining residues over time. Unlike [`cavity-trajectory`](cavity-trajectory.md),
it needs no reference map and follows all cavities, not one.

## Quick start

```bash
crevice track-cavities glut1_excerpt.gro glut1_excerpt.xtc \
    --region glut1_site_candidate_region.json --radius 8 --workers 4 -o tracking
```

```python
from crevice.tracking import track_cavities, write_tracking_outputs

result = track_cavities("glut1_excerpt.gro", "glut1_excerpt.xtc",
                        region="glut1_site_candidate_region.json", radius=8.0, workers=4)
write_tracking_outputs(result, "tracking", prefix="glut1")
print(len(result["lifetimes"]), [e["event"] for e in result["events"]][:5])
```

## How it works

1. Each frame is aligned to the first on matching CA atoms and cast with a
   rolling probe (1.4 Å inside a 6 Å bulk envelope) on one fixed lattice,
   so voxels of different frames line up exactly.
2. Cavities smaller than `--min-volume` are dropped; with `--centre` or
   `--region`, only cavities within `--radius` of the site are kept.
3. A cavity in one frame is **linked** to a cavity in the next when their
   voxels overlap enough: Jaccard index at least `--link-jaccard` (0.2), or
   at least half of the smaller one inside the larger (`--link-containment`).
4. Each identity continues into its best-overlapping linked cavity. A cavity
   with no link is a birth; one with none in the next frame dies; one linked
   to several cavities splits; several linked to one merge.

Only the previous frame's cavities are kept in memory and the trajectory is
read in chunks, so a 1001-frame run uses no more memory than a short one. The
full definitions are in the [method note](../methods/cavity-tracking.md).

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--region` or `--centre` + `--radius` | none, 8 Å | which cavities to follow |
| `--probe-radius`, `--outer-radius` | 1.4, 6.0 Å | the cast probes |
| `--spacing` | 0.5 Å | lattice spacing |
| `--min-volume` | 5 Å³ | smallest cavity kept |
| `--link-jaccard`, `--link-containment` | 0.2, 0.5 | how much overlap links two cavities |
| `--start`, `--stop`, `--stride` | all frames | frame range |
| `--workers`, `--chunk` | 1, 50 | parallel processes and frames read at a time |
| `--md-selection`, `--include-hetero` | protein, off | lipids and ions as walls |

All options are listed under [`crevice track-cavities`](../reference/cli/track-cavities.md).

## What you get

- **`*_tracking_events.csv`**: every birth, death, split and merge with its
  frame, time, the identities before and after, and the overlaps.
- **`*_tracking_identities.csv`**: one row per identity and frame: volume,
  voxel count, centroid and lining residues.
- **`*_tracking_lifetimes.csv`**: one row per identity: first and last frame,
  lifetime, frames present, how it started (first frame, birth, split) and
  ended (still present, death, merge), mean and maximum volume.
- **`*_tracking_frames.csv`**: per frame, the number and total volume of the
  tracked cavities and that frame's events.
- **`*_tracking.json`**: settings, definitions and event counts.
- **`*_tracking.png`**: volumes over time (log scale; the six most persistent
  identities in colour, the rest grey) above one lane per identity, with
  births, deaths, splits and merges marked.

:::{admonition} Interpreting results
:class: crevice-interpret

Identities follow measured geometry, not biology: a split or merge can come
from a probe-sized neck opening or closing, or from noise in a small cavity,
and frames far apart in time compress many changes into one event. Check how
the event counts change with `--link-jaccard`, `--probe-radius` and
`--spacing` before reading meaning into them.
:::

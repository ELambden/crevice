# Cavity tracking across frames

This note defines what `crevice track-cavities` (Python:
{func}`crevice.tracking.track_cavities` and {class}`crevice.tracking.CavityTracker`)
measures. The tool page is [Cavity tracking](../tools/cavity-tracking.md).

The tracker measures every cavity component of each trajectory frame on one
fixed lattice and links components of consecutive measured frames by voxel
overlap. A chain of links is a persistent **identity**; changes in topology
between frames are recorded as **events**.

```{admonition} What the numbers are
:class: note
Identities are geometric continuity of measured cavity components under the
stated probes, lattice and thresholds. They are not biological states,
transport events or proof that two cavities are connected, and an event can
come from noise in a small or barely connected component.
```

## Frames, lattice and components

| Item | Definition |
|---|---|
| Alignment | Every frame is rigidly fitted to the first frame of the range on matching C-alpha atoms ({func}`crevice.trajectory.align_frame_report`). |
| Cast | {func}`crevice.rolling.rolling_probe_cast` with `selection_mode="all"`: every inner-probe (`--probe-radius`, default 1.4 Å) component inside the rolling bulk-probe envelope (`--outer-radius`, default 6 Å). |
| Lattice | World axes, spacing `--spacing` (default 0.5 Å), with a node at the origin, so the voxels of different (aligned) frames coincide exactly. |
| Component | A cast component with volume ≥ `--min-volume` (default 5 Å³). With `--centre` or `--region`, only components with a voxel within `--radius` (default 8 Å) of the centre (for `--region`, the mean of the landmark centroids in the first frame). |
| Voxel set | The component's measured swept samples. |
| Volume, centroid | The cast volume (Å³); the mean of the voxel centres (Å). |
| Lining residues | Residues with a heavy-atom centre within `--lining-cutoff` (default 3.3 Å) of a voxel of the component (the cast lining definition). |
| Obstacles | Protein heavy atoms; with `--include-hetero`, also hetero atoms (lipids, ions) within the protein box plus 10 Å. |

## Links and identities

Component *A* of the previous measured frame and component *B* of the
current frame are **linked** when

* the Jaccard index `|A ∩ B| / |A ∪ B| ≥ --link-jaccard` (default 0.2), or
* the containment `|A ∩ B| / min(|A|, |B|) ≥ --link-containment` (default
  0.5; the smaller set lies mostly inside the larger one, as for a small piece
  budding off a large cavity; 0 disables this rule).

Links are then taken in order of decreasing Jaccard index (ties: larger
overlap, then component order). A link **continues** an identity when neither
of its components has been used yet, so each identity continues into at most
one component per frame.

## Events

| Event | Rule | Identities |
|---|---|---|
| `birth` | a current component with no link | new identity |
| `death` | a previous component with no link | its identity ends (`end_reason = death`) |
| `split` | a previous component linked to two or more current components | the continued child keeps the parent identity; each other child gets a new identity (`start_reason = split`, `parents`) |
| `merge` | a current component linked to two or more previous components | the continued parent's identity survives; each other parent ends (`end_reason = merge`, `merged_into`) |

Frames whose cast fails are recorded with `status = unresolved` and skipped;
the next measured frame is compared with the last measured one, and the
events record the number of skipped frames (`gap_frames`).

## Memory

Only the previous measured frame's voxel sets are kept, and the trajectory is
read `--chunk` frames at a time (default 50), so memory does not grow with the
number of frames. The tables grow by one row per identity and frame.

## Choosing thresholds

The defaults link a cavity that moves by a fraction of its own size between
frames. With frames far apart in time (for example every 5 ns), a cavity can
move or change shape enough to break the link; the result is a death and a
birth, not an error. Check sensitivity by rerunning with `--link-jaccard 0.1`
and `--link-containment 0`, and with a different `--probe-radius`: small
components near the probe size appear and vanish with small changes.

## Limitations

* Rigid C-alpha alignment: domain motions move cavities relative to the
  fitted frame and lower overlaps.
* Identities depend on the probe, the lattice and the neighbourhood; a
  component that grows beyond the neighbourhood is still tracked, but one
  that appears outside it is not.
* Events between widely spaced frames summarise everything that happened in
  between; they are not instantaneous transitions.

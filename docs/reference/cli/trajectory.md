# crevice trajectory

Frame series and MD trajectories: channel profiles (default), measured cavities (`--cavities`), a prepared named region (`--region DEFINITION`), water and hydration (`--water`) and cavity tracking (`--track`).

Choose one mode flag. `--track` and `--region` may appear together: cavity tracking has its own `--region DEFINITION` option (track around the region's landmark anchor), so `--track` wins and `--region` is passed to it.

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice trajectory`

Analyze a sequence of static structure frames. For the method, worked commands and outputs, see
[Channel profiles over a trajectory](../../tools/trajectory.md).

```{crevice-cli} trajectory
```

## `crevice trajectory --track`

*CREVICE 0.3: `crevice track-cavities`.*

Track cavity components across trajectory frames by voxel overlap, with births, deaths, splits and merges. For the method, worked commands and outputs, see
[Cavity tracking](../../tools/cavity-tracking.md).

```{crevice-cli} trajectory --track
```

## `crevice trajectory --cavities`

*CREVICE 0.3: `crevice cavity-trajectory`.*

Track measured cavities in every aligned trajectory frame; volumes, changing residue evidence and sectional width statistics. For the method, worked commands and outputs, see
[Cavities over a trajectory](../../tools/cavity-trajectory.md).

```{crevice-cli} trajectory --cavities
```

## `crevice trajectory --water`

*CREVICE 0.3: `crevice hydration`.*

Observed residue water contacts, geometric SASA and sampled trajectory hydration. For the method, worked commands and outputs, see
[Hydration](../../tools/hydration.md).

```{crevice-cli} trajectory --water
```

## `crevice trajectory --region DEFINITION`

*CREVICE 0.3: `crevice region-trajectory`.*

Run geometry and hydration using a prepared, versioned region definition. For the method, worked commands and outputs, see
[Named regions](../../tools/regions.md).

```{crevice-cli} trajectory --region
```

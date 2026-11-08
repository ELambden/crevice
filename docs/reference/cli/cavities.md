# crevice cavities

Cavities and casts in one structure. The default mode finds buried cavities on a clearance grid; `--cast` fills pockets, cavities and channels with a 3D rolling probe. Routes from a site to bulk solvent are found by [`crevice paths`](paths.md).

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice cavities`

Find buried cavities on a clearance grid. For the method, worked commands and outputs, see
[Buried cavities](../../tools/cavities.md).

```{crevice-cli} cavities
```

## `crevice cavities --cast`

*CREVICE 0.3: `crevice cast`.*

Fill protein-interior voids with a 3D rolling probe, including open pockets and multiple regions. For the method, worked commands and outputs, see
[Cavity and channel casts](../../tools/cast.md).

```{crevice-cli} cavities --cast
```

# crevice fetch

Download structures, or (with `--example`) a hash-checked example dataset.

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice fetch`

Fetch structures from RCSB by accession, or AlphaFold DB models by AF- identifier. For the method, worked commands and outputs, see
[Utility and batch commands](../../tools/utilities.md).

```{crevice-cli} fetch
```

## `crevice fetch --example NAME`

*CREVICE 0.3: `crevice fetch-example`.*

Download a named example dataset. Every file is checked against its SHA-256 hash and
cached locally. The datasets are excerpts of membrane simulations:

| Name | Contents | Walk-through |
| --- | --- | --- |
| `glut1-excerpt` | GLUT1 transporter, 21 frames over 100 ns, with a prepared region | [A membrane MD trajectory](../../examples/glut1-trajectory.md) |
| `aqp4-excerpt` | AQP4 water channel, 200 consecutive frames, waters near the pore | [Water through an aquaporin](../../examples/aqp4-water-permeation.md) |
| `ompf-excerpt` | OmpF porin trimer in 1 M KCl under a field, 100 consecutive frames, every ion | [Ions through a porin](../../examples/ompf-ion-permeation.md) |

```{crevice-cli} fetch --example
```

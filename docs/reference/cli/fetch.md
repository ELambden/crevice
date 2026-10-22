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

Download a named example dataset, such as the GLUT1 trajectory excerpt used in the
tutorials. Every file is checked against its SHA-256 hash and cached locally. For a
walk-through, see [A membrane MD trajectory](../../examples/glut1-trajectory.md).

```{crevice-cli} fetch --example
```

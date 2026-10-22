# crevice report

Complete outputs for one structure (default: the publication bundle), selected analyses (`--analyze`), or the benchmark systems (`--suite`, `--benchmark`).

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice report`

*CREVICE 0.3: `crevice publish`.*

Generate a publication-ready static-structure output bundle. For the method, worked commands and outputs, see
[One-command output bundle](../../tools/publish.md).

```{crevice-cli} report
```

## `crevice report --analyze`

*CREVICE 0.3: `crevice analyze`.*

Run selected analyses on one structure. For the method, worked commands and outputs, see
[Utility and batch commands](../../tools/utilities.md).

```{crevice-cli} report --analyze
```

## `crevice report --suite`

*CREVICE 0.3: `crevice static-suite`.*

Run publication-output static analyses across benchmark systems. For the method, worked commands and outputs, see
[Utility and batch commands](../../tools/utilities.md).

```{crevice-cli} report --suite
```

## `crevice report --benchmark`

*CREVICE 0.3: `crevice benchmark`.*

Run static checks on configured CREVICE PDB benchmarks. For the method, worked commands and outputs, see
[Utility and batch commands](../../tools/utilities.md).

```{crevice-cli} report --benchmark
```

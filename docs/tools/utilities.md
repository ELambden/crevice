# Utility and batch commands

A handful of commands help you get structures, run several analyses at once,
export features and run CREVICE across a set of structures.

## Downloading structures: `fetch`

```bash
crevice fetch 1GRM 4PYP AF-P02920-F1 --cache-dir .crevice/pdb
```

`fetch` downloads entries from the RCSB PDB, and models from AlphaFold DB by
their `AF-` identifier, into a local cache. Each download is checked and gets
a provenance record next to it. You rarely need it on its own, because every
command fetches what it needs, but it is handy before working offline. See
[fetching structures](../getting-started/inputs-and-selections.md#fetching-structures).
`crevice fetch --list` prints the bundled benchmark systems.
[All options](../reference/cli/fetch.md).

## Several analyses at once: `analyze`

```bash
crevice analyze 1GRM --analysis profile --analysis residues --analysis network \
    --out-dir results/1GRM_analyze --png
```

`analyze` runs any mix of `profile`, `residues`, `cavities`, `tunnels`,
`network` and `features` on one structure (repeat `--analysis` for each; the
default is profile and residues), writing each as CSV (with the JSON
records beside them) plus a manifest. `--png` adds figures for the profile,
residues, cavities and network. Unlike the standalone commands, its cavity
and tunnel grid caps (`--cavity-max-grid-points`, `--tunnel-max-grid-points`)
coarsen the spacing to stay under the cap, so raise them if you need a fine
grid. [All options](../reference/cli/analyze.md).

## Feature tables: `features`

```bash
crevice features 1GRM -o 1GRM_features.json
```

`features` exports numerical descriptors of a profile and its residues as JSON
and as `<stem>.csv` (`feature`, `value`, `unit`), ready for clustering or
your own machine-learning experiments (see {mod}`crevice.ml`). The features
describe geometry; CREVICE ships no trained model.
[All options](../reference/cli/features.md).

## Batch runs: `benchmark` and `static-suite`

`benchmark` runs quick checks, and `static-suite` full `publish`-style
bundles, on the bundled benchmark systems (1GRM, 6MVY, 1AF6, 4PYP, 2CHB) or the
ones you name.

```bash
crevice benchmark --fetch --output benchmark.json
crevice static-suite 1GRM 4PYP --fetch --out-dir suite
```

`benchmark --output b.json` also writes `b.csv`, one row per system with the
profile radii, top cavity volume and any error. `static-suite` writes
`static_benchmark_summary.csv` and `.json` and one bundle per system, with the
same tables as `publish` plus `*_profile_sensitivity.csv`. Both choose the
enclosure probe automatically for each system and record it
(`enclosure_radius_A`, `enclosure_probe_mode`); a system where no probe gives
a stable channel becomes an error row naming the probes tried. Pass a number,
such as `--enclosure-radius 0.8`, to use one probe everywhere.

The sensitivity rows of `static-suite` vary the *measurement* probe
(`--probe-radius`: 0, 1.0 and 1.4 Å by default) and keep each system's
enclosure probe. A measurement probe too wide for the channel (1.4 Å in 1GRM,
whose narrowest radius is about 1.33 Å) gives an `unresolved` row rather than
failing the system. Options: [`benchmark`](../reference/cli/benchmark.md),
[`static-suite`](../reference/cli/static-suite.md).

:::{admonition} Interpreting results
:class: crevice-interpret

The benchmark labels are provisional settings, not curated biological
classifications: `crevice.benchmarks.curated_benchmark_ids()` is empty, every
summary has a `curation_status` column, 4PYP (human GLUT1) has no assigned
cavity mode and 2CHB is set aside. A suite run checks that the software and
its outputs behave consistently; it isn't validation against reference
geometry.
:::

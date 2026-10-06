# Utility and batch commands

## `fetch`

Downloads structures from RCSB into a local cache. Each download is checked and
gets a provenance sidecar. See
[Inputs and selections](../getting-started/inputs-and-selections.md#fetching-structures).
`crevice fetch --list` prints the configured benchmark systems with their
`curation_status`.

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: fetch
```

## `analyze`

Runs any mix of `profile`, `residues`, `cavities`, `tunnels`, `network` and
`features` on one structure. It writes a manifest naming the outputs, the input
provenance and the parser report. Tables are CSV: `*_profile.csv`,
`*_residues.csv`, `*_cavities.csv`, `*_tunnels.csv`, `*_tunnel_points.csv`,
`*_network_nodes.csv`, `*_network_edges.csv` and `*_features.csv` (`feature`,
`value`, `unit`); the profile, cavity, tunnel and network JSON records are kept
beside them. `--png` adds figures for the profile, residues, cavities and
network (tunnels have tables only). The cavity and tunnel grid caps
(`--cavity-max-grid-points`, `--tunnel-max-grid-points`) **coarsen the requested
spacing** to stay under the cap. Raise the caps to keep a fine spacing.

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: analyze
```

## `features`

Exports numerical profile and residue features as JSON (with the profile and
top residues) and as `<stem>.csv` (`feature`, `value`, `unit`) for clustering
or machine-learning experiments (see {mod}`crevice.ml`). The features are
descriptive. No trained model ships with CREVICE.

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: features
```

`analyze` and `features` accept `--radii` (atomic radius set; default the standard table, see [Atomic radii](../methods/atomic-radii.md)).

## `benchmark` and `static-suite`

Run quick checks (`benchmark`) or full publication bundles (`static-suite`) on
the configured benchmark systems (1GRM, 6MVY, 1AF6, 4PYP, 2CHB).

:::{warning}
The benchmark labels are provisional configuration, not curated biological
classifications. `crevice.benchmarks.curated_benchmark_ids()` returns an empty
tuple, and every summary carries a `curation_status` column. 4PYP is human
GLUT1 and has no assigned cavity mode (`unassigned`). 2CHB (cholera toxin
B-pentamer) is set aside. A suite run is a
software regression and presentation check. It is not validation against
reference geometry.
:::

`benchmark --output b.json` also writes `b.csv` (one row per system:
`profile_min_radius_A`, `profile_mean_radius_A`, `top_cavity_volume_A3`, ...,
`error`). `static-suite` writes `static_benchmark_summary.csv` (units in the
column names, an `error` column for failed systems), the suite record
`static_benchmark_summary.json` and one bundle per system with the same CSV
tables as `publish` plus `*_profile_sensitivity.csv` (one row per probe
radius). No tunnel figure is written.
Both accept `--radii` (atomic radius set); a chosen set is recorded in the
summary JSON and in each system manifest.

Both commands use the automatic enclosure probe by default
(`--enclosure-radius auto`; see
[`profile`](profile.md#choosing-the-enclosure-probe---enclosure-radius-auto)),
chosen separately for each system. The chosen probe is recorded per system:
`enclosure_radius_A` and `enclosure_probe_mode` (`auto`, `explicit`, or
`not_used` with `--search-radius 0`) in both CSVs, plus
`enclosure_probe_reason` in the `benchmark` rows and the full `enclosure_probe`
record (with the selection) in the `static-suite` system summaries. If no
probe gives a stable channel, the system is an error row naming the probes
tried. A number (`--enclosure-radius 0.8` reproduces the outputs recorded
before 1 October 2026) runs that probe for every system.

In `static-suite`, the sensitivity rows vary the **measurement** probe
(`--probe-radius`, default 0, 1.0 and 1.4 Å, subtracted from each radius)
and keep the system's enclosure probe (`enclosure_radius_A` column). A
measurement probe that does not fit through the channel (for example 1.4 Å in
1GRM, whose narrowest radius is about 1.33 Å) gives a row with
`status: unresolved` and the reason, listed in
`profile_sensitivity_unresolved_probe_radii`; before 1 October 2026 such a row
made the whole system fail.

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: benchmark
```

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: static-suite
```

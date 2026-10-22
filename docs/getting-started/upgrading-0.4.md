# Upgrading to 0.4

CREVICE 0.4 groups the 25 commands of 0.3 into nine: `fetch`, `profile`,
`cavities`, `paths`, `residues`, `report`, `trajectory`, `region` and `view`.
Commands that answer related questions became **modes** of one command,
chosen with a flag (`crevice cavities --cast`) or, for `region`, an action word
(`crevice region prepare`).

The new command layout changes no result. Each mode runs the 0.3 command's own
code with the same options and the same defaults, and writes the same files.
Only `--legacy-outputs` is gone (see [below](#legacy-outputs-removed)). Two
measurement changes in 0.4 do change some numbers; see
[Results that change](#results-that-change).

## Old command names still work

Every 0.3 command name is a **deprecated alias**: it runs exactly the same
code as its new form and prints one note on standard error, for example

```text
crevice: note: 'crevice cast' is deprecated; use 'crevice cavities --cast' (same options and outputs; the 0.3 names are removed in 1.0)
```

Standard output and every output file are byte-for-byte the same for the old
name and the new form; only the recorded command line (`run.argv` in the
manifest) differs. This is checked for every command by
`tests/test_command_aliases.py`. The aliases will be removed in CREVICE 1.0,
so update scripts when convenient.

## Command mapping

| CREVICE 0.3 | CREVICE 0.4 |
|---|---|
| `crevice fetch IDS` | `crevice fetch IDS` (unchanged) |
| `crevice fetch-example NAME` | `crevice fetch --example NAME` |
| `crevice fetch-example --list` | `crevice fetch --example --list` |
| `crevice profile` | `crevice profile` (unchanged) |
| `crevice fingerprint` | `crevice profile --fingerprint` |
| `crevice cavities` | `crevice cavities` (unchanged) |
| `crevice cast` | `crevice cavities --cast` |
| `crevice tunnels` | `crevice cavities --tunnels` |
| `crevice paths` | `crevice paths` (unchanged) |
| `crevice residues` | `crevice residues` (unchanged) |
| `crevice network` | `crevice residues --network` |
| `crevice residue-evidence` | `crevice residues --evidence` |
| `crevice features` | `crevice residues --features` |
| `crevice publish` | `crevice report` |
| `crevice analyze` | `crevice report --analyze` |
| `crevice static-suite IDS` | `crevice report --suite IDS` |
| `crevice benchmark IDS` | `crevice report --benchmark IDS` |
| `crevice trajectory` | `crevice trajectory` (unchanged) |
| `crevice cavity-trajectory TOP TRAJ` | `crevice trajectory --cavities TOP TRAJ` |
| `crevice region-trajectory DEFINITION` | `crevice trajectory --region DEFINITION` |
| `crevice hydration STRUCTURE` | `crevice trajectory --water STRUCTURE` (a single structure without `--trajectory`, as before) |
| `crevice track-cavities TOP [TRAJ]` | `crevice trajectory --track TOP [TRAJ]` |
| `crevice region-init` | `crevice region init` |
| `crevice region-prepare` | `crevice region prepare` |
| `crevice region-compare` | `crevice region compare` |
| `crevice view` | `crevice view` (unchanged) |

All other arguments stay as they were, in any order; a mode flag may appear
anywhere after the command name (`crevice cavities 4PYP --cast --out-dir out`
works too).

## How modes work

- **Each mode keeps its own options and defaults.** For example
  `crevice cavities` (clearance grid) uses `--spacing 2.0`, while
  `crevice cavities --cast` (rolling-probe cast) uses `--spacing 0.5`;
  `--max-cavities` is 10 for the grid search and 1 for the cast. Nothing was
  merged or renamed, so no default changed. `crevice cavities --help` shows
  the default mode and lists the others; `crevice cavities --cast --help`
  shows the cast options.
- **One mode at a time.** Two mode flags of the same command
  (`crevice cavities --cast --tunnels`) are an error.
- **`--example NAME` and `--region DEFINITION` take a value**: it is the
  dataset name or the prepared region definition, the first input of that
  mode.
- **`--track` with `--region`.** Cavity tracking has its own `--region
  DEFINITION` option (track around the region's landmark anchor). When both
  appear, `--track` selects the mode and `--region` is passed to it, as in
  0.3: `crevice trajectory --track top.gro run.xtc --region site.json -o out`.
- **`--fingerprint` in two places.** `crevice profile --fingerprint` is the
  standalone fingerprint (0.3 `fingerprint`); `crevice report --fingerprint`
  still adds the fingerprint to the report bundle (0.3 `publish
  --fingerprint`).
- **Flag names.** `crevice residues --features` exports the numeric profile
  and residue descriptors (0.3 `features`). It is not called `--interactions`,
  because "interactions" already names the typed residue–residue interaction
  analysis of the bundle commands (`--outputs interactions`,
  `--skip-interactions`). `crevice report --suite` takes benchmark IDs, as
  `static-suite` did; it does not read a file of IDs.

## Settings files and Python

- Settings files written by `--print-config` in 0.3 keep working. Their
  `"command"` is the 0.3 name of the analysis (`"publish"`, `"cast"`,
  `"region-trajectory"` ...) and 0.4 still writes it that way, because the
  same name keys the output policy and the run record. The 0.4 spelling
  (`"report"`, `"cavities --cast"`, `"trajectory --region"`) is accepted too.
  No settings key was renamed.
- The only removed key is `outputs.legacy_outputs`: a 0.3 file that records it
  as `false` loads unchanged (the key is ignored); `true` is refused with a
  message pointing to `--all`.
- {func}`crevice.config.settings` and {class}`crevice.config.RunSettings`
  accept both spellings: `settings("cavities --cast", "4PYP", out_dir="out")`
  and `settings("cast", "4PYP", out_dir="out")` are the same.
- Manifests still record `"command"` as the analysis name, so results of 0.3
  and 0.4 runs can be compared field by field.

(legacy-outputs-removed)=
## Results that change

- **Channel profile radii** are measured at the exact top of each slice's
  clearance hill, so radii, minimum radius and bottleneck increase slightly
  (1GRM minimum radius 1.337 Å instead of 1.330 Å with `report`, 1.337 Å
  instead of 1.324 Å with `crevice profile 1GRM`). Results built on the
  profile follow: the fingerprint conductance estimate, hydration regions
  that follow the profile, and lateral exit legs with their cast segments.
  Channel casts are unchanged. See [Slice centres](../tools/profile.md#slice-centres).
- **HETATM residues in a polymer chain** (D-amino acids, selenomethionine,
  modified residues) are protein for PDB files, as they already were for
  mmCIF files. PDB files with such residues give different atom sets and
  results (1GRM: `cavities --cast 1GRM.pdb` gives 257.000 Å³ like the mmCIF
  file, 181.875 Å³ before). `--hetero-polymer exclude` restores the 0.3
  behaviour.

## `--legacy-outputs` removed

`--legacy-outputs` reproduced the CREVICE 0.2 file set during the 0.3 series
and was announced for removal. It is gone in 0.4, together with the Python
`outputs="legacy"` value. Use `--all` (Python `outputs="all"`) for every file
kind and optional analysis, or list what you need with `--outputs`
([Outputs](outputs.md)). Measurements were never affected by it.

## Updating scripts

Replace the command names using the table above. The commands' outputs are
unchanged, so results of earlier runs remain valid and comparable. The
worked-examples script (`examples/run_worked_examples.sh`) and every tutorial
in this documentation use the 0.4 forms.

# Shared options

These option families work the same way in every command that accepts them,
so they are described here once. Each command page lists which families it
takes.

(cli-structure-input)=
## Structure input

The structure argument can be a local PDB, mmCIF or GRO file, an RCSB
accession such as `1GRM`, or an AlphaFold DB identifier such as
`AF-P02920-F1`. Accessions are downloaded once and cached; if a file with the
same name exists locally, the file wins. See
[Inputs and selections](../../getting-started/inputs-and-selections.md) for
how atoms are chosen.

```{crevice-cli-options} structure-input
```

(cli-radii)=
## Atomic radii

One radius set is used for every clearance, cast, contact and surface-area
calculation in a command, including in worker processes, and it is recorded in
the output metadata. See [Atomic radii](../../methods/atomic-radii.md).

```{crevice-cli-options} radii
```

(cli-figures)=
## Figures and viewer scenes

Figures and PyMOL/VMD/ChimeraX scenes are clean by default: no titles,
captions or labels. The text is still stored in each PNG's metadata and in the
scene JSON, so `--annotate` only changes what is drawn. Display smoothing
changes how a surface looks, never what was measured
([Surface smoothing](../../SURFACE_SMOOTHING.md)).

```{crevice-cli-options} figures
```

(cli-channel-profile)=
## Channel profile settings

Commands that measure a through-channel profile share these settings. They
are explained with examples on the [`profile` page](../../tools/profile.md).

```{crevice-cli-options} channel-profile
```

(cli-hydration)=
## Hydration settings

Commands that analyse explicit water share these settings; see
[Hydration](../../tools/hydration.md).

```{crevice-cli-options} hydration
```

(cli-outputs)=
## Output selection

The bundle commands (`report`, `cavities --cast`, `report --analyze`, `residues --evidence`,
`trajectory --water`, `trajectory`, `trajectory --cavities`, `trajectory --region`) write a
compact default set: data, the manifest, key figures and PyMOL, VMD and
ChimeraX scenes, plus the typed residue interactions. `--outputs` chooses file
kinds and optional analyses (`--outputs help` lists the values and defaults of
a command), `--all` writes everything, `--dry-run` prints the plan without
reading or writing anything (`--legacy-outputs` was removed in 0.4). See
{ref}`Outputs <output-policy>`
and [Upgrading to 0.4](../../getting-started/upgrading-0.4.md).

```{crevice-cli-options} outputs
```

(cli-settings)=
## Settings files

Every command accepts `--print-config`, which prints its effective settings as
JSON (grouped into `input`, `geometry`, `radii`, `outputs`, `display`,
`trajectory`, `hydration` and `method`) and exits, and `--config FILE`, which
reads such a file; flags on the command line override it. Settings are
validated (units, ranges, files, combinations) before any input is read, and
the manifest records each value with its source. Python:
{class}`crevice.config.RunSettings`, {func}`crevice.config.settings`,
{func}`crevice.config.run`.

```{crevice-cli-options} settings
```

(cli-progress-errors)=
## Progress and errors

Every command accepts `-q/--quiet` and `--debug`, before or after the command
name.

- Analyses print a few status lines on standard error that say what they are
  doing and what they found, for example:

  ```text
  $ crevice profile kcsa_pore.pdb --axis z --origin 155.33,155.33,-43.614 --section-spacing 0.25 --radii bondi -o k.json
  crevice: read kcsa_pore.pdb: 3060 atoms in 0.0 s
  crevice: initialising the enclosure probe for kcsa_pore.pdb (per-plane probe, 0.30-3.00 A, section grid 0.25 A)
  crevice: pore radius found in 0.8 s: minimum 0.747 A at s = +5.0 A (no lower than 0.747 A between samples); mean 2.423 A over 81 samples, s = -36.9 to +22.9 A; per-plane enclosure probes 0.30-3.00 A
  crevice: profile: finished in 0.8 s
  ```

  (`kcsa_pore.pdb` is the KcsA tetramer of PDB 1K4C, protein atoms only.)

  A profile reports its sampled minimum and, in brackets, the exact minimum
  along the path between samples. A `trajectory` reports the mean of the
  frame minima with its blocking standard error, and the effective sample
  size, when the frames allow an estimate (see
  [Uncertainty](../../methods/uncertainty.md)); otherwise it says why it gives
  no uncertainty. Status lines are also written when standard error goes to a
  file, so batch logs keep them. `fetch` and `view` print none.
- Long runs (the frames of `trajectory`, `trajectory --cavities`,
  `trajectory --region` and `trajectory --water --trajectory`) also show a
  progress bar on standard error while it is a terminal. The bar is not drawn when standard error is redirected to a
  file or a pipe; set `CREVICE_PROGRESS=1` to log plain progress lines (every
  10%) instead.
- `--quiet` turns off the status lines and the progress bar;
  `CREVICE_PROGRESS=0` does the same everywhere. Neither ever changes standard
  output or any output file.
- An expected failure (a missing file, an unknown accession, a failed
  download, a channel that cannot be resolved) is reported as one
  `crevice: error: ...` line, followed for common cases by a
  `crevice: hint: ...` line with what to try, and the command exits with
  status 2. `--debug` (or `CREVICE_DEBUG=1`) also prints the full Python
  traceback, which is useful in a bug report.

(cli-output-files)=
## Output files

Output options differ from command to command, so they are listed on each
command page, but they follow the same pattern:

- `-o/--output` names the main JSON result. Tables are written as CSV beside
  it (`<stem>.csv`) unless the command says otherwise.
- `--png`, `--pdb` and `--pymol` add a figure, a dummy-atom trace and a PyMOL
  script. `--dpi` sets the figure resolution: 200 dots per inch by default,
  which is quick and keeps files small; use `--dpi 400` (Python `dpi=400`)
  for print or publication figures. Native-viewer renders are sized
  separately (see the `crevice_render` helpers).
- Commands that write a whole set of files, such as `report`, `cavities --cast` and the
  trajectory commands, take `--out-dir` and record what they wrote, with the
  settings used, in a `*_manifest.json`; [Output selection](#cli-outputs)
  chooses which files they write.

[Outputs](../../getting-started/outputs.md) describes the files in more detail.

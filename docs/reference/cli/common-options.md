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

- Long runs (the frames of `trajectory`, `trajectory --cavities`,
  `trajectory --region` and `trajectory --water --trajectory`, and the automatic
  enclosure-probe ladder) show progress on standard error while it is a
  terminal. Nothing is printed when standard error is redirected to a file or
  a pipe, or with `--quiet`; set `CREVICE_PROGRESS=1` to log plain progress
  lines (every 10%) to a file, or `CREVICE_PROGRESS=0` to silence progress
  everywhere. Progress never changes standard output or any output file.
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

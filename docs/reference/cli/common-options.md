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

(cli-output-files)=
## Output files

Output options differ from command to command, so they are listed on each
command page, but they follow the same pattern:

- `-o/--output` names the main JSON result. Tables are written as CSV beside
  it (`<stem>.csv`) unless the command says otherwise.
- `--png`, `--pdb` and `--pymol` add a figure, a dummy-atom trace and a PyMOL
  script. `--dpi` sets the figure resolution.
- Commands that write a whole set of files, such as `publish`, `cast` and the
  trajectory commands, take `--out-dir` and record what they wrote, with the
  settings used, in a `*_manifest.json`.

[Outputs](../../getting-started/outputs.md) describes the files in more detail.

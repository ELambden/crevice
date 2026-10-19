# Changelog

What changed in each release. The complete list, including known limitations,
is in [CHANGELOG.md](../../CHANGELOG.md) in the repository.

## 0.3.3 (19 October 2026)

- **Windows.** The test suite passes on Windows (two tests compared paths with
  backslashes incorrectly; CREVICE itself was unaffected).
- **`crevice view`** opens an output directory written on another operating
  system, for example a Windows bundle copied to Linux or macOS.

## 0.3.2 (19 October 2026)

**Faster, more portable, more transparent; no change to any result.**

- **Faster.** `region-trajectory` is about 1.8× faster, `paths` 1.2–1.6×,
  `publish` 1.3–1.4× and `profile --lateral-exits` 1.25× on one core, with
  every output file byte for byte the same. `scripts/timing.py` in the
  repository reproduces the timings.
- **Routes left out by the automatic probe.** `crevice paths` now says how
  many routes the automatic probe left out and which `--probe-radius` shows
  them, and lists them in `*_path_left_out.csv`. See
  [Curved and branched pathways](../tools/paths.md).
- **macOS and Windows** are tested automatically, parallel workers use the
  start method those systems need, and text files are UTF-8 everywhere.
- **[Limitations](../getting-started/limitations.md)**: one page with the
  caveats that apply across the tools.
- **Fixes.** The fingerprint's nearest lining residue no longer changes between
  runs when two residues are equally close; `--help` no longer fails on
  consoles that cannot print symbols such as ≤; `region-prepare` handles files
  on another Windows drive.

## 0.3.1 (16 October 2026)

- **Automatic probe for `crevice paths`.** `--probe-radius auto` (the default)
  tries the probes of the enclosure ladder (0.5–3.0 Å) in parallel and keeps
  the smallest probe that gives the most stable set of routes, per side with
  `--axis`. The choice and every candidate are printed and recorded in
  `*_path_probes.csv`; a number such as `--probe-radius 0.8` gives the fixed
  probe exactly as before. See [Curved and branched pathways](../tools/paths.md).
- **Fixes.** `crevice view` on a data-only trajectory directory explains how to
  write a scene; the `publish` summary counts the files actually written.
- **Documentation** describes the 0.3 defaults on every tool page and tutorial.

## 0.3.0 (15 October 2026)

**Commands write fewer files by default; no measurement changes.**

- **Compact default outputs.** `publish`, `cast`, `analyze`,
  `residue-evidence`, `hydration`, `trajectory`, `cavity-trajectory` and
  `region-trajectory` write the data, the key figures and one scene each for
  PyMOL, VMD and ChimeraX. A 1GRM `publish` writes 40 files instead of 96,
  with identical numbers. See [Outputs](../getting-started/outputs.md).
- **Water analyses are opt-in** in those commands (`--outputs default,hydration`
  or `--all`), except in `hydration` itself and with `--water-membership`. The
  typed residue interactions still run by default; `--skip-hydration` no
  longer removes them (use `--skip-interactions`).
- **`publish` extras are opt-in**: the cavity and tunnel searches, the residue
  network, the annotated radius figure, input copies and alias scenes.
- **Choosing outputs.** `--outputs LIST`, `--outputs help`, `--all` and
  `--dry-run` on those commands; `--legacy-outputs` (deprecated) writes the
  exact 0.2 file set during the 0.3 series. See
  [Upgrading from 0.2 to 0.3](../getting-started/upgrading-0.3.md).
- **Settings files.** `--print-config` and `--config FILE` (JSON) on every
  command, settings validated before any input is read, and a `run` record
  (version, command line, settings and their sources, input hashes, output
  plan) in each manifest. Python: {mod}`crevice.config`, {mod}`crevice.outputs`
  and {mod}`crevice.static_workflow`.
- **Curved and branched pathways.** `crevice paths` finds the widest routes
  from a seed, residues or a region to bulk solvent without an assumed axis,
  with junctions, lining residues and per-branch scenes. See
  [Curved and branched pathways](../tools/paths.md).
- **Cavity tracking.** `crevice track-cavities` gives cavities persistent
  identities across trajectory frames and records births, deaths, splits and
  merges. See [Cavity tracking](../tools/cavity-tracking.md).
- **Documentation** updated throughout for the new defaults, with a GLUT1
  tutorial for paths and tracking.

## 0.2.0 (14 October 2026)

- **Interactive 3D views.** `crevice view OUTDIR` writes one self-contained,
  offline HTML file with the protein, cast segments, lining residues, centre
  line and member waters as layers with checkboxes, and {func}`crevice.view`
  shows the same scene in Jupyter (new extra `crevice[view]`). The surfaces
  are the native display surfaces. See [Interactive 3D views](../tools/view.md).
- **Pore fingerprints.** `crevice fingerprint`, `publish --fingerprint` and
  {mod}`crevice.fingerprint` put the radius, Wimley–White lining
  hydrophobicity, segment chemistry, water along the path (trajectories) and a
  HOLE-style geometric conductance estimate on one coordinate along a channel.
  See [Pore fingerprints](../tools/fingerprint.md).
- **Faster, with identical results.** Channel profiles, the automatic
  enclosure-probe ladder and the side-exit search are several times faster,
  with byte-identical outputs. `cavity-trajectory` and `region-trajectory` use
  up to four worker processes by default.
- **Testing.** `slow`, `native` and `realdata` test markers, a fast subset
  (`pytest -m "not slow"`), parallel runs with pytest-xdist (`-n auto`) and a
  fast pull-request CI job.

## 0.1.2 (13 October 2026)

- Figures default to 200 dpi everywhere (use `--dpi 400` for print).
- Progress on standard error for long runs, `-q/--quiet`, and `--debug` for
  full tracebacks; expected errors print as one `crevice: error:` line with a
  hint.
- Clearer messages when a download fails.

## 0.1.1 (7 October 2026)

- Example data with `crevice fetch-example`.
- Channel casts split into ENTRY, LUMEN and EXIT segments, with lining
  residues per segment, and side-exit legs cast with `--lateral-exits`.
- An automatic enclosure probe for `profile` and every command that builds a
  channel.
- Atomic radius sets (`--radii`) throughout, recorded in every output.
- AlphaFold DB identifiers as structure input.
- CSV-first outputs; the reorganised documentation site.

## 0.1.0

The first public release, with:

- reading PDB, mmCIF and GRO files, downloading from the RCSB PDB and
  AlphaFold DB with provenance and cache checks, and reading MD trajectories
  through MDAnalysis;
- channel radius profiles, 3D cavity casts, cavity and tunnel searches,
  lining residues and residue interaction networks;
- the residues that form a measured cavity wall and their partners, with a
  choice of which residues are drawn as sticks (`residue-evidence
  --stick-residues`);
- channel profiles and cavity statistics over trajectories, and named,
  versioned regions;
- residue hydration, water density, typed interactions and water inside
  cavities frame by frame (`--water-membership`);
- the one-command `publish` bundle, including a clear fallback when no channel
  resolves;
- display smoothing (`--smooth X`) that never changes measurements;
- PyMOL, VMD and ChimeraX scenes, and this documentation site.

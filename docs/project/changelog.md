# Changelog

What changed in each release. The complete list, including known limitations,
is in [CHANGELOG.md](../../CHANGELOG.md) in the repository.

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

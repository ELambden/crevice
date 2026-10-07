# Changelog

What changed in each release. The complete list, including known limitations,
is in [CHANGELOG.md](../../CHANGELOG.md) in the repository.

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

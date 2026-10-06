# Changelog

The authoritative changelog is [CHANGELOG.md](../../CHANGELOG.md) at the
repository root, in [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
format.

## 0.1.0

The first public release. It provides:

- structure intake (PDB, mmCIF, GRO) with RCSB download, provenance and cache
  verification, and MD trajectory intake through MDAnalysis;
- channel radius profiles, 3D cavity casts, cavity and tunnel searches, residue
  annotation and residue interaction networks;
- measured-boundary residue evidence with configurable viewer sticks
  (`residue-evidence --stick-residues`);
- trajectory profile statistics, all-frame cavity statistics and named-region
  workflows;
- residue hydration, water density, typed interactions and instantaneous-cavity
  water membership (`--water-membership`);
- a one-command `publish` bundle, including an explicit fallback when no
  through-profile resolves;
- constrained display smoothing (`--smooth X`) that never changes measurements;
- PyMOL, VMD and ChimeraX scene scripts, and this documentation site.

See the repository changelog for the complete list and known limitations.

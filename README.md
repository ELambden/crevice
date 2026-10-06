# CREVICE

**Cavity and Residue Environment Visualisation, Interaction and Connectivity
Evaluation**

CREVICE is a Python package and command-line tool for measuring pores,
channels, tunnels and cavities in protein structures and molecular dynamics
(MD) trajectories, and for describing the residues and water that surround
them. It combines ideas from HOLE, HOLLOW, CAVER and 3V in one workflow and
writes numerical tables, figures and matching PyMOL, VMD and ChimeraX scenes.

With CREVICE you can:

- measure **through-channel radius profiles** with continuous atom-clearance
  checks and an automatically proposed or user-supplied axis;
- fill **pockets, cavities and channels** with atom-clear probe spheres and
  export the measured binary map and a smoothed display surface;
- search for **enclosed cavities** and **widest-path tunnels** on clearance
  grids;
- identify the **residues that form a measured boundary**, their nonlocal
  contact partners and residue interaction networks;
- follow a cavity through an **MD trajectory**: aligned per-frame volumes,
  sectional widths, changing boundary residues and contacts, with approximate
  uncertainty on means;
- measure **hydration**: explicit water contacts, solvent-accessible area,
  water density, typed interactions and water inside the instantaneous cavity;
- define, version and compare **named regions** over identical frames.

Every analysis runs from the `crevice` command or from Python.

## Installation

CREVICE needs Python 3.10 or newer.

```bash
python -m pip install crevice
```

NumPy, SciPy, MDAnalysis, Matplotlib and scikit-image are installed
automatically. Optional extras add libraries used by specific features:

| Extra | Adds | Used for | Without it |
| --- | --- | --- | --- |
| `structures` | Gemmi | entity-aware mmCIF reading | built-in mmCIF parser |
| `hydration` | MDTraj | DSSP secondary structure in hydration reports | DSSP reported as unavailable |
| `networks` | NetworkX | weighted residue paths (`crevice.connectivity`) | those functions raise an error naming the extra |
| `ml` | scikit-learn | PCA in `crevice.ml.reduce_dimensions` | column truncation, which is not PCA |
| `science` | all of the above | the complete optional feature set | |
| `test` | pytest and the optional libraries the tests use | running the test suite | |

For example: `python -m pip install "crevice[science]"`.

PyMOL, VMD and ChimeraX are separate programs and are not installed by pip.
CREVICE writes scene scripts for them.

## Quick start

Structures can be given as a local file (PDB, mmCIF or GRO), an RCSB
accession or an AlphaFold DB identifier (for example `AF-P02920-F1`); remote
entries are downloaded once and cached with a provenance record.

```bash
# Radius profile of the gramicidin A channel (table also written as 1GRM_profile.csv)
crevice profile 1GRM --axis auto -o 1GRM_profile.json --png 1GRM_profile_radius.png

# 3D cast of the dominant interior cavity of GLUT1, with viewer scenes
crevice cast 4PYP --out-dir 4PYP_cast

# Everything for one structure: profile, cast, residues, cavities, tunnels,
# network, hydration, figures and PyMOL/VMD/ChimeraX scenes
crevice publish 1GRM --out-dir results/1GRM --prefix 1GRM
```

The same profile from Python:

```python
from crevice import annotate_residues, load_structure, pore_profile

frame = load_structure(".crevice/pdb/1GRM.cif")  # cached by the command above, or any local file
profile = pore_profile(frame, axis="auto", samples=81)
contacts = annotate_residues(frame, profile.points)

print(profile.min_radius, profile.bottleneck.nearest_residue)
for contact in contacts[:5]:
    print(contact.residue, contact.role, round(contact.min_distance, 2))
```

Python functions return result objects and do not write files; writers such as
`crevice.write_profile_json` and plots such as `crevice.plot_profile_radius`
are called separately.

## Tools

| Command | Question it answers |
| --- | --- |
| `profile` | How wide is a through-channel along its axis? |
| `cast` | What 3D space does a pocket, cavity or channel occupy? |
| `cavities` | Are there enclosed grid voids? |
| `tunnels` | Which widest grid paths lead from a point to the outside? |
| `residues` | Which residues line a sampled profile? |
| `network` | Which residues touch each other and the channel? |
| `residue-evidence` | Which residues form a measured cavity boundary, and which residues contact them? |
| `hydration` | Where are explicit waters, and how exposed are residues? |
| `trajectory` | How do through-channel profiles and contacts vary across frames? |
| `cavity-trajectory` | How do a cavity's volume, widths and boundary residues change over an MD run? |
| `region-init`, `region-prepare`, `region-trajectory`, `region-compare` | How do named, versioned regions behave over identical frames? |
| `publish` | A complete static-structure bundle in one run |
| `fetch`, `analyze`, `features`, `benchmark`, `static-suite` | Downloads, selected analyses and batch runs |

Run `crevice <command> --help` for every option.

A trajectory example (GRO/XTC or any format MDAnalysis reads):

```bash
crevice cavity-trajectory system.gro run.xtc \
  --reference-volume-dx pocket_volume.dx --geometry-mode reference-region \
  --probe-radius 1.4 --region-margin 2 --spacing 0.2 --axis 0,0,1 \
  --max-frames 1100 --workers 4 --out-dir trajectory-results --prefix system
```

The reference map must be a measured binary map (`*_volume.dx` from `cast` or
`publish`) in the first frame's coordinates. `--max-frames` is a guard against
loading too many frames by accident, not a stride.

## Outputs

Workflow commands write a directory with a `*_manifest.json` that lists every
file, the input provenance and the parser that ran. Typical outputs are:

- CSV tables for every dataset you would analyse (profiles, cast regions,
  residues, cavities, tunnels, networks, per-frame trajectory and hydration
  series, statistics), with units in the column names (`radius_A`,
  `volume_A3`, `time_ps`); JSON only for manifests, provenance, full
  structured records read by other commands and viewer scene definitions;
- a measured binary OpenDX map (`*_volume.dx`), from which every volume, radius
  and residue attribution is computed, and separate smoothed display maps and
  meshes (`--smooth X`, default 0.4 Å, changes the picture only);
- figures (radius profiles, residue contacts with a contact-role legend,
  network summaries and chords, trajectory distributions; 400 dpi in
  publication bundles), clean by default (`--annotate` adds titles and labels);
- PyMOL (`.pml`), VMD (`.vmd`/`.tcl`) and ChimeraX (`.cxc`) scenes that share
  one viewer PDB, with `crevice_render` / `*_render.cxc` for high-resolution
  images. Keep the whole output directory together. No HTML is written.

Unresolved, unavailable and zero are kept distinct: a method that cannot give
a defensible answer records its status and reason instead of a number.

## Scope and interpretation

CREVICE measures geometry. Please read its results with these limits in mind:

- **Geometric output is not biological validation.** A cavity, path, volume,
  radius or residue ranking is a measurement under stated settings (grid
  spacing, probe radii, atom selection, alignment). None of the bundled
  benchmark entries is curated, and CREVICE does not decide which space is
  biologically relevant.
- **Results depend on parameters.** Check sensitivity to grid spacing, probe
  radius and region definition before interpreting a number, and report the
  settings with it. Requested grids are never coarsened silently.
- **Unresolved is not closed.** If no unique two-ended channel is found, the
  profile is reported as unresolved; that does not show the protein is closed.
  Missing water in a static structure is not evidence of a dry pocket.
- **Contacts and correlations are associations.** Residue networks, motion
  correlations and hydration statistics identify candidates for follow-up,
  not causal pathways or residues that control transport.
- **Trajectory uncertainty is approximate.** Mean intervals assume stationary,
  representative sampling and exclude force-field, grid and probe errors.
  Frame-range bands describe fluctuations, not confidence in the mean.
- **Viewer scenes should be inspected.** A generated scene is not checked just
  because CREVICE wrote it; open it in the viewer before relying on it.

Each command's documentation page states its method, assumptions, parameters,
outputs, testing status and known limitations.

## Documentation

The full documentation (installation, quick start, one page per tool, method
notes, CLI and Python API reference) is built from [`docs/`](docs/index.md).
To build it locally:

```bash
python -m pip install . -r docs/requirements.txt
python -m sphinx -b html docs docs/_build/html
```

## Contributing

Contributions, bug reports and careful reports of results that look wrong are
welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) and the
[code of conduct](CODE_OF_CONDUCT.md). To run the tests:

```bash
python -m pip install -e ".[test]"
python -m pytest -q -rs tests
```

## Citing CREVICE

A paper is in preparation:

> Lambden, E., Findlay, H. and Booth, P. J. *CREVICE: Analysis and
> Visualisation of Protein Cavities, Channels, and Residue Interaction
> Networks.* In preparation.

Until it is published, please cite the software using
[CITATION.cff](CITATION.cff) and give the version you used. Please also cite
the methods and tools your analysis relies on (for example HOLE, HOLLOW, CAVER,
3V and MDAnalysis).

## Licence

CREVICE is released under the [MIT licence](LICENSE).

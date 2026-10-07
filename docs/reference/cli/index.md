# Command-line reference

Everything CREVICE does is available from one command, `crevice`, followed by
a sub-command. `python -m crevice` does the same thing.

```bash
crevice profile 1GRM -o 1GRM_profile.json --png 1GRM_profile.png
crevice <command> --help      # every option, with units and defaults
```

Pick a command below to see its arguments and options. Options that many
commands share, such as how a structure is read, which atomic radii are used
and whether figures carry text, are described once on
[Shared options](common-options.md).

## Pores and channels

| Command | What it does |
|---|---|
| [`profile`](profile.md) | Radius profile of a channel that runs right through the protein |
| [`cast`](cast.md) | 3D cast of a pocket, cavity or channel, filled with atom-clear probe spheres |

## Cavities and tunnels

| Command | What it does |
|---|---|
| [`cavities`](cavities.md) | Buried cavities on a clearance grid |
| [`tunnels`](tunnels.md) | Widest routes from a start point to the outside |

## Residues and networks

| Command | What it does |
|---|---|
| [`residues`](residues.md) | Residues lining a channel profile |
| [`residue-evidence`](residue-evidence.md) | Residues that form a measured cast boundary, and the residues they touch |
| [`network`](network.md) | Residue contact networks, with the channel or cavity as a node |
| [`features`](features.md) | Profile and residue descriptors as JSON and CSV |

## Hydration

| Command | What it does |
|---|---|
| [`hydration`](hydration.md) | Water contacts, solvent-accessible area, water density and typed interactions |

## Trajectories and regions

| Command | What it does |
|---|---|
| [`trajectory`](trajectory.md) | Channel profiles and contacts across many frames |
| [`cavity-trajectory`](cavity-trajectory.md) | A cavity's volume, widths and boundary residues in every frame |
| [`region-init`](region-init.md) | Start a named region definition from landmarks |
| [`region-prepare`](region-prepare.md) | Check a region definition and write review scenes |
| [`region-trajectory`](region-trajectory.md) | Run geometry and hydration on a prepared region |
| [`region-compare`](region-compare.md) | Compare two regions on the same frames |

## Visualisation and outputs

| Command | What it does |
|---|---|
| [`publish`](publish.md) | Tables, figures and PyMOL/VMD/ChimeraX scenes for one structure in one run |
| [`analyze`](analyze.md) | Several analyses on one structure, written to one directory |

## Structures and batch runs

| Command | What it does |
|---|---|
| [`fetch`](fetch.md) | Download structures from the RCSB PDB or AlphaFold DB |
| [`fetch-example`](fetch-example.md) | Download a hash-checked example dataset for the tutorials |
| [`benchmark`](benchmark.md) | Quick checks on the bundled benchmark structures |
| [`static-suite`](static-suite.md) | `publish`-style analyses across several structures |

```{toctree}
:hidden:

common-options
profile
cast
cavities
tunnels
residues
residue-evidence
network
features
hydration
trajectory
cavity-trajectory
region-init
region-prepare
region-trajectory
region-compare
publish
analyze
fetch
fetch-example
benchmark
static-suite
```

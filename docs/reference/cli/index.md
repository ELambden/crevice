# Command-line reference

Everything CREVICE does is available from one command, `crevice`, followed by
one of nine commands. `python -m crevice` does the same thing.

```bash
crevice profile 1GRM -o 1GRM_profile.json --png 1GRM_profile.png
crevice cavities --cast 1GRM --out-dir 1GRM_cast
crevice <command> --help          # the command, its modes and the default mode's options
crevice <command> --MODE --help   # the options of one mode, e.g. crevice cavities --cast --help
```

Several commands have **modes**, chosen with one flag (`crevice cavities
--cast`, `crevice residues --evidence`) or, for `crevice region`, an action word
(`crevice region prepare`). Each mode keeps its own options and defaults. The
25 command names of CREVICE 0.3 still work as deprecated aliases that print a
one-line note; see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md).

Options that many commands share, such as how a structure is read, which
atomic radii are used and whether figures carry text, are described once on
[Shared options](common-options.md).

| Command | Modes | What it does |
|---|---|---|
| [`fetch`](fetch.md) | `--example NAME` | Download structures from the RCSB PDB or AlphaFold DB; a hash-checked example dataset |
| [`profile`](profile.md) | `--fingerprint` | Radius profile of a channel that runs right through the protein; the pore fingerprint (lining chemistry, water, geometric conductance estimate) |
| [`cavities`](cavities.md) | `--cast`, `--tunnels` | Buried cavities on a clearance grid; 3D rolling-probe casts of pockets, cavities and channels; widest access tunnels |
| [`paths`](paths.md) | | Curved and branched widest routes from a seed, residues or region to bulk |
| [`residues`](residues.md) | `--network`, `--evidence`, `--features` | Residues lining a channel; contact networks; the residues that form a measured cast boundary; profile and residue descriptors |
| [`report`](report.md) | `--analyze`, `--suite`, `--benchmark` | Tables, figures and PyMOL/VMD/ChimeraX scenes for one structure; selected analyses; benchmark systems |
| [`trajectory`](trajectory.md) | `--cavities`, `--region DEFINITION`, `--water`, `--track` | Channel profiles across frames; a cavity, a named region or water in every frame; persistent cavity identities |
| [`region`](region.md) | `init`, `prepare`, `compare` | Define, check and compare named regions |
| [`view`](view.md) | | A self-contained interactive 3D HTML view of an output directory |

```{toctree}
:hidden:

common-options
fetch
profile
cavities
paths
residues
report
trajectory
region
view
```

# Worked examples

Each example downloads public structures from the RCSB PDB (or AlphaFold DB) by accession, runs
CREVICE commands with stated settings, and shows the printed output, the files
written and the figures. Every command and Python block on these pages was run
with CREVICE {{ version }} and the numbers were copied from that run; the
figures are the files those runs wrote (viewer images were rendered by the
viewer named in the caption). The pages are not executed when the documentation
is built, so the site never needs network access; to rerun them, use
`examples/run_worked_examples.sh` from the source repository (below).

Numbers in an example are observations made with the stated settings, not
validated biological results. Several examples show how a setting (probe,
grid, radius set) changes them.

| Example | Structures | Commands and options shown | Time |
|---|---|---|---|
| [1GRM: pore radius profile and channel cast](1grm-profile-cast.md) | 1GRM | `fetch`, `profile --axis auto` (automatic enclosure probe), `publish`, cast scene in ChimeraX | ~1 min |
| [3UKM: a capped channel with lateral exits](3ukm-lateral-exits.md) | 3UKM, biological assembly 1 | `profile --lateral-exits`, `--assembly`, `--exclude-hetero`, exit legs in the plot and the scenes | ~10 min |
| [1OED: a wide pore with gaps in its wall](1oed-pore-domain.md) | 1OED | the automatic enclosure probe choosing 2.0 Å, `--enclosure-radius`, `residues` | ~2 min |
| [4PYP: interior cavity cast and its boundary residues](4pyp-cavity-residues.md) | 4PYP | an unresolved `profile`, `cast`, `residue-evidence`, `--stick-residues`, evidence scene in ChimeraX | ~2 min |
| [Atomic radii and figure text](custom-radii-annotate.md) | 1GRM | `--radii hole`, a custom radius file, `--annotate`, figure metadata | <1 min |
| [An AlphaFold DB model: fetch by identifier and cast](alphafold-model.md) | AF-P02920-F1 (AlphaFold DB) | `fetch AF-...`, `cast`, the predicted-model flags in the input report | <1 min |
| [A small public trajectory: the 1GRM NMR ensemble](trajectory-nmr-ensemble.md) | 1GRM (5 NMR models) | `trajectory --inspect`, per-frame profiles, distribution and contact tables, Python equivalent | <1 min |

Times are for one CPU core of a desktop computer.

## Rerunning the examples

From a checkout of the source repository, with CREVICE installed:

```bash
bash examples/run_worked_examples.sh results/worked-examples
```

The script downloads the structures into `results/worked-examples/.crevice/pdb`
and runs every command of these pages in its own subdirectory, writing each
command's output to `commands.log`. It needs network access for the
downloads, about 15 minutes and about 50 MB of disk space.

## Synthetic examples in the source repository

The `examples/` directory of the source repository also contains scripts that
use mathematical atom walls with known answers rather than proteins, so they
need no downloads:

| Script | What it shows |
| --- | --- |
| `examples/geometry_demo.py` | a synthetic channel analysed as a short trajectory: profile distribution, void cast and viewer bundle, with a manifest of parameters and file hashes |
| `examples/connectivity_validation.py` | continuous-path connectivity across grid spacings and offsets on synthetic walls (a 48-case sweep with JSON and figures) |

```bash
python examples/geometry_demo.py --output results/geometry-demo
python examples/connectivity_validation.py --output results/connectivity-check
```

```{toctree}
:hidden:

1grm-profile-cast
3ukm-lateral-exits
1oed-pore-domain
4pyp-cavity-residues
custom-radii-annotate
trajectory-nmr-ensemble
alphafold-model
```

# Installation

CREVICE needs Python 3.10 or newer. Continuous integration runs the test suite
on Python 3.10–3.13, including a run with the oldest supported versions of the
required libraries.

```bash
python -m pip install crevice
```

## Required dependencies

Installing CREVICE installs these libraries automatically:

| Library | Minimum | Used for |
|---|---|---|
| NumPy | 1.26 | arrays, grids, statistics |
| SciPy | 1.11 | numerical methods, spatial queries |
| MDAnalysis | 2.7 | GRO/XTC and other MD formats, atom selections, periodic distances |
| Matplotlib | 3.8 | figures |
| scikit-image | 0.22 | marching-cubes meshes for viewer surfaces and residue boundary attribution |

The minimum versions are pinned in `.github/constraints/minimum.txt` for the CI
job that tests the oldest supported versions.

## Optional extras

The extras add only libraries that are **not** already required:

| Extra | Adds | Used for | Without it |
|---|---|---|---|
| `structures` | Gemmi | validated, entity-aware mmCIF reading | built-in mmCIF parser; `--parser gemmi` fails and names the extra |
| `hydration` | MDTraj | DSSP secondary structure in hydration reports | DSSP reported as explicitly unavailable |
| `networks` | NetworkX | weighted residue paths (`network --connectivity-json`, {mod}`crevice.connectivity`) | those functions raise an error naming the extra |
| `ml` | scikit-learn | PCA in {func}`crevice.ml.reduce_dimensions` | column truncation, which is not PCA |
| `science` | all of the above | complete optional feature set | — |
| `test` | pytest, Gemmi, NetworkX, MDTraj (and `tomli` on Python 3.10) | running the regression suite | — |

For example, `python -m pip install "crevice[science]"` installs every optional
library. MDAnalysis and scikit-image are required. scikit-image supplies the
`marching_cubes` surface extraction that default `crevice cast` and
`crevice publish` runs use for viewer cast meshes, boundary-residue surface
attribution and hydration viewer meshes.

## Install from a source checkout

To work on CREVICE itself, or to run the test suite, install an editable copy
from a clone of the source repository:

```bash
git clone https://github.com/ELambden/crevice.git
cd crevice
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .                   # required libraries only
python -m pip install -e ".[science,test]"   # plus every optional library and pytest
```

This installs the `crevice` console script. You can also run
`python -m crevice`, and Python code imports the `crevice` package.

## Check the installation

```bash
crevice --help
python -c "import crevice; print(crevice.__version__)"
```

To run the regression suite, see
[Contributing](../development/contributing.md#running-the-tests).

## Molecular viewers (optional)

CREVICE writes scene scripts for PyMOL (`.pml`), VMD (`.vmd`/`.tcl`) and ChimeraX
(`.cxc`). The viewers are separate programs that pip never installs. A
generated script only counts as checked in a viewer after that viewer has been
run and its output inspected. The `scripts/check_pymol.py`, `check_vmd.py` and
`check_chimerax.py` helpers do those native checks. See
[Outputs](outputs.md#viewer-scripts).

## Releases

Released versions are listed in the [Changelog](../project/changelog.md).
CREVICE follows semantic versioning; before 1.0.0 a minor release may change the
command-line or Python interface or default outputs.

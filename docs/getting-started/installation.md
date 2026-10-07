# Installation

CREVICE runs on Python 3.10 or newer and installs with pip:

```bash
python -m pip install crevice
```

That's all you need for every analysis. Check that it worked:

```bash
crevice --help
python -c "import crevice; print(crevice.__version__)"
```

We recommend installing into a fresh virtual environment or conda environment,
so CREVICE's scientific libraries don't clash with other projects. The test
suite runs on Python 3.10 to 3.13, including a run with the oldest supported
library versions.

## What gets installed

pip brings in these libraries automatically:

| Library | Minimum | Used for |
|---|---|---|
| NumPy | 1.26 | arrays, grids, statistics |
| SciPy | 1.11 | numerical methods, spatial queries |
| MDAnalysis | 2.7 | GRO/XTC and other MD formats, atom selections, periodic distances |
| Matplotlib | 3.8 | figures |
| scikit-image | 0.22 | marching-cubes meshes for viewer surfaces and residue boundary attribution |

## Optional extras

A few features use extra libraries. Each extra adds only what is not already
installed:

| Extra | Adds | Used for | Without it |
|---|---|---|---|
| `structures` | Gemmi | validated, entity-aware mmCIF reading | built-in mmCIF parser; `--parser gemmi` fails and names the extra |
| `hydration` | MDTraj | DSSP secondary structure in hydration reports | DSSP reported as explicitly unavailable |
| `networks` | NetworkX | weighted residue paths (`network --connectivity-json`, {mod}`crevice.connectivity`) | those functions raise an error naming the extra |
| `ml` | scikit-learn | PCA in {func}`crevice.ml.reduce_dimensions` | column truncation, which is not PCA |
| `science` | all of the above | complete optional feature set | — |
| `test` | pytest, Gemmi, NetworkX, MDTraj (and `tomli` on Python 3.10) | running the regression suite | — |

To get everything at once:

```bash
python -m pip install "crevice[science]"
```

## Installing from source

If you want the latest development version, or to change CREVICE yourself,
install an editable copy from a clone of the repository:

```bash
git clone https://github.com/ELambden/crevice.git
cd crevice
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .                   # required libraries only
python -m pip install -e ".[science,test]"   # plus every optional library and pytest
```

Either way you get the `crevice` command (`python -m crevice` works too) and
the `crevice` Python package. `CONTRIBUTING.md` in the repository explains how
to run the tests.

## Molecular viewers (optional)

CREVICE writes ready-to-open scenes for PyMOL (`.pml`), VMD (`.vmd`/`.tcl`) and
ChimeraX (`.cxc`). The viewers themselves are separate programs, so install
whichever you like to use from its own website; pip never installs them. See
[Viewer scripts](outputs.md#viewer-scripts) for what each scene contains.

## Upgrading

```bash
python -m pip install --upgrade crevice
```

The [Changelog](../project/changelog.md) lists what changed in each release.
Until version 1.0, a minor release may change options or default outputs, so
note the version you used with your results.

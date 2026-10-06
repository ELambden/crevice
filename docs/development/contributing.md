# Contributing and developer guide

CREVICE welcomes contributions to its methods, tests and documentation. It also
welcomes careful reports of results that look wrong.

The repository-level [CONTRIBUTING.md](../../CONTRIBUTING.md) has the version
policy, pull-request checklist and dependency rules.
[CODE_OF_CONDUCT.md](../../CODE_OF_CONDUCT.md) and
[SECURITY.md](../../SECURITY.md) are also at the repository root.

## Source layout

| Path | Contents |
|---|---|
| `src/crevice/` | the package |
| `tests/` | regression, scientific-contract and optional native-viewer tests |
| `examples/` | synthetic, dependency-light examples (mathematical walls, not proteins) |
| `scripts/` | native-viewer checkers (`check_pymol.py`, `check_vmd.py`, `check_chimerax.py`, `check_analysis_viewers.py`, `check_region_viewers.py`) and the wheel install check (`check_wheel.py`) |
| `docs/` | this site; the upper-case `*.md` method notes are built in place |

(running-the-tests)=
## Running the tests

```bash
python -m pip install -e ".[test]"
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q -rs tests
```

Use `-rs` so that every skipped test is listed with its reason. Some tests skip
when an optional prerequisite is missing: the native PyMOL, VMD and ChimeraX
checks need those programs (PyMOL must be importable from the test
interpreter; ChimeraX is looked up under `$CREVICE_WORKSPACE_DIR/tools/bin` or
`.crevice/tools/bin` and also needs `xvfb-run`); the real-structure registry
test needs `CREVICE_RUN_PDB_INTEGRATION=1` and either `CREVICE_PDB_DIR` (a local
structure cache) or `CREVICE_FETCH_PDB=1` (network access); the Tcl parse check
of VMD scripts needs `tkinter`; and one trajectory-reader test skips when the
installed MDAnalysis writer records no time step. When you report a run, give
the exact command, the counts and every skip reason. Do not report only the
pass count.

## Native viewer checks

The `scripts/check_*.py` helpers run a generated scene in the real viewer
without a display, save a PNG and a log, and write a machine-readable report of
what they verified:

```bash
python scripts/check_pymol.py results/protein/protein_volume.pml --output results/pymol-check
python scripts/check_vmd.py results/protein/protein_volume.tcl --output results/vmd-check
python scripts/check_chimerax.py results/protein/protein_volume.cxc --output results/chimerax-check --xvfb
```

A change to scene generation should be checked this way in every affected
viewer, and the pull request should say which viewers and versions were run.

## Scientific standards for changes

- **Keep validation levels apart.** Unit and synthetic tests, observations on
  real inputs, native viewer checks, and biological or functional validation are
  different kinds of evidence (see
  [Validation levels](../tools/index.md#validation-levels)). A generated viewer
  script is not "tested in PyMOL" unless PyMOL actually ran and its output was
  checked.
- **Don't fill gaps.** Unresolved methods must stay unresolved, and missing
  observations must stay missing, not zero.
- **Keep numbers stable.** A change that alters numerical results needs its own
  reviewed change, with before-and-after values.
- **Be exact about grids.** Never coarsen a requested grid silently. Record the
  effective settings and the input provenance.

## Documentation changes

See [Documentation](documentation.md) for how to build these pages.

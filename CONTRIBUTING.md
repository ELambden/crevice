# Contributing to CREVICE

Thank you for helping. CREVICE is early-stage research software, so the most
useful contributions are careful ones: small changes, exact commands and an
honest statement of what was and was not checked.

Please follow the [Code of Conduct](CODE_OF_CONDUCT.md). Questions and ideas
are welcome in the repository's GitHub Discussions; use issues for reproducible
bugs and concrete feature requests. Report security problems privately as
described in [SECURITY.md](SECURITY.md).

## Development install

CREVICE requires Python 3.10 or newer. NumPy, SciPy, MDAnalysis, Matplotlib and
scikit-image are installed automatically as required dependencies.

```bash
git clone https://github.com/ELambden/crevice.git
cd crevice
python -m venv .venv
. .venv/bin/activate
python -m pip install -e ".[test]"
```

The `test` extra adds pytest and the optional libraries the regression suite
exercises (Gemmi, NetworkX, MDTraj). Native PyMOL, VMD and
ChimeraX are never installed by pip; tests that need them skip when absent.

## Running the tests

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q -rs tests
```

`-rs` prints every skip reason. Record the pass/skip/fail counts **and** the
skip reasons in your pull request: optional tests skip when their
prerequisites are missing, and a skip is not a pass. The optional prerequisites
are:

- native viewers: PyMOL importable from the test interpreter, `vmd` on `PATH`,
  and ChimeraX under `$CREVICE_WORKSPACE_DIR/tools/bin/chimerax` or
  `.crevice/tools/bin/chimerax` together with `xvfb-run`;
- `tkinter`, for the Tcl parse check of generated VMD scripts;
- the real-structure registry test: `CREVICE_RUN_PDB_INTEGRATION=1` plus either
  `CREVICE_PDB_DIR` (a local structure cache) or `CREVICE_FETCH_PDB=1`
  (network access).

Native viewer scenes can also be checked directly:

```bash
python scripts/check_pymol.py results/protein/protein_volume.pml --output results/pymol-check
python scripts/check_vmd.py results/protein/protein_volume.tcl --output results/vmd-check
python scripts/check_chimerax.py results/protein/protein_volume.cxc --output results/chimerax-check --xvfb
```

To check packaging locally:

```bash
python -m pip install build twine
python -m build --outdir dist            # sdist, then a wheel built from it
python -m twine check --strict dist/*
python scripts/check_wheel.py dist/crevice-*.whl --report dist/wheel-check.json
```

## Dependency policy

- Required: NumPy, SciPy, MDAnalysis, Matplotlib and scikit-image. Code may
  import them unconditionally, and tests must not `importorskip` them.
- Optional: everything else (Gemmi, NetworkX, MDTraj, scikit-learn, native
  viewers). Import them lazily, and when one is missing
  either fall back to a documented, tested path or fail with a message that
  names the install extra (for example `pip install 'crevice[structures]'`).
- Optional extras in `pyproject.toml` must not repeat required libraries.
  Lower bounds of the required libraries are pinned in
  `.github/constraints/minimum.txt` and tested in CI; keep both in sync
  (`tests/test_packaging.py` enforces this).

## Branches and pull requests

- `main` is the integration branch; do not push to it directly.
- Work on a short-lived branch named for its purpose, e.g. `fix/gro-units`,
  `feature/smooth-surfaces` or `docs/install`.
- Keep a pull request to one topic. Update tests, user documentation and the
  `[Unreleased]` section of `CHANGELOG.md` in the same pull request.
- Preserve existing behaviour and numerical results unless the change is meant
  to alter them; if it is, say so and show before/after values.
- CI must pass. Describe any check you could not run.

## Reporting evidence for a change

Describe the evidence for a change in the pull request itself:

1. the exact commands, CREVICE version or commit and environment
   (`pip freeze`);
2. inputs with their source (for example an RCSB accession) and SHA-256 hashes;
3. results, including test counts and every skip reason;
4. which outputs were inspected and how (for example, which viewer was run);
5. limitations and anything left unfinished.

Do not commit large inputs such as trajectories, or generated result
directories; attach small outputs or summaries to the pull request instead.

## Scientific honesty conventions

These rules apply to code comments, documentation, changelog entries, issues
and pull requests.

- **Say what kind of evidence you have.** Distinguish software or synthetic
  tests, real-system observations and functional or biological inference. A
  passing test shows the software behaves as specified, not that a biological
  conclusion is correct.
- **Benchmarks are not curated.** Pore axes and regions for the bundled
  structures are assumed, not determined. Never present a cavity count,
  volume, radius or residue ranking as a validated biological result.
- **Viewer scripts.** Do not describe a generated PyMOL, VMD or ChimeraX
  script as natively tested unless that viewer was actually run and its
  output checked.
- **Report uncertainty and sensitivity.** State grid spacing, probe radius,
  selections and frame ranges. When a result depends strongly on a parameter,
  say so rather than choosing a favourable value.
- **Mark unfinished work as unfinished.** Do not describe planned features as
  available.

## Versioning and releases

CREVICE follows [Semantic Versioning 2.0.0](https://semver.org/).

- The version has a single source: `__version__` in `src/crevice/__init__.py`.
  Hatchling reads it for the package metadata (`dynamic = ["version"]`); do
  not add a version anywhere else. `CITATION.cff` repeats it and must be
  updated at release time.
- Before 1.0.0, a MINOR bump (0.**y**.0) may change the public CLI or Python
  API or default outputs; a PATCH bump (0.y.**z**) is for fixes that do not
  change documented behaviour. Any change to numerical results of an existing
  analysis must be called out in the changelog regardless of the bump.
- Release steps are in `docs/development/releasing.md`: bump the version,
  update the changelog, push to `main`, then publish a GitHub Release tagged
  `vX.Y.Z`; the release workflow tests, builds and uploads it to PyPI.
- Pre-releases use SemVer suffixes in the tag (e.g. `v0.2.0-rc.1`) and the
  matching PEP 440 form in `__version__` (`0.2.0rc1`).

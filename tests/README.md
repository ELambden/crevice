# CREVICE tests

Run the suite from the repository root after installing the test extra:

```bash
python -m pip install -e ".[test]"
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q -rs tests
```

`-rs` lists every skipped test with its reason. Report the pass, skip and fail
counts together with every skip reason; a skip is not a pass.

The tests cover input identities and provenance, geometry and connectivity,
channel profiles and their uncertainty, cavity trajectories and named regions,
residue contacts and evidence, hydration, water density and membership, typed
chemistry, packaging policy and viewer exports. `test_scientific_contracts.py`
holds independent contract tests with known geometric answers.

These are software and synthetic checks. They do not validate any biological
conclusion, and a generated viewer script is not checked in that viewer unless
a native test actually ran it.

Optional tests skip when their prerequisites are absent:

| Tests | Prerequisite |
| --- | --- |
| native PyMOL checks | `pymol` importable from the test interpreter |
| native VMD checks | `vmd` on `PATH` |
| native ChimeraX checks | ChimeraX at `$CREVICE_WORKSPACE_DIR/tools/bin/chimerax` or `.crevice/tools/bin/chimerax`, and `xvfb-run` |
| VMD script Tcl parse check | `tkinter` |
| real-structure registry load | `CREVICE_RUN_PDB_INTEGRATION=1` and either `CREVICE_PDB_DIR` or `CREVICE_FETCH_PDB=1` |
| trajectory time-step case | an MDAnalysis writer that records a time step |

# Documentation

## Layout

The Sphinx source directory is `docs/` itself, so the upper-case method notes
(`docs/HYDRATION_METHODS.md` and the others) are built in place by
[MyST-Parser](https://myst-parser.readthedocs.io/). Site pages live in
subdirectories:

| Directory | Contents |
|---|---|
| `getting-started/` | installation, quick start, CLI or Python, inputs, outputs |
| `tools/` | one page per command, each with a generated option list |
| `methods/` | index of the in-place method notes |
| `examples/` | worked examples; add one page per example and list it in the index page's toctree. Figures are pre-rendered files in `examples/images/`; the pages are not executed during the build |
| `reference/cli.md` | complete CLI reference generated from `crevice.cli.build_parser()` by sphinx-argparse |
| `api/` | grouped Python API index and a module index (autodoc/autosummary) |
| `development/`, `project/` | contributor pages, changelog, citation, community |
| `_ext/crevice_docs.py` | local extension: repository-file links, the API-coverage check, resolution of bare type names in docstrings, and tuple defaults in type fields |

`conf.py` builds only the pages matched by its `include_patterns`; add a new
directory or method note there as well as to a toctree.

## Building locally

Use a separate environment so that the documentation tools do not change the
test environment:

```bash
python -m venv .venv-docs
.venv-docs/bin/python -m pip install . -r docs/requirements.txt
.venv-docs/bin/python -m sphinx -W --keep-going -b html docs docs/_build/html
.venv-docs/bin/python -m sphinx -n -W --keep-going -b html docs docs/_build/html-nitpicky
.venv-docs/bin/python -m sphinx -b linkcheck docs docs/_build/linkcheck
```

Delete `docs/api/generated/` first when modules or docstrings have been
removed or renamed: autosummary does not delete stale stub pages.

Installing `.` brings CREVICE's required scientific libraries, which autodoc
needs to import every module. `docs/requirements.txt` lists only the
documentation tools.

`-W` turns warnings into errors. Intersphinx needs network access to fetch the
Python, NumPy, SciPy, Matplotlib and MDAnalysis inventories. Without it, those
cross-links fail and so does the `-W` build. The link check also needs network
access, and some publisher sites may refuse automated requests.

API pages import `crevice` from `src/`, so the documentation always describes the
checked-out source. Optional heavy dependencies (Gemmi, MDTraj, scikit-image,
NetworkX, ...) are mocked for autodoc.

## Build checks

- **API coverage.** `api/modules.rst` is an explicit list of submodules. The
  build fails if a `crevice` submodule is missing from that list, or if the list
  names a module that does not exist.
- **Nitpicky references.** The `-n` build checks every cross-reference,
  including the types in NumPy-style docstrings. `conf.py` preprocesses those
  types (`napoleon_preprocess_types`, with aliases such as `array_like` and
  `sequence`), and lists in `nitpick_ignore_regex` the words that describe
  rather than name a type (output kinds such as `table` or `figure`, array
  shapes) and a few unparsable or uninventoried targets. A new unresolved name
  is an error under `-n -W`; add a real target rather than an ignore where
  possible.
- **Docstring coverage.** `tests/test_docstring_coverage.py` fails when a
  public module, class or function of `crevice` lacks a docstring, or when a
  parameter or dataclass field is not mentioned in it.
- **Links outside `docs/`.** A page may link to a file elsewhere in the
  repository (for example `CHANGELOG.md`). The link must point to a file that
  exists; otherwise the build fails. Such links become repository links when
  `CREVICE_REPOSITORY_URL` is set, and labelled repository paths otherwise. The
  files are never copied into the build.

## Keeping pages in step with the interface

The option lists on the tool pages and on the CLI reference are generated, so
they update automatically. The descriptions are written by hand and must be
updated when behaviour changes: the tool page's method, parameter and output
sections, the relevant method note, and `getting-started/cli-vs-python.md` when
a Python signature or default changes.

## Continuous integration and Read the Docs

The `documentation` job in `.github/workflows/tests.yml` installs `.` and
`docs/requirements.txt` on Python 3.12, and runs the `-W` HTML build as a
required step. It then runs linkcheck with `continue-on-error`, because
publisher sites often refuse automated requests, and uploads the HTML and the
linkcheck report as an artifact.

`.readthedocs.yaml` at the repository root installs the same two things and
builds with `fail_on_warning: true`. Repository links point to
`https://github.com/ELambden/crevice` by default; set the
`CREVICE_REPOSITORY_URL` environment variable to point them elsewhere (an empty
value turns them into labelled paths). The worked examples are not executed on
Read the Docs, so the build needs no data downloads; only the intersphinx
inventories are fetched.

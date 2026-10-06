# Releasing CREVICE

Releases are published from GitHub. Publishing a GitHub Release whose tag is
`vX.Y.Z` runs the `Release to PyPI` workflow
(`.github/workflows/release.yml`), which:

1. checks that the tag matches `__version__` in `src/crevice/__init__.py`;
2. installs the package with its test dependencies and runs the test suite;
3. builds the sdist and wheel and checks their metadata with `twine`;
4. uploads both files to PyPI using
   [Trusted Publishing](https://docs.pypi.org/trusted-publishers/), so no
   API token is stored in the repository.

If any step fails, nothing is uploaded. A PyPI version can never be replaced,
so every fix needs a new version number.

Pushing commits to `main` does **not** publish anything. It runs the test,
build and documentation checks (`.github/workflows/tests.yml`) and, once the
Read the Docs project is connected, rebuilds the development ("latest")
documentation.

## One-time setup

These steps need the maintainer's accounts and are done once.

1. **PyPI trusted publisher.** Sign in to <https://pypi.org>, open
   *Your account → Publishing*, and add a *pending publisher* with:
   - PyPI project name: `crevice`
   - Owner: `ELambden`
   - Repository name: `crevice`
   - Workflow name: `release.yml`
   - Environment name: `pypi`

   The first successful release creates the `crevice` project and binds it to
   this repository.
2. **TestPyPI trusted publisher (optional dry runs).** Repeat step 1 on
   <https://test.pypi.org> with environment name `testpypi`.
3. **GitHub environments.** In the repository's *Settings → Environments*,
   create `pypi` and `testpypi`. Adding yourself as a required reviewer on
   `pypi` makes every upload wait for a one-click approval.
4. **Read the Docs.** Import the repository at <https://readthedocs.org>. The
   `.readthedocs.yaml` file configures the build; enable builds for tags so
   that each release gets versioned documentation.
5. **Zenodo (archive DOI).** Enable the repository at
   <https://zenodo.org/account/settings/github/>. Each GitHub Release is then
   archived with its own DOI, and a concept DOI covers all versions.
   `CITATION.cff` supplies the metadata.

## Making a release

1. Update `__version__` in `src/crevice/__init__.py` and `version` in
   `CITATION.cff`, and set `date-released` in `CITATION.cff`.
2. In `CHANGELOG.md`, move the `[Unreleased]` entries under a new
   `[X.Y.Z] - YYYY-MM-DD` heading.
3. Commit, push to `main`, and wait for the checks to pass.
4. On GitHub choose *Releases → Draft a new release*, create the tag `vX.Y.Z`
   on that commit, paste the changelog entries as release notes and publish.
5. Watch the *Release to PyPI* workflow. When it finishes, check
   <https://pypi.org/project/crevice/> and install the new version in a fresh
   environment: `pip install --upgrade crevice`.

For a dry run, open *Actions → Release to PyPI → Run workflow*. This builds
and tests the current commit and uploads it to TestPyPI; a version already on
TestPyPI is skipped. Install it with
`pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ crevice`.

Pre-releases use a SemVer suffix in the tag (`v0.2.0-rc.1`) and the matching
PEP 440 form in `__version__` (`0.2.0rc1`); mark the GitHub Release as a
pre-release. PyPI installs pre-releases only when asked (`pip install --pre`).

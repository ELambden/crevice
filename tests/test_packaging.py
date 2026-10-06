"""Packaging and dependency-policy contracts.

These are software/metadata checks only. They pin the policy that NumPy, SciPy,
MDAnalysis and Matplotlib plus scikit-image (needed by
the default cast/publish viewer bundles) are required runtime dependencies,
that optional extras do not repeat them, and that the package version has a
single source. They make no scientific claim.
"""
from __future__ import annotations

import importlib
import importlib.metadata
from pathlib import Path
import re

import pytest

import crevice

try:  # Python >= 3.11
    import tomllib
except ModuleNotFoundError:  # Python 3.10: the test extra installs tomli
    try:
        import tomli as tomllib
    except ModuleNotFoundError:  # pragma: no cover
        tomllib = None

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = {"numpy": "numpy", "scipy": "scipy", "mdanalysis": "MDAnalysis",
            "matplotlib": "matplotlib", "scikit-image": "skimage"}
DISTRIBUTION = {"scikit-image": "scikit-image"}
needs_toml = pytest.mark.skipif(tomllib is None, reason="needs tomllib (Python >= 3.11) or tomli")


def _name(requirement: str) -> str:
    return re.split(r"[\s;<>=!~\[(]", requirement.strip(), maxsplit=1)[0].lower()


def _floor(requirement: str) -> str:
    match = re.search(r">=\s*([0-9][0-9.]*)", requirement)
    assert match, f"{requirement!r} has no >= lower bound"
    return match.group(1)


def _numeric(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in re.findall(r"\d+", version)[:3])


def _pyproject() -> dict:
    return tomllib.loads((ROOT / "pyproject.toml").read_text())


@pytest.mark.parametrize("module", sorted(REQUIRED.values()))
def test_required_runtime_library_imports(module):
    assert importlib.import_module(module).__name__ == module


@needs_toml
def test_required_dependencies_are_declared_with_lower_bounds():
    dependencies = _pyproject()["project"]["dependencies"]
    assert {_name(r) for r in dependencies} == set(REQUIRED)
    for requirement in dependencies:
        assert ";" not in requirement, "required dependencies must be unconditional"
        _floor(requirement)


@needs_toml
def test_installed_required_libraries_meet_declared_floors():
    for requirement in _pyproject()["project"]["dependencies"]:
        name = _name(requirement)
        installed = importlib.metadata.version(DISTRIBUTION.get(name, REQUIRED[name]))
        assert _numeric(installed) >= _numeric(_floor(requirement)), (requirement, installed)


@needs_toml
def test_minimum_constraints_pin_the_declared_floors():
    pins = {}
    for line in (ROOT / ".github/constraints/minimum.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            name, version = line.split("==")
            pins[name.strip().lower()] = version.strip()
    assert set(pins) == set(REQUIRED)
    for requirement in _pyproject()["project"]["dependencies"]:
        floor, pinned = _numeric(_floor(requirement)), _numeric(pins[_name(requirement)])
        assert pinned[:len(floor)] == floor, (requirement, pins[_name(requirement)])


@needs_toml
def test_optional_extras_do_not_repeat_required_dependencies():
    extras = _pyproject()["project"]["optional-dependencies"]
    for extra, requirements in extras.items():
        repeated = {_name(r) for r in requirements} & set(REQUIRED)
        assert not repeated, f"extra {extra!r} repeats required {sorted(repeated)}"
    # Optional libraries stay optional: every one used by CREVICE is in `science`.
    science = {_name(r) for r in extras["science"]}
    assert {"gemmi", "networkx", "mdtraj", "scikit-learn"} <= science


@needs_toml
def test_version_is_single_sourced_from_package():
    project = _pyproject()
    assert "version" not in project["project"]
    assert "version" in project["project"]["dynamic"]
    source = project["tool"]["hatch"]["version"]["path"]
    text = (ROOT / source).read_text()
    assert re.findall(r'^__version__ = "([^"]+)"$', text, re.M) == [crevice.__version__]


def test_version_is_semantic():
    # MAJOR.MINOR.PATCH, with PEP 440 pre/post/dev suffixes for pre-releases
    # (tag v0.2.0-rc.1 <-> __version__ 0.2.0rc1; see CONTRIBUTING.md).
    semver = r"(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)((a|b|rc)\d+)?(\.post\d+)?(\.dev\d+)?"
    assert re.fullmatch(semver, crevice.__version__)


def test_ci_does_not_hard_code_the_wheel_version():
    workflow = (ROOT / ".github/workflows/tests.yml").read_text()
    assert f"crevice-{crevice.__version__}" not in workflow


def test_missing_required_library_hint_asks_for_reinstall_not_an_extra(tmp_path, monkeypatch):
    # NumPy, SciPy and Matplotlib are required, so a missing one means a broken
    # install; the recorded hint must not point to an optional extra.
    import argparse
    import json

    from crevice import hydration_workflow

    monkeypatch.setattr(hydration_workflow, "hydration_dependencies", lambda: ["scipy"])
    args = argparse.Namespace(skip_hydration=False, prefix="demo", out_dir=str(tmp_path))
    files = {}
    hydration_workflow.append_static_hydration(args, None, tmp_path / "demo.pdb", files)
    record = json.loads((tmp_path / "demo_hydration.json").read_text())
    assert record["status"] == "unavailable_missing_required_dependencies"
    assert record["missing"] == ["scipy"]
    assert "reinstall" in record["install"] and "[" not in record["install"]


def test_install_hints_only_name_existing_optional_extras():
    extras = set(_pyproject()["project"]["optional-dependencies"]) if tomllib else None
    named = set()
    for path in (ROOT / "src/crevice").glob("*.py"):
        named |= set(re.findall(r"crevice\[([A-Za-z0-9_-]+)\]", path.read_text()))
    assert named, "expected at least one install-extra hint"
    assert "science" not in named, "required libraries must not be attributed to an extra"
    if extras is not None:
        assert named <= extras - {"science", "test"}

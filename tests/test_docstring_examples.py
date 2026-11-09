"""The ``Examples`` sections of the package docstrings run and print what they show.

Software checks: each module's doctests run with :mod:`doctest`. A module whose
optional dependency is not installed is skipped.
"""
from __future__ import annotations

import doctest
import importlib
import pkgutil

import pytest

import crevice

# __main__ runs the command line on import.
MODULES = sorted(info.name for info in pkgutil.walk_packages(crevice.__path__, "crevice.")
                 if not info.name.endswith(".__main__"))


@pytest.mark.parametrize("name", MODULES)
def test_docstring_examples(name):
    try:
        module = importlib.import_module(name)
    except ImportError as error:
        pytest.skip(f"{name} needs an optional dependency: {error}")
    result = doctest.testmod(module, optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE, verbose=False)
    assert result.failed == 0, f"{result.failed} of {result.attempted} examples failed in {name}"

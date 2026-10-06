"""Every public module, class, function and method of ``crevice`` has a substantive docstring.

The rule is checked on the source (AST), so it covers code that is never
imported by the tests. An object is *public* when its name has no leading
underscore and it is defined at module level or inside a public class
(nested functions are implementation details). A docstring is *substantive*
when:

* module: at least 40 characters;
* function or method: a summary line of at least 10 characters, and every
  named parameter (except ``self``/``cls``, including ``*args``/``**kwargs``
  names) appears in the docstring as a word, so each argument and option is
  explained;
* class: a summary line of at least 10 characters, and every annotated class
  attribute (dataclass field) and every ``__init__`` parameter appears in the
  class docstring (``__init__`` parameters may instead be in the
  ``__init__`` docstring).

Every ``scripts/*.py`` file present must also have a module docstring.
"""
import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "crevice"

# No exemptions: every public object in src/crevice must meet the rule.
ALLOWLIST: dict[str, set[str]] = {}


def _params(fn):
    a = fn.args
    names = [x.arg for x in a.posonlyargs + a.args + a.kwonlyargs]
    names += [x.arg for x in (a.vararg, a.kwarg) if x is not None]
    return [n for n in names if n not in ("self", "cls")]


def _mentions(doc, name):
    return re.search(rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])", doc) is not None


def _summary_ok(doc):
    lines = [line.strip() for line in (doc or "").strip().splitlines() if line.strip()]
    return bool(lines) and len(lines[0]) >= 10


def _function_problem(node, extra=""):
    doc = ast.get_docstring(node) or ""
    if not _summary_ok(doc):
        return "no docstring" if not doc else "summary line shorter than 10 characters"
    missing = [p for p in _params(node) if not _mentions(doc + "\n" + extra, p)]
    return f"parameters not described: {', '.join(missing)}" if missing else None


def _class_problems(node, qualname):
    problems = {}
    doc = ast.get_docstring(node) or ""
    if not _summary_ok(doc):
        problems[qualname] = "no docstring" if not doc else "summary line shorter than 10 characters"
    else:
        fields = [s.target.id for s in node.body if isinstance(s, ast.AnnAssign)
                  and isinstance(s.target, ast.Name) and not s.target.id.startswith("_")]
        init = next((s for s in node.body if isinstance(s, ast.FunctionDef) and s.name == "__init__"), None)
        init_doc = (ast.get_docstring(init) or "") if init else ""
        missing = [f for f in fields if not _mentions(doc, f)]
        missing += [p for p in (_params(init) if init else []) if p not in fields and not _mentions(doc + "\n" + init_doc, p)]
        if missing:
            problems[qualname] = f"fields/parameters not described: {', '.join(missing)}"
    for member in node.body:
        if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)) and not member.name.startswith("_"):
            problem = _function_problem(member)
            if problem:
                problems[f"{qualname}.{member.name}"] = problem
        elif isinstance(member, ast.ClassDef) and not member.name.startswith("_"):
            problems.update(_class_problems(member, f"{qualname}.{member.name}"))
    return problems


def docstring_problems(path):
    """Map each public object of a module that lacks a substantive docstring to the reason."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    problems = {}
    if len((ast.get_docstring(tree) or "").strip()) < 40:
        problems["<module>"] = "module docstring missing or shorter than 40 characters"
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
            problem = _function_problem(node)
            if problem:
                problems[node.name] = problem
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            problems.update(_class_problems(node, node.name))
    return problems


MODULES = sorted(PACKAGE.glob("*.py"))


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.stem)
def test_public_api_has_substantive_docstrings(path):
    problems = {name: why for name, why in docstring_problems(path).items()
                if name not in ALLOWLIST.get(path.stem, set())}
    assert not problems, "\n".join(f"crevice.{path.stem}.{name}: {why}" for name, why in sorted(problems.items()))


def test_allowlist_is_current():
    stale = [f"{module}.{name}" for module, names in ALLOWLIST.items()
             for name in names if name not in docstring_problems(PACKAGE / f"{module}.py")]
    assert not stale, "documented now; remove from ALLOWLIST: " + ", ".join(sorted(stale))


def test_rule_detects_missing_parameter_and_field_descriptions(tmp_path):
    sample = tmp_path / "sample.py"
    sample.write_text('"""A module docstring that is long enough for the rule."""\n'
                      'def f(a, b):\n    """Do something with a only."""\n'
                      'class C:\n    """A class with one documented field: x."""\n    x: int\n    y: int\n'
                      'def _private(z):\n    pass\n')
    assert docstring_problems(sample) == {"f": "parameters not described: b",
                                          "C": "fields/parameters not described: y"}


@pytest.mark.parametrize("path", sorted((ROOT / "scripts").glob("*.py")), ids=lambda p: p.stem)
def test_scripts_have_module_docstrings(path):
    assert len((ast.get_docstring(ast.parse(path.read_text(encoding="utf-8"))) or "").strip()) >= 40

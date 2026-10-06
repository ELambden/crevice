"""Local Sphinx helpers for the CREVICE documentation site.

1. Links from pages to repository files outside ``docs/`` (for example
   ``CHANGELOG.md`` or ``LICENSE``) cannot become site links.  This extension
   rewrites them into repository links when ``crevice_repository_url`` is
   configured, and into labelled repository paths otherwise, so such files are
   never copied into ``_downloads``.  A link to a file that does not exist in
   the repository is reported as a warning (an error under ``-W``).

2. Pages refer to the public namespace (``{func}`crevice.pore_profile```), but
   autodoc documents each object in its defining module
   (``crevice.channels.pore_profile``).  A ``missing-reference`` handler maps
   top-level ``crevice.NAME`` targets that are re-exported from a submodule to that
   canonical location, so both spellings link to the same entry.

3. Docstring type fields and signatures use bare CREVICE names (``Coord``,
   ``StructureFrame``); a second handler resolves an unambiguous bare name to
   its documented definition.

4. Tuple and list defaults in NumPy type fields (``default (0.1, 0.9)``) are
   moved into the description before napoleon parses the type (see
   ``_move_tuple_defaults``).

5. The API module index (``api/modules.rst``) is an explicit list.  At build time
   every importable ``crevice`` submodule is compared with that list and a warning
   (an error under ``-W``) is emitted for any undocumented or stale module, so a
   new module cannot silently drop out of the API reference.
"""

from __future__ import annotations

import os
import pkgutil
import re
from pathlib import Path, PurePosixPath

from docutils import nodes
from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.util import logging

logger = logging.getLogger(__name__)

#: Submodules deliberately absent from the API module index.
API_EXCLUDED_MODULES = {"__main__"}


def _outside_target(app: Sphinx, node: nodes.Element, docname: str) -> str | None:
    """Return the repository-relative path if ``node`` points outside ``docs/``.

    Paths are normalised lexically (symlinks are not followed), so the result is
    the same in any checkout.
    """
    if node.get("reftype") != "myst":
        return None
    srcdir = Path(app.srcdir).resolve()
    target = str(node.get("reftarget", ""))
    if not target:
        return None
    target = target.split("#", 1)[0]
    if node.get("refdomain") == "doc":
        resolved = Path(os.path.normpath(srcdir / target))
    else:
        resolved = Path(os.path.normpath(srcdir / PurePosixPath(docname).parent / target))
    try:
        resolved.relative_to(srcdir)
        return None  # inside the site; let MyST/Sphinx resolve it normally
    except ValueError:
        pass
    repo_root = srcdir.parent
    try:
        rel = resolved.relative_to(repo_root).as_posix()
    except ValueError:
        return target
    if node.get("refdomain") == "doc" and not rel.endswith((".md", ".html")):
        rel += ".md"
    return rel


def _rewrite_outside_links(app: Sphinx, doctree: nodes.document) -> None:
    docname = app.env.docname
    repo_url = (app.config.crevice_repository_url or "").rstrip("/")
    repo_root = Path(app.srcdir).resolve().parent
    for node in list(doctree.findall(lambda n: isinstance(n, (addnodes.pending_xref, addnodes.download_reference)))):
        rel = _outside_target(app, node, docname)
        if rel is None:
            continue
        text = node.astext()
        if not (repo_root / rel).exists():
            logger.warning(f"link to {rel!r}, which is not a file in the repository",
                           location=node, type="crevice", subtype="link")
        if repo_url:
            replacement: nodes.Node = nodes.reference(
                "", text, refuri=f"{repo_url}/blob/main/{rel}", classes=["crevice-repo-link"]
            )
        else:
            label = "repository file"
            replacement = nodes.inline(
                "",
                "",
                nodes.Text(text + " "),
                nodes.inline("", f"[{label}: ", classes=["crevice-path-label"]),
                nodes.literal("", rel),
                nodes.inline("", "]", classes=["crevice-path-label"]),
                classes=["crevice-repo-path"],
            )
        node.replace_self(replacement)


def _check_api_module_index(app: Sphinx) -> None:
    import crevice

    index = Path(app.srcdir) / "api" / "modules.rst"
    if not index.exists():
        logger.warning("API module index api/modules.rst is missing", type="crevice", subtype="api")
        return
    listed = set(re.findall(r"^\s*crevice\.([A-Za-z_][A-Za-z0-9_]*)\s*$", index.read_text(), flags=re.M))
    actual = {m.name for m in pkgutil.iter_modules(crevice.__path__)} - API_EXCLUDED_MODULES
    for name in sorted(actual - listed):
        logger.warning(f"crevice.{name} is not listed in api/modules.rst", type="crevice", subtype="api")
    for name in sorted(listed - actual):
        logger.warning(f"api/modules.rst lists crevice.{name}, which does not exist", type="crevice", subtype="api")


def _resolve_public_alias(app: Sphinx, env, node: nodes.Element, contnode: nodes.Element):
    """Resolve ``crevice.NAME`` to the submodule where NAME is documented."""
    if node.get("refdomain") != "py":
        return None
    target = node.get("reftarget", "")
    if not target.startswith("crevice.") or target.count(".") != 1:
        return None
    import crevice

    name = target.split(".", 1)[1]
    if name not in getattr(crevice, "__all__", ()):
        return None
    module = getattr(getattr(crevice, name), "__module__", None)
    if not module or not module.startswith("crevice.") or module == "crevice":
        return None
    domain = env.get_domain("py")
    alias = node.deepcopy()
    alias["reftarget"] = f"{module}.{name}"
    return domain.resolve_xref(env, node["refdoc"], app.builder, node["reftype"],
                               alias["reftarget"], alias, contnode)


_SHORT_NAMES: dict[str, str] | None = None


def _short_name_index() -> dict[str, str]:
    """Map each unambiguous public CREVICE name to its documented location.

    Docstring type fields and signatures often use bare names (``Coord``,
    ``StructureFrame``, ``RadiusSet``). A name defined in exactly one
    ``crevice`` submodule (classes and functions by ``__module__``; other
    module-level names by where they are assigned) maps to ``module.name``.
    Ambiguous names are left unresolved.
    """
    global _SHORT_NAMES
    if _SHORT_NAMES is None:
        import importlib
        import inspect

        import crevice

        found: dict[str, set[str]] = {}
        for info in pkgutil.iter_modules(crevice.__path__):
            if info.name in API_EXCLUDED_MODULES:
                continue
            module = importlib.import_module(f"crevice.{info.name}")
            for name, value in vars(module).items():
                if name.startswith("_"):
                    continue
                home = getattr(value, "__module__", None) if inspect.isclass(value) or inspect.isfunction(value) else None
                if home is None and not inspect.ismodule(value):
                    source = Path(getattr(module, "__file__", "") or "")
                    text = source.read_text(encoding="utf-8") if source.is_file() else ""
                    home = module.__name__ if re.search(rf"^{re.escape(name)}\s*[:=]", text, flags=re.M) else None
                if home and home.startswith("crevice."):
                    found.setdefault(name, set()).add(f"{home}.{name}")
        _SHORT_NAMES = {name: next(iter(places)) for name, places in found.items() if len(places) == 1}
    return _SHORT_NAMES


#: Bare names of standard-library and NumPy types used in docstrings and
#: annotations, resolved through intersphinx.
_EXTERNAL_SHORT_NAMES = {"Path": "pathlib.Path", "ndarray": "numpy.ndarray",
                         "MDAnalysis.Universe": "MDAnalysis.core.universe.Universe",
                         "MDAnalysis.AtomGroup": "MDAnalysis.core.groups.AtomGroup"}


def _resolve_short_name(app: Sphinx, env, node: nodes.Element, contnode: nodes.Element):
    """Resolve a bare ``Name`` (or ``~Name``) used as a Python type to its CREVICE definition."""
    if node.get("refdomain") != "py" or node.get("reftype") not in {"class", "obj", "data"}:
        return None
    target = str(node.get("reftarget", "")).lstrip("~")
    if target not in _EXTERNAL_SHORT_NAMES and not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", target):
        return None
    if target in _EXTERNAL_SHORT_NAMES:
        from sphinx.ext.intersphinx import missing_reference

        alias = node.deepcopy()
        alias["reftarget"] = _EXTERNAL_SHORT_NAMES[target]
        return missing_reference(app, env, alias, contnode)
    full = _short_name_index().get(target)
    if full is None:
        return None
    domain = env.get_domain("py")
    alias = node.deepcopy()
    alias["reftarget"] = full
    return domain.resolve_xref(env, node["refdoc"], app.builder, "obj", full, alias, contnode)


_TUPLE_DEFAULT = re.compile(r"^(?P<indent>\s*)(?P<head>[*\w][\w*, ]* : .*?), default (?P<value>[(\[].*[)\]])\s*$")


def _move_tuple_defaults(app: Sphinx, what, name, obj, options, lines: list[str]) -> None:
    """Move tuple or list defaults out of NumPy type fields, before napoleon reads them.

    With ``napoleon_preprocess_types`` a type such as ``tuple of float,
    default (0.1, 0.9)`` is split at its commas and parentheses into broken
    references and literals. The field becomes ``..., optional`` and the
    description gains a first line stating the default as a literal; the
    source docstring is unchanged.
    """
    result: list[str] = []
    for line in lines:
        match = _TUPLE_DEFAULT.match(line)
        if not match:
            result.append(line)
            continue
        indent = match["indent"]
        result.append(f"{indent}{match['head']}, optional")
        result.append(f"{indent}    Default: ``{match['value']}``.")
    lines[:] = result


def setup(app: Sphinx) -> dict:
    app.add_config_value("crevice_repository_url", "", "env")
    # Run before Sphinx's download-file collector (default priority 500) so
    # repository files outside docs/ are never copied into the build.
    app.connect("doctree-read", _rewrite_outside_links, priority=400)
    app.connect("builder-inited", _check_api_module_index)
    app.connect("missing-reference", _resolve_public_alias)
    app.connect("missing-reference", _resolve_short_name)
    # Before sphinx.ext.napoleon (default priority 500) converts the docstring.
    app.connect("autodoc-process-docstring", _move_tuple_defaults, priority=400)
    return {"version": "0.1", "parallel_read_safe": True, "parallel_write_safe": True}

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

6. ``crevice-cli`` and ``crevice-cli-options`` directives render the command
   line reference from ``crevice.cli.build_parser()``: one compact page per
   command (its own arguments and options as tables) plus one shared page for
   the option families that many commands accept (structure input, atomic
   radii, figure text, channel-profile and hydration settings).  A build check
   warns when a command has no page under ``reference/cli/``.
"""

from __future__ import annotations

import os
import pkgutil
import re
from pathlib import Path, PurePosixPath

from docutils import nodes
from docutils.parsers.rst import Directive
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


# -- Command-line reference ---------------------------------------------------

#: Option families documented once on the shared-options page. A command's
#: option is shown there instead of on the command's own page when its help
#: text is the family's canonical text and the command has at least
#: ``min_members`` options of the family.
CLI_FAMILIES = {
    "structure-input": ("Structure input", 3, ["cache_dir", "input_format", "assembly", "offline", "refetch",
                                              "timeout", "model", "parser", "md_selection", "chain_ids"]),
    "radii": ("Atomic radii", 1, ["radii", "allow_radii_mismatch"]),
    "figures": ("Figures and viewer scenes", 1, ["annotate", "mouth_guides", "exit_centre_lines", "smooth",
                                                 "surface_smoothing"]),
    "channel-profile": ("Channel profile settings", 5, ["axis", "origin", "section_spacing", "enclosure_radius",
                                                         "lateral_exits", "exit_bulk_radius", "exit_spacing",
                                                         "samples", "padding", "search_radius", "refinement_steps",
                                                         "probe_radius", "include_hydrogen", "exclude_hetero"]),
    "hydration": ("Hydration settings", 3, ["hydration_cutoff", "hydration_probe", "hydration_sasa_points",
                                            "water_density_spacing", "water_density_smoothing",
                                            "water_density_level", "skip_water_density", "skip_interactions",
                                            "skip_analysis_views", "skip_hydration"]),
}


def _cli_commands():
    """Return ``{name: (subparser, one-line help)}`` for every ``crevice`` subcommand, in parser order."""
    import argparse

    from crevice.cli import build_parser

    parser = build_parser()
    sub = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    helps = {a.dest: a.help for a in sub._choices_actions}
    return {name: (p, helps.get(name, "")) for name, p in sub.choices.items()}


def _cli_options(subparser):
    return [a for a in subparser._actions if a.option_strings and a.dest != "help"]


def _canonical_help() -> dict[str, str]:
    """The most common help text of each family option across all commands."""
    from collections import Counter

    counts: dict[str, Counter] = {}
    for parser, _ in _cli_commands().values():
        for action in _cli_options(parser):
            counts.setdefault(action.dest, Counter())[action.help or ""] += 1
    return {dest: c.most_common(1)[0][0] for dest, c in counts.items()}


def _shared_families(subparser) -> dict[str, list]:
    """Map family key to the command's actions documented on the shared page."""
    canonical = _canonical_help()
    options = _cli_options(subparser)
    found: dict[str, list] = {}
    for key, (_, min_members, dests) in CLI_FAMILIES.items():
        members = [a for a in options if a.dest in dests and (a.help or "") == canonical.get(a.dest)]
        if len(members) >= min_members:
            found[key] = members
    return found


def _option_label(action) -> str:
    import argparse

    label = ", ".join(action.option_strings)
    takes_value = not isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction,
                                          argparse._StoreConstAction, argparse._CountAction))
    if takes_value:
        metavar = action.metavar or (None if action.choices else _plain_metavar(action))
        if metavar is None:
            metavar = "{" + ",".join(str(c) for c in action.choices) + "}"
        label += f" {metavar}"
    return label


_PATH_DEST = re.compile(r"(output|csv|pdb|png|json|pymol|dx|dir|npz|file|path|definition|scene)$")


def _plain_metavar(action) -> str:
    """A short, readable placeholder for an option value without a declared metavar."""
    if action.type is int:
        return "N"
    if action.type is float:
        return "X"
    if _PATH_DEST.search(action.dest):
        return "PATH"
    return action.dest.upper()


_OPTION_TOKEN = re.compile(r"(?<![\w-])(--?[A-Za-z][\w-]*)")


def _help_nodes(text: str) -> list[nodes.Node]:
    """Help text with option names as literals (also keeps ``--`` out of smart-dash conversion)."""
    result: list[nodes.Node] = []
    for i, part in enumerate(_OPTION_TOKEN.split(text)):
        if part:
            result.append(nodes.literal("", part) if i % 2 else nodes.Text(part))
    return result


def _default_text(action) -> str:
    import argparse

    if action.required:
        return "required"
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "off"
    if action.default is None or action.default is argparse.SUPPRESS:
        return "\u2013"
    return str(action.default)


def _doc_link(target: str, text: str, *, literal: bool = False, ref: bool = False) -> nodes.Node:
    content = nodes.literal("", text) if literal else nodes.inline("", text)
    return addnodes.pending_xref("", content, refdomain="std", reftype="ref" if ref else "doc",
                                 reftarget=target, refexplicit=True, refwarn=True)


def _table(headers: list[str], rows: list[list[nodes.Node]], widths: list[int], classes: list[str]) -> nodes.table:
    table = nodes.table(classes=["crevice-cli-table", *classes])
    group = nodes.tgroup(cols=len(headers))
    table += group
    for width in widths:
        group += nodes.colspec(colwidth=width)
    head = nodes.thead()
    group += head
    row = nodes.row()
    for title in headers:
        row += nodes.entry("", nodes.paragraph("", title))
    head += row
    body = nodes.tbody()
    group += body
    for cells in rows:
        row = nodes.row()
        for cell in cells:
            row += nodes.entry("", cell if isinstance(cell, nodes.paragraph) else nodes.paragraph("", "", cell))
        body += row
    return table


def _option_rows(actions) -> list[list[nodes.Node]]:
    rows = []
    for action in actions:
        help_text = (action.help or "").strip()
        if action.choices and action.metavar is None and not help_text.lower().startswith("choices"):
            help_text = f"{help_text} (choices: {', '.join(str(c) for c in action.choices)})".strip()
        rows.append([nodes.literal("", _option_label(action)), nodes.Text(_default_text(action)),
                     nodes.paragraph("", "", *_help_nodes(help_text))])
    return rows


class CreviceCliDirective(Directive):
    """Render one ``crevice`` subcommand: usage, arguments and its own options."""

    required_arguments = 1
    has_content = False

    def run(self) -> list[nodes.Node]:
        name = self.arguments[0]
        commands = _cli_commands()
        if name not in commands:
            raise self.error(f"crevice has no command {name!r}")
        parser, _summary = commands[name]
        positionals = [a for a in parser._actions if not a.option_strings]
        shared = _shared_families(parser)
        shared_ids = {id(a) for members in shared.values() for a in members}
        own = [a for a in _cli_options(parser) if id(a) not in shared_ids]

        usage_parts = ["crevice", name]
        for action in positionals:
            meta = (action.metavar or action.dest).upper()
            usage_parts.append(f"{meta} ..." if action.nargs in ("+", "*") else meta)
        usage_parts.append("[options]")
        result: list[nodes.Node] = [nodes.literal_block("", " ".join(usage_parts), language="text",
                                                        classes=["crevice-cli-usage"])]
        if positionals:
            rows = [[nodes.literal("", (a.metavar or a.dest).upper()),
                     nodes.paragraph("", "", *_help_nodes((a.help or "").strip()))]
                    for a in positionals]
            result.append(nodes.rubric("", "Arguments"))
            result.append(_table(["Argument", "Description"], rows, [25, 75], ["crevice-cli-args"]))
        if own:
            result.append(nodes.rubric("", "Options"))
            result.append(_table(["Option", "Default", "Description"], _option_rows(own), [30, 12, 58],
                                 ["crevice-cli-options"]))
        if shared:
            para = nodes.paragraph(classes=["crevice-cli-shared"])
            para += nodes.strong("", "Shared options: ")
            for i, key in enumerate(shared):
                if i:
                    para += nodes.Text(", ")
                para += _doc_link(f"cli-{key}", CLI_FAMILIES[key][0].lower(), ref=True)
            para += nodes.Text(". These are described once, on the shared options page.")
            result.append(para)
        return result


class CreviceCliOptionsDirective(Directive):
    """Render one shared option family and the commands that accept it."""

    required_arguments = 1
    has_content = False

    def run(self) -> list[nodes.Node]:
        key = self.arguments[0]
        if key not in CLI_FAMILIES:
            raise self.error(f"unknown option family {key!r}")
        dests = CLI_FAMILIES[key][2]
        users: list[str] = []
        actions: dict[str, object] = {}
        for name, (parser, _) in _cli_commands().items():
            members = _shared_families(parser).get(key)
            if members:
                users.append(name)
                for action in members:
                    actions.setdefault(action.dest, action)
        ordered = [actions[d] for d in dests if d in actions]
        para = nodes.paragraph(classes=["crevice-cli-users"])
        para += nodes.strong("", "Accepted by: ")
        for i, name in enumerate(users):
            if i:
                para += nodes.Text(", ")
            para += _doc_link(f"/reference/cli/{name}", name, literal=True)
        return [para, _table(["Option", "Default", "Description"], _option_rows(ordered), [30, 12, 58],
                             ["crevice-cli-options"])]


def _check_cli_pages(app: Sphinx) -> None:
    folder = Path(app.srcdir) / "reference" / "cli"
    pages = {p.stem for p in folder.glob("*.md")} - {"index", "common-options"}
    names = set(_cli_commands())
    for name in sorted(names - pages):
        logger.warning(f"crevice command {name!r} has no page reference/cli/{name}.md", type="crevice", subtype="cli")
    for name in sorted(pages - names):
        logger.warning(f"reference/cli/{name}.md documents no crevice command", type="crevice", subtype="cli")


def setup(app: Sphinx) -> dict:
    app.add_directive("crevice-cli", CreviceCliDirective)
    app.add_directive("crevice-cli-options", CreviceCliOptionsDirective)
    app.connect("builder-inited", _check_cli_pages)
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

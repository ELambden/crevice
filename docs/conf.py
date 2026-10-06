"""Sphinx configuration for the CREVICE documentation site.

The Sphinx source directory is ``docs/`` itself, so the Markdown method notes
(``docs/HYDRATION_METHODS.md`` and friends) are built in place by MyST rather
than copied.  Site-specific pages live in subdirectories.  Only the files matched
by ``include_patterns`` below are part of the site.

Build locally (see ``docs/development/documentation.md``)::

    python -m sphinx -W --keep-going -b html docs docs/_build/html
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent

# Import CREVICE from the source tree so API pages always describe this checkout,
# even when an older installed copy is present in the environment.
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(DOCS_DIR / "_ext"))

import crevice  # noqa: E402  (import must follow the sys.path change)

# -- Project information -----------------------------------------------------

project = "CREVICE"
# The acronym expansion is written here, in README.md, CITATION.cff, pyproject.toml,
# src/crevice/__init__.py and the CLI description; keep them identical.
full_name = "Cavity and Residue Environment Visualisation, Interaction and Connectivity Evaluation"
author = "CREVICE contributors"
copyright = "2026, Edward Lambden"
release = crevice.__version__
version = release

# -- Repository link ---------------------------------------------------------
# Set CREVICE_REPOSITORY_URL (for example in the Read the Docs project settings)
# to enable "view source" links and repository links for files outside docs/.
CREVICE_REPOSITORY_URL = os.environ.get("CREVICE_REPOSITORY_URL", "https://github.com/ELambden/crevice")

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinxarg.ext",
    "sphinx_copybutton",
    "crevice_docs",  # local helpers: docs/_ext/crevice_docs.py
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
# The site is an explicit allowlist: a Markdown file in docs/ that is not
# matched here is not built (and cannot be published by accident).
include_patterns = [
    "index.md",
    "getting-started/*.md",
    "tools/*.md",
    "methods/*.md",
    "examples/*.md",
    "reference/*.md",
    "api/*.md",
    "api/*.rst",
    "api/generated/*.rst",
    "development/*.md",
    "project/*.md",
    "BIOLOGICAL_REGIONS.md",
    "CAVITY_OBSTACLES.md",
    "CAVITY_TRAJECTORY_METHODS.md",
    "HYDRATION_METHODS.md",
    "HYDRATION_NETWORK_VIEWS.md",
    "REGION_WORKFLOW.md",
    "SURFACE_SMOOTHING.md",
    "TRAJECTORY_RESIDUE_METHODS.md",
]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "requirements.txt"]
templates_path = ["_templates"]

# MyST: allow ::: fences, definition lists and field lists, and generate heading
# anchors so the existing method notes can be linked section by section.
myst_enable_extensions = ["colon_fence", "deflist", "fieldlist", "substitution"]
myst_heading_anchors = 3
# sphinx-argparse registers a domain without `resolve_any_xref`; restrict the
# domains MyST searches for bare Markdown links to avoid spurious warnings.
myst_ref_domains = ["std", "py"]
myst_substitutions = {
    "version": release,
    "full_name": full_name,
}

# -- Autodoc / autosummary ---------------------------------------------------

autosummary_generate = True
autosummary_imported_members = False
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autodoc_preserve_defaults = True
# Optional heavy dependencies are imported lazily inside CREVICE functions; mock
# them here so a docs-only environment never needs them at import time.
autodoc_mock_imports = ["gemmi", "mdtraj", "skimage", "networkx", "sklearn"]
napoleon_numpy_docstring = True
napoleon_google_docstring = True
# Figure and scene writers document every file and every colour/representation.
napoleon_custom_sections = [("Outputs", "params_style"), "Colours and representation"]
# Render "Attributes" sections as :ivar: fields. Dataclass fields are also
# collected by autodoc (undoc-members), so ".. attribute::" entries would be
# duplicate object descriptions and fail the -W build.
napoleon_use_ivar = True
# Nitpicky builds (sphinx -n): parse NumPy-style type fields so that
# "optional", "default 0.5", literal choices such as {"pdb", "cif"} and numbers
# are rendered as text rather than as cross-references, and map the informal
# type words used throughout the docstrings to real targets. (``Path`` is
# resolved by docs/_ext/crevice_docs.py instead: a napoleon alias for it would
# also rewrite the ``str | Path`` annotations in signatures.)
napoleon_preprocess_types = True
napoleon_type_aliases = {
    "array_like": ":term:`numpy:array_like`",
    "array-like": ":term:`numpy:array_like`",
    "ndarray": "numpy.ndarray",
    "sequence": "collections.abc.Sequence",
    "iterable": "collections.abc.Iterable",
    "mapping": "collections.abc.Mapping",
    "callable": "collections.abc.Callable",
}

# Words that are not Python objects but appear where napoleon expects a type:
# the file kinds of "Outputs" sections ("PREFIX_x.csv : table") and the array
# shapes of NumPy-style types ("array_like, shape (n, 3)"). They are plain
# descriptions, so they are exempt from nitpicky reference checking.
nitpick_ignore_regex = [
    (r"py:(class|obj)", r"(JSON|DX|PNG|BILD|VMD|ChimeraX|ChimeraX script|PyMOL|PyMOL script|image|figures?|tables?|"
                        r"arrays|grids|scenes?|scene data|structures?|network|tunnels|cavities|connectivity|profile|"
                        r"density|interactions|trajectory tables|volume scene|display|regions|data|text|CSV|NPZ|PDB|"
                        r"script|bundle|evidence table|member waters|snapshots|per frame|definition|review)"),
    (r"py:(class|obj)", r"(shape .*|\(.*|.*\)|\d.*|[a-z]|n[a-z]|three .*|float32|int32|int64|int8|bool_)"),
    # Dataclass signatures with ``field(default_factory=...)`` defaults cannot
    # be parsed by Sphinx, which then splits the generic annotations at commas.
    (r"py:(class|obj)", r"(dict|tuple|Mapping|typing\.Mapping)\[.*"),
    # Optional third-party types without an intersphinx inventory, and a
    # private helper named in a docstring.
    (r"py:(class|obj)", r"(gemmi\..*|sklearn\..*|_SphereQueries)"),
]

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "mdanalysis": ("https://docs.mdanalysis.org/stable/", None),
}
intersphinx_timeout = 20

# -- HTML output -------------------------------------------------------------

html_theme = "furo"
html_title = f"CREVICE {release}"
html_static_path = ["_static"]
html_css_files = ["crevice.css"]
html_theme_options = {}
if CREVICE_REPOSITORY_URL:
    html_theme_options.update(
        {"source_repository": CREVICE_REPOSITORY_URL, "source_branch": "main", "source_directory": "docs/"}
    )

# -- Link checking -----------------------------------------------------------

linkcheck_timeout = 20
linkcheck_retries = 2
linkcheck_anchors = False  # many publisher pages build anchors with JavaScript

# -- crevice_docs extension settings ---------------------------------------------

# Links to repository files outside docs/ (for example CHANGELOG.md) become
# repository links when a URL is configured, and labelled paths otherwise.
crevice_repository_url = CREVICE_REPOSITORY_URL

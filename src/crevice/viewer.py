"""Interactive 3D views of CREVICE results in Jupyter and in a standalone HTML file.

:func:`view` builds a :class:`CreviceView` from a ``publish``/``cast`` output
directory, a ``region-trajectory``/``cavity-trajectory`` hydration output
directory, or in-memory results (a structure frame plus a cast, its ENTRY/EXIT
segments and a channel profile). The view holds a plain scene description:
the protein as a cartoon, each displayed cast surface as its own triangle
mesh, lining residues, the profile centre line and member waters, each a layer
that can be shown or hidden.

Surfaces are extracted with ``skimage.measure.marching_cubes`` at the 0.5
isovalue of the same display fields that CREVICE writes for PyMOL, VMD and
ChimeraX (``PREFIX_<segment>_display.dx``, ``PREFIX_display.dx``), so the
interactive surfaces are the native-viewer display surfaces. They are display
geometry: measured volumes stay those of the binary ``*_volume.dx`` grids.

Two renderers share one scene loader (``_viewer_assets/crevice_view.js``) on
top of 3Dmol.js (BSD-3-Clause, Rego and Koes, Bioinformatics 2015):

* :meth:`CreviceView.show` displays the view in Jupyter/IPython through
  py3Dmol (the optional ``crevice[view]`` extra), with a checkbox per layer;
* :meth:`CreviceView.save_html` writes one self-contained HTML file with the
  vendored 3Dmol.js build inlined, which opens offline in a browser
  (``crevice view OUTDIR``).

An interactive file is an explicit opt-in output; no analysis command writes
HTML by default.

Examples
--------
>>> import numpy as np
>>> field = np.zeros((4, 4, 4)); field[1:3, 1:3, 1:3] = 1.0
>>> mesh = field_mesh(field, origin=(0.0, 0.0, 0.0), deltas=np.eye(3))
>>> mesh["vertices"].shape[1], mesh["faces"].shape[1]
(3, 3)
"""
from __future__ import annotations

import base64
import csv
import gzip
import html
import json
import math
from dataclasses import dataclass, field as _dc_field
from pathlib import Path

import numpy as np

#: Scene description schema identifier.
SCENE_SCHEMA = "crevice.interactive_view.v1"
#: Protein cartoon colour (grey-blue, as in the native scenes).
PROTEIN_HEX = "#8fa1b3"
#: Protein cartoon opacity (25% transparent).
PROTEIN_OPACITY = 0.75
#: ENTRY cast segment colour (yellow).
ENTRY_HEX = "#f0d43a"
#: LUMEN and unsegmented channel cast colour (teal).
LUMEN_HEX = "#1494a1"
#: EXIT cast segment colour (rust).
EXIT_HEX = "#c2410c"
#: Colour of every cast that is not a channel (violet).
NON_CHANNEL_HEX = "#a855f7"
#: Profile centre-line colour (blue, as the exit centre lines).
CENTRE_LINE_HEX = "#2978d6"
#: Profile centre-line tube radius (Å).
CENTRE_LINE_RADIUS_A = 0.12
#: Member-water sphere colour (orange, as the native hydration scenes).
WATER_HEX = "#ff8000"
#: Member-water sphere radius (Å).
WATER_RADIUS_A = 0.7
#: Lining-residue stick radius (Å), as in the native scenes.
STICK_RADIUS_A = 0.2

_ASSETS = Path(__file__).with_name("_viewer_assets")
_THREEDMOL = "3Dmol-min.js"
_THREEDMOL_LICENSE = "3Dmol-LICENSE.txt"
_THREEDMOL_VERSION = "2.5.5"
#: CDN location of the same 3Dmol.js build, used only with ``inline_js=False``.
THREEDMOL_CDN = f"https://cdn.jsdelivr.net/npm/3dmol@{_THREEDMOL_VERSION}/build/3Dmol-min.js"
#: Region-output surface opacities in the interactive view (the native scenes use
#: 0.18 and 0.42, which are nearly invisible behind a WebGL transparent cartoon).
_INTERACTIVE_OPACITY = {"cavity": 0.35, "water": 0.6}
_SEGMENT_HEX = {"entry": ENTRY_HEX, "lumen": LUMEN_HEX, "exit": EXIT_HEX}


def _hex(rgb) -> str:
    return "#" + "".join(f"{round(float(v) * 255):02x}" for v in rgb)


def read_dx_field(path):
    """Read a regular OpenDX scalar field (binary grid, display field or density).

    Parameters
    ----------
    path : str or Path
        OpenDX file with one ``gridpositions`` object, three ``delta`` lines
        and a ``data follows`` array in C order (z fastest), as CREVICE writes.

    Returns
    -------
    field : numpy.ndarray
        Values with the file's ``counts`` shape.
    origin : numpy.ndarray, shape (3,)
        Position of grid point (0, 0, 0), Å.
    deltas : numpy.ndarray, shape (3, 3)
        Grid step vectors as rows, Å.

    Raises
    ------
    ValueError
        When the header is incomplete or the value count does not match.
    """
    text = Path(path).read_text()
    head, sep, body = text.partition("data follows")
    if not sep:
        raise ValueError(f"{path}: no 'data follows' array")
    shape = origin = None
    deltas = []
    for line in head.splitlines():
        words = line.split()
        if "gridpositions" in line and "counts" in words:
            shape = tuple(int(v) for v in words[-3:])
        elif words[:1] == ["origin"]:
            origin = np.asarray([float(v) for v in words[1:4]])
        elif words[:1] == ["delta"]:
            deltas.append([float(v) for v in words[1:4]])
    if shape is None or origin is None or len(deltas) != 3:
        raise ValueError(f"{path}: incomplete OpenDX header")
    size = int(np.prod(shape))
    values = body.split("\n", 1)[1] if "\n" in body else ""
    cut = min((i for i in (values.find("attribute"), values.find("object"), values.find("component")) if i >= 0),
              default=len(values))
    data = np.array(values[:cut].split(), dtype=float)
    if data.size != size:
        raise ValueError(f"{path}: {data.size} values for counts {shape}")
    return data.reshape(shape), origin, np.asarray(deltas)


def field_mesh(field, origin, deltas, level=0.5):
    """Triangle mesh of one isosurface of a regular 3D field, in world coordinates.

    Parameters
    ----------
    field : array_like, shape (nx, ny, nz)
        Scalar field (CREVICE display field or binary grid).
    origin : array_like, shape (3,)
        World position of grid point (0, 0, 0), Å.
    deltas : array_like, shape (3, 3)
        Grid step vectors as rows, Å (rotated lattices allowed).
    level : float, default 0.5
        Isovalue; 0.5 is the CREVICE display surface.

    Returns
    -------
    dict or None
        ``vertices`` (N, 3) float32 Å, ``faces`` (M, 3) uint32 and unit
        ``normals`` (N, 3) float32 pointing out of the enclosed volume; ``None``
        when ``level`` is not crossed by the field.

    Notes
    -----
    Vertices are ``origin + ijk @ deltas`` for the marching-cubes grid
    coordinates ``ijk``, the transform the native PyMOL/VMD/ChimeraX scenes
    apply to the same field. Faces are re-wound when ``deltas`` is a reflection.
    """
    from skimage.measure import marching_cubes
    values = np.asarray(field, dtype=float)
    if values.ndim != 3 or not float(values.min()) < level < float(values.max()):
        return None
    deltas = np.asarray(deltas, dtype=float)
    padded = np.pad(values, 1, constant_values=min(float(values.min()), level - 1.0))
    ijk, faces, normals, _ = marching_cubes(padded, level=level)
    vertices = np.asarray(origin, float) + (ijk - 1.0) @ deltas
    normals = normals @ np.linalg.inv(deltas).T  # skimage normals point down the gradient (outward)
    normals /= np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-12)
    if np.linalg.det(deltas) > 0:
        faces = faces[:, ::-1]  # counter-clockwise seen from outside, the WebGL front face
    return {"vertices": vertices.astype(np.float32), "faces": faces.astype(np.uint32),
            "normals": normals.astype(np.float32)}


def cast_mesh(cast, smooth=None, smooth_supersample=None):
    """Display mesh of an in-memory cast, computed as the native scenes compute it.

    Parameters
    ----------
    cast : VoidCast
        Cast (or one segment wrapped as a cast) on a CREVICE lattice.
    smooth : float, optional
        Display smoothing length (Å), as ``--smooth``; default
        :data:`crevice.volume_export.DEFAULT_SMOOTH_A`.
    smooth_supersample : int, optional
        Supersampling factor for ``smooth`` (default 2).

    Returns
    -------
    dict or None
        As :func:`field_mesh`.
    """
    from .volume_export import cast_grid, resolve_display_smoothing
    smoothing = resolve_display_smoothing(smooth, smooth_supersample=smooth_supersample)
    binary, origin = cast_grid(cast)
    display, spacing, _ = smoothing.field(binary, cast.spacing)
    basis = np.asarray(cast.metadata.get("grid_basis", np.eye(3)), dtype=float)
    return field_mesh(display, origin, basis * spacing)


@dataclass
class Layer:
    """One toggleable part of an interactive scene.

    Parameters
    ----------
    name : str
        Identifier used by :meth:`CreviceView.toggle` (for example ``"lumen"``,
        ``"entry_1"``, ``"lining_lumen"``, ``"centre_line"``, ``"waters"``).
    label : str
        Checkbox label in the viewer controls (never drawn into the canvas).
    kind : str
        ``"cartoon"``, ``"mesh"``, ``"sticks"``, ``"line"`` or ``"spheres"``.
    color : str
        Hex colour (carbon colour for sticks).
    visible : bool, default True
        Shown when the viewer opens.
    opacity : float, default 1.0
        Surface or sphere opacity.
    data : dict
        Kind-specific geometry: ``mesh`` arrays, ``residues`` (list of
        ``{"chain", "resi"}``), ``paths`` (lists of points) or ``points``.
    description : str, default ""
        Tooltip text in the controls.
    """

    name: str
    label: str
    kind: str
    color: str
    visible: bool = True
    opacity: float = 1.0
    data: dict = _dc_field(default_factory=dict)
    description: str = ""


def _encode(array, dtype, compress):
    raw = np.ascontiguousarray(array, dtype=dtype).tobytes()
    if compress:
        raw = gzip.compress(raw, compresslevel=9, mtime=0)
    return {"dtype": "uint32" if dtype == np.uint32 else "float32",
            "encoding": "gzip+base64" if compress else "base64",
            "data": base64.b64encode(raw).decode("ascii")}


def _rotation_quaternion(matrix):
    """Unit quaternion (x, y, z, w) of a proper rotation matrix."""
    m = np.asarray(matrix, dtype=float)
    w = math.sqrt(max(0.0, 1.0 + m[0, 0] + m[1, 1] + m[2, 2])) / 2
    x = math.copysign(math.sqrt(max(0.0, 1.0 + m[0, 0] - m[1, 1] - m[2, 2])) / 2, m[2, 1] - m[1, 2])
    y = math.copysign(math.sqrt(max(0.0, 1.0 - m[0, 0] + m[1, 1] - m[2, 2])) / 2, m[0, 2] - m[2, 0])
    z = math.copysign(math.sqrt(max(0.0, 1.0 - m[0, 0] - m[1, 1] + m[2, 2])) / 2, m[1, 0] - m[0, 1])
    q = np.asarray([x, y, z, w])
    return (q / np.linalg.norm(q)).tolist()


class CreviceView:
    """Interactive scene of CREVICE results: protein, cast surfaces and overlays.

    Build one with :func:`view`. Every part is a :class:`Layer` that can be
    switched on or off with :meth:`toggle` (before display) or with its checkbox
    (in the displayed viewer).

    Parameters
    ----------
    layers : list of Layer
        Scene layers in display order.
    structure : str, optional
        Protein coordinates as PDB text (cartoon and sticks are drawn from it).
    title : str, default "CREVICE view"
        HTML page title; it is not drawn into the scene.
    camera : dict, optional
        ``center`` (Å), ``rotation`` (3x3, model to view) and ``zoom``; the
        default fits the whole protein. Directory views use
        the native scene rotation and frame the displayed surfaces.
    source : str, optional
        Where the results came from (recorded in the scene description).
    notes : list of str, optional
        Interpretation notes recorded in the scene description.

    Colours and representation
    --------------------------
    * Protein: cartoon, grey-blue ``#8fa1b3``, 25% transparent.
    * Channel cast segments: ENTRY yellow ``#f0d43a``, LUMEN teal ``#1494a1``,
      EXIT rust ``#c2410c``, one opaque surface per segment; an unsegmented
      channel cast is teal and every other cast violet ``#a855f7``.
      Hydration (region) outputs keep their native colours: reference cavity
      pale grey at 35% opacity, mean water density teal at 60% (the
      native scenes use 18% and 42%).
    * Lining residues: sticks of radius 0.2 Å, carbons in the segment colour
      (or the residue's hydration colour), other atoms by element; hidden by
      default.
    * Profile centre line: blue ``#2978d6`` tubes of 0.12 Å through the axial
      profile samples and any lateral exit legs; hidden by default.
    * Member waters: orange ``#ff8000`` spheres of 0.7 Å, shown by default.
    * White background; no text, labels or legends inside the 3D canvas. The
      layer names appear only as checkbox labels beside it.
    """

    def __init__(self, layers, structure=None, title="CREVICE view", camera=None, source=None, notes=None):
        self.layers = list(layers)
        self.structure = structure
        self.title = title
        self.camera = camera
        self.source = source
        self.notes = list(notes or ())
        names = [layer.name for layer in self.layers]
        if len(set(names)) != len(names):
            raise ValueError(f"layer names must be unique: {names}")

    @property
    def layer_names(self):
        """Names of the layers in display order (the arguments of :meth:`toggle`).

        Returns
        -------
        list of str
            Layer names.
        """
        return [layer.name for layer in self.layers]

    def layer(self, name):
        """Return one layer by name.

        Parameters
        ----------
        name : str
            Layer name (see :attr:`layer_names`).

        Returns
        -------
        Layer
            The named layer.

        Raises
        ------
        KeyError
            For an unknown name; the message lists the valid names.
        """
        for layer in self.layers:
            if layer.name == name:
                return layer
        raise KeyError(f"no layer {name!r}; layers are {', '.join(self.layer_names)}")

    def toggle(self, name, visible=None):
        """Show, hide or flip one layer (or every layer of one kind of segment).

        Parameters
        ----------
        name : str
            Layer name, or a group: ``"entry"``, ``"exit"`` or ``"lining"``
            match every ``entry_N``, ``exit_N`` or ``lining_*`` layer.
        visible : bool, optional
            New state; default flips the current state.

        Returns
        -------
        CreviceView
            ``self``, so calls can be chained (``v.toggle("lining").show()``).
        """
        matches = [layer for layer in self.layers
                   if layer.name == name or layer.name.startswith(name + "_")]
        if not matches:
            self.layer(name)  # raises KeyError with the valid names
        for layer in matches:
            layer.visible = (not layer.visible) if visible is None else bool(visible)
        return self

    def scene(self, compress=True):
        """JSON-ready scene description consumed by the JavaScript loader.

        Parameters
        ----------
        compress : bool, default True
            Gzip mesh arrays before base64 encoding (decoded in the browser with
            ``DecompressionStream``); ``False`` writes plain base64.

        Returns
        -------
        dict
            ``schema``, ``title``, ``source``, ``notes``, ``structure``
            (``format``/``data``), ``camera`` and ``layers``.
        """
        layers = []
        for layer in self.layers:
            row = {"name": layer.name, "label": layer.label, "kind": layer.kind, "color": layer.color,
                   "visible": bool(layer.visible), "opacity": float(layer.opacity),
                   "description": layer.description}
            if layer.kind == "mesh":
                mesh = layer.data
                row.update(vertices=_encode(mesh["vertices"], np.float32, compress),
                           faces=_encode(mesh["faces"], np.uint32, compress),
                           normals=_encode(mesh["normals"], np.float32, compress),
                           vertex_count=int(len(mesh["vertices"])), face_count=int(len(mesh["faces"])))
            elif layer.kind == "sticks":
                row.update(residues=layer.data["residues"], radius=STICK_RADIUS_A)
            elif layer.kind == "line":
                row.update(paths=[[[round(float(v), 3) for v in p] for p in path] for path in layer.data["paths"]],
                           radius=CENTRE_LINE_RADIUS_A)
            elif layer.kind == "spheres":
                row.update(points=[[round(float(v), 3) for v in p] for p in layer.data["points"]],
                           radius=float(layer.data.get("radius", WATER_RADIUS_A)))
            layers.append(row)
        return {"schema": SCENE_SCHEMA, "title": self.title, "source": self.source, "notes": self.notes,
                "structure": {"format": "pdb", "data": self.structure} if self.structure else None,
                "camera": self.camera, "layers": layers}

    def to_html(self, inline_js=True, width="100%", height="100vh"):
        """Standalone HTML page of the view, with layer checkboxes beside the canvas.

        Parameters
        ----------
        inline_js : bool, default True
            Inline the vendored 3Dmol.js build so the page works offline;
            ``False`` loads it from :data:`THREEDMOL_CDN` (network needed).
        width : str, default "100%"
            CSS width of the page body.
        height : str, default "100vh"
            CSS height of the viewer.

        Returns
        -------
        str
            Complete HTML document.
        """
        if inline_js:
            library = (_ASSETS / _THREEDMOL).read_text(encoding="utf-8").replace("</script", "<\\/script")
            notice = (_ASSETS / _THREEDMOL_LICENSE).read_text(encoding="utf-8").replace("-->", "- ->")
            library_tag = f"<!--\n{notice}\n-->\n<script>\n{library}\n</script>"
        else:
            library_tag = f'<script src="{THREEDMOL_CDN}"></script>'
        loader = (_ASSETS / "crevice_view.js").read_text(encoding="utf-8")
        scene = json.dumps(self.scene(), separators=(",", ":")).replace("</", "<\\/")
        return _PAGE.format(title=html.escape(self.title), library=library_tag, loader=loader, scene=scene,
                            width=width, height=height, css=_CSS)

    def save_html(self, path, inline_js=True):
        """Write the view as one self-contained interactive HTML file.

        Parameters
        ----------
        path : str or Path
            Output ``.html`` file (parent directories are created).
        inline_js : bool, default True
            Inline 3Dmol.js (works offline, adds about 0.54 MB); ``False``
            loads it from a CDN when the page is opened.

        Returns
        -------
        Path
            The written file.

        Outputs
        -------
        path : HTML
            Page with the 3D canvas, one checkbox per layer, *Reset view* and
            *Save PNG* buttons. Opening it runs JavaScript in the browser;
            nothing is sent over the network when ``inline_js`` is true.
        """
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(self.to_html(inline_js=inline_js), encoding="utf-8")
        return out

    def jupyter_html(self, width=720, height=520, offline=False):
        """HTML fragment for a notebook cell, rendered through py3Dmol.

        Parameters
        ----------
        width : int, default 720
            Canvas width in pixels.
        height : int, default 520
            Canvas height in pixels.
        offline : bool, default False
            Load the vendored 3Dmol.js as a data URI instead of py3Dmol's CDN
            (larger notebook, no network needed).

        Returns
        -------
        str
            HTML with py3Dmol's viewer and the layer checkboxes.

        Raises
        ------
        ImportError
            When py3Dmol is not installed (``pip install 'crevice[view]'``).
        """
        py3Dmol = _import_py3dmol()
        options = {"width": width, "height": height}
        if offline:
            data = base64.b64encode((_ASSETS / _THREEDMOL).read_bytes()).decode("ascii")
            options["js"] = "data:text/javascript;base64," + data
        widget = py3Dmol.view(**options)
        fragment = widget.write_html()
        uid = widget.uniqueid
        loader = (_ASSETS / "crevice_view.js").read_text(encoding="utf-8")
        scene = json.dumps(self.scene(), separators=(",", ":")).replace("</", "<\\/")
        controls = f"crevice_controls_{uid}"
        return (f'<style>{_CSS}</style><div class="crevice-row">{fragment}'
                f'<div class="crevice-controls" id="{controls}"></div></div>'
                f"<script>\n{loader}\n$3Dmolpromise.then(function() {{ CreviceView.build(viewer_{uid}, "
                f"document.getElementById('{controls}'), {scene}); }});\n</script>")

    def show(self, width=720, height=520, offline=False):
        """Display the view in Jupyter/IPython with a checkbox per layer.

        Parameters
        ----------
        width : int, default 720
            Canvas width in pixels.
        height : int, default 520
            Canvas height in pixels.
        offline : bool, default False
            Embed 3Dmol.js instead of loading it from py3Dmol's CDN.

        Returns
        -------
        None
            The viewer is published to the notebook output.

        Raises
        ------
        ImportError
            Without py3Dmol (``pip install 'crevice[view]'``) or IPython.
        """
        fragment = self.jupyter_html(width=width, height=height, offline=offline)
        try:
            from IPython.display import publish_display_data
        except ImportError as exc:  # pragma: no cover - IPython ships with py3Dmol
            raise ImportError("CreviceView.show() needs IPython/Jupyter; use save_html() elsewhere") from exc
        publish_display_data({"application/3dmoljs_load.v0": fragment, "text/html": fragment}, metadata={})

    def _ipython_display_(self):
        self.show()

    def __repr__(self):
        shown = ", ".join(f"{layer.name}{'' if layer.visible else ' (hidden)'}" for layer in self.layers)
        return f"CreviceView({shown})"


def _import_py3dmol():
    try:
        import py3Dmol
    except ImportError as exc:
        raise ImportError("Showing a CREVICE view in Jupyter needs py3Dmol: pip install 'crevice[view]' "
                          "(or pip install py3Dmol). CreviceView.save_html() and 'crevice view' "
                          "work without it.") from exc
    return py3Dmol


_CSS = """
.crevice-row{display:flex;gap:12px;align-items:flex-start;font:13px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif;color:#1f2933}
.crevice-controls{min-width:170px;display:flex;flex-direction:column;gap:4px}
.crevice-layer{display:flex;align-items:center;gap:6px;cursor:pointer;user-select:none}
.crevice-swatch{width:12px;height:12px;border-radius:2px;border:1px solid #0003;flex:none}
.crevice-buttons{display:flex;gap:6px;margin-top:8px}
.crevice-buttons button{font:inherit;padding:3px 8px;border:1px solid #9aa5b1;border-radius:4px;background:#f5f7fa;cursor:pointer}
"""

_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
html,body{{margin:0;background:#ffffff;width:{width}}}
{css}
.crevice-page{{display:flex;flex-wrap:wrap;gap:0;align-items:stretch}}
#crevice-viewer{{position:relative;flex:1 1 480px;height:{height};min-height:320px}}
.crevice-page .crevice-controls{{padding:12px 16px;box-sizing:border-box;font:13px/1.4 system-ui,-apple-system,"Segoe UI",sans-serif;color:#1f2933}}
.crevice-credit{{margin-top:10px;color:#52606d;font-size:11px}}
</style>
{library}
</head>
<body>
<div class="crevice-page">
<div id="crevice-viewer"></div>
<div><div class="crevice-controls" id="crevice-controls"></div>
<div class="crevice-controls crevice-credit">CREVICE interactive view; rendered with 3Dmol.js (BSD-3-Clause).</div></div>
</div>
<script>
{loader}
var creviceScene = {scene};
var creviceViewer = $3Dmol.createViewer(document.getElementById("crevice-viewer"), {{backgroundColor: "white"}});
window.creviceReady = CreviceView.build(creviceViewer, document.getElementById("crevice-controls"), creviceScene)
  .then(function (handle) {{ window.creviceHandle = handle; document.body.dataset.creviceReady = "1"; return handle; }});
</script>
</body>
</html>
"""


# --------------------------------------------------------------------------
# Scene construction
# --------------------------------------------------------------------------

def _resolve(path_text, directory):
    """Manifest paths are absolute at write time; fall back to the bundle directory."""
    if not path_text:
        return None
    path = Path(path_text)
    if path.is_file():
        return path
    local = Path(directory) / path.name
    return local if local.is_file() else None


def _mesh_layer(name, label, color, mesh, opacity=1.0, visible=True, description=""):
    return Layer(name, label, "mesh", color, visible=visible, opacity=opacity, data=mesh, description=description)


def _segment_label(name):
    kind, _, index = name.partition("_")
    return kind.upper() + (f" {index}" if index else "")


def _profile_paths(profile):
    """Axial profile samples plus lateral exit leg samples, as point lists."""
    if profile is None:
        return []
    if isinstance(profile, dict):
        points = [[p["x"], p["y"], p["z"]] for p in profile.get("points", ())]
        metadata = profile.get("metadata") or {}
    else:
        points = [list(p.position) for p in profile.points]
        metadata = profile.metadata or {}
    paths = [points] if len(points) > 1 else []
    for end in (metadata.get("exits") or {}).values():
        if not isinstance(end, dict):
            continue
        for leg in end.get("legs") or ():
            samples = leg.get("samples") or ()
            leg_points = [s.get("position_A") or [s.get("x"), s.get("y"), s.get("z")] for s in samples]
            leg_points = [p for p in leg_points if p is not None and None not in p]
            if len(leg_points) > 1:
                paths.append(leg_points)
    return paths


def _camera(rotation=None, center=None, span=None, coords=None, focus=None):
    """Camera dict for the loader: native rotation, centre and a zoom that frames ``focus`` points."""
    if rotation is None or center is None:
        return None
    zoom = 1.0
    if coords is not None and len(coords):
        extent = float(np.max(np.ptp(np.asarray(coords, float), axis=0)))
        target = span
        if focus is not None and len(focus):
            focus = np.asarray(focus, float)
            target = 1.6 * float(np.max(np.ptp(focus, axis=0)))
            center = focus.min(axis=0) / 2 + focus.max(axis=0) / 2
        if target:
            zoom = float(np.clip(extent / float(target), 1.0, 6.0))
    return {"center": [float(v) for v in center], "rotation": np.asarray(rotation, float).tolist(),
            "quaternion": _rotation_quaternion(rotation), "zoom": zoom}


def _mesh_points(layers):
    meshes = [layer.data["vertices"] for layer in layers if layer.kind == "mesh"]
    return np.concatenate(meshes) if meshes else None


def _pdb_coords(text):
    rows = [line for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))]
    return np.asarray([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in rows]) \
        if rows else np.empty((0, 3))


def _lining_layers(rows, colors):
    groups = {}
    for row in rows:
        group = groups.setdefault(row["segment"], {})
        group.setdefault(row.get("viewer_chain") or "", []).append(int(row["viewer_resid"]))
    layers = []
    for segment, by_chain in groups.items():
        kind = segment.split("_")[0]
        residues = [{"chain": chain, "resi": sorted(set(resids))} for chain, resids in by_chain.items()]
        count = sum(len(r["resi"]) for r in residues)
        layers.append(Layer(f"lining_{segment}", f"{_segment_label(segment)} lining ({count})", "sticks",
                            colors.get(segment, _SEGMENT_HEX.get(kind, NON_CHANNEL_HEX)), visible=False,
                            data={"residues": residues},
                            description="residues with an atom centre within the lining cutoff of the segment's "
                                        "measured grid nodes"))
    return layers


def _from_cast_manifest(directory, manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    files = manifest.get("files", manifest)
    structure_path = _resolve(files.get("structure_file"), directory)
    structure = structure_path.read_text() if structure_path else None
    scene_path = _resolve(files.get("scene_json"), directory)
    scene = json.loads(scene_path.read_text()) if scene_path else {}
    layers = [Layer("protein", "Protein", "cartoon", PROTEIN_HEX, opacity=PROTEIN_OPACITY)] if structure else []
    colors = {}
    segments = scene.get("cast_segments") or []
    if segments:
        order = sorted(segments, key=lambda r: ({"entry": 0, "lumen": 1, "exit": 2}.get(r["segment_type"], 3),
                                                r["segment"]))
        for row in order:
            dx = _resolve(files.get(f"segment_{row['segment']}_display_dx"), directory)
            if dx is None:
                continue
            mesh = field_mesh(*read_dx_field(dx))
            if mesh is None:
                continue
            color = row.get("color_hex") or _hex(row["color_rgb"])
            colors[row["segment"]] = color
            layers.append(_mesh_layer(row["segment"], _segment_label(row["segment"]), color, mesh,
                                      description=f"{row['segment_type']} segment display surface; measured "
                                                  f"volume {row.get('measured_volume_A3', float('nan')):.1f} Å³"))
    else:
        dx = _resolve(files.get("display_dx"), directory)
        if dx is not None:
            mesh = field_mesh(*read_dx_field(dx))
            if mesh is not None:
                color = _hex(scene["cast_color_rgb"]) if scene.get("cast_color_rgb") else NON_CHANNEL_HEX
                channel = color.lower() == LUMEN_HEX
                colors["lumen"] = color
                layers.append(_mesh_layer("cast", "Channel cast" if channel else "Cast", color, mesh,
                                          description="display surface of the measured cast"))
    lining_csv = _resolve(files.get("lining_residues_csv"), directory)
    if lining_csv:
        with open(lining_csv, newline="") as handle:
            layers += _lining_layers(list(csv.DictReader(handle)), colors)
    profile_path = _resolve(files.get("profile_json"), directory)
    paths = _profile_paths(json.loads(profile_path.read_text())) if profile_path else []
    if not any(len(p) for p in paths[1:]) and scene.get("exit_paths"):
        paths += [p.get("points") if isinstance(p, dict) else p for p in scene["exit_paths"]]
        paths = [p for p in paths if p and len(p) > 1]
    if paths:
        layers.append(Layer("centre_line", "Profile centre line", "line", CENTRE_LINE_HEX, visible=False,
                            data={"paths": paths}, description="axial profile samples and lateral exit legs"))
    camera = _camera(scene.get("camera_rotation"), scene.get("center"), scene.get("vertical_span"),
                     _pdb_coords(structure) if structure else None, focus=_mesh_points(layers))
    title = f"{manifest.get('prefix', Path(directory).name)} - CREVICE view"
    return CreviceView(layers, structure=structure, title=title, camera=camera, source=str(manifest_path),
                       notes=["Surfaces are display geometry; measured volumes are in the bundle's tables."])


def _from_analysis_scene(directory, scene_path, snapshot=0):
    scene = json.loads(Path(scene_path).read_text())
    root = Path(directory)
    snapshots = scene.get("snapshots") or []
    index = int(snapshot)
    if snapshots and not -len(snapshots) <= index < len(snapshots):
        raise ValueError(f"snapshot {index} out of range: {len(snapshots)} snapshots")
    chosen = snapshots[index] if snapshots else None
    protein = root / (chosen["pdb"] if chosen else scene["protein_pdb"])
    structure = protein.read_text() if protein.is_file() else None
    layers = [Layer("protein", "Protein", "cartoon", PROTEIN_HEX, opacity=PROTEIN_OPACITY)] if structure else []
    names = {"cavity": ("region", "Region cast"), "water": ("water_density", "Water density")}
    for surface in scene.get("surfaces") or ():
        dx = _resolve(surface.get("dx"), root)
        if dx is None:
            continue
        mesh = field_mesh(*read_dx_field(dx), level=float(surface.get("level", .5)))
        if mesh is None:
            continue
        name, label = names.get(surface["key"], (surface["key"], surface["key"]))
        layers.append(_mesh_layer(name, label, _hex(surface["color_rgb"]), mesh,
                                  opacity=_INTERACTIVE_OPACITY.get(surface["key"], float(surface.get("opacity", 1.))),
                                  description=surface.get("meaning", "")))
    residues = [r for r in scene.get("residues") or () if r.get("residue") in set(scene.get("display_residues") or ())
                or r.get("key") in set(scene.get("display_residues") or ())]
    sticks = {}
    for row in residues:
        chain, resid = row.get("viewer_chain"), row.get("viewer_resid")
        if resid is None:
            continue
        sticks.setdefault(chain or "", []).append(int(resid))
    if sticks:
        layers.append(Layer("lining_region", f"Region boundary residues ({sum(map(len, sticks.values()))})",
                            "sticks", "#7b8794", visible=False,
                            data={"residues": [{"chain": c, "resi": sorted(set(r))} for c, r in sticks.items()]},
                            description="direct boundary contributors of the reference region"))
    waters = scene.get("member_waters")
    if waters and chosen is not None:
        record = next((s for s in waters.get("snapshots") or ()
                       if s.get("source_frame_index") == chosen.get("source_frame_index")), None)
        if record and record.get("pdb") and (root / record["pdb"]).is_file():
            points = _pdb_coords((root / record["pdb"]).read_text())
            layers.append(Layer("waters", f"Member waters ({len(points)})", "spheres",
                                _hex(waters.get("color_rgb", (1., .5, 0.))), data={
                                    "points": points.tolist(), "radius": float(waters.get("sphere_radius_A", .7))},
                                description=waters.get("definition", "")))
    camera = scene.get("camera") or {}
    view = _camera(camera.get("camera_rotation"), camera.get("center"), camera.get("vertical_span"),
                   _pdb_coords(structure) if structure else None, focus=_mesh_points(layers))
    when = "" if chosen is None else f" frame {chosen['source_frame_index']}"
    return CreviceView(layers, structure=structure, title=f"{scene.get('prefix', root.name)}{when} - CREVICE view",
                       camera=view, source=str(scene_path),
                       notes=["Member waters are geometric membership (oxygen centre inside the measured region), "
                              "not binding."])


def _find_one(directory, pattern):
    found = sorted(Path(directory).glob(pattern))
    return found[0] if found else None


def view(source=None, *, frame=None, cast=None, segments=None, profile=None, lining=None, waters=None,
         kind="auto", snapshot=0, smooth=None, smooth_supersample=None, title=None):
    """Build an interactive view of CREVICE results for Jupyter or a browser.

    Parameters
    ----------
    source : str or Path, optional
        Output directory (or its manifest/scene JSON) of ``crevice publish``,
        ``crevice cast``, ``crevice region-trajectory``, ``crevice
        cavity-trajectory`` or ``crevice hydration``, or of ``crevice paths``
        (a ``*_paths.json``: one line layer per branch). Omit it to view
        in-memory results instead.
    frame : StructureFrame, optional
        In-memory protein (cartoon and lining sticks).
    cast : VoidCast, optional
        In-memory cast; drawn as the LUMEN when ``segments`` are given or the
        cast is a channel cast, otherwise as a violet cast.
    segments : sequence of VoidComponent, optional
        ENTRY/EXIT segment casts (:func:`crevice.cast_segments.segment_casts`).
    profile : PoreProfile, optional
        Channel profile for the centre-line layer.
    lining : sequence of dict, optional
        Rows of :func:`crevice.cast_segments.lining_residues` (``segment``,
        ``viewer_chain``, ``viewer_resid``) for the stick layers.
    waters : array_like, shape (n, 3), optional
        Member-water oxygen positions (Å) for the sphere layer.
    kind : {"auto", "cast", "region"}, default "auto"
        Which output to read from a directory: the cast/publish manifest or the
        hydration analysis scene. ``auto`` prefers the cast manifest.
    snapshot : int, default 0
        Hydration outputs: which saved snapshot (first, middle, last, ...) to
        show with its member waters.
    smooth : float, optional
        In-memory casts: display smoothing length (Å), as ``--smooth``.
    smooth_supersample : int, optional
        In-memory casts: display supersampling factor.
    title : str, optional
        Page title (not drawn into the scene).

    Returns
    -------
    CreviceView
        Call :meth:`~CreviceView.show` in Jupyter or
        :meth:`~CreviceView.save_html` anywhere.

    Raises
    ------
    FileNotFoundError
        When the directory holds no recognised CREVICE output.
    ValueError
        When neither ``source`` nor in-memory results are given, or when the
        directory is a data-only trajectory run (no scene was written; rerun with
        ``--outputs default,hydration`` or ``--all``).

    Colours and representation
    --------------------------
    As :class:`CreviceView`: grey-blue cartoon; ENTRY ``#f0d43a``, LUMEN
    ``#1494a1``, EXIT ``#c2410c`` and non-channel ``#a855f7`` surfaces; hidden
    lining sticks and centre line; orange member-water spheres; no text in the
    3D canvas.

    Examples
    --------
    >>> import crevice
    >>> v = crevice.view("3UKM_bundle")            # doctest: +SKIP
    >>> v.toggle("lining_lumen").show()            # doctest: +SKIP
    >>> v.save_html("3UKM_view.html")              # doctest: +SKIP
    """
    if source is not None:
        path = Path(source)
        if path.is_file():
            directory, target = path.parent, path
        elif path.is_dir():
            directory, target = path, None
        else:
            raise FileNotFoundError(f"no such output directory or file: {path}")
        if target is None:
            manifest = None
            if kind in ("auto", "cast"):
                manifest = next((p for p in sorted(directory.glob("*_manifest.json"))
                                 if not p.name.endswith(("_hydration_manifest.json", "_cavity_manifest.json"))), None)
            analysis = _find_one(directory, "*_analysis_scene.json") if kind in ("auto", "region") else None
            target = manifest or analysis
            if target is None and kind == "auto" and _find_one(directory, "*_paths.json") is not None:
                # crevice paths output: one line layer per branch (crevice.paths).
                from .paths import paths_view_from_json
                return paths_view_from_json(_find_one(directory, "*_paths.json"), title=title)
            if target is None and any(directory.glob("*_cavity_manifest.json")):
                raise ValueError(f"{directory}: this trajectory run wrote data only (the 0.3 default), so there "
                                 "is no scene to view; rerun it with --outputs default,hydration (or --all) to "
                                 "write the analysis scene")
            if target is None:
                raise FileNotFoundError(f"{directory}: no CREVICE cast/publish manifest (*_manifest.json) or "
                                        "hydration scene (*_analysis_scene.json)")
        if target.name.endswith("_paths.json"):
            from .paths import paths_view_from_json
            return paths_view_from_json(target, title=title)
        result = (_from_analysis_scene(directory, target, snapshot) if target.name.endswith("_analysis_scene.json")
                  else _from_cast_manifest(directory, target))
        if title:
            result.title = title
        return result
    if frame is None and cast is None:
        raise ValueError("give an output directory or in-memory results (frame, cast, segments, profile)")
    return _from_memory(frame, cast, segments, profile, lining, waters, smooth, smooth_supersample, title)


def _frame_pdb(frame):
    import tempfile
    from .presentation import write_viewer_structure
    with tempfile.TemporaryDirectory() as tmp:
        path = write_viewer_structure(frame, Path(tmp) / "frame.pdb")
        return path.read_text()


def _from_memory(frame, cast, segments, profile, lining, waters, smooth, smooth_supersample, title):
    from dataclasses import replace
    from .presentation import cast_rgb
    structure = _frame_pdb(frame) if frame is not None else None
    layers = [Layer("protein", "Protein", "cartoon", PROTEIN_HEX, opacity=PROTEIN_OPACITY)] if structure else []
    colors = {}
    if cast is not None:
        segments = tuple(segments or ())
        parts = []
        if segments:
            parts.append(("lumen", cast))
            parts += [(c.metadata["segment"], replace(cast, components=(c,), points=tuple(c.points)))
                      for c in segments]
        else:
            parts.append(("cast", cast))
        cast_color = _hex(cast_rgb(cast))
        for name, part in parts:
            if not part.points:
                continue
            mesh = cast_mesh(part, smooth=smooth, smooth_supersample=smooth_supersample)
            if mesh is None:
                continue
            kind = name.split("_")[0]
            color = cast_color if kind in ("lumen", "cast") else _SEGMENT_HEX.get(kind, NON_CHANNEL_HEX)
            colors["lumen" if kind == "cast" else name] = color
            label = _segment_label(name) if segments else ("Channel cast" if color == LUMEN_HEX else "Cast")
            layers.append(_mesh_layer(name, label, color, mesh))
        order = {"entry": 0, "lumen": 1, "exit": 2}
        meshes = sorted((layer for layer in layers if layer.kind == "mesh"),
                        key=lambda layer: (order.get(layer.name.split("_")[0], 3), layer.name))
        layers = [layer for layer in layers if layer.kind != "mesh"] + meshes
    if lining:
        layers += _lining_layers(lining, colors)
    paths = _profile_paths(profile)
    if paths:
        layers.append(Layer("centre_line", "Profile centre line", "line", CENTRE_LINE_HEX, visible=False,
                            data={"paths": paths}))
    if waters is not None and len(waters):
        layers.append(Layer("waters", f"Member waters ({len(waters)})", "spheres", WATER_HEX,
                            data={"points": np.asarray(waters, float).tolist(), "radius": WATER_RADIUS_A}))
    camera = None
    if frame is not None and cast is not None:
        from .presentation import scene_geometry
        try:
            geometry = scene_geometry(cast, frame, profile)
            camera = _camera(geometry.get("camera_rotation"), geometry.get("center"), geometry.get("vertical_span"),
                             np.asarray([a.coord for a in frame.atoms]), focus=_mesh_points(layers))
        except (ValueError, KeyError, TypeError):
            camera = None
    return CreviceView(layers, structure=structure, title=title or "CREVICE view", camera=camera, source="memory")

"""Display geometry and figure-text control, kept separate from measurements.

Two concerns live here:

* display geometry (cast extensions, viewer structures, camera geometry) that
  never changes a measured volume, radius or count;
* the single switch that decides whether CREVICE draws explanatory text
  (titles, captions, value call-outs, residue/landmark labels, viewer legends
  and 2D labels) into its figures and molecular scenes. Figures and scenes are
  clean by default; :func:`figure_annotations`, the ``annotate=True`` argument
  of every figure/scene writer, or the CLI ``--annotate`` flag turn the text on.
  Whatever the setting, the text is recorded with the output (PNG metadata, or
  the scene JSON ``annotations`` entry), so no information is lost.
"""
from __future__ import annotations

import contextlib
import contextvars
from dataclasses import replace
import functools
import json
import math
from pathlib import Path

from .geometry import add, scale, dot
from .models import StructureFrame
from .spatial import SpatialIndex


def extend_channel_cast(cast, frame, length=2.0):
    """Extrude terminal footprints, retaining only atom-clear attached columns.

    For each end of a connected channel cast (``export_mode ==
    "connected_section_fill"``), the samples in the terminal mouth plane are
    stepped outward along the axis one grid spacing at a time; a new sample is kept
    only if it and the segment to it clear every atom by at least the cast's
    minimum radius. This bounded continuation is a presentation aid, not a solvent
    volume or a measurement of the exterior lumen. The original cast is never
    modified, and the measured volume is recorded unchanged in the metadata.

    Parameters
    ----------
    cast : VoidCast
        Channel cast from the profile-guided fill.
    frame : StructureFrame
        Structure whose atoms clip the extension (the cast's own atom selection).
    length : float, default 2.0
        Requested extension beyond each mouth, Å (rounded down to whole grid
        steps). 0 returns ``cast`` unchanged.

    Returns
    -------
    VoidCast
        A copy with the extra samples, ``display_only=True``,
        ``extension_requested_A``, ``extension_actual_A`` (per end),
        ``extension_added_points`` and ``measured_volume_A3`` in its metadata; or
        ``cast`` itself for zero length, a non-channel cast or a focus-cropped cast
        (a crop boundary is not a channel mouth).

    Raises
    ------
    ValueError
        If ``length`` is negative or not finite.
    """
    if not math.isfinite(length) or length < 0:
        raise ValueError("cast_extension must be finite and non-negative")
    if not length or cast.metadata.get("export_mode") != "connected_section_fill":
        return cast
    if cast.metadata.get("focus_points"):
        return cast  # a crop boundary must not be advertised as a channel mouth
    spatial = SpatialIndex(frame.selected_atoms(
        include_hydrogen=cast.metadata.get("include_hydrogen", False),
        include_hetero=cast.metadata.get("include_hetero", True)))
    axis = cast.metadata["axis_direction"]
    planes = cast.metadata["mouth_planes_A"]
    steps = int(math.floor(length/cast.spacing+1e-9))
    components, added, extents = [], 0, [0.0, 0.0]
    for component in cast.components:
        points = list(component.points)
        for end, sign in enumerate((-1, 1)):
            terminal = [p for p in component.points if abs(p.t-planes[end]) < 1e-7]
            for k in range(1, steps+1):
                layer = []
                for point in terminal:
                    position = add(point.position, scale(axis, sign*cast.spacing))
                    raw, nearest = spatial.nearest_surface(position)
                    if raw >= cast.min_radius and spatial.segment_clearance(point.position, position) >= cast.min_radius:
                        layer.append(replace(point, index=len(points)+len(layer), position=position,
                                             t=point.t+sign*cast.spacing, radius=raw, raw_clearance=raw,
                                             nearest_atom_serial=nearest.serial, nearest_residue=nearest.residue_key.label))
                points.extend(layer)
                added += len(layer)
                terminal = layer
                if layer:
                    extents[end] = k*cast.spacing
        components.append(replace(component, points=tuple(points), volume=len(points)*cast.spacing**3))
    all_points = tuple(p for c in components for p in c.points)
    return replace(cast, components=tuple(components), points=all_points,
                   metadata={**cast.metadata, "display_only": True,
                             "extension_requested_A": length, "extension_actual_A": extents,
                             "extension_added_points": added,
                             "measured_volume_A3": cast.total_volume,
                             "extension_definition": "constant terminal footprint, clipped by atom-clear axial segments; excluded from measured volume"})


def orthogonal_map_cast(cast):
    """Same samples in a positive Cartesian lattice for PyMOL's DX reader.

    A cast measured on a rotated lattice (``metadata["grid_basis"]``) is
    re-expressed in that lattice's own coordinates, so its DX map has axis-aligned
    steps. Only used for display maps; positions are rotated, not resampled.

    Parameters
    ----------
    cast : VoidCast
        Cast whose metadata may hold ``grid_basis`` (rows = lattice axes).

    Returns
    -------
    VoidCast
        A copy with every sample position projected onto the basis rows and
        ``grid_basis`` set to the identity.
    """
    basis = cast.metadata.get("grid_basis", ((1,0,0),(0,1,0),(0,0,1)))
    def local(p):
        return replace(p, position=tuple(dot(p.position, row) for row in basis))
    return replace(cast, components=tuple(replace(c, points=tuple(local(p) for p in c.points))
                                          for c in cast.components),
                   points=tuple(local(p) for p in cast.points),
                   metadata={**cast.metadata, "grid_basis": ((1,0,0),(0,1,0),(0,0,1))})


def write_viewer_structure(frame, path):
    """Portable single-model PDB with explicit reversible identity remapping.

    Molecular viewers need classic PDB limits: one-character chains, residue
    numbers <= 9999 and atom serials <= 99999. Chains are kept when they are
    already single characters and otherwise mapped to unused characters;
    residues are renumbered 1, 2, ... per chain in input order; atoms are
    renumbered in order. Analysis retains original identities and coordinates.
    Viewer PDB precision is 0.001 Å; the sidecar records the original labels and
    renumbered atom serials.

    Parameters
    ----------
    frame : StructureFrame
        Structure to copy (all of its atoms).
    path : str or Path
        Output PDB path; the sidecar is written beside it with the suffix
        replaced by ``.identities.json``.

    Returns
    -------
    Path
        The written PDB.

    Outputs
    -------
    path : PDB
        ``ATOM``/``HETATM`` records, occupancy 1.00, B-factor 0.00, elements.
    <stem>.identities.json : JSON
        ``chain_mapping``, ``residue_mapping`` (original label to viewer chain
        and residue number), ``original_serials_in_viewer_order``, the source,
        model index, atom count and coordinate rounding.

    Raises
    ------
    ValueError
        Beyond 62 chains, 99999 atoms or 9999 residues per chain, or for
        coordinates that do not fit the PDB columns.
    """
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    chains = list(dict.fromkeys(a.chain_id for a in frame.atoms))
    if len(chains) > len(alphabet) or len(frame.atoms) > 99999:
        raise ValueError("Viewer PDB capacity exceeded (62 chains / 99999 atoms)")
    chain_map = {c: (c if len(c)==1 and c in alphabet else None) for c in chains}
    available = [c for c in alphabet if c not in chain_map.values()]
    for c in chains:
        if chain_map[c] is None:
            chain_map[c] = available.pop(0)
    residues = {}
    counts = {}
    for a in frame.atoms:
        if a.residue_key not in residues:
            counts[a.chain_id] = counts.get(a.chain_id, 0)+1
            residues[a.residue_key] = counts[a.chain_id]
    if max(counts.values()) > 9999:
        raise ValueError("Viewer PDB residue capacity exceeded")
    lines = ['REMARK CREVICE selected single-model viewer copy; identities in sidecar\n']
    for serial, a in enumerate(frame.atoms, 1):
        if any(not math.isfinite(v) or len(f'{v:8.3f}') > 8 for v in a.coord):
            raise ValueError("Coordinates exceed viewer PDB format")
        record = 'HETATM' if a.hetero else 'ATOM  '
        name = f' {a.name:<3}' if len(a.element)==1 and len(a.name)<4 else f'{a.name:<4}'
        lines.append(f'{record}{serial:5d} {name[:4]} {a.resname[:3]:>3} {chain_map[a.chain_id]}'
                     f'{residues[a.residue_key]:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}'
                     f'{1.0:6.2f}{0.0:6.2f}          {a.element:>2}\n')
    lines.append('END\n')
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(''.join(lines))
    mapping = {"source": frame.source, "model_index": frame.model_index,
               "atom_count": len(frame.atoms), "coordinate_rounding_A": 0.001,
               "chain_mapping": chain_map,
               "residue_mapping": [{"original": k.label, "viewer_chain": chain_map[k.chain_id],
                                     "viewer_resid": v} for k, v in residues.items()],
               "original_serials_in_viewer_order": [a.serial for a in frame.atoms]}
    output.with_suffix('.identities.json').write_text(json.dumps(mapping, indent=2)+'\n')
    return output


def scene_geometry(cast, frame=None, profile=None):
    """Shared camera and mouth rings in input coordinates; no geometry smoothing.

    The camera looks across the cast's third grid direction (for a rolling cast
    that is clearly elongated, its dominant direction instead), tilted by 16
    degrees, or 35 degrees for multimers (at least three chains) with at least
    three substantial regions (each >= 10% of the largest volume). The same
    camera is written into every viewer scene so PyMOL, VMD and ChimeraX show the
    same view.

    Parameters
    ----------
    cast : VoidCast
        Displayed cast.
    frame : StructureFrame, optional
        Structure, included in the framing when given.
    profile : PoreProfile, optional
        Channel profile; its mouths give the mouth-plane guide rings and its
        lateral exit legs (``metadata["exits"]``) give the exit centre lines.
        Ignored for focus-cropped casts.

    Returns
    -------
    dict
        ``camera_rotation`` (3 x 3 rows), ``camera_axis``, ``camera_tilt_degrees``,
        ``end_view_turn_degrees``, ``framing_margin``, ``center`` and
        ``vertical_span``/``radius`` (Å) for framing, ``mouth_rings`` (65-point
        circles 0.3 Å beyond the widest terminal sample, guides only, not
        measured apertures) and, when present, ``exit_paths`` (one centre line
        per lateral exit leg, drawn as blue tubes in the scenes).
    """
    import numpy as np
    basis = np.asarray(cast.metadata.get('grid_basis', np.eye(3)), dtype=float)
    u, v, axis = basis
    camera_axis_definition = 'third grid direction'
    if cast.mode == 'rolling':
        cloud = np.asarray([p.position for c in cast.components for p in c.points])
        _, singular, directions = np.linalg.svd(cloud-cloud.mean(0), full_matrices=False)
        if len(singular) >= 2 and singular[0] > 2.5*max(singular[1],1e-12):
            axis = directions[0]
            reference = basis[np.argmax(np.abs(basis@axis))]
            if axis@reference < 0:
                axis = -axis
            u = basis[np.argmin(np.abs(basis@axis))]
            u = u-(u@axis)*axis
            u /= np.linalg.norm(u)
            v = np.cross(axis,u)
            camera_axis_definition = 'dominant elongated void direction for display; not a channel assignment'
    camera = np.asarray([u, axis, np.cross(u, axis)])
    major_regions = sum(c.volume >= .1*max(x.volume for x in cast.components) for c in cast.components)
    multichain = frame is not None and len({a.chain_id for a in frame.atoms}) >= 3
    tilt_degrees = 35 if multichain and major_regions >= 3 else 16
    angle = math.radians(tilt_degrees)
    tilt = np.array([[1,0,0],[0,math.cos(angle),-math.sin(angle)],[0,math.sin(angle),math.cos(angle)]])
    camera = tilt @ camera
    points = [p.position for c in cast.components for p in c.points]
    if frame:
        points += [a.coord for a in frame.atoms]
    coords = np.asarray(points)
    low, high = coords.min(axis=0), coords.max(axis=0)
    center = (low+high)/2
    rings = []
    if profile and profile.metadata.get('channel_mouths') and not cast.metadata.get('focus_points'):
        for label, point in [('lower',profile.points[0]),('upper',profile.points[-1])]:
            t = profile.metadata['channel_mouths'][label]['coordinate_A']
            c = np.asarray(point.position)+(t-point.t)*np.asarray(profile.axis_direction)
            terminal = [p for comp in cast.components for p in comp.points
                        if abs(p.t-point.t) < cast.spacing+1e-6]
            radius = max([point.raw_clearance]+[float(np.linalg.norm(
                         (np.asarray(p.position)-c) @ np.array([u,v]).T)) for p in terminal])+0.3
            ring = [list(c + radius*(math.cos(a)*u+math.sin(a)*v))
                    for a in np.linspace(0,2*math.pi,65)]
            rings.append({'end':label, 'coordinate_A':t, 'points':ring,
                          'definition':'mouth-plane guide, not a measured aperture contour'})
    # Lateral exit legs of capped channel ends: centre lines only, one entry
    # per distinct leg (every leg of every capped end, widest first).
    exits = []
    if profile and not cast.metadata.get('focus_points'):
        for label, record in (profile.metadata.get('exits') or {}).items():
            if record.get('type') == 'lateral':
                for leg in record.get('legs', ()):
                    exits.append({'end': label, 'leg': leg.get('leg', 1),
                                  'points': [list(map(float, p)) for p in leg['vertices']],
                                  'definition': 'lateral exit leg centre line (capped end); not a measured surface'})
    view_span = np.ptp((coords-center) @ camera.T, axis=0)
    vertical_span = float(max(view_span[1], view_span[0]*.75)+4)*1.10
    return {'camera_rotation':camera.tolist(), 'camera_axis':axis.tolist(),
            'camera_tilt_degrees':tilt_degrees, 'end_view_turn_degrees':90-tilt_degrees,
            'framing_margin':1.10, 'overview_definition':'35-degree elevation for at least three substantial regions in a multimer; otherwise 16 degrees',
            'camera_axis_definition':camera_axis_definition, 'center':center.tolist(),
            'vertical_span':vertical_span,
            'radius':float(np.linalg.norm(high-low)/2+1), 'mouth_rings':rings,
            **({'exit_paths':exits} if exits else {})}


# ---------------------------------------------------------------------------
# Default cast colours and opt-in display guides
# ---------------------------------------------------------------------------

#: Channel cast surface colour (teal), RGB in 0-1; also ``#1494a1``.
CHANNEL_CAST_RGB = (0.08, 0.58, 0.63)

#: Non-channel cast colour (violet ``#a855f7``), RGB in 0-1. Used for every
#: cast that is not a resolved channel lumen (cavity, rolling-probe and
#: dominant-region casts, the unresolved-profile fallback) and for the lines of
#: their width-profile plots. Chosen to stay separable, including under
#: simulated protanopia and deuteranopia (OKLab ΔE×100 >= 10), from the channel
#: teal and profile blue, the exit-leg green ``#009e73``, the orange member-water
#: and interval colours and both ends of the red-to-blue hydration scale.
NON_CHANNEL_CAST_RGB = (0.659, 0.333, 0.969)

#: :data:`NON_CHANNEL_CAST_RGB` as a hex string for Matplotlib and ChimeraX.
NON_CHANNEL_CAST_HEX = "#a855f7"

_MOUTH_GUIDES = contextvars.ContextVar("crevice_mouth_guides", default=False)
_EXIT_CENTRE_LINES = contextvars.ContextVar("crevice_exit_centre_lines", default=False)


def cast_rgb(cast) -> tuple[float, float, float]:
    """Default surface colour of a cast: teal for channels, violet otherwise.

    Parameters
    ----------
    cast : VoidCast
        Any cast; ``cast.mode == "channel"`` is a resolved channel lumen.

    Returns
    -------
    tuple of float
        :data:`CHANNEL_CAST_RGB` or :data:`NON_CHANNEL_CAST_RGB`.

    Examples
    --------
    >>> from crevice.models import VoidCast
    >>> cast_rgb(VoidCast((), (), 0.5, 0.0, "rolling"))
    (0.659, 0.333, 0.969)
    """
    return CHANNEL_CAST_RGB if getattr(cast, "mode", None) == "channel" else NON_CHANNEL_CAST_RGB


def mouth_guides_enabled() -> bool:
    """Return whether channel-mouth guides are drawn in scenes and profile plots.

    Off by default: no mouth rings in the PyMOL, VMD and ChimeraX scenes and no
    orange mouth brackets or mouth lines on profile plots. Mouth positions are
    always written to the profile JSON, scene JSON and manifests.

    Returns
    -------
    bool
        ``False`` unless enabled by :func:`display_guides`, a writer's
        ``mouth_guides=True`` or the CLI ``--mouth-guides`` flag.
    """
    return _MOUTH_GUIDES.get()


def exit_centre_lines_enabled() -> bool:
    """Return whether lateral exit legs are also drawn as thin centre lines.

    Off by default: lateral exit legs are shown as cast surfaces
    (:func:`crevice.channel_exits.lateral_exit_casts`); the thin centre-line
    tubes are an option. The legs' vertices are always in the scene JSON.

    Returns
    -------
    bool
        ``False`` unless enabled by :func:`display_guides`, a writer's
        ``exit_centre_lines=True`` or the CLI ``--exit-centre-lines`` flag.
    """
    return _EXIT_CENTRE_LINES.get()


@contextlib.contextmanager
def display_guides(*, mouth_guides: bool | None = None, exit_centre_lines: bool | None = None):
    """Turn the optional display guides on or off inside a ``with`` block.

    Parameters
    ----------
    mouth_guides : bool, optional
        Draw channel-mouth rings in viewer scenes and mouth brackets/lines on
        profile plots. ``None`` keeps the surrounding setting.
    exit_centre_lines : bool, optional
        Draw lateral exit legs as thin centre-line tubes as well as casts.
        ``None`` keeps the surrounding setting.

    Examples
    --------
    >>> with display_guides(mouth_guides=True):
    ...     mouth_guides_enabled()
    True
    >>> mouth_guides_enabled()
    False
    """
    tokens = []
    if mouth_guides is not None:
        tokens.append((_MOUTH_GUIDES, _MOUTH_GUIDES.set(bool(mouth_guides))))
    if exit_centre_lines is not None:
        tokens.append((_EXIT_CENTRE_LINES, _EXIT_CENTRE_LINES.set(bool(exit_centre_lines))))
    try:
        yield
    finally:
        for var, token in reversed(tokens):
            var.reset(token)


# ---------------------------------------------------------------------------
# Figure and scene annotation control
# ---------------------------------------------------------------------------

_ANNOTATE = contextvars.ContextVar("crevice_annotate", default=False)

#: PNG metadata key recording whether explanatory text was drawn.
ANNOTATION_METADATA_KEY = "CREVICE annotations"


def annotations_enabled() -> bool:
    """Return whether figures and scenes written now will draw explanatory text.

    Returns
    -------
    bool
        ``False`` unless enabled by :func:`figure_annotations`, an
        ``annotate=True`` argument of an enclosing writer, or the CLI
        ``--annotate`` flag.
    """
    return _ANNOTATE.get()


@contextlib.contextmanager
def figure_annotations(enabled: bool = True):
    """Turn drawn figure and scene text on or off inside a ``with`` block.

    CREVICE figures and molecular scenes are clean by default: no titles,
    suptitles, captions, value call-outs, residue or landmark labels, viewer
    legends/colour keys or 2D labels. Axis labels with units, tick labels,
    colour bars of quantitative colour scales and minimal legends for panels
    with several data series are always drawn. Inside this context every figure
    and scene writer draws the explanatory text as well. Nested writers that are
    given an explicit ``annotate`` argument override the context for their own
    output.

    Parameters
    ----------
    enabled : bool, default True
        ``True`` draws the explanatory text; ``False`` forces clean output.

    Examples
    --------
    >>> from crevice.presentation import figure_annotations
    >>> with figure_annotations():          # doctest: +SKIP
    ...     write_static_publication_bundle(frame, structure_path=path, output_dir="out")
    """
    token = _ANNOTATE.set(bool(enabled))
    try:
        yield
    finally:
        _ANNOTATE.reset(token)


def annotate_option(func):
    """Decorator applying a writer's keyword-only ``annotate`` argument.

    ``annotate=None`` (the default of every decorated writer) inherits the
    surrounding setting; ``True``/``False`` applies to that call and every
    figure or scene it writes.

    Parameters
    ----------
    func : callable
        Writer with a keyword-only ``annotate`` parameter.

    Returns
    -------
    callable
        The wrapped writer; an explicit ``annotate`` runs it inside
        :func:`figure_annotations`.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        annotate = kwargs.get("annotate")
        if annotate is None:
            return func(*args, **kwargs)
        with figure_annotations(annotate):
            return func(*args, **kwargs)
    return wrapper


def _one_line(text) -> str:
    return " ".join(str(text).split())


class FigureText:
    """Explanatory text of one figure: always recorded, drawn only when annotating.

    Whether text is drawn is fixed when the object is created, from
    :func:`annotations_enabled`. Every title, note, label and axis detail passed
    through it is recorded and embedded in the saved file's metadata (see
    :meth:`metadata` and :func:`figure_description`).

    Parameters
    ----------
    shows : str
        One-sentence description of what the figure shows. It is written to the
        output's metadata (``Description``) whether or not text is drawn.

    Attributes
    ----------
    draw : bool
        Whether explanatory text is drawn.
    shows : str
        The description, on one line.
    records : list of str
        Recorded text, as ``"Title: ..."``, ``"Note: ..."``, ``"Label: ..."``,
        ``"Notice: ..."`` or ``"Axis detail: ..."``.

    Notes
    -----
    Methods mirror the Matplotlib calls they replace. ``notice`` is the only
    method that always draws: it states that a panel has no data at all, because
    an empty axis is otherwise ambiguous.

    Examples
    --------
    >>> from crevice.presentation import FigureText, figure_annotations
    >>> text = FigureText("Radius against position")
    >>> text.note("probe 0.8 A")
    >>> text.draw, text.metadata()["Description"].splitlines()
    (False, ['Radius against position', 'Note: probe 0.8 A'])
    >>> with figure_annotations():
    ...     FigureText("x").draw
    True
    """

    def __init__(self, shows: str):
        self.draw = annotations_enabled()
        self.shows = _one_line(shows)
        self.records: list[str] = []

    def _record(self, kind: str, text) -> None:
        text = _one_line(text)
        if text:
            self.records.append(f"{kind}: {text}")

    def title(self, ax, text, **kwargs) -> None:
        """Record a panel title; draw it only when annotating.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Panel.
        text : str
            Title text.
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.set_title`.
        """
        self._record("Title", text)
        if self.draw:
            ax.set_title(text, **kwargs)

    def suptitle(self, fig, text, **kwargs) -> None:
        """Record a figure title; draw it only when annotating.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure.
        text : str
            Title text.
        **kwargs
            Passed to :meth:`matplotlib.figure.Figure.suptitle`.
        """
        self._record("Title", text)
        if self.draw:
            fig.suptitle(text, **kwargs)

    def text(self, ax, x, y, text, **kwargs) -> None:
        """Record an in-plot note; draw it only when annotating.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Panel.
        x, y : float
            Position in the coordinates given by ``kwargs`` (data by default).
        text : str
            Note text.
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.text`.
        """
        self._record("Note", text)
        if self.draw:
            ax.text(x, y, text, **kwargs)

    def annotate(self, ax, text, xy, **kwargs) -> None:
        """Record a value call-out or point label; draw it only when annotating.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Panel.
        text : str
            Label text.
        xy : tuple of float
            Annotated point.
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.annotate`.
        """
        self._record("Label", text)
        if self.draw:
            ax.annotate(text, xy, **kwargs)

    def note(self, text) -> None:
        """Record a note that is never drawn.

        Parameters
        ----------
        text : str
            Note text, stored in the figure metadata only.
        """
        self._record("Note", text)

    def notice(self, ax, x, y, text, **kwargs) -> None:
        """Draw and record a no-data notice (always drawn).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Empty panel.
        x, y : float
            Position of the notice.
        text : str
            Notice, for example "No network nodes detected".
        **kwargs
            Passed to :meth:`matplotlib.axes.Axes.text`.
        """
        self._record("Notice", text)
        ax.text(x, y, text, **kwargs)

    def axis_label(self, setter, label: str, detail: str | None = None, **kwargs) -> None:
        """Set an axis label; an explanatory second line is drawn only when annotating.

        The label itself (with its unit) is always drawn.

        Parameters
        ----------
        setter : callable
            Label setter, for example ``ax.set_xlabel``.
        label : str
            Axis label with unit, for example ``"Radius (Å)"``.
        detail : str, optional
            Explanatory second line; recorded always, drawn only when annotating.
        **kwargs
            Passed to ``setter``.
        """
        if detail:
            self._record("Axis detail", f"{label} - {detail}")
        setter(label + ("\n" + detail if detail and self.draw else ""), **kwargs)

    def metadata(self) -> dict[str, str]:
        """Metadata fields describing the figure, for embedding in the file.

        Returns
        -------
        dict of str to str
            ``Title`` (the first recorded title, else the description),
            ``Description`` (the description followed by every record, one per
            line) and ``CREVICE annotations`` (``drawn`` or ``omitted ...``).
        """
        titles = [r.split(": ", 1)[1] for r in self.records if r.startswith("Title: ")]
        return {"Title": titles[0] if titles else self.shows,
                "Description": "\n".join([self.shows, *self.records]),
                ANNOTATION_METADATA_KEY: "drawn" if self.draw else
                    "omitted (recorded here; draw with annotate=True or --annotate)"}


def figure_description(path) -> dict[str, str]:
    """Read the description CREVICE embeds in a PNG figure.

    Parameters
    ----------
    path : str or Path
        A PNG written by CREVICE.

    Returns
    -------
    dict of str to str
        The PNG text chunks: ``Title``, ``Description`` (what the figure shows,
        followed by every title, note, label and axis detail, one per line),
        ``CREVICE annotations`` (``drawn`` or ``omitted ...``) and ``Software``.
    """
    from PIL import Image
    with Image.open(path) as image:
        return {str(k): str(v) for k, v in getattr(image, "text", {}).items()}


def scene_annotation_record(texts: dict[str, str], *, drawn: bool | None = None) -> dict:
    """Describe the text of a molecular scene for its JSON definition.

    Parameters
    ----------
    texts : dict of str to str
        Each text element of the scene (legend caption, colour-scale ticks,
        status banner, residue labels) keyed by role.
    drawn : bool, optional
        Whether the scripts draw the text; defaults to :func:`annotations_enabled`.

    Returns
    -------
    dict
        ``{"drawn": bool, "texts": texts, "how_to_draw": ...}``.
    """
    drawn = annotations_enabled() if drawn is None else bool(drawn)
    return {"drawn": drawn, "texts": {k: v for k, v in texts.items() if v},
            "how_to_draw": "Text is omitted from the viewer scripts by default; rerun with --annotate (CLI) or annotate=True (Python) to draw it."}

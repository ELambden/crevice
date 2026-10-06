"""Publication-ready CREVICE figures and visualization bundles."""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from dataclasses import replace
from typing import Sequence

from .analysis import annotate_residues, detect_cavities, pore_profile
from .tunnels import find_tunnels
from .voids import void_cast
from .io import (
    write_cavities_csv,
    write_cavities_json,
    write_cavities_pdb,
    write_connectivity_csv,
    write_network_csv,
    write_profile_csv,
    write_profile_json,
    write_profile_pdb,
    write_residue_contacts_csv,
    write_tunnel_points_csv,
    write_tunnels_csv,
    write_tunnels_json,
    write_tunnels_pdb,
    write_void_cast_csv,
    write_void_cast_json,
    write_void_cast_pdb,
)
from .models import Cavity, FrameAnalysis, PoreProfile, ResidueContact, ResidueInteractionNetwork, StructureFrame, Tunnel, TrajectoryAnalysis
from .networks import build_cavity_network, network_metrics
from .presentation import (FigureText, annotate_option, annotations_enabled, figure_annotations, figure_description,
                           mouth_guides_enabled)
from .radii import RadiusSet, radii_fields, radii_option

CREVICE_BLUE = "#1f77b4"
CREVICE_TEAL = "#008b8b"
CREVICE_RED = "#b22222"
CREVICE_ORANGE = "#d9822b"
CREVICE_GRAY = "#4d4d4d"
CREVICE_GREEN = "#2e7d32"
CREVICE_PURPLE = "#6f4aa8"
CREVICE_LIGHT = "#eef3f7"


def _draw_channel_mouths(ax, profile, text: FigureText | None = None):
    mouths = profile.metadata.get("channel_mouths")
    if not mouths:
        return
    context = profile.metadata.get("display_context", [])
    exits = profile.metadata.get("exits") or {}
    labelled = False
    for label, segment in zip(("lower", "upper"), context):
        # At a capped end the straight continuation runs into the cap; the
        # lateral exit leg (``_draw_exit_legs``) is drawn there instead.
        if (exits.get(label) or {}).get("type") == "lateral":
            continue
        ax.plot([p["t"] for p in segment], [p["radius"] for p in segment],
                color="#8d969e", linestyle="--", linewidth=1.5,
                label=None if labelled else "Exterior clearance")
        labelled = True
    for label in ("lower", "upper"):
        mouth = mouths[label]
        t = mouth["coordinate_A"]
        if not mouth_guides_enabled():
            # Mouth brackets are opt-in (--mouth-guides); positions stay in metadata.
            if text is not None:
                text.note(f"{label.title()} mouth at {t:.2f} Å (bracket {mouth['bracket_A'][0]:.2f} to "
                          f"{mouth['bracket_A'][1]:.2f} Å); mouth guides not drawn")
            continue
        ax.axvspan(*mouth["bracket_A"], color=CREVICE_ORANGE, alpha=0.25)
        ax.axvline(t, color=CREVICE_ORANGE, linewidth=1.2, linestyle=":")
        message = f"{label.title()} mouth\n{t:.2f} Å"
        kwargs = dict(transform=ax.get_xaxis_transform(), ha="right" if label == "lower" else "left",
                      va="top", fontsize=8, color="#996018")
        if text is None:
            ax.text(t, 0.98, message, **kwargs)
        else:
            text.text(ax, t, 0.98, message, **kwargs)
    lower, upper = (mouths[k]["coordinate_A"] for k in ("lower", "upper"))
    left, right = ax.get_xlim()
    ax.axvspan(left, lower, color="#e6e9ec", alpha=0.45, zorder=0)
    ax.axvspan(upper, right, color="#e6e9ec", alpha=0.45, zorder=0)
    ax.set_xlim(left, right)


CREVICE_EXIT = "#009e73"


def _draw_exit_legs(ax, profile, text: FigureText | None = None) -> int:
    """Draw every lateral exit leg as a continuation beyond its capped mouth.

    Each leg of ``profile.metadata["exits"][end]["legs"]`` (present only for
    ``lateral_exits`` profiles with a capped end) is drawn against
    x = axial coordinate of the profile sample at that end ± the leg's arc
    length (Å, + beyond the upper end, − beyond the lower end), y = the leg's
    free-sphere radius (``radii_A``). The x values beyond the mouth are
    therefore path length along the leg, not an axial coordinate; the legend
    says so. Bluish-green (``#009e73``) dash-dot lines, one per leg; the
    legend entry is added once. Returns the number of legs drawn.
    """
    exits = profile.metadata.get("exits") or {}
    drawn = 0
    for label, sign, point in (("lower", -1, profile.points[0]), ("upper", 1, profile.points[-1])):
        record = exits.get(label) or {}
        if record.get("type") != "lateral":
            continue
        for leg in record.get("legs", ()):
            arc = leg.get("arc_lengths_A") or []
            radii = leg.get("radii_A") or []
            if len(arc) != len(radii) or len(arc) < 2:
                continue
            ax.plot([point.t + sign*a for a in arc], radii, color=CREVICE_EXIT, linestyle="-.",
                    linewidth=1.8, label=None if drawn else "Lateral exit (path length beyond mouth)")
            drawn += 1
            if text is not None:
                text.note(f"Lateral exit leg {leg.get('leg', 1)} at the {label} end: bottleneck "
                          f"{leg['bottleneck_radius_A']:.3f} Å, length {leg['length_A']:.2f} Å; drawn against "
                          "path length beyond the mouth")
    return drawn


def _set_title(ax, text: FigureText, title: str | None, default: str, **kwargs) -> None:
    """A caller-supplied title is the caller's own text and is always drawn."""
    if title:
        text._record("Title", title)
        ax.set_title(title, **kwargs)
    else:
        text.title(ax, default, **kwargs)


def _profile_summary(profile: PoreProfile) -> str:
    b = profile.bottleneck
    return (f"minimum radius {profile.min_radius:.3f} Å at the bottleneck (sample {b.index}, "
            f"t = {b.t:.2f} Å, nearest residue {b.nearest_residue or 'unassigned'}); "
            f"mean radius {profile.mean_radius:.3f} Å over {len(profile.points)} samples")


@annotate_option
def plot_profile_radius(
    profile: PoreProfile,
    path: str | Path,
    *,
    dpi: int = 400,
    title: str | None = None,
    show_bottleneck: bool = True,
    annotate: bool | None = None,
) -> Path:
    """Plot pore radius against the axial coordinate of a resolved channel.

    Each profile sample is drawn at its axial coordinate *t* (Å, along the
    profile axis from the profile origin) with its measured radius: the largest
    atom-surface clearance in that cross-section after any probe-radius
    subtraction (see :func:`crevice.channels.pore_profile`). The minimum-radius
    sample is marked as the bottleneck. When the profile carries channel-mouth
    metadata (``publish`` profiles), the mouths and the display-only exterior
    clearance beyond them are also drawn. For a capped channel resolved with
    ``lateral_exits``, every lateral exit leg is drawn as a continuation
    beyond its mouth against path length along the leg (see Colours). Nothing
    is recalculated; the figure only displays ``profile``.

    Parameters
    ----------
    profile : PoreProfile
        A resolved channel profile.
    path : str or Path
        Output image; parent directories are created. The suffix selects the
        format (PNG recommended; PNG, PDF and SVG carry the description metadata).
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title is
        drawn unless annotating, when "CREVICE pore radius profile" is used.
    show_bottleneck : bool, default True
        Mark the minimum-radius sample.
    annotate : bool, optional
        Also draw the default title, the bottleneck radius call-out ("1.33 Å")
        and (with mouth guides) the mouth labels. ``None`` (default) inherits the
        surrounding setting, which is off unless enabled with
        :func:`crevice.presentation.figure_annotations` or the CLI ``--annotate`` flag.

    Returns
    -------
    Path
        The written image.

    Outputs
    -------
    path : image
        6.8 x 4.2 inch figure. Axes: "Axial coordinate (Å)" and "Pore radius (Å)",
        radius axis starting at zero. The PNG ``Description`` metadata records the
        minimum and mean radius, bottleneck sample and nearest residue, and every
        omitted title or label.

    Colours and representation
    --------------------------
    * Blue line (``#1f77b4``) through the samples, without per-sample
      markers, and a light-blue fill to zero: pore radius.
    * Red dashed vertical line and red dot: bottleneck (minimum radius).
    * Light-grey shading outside the channel mouths. Orange shaded bands and
      dotted orange lines (lower and upper mouth brackets and positions) only
      with mouth guides (:func:`crevice.presentation.display_guides`, CLI
      ``--mouth-guides``); off by default, when the mouth positions are listed
      in the image metadata instead.
    * Grey dashed line: exterior clearance straight beyond each axial (open)
      mouth (display only, not part of the channel statistics); omitted at a
      capped end.
    * Bluish-green (``#009e73``) dash-dot line, one per lateral exit leg:
      the leg's free-sphere radius (distance to the nearest atom surface minus
      the probe radius) against the capped mouth's axial coordinate plus
      (upper end) or minus (lower end) the path length along the leg. It is a
      tunnel-like clearance along a 3D path, not a planar section radius,
      and is not part of ``min_radius``. Leg bottlenecks and lengths are
      listed in the image metadata.
    * Light-grey grid; a legend names the pore radius, bottleneck, exterior
      clearance and lateral exit lines that are present.
    """

    plt = _pyplot()
    output = _ensure_parent(path)
    x = [point.t for point in profile.points]
    y = [point.radius for point in profile.points]
    text = FigureText("Pore radius against axial coordinate for one resolved channel profile; " + _profile_summary(profile))

    fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
    ax.plot(x, y, color=CREVICE_BLUE, linewidth=2.2, label="Pore radius")
    ax.fill_between(x, y, [0.0] * len(y), color=CREVICE_BLUE, alpha=0.14, linewidth=0)
    _draw_exit_legs(ax, profile, text)
    if show_bottleneck:
        bottleneck = profile.bottleneck
        ax.axvline(bottleneck.t, color=CREVICE_RED, linestyle="--", linewidth=1.25, label="Bottleneck")
        ax.scatter([bottleneck.t], [bottleneck.radius], color=CREVICE_RED, s=34, zorder=4)
        text.annotate(
            ax,
            f"{bottleneck.radius:.2f} Å",
            (bottleneck.t, bottleneck.radius),
            xytext=(6, 10),
            textcoords="offset points",
            fontsize=8.5,
            color=CREVICE_RED,
        )
    _draw_channel_mouths(ax, profile, text)
    ax.set_xlabel("Axial coordinate (Å)")
    ax.set_ylabel("Pore radius (Å)")
    _set_title(ax, text, title, "CREVICE pore radius profile")
    ax.set_ylim(bottom=0.0)
    ax.grid(True, color="#d8d8d8", linewidth=0.65, alpha=0.75)
    ax.legend(frameon=False, loc="best")  # several line series: a minimal legend is kept
    _style_axes(ax)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


@annotate_option
def plot_profile_radius_with_residues(
    profile: PoreProfile,
    path: str | Path,
    *,
    dpi: int = 400,
    title: str | None = None,
    label_every: int = 5,
    max_labels: int = 28,
    annotate: bool | None = None,
) -> Path:
    """Plot the pore-radius profile above a lane of nearest-residue landmarks.

    The upper panel is the radius profile of :func:`plot_profile_radius` (the
    bottleneck is always marked). The lower lane places one landmark per
    selected profile sample at its axial coordinate: every ``label_every``-th
    sample whose nearest pore-lining residue differs from the previous
    landmark, plus the bottleneck sample, thinned evenly to ``max_labels``
    (first, last and bottleneck kept). The nearest residue is the one owning the
    atom closest to the sample centre (see
    :func:`crevice.analysis.annotate_residues`). Residue names are drawn only
    when annotating; they are always listed in the file metadata.

    Parameters
    ----------
    profile : PoreProfile
        A resolved channel profile with nearest-residue assignments.
    path : str or Path
        Output image; parent directories are created.
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title is
        drawn unless annotating ("CREVICE pore radius profile with residue landmarks").
    label_every : int, default 5
        Landmark spacing in profile samples (CLI ``--residue-label-every``);
        must be at least 1.
    max_labels : int, default 28
        Maximum number of landmarks.
    annotate : bool, optional
        Also draw the default and lane titles, the bottleneck call-out, the mouth
        labels and the residue names (for example ``A:TRP11``) with
        connector lines. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Raises
    ------
    ValueError
        If ``label_every`` is less than 1.

    Outputs
    -------
    path : image
        7.4 x 5.2 inch figure: radius panel ("Pore radius (Å)") over the landmark
        lane ("Reaction coordinate (Å)"). The PNG ``Description`` metadata lists
        every landmark as residue and axial coordinate.

    Colours and representation
    --------------------------
    * Radius panel: as :func:`plot_profile_radius` (blue profile line without
      sample markers, red bottleneck line and dot, orange mouth brackets only
      with mouth guides, grey
      dashed exterior clearance at open ends, bluish-green dash-dot lateral
      exit legs against path length beyond a capped mouth).
    * Landmark dots on a light-grey baseline, coloured by the nearest residue's
      chemistry: purple (``#6f4aa8``) charged (Asp, Glu, Lys, Arg, His), teal
      (``#008b8b``) polar (Ser, Thr, Asn, Gln, Tyr, Cys), blue (``#1f77b4``) other.
    * When annotating: grey residue names at alternating heights, rotated 34°,
      joined to their dots by thin grey connectors.
    """

    if label_every < 1:
        raise ValueError("label_every must be at least 1")
    plt = _pyplot()
    output = _ensure_parent(path)
    x = [point.t for point in profile.points]
    y = [point.radius for point in profile.points]
    labels = _profile_residue_landmarks(profile, label_every=label_every, max_labels=max_labels)
    text = FigureText("Pore radius against axial coordinate with a lane of nearest pore-lining residue landmarks; "
                      + _profile_summary(profile))
    text.note("Residue landmarks (residue at axial coordinate): "
              + "; ".join(f"{_short_residue_label(label)} at {point.t:.2f} Å" for point, label in labels))

    fig, (ax, lane) = plt.subplots(
        2,
        1,
        figsize=(7.4, 5.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.05], "hspace": 0.08},
        constrained_layout=True,
    )
    ax.plot(x, y, color=CREVICE_BLUE, linewidth=2.2, label="Pore radius")
    ax.fill_between(x, y, [0.0] * len(y), color=CREVICE_BLUE, alpha=0.13, linewidth=0)
    _draw_exit_legs(ax, profile, text)
    bottleneck = profile.bottleneck
    ax.axvline(bottleneck.t, color=CREVICE_RED, linestyle="--", linewidth=1.2, label="Bottleneck")
    ax.scatter([bottleneck.t], [bottleneck.radius], color=CREVICE_RED, s=38, zorder=4)
    text.annotate(
        ax,
        f"Bottleneck\n{bottleneck.radius:.2f} Å",
        (bottleneck.t, bottleneck.radius),
        xytext=(8, 12),
        textcoords="offset points",
        fontsize=8.4,
        color=CREVICE_RED,
        ha="left",
        va="bottom",
    )
    _draw_channel_mouths(ax, profile, text)
    ax.set_ylabel("Pore radius (Å)")
    _set_title(ax, text, title, "CREVICE pore radius profile with residue landmarks")
    ax.set_ylim(bottom=0.0)
    ax.grid(True, color="#d8d8d8", linewidth=0.65, alpha=0.75)
    ax.legend(frameon=False, loc="best")  # several line series: a minimal legend is kept
    _style_axes(ax)

    lane.set_ylim(0.0, 1.0)
    lane.set_yticks([])
    lane.set_xlabel("Reaction coordinate (Å)")
    lane.axhline(0.62, color="#c6cbd1", linewidth=0.9)
    for rank, (point, label) in enumerate(labels):
        y_text = 0.78 if rank % 2 == 0 else 0.28
        if text.draw:
            lane.plot([point.t, point.t], [0.55, y_text - 0.08 if y_text > 0.5 else y_text + 0.08], color="#a9b2ba", linewidth=0.7)
        lane.scatter([point.t], [0.62], s=18, color=_landmark_color(label), zorder=3)
        if text.draw:
            lane.text(
                point.t,
                y_text,
                _short_residue_label(label),
                rotation=34,
                ha="left",
                va="center",
                fontsize=7.0,
                color=CREVICE_GRAY,
            )
    text.title(lane, "Nearest pore-lining residue landmarks", fontsize=8.8, color=CREVICE_GRAY, pad=2)
    _style_axes(lane)
    lane.spines["left"].set_visible(False)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


@annotate_option
def plot_residue_contacts(
    contacts: Sequence[ResidueContact],
    path: str | Path,
    *,
    top_n: int = 20,
    dpi: int = 400,
    title: str | None = None,
    annotate: bool | None = None,
) -> Path:
    """Rank residues contacting a profile or region by influence score.

    Contacts are sorted by decreasing influence score, then increasing minimum
    distance, and the first ``top_n`` are drawn as horizontal bars (highest at
    the top). The influence score and role come from
    :func:`crevice.analysis.annotate_residues`; the figure does not recompute them.

    Parameters
    ----------
    contacts : sequence of ResidueContact
        Residue contacts, for example from :func:`crevice.analysis.annotate_residues`.
    path : str or Path
        Output image; parent directories are created.
    top_n : int, default 20
        Number of residues shown.
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title unless
        annotating ("CREVICE residue contacts").
    annotate : bool, optional
        Also draw the default title. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Outputs
    -------
    path : image
        Horizontal bar chart, height scaled with the number of residues. Axes:
        "Influence score (higher = closer contact along more of the path)" and
        "Residue" (residue IDs as tick labels). With no contacts, the panel
        states "No residue contacts detected".

    Colours and representation
    --------------------------
    Bar colour encodes the contact role assigned by
    :func:`crevice.analysis.annotate_residues`, where the gap is the residue's
    smallest atom-surface distance to the region surface (profile sphere or
    region point) and "at the bottleneck" means a contact within two profile
    samples of the narrowest sample:

    * red (``#b22222``) ``bottleneck``: gap ≤ 1 Å at the bottleneck;
    * orange (``#d9822b``) ``bottleneck-nearby``: gap > 1 Å at the bottleneck;
    * blue (``#1f77b4``) ``lining``: gap ≤ 1 Å elsewhere;
    * grey (``#4d4d4d``) ``nearby``: gap > 1 Å elsewhere, within the contact cutoff.

    A legend below the axes names the roles present among the drawn bars
    (a categorical colour key, so it is drawn by default). White bar edges;
    light-grey vertical grid.
    """

    plt = _pyplot()
    output = _ensure_parent(path)
    ranked = sorted(contacts, key=lambda contact: (-contact.influence_score, contact.min_distance))[:top_n]
    text = FigureText(f"Top {len(ranked)} residues ranked by influence score; bar colour is the contact role "
                      "(gap = smallest residue atom-surface distance to the region surface; "
                      "bottleneck = within two profile samples of the narrowest sample)")
    height = max(3.6, 0.32 * max(len(ranked), 1) + 1.2)
    fig, ax = plt.subplots(figsize=(7.2, height), constrained_layout=True)
    if not ranked:
        _empty_panel(ax, "No residue contacts detected", text)
    else:
        from matplotlib.patches import Patch
        labels = [contact.residue for contact in reversed(ranked)]
        scores = [contact.influence_score for contact in reversed(ranked)]
        colors = [_role_color(contact.role) for contact in reversed(ranked)]
        ax.barh(labels, scores, color=colors, edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Influence score (higher = closer contact along more of the path)")
        ax.set_ylabel("Residue")
        _set_title(ax, text, title, "CREVICE residue contacts")
        ax.grid(True, axis="x", color="#d8d8d8", linewidth=0.65, alpha=0.75)
        present = {_role_key(contact.role) for contact in ranked}
        handles = [Patch(facecolor=_role_color(role), edgecolor="white", label=label)
                   for role, label in ROLE_LEGEND if role in present]
        fig.legend(handles=handles, title="Contact role", loc="outside lower center", ncol=2,
                   fontsize=7.5, title_fontsize=8, frameon=False)
    _style_axes(ax)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


@annotate_option
def plot_cavity_summary(
    cavities: Sequence[Cavity],
    path: str | Path,
    *,
    top_n: int = 10,
    dpi: int = 400,
    title: str | None = None,
    annotate: bool | None = None,
) -> Path:
    """Summarise detected cavities by volume and maximum radius.

    Cavities are sorted by their detection score (highest first) and the first
    ``top_n`` are drawn: volume as bars on the left axis and maximum clearance
    radius as a connected line on a second, right-hand axis. Values come from
    :func:`crevice.analysis.detect_cavities`.

    Parameters
    ----------
    cavities : sequence of Cavity
        Detected cavities.
    path : str or Path
        Output image; parent directories are created.
    top_n : int, default 10
        Number of cavities shown.
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title unless
        annotating ("CREVICE cavity summary").
    annotate : bool, optional
        Also draw the default title. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Outputs
    -------
    path : image
        7.0 x 4.4 inch figure. Axes: "Cavity" (IDs ``C1``, ``C2``, ...),
        "Volume (Å³)" (left) and "Max radius (Å)" (right). With no cavities the
        panel states "No cavities detected".

    Colours and representation
    --------------------------
    Teal bars (``#008b8b``): cavity volume. Red line with circular markers
    (``#b22222``): maximum radius, read on the right axis. Light-grey horizontal grid.
    """

    plt = _pyplot()
    output = _ensure_parent(path)
    ranked = sorted(cavities, key=lambda cavity: cavity.score, reverse=True)[:top_n]
    text = FigureText(f"Volume (bars, left axis) and maximum radius (line, right axis) of the top {len(ranked)} detected cavities by score")
    fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
    if not ranked:
        _empty_panel(ax, "No cavities detected", text)
    else:
        labels = [f"C{cavity.id}" for cavity in ranked]
        volumes = [cavity.volume for cavity in ranked]
        radii = [cavity.radius for cavity in ranked]
        ax.bar(labels, volumes, color=CREVICE_TEAL, alpha=0.86, label="Volume")
        ax.set_ylabel("Volume (Å³)")
        ax.set_xlabel("Cavity")
        twin = ax.twinx()
        twin.plot(labels, radii, color=CREVICE_RED, marker="o", linewidth=1.8, label="Max radius")
        twin.set_ylabel("Max radius (Å)")
        twin.spines["top"].set_visible(False)
        _set_title(ax, text, title, "CREVICE cavity summary")
        ax.grid(True, axis="y", color="#d8d8d8", linewidth=0.65, alpha=0.75)
    _style_axes(ax)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


@annotate_option
def plot_network_summary(
    network: ResidueInteractionNetwork,
    path: str | Path,
    *,
    top_n: int = 20,
    dpi: int = 400,
    title: str | None = None,
    annotate: bool | None = None,
) -> Path:
    """Rank residues by their number of residue contacts, split by contact class.

    For every residue node the residue-residue edges of ``network`` (from
    :func:`crevice.networks.build_cavity_network` or
    :func:`crevice.networks.build_residue_network`) are counted; this count is
    the residue's degree in the residue-only graph, i.e. the number of other
    residues within the network's contact cutoff under its distance metric
    (``network.metadata["cutoff"]``/``["distance_metric"]``). The ``top_n``
    residues with most contacts (ties by residue ID) are drawn as horizontal
    bars stacked by edge interaction class. Region pseudo-nodes and
    region-residue edges are not counted: which residues line the region and how
    strongly is shown by :func:`plot_residue_contacts`, and the full degree
    including region edges stays in the network table
    (``PREFIX_network_nodes.csv``).

    Parameters
    ----------
    network : ResidueInteractionNetwork
        Residue (and optional region) contact network.
    path : str or Path
        Output image; parent directories are created.
    top_n : int, default 20
        Number of residues shown.
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title unless
        annotating ("CREVICE residue-interaction network").
    annotate : bool, optional
        Also draw the default title. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Outputs
    -------
    path : image
        Horizontal stacked bar chart, height scaled with the number of
        residues. Axes: "Residue contacts (degree: residues with <metric> ≤
        <cutoff> Å)", for example "atom-centre distance ≤ 4.5 Å" or "atom-surface
        gap ≤ 4.5 Å", and "Residue" (residue IDs as tick labels). With no
        residue-residue edges the panel states "No residue contacts in the network".

    Colours and representation
    --------------------------
    Bar segments are coloured by the edge interaction class written by
    :mod:`crevice.networks`: purple (``#6f4aa8``) salt-bridge candidate, orange
    (``#d9822b``) oppositely charged residue contact, green (``#2e7d32``)
    hydrophobic, teal (``#008b8b``) polar, grey (``#4d4d4d``) other distance
    contact. A legend below the axes names the classes present (a
    categorical colour key). White segment edges; light-grey vertical grid.
    """

    plt = _pyplot()
    output = _ensure_parent(path)
    residue_ids = [node.id for node in network.nodes if node.kind != "region"]
    region_ids = {node.id for node in network.nodes if node.kind == "region"}
    counts: dict[str, dict[str, int]] = {node_id: {} for node_id in residue_ids}
    for edge in network.edges:
        if edge.source in region_ids or edge.target in region_ids or edge.interaction.startswith("region_"):
            continue
        key = _interaction_key(edge.interaction)
        for node_id in (edge.source, edge.target):
            row = counts.setdefault(node_id, {})
            row[key] = row.get(key, 0) + 1
    totals = {node_id: sum(row.values()) for node_id, row in counts.items()}
    top = [item for item in sorted(totals.items(), key=lambda item: (-item[1], item[0])) if item[1] > 0][:top_n]
    criterion = _contact_criterion(network)
    xlabel = f"Residue contacts (degree: residues with {criterion})"
    text = FigureText(f"Number of residue-residue contacts (degree in the residue contact network, {criterion}) of the "
                      f"{len(top)} residues with most contacts, split by contact class; region pseudo-node edges are not counted")
    height = max(3.6, 0.32 * max(len(top), 1) + 1.2)
    fig, ax = plt.subplots(figsize=(7.2, height), constrained_layout=True)
    if not top:
        _empty_panel(ax, "No residue contacts in the network", text)
    else:
        from matplotlib.patches import Patch
        labels = [item[0] for item in reversed(top)]
        left = [0] * len(labels)
        present = []
        for key, _label in INTERACTION_LEGEND:
            if key == "region":
                continue
            values = [counts[node_id].get(key, 0) for node_id in labels]
            if not any(values):
                continue
            present.append(key)
            ax.barh(labels, values, left=left, color=_interaction_color(key), edgecolor="white", linewidth=0.5)
            left = [a + b for a, b in zip(left, values)]
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Residue")
        ax.xaxis.get_major_locator().set_params(integer=True)
        _set_title(ax, text, title, "CREVICE residue-interaction network")
        ax.grid(True, axis="x", color="#d8d8d8", linewidth=0.65, alpha=0.75)
        handles = [Patch(facecolor=_interaction_color(key), edgecolor="white", label=label)
                   for key, label in INTERACTION_LEGEND if key in present]
        fig.legend(handles=handles, title="Contact class", loc="outside lower center", ncol=3,
                   fontsize=7.5, title_fontsize=8, frameon=False)
    _style_axes(ax)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


@annotate_option
def plot_network_chord(
    network: ResidueInteractionNetwork,
    path: str | Path,
    *,
    top_n: int = 24,
    max_edges: int = 90,
    dpi: int = 400,
    title: str | None = None,
    annotate: bool | None = None,
) -> Path:
    """Draw a chord-style diagram of the most connected network nodes.

    The ``top_n`` highest-degree nodes are placed on a circle (a region node,
    if any, first at the top, then residues in sorted order). Edges between
    shown nodes are drawn as quadratic curves through the centre, strongest
    first, up to ``max_edges``; line width scales with edge weight. Node names
    are the circle's category labels and are always drawn, like the tick labels
    of a bar chart.

    Parameters
    ----------
    network : ResidueInteractionNetwork
        Residue (and optional region) contact network.
    path : str or Path
        Output image; parent directories are created.
    top_n : int, default 24
        Number of nodes placed on the circle.
    max_edges : int, default 90
        Maximum number of edges drawn.
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title unless
        annotating ("CREVICE residue interaction chord diagram").
    annotate : bool, optional
        Also draw the default title. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Outputs
    -------
    path : image
        7.2 x 7.2 inch figure without axes. With no nodes it states "No network
        nodes detected".

    Colours and representation
    --------------------------
    * Nodes: red circles (``#b22222``, larger) for the region node, blue
      (``#1f77b4``) for residues; grey node names outside the circle.
    * Edges (28% opacity, width 0.55-2.85 pt by relative weight), coloured by
      the interaction class written by :mod:`crevice.networks`: red
      (``#b22222``) channel/cavity (region) contact, purple (``#6f4aa8``)
      salt-bridge candidate, orange (``#d9822b``) oppositely charged residue
      contact, green (``#2e7d32``) hydrophobic, teal (``#008b8b``) polar, grey
      (``#4d4d4d``) other distance contact.
    * A figure legend (up to three columns) below the whole circle, node names
      included, names the edge classes that are drawn (several series
      distinguished only by colour). The layout reserves its own band for it,
      so it cannot overlap a residue name whatever the number of nodes or the
      length of their names.
    """

    plt = _pyplot()
    from matplotlib.lines import Line2D
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import Circle, PathPatch

    output = _ensure_parent(path)
    degrees = network.degree()
    nodes_by_id = {node.id: node for node in network.nodes}
    selected = [node_id for node_id, _degree in sorted(degrees.items(), key=lambda item: (-item[1], item[0]))[:top_n]]
    if any(node.kind == "region" for node in network.nodes):
        for node in network.nodes:
            if node.kind == "region" and node.id not in selected:
                selected = [node.id] + selected[:-1]
                break
    selected_set = set(selected)
    edges = [edge for edge in network.edges if edge.source in selected_set and edge.target in selected_set]
    edges.sort(key=lambda edge: (-edge.weight, edge.distance, edge.source, edge.target))
    edges = edges[:max_edges]
    text = FigureText(f"Chord diagram of {len(selected)} highest-degree network nodes and {len(edges)} strongest edges among them; "
                      "edge colour is the interaction type and width the relative weight")

    fig, ax = plt.subplots(figsize=(7.2, 7.2), constrained_layout=True)
    ax.set_aspect("equal")
    ax.axis("off")
    if not selected:
        _empty_panel(ax, "No network nodes detected", text)
        _save_figure(fig, output, dpi=dpi, text=text)
        return output

    radius = 1.0
    angles = _node_angles(selected, nodes_by_id)
    positions = {node_id: (math.cos(angle) * radius, math.sin(angle) * radius) for node_id, angle in angles.items()}
    max_weight = max((edge.weight for edge in edges), default=1.0)
    for edge in edges:
        x1, y1 = positions[edge.source]
        x2, y2 = positions[edge.target]
        path_data = [
            (MplPath.MOVETO, (x1, y1)),
            (MplPath.CURVE3, (0.0, 0.0)),
            (MplPath.CURVE3, (x2, y2)),
        ]
        codes, verts = zip(*path_data)
        patch = PathPatch(
            MplPath(verts, codes),
            facecolor="none",
            edgecolor=_interaction_color(edge.interaction),
            linewidth=0.55 + 2.3 * (edge.weight / max_weight),
            alpha=0.28,
            zorder=1,
        )
        ax.add_patch(patch)

    for node_id in selected:
        node = nodes_by_id[node_id]
        x, y = positions[node_id]
        color = CREVICE_RED if node.kind == "region" else CREVICE_BLUE
        size = 0.043 if node.kind == "region" else 0.032
        ax.add_patch(Circle((x, y), size, facecolor=color, edgecolor="white", linewidth=0.8, zorder=3))
        angle = angles[node_id]
        ha = "left" if math.cos(angle) >= -1e-9 else "right"  # bottom label runs outward, not over its node
        rotation = math.degrees(angle)
        if 90 < rotation < 270:
            rotation += 180
        label_radius = 1.13
        # Node names are this chart's category (tick) labels.
        ax.text(
            math.cos(angle) * label_radius,
            math.sin(angle) * label_radius,
            _short_residue_label(node.label),
            ha=ha,
            va="center",
            rotation=rotation,
            rotation_mode="anchor",
            fontsize=7.2 if node.kind != "region" else 8.2,
            color=CREVICE_GRAY,
        )

    _set_title(ax, text, title, "CREVICE residue interaction chord diagram", fontsize=12, color="#222222", pad=18)
    drawn = {_interaction_key(edge.interaction) for edge in edges}
    legend_items = [Line2D([0], [0], color=_interaction_color(key), lw=2, label=label)
                    for key, label in INTERACTION_LEGEND if key in drawn]
    ax.set_xlim(-1.38, 1.38)
    ax.set_ylim(-1.32, 1.35)
    if legend_items:
        # A figure legend "outside" the axes: constrained layout reserves a band
        # below the axes' full extent, residue names included, so the legend can
        # never overlap a label however many nodes or how long their names are.
        fig.legend(handles=legend_items, frameon=False, loc="outside lower center",
                   ncol=min(3, len(legend_items)), fontsize=8)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output


@annotate_option
def plot_trajectory_profiles(
    trajectory, path, *, dpi=400,title=None,view="distribution",quantiles=(.1,.9),
    coordinate_grid=None,samples=101,confidence=None,block_length=None,
    bootstrap_replicates=2000,seed=20260911,distribution=None,quantity="radius",
    annotate=None,
):
    """Plot aligned per-frame channel profiles as a distribution or a timeseries.

    ``view="distribution"`` resamples every resolved frame profile onto a common
    axial coordinate (see :func:`crevice.ensemble.profile_distribution`) and
    draws the per-coordinate frame range (fluctuations), the mean and, when
    estimable, an approximate simultaneous confidence band for the mean from a
    block bootstrap over frames. A lower panel shows the fraction of frames
    whose profile covers each coordinate. ``view="timeseries"`` draws each
    frame's minimum, mean and maximum radius against time (ps), or frame index
    when times are missing; unresolved frames are gaps. Frame fluctuation bands
    and mean confidence bands are different quantities and are drawn separately.

    Parameters
    ----------
    trajectory : TrajectoryAnalysis or sequence of FrameAnalysis
        Per-frame results, for example from :func:`crevice.trajectory.analyze_trajectory`.
    path : str or Path
        Output image; parent directories are created.
    dpi : int, default 400
        Raster resolution in dots per inch.
    title : str, optional
        Caller-supplied title; always drawn when given. Otherwise no title unless
        annotating, when a status title with the resolved/total frame count is used.
    view : {"distribution", "timeseries"}, default "distribution"
        Figure type (CLI ``trajectory --view``).
    quantiles : tuple of float, default (0.1, 0.9)
        Lower and upper per-coordinate frame quantiles of the fluctuation band.
    coordinate_grid : sequence of float, optional
        Explicit axial coordinates (Å) for resampling; default spans the profiles.
    samples : int, default 101
        Number of resampled coordinates when ``coordinate_grid`` is not given.
    confidence : float, optional
        Nominal level of the mean band (for example 0.95); ``None`` draws no band.
    block_length : int, optional
        Bootstrap block length in frames; default estimated from autocorrelation.
    bootstrap_replicates : int, default 2000
        Bootstrap replicates for the mean band.
    seed : int, default 20260911
        Bootstrap random seed.
    distribution : dict, optional
        A precomputed :func:`~crevice.ensemble.profile_distribution` result, reused as is.
    quantity : {"radius", "diameter"}, default "radius"
        Plot radii or diameters (2 x radius); CLI ``--profile-quantity``.
    annotate : bool, optional
        Also draw the status title and, when the mean band is unavailable, the
        in-plot note giving the reason. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Raises
    ------
    ValueError
        For an unknown ``view`` or ``quantity``.

    Outputs
    -------
    path : image
        Distribution view: 7.8 x 6 inch, "Pore radius (Å)" (or "Clearance
        radius (Å)" for an unvalidated fixed-axis scan; "diameter" when
        requested) over "Frame coverage", x axis "Aligned axial coordinate (Å)".
        Timeseries view: 7.8 x 4.8 inch, "Time (ps)" or "Frame". The PNG
        ``Description`` metadata records the resolved/total frame count and
        confidence-band status.

    Colours and representation
    --------------------------
    * Distribution: light-teal band (17% opacity) = frame quantile range;
      darker teal band (46%) = approximate simultaneous mean confidence band;
      teal line = mean; grey line in the lower panel = frame coverage. A legend
      names the bands and mean (several series).
    * Timeseries: red = minimum, blue = mean, grey = maximum radius per frame,
      with a legend.
    * Light-grey grid; the radius axis starts at zero.
    """
    import numpy as np
    if view not in {"distribution","timeseries"}:raise ValueError("view must be distribution or timeseries")
    if quantity not in {"radius","diameter"}:raise ValueError("quantity must be radius or diameter")
    plt=_pyplot();output=_ensure_parent(path)
    frames=tuple(trajectory.frames if isinstance(trajectory,TrajectoryAnalysis) else trajectory)
    profiles=[f.profile for f in frames if f.profile is not None]
    scale=2 if quantity=="diameter" else 1
    text=FigureText(f"Aligned channel {quantity} {'distribution along the axis' if view=='distribution' else 'timeseries'}; "
                    f"{len(profiles)} of {len(frames)} frames have a resolved profile")
    if view=="distribution" and profiles:
        from .ensemble import profile_distribution
        distribution=distribution or profile_distribution(trajectory,quantiles=quantiles,coordinate_grid=coordinate_grid,
                      samples=samples,confidence=confidence,block_length=block_length,
                      bootstrap_replicates=bootstrap_replicates,seed=seed)
        rows=distribution['rows'];x=[r['coordinate'] for r in rows]
        values=lambda key:np.asarray([np.nan if r.get(key) is None else r[key]*scale for r in rows])
        fig,(ax,coverage)=plt.subplots(2,1,figsize=(7.8,6),sharex=True,gridspec_kw={'height_ratios':[4,1]},layout='constrained')
        ax.fill_between(x,values('lower_radius'),values('upper_radius'),color=CREVICE_TEAL,alpha=.17,
                        label=f'{100*quantiles[0]:g}–{100*quantiles[1]:g}% frame range (fluctuations)')
        interval=distribution.get('confidence_interval') or {}
        if interval.get('status')=='estimated':
            ax.fill_between(x,values('mean_band_lower'),values('mean_band_upper'),color=CREVICE_TEAL,alpha=.46,
                            label=f"{100*interval['confidence']:g}% approximate simultaneous mean band")
        elif interval:
            text.text(ax,.02,.02,'Mean confidence unavailable: '+interval.get('status','unknown').replace('_',' '),
                      transform=ax.transAxes,fontsize=8,wrap=True)
        ax.plot(x,values('mean_radius'),color=CREVICE_TEAL,lw=2,label='Mean')
        coverage.plot(x,[r['n_frames']/len(frames) for r in rows],color=CREVICE_GRAY,lw=1.5)
        coverage.set_ylim(-.05,1.05);coverage.set_ylabel('Frame\ncoverage');coverage.set_xlabel('Aligned axial coordinate (Å)')
        coverage.set_yticks([0,.5,1]);_style_axes(coverage)
        unvalidated=any(p.metadata.get('status')=='unvalidated_axis_scan' for p in profiles)
        _set_title(ax,text,title,('Unvalidated fixed-axis clearance' if unvalidated else 'Aligned pore profile')+f' · {len(profiles)}/{len(frames)} frames')
        ax.set_ylabel(('Clearance ' if unvalidated else 'Pore ')+quantity+' (Å)');ax.legend(frameon=False,fontsize=8)
    else:
        fig,ax=plt.subplots(figsize=(7.8,4.8),layout='constrained')
        if not profiles:
            _empty_panel(ax,'No reference through-channel profile resolved\nNo pore-width confidence band is available',text)
            _set_title(ax,text,title,f'Profile status · {len(frames)} frames retained')
        else:
            x=[f.time if f.time is not None else f.frame_index for f in frames]
            for attr,color,label in [('min_radius',CREVICE_RED,'Minimum'),('mean_radius',CREVICE_BLUE,'Mean'),('max_radius',CREVICE_GRAY,'Maximum')]:
                y=[getattr(f.profile,attr)*scale if f.profile is not None else np.nan for f in frames]
                ax.plot(x,y,color=color,lw=1.6,label=label)
            ax.set_xlabel('Time (ps)' if all(f.time is not None for f in frames) else 'Frame')
            unvalidated=any(p.metadata.get('status')=='unvalidated_axis_scan' for p in profiles)
            ax.set_ylabel(('Clearance ' if unvalidated else 'Pore ')+quantity+' (Å)');_set_title(ax,text,title,'Unvalidated fixed-axis timeseries' if unvalidated else 'Profile timeseries · unresolved frames remain gaps');ax.legend(frameon=False)
    ax.set_ylim(bottom=0);ax.grid(True,color='#d8d8d8',lw=.65,alpha=.75);_style_axes(ax)
    _save_figure(fig,output,dpi=dpi,text=text)
    return output


def write_publication_pymol_script(
    *,
    structure_path: str | Path,
    pore_cast_path: str | Path,
    script_path: str | Path,
    profile: PoreProfile | None = None,
    output_png: str | Path | None = None,
    width: int = 2400,
    height: int = 1800,
    dpi: int = 400,
    structure_object: str = "crevice_cartoon",
    pore_object: str = "crevice_pore_cast",
    pore_fill_object: str = "crevice_pore_fill",
    structure_local_name: str | None = None,
    pore_cast_local_name: str | None = None,
    output_png_local_name: str | None = None,
    bottleneck_resi: int | None = None,
) -> Path:
    """Write a PyMOL script that renders a pore cast inside its protein.

    The script reloads the structure and the pore-cast dummy-atom PDB (whose
    B-factor holds each sample's clearance radius), sets each dummy atom's van
    der Waals radius to ``max(0.25, B)`` and shows them as translucent spheres
    coloured by radius, together with a molecular surface of the same spheres.
    Files are read from the current directory when present there, otherwise
    from their absolute paths, so the bundle can be moved. The scene contains
    no text, labels or legends.

    Parameters
    ----------
    structure_path : str or Path
        Protein structure loaded as ``structure_object``.
    pore_cast_path : str or Path
        Pore-cast dummy-atom PDB (B-factor = clearance radius, Å).
    script_path : str or Path
        Output ``.pml`` path.
    profile : PoreProfile, optional
        Used only to locate the bottleneck when ``bottleneck_resi`` is not given.
    output_png : str or Path, optional
        When given, the script ray-traces a ``width`` x ``height`` PNG at ``dpi``.
    width, height : int, default 2400, 1800
        Rendered image size in pixels.
    dpi : int, default 400
        DPI recorded in the rendered PNG.
    structure_object, pore_object, pore_fill_object : str
        PyMOL object names (defaults ``crevice_cartoon``, ``crevice_pore_cast``,
        ``crevice_pore_fill``).
    structure_local_name, pore_cast_local_name, output_png_local_name : str, optional
        File names looked up in the working directory before the absolute paths.
    bottleneck_resi : int, optional
        Dummy-atom residue number of the bottleneck sample, highlighted separately.

    Returns
    -------
    Path
        The written script.

    Outputs
    -------
    script_path : PyMOL script
        Run with ``pymol script.pml`` from the bundle directory.
    output_png : PNG, optional
        Written by PyMOL when the script runs, not by this function.

    Colours and representation
    --------------------------
    * Protein: cyan cartoon (``[0.00, 0.82, 0.92]``), 12% transparent.
    * Pore samples: spheres, 55% transparent, coloured by radius on PyMOL's
      ``blue_white_red`` spectrum (narrowest blue, widest red, white between).
    * Pore fill: pink molecular surface (``[1.00, 0.05, 0.55]``) of the
      samples, 18% transparent.
    * Bottleneck sample: opaque firebrick sphere.
    * White background, no depth cue or shadows.
    """

    output = _ensure_parent(script_path)
    structure_abs = Path(structure_path).resolve()
    pore_cast_abs = Path(pore_cast_path).resolve()
    structure_local = structure_local_name or structure_abs.name
    pore_cast_local = pore_cast_local_name or pore_cast_abs.name
    if bottleneck_resi is None and profile is not None:
        bottleneck_resi = profile.bottleneck.index + 1
    png_block = ""
    if output_png is not None:
        png_abs = Path(output_png).resolve()
        png_local = output_png_local_name or png_abs.name
        png_block = f"""python
_required = [{structure_local!r}, {pore_cast_local!r}]
_png_target = _crevice_output({png_local!r}, {str(png_abs)!r}, _required)
cmd.ray({width}, {height})
cmd.png(_png_target, dpi={dpi})
python end
"""
    bottleneck_lines = ""
    if bottleneck_resi is not None:
        bottleneck_lines = f"""select crevice_bottleneck, {pore_object} and resi {bottleneck_resi}
color firebrick, crevice_bottleneck
set sphere_transparency, 0.0, crevice_bottleneck
"""
    script = f"""reinitialize
set retain_order, 1
bg_color white
set ray_opaque_background, on
set antialias, 4
set ray_trace_mode, 1
set ray_shadows, off
set ambient, 0.36
set direct, 0.72
set specular, 0.18
set shininess, 24
set depth_cue, 0
set cartoon_fancy_helices, on
set cartoon_sampling, 14
set_color crevice_cyan, [0.00, 0.82, 0.92]
set_color crevice_pink, [1.00, 0.05, 0.55]
set_color crevice_deep_blue, [0.08, 0.18, 0.55]
python
from pathlib import Path
from pymol import cmd

def _crevice_input(local_name, fallback):
    local = Path.cwd() / local_name
    return str(local if local.exists() else Path(fallback))

def _crevice_output(local_name, fallback, required_local_inputs):
    if all((Path.cwd() / name).exists() for name in required_local_inputs):
        return str(Path.cwd() / local_name)
    return str(Path(fallback))

cmd.load(_crevice_input({structure_local!r}, {str(structure_abs)!r}), {structure_object!r})
cmd.load(_crevice_input({pore_cast_local!r}, {str(pore_cast_abs)!r}), {pore_object!r})
python end
hide everything, {structure_object}
show cartoon, {structure_object}
color crevice_cyan, {structure_object}
set cartoon_transparency, 0.12, {structure_object}
hide everything, {pore_object}
alter {pore_object}, vdw=max(0.25, b)
rebuild
create {pore_fill_object}, {pore_object}
hide everything, {pore_fill_object}
show surface, {pore_fill_object}
color crevice_pink, {pore_fill_object}
set surface_color, crevice_pink, {pore_fill_object}
set surface_quality, 2, {pore_fill_object}
set transparency, 0.18, {pore_fill_object}
show spheres, {pore_object}
set sphere_quality, 4, {pore_object}
set sphere_transparency, 0.55, {pore_object}
spectrum b, blue_white_red, {pore_object}
{bottleneck_lines}orient {structure_object} or {pore_object} or {pore_fill_object}
zoom {structure_object} or {pore_object} or {pore_fill_object}, 6
{png_block}"""
    output.write_text(script)
    return output


@radii_option
@annotate_option
def write_static_publication_bundle(
    frame: StructureFrame,
    *,
    structure_path: str | Path,
    output_dir: str | Path,
    prefix: str = "crevice",
    axis: str | tuple[float, float, float] = "auto",
    origin: tuple[float, float, float] | None = None,
    section_spacing: float = 0.5,
    enclosure_radius: float | None = None,
    include_hydrogen: bool = False,
    include_hetero: bool = True,
    samples: int = 101,
    lateral_exits: bool = False,
    exit_bulk_radius: float = 6.0,
    exit_spacing: float = 0.5,
    search_radius: float = 6.0,
    refinement_steps: int = 4,
    probe_radius: float = 0.0,
    dpi: int = 400,
    include_cavities: bool = True,
    include_tunnels: bool = True,
    include_network: bool = True,
    residue_label_every: int = 5,
    cast_mode: str = "auto",
    cast_outer_radius: float = 6.0,
    cast_enclosure_fraction: float = .9,
    cast_extension: float = 2.0,
    cast_spacing: float = 0.5,
    cast_min_radius: float = 0.0,
    cast_max_components: int = 1,
    cast_selection: str = "dominant",
    cast_core_fraction: float = .08,
    surface_smoothing: float | None = None,
    cast_max_grid_points: int = 16_000_000,
    cast_focus_points: Sequence[tuple[float, float, float]] | None = None,
    cast_focus_radius: float = 8.0,
    cast_min_component_volume: float = 0.0,
    cast_max_component_volume: float | None = None,
    residue_groups: dict[str, str] | None = None,
    smooth: float | None = None,
    annotate: bool | None = None,
    radii: RadiusSet | str | None = None,
) -> dict[str, str]:
    """Write the complete static-structure bundle of ``crevice publish``.

    The frame is reduced to the requested atom selection, a single
    through-channel profile is resolved (:func:`crevice.channels.pore_profile`),
    residues within the profile's reach are annotated and a void cast is filled
    (:func:`crevice.voids.void_cast`). Tables, figures and native viewer scenes
    (PyMOL, VMD, ChimeraX) are then written, followed by the requested cavity,
    tunnel and network analyses. When no unique through-channel resolves and
    ``cast_mode`` is ``"auto"`` with default crops, the bundle falls back to a
    3D rolling-probe cast (:func:`write_rolling_probe_bundle`), still runs the
    requested extras on that measured region, writes
    ``PREFIX_profile_status.json`` and never invents a radius profile.
    ``cast_mode="rolling"`` delegates to :func:`write_rolling_probe_bundle`
    directly. Display smoothing changes display surfaces only; measured grids,
    volumes and profiles are identical for every ``smooth`` setting.

    Parameters
    ----------
    frame : StructureFrame
        Parsed structure.
    structure_path : str or Path
        Source file, copied into the bundle as ``PREFIX_structure.<ext>``.
    output_dir : str or Path
        Bundle directory (created).
    prefix : str, default "crevice"
        File-name prefix (CLI default: the structure file stem).
    axis : {"auto", "x", "y", "z"} or tuple of float, default "auto"
        Channel direction; ``auto`` uses the principal axis.
    origin : tuple of float, optional
        Seed point inside the intended channel (Å).
    section_spacing : float, default 0.5
        Cross-section grid spacing (Å).
    enclosure_radius : float, optional
        Probe radius (Å) deciding which cross-section samples are enclosed
        channel. ``None`` (default; CLI ``--enclosure-radius auto``) chooses
        it automatically (:func:`crevice.sections.connected_profile`); the
        choice is recorded in the profile metadata and in the manifest's
        ``profile_status.enclosure_probe``. Also the rolling-probe radius of
        the ``rolling`` cast mode and of the fallback cast written when no
        channel resolves; there ``None`` means the legacy 0.8 Å (recorded as
        ``fallback_cast_probe_A``).
    include_hydrogen : bool, default False
        Treat hydrogens as obstacles.
    include_hetero : bool, default True
        Treat HETATM records (ligands, ions, waters) as obstacles.
    samples : int, default 101
        Number of profile samples along the axis.
    search_radius : float, default 6.0
        Transverse centre-refinement half-width (Å).
    refinement_steps : int, default 4
        Centre-refinement iterations.
    probe_radius : float, default 0.0
        Radius (Å) subtracted from reported clearances.
    dpi : int, default 400
        Resolution of every Matplotlib figure and of rendered viewer images.
    include_cavities, include_tunnels, include_network : bool, default True
        Run the cavity search, tunnel search and residue network (CLI
        ``--skip-cavities``/``--skip-tunnels``/``--skip-network``).
    residue_label_every : int, default 5
        Landmark spacing (profile samples) in ``PREFIX_profile_radius_annotated.png``.
    cast_mode : {"auto", "channel", "cavity", "all", "rolling"}, default "auto"
        Which void components form the cast.
    cast_outer_radius : float, default 6.0
        Outer envelope probe of rolling casts (Å).
    cast_enclosure_fraction : float, default 0.9
        Burial threshold for rolling-cast cores.
    cast_extension : float, default 2.0
        Display-only cast continuation beyond each channel mouth (Å); excluded
        from measured volume.
    cast_spacing : float, default 0.5
        Cast grid spacing (Å); never silently coarsened.
    cast_min_radius : float, default 0.0
        Minimum atom-surface clearance of cast samples (Å).
    cast_max_components : int, default 1
        Maximum retained cast regions.
    cast_selection : {"dominant", "all"}, default "dominant"
        Region selection for rolling casts.
    cast_core_fraction : float, default 0.08
        Minimum core volume relative to the strongest core (dominant selection).
    surface_smoothing : float, optional
        Legacy display interpolation width (Å); cannot be combined with ``smooth``.
    cast_max_grid_points : int, default 16000000
        Grid allocation guard for the cast.
    cast_focus_points : sequence of tuple of float, optional
        Keep only cast samples within ``cast_focus_radius`` of these points (Å).
    cast_focus_radius : float, default 8.0
        Focus radius (Å).
    cast_min_component_volume, cast_max_component_volume : float, optional
        Inclusive volume limits (Å³) on final components after cropping.
    residue_groups : dict of str to str, optional
        Residue ID to helix/strand name; adds ``PREFIX_connectivity.{json,png}``.
    smooth : float, optional
        Display-surface smoothing length (Å, default 0.4; 0 = raw voxel
        boundary). See :func:`crevice.volume_export.resolve_display_smoothing`.
    lateral_exits : bool, default False
        Allow a capped channel end (blocked straight exit) to leave through the
        widest free lateral path to rolling-probe bulk; passed to
        :func:`crevice.analysis.pore_profile`. Off by default; without it the
        output is identical to earlier versions.
    exit_bulk_radius : float, default 6.0
        Rolling-probe radius (Å) defining bulk solvent for lateral exits.
    exit_spacing : float, default 0.5
        Lattice spacing (Å) of the lateral-exit path search.
    annotate : bool, optional
        Draw titles, value call-outs, mouth labels and residue landmark names in
        the figures (the viewer scenes of this bundle contain no text).
        ``None`` (default) inherits the surrounding setting, which is off unless
        enabled with :func:`crevice.presentation.figure_annotations` or ``--annotate``.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call (profiles, casts, cavities, tunnels, networks,
        hydration SASA) comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict of str to str
        Manifest mapping each output key to its path, including ``manifest_json``.

    Raises
    ------
    ChannelResolutionError
        When no channel resolves and a crop or non-auto mode forbids the fallback.
    ValueError
        For invalid cast or smoothing settings.

    Outputs
    -------
    PREFIX_manifest.json : JSON
        File index, DPI, void-cast summary and ``profile_status``.
    PREFIX_structure.<ext>, PREFIX_viewer.pdb, PREFIX_viewer.identities.json : structure
        Copy of the input; renumbered viewer PDB and its identity map.
    PREFIX_profile.csv : table
        One row per profile sample (:func:`crevice.io.write_profile_csv`:
        ``x_A``/``y_A``/``z_A``, ``radius_A``, ``raw_clearance_A``,
        ``axis_position_A``, nearest atom and residue).
    PREFIX_profile.json : JSON
        The full profile record (samples plus method, settings, bottleneck,
        channel mouths and exits), kept for machine-readable reuse.
    PREFIX_residue_contacts.csv : table
        Residues contacting the profile (or cast) with gaps (Å), roles and
        influence scores (:func:`crevice.io.write_residue_contacts_csv`).
    PREFIX_void_cast.csv : table
        One row per cast region (:func:`crevice.io.write_void_cast_csv`).
    PREFIX_void_cast.json, PREFIX_pore_cast.pdb : cast
        Full cast record (settings, diagnostics); dummy atoms with clearance in
        the B-factor.
    PREFIX_volume.dx, PREFIX_display.dx, PREFIX_pymol_local.dx, PREFIX_pymol_mesh.npz : grids
        Measured binary grid; display-only (smoothed, extended) field and mesh.
    PREFIX_volume.{pml,vmd,tcl,cxc}, PREFIX_pore_cast.{pml,vmd,tcl,cxc}, PREFIX_render.cxc : scenes
        Native viewer scenes of the cast in the protein (see
        :func:`crevice.volume_export.write_volume_viewer_bundle`);
        ``PREFIX_pore_spheres.pml`` keeps the dummy-sphere view
        (:func:`write_publication_pymol_script`).
    PREFIX_scene.json, PREFIX_volume_metadata.json, PREFIX_geometry_reference.npz, PREFIX_mouths.bild : scene data
        Camera, display definitions, reference points and mouth guides read by
        the viewer scripts.
    PREFIX_profile_radius.png : figure
        :func:`plot_profile_radius`.
    PREFIX_profile_radius_annotated.png : figure
        :func:`plot_profile_radius_with_residues` (landmark names only when annotating).
    PREFIX_residue_contacts.png : figure
        :func:`plot_residue_contacts` (with a contact-role legend).
    PREFIX_cavities.csv, PREFIX_cavities.json, PREFIX_cavities.pdb, PREFIX_cavity_summary.png : cavities
        :func:`crevice.analysis.detect_cavities` results as a table
        (:func:`crevice.io.write_cavities_csv`), JSON record and dummy atoms,
        and :func:`plot_cavity_summary`.
    PREFIX_tunnels.csv, PREFIX_tunnel_points.csv, PREFIX_tunnels.json, PREFIX_tunnels.pdb : tunnels
        :func:`crevice.tunnels.find_tunnels` results: one row per tunnel, one
        row per path node, the JSON record (with per-segment metadata) and dummy
        atoms. No tunnel figure is written.
    PREFIX_network_nodes.csv, PREFIX_network_edges.csv, PREFIX_network.json : network
        Contact network as node and edge tables (:func:`crevice.io.write_network_csv`)
        and the JSON record with metrics.
    PREFIX_network_summary.png, PREFIX_network_chord.png : figures
        :func:`plot_network_summary` and :func:`plot_network_chord`.
    PREFIX_connectivity.csv, PREFIX_connectivity_groups.csv, PREFIX_connectivity.json, PREFIX_connectivity.png : connectivity
        Only with ``residue_groups``: paths to the lining and group contacts
        (:func:`crevice.io.write_connectivity_csv`), the JSON record and
        :func:`crevice.connectivity.plot_lining_connectivity`.
    PREFIX_profile_status.json : JSON
        Only after the unresolved-profile fallback (whose other files are those
        of :func:`write_rolling_probe_bundle`).

    The CLI adds ``PREFIX_input_report.json`` (input provenance) and the
    standard hydration outputs (:func:`crevice.hydration_export.write_hydration_bundle`).
    No HTML, 3D preview image or tunnel summary figure is written.

    With ``lateral_exits``, the profile metadata gains ``exits``, ``path_type``
    and ``path_bottleneck_radius_A``, and the manifest's ``profile_status``
    records ``path_type``/``exit_types``. The measured cast (``total_volume``,
    ``PREFIX_volume.dx``, ``PREFIX_pore_cast.pdb``) covers only the axial lumen
    between the mouths. Each lateral exit leg also gets a probe-swept cast
    (:func:`crevice.channel_exits.lateral_exit_casts`), reported separately:
    ``metadata["lateral_exit_casts"]`` in ``PREFIX_void_cast.json`` and the
    manifest, one ``lateral_exit`` row per leg in ``PREFIX_void_cast.csv``,
    ``profile_status.lateral_exit_cast_volumes_A3`` and ``PREFIX_exit_casts.dx``.
    The scenes draw them joined to the channel cast.

    Every PNG embeds its description (``Title``, ``Description``) as metadata;
    read it with :func:`crevice.presentation.figure_description`.

    Colours and representation
    --------------------------
    * Figures: see each plotting function.
    * Viewer scenes: grey-blue protein cartoon (25% transparent) around an
      opaque cast surface: teal (``[0.08, 0.58, 0.63]``) for a channel cast,
      violet (``[0.659, 0.333, 0.969]``, ``#a855f7``) for any other cast,
      including the unresolved-profile fallback; white background. Multiple
      display regions use the cast colour, purple, blue and mauve, orange for a
      shared core and translucent grey for other voids.
    * Channel mouths: no guide by default; orange rings (``[0.90, 0.52, 0.12]``)
      with mouth guides (:func:`crevice.presentation.display_guides`, CLI
      ``--mouth-guides``). Positions are always in the JSON outputs.
    * Lateral exit legs (``lateral_exits``): drawn as cast surfaces joined to
      the channel cast at the capped mouth, in the channel colour. Thin blue
      centre-line tubes (RGB ``[0.16, 0.47, 0.84]``, radius 0.12 Å, PyMOL object
      ``crevice_exits``; VMD colour 23; ChimeraX cylinders in
      ``PREFIX_mouths.bild``) only with ``--exit-centre-lines``.
    """
    from .volume_export import resolve_display_smoothing
    from .sections import LEGACY_ENCLOSURE_RADIUS
    resolve_display_smoothing(smooth, surface_smoothing)
    # Rolling-probe casts need a number; the automatic choice applies to the
    # channel profile only.
    cast_probe = LEGACY_ENCLOSURE_RADIUS if enclosure_radius is None else enclosure_radius

    if cast_mode == 'rolling':
        if cast_focus_points is not None or cast_max_component_volume is not None or cast_min_radius != 0:
            raise ValueError('Rolling publication mode does not support focus, maximum-volume or minimum-swept-clearance crops; use the explicit rolling probe and component settings')
        return write_rolling_probe_bundle(frame, structure_path=structure_path,
                   output_dir=output_dir, prefix=prefix, spacing=cast_spacing,
                   probe_radius=cast_probe, outer_radius=cast_outer_radius,
                   enclosure_fraction=cast_enclosure_fraction,
                   min_component_volume=cast_min_component_volume,
                   max_components=cast_max_components, max_grid_points=cast_max_grid_points,
                   selection='all' if include_hetero else 'protein', dpi=dpi,
                   selection_mode=cast_selection, core_fraction=cast_core_fraction, surface_smoothing=surface_smoothing, smooth=smooth)
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, str] = {}

    if not math.isfinite(cast_extension) or cast_extension < 0:
        raise ValueError("cast_extension must be finite and non-negative")
    frame = replace(frame, atoms=frame.selected_atoms(include_hydrogen=include_hydrogen,
                                                      include_hetero=include_hetero))
    original_copy = _copy_structure_for_bundle(structure_path, root, prefix)
    from .presentation import write_viewer_structure
    structure_copy = write_viewer_structure(frame, root / f"{prefix}_viewer.pdb")
    from .sections import ChannelResolutionError
    try:
        profile = pore_profile(
            frame,
            axis=axis,
            origin=origin,
            section_spacing=section_spacing,
            enclosure_radius=enclosure_radius,
            include_hydrogen=include_hydrogen,
            include_hetero=include_hetero,
            samples=samples,
            search_radius=search_radius,
            refinement_steps=refinement_steps,
            probe_radius=probe_radius,
            lateral_exits=lateral_exits,
            exit_bulk_radius=exit_bulk_radius,
            exit_spacing=exit_spacing,
        )
    except ChannelResolutionError as exc:
        if (cast_mode != 'auto' or cast_focus_points is not None or cast_max_component_volume is not None
                or cast_min_radius != 0 or include_hydrogen):
            raise
        # General voids remain useful when no unique through-channel exists.
        # Preserve the explicit selected obstacle set and never fabricate a
        # radius profile for a multichannel, closed or one-ended pocket system.
        manifest = write_rolling_probe_bundle(frame, structure_path=structure_path,
            output_dir=root, prefix=prefix, spacing=cast_spacing, probe_radius=cast_probe,
            outer_radius=cast_outer_radius, enclosure_fraction=cast_enclosure_fraction,
            min_component_volume=cast_min_component_volume, max_components=cast_max_components,
            max_grid_points=cast_max_grid_points, selection='all', dpi=dpi,
            selection_mode=cast_selection, core_fraction=cast_core_fraction, surface_smoothing=surface_smoothing, smooth=smooth)
        saved = json.loads(Path(manifest['manifest_json']).read_text())
        saved['single_channel_unresolved_reason'] = str(exc)
        saved['auto_fallback'] = '3D rolling-probe cast; no single-channel profile'
        saved['publication_input_selection'] = {'include_hydrogen':include_hydrogen, 'include_hetero':include_hetero}
        # Requested extras do not depend on a profile; complete them explicitly.
        fallback = _complete_unresolved_profile_extras(frame, manifest, root, prefix, dpi=dpi,
            include_cavities=include_cavities, include_tunnels=include_tunnels, include_network=include_network,
            include_hydrogen=include_hydrogen, include_hetero=include_hetero, residue_groups=residue_groups)
        status = {'status': 'unresolved', 'reason': str(exc),
                  'interpretation': 'No unique through-channel resolved; no pore-radius profile or profile-derived features assigned.',
                  'fallback_analyses': fallback}
        selection = getattr(exc, 'probe_selection', None)
        if selection is not None:
            # Automatic enclosure probe: every candidate and why none was used.
            status['enclosure_probe_selection'] = selection
            saved['enclosure_probe'] = {'mode': 'auto', 'chosen_A': None, 'reason': selection['reason'],
                                        'fallback_cast_probe_A': cast_probe}
        manifest['profile_status_json'] = str(root / f'{prefix}_profile_status.json')
        _write_json(status, manifest['profile_status_json'])
        saved['profile_status'] = {k: status[k] for k in ('status', 'reason', 'interpretation')}
        if selection is not None:
            saved['profile_status']['enclosure_probe'] = saved.pop('enclosure_probe')
        saved['fallback_analyses'] = fallback
        saved.update(radii_fields())
        saved['files'].update({k: v for k, v in manifest.items() if k != 'manifest_json'})
        _write_json(saved, manifest['manifest_json'])
        import warnings
        warnings.warn('No unique through-channel profile resolved; wrote a 3D rolling-probe cast. '
                      'See single_channel_unresolved_reason in the manifest.', RuntimeWarning, stacklevel=3)  # caller of the annotate_option wrapper
        return manifest
    if profile.method == "axial-connected" and cast_extension:
        from .spatial import SpatialIndex
        from .geometry import add, scale
        spatial = SpatialIndex(frame.atoms)
        context = []
        for point, sign in ((profile.points[0], -1), (profile.points[-1], 1)):
            segment = []
            count = max(2, int(math.ceil(cast_extension / 0.25)))
            for k in range(count+1):
                dt = sign * cast_extension * k/count
                position = add(point.position, scale(profile.axis_direction, dt))
                raw = spatial.nearest_surface(position)[0]
                segment.append({"t": point.t+dt, "radius": max(0.0, raw-probe_radius)})
            context.append(segment)
        profile = replace(profile, metadata={**profile.metadata, "display_context": context,
                          "display_context_definition": "clearance along axial exit; outside channel statistics"})
    contacts = annotate_residues(frame, profile.points)
    cast = void_cast(
        frame,
        mode=cast_mode,
        axis=axis,
        spacing=cast_spacing,
        min_radius=cast_min_radius,
        max_components=cast_max_components,
        max_grid_points=cast_max_grid_points,
        focus_points=cast_focus_points,
        focus_radius=cast_focus_radius,
        min_component_volume=cast_min_component_volume,
        max_component_volume=cast_max_component_volume,
        profile=profile,
    )
    exit_casts = ()
    if (profile.metadata.get("exits") and cast.mode == "channel"
            and any(e.get("type") == "lateral" for e in profile.metadata["exits"].values())):
        # Lateral exit legs as probe-swept casts, reported separately; the
        # axial cast and its volume are unchanged.
        from .channel_exits import EXIT_CAST_DEFINITION, exit_cast_summary, lateral_exit_casts
        exit_casts = lateral_exit_casts(frame, profile, cast)
        if exit_casts:
            cast = replace(cast, metadata={**cast.metadata,
                                           "lateral_exit_casts": exit_cast_summary(exit_casts),
                                           "lateral_exit_cast_definition": EXIT_CAST_DEFINITION,
                                           "lateral_exit_casts_note": "separate volumes; not part of total_volume"})
    if cast_focus_points is not None or cast_mode in {"cavity", "all"}:
        from .models import VoidComponent
        cast_points = tuple(p for c in cast.components for p in c.points)
        region = VoidComponent(0, "selected_void", (0.0, 0.0, 0.0), cast.total_volume, cast_points)
        contacts = annotate_residues(frame, cast_points, identify_bottleneck=False)
    else:
        region = profile

    manifest["structure_file"] = str(structure_copy)
    manifest["original_structure_file"] = str(original_copy)
    manifest["viewer_identity_json"] = str(structure_copy.with_suffix(".identities.json"))
    manifest["profile_json"] = str(root / f"{prefix}_profile.json")
    manifest["profile_csv"] = str(root / f"{prefix}_profile.csv")
    manifest["pore_cast_pdb"] = str(root / f"{prefix}_pore_cast.pdb")
    manifest["pore_cast_pml"] = str(root / f"{prefix}_pore_cast.pml")
    manifest["void_cast_json"] = str(root / f"{prefix}_void_cast.json")
    manifest["pymol_render_png"] = str(root / f"{prefix}_pore_cast.png")
    manifest["profile_radius_png"] = str(root / f"{prefix}_profile_radius.png")
    manifest["profile_radius_annotated_png"] = str(root / f"{prefix}_profile_radius_annotated.png")
    manifest["void_cast_csv"] = str(root / f"{prefix}_void_cast.csv")
    manifest["residue_contacts_csv"] = str(root / f"{prefix}_residue_contacts.csv")
    manifest["residue_contacts_png"] = str(root / f"{prefix}_residue_contacts.png")

    write_profile_json(profile, manifest["profile_json"])
    write_profile_csv(profile, manifest["profile_csv"])
    write_void_cast_json(cast, manifest["void_cast_json"])
    write_void_cast_csv(cast, manifest["void_cast_csv"], extra_components=exit_casts)
    write_void_cast_pdb(cast, manifest["pore_cast_pdb"])
    write_residue_contacts_csv(contacts, manifest["residue_contacts_csv"])
    plot_profile_radius(profile, manifest["profile_radius_png"], dpi=dpi)
    plot_profile_radius_with_residues(
        profile,
        manifest["profile_radius_annotated_png"],
        dpi=dpi,
        label_every=residue_label_every,
    )
    plot_residue_contacts(contacts, manifest["residue_contacts_png"], dpi=dpi)
    write_publication_pymol_script(
        structure_path=structure_copy,
        pore_cast_path=manifest["pore_cast_pdb"],
        script_path=manifest["pore_cast_pml"],
        profile=profile,
        output_png=manifest["pymol_render_png"],
        dpi=dpi,
        structure_local_name=Path(structure_copy).name,
        pore_cast_local_name=Path(manifest["pore_cast_pdb"]).name,
        output_png_local_name=Path(manifest["pymol_render_png"]).name,
        bottleneck_resi=_pdb_resi_for_min_radius(cast.points),
    )

    from .volume_export import write_volume_viewer_bundle
    if cast.components:
        manifest.update(write_volume_viewer_bundle(cast, structure_path=structure_copy,
                        output_dir=root, prefix=prefix, frame=frame, profile=profile,
                        cast_extension=cast_extension, surface_smoothing=surface_smoothing, smooth=smooth,
                        exit_casts=exit_casts))
        if profile.method == "axial-connected":
            # The primary pore-cast scene must show the measured full grid.
            # Retain the dummy-atom representation separately for compatibility.
            legacy = root / f"{prefix}_pore_spheres.pml"
            shutil.copyfile(manifest["pore_cast_pml"], legacy)
            manifest["pore_spheres_pml"] = str(legacy)
            shutil.copyfile(manifest["volume_pml"], manifest["pore_cast_pml"])
            manifest["pymol_render_png"] = str(root / f"{prefix}_pymol.png")
        for extension in ['vmd','tcl','cxc']:
            target=root/f'{prefix}_pore_cast.{extension}'
            shutil.copyfile(manifest[f'volume_{extension}'],target)
            manifest[f'pore_cast_{extension}']=str(target)


    if include_cavities:
        _write_cavity_outputs(frame, root, prefix, manifest, dpi=dpi)

    if include_tunnels:
        _write_tunnel_outputs(frame, root, prefix, manifest, dpi=dpi)

    if include_network:
        _write_network_outputs(frame, region, contacts, root, prefix, manifest,
                               dpi=dpi, residue_groups=residue_groups)

    manifest_path = root / f"{prefix}_manifest.json"
    _write_json({"files": manifest, "dpi": dpi, "prefix": prefix, "void_cast": cast.to_dict(), **radii_fields(),
                 "profile_status": {"status": "resolved", "method": profile.method,
                                    "profile_json": manifest["profile_json"],
                                    **({"path_type": profile.metadata["path_type"],
                                        "exit_types": {k: e["type"] for k, e in profile.metadata["exits"].items()},
                                        "exit_leg_counts": {k: e.get("leg_count", 0) for k, e in profile.metadata["exits"].items()}}
                                       if "path_type" in profile.metadata else {}),
                                    **({"lateral_exit_cast_volumes_A3": {f"{c.metadata['end']}_{c.metadata['leg']}": c.volume
                                                                         for c in exit_casts},
                                        "lateral_exit_cast_note": "separate from void_cast.total_volume (axial channel only)"}
                                       if exit_casts else {}),
                                    **({"enclosure_probe": {
                                        "mode": "auto",
                                        "chosen_A": profile.metadata["enclosure_probe_selection"]["chosen_A"],
                                        "reason": profile.metadata["enclosure_probe_selection"]["reason"],
                                        "candidates_tried": [r["probe_A"] for r in profile.metadata["enclosure_probe_selection"]["candidates"]
                                                             if r["status"] != "not_tried"]}}
                                       if "enclosure_probe_selection" in profile.metadata else {})}}, manifest_path)
    manifest["manifest_json"] = str(manifest_path)
    return manifest


def _write_cavity_outputs(frame, root, prefix, manifest, *, dpi, **options):
    cavities = detect_cavities(frame, **options)
    manifest["cavities_json"] = str(root / f"{prefix}_cavities.json")
    manifest["cavities_csv"] = str(root / f"{prefix}_cavities.csv")
    manifest["cavities_pdb"] = str(root / f"{prefix}_cavities.pdb")
    manifest["cavity_summary_png"] = str(root / f"{prefix}_cavity_summary.png")
    write_cavities_json(cavities, manifest["cavities_json"])
    write_cavities_csv(cavities, manifest["cavities_csv"])
    write_cavities_pdb(cavities, manifest["cavities_pdb"])
    plot_cavity_summary(cavities, manifest["cavity_summary_png"], dpi=dpi)
    return cavities


def _write_tunnel_outputs(frame, root, prefix, manifest, *, dpi, **options):
    tunnels = find_tunnels(frame, **options)
    manifest["tunnels_json"] = str(root / f"{prefix}_tunnels.json")
    manifest["tunnels_csv"] = str(root / f"{prefix}_tunnels.csv")
    manifest["tunnel_points_csv"] = str(root / f"{prefix}_tunnel_points.csv")
    manifest["tunnels_pdb"] = str(root / f"{prefix}_tunnels.pdb")
    write_tunnels_json(tunnels, manifest["tunnels_json"])
    write_tunnels_csv(tunnels, manifest["tunnels_csv"])
    write_tunnel_points_csv(tunnels, manifest["tunnel_points_csv"])
    write_tunnels_pdb(tunnels, manifest["tunnels_pdb"])
    return tunnels


def _write_network_outputs(frame, region, contacts, root, prefix, manifest, *, dpi, residue_groups=None):
    network = build_cavity_network(frame, region, contacts=contacts, residue_groups=residue_groups)
    if residue_groups is not None:
        from .connectivity import lining_connectivity, plot_lining_connectivity
        manifest["connectivity_json"] = str(root / f"{prefix}_connectivity.json")
        manifest["connectivity_csv"] = str(root / f"{prefix}_connectivity.csv")
        manifest["connectivity_groups_csv"] = str(root / f"{prefix}_connectivity_groups.csv")
        manifest["connectivity_png"] = str(root / f"{prefix}_connectivity.png")
        report = lining_connectivity(network)
        _write_json(report, manifest["connectivity_json"])
        write_connectivity_csv(report, manifest["connectivity_csv"], manifest["connectivity_groups_csv"])
        plot_lining_connectivity(network, manifest["connectivity_png"], dpi=dpi)
    manifest["network_json"] = str(root / f"{prefix}_network.json")
    manifest["network_nodes_csv"] = str(root / f"{prefix}_network_nodes.csv")
    manifest["network_edges_csv"] = str(root / f"{prefix}_network_edges.csv")
    manifest["network_summary_png"] = str(root / f"{prefix}_network_summary.png")
    manifest["network_chord_png"] = str(root / f"{prefix}_network_chord.png")
    _write_json({"network": network.to_dict(), "metrics": network_metrics(network)}, manifest["network_json"])
    write_network_csv(network, manifest["network_nodes_csv"], manifest["network_edges_csv"])
    plot_network_summary(network, manifest["network_summary_png"], dpi=dpi)
    plot_network_chord(network, manifest["network_chord_png"], dpi=dpi)
    return network


def _complete_unresolved_profile_extras(frame, manifest, root, prefix, *, dpi, include_cavities,
                                        include_tunnels, include_network, include_hydrogen,
                                        include_hetero, residue_groups=None):
    """Run requested cavity/tunnel/network extras after an unresolved profile.

    No radius profile or profile-derived feature is assigned. The measured
    fallback cast (complete binary DX samples) is the only region used: it seeds
    the tunnel search at its maximum atom-clearance sample and defines the
    network region node. Obstacles are the publication atom selection, as for
    the fallback cast itself.
    """
    import numpy as np
    obstacles = {"include_hydrogen": include_hydrogen, "include_hetero": include_hetero,
                 "definition": "Publication atom selection; the same obstacle atoms as the fallback rolling-probe cast"}
    records: dict[str, object] = {}
    region_xyz = None
    source_dx = manifest.get("volume_dx")
    if source_dx:
        from .residue_evidence import read_binary_dx
        grid, origin, deltas = read_binary_dx(source_dx)
        region_xyz = np.argwhere(grid) @ deltas + origin
        if not len(region_xyz):
            region_xyz = None
    no_region = "The fallback cast has no measured region samples"

    if not include_cavities:
        records["cavities"] = {"status": "skipped_by_request"}
    else:
        cavities = _write_cavity_outputs(frame, root, prefix, manifest, dpi=dpi,
                                         include_hydrogen=include_hydrogen, include_hetero=include_hetero)
        records["cavities"] = {"status": "completed", "cavity_count": len(cavities), "obstacles": obstacles,
                               "definition": "Separate enclosed-cavity search (detect_cavities defaults); independent of the rolling-probe cast"}

    if not include_tunnels:
        records["tunnels"] = {"status": "skipped_by_request"}
    else:
        start = None
        start_definition = "Structure centroid (find_tunnels default); " + no_region.lower()
        if region_xyz is not None:
            from .rolling import _SphereQueries
            from .radii import atom_vdw_radius
            query = _SphereQueries(np.asarray([a.coord for a in frame.atoms]),
                                   np.asarray([atom_vdw_radius(a) for a in frame.atoms]))
            clearance, _ = query.points(region_xyz)
            start = tuple(float(v) for v in region_xyz[int(np.argmax(clearance))])
            start_definition = "Maximum atom-surface-clearance sample of the measured fallback cast"
        record = {"start": list(start) if start is not None else None, "start_definition": start_definition,
                  "obstacles": obstacles}
        try:
            tunnels = _write_tunnel_outputs(frame, root, prefix, manifest, dpi=dpi, start=start,
                                            include_hydrogen=include_hydrogen, include_hetero=include_hetero)
        except ValueError as exc:
            record.update(status="not_completed", reason=str(exc))
        else:
            record.update(status="completed", tunnel_count=len(tunnels))
        records["tunnels"] = record

    if not include_network:
        records["network"] = {"status": "skipped_by_request"}
    elif region_xyz is None:
        records["network"] = {"status": "not_completed", "reason": no_region}
    else:
        from .models import ChannelPoint, VoidComponent
        volume = float(len(region_xyz) * abs(np.linalg.det(deltas)))
        points = tuple(ChannelPoint(i, tuple(map(float, p)), 0.0, 0.0, float(p[2])) for i, p in enumerate(region_xyz))
        region = VoidComponent(1, "measured_void", tuple(map(float, region_xyz.mean(0))), volume, points)
        coords = tuple(p.position for p in points)
        contacts = annotate_residues(frame, coords, identify_bottleneck=False)
        network = _write_network_outputs(frame, region, contacts, root, prefix, manifest,
                                         dpi=dpi, residue_groups=residue_groups)
        records["network"] = {"status": "completed", "node_count": len(network.nodes), "edge_count": len(network.edges),
                              "region_source": str(source_dx), "region_sample_count": len(points),
                              "definition": "Standard publish network; the region node is every measured fallback-cast sample centre (point-to-atom-surface contacts), not a profile"}
    return records



@radii_option
@annotate_option
def write_rolling_probe_bundle(frame, *, structure_path, output_dir, prefix="crevice",
                               spacing=.5, probe_radius=.8, outer_radius=6.,
                               enclosure_fraction=.9, min_depth=0.,
                               min_component_volume=50., max_components=1,
                               selection_mode="dominant", core_fraction=.08, surface_smoothing=None,
                               max_grid_points=16_000_000, selection="protein", dpi=250, seed=None, seed_tolerance=.75,
                               smooth=None, annotate=None,
                               radii: RadiusSet | str | None = None):
    """Write the 3D rolling-probe cast bundle of ``crevice cast``.

    A probe of radius ``probe_radius`` is rolled over the selected atoms on a
    grid of ``spacing``; grid points it can reach inside an outer envelope
    (``outer_radius`` probe) and buried by at least ``enclosure_fraction`` form
    cavity cores, which are then grown into atom-clear regions
    (:func:`crevice.rolling.rolling_probe_cast`). No single pore profile is
    assigned. Tables, grids and native viewer scenes (one per region and
    combined, plus an end-on view) are written; no Matplotlib figure is
    written, because the native viewer renders are the quality view of a cast.

    Parameters
    ----------
    frame : StructureFrame
        Parsed structure.
    structure_path : str or Path
        Source file, copied into the bundle.
    output_dir : str or Path
        Bundle directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    spacing : float, default 0.5
        Grid spacing (Å); never silently coarsened.
    probe_radius : float, default 0.8
        Inner rolling-probe radius (Å).
    outer_radius : float, default 6.0
        Outer envelope probe radius (Å).
    enclosure_fraction : float, default 0.9
        Burial threshold for cavity cores.
    min_depth : float, default 0.0
        Extra depth inside the outer envelope (Å).
    min_component_volume : float, default 50.0
        Minimum region volume (Å³).
    max_components : int, default 1
        Maximum qualifying regions.
    selection_mode : {"dominant", "all"}, default "dominant"
        Substantial interior cores, or legacy enumeration of every void.
    core_fraction : float, default 0.08
        Minimum core volume relative to the strongest core (dominant mode).
    surface_smoothing : float, optional
        Legacy display interpolation width (Å); cannot be combined with ``smooth``.
    max_grid_points : int, default 16000000
        Grid allocation guard.
    selection : {"protein", "all"}, default "protein"
        Obstacles: protein heavy atoms, or all atoms including HETATM records.
    dpi : int, default 250
        Accepted for a uniform bundle signature; this bundle writes no
        Matplotlib figure (the viewer scenes are the visual output of a cast).
    seed : tuple of float, optional
        Interior seed (Å); keep only its probe-connected region.
    seed_tolerance : float, default 0.75
        Maximum seed-to-grid snap distance (Å).
    smooth : float, optional
        Display-surface smoothing length (Å, default 0.4).
    annotate : bool, optional
        Accepted for a uniform signature; the viewer scenes of this bundle
        contain no text either way. ``None`` (default) inherits the surrounding
        setting.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call (profiles, casts, cavities, tunnels, networks,
        hydration SASA) comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    dict of str to str
        Manifest mapping each output key to its path, including ``manifest_json``.

    Outputs
    -------
    PREFIX_manifest.json : JSON
        File index, atom selection, cast summary and interpretation.
    PREFIX_void_cast.csv : table
        One row per retained cast region: ``region_id``, ``label``,
        ``volume_A3``, ``point_count``, ``centroid_{x,y,z}_A``,
        ``max_clearance_A`` (see :func:`crevice.io.write_void_cast_csv`).
    PREFIX_void_cast.json, PREFIX_selection.json : JSON
        Full cast record (components, settings, dominance diagnostics) read by
        downstream tools, and the atom-selection report.
    PREFIX_region_annotations.json : JSON
        Display-region labels used by the viewer scenes (multi-region casts only).
    PREFIX_viewer.pdb, PREFIX_structure.<ext> : structure
        Renumbered viewer PDB and a copy of the input.
    PREFIX_volume.dx, PREFIX_display.dx, regions/region_NNN.dx : grids
        Measured binary grid, display field and each region alone.
    PREFIX_volume.{pml,vmd,tcl,cxc}, PREFIX_pore_cast.{pml,vmd,tcl,cxc} : scenes
        Native scenes (see :func:`crevice.volume_export.write_volume_viewer_bundle`).
    PREFIX_end_view.{pml,vmd,tcl,cxc} : scenes
        The same scene turned about x to look along the cast.

    Colours and representation
    --------------------------
    Scenes: grey-blue protein cartoon (25% transparent); cast regions as
    opaque surfaces in violet (``#a855f7``, the non-channel cast colour
    :data:`crevice.presentation.NON_CHANNEL_CAST_RGB`), purple, blue and mauve
    (orange for a shared core, translucent grey for other voids); white
    background; no text.
    """
    from .rolling import rolling_probe_cast, select_cast_atoms, chain_enclosed_regions
    from .presentation import write_viewer_structure
    from .volume_export import write_volume_viewer_bundle, write_void_cast_dx
    if Path(prefix).name != prefix or not prefix or any(c in prefix for c in '\n\r"'):
        raise ValueError("prefix must be a simple filename stem")
    from .volume_export import resolve_display_smoothing
    resolve_display_smoothing(smooth, surface_smoothing)
    frame, intake = select_cast_atoms(frame, selection)
    cast = rolling_probe_cast(frame, spacing=spacing, probe_radius=probe_radius,
                              outer_radius=outer_radius, enclosure_fraction=enclosure_fraction,
                              min_depth=min_depth, min_component_volume=min_component_volume,
                              max_components=max_components, max_grid_points=max_grid_points,
                              include_hetero=selection == 'all', selection_mode=selection_mode, core_fraction=core_fraction,
                              seed=seed,seed_tolerance=seed_tolerance)
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    original = _copy_structure_for_bundle(structure_path, root, prefix)
    viewer = write_viewer_structure(frame, root / f"{prefix}_viewer.pdb")
    manifest = {'structure_file': str(viewer), 'original_structure_file': str(original),
                'viewer_identity_json': str(viewer.with_suffix('.identities.json')),
                'void_cast_json': str(root / f'{prefix}_void_cast.json'),
                'void_cast_csv': str(root / f'{prefix}_void_cast.csv'),
                'selection_json': str(root / f'{prefix}_selection.json')}
    write_void_cast_json(cast, manifest['void_cast_json'])
    _write_json(intake, manifest['selection_json'])
    labels = {}
    if cast.components:
        if selection_mode == 'dominant':
            rows = {r['region_id']:r for r in cast.metadata['dominance'].get('regions', [])}
            regions, region_summary = {}, []
            for component in cast.components:
                row = rows[component.metadata['source_component_id']]
                name = 'shared' if row['kind'] == 'shared_interface_core' else f'cavity_{component.id}'
                regions[name] = component.points
                labels[component.id] = name
                region_summary.append({'label':name, 'volume_A3':component.volume, 'selection':row})
            if len(regions) == 1:
                regions = {}
        else:
            regions, region_summary = chain_enclosed_regions(frame, cast)
        manifest['region_annotations_json'] = str(root/f'{prefix}_region_annotations.json')
        _write_json({'regions':region_summary, 'interpretation':'Disjoint display annotation; does not change connected components or imply a biological channel count.'}, manifest['region_annotations_json'])
        manifest.update(write_volume_viewer_bundle(cast, structure_path=viewer, output_dir=root,
                                                   prefix=prefix, frame=frame, cast_extension=0, regions=regions, surface_smoothing=surface_smoothing, smooth=smooth))
        for extension, key in [('pml','volume_pml'), ('tcl','volume_tcl'), ('vmd','volume_vmd'), ('cxc','volume_cxc')]:
            alias = root / f'{prefix}_pore_cast.{extension}'
            shutil.copyfile(manifest[key], alias)
            manifest[f'pore_cast_{extension}'] = str(alias)
        end_turn = json.loads(Path(manifest['scene_json']).read_text()).get('end_view_turn_degrees',74)
        for extension,command in [('pml',f'turn x, {end_turn:g}'),('tcl',f'rotate x by {end_turn:g}'),('vmd',f'rotate x by {end_turn:g}'),('cxc',f'turn x {end_turn:g}')]:
            end_view = root/f'{prefix}_end_view.{extension}'
            end_view.write_text(Path(manifest[f'volume_{extension}']).read_text()+command+'\n')
            manifest[f'end_view_{extension}'] = str(end_view)
        # Each retained region can also be inspected independently in a viewer.
        region_dir = root / 'regions'
        region_dir.mkdir(exist_ok=True)
        for component in cast.components:
            write_void_cast_dx(cast, region_dir/f'region_{component.id:03d}.dx', component_id=component.id)
    write_void_cast_csv(cast, manifest['void_cast_csv'], labels=labels)
    manifest_path = root / f'{prefix}_manifest.json'
    _write_json({'files': manifest, 'prefix': prefix, 'atom_selection': intake,
                 'void_cast': cast.to_dict(), **radii_fields(),
                 'channel_profile_status': 'not_assigned_by_3D_void_detection',
                 'interpretation': 'Geometric void regions. Component count is not a biological channel count.',
                 'native_viewers': 'Scripts generated; running a viewer is a separate validation step.'}, manifest_path)
    manifest['manifest_json'] = str(manifest_path)
    return manifest


def _pdb_resi_for_min_radius(points: Sequence) -> int | None:
    if not points:
        return None
    return min(enumerate(points, start=1), key=lambda item: item[1].radius)[0]

def _pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("Publication figures require matplotlib. Matplotlib is a required CREVICE dependency; reinstall CREVICE or install matplotlib.") from exc
    return plt


def _save_figure(fig, path: Path, *, dpi: int, text: FigureText | None = None) -> None:
    """Save a figure with its description embedded in the file metadata.

    PNG files receive ``Software``, ``Title``, ``Description`` and
    ``CREVICE annotations`` text chunks (read them back with
    :func:`crevice.presentation.figure_description`); PDF and SVG files receive
    the title and description in their own metadata fields.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    described = text.metadata() if text is not None else {}
    if suffix == ".png":
        metadata = {"Software": "CREVICE", **described}
    elif suffix == ".pdf" and described:
        metadata = {"Creator": "CREVICE", "Title": described["Title"], "Subject": described["Description"]}
    elif suffix == ".svg" and described:
        metadata = {"Creator": "CREVICE", "Title": described["Title"], "Description": described["Description"]}
    else:
        metadata = {"Creator": "CREVICE"}
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white", metadata=metadata)
    fig.clf()
    _pyplot().close(fig)


def _style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)
    ax.title.set_fontsize(11)
    ax.xaxis.label.set_fontsize(10)
    ax.yaxis.label.set_fontsize(10)


def _empty_panel(ax, message: str, text: FigureText | None = None) -> None:
    # A no-data notice is always drawn: a blank axis would otherwise be ambiguous.
    if text is not None:
        text.notice(ax, 0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color=CREVICE_GRAY)
    else:
        ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color=CREVICE_GRAY)
    ax.set_xticks([])
    ax.set_yticks([])


def _role_color(role: str) -> str:
    if role == "bottleneck":
        return CREVICE_RED
    if role == "bottleneck-nearby":
        return CREVICE_ORANGE
    if role == "lining":
        return CREVICE_BLUE
    return CREVICE_GRAY


def _role_key(role: str) -> str:
    """Legend category of a contact role; unknown roles share the grey "nearby" colour."""
    return role if role in {"bottleneck", "bottleneck-nearby", "lining"} else "nearby"


# Contact roles of crevice.analysis._contact_role, in legend order. "Gap" is the
# residue's smallest atom-surface distance to the region surface.
ROLE_LEGEND = (
    ("bottleneck", "Bottleneck (gap ≤ 1 Å at the bottleneck)"),
    ("bottleneck-nearby", "Bottleneck-nearby (gap > 1 Å at the bottleneck)"),
    ("lining", "Lining (gap ≤ 1 Å elsewhere)"),
    ("nearby", "Nearby (gap > 1 Å, within cutoff)"),
)

# Edge interaction classes written by crevice.networks, in legend order.
INTERACTION_LEGEND = (
    ("region", "Channel/cavity contact"),
    ("salt_bridge_candidate", "Salt-bridge candidate"),
    ("oppositely_charged_residue_contact", "Oppositely charged"),
    ("hydrophobic_residue_contact", "Hydrophobic"),
    ("polar_residue_contact", "Polar"),
    ("distance_contact", "Other distance contact"),
)


def _interaction_key(interaction: str) -> str:
    """Legend category of a network edge interaction class."""
    if interaction.startswith("region_"):
        return "region"
    known = {key for key, _label in INTERACTION_LEGEND}
    return interaction if interaction in known else "distance_contact"


def _interaction_color(interaction: str) -> str:
    """Colour of an edge interaction class, keyed on the names crevice.networks writes."""
    return {
        "region": CREVICE_RED,
        "salt_bridge_candidate": CREVICE_PURPLE,
        "oppositely_charged_residue_contact": CREVICE_ORANGE,
        "hydrophobic_residue_contact": CREVICE_GREEN,
        "polar_residue_contact": CREVICE_TEAL,
    }.get(_interaction_key(interaction), CREVICE_GRAY)


def _contact_criterion(network: ResidueInteractionNetwork) -> str:
    """Plain-language residue contact rule of a network, from its metadata."""
    cutoff = network.metadata.get("cutoff")
    metric = network.metadata.get("distance_metric", "vdw_surface_gap")
    distance = "atom-centre distance" if metric in {"center", "atom_center_distance"} else "atom-surface gap"
    return f"{distance} ≤ {cutoff:g} Å" if isinstance(cutoff, (int, float)) else distance


def _landmark_color(label: str) -> str:
    residue = label.split(":")[-1][:3].upper()
    if residue in {"ASP", "GLU", "LYS", "ARG", "HIS"}:
        return CREVICE_PURPLE
    if residue in {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}:
        return CREVICE_TEAL
    return CREVICE_BLUE


def _profile_residue_landmarks(
    profile: PoreProfile,
    *,
    label_every: int,
    max_labels: int,
):
    landmarks = []
    last_label = None
    for index, point in enumerate(profile.points):
        label = point.nearest_residue
        if not label:
            continue
        if index != profile.bottleneck.index and index % label_every != 0:
            continue
        if label == last_label and index != profile.bottleneck.index:
            continue
        landmarks.append((point, label))
        last_label = label
    if not any(point.index == profile.bottleneck.index for point, _label in landmarks):
        landmarks.append((profile.bottleneck, profile.bottleneck.nearest_residue or "bottleneck"))
    landmarks.sort(key=lambda item: item[0].t)
    if len(landmarks) <= max_labels:
        return landmarks
    keep = {0, len(landmarks) - 1}
    bottleneck_pos = next((idx for idx, item in enumerate(landmarks) if item[0].index == profile.bottleneck.index), None)
    if bottleneck_pos is not None:
        keep.add(bottleneck_pos)
    slots = max(max_labels - len(keep), 0)
    if slots:
        step = max(1, len(landmarks) // slots)
        keep.update(range(0, len(landmarks), step))
    return [item for idx, item in enumerate(landmarks) if idx in sorted(keep)][:max_labels]


def _node_angles(selected: Sequence[str], nodes_by_id) -> dict[str, float]:
    region = [node_id for node_id in selected if nodes_by_id[node_id].kind == "region"]
    residues = sorted([node_id for node_id in selected if nodes_by_id[node_id].kind != "region"])
    ordered = region + residues
    return {node_id: (2.0 * math.pi * index / len(ordered)) + math.pi / 2.0 for index, node_id in enumerate(ordered)}


def _short_residue_label(label: str) -> str:
    if ":" not in label:
        return label.replace("region:", "")
    chain, residue = label.split(":", 1)
    return f"{chain}:{residue}"


def _copy_structure_for_bundle(structure_path: str | Path, root: Path, prefix: str) -> Path:
    source = Path(structure_path)
    suffix = "".join(source.suffixes) or ".pdb"
    destination = root / f"{prefix}_structure{suffix}"
    if source.exists() and source.resolve() != destination.resolve():
        shutil.copyfile(source, destination)
        return destination
    if source.exists():
        return source
    raise FileNotFoundError(f"Cannot copy missing structure file {source}")


def _ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    if output.parent != Path("."):
        output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _write_json(data: object, path: str | Path) -> None:
    output = _ensure_parent(path)
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

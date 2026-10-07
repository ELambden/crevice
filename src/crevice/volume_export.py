"""Measured occupancy maps and shared presentation fields for native viewers.

A void cast (:class:`crevice.models.VoidCast`) is a set of samples on a
regular lattice. :func:`cast_grid` turns it into the measured binary
occupancy grid, written unchanged as OpenDX by :func:`write_void_cast_dx`
(``*_volume.dx``; every volume and boundary measurement uses this grid).
Surfaces shown by PyMOL, VMD and ChimeraX come from a separate *display*
field contoured at 0.5: by default constrained Gaussian twicing on a 2x
supersampled grid (``--smooth X``, default 0.4 Å; see
``SMOOTH_DEFINITION`` and :func:`constrained_display_field`), or the
legacy sample-clamped Gaussian (``--surface-smoothing W``,
:func:`presentation_grid`). Display fields never change a measured sample's
inside/outside state. :func:`write_volume_viewer_bundle` writes the maps and
the three viewer scenes.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

from .models import VoidCast


def cast_grid(cast: VoidCast, *, component_id: int | None = None, max_points: int = 16_000_000):
    """Binary grid from full measured components, with a one-cell empty border.

    The 0.5 isosurface is a voxel-boundary approximation, with resolution set by
    ``cast.spacing``. No solvent probe, Gaussian broadening or smoothing is added.

    Parameters
    ----------
    cast : VoidCast
        Cast whose samples share one regular lattice (``metadata["grid_basis"]``
        rows, identity by default).
    component_id : int, optional
        Only this component; default all components.
    max_points : int, default 16000000
        Largest allowed grid (number of voxels).

    Returns
    -------
    grid : numpy.ndarray of float32
        1 for occupied samples, 0 elsewhere, indexed along the basis axes.
    origin : numpy.ndarray, shape (3,)
        Cartesian position (Å) of voxel (0, 0, 0), one spacing outside the
        occupied samples.

    Raises
    ------
    ValueError
        For an unknown component, an empty cast, a nonpositive spacing, a
        non-orthonormal basis, samples off one lattice, or a grid larger than
        ``max_points``.

    Examples
    --------
    >>> from crevice.models import ChannelPoint, VoidComponent, VoidCast
    >>> from crevice.volume_export import cast_grid
    >>> pts = tuple(ChannelPoint(index=i, position=(0.5 * i, 0.0, 0.0), radius=1.0,
    ...                          raw_clearance=1.0, t=0.0) for i in range(3))
    >>> comp = VoidComponent(id=1, kind="cavity", center=(0.5, 0.0, 0.0), volume=0.375, points=pts)
    >>> grid, origin = cast_grid(VoidCast(mode="cavity", spacing=0.5, min_radius=0.8,
    ...                                   components=(comp,), points=pts))
    >>> grid.shape, int(grid.sum()), origin.tolist()
    ((5, 3, 3), 3, [-0.5, -0.5, -0.5])
    """
    import numpy as np

    components = cast.components
    if component_id is not None:
        components = tuple(c for c in components if c.id == component_id)
        if not components:
            raise ValueError("Unknown component_id")
    points = [p.position for c in components for p in c.points]
    if not points:
        raise ValueError("Cannot render an empty cast")
    if not math.isfinite(cast.spacing) or cast.spacing <= 0:
        raise ValueError("Cast spacing must be finite and positive")
    coords = np.asarray(points)
    basis = np.asarray(cast.metadata.get("grid_basis", np.eye(3)), dtype=float)
    if basis.shape != (3, 3) or not np.allclose(basis @ basis.T, np.eye(3), atol=1e-8):
        raise ValueError("Cast grid_basis must be orthonormal")
    local = coords @ basis.T
    local_origin = local.min(axis=0) - cast.spacing
    origin = local_origin @ basis
    indices_float = (local - local_origin) / cast.spacing
    if not np.allclose(indices_float, np.rint(indices_float), atol=1e-5, rtol=0):
        raise ValueError("Cast points do not share a regular grid")
    indices = np.rint(indices_float).astype(int)
    shape = tuple(int(i) + 2 for i in indices.max(axis=0))
    if math.prod(shape) > max_points:
        raise ValueError("Volume map exceeds max_points; select a smaller region or increase the limit")
    grid = np.zeros(shape, dtype=np.float32)
    grid[tuple(indices.T)] = 1
    return grid, origin


def presentation_grid(grid, spacing, smoothing_A=.6):
    """Smooth interpolation while preserving every occupied/unoccupied sample.

    The legacy display method (``--surface-smoothing W``). Clamping around 0.5
    prevents a blur from filling an empty sample or deleting a narrow
    connection. Measurements always use the original binary map. Surface
    vertices remain within the original grid cells; this is display
    interpolation, not a new continuous atom-exclusion guarantee.

    Parameters
    ----------
    grid : numpy.ndarray
        Binary occupancy grid (:func:`cast_grid`).
    spacing : float
        Grid spacing, Å.
    smoothing_A : float, default 0.6
        Gaussian standard deviation, Å; 0 returns a copy of ``grid``.

    Returns
    -------
    numpy.ndarray of float32
        Gaussian-filtered field, set to at least 0.501 on occupied samples and at
        most 0.499 on empty ones.

    Raises
    ------
    ValueError
        If ``smoothing_A`` is negative or not finite.
    """
    import numpy as np
    if not math.isfinite(smoothing_A) or smoothing_A < 0:
        raise ValueError('surface_smoothing must be finite and non-negative')
    if smoothing_A == 0:
        return grid.copy()
    from scipy.ndimage import gaussian_filter
    smooth = gaussian_filter(grid, sigma=smoothing_A/spacing, mode='constant')
    return np.where(grid > .5, np.maximum(smooth, .501), np.minimum(smooth, .499)).astype(np.float32)


#: Default ``smooth`` width in Angstrom used when neither ``smooth`` nor the
#: legacy ``surface_smoothing`` option is supplied.
DEFAULT_SMOOTH_A = .4
#: Display-grid budget for the supersampled smoothing field.
DEFAULT_MAX_DISPLAY_POINTS = 8_000_000
SMOOTH_METHOD = 'constrained_gaussian_twicing'
SMOOTH_DEFINITION = (
    'Display only. The measured binary occupancy grid B (spacing h) is linearly '
    'resampled onto a display grid of spacing h/k (default k=2) that contains every '
    'measured sample, then smoothed by Gaussian twicing, G*(2B-G*B), where G is a '
    'Gaussian with standard deviation X Angstrom; twicing removes the leading-order '
    'shrinkage of a single Gaussian. The field is contoured at 0.5. Every display '
    'node lying on a measured sample, or inside a '
    'measured grid cell/face/edge whose corner samples all agree, is clamped to '
    'that occupancy (>=0.501 occupied, <=0.499 empty). The surface therefore '
    'moves only inside grid cells that contain the unsmoothed voxel boundary, every '
    'measured sample keeps its inside/outside classification, face-adjacent '
    'occupied samples stay connected and empty samples between regions stay empty. '
    'X=0 shows the raw 0.5 voxel boundary. Measured grids, volumes and profiles '
    'are unchanged.')
LEGACY_SMOOTH_DEFINITION = ('Gaussian display interpolation clamped to preserve every occupancy '
                            'sign at isovalue 0.5; binary measured DX unchanged')


def validate_display_width(value, name: str = 'smooth') -> float:
    """Return a finite non-negative display-smoothing length in Å.

    Parameters
    ----------
    value : float or str
        Length to check (booleans are refused).
    name : str, default "smooth"
        Option name used in the error message.

    Returns
    -------
    float
        The length in Å.

    Raises
    ------
    ValueError
        If ``value`` is not a finite non-negative number.

    Examples
    --------
    >>> from crevice.volume_export import validate_display_width
    >>> validate_display_width("0.4")
    0.4
    """
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a finite non-negative length in Angstrom')
    try:
        width = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{name} must be a finite non-negative length in Angstrom') from None
    if not math.isfinite(width) or width < 0:
        raise ValueError(f'{name} must be a finite non-negative length in Angstrom')
    return width


class DisplaySmoothing:
    """Resolved display-surface smoothing; never used by measurements.

    Create it with :func:`resolve_display_smoothing`.

    Parameters
    ----------
    method : {"constrained_gaussian_twicing", "legacy_clamped_gaussian", "none"}
        ``constrained_gaussian_twicing`` (``smooth``), ``legacy_clamped_gaussian``
        (``surface_smoothing``) or ``none``.
    width_A : float
        Gaussian standard deviation X, Å.
    option : str
        How it was selected (``"default"``, ``"smooth"`` or
        ``"surface_smoothing"``), recorded in the metadata.
    supersample : int, optional
        Display supersampling factor k from 1 to 8 (default 2) for the
        constrained method.
    max_display_points : int, default 8000000
        Display-grid budget; k is reduced until the display grid fits.

    Raises
    ------
    ValueError
        For an unknown method or an invalid ``supersample``.
    """

    __slots__ = ('method', 'width_A', 'option', 'supersample', 'max_display_points')

    def __init__(self, method, width_A, option, supersample=None,
                 max_display_points=DEFAULT_MAX_DISPLAY_POINTS):
        if method not in {SMOOTH_METHOD, 'legacy_clamped_gaussian', 'none'}:
            raise ValueError('Unknown display smoothing method')
        if supersample is not None and (isinstance(supersample, bool) or int(supersample) != supersample
                                        or not 1 <= int(supersample) <= 8):
            raise ValueError('smooth_supersample must be an integer from 1 to 8')
        self.method, self.width_A, self.option = method, float(width_A), option
        self.supersample = None if supersample is None else int(supersample)
        self.max_display_points = int(max_display_points)

    def __repr__(self):
        return f'DisplaySmoothing({self.method!r}, {self.width_A!r}, option={self.option!r})'

    def factor(self, shape, spacing) -> int:
        """Display supersampling factor k (display spacing = spacing/k).

        Parameters
        ----------
        shape : tuple of int
            Shape of the measured grid.
        spacing : float
            Measured grid spacing, Å (unused; kept for a uniform signature).

        Returns
        -------
        int
            1 unless the constrained method with a nonzero width is selected; then
            the requested factor (default 2), reduced until
            ``prod(k*(n-1)+1) <= max_display_points``.
        """
        if self.method != SMOOTH_METHOD or self.width_A == 0:
            return 1
        k = self.supersample or 2
        while k > 1 and math.prod(k*(n-1)+1 for n in shape) > self.max_display_points:
            k -= 1
        return k

    def field(self, grid, spacing):
        """Return ``(field, display_spacing, k)`` contoured at 0.5 by every viewer.

        Parameters
        ----------
        grid : numpy.ndarray
            Measured binary occupancy grid.
        spacing : float
            Measured grid spacing, Å.

        Returns
        -------
        field : numpy.ndarray of float32
            Display field (a copy of ``grid`` for ``none`` or zero width).
        display_spacing : float
            Spacing of ``field``, Å (``spacing / k``).
        k : int
            Supersampling factor used.
        """
        if self.method == 'legacy_clamped_gaussian':
            return presentation_grid(grid, spacing, self.width_A), spacing, 1
        if self.method == 'none' or self.width_A == 0:
            return grid.astype('float32', copy=True), spacing, 1
        k = self.factor(grid.shape, spacing)
        return constrained_display_field(grid, spacing, self.width_A, k), spacing/k, k

    def metadata(self, spacing, shape=None) -> dict:
        """Describe the display smoothing for manifests and scene JSON.

        Parameters
        ----------
        spacing : float
            Measured grid spacing, Å.
        shape : tuple of int, optional
            Measured grid shape; when given the supersampling factor and display
            spacing are reported, else they are ``None``.

        Returns
        -------
        dict
            ``method``, ``width_A``, ``gaussian_sigma_A``, ``selected_by``,
            ``measured_spacing_A``, ``supersample_factor``, ``display_spacing_A``,
            ``isovalue`` (0.5), ``definition`` and ``measurement_effect`` (none).
        """
        k = self.factor(shape, spacing) if shape is not None else None
        definition = {SMOOTH_METHOD: SMOOTH_DEFINITION, 'legacy_clamped_gaussian': LEGACY_SMOOTH_DEFINITION,
                      'none': 'Raw 0.5 isosurface of the measured binary occupancy grid'}[self.method]
        return {'method': self.method if self.width_A else 'none', 'width_A': self.width_A,
                'gaussian_sigma_A': self.width_A if self.width_A else 0.0,
                'selected_by': self.option, 'measured_spacing_A': spacing,
                'supersample_factor': k, 'display_spacing_A': None if k is None else spacing/k,
                'isovalue': .5, 'definition': definition,
                'measurement_effect': 'none; measured binary DX, volumes, profiles and attribution use the unsmoothed grid'}


def resolve_display_smoothing(smooth=None, surface_smoothing=None, *, smooth_supersample=None,
                              max_display_points=DEFAULT_MAX_DISPLAY_POINTS) -> DisplaySmoothing:
    """Resolve ``smooth``/legacy ``surface_smoothing`` into one display method.

    Supplying both is an error. ``surface_smoothing`` keeps the previous
    sample-clamped Gaussian on the measured grid. Otherwise ``smooth`` (default
    ``DEFAULT_SMOOTH_A`` = 0.4 Å) selects the constrained supersampled Gaussian; 0
    disables smoothing.

    Parameters
    ----------
    smooth : float or DisplaySmoothing, optional
        Display smoothing length X, Å (``--smooth``); an existing
        :class:`DisplaySmoothing` is returned unchanged.
    surface_smoothing : float, optional
        Legacy width W, Å (``--surface-smoothing``).
    smooth_supersample : int, optional
        Supersampling factor for ``smooth`` (default 2).
    max_display_points : int, default 8000000
        Display-grid budget.

    Returns
    -------
    DisplaySmoothing

    Raises
    ------
    ValueError
        If both options are given or a width is invalid.

    Examples
    --------
    >>> from crevice.volume_export import resolve_display_smoothing
    >>> resolve_display_smoothing()
    DisplaySmoothing('constrained_gaussian_twicing', 0.4, option='default')
    >>> resolve_display_smoothing(surface_smoothing=0.6).method
    'legacy_clamped_gaussian'
    """
    if isinstance(smooth, DisplaySmoothing):
        if surface_smoothing is not None:
            raise ValueError('Specify either smooth or surface_smoothing, not both')
        return smooth
    if smooth is not None and surface_smoothing is not None:
        raise ValueError('Specify either smooth or surface_smoothing, not both')
    if surface_smoothing is not None:
        width = validate_display_width(surface_smoothing, 'surface_smoothing')
        return DisplaySmoothing('legacy_clamped_gaussian' if width else 'none', width, 'surface_smoothing')
    width = validate_display_width(DEFAULT_SMOOTH_A if smooth is None else smooth, 'smooth')
    return DisplaySmoothing(SMOOTH_METHOD if width else 'none', width,
                            'default' if smooth is None else 'smooth', smooth_supersample, max_display_points)


def _upsample_linear(grid, k: int):
    """Separable linear resampling; node k*i equals measured sample i exactly."""
    import numpy as np
    out = np.asarray(grid, dtype=np.float32)
    for axis in range(3):
        n = out.shape[axis]
        position = np.arange(k*(n-1)+1)/k
        low = np.minimum(np.floor(position).astype(int), n-1)
        high = np.minimum(low+1, n-1)
        weight = (position-low).astype(np.float32)
        shape = [1, 1, 1]
        shape[axis] = -1
        weight = weight.reshape(shape)
        out = np.take(out, low, axis=axis)*(1-weight) + np.take(out, high, axis=axis)*weight
    return out


def _cell_agreement(grid, k: int):
    """Min/max measured occupancy over the cell, face, edge or sample containing each display node."""
    import numpy as np
    index = []
    for n in grid.shape:
        f = np.arange(k*(n-1)+1)
        low = f//k
        index.append((low, np.minimum(low+(f % k != 0), n-1)))
    lo = hi = None
    for choice in itertools.product((0, 1), repeat=3):
        values = grid[np.ix_(*(index[a][c] for a, c in enumerate(choice)))]
        lo = values if lo is None else np.minimum(lo, values)
        hi = values if hi is None else np.maximum(hi, values)
    return lo, hi


def constrained_display_field(grid, spacing, smooth_A, k: int = 2):
    """Constrained supersampled Gaussian display field (see ``SMOOTH_DEFINITION``).

    1. Linearly resample the binary grid onto spacing ``spacing/k`` (display node
       ``k*i`` equals measured sample ``i`` exactly).
    2. Gaussian twicing ``G*(2B - G*B)`` with standard deviation ``smooth_A``,
       which removes the leading-order shrinkage of a single blur.
    3. Clamp every display node inside a measured cell, face, edge or sample
       whose corner samples all agree: >= 0.501 if all occupied, <= 0.499 if all
       empty. The 0.5 surface therefore moves only inside cells that contain the
       unsmoothed voxel boundary.
    4. Round to four decimals; nodes exactly at 0.5 move 1e-4 towards their
       unrounded side so contours never pass through a node.

    Parameters
    ----------
    grid : array_like
        Binary 3D occupancy grid.
    spacing : float
        Measured grid spacing, Å.
    smooth_A : float
        Gaussian standard deviation X, Å; 0 returns a copy of ``grid``.
    k : int, default 2
        Supersampling factor.

    Returns
    -------
    numpy.ndarray of float32
        Field with spacing ``spacing/k`` and shape ``k*(n-1)+1`` per axis,
        sharing the measured origin and basis.

    Raises
    ------
    ValueError
        For a negative width, a nonpositive spacing, an invalid ``k`` or a grid
        that is not binary and 3D.
    """
    import numpy as np
    width = validate_display_width(smooth_A, 'smooth')
    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError('spacing must be finite and positive')
    if isinstance(k, bool) or int(k) != k or k < 1:
        raise ValueError('supersample factor must be a positive integer')
    k = int(k)
    grid = np.asarray(grid, dtype=np.float32)
    if grid.ndim != 3 or not np.isin(grid, (0, 1)).all():
        raise ValueError('constrained smoothing requires a binary 3D occupancy grid')
    if width == 0:
        return grid.copy()
    from scipy.ndimage import gaussian_filter
    fine = _upsample_linear(grid, k)
    sigma = width*k/spacing
    # Twicing G*(2B-G*B): unit gain at low spatial frequency to second order,
    # so constrictions and convex parts do not shrink like a single blur.
    field = gaussian_filter(2*fine - gaussian_filter(fine, sigma, mode='constant', truncate=4.),
                            sigma, mode='constant', truncate=4.)
    lo, hi = _cell_agreement(grid, k)
    field = np.where(lo > .5, np.maximum(field, .501), np.where(hi < .5, np.minimum(field, .499), field))
    # Four decimals keep portable DX files compact; every viewer and the PyMOL
    # mesh contour this same rounded field. A node exactly at the 0.5 level
    # would give degenerate/open contour triangles, so ties move 1e-4 towards
    # the unrounded side. Clamped nodes (0.501/0.499) are unaffected.
    rounded = np.round(field, 4)
    tie = rounded == .5
    rounded[tie] = np.where(field[tie] >= .5, .5001, .4999)
    return rounded.astype(np.float32)


def write_void_cast_dx(cast: VoidCast, path, *, component_id: int | None = None,
                       surface_smoothing: float = 0., smoothing: DisplaySmoothing | None = None) -> Path:
    """Write OpenDX scalar data in x/y/z order (z varies fastest).

    Without ``smoothing`` and with ``surface_smoothing=0`` this is the measured
    binary grid. A ``DisplaySmoothing`` writes its display field instead. The
    grid axes are the cast's lattice axes (``delta`` lines), so rotated lattices
    are preserved; the header comment states whether the map is binary or a
    presentation field. Read it back with
    :func:`crevice.residue_evidence.read_binary_dx`.

    Parameters
    ----------
    cast : VoidCast
        Cast to export.
    path : str or Path
        Output ``.dx`` path (parent directories are created).
    component_id : int, optional
        Only this component; default all.
    surface_smoothing : float, default 0.0
        Legacy display width, Å, used only when ``smoothing`` is ``None``.
    smoothing : DisplaySmoothing, optional
        Display method whose field is written instead of the binary grid.

    Returns
    -------
    Path
        The written file.
    """
    grid, origin = cast_grid(cast, component_id=component_id)
    if smoothing is None:
        grid, spacing = presentation_grid(grid, cast.spacing, surface_smoothing), cast.spacing
        shown = bool(surface_smoothing)
    else:
        grid, spacing, _ = smoothing.field(grid, cast.spacing)
        shown = smoothing.width_A > 0
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = " ".join(str(n) for n in grid.shape)
    with output.open("w") as handle:
        handle.write("# CREVICE selected-void "+("presentation" if shown else "binary")+" grid; isovalue 0.5; units Angstrom\n")
        handle.write(f"object 1 class gridpositions counts {counts}\n")
        handle.write("origin " + " ".join(f"{v:.9g}" for v in origin) + "\n")
        basis = cast.metadata.get("grid_basis", ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
        for direction in basis:
            handle.write("delta " + " ".join(f"{spacing*v:.12g}" for v in direction) + "\n")
        handle.write(f"object 2 class gridconnections counts {counts}\n")
        handle.write(f"object 3 class array type double rank 0 items {grid.size} data follows\n")
        values = grid.ravel(order="C")
        for i in range(0, len(values), 3):
            handle.write(" ".join(f"{v:.6g}" for v in values[i:i + 3]) + "\n")
        handle.write('attribute "dep" string "positions"\nobject "crevice_cast" class field\n')
        for name, number in (("positions", 1), ("connections", 2), ("data", 3)):
            handle.write(f'component "{name}" value {number}\n')
    return output


def _tcl_word(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"').replace("$", "\\$").replace("[", "\\[") + '"'


def write_volume_viewer_bundle(cast: VoidCast, *, structure_path, output_dir,
                              prefix: str = "crevice", frame=None, profile=None,
                              cast_extension: float = 2.0, regions=None, surface_smoothing: float | None = None,
                              smooth: float | None = None, smooth_supersample: int | None = None,
                              exit_casts=None, mouth_guides: bool | None = None,
                              exit_centre_lines: bool | None = None, cast_color=None,
                              segments=None, lining=None) -> dict[str, str]:
    """Write PyMOL, VMD and ChimeraX scenes of a cast inside its protein.

    The measured binary grid is written unchanged to ``PREFIX_volume.dx``. A
    display field, optionally extended beyond channel mouths by
    ``cast_extension`` and smoothed, is written separately and contoured at 0.5
    for every viewer; PyMOL receives the same triangles as a standalone mesh.
    Measured and display-only grids are never mixed. The scenes contain no
    text, labels or legends.

    Parameters
    ----------
    cast : VoidCast
        Measured cast.
    structure_path : str or Path
        Protein structure shown as cartoon (normally the bundle's viewer PDB).
    output_dir : str or Path
        Output directory (created).
    prefix : str, default "crevice"
        File-name prefix; must be a simple file stem.
    frame : StructureFrame, optional
        Structure used for display extensions, camera and reference coordinates.
    profile : PoreProfile, optional
        Channel profile; its mouths give the optional ring guides and its
        lateral exit legs the optional centre lines.
    cast_extension : float, default 2.0
        Display-only continuation beyond each mouth (Å); excluded from volume.
    regions : dict of str to sequence of ChannelPoint, optional
        Disjoint display regions (must partition the cast); each gets its own
        colour and objects.
    surface_smoothing : float, optional
        Legacy sample-clamped Gaussian width (Å) on the measured grid.
    smooth : float, optional
        Display smoothing length (Å, default ``DEFAULT_SMOOTH_A`` = 0.4; 0 = raw
        voxel boundary): constrained Gaussian twicing on a supersampled grid,
        described by ``SMOOTH_DEFINITION``. Supplying both options raises
        ``ValueError``.
    smooth_supersample : int, optional
        Display supersampling factor for ``smooth`` (default 2).
    exit_casts : sequence of VoidComponent, optional
        Lateral exit leg casts on the cast's lattice
        (:func:`crevice.channel_exits.lateral_exit_casts`), written as
        ``PREFIX_exit_casts.dx``. Without ``segments`` each leg is drawn as an
        EXIT segment (``crevice_exit_N``); ``PREFIX_volume.dx`` and the
        measured volume stay those of ``cast``.
    mouth_guides : bool, optional
        Draw the channel-mouth rings. ``None`` (default) follows
        :func:`crevice.presentation.mouth_guides_enabled` (off unless
        ``--mouth-guides``). Ring coordinates are always in the scene JSON.
    exit_centre_lines : bool, optional
        Also draw lateral exit legs as thin centre-line tubes. ``None``
        (default) follows :func:`crevice.presentation.exit_centre_lines_enabled`
        (off unless ``--exit-centre-lines``).
    cast_color : tuple of float, optional
        Cast surface RGB (0-1). Default
        :func:`crevice.presentation.cast_rgb`: teal for a channel cast, violet
        for every other cast.
    segments : sequence of VoidComponent, optional
        ENTRY/EXIT segment casts on the cast's lattice
        (:func:`crevice.cast_segments.segment_casts`; metadata ``segment``).
        When given (even empty) the scene draws ``cast`` as ``crevice_lumen``
        and each segment as its own object instead of ``crevice_volume``;
        display extension nodes beyond a mouth that has a segment are dropped.
        Cannot be combined with ``regions``.
    lining : sequence of dict, optional
        Rows of :func:`crevice.cast_segments.lining_residues` (needs
        ``structure_path`` to be the :func:`crevice.presentation.write_viewer_structure`
        copy of ``frame``); adds hidden lining-residue objects and selections.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_volume.dx : DX
        Measured binary grid (1 inside, 0 outside).
    PREFIX_display.dx, PREFIX_pymol_local.dx, PREFIX_pymol_mesh.npz : display
        Display field in world and orthogonal PyMOL-local frames; PyMOL mesh.
    PREFIX_volume.pml, PREFIX_volume.vmd, PREFIX_volume.tcl, PREFIX_volume.cxc : scenes
        Open from the bundle directory.
    PREFIX_render.cxc : ChimeraX script
        Saves a supersampled PNG.
    PREFIX_scene.json, PREFIX_volume_metadata.json, PREFIX_geometry_reference.npz : data
        Camera and mouth geometry; display definitions; expected surfaces for checks.
    PREFIX_mouths.bild : BILD
        Mouth ring guides and exit centre lines when drawn (otherwise only a
        colour line); opened by the ChimeraX scene only when it has guides.
    PREFIX_exit_casts.dx : DX
        Only with ``exit_casts``: the measured binary grid of the exit leg casts.
    PREFIX_<segment>_cast.dx, PREFIX_<segment>_display.dx, PREFIX_<segment>_display_local.dx : DX
        Segmented scenes only: the measured grid of each entry/exit segment (the
        lumen's is ``PREFIX_volume.dx``) and each object's display field in the
        world and orthogonal local frames.
    PREFIX_region_NNN.dx, PREFIX_region_NNN_local.dx, PREFIX_display_regions.json : regions
        Only with ``regions``.

    Colours and representation
    --------------------------
    * Protein: cartoon only, grey-blue (``[0.56, 0.63, 0.70]``, ChimeraX
      ``#8fa1b3``), 25% transparent.
    * Cast: opaque surface, teal (``[0.08, 0.58, 0.63]``, ``#1494a1``) for a
      channel cast (``cast.mode == "channel"``), violet
      (``[0.659, 0.333, 0.969]``, ``#a855f7``) for every other cast (cavity,
      rolling-probe and dominant-region casts, the unresolved-profile
      fallback); see :data:`crevice.presentation.NON_CHANNEL_CAST_RGB`.
    * Channel cast segments (``segments``, or ``exit_casts`` alone): one opaque
      surface object each, PyMOL object / VMD molecule / ChimeraX model
      ``crevice_entry_N`` yellow (``[0.941, 0.831, 0.227]``, ``#f0d43a``),
      ``crevice_lumen`` in the cast colour (teal) and ``crevice_exit_N`` rust
      (``[0.761, 0.255, 0.047]``, ``#c2410c``); each can be recoloured or
      hidden on its own. Separately smoothed surfaces meet at the mouth planes.
    * Lining residues (``lining``): stick objects ``crevice_<segment>_lining``
      (carbons in the segment colour, heteroatoms by element, radius 0.2 Å),
      disabled by default, and selections ``crevice_<segment>_lining_sel``
      (PyMOL; ``crevice_lining SEGMENT [on|off|toggle]``); VMD Licorice
      representations hidden by default (colour ids 10 lumen, 26 entry, 27
      exit; ``crevice_lining SEGMENT on|off``); ChimeraX named selections
      ``crevice_<segment>_lining``. No labels.
    * Channel mouths: no guide by default. With ``mouth_guides``: orange rings
      of 0.065 Å cylinders (``[0.90, 0.52, 0.12]``), PyMOL object
      ``crevice_mouths``.
    * Lateral exit centre lines: only with ``exit_centre_lines``; blue tubes
      of 0.12 Å cylinders (``[0.16, 0.47, 0.84]``), one per distinct leg;
      PyMOL object ``crevice_exits``, also in ``PREFIX_mouths.bild``.
    * Regions: the cast colour (teal or violet, as above), purple
      ``[0.42, 0.38, 0.72]``, blue ``[0.20, 0.53, 0.86]``,
      mauve ``[0.67, 0.40, 0.62]`` in turn; a shared core orange
      ``[0.90, 0.58, 0.18]``; other voids grey ``[0.55, 0.64, 0.69]`` at 25% opacity.
    * White background, orthographic camera.
    """
    smoothing = resolve_display_smoothing(smooth, surface_smoothing, smooth_supersample=smooth_supersample)
    from .presentation import (cast_rgb, exit_centre_lines_enabled, extend_channel_cast, mouth_guides_enabled,
                               orthogonal_map_cast, scene_geometry)
    draw_rings = mouth_guides_enabled() if mouth_guides is None else bool(mouth_guides)
    draw_lines = exit_centre_lines_enabled() if exit_centre_lines is None else bool(exit_centre_lines)
    color = tuple(float(v) for v in (cast_color if cast_color is not None else cast_rgb(cast)))
    color_hex = '#'+''.join(f'{round(v*255):02x}' for v in color)
    from .channel_exits import exit_cast_summary
    exit_casts = tuple(exit_casts or ())
    if Path(prefix).name != prefix or not prefix or any(c in prefix for c in '\n\r"'):
        raise ValueError("prefix must be a simple filename stem")
    if not math.isfinite(cast_extension) or cast_extension < 0:
        raise ValueError("cast_extension must be finite and non-negative")
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    structure = Path(structure_path).resolve()
    dx = write_void_cast_dx(cast, root / f"{prefix}_volume.dx")
    display = extend_channel_cast(cast, frame, cast_extension) if frame is not None else cast
    channel_display = display
    from dataclasses import replace as _replace
    if segments is None and exit_casts:
        # Older callers pass only the lateral exit legs: each is an EXIT segment.
        segments = tuple(_replace(c, metadata={**c.metadata, "segment": f"exit_{i}", "segment_type": "exit",
                                               "segment_index": i})
                         for i, c in enumerate(exit_casts, 1))
    segment_mode = segments is not None
    segments = tuple(segments or ())
    if segment_mode and regions:
        raise ValueError("Cast segments and display regions cannot be combined")
    if exit_casts:
        exit_cast = _replace(cast, components=exit_casts, points=tuple(p for c in exit_casts for p in c.points))
        manifest_exit_dx = write_void_cast_dx(exit_cast, root / f"{prefix}_exit_casts.dx")
    lumen_display = display
    if segments:
        # Segment nodes beyond a mouth replace any display extension there; the
        # lumen object keeps the rest of the display continuation.
        def _key(point):
            return tuple(round(v, 6) for v in point.position)
        taken = {_key(p) for c in segments for p in c.points}
        ends = {c.metadata.get("end") for c in segments}
        meta = cast.metadata
        if {"grid_origin", "grid_basis", "grid_dimensions"} <= set(meta):
            import numpy as _np
            _origin, _basis = _np.asarray(meta["grid_origin"], float), _np.asarray(meta["grid_basis"], float)
            _last = int(meta["grid_dimensions"][2])-1

            def _beyond_segment_mouth(point):
                k = int(round(float((_np.asarray(point.position)-_origin) @ _basis[2])/cast.spacing))
                return (k < 0 and "lower" in ends) or (k > _last and "upper" in ends)
        else:
            def _beyond_segment_mouth(point):
                return False
        # The display continuation beyond a mouth is dropped where a segment
        # fills that end, so the segment object is not hidden inside it.
        kept = tuple(c for c in (_replace(c, points=tuple(p for p in c.points if _key(p) not in taken
                                                          and not _beyond_segment_mouth(p)))
                                 for c in display.components) if c.points)
        lumen_display = _replace(display, components=kept, points=tuple(p for c in kept for p in c.points))
        display = _replace(display, components=kept+segments,
                           points=lumen_display.points+tuple(p for c in segments for p in c.points),
                           metadata={**display.metadata, "display_only": True})
    shown_dx = write_void_cast_dx(display, root / f"{prefix}_display.dx", smoothing=smoothing)
    local_dx = write_void_cast_dx(orthogonal_map_cast(display), root / f"{prefix}_pymol_local.dx", smoothing=smoothing)
    # Native validation compares actual triangles with original point positions,
    # independently of the DX basis and viewer object transform.
    import numpy as np
    reference_path = root / f"{prefix}_geometry_reference.npz"
    reference_points = np.asarray([p.position for c in display.components for p in c.points])
    grid, reference_origin = cast_grid(display)
    mask = grid > .5
    interior = mask.copy()
    for a in range(3):
        interior &= np.roll(mask, 1, axis=a) & np.roll(mask, -1, axis=a)
    reference_basis = np.asarray(display.metadata.get('grid_basis', np.eye(3)))
    indices = np.rint(((reference_points-reference_origin) @ reference_basis.T)/display.spacing).astype(int)
    boundary_points = reference_points[(mask & ~interior)[tuple(indices.T)]]
    np.savez_compressed(reference_path, points=reference_points, boundary_points=boundary_points,
                        protein_coords=np.asarray([a.coord for a in frame.atoms]) if frame else np.empty((0,3)))
    geometry = scene_geometry(display, frame, profile)
    if segments:
        # Mouth rings come from the channel cast alone; the camera frames everything.
        geometry['mouth_rings'] = scene_geometry(channel_display, frame, profile)['mouth_rings']
    geometry['mouth_rings_drawn'] = draw_rings and bool(geometry['mouth_rings'])
    geometry['exit_paths_drawn'] = draw_lines and bool(geometry.get('exit_paths'))
    geometry['cast_color_rgb'] = list(color)
    geometry_path = root / f"{prefix}_scene.json"
    geometry_path.write_text(json.dumps(geometry, indent=2)+'\n')
    basis = cast.metadata.get('grid_basis', ((1,0,0),(0,1,0),(0,0,1)))
    matrix = [basis[j][i] if j < 3 else 0 for i in range(3) for j in range(4)] + [0,0,0,1]
    manifest = {"volume_dx": str(dx), "display_dx": str(shown_dx),
                "pymol_local_dx": str(local_dx), "scene_json": str(geometry_path),
                "geometry_reference_npz": str(reference_path)}
    if exit_casts:
        manifest["exit_casts_dx"] = str(manifest_exit_dx)
    from skimage.measure import marching_cubes
    mesh_arrays = {}
    def expected_surface(part, key):
        binary, origin = cast_grid(part)
        field, display_spacing, _ = smoothing.field(binary, part.spacing)
        verts, faces, normals, _ = marching_cubes(field, level=.5, spacing=(display_spacing,)*3)
        part_basis = np.asarray(part.metadata.get('grid_basis', np.eye(3)))
        mesh_arrays[key+'_vertices'] = verts + origin @ part_basis.T
        mesh_arrays[key+'_faces'] = faces
        mesh_arrays[key+'_normals'] = normals
        return verts @ part_basis + origin
    union_vertices = expected_surface(display, 'volume')
    expected_surfaces = [] if (regions or segment_mode) else [union_vertices]
    region_rows, region_pml, region_tcl = [], [], ['set crevice_region_molecules {}']
    if regions:
        from dataclasses import replace
        from collections import Counter
        original_points = Counter(p.position for c in display.components for p in c.points)
        region_points = Counter(p.position for points in regions.values() for p in points)
        if original_points != region_points:
            raise ValueError('Display regions must partition the complete displayed cast exactly once')
        palette = [color, (.42,.38,.72), (.20,.53,.86), (.67,.40,.62)]
        region_pml += ['cmd.disable("crevice_volume")']
        region_tcl += ['mol showrep $crevice_volume 0 off',
                       'if {[lsearch -exact [material list] CREVICEOtherVoid] < 0} {material add CREVICEOtherVoid}',
                       'material change opacity CREVICEOtherVoid 0.25']
        for index,(name,points) in enumerate(regions.items(), 1):
            color = ((.90,.58,.18) if name == 'shared' and len(regions)>1 else
                     (.55,.64,.69) if name == 'other_void' else palette[(index-1)%len(palette)])
            opacity = .25 if name == 'other_void' else 1.0
            part = replace(display, components=(replace(display.components[0], points=tuple(points), volume=len(points)*display.spacing**3),), points=tuple(points))
            mesh_key = f'region_{index:03d}'
            expected_surfaces.append(expected_surface(part, mesh_key))
            world_path = write_void_cast_dx(part, root/f'{prefix}_region_{index:03d}.dx', smoothing=smoothing)
            local_path = write_void_cast_dx(orthogonal_map_cast(part), root/f'{prefix}_region_{index:03d}_local.dx', smoothing=smoothing)
            obj, map_name = f'crevice_volume_region_{index:03d}', f'crevice_region_grid_{index:03d}'
            region_pml += [f'crevice_surface({obj!r}, {mesh_key!r})',
                f'cmd.set_color("crevice_region_color_{index}", {color!r})',
                f'cmd.color("crevice_region_color_{index}", {obj!r})',
                f'cmd.set("cgo_transparency", {1-opacity}, {obj!r})']
            color_id = 17+(index-1)%8
            region_tcl += [f'set crevice_region [mol new [file join $crevice_root {_tcl_word(world_path.name)}] type dx waitfor all]',
                'mol delrep 0 $crevice_region', 'mol representation Isosurface 0.5 0 0 0 1 1',
                f'mol color ColorID {color_id}', 'mol material '+('CREVICEOtherVoid' if opacity<1 else 'CREVICECast'), 'mol addrep $crevice_region',
                f'color change rgb {color_id} '+' '.join(str(x) for x in color),
                'lappend crevice_region_molecules $crevice_region']
            region_rows.append({'name':name,'pymol_object':obj,'pymol_mesh_key':mesh_key, 'color_rgb':color,'opacity':opacity,'points':len(points),
                                'volume_dx':str(world_path), 'pymol_local_dx':str(local_path)})
        region_manifest = root/f'{prefix}_display_regions.json'
        region_manifest.write_text(json.dumps(region_rows,indent=2)+'\n')
        manifest['display_regions_json'] = str(region_manifest)
    segment_rows = []
    from .cast_segments import SEGMENT_RGB
    if segment_mode:
        parts = [("lumen", lumen_display, None)] + [
            (c.metadata["segment"], _replace(display, components=(c,), points=tuple(c.points)), c) for c in segments]
        for name, part, component in parts:
            if not part.points:
                continue
            kind = name.split("_")[0]
            rgb = color if kind == "lumen" else tuple(SEGMENT_RGB[kind])
            expected_surfaces.append(expected_surface(part, name))
            world_path = write_void_cast_dx(part, root/f'{prefix}_{name}_display.dx', smoothing=smoothing)
            local_path = write_void_cast_dx(orthogonal_map_cast(part), root/f'{prefix}_{name}_display_local.dx',
                                            smoothing=smoothing)
            if component is None:
                measured_path, volume = dx, cast.total_volume
            else:
                measured_path = write_void_cast_dx(_replace(cast, components=(component,),
                                                            points=tuple(component.points)),
                                                   root/f'{prefix}_{name}_cast.dx')
                volume = component.volume
            segment_rows.append({'segment': name, 'segment_type': kind, 'pymol_object': f'crevice_{name}',
                                 'vmd_molecule': f'crevice_{name}', 'chimerax_model': f'crevice_{name}',
                                 'pymol_mesh_key': name, 'color_rgb': list(rgb),
                                 'color_hex': '#'+''.join(f'{round(v*255):02x}' for v in rgb),
                                 'display_points': len(part.points), 'measured_volume_A3': volume,
                                 'measured_dx': str(measured_path), 'display_dx': str(world_path),
                                 'display_local_dx': str(local_path)})
            manifest[f'segment_{name}_display_dx'] = str(world_path)
            manifest[f'segment_{name}_display_local_dx'] = str(local_path)
            if component is not None:
                manifest[f'segment_{name}_dx'] = str(measured_path)
    lining_groups = {}
    for row in (lining or ()):
        group = lining_groups.setdefault(row['segment'], {'residues': [], 'by_chain': {}})
        group['residues'].append(row['residue'])
        group['by_chain'].setdefault(row['viewer_chain'], []).append(int(row['viewer_resid']))
    lining_rows = []
    for name, group in lining_groups.items():
        kind = name.split('_')[0]
        rgb = color if kind == 'lumen' else tuple(SEGMENT_RGB[kind])
        lining_rows.append({'segment': name, 'residue_count': len(group['residues']),
                            'residues': group['residues'],
                            'pymol_object': f'crevice_{name}_lining', 'pymol_selection': f'crevice_{name}_lining_sel',
                            'chimerax_name': f'crevice_{name}_lining', 'vmd_representation': f'crevice_lining {name}',
                            'pymol_spec': ' or '.join(f'(chain {c} and resi {"+".join(map(str, r))})'
                                                      for c, r in group['by_chain'].items()),
                            'vmd_spec': ' or '.join(f'(chain {c} and resid {" ".join(map(str, r))})'
                                                    for c, r in group['by_chain'].items()),
                            'chimerax_spec': '#1'+''.join(f'/{c}:{",".join(map(str, r))}'
                                                         for c, r in group['by_chain'].items()),
                            'color_rgb': list(rgb), 'shown_by_default': False})
    if segment_mode or lining_rows:
        geometry['cast_segments'] = [{k: r[k] for k in ('segment', 'segment_type', 'pymol_object', 'color_rgb',
                                                     'color_hex', 'display_points', 'measured_volume_A3')}
                                         for r in segment_rows]
        geometry['lining_selections'] = [{k: r[k] for k in ('segment', 'residue_count', 'pymol_object',
                                                             'pymol_selection', 'chimerax_name', 'shown_by_default')}
                                         for r in lining_rows]
        geometry_path.write_text(json.dumps(geometry, indent=2)+'\n')
    np.savez_compressed(reference_path, points=reference_points, boundary_points=boundary_points,
                        surface_vertices=np.concatenate(expected_surfaces),
                        protein_coords=np.asarray([a.coord for a in frame.atoms]) if frame else np.empty((0,3)))
    mesh_path = root/f'{prefix}_pymol_mesh.npz'
    np.savez_compressed(mesh_path, **mesh_arrays)
    manifest['pymol_mesh_npz'] = str(mesh_path)
    pml = root / f"{prefix}_volume.pml"
    if segment_mode:
        cast_load_pml = ['cmd.delete("crevice_lumen*")', 'cmd.delete("crevice_entry_*")', 'cmd.delete("crevice_exit_*")']
        cast_color_pml = []
        for row in segment_rows:
            obj = row['pymol_object']
            cast_load_pml.append(f'crevice_surface({obj!r}, {row["pymol_mesh_key"]!r})')
            cast_color_pml += [f'cmd.set_color("{obj}_color", {row["color_rgb"]!r})', f'cmd.color("{obj}_color", {obj!r})',
                               f'cmd.set("cgo_transparency", 0.0, {obj!r})']
    else:
        cast_load_pml = ['crevice_surface("crevice_volume", "volume")']
        cast_color_pml = ['cmd.set_color("crevice_cast_color", '+repr(list(color))+')',
                          'cmd.color("crevice_cast_color", "crevice_volume")',
                          'cmd.set("cgo_transparency", 0.0, "crevice_volume")']
    lining_pml = []
    if lining_rows:
        lining_pml = ['# Lining residues: one stick object per segment, hidden by default (no labels).',
                      '# Toggle: crevice_lining lumen (or entry_1, exit_1, ...; on|off|toggle);',
                      '# highlight in place: color red, crevice_lumen_lining_sel; show sticks, crevice_lumen_lining_sel']
        for row in lining_rows:
            obj, sel = row['pymol_object'], row['pymol_selection']
            lining_pml += [f'cmd.select({sel!r}, "crevice_protein and ({row["pymol_spec"]})", enable=0)',
                           f'cmd.create({obj!r}, {sel!r}, zoom=0)',
                           f'cmd.hide("everything", {obj!r})', f'cmd.show("sticks", {obj!r} + " and not hydro")',
                           f'cmd.set_color("{obj}_color", {row["color_rgb"]!r})',
                           f'cmd.color("{obj}_color", {obj!r} + " and elem C")',
                           f'cmd.color("atomic", {obj!r} + " and not elem C")',
                           f'cmd.set("stick_radius", 0.2, {obj!r})', f'cmd.disable({obj!r})']
        lining_pml += ['def crevice_lining(segment="lumen", state="toggle"):',
                       '    name = "crevice_%s_lining" % segment',
                       '    shown = name in cmd.get_names("objects", enabled_only=1)',
                       '    (cmd.enable if state == "on" or (state == "toggle" and not shown) else cmd.disable)(name)',
                       'cmd.extend("crevice_lining", crevice_lining)']
    # Capped channel ends: lateral exit legs as thin blue centre-line tubes.
    exit_pml = ''
    if geometry['exit_paths_drawn']:
        exit_pml = '''cmd.delete("crevice_exits")
crevice_exit_cgo = []
for leg in crevice_scene.get("exit_paths", []):
    for a, b in zip(leg["points"], leg["points"][1:]):
        crevice_exit_cgo.extend([CYLINDER, *a, *b, 0.12, 0.16, 0.47, 0.84, 0.16, 0.47, 0.84])
cmd.load_cgo(crevice_exit_cgo, "crevice_exits")
'''
    pml.write_text(f'''# CREVICE display extensions are excluded from measured volume.
python
# Measured grid file: {dx.name!r}
from pathlib import Path
import inspect
import json
from pymol import cmd
from pymol.cgo import CYLINDER, BEGIN, END, TRIANGLES, NORMAL, VERTEX
import numpy as np
_crevice_script_root = Path(getattr(cmd._pymol, "__script__",
    inspect.currentframe().f_code.co_filename)).resolve().parent
def crevice_input(name, fallback):
    for directory in (_crevice_script_root, Path.cwd(), Path(fallback).parent):
        path = directory / name
        if path.is_file():
            return str(path)
    raise FileNotFoundError(name)
for name in ("crevice_protein", "crevice_volume", "crevice_mouths"):
    cmd.delete(name)
cmd.delete("crevice_volume_region_*")
cmd.delete("crevice_grid")
cmd.delete("crevice_region_grid*")
cmd.load(crevice_input({structure.name!r}, {str(structure)!r}), "crevice_protein", state=1)
# Standalone triangle surfaces survive map deletion and later scene rebuilds.
crevice_mesh = np.load(crevice_input({mesh_path.name!r}, {str(mesh_path)!r}), allow_pickle=False)
def crevice_surface(name, key):
    faces = crevice_mesh[key+"_faces"]
    body = np.empty((len(faces), 3, 8), dtype=float)
    body[:, :, 0] = NORMAL
    body[:, :, 1:4] = crevice_mesh[key+"_normals"][faces]
    body[:, :, 4] = VERTEX
    body[:, :, 5:8] = crevice_mesh[key+"_vertices"][faces]
    cmd.load_cgo([BEGIN, TRIANGLES, *body.ravel().tolist(), END], name, state=1, zoom=0)
    cmd.set_object_ttt(name, {matrix!r})
{chr(10).join(cast_load_pml)}
cmd.hide("everything", "crevice_protein")
cmd.dss("crevice_protein")
cmd.show("cartoon", "crevice_protein")
cmd.set_color("crevice_context", [0.56, 0.63, 0.70])
cmd.color("crevice_context", "crevice_protein")
cmd.set("cartoon_transparency", 0.25, "crevice_protein")
cmd.set("cartoon_sampling", 32, "crevice_protein")
cmd.set("cartoon_loop_quality", 24, "crevice_protein")
cmd.set("cartoon_oval_quality", 24, "crevice_protein")
cmd.set("cartoon_tube_quality", 24, "crevice_protein")
cmd.set("cartoon_fancy_helices", 1, "crevice_protein")
cmd.set("cartoon_fancy_sheets", 1, "crevice_protein")
cmd.set("cartoon_smooth_loops", 1, "crevice_protein")
cmd.set("cartoon_flat_sheets", 1, "crevice_protein")
cmd.set("cartoon_loop_radius", 0.18, "crevice_protein")
cmd.set("cartoon_oval_length", 1.25, "crevice_protein")
cmd.set("cartoon_oval_width", 0.24, "crevice_protein")
{chr(10).join(cast_color_pml)}
{chr(10).join(region_pml + lining_pml)}
crevice_mesh.close()
cmd.set("two_sided_lighting", 1)
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("orthoscopic", 1)
cmd.set("depth_cue", 0)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.35)
cmd.set("direct", 0.65)
cmd.set("light_count", 3)
cmd.set("transparency_mode", 2)
cmd.set("ray_transparency_oblique", 0.2)
cmd.set("ray_trace_mode", 0)
cmd.set("ambient_occlusion_mode", 2)
cmd.set("ambient_occlusion_scale", 12)
cmd.set("ambient_occlusion_smooth", 10)
cmd.set("specular", 0.08)
cmd.set("shininess", 18)
cmd.set("antialias", 2)
with open(crevice_input({geometry_path.name!r}, {str(geometry_path)!r})) as handle:
    crevice_scene = json.load(handle)
crevice_cgo = []
for ring in (crevice_scene["mouth_rings"] if crevice_scene.get("mouth_rings_drawn") else []):
    for a, b in zip(ring["points"], ring["points"][1:]):
        crevice_cgo.extend([CYLINDER, *a, *b, 0.065, 0.90, 0.52, 0.12, 0.90, 0.52, 0.12])
if crevice_cgo:
    cmd.load_cgo(crevice_cgo, "crevice_mouths")
{exit_pml}cmd.set("field_of_view", 25)
crevice_distance = crevice_scene["vertical_span"] / (2 * 0.22169466264293988)
# PyMOL get_view/set_view uses a column-major model-to-camera matrix.
crevice_rotation = [crevice_scene["camera_rotation"][i][j] for j in range(3) for i in range(3)]
cmd.set_view([*crevice_rotation, 0, 0, -crevice_distance, *crevice_scene["center"],
              crevice_distance-crevice_scene["radius"]*1.5,
              crevice_distance+crevice_scene["radius"]*1.5, 1])
cmd.viewport(1200, 900)
def crevice_render(path=None):
    target = path or str(_crevice_script_root / {str(prefix + '_pymol.png')!r})
    cmd.png(target, width=3200, height=2400, dpi=400, ray=1)
cmd.extend("crevice_render", crevice_render)
# Type crevice_render for the final 3200 x 2400 image.
python end
''')
    manifest["volume_pml"] = str(pml)
    rotation = geometry['camera_rotation']
    center = geometry['center']
    rot_tcl = '{' + ' '.join('{'+' '.join(f'{v:.12g}' for v in row)+' 0}' for row in rotation) + ' {0 0 0 1}}'
    center_tcl = '{'+' '.join('{'+ ' '.join(str(int(i==j)) for j in range(3)) +f' {-center[i]:.12g}'+'}' for i in range(3))+' {0 0 0 1}}'
    scale_factor = 2.0 / geometry['vertical_span']
    scale_tcl = '{'+' '.join('{'+ ' '.join(f'{scale_factor if i==j else 0:.12g}' for j in range(3))+' 0}' for i in range(3))+' {0 0 0 1}}'
    ring_tcl = []
    for ring in (geometry['mouth_rings'] if draw_rings else []):
        ring_tcl.append('graphics $crevice_protein color 3')
        for a, b in zip(ring['points'], ring['points'][1:]):
            aa, bb = ('{'+' '.join(f'{v:.9g}' for v in p)+'}' for p in (a,b))
            ring_tcl.append(f'graphics $crevice_protein cylinder {aa} {bb} radius 0.065 resolution 16 filled yes')
    if geometry['exit_paths_drawn']:
        ring_tcl.append('color change rgb 23 0.16 0.47 0.84')
        ring_tcl.append('graphics $crevice_protein color 23')
        for leg in geometry['exit_paths']:
            for a, b in zip(leg['points'], leg['points'][1:]):
                aa, bb = ('{'+' '.join(f'{v:.9g}' for v in p)+'}' for p in (a,b))
                ring_tcl.append(f'graphics $crevice_protein cylinder {aa} {bb} radius 0.12 resolution 16 filled yes')
    tcl = root / f"{prefix}_volume.tcl"
    vmd_main_dx = Path(segment_rows[0]['display_dx']) if segment_rows and segment_rows[0]['segment'] == 'lumen' else shown_dx
    segment_tcl, segment_list = [], ''
    if segment_mode:
        # Colour ids: 10 lumen (cast colour), 26 entry, 27 exit.
        segment_tcl = ['set crevice_segment_molecules {}', 'array set crevice_lining_rep {}',
                       'color change rgb 26 '+' '.join(f'{v:.6g}' for v in SEGMENT_RGB['entry']),
                       'color change rgb 27 '+' '.join(f'{v:.6g}' for v in SEGMENT_RGB['exit'])]
        for row in segment_rows:
            if row['segment'] == 'lumen':
                segment_tcl.append('mol rename $crevice_volume crevice_lumen')
                continue
            color_id = 26 if row['segment_type'] == 'entry' else 27
            segment_tcl += [f'set crevice_segment [mol new [file join $crevice_root {_tcl_word(Path(row["display_dx"]).name)}] type dx waitfor all]',
                            'mol delrep 0 $crevice_segment', 'mol representation Isosurface 0.5 0 0 0 1 1',
                            f'mol color ColorID {color_id}', 'mol material CREVICECast', 'mol addrep $crevice_segment',
                            f'mol rename $crevice_segment {row["vmd_molecule"]}',
                            'lappend crevice_segment_molecules $crevice_segment']
        for row in lining_rows:
            kind = row['segment'].split('_')[0]
            color_id = {'lumen': 10, 'entry': 26, 'exit': 27}[kind]
            segment_tcl += ['mol representation Licorice 0.2 12 12', f'mol selection {{{row["vmd_spec"]}}}',
                            f'mol color ColorID {color_id}', 'mol material Opaque', 'mol addrep $crevice_protein',
                            f'set crevice_lining_rep({row["segment"]}) [expr {{[molinfo $crevice_protein get numreps]-1}}]',
                            f'mol showrep $crevice_protein $crevice_lining_rep({row["segment"]}) off']
        if lining_rows:
            segment_tcl += ['# Lining residues (hidden by default): crevice_lining lumen on|off (or entry_1, exit_1, ...)',
                            'proc crevice_lining {{segment lumen} {state on}} {',
                            '    global crevice_protein crevice_lining_rep',
                            '    mol showrep $crevice_protein $crevice_lining_rep($segment) [expr {$state eq "on"}]',
                            '}']
        segment_list = ' $crevice_segment_molecules'
    tcl.write_text(f'''# Measured grid: {dx.name}; display-only continuation: {shown_dx.name}
set crevice_root [file dirname [file normalize [info script]]]
set crevice_structure [file join $crevice_root {_tcl_word(structure.name)}]
if {{![file exists $crevice_structure]}} {{set crevice_structure {_tcl_word(str(structure))}}}
set crevice_protein [mol new $crevice_structure waitfor all]
mol delrep 0 $crevice_protein
if {{[lsearch -exact [material list] CREVICEContext] < 0}} {{material add CREVICEContext}}
material change opacity CREVICEContext 0.75
material change ambient CREVICEContext 0.30
material change diffuse CREVICEContext 0.65
material change specular CREVICEContext 0.08
mol representation NewCartoon 0.25 40 4.0 1
mol selection all
mol color ColorID 2
mol material CREVICEContext
mol addrep $crevice_protein
set crevice_volume [mol new [file join $crevice_root {_tcl_word(vmd_main_dx.name)}] type dx waitfor all]
mol delrep 0 $crevice_volume
mol selection all
mol representation Isosurface 0.5 0 0 0 1 1
mol color ColorID 10
if {{[lsearch -exact [material list] CREVICECast] < 0}} {{material add CREVICECast}}
material change ambient CREVICECast 0.30
material change diffuse CREVICECast 0.70
material change specular CREVICECast 0.08
material change shininess CREVICECast 0.12
material change opacity CREVICECast 1.0
mol material CREVICECast
mol addrep $crevice_volume
color change rgb 10 {' '.join(f'{v:.6g}' for v in color)}
color change rgb 2 0.56 0.63 0.70
color change rgb 3 0.90 0.52 0.12
{chr(10).join(ring_tcl)}
{chr(10).join(region_tcl + segment_tcl)}
color Display Background white
display projection Orthographic
display depthcue off
display antialias on
display shadows off
display ambientocclusion on
display aoambient 0.8
display aodirect 0.4
axes location Off
display resetview
foreach crevice_molecule [concat [list $crevice_protein $crevice_volume] $crevice_region_molecules{segment_list}] {{
    molinfo $crevice_molecule set center_matrix [list {center_tcl}]
    molinfo $crevice_molecule set rotate_matrix [list {rot_tcl}]
    molinfo $crevice_molecule set scale_matrix [list {scale_tcl}]
}}
display height 4.0
display resize 1200 900
display update
proc crevice_render {{}} {{
    global crevice_root
    display resize 2400 1800
    display redraw
    display update
    render TachyonInternal [file join $crevice_root {_tcl_word(prefix + '_vmd.tga')}]
}}
# Type crevice_render for the final 2400 x 1800 image.
''')
    manifest["volume_tcl"] = str(tcl)
    vmd=root/f"{prefix}_volume.vmd"
    vmd.write_text(tcl.read_text())
    manifest["volume_vmd"]=str(vmd)
    bild = root / f"{prefix}_mouths.bild"
    bild_lines = ['.color 0.90 0.52 0.12']
    for ring in (geometry['mouth_rings'] if draw_rings else []):
        for a, b in zip(ring['points'], ring['points'][1:]):
            bild_lines.append('.cylinder '+' '.join(f'{v:.9g}' for v in (*a,*b))+ ' 0.065')
    if geometry['exit_paths_drawn']:
        bild_lines.append('.color 0.16 0.47 0.84')
        for leg in geometry['exit_paths']:
            for a, b in zip(leg['points'], leg['points'][1:]):
                bild_lines.append('.cylinder '+' '.join(f'{v:.9g}' for v in (*a,*b))+ ' 0.12')
    bild.write_text('\n'.join(bild_lines)+'\n')
    manifest['mouth_guides_bild'] = str(bild)
    camera_matrix = [rotation[j][i] if j<3 else center[i]+rotation[2][i]*geometry['radius']*3
                     for i in range(3) for j in range(4)]
    cxc = root / f"{prefix}_volume.cxc"
    # ChimeraX changes to the command script's directory while opening .cxc.
    # Its APBS/OpenDX reader ignores off-diagonal delta directions. Load the
    # orthogonal local map and explicitly place each model in input coordinates.
    commands = [f"open {json.dumps(structure.name if structure.parent == root else str(structure))} id #1",
                f"open {json.dumps(Path(segment_rows[0]['display_local_dx']).name if segment_rows and segment_rows[0]['segment'] == 'lumen' else local_dx.name)} id #2",
                "view matrix models #2,"+",".join(f"{v:.12g}" for v in matrix[:12]),
                "hide #1 atoms", "cartoon #1", "color #1 #8fa1b3",
                "transparency #1 25 target c",
                "volume #2 style surface level 0.5 step 1 surfaceSmoothing false",
                f"color #2 {color_hex}", "transparency #2 0",
                "set bgColor white", "camera ortho", "lighting soft", "material dull", "graphics quality 3",
                "graphics silhouettes false"]
    guides_drawn = geometry['mouth_rings_drawn'] or geometry['exit_paths_drawn']
    if guides_drawn:
        commands.append(f"open {json.dumps(bild.name)} id #3")
    if region_rows:
        commands.append('hide #2 models')
        for index,row in enumerate(region_rows, 4 if guides_drawn else 3):
            color = '#'+''.join(f'{round(v*255):02x}' for v in row['color_rgb'])
            commands += [f'open {json.dumps(Path(row["pymol_local_dx"]).name)} id #{index}',
                         f'view matrix models #{index},'+",".join(f"{v:.12g}" for v in matrix[:12]),
                         f'volume #{index} style surface level 0.5 step 1 surfaceSmoothing false',
                         f'color #{index} {color}', f'transparency #{index} {100*(1-row["opacity"]):g}']
    for index, row in enumerate(segment_rows, 10):
        if row['segment'] == 'lumen':
            commands.append('rename #2 crevice_lumen')
            continue
        commands += [f'open {json.dumps(Path(row["display_local_dx"]).name)} id #{index}',
                     f'view matrix models #{index},'+",".join(f"{v:.12g}" for v in matrix[:12]),
                     f'volume #{index} style surface level 0.5 step 1 surfaceSmoothing false',
                     f'color #{index} {row["color_hex"]}', f'transparency #{index} 0',
                     f'rename #{index} {row["chimerax_model"]}']
    for row in lining_rows:
        commands.append(f'name frozen {row["chimerax_name"]} {row["chimerax_spec"]}')
    if lining_rows:
        commands.append('# Lining residues (hidden by default): show crevice_lumen_lining atoms; '
                        'style crevice_lumen_lining stick; hide crevice_lumen_lining atoms')
    commands += ["windowsize 900 675", "view matrix camera "+','.join(f'{v:.12g}' for v in camera_matrix), "view", "zoom 0.9"]
    cxc.write_text("\n".join(commands) + "\n")
    manifest["volume_cxc"] = str(cxc)
    render_cxc = root / f'{prefix}_render.cxc'
    render_cxc.write_text(f'open {json.dumps(cxc.name)}\nsave {json.dumps(prefix+"_chimerax.png")} width 2400 height 1800 supersample 3\n')
    manifest['chimerax_render_cxc'] = str(render_cxc)
    display_smoothing = smoothing.metadata(display.spacing, cast_grid(display)[0].shape)
    info = root / f"{prefix}_volume_metadata.json"
    info.write_text(json.dumps({"spacing_A": cast.spacing, "isovalue": 0.5,
                  "grid_value": "selected_void_occupancy", "smoothing": smoothing.width_A > 0,
                  "surface_smoothing_A": smoothing.width_A if smoothing.option == 'surface_smoothing' else None,
                  "smooth_A": smoothing.width_A if smoothing.option != 'surface_smoothing' else None,
                  "smoothing_definition": display_smoothing['definition'],
                  "display_smoothing": display_smoothing,
                  "unsmoothed_reference_dx": str(dx),
                  "material": "matte",
                  "protein_representation": "cartoon_only",
                  "cartoon_transparency": .25,
                  "pymol_render_size": [3200, 2400],
                  "pymol_representation": "standalone_triangle_mesh",
                  "pymol_helper_grids": "not loaded; stale helper objects deleted",
                  "pymol_surface_definition": "triangles and interpolated normals from the same 0.5 display field, independent of a live map",
                  "measurement_points": sum(c.point_count for c in cast.components),
                  "measured_volume_A3": cast.total_volume,
                  "subsampled_display_points": cast.point_count,
                  "display_points": sum(c.point_count for c in display.components),
                  "display_extension": {k:v for k,v in display.metadata.items() if k.startswith('extension_')},
                  "display_volume_is_measurement": False,
                  "pymol_local_to_input_matrix": matrix,
                  "channel_mouths": cast.metadata.get('channel_mouths'),
                  "cast_color_rgb": list(color),
                  "mouth_guides_drawn": geometry['mouth_rings_drawn'],
                  "exit_centre_lines_drawn": geometry['exit_paths_drawn'],
                  **({"lateral_exit_casts": exit_cast_summary(exit_casts),
                      "lateral_exit_casts_display": "separate EXIT segment objects (crevice_exit_N); volumes separate from measured_volume_A3"}
                     if exit_casts else {}),
                  **({"cast_segments": segment_rows,
                      "cast_segments_display": "one object per segment (crevice_entry_N, crevice_lumen, crevice_exit_N); "
                                               "measured_volume_A3 is the LUMEN only"}
                     if segment_mode else {}),
                  **({"lining_selections": lining_rows} if lining_rows else {}),
                  "representation": "voxel boundary approximation, not electron density",
                  "boundary_validation": "not_continuously_validated",
                  "connectivity": cast.metadata.get("connectivity"),
                  "component_ids": [c.id for c in cast.components],
                  "display_regions": region_rows,
                  "files": manifest}, indent=2) + "\n")
    manifest["volume_metadata_json"] = str(info)
    return manifest

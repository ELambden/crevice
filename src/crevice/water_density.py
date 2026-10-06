"""Aligned water-oxygen number density with auditable counts and periodic images.

Used by ``crevice hydration`` and the hydration pass of ``publish``,
``cavity-trajectory`` and ``region-trajectory``. :class:`WaterDensity`
accumulates aligned water-oxygen positions frame by frame on a fixed grid
around a reference neighbourhood (every frame counts in the denominator, so
absent water is a zero, not a gap); :func:`write_density_bundle` writes the
measured and display maps, projections and metadata; :func:`write_scalar_dx`
writes any scalar grid as OpenDX. Densities are observations in a fixed
region, not hydration free energies or accessibility tests.
"""
from itertools import product
from pathlib import Path
import json
import math

from .presentation import FigureText, annotate_option


def write_scalar_dx(values, origin, spacing, path, *, basis=None):
    """Write a 3D scalar field as OpenDX (samples at voxel centres; last index fastest).

    Parameters
    ----------
    values : array_like, shape (nx, ny, nz)
        Finite scalar values.
    origin : sequence of float
        Position of sample (0, 0, 0), Å (a voxel centre).
    spacing : float
        Grid spacing, Å.
    path : str or Path
        Output ``.dx`` path.
    basis : array_like, shape (3, 3), optional
        Grid axes as rows (unit vectors); default the Cartesian axes.

    Returns
    -------
    str
        ``path`` as a string.

    Raises
    ------
    ValueError
        If ``values`` is not a finite 3D array.
    """
    import numpy as np
    values = np.asarray(values, dtype=float)
    if values.ndim != 3 or not np.isfinite(values).all():
        raise ValueError('DX requires a finite three-dimensional scalar field')
    nx, ny, nz = values.shape
    lines = [f'object 1 class gridpositions counts {nx} {ny} {nz}',
             'origin ' + ' '.join(f'{v:.12g}' for v in origin)]
    lines += ['delta ' + ' '.join(f'{v:.12g}' for v in row) for row in (np.eye(3) if basis is None else np.asarray(basis))*spacing]
    lines += [f'object 2 class gridconnections counts {nx} {ny} {nz}',
              f'object 3 class array type double rank 0 items {values.size} data follows']
    flat = values.ravel(order='C')
    lines += [' '.join(f'{v:.12g}' for v in flat[i:i+3]) for i in range(0, len(flat), 3)]
    lines += ['attribute "dep" string "positions"', 'object "water density" class field',
              'component "positions" value 1', 'component "connections" value 2', 'component "data" value 3']
    Path(path).write_text('\n'.join(lines)+'\n')
    return str(path)


class WaterDensity:
    """Fixed reference neighbourhood; every frame contributes to the denominator.

    The grid covers the reference points padded by ``margin + 4*smoothing +
    2*spacing``; only oxygens within ``margin`` of a reference point are binned
    (``definition`` names this domain). Call :meth:`add` once per frame, exactly
    ``expected_frames`` times, then :meth:`result`.

    Parameters
    ----------
    reference_points : array_like, shape (n, 3)
        Points defining the neighbourhood (for example cavity samples or focus
        residue atoms), aligned frame, Å.
    margin : float, default 2.0
        Neighbourhood margin, Å.
    spacing : float, default 0.5
        Voxel edge, Å.
    smoothing : float, default 0.5
        Gaussian standard deviation (Å) of the display map only; 0 = none.
    expected_frames : int, default 1
        Number of frames that will be added.
    max_points : int, default 2000000
        Largest allowed grid.
    definition : str, default "reference_neighbourhood"
        Name of the domain, recorded in the metadata.

    Raises
    ------
    ValueError
        For non-finite or empty reference points, a nonpositive margin or
        spacing, negative smoothing, a nonpositive frame count or a grid above
        ``max_points``.
    """

    def __init__(self, reference_points, *, margin=2., spacing=.5, smoothing=.5,
                 expected_frames=1, max_points=2_000_000, definition='reference_neighbourhood'):
        import numpy as np
        from scipy.spatial import cKDTree
        points = np.asarray(reference_points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not len(points) or not np.isfinite(points).all():
            raise ValueError('Density requires finite reference points')
        if not all(math.isfinite(v) and v > 0 for v in [margin, spacing]):
            raise ValueError('Density margin and spacing must be positive')
        if not math.isfinite(smoothing) or smoothing < 0:
            raise ValueError('Density smoothing must be nonnegative')
        if int(expected_frames) != expected_frames or expected_frames < 1:
            raise ValueError('Density requires a positive frame count')
        self.spacing, self.margin, self.smoothing = float(spacing), float(margin), float(smoothing)
        self.expected_frames = int(expected_frames)
        self.definition = definition
        self.tree = cKDTree(points)
        self.anchor = points.mean(0)
        padding = margin + 4*smoothing + 2*spacing
        self.edge = np.floor((points.min(0)-padding)/spacing)*spacing
        self.shape = np.ceil((points.max(0)+padding-self.edge)/spacing).astype(int)
        if int(np.prod(self.shape)) > max_points:
            raise ValueError('Water density grid exceeds max_points; choose a coarser explicit spacing or smaller residue focus')
        self.origin = self.edge + spacing/2
        self.upper = self.edge + self.shape*spacing
        self.counts = np.zeros(tuple(self.shape), dtype=np.int64)
        self.occupied_frames = np.zeros_like(self.counts)
        self.first_half_counts = np.zeros_like(self.counts)
        self.frame_counts = []
        self.image_counts = []

    def aligned_images(self, water_xyz, rotation=None, translation=None, box=None):
        """Enumerate exactly the lattice images that can intersect the map box.

        Parameters
        ----------
        water_xyz : array_like, shape (n, 3)
            Oxygen coordinates in the original frame, Å.
        rotation : array_like, shape (3, 3), optional
            Proper rotation of the frame's fit (``aligned = xyz @ rotation +
            translation``); default identity.
        translation : array_like, shape (3,), optional
            Translation of the fit, Å; default zero.
        box : array_like, optional
            Periodic cell ``[a, b, c, alpha, beta, gamma]``; ``None`` yields the
            aligned coordinates once.

        Yields
        ------
        numpy.ndarray, shape (n, 3)
            Aligned coordinates of one periodic image of every oxygen (central image
            about the anchor plus every lattice offset whose image can reach the
            grid).

        Raises
        ------
        ValueError
            For non-finite coordinates, an improper rotation or a domain spanning
            more than 125 periodic cells.
        """
        import numpy as np
        xyz = np.asarray(water_xyz, dtype=float).reshape((-1, 3))
        if not np.isfinite(xyz).all():
            raise ValueError('Nonfinite water coordinates')
        rotation = np.eye(3) if rotation is None else np.asarray(rotation, dtype=float)
        translation = np.zeros(3) if translation is None else np.asarray(translation, dtype=float)
        if not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-7) or not np.isclose(np.linalg.det(rotation), 1):
            raise ValueError('Density alignment requires a proper rigid rotation')
        if box is None:
            yield xyz @ rotation + translation
            return
        from .hydration import checked_box
        from MDAnalysis.lib.mdamath import triclinic_vectors
        checked_box(box, self.spacing)
        cell = triclinic_vectors(box).astype(float)
        inverse_cell = np.linalg.inv(cell)
        anchor_original = (self.anchor-translation) @ rotation.T
        fractional = (xyz-anchor_original) @ inverse_cell
        central = (fractional-np.floor(fractional+.5)) @ cell
        corners = np.asarray(list(product(*zip(self.edge, self.upper))))
        fractions = ((corners-self.anchor) @ rotation.T) @ inverse_cell
        low, high = np.floor(fractions.min(0)+.5).astype(int), np.floor(fractions.max(0)+.5).astype(int)
        if int(np.prod(high-low+1)) > 125:
            raise ValueError('Density domain spans too many periodic cells')
        for offset in product(*(range(a, b+1) for a, b in zip(low, high))):
            yield (central + np.asarray(offset) @ cell) @ rotation + self.anchor

    def add(self, waters, *, rotation=None, translation=None, box=None):
        """Bin one frame's aligned oxygens (all periodic images) into the grid.

        Each image inside the grid and within ``margin`` of a reference point adds
        one count to its voxel; the voxel's occupied-frame count rises by one. The
        first half of the declared frames is also accumulated separately for a
        split-half comparison.

        Parameters
        ----------
        waters : array_like, shape (n, 3)
            Oxygen coordinates in the original frame, Å.
        rotation, translation, box
            As in :meth:`aligned_images`.

        Raises
        ------
        ValueError
            If more frames are added than declared.
        """
        import numpy as np
        if len(self.frame_counts) >= self.expected_frames:
            raise ValueError('More density frames than declared')
        chunks = []
        images = 0
        for xyz in self.aligned_images(waters, rotation, translation, box):
            images += 1
            inside = np.all((xyz >= self.edge) & (xyz < self.upper), axis=1)
            xyz = xyz[inside]
            if not len(xyz):
                continue
            xyz = xyz[self.tree.query(xyz)[0] <= self.margin+1e-9]
            if len(xyz):
                scaled = (xyz-self.edge)/self.spacing
                # Lattice round trips can displace an exact bin boundary by a
                # few floating-point ulps. Snap only at numerical precision.
                nearest = np.rint(scaled)
                tolerance = 64*np.finfo(float).eps*np.maximum(1, np.abs(scaled))
                scaled = np.where(np.abs(scaled-nearest) <= tolerance, nearest, scaled)
                indices = np.floor(scaled).astype(int)
                chunks.append(np.ravel_multi_index(indices.T, tuple(self.shape)))
        flat = np.concatenate(chunks) if chunks else np.empty(0, dtype=int)
        unique, counts = np.unique(flat, return_counts=True)
        self.counts.ravel()[unique] += counts
        self.occupied_frames.ravel()[unique] += 1
        if len(self.frame_counts) < self.expected_frames//2:
            self.first_half_counts.ravel()[unique] += counts
        self.frame_counts.append(len(flat))
        self.image_counts.append(images)

    def result(self, water_population):
        """Densities and audit counts after every declared frame was added.

        Parameters
        ----------
        water_population : int
            Number of water oxygens in the system; 0 marks the density as
            unavailable (NaN), never zero.

        Returns
        -------
        dict
            ``counts``, ``occupied_frames``, ``first_half_counts`` (int grids),
            ``density_A3`` (counts / (frames x voxel volume), oxygens per Å³),
            ``voxel_occupancy`` (fraction of frames with an oxygen in the voxel),
            ``display_density_A3`` (Gaussian-smoothed copy for display), ``origin``
            (first voxel centre, Å), ``spacing_A``, ``frame_counts`` and
            ``metadata`` (status, normalisation, integrals, periodic-image counts,
            interpretation).

        Raises
        ------
        ValueError
            If fewer frames were added than declared.
        """
        import numpy as np
        from scipy.ndimage import gaussian_filter
        n = len(self.frame_counts)
        if n != self.expected_frames:
            raise ValueError('Density frame coverage is incomplete')
        observed = bool(water_population)
        density = self.counts.astype(float)/(n*self.spacing**3) if observed else np.full(self.counts.shape, np.nan)
        occupancy = self.occupied_frames.astype(float)/n if observed else np.full(self.counts.shape, np.nan)
        shown = gaussian_filter(density, self.smoothing/self.spacing, mode='constant', truncate=4) if self.smoothing and observed else density.copy()
        return {'counts': self.counts, 'occupied_frames': self.occupied_frames,
                'first_half_counts': self.first_half_counts, 'density_A3': density,
                'voxel_occupancy': occupancy, 'display_density_A3': shown,
                'origin': self.origin, 'spacing_A': self.spacing,
                'frame_counts': np.asarray(self.frame_counts),
                'metadata': {'status': 'observed' if observed else 'unavailable_no_explicit_water',
                    'frame_count': n, 'water_population': int(water_population), 'grid_shape': self.shape.tolist(),
                    'origin_A': self.origin.tolist(), 'bin_lower_edge_A': self.edge.tolist(),
                    'spacing_A': self.spacing, 'voxel_volume_A3': self.spacing**3,
                    'units': 'water oxygens per Å³', 'domain': self.definition,
                    'region_margin_A': self.margin, 'normalization': 'sum oxygen counts / (all observed frames * voxel volume)',
                    'voxel_occupancy_definition': 'fraction of all frames with one or more water oxygens in the voxel; depends on voxel size',
                    'display_gaussian_sigma_A': self.smoothing,
                    'mean_domain_water_count': float(np.mean(self.frame_counts)) if observed else None,
                    'density_integral_water_count': float(np.sum(density)*self.spacing**3) if observed else None,
                    'display_integral_water_count': float(np.sum(shown)*self.spacing**3) if observed else None,
                    'maximum_periodic_images_enumerated': max(self.image_counts),
                    'periodicity': 'lattice images enumerated around the aligned reference; a domain spanning cells can contain multiple images of one water',
                    'interpretation': 'fixed-region oxygen observations, not a hydration free energy or ligand-accessibility test; single-frame fields are snapshot observations'}}


@annotate_option
def write_density_bundle(result, output_dir, prefix, *, level=.05, dpi=240, annotate=None):
    """Write aligned water-oxygen density maps, projections and metadata.

    ``result`` (from :meth:`WaterDensity.result`) bins aligned water-oxygen
    positions from every analysed frame into cubic voxels. ``density_A3`` is
    the mean number density (oxygens per Å³), ``voxel_occupancy`` the fraction
    of frames with at least one oxygen in the voxel, and
    ``display_density_A3`` a display-only Gaussian-smoothed copy used for the
    viewer surface at ``level``. A three-panel figure projects the unsmoothed
    density along z, y and x. Nothing is written except the JSON when no water
    was observed.

    Parameters
    ----------
    result : dict
        Density result with ``density_A3``, ``voxel_occupancy``,
        ``display_density_A3``, ``origin``, ``spacing_A`` and ``metadata``.
    output_dir : str or Path
        Output directory.
    prefix : str
        File-name prefix.
    level : float, default 0.05
        Display isovalue (water oxygens per Å³) of the viewer surface; must be positive.
    dpi : int, default 240
        Figure resolution.
    annotate : bool, optional
        Draw the figure title. ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    dict of str to str
        Output keys and paths.

    Outputs
    -------
    PREFIX_water_density.dx, PREFIX_water_occupancy.dx, PREFIX_water_density_display.dx : DX
        Mean number density (Å⁻³), voxel occupancy fraction, display-smoothed density.
    PREFIX_water_density.npz : arrays
        The same fields with grid origin and spacing.
    PREFIX_water_density.png : figure
        Three projections (X-Y, X-Z, Y-Z, axes in Å) of the density summed along
        the third axis ("Projected oxygen density (Å⁻²)").
    PREFIX_water_density.json : JSON
        Metadata, including the display level and whether a surface exists at it.

    Colours and representation
    --------------------------
    Each projection uses the ``magma`` colour scale (black = none, pale yellow =
    highest) with its own colour bar "Projected oxygen density (Å⁻²)". In the
    viewer scenes the display density is a teal translucent isosurface (see
    :func:`crevice.analysis_viewers.write_analysis_assets`).
    """
    import numpy as np
    from .figures import _pyplot, _save_figure
    if not math.isfinite(level) or level <= 0:
        raise ValueError('Water density display level must be positive')
    root = Path(output_dir)
    files = {}
    def path(key, suffix):
        p = root/(prefix+suffix)
        files[key] = str(p)
        return p
    metadata = dict(result['metadata'], display_level_A3=float(level))
    metadata['display_surface_available'] = bool(np.nanmax(result['display_density_A3']) > level) if metadata['status']=='observed' else False
    if metadata['status']=='observed':
        for key, field, ending in [('water_density_dx', 'density_A3', '_water_density.dx'),
                                   ('water_occupancy_dx', 'voxel_occupancy', '_water_occupancy.dx'),
                                   ('water_density_display_dx', 'display_density_A3', '_water_density_display.dx')]:
            write_scalar_dx(result[field], result['origin'], result['spacing_A'], path(key, ending))
        np.savez_compressed(path('water_density_npz','_water_density.npz'), **{k: v for k,v in result.items() if k!='metadata'})
        plt = _pyplot()
        text = FigureText(('Mean aligned water-oxygen density over all analysed frames' if metadata['frame_count']>1 else
                           'Observed water-oxygen positions in one structure')+', projected along z, y and x')
        fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
        for ax, axis, names in zip(axes, [2, 1, 0], [('X','Y'),('X','Z'),('Y','Z')]):
            field = result['density_A3'].sum(axis=axis)*result['spacing_A']
            remaining = [i for i in range(3) if i!=axis]
            low = result['origin'][remaining]-result['spacing_A']/2
            high = low+np.asarray(field.shape)*result['spacing_A']
            im = ax.imshow(field.T, origin='lower', extent=[low[0],high[0],low[1],high[1]], cmap='magma')
            ax.set_xlabel(names[0]+' (Å)'); ax.set_ylabel(names[1]+' (Å)')
            fig.colorbar(im, ax=ax, label='Projected oxygen density (Å⁻²)', shrink=.8)
        text.suptitle(fig, 'Water oxygen density · all-frame mean' if metadata['frame_count']>1 else 'Observed water positions · static snapshot')
        _save_figure(fig, path('water_density_png','_water_density.png'), dpi=dpi, text=text)
    path('water_density_json','_water_density.json').write_text(json.dumps(metadata,indent=2,allow_nan=False,ensure_ascii=False)+'\n',encoding='utf-8')
    return files

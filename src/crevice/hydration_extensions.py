"""Shared spatial/chemical extensions for standard hydration workflows.

Both the static (single-structure) and trajectory hydration workflows can
add two extra analyses to their results:

* a water-oxygen number density on a grid around a fixed reference region
  (:class:`crevice.water_density.WaterDensity`), and
* chemically typed residue interactions with explicit hydrogens
  (:class:`crevice.interaction_analysis.InteractionAnalysis`).

This module holds their default settings and attaches them to a static
result.
"""
from pathlib import Path


def defaults(options=None):
    """Extension settings with defaults filled in.

    Parameters
    ----------
    options : dict, optional
        Overrides; unknown keys are kept unchanged.

    Returns
    -------
    dict
        ``density_spacing`` (0.5 Å grid spacing), ``density_smoothing`` (0.5 Å
        display-only Gaussian sigma; 0 keeps the raw density),
        ``density_level`` (0.05 water oxygens per Å³, the display isovalue),
        ``density_max_points`` (8,000,000 grid points), and the flags
        ``skip_density``, ``skip_interactions`` and ``skip_views`` (all
        ``False``).

    Examples
    --------
    >>> from crevice.hydration_extensions import defaults
    >>> settings = defaults({"skip_views": True})
    >>> settings["density_spacing"], settings["skip_views"]
    (0.5, True)
    """
    result = {'density_spacing':.5,'density_smoothing':.5,'density_level':.05,
              'density_max_points':8_000_000,'skip_density':False,'skip_interactions':False,
              'skip_views':False}
    result.update(options or {})
    return result


def reference_points(frame, focus, volume_path=None):
    """Points that define the fixed density neighbourhood.

    Parameters
    ----------
    frame : StructureFrame
    focus : iterable of str
        Residue labels; their atoms are used when no volume is given.
    volume_path : str or pathlib.Path, optional
        Binary OpenDX map of a measured cavity; its occupied voxel positions
        are used instead.

    Returns
    -------
    points : ndarray, shape (n, 3)
        Reference positions in Å.
    definition : str
        ``"fixed_measured_cavity_neighbourhood"`` or
        ``"fixed_reference_residue_neighbourhood"``.
    """
    import numpy as np
    if volume_path:
        from .residue_evidence import read_binary_dx
        grid,origin,deltas = read_binary_dx(volume_path)
        return np.argwhere(grid) @ deltas + origin, 'fixed_measured_cavity_neighbourhood'
    return np.asarray([a.coord for a in frame.atoms if a.residue_key.label in set(focus)]), 'fixed_reference_residue_neighbourhood'


def attach_static_extensions(analysis, frame, *, box=None, volume_path=None, options=None):
    """Add the water density and interaction analyses to a static hydration result.

    The protein reference used for the density is the frame's non-HETATM,
    non-hydrogen atoms. The density margin is 2 Å around a cavity map, or the
    contact cutoff around the focus residues. Waters are the occupied oxygens
    of residues in :data:`crevice.hydration.WATER_RESNAMES`.

    Parameters
    ----------
    analysis : dict
        Result of :func:`crevice.hydration.static_hydration`; updated in place.
    frame : StructureFrame
        The structure the result was computed from.
    box : array_like, optional
        Unit cell for periodic distances.
    volume_path : str or pathlib.Path, optional
        Cavity map defining the density neighbourhood.
    options : dict, optional
        Settings for :func:`defaults`.

    Returns
    -------
    dict
        ``analysis`` with ``water_density`` and ``interactions`` (unless
        skipped), plus reference-frame entries used by the viewer writers.
    """
    import numpy as np
    from .models import StructureFrame
    from .hydration import WATER_RESNAMES
    from .water_density import WaterDensity
    from .interaction_analysis import chemical_topology_from_frame,InteractionAnalysis
    settings = defaults(options)
    protein_indices = [i for i,a in enumerate(frame.atoms) if not a.hetero and a.element.upper() not in {'H','D'}]
    protein = StructureFrame(tuple(frame.atoms[i] for i in protein_indices),frame.source,frame.model_index,frame.frame_index)
    waters = [i for i,a in enumerate(frame.atoms) if a.resname.upper() in WATER_RESNAMES and a.element.upper()=='O' and (a.occupancy is None or a.occupancy>0)]
    xyz = np.asarray([a.coord for a in frame.atoms])
    if not settings['skip_density']:
        points,definition = reference_points(protein,analysis['residue_ids'],volume_path)
        density = WaterDensity(points,margin=2. if volume_path else analysis['settings']['contact_cutoff_A'],
                     spacing=settings['density_spacing'],smoothing=settings['density_smoothing'],
                     max_points=settings['density_max_points'],definition=definition)
        density.add(xyz[waters],box=box)
        analysis['water_density'] = density.result(len(waters))
    if not settings['skip_interactions']:
        chemistry = chemical_topology_from_frame(frame,waters,box=box)
        interactions = InteractionAnalysis(chemistry,analysis['residue_ids'],1)
        interactions.add(xyz,box=box)
        analysis['interactions'] = interactions.result()
    analysis.update(reference_frame=protein,reference_atom_indices=np.asarray(protein_indices),
                    reference_full_coordinates=xyz,reference_rotation=np.eye(3),reference_box=box,
                    snapshot_frames=[protein],snapshot_times_ps=[None],analysis_options=settings,
                    reference_volume_path=str(volume_path) if volume_path else None)
    return analysis

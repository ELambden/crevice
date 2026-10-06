"""Explicit physical obstacles, independent of protein display/network selection.

In trajectory cavity analysis (``crevice cavity-trajectory`` with
reference-region geometry), the space measured in a fixed reference
neighbourhood can be bounded by more than the protein: lipids, ligands, ions
or waters can be declared as van der Waals obstacles with an MDAnalysis
selection. This module reads those atoms frame by frame, moves them into the
aligned reference frame, and enumerates periodic images so that no obstacle
near the region is missed.

Main entry points: :class:`ObstacleSource` and
:func:`aligned_obstacle_images`.

Limitations
    Obstacles are hard spheres with element radii. Treating lipids or
    waters as obstacles is a modelling choice; the result is not a
    substrate-accessibility or permeation model.
"""
from dataclasses import replace
from itertools import product
import os


def aligned_obstacle_images(atoms, coordinates, lower, upper, *, rotation, translation,
                            padding, box=None):
    """All atom images in a padded aligned box; retain source atom identities.

    The caller supplies a conservative clearance bound as padding. No nearest-
    image approximation is made when the box spans more than one periodic cell.

    Coordinates are transformed with the alignment convention
    ``x @ rotation + translation``. Without a periodic box, atoms inside the
    padded domain ``[lower - padding, upper + padding]`` are kept. With a box,
    atoms are first wrapped to the cell centred on the domain and then every
    lattice translation that can reach the domain is enumerated, so an atom can
    appear more than once.

    Parameters
    ----------
    atoms : sequence of Atom
        Topology atoms (identities are copied to the images).
    coordinates : array_like, shape (n, 3)
        Current positions (Å) in the simulation frame.
    lower, upper : array_like, shape (3,)
        Corners of the analysis domain in the aligned frame (Å).
    rotation : array_like, shape (3, 3)
        Proper rotation of the alignment (row convention).
    translation : array_like, shape (3,)
        Translation of the alignment (Å).
    padding : float
        Margin added to the domain (Å); must be finite and non-negative.
    box : array_like, shape (6,), optional
        MDAnalysis unit cell ``[a, b, c, alpha, beta, gamma]`` (Å, degrees).

    Returns
    -------
    atoms : tuple of Atom
        Retained images with aligned coordinates.
    report : dict
        Number of periodic images enumerated, atoms retained, whether the box
        was periodic and the padding used.

    Raises
    ------
    ValueError
        For mismatched or non-finite coordinates, an invalid domain or padding,
        an improper transform, an invalid cell, or a domain spanning more than
        343 periodic cells.
    """
    import numpy as np
    xyz = np.asarray(coordinates, dtype=float)
    rotation, translation = np.asarray(rotation), np.asarray(translation)
    lower, upper = np.asarray(lower)-padding, np.asarray(upper)+padding
    if xyz.shape != (len(atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError('Obstacle coordinates must be finite and match the topology')
    if not np.isfinite(padding) or padding < 0 or np.any(upper < lower):
        raise ValueError('Invalid obstacle domain or padding')
    if rotation.shape != (3,3) or translation.shape != (3,) or not np.isfinite(translation).all() or not np.allclose(rotation@rotation.T,np.eye(3),atol=1e-7) or not np.isclose(np.linalg.det(rotation),1):
        raise ValueError('Obstacle alignment requires a finite proper rigid transform')
    if box is None:
        chunks = [xyz@rotation+translation]
    else:
        from .hydration import checked_box
        from MDAnalysis.lib.mdamath import triclinic_vectors
        checked_box(box, 0.)
        cell = triclinic_vectors(box).astype(float); inverse = np.linalg.inv(cell)
        anchor = (lower+upper)/2
        fractional = (xyz-(anchor-translation)@rotation.T)@inverse
        central = (fractional-np.floor(fractional+.5))@cell
        corners = np.asarray(list(product(*zip(lower,upper))))
        fractions = ((corners-anchor)@rotation.T)@inverse
        lo, hi = np.floor(fractions.min(0)+.5).astype(int), np.floor(fractions.max(0)+.5).astype(int)
        if int(np.prod(hi-lo+1)) > 343:
            raise ValueError('Obstacle domain spans too many periodic cells')
        chunks = ((central+np.asarray(offset)@cell)@rotation+anchor
                  for offset in product(*(range(a,b+1) for a,b in zip(lo,hi))))
    result = []; images = 0
    for coords in chunks:
        images += 1
        keep = np.flatnonzero(np.all((coords>=lower)&(coords<=upper),axis=1))
        result.extend(replace(atoms[int(i)],x=float(coords[i,0]),y=float(coords[i,1]),z=float(coords[i,2])) for i in keep)
    return tuple(result), {'periodic_images_enumerated':images,'retained_atom_images':len(result),
                          'periodic':box is not None,'clearance_bound_padding_A':float(padding)}


class ObstacleSource:
    """Read only one physical-obstacle snapshot at a time, reopening after fork.

    Heavy atoms (not H or D) of ``obstacle_selection`` become obstacles. The
    selection must contain every analysed protein heavy atom; those reuse the
    reference protein's :class:`~crevice.models.Atom` records, and every other atom gets a
    separate chain namespace (``"<chain>/obstacle/<residue index>"``,
    ``hetero=True``) so solvent or ligand residue numbers cannot collide with
    protein residues. The MDAnalysis Universe is dropped when pickled and
    reopened in a new process.

    Parameters
    ----------
    topology, trajectory : str or pathlib.Path
        The same MD files as the protein analysis.
    protein_selection : str
        MDAnalysis selection used for the protein frames.
    obstacle_selection : str
        MDAnalysis selection of obstacle atoms (for example
        ``"protein or resname POPC"``).
    reference_frame : StructureFrame
        The protein reference frame; its atom count must match
        ``protein_selection``.
    frame_times : mapping of int to float
        Source frame index to time (ps) of the analysed protein frames.
    pbc : str, default "check"
        ``"none"`` disables periodic images; any other value uses the box.

    Raises
    ------
    ValueError
        If the selection has no heavy atoms, atom counts differ, protein heavy
        atoms are missing from the obstacle selection, or identities differ.
    """
    def __init__(self, topology, trajectory, *, protein_selection, obstacle_selection,
                 reference_frame, frame_times, pbc='check'):
        import numpy as np
        from collections import Counter
        from .trajectory import _open_universe, _element_names, _chain_labels
        from .models import Atom
        from .radii import atom_vdw_radius
        self.topology, self.trajectory = str(topology), str(trajectory)
        self.selection, self.pbc = obstacle_selection, pbc
        self.frame_times = dict(frame_times)
        self._universe = _open_universe(topology,trajectory); self._pid = os.getpid()
        u = self._universe
        display = u.select_atoms(protein_selection)
        chosen = u.select_atoms(obstacle_selection)
        elements = np.asarray(_element_names(chosen))
        chosen = chosen[~np.isin(elements,['H','D'])]; elements = elements[~np.isin(elements,['H','D'])]
        from .radii import UnrecognisedElementWarning, effective_radii, unrecognised_element_message
        import warnings
        active = effective_radii()
        known = np.asarray([active.has_radius(str(a.name), str(a.resname), str(e)) for a, e in zip(chosen, elements)], dtype=bool)
        if not known.all():
            missing_e = dict(sorted(Counter(str(e) for e in elements[~known]).items()))
            missing_n = dict(sorted(Counter(f'{a.resname}:{a.name}' for a in chosen[~known]).items()))
            warnings.warn(f'{trajectory} (obstacle selection): ' + unrecognised_element_message(missing_e, missing_n),
                          UnrecognisedElementWarning, stacklevel=2)
            chosen = chosen[known]; elements = elements[known]
        if not len(chosen):
            raise ValueError('Obstacle selection contains no heavy atoms')
        if len(display)!=len(reference_frame.atoms):
            raise ValueError('Obstacle source and protein reference atom counts differ')
        protein = {int(i):a for i,a in zip(display.indices,reference_frame.atoms) if not a.hetero and a.element.upper() not in {'H','D'}}
        if not set(protein).issubset(set(map(int,chosen.indices))):
            raise ValueError('Obstacle selection must retain every analysed protein heavy atom')
        chains = _chain_labels(chosen)[1]
        self.indices = chosen.indices.copy(); self.atoms = []
        for i,(a,element,chain) in enumerate(zip(chosen,elements,chains)):
            if int(a.index) in protein:
                atom = protein[int(a.index)]
                if atom.name != str(a.name) or atom.resid != int(a.resid):
                    raise ValueError('Obstacle topology and reference protein identities differ')
            else:
                # Separate namespace prevents solvent/ligand residue-number collisions.
                atom = Atom(int(a.index)+1,str(a.name),str(a.resname),f'{chain}/obstacle/{int(a.resindex)}',int(a.resid),0.,0.,0.,str(element),hetero=True)
            self.atoms.append(atom)
        self.atoms = tuple(self.atoms)
        self.radii = np.asarray([atom_vdw_radius(a) for a in self.atoms])
        self.metadata = {'selection':obstacle_selection,'protein_selection':protein_selection,
            'selected_heavy_atom_count':len(chosen),'analysed_protein_heavy_atom_count':len(protein),
            'other_heavy_atom_count':len(chosen)-len(protein),
            'residue_atom_counts':dict(sorted(Counter(map(str,chosen.resnames)).items())),
            'topology_atom_indices':list(map(int,self.indices)),
            'storage':'one obstacle snapshot per process; same source indices/times as protein',
            'periodicity':'enumerated lattice images in the aligned regional domain' if pbc!='none' else 'disabled explicitly',
            'scope':'reference-region geometry only; selected atoms are van der Waals obstacles, not a substrate-accessibility model'}

    def __getstate__(self):
        state = self.__dict__.copy(); state['_universe'] = None; state['_pid'] = None
        return state

    def snapshot(self, index):
        """Obstacle coordinates at one source frame.

        Parameters
        ----------
        index : int
            Source trajectory frame index; must be one of the analysed frames.

        Returns
        -------
        coordinates : ndarray, shape (n, 3)
            Obstacle positions in Å (simulation frame).
        box : ndarray or None
            Unit cell, or ``None`` with ``pbc="none"`` or no box.

        Raises
        ------
        ValueError
            If the frame was not analysed or its time differs from the protein
            frame by more than 1e-5 ps.
        """
        import numpy as np
        from .trajectory import _open_universe
        if self._pid != os.getpid() or self._universe is None:
            self._universe = _open_universe(self.topology,self.trajectory); self._pid = os.getpid()
        if index not in self.frame_times:
            raise ValueError('Obstacle frame is absent from the declared protein analysis')
        ts = self._universe.trajectory[index]
        if not np.isclose(float(ts.time),self.frame_times[index],rtol=0,atol=1e-5):
            raise ValueError('Obstacle and protein timestamps differ')
        box = ts.dimensions.copy() if self.pbc!='none' and ts.dimensions is not None else None
        return self._universe.atoms[self.indices].positions.copy(), box

    def prepare(self, frame, aligned, fit, reference, *, margin, probe_radius, max_grid_points=2_000_000):
        """Aligned obstacle images for one frame's reference-region grid.

        The padding is a conservative bound: the largest protein-only clearance at
        any regional grid sample (at least ``probe_radius + spacing`` and 1.5 Å)
        plus the largest obstacle radius. An atom further away than that cannot be
        the nearest sphere at any sample or block a segment between samples.

        Parameters
        ----------
        frame : StructureFrame
            Unaligned protein frame (supplies ``frame_index``).
        aligned : StructureFrame
            The same frame after alignment to the reference.
        fit : dict
            Alignment report with ``rotation_rows`` and ``translation_A``.
        reference : CavityReference
            Fixed reference region.
        margin : float
            Reference-neighbourhood margin (Å).
        probe_radius : float
            Probe radius (Å).
        max_grid_points : int, default 2000000

        Returns
        -------
        atoms : tuple of Atom
        report : dict
            See :func:`aligned_obstacle_images`.

        Raises
        ------
        ValueError
            If no obstacle image intersects the regional domain, or as for
            :meth:`snapshot`.
        """
        import numpy as np
        from .regional_volume import regional_grid
        from .rolling import _SphereQueries
        from .radii import atom_vdw_radius
        _,xyz,_,_,_,_ = regional_grid(reference,margin,probe_radius,max_grid_points)
        protein = aligned.selected_atoms(include_hydrogen=False,include_hetero=False)
        q = _SphereQueries(np.asarray([a.coord for a in protein]),np.asarray([atom_vdw_radius(a) for a in protein]))
        # Upper bound from protein alone: an excluded atom cannot be the nearest
        # sphere at any grid sample or block an eligible segment between samples.
        maximum = max(float(q.points(xyz)[0].max()),probe_radius+reference.spacing,1.5)
        coords, box = self.snapshot(frame.frame_index)
        atoms, report = aligned_obstacle_images(self.atoms,coords,xyz.min(0),xyz.max(0),
            rotation=fit['rotation_rows'],translation=fit['translation_A'],
            padding=maximum+float(self.radii.max())+1e-6,box=box)
        if not atoms:
            raise ValueError('No physical obstacle images intersect the regional domain')
        return atoms, report

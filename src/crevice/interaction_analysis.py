"""All-frame chemically typed residue networks and conformation summaries.

:class:`InteractionAnalysis` accumulates the per-frame observations of
:class:`crevice.interaction_chemistry.ChemicalTopology` (hydrogen bonds, salt
bridge candidates, aromatic contacts, water bridges, torsions, secondary
structure) over a trajectory. :func:`write_interaction_bundle` turns them
into occupancies (fraction of frames in which each typed residue pair was
observed), circular torsion statistics and output files.

Interpretation
    Occupancy is the fraction of analysed frames with a detected,
    geometry-qualified observation. Missing chemistry (for example no
    explicit hydrogens) is not a negative observation, and occupancies are
    descriptive, not interaction energies or evidence of functional control.
"""
from pathlib import Path
from collections import defaultdict
import csv
import json


def chemical_topology_from_frame(frame, water_indices, *, box=None):
    """Chemical topology for a single static structure.

    Protein atoms are the non-HETATM atoms; no topology bonds are used, so
    hydrogens are attached by distance.

    Parameters
    ----------
    frame : StructureFrame
    water_indices : sequence of int
        Indices (into ``frame.atoms``) of water oxygens.
    box : array_like, optional
        Unit cell for periodic distances.

    Returns
    -------
    ChemicalTopology
    """
    from .interaction_chemistry import ChemicalTopology
    protein = [i for i,a in enumerate(frame.atoms) if not a.hetero]
    return ChemicalTopology([a.name for a in frame.atoms], [a.element.upper() for a in frame.atoms],
        [a.residue_key.label for a in frame.atoms], [a.resname for a in frame.atoms], protein,
        water_indices, [a.coord for a in frame.atoms], box=box, chains=[a.chain_id for a in frame.atoms])


def chemical_topology_from_universe(universe, selected, reference_frame, water_indices):
    """Chemical topology for all atoms of an MDAnalysis Universe.

    Protein residues take their labels from ``reference_frame``; every other
    residue is labelled ``"water-resindex:<index>"``. Topology bonds are used
    when available (GRO files usually have none), and segment IDs serve as
    chains.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Positioned at the reference frame.
    selected : MDAnalysis.AtomGroup
        The protein atoms of ``reference_frame``, in the same order.
    reference_frame : StructureFrame
    water_indices : sequence of int
        Universe indices of water oxygens.

    Returns
    -------
    ChemicalTopology

    Raises
    ------
    ValueError
        If ``selected`` does not match ``reference_frame`` atom for atom.
    """
    import numpy as np
    from .trajectory import _element_names
    from MDAnalysis.exceptions import NoDataError
    from .interaction_chemistry import ChemicalTopology
    atoms = universe.atoms
    if len(selected)!=len(reference_frame.atoms) or list(selected.names)!=[a.name for a in reference_frame.atoms]:
        raise ValueError('Chemical source atom selection differs from the protein frame')
    labels = np.asarray(['water-resindex:'+str(i) for i in atoms.resindices],dtype=object)
    mapping = {int(atoms.resindices[index]):a.residue_key.label for index,a in zip(selected.indices,reference_frame.atoms)}
    for resindex,key in mapping.items():
        labels[atoms.resindices==resindex] = key
    protein_indices = np.flatnonzero(np.isin(atoms.resindices,list(mapping)))
    try:
        bonds = atoms.bonds.indices
    except (AttributeError, NoDataError):
        # GRO/PDB topologies commonly have no bonds; fallback is explicit in metadata.
        bonds = None
    return ChemicalTopology(atoms.names,_element_names(atoms),labels,atoms.resnames,protein_indices,
                            water_indices,atoms.positions,box=universe.dimensions,
                            chains=atoms.segids,bonds=bonds)


class InteractionAnalysis:
    """Accumulate typed interactions and torsions frame by frame.

    Parameters
    ----------
    topology : ChemicalTopology
    residue_ids : iterable of str
        Focus residues (sorted internally); all must be protein residues of the
        topology.
    expected_frames : int
        Number of frames that will be added; :meth:`result` requires exactly
        this many.

    Raises
    ------
    ValueError
        If a focus residue is unknown.
    """
    def __init__(self, topology, residue_ids, expected_frames):
        import numpy as np
        self.topology = topology
        self.ids = sorted(residue_ids)
        self.lookup = {key:i for i,key in enumerate(self.ids)}
        if set(self.ids)-set(topology.residues):
            raise ValueError('Interaction focus contains unknown protein residues')
        self.expected_frames = int(expected_frames)
        self.frame = 0
        self.torsion_names = ['phi','psi','chi1','chi2','chi3','chi4']
        self.torsions = np.full((expected_frames,len(self.ids),6),np.nan)
        self.secondary = np.full((expected_frames,len(self.ids)),'NA',dtype='U2')
        self.water_counts = np.zeros((expected_frames,len(self.ids)),dtype=int)
        self.water_contacts = []
        self.edges = {}

    def add(self, coordinates, *, box=None):
        """Measure one frame and add it to the accumulators.

        Parameters
        ----------
        coordinates : array_like, shape (n, 3)
            Coordinates of all topology atoms (Å).
        box : array_like, optional
            Unit cell for periodic distances.

        Raises
        ------
        ValueError
            If more than ``expected_frames`` frames are added.
        """
        from .interaction_chemistry import ChemicalTopology
        if self.frame>=self.expected_frames:
            raise ValueError('Too many chemical interaction frames')
        measured = self.topology.measure(coordinates,self.ids,box=box)
        grouped = defaultdict(list)
        for observation in measured['observations']:
            grouped[tuple(observation[k] for k in ['kind','source','target'])].append(observation)
        for key, observations in grouped.items():
            if key not in self.edges:
                self.edges[key] = {'frames':[],'counts':[],'distances':[],'angles':[],
                                  'example_observation':observations[0], 'example_frame_row':self.frame,
                                  'reference_observation':observations[0] if self.frame==0 else None}
            edge = self.edges[key]
            edge['frames'].append(self.frame)
            edge['counts'].append(len(observations))
            edge['distances'].append(min(o['distance_A'] for o in observations))
            angles = [o['angle_deg'] for o in observations if o['angle_deg'] is not None]
            edge['angles'].extend(angles)
        for key,row in measured['torsions'].items():
            j = self.lookup[key]
            self.torsions[self.frame,j] = [row[k] for k in self.torsion_names]
            self.secondary[self.frame,j] = measured['secondary'][key]
        waters = defaultdict(set)
        for key,water,bond in measured['water_bonds']:
            waters[key].add(water)
        for key, observed in waters.items():
            j = self.lookup[key]
            self.water_counts[self.frame,j] = len(observed)
            self.water_contacts.extend((self.frame,j,int(w)) for w in sorted(observed))
        self.frame += 1

    def result(self):
        """Collected per-frame results.

        Returns
        -------
        dict
            ``residue_ids``, ``frame_count``, topology ``metadata``, ``edges``
            (keyed by ``(kind, source, target)``: frame rows, per-frame counts and
            minimum distances, observed angles and example observations),
            ``torsions_degrees`` (frames x residues x 6), ``torsion_names``,
            ``secondary_structure``, ``hydrogen_bonded_water_count`` (NaN when no
            waters or no explicit donor hydrogens), ``hydrogen_bonded_water_contacts``
            and per-atom residue IDs and names.

        Raises
        ------
        ValueError
            If fewer frames than ``expected_frames`` were added.
        """
        import numpy as np
        if self.frame!=self.expected_frames:
            raise ValueError('Chemical interaction frame coverage is incomplete')
        metadata=self.topology.metadata()
        water_counts=self.water_counts.astype(float)
        if not metadata['water_oxygen_population'] or metadata['hydrogen_bond_status']=='unavailable_no_explicit_donor_hydrogens':water_counts[:]=np.nan
        return {'residue_ids':self.ids,'frame_count':self.frame,'metadata':metadata,
                'edges':self.edges,'torsions_degrees':self.torsions,'torsion_names':self.torsion_names,
                'secondary_structure':self.secondary,'hydrogen_bonded_water_count':water_counts,
                'hydrogen_bonded_water_contacts':np.asarray(self.water_contacts,dtype=int).reshape((-1,3)),
                'atom_residue_ids':self.topology.residue_ids,'atom_names':self.topology.names}


def write_interaction_bundle(result, output_dir, prefix, *, times=None, hydration_arrays=None,
                             confidence=.95, block_length=None, replicates=2000, seed=20260914):
    """Write interaction occupancies and conformation summaries to files.

    Parameters
    ----------
    result : dict
        From :meth:`InteractionAnalysis.result`.
    output_dir : str or pathlib.Path
        Existing output directory.
    prefix : str
        File-name prefix.
    times : array_like, optional
        Frame times (ps). Confidence intervals are computed only when times are
        given and regularly spaced.
    hydration_arrays : dict, optional
        Hydration ``arrays`` (with ``residue_ids`` and ``water_count``); adds
        chi1 step statistics split by wet and dry frame pairs.
    confidence : float, default 0.95
        Confidence level for occupancy intervals; a false value disables them.
    block_length : int, optional
        Batch length in frames for the bootstrap (automatic when omitted).
    replicates : int, default 2000
        Bootstrap replicates.
    seed : int, default 20260914
        Bootstrap random seed.

    Returns
    -------
    dict
        Paths of the files written (see Outputs).

    Outputs
    -------
    <prefix>_interactions.json : JSON
        Full structured report (edge occupancy statistics with example and
        reference observations, conformation summaries, chemistry metadata),
        read by the analysis scenes and ``region-compare``.
    <prefix>_interaction_edges.csv : table
        One row per typed residue pair: ``kind``, ``source``, ``target``,
        observed and denominator frames, ``observed_occupancy``,
        ``mean_minimum_contact_distance_A``, ``mean_observed_angle_degrees``,
        ``example_frame_row`` and the occupancy interval (header only when no
        typed contact is observed).
    <prefix>_interaction_edge_frames.csv : table
        One row per frame and observed edge: ``frame_row`` (0-based analysed
        frame), ``time_ps`` (empty for a static structure), ``kind``,
        ``source``, ``target`` and ``observation_count``.
    <prefix>_residue_conformations.csv : table
        Per-residue torsion circular means, chi1 rotamer and secondary-structure
        fractions.
    <prefix>_residue_conformation_frames.csv : table
        One row per frame and residue: ``frame_row``, ``time_ps``, one
        ``<torsion>_degrees`` column per torsion (phi, psi, chi1-chi4; empty
        where undefined), ``secondary_structure`` (``H``, ``E``, ``C`` or
        ``NA``) and ``hydrogen_bonded_water_count``.
    <prefix>_interaction_frames.npz : arrays
        The same per-frame arrays plus atom identities and sparse edge frames.

    Notes
    -----
    Per residue, torsions are summarised by the circular mean and mean
    resultant length (1 = no spread). chi1 is binned as
    ``gauche_plus`` [0, 120), ``trans`` [120, 240) and ``gauche_minus``
    [240, 360) degrees; these coarse bins are not a rotamer library.
    Secondary-structure fractions use only frames with a DSSP label. Wet/dry
    chi1 steps are descriptive associations with no causal or kinetic claim.
    """
    import numpy as np
    from .cavity_trajectory import describe_series
    root = Path(output_dir)
    n, ids = result['frame_count'], result['residue_ids']
    regular = times is not None and (n<3 or np.allclose(np.diff(times),np.diff(times)[0]))
    files = {}
    def path(key,suffix):
        p = root/(prefix+suffix); files[key] = str(p); return p
    def stats(series, bounded=False):
        stat = describe_series(series,confidence=confidence or .95,block_length=block_length,
                               replicates=replicates,seed=seed,regular=regular and bool(confidence))
        if not confidence:
            stat['confidence_interval'] = {'status':'not_requested'}
        for key in ['pointwise_lower','pointwise_upper','simultaneous_lower','simultaneous_upper']:
            ci = stat['confidence_interval']
            if key in ci:
                ci[key] = [float(np.clip(v,0,1)) if bounded and v is not None else max(0.,v) if v is not None else None for v in ci[key]]
        return stat
    edges, edge_rows, sparse = [], [], []
    for edge_number,(key,raw) in enumerate(sorted(result['edges'].items())):
        kind,left,right = key
        occupancy = np.zeros(n); occupancy[raw['frames']] = 1
        stat = stats(occupancy,True)
        record = {'kind':kind,'source':left,'target':right,'observed_frames':len(raw['frames']),
                  'denominator_frames':n,'observed_occupancy':len(raw['frames'])/n,'occupancy_statistics':stat,
                  'mean_minimum_contact_distance_A':float(np.mean(raw['distances'])),
                  'mean_observed_angle_degrees':float(np.mean(raw['angles'])) if raw['angles'] else None,
                  'example_observation':raw['example_observation'],'example_frame_row':raw['example_frame_row'],
                  'reference_observation':raw['reference_observation'],
                  'interpretation':'occupancy of detected template/geometry-qualified observations; missing chemistry is not a negative observation'}
        edges.append(record)
        row = {k:v for k,v in record.items() if k not in {'occupancy_statistics','example_observation','reference_observation','interpretation'}}
        ci = stat['confidence_interval']; row.update(CI_status=ci['status'],
             mean_occupancy_lower=ci.get('pointwise_lower',[None])[0],mean_occupancy_upper=ci.get('pointwise_upper',[None])[0])
        edge_rows.append(row)
        sparse.extend((frame,edge_number,count) for frame,count in zip(raw['frames'],raw['counts']))
    conformations = []
    for j,key in enumerate(ids):
        row = {'residue':key}
        for k,name in enumerate(result['torsion_names']):
            values = result['torsions_degrees'][:,j,k]
            good = np.isfinite(values)
            if good.any():
                z = np.exp(1j*np.radians(values[good])).mean()
                row[name+'_circular_mean_degrees'] = float(np.degrees(np.angle(z))) if abs(z)>1e-8 else None
                row[name+'_mean_resultant_length'] = float(abs(z))
            else:
                row[name+'_circular_mean_degrees'] = row[name+'_mean_resultant_length'] = None
            row[name+'_observed_frames'] = int(good.sum())
        chi = result['torsions_degrees'][:,j,2]
        valid = np.isfinite(chi)
        bins = np.floor(np.mod(chi[valid],360)/120).astype(int)
        for k,label in enumerate(['gauche_plus','trans','gauche_minus']):
            row['chi1_'+label+'_fraction'] = float(np.mean(bins==k)) if len(bins) else None
        secondary = result['secondary_structure'][:,j]
        eligible = secondary!='NA'
        row['secondary_structure_observed_frames'] = int(eligible.sum())
        for label,code in [('helix','H'),('sheet','E'),('coil','C')]:
            row[label+'_fraction'] = float(np.mean(secondary[eligible]==code)) if eligible.any() else None
        row['mean_hydrogen_bonded_waters'] = float(np.nanmean(result['hydrogen_bonded_water_count'][:,j])) if np.isfinite(result['hydrogen_bonded_water_count'][:,j]).any() else None
        if hydration_arrays is not None and key in list(hydration_arrays['residue_ids']):
            hi = list(hydration_arrays['residue_ids']).index(key)
            water = hydration_arrays['water_count'][:,hi]
            delta = np.abs(np.degrees(np.angle(np.exp(1j*np.radians(np.diff(chi))))))
            good = np.isfinite(delta)&np.isfinite(water[:-1])&np.isfinite(water[1:])
            wet = good&(water[:-1]>0)&(water[1:]>0)
            dry = good&(water[:-1]==0)&(water[1:]==0)
            row['chi1_wet_endpoint_intervals'] = int(wet.sum()); row['chi1_dry_endpoint_intervals'] = int(dry.sum())
            row['mean_chi1_step_wet_endpoints_degrees'] = float(delta[wet].mean()) if wet.any() else None
            row['mean_chi1_step_dry_endpoints_degrees'] = float(delta[dry].mean()) if dry.any() else None
        conformations.append(row)
    report = {'frame_count':n,'residue_count':len(ids),'metadata':result['metadata'],'edges':edges,
              'residues':conformations,'residue_ids':ids,
              'torsion_interpretation':'circular summaries and coarse chi1 bins; endpoint hydration/step associations are descriptive, with no causal or kinetic claim',
              'confidence_interpretation':'pointwise approximate occupancy intervals only; no simultaneous family-wide confidence or chemistry/systematic uncertainty'}
    path('interaction_json','_interactions.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    edge_fields = ['kind','source','target','observed_frames','denominator_frames','observed_occupancy',
                   'mean_minimum_contact_distance_A','mean_observed_angle_degrees','example_frame_row',
                   'CI_status','mean_occupancy_lower','mean_occupancy_upper']
    time_of = (lambda row: float(times[row])) if times is not None else (lambda row: '')
    edge_keys = sorted(result['edges'])
    edge_frame_rows = [{'frame_row':int(frame),'time_ps':time_of(frame),'kind':edge_keys[number][0],
                        'source':edge_keys[number][1],'target':edge_keys[number][2],'observation_count':int(count)}
                       for frame,number,count in sparse]
    torsion_names = list(result['torsion_names'])
    conformation_frame_rows = []
    for i in range(n):
        for j,key in enumerate(ids):
            row = {'frame_row':i,'time_ps':time_of(i),'residue':key}
            for k,name in enumerate(torsion_names):
                value = float(result['torsions_degrees'][i,j,k])
                row[name+'_degrees'] = value if np.isfinite(value) else ''
            row['secondary_structure'] = str(result['secondary_structure'][i,j])
            water = float(result['hydrogen_bonded_water_count'][i,j])
            row['hydrogen_bonded_water_count'] = water if np.isfinite(water) else ''
            conformation_frame_rows.append(row)
    for key,suffix,rows,fields in [
            ('interaction_edges_csv','_interaction_edges.csv',edge_rows,edge_fields),
            ('residue_conformations_csv','_residue_conformations.csv',conformations,None),
            ('interaction_edge_frames_csv','_interaction_edge_frames.csv',edge_frame_rows,
             ['frame_row','time_ps','kind','source','target','observation_count']),
            ('residue_conformation_frames_csv','_residue_conformation_frames.csv',conformation_frame_rows,None)]:
        p = path(key,suffix)
        with p.open('w',newline='',encoding='utf-8') as f:
            if rows or fields:
                fields = fields or list(dict.fromkeys(k for r in rows for k in r))
                writer = csv.DictWriter(f,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    np.savez_compressed(path('interaction_frames_npz','_interaction_frames.npz'),
                        residue_ids=np.asarray(ids),torsion_names=np.asarray(result['torsion_names']),
                        torsions_degrees=result['torsions_degrees'],secondary_structure=result['secondary_structure'],
                        hydrogen_bonded_water_count=result['hydrogen_bonded_water_count'],
                        hydrogen_bonded_water_contacts=result['hydrogen_bonded_water_contacts'],
                        atom_residue_ids=result['atom_residue_ids'],atom_names=result['atom_names'],
                        edge_frame_indices=np.asarray(sparse,dtype=int).reshape((-1,3)),
                        edge_keys=np.asarray(list(sorted(result['edges'])),dtype=str).reshape((-1,3)))
    return files

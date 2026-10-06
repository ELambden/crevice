"""Conservative residue chemistry and explicit-H directional interaction evidence.

Templates describe standard amino acids and named protonation variants. No
hydrogens, protonation states, ligand chemistry or interaction energies are inferred.

What is detected (between residues of a focus set, per frame)
    ``hydrogen_bond``
        Donor-H...acceptor with an explicit hydrogen: D-H bond 0.4-1.25 Å
        (1.65 Å for S), D...A at most 3.5 Å, H...A at most 2.5 Å and angle
        D-H...A at least 150 degrees, between different residues.
    ``salt_bridge_candidate``
        Oppositely charged side-chain atoms (ARG NE/NH1/NH2, LYS NZ,
        doubly protonated HIS vs ASP/GLU carboxylate O, CYM SG) within
        4 Å. Formal-charge candidates only.
    ``aromatic_parallel_candidate`` / ``aromatic_edge_face_candidate``
        Complete, planar (RMS at most 0.15 Å) rings with centroids 2.5-5.5 Å
        apart; parallel if the normals are within 30 degrees and the larger
        lateral offset is at most 2 Å; edge-face if at least 60 degrees and
        the smaller offset is at most 2 Å.
    ``water_bridge_hbond``
        Two residues each hydrogen-bonded (as above, with explicit
        hydrogens) to the same water.

Also measured: backbone phi/psi and side-chain chi1-chi4 torsions, and
simplified DSSP secondary structure (H/E/C) through MDTraj when installed.

Residue templates cover the 20 standard amino acids, the protonation
variants in :data:`PARENTS` (for example HSD/HSE/HSP, HID/HIE/HIP, ASH, GLH,
LYN, CYM, CYX) and the D-amino acids DLE and DVA. Residues outside the
templates are listed as unknown and get no typed atoms.

Limitations
    Without explicit hydrogens no hydrogen bond can be observed; the missing
    observation is not evidence of absence. Protonation states are taken
    from residue names, never guessed. Geometry criteria are conventional
    cut-offs, not energies.
"""
from collections import defaultdict
from itertools import combinations
import math

#: Protonation-variant and D-amino-acid residue names mapped to their standard parent.
PARENTS = {'HSD':'HIS', 'HSE':'HIS', 'HSP':'HIS', 'HID':'HIS', 'HIE':'HIS',
           'HIP':'HIS', 'ASH':'ASP', 'GLH':'GLU', 'LYN':'LYS', 'CYM':'CYS', 'CYX':'CYS'}
# CCD D-PEPTIDE LINKING components; only chemical atom templates are shared.
D_PARENTS = {'DLE':'LEU','DVA':'VAL'}
PARENTS.update(D_PARENTS)
STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
DONORS = {'ARG':['NE','NH1','NH2'], 'ASN':['ND2'], 'CYS':['SG'], 'GLN':['NE2'],
          'HIS':['ND1','NE2'], 'LYS':['NZ'], 'SER':['OG'], 'THR':['OG1'],
          'TRP':['NE1'], 'TYR':['OH'], 'ASP':['OD1','OD2'], 'GLU':['OE1','OE2']}
ACCEPTORS = {'ASN':['OD1'], 'ASP':['OD1','OD2'], 'CYS':['SG'], 'GLN':['OE1'],
             'GLU':['OE1','OE2'], 'HIS':['ND1','NE2'], 'MET':['SD'],
             'SER':['OG'], 'THR':['OG1'], 'TYR':['OH']}
RINGS = {'PHE':[['CG','CD1','CE1','CZ','CE2','CD2']],
         'TYR':[['CG','CD1','CE1','CZ','CE2','CD2']],
         'HIS':[['CG','ND1','CE1','NE2','CD2']],
         'TRP':[['CG','CD1','NE1','CE2','CD2'], ['CD2','CE2','CZ2','CH2','CZ3','CE3']]}
CHI = {
 'ARG': [('N','CA','CB','CG'),('CA','CB','CG','CD'),('CB','CG','CD','NE'),('CG','CD','NE','CZ')],
 'ASN': [('N','CA','CB','CG'),('CA','CB','CG','OD1')],
 'ASP': [('N','CA','CB','CG'),('CA','CB','CG','OD1')],
 'CYS': [('N','CA','CB','SG')],
 'GLN': [('N','CA','CB','CG'),('CA','CB','CG','CD'),('CB','CG','CD','OE1')],
 'GLU': [('N','CA','CB','CG'),('CA','CB','CG','CD'),('CB','CG','CD','OE1')],
 'HIS': [('N','CA','CB','CG'),('CA','CB','CG','ND1')],
 'ILE': [('N','CA','CB','CG1'),('CA','CB','CG1','CD1')],
 'LEU': [('N','CA','CB','CG'),('CA','CB','CG','CD1')],
 'LYS': [('N','CA','CB','CG'),('CA','CB','CG','CD'),('CB','CG','CD','CE'),('CG','CD','CE','NZ')],
 'MET': [('N','CA','CB','CG'),('CA','CB','CG','SD'),('CB','CG','SD','CE')],
 'PHE': [('N','CA','CB','CG'),('CA','CB','CG','CD1')],
 'PRO': [('N','CA','CB','CG'),('CA','CB','CG','CD')],
 'SER': [('N','CA','CB','OG')], 'THR':[('N','CA','CB','OG1')],
 'TRP': [('N','CA','CB','CG'),('CA','CB','CG','CD1')],
 'TYR': [('N','CA','CB','CG'),('CA','CB','CG','CD1')],
 'VAL': [('N','CA','CB','CG1')],
}


def periodic_vectors(vectors, box):
    """Minimum-image displacement vectors.

    Parameters
    ----------
    vectors : array_like, shape (n, 3)
        Displacements in Å.
    box : array_like or None
        Unit cell; ``None`` returns the vectors unchanged.

    Returns
    -------
    ndarray

    See Also
    --------
    crevice.hydration.minimum_vectors
    """
    from .hydration import minimum_vectors
    return minimum_vectors(vectors, box)


def torsion_degrees(points, box=None):
    """Dihedral angle of four points, in degrees.

    Parameters
    ----------
    points : array_like, shape (4, 3)
        Atom positions (Å) in bond order.
    box : array_like, optional
        Unit cell; bond vectors use the minimum image.

    Returns
    -------
    float
        Angle in ``(-180, 180]`` degrees (IUPAC sign convention), or NaN if any
        of the three bond vectors is shorter than 0.4 Å or longer than 2.2 Å
        (a broken or missing link) or the geometry is degenerate.

    Examples
    --------
    >>> from crevice.interaction_chemistry import torsion_degrees
    >>> round(torsion_degrees([(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.5, 0.0), (-1.0, 1.5, 0.0)]), 1)
    180.0
    """
    import numpy as np
    p = np.asarray(points, dtype=float)
    b = periodic_vectors(np.diff(p, axis=0), box)
    lengths = np.linalg.norm(b, axis=1)
    if not np.isfinite(b).all() or np.any((lengths < .4) | (lengths > 2.2)):
        return float('nan')
    axis = b[1]/lengths[1]
    left = -b[0]-np.dot(-b[0], axis)*axis
    right = b[2]-np.dot(b[2], axis)*axis
    if min(np.linalg.norm(left), np.linalg.norm(right)) < 1e-8:
        return float('nan')
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis,left),right), np.dot(left,right))))


class ChemicalTopology:
    """Typed donors, acceptors, charges and rings for a fixed set of atoms.

    Built once from topology and reference coordinates; :meth:`measure` is then
    called with the coordinates of each frame.

    Hydrogen attachment
        With ``bonds``, a hydrogen is attached to its single bonded heavy atom
        in the same residue. Without bonds, it is attached to the nearest heavy
        atom of its residue, unless the second nearest is within 0.15 Å of the
        same distance (ambiguous, skipped). Either way the H-heavy distance must
        be 0.4-1.25 Å (1.65 Å for S).

    Typing rules (standard residues only)
        Donors are template side-chain donors plus backbone N (proline only if
        an H is attached). Acidic OH and histidine donors count only where an
        explicit H is observed. Acceptors are backbone O/OXT/OT1/OT2 plus
        template acceptors, excluding atoms that carry an H, histidine ring N
        atoms whose tautomer is unresolved (HSD/HID accept at NE2 only,
        HSE/HIE at ND1 only, HSP/HIP at neither), protonated acids without an
        observed acid H, and CYX SG.

    Parameters
    ----------
    names, elements, residue_ids, resnames : sequence of str
        Per-atom atom name, element (upper case), residue label and residue
        name, for *all* atoms (protein, water and others).
    protein_indices : sequence of int
        Atoms belonging to protein residues.
    water_oxygen_indices : sequence of int
        Water oxygen atoms.
    reference_coordinates : array_like, shape (n, 3)
        Coordinates (Å) used for hydrogen attachment and peptide links.
    box : array_like, optional
        Unit cell for periodic distances.
    chains : sequence of str, optional
        Chain ID per atom; consecutive residues are linked only within a chain
        and when their C-N distance is at most 1.9 Å.
    bonds : sequence of (int, int), optional
        Topology bonds; used for hydrogen attachment when given.

    Attributes
    ----------
    residues : list of str
        Protein residue labels in input order.
    unknown_residues : list of str
        Residues without a template.
    ambiguous_residues : list of str
        Residues with duplicate atom names (excluded from typing).
    ambiguous_hydrogens : int
        Hydrogens that could not be attached unambiguously.
    donor_pairs : ndarray, shape (k, 2)
        ``(heavy atom, hydrogen)`` index pairs.
    acceptors : ndarray of int
        Indices of hydrogen-bond acceptor atoms.
    positive : ndarray of int
        Indices of positively charged atoms.
    negative : ndarray of int
        Indices of negatively charged atoms.
    rings : list
        ``(residue, ring number, atom indices)`` for complete rings.
    dssp_status : str
        Whether simplified DSSP is available.

    Raises
    ------
    ValueError
        If the per-atom arrays and coordinates differ in length.
    """
    def __init__(self, names, elements, residue_ids, resnames, protein_indices,
                 water_oxygen_indices, reference_coordinates, *, box=None, chains=None, bonds=None):
        import numpy as np
        self.names = np.asarray(names, dtype=str)
        self.elements = np.asarray(elements, dtype=str)
        self.residue_ids = np.asarray(residue_ids, dtype=str)
        self.resnames = np.asarray(resnames, dtype=str)
        self.protein_indices = np.asarray(protein_indices, dtype=int)
        self.water_oxygens = np.asarray(water_oxygen_indices, dtype=int)
        self.chains = np.asarray(chains if chains is not None else ['']*len(names), dtype=str)
        xyz = np.asarray(reference_coordinates, dtype=float)
        if any(len(x)!=len(xyz) for x in [self.names,self.elements,self.residue_ids,self.resnames,self.chains]):
            raise ValueError('Chemical topology and coordinate sizes differ')
        self.groups = defaultdict(list)
        for i in self.protein_indices:
            self.groups[str(self.residue_ids[i])].append(int(i))
        self.residues = list(self.groups)
        self.atom_names = {}
        self.unknown_residues = []
        self.ambiguous_residues = []
        for key, indices in self.groups.items():
            names_here = [str(self.names[i]) for i in indices]
            if len(names_here)!=len(set(names_here)):
                self.ambiguous_residues.append(key)
                self.atom_names[key] = {}
                continue
            self.atom_names[key] = dict(zip(names_here, indices))
            name = self.resnames[indices[0]]
            if PARENTS.get(name,name) not in STANDARD:
                self.unknown_residues.append(key)
        all_groups = defaultdict(list)
        for i, key in enumerate(self.residue_ids):
            all_groups[str(key)].append(i)
        relevant = set(self.residues) | set(self.residue_ids[self.water_oxygens])
        attached = defaultdict(list)
        ambiguous = 0
        self.attachment_method = 'supplied bonds' if bonds is not None else 'nearest heavy atom within the same residue, with covalent-distance and ambiguity guards'
        if bonds is not None:
            partners = defaultdict(list)
            for a,b in bonds:
                partners[int(a)].append(int(b)); partners[int(b)].append(int(a))
        for key in relevant:
            indices = all_groups[key]
            hs = [i for i in indices if self.elements[i] in {'H','D'}]
            heavy = [i for i in indices if self.elements[i] not in {'H','D'}]
            if not hs or not heavy:
                continue
            distances = np.linalg.norm(periodic_vectors((xyz[hs,None]-xyz[heavy]).reshape(-1,3),box).reshape(len(hs),len(heavy),3),axis=2)
            for row,h in enumerate(hs):
                if bonds is not None:
                    candidates = [a for a in partners[h] if a in heavy]
                    if len(candidates)!=1:
                        ambiguous += 1; continue
                    a = candidates[0]; distance = distances[row,heavy.index(a)]
                else:
                    order = np.argsort(distances[row]); a = heavy[order[0]]; distance = distances[row,order[0]]
                    if len(order)>1 and distances[row,order[1]]-distance < .15:
                        ambiguous += 1; continue
                maximum = 1.65 if self.elements[a]=='S' else 1.25
                if .4 < distance <= maximum:
                    attached[a].append(h)
        self.attached = dict(attached)
        self.ambiguous_hydrogens = ambiguous
        donor_pairs, acceptors, positive, negative, rings = [], [], [], [], []
        expected_donors = []
        for key in self.residues:
            mapping = self.atom_names[key]; name = str(self.resnames[self.groups[key][0]]); parent = PARENTS.get(name,name)
            if parent not in STANDARD or key in self.ambiguous_residues:
                continue
            dnames = list(DONORS.get(parent,[]))
            if name in {'CYM','CYX'}:
                dnames = []
            if parent!='PRO' or attached.get(mapping.get('N',-1)):
                dnames += ['N']
            anames = ['O','OXT','OT1','OT2'] + ACCEPTORS.get(parent,[])
            for atom_name in dnames:
                if atom_name not in mapping:
                    continue
                atom = mapping[atom_name]
                # Acidic OH and histidine donors exist only where an explicit H is observed.
                if parent not in {'ASP','GLU','HIS'} or attached.get(atom):
                    expected_donors.append(atom)
                donor_pairs += [(atom,h) for h in attached.get(atom,[])]
            for atom_name in anames:
                if atom_name not in mapping:
                    continue
                atom = mapping[atom_name]
                if name in {'ASH','GLH'} and atom_name in ACCEPTORS.get(parent,[]) and not any(attached.get(mapping.get(n,-1)) for n in ACCEPTORS.get(parent,[])):
                    continue  # Protonated acid without its H has unresolved OH identity.
                if parent=='HIS' and atom_name in {'ND1','NE2'}:
                    if attached.get(atom):
                        continue
                    # A generic HIS without any imidazole H has unresolved tautomerism.
                    fixed = {'HSD':'NE2','HID':'NE2','HSE':'ND1','HIE':'ND1'}.get(name)
                    other_h = any(attached.get(mapping.get(n,-1)) for n in ['ND1','NE2'])
                    if name in {'HSP','HIP'} or (fixed is None and not other_h) or (fixed and atom_name!=fixed):
                        continue
                if (parent in {'ASP','GLU'} or atom_name in {'OXT','OT1','OT2'}) and attached.get(atom):
                    continue
                if name=='CYX' and atom_name=='SG':
                    continue
                acceptors.append(atom)
            if parent=='ARG':
                positive += [mapping[n] for n in ['NE','NH1','NH2'] if n in mapping]
            elif parent=='LYS' and name!='LYN':
                positive += [mapping['NZ']] if 'NZ' in mapping else []
            elif name in {'HSP','HIP'}:
                positive += [mapping[n] for n in ['ND1','NE2'] if n in mapping]
            acid_names = {'ASP':['OD1','OD2'],'GLU':['OE1','OE2']}.get(parent,[])
            if name not in {'ASH','GLH'} and not any(attached.get(mapping.get(n,-1)) for n in acid_names):
                negative += [mapping[n] for n in acid_names if n in mapping]
            if name=='CYM' and 'SG' in mapping:
                negative.append(mapping['SG'])
            for ring_number, atom_names in enumerate(RINGS.get(parent,[])):
                if all(n in mapping for n in atom_names):
                    rings.append((key,ring_number,[mapping[n] for n in atom_names]))
        self.donor_pairs = np.asarray(donor_pairs,dtype=int).reshape(-1,2)
        self.acceptors = np.asarray(sorted(set(acceptors)),dtype=int)
        self.positive = np.asarray(positive,dtype=int)
        self.negative = np.asarray(negative,dtype=int)
        self.rings = rings
        self.expected_donor_sites = sorted(set(expected_donors))
        self.water_donor_pairs = np.asarray([(o,h) for o in self.water_oxygens for h in attached.get(int(o),[])],dtype=int).reshape(-1,2)
        self.previous, self.following = {}, {}
        for left,right in zip(self.residues,self.residues[1:]):
            a,b = self.atom_names[left],self.atom_names[right]
            if 'C' in a and 'N' in b and self.chains[a['C']]==self.chains[b['N']]:
                if np.linalg.norm(periodic_vectors((xyz[b['N']]-xyz[a['C']])[None,:],box)[0]) <= 1.9:
                    self.following[left] = right; self.previous[right] = left
        self._prepare_dssp(xyz)

    def _prepare_dssp(self, xyz):
        """Build an MDTraj backbone topology (N, CA, C, O of standard L-residues) for
        simplified DSSP; chains break wherever no peptide link was found.
        """
        self.dssp_status = 'unavailable_optional_mdtraj'
        self.md_topology = None
        try:
            import mdtraj as md
        except ImportError:
            return
        topology = md.Topology(); indices = []; residues = []; chain = None
        for key in self.residues:
            mapping = self.atom_names[key]
            if not all(n in mapping for n in ['N','CA','C','O']) or key in self.unknown_residues or str(self.resnames[self.groups[key][0]]) in D_PARENTS:
                continue
            if key not in self.previous or not residues or self.previous[key]!=residues[-1]:
                chain = topology.add_chain()
            name = str(self.resnames[self.groups[key][0]])
            residue = topology.add_residue(PARENTS.get(name,name),chain)
            for atom_name in ['N','CA','C','O']:
                element = md.element.nitrogen if atom_name=='N' else md.element.oxygen if atom_name=='O' else md.element.carbon
                topology.add_atom(atom_name,element,residue)
                indices.append(mapping[atom_name])
            residues.append(key)
        if residues:
            self.md_topology, self.md_indices, self.md_residues = topology,indices,residues
            self.dssp_status = 'mdtraj_simplified_DSSP_backbone_geometry'
        else:
            self.dssp_status = 'unavailable_complete_standard_backbones'

    def hydrogen_bonds(self, xyz, donors, acceptors, box):
        """Geometric hydrogen bonds between donor-H pairs and acceptors.

        Parameters
        ----------
        xyz : ndarray, shape (n, 3)
            Coordinates of all atoms (Å).
        donors : ndarray, shape (k, 2)
            ``(donor, hydrogen)`` index pairs.
        acceptors : ndarray of int
        box : array_like, optional

        Returns
        -------
        list of dict
            ``donor``, ``hydrogen``, ``acceptor``, ``distance_A`` (donor-acceptor
            distance) and ``angle_deg`` (D-H...A) for every pair passing the
            criteria in the module description.
        """
        import numpy as np
        from .hydration import distance_pairs
        if not len(donors) or not len(acceptors):
            return []
        dh = periodic_vectors(xyz[donors[:,1]]-xyz[donors[:,0]],box)
        maximum = np.where(self.elements[donors[:,0]]=='S',1.65,1.25)
        valid = (np.linalg.norm(dh,axis=1)>.4)&(np.linalg.norm(dh,axis=1)<=maximum)
        donors = donors[valid]
        pairs, distances = distance_pairs(xyz[donors[:,0]],xyz[acceptors],3.5,box)
        if not len(pairs):
            return []
        d,h = donors[pairs[:,0]].T; a = acceptors[pairs[:,1]]
        hd = periodic_vectors(xyz[d]-xyz[h],box); ha = periodic_vectors(xyz[a]-xyz[h],box)
        norms = np.linalg.norm(hd,axis=1)*np.linalg.norm(ha,axis=1)
        angle = np.degrees(np.arccos(np.clip(np.sum(hd*ha,axis=1)/np.maximum(norms,1e-12),-1,1)))
        valid = (distances>1.2)&(np.linalg.norm(ha,axis=1)<=2.5)&(angle>=150)&(self.residue_ids[d]!=self.residue_ids[a])
        return [{'donor':int(di),'hydrogen':int(hi),'acceptor':int(ai),'distance_A':float(dd),'angle_deg':float(aa)}
                for di,hi,ai,dd,aa in zip(d[valid],h[valid],a[valid],distances[valid],angle[valid])]

    def measure(self, coordinates, focus, *, box=None):
        """All interaction observations, torsions and secondary structure for a frame.

        Parameters
        ----------
        coordinates : array_like, shape (n, 3)
            Coordinates of all atoms (Å), in the topology's atom order.
        focus : iterable of str
            Residue labels to report. Observations are kept only when both
            residues are in the focus set.
        box : array_like, optional

        Returns
        -------
        dict
            ``observations`` (one dict per detected interaction with ``kind``,
            residue labels ``source`` < ``target``, atom indices, ``distance_A``
            and ``angle_deg``), ``water_bonds`` (residue-water hydrogen bonds),
            ``torsions`` (per residue: ``phi``, ``psi``, ``chi1``-``chi4`` in
            degrees, NaN when unavailable) and ``secondary`` (``"H"``, ``"E"``,
            ``"C"`` or ``"NA"``).
        """
        import numpy as np
        from .hydration import distance_pairs
        xyz = np.asarray(coordinates,dtype=float); focus = set(focus)
        observations = []
        def record(kind, left, right, atom_indices, distance, angle=None, **extra):
            if left==right or left not in focus or right not in focus:
                return
            observations.append({'kind':kind,'source':min(left,right),'target':max(left,right),
                                 'atom_indices':[int(i) for i in atom_indices], 'distance_A':float(distance),
                                 'angle_deg':float(angle) if angle is not None else None,**extra})
        for bond in self.hydrogen_bonds(xyz,self.donor_pairs,self.acceptors,box):
            d,h,a = (bond[k] for k in ['donor','hydrogen','acceptor'])
            record('hydrogen_bond',str(self.residue_ids[d]),str(self.residue_ids[a]),[d,h,a],bond['distance_A'],bond['angle_deg'])
        pairs, distances = distance_pairs(xyz[self.positive],xyz[self.negative],4.,box) if len(self.positive) and len(self.negative) else ([],[])
        for (i,j), distance in zip(pairs,distances):
            a,b = int(self.positive[i]),int(self.negative[j])
            record('salt_bridge_candidate',str(self.residue_ids[a]),str(self.residue_ids[b]),[a,b],distance)
        rings = []
        for key,number,indices in self.rings:
            if key not in focus:
                continue
            points = xyz[indices[0]]+periodic_vectors(xyz[indices]-xyz[indices[0]],box)
            center = points.mean(0); _,_,vh = np.linalg.svd(points-center,full_matrices=False)
            rms = float(np.sqrt(np.mean(((points-center)@vh[-1])**2)))
            if rms<=.15:
                rings.append((key,number,indices,center,vh[-1],rms))
        for left,right in combinations(rings,2):
            if left[0]==right[0]:
                continue
            delta = periodic_vectors((right[3]-left[3])[None,:],box)[0]; distance = np.linalg.norm(delta)
            if not 2.5<=distance<=5.5:
                continue
            angle = np.degrees(np.arccos(np.clip(abs(np.dot(left[4],right[4])),0,1)))
            offsets = [float(np.linalg.norm(delta-np.dot(delta,n)*n)) for n in [left[4],right[4]]]
            kind = 'aromatic_parallel_candidate' if angle<=30 and max(offsets)<=2 else 'aromatic_edge_face_candidate' if angle>=60 and min(offsets)<=2 else None
            if kind:
                record(kind,left[0],right[0],[left[2][0],right[2][0]],distance,angle,ring_offsets_A=offsets)
        # Each bridge requires two explicit-H, directional protein-water bonds.
        owners = defaultdict(set); water_bonds = []
        for bond in self.hydrogen_bonds(xyz,self.donor_pairs,self.water_oxygens,box):
            key = str(self.residue_ids[bond['donor']]); water = bond['acceptor']
            if key in focus:
                owners[water].add(key); water_bonds.append((key,water,bond))
        for bond in self.hydrogen_bonds(xyz,self.water_donor_pairs,self.acceptors,box):
            key = str(self.residue_ids[bond['acceptor']]); water = bond['donor']
            if key in focus:
                owners[water].add(key); water_bonds.append((key,water,bond))
        water_examples = {(water,key):bond for key,water,bond in water_bonds}
        for water,keys in owners.items():
            for left,right in combinations(sorted(keys),2):
                first,second = water_examples[water,left],water_examples[water,right]
                a = first['acceptor'] if first['donor']==water else first['donor']
                b = second['acceptor'] if second['donor']==water else second['donor']
                record('water_bridge_hbond',left,right,[a,water,b],max(first['distance_A'],second['distance_A']),
                       min(first['angle_deg'],second['angle_deg']),water_id=str(self.residue_ids[water]))
        torsions = {}
        for key in sorted(focus):
            mapping = self.atom_names[key]; name = str(self.resnames[self.groups[key][0]]); parent = PARENTS.get(name,name)
            row = {'phi':float('nan'),'psi':float('nan'), **{'chi'+str(i):float('nan') for i in range(1,5)}}
            for number,names in enumerate(CHI.get(parent,[]),1):
                if all(n in mapping for n in names):
                    row['chi'+str(number)] = torsion_degrees(xyz[[mapping[n] for n in names]],box)
            if key in self.previous and all(n in mapping for n in ['N','CA','C']):
                previous = self.atom_names[self.previous[key]]
                row['phi'] = torsion_degrees(xyz[[previous['C'],mapping['N'],mapping['CA'],mapping['C']]],box)
            if key in self.following and all(n in mapping for n in ['N','CA','C']):
                following = self.atom_names[self.following[key]]
                row['psi'] = torsion_degrees(xyz[[mapping['N'],mapping['CA'],mapping['C'],following['N']]],box)
            torsions[key] = row
        secondary = {key:'NA' for key in focus}
        if self.md_topology is not None:
            import mdtraj as md
            trajectory = md.Trajectory((xyz[self.md_indices]/10)[None,:,:],self.md_topology)
            labels = md.compute_dssp(trajectory,simplified=True)[0]
            secondary.update({key:str(label) for key,label in zip(self.md_residues,labels) if key in focus})
        return {'observations':observations,'water_bonds':water_bonds,'torsions':torsions,'secondary':secondary}

    def metadata(self):
        """Definitions, counts and exclusions for reports.

        Returns
        -------
        dict
            Template scope, unknown and ambiguous residues, hydrogen attachment
            method, typed site counts, the hydrogen-bond, salt-bridge, aromatic,
            torsion and secondary-structure definitions, and
            ``hydrogen_bond_status`` (``"unavailable_no_explicit_donor_hydrogens"``
            when no donor carries a hydrogen).
        """
        donor_sites = set(int(i) for i in self.donor_pairs[:,0])
        return {'chemistry':'standard amino-acid atom templates plus named protonation variants and CCD DLE/DVA; unsupported residues excluded',
                'unknown_residues':self.unknown_residues,
                'modified_residue_templates':[{'residue':key,'component':str(self.resnames[self.groups[key][0]]),'chemical_parent':D_PARENTS[str(self.resnames[self.groups[key][0]])],'chirality':'D; coordinates and signed torsions retained'} for key in self.residues if str(self.resnames[self.groups[key][0]]) in D_PARENTS],
                'D_residue_secondary_policy':'DLE/DVA omitted from unvalidated DSSP classification; chemical templates and raw torsions remain available',
                'ambiguous_residues_excluded':self.ambiguous_residues,
                'water_oxygen_population':len(self.water_oxygens),'hydrogen_attachment':self.attachment_method,
                'ambiguous_hydrogens_excluded':self.ambiguous_hydrogens,
                'typed_donor_sites':len(self.expected_donor_sites),'typed_donor_sites_with_explicit_H':len(donor_sites),
                'typed_acceptor_atoms':len(self.acceptors),'water_oxygens_with_explicit_H':len(set(self.water_donor_pairs[:,0].tolist())),
                'hydrogen_bond_definition':{'maximum_DA_A':3.5,'maximum_HA_A':2.5,'minimum_DHA_degrees':150},
                'hydrogen_bond_status':'observed_explicit_hydrogen_geometry' if len(donor_sites) or len(self.water_donor_pairs) else 'unavailable_no_explicit_donor_hydrogens',
                'missing_hydrogen_policy':'no hydrogen construction; absent/untyped donor observations are not evidence of absent bonds',
                'salt_bridge_definition':'oppositely charged named sidechain groups within 4 A; formal charge/protonation candidates, no energetic assignment',
                'aromatic_definition':'complete template rings, plane RMS <= 0.15 A, centroid distance 2.5–5.5 A; parallel <=30 degrees/max offset <=2 A; edge-face >=60 degrees/min offset <=2 A',
                'secondary_structure':self.dssp_status,'secondary_definition':'MDTraj simplified DSSP H/E/C with backbone-based virtual-H geometry; distinct from explicit-H interaction detection',
                'torsion_definition':'phi/psi and available chi1–chi4 in degrees; broken/missing quadruplets are NaN; coarse chi bins are not a residue-specific rotamer library',
                'interpretation':'geometric chemical evidence and descriptive occupancy, not interaction energies or functional control'}

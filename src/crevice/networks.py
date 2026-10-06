"""Residue contact networks, optionally with a channel/cavity pseudo-node.

Nodes are residues; an edge joins two residues whose closest atoms are within
a distance cutoff (4.5 Å by default). Two distance metrics are available:

``"surface"`` (default of the Python functions)
    Van der Waals surface gap: centre distance minus both atomic radii.
    Negative for overlapping spheres, such as covalently bonded neighbours.
``"center"`` (default of ``crevice network``)
    Plain atom-centre distance, the usual definition for residue contact
    networks.

Each edge is labelled with a coarse class from residue names and atom pairs
(``salt_bridge_candidate``, ``oppositely_charged_residue_contact``,
``hydrophobic_residue_contact``, ``polar_residue_contact`` or
``distance_contact``). :func:`build_cavity_network` adds a pseudo-node for a
profile, tunnel, cavity or cast component, linked to the residues that line
it.

Main entry points: :func:`build_residue_network`,
:func:`build_cavity_network`, :func:`network_metrics`,
:func:`compare_networks`. Command: ``crevice network``.

Limitations
    Edges are geometric proximity with chemical-class labels. They are not
    hydrogen bonds, energies or evidence of allosteric coupling. Sequence
    neighbours are not excluded, so the backbone chain always appears as
    edges. For explicit-hydrogen hydrogen bonds over trajectories, see
    :mod:`crevice.interaction_chemistry`.
"""

from __future__ import annotations

import itertools
import math
from typing import Mapping, Sequence

from .analysis import annotate_residues
from .geometry import distance
from .models import (
    Cavity,
    NetworkEdge,
    NetworkNode,
    PoreProfile,
    ResidueContact,
    ResidueInteractionNetwork,
    StructureFrame,
    Tunnel,
    VoidComponent,
)
from .radii import RadiusSet, atom_vdw_radius, radii_option, residue_properties


@radii_option
def build_residue_network(
    frame: StructureFrame,
    *,
    cutoff: float = 4.5,
    include_hydrogen: bool = False,
    include_hetero: bool = False,
    residue_groups: Mapping[str, str] | None = None,
    distance_metric: str = "surface",
    radii: RadiusSet | str | None = None,
) -> ResidueInteractionNetwork:
    """Build residue contacts. The legacy default uses van der Waals gaps.

    ``distance_metric="center"`` uses the minimum heavy-atom centre distance,
    matching the usual 4.5 Å proximity definition for residue contact networks.

    Parameters
    ----------
    frame : StructureFrame
    cutoff : float, default 4.5
        Maximum distance in Å (inclusive) under the chosen metric.
    include_hydrogen : bool, default False
        Include hydrogen atoms in the distance search.
    include_hetero : bool, default False
        Include HETATM residues (ligands, waters, ions) as nodes. Note that
        ``crevice network`` keeps heteroatoms, while the networks written by
        ``publish`` and ``analyze`` use this default.
    residue_groups : mapping of str to str, optional
        Residue label to group name (for example a helix name); stored in each
        node's ``metadata["group"]`` (``"unassigned"`` otherwise). Every key
        must be a residue in the selection.
    distance_metric : {"surface", "center"}, default "surface"
        See the module description. ``crevice network`` and ``crevice
        trajectory`` default to ``"center"``; ``publish`` and ``analyze`` use
        this default. Unifying the defaults is an open interface decision; the
        metric used is recorded in ``metadata["distance_metric"]``.
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    ResidueInteractionNetwork
        Nodes for every residue with a selected atom (in structure order) and
        one edge per residue pair, with ``distance`` the minimum over atom
        pairs (Å). Edge ``weight`` is ``1/(1 + max(gap, 0))`` for the surface
        metric and ``1/(1 + distance)`` for the centre metric.

    Raises
    ------
    ValueError
        For an unknown metric, a negative or non-finite cutoff, or invalid
        ``residue_groups``.
    ImportError
        For the centre metric if NumPy/SciPy cannot be imported.

    Notes
    -----
    Edge classes, in order of precedence:

    1. ``salt_bridge_candidate``: a LYS NZ or ARG NE/NH1/NH2 nitrogen within
       4.0 Å (centre distance) of an ASP OD1/OD2 or GLU OE1/OE2 oxygen.
    2. ``oppositely_charged_residue_contact``: one residue positive and the
       other negative by name (HIS counts as positive).
    3. ``hydrophobic_residue_contact``: both residues hydrophobic.
    4. ``polar_residue_contact``: surface gap of the closest pair at most
       0.5 Å and at least one residue polar.
    5. ``distance_contact``: anything else.

    With the centre metric, edges also record the observed sequence separation
    within a chain and a ``nonlocal_contact`` flag (different chains or
    separation greater than 2).

    Examples
    --------
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.networks import build_residue_network
    >>> frame = StructureFrame((
    ...     Atom(1, "NZ", "LYS", "A", 1, 0.0, 0.0, 0.0, "N"),
    ...     Atom(2, "OD1", "ASP", "A", 5, 3.0, 0.0, 0.0, "O"),
    ...     Atom(3, "CB", "ALA", "A", 9, 20.0, 0.0, 0.0, "C"),
    ... ))
    >>> network = build_residue_network(frame, distance_metric="center")
    >>> [(e.source, e.target, e.interaction, e.distance) for e in network.edges]
    [('A:ASP5', 'A:LYS1', 'salt_bridge_candidate', 3.0)]
    """
    if distance_metric not in {"surface","center"}:
        raise ValueError("distance_metric must be surface or center")

    if not math.isfinite(cutoff) or cutoff < 0:
        raise ValueError("cutoff must be non-negative")
    atoms = frame.selected_atoms(include_hydrogen=include_hydrogen, include_hetero=include_hetero)
    selected = set(atoms)
    residues = tuple(residue for residue in frame.residues() if any(atom in selected for atom in residue.atoms))
    residue_by_label = {residue.label: residue for residue in residues}
    groups = dict(residue_groups or {})
    if set(groups) - set(residue_by_label):
        raise ValueError("residue_groups contains residues outside the selected structure")
    if any(not isinstance(value, str) or not value.strip() for value in groups.values()):
        raise ValueError("Residue group names must be non-empty strings")
    nodes = tuple(
        NetworkNode(
            id=residue.label,
            kind="residue",
            label=residue.label,
            residue=residue.label,
            properties=residue_properties(residue.key.resname),
            metadata={"resname": residue.key.resname, "chain_id": residue.key.chain_id,
                      "resid": residue.key.resid, "icode": residue.key.icode,
                      "group": groups.get(residue.label, "unassigned"),
                      "position": tuple(sum(a.coord[i] for a in residue.atoms if a in selected)
                                        / sum(a in selected for a in residue.atoms) for i in range(3))},
        )
        for residue in residues
    )

    if distance_metric=="center":
        return _center_contact_network(atoms,nodes,cutoff,include_hydrogen,include_hetero,groups)

    radii = {index: atom_vdw_radius(atom) for index, atom in enumerate(atoms)}
    max_radius = max(radii.values(), default=1.7)
    center_cutoff = cutoff + 2.0 * max_radius
    cell_size = max(center_cutoff, 0.1)
    cells: dict[tuple[int, int, int], list[tuple[int, object]]] = {}
    for index, atom in enumerate(atoms):
        cells.setdefault(_cell_index(atom.coord, cell_size), []).append((index, atom))

    best_by_pair: dict[tuple[str, str], float] = {}
    evidence = {}
    ionic_pairs = {}
    for cell_idx in sorted(cells):
        for neighbor_idx in _neighbor_cell_indices(cell_idx):
            if neighbor_idx not in cells or neighbor_idx < cell_idx:
                continue
            left_items = cells[cell_idx]
            right_items = cells[neighbor_idx]
            if neighbor_idx == cell_idx:
                pairs = itertools.combinations(left_items, 2)
            else:
                pairs = itertools.product(left_items, right_items)
            for (left_index, left_atom), (right_index, right_atom) in pairs:
                left_label = left_atom.residue_key.label
                right_label = right_atom.residue_key.label
                if left_label == right_label:
                    continue
                center_distance = distance(left_atom.coord, right_atom.coord)
                max_pair_center = cutoff + radii[left_index] + radii[right_index]
                if center_distance > max_pair_center:
                    continue
                surface_distance = center_distance - radii[left_index] - radii[right_index]
                pair = tuple(sorted((left_label, right_label)))
                if _charged_pair(left_atom, right_atom) and center_distance <= 4.0:
                    ionic_pairs.setdefault(pair, []).append((left_atom.name, right_atom.name, center_distance))
                previous = best_by_pair.get(pair)
                if previous is None or surface_distance < previous:
                    best_by_pair[pair] = surface_distance
                    evidence[pair] = {"atoms": [left_atom.name, right_atom.name],
                                      "residues": [left_label, right_label],
                                      "center_distance": center_distance}

    edges: list[NetworkEdge] = []
    for (left_label, right_label), surface_distance in sorted(best_by_pair.items()):
        left = residue_by_label[left_label]
        right = residue_by_label[right_label]
        interaction = _interaction_type(left.key.resname, right.key.resname, surface_distance)
        if (left_label, right_label) in ionic_pairs:
            interaction = "salt_bridge_candidate"
        edges.append(
            NetworkEdge(
                source=left_label,
                target=right_label,
                interaction=interaction,
                distance=surface_distance,
                weight=1.0 / (1.0 + max(surface_distance, 0.0)),
                metadata={"closest_pair": evidence[(left_label, right_label)],
                          "charged_atom_pairs": ionic_pairs.get((left_label, right_label), []),
                          "distance_units": "angstrom"},
            )
        )
    return ResidueInteractionNetwork(
        nodes=nodes,
        edges=tuple(edges),
        metadata={
            "kind": "residue_interaction_network",
            "cutoff": cutoff,
            "include_hydrogen": include_hydrogen,
            "include_hetero": include_hetero,
            "method": "atom_cell_list",
            "interaction_rules": "contact_v2_charged_sidechain_atoms_4A",
            "distance_metric":"vdw_surface_gap",
            "weight_kind": "inverse_surface_distance",
            "group_source": "user_supplied" if groups else "unassigned",
        },
    )

@radii_option
def build_cavity_network(
    frame: StructureFrame,
    region: PoreProfile | Cavity | Tunnel | VoidComponent,
    *,
    contacts: Sequence[ResidueContact] | None = None,
    cutoff: float = 4.5,
    include_hydrogen: bool = False,
    include_hetero: bool = False,
    residue_groups: Mapping[str, str] | None = None,
    distance_metric: str = "surface",
    radii: RadiusSet | str | None = None,
) -> ResidueInteractionNetwork:
    """Build a residue network augmented with a channel/cavity/tunnel pseudo-node.

    The residue network is built with :func:`build_residue_network`; then a
    node of kind ``"region"`` is added and connected to every residue contact
    of the region.

    Parameters
    ----------
    frame : StructureFrame
    region : PoreProfile, Cavity, Tunnel or VoidComponent
        The region. Profiles and tunnels contribute their points (with their
        clearance spheres), and their bottleneck is identified; void components
        contribute their points; a :class:`~crevice.models.Cavity` contributes only its centre
        point.
    contacts : sequence of ResidueContact, optional
        Pre-computed region contacts. Computed with
        :func:`crevice.analysis.annotate_residues` (same ``cutoff`` and atom
        selection) when omitted.
    cutoff, include_hydrogen, include_hetero, residue_groups, distance_metric
        As for :func:`build_residue_network`; ``cutoff`` also limits the region
        contacts (Å, surface gap to the region).
    radii : RadiusSet, str or None, optional
        Atomic radius set: a preset name (``"default"``, ``"bondi"``,
        ``"hole"``, ``"charmm_like"``), the path of a JSON, CSV or HOLE ``.rad`` radius
        file, or a :class:`~crevice.radii.RadiusSet`. ``None`` (default) uses
        the set in effect, which is
        :data:`~crevice.radii.DEFAULT_RADII` (standard table plus CHARMM36 ion radii) unless a caller chose another
        (:func:`~crevice.radii.use_radii`, ``--radii``). Every atomic radius
        used by this call comes from that set, and results with a ``metadata``
        dict record it as ``metadata["radii"]``; see
        :doc:`/methods/atomic-radii`.

    Returns
    -------
    ResidueInteractionNetwork
        The region node first, then residue nodes. Region edges have
        ``interaction`` ``"region_<role>"`` (for example ``"region_lining"``),
        ``distance`` equal to the contact's ``min_distance`` and ``weight`` equal
        to its influence score. Contacts with residues that are not nodes (for
        example HETATM residues when ``include_hetero=False``) get no edge.

    Raises
    ------
    TypeError
        For an unsupported region type.
    ValueError
        As for :func:`build_residue_network`.

    Notes
    -----
    For a :class:`~crevice.models.Cavity`, contacts are distances to the centre point, not to
    the cavity surface, so only residues within ``cutoff`` of the centre are
    linked. Use the cast's :class:`~crevice.models.VoidComponent` for a surface-based network.
    """

    network = build_residue_network(
        frame,
        cutoff=cutoff,
        include_hydrogen=include_hydrogen,
        include_hetero=include_hetero,
        residue_groups=residue_groups,
        distance_metric=distance_metric,
    )
    region_id, region_label, points = _region_descriptor(region)
    if contacts is None:
        contacts = annotate_residues(
            frame,
            points,
            identify_bottleneck=isinstance(region, (PoreProfile, Tunnel)),
            cutoff=cutoff,
            include_hydrogen=include_hydrogen,
            include_hetero=include_hetero,
        )
    region_node = NetworkNode(
        id=region_id,
        kind="region",
        label=region_label,
        metadata={"contact_count": len(contacts)},
    )
    residue_ids = {node.id for node in network.nodes}
    region_edges = [
        NetworkEdge(
            source=region_id,
            target=contact.residue,
            interaction=f"region_{contact.role}",
            distance=contact.min_distance,
            weight=contact.influence_score,
            metadata={"point_count": contact.point_count, "properties": list(contact.properties)},
        )
        for contact in contacts
        if contact.residue in residue_ids
    ]
    return ResidueInteractionNetwork(
        nodes=(region_node,) + network.nodes,
        edges=tuple(region_edges) + network.edges,
        metadata={**network.metadata, "kind": "region_residue_network", "region": region_label},
    )


def network_metrics(network: ResidueInteractionNetwork) -> dict[str, object]:
    """Small graph summaries without requiring NetworkX.

    Parameters
    ----------
    network : ResidueInteractionNetwork

    Returns
    -------
    dict
        ``node_count``, ``edge_count``, ``density`` (edges divided by
        ``n(n-1)/2``; divided by 1 for fewer than two nodes), ``degree`` (every
        node) and ``top_degree`` (up to ten ``(node, degree)`` pairs, highest
        first, ties by node id). Region pseudo-nodes and their edges are
        included.
    """

    degrees = network.degree()
    node_count = len(network.nodes)
    possible_edges = node_count * (node_count - 1) / 2 if node_count > 1 else 1.0
    top_degree = sorted(degrees.items(), key=lambda item: (-item[1], item[0]))[:10]
    return {
        "node_count": node_count,
        "edge_count": len(network.edges),
        "density": len(network.edges) / possible_edges,
        "degree": degrees,
        "top_degree": top_degree,
    }


def compare_networks(networks: Sequence[ResidueInteractionNetwork]) -> list[dict[str, object]]:
    """Pairwise Jaccard similarity of node sets and edge sets.

    Edges are compared as ``(min(source, target), max(source, target),
    interaction)``, so the same pair with a different class counts as a
    different edge. Two empty sets have similarity 1.0.

    Parameters
    ----------
    networks : sequence of ResidueInteractionNetwork

    Returns
    -------
    list of dict
        One row per pair ``i < j`` with ``left``, ``right``,
        ``node_jaccard`` and ``edge_jaccard``.
    """

    comparisons: list[dict[str, object]] = []
    for left_index, right_index in itertools.combinations(range(len(networks)), 2):
        left = networks[left_index]
        right = networks[right_index]
        left_nodes = {node.id for node in left.nodes}
        right_nodes = {node.id for node in right.nodes}
        left_edges = _edge_set(left)
        right_edges = _edge_set(right)
        comparisons.append(
            {
                "left": left_index,
                "right": right_index,
                "node_jaccard": _jaccard(left_nodes, right_nodes),
                "edge_jaccard": _jaccard(left_edges, right_edges),
            }
        )
    return comparisons



def _cell_index(coord, cell_size: float) -> tuple[int, int, int]:
    """Integer cell of a coordinate in a cubic cell list."""
    return (
        math.floor(coord[0] / cell_size),
        math.floor(coord[1] / cell_size),
        math.floor(coord[2] / cell_size),
    )


def _neighbor_cell_indices(idx: tuple[int, int, int]):
    """The 27 cells around (and including) a cell."""
    for offset in itertools.product((-1, 0, 1), repeat=3):
        yield (idx[0] + offset[0], idx[1] + offset[1], idx[2] + offset[2])


def _min_residue_surface_distance(left_atoms, right_atoms) -> float:
    """Minimum van der Waals surface gap between two atom sets (Å)."""
    best = float("inf")
    for left in left_atoms:
        left_radius = atom_vdw_radius(left)
        for right in right_atoms:
            surface_distance = distance(left.coord, right.coord) - left_radius - atom_vdw_radius(right)
            if surface_distance < best:
                best = surface_distance
    return best


def _charged_pair(left, right) -> bool:
    """True for a LYS NZ / ARG NE, NH1, NH2 nitrogen paired with an ASP OD1, OD2 /
    GLU OE1, OE2 oxygen (names and elements must match).
    """
    positive = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}}
    negative = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
    def charged(atom, names, element):
        return atom.name.upper() in names.get(atom.resname.upper(), set()) and atom.element.upper() == element
    return ((charged(left, positive, "N") and charged(right, negative, "O"))
            or (charged(right, positive, "N") and charged(left, negative, "O")))


def _interaction_type(left_resname: str, right_resname: str, surface_distance: float) -> str:
    """Residue-name class of a contact; see :func:`build_residue_network`."""
    left_props = set(residue_properties(left_resname))
    right_props = set(residue_properties(right_resname))
    if {"positive", "negative"} <= (left_props | right_props) and left_props != right_props:
        return "oppositely_charged_residue_contact"
    if "hydrophobic" in left_props and "hydrophobic" in right_props:
        return "hydrophobic_residue_contact"
    if surface_distance <= 0.5 and ("polar" in left_props or "polar" in right_props):
        return "polar_residue_contact"
    return "distance_contact"


def _region_descriptor(region: PoreProfile | Cavity | Tunnel):
    """Node id, label and contact points for a region object."""
    if isinstance(region, PoreProfile):
        return "region:channel", "channel", region.points
    if isinstance(region, Cavity):
        return f"region:cavity:{region.id}", f"cavity:{region.id}", (region.center,)
    if isinstance(region, Tunnel):
        return f"region:tunnel:{region.id}", f"tunnel:{region.id}", region.points
    if isinstance(region, VoidComponent):
        return f"region:void:{region.id}", f"void:{region.id}", region.points
    raise TypeError("region must be a PoreProfile, Cavity, Tunnel, or VoidComponent")


def _edge_set(network: ResidueInteractionNetwork) -> set[tuple[str, str, str]]:
    """Edges as ``(min id, max id, interaction)`` tuples."""
    return {
        (min(edge.source, edge.target), max(edge.source, edge.target), edge.interaction)
        for edge in network.edges
    }


def _jaccard(left: set, right: set) -> float:
    """Jaccard index; 1.0 when both sets are empty."""
    if not left and not right:
        return 1.0
    return len(left & right) / len(left | right)


def _center_contact_network(atoms,nodes,cutoff,include_hydrogen,include_hetero,groups):
    """Exact cutoff pairs with spatial acceleration and atom-level evidence.

    Uses a KD-tree to find all atom pairs within ``cutoff`` (centre distance)
    and keeps the closest pair per residue pair.
    """
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise ImportError("Center-distance contact networks require NumPy and SciPy, which are required CREVICE dependencies but could not be imported; reinstall CREVICE (pip install --force-reinstall crevice)") from exc
    xyz=np.asarray([a.coord for a in atoms])
    pairs=cKDTree(xyz).query_pairs(cutoff,output_type='ndarray')
    best={};charged={}
    if len(pairs):
        ds=np.linalg.norm(xyz[pairs[:,0]]-xyz[pairs[:,1]],axis=1)
        for (i,j),d in zip(pairs,ds):
            a,b=atoms[int(i)],atoms[int(j)];left,right=a.residue_key.label,b.residue_key.label
            if left==right:continue
            key=tuple(sorted((left,right)))
            value=(float(d),int(i),int(j))
            if key not in best or value<best[key]:best[key]=value
            if d<=4 and _charged_pair(a,b):
                charged.setdefault(key,[]).append((a.name,b.name,float(d)))
    order={};chain_counts={}
    for node in nodes:
        chain=node.metadata['chain_id'];order[node.id]=chain_counts.get(chain,0);chain_counts[chain]=order[node.id]+1
    lookup={n.id:n for n in nodes};edges=[]
    for (left,right),(d,i,j) in sorted(best.items()):
        a,b=atoms[i],atoms[j];gap=d-atom_vdw_radius(a)-atom_vdw_radius(b)
        same_chain=lookup[left].metadata['chain_id']==lookup[right].metadata['chain_id']
        separation=abs(order[left]-order[right]) if same_chain else None
        edges.append(NetworkEdge(left,right,'salt_bridge_candidate' if (left,right) in charged else _interaction_type(a.resname,b.resname,gap),
                    d,1/(1+d),{'closest_pair':{'atoms':[a.name,b.name],'residues':[a.residue_key.label,b.residue_key.label],
                    'atom_serials':[a.serial,b.serial],'center_distance':d,'surface_distance':gap},
                    'charged_atom_pairs':charged.get((left,right),[]),'distance_units':'angstrom',
                    'sequence_separation_observed':separation,'different_chains':not same_chain,
                    'nonlocal_contact':not same_chain or separation>2,
                    'sequence_separation_definition':'index difference in observed residue order within chain'}))
    return ResidueInteractionNetwork(nodes,tuple(edges),{'kind':'residue_interaction_network','cutoff':cutoff,
            'distance_metric':'atom_center_distance','include_hydrogen':include_hydrogen,'include_hetero':include_hetero,
            'method':'exact_kdtree_atom_pairs','interaction_rules':'contact_v2_charged_sidechain_atoms_4A',
            'weight_kind':'inverse_center_distance','group_source':'user_supplied' if groups else 'unassigned',
            'interpretation':'Geometric proximity and chemical-class candidates; no hydrogen-bond or energy assignment'})

"""Residue contact networks, optionally with a channel/cavity pseudo-node.

Nodes are residues; an edge joins two residues whose closest atoms are within
a distance cutoff (4.5 Å by default). Two distance metrics are available:

``"center"`` (default)
    Minimum heavy-atom centre distance, the usual definition for residue
    contact networks (protein structure networks and contact maps use
    4-5 Å).
``"surface"``
    Van der Waals surface gap: centre distance minus both atomic radii.
    Negative for overlapping spheres, such as covalently bonded neighbours.
    Much more permissive at the same cutoff (4.5 Å of gap is about 8 Å
    between carbon centres).

Each edge is labelled with a coarse class from residue names and atom pairs
(``salt_bridge_candidate``, ``oppositely_charged_residue_contact``,
``hydrophobic_residue_contact``, ``polar_residue_contact`` or
``distance_contact``). :func:`build_cavity_network` adds a pseudo-node for a
profile, cavity or cast component, linked to the residues that line
it.

Main entry points: :func:`build_residue_network`,
:func:`build_cavity_network`, :func:`network_metrics`,
:func:`compare_networks`. Command: ``crevice residues --network``.

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

import numpy as np
from scipy.spatial import cKDTree

from .analysis import annotate_residues
from .interaction_chemistry import SALT_BRIDGE_CUTOFF, formal_charge
from .models import (
    Cavity,
    NetworkEdge,
    NetworkNode,
    PoreProfile,
    ResidueContact,
    ResidueInteractionNetwork,
    StructureFrame,
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
    distance_metric: str = "center",
    radii: RadiusSet | str | None = None,
) -> ResidueInteractionNetwork:
    """Build residue contacts from the closest atom pair of every residue pair.

    Parameters
    ----------
    frame : StructureFrame
    cutoff : float, default 4.5
        Maximum distance in Å (inclusive) under the chosen metric.
    include_hydrogen : bool, default False
        Include hydrogen atoms in the distance search.
    include_hetero : bool, default False
        Include HETATM residues (ligands, waters, ions) as nodes. Note that
        ``crevice residues --network`` keeps heteroatoms, while the networks written by
        ``report`` and ``report --analyze`` use this default.
    residue_groups : mapping of str to str, optional
        Residue label to group name (for example a helix name); stored in each
        node's ``metadata["group"]`` (``"unassigned"`` otherwise). Every key
        must be a residue in the selection.
    distance_metric : {"center", "surface"}, default "center"
        See the module description; every CREVICE command uses the default
        unless ``--contact-metric surface`` is given. The metric is recorded in
        ``metadata["distance_metric"]``.
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
        pairs (Å) under the metric. Edge ``weight`` is ``1/(1 + distance)`` for
        the centre metric and ``1/(1 + max(gap, 0))`` for the surface metric.

    Raises
    ------
    ValueError
        For an unknown metric, a negative or non-finite cutoff, or invalid
        ``residue_groups``.

    Notes
    -----
    Edge classes, in order of precedence:

    1. ``salt_bridge_candidate``: oppositely charged side-chain atoms of
       :data:`crevice.interaction_chemistry.CHARGED_SITES` (ARG NE/NH1/NH2,
       LYS NZ, HSP/HIP ND1/NE2 vs ASP OD1/OD2, GLU OE1/OE2, CYM SG) within
       4.0 Å centre distance, the definition of the typed interactions.
       Protonation comes from residue names (LYN, ASH, GLH are neutral).
    2. ``oppositely_charged_residue_contact``: one residue positive and the
       other negative by name (HIS counts as positive).
    3. ``hydrophobic_residue_contact``: both residues hydrophobic.
    4. ``polar_residue_contact``: surface gap of the closest pair at most
       0.5 Å and at least one residue polar.
    5. ``distance_contact``: anything else.

    Edges also record the closest atom pair (names, serials, centre distance
    and surface gap), the observed sequence separation within a chain and a
    ``nonlocal_contact`` flag (different chains or separation greater than 2).

    Examples
    --------
    >>> from crevice.models import Atom, StructureFrame
    >>> from crevice.networks import build_residue_network
    >>> frame = StructureFrame((
    ...     Atom(1, "NZ", "LYS", "A", 1, 0.0, 0.0, 0.0, "N"),
    ...     Atom(2, "OD1", "ASP", "A", 5, 3.0, 0.0, 0.0, "O"),
    ...     Atom(3, "CB", "ALA", "A", 9, 20.0, 0.0, 0.0, "C"),
    ... ))
    >>> network = build_residue_network(frame)
    >>> [(e.source, e.target, e.interaction, e.distance) for e in network.edges]
    [('A:ASP5', 'A:LYS1', 'salt_bridge_candidate', 3.0)]
    """
    if distance_metric not in {"surface", "center"}:
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
            metadata={
                "resname": residue.key.resname,
                "chain_id": residue.key.chain_id,
                "resid": residue.key.resid,
                "icode": residue.key.icode,
                "group": groups.get(residue.label, "unassigned"),
                "position": tuple(
                    sum(a.coord[i] for a in residue.atoms if a in selected) / sum(a in selected for a in residue.atoms)
                    for i in range(3)
                ),
            },
        )
        for residue in residues
    )
    return _contact_network(atoms, nodes, cutoff, distance_metric == "center", include_hydrogen, include_hetero, groups)


@radii_option
def build_cavity_network(
    frame: StructureFrame,
    region: PoreProfile | Cavity | VoidComponent,
    *,
    contacts: Sequence[ResidueContact] | None = None,
    cutoff: float = 4.5,
    include_hydrogen: bool = False,
    include_hetero: bool = False,
    residue_groups: Mapping[str, str] | None = None,
    distance_metric: str = "center",
    radii: RadiusSet | str | None = None,
) -> ResidueInteractionNetwork:
    """Build a residue network augmented with a channel/cavity pseudo-node.

    The residue network is built with :func:`build_residue_network`; then a
    node of kind ``"region"`` is added and connected to every residue contact
    of the region.

    Parameters
    ----------
    frame : StructureFrame
    region : PoreProfile, Cavity or VoidComponent
        The region. Profiles contribute their points (with their
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
            identify_bottleneck=isinstance(region, PoreProfile),
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


def _charge(atom) -> int:
    """Formal side-chain charge of an atom whose element matches its name (see :data:`CHARGED_SITES`)."""
    return formal_charge(atom.resname, atom.name) if atom.element.upper() == atom.name[:1].upper() else 0


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


def _region_descriptor(region: PoreProfile | Cavity | VoidComponent):
    """Node id, label and contact points for a region object."""
    if isinstance(region, PoreProfile):
        return "region:channel", "channel", region.points
    if isinstance(region, Cavity):
        return f"region:cavity:{region.id}", f"cavity:{region.id}", (region.center,)
    if isinstance(region, VoidComponent):
        return f"region:void:{region.id}", f"void:{region.id}", region.points
    raise TypeError("region must be a PoreProfile, Cavity or VoidComponent")


def closest_edges(edges) -> dict[tuple[str, str], NetworkEdge]:
    """The shortest edge of every node pair: ``{(min id, max id): edge}``.

    The pooling unit of contact occupancy (:func:`crevice.connectivity.persistent_network`,
    :class:`crevice.residue_dynamics.ResidueDynamics`): a pair counts once per
    frame whatever its interaction label.

    Parameters
    ----------
    edges : iterable of NetworkEdge
        Edges of one network (one frame).

    Returns
    -------
    dict
        ``{(min id, max id): edge}`` with the smallest ``distance`` per pair.
    """
    best = {}
    for edge in edges:
        pair = tuple(sorted((edge.source, edge.target)))
        if pair not in best or edge.distance < best[pair].distance:
            best[pair] = edge
    return best


def _edge_set(network: ResidueInteractionNetwork) -> set[tuple[str, str, str]]:
    """Edges as ``(min id, max id, interaction)`` tuples."""
    return {(min(edge.source, edge.target), max(edge.source, edge.target), edge.interaction) for edge in network.edges}


def _jaccard(left: set, right: set) -> float:
    """Jaccard index; 1.0 when both sets are empty."""
    if not left and not right:
        return 1.0
    return len(left & right) / len(left | right)


def _contact_network(atoms, nodes, cutoff, center, include_hydrogen, include_hetero, groups):
    """Closest atom pair of every residue pair within ``cutoff`` (KD-tree, exact).

    ``center`` selects the centre-distance metric; otherwise the surface gap is
    the edge distance and the search reaches ``cutoff`` plus the two largest radii.
    """
    xyz = np.asarray([a.coord for a in atoms], dtype=float).reshape(-1, 3)
    radii = np.asarray([atom_vdw_radius(a) for a in atoms], dtype=float)
    labels = [a.residue_key.label for a in atoms]
    codes = {label: k for k, label in enumerate(dict.fromkeys(labels))}
    residue = np.asarray([codes[label] for label in labels], dtype=np.int64)
    charge = np.asarray([_charge(a) for a in atoms], dtype=np.int64)
    reach = cutoff if center else cutoff + 2.0 * float(radii.max(initial=0.0))
    i, j = cKDTree(xyz).query_pairs(reach, output_type='ndarray').reshape(-1, 2).T
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    gap = d - radii[i] - radii[j]
    measure = d if center else gap
    keep = (residue[i] != residue[j]) & (measure <= cutoff)
    i, j, d, gap, measure = i[keep], j[keep], d[keep], gap[keep], measure[keep]
    pair = np.minimum(residue[i], residue[j]) * len(codes) + np.maximum(residue[i], residue[j])
    order = np.lexsort((j, i, measure, pair))
    _, first = np.unique(pair[order], return_index=True)

    def key(k):
        return tuple(sorted((labels[i[k]], labels[j[k]])))

    best = {key(k): k for k in order[first]}
    charged = {}
    salt = np.flatnonzero((d <= SALT_BRIDGE_CUTOFF) & (charge[i] * charge[j] < 0))
    for k in salt[np.lexsort((j[salt], i[salt]))]:
        charged.setdefault(key(k), []).append((atoms[i[k]].name, atoms[j[k]].name, float(d[k])))
    position = {}
    chain_counts = {}
    for node in nodes:
        chain = node.metadata['chain_id']
        position[node.id] = chain_counts.get(chain, 0)
        chain_counts[chain] = position[node.id] + 1
    lookup = {n.id: n for n in nodes}
    edges = []
    for (left, right), k in sorted(best.items()):
        a, b = atoms[i[k]], atoms[j[k]]
        distance, center_distance, surface_distance = float(measure[k]), float(d[k]), float(gap[k])
        same_chain = lookup[left].metadata['chain_id'] == lookup[right].metadata['chain_id']
        separation = abs(position[left] - position[right]) if same_chain else None
        edges.append(
            NetworkEdge(
                left,
                right,
                'salt_bridge_candidate'
                if (left, right) in charged
                else _interaction_type(a.resname, b.resname, surface_distance),
                distance,
                1 / (1 + max(distance, 0.0)),
                {
                    'closest_pair': {
                        'atoms': [a.name, b.name],
                        'residues': [a.residue_key.label, b.residue_key.label],
                        'atom_serials': [a.serial, b.serial],
                        'center_distance': center_distance,
                        'surface_distance': surface_distance,
                    },
                    'charged_atom_pairs': charged.get((left, right), []),
                    'distance_units': 'angstrom',
                    'sequence_separation_observed': separation,
                    'different_chains': not same_chain,
                    'nonlocal_contact': not same_chain or separation > 2,
                    'sequence_separation_definition': 'index difference in observed residue order within chain',
                },
            )
        )
    return ResidueInteractionNetwork(
        nodes,
        tuple(edges),
        {
            'kind': 'residue_interaction_network',
            'cutoff': cutoff,
            'distance_metric': 'atom_center_distance' if center else 'vdw_surface_gap',
            'include_hydrogen': include_hydrogen,
            'include_hetero': include_hetero,
            'method': 'exact_kdtree_atom_pairs',
            'interaction_rules': 'contact_v3_formal_charge_sites_4A',
            'weight_kind': 'inverse_center_distance' if center else 'inverse_surface_distance',
            'group_source': 'user_supplied' if groups else 'unassigned',
            'interpretation': 'Geometric proximity and chemical-class candidates; no hydrogen-bond or energy assignment',
        },
    )

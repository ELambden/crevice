"""One contact engine and one salt-bridge definition across CREVICE.

Software tests on synthetic atoms: both network metrics come from the same
KD-tree search, the network ``salt_bridge_candidate`` label uses the charged
sites of the typed interactions, and per-frame networks pool per residue pair
in both ``persistent_network`` and ``ResidueDynamics``.
"""
import itertools

import pytest

from crevice import Atom, StructureFrame, ResidueDynamics, build_residue_network, persistent_network
from crevice.interaction_chemistry import CHARGED_SITES, formal_charge
from crevice.models import NetworkEdge, NetworkNode, ResidueInteractionNetwork
from crevice.interaction_analysis import chemical_topology_from_frame

SITES = [(resname, name) for resname, (_, names) in CHARGED_SITES.items() for name in names[:1]]
NEUTRAL = [("LYN", "NZ"), ("ASH", "OD1"), ("GLH", "OE1"), ("HSD", "NE2"), ("HIS", "NE2"), ("CYS", "SG")]


def pair_frame(left, right, separation=3.0):
    (lres, lname), (rres, rname) = left, right
    return StructureFrame((Atom(1, lname, lres, "A", 1, 0.0, 0.0, 0.0, lname[0]),
                           Atom(2, rname, rres, "A", 2, separation, 0.0, 0.0, rname[0])))


def typed_salt_bridge(frame):
    model = chemical_topology_from_frame(frame, [])
    found = model.measure([a.coord for a in frame.atoms], ["A:%s1" % frame.atoms[0].resname,
                                                             "A:%s2" % frame.atoms[1].resname])
    return any(o["kind"] == "salt_bridge_candidate" for o in found["observations"])


@pytest.mark.parametrize("left,right", list(itertools.combinations(SITES + NEUTRAL, 2)))
def test_network_label_and_typed_interaction_agree(left, right):
    frame = pair_frame(left, right)
    expected = formal_charge(*left) * formal_charge(*right) < 0
    network = build_residue_network(frame)
    assert (network.edges[0].interaction == "salt_bridge_candidate") is expected
    assert typed_salt_bridge(frame) is expected


def test_protonation_variants_by_residue_name():
    assert formal_charge("HIP", "ND1") == formal_charge("HSP", "NE2") == 1
    assert formal_charge("CYM", "SG") == -1
    assert formal_charge("LYN", "NZ") == formal_charge("ASH", "OD1") == formal_charge("GLH", "OE2") == 0
    assert formal_charge("LYS", "CA") == 0


def test_salt_bridge_cutoff_is_inclusive_4_angstrom():
    assert build_residue_network(pair_frame(("LYS", "NZ"), ("GLU", "OE1"), 4.0)).edges[0].interaction \
        == "salt_bridge_candidate"
    assert build_residue_network(pair_frame(("LYS", "NZ"), ("GLU", "OE1"), 4.01)).edges[0].interaction \
        != "salt_bridge_candidate"


def test_centre_is_the_default_and_surface_is_wider():
    frame = StructureFrame(tuple(Atom(i, "CA", "ALA", "A", i, 3.9 * i, 0.0, 0.0, "C") for i in range(1, 5)))
    centre = build_residue_network(frame)
    surface = build_residue_network(frame, distance_metric="surface")
    assert centre.metadata["distance_metric"] == "atom_center_distance"
    assert surface.metadata["distance_metric"] == "vdw_surface_gap"
    assert [(e.source, e.target) for e in centre.edges] == [("A:ALA1", "A:ALA2"), ("A:ALA2", "A:ALA3"),
                                                            ("A:ALA3", "A:ALA4")]
    # Surface gap 3.9*k - 2*1.7: every pair up to 7.8 Å apart (gap 4.4 Å) is a contact.
    assert len(surface.edges) == 5
    for edge in centre.edges + surface.edges:
        pair = edge.metadata["closest_pair"]
        assert pair["surface_distance"] == pytest.approx(pair["center_distance"] - 3.4)
    assert all(e.distance == e.metadata["closest_pair"]["center_distance"] for e in centre.edges)
    assert all(e.distance == e.metadata["closest_pair"]["surface_distance"] for e in surface.edges)


def frames_with_changing_label():
    nodes = tuple(NetworkNode(n, "residue", n, metadata={"position": (float(i), 0.0, 0.0)})
                  for i, n in enumerate("ab"))
    meta = {"cutoff": 4.5, "distance_metric": "atom_center_distance"}
    return [ResidueInteractionNetwork(nodes, (NetworkEdge("a", "b", label, d),), meta)
            for label, d in (("salt_bridge_candidate", 3.0), ("polar_residue_contact", 3.5),
                             ("salt_bridge_candidate", 3.1))]


def test_persistent_network_pools_per_pair_like_residue_dynamics():
    frames = frames_with_changing_label()
    pooled = persistent_network(frames, min_occupancy=0.75)
    (edge,) = pooled.edges
    assert edge.weight == 1.0 and edge.interaction == "salt_bridge_candidate"
    assert edge.metadata["contact_type_fractions"] == {"polar_residue_contact": pytest.approx(1 / 3),
                                                       "salt_bridge_candidate": pytest.approx(2 / 3)}
    dynamics = ResidueDynamics(min_occupancy=0.75)
    for k, network in enumerate(frames):
        dynamics.add(None, network, frame_index=k, time=float(k))
    (row,) = dynamics.finish(aligned=False)["contacts"]
    assert row["occupancy"] == edge.weight
    assert row["mean_distance_when_present_A"] == pytest.approx(edge.distance)
    assert row["contact_type_fractions"] == pytest.approx(edge.metadata["contact_type_fractions"])

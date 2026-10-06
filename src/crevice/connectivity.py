"""Persistent residue contacts and paths linking lining residues to structural groups."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import replace
from typing import Sequence

from .models import NetworkEdge, ResidueInteractionNetwork
from .presentation import FigureText, annotate_option


def persistent_network(networks: Sequence[ResidueInteractionNetwork], *,
                       min_occupancy: float = 0.5) -> ResidueInteractionNetwork:
    """Pool contacts conditional on both nodes being observed in a frame.

    Weights are empirical occupancy. Distances average only contact-present
    frames. Correlated frames do not give independent uncertainty estimates.

    A typed contact (pair plus interaction class) has occupancy = frames with the
    contact / frames in which both residues are present; contacts below
    ``min_occupancy`` are dropped. Each pooled node records ``observed_frames``
    and ``lining_occupancy`` (fraction of its observed frames with a region-lining
    contact).

    Parameters
    ----------
    networks : sequence of ResidueInteractionNetwork
        One per-frame network per frame, built with identical settings (cutoff,
        metric, hydrogen/hetero selection, interaction rules).
    min_occupancy : float, default 0.5
        Smallest kept occupancy, in [0, 1].

    Returns
    -------
    ResidueInteractionNetwork
        Pooled network: edge ``weight`` = occupancy, ``distance`` = mean
        distance (Å) over contact frames, edge metadata with the frame counts;
        network metadata ``weight_kind = "contact_occupancy"``.

    Raises
    ------
    ValueError
        For no networks, an invalid ``min_occupancy``, inconsistent settings or
        node identities, duplicate contacts, or already pooled input.
    """
    if not networks:
        raise ValueError("At least one network is required")
    if not math.isfinite(min_occupancy) or not 0 <= min_occupancy <= 1:
        raise ValueError("min_occupancy must be in [0, 1]")
    for key in ("cutoff", "include_hydrogen", "include_hetero", "interaction_rules", "distance_metric"):
        if any(n.metadata.get(key) != networks[0].metadata.get(key) for n in networks):
            raise ValueError(f"Incompatible per-frame network setting: {key}")
    if any(n.metadata.get("weight_kind") == "contact_occupancy" for n in networks):
        raise ValueError("Provide per-frame networks, not already pooled networks")
    nodes, observed, contacts = {}, defaultdict(set), defaultdict(dict)
    lining = defaultdict(set)
    for frame_index, network in enumerate(networks):
        ids = {n.id for n in network.nodes}
        if len(ids) != len(network.nodes):
            raise ValueError("Duplicate node identities")
        for node in network.nodes:
            if node.id in nodes and (node.kind != nodes[node.id].kind or
                                    node.metadata.get("group") != nodes[node.id].metadata.get("group")):
                raise ValueError("Node kind and structural group must remain consistent")
            nodes.setdefault(node.id, node)
            observed[node.id].add(frame_index)
        for edge in network.edges:
            if edge.source not in ids or edge.target not in ids:
                raise ValueError("An edge references an absent node")
            key = (*sorted((edge.source, edge.target)), edge.interaction)
            if frame_index in contacts[key]:
                raise ValueError("Duplicate typed contact in a frame")
            contacts[key][frame_index] = edge
            if edge.interaction in {"region_lining", "region_bottleneck"}:
                for node_id in (edge.source, edge.target):
                    if nodes[node_id].kind == "residue":
                        lining[node_id].add(frame_index)
    edges = []
    for (left, right, interaction), frames in sorted(contacts.items()):
        eligible = len(observed[left] & observed[right])
        occupancy = len(frames) / eligible
        if occupancy < min_occupancy:
            continue
        edges.append(NetworkEdge(left, right, interaction,
                     sum(e.distance for e in frames.values()) / len(frames), occupancy,
                     {"occupancy": occupancy, "contact_frames": len(frames),
                      "eligible_frames": eligible, "total_frames": len(networks)}))
    pooled_nodes = tuple(replace(node, metadata={**node.metadata,
                         "observed_frames": len(observed[node.id]),
                         "lining_occupancy": len(lining[node.id]) / len(observed[node.id])})
                         for _, node in sorted(nodes.items()))
    return ResidueInteractionNetwork(pooled_nodes, tuple(edges),
               {**networks[0].metadata, "weight_kind": "contact_occupancy",
                "frame_count": len(networks), "min_occupancy": min_occupancy,
                "missing_nodes": "pairwise_observed_denominator",
                "uncertainty": "not_estimated", "position_source": "first_observed_frame"})


def lining_connectivity(network: ResidueInteractionNetwork, *,
                        lining_residues: Sequence[str] | None = None,
                        removed_residues: Sequence[str] = ()) -> dict:
    """Find weighted graph paths from pore-lining residues using NetworkX.

    Region pseudo-nodes are excluded from traversal. Each edge cost is 1/weight;
    stronger contacts are cheaper. Parallel interaction types use the strongest
    weight without adding their probabilities. Paths are associations, not
    transport paths or evidence of physical signal transmission.

    Parameters
    ----------
    network : ResidueInteractionNetwork
        Static or pooled (:func:`persistent_network`) residue network.
    lining_residues : sequence of str, optional
        Source residues. Default: residues with a ``region_lining`` or
        ``region_bottleneck`` edge, plus pooled residues whose
        ``lining_occupancy`` reaches the network's ``min_occupancy``.
    removed_residues : sequence of str, default ()
        Residues removed from the graph before the search (a graph sensitivity
        experiment, as in ``network --remove-residue``).

    Returns
    -------
    dict
        ``residues`` (per residue: group, ``is_lining``, ``removed``,
        ``reachable``, ``cost_to_lining``, ``hops_to_lining``, ``path``,
        ``path_hops``, ``lining_occupancy``, ``bridge_endpoint``,
        ``articulation_point``), ``lining_residues``, ``removed_residues``,
        ``group_contacts`` (contact count and total strength between structural
        groups) and definitions. Written as CSV by
        :func:`crevice.io.write_connectivity_csv`.

    Raises
    ------
    ImportError
        Without NetworkX (``pip install 'crevice[networks]'``).
    ValueError
        For unknown residue IDs or negative/non-finite weights.
    """
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError("Connectivity analysis requires NetworkX: pip install 'crevice[networks]'") from exc
    nodes = {n.id: n for n in network.nodes if n.kind == "residue"}
    removed = set(removed_residues)
    if removed - set(nodes):
        raise ValueError("removed_residues contains unknown residue IDs")
    if lining_residues is None:
        sources = {node_id for e in network.edges
                   if e.interaction in {"region_lining", "region_bottleneck"}
                   for node_id in (e.source, e.target) if node_id in nodes}
        sources.update(n.id for n in nodes.values()
                       if n.metadata.get("lining_occupancy", 0) > 0
                       and n.metadata["lining_occupancy"] >= network.metadata.get("min_occupancy", 0.5))
    else:
        sources = set(lining_residues)
        if sources - set(nodes):
            raise ValueError("lining_residues contains unknown residue IDs")
    graph = nx.Graph()
    graph.add_nodes_from(sorted(set(nodes) - removed))
    for edge in network.edges:
        if edge.source not in graph or edge.target not in graph or edge.source == edge.target:
            continue
        if not math.isfinite(edge.weight) or edge.weight < 0:
            raise ValueError("Path weights must be finite and non-negative")
        if edge.weight == 0:
            continue
        previous = graph.get_edge_data(edge.source, edge.target)
        if previous is None or edge.weight > previous["strength"]:
            graph.add_edge(edge.source, edge.target, strength=edge.weight,
                           cost=1.0 / edge.weight, interaction=edge.interaction)
    active_sources = sorted(sources - removed)
    costs, paths = nx.multi_source_dijkstra(graph, active_sources, weight="cost") if active_sources else ({}, {})
    hops = nx.multi_source_dijkstra_path_length(graph, active_sources, weight=lambda u, v, d: 1) if active_sources else {}
    bridges = {n for pair in nx.bridges(graph) for n in pair}
    articulations = set(nx.articulation_points(graph))
    rows = []
    for node_id, node in sorted(nodes.items()):
        path = paths.get(node_id, [])
        rows.append({"residue": node_id, "group": node.metadata.get("group", "unassigned"),
                     "is_lining": node_id in sources, "removed": node_id in removed,
                     "reachable": node_id in costs, "cost_to_lining": costs.get(node_id),
                     "hops_to_lining": hops.get(node_id), "path": path,
                     "path_hops": len(path) - 1 if path else None,
                     "lining_occupancy": node.metadata.get("lining_occupancy"),
                     "bridge_endpoint": node_id in bridges,
                     "articulation_point": node_id in articulations})
    group_edges = defaultdict(lambda: {"contact_count": 0, "total_strength": 0.0})
    for left, right, edge in graph.edges(data=True):
        a, b = (nodes[n].metadata.get("group", "unassigned") for n in (left, right))
        if a != b:
            entry = group_edges[tuple(sorted((a, b)))]
            entry["contact_count"] += 1
            entry["total_strength"] += edge["strength"]
    return {"residues": rows, "lining_residues": sorted(sources),
            "removed_residues": sorted(removed),
            "group_contacts": [{"source": a, "target": b, **v} for (a, b), v in sorted(group_edges.items())],
            "weight_transform": "1/positive_weight", "parallel_types": "maximum_strength",
            "region_nodes_used_for_paths": False,
            "interpretation": "Contact-network association; removal is a graph sensitivity experiment"}


@annotate_option
def plot_lining_connectivity(network: ResidueInteractionNetwork, path, *,
                             lining_residues=None, max_hops: int = 3, max_nodes: int = 50,
                             dpi: int = 400, annotate: bool | None = None):
    """Draw residues in columns by contact-graph distance from the pore lining.

    :func:`lining_connectivity` computes, for every residue, the minimum number
    of residue-residue contact edges to a lining residue (region nodes are not
    used as shortcuts). Residues within ``max_hops`` are placed in columns
    (lining = column 0), sorted by hop count, group and ID, up to ``max_nodes``;
    grey lines join shown residues that are in contact. Residue names are the
    chart's category labels and are always drawn next to their markers.

    Parameters
    ----------
    network : ResidueInteractionNetwork
        Contact network with ``residue_groups`` (helix/strand names).
    path : str or Path
        Output image; parent directories are created.
    lining_residues : iterable of str, optional
        Lining residue IDs; default the residues contacting the region node.
    max_hops : int, default 3
        Largest graph distance shown (columns 0..max_hops).
    max_nodes : int, default 50
        Maximum residues shown.
    dpi : int, default 400
        Raster resolution in dots per inch.
    annotate : bool, optional
        Also draw the title "Pore-lining connectivity (N residues shown)".
        ``None`` (default) inherits the surrounding setting.

    Returns
    -------
    Path
        The written image.

    Raises
    ------
    ValueError
        If ``max_hops`` < 0 or ``max_nodes`` < 1.

    Outputs
    -------
    path : image
        Width 9 inch, height scaled with the largest column. x axis "Minimum
        number of residue contacts from the lining" (ticks "Lining", 1, 2, ...);
        no y axis. With no connected residues the panel says so.

    Colours and representation
    --------------------------
    Markers coloured by residue group (Matplotlib ``tab20`` cycle): squares for
    lining residues, circles for others; residue names to the right of each
    marker; grey (``#b5b5b5``) contact lines at 60% opacity. A legend names the
    groups when there are at least two (or when annotating).
    """
    from pathlib import Path
    from matplotlib.lines import Line2D
    from .figures import _pyplot, _save_figure

    if max_hops < 0 or max_nodes < 1:
        raise ValueError("max_hops must be >= 0 and max_nodes >= 1")
    report = lining_connectivity(network, lining_residues=lining_residues)
    rows = sorted((r for r in report["residues"] if r["reachable"] and r["hops_to_lining"] <= max_hops),
                  key=lambda r: (r["hops_to_lining"], r["group"], r["residue"]))[:max_nodes]
    text = FigureText(f"{len(rows)} residues within {max_hops} contact hops of the pore lining, in columns by graph distance; "
                      "marker colour is the residue group, squares are lining residues")
    plt = _pyplot()
    columns = defaultdict(list)
    for row in rows:
        columns[row["hops_to_lining"]].append(row)
    height = max(4.0, max((len(c) for c in columns.values()), default=0) * 0.32 + 1.8)
    fig, ax = plt.subplots(figsize=(9, height), constrained_layout=True)
    positions = {r["residue"]: (hop, i - (len(column) - 1) / 2)
                 for hop, column in columns.items() for i, r in enumerate(column)}
    groups = sorted({r["group"] for r in rows})
    colors = {group: plt.get_cmap("tab20")(i % 20) for i, group in enumerate(groups)}
    shown_edges = set()
    for edge in network.edges:
        pair = tuple(sorted((edge.source, edge.target)))
        if pair in shown_edges or not all(n in positions for n in pair):
            continue
        shown_edges.add(pair)
        x, y = zip(*(positions[n] for n in pair))
        ax.plot(x, y, color="#b5b5b5", lw=0.6, alpha=0.6, zorder=1)
    for row in rows:
        x, y = positions[row["residue"]]
        ax.scatter(x, y, s=55, color=colors[row["group"]], marker="s" if row["is_lining"] else "o", zorder=3)
        # Residue names are this chart's category labels.
        ax.annotate(row["residue"], (x, y), xytext=(7, 0), textcoords="offset points", va="center", fontsize=7)
    if not rows:
        text.notice(ax, 0.5, 0.5, "No residues connected to the selected lining", transform=ax.transAxes, ha="center")
    ax.set_xticks(range(max_hops + 1), ["Lining"] + [str(i) for i in range(1, max_hops + 1)])
    ax.set_xlabel("Minimum number of residue contacts from the lining")
    ax.set_yticks([])
    ax.set_xlim(-0.2, max_hops + 0.65)
    text.title(ax, f"Pore-lining connectivity ({len(rows)} residues shown)")
    for spine in ax.spines.values():
        spine.set_visible(False)
    if len(groups) > 1 or (groups and text.draw):
        ax.legend(handles=[Line2D([], [], marker="o", linestyle="", color=colors[g], label=g) for g in groups],
                  frameon=False, loc="upper left", bbox_to_anchor=(1, 1), fontsize=8)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    _save_figure(fig, output, dpi=dpi, text=text)
    return output

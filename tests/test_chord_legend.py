"""The chord-diagram legend sits below the whole circle and never overlaps a node name."""
import itertools

import pytest

import crevice.figures as figures
from crevice.figures import INTERACTION_LEGEND, plot_network_chord
from crevice.models import NetworkEdge, NetworkNode, ResidueInteractionNetwork


def network(count, name_length, with_region=True):
    kinds = [key for key, _ in INTERACTION_LEGEND if key != "region"]
    nodes = [NetworkNode(id=f"A:{'X' * name_length}{i}", kind="residue", label=f"A:{'X' * name_length}{i}")
             for i in range(count)]
    edges = [NetworkEdge(source=a.id, target=b.id, interaction=kinds[(i + j) % len(kinds)], distance=4.0,
                         weight=1.0 + (i * j) % 3)
             for (i, a), (j, b) in itertools.combinations(enumerate(nodes), 2) if (i + j) % 3 == 0]
    if with_region:
        region = NetworkNode(id="region:channel", kind="region", label="region:channel")
        edges += [NetworkEdge(source=region.id, target=n.id, interaction="channel_contact", distance=3.0)
                  for n in nodes[:5]]
        nodes = [region, *nodes]
    return ResidueInteractionNetwork(nodes=tuple(nodes), edges=tuple(edges))


@pytest.mark.parametrize("count,name_length,top_n", [(8, 3, 24), (40, 3, 24), (80, 12, 60)])
def test_legend_is_below_every_label(tmp_path, monkeypatch, count, name_length, top_n):
    captured = {}
    original = figures._save_figure

    def capture(fig, path, **kwargs):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        ax = fig.axes[0]
        captured["labels"] = [t.get_window_extent(renderer) for t in ax.texts if t.get_text()]
        legends = list(fig.legends) + ([ax.get_legend()] if ax.get_legend() is not None else [])
        captured["legends"] = [legend.get_window_extent(renderer) for legend in legends]
        original(fig, path, **kwargs)

    monkeypatch.setattr(figures, "_save_figure", capture)
    out = plot_network_chord(network(count, name_length), tmp_path / "chord.png", top_n=top_n, dpi=60)
    assert out.is_file()
    assert len(captured["legends"]) == 1
    legend = captured["legends"][0]
    assert len(captured["labels"]) == min(top_n, count + 1)
    lowest_label = min(box.y0 for box in captured["labels"])
    assert legend.y1 <= lowest_label + 1e-6, "legend overlaps or rises above the lowest node name"
    assert not any(legend.overlaps(box) for box in captured["labels"])

"""Length units in figures are written as Å, never as an ASCII "A".

Software test on a synthetic channel: every axis label, in-plot text and PNG
Description of the static publication bundle (with annotations drawn) is
collected and searched for ASCII ångström units such as "(A)", "(A^3)" or
"1.33 A".
"""
# Tests here exercise other behaviour with the legacy fixed 0.8 A enclosure probe;
# the automatic probe is tested in test_auto_enclosure.py.
import re
from pathlib import Path

import pytest

import crevice.figures as figures
from crevice.presentation import figure_description, write_viewer_structure
from test_channel_coordinate import rectangular_channel

ASCII_UNIT = re.compile(r"\((?:A|A\^?[23]|A-?[123])\)|\bA\^\d|\d A\b|\d A[,;)]")


def _labels(fig):
    texts = []
    for ax in fig.axes:
        texts += [ax.get_xlabel(), ax.get_ylabel(), ax.get_title()]
        texts += [t.get_text() for t in ax.texts]
    texts += [t.get_text() for t in fig.texts]
    return [t for t in texts if t]


@pytest.fixture
def labels(monkeypatch):
    seen = {}
    original = figures._save_figure

    def spy(fig, path, *, dpi, text=None):
        seen[Path(path).name] = _labels(fig)
        return original(fig, path, dpi=dpi, text=text)
    monkeypatch.setattr(figures, "_save_figure", spy)
    return seen


def test_publication_figures_use_the_angstrom_sign(tmp_path, labels):
    frame = rectangular_channel()
    source = write_viewer_structure(frame, tmp_path / "channel.pdb")
    root = tmp_path / "bundle"
    figures.write_static_publication_bundle(frame, enclosure_radius=0.8, structure_path=source, output_dir=root, prefix="c",
                                            samples=41, search_radius=8, cast_spacing=.5, dpi=60,
                                            annotate=True)
    assert len(labels) >= 6
    offending = {name: [t for t in texts if ASCII_UNIT.search(t)] for name, texts in labels.items()}
    assert not any(offending.values()), offending
    for png in root.glob("*.png"):
        description = figure_description(png).get("Description", "")
        assert not ASCII_UNIT.search(description), (png.name, description)
    profile = labels["c_profile_radius.png"]
    assert "Axial coordinate (Å)" in profile and "Pore radius (Å)" in profile
    assert any(t.endswith(" Å") for t in profile)
    assert "Volume (Å³)" in labels["c_cavity_summary.png"] or "No cavities detected" in " ".join(labels["c_cavity_summary.png"])

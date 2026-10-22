"""Public plot functions return the matplotlib Figure and still save when given a path.

Software tests on a synthetic ring channel. Saved images are unchanged by this
API change (byte-identical PNGs before/after on 1GRM and a synthetic
trajectory, checked when the change was made).
"""
import math

import numpy as np
import pytest

from crevice import (analyze_trajectory, annotate_residues, build_cavity_network, detect_cavities,
                     plot_cavity_summary, plot_lining_connectivity, plot_network_chord, plot_network_summary,
                     plot_profile_radius, plot_profile_radius_with_residues, plot_residue_contacts,
                     plot_trajectory_profiles, pore_profile, profile_distribution, trajectory_from_frames)
from crevice.models import Atom, StructureFrame

matplotlib = pytest.importorskip("matplotlib")
from matplotlib.figure import Figure


def _ring(noise=None):
    atoms = []
    for ring, z in enumerate(np.arange(-9.0, 9.1, 1.5)):
        for k in range(10):
            a = 2 * math.pi * k / 10
            dx = noise[len(atoms)] if noise is not None else (0, 0, 0)
            atoms.append(Atom(len(atoms) + 1, "CA", ["ALA", "SER", "LEU"][k % 3], "A", len(atoms) + 1,
                              4.2 * math.cos(a) + dx[0], 4.2 * math.sin(a) + dx[1], float(z) + dx[2], "C"))
    return StructureFrame(tuple(atoms))


@pytest.fixture(scope="module")
def static_inputs():
    frame = _ring()
    profile = pore_profile(frame, axis="z", samples=31)
    contacts = annotate_residues(frame, profile.points)
    network = build_cavity_network(frame, profile, contacts=contacts)
    return frame, profile, contacts, network


def _calls(static_inputs):
    frame, profile, contacts, network = static_inputs
    return {
        "profile_radius": lambda p: plot_profile_radius(profile, p, dpi=50),
        "profile_radius_with_residues": lambda p: plot_profile_radius_with_residues(profile, p, dpi=50),
        "residue_contacts": lambda p: plot_residue_contacts(contacts, p, dpi=50),
        "cavity_summary": lambda p: plot_cavity_summary(detect_cavities(frame), p, dpi=50),
        "network_summary": lambda p: plot_network_summary(network, p, dpi=50),
        "network_chord": lambda p: plot_network_chord(network, p, dpi=50),
        "lining_connectivity": lambda p: plot_lining_connectivity(network, p, dpi=50),
    }


@pytest.mark.parametrize("name", ["profile_radius", "profile_radius_with_residues", "residue_contacts",
                                  "cavity_summary", "network_summary", "network_chord", "lining_connectivity"])
def test_plot_returns_figure_and_saves(static_inputs, tmp_path, name):
    call = _calls(static_inputs)[name]
    target = tmp_path / "nested" / f"{name}.png"
    fig = call(target)
    assert isinstance(fig, Figure) and fig.axes, "the returned figure keeps its drawing"
    assert target.is_file() and target.stat().st_size > 0
    assert "Description" in fig.crevice_description
    again = tmp_path / f"{name}_again.png"
    fig.savefig(again)
    assert again.stat().st_size > 0
    unsaved = call(None)
    assert isinstance(unsaved, Figure) and unsaved.axes
    assert sorted(p.name for p in tmp_path.rglob("*.png")) == sorted([f"{name}.png", f"{name}_again.png"])


def test_plot_trajectory_profiles_returns_figure(tmp_path):
    rng = np.random.default_rng(1)
    frames = [StructureFrame(_ring(rng.normal(scale=0.1, size=(130, 3))).atoms, frame_index=i) for i in range(4)]
    result = analyze_trajectory(trajectory_from_frames(frames, time_step=10.0), analyses=("profile",),
                                profile_kwargs={"axis": "z", "samples": 21})
    fig = plot_trajectory_profiles(result, tmp_path / "traj.png", dpi=50, distribution=profile_distribution(result))
    assert isinstance(fig, Figure) and (tmp_path / "traj.png").is_file()
    assert isinstance(plot_trajectory_profiles(result, dpi=50), Figure)


def test_module_plot_functions_return_figures(tmp_path):
    from crevice.paths import find_paths, plot_paths
    frame = _ring()
    result = find_paths(frame, seed=(0.0, 0.0, 0.0), spacing=0.7)
    fig = plot_paths(result, tmp_path / "paths.png")
    assert isinstance(fig, Figure) and (tmp_path / "paths.png").is_file()
    assert isinstance(plot_paths(result), Figure)

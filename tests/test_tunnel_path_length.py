"""Tunnel points carry their cumulative path length in Å as ``t``.

Synthetic cylinder only (software test). Previously ``t`` was the step number,
unlike every other profile; the grid engine ``crevice.analysis.find_tunnels``
kept the step number until 30 Sep 2026. The step number stays available as
``ChannelPoint.index``; nothing else about a tunnel may change.
"""
import json
import math
from dataclasses import replace

import pytest

from crevice.analysis import find_tunnels as grid_engine
from crevice.cli import main
from crevice.ensemble import profile_distribution
from crevice.models import Atom, FrameAnalysis, StructureFrame
from crevice.presentation import write_viewer_structure
from crevice.tunnels import find_tunnels, tunnel_profile


def cylinder():
    atoms = [Atom(12 * k + j + 1, "C", "ALA", "A", k + 1,
                  6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6), -10.5 + 1.5 * k, "C")
             for k in range(15) for j in range(12)]
    return StructureFrame(tuple(atoms))


OPTIONS = dict(start=(0.3, 0.2, 0.1), spacing=2.0, min_radius=0.8, max_tunnels=3)


def test_t_is_cumulative_path_length_from_the_start():
    tunnels = find_tunnels(cylinder(), **OPTIONS)
    assert tunnels
    for tunnel in tunnels:
        points = tunnel.points
        assert points[0].t == 0.0
        assert [p.index for p in points] == list(range(len(points)))
        running = 0.0
        for a, b in zip(points, points[1:]):
            running += math.dist(a.position, b.position)
            assert b.t == pytest.approx(running, abs=1e-12)
        assert points[-1].t == pytest.approx(tunnel.length, abs=1e-9)


def test_grid_engine_returns_the_same_path_length_t():
    # crevice.analysis.find_tunnels (the grid engine) now also reports t in
    # Å, so both entry points return identical tunnels.
    frame = cylinder()
    raw = grid_engine(frame, **OPTIONS)
    converted = find_tunnels(frame, **OPTIONS)
    assert raw and raw == converted
    for tunnel in raw:
        assert tunnel.points[-1].t == pytest.approx(tunnel.length, abs=1e-9)
        assert [p.index for p in tunnel.points] == list(range(len(tunnel.points)))


def test_tunnel_profile_is_in_angstrom():
    tunnel = find_tunnels(cylinder(), **OPTIONS)[0]
    profile = tunnel_profile(tunnel)
    assert profile.method == "tunnel-path"
    assert profile.points[-1].t == pytest.approx(tunnel.length)
    assert profile.metadata["t_definition"].startswith("cumulative path length")


def test_axial_statistics_still_refuse_path_profiles_with_a_clear_reason():
    profile = tunnel_profile(find_tunnels(cylinder(), **OPTIONS)[0])
    with pytest.raises(ValueError, match="measured along its own path"):
        profile_distribution([FrameAnalysis(0, profile=profile), FrameAnalysis(1, profile=profile)])


def test_cli_tunnel_json_reports_path_length(tmp_path):
    source = write_viewer_structure(cylinder(), tmp_path / "cyl.pdb")
    out = tmp_path / "tunnels.json"
    assert main(["tunnels", str(source), "-o", str(out), "--start", "0.3,0.2,0.1", "--spacing", "2",
                 "--min-radius", "0.8", "--max-tunnels", "1", "--include-hetero"]) == 0
    tunnel = json.loads(out.read_text())["tunnels"][0]
    positions = [p["axis_position"] for p in tunnel["points"]]
    assert positions[0] == 0.0
    assert positions[-1] == pytest.approx(tunnel["length"])

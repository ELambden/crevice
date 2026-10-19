"""Worker pools under the ``spawn`` start method (the default on macOS and Windows).

Linux runs the pools with ``fork``; these tests force ``spawn`` through
``CREVICE_MP_START_METHOD`` so that every platform's start method is exercised
on Linux too, and check that parallel results equal serial ones. Software tests
on synthetic geometry only; no biological result is implied.
"""
import json
import sys
from dataclasses import replace

import numpy as np
import pytest

from crevice.cavity_trajectory import _process_context
from test_paths_autoprobe import _strip_timing, two_sided  # noqa: F401  (fixture)
from test_tracking import _shell_frame


@pytest.fixture
def spawn(monkeypatch):
    monkeypatch.setenv("CREVICE_MP_START_METHOD", "spawn")
    assert _process_context().get_start_method() == "spawn"


def test_default_start_method_is_fork_only_on_linux(monkeypatch):
    monkeypatch.delenv("CREVICE_MP_START_METHOD", raising=False)
    expected = "fork" if sys.platform.startswith("linux") else "spawn"
    assert _process_context().get_start_method() == expected


def test_parallel_paths_equal_serial_under_spawn(spawn, two_sided):  # noqa: F811
    from crevice.paths import find_paths
    serial = _strip_timing(find_paths(two_sided, seed=(0.0, 0.0, 3.75), axis="z", workers=1).to_dict())
    parallel = _strip_timing(find_paths(two_sided, seed=(0.0, 0.0, 3.75), axis="z", workers=2).to_dict())
    assert json.dumps(serial, sort_keys=True, default=float) == json.dumps(parallel, sort_keys=True, default=float)


def test_parallel_cavity_trajectory_equals_serial_under_spawn(spawn, tmp_path):
    from crevice.cavity_trajectory import CavityReference, analyze_cavity_trajectory
    from crevice.rolling import rolling_probe_cast
    from crevice.trajectory import Trajectory
    from crevice.volume_export import write_void_cast_dx
    from test_channel_coordinate import rectangular_channel, transform
    frame = rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True)
    cast = rolling_probe_cast(frame, spacing=.75, selection_mode='all', enclosure_fraction=0, min_component_volume=1)
    main = max(cast.components, key=lambda c: c.volume)
    cast = replace(cast, components=(main,), points=main.points)
    ref = CavityReference.from_dx(write_void_cast_dx(cast, tmp_path / 'ref.dx'), axis=(0, 0, 1))
    moved, _, _ = transform(frame)
    traj = Trajectory([replace(frame, frame_index=0), replace(moved, frame_index=1)], time_step=1.)
    options = dict(geometry_mode='reference-region', geometry={'probe_radius': 1.4}, radii="charmm_like")
    serial = analyze_cavity_trajectory(traj, ref, workers=1, **options)
    parallel = analyze_cavity_trajectory(traj, ref, workers=2, **options)
    assert [r['volume_A3'] for r in serial['frames']] == [r['volume_A3'] for r in parallel['frames']]
    assert parallel['settings']['radii'] == serial['settings']['radii']


def _write_models(frames, path):
    with open(path, "w", encoding="utf-8") as handle:
        for model, frame in enumerate(frames, 1):
            handle.write(f"MODEL     {model:4d}\n")
            for a in frame.atoms:
                handle.write(f"ATOM  {a.serial:5d}  CA  ALA A{a.resid % 10000:4d}    "
                             f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C\n")
            handle.write("ENDMDL\n")
        handle.write("END\n")


@pytest.mark.filterwarnings("ignore:Reader has no dt information")
def test_parallel_cavity_tracking_equals_serial_under_spawn(spawn, tmp_path):
    from crevice.tracking import track_cavities
    frame = _shell_frame([(0, 0, 0)])
    rng = np.random.default_rng(1)
    frames = [replace(frame, atoms=tuple(replace(a, x=a.x + d[0], y=a.y + d[1], z=a.z + d[2])
                                         for a, d in zip(frame.atoms, rng.normal(0, s, (len(frame.atoms), 3)))))
              for s in (0.0, 0.2, 0.4)]
    _write_models(frames, tmp_path / "models.pdb")
    options = dict(spacing=0.5, probe_radius=1.0, min_volume=1.0, centre=(0, 0, 0), radius=6, selection="all")
    serial = track_cavities(str(tmp_path / "models.pdb"), workers=1, **options)
    parallel = track_cavities(str(tmp_path / "models.pdb"), workers=2, **options)
    assert [r["status"] for r in serial["frames"]] == ["measured"] * 3
    for result in (serial, parallel):
        result["settings"].pop("workers")
    assert json.dumps(serial, sort_keys=True, default=str) == json.dumps(parallel, sort_keys=True, default=str)

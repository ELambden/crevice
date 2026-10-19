"""Cavity identity tracking (crevice.tracking) on synthetic voxel blobs and frames.

Software tests only: hand-built voxel sets with known motion, splitting,
merging, appearance and disappearance, and one synthetic structure series. No
biological result is implied.
"""
import csv

import numpy as np
import pytest

from crevice.models import Atom, StructureFrame
from crevice.tracking import CavityTracker, TrackedComponent, measure_frame, voxel_keys, write_tracking_outputs


def ball(centre, radius=3.0, spacing=0.5):
    """Voxel component: lattice nodes within ``radius`` of ``centre``."""
    r = int(radius / spacing) + 1
    c = np.rint(np.asarray(centre) / spacing).astype(int)
    grid = np.array(np.meshgrid(*[np.arange(v - r, v + r + 1) for v in c], indexing="ij")).reshape(3, -1).T * spacing
    pts = grid[np.linalg.norm(grid - np.asarray(centre), axis=1) <= radius]
    return TrackedComponent(voxel_keys(pts, spacing), len(pts) * spacing ** 3, tuple(pts.mean(axis=0)), ["A:ALA1"])


def events(tracker):
    return [(e["frame_index"], e["event"]) for e in tracker.events]


def test_moving_blob_keeps_identity():
    tracker = CavityTracker()
    for k in range(6):
        tracker.update(k, 10.0 * k, [ball((0.6 * k, 0.0, 0.0))])
    assert events(tracker) == []
    assert {o["identity"] for o in tracker.observations} == {1}
    life = tracker.lifetimes()
    assert life[0]["frames_present"] == 6 and life[0]["end_reason"] == "present_at_end"


def test_split_gives_split_event_and_new_identity():
    tracker = CavityTracker()
    tracker.update(0, 0.0, [ball((0, 0, 0), 5.0)])
    tracker.update(1, 1.0, [ball((-2.5, 0, 0), 2.6), ball((2.5, 0, 0), 2.4)])
    assert events(tracker) == [(1, "split")]
    split = tracker.events[0]
    assert split["parents"] == [1] and sorted(split["children"]) == [1, 2]
    rows = {r["identity"]: r for r in tracker.lifetimes()}
    assert rows[2]["start_reason"] == "split" and rows[2]["parents"] == "1"


def test_merge_gives_merge_event():
    tracker = CavityTracker()
    tracker.update(0, 0.0, [ball((-2.5, 0, 0), 2.6), ball((2.5, 0, 0), 2.4)])
    tracker.update(1, 1.0, [ball((0, 0, 0), 5.0)])
    assert events(tracker) == [(1, "merge")]
    merge = tracker.events[0]
    assert sorted(merge["parents"]) == [1, 2] and merge["children"] == [1]
    rows = {r["identity"]: r for r in tracker.lifetimes()}
    assert rows[2]["end_reason"] == "merge" and rows[2]["merged_into"] == 1


def test_birth_and_death():
    tracker = CavityTracker()
    tracker.update(0, 0.0, [ball((0, 0, 0))])
    tracker.update(1, 1.0, [ball((0, 0, 0)), ball((20, 0, 0))])
    tracker.update(2, 2.0, [ball((20, 0, 0))])
    assert events(tracker) == [(1, "birth"), (2, "death")]
    rows = {r["identity"]: r for r in tracker.lifetimes()}
    assert rows[1]["end_reason"] == "death" and rows[2]["start_reason"] == "birth"
    assert rows[1]["lifetime_ps"] == 1.0


def test_jump_beyond_overlap_is_death_plus_birth():
    tracker = CavityTracker()
    tracker.update(0, 0.0, [ball((0, 0, 0))])
    tracker.update(1, 1.0, [ball((7, 0, 0))])
    assert sorted(e for _, e in events(tracker)) == ["birth", "death"]


def test_thresholds_validated():
    with pytest.raises(ValueError):
        CavityTracker(link_jaccard=0)


def _shell_frame(cavities, serial=0):
    """Block of atoms with spherical cavities carved at the given centres (radius 4 Å)."""
    g = np.arange(-12, 12.01, 1.5)
    grid = np.array(np.meshgrid(g, g, g, indexing="ij")).reshape(3, -1).T
    keep = np.ones(len(grid), bool)
    for c in cavities:
        keep &= np.linalg.norm(grid - np.asarray(c), axis=1) > 4.0
    atoms = tuple(Atom(i + 1, "CA", "ALA", "A", i + 1, *map(float, p), "C") for i, p in enumerate(grid[keep]))
    return StructureFrame(atoms, frame_index=serial)


def test_measure_frame_and_outputs(tmp_path):
    reference = _shell_frame([(0, 0, 0)])
    tracker = CavityTracker()
    frames = [reference, _shell_frame([(0, 0, 0)], 1)]
    for frame in frames:
        rec = measure_frame(frame, reference, spacing=0.5, probe_radius=1.0, min_volume=1.0, centre=(0, 0, 0), radius=6)
        assert rec["status"] == "measured" and len(rec["components"]) == 1
        tracker.update(frame.frame_index, float(frame.frame_index), rec["components"])
    assert tracker.events == []
    result = {"frames": [{"frame_index": f.frame_index, "time_ps": float(f.frame_index), "status": "measured",
                          "components": 1} for f in frames],
              "events": tracker.events, "observations": tracker.observations, "lifetimes": tracker.lifetimes(),
              "settings": {}, "interpretation": "test"}
    files = write_tracking_outputs(result, tmp_path, prefix="t")
    assert {p.name for p in files} >= {"t_tracking_events.csv", "t_tracking_identities.csv",
                                        "t_tracking_lifetimes.csv", "t_tracking_frames.csv", "t_tracking.json",
                                        "t_tracking.png"}
    rows = list(csv.DictReader((tmp_path / "t_tracking_identities.csv").open(encoding="utf-8")))
    assert len(rows) == 2 and rows[0]["lining_residues"]

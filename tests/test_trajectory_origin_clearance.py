"""A fixed channel origin that loses its clearance in one trajectory frame.

Synthetic software tests: a 12-atom-ring cylinder (radius 6 A, along z) whose
second frame has one wall atom moved onto the origin. That frame must be
recorded as unresolved (with the reason) while the other frames keep exactly
the profiles they would have without it.
"""

import csv
import json
import math
import warnings

from dataclasses import replace

import pytest

from crevice.models import Atom, StructureFrame
from crevice.sections import OriginClearanceError
from crevice.trajectory import OriginClearanceWarning, analyze_trajectory

ORIGIN = (0.0, 0.0, 0.0)


@pytest.fixture(autouse=True)
def _short_ladder(monkeypatch):
    """Two automatic enclosure-probe candidates keep these tests fast."""
    import crevice.sections as sections
    monkeypatch.setattr(sections, "AUTO_ENCLOSURE_LADDER", (0.5, 0.8))


def _cylinder(*, blocked=False):
    """Ring-stacked cylinder; with ``blocked`` the wall atom at z = 0, angle 0 sits on the origin."""
    atoms = []
    for k in range(15):
        for j in range(12):
            x, y, z = 6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6), -10.5 + 1.5 * k
            if blocked and k == 7 and j == 0:
                x, y, z = ORIGIN
            atoms.append(Atom(12 * k + j + 1, f"C{j}", "ALA", "A", k + 1, x, y, z, "C"))
    return StructureFrame(tuple(atoms))


def _intruded(distance):
    """The cylinder with its wall atom at z = 0, angle 0 moved to ``distance`` A from the axis."""
    frame = _cylinder()
    atoms = list(frame.atoms)
    atoms[7 * 12] = replace(atoms[7 * 12], x=distance, y=0.0, z=0.0)
    return StructureFrame(tuple(atoms))


def _profile_dict(frame):
    profile = frame.profile
    return json.dumps([profile.to_dict(), [p.to_dict() for p in profile.points]], sort_keys=True)


@pytest.mark.parametrize("enclosure_radius", [None, 0.8])
def test_origin_without_clearance_in_one_frame_is_recorded_and_the_run_continues(enclosure_radius):
    open_, blocked = _cylinder(), _cylinder(blocked=True)
    kwargs = {"axis": "z", "origin": ORIGIN, "enclosure_radius": enclosure_radius}
    with pytest.warns(OriginClearanceWarning, match=r"1 of 3 frames are unresolved.*frames 1\)"):
        result = analyze_trajectory([open_, blocked, open_], align=False, profile_kwargs=kwargs)
    first, middle, last = result.frames
    assert first.profile is not None and last.profile is not None
    assert middle.profile is None
    assert middle.metadata["profile_status"] == "unresolved"
    assert "Channel origin must have clearance" in middle.metadata["profile_reason"]
    probe = middle.metadata["enclosure_probe"]
    assert probe["source"] == "none" and probe["origin_clearance_A"] < 0
    # Automatic probe: the per-frame fallback is tried and fails the same way.
    assert ("fallback_failure" in probe) == (enclosure_radius is None)
    assert result.metadata["resolved_profile_count"] == 2
    assert result.metadata["origin_clearance"]["unresolved_frames"] == 1
    assert result.metadata["origin_clearance"]["frame_indices"] == [1]
    assert result.metadata["enclosure_probe"]["fallback"]["unresolved_frames"] == 1

    # Frames that resolve are unchanged by the unresolved frame between them.
    with warnings.catch_warnings():
        warnings.simplefilter("error", OriginClearanceWarning)
        reference = analyze_trajectory([open_, open_], align=False, profile_kwargs=kwargs)
    assert "origin_clearance" not in reference.metadata
    for frame in (first, last):
        assert _profile_dict(frame) == _profile_dict(reference.frames[0])


def test_origin_without_clearance_raises_when_asked_including_on_the_reference_frame():
    open_, blocked = _cylinder(), _cylinder(blocked=True)
    kwargs = {"axis": "z", "origin": ORIGIN}
    with pytest.raises(OriginClearanceError, match="clearance at the origin"):
        analyze_trajectory([open_, blocked], align=False, profile_kwargs=kwargs, on_profile_error="raise")
    with pytest.raises(OriginClearanceError) as caught:
        analyze_trajectory([blocked, open_], align=False, profile_kwargs=kwargs, on_profile_error="raise")
    assert isinstance(caught.value, ValueError) and caught.value.clearance < 0


def _no_reference_axis(result):
    """The per-frame fields that an unresolved reference frame leaves (status, probe record)."""
    return [(f.profile, f.metadata["profile_status"], f.metadata["enclosure_probe"]["source"]) for f in result.frames]


@pytest.mark.parametrize("enclosure_radius", [None, 0.8])
def test_origin_without_clearance_on_the_reference_frame_is_recorded_like_an_unresolved_reference(enclosure_radius):
    # On the reference frame no axis or probe can be fixed, so (as for a reference
    # frame whose channel is unresolved) no frame is profiled; later frames are not attempted.
    open_, blocked = _cylinder(), _cylinder(blocked=True)
    kwargs = {"axis": "z", "origin": ORIGIN, "enclosure_radius": enclosure_radius}
    with pytest.warns(OriginClearanceWarning, match=r"reference frame \(frame 0; clearance -.*no frame is profiled") as caught:
        result = analyze_trajectory([blocked, open_, open_], align=False, profile_kwargs=kwargs)
    assert len([w for w in caught if issubclass(w.category, OriginClearanceWarning)]) == 1
    assert _no_reference_axis(result) == [(None, "not_attempted_reference_unresolved", "not_attempted")] * 3
    reason = result.metadata["reference_profile_reason"]
    assert "Channel origin must have clearance" in reason
    assert result.metadata["reference_profile_status"] == "unresolved"
    assert all(f.metadata["profile_reason"] == reason for f in result.frames)
    assert result.frames[0].metadata["enclosure_probe"]["origin_clearance_A"] < 0
    assert all("origin_clearance_A" not in f.metadata["enclosure_probe"] for f in result.frames[1:])
    record = result.metadata["origin_clearance"]
    assert (record["unresolved_frames"], record["frame_indices"], record["reference_frame"]) == (1, [0], True)
    assert result.metadata["resolved_profile_count"] == 0
    assert result.metadata["enclosure_probe"]["chosen_A"] == enclosure_radius
    assert result.metadata["enclosure_probe"]["fallback"]["enabled"] is False

    # Same per-frame outcome as a reference frame whose channel is unresolved
    # (here: an axis across the cylinder, along which no channel has two open ends).
    with warnings.catch_warnings():
        warnings.simplefilter("error", OriginClearanceWarning)
        unresolved = analyze_trajectory([open_, open_, open_], align=False,
                                        profile_kwargs={**kwargs, "axis": "x"})
    assert result.metadata["reference_profile_status"] == unresolved.metadata["reference_profile_status"]
    assert _no_reference_axis(result) == _no_reference_axis(unresolved)
    assert "origin_clearance" not in unresolved.metadata


def test_a_fallback_probe_can_still_resolve_a_frame_whose_origin_is_too_tight_for_the_pinned_probe():
    # Clearance at the origin 0.9 A: too little for the pinned 1.5 A probe, enough for 0.5 A.
    with warnings.catch_warnings():
        warnings.simplefilter("error", OriginClearanceWarning)
        result = analyze_trajectory([_cylinder(), _intruded(2.6)], align=False, probe_fallback=True,
                                    profile_kwargs={"axis": "z", "origin": ORIGIN, "enclosure_radius": 1.5})
    frame = result.frames[1]
    probe = frame.metadata["enclosure_probe"]
    assert frame.profile is not None and probe["source"] == "fallback" and probe["radius_A"] == 0.5
    assert "clearance at the origin 0.90 A" in probe["pinned_failure"]
    assert probe["origin_clearance_A"] == pytest.approx(0.9)
    assert "origin_clearance" not in result.metadata


def _write_pdb(frame, path):
    path.write_text("\n".join(
        f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00  0.00           C"
        for a in frame.atoms) + "\nEND\n", encoding="utf-8", newline="\n")
    return str(path)


def test_trajectory_cli_records_the_frame_and_warns(tmp_path, capsys):
    from crevice.cli import main
    open_, blocked = _cylinder(), _cylinder(blocked=True)
    paths = [_write_pdb(f, tmp_path / f"f{i}.pdb") for i, f in enumerate((open_, blocked, open_))]
    output = tmp_path / "run.json"
    with warnings.catch_warnings():
        warnings.simplefilter("error", OriginClearanceWarning)
        assert main(["trajectory", *paths, "-o", str(output), "--axis", "z", "--origin", "0,0,0",
                     "--no-align", "--confidence", "0", "--skip-hydration", "--skip-analysis-views"]) == 0
    err = capsys.readouterr().err
    assert "crevice: warning: 1 of 3 frames are unresolved where --origin has too little clearance" in err
    with (tmp_path / "run_frames.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [r["profile_status"] for r in rows] == ["resolved", "unresolved", "resolved"]
    assert "Channel origin must have clearance" in rows[1]["profile_reason"]
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["metadata"]["origin_clearance"]["frame_indices"] == [1]


def test_trajectory_cli_records_a_reference_frame_without_clearance_and_exits_0(tmp_path, capsys):
    from crevice.cli import main
    open_, blocked = _cylinder(), _cylinder(blocked=True)
    paths = [_write_pdb(f, tmp_path / f"f{i}.pdb") for i, f in enumerate((blocked, open_, open_))]
    output = tmp_path / "run.json"
    args = ["trajectory", *paths, "-o", str(output), "--axis", "z", "--origin", "0,0,0",
            "--no-align", "--confidence", "0", "--skip-hydration", "--skip-analysis-views"]
    with warnings.catch_warnings():
        warnings.simplefilter("error", OriginClearanceWarning)
        assert main(args) == 0
    err = capsys.readouterr().err
    assert err.count("crevice: warning:") == 1
    assert "--origin has too little clearance for the enclosure probe on the reference frame (source frame 0)" in err
    with (tmp_path / "run_frames.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [r["profile_status"] for r in rows] == ["not_attempted_reference_unresolved"] * 3
    assert all("Channel origin must have clearance" in r["profile_reason"] for r in rows)
    data = json.loads(output.read_text(encoding="utf-8"))
    assert data["metadata"]["origin_clearance"]["frame_indices"] == [0]
    assert data["metadata"]["origin_clearance"]["reference_frame"] is True
    # --profile-errors raise still stops the run.
    assert main(args + ["--profile-errors", "raise"]) != 0
    assert "Channel origin must have clearance" in capsys.readouterr().err

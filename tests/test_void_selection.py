"""Final cast selection contracts, distinct from continuous-path validation."""

from dataclasses import replace

import pytest

from crevice import StructureFrame, detect_cavities, void_cast
from crevice.models import ChannelPoint, PoreProfile
from crevice.volume_export import cast_grid
from crevice import voids
from test_ensemble_cast_network import atom, cage


def candidate_frame(monkeypatch, indices):
    """Inject a known candidate graph to isolate selection from atom clearance."""
    frame = StructureFrame((atom(1, (0, 0, 0)), atom(2, (20, 20, 20))))

    def collect(_spatial, *, min_corner, spacing, **_kwargs):
        return {idx: voids._GridNode(voids._grid_point(min_corner, spacing, idx),
                                    2.0, frame.atoms[0], 8) for idx in indices}

    monkeypatch.setattr(voids, "_collect_void_nodes", collect)
    # These injected graphs test selection only; physical edges have separate tests.
    monkeypatch.setattr(voids.GridConnectivity, "allows", lambda *args: True)
    return frame


def line_profile(x=1, blocked=()):
    return PoreProfile((0, 0, 0), (0, 0, 1), tuple(
        ChannelPoint(z, (x, 1, z), 0 if z in blocked else 1.1,
                     -1 if z in blocked else 1.1, z) for z in range(21)))


def test_volume_rejection_does_not_hide_smaller_atomistic_cavity():
    smaller = tuple(replace(a, x=a.x * 0.8 + 16, y=a.y * 0.8, z=a.z * 0.8,
                            serial=a.serial + 1000, resid=a.resid + 1000) for a in cage())
    frame = StructureFrame(cage() + smaller)
    options = dict(mode="cavity", axis="z", spacing=1, min_radius=0.5)
    full = void_cast(frame, **options)
    assert full.component_count == 2
    large, small = full.components
    assert large.volume > small.volume
    ceiling = (large.volume + small.volume) / 2
    selected = void_cast(frame, **options, max_components=1, max_component_volume=ceiling)
    assert selected.component_count == 1
    assert selected.total_volume == small.volume
    assert selected.components[0].center == small.center
    cavities = detect_cavities(frame, spacing=1, min_radius=0.5,
                               max_cavities=1, max_component_volume=ceiling)
    assert len(cavities) == 1
    assert cavities[0].volume == small.volume


def test_centerline_miss_does_not_hide_another_source(monkeypatch):
    frame = candidate_frame(monkeypatch, [(1, 1, z) for z in range(16)] +
                            [(7, 1, z) for z in range(9)])
    cast = void_cast(frame, mode="channel", axis="z", spacing=1,
                     profile=line_profile(x=7), max_components=1)
    assert cast.component_count == 1
    assert cast.total_volume == 9
    assert all(p.position[0] == 7 for p in cast.components[0].points)


def test_final_limit_and_volume_bounds_apply_after_centerline_split(monkeypatch):
    frame = candidate_frame(monkeypatch, [(1, 1, z) for z in range(21)])
    options = dict(mode="channel", axis="z", spacing=1,
                   profile=line_profile(blocked=range(8, 13)),
                   min_component_volume=8, max_component_volume=8)
    full = void_cast(frame, **options)
    assert full.component_count == 2
    limited = void_cast(frame, **options, max_components=1, max_export_points=1)
    assert limited.component_count == 1
    assert limited.total_volume == 8
    assert len(limited.points) == 1
    assert cast_grid(limited)[0].sum() == 8
    assert limited.components[0].id == 1
    assert limited.metadata["selection"]["eligible_component_count"] == 2
    assert limited.metadata["selection"]["omitted_by_limit_count"] == 1


@pytest.mark.parametrize("mode", ["all", "auto"])
def test_focus_preserves_original_open_source_topology(monkeypatch, mode):
    frame = candidate_frame(monkeypatch, [(1, 1, z) for z in range(21)])
    cast = void_cast(frame, mode=mode, axis="z", spacing=1, profile=line_profile(),
                     focus_points=[(1, 1, 10)], focus_radius=5)
    assert cast.component_count == 1
    component = cast.components[0]
    assert component.boundary_faces == ()
    assert component.kind != "buried_cavity"
    assert component.metadata["source_kind"] == "through_channel"
    assert component.metadata["source_boundary_faces"] == ("z+", "z-")
    assert component.metadata["selection_cropped"] is True
    assert component.volume == 11
    if mode == "auto":
        assert cast.metadata["export_mode"] == "centerline_fill"


def test_selection_decisions_and_contiguous_output_ids(monkeypatch):
    frame = candidate_frame(monkeypatch, [(x, 1, z) for x, n in ((1, 12), (7, 10), (13, 8))
                                         for z in range(n)])
    cast = void_cast(frame, mode="all", axis="z", spacing=1, max_components=1,
                     min_component_volume=8, max_component_volume=10)
    assert cast.component_count == 1
    assert cast.components[0].id == 1
    assert cast.total_volume == 10
    candidates = cast.metadata["selection"]["candidates"]
    assert [row["volume"] for row in candidates] == [12, 10, 8]
    assert [row["decision"] for row in candidates] == ["above_max_volume", "selected", "component_limit"]
    assert [row["output_component_id"] for row in candidates] == [None, 1, None]
    assert cast.components[0].metadata["candidate_component_id"] == candidates[1]["candidate_component_id"]
    assert len(cast.metadata["selected_source_components"]) == 1


@pytest.mark.parametrize("limits, decision", [
    ({"min_component_volume": 9}, "below_min_volume"),
    ({"max_component_volume": 7}, "above_max_volume"),
])
def test_all_rejected_is_empty_with_explanation(monkeypatch, limits, decision):
    frame = candidate_frame(monkeypatch, [(1, 1, z) for z in range(8)])
    cast = void_cast(frame, mode="all", axis="z", spacing=1, **limits)
    assert cast.component_count == cast.point_count == 0
    assert cast.total_volume == 0
    assert cast.metadata["selected_source_components"] == []
    assert cast.metadata["selection"]["candidates"][0]["decision"] == decision


def test_equal_volume_selection_is_deterministic(monkeypatch):
    indices = [(x, 1, z) for x in (1, 9) for z in range(8)]
    outputs = []
    for order in (indices, list(reversed(indices)), indices[::2] + indices[1::2]):
        frame = candidate_frame(monkeypatch, order)
        cast = void_cast(frame, mode="all", axis="z", spacing=1, max_components=1)
        outputs.append(cast.to_dict(include_points=True))
        assert all(p.position[0] == 1 for p in cast.components[0].points)
    assert outputs[0] == outputs[1] == outputs[2]


@pytest.mark.parametrize("limit", [0, -1, 1.5, float("nan"), float("inf"), True])
def test_component_limit_requires_positive_integer(limit):
    frame = StructureFrame((atom(1, (0, 0, 0)),))
    with pytest.raises(ValueError, match="max_components"):
        void_cast(frame, max_components=limit)


def test_small_crop_reports_point_threshold_exclusion(monkeypatch):
    frame = candidate_frame(monkeypatch, [(1, 1, z) for z in range(21)])
    cast = void_cast(frame, mode="channel", axis="z", spacing=1,
                     profile=line_profile(blocked=range(3, 21)))
    assert cast.component_count == 0
    assert cast.metadata["cropped_point_count"] == 3
    assert cast.metadata["selection"]["discarded_small_component_points"] == 3
    assert cast.metadata["selection"]["candidates"] == []


def test_metadata_records_the_spacing_actually_used():
    """A grid-point cap coarsens the requested spacing; both must be reported."""
    frame = StructureFrame(cage())
    fine = voids.void_cast(frame, mode="all", spacing=1.0, min_radius=0.4,
                           max_grid_points=10_000_000)
    assert fine.metadata["effective_spacing_A"] == 1.0
    assert fine.metadata["spacing_adjusted"] is False

    capped = voids.void_cast(frame, mode="all", spacing=1.0, min_radius=0.4,
                             max_grid_points=64)
    assert capped.metadata["spacing_adjusted"] is True
    assert capped.metadata["requested_spacing_A"] == 1.0
    assert capped.metadata["effective_spacing_A"] > 1.0
    # Every volume and radius downstream derives from the effective spacing.
    assert capped.metadata["grid_point_count"] <= fine.metadata["grid_point_count"]

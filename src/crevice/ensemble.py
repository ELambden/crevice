"""Pool channel profiles from many frames into radius distributions.

Each frame's radius profile is linearly interpolated onto one shared axial
coordinate grid (Å, measured along the common aligned axis), and per-position
statistics are computed across frames: mean, median and empirical quantile
bands (10th/90th percentile by default). Optionally a confidence band for the
mean is added with a serial-correlation-aware bootstrap
(:func:`crevice.uncertainty.block_mean_confidence`).

Interpretation
    Quantile bands describe how much the radius *fluctuates* between frames.
    They are not confidence intervals. Confidence bands (when requested)
    describe statistical uncertainty of the mean for a stationary trajectory,
    not force-field, sampling or geometric-model error.

Main entry points: :func:`profile_distribution` and
:func:`write_profile_distribution_csv`.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Sequence

from .models import FrameAnalysis, TrajectoryAnalysis


def profile_distribution(
    trajectory: TrajectoryAnalysis | Sequence[FrameAnalysis],
    *,
    quantiles: tuple[float, float] = (0.1, 0.9),
    samples: int = 101,
    coordinate_grid: Sequence[float] | None = None,
    confidence: float | None = None,
    block_length: int | None = None,
    bootstrap_replicates: int = 2000,
    seed: int = 20260911,
    assume_aligned: bool = False,
) -> dict:
    """Interpolate radii on a shared axial coordinate without extrapolation.

    Inputs must share a physical reference frame and axis orientation (for
    example the output of :func:`crevice.trajectory.analyze_trajectory` with
    alignment). Origins may differ along the axis: their axial offsets are
    included. The default grid covers only the axial interval that every
    profile spans. Explicit grids retain missing coverage as null values.

    Parameters
    ----------
    trajectory : TrajectoryAnalysis or sequence of FrameAnalysis
        Frames without a profile count towards the total but contribute no
        radii (see ``summary_population`` in the result).
    quantiles : tuple of float, default (0.1, 0.9)
        Lower and upper quantiles of the band, ``0 <= lower < upper <= 1``.
    samples : int, default 101
        Number of grid positions for the default common-interval grid.
    coordinate_grid : sequence of float, optional
        Explicit strictly increasing axial positions (Å, relative to the first
        profile's origin). Positions outside a frame's profile are missing for
        that frame.
    confidence : float, optional
        If given (for example 0.95), add a bootstrap confidence band for the
        mean at each position. Requires aligned input (see ``assume_aligned``)
        and regularly spaced, increasing frame times.
    block_length : int, optional
        Batch length in frames for the bootstrap; chosen automatically when
        omitted (see :func:`crevice.uncertainty.block_mean_confidence`).
    bootstrap_replicates : int, default 2000
        Bootstrap replicates (at least 200).
    seed : int, default 20260911
        Random seed for the bootstrap.
    assume_aligned : bool, default False
        Declare that the frames were aligned. Not needed for a
        :class:`~crevice.models.TrajectoryAnalysis` whose metadata records
        ``aligned_to_first_frame``.

    Returns
    -------
    dict
        ``rows``: one dictionary per grid position with ``coordinate`` (Å),
        ``n_frames`` (frames covering it), ``coverage_fraction`` (of all
        frames), ``mean_radius``, ``median_radius``, ``lower_radius`` and
        ``upper_radius`` (Å; ``None`` without coverage) and, with
        ``confidence``, ``mean_ci_lower``/``mean_ci_upper`` (pointwise) and
        ``mean_band_lower``/``mean_band_upper`` (simultaneous). Top-level keys
        record the axis, probe radius, frame counts, the interval result and
        the definitions used (``band_kind="empirical_frame_quantiles"``,
        ``interpolation="linear_no_extrapolation"``), plus
        ``enclosure_radii_A`` (the distinct enclosure probes pooled) and
        ``enclosure_probe_fallback_frame_count`` (frames profiled with the
        per-frame fallback probe of
        :func:`crevice.trajectory.analyze_trajectory`; only those may differ
        from the reference's enclosure probe).

    Raises
    ------
    ValueError
        For invalid quantiles or samples; no profiles; methods other than
        ``"axial-refined"``/``"axial-connected"``; profiles with different
        methods, probe radii, atom selections, section settings or axis
        directions, or origins offset perpendicular to the axis; non-increasing
        coordinates or negative radii; no common interval; an invalid
        ``coordinate_grid``; or a confidence band requested for unaligned or
        irregularly sampled frames.

    Notes
    -----
    Profiles from fixed-axis scans (``search_radius=0``) get no confidence band;
    the interval status is ``"unavailable_unvalidated_axis_scan"``.

    Examples
    --------
    >>> from crevice.models import ChannelPoint, FrameAnalysis, PoreProfile
    >>> from crevice.ensemble import profile_distribution
    >>> def frame(i, radii):
    ...     points = tuple(ChannelPoint(k, (0.0, 0.0, float(k)), r, r, float(k))
    ...                    for k, r in enumerate(radii))
    ...     return FrameAnalysis(i, profile=PoreProfile((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), points))
    >>> result = profile_distribution([frame(0, [2.0, 1.0, 2.0]), frame(1, [3.0, 2.0, 3.0])],
    ...                               samples=3)
    >>> [(row["coordinate"], row["mean_radius"]) for row in result["rows"]]
    [(0.0, 2.5), (1.0, 1.5), (2.0, 2.5)]
    """
    import numpy as np

    if len(quantiles) != 2 or not 0 <= quantiles[0] < quantiles[1] <= 1:
        raise ValueError("quantiles must satisfy 0 <= lower < upper <= 1")
    if samples < 2:
        raise ValueError("samples must be at least 2")
    frames = tuple(trajectory.frames if isinstance(trajectory, TrajectoryAnalysis) else trajectory)
    profiles = [frame.profile for frame in frames if frame.profile is not None]
    profile_indices = [i for i,frame in enumerate(frames) if frame.profile is not None]
    if not profiles:
        raise ValueError("No profile data available")
    reference = profiles[0]
    if reference.method not in {"axial-refined", "axial-connected"}:
        raise ValueError(
            f"An axial distribution requires axial-refined or axial-connected profiles; a "
            f"{reference.method!r} profile is measured along its own path (t in Å from its start), "
            "which is not a shared axis across frames")
    direction = np.asarray(reference.axis_direction, dtype=float)
    if not np.isfinite(direction).all() or not np.isclose(np.linalg.norm(direction), 1):
        raise ValueError("Profile axes must be unit vectors")
    coordinates, radii = [], []
    # Frames profiled by the per-frame enclosure-probe fallback of
    # analyze_trajectory carry their own probe; only their probe may differ.
    fallback = [((frames[i].metadata or {}).get("enclosure_probe") or {}).get("source") == "fallback"
                for i in profile_indices]
    for profile, is_fallback in zip(profiles, fallback):
        if profile.method != reference.method or profile.probe_radius != reference.probe_radius:
            raise ValueError("Cannot pool profiles with different methods or probe radii")
        keys = ("include_hydrogen", "include_hetero", "section_spacing_A") + (() if is_fallback else ("enclosure_radius_A",))
        if any(profile.metadata.get(k) != reference.metadata.get(k) for k in keys):
            raise ValueError("Cannot pool profiles with different atom selections")
        if not np.allclose(profile.axis_direction, direction, atol=1e-6, rtol=0):
            raise ValueError("Profiles must have a common aligned axis orientation")
        delta_origin=np.asarray(profile.axis_origin)-reference.axis_origin
        offset=np.dot(delta_origin,direction)
        if np.linalg.norm(delta_origin-offset*direction)>1e-5:
            raise ValueError("Profiles must share the same aligned channel origin perpendicular to the axis")
        x = np.array([p.t + offset for p in profile.points], dtype=float)
        y = np.array([p.radius for p in profile.points], dtype=float)
        if len(x) < 2 or not np.isfinite(x).all() or not np.all(np.diff(x) > 0):
            raise ValueError("Profile coordinates must be finite and strictly increasing")
        if not np.isfinite(y).all() or np.any(y < 0):
            raise ValueError("Profile radii must be finite and non-negative")
        coordinates.append(x)
        radii.append(y)
    if coordinate_grid is None:
        start, end = max(x[0] for x in coordinates), min(x[-1] for x in coordinates)
        if start >= end:
            raise ValueError("Profiles have no common axial interval")
        grid = np.linspace(start, end, samples)
    else:
        grid = np.asarray(coordinate_grid, dtype=float)
        if grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all() or not np.all(np.diff(grid) > 0):
            raise ValueError("coordinate_grid must be finite and strictly increasing")
    matrix = np.full((len(frames),len(grid)),np.nan)
    for i,x,y in zip(profile_indices,coordinates,radii):
        matrix[i]=np.interp(grid,x,y,left=np.nan,right=np.nan)
    interval = None
    if confidence is not None:
        from .uncertainty import block_mean_confidence
        aligned = assume_aligned or (isinstance(trajectory,TrajectoryAnalysis) and trajectory.metadata.get("aligned_to_first_frame",False))
        if not aligned:
            raise ValueError("Confidence bands require verified alignment or explicit assume_aligned=True")
        times = [frame.time if frame.time is not None else frame.frame_index for frame in frames]
        delta = np.diff(times)
        if len(delta) and (np.any(delta<=0) or not np.allclose(delta,delta[0],rtol=1e-5,atol=1e-8)):
            raise ValueError("Block confidence bands require ordered, regularly sampled frames")
        if any(p.metadata.get("status")=="unvalidated_axis_scan" for p in profiles):
            interval={"status":"unavailable_unvalidated_axis_scan","confidence":confidence}
        else:
            interval=block_mean_confidence(matrix,confidence=confidence,block_length=block_length,
                                          replicates=bootstrap_replicates,seed=seed)
    rows = []
    for index, coordinate in enumerate(grid):
        values = matrix[:, index]
        values = values[np.isfinite(values)]
        lower, median, upper = np.quantile(values, [quantiles[0], 0.5, quantiles[1]]) if len(values) else (None,) * 3
        rows.append({"coordinate": float(coordinate), "n_frames": int(len(values)),
                     "coverage_fraction":len(values)/len(frames),
                     "mean_radius": float(values.mean()) if len(values) else None,
                     "median_radius": float(median) if median is not None else None,
                     "lower_radius": float(lower) if lower is not None else None,
                     "upper_radius": float(upper) if upper is not None else None})
        if interval is not None:
            for source,target in [("pointwise_lower","mean_ci_lower"),("pointwise_upper","mean_ci_upper"),
                                  ("simultaneous_lower","mean_band_lower"),("simultaneous_upper","mean_band_upper")]:
                rows[-1][target]=interval.get(source,[None]*len(grid))[index]
    return {"rows": rows, "total_frame_count":len(frames),"confidence_interval":interval,
            "summary_population":"all frames" if len(profiles)==len(frames) else "conditional on resolved profiles", "quantiles": list(quantiles), "frame_count": len(profiles),
            "missing_profile_count": len(frames) - len(profiles), "units": "angstrom",
            "band_kind": "empirical_frame_quantiles", "interpolation": "linear_no_extrapolation",
            "axis_origin": reference.axis_origin, "axis_direction": reference.axis_direction,
            "probe_radius": reference.probe_radius,
            "enclosure_probe_fallback_frame_count": int(sum(fallback)),
            "enclosure_radii_A": sorted({p.metadata.get("enclosure_radius_A") for p in profiles
                                         if p.metadata.get("enclosure_radius_A") is not None})}


DISTRIBUTION_COLUMNS = {"coordinate": "coordinate_A", "mean_radius": "mean_radius_A", "median_radius": "median_radius_A",
                        "lower_radius": "lower_radius_A", "upper_radius": "upper_radius_A",
                        "mean_ci_lower": "mean_ci_lower_A", "mean_ci_upper": "mean_ci_upper_A",
                        "mean_band_lower": "mean_band_lower_A", "mean_band_upper": "mean_band_upper_A"}


def write_profile_distribution_csv(distribution: dict, path: str | Path) -> None:
    """Write the ``rows`` of a :func:`~crevice.ensemble.profile_distribution` result to CSV.

    Columns follow the keys of the first row, with the unit appended to every
    length (``coordinate_A``, ``mean_radius_A``, ``median_radius_A``,
    ``lower_radius_A``/``upper_radius_A`` (frame quantiles),
    ``mean_ci_lower_A``/``mean_ci_upper_A`` and
    ``mean_band_lower_A``/``mean_band_upper_A``); ``n_frames`` and
    ``coverage_fraction`` are unitless. A header of ``coordinate_A, n_frames,
    mean_radius_A`` is written when there are no rows. Parent directories are
    created.

    Parameters
    ----------
    distribution : dict
        Result of :func:`profile_distribution`.
    path : str or pathlib.Path
        Output file (overwritten).
    """
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        fields=list(distribution["rows"][0]) if distribution["rows"] else ["coordinate","n_frames","mean_radius"]
        writer = csv.DictWriter(handle, fieldnames=[DISTRIBUTION_COLUMNS.get(field, field) for field in fields])
        writer.writeheader()
        writer.writerows({DISTRIBUTION_COLUMNS.get(k, k): v for k, v in row.items()} for row in distribution["rows"])

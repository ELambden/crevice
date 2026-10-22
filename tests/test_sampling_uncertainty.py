"""Effective sample size, blocking and block bootstrap on series with known autocorrelation.

Synthetic software tests: AR(1) series x_t = phi x_(t-1) + e_t have integrated
autocorrelation time tau = (1 + phi) / (1 - phi) and variance 1 / (1 - phi^2),
so the standard error of the mean of N frames is sqrt(tau / (N (1 - phi^2))).
"""
import math

import numpy as np
import pytest

from crevice.uncertainty import (SAMPLING_METHODS, blocking_standard_error, integrated_autocorrelation_time,
                                 moving_block_bootstrap, sampling_uncertainty, trajectory_frame_statistics)


def ar1(phi, n, seed):
    rng = np.random.default_rng(seed)
    e = rng.normal(size=n + 1000)
    x = np.empty_like(e)
    x[0] = e[0]
    for i in range(1, len(e)):
        x[i] = phi * x[i - 1] + e[i]
    return x[1000:]


def true_se(phi, n):
    return math.sqrt((1 + phi) / (1 - phi) / n / (1 - phi ** 2))


@pytest.mark.parametrize("phi", [0.0, 0.5, 0.8])
def test_integrated_autocorrelation_time_recovers_ar1_tau(phi):
    tau = (1 + phi) / (1 - phi)
    estimates = [integrated_autocorrelation_time(ar1(phi, 5000, seed))["tau_int_frames"] for seed in range(20)]
    assert np.mean(estimates) == pytest.approx(tau, rel=0.15)
    assert integrated_autocorrelation_time(ar1(phi, 5000, 1))["window_converged"]


@pytest.mark.parametrize("phi", [0.0, 0.5, 0.8])
def test_blocking_and_bootstrap_standard_errors_match_the_ar1_value(phi):
    n = 4096
    blocking, boot = [], []
    for seed in range(20):
        r = sampling_uncertainty(ar1(phi, n, seed), seed=seed, replicates=1000)
        assert r["status"] == "estimated" and r["block_bootstrap_status"] == "estimated"
        blocking.append(r["blocking_standard_error"])
        boot.append(r["block_bootstrap_standard_error"])
        assert r["block_bootstrap_block_length_frames"] == math.ceil(5 * r["integrated_autocorrelation_time_frames"])
        assert r["blocking_block_size_frames"] & (r["blocking_block_size_frames"] - 1) == 0  # power of two
    assert np.mean(blocking) == pytest.approx(true_se(phi, n), rel=0.2)
    assert np.mean(boot) == pytest.approx(true_se(phi, n), rel=0.2)


def test_block_bootstrap_interval_covers_the_true_mean_at_about_the_nominal_rate():
    hits = 0
    for seed in range(100):
        r = sampling_uncertainty(ar1(0.6, 2000, 100 + seed), seed=seed, replicates=500)
        hits += r["block_bootstrap_lower"] <= 0.0 <= r["block_bootstrap_upper"]
    assert 0.85 <= hits / 100 <= 1.0


def test_effective_sample_size_of_white_noise_is_about_n():
    r = sampling_uncertainty(np.random.default_rng(5).normal(size=1000))
    assert r["effective_sample_size"] == pytest.approx(1000, rel=0.15)
    assert r["blocking_standard_error"] == pytest.approx(1 / math.sqrt(1000), rel=0.25)


def test_statuses_for_unusable_series():
    assert sampling_uncertainty([1.0] * 50)["status"] == "unavailable_observed_constant"
    assert sampling_uncertainty([1.0, np.nan] * 20)["status"] == "unavailable_missing_frames"
    assert sampling_uncertainty(np.arange(10.0))["status"] == "unavailable_too_few_frames"
    assert sampling_uncertainty(np.arange(50.0), regular=False)["status"] == "unavailable_irregular_timestamps"
    short = sampling_uncertainty(ar1(0.97, 300, 3))
    assert short["block_bootstrap_status"] == "insufficient_blocks" and short["effective_sample_size"] is not None
    assert sampling_uncertainty(ar1(0.5, 400, 3), confidence=0)["block_bootstrap_status"] == "not_requested"
    user = sampling_uncertainty(ar1(0.5, 400, 3), block_length=10)
    assert user["block_bootstrap_block_length_frames"] == 10 and user["block_bootstrap_block_length_source"] == "user"
    assert set(SAMPLING_METHODS) == {"effective_sample_size", "blocking", "block_bootstrap"}


def test_bootstrap_is_reproducible_and_validates_input():
    x = ar1(0.3, 300, 9)
    assert moving_block_bootstrap(x, block_length=5, seed=1) == moving_block_bootstrap(x, block_length=5, seed=1)
    with pytest.raises(ValueError):
        moving_block_bootstrap(x, block_length=0)
    with pytest.raises(ValueError):
        blocking_standard_error([1.0, float("nan")])


def test_describe_series_and_summary_rows_carry_the_new_fields_without_changing_old_ones():
    from crevice.cavity_trajectory import describe_series
    from crevice.io import SUMMARY_STATISTICS_FIELDS, summary_statistics_rows
    x = ar1(0.5, 400, 11)
    record = describe_series(x)
    assert record["sampling_uncertainty"]["status"] == "estimated"
    old = {k: v for k, v in record.items() if k != "sampling_uncertainty"}
    assert old["mean"] == pytest.approx(float(np.mean(x))) and "confidence_interval" in old
    row = summary_statistics_rows({"volume_A3": record})[0]
    assert SUMMARY_STATISTICS_FIELDS[:27][-1] == "replicates"
    assert row["effective_sample_size"] == record["sampling_uncertainty"]["effective_sample_size"]
    assert row["block_bootstrap_block_length_frames"] == record["sampling_uncertainty"]["block_bootstrap_block_length_frames"]


def test_trajectory_frame_statistics_on_a_synthetic_channel():
    from crevice import analyze_trajectory, trajectory_from_frames
    from crevice.models import Atom, StructureFrame
    rng = np.random.default_rng(2)
    atoms = []
    for ring, z in enumerate(np.arange(-6.0, 6.1, 1.5)):
        for k in range(8):
            a = 2 * math.pi * k / 8
            atoms.append((4.0 * math.cos(a), 4.0 * math.sin(a), float(z)))
    base = np.asarray(atoms)
    frames = [StructureFrame(tuple(Atom(i + 1, "CA", "ALA", "A", i + 1, *p, "C")
                                   for i, p in enumerate(base + rng.normal(scale=0.1, size=base.shape))), frame_index=j)
              for j in range(24)]
    result = analyze_trajectory(trajectory_from_frames(frames, time_step=5.0), analyses=("profile",),
                                profile_kwargs={"axis": "z", "samples": 15})
    stats = trajectory_frame_statistics(result, replicates=200)
    assert {"min_radius_A", "mean_radius_A", "max_radius_A", "bottleneck_axis_position_A"} <= set(stats["quantities"])
    sampling = stats["quantities"]["min_radius_A"]["sampling_uncertainty"]
    assert sampling["frame_count"] == 24 and sampling["effective_sample_size"] is not None

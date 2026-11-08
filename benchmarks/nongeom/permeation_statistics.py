"""Calibration of the permeation estimators: D_n (pf), the pd interval and the ion current.

Software verification against known answers (synthetic data only):

A. ``dn_estimator.csv``: :func:`crevice.permeation.collective_diffusion` on
   collective coordinates n(t) = random walk with known D_n, optionally plus
   a bounded Ornstein-Uhlenbeck term (fast local rattling of the water column
   with correlation time tau), for trajectories of 201 to 20001 frames.
   Reported: mean estimate / true value, relative SD, coverage of the
   reported interval (5 segments, Student t), how often the lower bound is
   clipped at zero.
B. ``dn_end_to_end.csv``: the same through :func:`analyze_permeation` on
   independent walkers on a periodic line (exact D_n = N D / Lambda^2),
   replicated, to check the interval of the whole pipeline.
C. ``pd_interval.csv``: coverage of the interval that CREVICE reports for
   pd (and q0) before (exact Poisson, Garwood) and after (Garwood widened by
   the dispersion index of the counts in 10 time blocks,
   :func:`crevice.permeation.dispersion_count_interval`), for independent,
   single-file and gated walkers. Counts come from the independent reference
   counter (CREVICE's counts are identical event by event; see
   permeation_counting.py); the true rate is the pooled rate over all
   replicates.
D. ``ion_current.csv``: cations and anions drifting in opposite directions on
   a periodic line (drift v, so the exact mean current is
   I = e (N+ + N-) v / Lambda), OmpF-like length scales. CREVICE's
   count-based current and a current from the ion collective coordinate
   (I_n = e sum q n(T) / T, read from the per-frame output) are compared with
   the exact value at saved-frame intervals of 1-500 ps.
E. ``ion_current_coverage.csv``: coverage of the Poisson and block-bootstrap
   intervals (before) and of the dispersion-widened interval (after) for
   constant drift and for drift switched by a common gate.

Run::

    PYTHONPATH=src python benchmarks/nongeom/permeation_statistics.py [--quick]
"""
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import numpy as np

from nongeom_common import (RESULTS, dump, ensure_src_on_path, expected_window_fit, machine_record, memory_universe,
                            ou_noise,
                            reference_axial_lateral, reference_events, reference_states, ring_walkers,
                            single_file_relabel, summarize_events, write_csv)

ensure_src_on_path()
from crevice.permeation import (CylinderChannel, analyze_permeation, collective_diffusion,  # noqa: E402
                                dispersion_count_interval)

N_PROTEIN = 6
E_PER_PS_IN_PA = 1.602176634e-19 / 1e-12 * 1e12   # e / ps in pA


def part_a(quick, rng):
    rows = []
    D = 0.5   # n^2 per frame (dt = 1)
    lengths = (201, 1001) if quick else (201, 1001, 5001, 20001)
    noises = ((0.0, 0.0), (1.0, 0.5), (1.0, 5.0)) if quick else \
        ((0.0, 0.0), (1.0, 0.5), (1.0, 2.0), (1.0, 5.0), (1.0, 20.0), (3.0, 5.0))
    reps = 100 if quick else 1000
    for frames in lengths:
        for sd, tau in noises:
            est, cover, clipped, status = [], 0, 0, 0
            for _ in range(reps):
                walk = np.concatenate([[0.0], np.cumsum(rng.normal(0, math.sqrt(2 * D), frames - 1))])
                n = walk + (ou_noise(frames, tau_frames=tau, sd=sd, rng=rng) if sd > 0 else 0.0)
                rec = collective_diffusion(n, 1.0, blocks=5)
                if rec["status"] != "estimated":
                    status += 1
                    continue
                est.append(rec["D_per_ps"])
                ci = rec.get("interval_per_s")
                if ci:
                    cover += ci[0] * 1e-12 <= D <= ci[1] * 1e-12
                    clipped += bool(rec.get("interval_clipped_at_zero"))
            est = np.asarray(est)
            rows.append({"frames": frames, "noise_sd": sd, "noise_tau_frames": tau, "replicates": reps,
                         "fit_lags_frames": rec["fit_lags_frames"], "estimated": len(est),
                         "mean_over_true": float(est.mean() / D) if len(est) else None,
                         "relative_sd": float(est.std(ddof=1) / D) if len(est) > 1 else None,
                         "interval_coverage": cover / len(est) if len(est) else None,
                         "lower_bound_clipped_fraction": clipped / len(est) if len(est) else None,
                         "not_estimated": status})
            print(f"Dn frames={frames} sd={sd} tau={tau}: mean/true={rows[-1]['mean_over_true']:.3f} "
                  f"rsd={rows[-1]['relative_sd']:.3f} cover={rows[-1]['interval_coverage']:.3f}", flush=True)
    return rows


def ring_universe(z, box_z, resnames):
    frames, n = z.shape
    pos = np.zeros((frames, n, 3))
    pos[:, :, 2] = z
    return memory_universe(resnames, pos, box=(30.0, 30.0, box_z), dt=1.0)


def part_b(quick, rng):
    box_z, half, n = 20.0, 5.0, 100
    frames = 2000 if quick else 5000
    reps = 5 if quick else 40
    sigma = 1.0
    D = sigma ** 2 / 2
    true = n * D / box_z ** 2
    rows = []
    for rep in range(reps):
        z = ring_walkers(n, frames, D=D, dt=1.0, box_z=box_z, rng=rng)
        u = ring_universe(z, box_z, ["TIP3"] * n)
        ch = CylinderChannel("pore", (0, 0, 1), 5.0, 15.0, point=(15.0, 15.0, 10.0), radius=3.0)
        res = analyze_permeation(u, species=["TIP3"], channels=[ch], align=False, msd_blocks=5)
        w = res.summary["channels"][0]["species"][0]["water"]
        ci = w["pf_interval_cm3_per_s"]
        v = w["v_w_cm3"]
        lo_lag, hi_lag = w["collective_diffusion"]["fit_lags_frames"]
        window = expected_window_fit(lo_lag, hi_lag, 1.0, n, D, 2 * half, box_z)
        rows.append({"replicate": rep, "frames": frames, "molecules": n, "D_n_true_per_ps": true,
                     "D_n_per_ps": w["D_n_per_s"] * 1e-12, "D_n_over_true": w["D_n_per_s"] * 1e-12 / true,
                     "pf_cm3_per_s": w["pf_cm3_per_s"], "pf_true_cm3_per_s": v * true * 1e12,
                     "pf_interval_covers_true": bool(ci and ci[0] <= v * true * 1e12 <= ci[1]),
                     "fit_lags_frames": [lo_lag, hi_lag], "D_n_window_expected_per_ps": window,
                     "D_n_over_window_expected": w["D_n_per_s"] * 1e-12 / window,
                     "pf_interval_covers_window_expected": bool(ci and ci[0] <= v * window * 1e12 <= ci[1])})
    print(f"Dn end-to-end: mean/true={np.mean([r['D_n_over_true'] for r in rows]):.3f} "
          f"mean/window={np.mean([r['D_n_over_window_expected'] for r in rows]):.3f} "
          f"cover(true)={np.mean([r['pf_interval_covers_true'] for r in rows]):.3f} "
          f"cover(window)={np.mean([r['pf_interval_covers_window_expected'] for r in rows]):.3f}", flush=True)
    return rows


def event_series(z, box_z, zmin, zmax):
    """Per-frame complete events (reference counter, compartment convention) for walkers on the axis."""
    s = z + box_z / 2
    img = (s - box_z / 2) - box_z * np.round((s - box_z / 2) / box_z) + box_z / 2
    states = reference_states(img, np.zeros_like(img), zmin, zmax, 3.0)
    x = np.zeros(len(z))
    for ev in reference_events(states, s, zmin, zmax, convention="compartment", s_image=img):
        if ev["status"] == "complete":
            x[ev["exit_frame"]] += 1
    return x


def gate(frames, mean_dwell, rng):
    """Two-state open (1) / closed (0) gate, stationary open probability 0.5, mean dwell in frames."""
    g = np.empty(frames)
    g[0] = rng.random() < 0.5
    flips = rng.random(frames) < 1.0 / mean_dwell
    for k in range(1, frames):
        g[k] = 1 - g[k - 1] if flips[k] else g[k - 1]
    return g


def gated_walkers(n, frames, *, D, box_z, rng, mean_dwell, drift=0.0):
    """Walkers whose steps (diffusion and drift) are switched off while a common gate is closed."""
    g = gate(frames - 1, mean_dwell, rng)
    x0 = rng.uniform(-box_z / 2, box_z / 2, n)
    steps = rng.normal(drift, math.sqrt(2 * D), size=(frames - 1, n)) * g[:, None]
    return np.vstack([x0[None], x0[None] + np.cumsum(steps, axis=0)]), g


def part_c(quick, rng):
    """Coverage of the Garwood interval (before) and of the dispersion-widened interval (after)."""
    box_z, half = 20.0, 5.0
    sigma = 0.5
    D = sigma ** 2 / 2
    reps = 30 if quick else 300
    rows = []
    models = (("independent", 20, 2000 if quick else 4000), ("single_file", 16, 4000 if quick else 20000),
              ("gated_fast_dwell_T/100", 20, 4000 if quick else 20000), ("gated_slow_dwell_T/20", 20, 4000 if quick else 20000))
    for model, n, frames in models:
        series = []
        for _ in range(reps):
            if model.startswith("gated"):
                dwell = frames / (100 if "fast" in model else 20)
                z, _g = gated_walkers(n, frames, D=D, box_z=box_z, rng=rng, mean_dwell=dwell)
            else:
                x = ring_walkers(n, frames, D=D, dt=1.0, box_z=box_z, rng=rng)
                z = single_file_relabel(x, box_z) if model == "single_file" else x
            series.append(event_series(z, box_z, box_z / 2 - half, box_z / 2 + half))
        totals = np.asarray([s.sum() for s in series])
        lam = totals.mean()   # pooled expected total count per replicate
        recs = [dispersion_count_interval(s) for s in series]
        cover_g = np.mean([r["garwood_interval"][0] <= lam <= r["garwood_interval"][1] for r in recs])
        cover_d = np.mean([r["interval"][0] <= lam <= r["interval"][1] for r in recs])
        width = np.mean([(r["interval"][1] - r["interval"][0]) / (r["garwood_interval"][1] - r["garwood_interval"][0])
                         for r in recs])
        phis = [r["dispersion_index"] for r in recs if r["dispersion_index"] is not None]
        rows.append({"model": model, "molecules": n, "frames": frames, "replicates": reps,
                     "mean_events_per_replicate": float(lam), "variance_of_events": float(totals.var(ddof=1)),
                     "index_of_dispersion_across_replicates": float(totals.var(ddof=1) / lam) if lam else None,
                     "mean_block_dispersion_index": float(np.mean(phis)) if phis else None,
                     "garwood_95_coverage_before": float(cover_g), "dispersion_95_coverage_after": float(cover_d),
                     "mean_width_after_over_before": float(width),
                     "note": "q0 and pd intervals are this count interval scaled; truth = pooled mean; "
                             "interval from crevice.permeation.dispersion_count_interval (10 blocks)"})
        print(f"pd interval {model}: mean={lam:.1f} dispersion={rows[-1]['index_of_dispersion_across_replicates']:.2f} "
              f"block phi={rows[-1]['mean_block_dispersion_index']:.2f} coverage before={cover_g:.3f} "
              f"after={cover_d:.3f} width ratio={width:.2f}", flush=True)
    return rows


def part_d(quick, rng):
    box_z, half = 100.0, 15.0              # OmpF-like: 30 A channel extent in a ~100 A periodic box
    D = 0.2                                 # A^2/ps, near bulk K+ / Cl-
    v = D * 0.005 / 0.0267                  # A/ps: drift D qE/kT for E = 0.05 V/nm at 310 K (illustrative)
    n_each = 20
    base_frames = 20000 if quick else 100000  # 1 ps steps: 20 or 100 ns
    rng_state = rng
    zp = ring_walkers(n_each, base_frames, D=D, dt=1.0, box_z=box_z, drift=+v, rng=rng_state)
    zm = ring_walkers(n_each, base_frames, D=D, dt=1.0, box_z=box_z, drift=-v, rng=rng_state)
    z = np.hstack([zp, zm])
    u = ring_universe(z, box_z, ["SOD"] * n_each + ["CLA"] * n_each)
    c = box_z / 2
    ch = CylinderChannel("pore", (0, 0, 1), c - half, c + half, point=(15.0, 15.0, c), radius=3.0)
    exact_pA = E_PER_PS_IN_PA * (2 * n_each) * v / box_z
    rows = []
    for stride in ((1, 100) if quick else (1, 10, 50, 100, 200, 500)):
        res = analyze_permeation(u, species=["SOD,CLA"], channels=[ch], align=False, stride=stride, msd_blocks=5)
        tot = res.summary["channels"][0]["ion_current_total"]
        T_ps = res.summary["frames"]["duration_ps"]
        series = res.frames["series"]
        n_sod = series["pore/SOD"]["n"][-1]
        n_cla = series["pore/CLA"]["n"][-1]
        i_n = E_PER_PS_IN_PA * (n_sod - n_cla) / T_ps
        unresolved = sum(r["unresolved_crossings"] for r in res.summary["channels"][0]["species"])
        events = sum(r["events_total"] for r in res.summary["channels"][0]["species"])
        p_lo, p_hi = tot["poisson"]["interval_pA"]
        boot = tot["block_bootstrap"].get("interval_pA")
        rows.append({"saved_interval_ps": stride, "sigma_per_saved_frame_A": math.sqrt(2 * D * stride),
                     "sigma_over_L": math.sqrt(2 * D * stride) / (2 * half), "frames_analysed": res.summary["frames"]["analysed"],
                     "exact_current_pA": exact_pA, "crevice_current_pA": tot["current_pA"],
                     "crevice_over_exact": tot["current_pA"] / exact_pA,
                     "collective_current_pA": i_n, "collective_over_exact": i_n / exact_pA,
                     "complete_events": events, "unresolved_crossings": unresolved,
                     "unresolved_fraction": unresolved / (events + unresolved) if events + unresolved else None,
                     "poisson_interval_pA": [p_lo, p_hi], "poisson_covers_exact": p_lo <= exact_pA <= p_hi,
                     "bootstrap_interval_pA": boot, "bootstrap_covers_exact": bool(boot and boot[0] <= exact_pA <= boot[1])})
        print(f"current stride={stride}: crevice/exact={rows[-1]['crevice_over_exact']:.3f} "
              f"collective/exact={rows[-1]['collective_over_exact']:.3f} unresolved={unresolved}/{events + unresolved}",
              flush=True)
    return rows


def part_d_coverage(quick, rng):
    """Interval coverage of the current at a 10 ps saving interval, over replicate trajectories.

    'drift': constant drift (events close to independent). 'gated_drift_*': the drift and diffusion of all ions
    are switched on and off by a common two-state gate (mean dwell 1/100 or 1/20 of the run, i.e. a tenth or half
    of a dispersion block), so the net charge comes in bursts; the expected current is half the open-gate current.
    """
    box_z, half, D = 100.0, 15.0, 0.2
    v = D * 0.005 / 0.0267
    n_each, frames, step = 20, (2000 if quick else 10000), 10.0   # 10 ps per saved frame
    reps = 5 if quick else 40
    rows = []
    for model in ("drift", "gated_drift_dwell_T/100", "gated_drift_dwell_T/20"):
        open_pA = E_PER_PS_IN_PA * (2 * n_each) * v / box_z
        exact_pA = open_pA * (0.5 if model.startswith("gated") else 1.0)
        cover = {"poisson": 0, "bootstrap": 0, "dispersion": 0}
        ratios, phis = [], []
        for _ in range(reps):
            if model.startswith("gated"):
                dwell = frames / (100 if model.endswith("T/100") else 20)
                zp, g = gated_walkers(n_each, frames, D=D * step, box_z=box_z, rng=rng, mean_dwell=dwell,
                                      drift=+v * step)
                zm = np.empty_like(zp)
                x0 = rng.uniform(-box_z / 2, box_z / 2, n_each)
                stepsm = rng.normal(-v * step, math.sqrt(2 * D * step), size=(frames - 1, n_each)) * g[:, None]
                zm[0] = x0
                zm[1:] = x0 + np.cumsum(stepsm, axis=0)
            else:
                zp = ring_walkers(n_each, frames, D=D, dt=step, box_z=box_z, drift=+v, rng=rng)
                zm = ring_walkers(n_each, frames, D=D, dt=step, box_z=box_z, drift=-v, rng=rng)
            pos = np.zeros((frames, 2 * n_each, 3))
            pos[:, :, 2] = np.hstack([zp, zm])
            u = memory_universe(["SOD"] * n_each + ["CLA"] * n_each, pos, box=(30.0, 30.0, box_z), dt=step)
            c = box_z / 2
            ch = CylinderChannel("pore", (0, 0, 1), c - half, c + half, point=(15.0, 15.0, c), radius=3.0)
            tot = analyze_permeation(u, species=["SOD,CLA"], channels=[ch], align=False)
            tot = tot.summary["channels"][0]["ion_current_total"]
            lo, hi = tot["poisson"]["interval_pA"]
            cover["poisson"] += lo <= exact_pA <= hi
            b = tot["block_bootstrap"].get("interval_pA")
            cover["bootstrap"] += bool(b and b[0] <= exact_pA <= b[1])
            lo, hi = tot["interval_pA"]
            cover["dispersion"] += lo <= exact_pA <= hi
            ratios.append(tot["current_pA"] / exact_pA)
            if tot.get("dispersion_index") is not None:
                phis.append(tot["dispersion_index"])
        row = {"model": model, "saved_interval_ps": step, "frames": frames, "replicates": reps,
               "exact_current_pA": exact_pA, "mean_current_over_exact": float(np.mean(ratios)),
               "sd_current_over_exact": float(np.std(ratios, ddof=1)),
               "mean_net_charge_dispersion_index": float(np.mean(phis)) if phis else None,
               "poisson_95_coverage_before": cover["poisson"] / reps,
               "block_bootstrap_95_coverage_before": cover["bootstrap"] / reps,
               "dispersion_95_coverage_after": cover["dispersion"] / reps}
        print(f"current coverage: {row}", flush=True)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--results", default=str(RESULTS))
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--parts", default="abcde", help="subset of parts a-e to run")
    args = parser.parse_args()
    out = Path(args.results)
    rng = np.random.default_rng(args.seed)
    tag = "_quick" if args.quick else ""
    t0 = time.perf_counter()
    summary = {"script": "benchmarks/nongeom/permeation_statistics.py", "quick": args.quick, "seed": args.seed,
               "machine": machine_record()}
    if "a" in args.parts:
        rows = part_a(args.quick, rng)
        write_csv(out / f"dn_estimator{tag}.csv", rows, comment="collective_diffusion on random-walk n(t) (+ OU "
                  "noise) with known D_n; dt = 1 frame, D = 0.5 per frame. Synthetic.")
        summary["dn_estimator"] = rows
    if "b" in args.parts:
        rows = part_b(args.quick, rng)
        write_csv(out / f"dn_end_to_end{tag}.csv", rows, comment="analyze_permeation on independent walkers; exact "
                  "D_n = N D / Lambda^2. Synthetic.")
        summary["dn_end_to_end"] = {"replicates": len(rows),
                                    "mean_D_n_over_true": float(np.mean([r["D_n_over_true"] for r in rows])),
                                    "sd_D_n_over_true": float(np.std([r["D_n_over_true"] for r in rows], ddof=1)),
                                    "pf_interval_coverage": float(np.mean([r["pf_interval_covers_true"] for r in rows])),
                                    "mean_D_n_over_window_expected": float(np.mean([r["D_n_over_window_expected"] for r in rows])),
                                    "pf_interval_coverage_of_window_expected":
                                        float(np.mean([r["pf_interval_covers_window_expected"] for r in rows]))}
    if "c" in args.parts:
        rows = part_c(args.quick, rng)
        write_csv(out / f"pd_interval{tag}.csv", rows, comment="Coverage of the Garwood interval used for q0/pd; "
                  "truth = pooled mean count. Synthetic.")
        summary["pd_interval"] = rows
    if "d" in args.parts:
        rows = part_d(args.quick, rng)
        write_csv(out / f"ion_current{tag}.csv", rows, comment="Drift-diffusion of cations and anions on a 100 A "
                  "period through a 30 A channel; exact I = e (N+ + N-) v / Lambda. Synthetic.")
        summary["ion_current"] = rows
    if "e" in args.parts:
        rows = part_d_coverage(args.quick, rng)
        write_csv(out / f"ion_current_coverage{tag}.csv", rows, comment="Replicate coverage of the current "
                  "intervals at a 10 ps saving interval. Synthetic.")
        summary["ion_current_coverage"] = rows
    summary["wall_seconds"] = time.perf_counter() - t0
    dump(out / f"permeation_statistics_summary{tag}{'' if args.parts == 'abcde' else '_' + args.parts}.json", summary)
    print(f"done in {summary['wall_seconds']:.0f} s")


if __name__ == "__main__":
    main()

"""Coverage of CREVICE's correlation-aware intervals on series with known autocorrelation.

Every per-frame summary of a trajectory is described by
:func:`crevice.cavity_trajectory.describe_series`, which reports

* ``confidence_interval``: the studentized batch bootstrap of
  :func:`crevice.uncertainty.block_mean_confidence` (primary interval), with
  the Geyer-sequence statistical inefficiency g;
* ``sampling_uncertainty``: tau_int (Sokal window), the effective sample size
  N / tau_int, the Flyvbjerg-Petersen blocking standard error and a
  moving-block bootstrap percentile interval.

This script draws stationary series with an exactly known integrated
autocorrelation time and standard error of the mean:

* Gaussian AR(1), ``x_t = phi x_(t-1) + e_t`` (tau = (1 + phi) / (1 - phi));
* a two-state Markov chain (binary occupancy, as for a water contact or a
  lining flag) with stationary probability p and lag-one correlation lambda,
  whose autocorrelation is lambda^t, so tau and the SE follow the same formulas.

For each (process, N) cell it reports how often each interval contains the
true mean (0, or p), how often an interval is withheld, and the bias of the
tau, g and SE estimates, as a function of N / tau.

Independent cross-checks (``--reference-python``, an interpreter with pymbar
and emcee, for example ``.venv-nongeom/bin/python``): the same series are
regenerated from their seeds in a separate process and analysed with
``pymbar.timeseries.statistical_inefficiency`` (Chodera et al., J. Chem.
Theory Comput. 3, 26 (2007)), ``emcee.autocorr.integrated_time`` (Sokal
window, c = 5) and a transcription of Jonsson's reference blocking code
(Phys. Rev. E 98, 043304 (2018); power-of-two series only). Results:
``uncertainty_coverage.csv`` and ``uncertainty_reference.csv``.

Run::

    PYTHONPATH=src python benchmarks/nongeom/uncertainty_coverage.py \\
        --reference-python .venv-nongeom/bin/python [--quick]
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from nongeom_common import NONGEOM, RESULTS, ar1, ar1_true_se, ar1_true_tau, dump, ensure_src_on_path, machine_record, \
    write_csv

ensure_src_on_path()


def markov_binary(n, p, lam, rng):
    """Two-state chain with P(x=1) = p and autocorrelation lam^t."""
    a = p * (1 - lam)          # P(0 -> 1)
    b = (1 - p) * (1 - lam)    # P(1 -> 0)
    x = np.empty(n)
    x[0] = rng.random() < p
    u = rng.random(n)
    for t in range(1, n):
        x[t] = (u[t] < 1 - b) if x[t - 1] else (u[t] < a)
    return x


def series_for(cell, rep, seed):
    """Deterministic series for a cell (also used by the reference worker)."""
    rng = np.random.default_rng([seed, cell["index"], rep])
    if cell["process"] == "ar1":
        return ar1(cell["n"], cell["phi"], rng)
    return markov_binary(cell["n"], cell["p"], cell["phi"], rng)


def cells(quick):
    out = []
    lengths = (200, 1000) if quick else (200, 1000, 5000)
    phis = (0.0, 0.9) if quick else (0.0, 0.5, 0.8, 0.9, 0.95, 0.98)
    for n in lengths:
        for phi in phis:
            out.append({"process": "ar1", "n": n, "phi": phi, "p": None})
        for lam in ((0.9,) if quick else (0.5, 0.9, 0.98)):
            out.append({"process": "markov_binary", "n": n, "phi": lam, "p": 0.2})
    for k, c in enumerate(out):
        c["index"] = k
        c["true_mean"] = 0.0 if c["process"] == "ar1" else c["p"]
        sd = 1.0 if c["process"] == "ar1" else math.sqrt(c["p"] * (1 - c["p"]))
        c["true_tau"] = ar1_true_tau(c["phi"])
        c["true_se"] = ar1_true_se(c["n"], c["phi"], sd)
    return out


def analyse_cell(cell, reps, seed):
    from crevice.cavity_trajectory import describe_series
    from scipy.stats import norm, t as student
    z = norm.ppf(0.975)
    mu = cell["true_mean"]
    acc = {k: [] for k in ("tau", "g", "ess", "blocking_se", "mbb_se", "ci_cover", "mbb_cover", "blocking_cover",
                           "ess_cover", "ci_estimated", "mbb_estimated", "ci_width", "mbb_width", "window_converged",
                           "blocking_converged")}
    for rep in range(reps):
        x = series_for(cell, rep, seed)
        d = describe_series(x)
        su = d["sampling_uncertainty"]
        ci = d["confidence_interval"]
        mean = d["mean"]
        if su["status"] == "estimated":
            acc["tau"].append(su["integrated_autocorrelation_time_frames"])
            acc["ess"].append(su["effective_sample_size"])
            acc["blocking_se"].append(su["blocking_standard_error"])
            acc["window_converged"].append(su["autocorrelation_window_converged"])
            acc["blocking_converged"].append(su["blocking_converged"])
            se_b = su["blocking_standard_error"]
            acc["blocking_cover"].append(abs(mean - mu) <= z * se_b)
            ess = su["effective_sample_size"]
            se_e = d["sample_sd"] / math.sqrt(ess)
            acc["ess_cover"].append(abs(mean - mu) <= student.ppf(0.975, max(ess - 1, 1)) * se_e)
            ok = su["block_bootstrap_status"] == "estimated"
            acc["mbb_estimated"].append(ok)
            if ok:
                acc["mbb_se"].append(su["block_bootstrap_standard_error"])
                acc["mbb_cover"].append(su["block_bootstrap_lower"] <= mu <= su["block_bootstrap_upper"])
                acc["mbb_width"].append(su["block_bootstrap_upper"] - su["block_bootstrap_lower"])
        if d.get("statistical_inefficiency") is not None:
            acc["g"].append(d["statistical_inefficiency"])
        ok = isinstance(ci, dict) and ci.get("status") == "estimated" and ci["pointwise_lower"][0] is not None
        acc["ci_estimated"].append(ok)
        if ok:
            acc["ci_cover"].append(ci["pointwise_lower"][0] <= mu <= ci["pointwise_upper"][0])
            acc["ci_width"].append(ci["pointwise_upper"][0] - ci["pointwise_lower"][0])
    mean_of = lambda k: float(np.mean(acc[k])) if acc[k] else None  # noqa: E731
    return {
        "process": cell["process"], "n_frames": cell["n"], "phi_or_lambda": cell["phi"], "p": cell["p"],
        "replicates": reps, "true_tau": cell["true_tau"], "n_over_true_tau": cell["n"] / cell["true_tau"],
        "mean_tau_int": mean_of("tau"), "mean_g_geyer": mean_of("g"), "mean_ess": mean_of("ess"),
        "true_se": cell["true_se"], "mean_blocking_se": mean_of("blocking_se"), "mean_mbb_se": mean_of("mbb_se"),
        "batch_ci_estimated_fraction": mean_of("ci_estimated"), "batch_ci_coverage": mean_of("ci_cover"),
        "batch_ci_mean_width_over_true": (mean_of("ci_width") / (2 * 1.959964 * cell["true_se"])) if acc["ci_width"] else None,
        "mbb_estimated_fraction": mean_of("mbb_estimated"), "mbb_coverage": mean_of("mbb_cover"),
        "blocking_1p96se_coverage": mean_of("blocking_cover"), "ess_t_coverage": mean_of("ess_cover"),
        "window_converged_fraction": mean_of("window_converged"), "blocking_converged_fraction": mean_of("blocking_converged"),
    }


def run_reference(python, cell_list, reps, seed, extra_n):
    """Run the pymbar/emcee/Jonsson worker in a separate interpreter; returns its JSON."""
    worker = NONGEOM / "uncertainty_reference_worker.py"
    payload = json.dumps({"cells": cell_list, "reps": reps, "seed": seed, "blocking_n": extra_n})
    proc = subprocess.run([python, "-I", str(worker)], input=payload, capture_output=True, text=True, check=True,
                          env={"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"})
    return json.loads(proc.stdout)


def crevice_estimates(cell, reps, seed):
    from crevice.uncertainty import blocking_standard_error, integrated_autocorrelation_time, statistical_inefficiency
    out = {"tau": [], "g": [], "blocking_se": []}
    for rep in range(reps):
        x = series_for(cell, rep, seed)
        out["tau"].append(integrated_autocorrelation_time(x)["tau_int_frames"])
        out["g"].append(statistical_inefficiency(x))
        out["blocking_se"].append(blocking_standard_error(x)["standard_error"])
    return {k: np.asarray(v) for k, v in out.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--reference-python", default=None,
                        help="interpreter with pymbar and emcee (cross-check); omitted = skipped")
    parser.add_argument("--results", default=str(RESULTS))
    parser.add_argument("--seed", type=int, default=20261027)
    args = parser.parse_args()
    out = Path(args.results)
    tag = "_quick" if args.quick else ""
    reps = args.replicates or (40 if args.quick else 400)
    t0 = time.perf_counter()
    cell_list = cells(args.quick)
    rows = []
    for cell in cell_list:
        row = analyse_cell(cell, reps, args.seed)
        rows.append(row)
        print(f"{row['process']} N={row['n_frames']} phi={row['phi_or_lambda']}: N/tau={row['n_over_true_tau']:.0f} "
              f"tau={row['mean_tau_int']} batchCI cov={row['batch_ci_coverage']} (est {row['batch_ci_estimated_fraction']}) "
              f"MBB cov={row['mbb_coverage']} (est {row['mbb_estimated_fraction']}) block cov={row['blocking_1p96se_coverage']}",
              flush=True)
    write_csv(out / f"uncertainty_coverage{tag}.csv", rows,
              comment="Coverage of 95% intervals for the mean on stationary series with known autocorrelation "
                      "(describe_series: batch_ci = studentized batch bootstrap, mbb = moving-block bootstrap, blocking = "
                      "mean +/- 1.96 x Flyvbjerg-Petersen SE, ess_t = mean +/- t x SD/sqrt(N/tau_int)). Coverage is "
                      "conditional on an interval being reported. Synthetic.")
    summary = {"script": "benchmarks/nongeom/uncertainty_coverage.py", "quick": args.quick, "seed": args.seed,
               "replicates": reps, "machine": machine_record(), "coverage": rows}
    if args.reference_python:
        ref_reps = min(reps, 100)
        blocking_n = [256, 1024] if args.quick else [1024, 4096]
        ref = run_reference(args.reference_python, cell_list, ref_reps, args.seed, blocking_n)
        ref_rows = []
        for cell in cell_list:
            mine = crevice_estimates(cell, ref_reps, args.seed)
            theirs = ref["cells"][str(cell["index"])]
            emcee_tau = np.asarray(theirs["emcee_tau"], float)
            pymbar_g = np.asarray(theirs["pymbar_g"], float)
            ref_rows.append({
                "process": cell["process"], "n_frames": cell["n"], "phi_or_lambda": cell["phi"], "replicates": ref_reps,
                "true_tau": cell["true_tau"],
                "crevice_tau_mean": float(mine["tau"].mean()), "emcee_tau_mean": float(np.nanmean(emcee_tau)),
                "tau_max_abs_rel_diff_crevice_vs_emcee": float(np.nanmax(np.abs(mine["tau"] - emcee_tau) / emcee_tau)),
                # CREVICE clips tau_int to [1, N]; emcee does not (white noise can give tau < 1).
                "tau_max_abs_rel_diff_vs_emcee_clipped_at_1": float(np.nanmax(
                    np.abs(mine["tau"] - np.maximum(emcee_tau, 1.0)) / np.maximum(emcee_tau, 1.0))),
                "emcee_tau_below_1_fraction": float(np.mean(emcee_tau < 1)),
                "crevice_g_geyer_mean": float(mine["g"].mean()), "pymbar_g_mean": float(np.nanmean(pymbar_g)),
                "pymbar_refused_constant_series": int(np.sum(np.isnan(pymbar_g))),
                "crevice_g_rmse_vs_true": float(np.sqrt(np.mean((mine["g"] - cell["true_tau"]) ** 2))),
                "pymbar_g_rmse_vs_true": float(np.sqrt(np.nanmean((pymbar_g - cell["true_tau"]) ** 2))),
                "g_max_abs_rel_diff_crevice_vs_pymbar": float(np.nanmax(np.abs(mine["g"] - pymbar_g) / pymbar_g)),
                "crevice_tau_rmse_vs_true": float(np.sqrt(np.mean((mine["tau"] - cell["true_tau"]) ** 2))),
            })
            print(f"ref {cell['process']} N={cell['n']} phi={cell['phi']}: tau crevice/emcee max rel diff "
                  f"{ref_rows[-1]['tau_max_abs_rel_diff_crevice_vs_emcee']:.2e} (clipped {ref_rows[-1]['tau_max_abs_rel_diff_vs_emcee_clipped_at_1']:.1e}); g crevice {ref_rows[-1]['crevice_g_geyer_mean']:.2f} "
                  f"pymbar {ref_rows[-1]['pymbar_g_mean']:.2f} true {cell['true_tau']:.2f}", flush=True)
        # Blocking: CREVICE vs a transcription of Jonsson's reference code on power-of-two series.
        from crevice.uncertainty import blocking_standard_error
        block_rows = []
        for item in ref["blocking"]:
            rng = np.random.default_rng([args.seed, 999, item["n"], int(item["phi"] * 1000), item["rep"]])
            x = ar1(item["n"], item["phi"], rng)
            mine = blocking_standard_error(x)["standard_error"]
            block_rows.append(abs(mine - item["jonsson_se"]) / item["jonsson_se"])
        write_csv(out / f"uncertainty_reference{tag}.csv", ref_rows,
                  comment="CREVICE tau_int / Geyer g vs emcee.autocorr.integrated_time and "
                          "pymbar.timeseries.statistical_inefficiency on identical series. Synthetic.")
        summary["reference"] = {"rows": ref_rows, "versions": ref["versions"],
                                "blocking_vs_jonsson_series": len(block_rows),
                                "blocking_vs_jonsson_max_rel_diff": float(max(block_rows)) if block_rows else None}
        print(f"blocking vs Jonsson reference: {len(block_rows)} series, max rel diff "
              f"{summary['reference']['blocking_vs_jonsson_max_rel_diff']:.2e}", flush=True)
    summary["wall_seconds"] = time.perf_counter() - t0
    dump(out / f"uncertainty_coverage_summary{tag}.json", summary)
    print(f"done in {summary['wall_seconds']:.0f} s")


if __name__ == "__main__":
    sys.exit(main())

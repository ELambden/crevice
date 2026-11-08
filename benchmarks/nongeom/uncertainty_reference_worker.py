"""Reference estimators for ``uncertainty_coverage.py`` (runs in a separate interpreter).

Reads a JSON request on stdin (cells, replicates, seed, power-of-two lengths
for the blocking comparison), regenerates the same series from their seeds and
writes JSON to stdout with, per series:

* ``pymbar.timeseries.statistical_inefficiency`` (fast=False, fft=False: the direct sum);
* ``emcee.autocorr.integrated_time`` (c=5, tol=0, quiet=True);
* the blocking standard error of a transcription of M. Jonsson's reference
  implementation (github.com/computative/block, python/tictoc.py; series of
  length 2^d), written here without CREVICE.

This file imports nothing from CREVICE. It needs numpy, scipy, pymbar, emcee.
"""
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


def ar1(n, phi, rng, sigma=1.0):
    e = rng.normal(0, sigma * math.sqrt(1 - phi ** 2), n)
    x = np.empty(n)
    x[0] = rng.normal(0, sigma)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x


def markov_binary(n, p, lam, rng):
    a = p * (1 - lam)
    b = (1 - p) * (1 - lam)
    x = np.empty(n)
    x[0] = rng.random() < p
    u = rng.random(n)
    for t in range(1, n):
        x[t] = (u[t] < 1 - b) if x[t - 1] else (u[t] < a)
    return x


def series_for(cell, rep, seed):
    rng = np.random.default_rng([seed, cell["index"], rep])
    if cell["process"] == "ar1":
        return ar1(cell["n"], cell["phi"], rng)
    return markov_binary(cell["n"], cell["p"], cell["phi"], rng)


def jonsson_block(x):
    """Transcription of Jonsson's tictoc.py ``block`` (length must be 2^d); returns the SE of the mean."""
    from scipy.stats import chi2
    n = len(x)
    d = int(math.log2(n))
    s, gamma = np.zeros(d), np.zeros(d)
    mu = np.mean(x)
    for i in range(d):
        n = len(x)
        gamma[i] = (n) ** (-1) * np.sum((x[0:(n - 1)] - mu) * (x[1:n] - mu))
        s[i] = np.var(x)
        x = 0.5 * (x[0::2] + x[1::2])
    M = (np.cumsum(((gamma / s) ** 2 * 2 ** np.arange(1, d + 1)[::-1])[::-1]))[::-1]
    q = chi2.ppf(0.99, np.arange(1, d + 1))
    k = d - 1
    for k in range(d):
        if M[k] < q[k]:
            break
    return math.sqrt(s[k] / 2 ** (d - k))


def main():
    warnings.filterwarnings("ignore")
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        import emcee
        import pymbar
        from pymbar import timeseries
    req = json.loads(sys.stdin.read())
    out = {"versions": {"pymbar": pymbar.__version__, "emcee": emcee.__version__, "numpy": np.__version__},
           "cells": {}, "blocking": []}
    for cell in req["cells"]:
        g, tau = [], []
        for rep in range(req["reps"]):
            x = series_for(cell, rep, req["seed"])
            with contextlib.redirect_stdout(io.StringIO()):
                try:   # pymbar refuses constant series (possible for short binary chains)
                    g.append(float(timeseries.statistical_inefficiency(x, fast=False, fft=False)))
                except Exception:
                    g.append(float("nan"))
                try:
                    tau.append(float(emcee.autocorr.integrated_time(x, c=5, tol=0, quiet=True)[0]))
                except Exception:
                    tau.append(float("nan"))
        out["cells"][str(cell["index"])] = {"pymbar_g": g, "emcee_tau": tau}
    for n in req["blocking_n"]:
        for phi in (0.0, 0.5, 0.9, 0.98):
            for rep in range(20):
                rng = np.random.default_rng([req["seed"], 999, n, int(phi * 1000), rep])
                out["blocking"].append({"n": n, "phi": phi, "rep": rep, "jonsson_se": jonsson_block(ar1(n, phi, rng))})
    sys.stdout.write(json.dumps(out))


if __name__ == "__main__":
    main()

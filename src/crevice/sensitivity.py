"""Definition sensitivity of headline numbers (``--sensitivity``).

A geometric measurement depends on its definition: the region margin, the
probe radius and the grid spacing. :func:`run_sensitivity` reruns a
measurement with one setting changed at a time and reports how far each
headline quantity moves from the baseline run:

* **margin** ±1 Å (region and cavity trajectories in reference-region mode);
* **probe** the neighbouring radii of the probe ladder
  (:data:`crevice.sections.AUTO_ENCLOSURE_LADDER`, 0.5–3.0 Å) around the
  probe actually used: the enclosure probe of a profile, the rolling probe of
  a cast or region;
* **grid** 0.25 and 0.5 Å (the cast or region lattice, or the profile's
  transverse section grid), skipping the baseline value.

The spread is descriptive: it shows how strongly a number depends on choices
that have no unique correct value. It is not a confidence interval and the
variants are not equally plausible alternatives. Each variant reruns the
whole measurement, so the cost is roughly (number of variants) x (baseline
run time); trajectory commands rerun only a few evenly spaced frames
(``--sensitivity-frames``, default 10). Nothing runs unless requested.

Main entry points: :func:`run_sensitivity`, :func:`ladder_neighbours`,
:func:`write_sensitivity`.
"""

from __future__ import annotations

import csv
import json
import math
import time
from pathlib import Path

#: Grid spacings (Å) of the grid variants.
GRID_VARIANTS_A = (0.25, 0.5)
#: Margin change (Å) of the margin variants.
MARGIN_DELTA_A = 1.0

#: Columns of :func:`write_sensitivity`'s CSV.
SENSITIVITY_COLUMNS = ("quantity", "variant", "parameter", "value", "baseline_value", "result", "baseline_result",
                       "difference", "relative_difference_percent", "status", "seconds")


def ladder_neighbours(value, ladder=None):
    """The ladder probes immediately below and above ``value``.

    Parameters
    ----------
    value : float
        The probe radius used (Å).
    ladder : sequence of float, optional
        Candidate probes; default :data:`crevice.sections.AUTO_ENCLOSURE_LADDER`.

    Returns
    -------
    tuple of float
        Zero, one or two neighbours, in increasing order.

    Examples
    --------
    >>> from crevice.sensitivity import ladder_neighbours
    >>> ladder_neighbours(0.8)
    (0.75, 0.9)
    >>> ladder_neighbours(0.5)
    (0.55,)
    """
    if ladder is None:
        from .sections import AUTO_ENCLOSURE_LADDER as ladder
    values = sorted(float(v) for v in ladder)
    below = [v for v in values if v < float(value) - 1e-9]
    above = [v for v in values if v > float(value) + 1e-9]
    return tuple(([below[-1]] if below else []) + ([above[0]] if above else []))


def standard_variants(*, margin=None, probe=None, grid=None, grid_variants=GRID_VARIANTS_A):
    """The variant list for the settings that apply to a command.

    Parameters
    ----------
    margin, probe, grid : float, optional
        Baseline margin, probe radius and grid spacing (Å); ``None`` when the
        setting does not apply.
    grid_variants : sequence of float
        Grid spacings to try (the baseline value is skipped).

    Returns
    -------
    list of tuple
        ``(parameter, value)`` pairs.

    Examples
    --------
    >>> from crevice.sensitivity import standard_variants
    >>> standard_variants(margin=2.0, probe=0.8, grid=0.5)
    [('margin_A', 1.0), ('margin_A', 3.0), ('probe_radius_A', 0.75), ('probe_radius_A', 0.9), ('grid_spacing_A', 0.25)]
    """
    variants = []
    if margin is not None:
        variants += [("margin_A", m) for m in (margin - MARGIN_DELTA_A, margin + MARGIN_DELTA_A) if m >= 0]
    if probe is not None:
        variants += [("probe_radius_A", p) for p in ladder_neighbours(probe)]
    if grid is not None:
        variants += [("grid_spacing_A", g) for g in grid_variants if not math.isclose(g, grid)]
    return variants


def run_sensitivity(compute, baseline, variants, *, command="", note=""):
    """Rerun ``compute`` with one setting changed at a time and tabulate the spread.

    Parameters
    ----------
    compute : callable
        ``compute(overrides)`` returns a dict of headline quantity to number
        (``None`` when not measured); ``overrides`` maps a parameter name of
        ``variants`` to its value (empty for the baseline).
    baseline : dict
        Baseline value of every varied parameter (recorded).
    variants : sequence of tuple
        ``(parameter, value)`` pairs, for example from :func:`standard_variants`.
    command : str, optional
        Command name, recorded.
    note : str, optional
        Scope note, recorded (for example the frame subset).

    Returns
    -------
    dict
        ``runs`` (one per variant, baseline first, with ``quantities``,
        ``status`` and ``seconds``), ``spread`` (per quantity: baseline,
        minimum, maximum, range and range relative to the baseline in
        percent, over the runs that succeeded), ``rows`` (long table for
        :func:`write_sensitivity`), ``cost_seconds``, ``baseline_settings``
        and ``interpretation``. A variant that fails is recorded with its
        error and left out of the spread.
    """
    runs = []
    for parameter, value in [(None, None), *variants]:
        start = time.perf_counter()
        overrides = {} if parameter is None else {parameter: value}
        try:
            quantities, status = dict(compute(overrides)), "ok"
        except Exception as exc:  # a variant may legitimately fail (for example no channel at that probe)
            quantities, status = {}, f"failed: {type(exc).__name__}: {exc}"
        runs.append({"variant": "baseline" if parameter is None else f"{parameter}={value:g}",
                     "parameter": parameter or "", "value": value, "quantities": quantities, "status": status,
                     "seconds": round(time.perf_counter() - start, 3)})
    base = runs[0]["quantities"]
    names = list(dict.fromkeys(name for run in runs for name in run["quantities"]))
    spread, rows = {}, []
    for name in names:
        values = [run["quantities"].get(name) for run in runs if run["status"] == "ok"]
        values = [float(v) for v in values if v is not None and math.isfinite(float(v))]
        reference = base.get(name)
        low, high = (min(values), max(values)) if values else (None, None)
        spread[name] = {"baseline": reference, "minimum": low, "maximum": high,
                        "range": high - low if values else None,
                        "relative_range_percent": (100 * (high - low) / abs(reference)
                                                   if values and reference not in (None, 0) else None),
                        "runs": len(values)}
        for run in runs:
            result = run["quantities"].get(name)
            difference = result - reference if result is not None and reference is not None else None
            rows.append({"quantity": name, "variant": run["variant"], "parameter": run["parameter"],
                         "value": run["value"], "baseline_value": baseline.get(run["parameter"]),
                         "result": result, "baseline_result": reference, "difference": difference,
                         "relative_difference_percent": (100 * difference / abs(reference)
                                                         if difference is not None and reference not in (None, 0) else None),
                         "status": run["status"], "seconds": run["seconds"]})
    return {"command": command, "baseline_settings": dict(baseline), "variants": [list(v) for v in variants],
            "runs": runs, "spread": spread, "rows": rows,
            "cost_seconds": round(sum(run["seconds"] for run in runs), 3), "note": note,
            "definition": "one setting changed at a time: margin +/-1 A, neighbouring probe-ladder radii, "
                          "grid 0.25/0.5 A (where the setting applies)",
            "interpretation": "spread of the headline numbers over measurement definitions; descriptive, "
                              "not a confidence interval and not a ranking of the variants"}


def write_sensitivity(result, csv_path, json_path=None):
    """Write the long table (CSV) and, optionally, the full record (JSON).

    Parameters
    ----------
    result : dict
        From :func:`run_sensitivity`.
    csv_path : str or Path
        Output CSV with :data:`SENSITIVITY_COLUMNS` (one row per quantity and variant).
    json_path : str or Path, optional
        Output JSON (the whole ``result``).

    Returns
    -------
    list of Path
        The files written.
    """
    written = []
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SENSITIVITY_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: ("" if row.get(k) is None else row.get(k)) for k in SENSITIVITY_COLUMNS}
                         for row in result["rows"])
    written.append(csv_path)
    if json_path is not None:
        json_path = Path(json_path)
        from ._portable_paths import dumps_record  # 0.4: recorded paths relative to this file
        json_path.write_text(dumps_record(result, json_path, indent=2, sort_keys=True, default=str) + "\n",
                             encoding="utf-8", newline="\n")
        written.append(json_path)
    return written


def describe_spread(result) -> str:
    """One line per headline quantity: baseline and the min–max over variants.

    Parameters
    ----------
    result : dict
        From :func:`run_sensitivity`.

    Returns
    -------
    str
        The printed summary.

    Examples
    --------
    >>> from crevice.sensitivity import describe_spread
    >>> print(describe_spread({"spread": {"v": {"baseline": 2.0, "minimum": 1.0, "maximum": 3.0,
    ...       "relative_range_percent": 100.0, "runs": 3}}, "cost_seconds": 1.5}))
    sensitivity v: baseline 2 range 1 to 3 (100% of baseline, 3 runs); cost 1.5 s
    """
    lines = []
    for name, s in result["spread"].items():
        if s["minimum"] is None:
            lines.append(f"sensitivity {name}: not measured")
            continue
        relative = f"{s['relative_range_percent']:.3g}% of baseline, " if s["relative_range_percent"] is not None else ""
        lines.append(f"sensitivity {name}: baseline {s['baseline']:.6g} range {s['minimum']:.6g} to "
                     f"{s['maximum']:.6g} ({relative}{s['runs']} runs)")
    if lines:
        lines[-1] += f"; cost {result['cost_seconds']:.3g} s"
    return "\n".join(lines)

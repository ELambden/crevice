"""Sensitivity of CREVICE's automatic enclosure probe to its grouping rule and ladder.

Usage (CREVICE environment; after ``curate_channels.py``)::

    python benchmarks/probe_sensitivity.py [--work-dir DIR] [--only 1GRM]

The automatic probe (``crevice profile --enclosure-radius auto``) resolves the
profile at every probe of a ladder (0.5-3.0 Å, 21 probes), groups resolved
results into channels (minimum radius within max(0.05 Å, 2 %), axial ranges
overlapping >= 60 % of their union, same axis and exit types) and uses the
smallest probe of the channel resolved at the most probes (at least 2).

This script reruns that rule on each benchmark channel with one setting
changed at a time:

* minimum-radius tolerance: absolute 0.025 / 0.1 Å (default 0.05); fraction
  1 % / 5 % (default 2 %);
* axial-overlap fraction 0.4 / 0.8 (default 0.6);
* minimum support 3 / 4 probes (default 2; applied to the recorded groups);
* ladder range 0.5-2.0, 0.5-4.0 (adds 3.3, 3.6, 4.0), 0.4-3.0 (adds 0.40,
  0.45), 0.6-3.0, and every other ladder probe.

It patches module constants of :mod:`crevice.sections` in this process only
(the package is not modified) and caches each explicit-probe profile, so a
variant reuses identical single-probe results. Settings otherwise as the
HOLE comparison (curated atoms, axis and seed; 0.25 Å section grid; 81
samples; Bondi radii, CREVICE's default for protein atoms), with
``--lateral-exits`` and a 12 Å search window where the comparison needed them.
Outputs ``results/probe_sensitivity/variants.csv`` (one row per channel and
variant) and ``results/probe_sensitivity/summary.csv``.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INPUTS, REPO, RESULTS, default_work_dir, read_csv, write_csv  # noqa: E402

sys.path.insert(0, str(REPO / "src"))
import crevice  # noqa: E402
from crevice import sections  # noqa: E402
from crevice.sections import ChannelResolutionError  # noqa: E402

DEFAULT_LADDER = tuple(sections.AUTO_ENCLOSURE_LADDER)
ORIGINAL = sections.connected_profile
_CACHE: dict = {}


def cached_connected_profile(atoms, **kwargs):
    if kwargs.get("enclosure_radius") is None:
        return ORIGINAL(atoms, **kwargs)
    # One structure at a time (the cache is cleared per channel), so atoms need not be in the key.
    key = (len(atoms), round(float(kwargs["enclosure_radius"]), 6), kwargs.get("lateral_exits"),
           kwargs.get("search_radius"), kwargs.get("samples"))
    if key not in _CACHE:
        try:
            _CACHE[key] = ("ok", ORIGINAL(atoms, **kwargs))
        except ChannelResolutionError as exc:
            _CACHE[key] = ("err", exc)
    status, value = _CACHE[key]
    if status == "err":
        raise value
    return value


sections.connected_profile = cached_connected_profile


def variants():
    d = DEFAULT_LADDER
    yield "default", {}, d, 2
    yield "min_radius_abs_0.025", {"_SAME_MIN_RADIUS_A": 0.025}, d, 2
    yield "min_radius_abs_0.1", {"_SAME_MIN_RADIUS_A": 0.1}, d, 2
    yield "min_radius_frac_0.01", {"_SAME_MIN_RADIUS_FRACTION": 0.01}, d, 2
    yield "min_radius_frac_0.05", {"_SAME_MIN_RADIUS_FRACTION": 0.05}, d, 2
    yield "overlap_0.4", {"_SAME_RANGE_FRACTION": 0.4}, d, 2
    yield "overlap_0.8", {"_SAME_RANGE_FRACTION": 0.8}, d, 2
    yield "min_support_3", {}, d, 3
    yield "min_support_4", {}, d, 4
    yield "ladder_0.5-2.0", {}, tuple(p for p in d if p <= 2.0), 2
    yield "ladder_0.5-4.0", {}, d + (3.3, 3.6, 4.0), 2
    yield "ladder_0.4-3.0", {}, (0.40, 0.45) + d, 2
    yield "ladder_0.6-3.0", {}, tuple(p for p in d if p >= 0.6), 2
    yield "ladder_every_other", {}, d[::2], 2


def apply_support(record: dict, k: int):
    """Re-apply the choice with a larger minimum support to the recorded groups."""
    groups = [g for g in record.get("groups", []) if g["support"] >= k]
    if not groups:
        return None, None
    best = max(groups, key=lambda g: (g["support"], -g["probes_A"][0]))
    return best["probes_A"][0], best["min_radius_A"]


def run_variant(frame, opts, consts, ladder, support):
    saved = {name: getattr(sections, name) for name in consts}
    saved_ladder = sections.AUTO_ENCLOSURE_LADDER
    try:
        for name, value in consts.items():
            setattr(sections, name, value)
        sections.AUTO_ENCLOSURE_LADDER = tuple(ladder)
        try:
            prof = crevice.pore_profile(frame, enclosure_radius=None, **opts)
        except ChannelResolutionError as exc:
            record = getattr(exc, "probe_selection", {}) or {}
            if support > 2:
                return {"status": "unresolved", "reason": "no channel at >= 2 probes"}
            return {"status": "unresolved", "reason": record.get("reason", str(exc)[:120]),
                    "groups": len(record.get("groups", []))}
        record = prof.metadata["enclosure_probe_selection"]
        if support > 2:
            chosen, min_r = apply_support(record, support)
            if chosen is None:
                return {"status": "unresolved", "reason": f"no channel at >= {support} probes"}
            return {"status": "resolved", "chosen_probe_A": chosen, "min_radius_A": min_r,
                    "support": next(g["support"] for g in record["groups"] if g["probes_A"][0] == chosen),
                    "groups": len(record["groups"])}
        group = next(g for g in record["groups"] if g["group"] == record["chosen_group"])
        return {"status": "resolved", "chosen_probe_A": record["chosen_A"], "min_radius_A": prof.min_radius,
                "support": group["support"], "groups": len(record["groups"])}
    finally:
        for name, value in saved.items():
            setattr(sections, name, value)
        sections.AUTO_ENCLOSURE_LADDER = saved_ladder


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args(argv)
    rows, summary = [], []
    for row in read_csv(INPUTS / "channels.csv"):
        acc = row["id"]
        if args.only and acc not in args.only:
            continue
        cur = args.work_dir / "channels" / acc
        axis = json.loads((cur / "axis.json").read_text())
        comp_path = cur / "comparison_bondi.json"
        used = json.loads(comp_path.read_text())["summary"].get("crevice_attempt_used") if comp_path.exists() else None
        opts = {"axis": tuple(axis["axis_unit"]), "origin": tuple(axis["origin"]), "samples": 81,
                "section_spacing": 0.25, "radii": "bondi"}
        if used in ("lateral_exits", "lateral_exits_search12"):
            opts["lateral_exits"] = True
        if used == "lateral_exits_search12":
            opts["search_radius"] = 12.0
        frame = crevice.load_structure(str(cur / "curated.pdb"))
        _CACHE.clear()
        base = None
        t0 = time.perf_counter()
        for name, consts, ladder, support in variants():
            res = run_variant(frame, opts, consts, ladder, support)
            if base is None:
                base = res
            changed_status = res["status"] != base["status"]
            changed_probe = res.get("chosen_probe_A") != base.get("chosen_probe_A")
            dmin = (res["min_radius_A"] - base["min_radius_A"]) if (res.get("min_radius_A") is not None
                                                                  and base.get("min_radius_A") is not None) else None
            rows.append({"id": acc, "options_used": used or "default", "variant": name, **res,
                         "status_changed": changed_status, "probe_changed": changed_probe,
                         "min_radius_change_A": dmin})
            print(acc, name, res)
        vr = [r for r in rows if r["id"] == acc and r["variant"] != "default"]
        summary.append({"id": acc, "options_used": used or "default", "default_status": base["status"],
                        "default_probe_A": base.get("chosen_probe_A"), "default_min_radius_A": base.get("min_radius_A"),
                        "variants": len(vr), "status_changed": sum(r["status_changed"] for r in vr),
                        "probe_changed": sum(r["probe_changed"] for r in vr),
                        "min_radius_changed_gt_0.01A": sum(1 for r in vr if r["min_radius_change_A"] is not None
                                                           and abs(r["min_radius_change_A"]) > 0.01),
                        "max_abs_min_radius_change_A": max((abs(r["min_radius_change_A"]) for r in vr
                                                            if r["min_radius_change_A"] is not None), default=None),
                        "wall_s": time.perf_counter() - t0})
    if not args.only:
        write_csv(args.results / "probe_sensitivity" / "variants.csv", rows,
                  comment="Automatic enclosure-probe choice under one-at-a-time changes of the grouping rule and "
                          "ladder (probe_sensitivity.py). Bondi radii; curated atoms, axis and seed.")
        write_csv(args.results / "probe_sensitivity" / "summary.csv", summary,
                  comment="Per channel: how many of the 13 variants change the status, the chosen probe or the "
                          "minimum radius (> 0.01 A).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

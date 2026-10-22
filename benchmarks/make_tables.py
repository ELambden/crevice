"""Collect the benchmark result CSVs into the summary tables used by the paper.

Usage::

    python benchmarks/make_tables.py [--results benchmarks/results]

Reads the per-comparison CSVs and writes

``validation_channels.csv``
    One row per channel and radius set: CREVICE status with default options
    and with the attempt used, enclosure probe, path type, RMSD / max |diff|
    against HOLE in the interior window, both minimum radii and positions,
    HOLE's global minimum and whether it lies inside CREVICE's resolved range
    (coverage), and the conductance ratio where a measurement exists.
``validation_summary.csv``
    One row per benchmark (channels vs HOLE, conductance, tunnels vs CAVER,
    volumes, trajectory pocket, probe sensitivity) with the headline numbers.

Only numbers already in the result CSVs are used; nothing is recomputed.
"""
from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS, read_csv, write_csv  # noqa: E402


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results", type=Path, default=RESULTS)
    args = ap.parse_args(argv)
    R = args.results
    rows, summary = [], []
    for tag in ("simple", "bondi"):
        path = R / "channels" / f"channels_vs_hole_{tag}.csv"
        if not path.exists():
            continue
        cond = {r["id"]: r for r in read_csv(R / "conductance" / f"conductance_vs_measured_{tag}.csv")} \
            if (R / "conductance" / f"conductance_vs_measured_{tag}.csv").exists() else {}
        chans = read_csv(path)
        for r in chans:
            lo, hi, gs = f(r.get("crevice_s_min_A")), f(r.get("crevice_s_max_A")), f(r.get("hole_min_s_global_A"))
            inside = (lo is not None and gs is not None and lo <= gs <= hi)
            c = cond.get(r["id"], {})
            rows.append({"id": r["id"], "name": r["name"], "state": r["state"], "curation": r["curation"],
                         "radius_set": tag, "crevice_status_default_options": r["crevice_status_default"],
                         "crevice_status": r["crevice_status"], "crevice_options_used": r["crevice_attempt_used"],
                         "enclosure_probe_A": r.get("crevice_chosen_probe_A"), "path_type": r.get("crevice_path_type"),
                         "compared_profile": r.get("compared_profile"), "n_compared": r.get("n_compared"),
                         "rmsd_A": r.get("rmsd_A"), "max_abs_diff_A": r.get("max_abs_diff_A"),
                         "crevice_min_radius_A": r.get("crevice_min_radius_A"), "crevice_min_s_A": r.get("crevice_min_s_A"),
                         "hole_min_in_range_A": r.get("hole_min_radius_in_range_A"),
                         "hole_min_in_range_s_A": r.get("hole_min_s_in_range_A"),
                         "crevice_range_A": f"{lo:.1f} to {hi:.1f}" if lo is not None else "",
                         "hole_global_min_A": r.get("hole_min_radius_global_A"), "hole_global_min_s_A": gs,
                         "hole_global_min_inside_crevice_range": inside if lo is not None else "",
                         "measured_pS": c.get("measured_pS_per_pore"), "crevice_g_macro_pS": c.get("crevice_g_macro_pS"),
                         "crevice_ratio": c.get("crevice_macro_ratio"), "hole_g_macro_full_pS": c.get("hole_g_macro_full_pS"),
                         "hole_full_ratio": c.get("hole_full_ratio"), "conductance_comparable": c.get("comparable")})
        res = [r for r in chans if r["crevice_status"] == "resolved" and f(r.get("rmsd_A")) is not None]
        dmin = [abs(f(r["crevice_min_radius_A"]) - f(r["hole_min_radius_in_range_A"])) for r in res]
        summary.append({"benchmark": f"channel profiles vs HOLE 2.3.1 ({tag}.rad in both)",
                        "entries": len(chans),
                        "result": (f"resolved with default options {sum(r['crevice_status_default'] == 'resolved' for r in chans)}"
                                   f"/{len(chans)}; with documented options (lateral exits, 12 A search) "
                                   f"{sum(r['crevice_status'] == 'resolved' for r in chans)}/{len(chans)}; "
                                   f"RMSD {min(f(r['rmsd_A']) for r in res):.3f}-{max(f(r['rmsd_A']) for r in res):.3f} A "
                                   f"(median {statistics.median(f(r['rmsd_A']) for r in res):.3f}); max |diff| "
                                   f"{max(f(r['max_abs_diff_A']) for r in res):.3f} A; |min radius diff| "
                                   f"{min(dmin):.3f}-{max(dmin):.3f} A") if res else "no resolved entry"})
    write_csv(R / "validation_channels.csv", rows,
              comment="Paper table: channel profiles vs HOLE and geometric conductance vs measurement (from the "
                      "per-comparison CSVs). Software comparison on curated static structures.")
    for name, path, text in (
            ("tunnels vs CAVER 3.0.3", R / "tunnels" / "summary.csv", None),
            ("pocket volumes vs pyKVFinder/POVME", R / "volumes" / "volumes.csv", None),
            ("GLUT1 trajectory pocket vs MDpocket", R / "trajectory_pocket" / "summary.csv", None),
            ("automatic-probe sensitivity", R / "probe_sensitivity" / "summary.csv", None)):
        if not path.exists():
            summary.append({"benchmark": name, "entries": 0, "result": "not run"})
            continue
        data = read_csv(path)
        if name.startswith("tunnels"):
            auto = [r for r in data if r["crevice_mode"] == "auto"]
            parts = [f"{r['id']}: CAVER top cluster {f(r['caver_top_cluster_bottleneck_A']):.2f} A / CREVICE "
                     f"{f(r['crevice_bottleneck_of_that_route_A']):.2f} A ({f(r['angle_to_top_cluster_deg']):.0f} deg)"
                     for r in auto if r.get("crevice_bottleneck_of_that_route_A")]
            summary.append({"benchmark": name, "entries": len(auto), "result": "; ".join(parts)})
        elif name.startswith("pocket"):
            ids = list(dict.fromkeys(r["id"] for r in data))
            parts = []
            for i in ids:
                v = {r["tool"]: f(r["volume_A3"]) for r in data if r["id"] == i}
                parts.append(f"{i}: cast {v.get('crevice_cast')}, pyKVFinder(0.8) {v.get('pykvfinder_probe_in_0.8')}, "
                             f"pyKVFinder(1.4) {v.get('pykvfinder_probe_in_1.4')}, POVME {v.get('povme')} A3")
            summary.append({"benchmark": name, "entries": len(ids), "result": "; ".join(parts)})
        elif name.startswith("GLUT1"):
            r = data[0]
            summary.append({"benchmark": name, "entries": r.get("frames"),
                            "result": f"mean CREVICE {r.get('crevice_mean_A3')} A3, MDpocket {r.get('mdpocket_mean_A3')} A3, "
                                      f"Pearson r {r.get('pearson_r')}"})
        else:
            res = [r for r in data if r["default_status"] == "resolved"]
            summary.append({"benchmark": name, "entries": len(data),
                            "result": (f"variants per channel {data[0]['variants']}; over the {len(res)} channels resolved by "
                                       f"default: probe changed in {sum(int(r['probe_changed']) for r in res)} of "
                                       f"{sum(int(r['variants']) for r in res)} variant runs, status changed in "
                                       f"{sum(int(r['status_changed']) for r in res)}, min radius changed > 0.01 A in "
                                       f"{sum(int(r['min_radius_changed_gt_0.01A']) for r in res)}; unresolved by default: "
                                       f"{', '.join(r['id'] for r in data if r['default_status'] != 'resolved')}")})
    write_csv(R / "validation_summary.csv", summary, comment="Headline numbers of each benchmark (make_tables.py).")
    for s in summary:
        print(s["benchmark"], "::", s["result"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Geometric conductance estimate vs measured single-channel conductance.

Usage (after ``channels_vs_hole.py``)::

    python benchmarks/conductance.py [--work-dir DIR] [--radius-set simple|bondi]

For every channel in ``inputs/conductance_literature.csv`` this computes, under
the bulk resistivity of the recording solution given there:

* ``crevice_g_macro_pS`` / ``crevice_g_access_pS`` /
  ``crevice_g_hole_corrected_pS``: :func:`crevice.fingerprint.conductance_estimate`
  of CREVICE's resolved profile (the dense pass written by
  ``channels_vs_hole.py``; mouth to mouth, which is the fingerprint's LUMEN
  segment of a profile without casts);
* ``hole_g_macro_same_range_pS``: the same formula on HOLE's sampled profile
  restricted to CREVICE's resolved range (like-for-like with CREVICE);
* ``hole_g_macro_full_pS``: HOLE's own ``Gmacro`` (its ``F`` integrates the
  whole HOLE profile out to ``ENDRAD``), rescaled from HOLE's 1 M KCl
  resistivity (1/12 Ohm m) to the recording resistivity.

Ratios are estimate / measured. A channel CREVICE does not resolve has no
CREVICE estimate (reported as unresolved, never as zero). The per-pore value
is used for OmpF (trimer value / 3) because one monomer pore is profiled.
Output: ``results/conductance/conductance_vs_measured_<set>.csv``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INPUTS, REPO, RESULTS, default_work_dir, read_csv, write_csv  # noqa: E402

sys.path.insert(0, str(REPO / "src"))
from crevice.fingerprint import conductance_estimate  # noqa: E402

HOLE_RHO_1M_KCL = 1.0 / 12.0
PER_PORE_DIVISOR = {"2OMF": 3.0}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--radius-set", choices=["simple", "bondi"], default="simple")
    args = ap.parse_args(argv)
    tag = args.radius_set
    rows = []
    for lit in read_csv(INPUTS / "conductance_literature.csv"):
        acc = lit["id"]
        cur = args.work_dir / "channels" / acc
        comp = json.loads((cur / f"comparison_{tag}.json").read_text())["summary"]
        axis = json.loads((cur / "axis.json").read_text())
        u, p0 = np.array(axis["axis_unit"]), np.array(axis["origin"])
        rho = float(lit["electrolyte_resistivity_ohm_m"])
        measured = float(lit["measured_pS"]) / PER_PORE_DIVISOR.get(acc, 1.0)
        row = {"id": acc, "channel": lit["channel"], "radius_set": f"HOLE {tag}.rad",
               "structure_state_matches": lit["structure_state_matches"], "comparable": lit["comparable"],
               "measured_pS_per_pore": measured, "conditions": lit["conditions"], "resistivity_ohm_m": rho,
               "source": lit["source"], "crevice_status": comp["crevice_status"]}
        hole = json.loads((cur / f"hole_{tag}" / "hole_profile.json").read_text())
        if hole.get("Gmacro"):
            g_full = hole["Gmacro"] * HOLE_RHO_1M_KCL / rho
            row.update(hole_g_macro_full_pS=g_full, hole_full_ratio=g_full / measured)
        used = comp.get("crevice_attempt_used")
        dense = cur / f"crevice_{tag}" / used / "dense.json" if used else None
        if dense is not None and not dense.exists():
            dense = dense.with_name("auto.json")   # dense pass failed: the compared profile is the automatic pass
        if dense and dense.exists():
            prof = json.loads(dense.read_text())
            pts = prof["points"]
            s = [(np.array([p["x"], p["y"], p["z"]]) - p0) @ u for p in pts]
            order = np.argsort(s)
            s = [float(s[i]) for i in order]
            r = [float(pts[i]["radius"]) for i in order]
            g = conductance_estimate(s, r, resistivity_ohm_m=rho, electrolyte=None)
            row.update(crevice_profile_pass=dense.stem, crevice_path_type=comp.get("crevice_path_type"), crevice_lumen_length_A=g["length_A"],
                       crevice_min_radius_A=g["min_radius_A"], crevice_g_macro_pS=g["g_macro_pS"],
                       crevice_g_access_pS=g["g_access_pS"], crevice_g_hole_corrected_pS=g["g_hole_corrected_pS"],
                       crevice_macro_ratio=g["g_macro_pS"] / measured, crevice_access_ratio=g["g_access_pS"] / measured,
                       crevice_hole_corrected_ratio=g["g_hole_corrected_pS"] / measured)
            hs, hr = np.array(hole["s"]), np.array(hole["radius"])
            o = np.argsort(hs)
            hs, hr = hs[o], hr[o]
            inside = (hs >= s[0]) & (hs <= s[-1])
            if inside.sum() >= 2:
                gh = conductance_estimate(hs[inside].tolist(), hr[inside].tolist(), resistivity_ohm_m=rho, electrolyte=None)
                row.update(hole_g_macro_same_range_pS=gh["g_macro_pS"], hole_same_range_ratio=gh["g_macro_pS"] / measured)
        rows.append(row)
        print({k: row.get(k) for k in ("id", "crevice_status", "measured_pS_per_pore", "crevice_g_macro_pS",
                                       "crevice_macro_ratio", "hole_g_macro_full_pS")})
    write_csv(args.results / "conductance" / f"conductance_vs_measured_{tag}.csv", rows,
              comment=("Geometric conductance estimates (CREVICE and HOLE) vs measured single-channel conductance. "
                       "Order-of-magnitude estimates; see inputs/conductance_literature.csv for sources and "
                       "conditions. Ratio = estimate / measured (per pore)."))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

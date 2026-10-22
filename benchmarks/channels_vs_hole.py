"""Matched pore-radius profiles: CREVICE ``profile`` vs HOLE 2.3.1 on the channel set.

Usage (CREVICE environment; HOLE located by ``common.find_tool``)::

    python benchmarks/curate_channels.py
    python benchmarks/channels_vs_hole.py [--work-dir DIR] [--only 1GRM 3UKM]

Matched settings
----------------
Both tools read the same curated PDB (``curate_channels.py``), the same
radius file (HOLE's own ``simple.rad``, passed to CREVICE as ``--radii PATH``),
the same axis (HOLE ``CVECT`` = CREVICE ``--axis``) and the same seed (HOLE
``CPOINT`` = CREVICE ``--origin``). The compared quantity is, for each plane
normal to the axis, the radius of the largest sphere centred in that plane
that touches no atom (HOLE "sampled" rows; CREVICE per-sample ``radius`` with
``--probe-radius 0``).

CREVICE attempts
----------------
CREVICE is run with its automatic enclosure probe and otherwise default
options except a 0.25 Å section grid. If the default does not resolve a
channel, two documented variants are tried in order, and every attempt is
recorded: ``--lateral-exits`` (capped ends may leave sideways) and then also
``--search-radius 12`` (transverse search half-width 12 Å instead of 6 Å, for
pores wider than the default window). The first resolved attempt is used; a
second, dense pass re-runs it with the chosen probe and ~0.2 Å axial steps
(the profile an explicit probe gives is exactly the automatic one). If the
dense pass fails, the automatic pass (81+ samples, steps <= 0.75 Å) is compared
and the failure is recorded (``compared_profile``, ``crevice_dense_error``).

Comparison
----------
Axial coordinate ``s = (centre - seed) . axis`` for the sphere centres of both
tools. HOLE is interpolated linearly at CREVICE's ``s`` values inside the
*interior window*: CREVICE's resolved range trimmed by 2 Å at each end
(by a quarter of the range when it is shorter than 8 Å, e.g. capped channels) (HOLE
keeps growing beyond CREVICE's mouths until ``ENDRAD``). Reported: RMSD, mean
and maximum absolute difference, the minimum radius of each tool (HOLE over
its sampled rows within CREVICE's resolved range; also its global minimum)
and their positions. Outputs: ``results/channels/<id>_profiles.csv`` (both
profiles on CREVICE's samples) and ``results/channels/channels_vs_hole.csv``.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (INPUTS, RESULTS, crevice_env, default_work_dir, find_tool, python_for_crevice,  # noqa: E402
                    read_csv, run_logged, write_csv, write_json)
from curate_channels import find_simple_rad  # noqa: E402

HOLE_SAMPLE_A = 0.2
HOLE_ENDRAD_A = 30.0
SECTION_SPACING_A = 0.25
INTERIOR_TRIM_A = 2.0
ATTEMPTS = (("default", []), ("lateral_exits", ["--lateral-exits"]),
            ("lateral_exits_search12", ["--lateral-exits", "--search-radius", "12"]))


# ----------------------------------------------------------------------------- HOLE
def _float(text: str):
    """HOLE prints ``*****`` for values that overflow its format."""
    try:
        return float(text)
    except ValueError:
        return None


def parse_hole_output(text: str) -> dict:
    rows, on = [], False
    result: dict = {}
    for line in text.splitlines():
        if line.startswith(" cenxyz.cvec"):
            on, rows = True, []
            continue
        if on:
            if not line.strip():
                on = False
                continue
            parts = line.split()
            rows.append({"cvec": float(parts[0]), "radius": float(parts[1]), "cen_line_D": _float(parts[2]),
                         "sum_s_over_area": _float(parts[3]), "tag": "sampled" if "sampled" in line else "mid-point"})
        s = line.strip()
        if s.startswith("(TAG"):
            toks = s.strip("()").replace("=", " ").split()
            for key in ("Rmin", "Gmacro", "Conn_Gmacro"):
                if key in toks:
                    result[key] = _float(toks[toks.index(key) + 1])
        if "F= sum(ds/area) along channel is" in s:
            result["F_per_A"] = _float(s.split("is")[1].split()[0])
        if s.startswith("Minimum radius found:"):
            result["min_radius_line"] = s
    result["rows"] = rows
    return result


def read_sph(path: Path):
    centres, radii = [], []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20] == "SPH" and line[22:26].strip() != "-888":
            centres.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            radii.append(float(line[54:60]))
    return np.array(centres), np.array(radii)


def run_hole(hole: Path, pdb: Path, rad: Path, origin, axis, run: Path, *, sample=HOLE_SAMPLE_A,
             endrad=HOLE_ENDRAD_A, seed=1) -> dict:
    """Run HOLE in ``run`` (short relative names: HOLE reads 80-column lines)."""
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)
    shutil.copy(pdb, run / "in.pdb")
    shutil.copy(rad, run / "m.rad")
    deck = ["COORD in.pdb", "RADIUS m.rad", "SPHPDB hole.sph", f"SAMPLE {sample}", f"ENDRAD {endrad}", "IGNORE",
            "SHORTO 0", f"RASEED {seed}", "CPOINT {:.6f} {:.6f} {:.6f}".format(*origin),
            "CVECT {:.8f} {:.8f} {:.8f}".format(*axis)]
    (run / "hole.inp").write_text("\n".join(deck) + "\n")
    import time
    t0 = time.perf_counter()
    with open(run / "hole.inp") as fin, open(run / "hole.out", "w") as fout:
        import subprocess
        code = subprocess.run([str(hole)], stdin=fin, stdout=fout, stderr=subprocess.STDOUT, cwd=run).returncode
    wall = time.perf_counter() - t0
    parsed = parse_hole_output((run / "hole.out").read_text())
    centres, sph_r = read_sph(run / "hole.sph")
    # Match each "sampled" row to its sphere by the projection on the axis (HOLE's
    # cenxyz.cvec = centre . cvect) and the radius; the .sph order is the run order.
    sampled = [r for r in parsed["rows"] if r["tag"] == "sampled"]
    u, p0 = np.asarray(axis, float), np.asarray(origin, float)
    u = u / np.linalg.norm(u)
    proj = centres @ u if len(centres) else np.array([])
    keep_c, keep_r, errs = [], [], []
    for r in sampled:
        if not len(proj):
            break
        k = int(np.argmin(np.abs(proj - r["cvec"]) + np.abs(sph_r - r["radius"])))
        errs.append(max(abs(proj[k] - r["cvec"]), abs(sph_r[k] - r["radius"])))
        keep_c.append(centres[k])
        keep_r.append(r["radius"])
    centres_m = np.array(keep_c).reshape(-1, 3)
    match_err = float(max(errs)) if errs else None
    s = ((centres_m - p0) @ u).tolist()
    n = len(keep_r)
    parsed.update({"exit_code": code, "wall_s": wall, "s": s, "radius": keep_r,
                   "centres": centres_m.tolist(), "sph_match_max_abs": match_err, "n_sampled": n})
    write_json(run / "hole_profile.json", {k: v for k, v in parsed.items() if k != "rows"})
    return parsed


# ----------------------------------------------------------------------------- CREVICE
def crevice_profile(pdb: Path, rad: Path, origin, axis, run: Path, extra: list[str], samples: int | None = None,
                    probe: float | None = None, tag: str = "pass") -> tuple[int, float, dict | None]:
    run.mkdir(parents=True, exist_ok=True)
    out = run / f"{tag}.json"
    cmd = [python_for_crevice(), "-m", "crevice", "profile", str(pdb), "--offline",
           "--axis=" + ",".join(f"{v:.10f}" for v in axis), "--origin=" + ",".join(f"{v:.6f}" for v in origin),
           "--section-spacing", str(SECTION_SPACING_A), "--radii", str(rad), *extra, "-o", str(out)]
    if samples:
        cmd += ["--samples", str(samples)]
    if probe is not None:
        cmd += ["--enclosure-radius", f"{probe:g}"]
    code, wall = run_logged(cmd, run / "commands.log", env=crevice_env())
    data = json.loads(out.read_text()) if code == 0 and out.exists() else None
    return code, wall, data


def last_error(log: Path) -> str:
    lines = [l for l in log.read_text().splitlines() if l.strip() and not l.startswith(("$", "["))]
    for l in reversed(lines):
        if "Error" in l or "error" in l or "No enclosure" in l or "refus" in l:
            return l.strip()[:400]
    return lines[-1].strip()[:400] if lines else ""


def run_crevice(pdb, rad, origin, axis, base: Path) -> dict:
    record = {"attempts": []}
    for name, extra in ATTEMPTS:
        run = base / name
        if run.exists():
            shutil.rmtree(run)
        code, wall, data = crevice_profile(pdb, rad, origin, axis, run, extra, tag="auto")
        att = {"attempt": name, "options": " ".join(extra), "exit_code": code, "wall_s": wall}
        if data is None:
            att["status"] = "unresolved"
            att["reason"] = last_error(run / "commands.log")
            record["attempts"].append(att)
            continue
        sel = data["metadata"].get("enclosure_probe_selection", {})
        att.update(status="resolved", chosen_probe_A=sel.get("chosen_A"), selection_reason=sel.get("reason"),
                   path_type=data["metadata"].get("path_type", "through_channel"), min_radius_A=data["min_radius"])
        record["attempts"].append(att)
        # Dense pass with the chosen probe, ~0.2 Å steps.
        t = [p["axis_position"] for p in data["points"]]
        samples = max(81, math.ceil((max(t) - min(t)) / 0.2) + 1)
        code2, wall2, dense = crevice_profile(pdb, rad, origin, axis, run, extra, samples=samples,
                                              probe=sel.get("chosen_A"), tag="dense")
        record.update(used=name, auto=data, dense=dense, dense_exit=code2, dense_wall_s=wall2, samples=samples)
        return record
    record["used"] = None
    return record


# ----------------------------------------------------------------------------- comparison
def compare(crev: dict, hole: dict, origin, axis) -> tuple[dict, list[dict]]:
    u = np.asarray(axis, float) / np.linalg.norm(axis)
    p0 = np.asarray(origin, float)
    pts = crev["points"]
    cs = np.array([(np.array([p["x"], p["y"], p["z"]]) - p0) @ u for p in pts])
    cr = np.array([p["radius"] for p in pts])
    hs, hr = np.array(hole["s"]), np.array(hole["radius"])
    order = np.argsort(hs)
    hs, hr = hs[order], hr[order]
    lo, hi = cs.min(), cs.max()
    trim = min(INTERIOR_TRIM_A, 0.25 * (hi - lo))
    window = (cs >= lo + trim) & (cs <= hi - trim) & (cs >= hs.min()) & (cs <= hs.max())
    h_at = np.interp(cs, hs, hr, left=np.nan, right=np.nan)
    diff = cr - h_at
    d = diff[window]
    in_range = (hs >= lo) & (hs <= hi)
    metrics = {
        "crevice_s_range_A": [float(lo), float(hi)],
        "interior_window_A": [float(lo + trim), float(hi - trim)], "trim_A": float(trim),
        "n_compared": int(window.sum()),
        "rmsd_A": float(np.sqrt(np.mean(d ** 2))) if len(d) else None,
        "mean_diff_A": float(np.mean(d)) if len(d) else None,
        "max_abs_diff_A": float(np.max(np.abs(d))) if len(d) else None,
        "s_of_max_abs_diff_A": float(cs[window][np.argmax(np.abs(d))]) if len(d) else None,
        "crevice_min_radius_A": float(cr.min()), "crevice_min_s_A": float(cs[np.argmin(cr)]),
        "crevice_continuous_bottleneck_A": crev["metadata"].get("continuous_bottleneck_radius_A"),
        "hole_min_radius_in_range_A": float(hr[in_range].min()) if in_range.any() else None,
        "hole_min_s_in_range_A": float(hs[in_range][np.argmin(hr[in_range])]) if in_range.any() else None,
        "hole_min_radius_global_A": float(hr.min()), "hole_min_s_global_A": float(hs[np.argmin(hr)]),
    }
    rows = [{"s_A": float(s), "crevice_radius_A": float(r), "hole_radius_interp_A": float(h),
             "diff_A": float(x), "in_interior_window": bool(w)} for s, r, h, x, w in zip(cs, cr, h_at, diff, window)]
    return metrics, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--radius-set", choices=["simple", "bondi"], default="simple",
                    help="HOLE distribution radius file given to both tools (default simple.rad, HOLE's default)")
    args = ap.parse_args(argv)
    hole = find_tool("hole")
    if hole is None:
        raise SystemExit("HOLE not found; set CREVICE_HOLE or run install_comparators.sh")
    rad = find_simple_rad(args.radius_set)
    tag = args.radius_set
    summary = []
    for row in read_csv(INPUTS / "channels.csv"):
        acc = row["id"]
        if args.only and acc not in args.only:
            continue
        cur = args.work_dir / "channels" / acc
        axis_rec = json.loads((cur / "axis.json").read_text())
        u, p0 = axis_rec["axis_unit"], axis_rec["origin"]
        pdb = cur / "curated.pdb"
        h = run_hole(hole, pdb, rad, p0, u, cur / f"hole_{tag}")
        c = run_crevice(pdb, rad, p0, u, cur / f"crevice_{tag}")
        rec = {"id": acc, "pdb_id": row["pdb_id"], "name": row["name"], "state": row["state"],
               "curation": row["curation"], "atoms": axis_rec["atoms_kept"], "radius_set": f"HOLE {tag}.rad (both tools)",
               "hole_exit_code": h["exit_code"], "hole_wall_s": h["wall_s"], "hole_Gmacro_1M_KCl_pS": h.get("Gmacro"),
               "crevice_status_default": c["attempts"][0]["status"],
               "crevice_status": "resolved" if c["used"] else "unresolved", "crevice_attempt_used": c["used"] or "",
               "crevice_attempts": "; ".join(f"{a['attempt']}={a['status']}" for a in c["attempts"]),
               "crevice_unresolved_reason": next((a.get("reason", "") for a in c["attempts"] if a["status"] != "resolved"), ""),
               "crevice_wall_s_auto": sum(a["wall_s"] for a in c["attempts"])}
        hs = np.array(h["s"])
        hr = np.array(h["radius"])
        rec.update(hole_min_radius_global_A=float(hr.min()) if len(hr) else None,
                   hole_min_s_global_A=float(hs[np.argmin(hr)]) if len(hr) else None)
        if c["used"]:
            att = next(a for a in c["attempts"] if a["attempt"] == c["used"])
            rec.update(crevice_chosen_probe_A=att["chosen_probe_A"], crevice_path_type=att["path_type"],
                       crevice_auto_min_radius_A=att["min_radius_A"], crevice_dense_samples=c["samples"],
                       crevice_dense_exit_code=c["dense_exit"])
            compared = c["dense"] if c["dense"] is not None else c["auto"]
            rec["compared_profile"] = "dense" if c["dense"] is not None else "auto (dense pass failed)"
            if c["dense"] is None:
                rec["crevice_dense_error"] = last_error(cur / f"crevice_{tag}" / c["used"] / "commands.log")
            if compared is not None:
                metrics, prof = compare(compared, h, p0, u)
                rec.update({k: v for k, v in metrics.items() if not isinstance(v, list)})
                rec["crevice_s_min_A"], rec["crevice_s_max_A"] = metrics["crevice_s_range_A"]
                write_csv(args.results / "channels" / f"{acc}_profiles_{tag}.csv", prof,
                          comment=f"{acc}: CREVICE samples ({rec['compared_profile']} pass, enclosure probe {att['chosen_probe_A']} A, "
                                  f"attempt {c['used']}) with HOLE 2.3.1 interpolated at the same s. "
                                  f"s = (centre - seed).axis, A; radius set HOLE {tag}.rad in both tools.")
        write_json(cur / f"comparison_{tag}.json", {"summary": rec, "crevice_attempts": c["attempts"]})
        summary.append(rec)
        print(json.dumps({k: rec.get(k) for k in ("id", "crevice_status", "crevice_attempt_used", "crevice_chosen_probe_A",
                                                   "rmsd_A", "max_abs_diff_A", "crevice_min_radius_A",
                                                   "hole_min_radius_in_range_A", "hole_min_radius_global_A")}))
    if not args.only:
        write_csv(args.results / "channels" / f"channels_vs_hole_{tag}.csv", summary,
                  comment=f"CREVICE profile vs HOLE 2.3.1, matched atoms, radii ({tag}.rad), axis and seed. "
                          "Software comparison on static structures; not a biological validation.")
    else:
        write_csv(args.work_dir / "channels" / (f"partial_{tag}_" + "_".join(args.only) + ".csv"), summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

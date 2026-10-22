"""Enzyme tunnels: CREVICE ``paths`` vs CAVER 3.0.3 on the tunnel set.

Usage (CREVICE environment; Java and CAVER located by ``common.find_tool``)::

    python benchmarks/fetch_inputs.py --set tunnels
    python benchmarks/tunnels_vs_caver.py [--work-dir DIR] [--only 1CQW]

Matched settings
----------------
Both tools read the same curated PDB (polymer heavy atoms) and start from the
same point (centroid of the catalytic atoms listed in ``inputs/tunnels.csv``;
CAVER may move it by up to ``max_distance`` 3 Å, CREVICE snaps it to the
widest free node within ``--seed-tolerance`` 3 Å). Atomic radii: CAVER's
default table (``bin/atom_radii.csv``: C 1.70, N 1.55, O 1.52, S 1.80 Å) is
the Bondi set, so CREVICE runs with ``--radii bondi``. CAVER runs with its
default configuration plus ``seed 1``. CREVICE runs twice: with its default
automatic path probe, and with ``--probe-radius 0.9`` (CAVER's default
probe). The tools define a tunnel differently: CAVER's tunnels are Voronoi
paths from the start to the molecular surface (``shell_radius`` 3 Å,
``shell_depth`` 4 Å), clustered, with length measured along the centre line
to the surface; CREVICE's routes are lattice paths to rolling-probe bulk
(6 Å) or to an exposed node, relaxed into a centre line. Lengths therefore
differ by construction; bottleneck radii are comparable (largest sphere that
passes).

Route agreement
---------------
Each tool's route is summarised by its end point. Every CAVER tunnel cluster
is assigned to the CREVICE route with the smallest angle between the
start-to-end directions (at most 45°); several CAVER clusters can fall into one
CREVICE route, because CREVICE merges routes whose exits are closer than 12 Å
while CAVER clusters at 3.5 Å. For reported CREVICE routes (which carry a centre line) the mean
distance from CREVICE's centre-line points to CAVER's centre line is also
given. Routes that CREVICE leaves out under the automatic probe
(``*_path_left_out.csv``) are matched as well and flagged.

Outputs: ``results/tunnels/tunnels_vs_caver.csv`` (one row per route pair or
unmatched route) and ``results/tunnels/summary.csv`` (one row per enzyme).
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (INPUTS, RESULTS, crevice_env, default_work_dir, find_tool, python_for_crevice,  # noqa: E402
                    read_csv, run_logged, write_csv, write_json)
from curate_channels import pattern_atoms, select, write_pdb, xyz  # noqa: E402
from fetch_inputs import structure_file  # noqa: E402

MATCH_ANGLE_DEG = 45.0
CAVER_PROBE_A = 0.9


def curate(row: dict, work: Path) -> tuple[Path, np.ndarray, dict]:
    src = structure_file(work, row["pdb_id"], row["assembly"])
    st = gemmi.read_structure(str(src))
    model, kept, dropped = select(st, row)
    out = work / "tunnels" / row["id"]
    (out / "pdb").mkdir(parents=True, exist_ok=True)
    chains = list(dict.fromkeys(c for c, _, _ in kept))
    pdb = out / "pdb" / f"{row['id']}.pdb"
    write_pdb(kept, {c: c[0] for c in chains}, pdb, f"{row['id']} curated for the CREVICE tunnel benchmark")
    atoms = pattern_atoms(model, row["seed_atoms"])
    seed = np.array([xyz(a) for _, _, a in atoms]).mean(0)
    coords = np.array([xyz(a) for _, _, a in kept])
    rec = {"id": row["id"], "source": src.name, "atoms_kept": len(kept), "dropped": dict(dropped),
           "seed": seed.tolist(), "seed_atoms": [f"{c}:{r.name}{r.seqid.num}:{a.name}" for c, r, a in atoms],
           "seed_nearest_atom_centre_A": float(np.linalg.norm(coords - seed, axis=1).min())}
    write_json(out / "curation.json", rec)
    return pdb, seed, rec


def run_caver(java: Path, caver_jar: Path, pdb: Path, seed, run: Path) -> dict:
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)
    conf = run / "config.txt"
    conf.write_text("starting_point_coordinates {:.3f} {:.3f} {:.3f}\nseed 1\n".format(*seed))
    home = caver_jar.parent
    code, wall = run_logged([str(java), "-Xmx4g", "-jar", str(caver_jar), "-home", str(home), "-pdb", str(pdb.parent),
                             "-conf", str(conf), "-out", str(run / "out")], run / "caver.log")
    tunnels = []
    chars = run / "out" / "analysis" / "tunnel_characteristics.csv"
    profiles = run / "out" / "analysis" / "tunnel_profiles.csv"
    if code == 0 and chars.exists():
        with open(chars) as fh:
            reader = csv.reader(fh, skipinitialspace=True)
            header = next(reader)
            for r in reader:
                d = dict(zip(header, r))
                tunnels.append({"cluster": int(d["Tunnel cluster"]), "tunnel": int(d["Tunnel"]),
                                "bottleneck_radius_A": float(d["Bottleneck radius"]), "length_A": float(d["Length"]),
                                "curvature": float(d["Curvature"]), "throughput": float(d["Throughput"])})
        coords: dict = {}
        with open(profiles) as fh:
            reader = csv.reader(fh, skipinitialspace=True)
            next(reader)
            for r in reader:
                key = (int(r[1]), int(r[2]))
                axis = r[12].strip()
                if axis in ("X", "Y", "Z", "R", "length"):
                    coords.setdefault(key, {})[axis] = [float(v) for v in r[13:] if v.strip()]
        for t in tunnels:
            c = coords.get((t["cluster"], t["tunnel"]), {})
            if all(k in c for k in "XYZ"):
                t["points"] = np.column_stack([c["X"], c["Y"], c["Z"]]).tolist()
                t["radii"] = c.get("R")
                t["start"], t["end"] = t["points"][0], t["points"][-1]
    return {"exit_code": code, "wall_s": wall, "tunnels": tunnels}


def run_paths(pdb: Path, seed, run: Path, probe: str) -> dict:
    if run.exists():
        shutil.rmtree(run)
    cmd = [python_for_crevice(), "-m", "crevice", "paths", str(pdb), "--seed=" + ",".join(f"{v:.3f}" for v in seed),
           "--radii", "bondi", "--no-scenes", "--probe-radius", probe, "--workers", "2", "-o", str(run),
           "--prefix", "paths"]
    code, wall = run_logged(cmd, run.parent / f"{run.name}.log", env=crevice_env())
    routes, left = [], []
    if code == 0 and (run / "paths_paths.json").exists():
        data = json.loads((run / "paths_paths.json").read_text())
        for b in data["branches"]:
            routes.append({"route": f"B{b['branch']}", "reported": True, "bottleneck_radius_A": b["bottleneck_radius_A"],
                           "length_A": b["length_A"], "start": b["start_point_A"], "end": b["exit_point_A"],
                           "points": b["points_A"], "end_type": b["end_type"], "probe_A": b.get("probe_radius_A")})
        lo = run / "paths_path_left_out.csv"
        if lo.exists():
            for r in read_csv(lo):
                left.append({"route": r["route"], "reported": False, "bottleneck_radius_A": float(r["width_A"]),
                             "length_A": float(r["lattice_length_A"]), "start": None,
                             "end": [float(r["mouth_x_A"]), float(r["mouth_y_A"]), float(r["mouth_z_A"])],
                             "points": None, "end_type": r["end_type"], "found_at_probes_A": r["found_at_probes_A"]})
        settings = data.get("settings", {})
        start = data.get("start", {})
    else:
        settings, start = {}, {}
    return {"exit_code": code, "wall_s": wall, "routes": routes, "left_out": left, "settings": settings,
            "start": start}


def angle(seed, a, b) -> float:
    u, v = np.asarray(a) - seed, np.asarray(b) - seed
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))


def mean_line_distance(points, other) -> float | None:
    if not points or not other:
        return None
    p, q = np.asarray(points), np.asarray(other)
    return float(np.mean(np.min(np.linalg.norm(p[:, None, :] - q[None, :, :], axis=2), axis=1)))


def match(seed, crevice_routes, caver_tunnels):
    """Assign every CAVER tunnel to the CREVICE route with the closest end direction (<= 45 deg).

    Several CAVER clusters may fall into one CREVICE route: CREVICE keeps routes
    whose exits are closer than ``--separation`` (12 Å by default) as one route,
    while CAVER clusters tunnels at 3.5 Å.
    """
    out, used_r, used_t = [], set(), set()
    for j, t in enumerate(caver_tunnels):
        if "end" not in t:
            continue
        best = None
        for i, r in enumerate(crevice_routes):
            ang = angle(seed, r["end"], t["end"])
            # Prefer reported routes over left-out ones at equal angle.
            key = (ang, not r["reported"])
            if ang <= MATCH_ANGLE_DEG and (best is None or key < best[0]):
                best = (key, i)
        if best is not None:
            out.append((best[1], j, best[0][0]))
            used_r.add(best[1])
            used_t.add(j)
    return out, used_r, used_t


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args(argv)
    java, jar = find_tool("java"), find_tool("caver")
    if java is None or jar is None:
        raise SystemExit("Java or CAVER not found; set CREVICE_JAVA / CREVICE_CAVER_JAR or run install_comparators.sh")
    pair_rows, summary = [], []
    for row in read_csv(INPUTS / "tunnels.csv"):
        if args.only and row["id"] not in args.only:
            continue
        pdb, seed, cur = curate(row, args.work_dir)
        base = args.work_dir / "tunnels" / row["id"]
        cav = run_caver(java, jar, pdb, seed, base / "caver")
        for mode, probe in (("auto", "auto"), ("probe0.9", f"{CAVER_PROBE_A}")):
            cre = run_paths(pdb, seed, base / f"crevice_{mode}", probe)
            routes = cre["routes"] + cre["left_out"]
            pairs, used_r, used_t = match(seed, routes, cav["tunnels"])
            for i, j, ang in pairs:
                r, t = routes[i], cav["tunnels"][j]
                pair_rows.append({"id": row["id"], "crevice_mode": mode, "crevice_route": r["route"],
                                  "crevice_reported": r["reported"], "caver_cluster": t["cluster"],
                                  "end_direction_angle_deg": ang,
                                  "crevice_bottleneck_A": r["bottleneck_radius_A"], "caver_bottleneck_A": t["bottleneck_radius_A"],
                                  "bottleneck_diff_A": r["bottleneck_radius_A"] - t["bottleneck_radius_A"],
                                  "crevice_length_A": r["length_A"], "caver_length_A": t["length_A"],
                                  "mean_centre_line_distance_A": mean_line_distance(r["points"], t.get("points")),
                                  "crevice_end_type": r["end_type"], "caver_throughput": t["throughput"]})
            for i, r in enumerate(routes):
                if i not in used_r:
                    pair_rows.append({"id": row["id"], "crevice_mode": mode, "crevice_route": r["route"],
                                      "crevice_reported": r["reported"], "caver_cluster": "unmatched",
                                      "crevice_bottleneck_A": r["bottleneck_radius_A"], "crevice_length_A": r["length_A"],
                                      "crevice_end_type": r["end_type"]})
            for j, t in enumerate(cav["tunnels"]):
                if j not in used_t:
                    pair_rows.append({"id": row["id"], "crevice_mode": mode, "crevice_route": "unmatched",
                                      "caver_cluster": t["cluster"], "caver_bottleneck_A": t["bottleneck_radius_A"],
                                      "caver_length_A": t["length_A"], "caver_throughput": t["throughput"]})
            widest_c = max((r["bottleneck_radius_A"] for r in cre["routes"]), default=None)
            widest_v = max((t["bottleneck_radius_A"] for t in cav["tunnels"]), default=None)
            top = cav["tunnels"][0] if cav["tunnels"] else None
            top_pair = next(((routes[i], a) for i, j, a in pairs if j == 0), (None, None))
            summary.append({"id": row["id"], "name": row["name"], "curation": row["curation"], "crevice_mode": mode,
                            "crevice_probe_A": cre["settings"].get("probe_selection", {}).get("chosen_A",
                                                                                               cre["settings"].get("probe_radius_A")),
                            "crevice_routes_reported": len(cre["routes"]), "crevice_routes_left_out": len(cre["left_out"]),
                            "caver_tunnel_clusters": len(cav["tunnels"]),
                            "caver_clusters_matched": len(pairs), "crevice_reported_routes_matched": len({i for i, _, _ in pairs if routes[i]["reported"]}),
                            "caver_clusters_matched_to_reported": sum(routes[i]["reported"] for i, _, _ in pairs),
                            "crevice_widest_bottleneck_A": widest_c, "caver_widest_bottleneck_A": widest_v,
                            "caver_top_cluster_bottleneck_A": top["bottleneck_radius_A"] if top else None,
                            "caver_top_cluster_length_A": top["length_A"] if top else None,
                            "crevice_route_matching_top_cluster": top_pair[0]["route"] if top_pair[0] else "",
                            "crevice_bottleneck_of_that_route_A": top_pair[0]["bottleneck_radius_A"] if top_pair[0] else None,
                            "crevice_length_of_that_route_A": top_pair[0]["length_A"] if top_pair[0] else None,
                            "angle_to_top_cluster_deg": top_pair[1],
                            "crevice_wall_s": cre["wall_s"], "caver_wall_s": cav["wall_s"],
                            "crevice_exit_code": cre["exit_code"], "caver_exit_code": cav["exit_code"]})
            write_json(base / f"comparison_{mode}.json", {"curation": cur, "caver": {k: v for k, v in cav.items()},
                                                          "crevice": cre})
        print(json.dumps(summary[-1]))
    if not args.only:
        write_csv(args.results / "tunnels" / "tunnels_vs_caver.csv", pair_rows,
                  comment="CREVICE paths vs CAVER 3.0.3: each CAVER cluster assigned to the CREVICE route with the closest end direction (<= 45 deg); "
                          "Bondi radii in both; same start point. Software comparison, not a validated result.")
        write_csv(args.results / "tunnels" / "summary.csv", summary,
                  comment="Per-enzyme summary of the CREVICE/CAVER tunnel comparison.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

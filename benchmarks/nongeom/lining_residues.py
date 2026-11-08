"""Boundary (lining) residue attribution, partners and residue contact dynamics against known answers.

A. ``lining_synthetic.csv`` (exact ground truth): a spherical cavity of
   radius R whose wall is a dense shell of atoms split into "residues" by
   polar-angle bands (cos theta intervals). On a sphere, a band between
   cos theta = a and b has area fraction (b - a) / 2 exactly, so the share of
   the boundary each lining residue should receive is known analytically. A
   second shell of residues (far in sequence) sits behind the wall within
   4.5 A of it (expected role ``nonlocal_lining_partner``), and one residue
   lies far away (``other``). :func:`crevice.residue_evidence.boundary_residue_evidence`
   is run on the voxelised cavity at several grid spacings; an independent
   reference assigns uniformly spread points of the exact sphere to the
   nearest atom surface (no marching cubes, no CREVICE code).
B. ``lining_4pyp.csv`` (real structure, cross-check of definitions): the
   documented 4PYP cast (``crevice cavities --cast``) and its boundary residues
   (``crevice residues --evidence``), compared with (i) a distance-based
   lining definition (any heavy-atom surface within lining_distance + h/2 of
   an occupied voxel of the cast) and (ii) the residues within 4 A of the
   bound detergent BNG, which sits in the cavity in the crystal.
C. ``contacts_glut1.csv`` (real trajectory, implementation cross-check):
   residue-pair contact occupancy (minimum heavy-atom centre distance <= 4.5 A)
   and centroid DCCM from :class:`crevice.residue_dynamics.ResidueDynamics`
   with the networks CREVICE builds per aligned frame, vs MDAnalysis
   (``capped_distance``; Kabsch fit on CA; centroid cross-correlation).

Run::

    PYTHONPATH=src python benchmarks/nongeom/lining_residues.py --parts ab
    PYTHONPATH=src python benchmarks/nongeom/lining_residues.py --parts c --topology system.gro --trajectory production.xtc
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from nongeom_common import RESULTS, dump, ensure_src_on_path, machine_record, write_csv

ensure_src_on_path()

#: Bondi radii (the CREVICE default set agrees for these elements).
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def pdb_line(serial, name, resname, chain, resid, xyz, element):
    return (f"ATOM  {serial:5d} {name:<4s} {resname:3s} {chain}{resid:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2s}\n")


def synthetic_cavity(path, *, radius, edges, wall_atoms=2400, gap=0.3):
    """Write a PDB with a band-partitioned wall, a partner shell and a distant residue; return the truth."""
    wall_r = radius + RADII["C"] + gap
    dirs = fibonacci_sphere(wall_atoms)
    band = np.clip(np.searchsorted(edges, dirs[:, 2], side="right") - 1, 0, len(edges) - 2)
    lines, serial = [], 1
    truth = {}
    for k in range(len(edges) - 1):
        resid = 1 + k
        for xyz in dirs[band == k] * wall_r:
            lines.append(pdb_line(serial, "C", "ALA", "A", resid, xyz, "C"))
            serial += 1
        truth[f"A:ALA{resid}"] = {"role": "boundary_lining", "fraction": (edges[k + 1] - edges[k]) / 2}
    # Partner shell: 6 residues of another chain, 3.8 A behind the wall (within the 4.5 A contact cutoff).
    # (Partners are "nonlocal": another chain or more than two positions apart in the *observed* residue
    # order, so a residue numbered 100 written right after residue 5 of the same chain would count as local.)
    pdirs = fibonacci_sphere(120)
    for j, xyz in enumerate(pdirs * (wall_r + 3.8)):
        resid = 100 + j // 20
        lines.append(pdb_line(serial, "C", "GLY", "B", resid, xyz, "C"))
        serial += 1
        truth[f"B:GLY{resid}"] = {"role": "nonlocal_lining_partner", "fraction": 0.0}
    lines.append(pdb_line(serial, "C", "LEU", "C", 500, (0.0, 0.0, wall_r + 20.0), "C"))
    truth["C:LEU500"] = {"role": "other", "fraction": 0.0}
    Path(path).write_text("".join(lines) + "END\n")
    return truth, wall_r


def reference_fractions(atom_xyz, atom_res, radius_surface, lining_distance, n_points=200000):
    """Nearest-atom-surface attribution of uniformly spread points on the exact sphere (area-weighted)."""
    from scipy.spatial import cKDTree
    pts = fibonacci_sphere(n_points) * radius_surface
    tree = cKDTree(atom_xyz)
    d, idx = tree.query(pts)
    gap = d - RADII["C"]    # every wall atom is carbon
    owner = np.asarray(atom_res)[idx]
    keep = gap <= lining_distance
    out = {}
    for r in np.unique(owner[keep]):
        out[str(r)] = float(np.sum(owner[keep] == r) / n_points)
    return out, float(keep.mean())


def part_a(args):
    from crevice import load_structure
    from crevice.residue_evidence import boundary_residue_evidence
    rows = []
    work = Path(args.work_dir) / "lining_synthetic"
    work.mkdir(parents=True, exist_ok=True)
    configs = [("equal_bands", [-1.0, -0.6, -0.2, 0.2, 0.6, 1.0]), ("unequal_bands", [-1.0, -0.6, 0.0, 0.3, 0.9, 1.0])]
    radius = 6.0
    for label, edges in configs:
        pdb = work / f"{label}.pdb"
        truth, wall_r = synthetic_cavity(pdb, radius=radius, edges=np.asarray(edges))
        frame = load_structure(str(pdb))
        atoms = frame.selected_atoms(include_hydrogen=False, include_hetero=False)
        xyz = np.asarray([a.coord for a in atoms])
        res = [a.residue_key.label for a in atoms]
        wall = np.asarray([r.split(":")[1].startswith("ALA") for r in res])
        for h in (0.25, 0.5, 0.7):
            n = int(math.ceil((radius + 2 * h) / h))
            ax = np.arange(-n, n + 1) * h
            gx, gy, gz = np.meshgrid(ax, ax, ax, indexing="ij")
            grid = (np.sqrt(gx ** 2 + gy ** 2 + gz ** 2) <= radius).astype(float)
            origin = np.array([ax[0]] * 3)
            t0 = time.perf_counter()
            report = boundary_residue_evidence(frame, grid, origin, np.eye(3) * h, lining_distance=1.5, contact_cutoff=4.5)
            seconds = time.perf_counter() - t0
            # The marching-cubes surface of a binary sphere lies about h/2 outside the last occupied voxel centre.
            ref, ref_assigned = reference_fractions(xyz[wall], np.asarray(res)[wall], radius + h / 2, 1.5)
            by = {r["residue"]: r for r in report["residues"]}
            for resid, t in truth.items():
                got = by.get(resid, {})
                rows.append({"configuration": label, "spacing_A": h, "residue": resid, "expected_role": t["role"],
                             "crevice_role": got.get("role"), "role_matches": got.get("role") == t["role"],
                             "analytic_fraction": t["fraction"], "reference_fraction": ref.get(resid, 0.0),
                             "crevice_fraction": got.get("boundary_fraction", 0.0),
                             "crevice_minus_analytic": got.get("boundary_fraction", 0.0) - t["fraction"],
                             "crevice_minus_reference": got.get("boundary_fraction", 0.0) - ref.get(resid, 0.0),
                             "crevice_assigned_fraction": report["assigned_boundary_fraction"],
                             "reference_assigned_fraction": ref_assigned,
                             "crevice_boundary_area_A2": report["boundary_area_A2"],
                             "sphere_area_A2": 4 * math.pi * (radius + h / 2) ** 2, "seconds": seconds})
            worst = max(abs(r["crevice_minus_analytic"]) for r in rows[-len(truth):])
            print(f"lining {label} h={h}: roles ok={all(r['role_matches'] for r in rows[-len(truth):])} "
                  f"max |crevice-analytic|={worst:.4f} area={report['boundary_area_A2']:.1f} vs "
                  f"{4 * math.pi * (radius + h / 2) ** 2:.1f}", flush=True)
    return rows


def run(cmd, log):
    with open(log, "a") as handle:
        handle.write("$ " + " ".join(map(str, cmd)) + "\n")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        handle.write(proc.stdout + proc.stderr + f"[exit {proc.returncode}]\n")
    if proc.returncode:
        raise RuntimeError(f"{cmd} failed; see {log}")
    return proc.stdout


def part_b(args):
    from scipy.spatial import cKDTree
    from crevice.residue_evidence import read_binary_dx
    work = Path(args.work_dir) / "lining_4pyp"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    log = work / "commands.log"
    cif = Path(args.cif_4pyp)
    py = [sys.executable, "-m", "crevice"]
    out_cast = run(py + ["cavities", "--cast", str(cif), "--out-dir", str(work / "cast")], log)
    out_ev = run(py + ["residues", "--evidence", str(work / "cast" / "4PYP_viewer.pdb"), "--volume-dx",
                       str(work / "cast" / "4PYP_volume.dx"), "--out-dir", str(work / "evidence"), "--prefix", "4PYP"], log)
    report = json.loads((work / "evidence" / "4PYP_residue_evidence.json").read_text())
    lining = {r["residue"]: r for r in report["residues"] if r["role"] == "boundary_lining"}
    grid, origin, deltas = read_binary_dx(work / "cast" / "4PYP_volume.dx")
    h = float(abs(np.linalg.det(deltas)) ** (1 / 3))
    occupied = np.argwhere(grid > 0.5) @ deltas + origin
    tree = cKDTree(occupied)
    # Independent read of the deposited structure (gemmi-free mmCIF parse via MDAnalysis is not available
    # for CIF; use CREVICE's loader only to read coordinates, then apply the reference definitions here).
    from crevice import load_structure
    frame = load_structure(str(cif))
    atoms = [a for a in frame.atoms if (a.element or "").upper() != "H"]
    protein = [a for a in atoms if not a.hetero]
    ligand = [a for a in atoms if a.hetero and a.resname.upper() == "BNG"]
    lab = lambda a: a.residue_key.label  # noqa: E731
    xyz = np.asarray([a.coord for a in protein])
    rad = np.asarray([RADII.get((a.element or "C").upper(), 1.7) for a in protein])
    d, _ = tree.query(xyz)
    gap = d - rad
    distance_based = {}
    for a, g in zip(protein, gap):
        if g <= 1.5 + h / 2:
            distance_based[lab(a)] = min(g, distance_based.get(lab(a), 9e9))
    lig_xyz = np.asarray([a.coord for a in ligand]) if ligand else np.empty((0, 3))
    lig_res = set()
    if len(lig_xyz):
        lt = cKDTree(lig_xyz)
        for a in protein:
            if lt.query(a.coord)[0] <= 4.0:
                lig_res.add(lab(a))
    names_lining = set(lining)
    # Viewer identities are restored to the original chain/numbering in the evidence report; align label styles.
    def norm(label):
        chain, rest = label.split(":")
        return f"{chain}:{rest}"
    lining_n = {norm(x) for x in names_lining}
    union = lining_n | set(distance_based)
    rows = []
    for r in sorted(union | lig_res):
        rows.append({"residue": r, "crevice_lining": r in lining_n,
                     "crevice_boundary_fraction": lining.get(r, {}).get("boundary_fraction"),
                     "distance_based_lining": r in distance_based,
                     "distance_based_min_gap_A": distance_based.get(r),
                     "within_4A_of_BNG": r in lig_res})
    inter = lining_n & set(distance_based)
    strip = lambda line: re.sub(r" out_dir=\S+", "", line)  # noqa: E731
    summary = {"cast_stdout": strip(out_cast.strip().splitlines()[-1]),
               "evidence_stdout": strip(out_ev.strip().splitlines()[-1]),
               "spacing_A": h, "crevice_lining": len(lining_n), "distance_based": len(distance_based),
               "both": len(inter), "jaccard": len(inter) / len(union) if union else None,
               "crevice_only": sorted(lining_n - set(distance_based)),
               "distance_only_count": len(set(distance_based) - lining_n),
               "bng_contact_residues": len(lig_res), "bng_contacts_in_crevice_lining": len(lig_res & lining_n),
               "bng_contacts_not_lining": sorted(lig_res - lining_n), "bng_heavy_atoms": len(ligand)}
    print(f"4PYP: {summary}", flush=True)
    return rows, summary


def part_c(args):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
    from crevice.networks import build_residue_network
    from crevice.residue_dynamics import ResidueDynamics
    from crevice.trajectory import align_frame_report, load_trajectory_report
    t0 = time.perf_counter()
    traj, _ = load_trajectory_report(args.topology, args.trajectory, selection="protein and not name H*",
                                     start=0, stop=args.stop, stride=args.stride, max_frames=None)
    frames = traj.frames
    first = frames[0]
    dyn = ResidueDynamics(min_occupancy=0.0)
    times = traj.times()
    for f, tt in zip(frames, times):
        aligned, _fit = align_frame_report(f, first)
        dyn.add(None, build_residue_network(aligned, cutoff=args.cutoff, distance_metric="center"),
                frame_index=f.frame_index, time=tt)
    report = dyn.finish(aligned=True)
    crevice_seconds = time.perf_counter() - t0
    # Reference with MDAnalysis.
    u = mda.Universe(args.topology, args.trajectory)
    heavy = u.select_atoms("protein and not name H*")
    ca = u.select_atoms("protein and name CA")
    resid_of = heavy.resids
    labels = {r.resid: f"SYSTEM:{r.resname}{r.resid}" for r in heavy.residues}
    idx = list(range(0, args.stop, args.stride))
    counts = {}
    centroids = []
    u.trajectory[idx[0]]
    ref_ca = ca.positions.astype(float).copy()
    for k in idx:
        u.trajectory[k]
        pairs = capped_distance(heavy.positions, heavy.positions, args.cutoff, return_distances=False)
        a, b = resid_of[pairs[:, 0]], resid_of[pairs[:, 1]]
        keep = a < b
        for pair in set(zip(a[keep].tolist(), b[keep].tolist())):
            counts[pair] = counts.get(pair, 0) + 1
        # Kabsch on CA onto the first frame, then residue heavy-atom centroids.
        x = ca.positions.astype(float)
        xc, yc = x.mean(0), ref_ca.mean(0)
        hmat = (x - xc).T @ (ref_ca - yc)
        uu, _, vt = np.linalg.svd(hmat)
        dsign = np.sign(np.linalg.det(uu @ vt))
        rot = uu @ np.diag([1, 1, dsign]) @ vt
        pos = (heavy.positions.astype(float) - xc) @ rot + yc
        centroids.append(np.asarray([pos[resid_of == r].mean(0) for r in heavy.residues.resids]))
    centroids = np.asarray(centroids)
    fl = centroids - centroids.mean(0)
    order = {r: i for i, r in enumerate(heavy.residues.resids)}
    n = len(idx)
    mine = {(row["source"], row["target"]): row for row in report["contacts"]}
    rows, mismatched, dccm_diff = [], 0, []
    keys = set()
    for (ra, rb), c in counts.items():
        la, lb = sorted((labels[ra], labels[rb]))
        keys.add((la, lb))
        got = mine.get((la, lb))
        occ_c = got["occupancy"] if got else 0.0
        i, j = order[ra], order[rb]
        num = np.mean(np.sum(fl[:, i] * fl[:, j], axis=1))
        den = math.sqrt(np.mean(np.sum(fl[:, i] ** 2, axis=1)) * np.mean(np.sum(fl[:, j] ** 2, axis=1)))
        dccm_ref = float(num / den) if den > 0 else None
        if got and got.get("dccm") is not None and dccm_ref is not None:
            dccm_diff.append(abs(got["dccm"] - dccm_ref))
        mismatched += abs(occ_c - c / n) > 1e-12
        rows.append({"source": la, "target": lb, "reference_occupancy": c / n, "crevice_occupancy": occ_c,
                     "crevice_dccm": got.get("dccm") if got else None, "reference_dccm": dccm_ref})
    extra = [k for k in mine if k not in keys]
    summary = {"frames": n, "stride": args.stride, "pairs_reference": len(counts), "pairs_crevice": len(mine),
               "pairs_only_in_crevice": len(extra), "occupancy_mismatches": int(mismatched),
               "dccm_max_abs_difference": float(max(dccm_diff)) if dccm_diff else None,
               "crevice_seconds": crevice_seconds}
    print(f"GLUT1 contacts: {summary}", flush=True)
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--parts", default="abc")
    parser.add_argument("--work-dir", default=None, help="scratch outputs (default <workspace>/benchmarks/nongeom)")
    parser.add_argument("--cif-4pyp", default=None, help="4PYP mmCIF (default: <workspace>/pdb/4PYP.cif)")
    parser.add_argument("--topology", help="MD topology (GRO, TPR, PSF, ...); needed for the trajectory parts")
    parser.add_argument("--trajectory", help="MD trajectory; needed for the trajectory parts")
    parser.add_argument("--stride", type=int, default=50)
    parser.add_argument("--stop", type=int, default=1001)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--results", default=str(RESULTS))
    args = parser.parse_args()
    from common import default_work_dir, workspace_dir
    if args.work_dir is None:
        args.work_dir = str(default_work_dir() / "nongeom")
    if args.cif_4pyp is None:
        args.cif_4pyp = str(workspace_dir() / "pdb" / "4PYP.cif")
    if "c" in args.parts and not (args.topology and args.trajectory):
        parser.error("part c needs --topology and --trajectory")
    out = Path(args.results)
    summary = {"script": "benchmarks/nongeom/lining_residues.py", "machine": machine_record()}
    if "a" in args.parts:
        rows = part_a(args)
        write_csv(out / "lining_synthetic.csv", rows, comment="Spherical cavity with polar-band residues: boundary "
                  "fraction per residue vs the exact band area fraction and vs an independent nearest-surface "
                  "attribution on the exact sphere. Synthetic.")
        summary["synthetic"] = {"rows": len(rows), "roles_all_match": all(r["role_matches"] for r in rows),
                                "max_abs_crevice_minus_analytic_by_spacing": {
                                    str(h): max(abs(r["crevice_minus_analytic"]) for r in rows if r["spacing_A"] == h)
                                    for h in sorted({r["spacing_A"] for r in rows})},
                                "max_abs_crevice_minus_reference_by_spacing": {
                                    str(h): max(abs(r["crevice_minus_reference"]) for r in rows if r["spacing_A"] == h)
                                    for h in sorted({r["spacing_A"] for r in rows})}}
    if "b" in args.parts:
        rows, s = part_b(args)
        write_csv(out / "lining_4pyp.csv", rows, comment="4PYP cast boundary residues (CREVICE nearest-surface "
                  "attribution) vs a distance-based lining definition and vs residues within 4 A of the crystal "
                  "BNG. Definition cross-check on a static structure; not a functional validation.")
        summary["4pyp"] = s
    if "c" in args.parts:
        rows, s = part_c(args)
        write_csv(out / "contacts_glut1.csv", rows, comment="Residue contact occupancy (min heavy-atom centre "
                  "distance <= 4.5 A) and centroid DCCM: CREVICE ResidueDynamics vs MDAnalysis. Implementation "
                  "cross-check on an MD trajectory that is not a validated reference system.")
        summary["contacts_glut1"] = s
    dump(out / f"lining_residues_summary{'' if args.parts == 'abc' else '_' + args.parts}.json", summary)


if __name__ == "__main__":
    main()

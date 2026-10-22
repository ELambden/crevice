"""Trajectory pocket volume: CREVICE ``trajectory --region`` vs MDpocket on the GLUT1 excerpt.

Usage (CREVICE environment with MDAnalysis; MDpocket located by
``common.find_tool``)::

    crevice fetch --example glut1-excerpt --cache-dir <work-dir>/examples
    python benchmarks/mdpocket_vs_region.py [--work-dir DIR]

Input: the public GLUT1 excerpt (21 frames, every 5 ns of a 100 ns membrane
simulation; CC-BY-4.0) and its prepared region definition for the candidate
sugar site. Both tools analyse the same region:

* CREVICE: ``crevice trajectory --region`` with the prepared definition
  (default settings): per-frame measured volume of the region (probe 1.4 Å,
  0.2 Å grid, frames aligned on protein CA).
* MDpocket (fpocket 4.1.5): frames aligned on the same protein CA atoms with
  MDAnalysis and written as one PDB per frame (all atoms except TIP3 water, as
  CREVICE's obstacle selection), then ``mdpocket -L frames.txt
  --selected_pocket region.pdb`` where ``region.pdb`` holds the reference
  region's grid points thinned to a 1 Å lattice (the selected-pocket file is
  a set of site points). Volume = MDpocket's ``pock_volume`` per snapshot
  (alpha-sphere based, its own definition).

Volumes are not expected to agree in magnitude (different definitions); the
comparison reports both series, the ratio of means and the Pearson
correlation over frames. Output: ``results/trajectory_pocket/glut1_volumes.csv``
(per frame) and ``results/trajectory_pocket/summary.csv``.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO, RESULTS, crevice_env, default_work_dir, find_tool, python_for_crevice, run_logged, write_csv  # noqa: E402

sys.path.insert(0, str(REPO / "src"))
from crevice.residue_evidence import read_binary_dx  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args(argv)
    mdpocket = find_tool("mdpocket")
    if mdpocket is None:
        raise SystemExit("mdpocket not found; set CREVICE_MDPOCKET or run install_comparators.sh")
    ex = args.work_dir / "examples" / "glut1-excerpt"
    region = ex / "glut1_site_candidate_region.json"
    if not region.exists():
        raise SystemExit(f"{region} missing; run crevice fetch --example glut1-excerpt --cache-dir {ex.parent}")
    base = args.work_dir / "mdpocket"
    base.mkdir(parents=True, exist_ok=True)

    # CREVICE
    out = base / "crevice_region"
    if out.exists():
        shutil.rmtree(out)
    code, wall_c = run_logged([python_for_crevice(), "-m", "crevice", "trajectory", "--region", str(region), "--out-dir", str(out),
                               "--workers", str(args.workers), "--skip-analysis-views", "--skip-water-density",
                               "--skip-interactions", "-q"], base / "crevice_region.log", env=crevice_env())
    if code:
        raise SystemExit(f"trajectory --region failed; see {base / 'crevice_region.log'}")
    frames_csv = next(out.glob("*_cavity_frames.csv"))
    with open(frames_csv) as fh:
        crows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    vol_key = "volume_A3"

    # MDpocket input: aligned per-frame PDBs (non-water atoms).
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    defn = json.loads(region.read_text())
    u = mda.Universe(str(ex / defn["system"]["topology"]["path"]), str(ex / defn["system"]["trajectory"]["path"]))
    ref = mda.Universe(str(ex / defn["system"]["topology"]["path"]), str(ex / defn["system"]["trajectory"]["path"]))
    ref.trajectory[defn["system"].get("reference_frame", 0)]
    sel = defn["system"].get("obstacle_selection", "not resname TIP3")
    atoms = u.select_atoms(sel)
    pdb_dir = base / "frames"
    if pdb_dir.exists():
        shutil.rmtree(pdb_dir)
    pdb_dir.mkdir()
    listing = []
    for ts in u.trajectory:
        align.alignto(u, ref, select="protein and name CA", weights=None)
        path = pdb_dir / f"f{ts.frame:04d}.pdb"
        atoms.write(str(path))
        listing.append(f"frames/{path.name}")   # relative: mdpocket aborts on long paths
    (base / "frames.txt").write_text("\n".join(listing) + "\n")
    # Selected pocket: reference-region grid points on a 1 Å lattice.
    grid, origin, deltas = read_binary_dx(ex / defn["reference"]["asset"]["path"])
    idx = np.argwhere(grid > 0.5)
    pts = origin + idx @ deltas
    lattice = np.unique(np.round(pts), axis=0)
    inside = []
    for p in lattice:
        k = np.round((np.linalg.solve(deltas.T, p - origin))).astype(int)
        if np.all(k >= 0) and np.all(k < grid.shape) and grid[tuple(k)] > 0.5:
            inside.append(p)
    lines = [f"ATOM  {i % 100000:5d}  C   STP     1    {x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00" for i, (x, y, z) in enumerate(inside, 1)]
    (base / "region.pdb").write_text("\n".join(lines) + "\nEND\n")
    # Short arguments, including argv[0]: mdpocket (4.1.5 and 4.2.3) aborts with "malloc(): corrupted top
    # size" when invoked through a long path, so it is run through a symlink in the work directory.
    link = base / "mdpocket"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(mdpocket)
    for old in base.glob("mdpout*"):
        old.unlink()
    code, wall_m = run_logged(["./mdpocket", "-L", "frames.txt", "--selected_pocket", "region.pdb", "-o", "mdpout"],
                              base / "mdpocket.log", cwd=base)
    desc = base / "mdpout_descriptors.txt"
    mvals = []
    if desc.exists():
        rows = [l.split() for l in desc.read_text().splitlines() if l.strip()]
        head = rows[0]
        mvals = [float(r[head.index("pock_volume")]) for r in rows[1:]]
    per_frame = []
    for i, r in enumerate(crows):
        per_frame.append({"frame": i, "crevice_volume_A3": float(r[vol_key]),
                          "mdpocket_pock_volume_A3": mvals[i] if i < len(mvals) else None})
    c = np.array([r["crevice_volume_A3"] for r in per_frame], float)
    m = np.array([r["mdpocket_pock_volume_A3"] for r in per_frame if r["mdpocket_pock_volume_A3"] is not None], float)
    summary = {"frames": len(c), "crevice_volume_column": vol_key, "crevice_mean_A3": float(c.mean()),
               "crevice_sd_A3": float(c.std(ddof=1)), "mdpocket_frames": len(m),
               "mdpocket_mean_A3": float(m.mean()) if len(m) else None,
               "mdpocket_sd_A3": float(m.std(ddof=1)) if len(m) > 1 else None,
               "ratio_of_means_mdpocket_over_crevice": float(m.mean() / c.mean()) if len(m) else None,
               "pearson_r": float(np.corrcoef(c[:len(m)], m)[0, 1]) if len(m) == len(c) and len(m) > 2 else None,
               "region_points_1A": len(inside), "crevice_wall_s": wall_c, "mdpocket_wall_s": wall_m,
               "mdpocket_exit_code": code}
    write_csv(args.results / "trajectory_pocket" / "glut1_volumes.csv", per_frame,
              comment="GLUT1 excerpt (21 frames): CREVICE region-trajectory volume vs MDpocket pock_volume in the "
                      "same reference region. Different volume definitions; compare trends, not magnitudes.")
    write_csv(args.results / "trajectory_pocket" / "summary.csv", [summary])
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

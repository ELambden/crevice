"""Comparator workers for ``timing_trajectory.py`` (run as separate, timed processes).

Usage::

    python trajectory_worker.py hole TOPOLOGY TRAJECTORY --select SEL --cpoint X,Y,Z --cvect X,Y,Z
        --stop N --hole PATH --radii PATH --out per_frame.csv
    python trajectory_worker.py mdpocket-input TOPOLOGY TRAJECTORY --select SEL --obstacles SEL
        --stop N --out-dir DIR

``hole`` runs HOLE 2.3.1 on every frame through mdahole2 (``HoleAnalysis``; one
PDB written and one HOLE process per frame, as mdahole2 does). Each frame is
first fitted on the CA atoms of ``--select`` to the first frame (as CREVICE
does), so the channel point ``--cpoint`` and direction ``--cvect`` are given
in first-frame coordinates. Writes the minimum radius of each frame. Needs an
environment with MDAnalysis and mdahole2.

``mdpocket-input`` writes the input MDpocket reads: the first frame of the
``--obstacles`` atoms as ``top.pdb`` and the frames fitted on the same CA atoms
as ``traj.xtc`` (MDpocket expects superimposed frames), and the frame count
as ``frames_written.txt``. Needs MDAnalysis.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
import tempfile
from pathlib import Path


def _vector(text: str) -> list[float]:
    values = [float(v) for v in text.split(",")]
    if len(values) != 3:
        raise argparse.ArgumentTypeError(f"expected x,y,z, got {text!r}")
    return values


def _fitted_universe(topology: str, trajectory: str, select: str):
    """Universe whose frames are fitted on the CA atoms of ``select`` to the first frame (on the fly)."""
    import MDAnalysis as mda
    from MDAnalysis import transformations

    u = mda.Universe(topology, trajectory)
    ref = mda.Universe(topology, trajectory)
    ref.trajectory[0]
    fit = f"({select}) and name CA"
    u.trajectory.add_transformations(transformations.fit_rot_trans(u.select_atoms(fit), ref.select_atoms(fit)))
    return u


def run_hole(args) -> int:
    from mdahole2.analysis import HoleAnalysis

    u = _fitted_universe(args.topology, args.trajectory, args.select)
    work = Path(tempfile.mkdtemp(prefix="hole_", dir=args.tmp or None))
    rows = []
    with HoleAnalysis(u, select=args.select, executable=args.hole, vdwradii_file=args.radii, cpoint=args.cpoint,
                      cvect=args.cvect, sample=args.sample, end_radius=args.end_radius, tmpdir=str(work)) as hole:
        hole.run(stop=args.stop)
        for frame, profile in hole.results.profiles.items():
            radii = profile.radius
            rows.append({"frame": frame, "points": len(radii), "min_radius_A": float(radii.min()) if len(radii) else ""})
    shutil.rmtree(work, ignore_errors=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["frame", "points", "min_radius_A"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"hole frames={len(rows)}")
    return 0


def run_mdpocket_input(args) -> int:
    import MDAnalysis as mda

    u = _fitted_universe(args.topology, args.trajectory, args.select)
    atoms = u.select_atoms(args.obstacles)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    u.trajectory[0]
    atoms.write(str(out / "top.pdb"))
    count = 0
    with mda.Writer(str(out / "traj.xtc"), atoms.n_atoms) as writer:
        for _ in u.trajectory[:args.stop]:
            writer.write(atoms)
            count += 1
    (out / "frames_written.txt").write_text(f"{count}\n")
    print(f"mdpocket input frames={count} atoms={atoms.n_atoms}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="task", required=True)
    for name in ("hole", "mdpocket-input"):
        p = sub.add_parser(name)
        p.add_argument("topology")
        p.add_argument("trajectory")
        p.add_argument("--select", default="protein", help="analysed atoms; their CA atoms define the fit")
        p.add_argument("--stop", type=int, default=None, help="stop before this frame (default: all)")
    hole = sub.choices["hole"]
    hole.add_argument("--cpoint", type=_vector, required=True)
    hole.add_argument("--cvect", type=_vector, required=True)
    hole.add_argument("--hole", required=True, help="HOLE executable")
    hole.add_argument("--radii", default=None, help="HOLE radius file (default: mdahole2's simple2.rad)")
    hole.add_argument("--sample", type=float, default=0.2)
    hole.add_argument("--end-radius", type=float, default=22.0)
    hole.add_argument("--tmp", default=None, help="directory for HOLE's per-frame files (short path)")
    hole.add_argument("--out", required=True)
    md = sub.choices["mdpocket-input"]
    md.add_argument("--obstacles", default="not resname TIP3 TIP3P SOL WAT HOH")
    md.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    return run_hole(args) if args.task == "hole" else run_mdpocket_input(args)


if __name__ == "__main__":
    sys.exit(main())

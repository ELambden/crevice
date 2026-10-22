"""Run pyKVFinder on one structure and report the cavity that contains a seed.

Called by ``volumes_vs_comparators.py`` with the interpreter of the
environment where pyKVFinder (GPL-3.0-or-later) is installed::

    python kvfinder_worker.py structure.pdb x,y,z out.json [--probe-in 1.4] [--step 0.6]

Atomic radii are replaced by the Bondi element radii (C 1.70, N 1.55,
O 1.52, S 1.80 Å) so that every tool in the volume benchmark uses the same
atoms and radii. The selected cavity is the one whose grid point is nearest
the seed (within 3 Å).
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pyKVFinder

BONDI = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20, "P": 1.80, "F": 1.47}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("seed")
    ap.add_argument("out")
    ap.add_argument("--probe-in", type=float, default=1.4)
    ap.add_argument("--probe-out", type=float, default=4.0)
    ap.add_argument("--step", type=float, default=0.6)
    args = ap.parse_args()
    seed = np.array([float(v) for v in args.seed.split(",")])
    t0 = time.perf_counter()
    atomic = pyKVFinder.read_pdb(args.pdb)
    elements = []
    for line in open(args.pdb):
        if line.startswith(("ATOM", "HETATM")):
            elements.append(line[76:78].strip() or line[12:16].strip()[0])
    atomic[:, -1] = [BONDI[e] for e in elements]
    vertices = pyKVFinder.get_vertices(atomic, probe_out=args.probe_out, step=args.step)
    ncav, cavities = pyKVFinder.detect(atomic, vertices, step=args.step, probe_in=args.probe_in,
                                       probe_out=args.probe_out)
    surface, volume, area = pyKVFinder.spatial(cavities, step=args.step)
    origin = np.asarray(vertices[0], float)
    idx = np.argwhere(cavities > 1)
    record = {"ncav": int(ncav), "probe_in_A": args.probe_in, "probe_out_A": args.probe_out, "step_A": args.step,
              "pyKVFinder_version": pyKVFinder.__version__, "volumes_A3": volume, "areas_A2": area}
    if len(idx):
        # Grid point (i, j, k) lies at origin + step * (i, j, k) when no box rotation is applied.
        pts = origin + args.step * idx
        d = np.linalg.norm(pts - seed, axis=1)
        k = int(np.argmin(d))
        label = int(cavities[tuple(idx[k])])
        # pyKVFinder names cavities KAA, KAB, ... for labels 2, 3, ...
        n = label - 2
        name = "K" + chr(65 + n // 26) + chr(65 + n % 26)
        record.update(selected=name, seed_to_cavity_point_A=float(d[k]), volume_A3=volume.get(name),
                      area_A2=area.get(name), selected_points=int((cavities == label).sum()))
    record["wall_s"] = time.perf_counter() - t0
    json.dump(record, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()

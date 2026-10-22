"""Pocket and cavity volumes: CREVICE vs pyKVFinder and POVME on the pocket set.

Usage (CREVICE environment)::

    python benchmarks/fetch_inputs.py --set pockets
    python benchmarks/volumes_vs_comparators.py [--work-dir DIR] [--only 1L83]

Comparators are optional and located by environment variables:
``CREVICE_KVFINDER_PYTHON`` (an interpreter with pyKVFinder installed) and
``CREVICE_POVME_PYTHON`` + ``CREVICE_POVME_SCRIPT`` (a Python 2.7 interpreter
and POVME 3's ``POVME3.py``; POVME 3.0 is Python 2 only). A missing comparator
is reported, not faked.

Same atoms (polymer heavy atoms, ligand removed), same seed (ligand centroid)
and the same Bondi radii (C 1.70, N 1.55, O 1.52, S 1.80 Å) in every tool.
Each tool keeps its own volume definition, which differs by construction:

``crevice_cast``
    ``crevice cavities --cast --seed`` with defaults (0.8 Å rolling probe, 0.5 Å grid,
    6 Å outer envelope probe): the union of probe spheres swept through the
    seed-connected interior.
``crevice_cavities``
    ``crevice cavities --spacing 0.5`` (1.4 Å probe): grid points where a
    1.4 Å probe *centre* fits in an enclosed void (probe-centre volume); the
    cavity nearest the seed. Open pockets are not cavities here.
``pykvfinder``
    pyKVFinder defaults (probe in 1.4 Å, probe out 4.0 Å, step 0.6 Å): grid
    points inside the probe-out envelope and outside the probe-in solvent-
    excluded surface; the cavity nearest the seed. Also run with probe in
    0.8 Å (CREVICE's cast probe).
``povme``
    POVME 3 with an inclusion sphere of 12 Å at the seed, a 2 Å seed sphere
    (contiguous points), 0.5 Å grid, distance cutoff 1.09 Å (default) and
    convex-hull exclusion: grid points farther than 1.09 Å from every atom
    surface. Bounded by the inclusion sphere by definition.

Output: ``results/volumes/volumes.csv`` (one row per structure and tool).
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (BENCH, INPUTS, RESULTS, crevice_env, default_work_dir, python_for_crevice,  # noqa: E402
                    read_csv, run_logged, write_csv, write_json)
from curate_channels import pattern_atoms, select, write_pdb, xyz  # noqa: E402
from fetch_inputs import structure_file  # noqa: E402

POVME_SPHERE_A = 12.0
GRID_A = 0.5


def curate(row, work):
    src = structure_file(work, row["pdb_id"], row["assembly"])
    st = gemmi.read_structure(str(src))
    model, kept, dropped = select(st, row)
    out = work / "pockets" / row["id"]
    out.mkdir(parents=True, exist_ok=True)
    pdb = out / f"{row['id']}.pdb"
    write_pdb(kept, {c: c[0] for c in dict.fromkeys(c for c, _, _ in kept)}, pdb, f"{row['id']} curated pocket benchmark")
    seed = np.array([xyz(a) for _, _, a in pattern_atoms(model, row["seed_atoms"])]).mean(0)
    write_json(out / "curation.json", {"id": row["id"], "atoms_kept": len(kept), "dropped": dict(dropped),
                                       "seed": seed.tolist()})
    return pdb, seed


def run_cast(pdb, seed, out):
    if out.exists():
        shutil.rmtree(out)
    cmd = [python_for_crevice(), "-m", "crevice", "cavities", "--cast", str(pdb), "--seed=" + ",".join(f"{v:.3f}" for v in seed),
           "--radii", "bondi", "--skip-analysis-views", "--out-dir", str(out), "--prefix", "cast"]
    code, wall = run_logged(cmd, out.parent / "cast.log", env=crevice_env())
    vol = None
    f = out / "cast_void_cast.json"
    if code == 0 and f.exists():
        vol = json.loads(f.read_text())["total_volume"]
    return {"tool": "crevice_cast", "settings": "probe 0.8 A, grid 0.5 A, outer probe 6 A (defaults)",
            "volume_A3": vol, "wall_s": wall, "exit_code": code}


def run_cavities(pdb, seed, out):
    out.mkdir(parents=True, exist_ok=True)
    cmd = [python_for_crevice(), "-m", "crevice", "cavities", str(pdb), "--radii", "bondi", "--spacing", str(GRID_A),
           "-o", str(out / "cavities.json")]
    code, wall = run_logged(cmd, out / "cavities.log", env=crevice_env())
    rec = {"tool": "crevice_cavities", "settings": "probe 1.4 A (centre volume), grid 0.5 A", "wall_s": wall,
           "exit_code": code, "volume_A3": None}
    if code == 0:
        data = json.loads((out / "cavities.json").read_text())
        cavs = data.get("cavities", data if isinstance(data, list) else [])
        best = None
        for c in cavs:
            centre = np.array(c.get("center") or c.get("centre"))
            d = float(np.linalg.norm(centre - seed))
            if best is None or d < best[0]:
                best = (d, c)
        if best is not None:
            rec.update(volume_A3=best[1]["volume"], seed_to_centre_A=best[0],
                       note="enclosed cavity whose centre is nearest the seed")
        else:
            rec["note"] = "no enclosed cavity"
    return rec


def run_kvfinder(pdb, seed, out, probe_in):
    exe = os.environ.get("CREVICE_KVFINDER_PYTHON")
    rec = {"tool": f"pykvfinder_probe_in_{probe_in:g}", "settings": f"probe in {probe_in:g} A, probe out 4 A, step 0.6 A",
           "volume_A3": None}
    if not exe:
        rec["note"] = "pyKVFinder not configured (CREVICE_KVFINDER_PYTHON)"
        return rec
    js = out / f"kvfinder_{probe_in:g}.json"
    code, wall = run_logged([exe, str(BENCH / "kvfinder_worker.py"), str(pdb), ",".join(f"{v:.3f}" for v in seed),
                             str(js), "--probe-in", str(probe_in)], out / "kvfinder.log")
    rec.update(exit_code=code, wall_s=wall)
    if code == 0:
        d = json.loads(js.read_text())
        rec.update(volume_A3=d.get("volume_A3"), seed_to_cavity_point_A=d.get("seed_to_cavity_point_A"),
                   note=f"cavity {d.get('selected')} of {d['ncav']}; pyKVFinder {d['pyKVFinder_version']}")
    return rec


def run_povme(pdb, seed, out):
    exe, script = os.environ.get("CREVICE_POVME_PYTHON"), os.environ.get("CREVICE_POVME_SCRIPT")
    rec = {"tool": "povme", "settings": f"inclusion sphere {POVME_SPHERE_A:g} A, seed sphere 2 A, grid {GRID_A} A, "
                                        "distance cutoff 1.09 A, convex hull", "volume_A3": None}
    if not (exe and script):
        rec["note"] = "POVME not configured (CREVICE_POVME_PYTHON, CREVICE_POVME_SCRIPT)"
        return rec
    run = out / "povme"
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)
    shutil.copy(pdb, run / "in.pdb")
    x, y, z = seed
    (run / "povme.ini").write_text("\n".join([
        f"GridSpacing {GRID_A}", f"InclusionSphere {x:.3f} {y:.3f} {z:.3f} {POVME_SPHERE_A}",
        f"SeedSphere {x:.3f} {y:.3f} {z:.3f} 2.0", "ContiguousPointsCriteria 3", "PDBFileName in.pdb",
        "DistanceCutoff 1.09", "ConvexHullExclusion first", "NumProcessors 1", "OutputFilenamePrefix ./out/povme_",
        "CompressOutput false"]) + "\n")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(script).resolve().parent.parent)
    code, wall = run_logged([exe, script, "povme.ini"], run / "povme.log", cwd=run, env=env)
    rec.update(exit_code=code, wall_s=wall)
    vols = run / "out" / "povme_volumes.tabbed.txt"
    if code == 0 and vols.exists():
        line = [l for l in vols.read_text().splitlines() if l.strip()][-1]
        rec["volume_A3"] = float(line.split()[-1])
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args(argv)
    rows = []
    for row in read_csv(INPUTS / "pockets.csv"):
        if args.only and row["id"] not in args.only:
            continue
        pdb, seed = curate(row, args.work_dir)
        out = args.work_dir / "pockets" / row["id"]
        recs = [run_cast(pdb, seed, out / "cast"), run_cavities(pdb, seed, out / "cavities"),
                run_kvfinder(pdb, seed, out, 1.4), run_kvfinder(pdb, seed, out, 0.8), run_povme(pdb, seed, out)]
        for r in recs:
            rows.append({"id": row["id"], "name": row["name"], "kind": row["kind"], "curation": row["curation"], **r})
            print(rows[-1])
    if not args.only:
        write_csv(args.results / "volumes" / "volumes.csv", rows,
                  fields=["id", "name", "kind", "curation", "tool", "settings", "volume_A3", "seed_to_centre_A",
                          "seed_to_cavity_point_A", "note", "wall_s", "exit_code"],
                  comment="Pocket/cavity volumes by tool; same atoms, seed and Bondi radii; each tool's own volume "
                          "definition (see volumes_vs_comparators.py). Not a validated biological volume.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

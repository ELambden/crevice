"""Water contact counts and region membership checked against MDAnalysis on an MD trajectory.

Implementation cross-checks on real coordinates (the system is not a
validated reference; only the bookkeeping is tested):

A. ``hydration_contacts_<system>.csv``: per frame and residue, the number of
   distinct water oxygens within the contact cutoff (3.5 Å) of any heavy atom
   of the residue, from :func:`crevice.hydration_trajectory.analyze_hydration_trajectory`
   and from an independent MDAnalysis computation
   (``MDAnalysis.lib.distances.capped_distance`` with the periodic box).
   Mismatches are listed with the distance of the closest pair to the cutoff.
B. ``membership_bulk_<system>.csv``: a spherical region of radius R on a
   lattice of spacing h, placed in bulk water so that it straddles a periodic
   face of the box, given to :func:`crevice.water_membership.water_membership`
   frame by frame. Compared with (i) an independent voxel count (minimum image
   about the sphere centre with NumPy, half-up voxel index, same membership
   rule) and (ii) the bulk number density of water oxygens measured with
   MDAnalysis in a slab far from the membrane (members / voxel volume should
   match it within sampling error).

The committed results are for a 100 ns GLUT1 membrane trajectory (CHARMM
naming, 1001 frames at 100 ps; ``--system glut1``); any MD system with
explicit water can be given with ``--topology/--trajectory``.

Run::

    PYTHONPATH=src python benchmarks/nongeom/water_occupancy.py --topology system.gro --trajectory production.xtc
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from nongeom_common import RESULTS, dump, ensure_src_on_path, machine_record, write_csv

ensure_src_on_path()


def reference_contacts(universe, frames, residue_resids, cutoff, protein_sel, water_sel):
    """Distinct water oxygens within ``cutoff`` of any heavy atom of each residue (MDAnalysis, periodic)."""
    from MDAnalysis.lib.distances import capped_distance
    protein = universe.select_atoms(protein_sel)
    water = universe.select_atoms(water_sel)
    res_of_atom = protein.resids
    columns = {r: k for k, r in enumerate(residue_resids)}
    counts = np.zeros((len(frames), len(residue_resids)))
    margins = []
    for row, f in enumerate(frames):
        ts = universe.trajectory[f]
        pairs, dist = capped_distance(protein.positions, water.positions, cutoff + 0.01, box=ts.dimensions)
        inside = dist <= cutoff
        seen = {}
        for (i, j), d, ok in zip(pairs, dist, inside):
            if ok:
                seen.setdefault(res_of_atom[i], set()).add(j)
        for r, waters in seen.items():
            if r in columns:
                counts[row, columns[r]] = len(waters)
        margins.append(float(np.min(np.abs(dist - cutoff))) if len(dist) else None)
    return counts, margins


def part_a(args):
    import MDAnalysis as mda
    from crevice.hydration_trajectory import analyze_hydration_trajectory
    frames = list(range(args.start, args.stop, args.stride))
    t0 = time.perf_counter()
    result = analyze_hydration_trajectory(args.topology, args.trajectory, start=args.start, stop=args.stop,
                                          stride=args.stride, sasa_samples=32, cutoff=args.cutoff,
                                          analysis_options={"skip_density": True, "skip_interactions": True,
                                                            "skip_views": True})
    crevice_seconds = time.perf_counter() - t0
    arrays = result["arrays"]
    labels = [str(x) for x in arrays["residue_ids"]]
    resids = [int("".join(ch for ch in lab.split(":")[-1] if ch.isdigit() or ch == "-")) for lab in labels]
    assert list(arrays["frame_indices"]) == frames
    u = mda.Universe(args.topology, args.trajectory)
    t0 = time.perf_counter()
    ref, margins = reference_contacts(u, frames, resids, args.cutoff, "protein and not name H*",
                                      f"resname {args.water_resname} and name {args.water_oxygen}")
    ref_seconds = time.perf_counter() - t0
    mine = np.asarray(arrays["water_count"], dtype=float)
    diff = mine - ref
    rows = []
    for k, lab in enumerate(labels):
        rows.append({"residue": lab, "frames": len(frames), "crevice_mean_water_count": float(mine[:, k].mean()),
                     "reference_mean_water_count": float(ref[:, k].mean()),
                     "crevice_occupancy": float(np.mean(mine[:, k] > 0)), "reference_occupancy": float(np.mean(ref[:, k] > 0)),
                     "frames_differing": int(np.sum(diff[:, k] != 0)), "max_abs_difference": float(np.abs(diff[:, k]).max())})
    summary = {"frames": len(frames), "frame_range": [args.start, args.stop, args.stride], "residues": len(labels),
               "residue_frame_pairs": int(mine.size), "pairs_differing": int(np.sum(diff != 0)),
               "max_abs_difference": float(np.abs(diff).max()), "total_contacts_crevice": float(mine.sum()),
               "total_contacts_reference": float(ref.sum()),
               "closest_pair_to_cutoff_A": min(m for m in margins if m is not None),
               "crevice_seconds": crevice_seconds, "reference_seconds": ref_seconds,
               "water_population": result["water_population"]}
    print(f"hydration contacts: {summary}", flush=True)
    return rows, summary


def part_b(args):
    import MDAnalysis as mda
    from crevice.cavity_trajectory import CavityReference
    from crevice.models import ChannelPoint, VoidComponent
    from crevice.water_membership import water_membership
    u = mda.Universe(args.topology, args.trajectory)
    oxygens = u.select_atoms(f"resname {args.water_resname} and name {args.water_oxygen}")
    membrane = u.select_atoms(args.membrane_selection)
    ids = oxygens.resindices.astype(np.int64)
    u.trajectory[0]
    box = u.dimensions[:3].astype(float)
    zmem = float(np.median(membrane.positions[:, 2]))
    # Bulk point: half a box away from the membrane centre along z, on the x periodic face.
    centre = np.array([0.0, box[1] / 2, (zmem + box[2] / 2) % box[2]])
    h, radius = args.spacing, args.radius
    r = np.arange(-int(np.ceil(radius / h)), int(np.ceil(radius / h)) + 1) * h
    lattice = np.array(np.meshgrid(r, r, r, indexing="ij")).reshape(3, -1).T
    occupied = lattice[np.linalg.norm(lattice, axis=1) <= radius + 1e-9] + centre
    axial = np.arange(int(np.floor((occupied[:, 2] - centre[2]).min() / h)) - 2,
                      int(np.ceil((occupied[:, 2] - centre[2]).max() / h)) + 3)
    reference = CavityReference(occupied, np.eye(3), centre, h, axial, None, None, None)
    points = tuple(ChannelPoint(i, tuple(map(float, p)), 1.0, 1.0, float(p[2])) for i, p in enumerate(occupied))
    component = VoidComponent(1, "bulk_sphere", tuple(centre), float(len(points)), points)
    identity = {"rotation_rows": np.eye(3).tolist(), "translation_A": [0.0, 0.0, 0.0]}
    volume = len(occupied) * h ** 3
    rows, densities = [], []
    frames = list(range(args.start, args.stop, args.membership_stride))
    for f in frames:
        ts = u.trajectory[f]
        x = oxygens.positions.astype(float)
        got = water_membership(x, ids, ts.dimensions.astype(float), identity, reference, component, region_margin=0.0,
                               geometry_mode="rolling", tolerances=())
        # Independent voxel count: minimum image about the centre (orthorhombic), half-up voxel index.
        L = ts.dimensions[:3].astype(float)
        d = x - centre
        d -= L * np.round(d / L)
        idx = np.floor(d / h + 0.5)
        ref_members = np.flatnonzero(np.linalg.norm(idx * h, axis=1) <= radius + 1e-9)
        same = sorted(ids[ref_members].tolist()) == sorted(np.asarray(got["member_water_ids"]).tolist())
        # Bulk density from a slab of thickness 10 A centred on the same z, MDAnalysis positions.
        dz = x[:, 2] - centre[2]
        dz -= L[2] * np.round(dz / L[2])
        slab = np.abs(dz) <= 5.0
        rho = slab.sum() / (L[0] * L[1] * 10.0)
        densities.append(rho)
        rows.append({"frame": f, "time_ps": float(ts.time), "crevice_members": got["cavity_member_count"],
                     "reference_voxel_members": len(ref_members), "identical_member_ids": same,
                     "region_volume_A3": volume, "crevice_density_per_A3": got["cavity_member_count"] / volume,
                     "slab_bulk_density_per_A3": rho, "box_A": L.tolist()})
    members = np.asarray([r["crevice_members"] for r in rows], float)
    rho = float(np.mean(densities))
    summary = {"frames": len(rows), "centre_A": centre.tolist(), "radius_A": radius, "spacing_A": h,
               "voxels": len(occupied), "volume_A3": volume, "sphere_volume_A3": 4 / 3 * np.pi * radius ** 3,
               "straddles_periodic_face": bool(centre[0] - radius < 0), "membrane_median_z_A": zmem,
               "all_member_ids_identical": all(r["identical_member_ids"] for r in rows),
               "mean_members": float(members.mean()), "sd_members": float(members.std(ddof=1)),
               "mean_member_density_per_A3": float(members.mean() / volume),
               "mean_slab_density_per_A3": rho, "density_ratio": float(members.mean() / volume / rho),
               "expected_members_from_slab_density": rho * volume,
               "note": "frames are correlated; the ratio's sampling error is roughly sd/sqrt(frames) / mean"}
    print(f"membership: {summary}", flush=True)
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--topology", help="MD topology (GRO, TPR, PSF, ...); needed for the trajectory parts")
    parser.add_argument("--trajectory", help="MD trajectory; needed for the trajectory parts")
    parser.add_argument("--system", default="glut1", help="label used in the result file names")
    parser.add_argument("--water-resname", default="TIP3")
    parser.add_argument("--water-oxygen", default="OH2")
    parser.add_argument("--membrane-selection", default="resname DOPC DOPE DOPS POPC DOPG CHL1 and name P")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=1001)
    parser.add_argument("--stride", type=int, default=50, help="frame stride of the hydration comparison")
    parser.add_argument("--membership-stride", type=int, default=10)
    parser.add_argument("--cutoff", type=float, default=3.5)
    parser.add_argument("--radius", type=float, default=8.0)
    parser.add_argument("--spacing", type=float, default=0.5)
    parser.add_argument("--parts", default="ab")
    parser.add_argument("--results", default=str(RESULTS))
    args = parser.parse_args()
    if not (args.topology and args.trajectory):
        parser.error("--topology and --trajectory are required")
    out = Path(args.results)
    summary = {"script": "benchmarks/nongeom/water_occupancy.py", "machine": machine_record(),
               "inputs": {"topology": Path(args.topology).name, "trajectory": Path(args.trajectory).name}}
    if "a" in args.parts:
        rows, s = part_a(args)
        write_csv(out / f"hydration_contacts_{args.system}.csv", rows,
                  comment="Per-residue water-oxygen contacts (3.5 A of any heavy atom): CREVICE vs MDAnalysis "
                          "capped_distance with the periodic box. Implementation cross-check on an MD trajectory "
                          "that is not a validated reference system.")
        summary["hydration_contacts"] = s
    if "b" in args.parts:
        rows, s = part_b(args)
        write_csv(out / f"membership_bulk_{args.system}.csv", rows,
                  comment="Water membership of a lattice sphere in bulk water straddling a periodic face: CREVICE vs "
                          "an independent voxel count and vs the bulk oxygen density from a slab. Implementation "
                          "cross-check; not a validated reference system.")
        summary["membership_bulk"] = s
    dump(out / f"water_occupancy_summary_{args.system}.json", summary)


if __name__ == "__main__":
    main()

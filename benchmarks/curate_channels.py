"""Curate the channel set: exact wall atoms, axis and seed for every entry.

Usage (CREVICE environment; needs gemmi and numpy)::

    python benchmarks/curate_channels.py [--work-dir DIR] [--only 1GRM 3UKM]

For each row of ``inputs/channels.csv`` this writes, in
``<work-dir>/channels/<id>/``:

``curated.pdb``
    The wall atoms both tools read: polymer-entity heavy atoms of model 1
    (altloc blank or A) of the listed chains and residues, as ATOM records
    with element columns and single-letter chain names.
``axis.json``
    Axis direction, seed point, how they were derived, what was removed, the
    seed clearance with the benchmark radius set, a check of the assembly
    file against the deposited entry and the angle/offset between the axis
    and RCSB's symmetry axis for the assembly (an independent cross-check).

and the summary table ``results/channels/curation.csv``. Cached inputs are
never modified.
"""
from __future__ import annotations

import argparse
import json
import re
import string
import sys
from collections import Counter
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INPUTS, REPO, RESULTS, default_work_dir, read_csv, sha256, write_csv, write_json  # noqa: E402
from fetch_inputs import structure_file  # noqa: E402

sys.path.insert(0, str(REPO / "src"))
from crevice.radii import load_radius_file  # noqa: E402

HOLE_SIMPLE_RAD = "simple.rad"
#: Seeds with less clearance (Å, benchmark radii) are moved by the recorded rule in ``curate``.
MIN_SEED_CLEARANCE = 1.0


def find_simple_rad(name: str = "simple") -> Path:
    """A radius file of the HOLE distribution (``simple`` or ``bondi``)."""
    from common import find_tool
    hole = find_tool("hole")
    candidates = []
    if hole is not None:
        candidates += [hole.parent.parent / f"share/hole2/rad/{name}.rad", hole.parent.parent / f"rad/{name}.rad"]
    for c in candidates:
        if c.exists():
            return c
    raise SystemExit(f"HOLE {name}.rad not found; install HOLE first (install_comparators.sh)")


# ----------------------------------------------------------------------------- selection
def residue_matches(spec: str, residue: gemmi.Residue) -> bool:
    if spec == "*":
        return True
    if re.fullmatch(r"-?\d+-\d+", spec):
        lo, hi = spec.rsplit("-", 1)
        return int(lo) <= residue.seqid.num <= int(hi)
    if re.fullmatch(r"-?\d+", spec):
        return residue.seqid.num == int(spec)
    return residue.name == spec


def pattern_atoms(model: gemmi.Model, pattern: str):
    """Atoms matching ``chain:residue:atom`` patterns separated by ';' (any entity)."""
    found = []
    for item in pattern.split(";"):
        chain_spec, res_spec, atom_spec = item.split(":")
        for chain in model:
            if chain_spec not in ("*", chain.name):
                continue
            for residue in chain:
                if not residue_matches(res_spec, residue):
                    continue
                for atom in residue:
                    if atom_spec in ("*", atom.name) and atom.altloc in ("\0", "A"):
                        found.append((chain.name, residue, atom))
    return found


def select(structure: gemmi.Structure, row: dict):
    structure.setup_entities()
    model = structure[0]
    chains = None if row["keep_chains"] == "all" else row["keep_chains"].split(";")
    lo = hi = None
    if row["keep_residues"] != "all":
        lo, hi = (int(v) for v in row["keep_residues"].split("-"))
    kept, dropped = [], Counter()
    for chain in model:
        for residue in chain:
            for atom in residue:
                if residue.entity_type != gemmi.EntityType.Polymer:
                    dropped[f"{residue.name} ({residue.entity_type.name.lower()})"] += 1
                elif chains is not None and chain.name not in chains:
                    dropped[f"chain {chain.name}"] += 1
                elif lo is not None and not (lo <= residue.seqid.num <= hi):
                    dropped[f"residues outside {lo}-{hi}"] += 1
                elif atom.is_hydrogen():
                    dropped["hydrogen"] += 1
                elif atom.altloc not in ("\0", "A"):
                    dropped["altloc not blank/A"] += 1
                else:
                    kept.append((chain.name, residue, atom))
    return model, kept, dropped


def write_pdb(atoms, chain_map: dict, path: Path, title: str) -> None:
    lines = [f"REMARK   1 {title}"[:80]]
    for serial, (chain, residue, atom) in enumerate(atoms, 1):
        name = atom.name
        name_field = f" {name:<3s}" if len(name) < 4 else name[:4]
        lines.append(
            "ATOM  {serial:5d} {name} {res:>3s} {chain:1s}{seq:4d}{icode:1s}   "
            "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {el:>2s}".format(
                serial=serial % 100000, name=name_field, res=residue.name[:3], chain=chain_map[chain],
                seq=residue.seqid.num, icode=residue.seqid.icode.strip() or " ",
                x=atom.pos.x, y=atom.pos.y, z=atom.pos.z, occ=1.0, b=atom.b_iso, el=atom.element.name.upper()))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------- geometry
def kabsch(mobile: np.ndarray, target: np.ndarray):
    mc, tc = mobile.mean(0), target.mean(0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1, 1, d]) @ u.T
    t = tc - rot @ mc
    rmsd = float(np.sqrt(((mobile @ rot.T + t - target) ** 2).sum(1).mean()))
    return rot, t, rmsd


def rotation_axis(rot: np.ndarray, t: np.ndarray):
    """Direction, angle (deg), point nearest the origin and screw (Å) of x -> R x + t."""
    w, v = np.linalg.eig(rot)
    direction = np.real(v[:, np.argmin(np.abs(w - 1))])
    direction /= np.linalg.norm(direction)
    angle = float(np.degrees(np.arccos(np.clip((np.trace(rot) - 1) / 2, -1, 1))))
    t_perp = t - np.dot(t, direction) * direction
    point, *_ = np.linalg.lstsq(np.eye(3) - rot, t_perp, rcond=None)
    point -= np.dot(point, direction) * direction
    return direction, angle, point, float(np.dot(t, direction))


def principal_axis(coords: np.ndarray):
    centred = coords - coords.mean(0)
    evals, evecs = np.linalg.eigh(centred.T @ centred / len(coords))
    return evecs[:, -1] / np.linalg.norm(evecs[:, -1]), evals[::-1]


def xyz(atom) -> list[float]:
    return [atom.pos.x, atom.pos.y, atom.pos.z]


def matched_ca(model, a: str, b: str):
    ca = {(r.seqid.num, r.seqid.icode): at.pos for r in model[a] for at in r if at.name == "CA"}
    cb = {(r.seqid.num, r.seqid.icode): at.pos for r in model[b] for at in r if at.name == "CA"}
    keys = sorted(set(ca) & set(cb))
    return (np.array([[ca[k].x, ca[k].y, ca[k].z] for k in keys]),
            np.array([[cb[k].x, cb[k].y, cb[k].z] for k in keys]))


def rcsb_axis_check(rcsb_assembly: Path, direction: np.ndarray, point: np.ndarray) -> dict:
    """Angle and offset between ``direction`` through ``point`` and each RCSB rotation axis."""
    if not rcsb_assembly.exists():
        return {"rcsb_axis": "metadata not found"}
    data = json.loads(rcsb_assembly.read_text())
    best = None
    for sym in data.get("rcsb_struct_symmetry", []):
        for axis in sym.get("rotation_axes", []) or []:
            start, end = np.array(axis["start"]), np.array(axis["end"])
            rdir = (end - start) / np.linalg.norm(end - start)
            angle = float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(rdir, direction)))))))
            off_vec = (start - point) - np.dot(start - point, direction) * direction
            rec = {"rcsb_symmetry": f"{sym.get('symbol')} ({sym.get('kind')})", "rcsb_axis_angle_deg": angle,
                   "rcsb_axis_offset_A": float(np.linalg.norm(off_vec)), "rcsb_axis_order": axis.get("order")}
            if best is None or angle < best["rcsb_axis_angle_deg"]:
                best = rec
    return best or {"rcsb_axis": "no rotation axis listed"}


def identity_check(assembly_path: Path, deposited_path: Path) -> dict:
    """Compare same-named chains of an assembly file with the deposited entry (max |dxyz|)."""
    if assembly_path == deposited_path:
        return {"identity": "deposited file", "max_abs_coordinate_difference_A": 0.0}
    if not deposited_path.exists():
        return {"identity": "deposited file missing"}
    asm, dep = gemmi.read_structure(str(assembly_path))[0], gemmi.read_structure(str(deposited_path))[0]
    diffs, n = [], 0
    for chain in asm:
        if "-" in chain.name or dep.find_chain(chain.name) is None:
            continue
        other = {(r.seqid.num, r.name, a.name, a.altloc): a.pos for r in dep[chain.name] for a in r}
        for r in chain:
            for a in r:
                p = other.get((r.seqid.num, r.name, a.name, a.altloc))
                if p is not None:
                    n += 1
                    diffs.append(max(abs(p.x - a.pos.x), abs(p.y - a.pos.y), abs(p.z - a.pos.z)))
    return {"identity": "assembly chains compared with the deposited entry", "atoms_compared": n,
            "max_abs_coordinate_difference_A": float(max(diffs)) if diffs else None}


def curate(row: dict, work: Path, radii) -> dict:
    acc = row["id"]
    src = structure_file(work, row["pdb_id"], row["assembly"])
    deposited = structure_file(work, row["pdb_id"], "deposited")
    st = gemmi.read_structure(str(src))
    model, kept, dropped = select(st, row)
    out = work / "channels" / acc
    out.mkdir(parents=True, exist_ok=True)
    names = [c for c in dict.fromkeys(ch for ch, _, _ in kept)]
    letters = string.ascii_uppercase + string.ascii_lowercase + string.digits
    chain_map = {name: (name if len(name) == 1 and name not in names[: names.index(name)] else letters[i])
                 for i, name in enumerate(names)}
    if len(set(chain_map.values())) != len(chain_map):
        chain_map = {name: letters[i] for i, name in enumerate(names)}
    pdb_path = out / "curated.pdb"
    write_pdb(kept, chain_map, pdb_path, f"{acc} curated for the CREVICE benchmark from {src.name}")
    coords = np.array([xyz(a) for _, _, a in kept])
    radius = np.array([radii.radius_for_names(a.name, r.name, a.element.name) for _, r, a in kept])
    centroid = coords.mean(0)

    method = row["axis_method"]
    record: dict = {"id": acc, "pdb_id": row["pdb_id"], "source": src.name, "source_sha256": sha256(src),
                    "assembly": row["assembly"], "atoms_kept": len(kept),
                    "element_counts": dict(Counter(a.element.name for _, _, a in kept)),
                    "dropped": dict(dropped), "chain_map": chain_map, "axis_method": method,
                    "curation": row["curation"]}
    if method == "principal_backbone":
        bb = np.array([xyz(a) for _, _, a in kept if a.name in ("N", "CA", "C", "O")])
        direction, evals = principal_axis(bb)
        record["backbone_covariance_eigenvalues_A2"] = evals.tolist()
        axis_point = bb.mean(0)
    elif method.startswith("cn_superposition"):
        _, a, b = method.split(":")
        mob, tgt = matched_ca(model, a, b)
        rot, t, rmsd = kabsch(mob, tgt)
        direction, angle, axis_point, screw = rotation_axis(rot, t)
        record["pair_operator"] = {"pair": f"{a}->{b}", "matched_CA": len(mob), "rmsd_A": rmsd,
                                   "rotation_deg": angle, "screw_A": screw}
        # Least-squares Cn operator over every chain: order the chains around the
        # provisional axis, then superpose chain i onto chain i+1 for all i at once.
        ring = [c for c in names if model.find_chain(c) is not None]
        cents = {c: np.array([xyz(at) for ch2, _, at in kept if ch2 == c]).mean(0) for c in ring}
        e1 = np.cross(direction, [1.0, 0, 0] if abs(direction[0]) < 0.9 else [0, 1.0, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(direction, e1)
        ang = {c: np.arctan2(np.dot(cents[c] - axis_point, e2), np.dot(cents[c] - axis_point, e1)) for c in ring}
        ring.sort(key=lambda c: ang[c])
        mobs, tgts = [], []
        for i, c in enumerate(ring):
            m_, t_ = matched_ca(model, c, ring[(i + 1) % len(ring)])
            mobs.append(m_)
            tgts.append(t_)
        rot, t, rmsd = kabsch(np.vstack(mobs), np.vstack(tgts))
        direction2, angle, axis_point, screw = rotation_axis(rot, t)
        direction = direction2 if np.dot(direction2, direction) > 0 else -direction2
        record["operator"] = {"method": "least-squares cyclic operator chain i -> i+1 over all kept chains (CA)",
                              "ring_order": ring, "matched_CA": int(sum(len(m_) for m_ in mobs)), "rmsd_A": rmsd,
                              "rotation_deg": angle, "expected_rotation_deg": 360.0 / len(ring), "screw_A": screw}
    elif method == "centroid_plane_normal":
        cents = np.array([np.array([xyz(at) for ch2, _, at in kept if ch2 == ch]).mean(0) for ch in names])
        centre = cents.mean(0)
        _, sv, vt = np.linalg.svd(cents - centre)
        direction = vt[-1] / np.linalg.norm(vt[-1])
        axis_point = centre
        record["chain_centroid_plane_singular_values"] = sv.tolist()
    else:
        raise SystemExit(f"{acc}: unknown axis_method {method}")

    sign = row["sign"]
    if sign == "+z":
        if direction[2] < 0:
            direction = -direction
        record["sign"] = "+z of the deposited frame (arbitrary)"
    else:
        target = pattern_atoms(model, sign.split(":", 1)[1])
        if not target:
            raise SystemExit(f"{acc}: sign pattern matched no atoms")
        tc = np.array([xyz(a) for _, _, a in target]).mean(0)
        if np.dot(tc - centroid, direction) < 0:
            direction = -direction
        record["sign"] = f"+s towards the centroid of {len(target)} atoms {sign.split(':', 1)[1]}"

    origin_method = row["origin_method"]
    if origin_method == "backbone_centroid":
        origin = axis_point
    elif origin_method == "axis_at_centroid":
        origin = axis_point + np.dot(centroid - axis_point, direction) * direction
    elif origin_method.startswith("axis_at_atoms:"):
        target = pattern_atoms(model, origin_method.split(":", 1)[1])
        if not target:
            raise SystemExit(f"{acc}: origin pattern matched no atoms")
        level = np.array([xyz(a) for _, _, a in target]).mean(0)
        origin = axis_point + np.dot(level - axis_point, direction) * direction
        record["origin_atoms"] = [f"{c}:{r.name}{r.seqid.num}:{a.name}" for c, r, a in target]
    elif origin_method.startswith("atoms:"):
        target = pattern_atoms(model, origin_method.split(":", 1)[1])
        if not target:
            raise SystemExit(f"{acc}: origin pattern matched no atoms")
        origin = np.array([xyz(a) for _, _, a in target]).mean(0)
        record["origin_atoms"] = [f"{c}:{r.name}{r.seqid.num}:{a.name}" for c, r, a in target]
    else:
        raise SystemExit(f"{acc}: unknown origin_method {origin_method}")

    def clear(p):
        return float((np.linalg.norm(coords - p, axis=1) - radius).min())

    record["origin_initial"] = origin.tolist()
    record["origin_initial_clearance_A"] = clear(origin)
    if clear(origin) < MIN_SEED_CLEARANCE:
        # Deterministic seed adjustment, recorded: monomer pores search the plane
        # normal to the axis within 3 Å; on-axis seeds search +-5 Å along the axis.
        if origin_method.startswith("atoms:"):
            e1 = np.cross(direction, [1.0, 0, 0] if abs(direction[0]) < 0.9 else [0, 1.0, 0])
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(direction, e1)
            grid = [origin + i * 0.1 * e1 + j * 0.1 * e2 for i in range(-30, 31) for j in range(-30, 31)
                    if i * i + j * j <= 900]
            rule = "max-clearance point within 3 A of the matched-atom centroid in the plane normal to the axis (0.1 A grid)"
        else:
            grid = [origin + k * 0.1 * direction for k in range(-50, 51)]
            rule = "max-clearance point on the axis within +-5 A of the centroid level (0.1 A steps)"
        values = [clear(p) for p in grid]
        origin = grid[int(np.argmax(values))]
        record["origin_adjustment"] = (f"initial seed clearance below {MIN_SEED_CLEARANCE} A; seed moved to the {rule}")

    clearance = np.linalg.norm(coords - origin, axis=1) - radius
    record.update({
        "axis_unit": direction.tolist(), "origin": origin.tolist(), "origin_method": origin_method,
        "origin_clearance_A": float(clearance.min()),
        "axial_extent_of_atoms_A": [float((coords - origin).dot(direction).min()),
                                    float((coords - origin).dot(direction).max())],
        "curated_pdb_sha256": sha256(pdb_path),
    })
    # The seed of a monomer pore is off the symmetry axis: compare the direction only.
    record.update(rcsb_axis_check(work / "rcsb" / f"{row['pdb_id']}_assembly{row['assembly'] if row['assembly'] != 'deposited' else '1'}.json",
                                  direction, axis_point))
    record.update(identity_check(src, deposited))
    write_json(out / "axis.json", record)
    return record


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args(argv)
    rad = find_simple_rad()
    radii = load_radius_file(rad)
    rows = []
    for row in read_csv(INPUTS / "channels.csv"):
        if args.only and row["id"] not in args.only:
            continue
        rec = curate(row, args.work_dir, radii)
        rows.append({"id": rec["id"], "pdb_id": rec["pdb_id"], "assembly": rec["assembly"], "curation": rec["curation"],
                     "atoms_kept": rec["atoms_kept"], "dropped": "; ".join(f"{k} {v}" for k, v in rec["dropped"].items()),
                     "axis_method": rec["axis_method"], "axis_unit": rec["axis_unit"], "origin_A": rec["origin"],
                     "origin_method": rec["origin_method"], "origin_clearance_A": rec["origin_clearance_A"],
                     "operator_rmsd_A": rec.get("operator", {}).get("rmsd_A"),
                     "operator_rotation_deg": rec.get("operator", {}).get("rotation_deg"),
                     "origin_adjustment": rec.get("origin_adjustment", ""),
                     "rcsb_symmetry": rec.get("rcsb_symmetry"), "rcsb_axis_angle_deg": rec.get("rcsb_axis_angle_deg"),
                     "rcsb_axis_offset_A": rec.get("rcsb_axis_offset_A"),
                     "assembly_vs_deposited_max_diff_A": rec.get("max_abs_coordinate_difference_A"),
                     "curated_pdb_sha256": rec["curated_pdb_sha256"][:16]})
        print(f"{rec['id']}: {rec['atoms_kept']} atoms, axis {np.round(rec['axis_unit'], 4).tolist()}, "
              f"seed clearance {rec['origin_clearance_A']:.2f} A, RCSB {rec.get('rcsb_axis_angle_deg')}")
    if not args.only:
        write_csv(args.results / "channels" / "curation.csv", rows,
                  comment="Channel curation summary written by curate_channels.py (radius set: HOLE simple.rad).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

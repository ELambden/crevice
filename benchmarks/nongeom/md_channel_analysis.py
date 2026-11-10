"""Water and ion permeation of a real MD trajectory: CREVICE vs independent references, with definition sensitivity.

For one system described in a JSON file (``inputs/md_systems.json``), this
script

1. defines each pore as an explicit cylinder (axis, a point from the centroid
   of named landmark atoms, an axial extent relative to that point and a
   lateral radius), so that the result does not depend on the automatic
   channel profile;
2. runs :func:`crevice.permeation.analyze_permeation` for the base definition
   and for variants of the radius and of the extent (review item: counts
   depend on the lumen definition);
3. repeats the counting with the independent reference counter of
   ``nongeom_common`` (its own Kabsch fit and minimum image; 'compartment',
   the CREVICE default, for every definition and 'strict' for the base
   definition) and compares the events one by one;
4. for ions, also computes the whole-box displacement current of
   Aksimentiev and Schulten (Biophys. J. 88, 3745 (2005)),
   I(t) = sum_i q_i [z_i(t + dt) - z_i(t)] / (dt L_z), which does not use any
   channel definition, and the current from the ion collective coordinate
   (I_n = e sum_s q_s n_s(T) / T, CREVICE per-frame output);
5. reports pf, pd, pf/pd per pore, the current and conductance, the same in
   consecutive windows (stationarity), and the literature values of
   ``inputs/nongeom_literature.csv`` next to them. Agreement with a
   literature value is reported as a ratio; it is a comparison, not a
   validation of either.

Example (AQP4 pore-water monomer, corrected definition; see ``md_systems/aqp4.md``)::

    PYTHONPATH=src python benchmarks/nongeom/md_channel_analysis.py aqp4_monomer_porewaters_v2 \\
        --topology /path/prod.tpr --trajectory /path/prod.xtc --out-dir /path/analysis

The synthetic smoke test is ``tests/test_nongeom_benchmarks.py``.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

from nongeom_common import (NONGEOM, STATE_BELOW, STATE_INSIDE, STATE_ABOVE, STATE_OUT, dump, ensure_src_on_path,
                            event_keys, machine_record, reference_events, summarize_events, write_csv)

ensure_src_on_path()

E_PER_PS_IN_PA = 1.602176634e-19 / 1e-12 * 1e12
V_W_CM3 = 18.0 / 6.02214076e23


def load_systems(path=None):
    path = Path(path) if path else NONGEOM / "inputs" / "md_systems.json"
    return json.loads(path.read_text())


def literature_rows(system_name, path=None):
    path = Path(path) if path else NONGEOM / "inputs" / "nongeom_literature.csv"
    with open(path, newline="") as handle:
        rows = [r for r in csv.DictReader(line for line in handle if not line.startswith("#"))]
    return [r for r in rows if system_name in r["compare_with"].split(";")]


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def axis_from_reference(universe, spec, reference_path):
    """Reference axis (for example the OPM membrane normal) mapped into the trajectory frame by a CA fit."""
    import MDAnalysis as mda
    if not reference_path:
        raise ValueError("this pore needs --axis-reference (a structure in the frame where 'reference_axis' holds)")
    ref = mda.Universe(str(reference_path))
    mobile = universe.select_atoms(spec["mobile_selection"]).positions.astype(float)
    target = ref.select_atoms(spec["reference_selection"]).positions.astype(float)
    if len(mobile) != len(target):
        raise ValueError(f"axis reference: {len(mobile)} mobile vs {len(target)} reference atoms")
    rot, _, _ = _kabsch(mobile, target)        # mobile @ rot ~ target (row vectors)
    axis = unit(np.asarray(spec.get("reference_axis", [0, 0, 1]), float) @ rot.T)
    return axis if axis[2] >= 0 else -axis


def define_pores(universe, system, axis_reference=None):
    """Axis, point and extent of each pore in the coordinates of the first analysed frame."""
    pores = []
    for pore in system["pores"]:
        point = universe.select_atoms(pore["point_selection"]).positions.mean(axis=0).astype(float)
        if pore.get("axis_from_reference"):
            axis = axis_from_reference(universe, pore["axis_from_reference"], axis_reference)
        elif pore.get("axis_selections"):
            lo = universe.select_atoms(pore["axis_selections"][0]).positions.mean(axis=0)
            hi = universe.select_atoms(pore["axis_selections"][1]).positions.mean(axis=0)
            axis = unit(hi - lo)
        else:
            axis = unit(pore.get("axis", [0, 0, 1]))
        if pore.get("point_offset_A") is not None:
            # Optional shift of the landmark centroid (frame of the first analysed frame), for example onto the
            # pore-water centre line when the landmark atoms sit to one side of the lumen.
            point = point + np.asarray(pore["point_offset_A"], float)
        s0 = float(point @ axis)
        diag = {}
        if pore.get("diagnostic_axis_selections") or pore.get("axis_selections"):
            sels = pore.get("diagnostic_axis_selections") or pore["axis_selections"]
            a = universe.select_atoms(sels[0]).positions.mean(axis=0)
            b = universe.select_atoms(sels[1]).positions.mean(axis=0)
            diag["landmark_axis_angle_to_z_deg"] = float(np.degrees(np.arccos(abs(unit(b - a)[2]))))
        diag["axis_angle_to_z_deg"] = float(np.degrees(np.arccos(abs(axis[2]))))
        pores.append({"name": pore["name"], "axis": axis, "point": point, **diag,
                      "zmin": s0 + pore["extent_relative_A"][0], "zmax": s0 + pore["extent_relative_A"][1],
                      "radius": pore.get("radius_A")})
    return pores


def variants(system):
    """(label, delta extent per end, radius override) for the definition sensitivity."""
    out = [("base", 0.0, None)]
    for dr in system.get("sensitivity", {}).get("radius_A", []):
        out.append((f"radius_{dr:g}A", 0.0, dr))
    for de in system.get("sensitivity", {}).get("extent_delta_A", []):
        out.append((f"extent_{de:+g}A", de, None))
    # 'no_radius' counts every molecule in a box-wide slab; a system may switch it off ("no_radius": false), for
    # example a large multi-pore box where it only adds membrane-wide transits and memory.
    if system.get("sensitivity", {}).get("no_radius", True):
        out.append(("no_radius", 0.0, "none"))
    return out


def channel_objects(pores, variant):
    from crevice.permeation import CylinderChannel
    label, de, radius = variant
    chans = []
    for p in pores:
        r = p["radius"] if radius is None else (None if radius == "none" else radius)
        chans.append(CylinderChannel(f"{p['name']}@{label}", tuple(p["axis"]), p["zmin"] - de, p["zmax"] + de,
                                     point=tuple(p["point"]), radius=r))
    return chans


def _kabsch(mobile, reference):
    xc, yc = mobile.mean(axis=0), reference.mean(axis=0)
    h = (mobile - xc).T @ (reference - yc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt)) or 1.0
    return u @ np.diag([1.0, 1.0, d]) @ vt, xc, yc


def reference_pass(universe, system, pores, var_list, frames, species_atoms, charges):
    """One pass over the trajectory: states per channel variant, unresolved flags, displacement current."""
    fit = universe.select_atoms(system.get("alignment_selection", "protein and name CA"))
    universe.trajectory[frames[0]]
    ref_fit = fit.positions.astype(float).copy()
    idx = np.concatenate(species_atoms)
    q = np.concatenate([np.full(len(a), charges[k]) for k, a in enumerate(species_atoms)])
    # Displacement current over every ion (lab frame, no channel definition).
    ions = np.flatnonzero(q != 0)
    ion_atoms = universe.atoms[idx[ions]]
    disp_charge, times, box_z = [], [], []
    prev_lab_z = None
    # Prefilter: molecules that ever come within (largest radius + 3 A) of a pore axis inside its widest extent + 3 A.
    r_max = max([p["radius"] or 0.0 for p in pores] + [v[2] for v in var_list if isinstance(v[2], float)]) + 3.0
    d_max = max([0.0] + [v[1] for v in var_list]) + 3.0
    # The 'no_radius' variant counts every molecule within the axial extent, so it needs no lateral prefilter.
    lat_limit = np.inf if any(var[2] == "none" for var in var_list) else r_max
    atoms_all = universe.atoms[idx]
    candidate = np.zeros(len(idx), bool)
    for f in frames:
        ts = universe.trajectory[f]
        box = np.asarray(ts.dimensions[:3], float)
        rot, xc, yc = _kabsch(fit.positions.astype(float), ref_fit)
        x = atoms_all.positions.astype(float)
        for p in pores:
            centre = p["point"] + ((p["zmin"] + p["zmax"]) / 2 - p["point"] @ p["axis"]) * p["axis"]
            v = x - ((centre - yc) @ rot.T + xc)
            v -= box * np.round(v / box)
            v = v @ rot
            s = v @ p["axis"]
            lat = np.linalg.norm(v - np.outer(s, p["axis"]), axis=1)
            half = (p["zmax"] - p["zmin"]) / 2 + d_max
            candidate |= (np.abs(s) <= half) & (lat <= lat_limit)
        if prev_lab_z is not None:
            zz = ion_atoms.positions[:, 2].astype(float)
            dz = zz - prev_lab_z
            dz -= box[2] * np.round(dz / box[2])
            disp_charge.append(float(np.sum(q[ions] * dz) / box[2]))   # e moved through the full box height
        prev_lab_z = ion_atoms.positions[:, 2].astype(float).copy()
        times.append(float(ts.time))
        box_z.append(box[2])
    keep = np.flatnonzero(candidate)
    atoms = universe.atoms[idx[keep]]
    defs = []
    for v in var_list:
        label, de, radius = v
        for p in pores:
            r = p["radius"] if radius is None else (None if radius == "none" else radius)
            defs.append({"label": label, "pore": p["name"], "axis": p["axis"], "point": p["point"],
                         "zmin": p["zmin"] - de, "zmax": p["zmax"] + de, "radius": r})
    n_mol = len(keep)
    states = [np.zeros((len(frames), n_mol), dtype=np.int8) for _ in defs]
    # Periodic copy of each channel that a frame's minimum image refers to: whole box heights (z, the membrane
    # normal of these orthorhombic boxes) between the unwrapped and the minimum-image lab position.
    copies = [np.zeros((len(frames), n_mol), dtype=np.int32) for _ in defs]
    unresolved = [[] for _ in defs]
    prev = [None] * len(defs)
    prev_x = None
    x_unwrapped = None
    for k, f in enumerate(frames):
        ts = universe.trajectory[f]
        box = np.asarray(ts.dimensions[:3], float)
        x = atoms.positions.astype(float)
        rot, xc, yc = _kabsch(fit.positions.astype(float), ref_fit)
        if prev_x is not None:
            # Unwrapped step: minimum-image lab displacement since the previous frame, in the reference frame
            # (the minimum image of x - centre jumps by a box length when a molecule crosses the periodic face).
            dx = x - prev_x
            dx -= box * np.round(dx / box)
            x_unwrapped = x_unwrapped + dx
            dx = dx @ rot
        else:
            x_unwrapped = x.copy()
        prev_x = x
        for j, d in enumerate(defs):
            axis = d["axis"]
            centre = d["point"] + ((d["zmin"] + d["zmax"]) / 2 - d["point"] @ axis) * axis
            centre_lab = (centre - yc) @ rot.T + xc
            v = x - centre_lab
            v -= box * np.round(v / box)          # orthorhombic minimum image (lab frame)
            copies[j][k] = np.round(((x_unwrapped - centre_lab) - v)[:, 2] / box[2])
            v = v @ rot                           # into the reference frame
            s = centre @ axis + v @ axis
            lat = np.linalg.norm(v - np.outer(v @ axis, axis), axis=1)
            st = np.where(s < d["zmin"], STATE_BELOW, np.where(s > d["zmax"], STATE_ABOVE, STATE_OUT)).astype(np.int8)
            slab = st == STATE_OUT
            ok = np.ones(n_mol, bool) if d["radius"] is None else lat <= d["radius"]
            st[slab & ok] = STATE_INSIDE
            states[j][k] = st
            if prev[j] is not None:
                p_state, p_s, p_v = prev[j]
                step = dx @ axis
                landed = p_s + step
                up = (p_state == STATE_BELOW) & (landed > d["zmax"])
                dn = (p_state == STATE_ABOVE) & (landed < d["zmin"])
                for i in np.flatnonzero(up | dn):
                    unresolved[j].append({"molecule": int(keep[i]), "direction": "up" if up[i] else "down",
                                          "entry_frame": k - 1, "exit_frame": k, "status": "unresolved_single_interval"})
            prev[j] = (st, s, v)
    return defs, states, copies, unresolved, keep, np.asarray(disp_charge), np.asarray(times), np.asarray(box_z)


def reference_counts(defs, states, copies, unresolved, keep):
    out = []
    for d, st, cp, unres in zip(defs, states, copies, unresolved):
        # Every molecule that entered the axial extent (inside the lumen or beside it): a transit outside the
        # lumen is listed under the compartment rule even for a molecule never in this lumen.
        local = np.flatnonzero(((st == STATE_INSIDE) | (st == STATE_OUT)).any(axis=0))
        ever = keep[local]
        st = st[:, local]
        sub, sub_copies = st, cp[:, local]
        # reference_events needs unwrapped s only for unresolved detection, done online above.
        strict = reference_events(sub, np.zeros(sub.shape), 0.0, 0.0, convention="strict", copies=sub_copies)
        comp = reference_events(sub, np.zeros(sub.shape), 0.0, 0.0, convention="compartment", copies=sub_copies)
        strict = [dict(e, molecule=int(ever[e["molecule"]])) for e in strict if e["status"] == "complete"] + unres
        comp = [dict(e, molecule=int(ever[e["molecule"]])) for e in comp
                if e["status"] in ("complete", "outside_lumen")] + unres
        out.append({"def": d, "strict": strict, "compartment": comp})
    return out


def windows_summary(result, pores_names, n_windows):
    """pf/pd and occupancy in consecutive windows from the per-frame output (stationarity check)."""
    from crevice.permeation import collective_diffusion
    times = np.asarray(result.frames["time_ps"])
    F = len(times)
    edges = np.linspace(0, F, n_windows + 1).astype(int)
    rows = []
    for key, series in result.frames["series"].items():
        n = np.asarray(series["n"])
        occ = np.asarray(series["occupancy"])
        net = np.asarray(series["net_crossings"])
        for w in range(n_windows):
            a, b = edges[w], edges[w + 1]
            if b - a < 20:
                continue
            dt = float(np.median(np.diff(times[a:b])))
            dn = collective_diffusion(n[a:b], dt, blocks=min(5, max(2, (b - a) // 40)))
            rows.append({"series": key, "window": w, "first_frame": int(a), "last_frame": int(b - 1),
                         "t0_ps": float(times[a]), "t1_ps": float(times[b - 1]), "mean_occupancy": float(occ[a:b].mean()),
                         "net_crossings": int(net[b - 1] - (net[a - 1] if a > 0 else 0)),
                         "D_n_per_s": dn.get("D_per_s"), "pf_cm3_per_s": V_W_CM3 * dn["D_per_s"] if dn.get("D_per_s") else None})
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("system", help="key in inputs/md_systems.json")
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--systems", default=None, help="systems JSON (default inputs/md_systems.json)")
    parser.add_argument("--literature", default=None)
    parser.add_argument("--start", type=int, default=None, help="first frame (default: the system's equilibration cut)")
    parser.add_argument("--stop", type=int, default=None)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--windows", type=int, default=5)
    parser.add_argument("--efield", type=float, default=None, help="V/nm along z (default from the system JSON)")
    parser.add_argument("--axis-reference", default=None,
                        help="structure for pores with 'axis_from_reference' (for example the OPM-oriented model)")
    args = parser.parse_args(argv)
    import MDAnalysis as mda
    from crevice.permeation import PermeationResult, analyze_permeation, parse_permeation_tokens, resolve_species
    system = load_systems(args.systems)[args.system]
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    u = mda.Universe(args.topology, args.trajectory)
    start = args.start if args.start is not None else system.get("start_frame", 0)
    frames = list(range(start, len(u.trajectory) if args.stop is None else min(args.stop, len(u.trajectory)), args.stride))
    u.trajectory[frames[0]]
    pores = define_pores(u, system, args.axis_reference)
    var_list = variants(system)
    efield = args.efield if args.efield is not None else system.get("efield_V_per_nm")
    species = system["species"]
    found, absent = resolve_species(u, parse_permeation_tokens(species))
    # CREVICE, one call per variant (all pores of a variant together).
    crevice = {}
    for v in var_list:
        for conv in (("compartment", "strict") if v[0] == "base" else ("compartment",)):
            res = analyze_permeation(u, species=species, channels=channel_objects(pores, v), start=frames[0],
                                     stop=frames[-1] + 1, stride=args.stride,
                                     alignment_selection=system.get("alignment_selection", "protein and name CA"),
                                     efield=[0, 0, efield] if efield else None,
                                     efield_source="md_systems.json" if efield else None, convention=conv)
            if v[0] != "base":
                # Keep only what the comparison needs: the per-lumen-visit rows of every variant of a
                # four-pore box held together exceed 18 GB for 1839 frames.
                res = PermeationResult(summary=res.summary, events=res.events, visits=[], frames={})
            crevice[(v[0], conv)] = res
            if v[0] == "base":
                from crevice.permeation import write_permeation_outputs
                write_permeation_outputs(res, out, prefix=f"{args.system}_{conv}")
    crevice_seconds = time.perf_counter() - t0
    # Independent reference pass.
    t1 = time.perf_counter()
    species_atoms = [s.atom_indices for s in found]
    charges = [s.charge for s in found]
    offsets = np.cumsum([0] + [len(a) for a in species_atoms])
    defs, states, copies, unresolved, keep, disp, times, box_z = reference_pass(u, system, pores, var_list, frames,
                                                                         species_atoms, charges)
    ref = reference_counts(defs, states, copies, unresolved, keep)
    reference_seconds = time.perf_counter() - t1
    duration_ps = float(times[-1] - times[0])
    rows = []
    all_atoms = np.concatenate(species_atoms)
    position = {int(a): i for i, a in enumerate(all_atoms)}
    index = {f: k for k, f in enumerate(frames)}
    for item in ref:
        d = item["def"]
        name = f"{d['pore']}@{d['label']}"
        for conv in ("compartment", "strict"):
            if (d["label"], conv) not in crevice:
                continue
            res = crevice[(d["label"], conv)]
            chan = next(c for c in res.summary["channels"] if c["definition"]["name"] == name)
            keys = sorted((position[int(e["atom_index"])], e["direction"], index[int(e["entry_frame"])],
                           index[int(e["exit_frame"])], e["status"]) for e in res.events if e["channel"] == name)
            for j, sp in enumerate(found):
                lo, hi = offsets[j], offsets[j + 1]
                refc = summarize_events([e for e in item[conv] if lo <= e["molecule"] < hi])
                row = chan["species"][j]
                rec = {"pore": d["pore"], "definition": d["label"], "convention": conv, "zmin_A": d["zmin"],
                       "zmax_A": d["zmax"], "radius_A": d["radius"], "species": sp.label,
                       "crevice_up": row["events_up"], "crevice_down": row["events_down"],
                       "crevice_unresolved": row["unresolved_crossings"],
                       "crevice_outside_lumen": row["transits_outside_lumen"],
                       "reference_up": refc["up"], "reference_down": refc["down"],
                       "reference_unresolved": refc["unresolved"], "reference_outside_lumen": refc["outside_lumen"],
                       "event_dispersion_index": row["event_dispersion"]["total"],
                       "mean_occupancy": row["occupancy"]["mean"]}
                if "water" in row:
                    w = row["water"]
                    rec.update(pf_cm3_per_s=w["pf_cm3_per_s"], pf_interval=w["pf_interval_cm3_per_s"],
                               pd_cm3_per_s=w["pd_cm3_per_s"], pd_interval=w["pd_interval_cm3_per_s"],
                               pd_garwood_interval=w["pd_garwood_interval_cm3_per_s"], pf_over_pd=w["pf_over_pd"])
                if "ion" in row:
                    rec.update(current_pA=row["ion"]["current_pA"])
                rows.append(rec)
            rows[-1]["events_identical_to_reference_all_species"] = keys == event_keys(item[conv])
            # The reference follows only molecules that come near a pore axis (prefilter), so transits outside the
            # lumen and unresolved crossings of molecules far from it are missing from its lists by construction;
            # the complete events must agree exactly.
            rows[-1]["complete_events_identical_to_reference"] = (
                [k for k in keys if k[-1] == "complete"] == [k for k in event_keys(item[conv]) if k[-1] == "complete"])
    write_csv(out / f"{args.system}_permeation_definitions.csv", rows,
              comment="CREVICE vs independent reference counts per pore, species and channel definition.")
    # Currents.
    summary = {"system": args.system, "description": system.get("description"), "frames": len(frames),
               "first_frame": frames[0], "stride": args.stride, "duration_ps": duration_ps,
               "pores": [{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in p.items()} for p in pores],
               "absent_species_tokens": absent, "crevice_seconds": crevice_seconds,
               "reference_seconds": reference_seconds, "machine": machine_record()}
    if any(c != 0 for c in charges):
        disp_current = E_PER_PS_IN_PA * float(np.sum(disp)) / duration_ps
        base = crevice[("base", "compartment")]
        per_pore = [c["ion_current_total"]["current_pA"] for c in base.summary["channels"] if c["ion_current_total"]]
        n_current = 0.0
        for key, series in base.frames["series"].items():
            sp = next(s for s in found if key.endswith("/" + s.label))
            n_current += sp.charge * series["n"][-1]
        n_current = E_PER_PS_IN_PA * n_current / duration_ps
        voltage_mV = 1000 * efield * float(np.mean(box_z)) / 10 if efield else None
        # Block intervals for the displacement current (10 blocks of frames).
        blocks = np.array_split(disp, 10)
        block_i = np.asarray([E_PER_PS_IN_PA * b.sum() / (len(b) * duration_ps / len(disp)) for b in blocks if len(b)])
        se = float(block_i.std(ddof=1) / math.sqrt(len(block_i)))
        summary["currents"] = {
            "displacement_current_pA": disp_current, "displacement_current_block_se_pA": se,
            "crevice_count_current_sum_over_pores_pA": float(sum(per_pore)),
            "crevice_collective_current_sum_over_pores_pA": n_current,
            "voltage_mV_from_field": voltage_mV,
            "conductance_displacement_nS": disp_current / voltage_mV if voltage_mV else None,
            "conductance_count_nS": float(sum(per_pore)) / voltage_mV if voltage_mV else None,
            "definitions": {"displacement": "Aksimentiev & Schulten 2005: sum q dz / (dt Lz) over all ions in the box",
                            "count": "CREVICE complete crossings of each pore extent, summed over pores",
                            "collective": "e sum_s q_s n_s(T) / T over pores (CREVICE frames output)"}}
    summary["windows"] = windows_summary(crevice[("base", "compartment")], [p["name"] for p in pores], args.windows)
    lit = literature_rows(args.system, args.literature)
    comparisons = []
    base_rows = [r for r in rows if r["definition"] == "base" and r["convention"] == "compartment"]
    for entry in lit:
        value = float(entry["value"])
        q = entry["quantity"]
        if q == "pf_per_monomer_cm3_per_s":
            ours = [r["pf_cm3_per_s"] for r in base_rows if r.get("pf_cm3_per_s") is not None]
            mine = float(np.mean(ours)) if ours else None
        elif q == "conductance_per_trimer_nS" and "currents" in summary:
            mine = summary["currents"]["conductance_displacement_nS"]
        else:
            mine = None
        comparisons.append({"quantity": q, "literature_value": value, "literature_uncertainty": entry["uncertainty"],
                            "conditions": entry["conditions"], "source": entry["source"], "doi": entry["doi"],
                            "this_run": mine, "ratio": mine / value if mine else None})
    summary["literature_comparison"] = comparisons
    summary["status"] = ("comparison with literature values only: a single trajectory, sampling intervals only; "
                         "not a validation of the force field or of CREVICE")
    dump(out / f"{args.system}_md_channel_summary.json", summary)
    print(json.dumps({k: summary[k] for k in ("system", "frames", "duration_ps")}, indent=1))
    for r in base_rows:
        print({k: r[k] for k in ("pore", "species", "crevice_up", "crevice_down", "reference_up",
                                 "reference_down", "crevice_outside_lumen") if k in r})
    if "currents" in summary:
        print(summary["currents"])
    return summary


if __name__ == "__main__":
    main()

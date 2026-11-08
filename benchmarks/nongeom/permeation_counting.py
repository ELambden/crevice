"""Permeation event counting and the collective coordinate n(t) against exact answers.

Software verification of :func:`crevice.permeation.analyze_permeation` on
synthetic trajectories whose answers are known by construction or
analytically, plus an event-by-event comparison with an independent reference
counter (``nongeom_common.reference_events``, written from the definition of
Zhu, Tajkhorshid and Schulten 2004 without calling CREVICE).

Parts (all written to ``benchmarks/results/nongeom/``):

A. ``permeation_cases.csv``: a catalogue of scripted single-molecule paths
   (full crossings, entry and exit through the same mouth, a molecule sitting
   in the mouth, re-crossings, periodic loops, lateral exits, a lateral
   approach at the mouth, whole-extent jumps, exact boundary positions, a
   channel straddling the periodic boundary, frame stride, a two-frame
   trajectory). Expected counts are written by hand from the definition.
B. ``permeation_ring.csv``: independent Brownian walkers on a periodic line
   through a channel of length L in a period Lambda. Exact continuum results:
   unidirectional rate q0 = N D / (L Lambda) and collective diffusion
   D_n = N D / Lambda^2 (so pf/pd = L / Lambda). Discrete frames see the
   channel ends shifted outward by beta*sigma (Siegmund's correction,
   beta = 0.5826, sigma = sqrt(2 D dt)), giving q0 ~ N D / ((L + 2 beta sigma) Lambda).
   The frame step is varied (sigma / L from 0.01 to 0.3).
C. ``permeation_single_file.csv``: the same walkers with hard-core
   reflection (single file). D_n is unchanged (exactly N D / Lambda^2),
   q0 drops; the same trajectory is analysed at frame strides 1-20 to show
   that pd depends on the saving interval while pf does not.
D. ``permeation_lumen_sensitivity.csv``: 3D Brownian particles in a box with
   an impenetrable membrane slab pierced by a cylindrical pore. Every
   membrane transit passes through the pore, so the physical transit count is
   known; CREVICE is run with channel definitions that match the pore or
   extend beyond it (axially into the bulk, laterally narrower or wider), and
   CREVICE and the reference counter are compared under both event
   conventions: 'compartment' (the default: bulk side remembered, event
   attributed to the channel if the lumen was visited on the way) and
   'strict' (uninterrupted stay in the lumen; a lateral exit forgets the side).

Run::

    PYTHONPATH=src python benchmarks/nongeom/permeation_counting.py            # full (about 10 min)
    PYTHONPATH=src python benchmarks/nongeom/permeation_counting.py --quick    # smoke test
"""
from __future__ import annotations

import argparse
import math
import time

import numpy as np

from nongeom_common import (RESULTS, SIEGMUND_BETA, crevice_event_keys, dump, ensure_src_on_path, event_keys,
                            expected_window_fit,
                            machine_record, memory_universe, reference_axial_lateral, reference_collective_coordinate,
                            reference_events, reference_states, ring_walkers, single_file_relabel, summarize_events,
                            write_csv)

ensure_src_on_path()
from crevice.permeation import CylinderChannel, analyze_permeation, collective_diffusion  # noqa: E402

N_PROTEIN = 6   # CA ring atoms come first in every synthetic universe


def crevice_counts(result, channel_index=0, species_index=0):
    row = result.summary["channels"][channel_index]["species"][species_index]
    return row


# ---------------------------------------------------------------------------
# A. Scripted catalogue
# ---------------------------------------------------------------------------

def scripted_cases():
    """(name, path [(x, z) per frame], expected strict (up, down, unresolved),
    expected compartment (up, down, outside_lumen), note). Unresolved crossings are the same in both conventions."""
    z = lambda seq, x=0.0: [(x, v) for v in seq]  # noqa: E731
    sitting = [-10] + [-5.5, -4.5] * 10 + [-10]
    loop = list(np.arange(-12, 69, 2.0))
    return [
        ("full_up", z([-10, -8, -4, 0, 4, 8, 10]), (1, 0, 0), (1, 0, 0), "strict entry frame 2, exit frame 5"),
        ("full_down", z([10, 8, 4, 0, -4, -8, -10]), (0, 1, 0), (0, 1, 0), ""),
        ("same_mouth_bottom", z([-10, -6, -3, 0, -3, -6, -10]), (0, 0, 0), (0, 0, 0), "enters and leaves below"),
        ("same_mouth_top", z([10, 6, 3, 0, 3, 6, 10]), (0, 0, 0), (0, 0, 0), "enters and leaves above"),
        ("sits_in_mouth", z(sitting), (0, 0, 0), (0, 0, 0), "10 entries through the lower mouth, never crosses"),
        ("recrossing_up_then_down", z([-10, -4, 0, 4, 8, 4, 0, -4, -8]), (1, 1, 0), (1, 1, 0),
         "exits above at frame 4 and comes straight back"),
        ("periodic_loop_up_twice", z(loop), (2, 0, 0), (2, 0, 0), "box 40 A: passes the channel image twice"),
        ("bulk_periodic_loop_outside_radius", z(loop, x=10.0), (0, 0, 0), (0, 0, 2),
         "passes the axial extent twice, never in the lumen: two outside_lumen transits in the compartment rule"),
        ("started_inside", z([0, 2, 4, 7, 9]), (0, 0, 0), (0, 0, 0), "censored at the start"),
        ("inside_at_end", z([-10, -6, -2, 0]), (0, 0, 0), (0, 0, 0), "censored at the end"),
        ("lateral_exit_mid_channel", [(0, -10), (0, -7), (0, -3), (8, -2), (0, 2), (0, 7), (0, 10)], (0, 0, 0),
         (1, 0, 0), "inside -> out (x = 8 A) -> inside -> above"),
        ("lateral_approach_at_mouth", [(8, -12), (8, -4), (0, -3), (0, 0), (0, 3), (0, 8)], (0, 0, 0), (1, 0, 0),
         "enters the axial extent outside the radius, then moves into the lumen and crosses"),
        ("whole_extent_jump", z([-8, -8, 8, 8]), (0, 0, 1), (0, 0, 0), "crosses 16 A between two frames (unresolved)"),
        ("exact_boundaries", z([-10, -5.0, 0, 5.0, 10]), (1, 0, 0), (1, 0, 0), "z = zmin and z = zmax count as inside"),
        ("down_then_lateral_exit", [(0, 10), (0, 3), (0, 0), (8, 0), (8, -8), (8, -10)], (0, 0, 0), (0, 1, 0),
         "enters from above, leaves sideways inside the extent, reaches below"),
        ("passes_extent_beside_lumen", z([-10, -4, 0, 4, 10], x=8.0), (0, 0, 0), (0, 0, 1),
         "crosses the axial extent 8 A from the axis (for example through another pore)"),
    ]


def run_cases():
    cases = scripted_cases()
    frames = max(len(p) for _, p, *_ in cases)
    pos = np.zeros((frames, len(cases), 3))
    for j, (_, path, *_rest) in enumerate(cases):
        for t in range(frames):
            x, zz = path[min(t, len(path) - 1)]
            pos[t, j] = (x, 0.0, zz)
    box = (30.0, 30.0, 40.0)
    u = memory_universe(["TIP3"] * len(cases), pos, box=box, dt=10.0)
    channel = CylinderChannel("pore", (0, 0, 1), 15.0, 25.0, point=(15.0, 15.0, 20.0), radius=3.0)
    res = {conv: analyze_permeation(u, species=["TIP3"], channels=[channel], msd_blocks=2, convention=conv)
           for conv in ("compartment", "strict")}
    s, img, lat, _ = reference_axial_lateral(u, np.arange(N_PROTEIN, N_PROTEIN + len(cases)), (15, 15, 20))
    states = reference_states(img, lat, 15.0, 25.0, 3.0)
    ref = {conv: reference_events(states, s, 15.0, 25.0, convention=conv, s_image=img)
           for conv in ("compartment", "strict")}
    rows = []
    for j, (name, path, strict, comp, note) in enumerate(cases):
        got = {}
        for conv in ("compartment", "strict"):
            mine = [e for e in res[conv].events if int(e["atom_index"]) == N_PROTEIN + j]
            got[conv] = summarize_events([{"direction": e["direction"], "status": e["status"]} for e in mine])
            got["ref_" + conv] = summarize_events([e for e in ref[conv] if e["molecule"] == j])
        cs, cc, rs, rc = got["strict"], got["compartment"], got["ref_strict"], got["ref_compartment"]
        transit = [e["transit_time_ps"] for e in res["compartment"].events
                   if int(e["atom_index"]) == N_PROTEIN + j and e["status"] == "complete"]
        rows.append({"case": name, "frames": len(path),
                     "expected_strict": list(strict), "crevice_strict": [cs["up"], cs["down"], cs["unresolved"]],
                     "reference_strict": [rs["up"], rs["down"], rs["unresolved"]],
                     "expected_compartment": list(comp),
                     "crevice_compartment": [cc["up"], cc["down"], cc["outside_lumen"]],
                     "reference_compartment": [rc["up"], rc["down"], rc["outside_lumen"]],
                     "crevice_unresolved_compartment": cc["unresolved"],
                     "crevice_compartment_transit_times_ps": transit,
                     "crevice_matches_expected": (cs["up"], cs["down"], cs["unresolved"]) == strict
                     and (cc["up"], cc["down"], cc["outside_lumen"]) == comp and cc["unresolved"] == strict[2],
                     "reference_matches_expected": (rs["up"], rs["down"], rs["unresolved"]) == strict
                     and (rc["up"], rc["down"], rc["outside_lumen"]) == comp, "note": note})
    extra = []
    # Channel straddling the periodic boundary: the channel centre is on the box face z = 0.
    zs = np.arange(-10, 10.1, 2.0)
    pos2 = np.zeros((len(zs), 1, 3))
    pos2[:, 0, 2] = zs
    u2 = memory_universe(["TIP3"], pos2, box=box, dt=10.0, centre=(15.0, 15.0, 0.0))
    ch2 = CylinderChannel("edge", (0, 0, 1), -5.0, 5.0, point=(15.0, 15.0, 0.0), radius=3.0)
    r2 = crevice_counts(analyze_permeation(u2, species=["TIP3"], channels=[ch2], msd_blocks=2, align=False))
    extra.append({"case": "channel_straddles_periodic_face", "frames": len(zs), "expected_up": 1, "expected_down": 0,
                  "expected_unresolved": 0, "crevice_up": r2["events_up"], "crevice_down": r2["events_down"],
                  "crevice_unresolved": r2["unresolved_crossings"],
                  "crevice_matches_expected": (r2["events_up"], r2["events_down"], r2["unresolved_crossings"]) == (1, 0, 0),
                  "note": "coordinates wrapped by the box; channel centred on z = 0 (alignment off)"})
    # Frame stride on the full_up path: stride 3 keeps one inside frame, stride 6 none.
    full = np.zeros((7, 1, 3))
    full[:, 0, 2] = [-8, -6, -3, 0, 3, 6, 8]   # stride 6 jumps 16 A (< half the 40 A box)
    u3 = memory_universe(["TIP3"], full, box=box, dt=10.0)
    for stride, expected in ((1, (1, 0, 0)), (3, (1, 0, 0)), (6, (0, 0, 1))):
        r3 = crevice_counts(analyze_permeation(u3, species=["TIP3"], channels=[channel], msd_blocks=2, stride=stride))
        got = (r3["events_up"], r3["events_down"], r3["unresolved_crossings"])
        extra.append({"case": f"full_up_stride_{stride}", "frames": len(range(0, 7, stride)), "expected_up": expected[0],
                      "expected_down": expected[1], "expected_unresolved": expected[2], "crevice_up": got[0],
                      "crevice_down": got[1], "crevice_unresolved": got[2], "crevice_matches_expected": got == expected,
                      "note": "frame stride: an event needs one saved frame inside the lumen"})
    # Two-frame trajectory: counting works, D_n is reported as insufficient.
    two = np.zeros((2, 1, 3))
    two[:, 0, 2] = [-10, 0]
    res4 = analyze_permeation(memory_universe(["TIP3"], two, box=box, dt=10.0), species=["TIP3"], channels=[channel],
                              msd_blocks=2)
    r4 = crevice_counts(res4)
    extra.append({"case": "two_frame_trajectory", "frames": 2, "expected_up": 0, "expected_down": 0,
                  "expected_unresolved": 0, "crevice_up": r4["events_up"], "crevice_down": r4["events_down"],
                  "crevice_unresolved": r4["unresolved_crossings"],
                  "crevice_matches_expected": (r4["events_up"], r4["events_down"]) == (0, 0)
                  and r4["water"]["collective_diffusion"]["status"] == "insufficient_frames",
                  "note": "D_n status " + r4["water"]["collective_diffusion"]["status"]})
    # n(t) when one frame step wraps through the periodic face opposite the channel and lands inside it
    # (period 20 A, channel 5..15 A): the 0.5 A travelled inside the channel image is not added.
    wrap = np.array([[9.5], [15.5], [18.0], [22.0]])
    u5 = ring_universe(wrap, 20.0)
    ch5 = CylinderChannel("pore", (0, 0, 1), 5.0, 15.0, point=(15.0, 15.0, 10.0), radius=3.0)
    r5 = analyze_permeation(u5, species=["TIP3"], channels=[ch5], align=False, msd_blocks=2)
    s5, _, _, _ = reference_axial_lateral(u5, [N_PROTEIN], (15.0, 15.0, 10.0))
    n_cv = r5.frames["series"]["pore/TIP3"]["n"]
    n_ref = reference_collective_coordinate(s5, 5.0, 15.0, 20.0).tolist()
    extra.append({"case": "n_entry_through_opposite_periodic_face", "frames": 4, "crevice_n": n_cv,
                  "reference_n": n_ref, "crevice_matches_expected": bool(np.allclose(n_cv, n_ref)),
                  "note": "expected n = 0, 0.05, 0.3, 0.7; CREVICE clips the step in the image of the previous "
                          "position, so a step that wraps into the channel loses the part inside (bounded, < 1 per "
                          "such step; needs a frame step longer than half the bulk path, (box - L) / 2)"})
    return rows + extra


# ---------------------------------------------------------------------------
# B / C. Ring walkers (independent and single file)
# ---------------------------------------------------------------------------

def ring_universe(z_unwrapped, box_z):
    frames, n = z_unwrapped.shape
    pos = np.zeros((frames, n, 3))
    pos[:, :, 2] = z_unwrapped
    return memory_universe(["TIP3"] * n, pos, box=(30.0, 30.0, box_z), dt=1.0)


def analyse_ring(z_unwrapped, box_z, half_length, *, stride=1, msd_max_lag=None):
    """CREVICE and reference results for walkers at x = y = 0 (always laterally inside)."""
    u = ring_universe(z_unwrapped, box_z)
    c = box_z / 2
    zmin, zmax = c - half_length, c + half_length
    channel = CylinderChannel("pore", (0, 0, 1), zmin, zmax, point=(15.0, 15.0, c), radius=3.0)
    t0 = time.perf_counter()
    result = analyze_permeation(u, species=["TIP3"], channels=[channel], stride=stride, msd_blocks=5,
                                msd_max_lag=msd_max_lag, align=False)
    seconds = time.perf_counter() - t0
    n = z_unwrapped.shape[1]
    s, img, lat, _ = reference_axial_lateral(u, np.arange(N_PROTEIN, N_PROTEIN + n), (15, 15, c), stride=stride)
    states = reference_states(img, lat, zmin, zmax, 3.0)
    ref = reference_events(states, s, zmin, zmax, convention="compartment", s_image=img)
    n_ref = reference_collective_coordinate(s, zmin, zmax, box_z)
    n_cv = np.asarray(result.frames["series"]["pore/TIP3"]["n"])
    frame_numbers = list(range(0, z_unwrapped.shape[0], stride))
    mine, theirs = set(crevice_event_keys(result, N_PROTEIN, frame_numbers)), set(event_keys(ref))
    same = mine == theirs
    # Differences caused by a molecule exactly half a period from the channel centre (float32 tie of the
    # minimum image: the two programs may round it to opposite images) are flagged, not hidden.
    ties = 0
    for i, _, entry, exit_, _ in mine ^ theirs:
        window = img[max(0, entry - 3):exit_ + 1, i]
        ties += bool(np.any(np.isclose(np.abs(window - c), box_z / 2, rtol=0, atol=1e-4)))
    return result, {"seconds": seconds, "events_identical_to_reference": same,
                    "differing_events": len(mine ^ theirs), "differing_events_at_half_period_ties": ties,
                    "n_max_abs_difference": float(np.max(np.abs(n_cv - n_ref))),
                    "reference": summarize_events(ref)}


def run_ring(quick, rng):
    box_z, half = 20.0, 5.0
    L = 2 * half
    n, frames = (40, 4000) if quick else (100, 20000)
    sigmas = (0.5, 2.0) if quick else (0.1, 0.25, 0.5, 1.0, 2.0, 3.0)
    reps = 1 if quick else 4
    rows = []
    for sigma in sigmas:
        D = sigma ** 2 / 2   # dt = 1 ps
        for rep in range(reps):
            z = ring_walkers(n, frames, D=D, dt=1.0, box_z=box_z, rng=rng)
            result, check = analyse_ring(z, box_z, half)
            row = crevice_counts(result)
            T = frames - 1
            q0_cont = n * D / (L * box_z)
            q0_disc = n * D / ((L + 2 * SIEGMUND_BETA * sigma) * box_z)
            q0 = (row["events_up"] + row["events_down"]) / (2 * T)
            dn_true = n * D / box_z ** 2
            dn = row["water"]["D_n_per_s"] * 1e-12
            lo_lag, hi_lag = row["water"]["collective_diffusion"]["fit_lags_frames"]
            dn_window = expected_window_fit(lo_lag, hi_lag, 1.0, n, D, L, box_z)
            ci = row["water"]["collective_diffusion"].get("interval_per_s")
            rows.append({"model": "independent", "sigma_A": sigma, "sigma_over_L": sigma / L, "replicate": rep,
                         "molecules": n, "frames": frames, "L_A": L, "period_A": box_z, "D_A2_per_ps": D,
                         "events_up": row["events_up"], "events_down": row["events_down"],
                         "unresolved": row["unresolved_crossings"],
                         "q0_per_ps": q0, "q0_continuum_per_ps": q0_cont, "q0_siegmund_per_ps": q0_disc,
                         "q0_over_continuum": q0 / q0_cont, "q0_over_siegmund": q0 / q0_disc,
                         "D_n_per_ps": dn, "D_n_true_per_ps": dn_true, "D_n_over_true": dn / dn_true,
                         "fit_lags_frames": [lo_lag, hi_lag], "D_n_window_expected_per_ps": dn_window,
                         "D_n_over_window_expected": dn / dn_window, "window_expected_over_true": dn_window / dn_true,
                         "D_n_interval_covers_true": bool(ci and ci[0] * 1e-12 <= dn_true <= ci[1] * 1e-12),
                         "D_n_interval_covers_window_expected": bool(ci and ci[0] * 1e-12 <= dn_window <= ci[1] * 1e-12),
                         "pf_over_pd": row["water"]["pf_over_pd"], "pf_over_pd_true": L / box_z,
                         "events_identical_to_reference": check["events_identical_to_reference"],
                         "differing_events": check["differing_events"],
                         "differing_events_at_half_period_ties": check["differing_events_at_half_period_ties"],
                         "n_max_abs_difference_vs_reference": check["n_max_abs_difference"],
                         "crevice_seconds": check["seconds"]})
            print(f"ring sigma={sigma} rep={rep}: q0/siegmund={q0 / q0_disc:.3f} q0/cont={q0 / q0_cont:.3f} "
                  f"Dn/true={dn / dn_true:.3f} Dn/window={dn / dn_window:.3f} ref-identical={check['events_identical_to_reference']} "
                  f"dn-diff={check['n_max_abs_difference']:.2e}", flush=True)
    return rows


def run_single_file(quick, rng):
    box_z, half = 20.0, 5.0
    L = 2 * half
    n, frames = (10, 4000) if quick else (16, 40000)
    sigma = 0.5
    D = sigma ** 2 / 2
    reps = 1 if quick else 3
    rows = []
    for rep in range(reps):
        x = ring_walkers(n, frames, D=D, dt=1.0, box_z=box_z, rng=rng)
        y = single_file_relabel(x, box_z)
        for model, z in (("independent", x), ("single_file", y)):
            for stride in ((1, 10) if quick else (1, 2, 5, 10, 20)):
                result, check = analyse_ring(z, box_z, half, stride=stride, msd_max_lag=max(2, 1000 // stride))
                row = crevice_counts(result)
                T = (len(range(0, frames, stride)) - 1) * stride
                q0 = (row["events_up"] + row["events_down"]) / (2 * T)
                dn = row["water"]["D_n_per_s"] * 1e-12
                dn_true = n * D / box_z ** 2
                lo_lag, hi_lag = row["water"]["collective_diffusion"]["fit_lags_frames"]
                dn_window = expected_window_fit(lo_lag, hi_lag, float(stride), n, D, L, box_z)
                rows.append({"model": model, "replicate": rep, "stride_frames": stride, "frame_interval_ps": stride,
                             "sigma_per_saved_frame_A": sigma * math.sqrt(stride), "molecules": n, "frames": frames,
                             "events_up": row["events_up"], "events_down": row["events_down"],
                             "unresolved": row["unresolved_crossings"], "q0_per_ps": q0,
                             "q0_independent_continuum_per_ps": n * D / (L * box_z),
                             "D_n_per_ps": dn, "D_n_true_per_ps": dn_true, "D_n_over_true": dn / dn_true,
                             "fit_lags_frames": [lo_lag, hi_lag], "D_n_over_window_expected": dn / dn_window,
                             "window_expected_over_true": dn_window / dn_true,
                             "pf_over_pd": row["water"]["pf_over_pd"],
                             "mean_occupancy": row["occupancy"]["mean"],
                             "events_identical_to_reference": check["events_identical_to_reference"],
                             "differing_events": check["differing_events"],
                             "n_max_abs_difference_vs_reference": check["n_max_abs_difference"]})
                print(f"{model} rep={rep} stride={stride}: q0={q0:.3e} Dn/true={dn / dn_true:.3f} "
                      f"pf/pd={row['water']['pf_over_pd']} identical={check['events_identical_to_reference']}", flush=True)
    return rows


# ---------------------------------------------------------------------------
# D. 3D Brownian dynamics through a pore in an impenetrable membrane
# ---------------------------------------------------------------------------

def membrane_bd(n, frames, *, D, dt, box, half_thickness, pore_radius, rng):
    """Rejection Brownian dynamics; positions relative to the box centre (pore axis through it)."""
    box = np.asarray(box, float)

    def blocked(p):
        return (np.abs(p[:, 2]) < half_thickness) & (np.hypot(p[:, 0], p[:, 1]) > pore_radius)

    x = rng.uniform(-box / 2, box / 2, size=(n, 3))
    bad = blocked(x)
    while bad.any():
        x[bad] = rng.uniform(-box / 2, box / 2, size=(bad.sum(), 3))
        bad = blocked(x)
    out = np.empty((frames, n, 3))
    out[0] = x
    sigma = math.sqrt(2 * D * dt)
    for t in range(1, frames):
        trial = x + rng.normal(0, sigma, size=x.shape)
        trial -= box * np.round(trial / box)
        ok = ~blocked(trial)
        x = np.where(ok[:, None], trial, x)
        out[t] = x
    return out


def run_lumen_sensitivity(quick, rng):
    box = (30.0, 30.0, 40.0)
    half_t, pore = 5.0, 4.0
    n, frames = (150, 3000) if quick else (600, 20000)
    D = 0.2
    pos = membrane_bd(n, frames, D=D, dt=1.0, box=box, half_thickness=half_t, pore_radius=pore, rng=rng)
    u = memory_universe(["TIP3"] * n, pos, box=box, dt=1.0, wrap=True)
    c = np.asarray(box) / 2
    definitions = [("membrane_slab_pore_radius", half_t, pore), ("membrane_slab_no_radius", half_t, None),
                   ("membrane_slab_radius_minus_1A", half_t, pore - 1.0),
                   ("membrane_slab_radius_plus_1A", half_t, pore + 1.0),
                   ("membrane_slab_radius_plus_2A", half_t, pore + 2.0),
                   ("extent_plus_2A_pore_radius", half_t + 2.0, pore), ("extent_plus_2A_radius_plus_2A", half_t + 2.0, pore + 2.0),
                   ("extent_plus_4A_pore_radius", half_t + 4.0, pore), ("extent_plus_4A_radius_plus_4A", half_t + 4.0, pore + 4.0)]
    channels = [CylinderChannel(name, (0, 0, 1), c[2] - h, c[2] + h, point=tuple(c), radius=r) for name, h, r in definitions]
    # Ground truth from every 1 ps step: every membrane transit passes through the pore.
    s1, img1, lat1, _ = reference_axial_lateral(u, np.arange(N_PROTEIN, N_PROTEIN + n), tuple(c))
    physical = summarize_events(reference_events(reference_states(img1, lat1, c[2] - half_t, c[2] + half_t, None), s1,
                                                 c[2] - half_t, c[2] + half_t, s_image=img1))
    phys_total = physical["up"] + physical["down"]
    rows, seconds = [], 0.0
    for stride in ((1, 10) if quick else (1, 10, 50)):
        s, img, lat, _ = reference_axial_lateral(u, np.arange(N_PROTEIN, N_PROTEIN + n), tuple(c), stride=stride)
        frame_numbers = list(range(0, frames, stride))
        index = {f: j for j, f in enumerate(frame_numbers)}
        results = {}
        for conv in ("compartment", "strict"):
            t0 = time.perf_counter()
            results[conv] = analyze_permeation(u, species=["TIP3"], channels=channels, msd_blocks=5, align=False,
                                               stride=stride, convention=conv)
            seconds += time.perf_counter() - t0
        for k, (name, h, r) in enumerate(definitions):
            states = reference_states(img, lat, c[2] - h, c[2] + h, r)
            axial = summarize_events(reference_events(reference_states(img, lat, c[2] - h, c[2] + h, None), s,
                                                      c[2] - h, c[2] + h, convention="compartment", s_image=img))
            row = {"definition": name, "saved_interval_ps": stride, "sigma_per_saved_frame_A": math.sqrt(2 * D * stride),
                   "half_extent_A": h, "lateral_radius_A": r, "membrane_half_thickness_A": half_t, "pore_radius_A": pore,
                   "molecules": n, "frames": frames, "physical_transits_1ps": phys_total,
                   "physical_net_up_minus_down": physical["up"] - physical["down"],
                   "axial_only_same_extent_and_interval": axial["up"] + axial["down"]}
            for conv in ("compartment", "strict"):
                cr = crevice_counts(results[conv], channel_index=k)
                ref = reference_events(states, s, c[2] - h, c[2] + h, convention=conv, s_image=img)
                keys = sorted((int(e["atom_index"]) - N_PROTEIN, e["direction"], index[int(e["entry_frame"])],
                               index[int(e["exit_frame"])], e["status"])
                              for e in results[conv].events if e["channel"] == name)
                total = cr["events_up"] + cr["events_down"]
                row.update({f"{conv}_events": total, f"{conv}_over_physical": total / phys_total if phys_total else None,
                            f"{conv}_net_up_minus_down": cr["events_up"] - cr["events_down"],
                            f"{conv}_unresolved": cr["unresolved_crossings"],
                            f"{conv}_events_identical_to_reference": keys == event_keys(ref)})
                if conv == "compartment":
                    row["compartment_outside_lumen"] = cr["transits_outside_lumen"]
                    row["compartment_over_axial_only"] = (total / row["axial_only_same_extent_and_interval"]
                                                          if row["axial_only_same_extent_and_interval"] else None)
                    row["compartment_pf_over_pd"] = cr["water"]["pf_over_pd"]
                    row["compartment_pd_dispersion_index"] = cr["water"]["event_dispersion_index"]
            row["mean_occupancy"] = crevice_counts(results["compartment"], channel_index=k)["occupancy"]["mean"]
            rows.append(row)
            print(f"lumen stride={stride} {name}: compartment={row['compartment_events']} "
                  f"(outside {row['compartment_outside_lumen']}) strict={row['strict_events']} physical={phys_total} "
                  f"identical={row['compartment_events_identical_to_reference']}/{row['strict_events_identical_to_reference']}",
                  flush=True)
    return rows, seconds


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--quick", action="store_true", help="small sizes (smoke test)")
    parser.add_argument("--results", default=str(RESULTS), help="output directory")
    parser.add_argument("--seed", type=int, default=20261025)
    args = parser.parse_args()
    out = __import__("pathlib").Path(args.results)
    rng = np.random.default_rng(args.seed)
    tag = "_quick" if args.quick else ""
    t0 = time.perf_counter()
    cases = run_cases()
    write_csv(out / f"permeation_cases{tag}.csv", cases,
              comment="Scripted permeation paths: expected counts written from the definition vs CREVICE and the "
                      "independent reference counter, for the strict convention [up, down, unresolved] and the "
                      "compartment convention [up, down, outside_lumen]. Software verification only.")
    ring = run_ring(args.quick, rng)
    write_csv(out / f"permeation_ring{tag}.csv", ring,
              comment="Independent Brownian walkers on a periodic line through a channel: exact q0 = N D/(L Lambda), "
                      "Siegmund-corrected q0 = N D/((L + 2 beta sigma) Lambda), D_n = N D/Lambda^2. Synthetic.")
    sf = run_single_file(args.quick, rng)
    write_csv(out / f"permeation_single_file{tag}.csv", sf,
              comment="Single-file (hard-core) vs independent walkers on the same period; D_n exact = N D/Lambda^2 "
                      "for both; frame stride varied on the same trajectory. Synthetic.")
    lumen, seconds = run_lumen_sensitivity(args.quick, rng)
    write_csv(out / f"permeation_lumen_sensitivity{tag}.csv", lumen,
              comment="3D Brownian particles, impenetrable membrane slab with a cylindrical pore: physical transits vs "
                      "counts under channel definitions that match or exceed the pore. Synthetic.")
    summary = {
        "script": "benchmarks/nongeom/permeation_counting.py", "quick": args.quick, "seed": args.seed,
        "machine": machine_record(), "wall_seconds": time.perf_counter() - t0,
        "cases_all_match_expected": all(r["crevice_matches_expected"] for r in cases
                                        if r["case"] != "n_entry_through_opposite_periodic_face"),
        "cases_mismatching": [r["case"] for r in cases if not r["crevice_matches_expected"]],
        "cases_reference_all_match_expected": all(r.get("reference_matches_expected", True) for r in cases),
        "ring_events_identical_to_reference": all(r["events_identical_to_reference"] for r in ring + sf),
        "ring_runs_identical": f'{sum(r["events_identical_to_reference"] for r in ring + sf)}/{len(ring + sf)}',
        "ring_n_max_abs_difference": max(r["n_max_abs_difference_vs_reference"] for r in ring + sf),
        "ring_q0_over_siegmund_by_sigma": {str(s): float(np.mean([r["q0_over_siegmund"] for r in ring if r["sigma_A"] == s]))
                                           for s in sorted({r["sigma_A"] for r in ring})},
        "ring_q0_over_continuum_by_sigma": {str(s): float(np.mean([r["q0_over_continuum"] for r in ring if r["sigma_A"] == s]))
                                            for s in sorted({r["sigma_A"] for r in ring})},
        "ring_D_n_over_true_mean": float(np.mean([r["D_n_over_true"] for r in ring])),
        "ring_D_n_over_window_expected_by_sigma": {
            str(s): [float(np.mean([r["D_n_over_window_expected"] for r in ring if r["sigma_A"] == s])),
                     float(np.std([r["D_n_over_window_expected"] for r in ring if r["sigma_A"] == s], ddof=1))
                     if sum(r["sigma_A"] == s for r in ring) > 1 else None]
            for s in sorted({r["sigma_A"] for r in ring})},
        "ring_window_expected_over_true_by_sigma": {str(s): next(r["window_expected_over_true"] for r in ring
                                                                 if r["sigma_A"] == s)
                                                    for s in sorted({r["sigma_A"] for r in ring})},
        "ring_differing_events_total": sum(r["differing_events"] for r in ring),
        "ring_differing_events_at_half_period_ties": sum(r["differing_events_at_half_period_ties"] for r in ring),
        "ring_D_n_over_true_sd": float(np.std([r["D_n_over_true"] for r in ring], ddof=1)) if len(ring) > 1 else None,
        "lumen_seconds": seconds,
        "lumen": {f'{r["definition"]}@{r["saved_interval_ps"]}ps': {
            "compartment_over_physical": r["compartment_over_physical"],
            "compartment_outside_lumen": r["compartment_outside_lumen"],
            "strict_over_physical": r["strict_over_physical"],
            "identical_to_reference": r["compartment_events_identical_to_reference"]
            and r["strict_events_identical_to_reference"]} for r in lumen},
    }
    dump(out / f"permeation_counting_summary{tag}.json", summary)
    print(f"done in {summary['wall_seconds']:.0f} s; cases match: {summary['cases_all_match_expected']}")


if __name__ == "__main__":
    main()

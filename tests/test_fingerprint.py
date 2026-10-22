"""Software-level tests of crevice.fingerprint on synthetic geometry (no real-system validation)."""
import csv
import json
import math

import pytest

from crevice.channels import profile_along_path
from crevice.fingerprint import (ELECTROLYTES, conductance_estimate, pore_fingerprint, resistance_integral,
                                 residue_classes, residue_hydrophobicity, write_fingerprint_bundle)
from crevice.models import Atom, ChannelPoint, PoreProfile, StructureFrame


def _ring_frame(resname_for_z, *, wall=4.5, zs=range(-10, 11, 2), per_ring=10):
    """Rings of carbon atoms around the z axis; one residue per ring, named by z."""
    atoms, serial = [], 1
    for resid, z in enumerate(zs, 1):
        for k in range(per_ring):
            a = 2*math.pi*k/per_ring
            atoms.append(Atom(serial, f"C{k}", resname_for_z(z), "A", resid, wall*math.cos(a), wall*math.sin(a),
                              float(z), "C"))
            serial += 1
    return StructureFrame(tuple(atoms))


def _cylinder_profile(radius, length, n=41):
    pts = tuple(ChannelPoint(i, (0.0, 0.0, -length/2 + length*i/(n-1)), radius, radius, -length/2 + length*i/(n-1))
                for i in range(n))
    return PoreProfile((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), pts, 0.0, method="axial-refined")


def test_cylinder_conductance_matches_analytic_integral():
    a, length = 2.0, 30.0
    rho = ELECTROLYTES["kcl-150mm"]
    s = [length*i/40 for i in range(41)]
    g = conductance_estimate(s, [a]*41, resistivity_ohm_m=rho)
    analytic_factor = length/(math.pi*a*a)                       # 1/Å
    assert g["geometric_factor_per_A"] == pytest.approx(analytic_factor, rel=1e-12)
    expected_pS = 1e12/(rho*analytic_factor*1e10)
    assert g["g_macro_pS"] == pytest.approx(expected_pS, rel=1e-12)
    access = rho*1e10*(2/(4*a))                                  # Hall access resistance, both mouths
    assert g["g_access_pS"] == pytest.approx(1e12/(rho*analytic_factor*1e10 + access), rel=1e-12)
    assert g["g_hole_corrected_pS"] == pytest.approx(expected_pS/5.594)
    # HOLE's own convention: Gmacro = 1200/F pS for rho = 1/12 Ohm m
    assert conductance_estimate(s, [a]*41, resistivity_ohm_m=1/12)["g_macro_pS"] == pytest.approx(1200/analytic_factor)


def test_conductance_occluded_and_invalid():
    g = conductance_estimate([0.0, 1.0, 2.0], [1.0, 0.0, 1.0])
    assert g["occluded"] and g["g_macro_pS"] == 0.0 and math.isinf(resistance_integral([0, 1], [0, 1]))
    with pytest.raises(ValueError):
        conductance_estimate([0.0, 1.0], [1.0, 1.0], resistivity_ohm_m=0)
    with pytest.raises(ValueError):
        conductance_estimate([0.0], [1.0])


def test_cylinder_profile_fingerprint_conductance_uses_lumen():
    frame = _ring_frame(lambda z: "LEU")
    fp = pore_fingerprint(frame, _cylinder_profile(2.0, 20.0))
    assert [r["segment"] for r in fp.samples] == ["lumen"]*41
    assert fp.samples[0]["s_A"] == 0.0 and fp.samples[-1]["s_A"] == pytest.approx(20.0)
    assert fp.conductance["geometric_factor_per_A"] == pytest.approx(20.0/(math.pi*4.0))


def test_hydrophobic_polar_wall_sign_pattern():
    frame = _ring_frame(lambda z: "LEU" if z < 0 else "ASP")
    path = [(0.0, 0.0, z/2) for z in range(-20, 21)]
    profile = profile_along_path(frame, path)
    fp = pore_fingerprint(frame, profile, hydrophobicity_bandwidth=2.0)
    h = {round(r["s_A"]): r["hydrophobicity_kcal_mol"] for r in fp.samples}
    assert h[2] > 0 and h[18] < 0                                 # s = z + 10
    assert h[2] == pytest.approx(1.25, abs=0.05) and h[18] == pytest.approx(-3.64, abs=0.05)
    lumen = fp.segments[0]
    assert lumen["n_hydrophobic"] == 5 and lumen["n_charged"] == 6
    assert {e["parent"] for e in fp.lining} == {"LEU", "ASP"}


def test_segment_boundaries_from_lateral_leg_and_labels():
    profile = _cylinder_profile(2.0, 10.0, n=11)
    leg1 = {"leg": 1, "points": [[0.0, 0.0, -5.0], [1.0, 0.0, -5.5], [2.0, 0.0, -6.0]],
            "raw_clearances_A": [2.0, 2.2, 2.4], "arc_lengths_A": [0.0, 1.118, 2.236], "length_A": 2.236,
            "bottleneck_radius_A": 2.0}
    leg2 = {**leg1, "leg": 2, "length_A": 3.0}
    profile = PoreProfile(profile.axis_origin, profile.axis_direction, profile.points, 0.0, method="axial-refined",
                          metadata={"exits": {"lower": {"type": "lateral", "legs": [leg2, leg1]},
                                              "upper": {"type": "axial"}}})
    frame = _ring_frame(lambda z: "SER")
    fp = pore_fingerprint(frame, profile)
    exit_rows = [r for r in fp.samples if r["segment"] == "exit_1"]
    assert [round(r["s_A"], 3) for r in exit_rows] == [-2.236, -1.118]   # mouth sample (arc 0) is the lumen end
    assert all(r["source"] == "lateral_exit_leg" for r in exit_rows)
    assert [s["segment"] for s in fp.segments] == ["lumen", "exit_1"]
    lumen = fp.segments[0]
    assert lumen["s_start_A"] == 0.0 and lumen["s_end_A"] == pytest.approx(10.0)
    assert [(leg["segment"], leg["on_path"]) for leg in fp.exit_legs] == [("exit_1", True), ("exit_2", False)]
    flipped = pore_fingerprint(frame, profile, entry_end="start")
    assert {r["segment"] for r in flipped.samples} == {"lumen", "entry_1"}


def test_residue_scale_and_classes():
    assert residue_hydrophobicity("TRP") == ("TRP", 2.09)
    assert residue_hydrophobicity("LYS") == ("LYS", -2.80)
    assert residue_hydrophobicity("ASH") == ("ASP", -0.43)
    assert residue_hydrophobicity("DLE") == ("LEU", 1.25)
    assert residue_hydrophobicity("TIP3") == (None, None)
    assert residue_classes("ARG") == ("charged",) and residue_classes("CLA") == ("other",)


def test_bundle_files(tmp_path):
    frame = _ring_frame(lambda z: "LEU" if z < 0 else "SER")
    fp = pore_fingerprint(frame, profile_along_path(frame, [(0.0, 0.0, z/2) for z in range(-20, 21)]))
    files = write_fingerprint_bundle(fp, tmp_path, "syn", dpi=60)
    rows = list(csv.DictReader(open(files["fingerprint_csv"], encoding="utf-8")))
    assert len(rows) == 41 and {"s_A", "radius_A", "hydrophobicity_kcal_mol", "segment"} <= set(rows[0])
    summary = json.loads(open(files["fingerprint_json"], encoding="utf-8").read())
    assert summary["conductance"]["electrolyte"] == "kcl-150mm" and "not a measurement" in summary["conductance"]["caveat"]
    assert open(files["fingerprint_png"], "rb").read(8) == b"\x89PNG\r\n\x1a\n"


def test_cli_parser_has_fingerprint():
    from crevice.cli import build_parser
    args = build_parser().parse_args(["profile", "--fingerprint", "x.pdb", "--out-dir", "o", "--electrolyte", "kcl-1m"])
    assert args.func.__name__ == "_cmd_fingerprint" and args.electrolyte == "kcl-1m"
    assert build_parser().parse_args(["report", "x.pdb", "--out-dir", "o", "--fingerprint"]).fingerprint


def test_count_path_waters_assigns_nearest_sample_inside_free_sphere():
    from crevice.fingerprint import count_path_waters
    centres = [(0.0, 0.0, z) for z in (0.0, 1.0, 2.0)]
    oxygens = [(0.0, 0.0, 0.1), (0.4, 0.0, 1.0), (1.2, 0.0, 1.0), (0.0, 0.0, 2.4), (0.0, 0.0, 9.0)]
    counts = count_path_waters(centres, [0.5, 1.0, -0.2], oxygens)
    assert counts.tolist() == [1.0, 1.0, 0.0]   # (1.2,0,1) outside r=1; sample 2 occluded (negative clearance)


def test_nearest_lining_residue_ties_do_not_depend_on_hash_order():
    """Samples midway between two rings are equally close to both; the choice must be reproducible."""
    import os
    import subprocess
    import sys
    code = ("import sys; sys.path.insert(0, 'tests'); import test_fingerprint as t; "
            "fp = t.pore_fingerprint(t._ring_frame(lambda z: 'LEU'), t._cylinder_profile(2.0, 20.0)); "
            "print([r['nearest_lining_residue'] for r in fp.samples])")
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    outputs = set()
    for seed in ("1", "2", "3", "4"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        outputs.add(subprocess.run([sys.executable, "-c", code], cwd=root, env=env, capture_output=True,
                                   text=True, encoding="utf-8", check=True).stdout)
    assert len(outputs) == 1

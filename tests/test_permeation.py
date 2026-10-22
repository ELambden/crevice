"""Permeation events, pf/pd and ion current on synthetic trajectories with scripted crossings.

These are software tests: every trajectory is generated here, so the expected
counts, directions, frames and transit times are known exactly. They do not
validate any real channel.
"""
import json
import math

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")
from MDAnalysis.coordinates.memory import MemoryReader  # noqa: E402

from crevice.permeation import (CylinderChannel, ProfileChannel, RegionChannel, analyze_permeation,  # noqa: E402
                                collective_diffusion, parse_charge_overrides, parse_permeation_tokens,
                                poisson_interval, read_applied_field, voltage_record, write_permeation_outputs,
                                WATER_MOLECULE_VOLUME_CM3)

BOX = (30.0, 30.0, 40.0)
RING = [(6 * math.cos(a), 6 * math.sin(a), 0.0) for a in np.linspace(0, 2 * math.pi, 6, endpoint=False)]


def build(molecules, frames, *, dt=10.0, box=BOX, transform=None, wrap=True):
    """Universe: a fixed 6-atom CA ring (protein) plus one residue per scripted molecule.

    ``molecules`` is a list of ``(resname, atom names, path)`` where ``path(t)``
    gives the representative position at frame ``t`` (other atoms are offset).
    """
    names, resnames, resindex, coords = [], [], [], []
    for k in range(len(RING)):
        names.append("CA"); resnames.append("ALA"); resindex.append(k)
    res = len(RING)
    for resname, atom_names, _ in molecules:
        for name in atom_names:
            names.append(name); resnames.append(resname); resindex.append(res)
        res += 1
    n_atoms = len(names)
    u = mda.Universe.empty(n_atoms, n_residues=res, atom_resindex=resindex, trajectory=False)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", [resnames[resindex.index(r)] for r in range(res)])
    u.add_TopologyAttr("resids", list(range(1, res + 1)))
    u.add_TopologyAttr("segids", ["SYS"])
    xyz = np.zeros((frames, n_atoms, 3))
    centre = np.asarray(box) / 2
    for t in range(frames):
        rows = [np.asarray(p) + centre for p in RING]
        for _, atom_names, path in molecules:
            p = np.asarray(path(t), float) + centre
            rows += [p + [0.3 * j, 0, 0] for j in range(len(atom_names))]
        frame = np.asarray(rows)
        if transform is not None:
            frame = transform(t, frame)
        if wrap:
            frame = np.mod(frame, box)
        xyz[t] = frame
    u.load_new(xyz.astype(np.float32), format=MemoryReader, dimensions=np.array([*box, 90, 90, 90], np.float32), dt=dt)
    return u


def path_from_z(zs, x=0.0, y=0.0):
    return lambda t: (x, y, zs[min(t, len(zs) - 1)])


def channel(**kw):
    return CylinderChannel("channel", (0, 0, 1), kw.pop("zmin", -5.0) + 20, kw.pop("zmax", 5.0) + 20,
                           point=(15.0, 15.0, 20.0), radius=kw.pop("radius", 3.0), **kw)


WATER = ("OH2", "H1", "H2")


def run(u, species=("TIP3",), channels=None, **kw):
    return analyze_permeation(u, species=species, channels=channels or [channel()], msd_blocks=2, **kw)


def test_scripted_crossings_counts_directions_frames_and_transit_times():
    up = [-10, -8, -4, 0, 4, 8, 10, 10, 10, 10]           # enters at frame 2, above at 5
    down = [12, 12, 6, 3, -1, -3, -6, -9, -9, -9]         # enters 3, below at 6
    same_side = [-10, -7, -3, 1, -2, -6, -9, -9, -9, -9]   # enters and leaves below: no event
    started_inside = [0, 2, 4, 7, 9, 9, 9, 9, 9, 9]       # starts inside: not an event
    u = build([("TIP3", WATER, path_from_z(up)), ("TIP3", WATER, path_from_z(down)),
               ("TIP3", WATER, path_from_z(same_side)), ("TIP3", WATER, path_from_z(started_inside))], 10)
    result = run(u)
    row = result.summary["channels"][0]["species"][0]
    assert (row["events_up"], row["events_down"], row["net_up_minus_down"]) == (1, 1, 0)
    events = sorted(result.events, key=lambda e: e["entry_frame"])
    assert [(e["direction"], e["entry_frame"], e["exit_frame"], e["transit_time_ps"]) for e in events] == [
        ("up", 2, 5, 30.0), ("down", 3, 6, 30.0)]
    assert [e["resid"] for e in events] == [7, 8]
    exits = {(v["resid"], v["entered_from"], v["exited_to"], v["censored"]) for v in result.visits}
    assert (9, "below", "below", False) in exits and (10, "unknown", "above", True) in exits
    occupancy = result.frames["series"]["channel/TIP3"]["occupancy"]
    assert occupancy[:6] == [1, 1, 3, 3, 3, 1]
    assert row["water"]["q0_per_ns"] == pytest.approx(1 / 90 * 1e3)   # (1 + 1) / 2 per 90 ps
    assert row["water"]["pd_cm3_per_s"] == pytest.approx(WATER_MOLECULE_VOLUME_CM3 * 1 / 90e-12)


def test_periodic_wrapping_counts_two_complete_crossings_and_not_the_bulk_return():
    # Up through the channel, on up through the periodic boundary (box 40 Å) and up through it again.
    zs = list(np.arange(-12, 61, 2.0))
    u = build([("TIP3", WATER, path_from_z(zs))], len(zs))
    result = run(u)
    row = result.summary["channels"][0]["species"][0]
    assert (row["events_up"], row["events_down"], row["unresolved_crossings"]) == (2, 0, 0)
    # A molecule that only travels through the bulk (outside the lumen radius) never permeates.
    u = build([("TIP3", WATER, path_from_z(zs, x=10.0))], len(zs))
    row = run(u).summary["channels"][0]["species"][0]
    assert row["events_total"] == 0 and row["lumen_visits"] == 0


def test_bulk_return_through_boundary_then_entry_from_the_far_side_is_not_an_event():
    # Starts above, goes up through the boundary to below (bulk only), enters the channel and returns below.
    zs = [8, 12, 16, 19, 22, 26, 30, 34, -8, -4, -2, -6, -10]   # 26..34 + wrap: 34 Å -> -6 Å image
    u = build([("TIP3", WATER, lambda t: (10.0 if 1 <= t <= 7 else 0.0, 0.0, zs[t]))], len(zs))
    row = run(u).summary["channels"][0]["species"][0]
    assert row["events_total"] == 0 and row["unresolved_crossings"] == 0


def test_leaving_laterally_clears_the_side_and_whole_extent_jumps_are_unresolved():
    lateral = [(0, -10), (0, -7), (0, -3), (8, -2), (0, 2), (0, 7), (0, 10)]   # inside -> out -> inside -> above
    jump = [-8, -8, 8, 8, 8, 8, 8]                                             # crosses 16 Å in one frame
    u = build([("TIP3", WATER, lambda t: (lateral[t][0], 0, lateral[t][1])),
               ("TIP3", WATER, path_from_z(jump))], 7)
    result = run(u)
    row = result.summary["channels"][0]["species"][0]
    assert row["events_total"] == 0
    assert row["unresolved_crossings"] == 1
    assert [e["status"] for e in result.events] == ["unresolved_single_interval"]
    assert any("between two saved frames" in w for w in result.summary["warnings"])


def test_alignment_removes_rigid_motion_of_the_whole_system():
    zs = [-10, -8, -4, 0, 4, 8, 10, 10]

    def moved(t, frame):
        angle = 0.2 * t
        rot = np.array([[1, 0, 0], [0, math.cos(angle), -math.sin(angle)], [0, math.sin(angle), math.cos(angle)]])
        c = np.asarray(BOX) / 2
        return (frame - c) @ rot.T + c + [1.5 * t, -t, 0.5 * t]
    reference = run(build([("TIP3", WATER, path_from_z(zs))], 8))
    rotated = run(build([("TIP3", WATER, path_from_z(zs))], 8, transform=moved, wrap=False))
    assert rotated.events[0]["direction"] == "up"
    assert [(e["entry_frame"], e["exit_frame"]) for e in rotated.events] == \
           [(e["entry_frame"], e["exit_frame"]) for e in reference.events]
    assert rotated.summary["reference_frame"]["fit_rmsd_A"]["max"] < 1e-3


def test_multiple_pores_are_counted_separately():
    zs = [-10, -6, -2, 2, 6, 10, 10]
    pores = [CylinderChannel(f"pore{i}", (0, 0, 1), 15, 25, point=(15 + x, 15, 20), radius=2.5)
             for i, x in ((1, -8), (2, 8))]
    u = build([("TIP3", WATER, path_from_z(zs, x=-8)), ("TIP3", WATER, path_from_z(zs[::-1], x=-8)),
               ("TIP3", WATER, path_from_z(zs, x=8))], 7)
    result = analyze_permeation(u, species=["TIP3"], channels=pores, msd_blocks=2)
    by = {c["definition"]["name"]: c["species"][0] for c in result.summary["channels"]}
    assert (by["pore1"]["events_up"], by["pore1"]["events_down"]) == (1, 1)
    assert (by["pore2"]["events_up"], by["pore2"]["events_down"]) == (1, 0)
    assert result.summary["totals_over_channels"][0]["events_up"] == 2


def test_martini_water_mapping_ion_charges_from_atom_names_and_conductance():
    up = [-10, -6, -2, 2, 6, 10, 10, 10]
    u = build([("W", ("W",), path_from_z(up)), ("ION", ("NA+",), path_from_z(up)),
               ("ION", ("CL-",), path_from_z(up[::-1])), ("BB_FAKE", ("BB",), lambda t: (12, 12, 0))], 8, dt=20.0)
    result = analyze_permeation(u, species=["W,ION"], channels=[channel()], voltage_mV=100, msd_blocks=2,
                                alignment_selection="name CA")
    s = result.summary
    assert s["coarse_grained"] is True
    labels = {sp["label"]: sp for sp in s["species"]}
    assert labels["ION:NA+"]["charge_e"] == 1 and labels["ION:CL-"]["charge_e"] == -1
    assert labels["W"]["water_mapping"] == 4 and labels["W"]["kind"] == "water"
    rows = {r["species"]: r for r in s["channels"][0]["species"]}
    duration_s = 7 * 20e-12
    assert rows["W"]["water"]["pd_cm3_per_s"] == pytest.approx(WATER_MOLECULE_VOLUME_CM3 * 4 * 0.5 / duration_s)
    current = s["channels"][0]["ion_current_total"]
    # +1 up and -1 down are both +1 e along +axis.
    expected_pA = 2 * 1.602176634e-19 / duration_s * 1e12
    assert current["current_pA"] == pytest.approx(expected_pA)
    assert current["conductance_pS"] == pytest.approx(expected_pA / 0.1)
    assert s["conductance_defined"] is True


def test_conductance_is_undefined_without_voltage_and_zero_events_give_a_poisson_bound():
    u = build([("SOD", ("SOD",), lambda t: (0, 0, -12)), ("CLA", ("CLA",), lambda t: (0, 0, 12))], 12)
    result = run(u, species=["SOD,CLA"])
    s = result.summary
    current = s["channels"][0]["ion_current_total"]
    assert current["current_pA"] == 0 and current["conductance_pS"] is None
    assert current["conductance_status"].startswith("undefined")
    assert s["conductance_defined"] is False
    low, high = current["poisson"]["interval_pA"]
    assert low < 0 < high
    assert current["block_bootstrap"]["status"] == "no_variation_between_blocks"


def test_charge_overrides_and_absent_tokens():
    u = build([("SOD", ("SOD",), lambda t: (0, 0, -12))], 4)
    result = run(u, species=["SOD,K,POT"], charges=parse_charge_overrides(["SOD=0.75"]))
    sp = result.summary["species"][0]
    assert (sp["charge_e"], sp["charge_source"]) == (0.75, "override")
    assert result.summary["absent_tokens"] == ["K", "POT"]
    with pytest.raises(ValueError, match="no --permeation species"):
        run(u, species=["XYZ"])
    with pytest.raises(ValueError, match="NAME=q"):
        parse_charge_overrides(["SOD"])
    assert parse_permeation_tokens(["W,NA", "w2=resname TIP3 and name OH2", "NA"]) == [
        "W", "NA", "w2=resname TIP3 and name OH2"]


def test_selection_species_and_profile_and_region_channels(tmp_path):
    from crevice.models import ChannelPoint, PoreProfile
    zs = [-10, -6, -2, 2, 6, 10, 10]
    u = build([("TIP3", WATER, path_from_z(zs))], 7)
    points = tuple(ChannelPoint(i, (15.0, 15.0, 20.0 + z), 2.0, 2.0, z) for i, z in enumerate(range(-5, 6)))
    profile = PoreProfile((15.0, 15.0, 20.0), (0.0, 0.0, 1.0), points)
    prof = ProfileChannel("channel", [profile, None, profile], [0, 3, 5], tolerance=0.5)
    result = analyze_permeation(u, species=["wat=resname TIP3 and name OH2"], channels=[prof], msd_blocks=2)
    assert result.summary["channels"][0]["species"][0]["events_up"] == 1
    assert result.summary["species"][0]["label"] == "wat"
    # Region map: a 4 x 4 x 10 Å box of voxels around the axis.
    dx = tmp_path / "region.dx"
    nx, ny, nz = 5, 5, 11
    values = np.zeros((nx, ny, nz))
    values[1:4, 1:4, :] = 1
    header = (f"object 1 class gridpositions counts {nx} {ny} {nz}\norigin 13 13 15\n"
              "delta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\n"
              f"object 2 class gridconnections counts {nx} {ny} {nz}\n"
              f"object 3 class array type double rank 0 items {values.size} data follows\n")
    dx.write_text(header + "\n".join(f"{v:g}" for v in values.ravel()) + "\n")
    region = RegionChannel("region", dx, (0, 0, 1))
    result = analyze_permeation(u, species=["TIP3"], channels=[region], msd_blocks=2)
    row = result.summary["channels"][0]["species"][0]
    assert row["events_up"] == 1
    assert result.summary["channels"][0]["definition"]["extent_A"] == [14.5, 25.5]


def test_collective_diffusion_recovers_a_known_random_walk():
    rng = np.random.default_rng(7)
    d, dt = 0.02, 2.0
    n = np.concatenate([[0], np.cumsum(rng.normal(0, math.sqrt(2 * d * dt), 50000))])
    record = collective_diffusion(n, dt, blocks=5, max_lag=20)
    assert record["D_per_ps"] == pytest.approx(d, rel=0.08)
    low, high = record["interval_per_s"]
    assert low < d * 1e12 < high


def test_collective_coordinate_of_a_rigid_single_file_column():
    """A rigid column of molecules (spacing a) doing one random walk X(t): n = X / a + bounded, so D_n = D / a^2.

    This checks the construction of n(t) (clipped displacements inside the
    extent divided by its length, unwrapped across the periodic boundary) and
    pf = v_w D_n against a value known in closed form.
    """
    rng = np.random.default_rng(3)
    frames, d, dt, a = 5000, 0.05, 1.0, 2.5
    x = np.concatenate([[0], np.cumsum(rng.normal(0, math.sqrt(2 * d * dt), frames - 1))])
    molecules = [("TIP3", ("OH2",), (lambda j: (lambda t: (0.0, 0.0, -20 + 0.3 + j * a + x[t])))(j)) for j in range(16)]
    u = build(molecules, frames, dt=dt)
    result = analyze_permeation(u, species=["TIP3"], channels=[channel()], msd_blocks=5, msd_max_lag=50)
    row = result.summary["channels"][0]["species"][0]
    expected = d / a ** 2 * 1e12
    assert row["water"]["D_n_per_s"] == pytest.approx(expected, rel=0.15)
    assert row["water"]["pf_cm3_per_s"] == pytest.approx(WATER_MOLECULE_VOLUME_CM3 * row["water"]["D_n_per_s"])
    n = np.asarray(result.frames["series"]["channel/TIP3"]["n"])
    assert np.max(np.abs(n - x / a)) < 1.0                       # n tracks X / a up to a bounded term
    assert row["occupancy"]["mean"] == pytest.approx(4.0, abs=0.6)


def test_poisson_interval_and_voltage_records(tmp_path):
    low, high = poisson_interval(4)
    assert low == pytest.approx(1.0899, abs=1e-3) and high == pytest.approx(10.2416, abs=1e-3)
    rec = voltage_record(efield=[0, 0, 0.05], box_lengths_A=[100, 100, 120], axis=(0, 0, -1))
    assert rec["status"] == "from_field" and rec["voltage_mV"] == pytest.approx(-600.0)
    assert voltage_record()["status"] == "undefined"
    assert voltage_record(efield=[0, 0, 0.05], box_lengths_A=[100, 100, 120], axis=(1, 0, 0))["status"] == "undefined"
    mdp = tmp_path / "run.mdp"
    mdp.write_text("; comment\nelectric-field-x = 0 0 0 0\nelectric-field-z = 0.05 0 0 0 ; static\n")
    field = read_applied_field(mdp)
    assert field["E_V_per_nm"] == [0, 0, 0.05] and field["static"] is True
    mdp.write_text("E-z = 1 0.065 0\n")
    assert read_applied_field(mdp)["E_V_per_nm"] == [0, 0, 0.065]
    mdp.write_text("electric-field-z = 0.05 150 0 0\n")
    assert read_applied_field(mdp)["static"] is False


def test_outputs_are_written_and_the_figure_is_clean_by_default(tmp_path):
    from crevice.presentation import figure_description
    zs = [-10, -8, -4, 0, 4, 8, 10, 10, 10, 10]
    u = build([("TIP3", WATER, path_from_z(zs)), ("SOD", ("SOD",), path_from_z(zs[::-1]))], 10)
    result = run(u, species=["TIP3,SOD"])
    files = write_permeation_outputs(result, tmp_path, "sys")
    summary = json.loads((tmp_path / "sys_permeation.json").read_text())
    assert summary["files"]["events_csv"] == "sys_permeation_events.csv"
    assert summary["units"]["permeability"] == "cm^3/s"
    lines = (tmp_path / "sys_permeation_events.csv").read_text().splitlines()
    assert len(lines) == 3 and lines[0].startswith("channel,species,resid")
    frames = (tmp_path / "sys_permeation_frames.csv").read_text().splitlines()
    assert frames[0] == ("frame,time_ps,channel/TIP3:occupancy,channel/TIP3:net_crossings,channel/TIP3:n,"
                         "channel/SOD:occupancy,channel/SOD:net_crossings,channel/SOD:n")
    assert len(frames) == 11
    meta = figure_description(files["figure_png"])
    assert meta["CREVICE annotations"].startswith("omitted")
    assert "Permeation events" in meta["Description"]


def _write_system(tmp_path, molecules, frames):
    u = build(molecules, frames)
    top = tmp_path / "system.pdb"
    traj = tmp_path / "system.xtc"
    u.atoms.write(str(top))
    with mda.Writer(str(traj), u.atoms.n_atoms) as writer:
        for _ in u.trajectory:
            writer.write(u.atoms)
    return top, traj


def test_cli_permeation_only_with_explicit_extent(tmp_path, capsys):
    from crevice.cli import main
    zs = [-10, -8, -4, 0, 4, 8, 10, 10, 10, 10]
    top, traj = _write_system(tmp_path, [("TIP3", WATER, path_from_z(zs)), ("SOD", ("SOD",), path_from_z(zs[::-1]))], 10)
    out = tmp_path / "run" / "aq.json"
    code = main(["trajectory", str(traj), "--topology", str(top), "--permeation", "TIP3,SOD,CLA", "--permeation-only",
                 "--axis", "z", "--channel-extent", "15,25", "--channel-radius", "3", "--origin", "15,15,20",
                 "--voltage", "100", "-o", str(out)])
    assert code == 0
    record = json.loads(out.read_text())
    assert record["mode"] == "permeation-only"
    summary = json.loads((tmp_path / "run" / "aq_permeation.json").read_text())
    rows = {r["species"]: r for r in summary["channels"][0]["species"]}
    assert (rows["TIP3"]["events_up"], rows["SOD"]["events_down"]) == (1, 1)
    assert summary["absent_tokens"] == ["CLA"]
    assert summary["channels"][0]["ion_current_total"]["conductance_pS"] < 0
    assert not (tmp_path / "run" / "aq_permeation.png").exists()   # trajectory writes no figures by default
    printed = capsys.readouterr().out
    assert "permeation channel TIP3: up=1 down=0" in printed


def test_cli_permeation_option_errors(tmp_path):
    from crevice.cli import main
    zs = [-10, 0, 10]
    top, traj = _write_system(tmp_path, [("TIP3", WATER, path_from_z(zs))], 3)
    base = ["trajectory", str(traj), "--topology", str(top), "-o", str(tmp_path / "x.json")]
    assert main(base + ["--permeation-only"]) != 0
    assert main(base + ["--permeation", "TIP3", "--permeation-only"]) != 0           # no extent
    assert main(base + ["--channel-extent", "1,2"]) != 0                              # extent without species
    assert main(base + ["--permeation", "TIP3", "--permeation-only", "--channel-extent", "15,25"]) != 0  # no axis


def test_cli_profile_defined_channel(tmp_path):
    """The default channel is the per-frame profile measured by the trajectory command."""
    from crevice.cli import main
    saved = RING[:]
    RING[:] = [(6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6), -10.5 + 1.5 * k)
               for k in range(15) for j in range(12)]
    try:
        zs = list(np.arange(-18, 19, 3.0))
        top, traj = _write_system(tmp_path, [("TIP3", WATER, path_from_z(zs))], len(zs))
    finally:
        RING[:] = saved
    out = tmp_path / "prof.json"
    code = main(["trajectory", str(traj), "--topology", str(top), "--selection", "resname ALA", "--axis", "z",
                 "--permeation", "TIP3", "--outputs", "none", "-o", str(out)])
    assert code == 0
    data = json.loads(out.read_text())
    assert data["metadata"]["permeation"]["channels"][0]["kind"] == "profile"
    summary = json.loads((tmp_path / "prof_permeation.json").read_text())
    row = summary["channels"][0]["species"][0]
    assert (row["events_up"], row["events_down"]) == (1, 0)


@pytest.mark.parametrize("reference", ["origin_without_clearance", "channel_unresolved"])
def test_cli_profile_defined_channel_needs_a_resolved_reference_frame(tmp_path, capsys, reference):
    """Without a reference profile (origin on a wall atom, or no channel along the axis) the run stops clearly.

    The structural outputs are written first (the trajectory command records the
    unresolved reference frame and exits 0 without ``--permeation``); the
    profile-defined permeation channel then has no extent, which is an error.
    """
    from crevice.cli import main
    saved = RING[:]
    RING[:] = [(6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6), -10.5 + 1.5 * k)
               for k in range(15) for j in range(12)]
    try:
        zs = list(np.arange(-18, 19, 6.0))
        top, traj = _write_system(tmp_path, [("TIP3", WATER, path_from_z(zs))], len(zs))
    finally:
        RING[:] = saved
    out = tmp_path / "prof.json"
    # The box centre is (15, 15, 20); a wall atom sits at (21, 15, 20).
    geometry = (["--axis", "z", "--origin", "21,15,20"] if reference == "origin_without_clearance"
                else ["--axis", "x", "--origin", "15,15,20"])
    base = ["trajectory", str(traj), "--topology", str(top), "--selection", "resname ALA", *geometry,
            "--outputs", "none", "--no-align", "-o", str(out)]
    assert main(base) == 0
    data = json.loads(out.read_text())
    assert data["metadata"]["reference_profile_status"] == "unresolved"
    assert ("origin_clearance" in data["metadata"]) == (reference == "origin_without_clearance")
    capsys.readouterr()
    assert main(base + ["--permeation", "TIP3"]) != 0
    assert "crevice: error: no channel profile was resolved, so the channel extent is unknown" in capsys.readouterr().err
    assert not (tmp_path / "prof_permeation.json").exists()


def test_d_n_interval_is_clipped_at_zero_and_json_serialisable():
    rng = np.random.default_rng(5)
    n = np.concatenate([[0], np.cumsum(rng.normal(0, 0.01, 400))]) + rng.normal(0, 1.0, 401)
    record = collective_diffusion(n, 1.0, blocks=5)
    json.dumps(record)
    assert record["interval_per_s"][0] >= 0
    assert record["interval_clipped_at_zero"] is (record["interval_unclipped_per_s"][0] < 0)

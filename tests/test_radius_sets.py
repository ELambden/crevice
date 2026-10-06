"""Atomic radius sets: presets, radius files, the active-set context and provenance.

Synthetic structures only (software tests). A 12-atom x 15-ring carbon
cylinder of radius 6 Å has its axial clearance at 6 Å minus the carbon radius,
so every radius set must move the measured minimum radius by exactly the change
in the carbon radius. With no set chosen, every non-ion atom must get exactly
its standard-table radius, and outputs must record the default set.
"""
from __future__ import annotations

import json
import math
import pickle
import threading

import pytest

from crevice import pore_profile, profile_along_path
from crevice.cli import main
from crevice.models import Atom, StructureFrame
from crevice.radii import (
    DEFAULT_RADII,
    ION_RADII,
    RADIUS_PRESETS,
    VDW_RADII,
    RadiusSet,
    active_radii,
    atom_vdw_radius,
    custom_radii_set,
    infer_element,
    load_radius_file,
    radii_preset,
    radii_provenance,
    read_hole_radius_file,
    resolve_radii,
    use_radii,
    vdw_radius,
)

# HOLE 2.3.1 rad/simple.rad, verbatim.
SIMPLE_RAD = """remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
remark: Simple - Only use one value for each element C O H etc.
remark: van der Waals radii
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR P??? ??? 2.10
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 1.00
VDWR D2?  ASN 1.00
remark: amber lone pairs on sulphurs
VDWR LP?? ??? 0.00
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85
"""


def cylinder(name="C", resname="ALA", element="C"):
    atoms = [Atom(12 * k + j + 1, name, resname, "A", k + 1,
                  6 * math.cos(math.pi * j / 6), 6 * math.sin(math.pi * j / 6), -10.5 + 1.5 * k, element)
             for k in range(15) for j in range(12)]
    return StructureFrame(tuple(atoms))


def atom(name, resname, element=""):
    return Atom(1, name, resname, "A", 1, 0.0, 0.0, 0.0, element)


# --- default behaviour is the standard table ---------------------------------

def test_no_set_means_the_default_set_and_it_is_recorded():
    assert active_radii() is None
    assert radii_provenance()["name"] == "default"
    assert radii_provenance() == DEFAULT_RADII.provenance()
    for name, resname, element in [("CA", "ALA", "C"), ("OG", "SER", ""), ("GD", "GD", "GD"), ("CA", "CA", "C")]:
        assert atom_vdw_radius(atom(name, resname, element)) == vdw_radius(infer_element(name, element, resname))
    profile = pore_profile(cylinder(), axis="z", enclosure_radius=0.8)
    assert profile.min_radius == pytest.approx(6 - 1.70, abs=1e-9)
    assert profile.metadata["radii"]["name"] == "default"
    assert profile.metadata["radii"]["table_sha256"] == DEFAULT_RADII.table_sha256


def test_bondi_preset_equals_the_default_except_for_ion_residues_and_is_recorded():
    frame = cylinder()
    default = pore_profile(frame, axis="z", enclosure_radius=0.8)
    bondi = pore_profile(frame, axis="z", enclosure_radius=0.8, radii="bondi")
    assert [p.radius for p in bondi.points] == [p.radius for p in default.points]
    assert bondi.metadata["radii"]["name"] == "bondi"
    assert bondi.metadata["radii"]["element_overrides_A"] == {}
    assert bondi.metadata["radii"]["ion_radii_A"] == {}
    assert dict(radii_preset("bondi").element_radii) == VDW_RADII
    assert dict(DEFAULT_RADII.element_radii) == VDW_RADII
    assert radii_preset("bondi").radius_for_names("SOD", "SOD") == 2.27
    assert DEFAULT_RADII.radius_for_names("SOD", "SOD") == ION_RADII["NA"] == 1.41075


# --- a non-default set changes the geometry by exactly the radius change ------

@pytest.mark.parametrize("radii, carbon", [("hole", 1.85), ("charmm_like", 2.00),
                                           (custom_radii_set(element_radii={"C": 1.5}), 1.5)])
def test_non_default_set_moves_the_cylinder_radius_by_the_radius_change(radii, carbon):
    frame = cylinder()
    profile = pore_profile(frame, axis="z", enclosure_radius=0.8, radii=radii)
    assert profile.min_radius == pytest.approx(6 - carbon, abs=1e-9)
    straight = profile_along_path(frame, [(0, 0, -3), (0, 0, 3)], radii=radii)
    assert straight.min_radius == pytest.approx(6 - carbon, abs=1e-9)
    assert "radii" in profile.metadata


def test_cli_records_the_radius_set_always(tmp_path):
    lines = [f"ATOM  {a.serial:5d}  C   ALA A{a.resid:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}  1.00 10.00           C"
             for a in cylinder().atoms]
    structure = tmp_path / "cyl.pdb"
    structure.write_text("\n".join(lines) + "\nEND\n")
    assert main(["profile", str(structure), "--axis", "z", "--enclosure-radius", "0.8", "-o", str(tmp_path / "default.json")]) == 0
    assert main(["profile", str(structure), "--axis", "z", "--enclosure-radius", "0.8", "-o", str(tmp_path / "hole.json"), "--radii", "hole"]) == 0
    default = json.loads((tmp_path / "default.json").read_text())
    hole = json.loads((tmp_path / "hole.json").read_text())
    assert default["metadata"]["radii"]["name"] == "default"
    assert hole["metadata"]["radii"]["name"] == "hole"
    assert hole["metadata"]["radii"]["source"] == "HOLE 2.3.1 rad/simple.rad"
    assert hole["min_radius"] == pytest.approx(default["min_radius"] - 0.15, abs=1e-9)
    assert main(["profile", str(structure), "-o", str(tmp_path / "x.json"), "--radii", "no-such-set"]) == 2


# --- the context --------------------------------------------------------------

def test_use_radii_nests_restores_and_none_inherits():
    carbon = atom("CB", "ALA", "C")
    with use_radii("hole") as outer:
        assert outer.name == "hole" and atom_vdw_radius(carbon) == 1.85
        with use_radii(None) as inherited:
            assert inherited is outer
        with use_radii("charmm_like"):
            assert atom_vdw_radius(carbon) == 2.00
        assert atom_vdw_radius(carbon) == 1.85
    assert active_radii() is None and atom_vdw_radius(carbon) == 1.70
    with pytest.raises(RuntimeError):
        with use_radii("hole"):
            raise RuntimeError("boom")
    assert active_radii() is None


def test_context_is_private_to_each_thread():
    seen = {}

    def worker():
        seen["thread"] = active_radii()

    with use_radii("hole"):
        thread = threading.Thread(target=worker)
        thread.start()
        thread.join()
    assert seen["thread"] is None


def test_radius_sets_pickle_for_worker_processes():
    for radii in [*RADIUS_PRESETS.values(), custom_radii_set(base="bondi", atom_radii={"SOD:SOD": 1.02})]:
        copy = pickle.loads(pickle.dumps(radii))
        assert copy == radii and copy.table_sha256 == radii.table_sha256


def test_parallel_cavity_trajectory_workers_use_the_chosen_set(tmp_path):
    from dataclasses import replace

    from crevice.cavity_trajectory import CavityReference, analyze_cavity_trajectory
    from crevice.rolling import rolling_probe_cast
    from crevice.trajectory import Trajectory
    from crevice.volume_export import write_void_cast_dx
    from test_channel_coordinate import rectangular_channel, transform

    frame = rectangular_channel(half_x=4, half_y=5, half_length=8, capped=True)
    cast = rolling_probe_cast(frame, spacing=.75, selection_mode='all', enclosure_fraction=0, min_component_volume=1)
    main_component = max(cast.components, key=lambda c: c.volume)
    cast = replace(cast, components=(main_component,), points=main_component.points)
    ref = CavityReference.from_dx(write_void_cast_dx(cast, tmp_path / 'ref.dx'), axis=(0, 0, 1))
    moved, _, _ = transform(frame)
    traj = Trajectory([replace(frame, frame_index=0), replace(moved, frame_index=1)], time_step=1.)
    options = dict(geometry_mode='reference-region', geometry={'probe_radius': 1.4})
    default = analyze_cavity_trajectory(traj, ref, workers=1, **options)
    serial = analyze_cavity_trajectory(traj, ref, workers=1, radii="charmm_like", **options)
    parallel = analyze_cavity_trajectory(traj, ref, workers=2, radii="charmm_like", **options)
    volumes = lambda result: [row['volume_A3'] for row in result['frames']]
    assert volumes(parallel) == volumes(serial)
    assert all(big < small for big, small in zip(volumes(serial), volumes(default)))


# --- presets -----------------------------------------------------------------

def test_hole_preset_is_holes_simple_rad_matched_by_atom_name():
    hole = radii_preset("hole")
    assert hole.hole_records == read_hole_radius_file_from_text(SIMPLE_RAD)
    assert hole.radius_for_names("CA", "ALA") == 1.85
    assert hole.radius_for_names("OG", "SER") == 1.65
    assert hole.radius_for_names("SD", "MET") == 2.00
    assert hole.radius_for_names("NZ", "LYS") == 1.75
    assert hole.radius_for_names("HA", "ALA") == 1.00
    assert hole.radius_for_names("P", "DOPC") == 2.10
    # Name, not element: a CHARMM sodium ion named SOD matches S??? as in HOLE.
    assert hole.radius_for_names("SOD", "SOD") == 2.00
    with pytest.raises(ValueError, match="no VDWR record"):
        hole.radius_for_names("ZN", "ZN")


def read_hole_radius_file_from_text(text):
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "simple.rad"
        path.write_text(text)
        return read_hole_radius_file(path)


def test_unknown_preset_or_path_is_an_error():
    with pytest.raises(ValueError, match="valid presets"):
        radii_preset("amber")
    with pytest.raises(ValueError, match="Unknown radius set"):
        resolve_radii("no-such-file.json")


# --- radius files ------------------------------------------------------------

def test_hole_rad_file_uses_first_match_fixed_columns_and_is_hashed(tmp_path):
    path = tmp_path / "mixed.rad"
    path.write_text("remark specific records before general\n"
                    "VDWR CB   ALA 2.00\n"
                    "VDWR C??? ??? 1.90\n"
                    "vdwr o??? ??? 1.40\n"
                    "BOND C??? 0.85\n"
                    "VDWR ???? ??? 1.70\n")
    radii = load_radius_file(path)
    assert radii.radius_for_names("CB", "ALA") == 2.00
    assert radii.radius_for_names("CB", "SER") == 1.90
    assert radii.radius_for_names("OG", "SER") == 1.40
    assert radii.radius_for_names("ZN", "ZN") == 1.70
    record = radii.provenance()
    assert record["source"] == str(path.resolve())
    assert len(record["source_sha256"]) == 64
    assert record["hole_records"][0] == ["CB  ", "ALA", 2.0]


@pytest.mark.parametrize("line, message", [("VDWR C??? ??? 2\n", "decimal point"),
                                           ("VDWR C??? ???\n", "decimal point"),
                                           ("VDWR\tC??? ??? 1.8\n", "tab")])
def test_hole_rad_values_that_hole_would_misread_are_refused(tmp_path, line, message):
    path = tmp_path / "bad.rad"
    path.write_text(line)
    with pytest.raises(ValueError, match=message):
        load_radius_file(path)


def test_json_radius_file_with_base_and_overrides(tmp_path):
    path = tmp_path / "ions.json"
    path.write_text(json.dumps({"name": "bondi-ionic", "base": "bondi",
                                "atom_radii": {"SOD:SOD": 1.02, "CLA:CLA": 1.81},
                                "element_radii": {"C": 1.8}, "citation": "test"}))
    radii = load_radius_file(path)
    assert radii.name == "bondi-ionic"
    assert radii.radius_for_names("SOD", "SOD") == 1.02
    assert radii.radius_for_names("OG", "SER") == 1.52          # from the base
    assert radii.radius_for_names("CA", "ALA", "C") == 1.8      # override
    record = radii.provenance()
    assert record["element_overrides_A"] == {"C": 1.8}
    assert record["atom_radii_A"] == {"CLA:CLA": 1.81, "SOD:SOD": 1.02}
    path.write_text(json.dumps({"element_radii": {"Xx": 1.0}}))
    with pytest.raises(ValueError, match="not an element"):
        load_radius_file(path)
    path.write_text(json.dumps({"radii": {}}))
    with pytest.raises(ValueError, match="unknown keys"):
        load_radius_file(path)


def test_csv_radius_file(tmp_path):
    path = tmp_path / "small.csv"
    path.write_text("# a comment\nkind,key,radius\nelement,C,1.9\natom,OG,1.4\ndefault,,2.0\n")
    radii = load_radius_file(path)
    assert radii.radius_for_names("CA", "ALA", "C") == 1.9
    assert radii.radius_for_names("OG", "SER") == 1.4
    assert radii.radius_for_names("NZ", "LYS", "N") == 2.0     # no base: default radius
    assert "N" in radii.provenance()["standard_elements_missing"]
    path.write_text("kind,key,radius\nbase,bondi,\natom,SOD:SOD,1.02\n")
    assert load_radius_file(path).radius_for_names("NZ", "LYS", "N") == 1.55
    path.write_text("kind,key,radius\nwhatever,C,1\n")
    with pytest.raises(ValueError, match="kind must be"):
        load_radius_file(path)


def test_precedence_exact_then_name_then_element_then_default():
    radii = custom_radii_set(element_radii={"O": 1.4}, atom_radii={"OG": 1.5, "SER:OG": 1.6}, default_radius=2.2)
    assert radii.radius_for_names("OG", "SER", "O") == 1.6
    assert radii.radius_for_names("OG", "THR", "O") == 1.5
    assert radii.radius_for_names("OD1", "ASP", "O") == 1.4
    assert radii.radius_for_names("NZ", "LYS", "N") == 2.2
    assert isinstance(radii, RadiusSet)

"""A small sealed synthetic channel for tests that need a resolvable profile.

Five rings of 16 atoms (radius 5 Å, 1.5 Å apart along z), split by quadrant
into four residues: LEU (+x), LYS (+y), ASP (-x) and PHE (-y). The ring is
dense enough that no enclosure probe leaks through the wall, so connected
profiles resolve along z with a clearance of about 3.3-3.5 Å. Software
fixture only.
"""
import math

from crevice.models import Atom, StructureFrame

#: (resid, resname, element) by quadrant, counter-clockwise from +x.
QUADRANTS = ((1, "LEU", "C"), (3, "LYS", "N"), (2, "ASP", "O"), (4, "PHE", "C"))
LAYERS = (-3.0, -1.5, 0.0, 1.5, 3.0)


def channel_rows():
    """``(serial, name, resname, resid, element, x, y, z)`` for every atom."""
    rows, serial = [], 1
    for layer, z in enumerate(LAYERS):
        for k in range(16):
            resid, resname, element = QUADRANTS[((k + 2) // 4) % 4]
            angle = math.pi * k / 8
            rows.append((serial, f"{element}{layer}{k % 4}", resname, resid, element,
                         5.0 * math.cos(angle), 5.0 * math.sin(angle), z))
            serial += 1
    return rows


def channel_frame(source="channel"):
    """The channel as a :class:`~crevice.models.StructureFrame`."""
    return StructureFrame(tuple(Atom(serial, name, resname, "A", resid, x, y, z, element)
                                for serial, name, resname, resid, element, x, y, z in channel_rows()), source=source)


def channel_pdb_text():
    """The channel as PDB text (element column set)."""
    lines = [f"ATOM  {serial:5d} {name:<4s} {resname:3s} A{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
             f"          {element:>2s}"
             for serial, name, resname, resid, element, x, y, z in channel_rows()]
    return "\n".join(lines + ["END"]) + "\n"

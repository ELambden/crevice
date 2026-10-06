"""Cast metadata says exactly what the code did (synthetic shell, software test)."""
import math

import pytest

from crevice.models import Atom, StructureFrame
from crevice.rolling import rolling_probe_cast


def shell():
    golden = math.pi * (3 - math.sqrt(5))
    atoms = []
    for i in range(250):
        y = 1 - (2 * i + 1) / 250
        r = math.sqrt(1 - y * y)
        atoms.append(Atom(i + 1, "C", "ALA", "A", i + 1, 7 * math.cos(golden * i) * r, 7 * y,
                          7 * math.sin(golden * i) * r, "C"))
    return StructureFrame(tuple(atoms))


@pytest.mark.parametrize("mode, fraction, phrase", [
    ("all", 0.0, "not applied"),
    ("all", 0.9, "hard crop of the probe centres"),
    ("dominant", 0.9, "qualifies buried cores only"),
])
def test_enclosure_definition_matches_the_mode(mode, fraction, phrase):
    cast = rolling_probe_cast(shell(), spacing=1.0, selection_mode=mode, enclosure_fraction=fraction)
    meta = cast.metadata
    assert phrase in meta["enclosure_definition"]
    assert meta["enclosure_required_rays"] == math.ceil(26 * fraction)
    assert "probe-inflated" in meta["enclosure_definition"]

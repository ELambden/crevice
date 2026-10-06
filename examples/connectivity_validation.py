"""Reproduce synthetic wall connectivity across grid spacings and offsets."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import time

import crevice
from crevice import Atom, StructureFrame
from crevice.grid import GridConnectivity
from crevice.spatial import SpatialIndex
from crevice.voids import _GridNode, _cluster_indices


def gate_frame(kind):
    atoms = []
    for x, y in itertools.product(range(-6, 7, 2), repeat=2):
        if kind == "wide" and abs(x) <= 2 and abs(y) <= 2:
            continue
        if kind == "narrow" and x == y == 0:
            continue
        serial = len(atoms) + 1
        atoms.append(Atom(serial, "CA", "ALA", "A", serial, x, y, 0, element="C"))
    return StructureFrame(tuple(atoms))


def measure_gate(kind, spacing, phase):
    """Look for a bottom-to-top path inside a fixed lateral region of a wall."""
    frame = gate_frame(kind)
    spatial = SpatialIndex(frame.atoms)
    low, high = (-1, -1, -3), (1, 1, 3)
    origin = tuple(value + phase * spacing for value in low)
    dims = tuple(math.floor((hi - lo) / spacing) + 1 for lo, hi in zip(origin, high))
    if min(dims) < 2 or math.prod(dims) > 100_000:
        raise ValueError("Demo grid must have 2+ points per dimension and at most 100000 points")
    connectivity = GridConnectivity(spatial, origin, spacing, 0.2)
    nodes = {}
    start = time.perf_counter()
    for index in itertools.product(*(range(n) for n in dims)):
        point = connectivity.position(index)
        clearance, nearest = spatial.nearest_surface(point)
        if clearance >= 0.2:
            nodes[index] = _GridNode(point, clearance, nearest, 0)
    components = _cluster_indices(nodes, edge_allowed=lambda a, b:
                                  connectivity.allows(a, b, nodes[a].clearance, nodes[b].clearance))
    through = [c for c in components if any(p[2] == 0 for p in c) and any(p[2] == dims[2] - 1 for p in c)]
    return {"wall": kind, "spacing_A": spacing, "phase_fraction": phase,
            "grid_origin": origin, "grid_dimensions": dims, "grid_points": math.prod(dims),
            "accessible_nodes": len(nodes), "components": len(components),
            "through_components": len(through), "connected": bool(through),
            "continuum_connection_expected": kind != "closed",
            "elapsed_seconds": time.perf_counter() - start,
            "connectivity": connectivity.to_dict()}


def plot_results(rows, spacings, phases, path):
    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from crevice.figures import _pyplot, _save_figure

    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.8), constrained_layout=True)
    for ax, kind in zip(axes, ("closed", "wide", "narrow")):
        values = np.array([[int(next(r["connected"] for r in rows
                           if (r["wall"], r["spacing_A"], r["phase_fraction"]) == (kind, s, p)))
                           for p in phases] for s in spacings])
        ax.imshow(values, vmin=0, vmax=1, cmap=ListedColormap(["#e9edef", "#178c94"]))
        ax.set(title=kind.title() + " wall", xlabel="Grid offset / spacing",
               xticks=range(len(phases)), xticklabels=[f"{p:g}" for p in phases],
               yticks=range(len(spacings)), yticklabels=[f"{s:g}" for s in spacings])
        if ax is axes[0]:
            ax.set_ylabel("Spacing (Angstrom)")
    fig.legend(handles=[Patch(color="#178c94", label="Grid connection found"),
                        Patch(color="#e9edef", label="No grid connection found")],
               loc="outside lower center", ncol=2, frameon=False)
    _save_figure(fig, path, dpi=200)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--spacings", type=float, nargs="+", default=[1, 0.5, 0.25, 0.125])
    parser.add_argument("--phases", type=float, nargs="+", default=[0, 0.25, 0.5, 0.75])
    args = parser.parse_args()
    if any(not math.isfinite(s) or s <= 0 for s in args.spacings):
        parser.error("spacings must be finite and positive")
    if any(not math.isfinite(p) or not 0 <= p < 1 for p in args.phases):
        parser.error("phases must be finite and in [0, 1)")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = [measure_gate(kind, spacing, phase) for kind in ("closed", "wide", "narrow")
            for spacing in args.spacings for phase in args.phases]
    unexpected = [r for r in rows if (r["wall"] == "closed" and r["connected"])
                  or (r["wall"] == "wide" and not r["connected"])
                  or (r["wall"] == "narrow" and r["spacing_A"] <= 0.125 and not r["connected"])]
    source_root = Path(crevice.__file__).parent
    sources = [source_root / name for name in ("geometry.py", "spatial.py", "grid.py", "voids.py", "analysis.py", "models.py")]
    sources.append(Path(__file__))
    report = {"status": "passed" if not unexpected else "failed", "python": platform.python_version(),
              "crevice_version": crevice.__version__, "source": str(source_root),
              "dataset": "synthetic carbon-sphere walls, not proteins or MD",
              "atom_spacing_A": 2, "atom_radius_A": 1.7, "required_clearance_A": 0.2,
              "bounds_A": [[-1, -1, -3], [1, 1, 3]],
              "interpretation": "Missing a grid path is not proof of closure; narrow-wall phase dependence is expected",
              "acceptance": "Closed disconnected; wide connected; narrow connected at <=0.125 A. Coarser narrow misses allowed.",
              "rows": rows, "unexpected_cases": unexpected,
              "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    path = output / "connectivity_sweep.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    plot_results(rows, args.spacings, args.phases, output / "connectivity_sweep.png")
    plot_results(rows, args.spacings, args.phases, output / "connectivity_sweep.svg")
    print(json.dumps({"status": report["status"], "cases": len(rows), "unexpected_cases": len(unexpected),
                     "narrow_connected": sum(r["connected"] for r in rows if r["wall"] == "narrow"),
                     "report": str(path)}, indent=2))
    if unexpected:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

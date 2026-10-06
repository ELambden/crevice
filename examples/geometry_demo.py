"""Reproducible mathematical channel demonstration, with no biological inference.

Run from a source checkout with CREVICE installed:
    python examples/geometry_demo.py --output results/geometry-demo
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
from pathlib import Path

from crevice import (Atom, StructureFrame, analyze_trajectory, profile_distribution,
                   plot_trajectory_profiles, write_profile_distribution_csv,
                   void_cast, write_volume_viewer_bundle,
                   write_void_cast_csv, write_void_cast_json, write_void_cast_pdb, write_trajectory_json,
                   write_trajectory_csv)


def wall_frame(phase: float) -> StructureFrame:
    atoms = []
    for z in range(-8, 9):
        radius = 4.6 + 1.4 * (z / 8) ** 2 + 0.55 * math.sin(phase) * math.exp(-(z / 4) ** 2)
        for j in range(32):
            angle = math.tau * j / 32
            atoms.append(Atom(len(atoms) + 1, "C", "ALA", "A", len(atoms) + 1,
                              radius * math.cos(angle), radius * math.sin(angle), float(z), element="C"))
    return StructureFrame(tuple(atoms), source="synthetic_carbon_wall")


def write_wall_pdb(frame: StructureFrame, path: Path) -> None:
    lines = ["REMARK SYNTHETIC CARBON WALL; NOT A MOLECULAR STRUCTURE"]
    for a in frame.atoms:
        lines.append(f"ATOM  {a.serial:5d} {a.name:>4s} {a.resname:>3s} {a.chain_id}{a.resid:4d}    "
                     f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{1.0:6.2f}{0.0:6.2f}          {a.element:>2s}")
    path.write_text("\n".join(lines) + "\nEND\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/geometry-demo"))
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    frames = [wall_frame(i * math.tau / 24) for i in range(24)]
    analysis = analyze_trajectory(frames, profile_kwargs={"axis": "z", "samples": 65, "search_radius": 0})
    distribution = profile_distribution(analysis)
    write_trajectory_json(analysis, root / "trajectory.json")
    write_trajectory_csv(analysis, root / "trajectory_frames.csv", root / "trajectory_profiles.csv")
    write_profile_distribution_csv(distribution, root / "profile_distribution.csv")
    (root / "profile_distribution.json").write_text(json.dumps(distribution, indent=2) + "\n")
    plot_trajectory_profiles(analysis, root / "profile_distribution.png", title="Synthetic channel: axial radius distribution")
    plot_trajectory_profiles(analysis, root / "profile_distribution.svg", title="Synthetic channel: axial radius distribution")
    structure = root / "wall.pdb"
    write_wall_pdb(frames[0], structure)
    write_wall_pdb(frames[6], root / "wall_expanded.pdb")
    cast = void_cast(frames[0], mode="channel", axis="z", spacing=0.75,
                     profile=analysis.frames[0].profile)
    write_void_cast_json(cast, root / "void_cast.json", include_points=True)
    write_void_cast_pdb(cast, root / "void_cast.pdb")
    write_void_cast_csv(cast, root / "void_cast.csv")
    files = write_volume_viewer_bundle(cast, structure_path=structure, output_dir=root, prefix="synthetic")
    manifest = {
        "dataset": "synthetic carbon-wall channel, not a protein or MD simulation",
        "frames": len(frames), "axis": "z", "profile_samples": 65,
        "profile_search_radius_A": 0, "cast_spacing_A": cast.spacing,
        "cast_volume_A3": cast.total_volume, "cast_components": cast.component_count,
        "python": platform.python_version(),
        "dependencies": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib", "scikit-image")},
        "files": files,
        "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(root.iterdir()) if p.is_file() and p.name != "manifest.json"},
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(root)


if __name__ == "__main__":
    main()


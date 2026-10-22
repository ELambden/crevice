"""Command-line options of the permeation analysis (``--permeation`` and related flags).

The options are added to a trajectory command by :func:`add_permeation_arguments`
and run by :func:`run_permeation` after (or, with ``--permeation-only``,
instead of) the command's structural analysis. The analysis itself is
:func:`crevice.permeation.analyze_permeation`.

Channel definition, in order of precedence:

1. ``--permeation-region FILE``: a region map (OpenDX, in the coordinates of
   the first analysed frame) or a prepared region definition JSON;
2. ``--channel-extent ZMIN,ZMAX`` along ``--axis`` (x, y, z or a vector),
   optionally limited laterally by ``--channel-radius`` around the axis line
   through ``--origin`` (or each ``--permeation-pore`` point, one channel per
   pore);
3. otherwise the per-frame profiles measured by the trajectory command.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _pair(value: str):
    try:
        lo, hi = (float(v) for v in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected ZMIN,ZMAX (Å)") from exc
    if not lo < hi:
        raise argparse.ArgumentTypeError("ZMIN must be smaller than ZMAX")
    return lo, hi


def _vector(value: str):
    try:
        values = tuple(float(v) for v in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected X,Y,Z") from exc
    if len(values) != 3:
        raise argparse.ArgumentTypeError("expected X,Y,Z")
    return values


def _field(value: str):
    try:
        values = [float(v) for v in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected E (along z) or EX,EY,EZ in V/nm") from exc
    if len(values) == 1:
        return [0.0, 0.0, values[0]]
    if len(values) != 3:
        raise argparse.ArgumentTypeError("expected E (along z) or EX,EY,EZ in V/nm")
    return values


def add_permeation_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the ``--permeation`` option group to a trajectory command parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The command's parser; it must already define ``--axis``, ``--origin``,
        ``--topology``, ``--selection``, ``--start``, ``--stop``, ``--align``,
        ``--confidence``, ``--bootstrap-replicates``, ``--bootstrap-seed`` and
        ``--dpi``.
    """
    group = parser.add_argument_group(
        "permeation", "Permeation events, water pf/pd and ion current through the channel (needs --topology)")
    group.add_argument("--permeation", action="append", metavar="SPECIES",
                       help="Species to follow: comma-separated residue or atom names (e.g. TIP3,SOD,CLA or W,ION for "
                            "Martini) or one MDAnalysis selection ('LABEL=resname TIP3 and name OH2'); repeatable. "
                            "Names that match nothing are listed and skipped")
    group.add_argument("--channel-extent", type=_pair, metavar="ZMIN,ZMAX",
                       help="Explicit channel extent (Å) along --axis, as x.axis in the coordinates of the first "
                            "analysed frame; default: each frame's measured profile")
    group.add_argument("--channel-radius", type=float, metavar="Å",
                       help="Lateral lumen radius around the axis line for --channel-extent; default no lateral limit")
    group.add_argument("--permeation-pore", action="append", type=_vector, metavar="X,Y,Z",
                       help="A point on the axis line of one pore (repeat for several pores, e.g. an aquaporin "
                            "tetramer); with --channel-extent. Default: --origin, else the alignment-atom centroid")
    group.add_argument("--permeation-region", metavar="FILE",
                       help="Channel from a region map (.dx, first-frame coordinates) or a prepared region "
                            "definition (.json); needs an axis (--axis, or the definition's)")
    group.add_argument("--lumen-tolerance", type=float, metavar="Å",
                       help="Added to the profile radius (profile channels) or used to dilate the region map; "
                            "default 1.0 for profiles, 0 for regions unless given")
    group.add_argument("--charge", action="append", metavar="NAME=q",
                       help="Charge override (e) by species label, residue name or atom name; repeatable")
    group.add_argument("--voltage", type=float, metavar="mV",
                       help="Known transmembrane voltage (positive drives cations along +axis); enables conductance")
    group.add_argument("--efield", type=_field, metavar="V/nm",
                       help="Applied static field (E along z, or EX,EY,EZ); V = E x mean box length")
    group.add_argument("--efield-from", metavar="FILE",
                       help="Read the applied field from a GROMACS .mdp or .tpr (the .tpr needs gmx on the PATH)")
    group.add_argument("--water-mapping", type=int, metavar="N",
                       help="Water molecules per Martini water bead (default 4)")
    group.add_argument("--permeation-stride", type=int, default=1, metavar="N",
                       help="Read every Nth trajectory frame for permeation (default 1: every saved frame)")
    group.add_argument("--msd-blocks", type=int, default=5, metavar="N",
                       help="Contiguous segments for the D_n (pf) interval (default 5)")
    group.add_argument("--permeation-only", action="store_true",
                       help="Run only the permeation analysis (needs --channel-extent or --permeation-region); "
                            "also the route for coarse-grained input, whose beads have no atomic radii")


def permeation_requested(args: argparse.Namespace) -> bool:
    """Whether the parsed arguments ask for the permeation analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command arguments.

    Returns
    -------
    bool
    """
    return bool(getattr(args, "permeation", None))


def validate_permeation_arguments(args: argparse.Namespace) -> None:
    """Check the permeation options before any trajectory is read.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command arguments.

    Raises
    ------
    ValueError
        For an inconsistent combination.
    """
    if getattr(args, "permeation_only", False) and not permeation_requested(args):
        raise ValueError("--permeation-only needs --permeation SPECIES")
    if not permeation_requested(args):
        for flag in ("channel_extent", "channel_radius", "permeation_pore", "permeation_region", "voltage", "efield",
                     "efield_from", "charge"):
            if getattr(args, flag, None):
                raise ValueError(f"--{flag.replace('_', '-')} needs --permeation SPECIES")
        return
    if not getattr(args, "topology", None):
        raise ValueError("--permeation needs an MD trajectory: pass --topology and one trajectory file")
    if args.permeation_only and not (args.channel_extent or args.permeation_region):
        raise ValueError("--permeation-only needs --channel-extent ZMIN,ZMAX or --permeation-region FILE")
    if (args.permeation_pore or args.channel_radius) and not args.channel_extent:
        raise ValueError("--permeation-pore and --channel-radius need --channel-extent")
    if args.channel_extent and args.permeation_region:
        raise ValueError("choose --channel-extent or --permeation-region")
    if args.voltage is not None and (args.efield or args.efield_from):
        raise ValueError("choose --voltage or an applied field (--efield/--efield-from)")
    if args.efield and args.efield_from:
        raise ValueError("choose --efield or --efield-from")
    if args.permeation_stride < 1:
        raise ValueError("--permeation-stride must be at least 1")


def _axis(args):
    value = getattr(args, "axis", "auto")
    if isinstance(value, str):
        named = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0)}
        if value.lower() in named:
            return named[value.lower()]
        if value.lower() == "auto":
            return None
        return _vector(value)
    return tuple(value)


def _channels(args, result):
    from .permeation import CylinderChannel, ProfileChannel, RegionChannel, region_channel_from_definition
    axis = _axis(args)
    if args.permeation_region:
        path = Path(args.permeation_region)
        tolerance = 0.0 if args.lumen_tolerance is None else args.lumen_tolerance
        if path.suffix.lower() == ".json":
            channel = region_channel_from_definition(path, tolerance=tolerance)
            if channel.region_definition["reference_frame"] != args.start:
                raise ValueError("the region is defined on frame "
                                 f"{channel.region_definition['reference_frame']}; use --start "
                                 f"{channel.region_definition['reference_frame']} so it is the reference frame")
            return [channel]
        if axis is None:
            raise ValueError("--permeation-region with a .dx map needs --axis x|y|z|X,Y,Z")
        return [RegionChannel(path.stem, path, axis, tolerance=tolerance)]
    if args.channel_extent:
        if axis is None:
            raise ValueError("--channel-extent needs --axis x|y|z|X,Y,Z (the extent is measured along it)")
        lo, hi = args.channel_extent
        points = args.permeation_pore or ([tuple(args.origin)] if getattr(args, "origin", None) else [None])
        return [CylinderChannel(f"pore{i + 1}" if len(points) > 1 else "channel", axis, lo, hi, point=p,
                                radius=args.channel_radius) for i, p in enumerate(points)]
    if result is None or not any(f.profile is not None for f in result.frames):
        raise ValueError("no channel profile was resolved, so the channel extent is unknown; give --channel-extent "
                         "ZMIN,ZMAX with --axis (and optionally --channel-radius)")
    profiles = [f.profile for f in result.frames]
    numbers = [f.frame_index for f in result.frames]
    tolerance = 1.0 if args.lumen_tolerance is None else args.lumen_tolerance
    return [ProfileChannel("channel", profiles, numbers, tolerance=tolerance)]


def _alignment_selection(args, universe):
    base = getattr(args, "selection", None) or "protein"
    for name in ("CA", "BB"):
        selection = f"({base}) and name {name}"
        try:
            if len(universe.select_atoms(selection)) >= 3:
                return selection
        except Exception:
            break
    return base


def run_permeation(args: argparse.Namespace, result=None) -> dict:
    """Run the permeation analysis for a trajectory command and write its files.

    Files are written next to ``-o`` with its stem as prefix
    (``STEM_permeation.json``, ``..._events.csv``, ``..._visits.csv``,
    ``..._frames.csv`` and, when figures are selected, ``..._permeation.png``).

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command arguments (see :func:`add_permeation_arguments`).
    result : crevice.models.TrajectoryAnalysis, optional
        The command's structural analysis; its profiles define the channel when
        no explicit extent or region is given.

    Returns
    -------
    dict
        The compact record stored in the trajectory JSON (``files``,
        ``coarse_grained``, per-species totals).
    """
    from .permeation import analyze_permeation, parse_charge_overrides, read_applied_field, write_permeation_outputs
    from .trajectory import _open_universe
    from .cli import _FrameProgress
    universe = _open_universe(args.topology, args.structures[0])
    channels = _channels(args, result)
    efield, source = args.efield, None
    if args.efield_from:
        record = read_applied_field(args.efield_from)
        if not record["static"]:
            raise ValueError(f"{args.efield_from}: the applied field is time dependent; conductance needs a static field")
        efield, source = record["E_V_per_nm"], record["source"]
    elif args.efield:
        source = "--efield"
    selection = _alignment_selection(args, universe)
    total = len(range(args.start, min(args.stop or len(universe.trajectory), len(universe.trajectory)),
                      args.permeation_stride))
    progress = _FrameProgress(total, "permeation frames")
    try:
        analysis = analyze_permeation(
            universe, species=args.permeation, channels=channels, start=args.start, stop=args.stop,
            stride=args.permeation_stride, align=args.align, alignment_selection=selection,
            charges=parse_charge_overrides(args.charge), water_mapping=args.water_mapping, voltage_mV=args.voltage,
            efield=efield, efield_source=source, msd_blocks=args.msd_blocks,
            bootstrap_replicates=args.bootstrap_replicates, bootstrap_seed=args.bootstrap_seed,
            confidence=args.confidence or 0.95, progress_callback=progress)
    finally:
        progress.close()
    output = Path(args.output)
    files = write_permeation_outputs(analysis, output.parent, output.stem, dpi=args.dpi)
    summary = analysis.summary
    for channel in summary["channels"]:
        name = channel["definition"]["name"]
        for row in channel["species"]:
            line = (f"permeation {name} {row['species']}: up={row['events_up']} down={row['events_down']} "
                    f"unresolved={row['unresolved_crossings']} mean_inside={row['occupancy']['mean']:.3g}")
            water = row.get("water")
            if water:
                pf = water["pf_cm3_per_s"]
                line += f" pd={water['pd_cm3_per_s']:.3g} pf={pf:.3g} cm3/s" if pf is not None else \
                        f" pd={water['pd_cm3_per_s']:.3g} cm3/s pf=unavailable"
            print(line)
        current = channel.get("ion_current_total")
        if current:
            g = current.get("conductance_pS")
            print(f"permeation {name} ion current={current['current_pA']:.3g} pA "
                  + (f"conductance={g:.3g} pS" if g is not None else "conductance undefined (no voltage)"))
    for warning in summary["warnings"]:
        print(f"crevice: warning: {warning}", file=sys.stderr)
    return {"files": {k: Path(v).name for k, v in files.items()}, "coarse_grained": summary["coarse_grained"],
            "totals_over_channels": summary["totals_over_channels"], "voltage": summary["voltage"],
            "channels": [c["definition"] for c in summary["channels"]]}


def write_permeation_only_record(args: argparse.Namespace, record: dict) -> None:
    """Write the ``-o`` JSON of a ``--permeation-only`` run.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments (``output``, ``topology``, ``structures``).
    record : dict
        From :func:`run_permeation`.
    """
    data = {"mode": "permeation-only", "topology": str(args.topology), "trajectory": str(args.structures[0]),
            "permeation": record, "coarse_grained": record["coarse_grained"]}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    from ._portable_paths import dumps_record  # 0.4: recorded paths relative to this file
    Path(args.output).write_text(dumps_record(data, args.output, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

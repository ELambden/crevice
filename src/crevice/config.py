"""Typed run settings shared by the ``crevice`` command line and Python.

Every option of every ``crevice`` command belongs to one section of
:class:`RunSettings`:

=================  ==========================================================
``input``          structures, trajectories, fetching, model/parser choice,
                   atom selections and residue lists
``geometry``       profile, probe, cast, cavity, contact and region settings
                   (lengths in Å, volumes in Å³)
``radii``          the atomic radius set (``--radii``)
``outputs``        output directory, prefix, explicit file paths and the
                   output policy (``--outputs``, ``--all``,
                   ``--legacy-outputs``, ``--dry-run``, ``--skip-*``)
``display``        figure resolution, annotation, guides and display
                   smoothing (display only; never changes a measurement)
``trajectory``     frame range and guards, alignment, periodic handling,
                   workers and statistics
``hydration``      water-contact, SASA and water-density settings
``method``         command-specific options without a shared meaning
=================  ==========================================================

A field name is the option's argparse destination (``--section-spacing`` is
``geometry.section_spacing``). Defaults are read from the parser itself
(:func:`crevice.cli.build_parser`), so the CLI and Python can never disagree
about a default and no numerical default is defined twice.

The same object drives both interfaces:

* CLI: ``crevice publish 1GRM --print-config`` prints the effective settings
  as JSON; ``crevice publish 1GRM --config run.json`` loads them (precedence:
  built-in default < ``--config`` file < command-line flag). JSON only: TOML
  would need an extra dependency on Python 3.10.
* Python: :func:`settings` builds a :class:`RunSettings`, :func:`run` executes
  it through exactly the code path of :func:`crevice.cli.main` and returns the
  manifest record.

:meth:`RunSettings.validate` checks units, ranges, input files and option
combinations before any structure is loaded or grid allocated. Every bundle
command records the effective settings, where each value came from and the
input file hashes in its manifest under ``run``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

#: Schema identifier written by :meth:`RunSettings.to_dict`.
SCHEMA = "crevice.settings/1"

#: Section names, in output order.
SECTIONS = ("input", "geometry", "radii", "outputs", "display", "trajectory", "hydration", "method")


def _f(unit: str | None = None, check: str | None = None, default: Any = None):
    """A settings field with its unit and range check recorded as metadata."""
    return field(default=default, metadata={"unit": unit, "check": check})


@dataclass
class InputSettings:
    """Inputs and atom selections.

    Attributes
    ----------
    structure, structures, topology, trajectory, definition : str or list of str
        Structure argument(s), MD topology/trajectory and region definition.
    left_statistics, right_statistics, left_definition, right_definition : str
        ``region-compare`` inputs.
    left_input_provenance, right_input_provenance : str
        ``region-compare`` historical input records.
    cache_dir, input_format, assembly, offline, refetch, timeout : object
        Fetching of accessions (``timeout`` in s).
    model, parser, md_selection, chain_ids : object
        Model index, mmCIF parser, GRO selection and mmCIF chain namespace.
    selection, obstacle_selection, water_selection : str
        MDAnalysis selections (``cast --selection`` is ``protein``/``all``).
    include_hydrogen, exclude_hetero, include_hetero : bool
        Obstacle atom policy of the profile family and cavity/tunnel searches.
    volume_dx, reference_volume_dx, cavity_results, scene : str
        Measured maps, earlier results and scenes the command builds on.
    residues_json, residue_groups, lining_evidence : str
        Residue lists and groups (JSON files).
    alignment_residue, alignment_residues_json, remove_residue : list of str or str
        Alignment scaffold and network-removal residues.
    landmarks_json, region_id : str
        ``region-init`` landmarks and region identifier.
    pdb_ids, name, source_dir, base_url : object
        Accessions and example datasets of ``fetch``/``benchmark``/``fetch-example``.
    outdir : str
        ``view`` input directory.
    """

    structure: Any = None
    structures: Any = None
    topology: Any = None
    trajectory: Any = None
    definition: Any = None
    left_statistics: Any = None
    right_statistics: Any = None
    left_definition: Any = None
    right_definition: Any = None
    left_input_provenance: Any = None
    right_input_provenance: Any = None
    cache_dir: Any = None
    input_format: Any = None
    assembly: Any = None
    offline: Any = None
    refetch: Any = None
    timeout: Any = _f("s", "positive")
    model: Any = _f(None, "nonnegative")
    parser: Any = None
    md_selection: Any = None
    chain_ids: Any = None
    selection: Any = None
    obstacle_selection: Any = None
    water_selection: Any = None
    include_hydrogen: Any = None
    exclude_hetero: Any = None
    include_hetero: Any = None
    volume_dx: Any = _f(None, "file")
    reference_volume_dx: Any = _f(None, "file")
    cavity_results: Any = _f(None, "file")
    scene: Any = _f(None, "file")
    residues_json: Any = _f(None, "file")
    residue_groups: Any = _f(None, "file")
    lining_evidence: Any = _f(None, "file")
    alignment_residue: Any = None
    alignment_residues_json: Any = _f(None, "file")
    remove_residue: Any = None
    landmarks_json: Any = _f(None, "file")
    region_id: Any = None
    pdb_ids: Any = None
    name: Any = None
    source_dir: Any = None
    base_url: Any = None
    outdir: Any = None


@dataclass
class GeometrySettings:
    """Geometric method settings: lengths in Å, volumes in Å³, fractions in [0, 1].

    Attributes
    ----------
    axis, origin, section_spacing, enclosure_radius, lateral_exits : object
        Channel direction and seed, cross-section spacing (Å) and enclosure probe (Å; ``None`` = automatic).
    exit_bulk_radius, exit_spacing, samples, padding, search_radius, refinement_steps : object
        Exit search (Å), profile sampling and centre refinement.
    probe_radius : float
        Probe radius (Å); subtracted from profile radii, or the rolling probe of casts.
    spacing, outer_radius, enclosure_fraction, min_depth, min_component_volume : float
        Cast and region grid settings.
    max_components, cavity_selection, core_fraction, max_grid_points, seed, seed_tolerance : object
        Cast component selection, grid guard and seed.
    cast_mode, cast_outer_radius, cast_enclosure_fraction, cast_extension, cast_spacing : object
        ``publish`` cast settings.
    cast_min_radius, cast_max_components, cast_selection, cast_core_fraction, cast_max_grid_points : object
        ``publish`` cast settings.
    focus, focus_radius, max_component_volume, entry_end, lining_cutoff : object
        ``publish`` crops, segment labels and lining-residue distance (Å).
    cutoff, contact_cutoff, contact_metric, lining_distance, region : object
        Residue contact settings (Å).
    geometry_mode, region_margin, grid_phase, profile_padding, minimum_overlap, ambiguity_ratio : object
        Cavity-trajectory geometry.
    min_radius, max_cavities, max_tunnels, tunnel_start : object
        Cavity and tunnel searches (``tunnel_start`` is ``--start`` of ``tunnels``/``analyze``).
    cavity_spacing, cavity_min_radius, cavity_max_grid_points : object
        ``analyze``/``static-suite`` cavity search.
    tunnel_spacing, tunnel_min_radius, tunnel_max_grid_points : object
        ``analyze``/``static-suite`` tunnel search.
    radius : float
        ``region-init`` landmark sphere radius (Å).
    hydrophobicity_bandwidth, lining_shell, electrolyte, resistivity : object
        ``fingerprint`` settings.
    """

    axis: Any = None
    origin: Any = None
    section_spacing: Any = _f("Å", "positive")
    enclosure_radius: Any = _f("Å", "nonnegative")
    lateral_exits: Any = None
    exit_bulk_radius: Any = _f("Å", "positive")
    exit_spacing: Any = _f("Å", "positive")
    samples: Any = _f(None, "count2")
    padding: Any = _f("Å", "nonnegative")
    search_radius: Any = _f("Å", "nonnegative")
    refinement_steps: Any = _f(None, "nonnegative")
    probe_radius: Any = _f("Å", "nonnegative")
    spacing: Any = _f("Å", "positive")
    outer_radius: Any = _f("Å", "positive")
    enclosure_fraction: Any = _f(None, "fraction")
    min_depth: Any = _f("Å", "nonnegative")
    min_component_volume: Any = _f("Å³", "nonnegative")
    max_components: Any = _f(None, "count1")
    cavity_selection: Any = None
    core_fraction: Any = _f(None, "fraction")
    max_grid_points: Any = _f(None, "count1")
    seed: Any = None
    seed_tolerance: Any = _f("Å", "nonnegative")
    cast_mode: Any = None
    cast_outer_radius: Any = _f("Å", "positive")
    cast_enclosure_fraction: Any = _f(None, "fraction")
    cast_extension: Any = _f("Å", "nonnegative")
    cast_spacing: Any = _f("Å", "positive")
    cast_min_radius: Any = _f("Å", "nonnegative")
    cast_max_components: Any = _f(None, "count1")
    cast_selection: Any = None
    cast_core_fraction: Any = _f(None, "fraction")
    cast_max_grid_points: Any = _f(None, "count1")
    focus: Any = None
    focus_radius: Any = _f("Å", "positive")
    max_component_volume: Any = _f("Å³", "nonnegative")
    entry_end: Any = None
    lining_cutoff: Any = _f("Å", "positive")
    cutoff: Any = _f("Å", "positive")
    contact_cutoff: Any = _f("Å", "positive")
    contact_metric: Any = None
    lining_distance: Any = _f("Å", "nonnegative")
    region: Any = None
    geometry_mode: Any = None
    region_margin: Any = _f("Å", "nonnegative")
    grid_phase: Any = None
    profile_padding: Any = _f("Å", "nonnegative")
    minimum_overlap: Any = _f(None, "fraction_open")
    ambiguity_ratio: Any = _f(None, "fraction")
    min_radius: Any = _f("Å", "nonnegative")
    max_cavities: Any = _f(None, "count1")
    max_tunnels: Any = _f(None, "count1")
    tunnel_start: Any = None
    cavity_spacing: Any = _f("Å", "positive")
    cavity_min_radius: Any = _f("Å", "nonnegative")
    cavity_max_grid_points: Any = _f(None, "count1")
    tunnel_spacing: Any = _f("Å", "positive")
    tunnel_min_radius: Any = _f("Å", "nonnegative")
    tunnel_max_grid_points: Any = _f(None, "count1")
    radius: Any = _f("Å", "positive")
    hydrophobicity_bandwidth: Any = _f("Å", "positive")
    lining_shell: Any = _f("Å", "nonnegative")
    electrolyte: Any = None
    resistivity: Any = _f("Ω m", "positive")


@dataclass
class RadiiSettings:
    """Atomic radius set.

    Attributes
    ----------
    radii : str or None
        Preset name or radius file (``None`` = the default set).
    allow_radii_mismatch : bool
        Proceed when the set differs from a region definition or earlier result.
    """

    radii: Any = None
    allow_radii_mismatch: Any = None


@dataclass
class OutputSettings:
    """Where and what to write.

    Attributes
    ----------
    out_dir, output, prefix : str
        Output directory, single output file and file-name prefix.
    outputs, all_outputs, legacy_outputs, dry_run : object
        Output policy (:mod:`crevice.outputs`) and plan-only runs.
    skip_hydration, skip_cavities, skip_tunnels, skip_network : bool
        Remove an optional analysis.
    fingerprint : bool
        ``publish --fingerprint``.
    csv, pdb, png, annotated_png, json, pymol, pymol_png, chord_png : object
        Explicit file paths of the single-file commands (``analyze --png`` is a switch).
    connectivity_json, connectivity_png, report_json, distribution_csv, network_json : str
        Explicit file paths.
    overwrite, list, open, cdn, show, hide, kind, snapshot, fetch, format : object
        Utility-command options.
    """

    out_dir: Any = None
    output: Any = None
    prefix: Any = None
    outputs: Any = None
    all_outputs: Any = None
    legacy_outputs: Any = None
    dry_run: Any = None
    skip_hydration: Any = None
    skip_cavities: Any = None
    skip_tunnels: Any = None
    skip_network: Any = None
    fingerprint: Any = None
    csv: Any = None
    pdb: Any = None
    png: Any = None
    annotated_png: Any = None
    json: Any = None
    pymol: Any = None
    pymol_png: Any = None
    chord_png: Any = None
    connectivity_json: Any = None
    connectivity_png: Any = None
    report_json: Any = None
    distribution_csv: Any = None
    network_json: Any = None
    overwrite: Any = None
    list: Any = None
    open: Any = None
    cdn: Any = None
    show: Any = None
    hide: Any = None
    kind: Any = None
    snapshot: Any = None
    fetch: Any = None
    format: Any = None


@dataclass
class DisplaySettings:
    """Display-only settings (never change a measurement).

    Attributes
    ----------
    dpi : int
        Figure resolution (dots per inch).
    annotate, mouth_guides, exit_centre_lines : bool
        Text in figures/scenes and optional guides.
    smooth, surface_smoothing : float
        Display-surface smoothing (Å).
    residue_label_every, stick_residues, view, profile_quantity : object
        Figure and scene content choices.
    """

    dpi: Any = _f("dpi", "count1")
    annotate: Any = None
    mouth_guides: Any = None
    exit_centre_lines: Any = None
    smooth: Any = _f("Å", "nonnegative")
    surface_smoothing: Any = _f("Å", "nonnegative")
    residue_label_every: Any = _f(None, "count1")
    stick_residues: Any = None
    view: Any = None
    profile_quantity: Any = None


@dataclass
class TrajectorySettings:
    """Frames, alignment, periodic handling, workers and statistics.

    Attributes
    ----------
    start, stop, stride, max_frames : int
        Frame range (``start``/``stop`` are source frame indices) and the frame guard.
    pbc, workers, align, inspect, profile_errors, probe_fallback : object
        Periodic handling, worker processes, alignment and error policy.
    confidence, block_length, bootstrap_replicates, bootstrap_seed : object
        Approximate mean intervals.
    occupancy, quantiles, analysis : object
        Contact occupancy threshold, frame quantile band and analysis list.
    """

    start: Any = _f(None, "nonnegative")
    stop: Any = _f(None, "nonnegative")
    stride: Any = _f(None, "count1")
    max_frames: Any = _f(None, "nonnegative")
    pbc: Any = None
    workers: Any = _f(None, "count1")
    align: Any = None
    inspect: Any = None
    profile_errors: Any = None
    probe_fallback: Any = None
    confidence: Any = _f(None, "confidence")
    block_length: Any = _f("frames", "count1")
    bootstrap_replicates: Any = _f(None, "count1")
    bootstrap_seed: Any = None
    occupancy: Any = _f(None, "fraction")
    quantiles: Any = None
    analysis: Any = None


@dataclass
class HydrationSettings:
    """Hydration, water density and interaction settings.

    Attributes
    ----------
    hydration_cutoff, hydration_probe : float
        Water-oxygen contact distance and SASA probe (Å).
    hydration_sasa_points : int
        Quadrature points per atom (at least 32).
    water_density_spacing, water_density_smoothing, water_density_level : float
        Density bin (Å), display smoothing (Å) and display level (oxygens per Å³).
    skip_water_density, skip_interactions, skip_analysis_views, water_membership, no_water : bool
        Optional parts of the hydration analysis.
    """

    hydration_cutoff: Any = _f("Å", "positive")
    hydration_probe: Any = _f("Å", "positive")
    hydration_sasa_points: Any = _f(None, "count32")
    water_density_spacing: Any = _f("Å", "positive")
    water_density_smoothing: Any = _f("Å", "nonnegative")
    water_density_level: Any = _f("1/Å³", "positive")
    skip_water_density: Any = None
    skip_interactions: Any = None
    skip_analysis_views: Any = None
    water_membership: Any = None
    no_water: Any = None


_SECTION_TYPES = {"input": InputSettings, "geometry": GeometrySettings, "radii": RadiiSettings,
                  "outputs": OutputSettings, "display": DisplaySettings, "trajectory": TrajectorySettings,
                  "hydration": HydrationSettings}
#: argparse destination -> (section, field name).
FIELD_OF: dict = {f.name: (section, f.name) for section, cls in _SECTION_TYPES.items() for f in fields(cls)}
#: ``--start`` is a coordinate for these commands, a frame index elsewhere.
_TUNNEL_START_COMMANDS = {"tunnels", "analyze"}
#: Never part of the settings.
_EXCLUDED = {"func", "command", "quiet", "debug", "config", "print_config", "help"}


def _section_of(command: str, dest: str) -> tuple:
    if dest == "start" and command in _TUNNEL_START_COMMANDS:
        return "geometry", "tunnel_start"
    return FIELD_OF.get(dest, ("method", dest))


def _dest_of(command: str, section: str, name: str) -> str:
    if section == "geometry" and name == "tunnel_start" and command in _TUNNEL_START_COMMANDS:
        return "start"
    return name


def _subparser(command: str) -> argparse.ArgumentParser:
    from .cli import build_parser
    parser = build_parser()
    action = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    if command not in action.choices:
        raise ValueError(f"unknown crevice command: {command}")
    return action.choices[command]


def command_actions(command: str) -> dict:
    """The argparse actions of ``command``, by destination (settings options only).

    Parameters
    ----------
    command : str
        Command name.

    Returns
    -------
    dict of str to argparse.Action
    """
    return {a.dest: a for a in _subparser(command)._actions
            if a.dest not in _EXCLUDED and a.dest != argparse.SUPPRESS}


def _jsonable(value):
    if isinstance(value, tuple):
        return [_jsonable(v) for v in value]
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


@dataclass
class RunSettings:
    """The complete, typed settings of one ``crevice`` command run.

    Attributes
    ----------
    command : str
        Command name (``publish``, ``cast``, ...).
    input, geometry, radii, outputs, display, trajectory, hydration : dataclass
        Section dataclasses (:class:`InputSettings` ... :class:`HydrationSettings`).
    method : dict
        Command-specific options without a shared section.
    applicable : frozenset of str
        Argparse destinations the command accepts; only these are serialised.
    sources : dict of str to str
        Where each value came from: ``default``, ``config``, ``cli`` or ``python``.
    config_file : str or None
        The ``--config`` file, when one was used.
    """

    command: str
    input: InputSettings = field(default_factory=InputSettings)
    geometry: GeometrySettings = field(default_factory=GeometrySettings)
    radii: RadiiSettings = field(default_factory=RadiiSettings)
    outputs: OutputSettings = field(default_factory=OutputSettings)
    display: DisplaySettings = field(default_factory=DisplaySettings)
    trajectory: TrajectorySettings = field(default_factory=TrajectorySettings)
    hydration: HydrationSettings = field(default_factory=HydrationSettings)
    method: dict = field(default_factory=dict)
    applicable: frozenset = frozenset()
    sources: dict = field(default_factory=dict)
    config_file: str | None = None

    # -- construction -------------------------------------------------------
    @classmethod
    def defaults(cls, command: str) -> "RunSettings":
        """Settings with every option at its built-in (parser) default.

        Parameters
        ----------
        command : str
            Command name.

        Returns
        -------
        RunSettings
        """
        values = {dest: action.default for dest, action in command_actions(command).items()}
        return cls.from_values(command, values, sources={dest: "default" for dest in values})

    @classmethod
    def from_values(cls, command: str, values: dict, *, sources: dict | None = None,
                    config_file: str | None = None) -> "RunSettings":
        """Build settings from a flat ``{argparse destination: value}`` mapping.

        Parameters
        ----------
        command : str
            Command name.
        values : dict
            Values by argparse destination (for example a parsed namespace's ``vars``).
        sources : dict, optional
            Source of each value.
        config_file : str, optional
            The configuration file used.

        Returns
        -------
        RunSettings
        """
        settings = cls(command)
        applicable = set()
        for dest, value in values.items():
            if dest in _EXCLUDED or dest.startswith("_"):
                continue
            section, name = _section_of(command, dest)
            if section == "method":
                settings.method[name] = value
            else:
                setattr(getattr(settings, section), name, value)
            applicable.add(dest)
        settings.applicable = frozenset(applicable)
        settings.sources = dict(sources or {})
        settings.config_file = config_file
        return settings

    @classmethod
    def from_namespace(cls, args: argparse.Namespace, *, sources: dict | None = None,
                       config_file: str | None = None) -> "RunSettings":
        """Settings of a parsed ``crevice`` command line.

        Parameters
        ----------
        args : argparse.Namespace
            Result of ``build_parser().parse_args(...)``.
        sources : dict, optional
            Source of each value (see :func:`setting_sources`).
        config_file : str, optional
            The ``--config`` file, when used.

        Returns
        -------
        RunSettings
        """
        actions = command_actions(args.command)
        values = {dest: getattr(args, dest) for dest in actions if hasattr(args, dest)}
        return cls.from_values(args.command, values, sources=sources, config_file=config_file)

    @classmethod
    def from_dict(cls, data: dict) -> "RunSettings":
        """Rebuild settings from :meth:`to_dict` output (or a ``--config`` file).

        Missing options take their defaults; unknown options are an error.

        Parameters
        ----------
        data : dict
            ``{"schema": ..., "command": ..., "<section>": {...}, ...}``.

        Returns
        -------
        RunSettings

        Raises
        ------
        ValueError
            Unknown command, section or option.
        """
        command = data.get("command")
        if not command:
            raise ValueError("settings need a 'command'")
        settings = cls.defaults(command)
        values = flatten_settings(data, command)
        for dest, value in values.items():
            settings.set(dest, value, source="config")
        return settings

    # -- access ---------------------------------------------------------------
    def get(self, dest: str):
        """Value of the option with argparse destination ``dest``.

        Parameters
        ----------
        dest : str
            Argparse destination.

        Returns
        -------
        object
        """
        section, name = _section_of(self.command, dest)
        return self.method.get(name) if section == "method" else getattr(getattr(self, section), name)

    def set(self, dest: str, value, *, source: str = "python") -> None:
        """Set an option by argparse destination, recording its ``source``.

        Parameters
        ----------
        dest : str
            Argparse destination.
        value : object
            New value (lists given for coordinate options become tuples).
        source : str, default "python"
            Recorded source.

        Raises
        ------
        ValueError
            When ``crevice <command>`` has no such option.
        """
        if dest not in self.applicable:
            raise ValueError(f"crevice {self.command} has no option '{dest}'")
        value = _normalise(self.command, dest, value)
        section, name = _section_of(self.command, dest)
        if section == "method":
            self.method[name] = value
        else:
            setattr(getattr(self, section), name, value)
        self.sources[dest] = source

    def values(self) -> dict:
        """Flat ``{argparse destination: value}`` of every applicable option.

        Returns
        -------
        dict
        """
        return {dest: self.get(dest) for dest in sorted(self.applicable)}

    def to_namespace(self) -> argparse.Namespace:
        """An :class:`argparse.Namespace` identical to what the parser would produce.

        Returns
        -------
        argparse.Namespace
            Includes ``command`` and the handler ``func``.
        """
        parser = _subparser(self.command)
        namespace = argparse.Namespace(command=self.command, quiet=False, debug=False, config=None,
                                       print_config=False, **self.values())
        namespace.func = parser.get_default("func")
        return namespace

    def to_dict(self) -> dict:
        """JSON-ready record: ``schema``, ``command`` and one mapping per section.

        Only the options the command accepts are included. This is the format
        of ``--print-config`` and ``--config``.

        Returns
        -------
        dict
        """
        out = {"schema": SCHEMA, "command": self.command}
        for section in SECTIONS:
            block = {}
            for dest in sorted(self.applicable):
                sec, name = _section_of(self.command, dest)
                if sec == section:
                    block[name] = _jsonable(self.get(dest))
            if block:
                out[section] = block
        return out

    def units(self) -> dict:
        """Unit of every applicable numeric option that has one (``section.name`` keys).

        Returns
        -------
        dict of str to str
        """
        out = {}
        for dest in sorted(self.applicable):
            section, name = _section_of(self.command, dest)
            if section in _SECTION_TYPES:
                meta = _SECTION_TYPES[section].__dataclass_fields__[name].metadata
                if meta.get("unit"):
                    out[f"{section}.{name}"] = meta["unit"]
        return out

    # -- validation -----------------------------------------------------------
    def validate(self) -> "RunSettings":
        """Check units, ranges, input files and option combinations.

        Runs before any structure is loaded, grid allocated or file written.

        Returns
        -------
        RunSettings
            ``self``, for chaining.

        Raises
        ------
        ValueError
            The first problem found, naming the option.
        FileNotFoundError
            A required input file (map, residue list, definition) is missing.
        """
        for dest in sorted(self.applicable):
            section, name = _section_of(self.command, dest)
            if section not in _SECTION_TYPES:
                continue
            check = _SECTION_TYPES[section].__dataclass_fields__[name].metadata.get("check")
            if check:
                _check(self.command, dest, self.get(dest), check)
        self._check_combinations()
        return self

    def _check_combinations(self) -> None:
        command, get, has = self.command, self.get, self.applicable.__contains__
        if has("hydration_sasa_points") and get("hydration_sasa_points") is not None and get("hydration_sasa_points") < 32:
            raise ValueError("--hydration-sasa-points: hydration SASA requires at least 32 points per atom")
        if command in {"hydration", "cavity-trajectory", "region-trajectory"} and get("max_frames") == 0:
            raise ValueError("--max-frames must be at least 1 for this command")
        if has("start") and has("stop") and command not in _TUNNEL_START_COMMANDS:
            start, stop = get("start"), get("stop")
            if stop is not None and start is not None and stop <= start:
                raise ValueError(f"--stop ({stop}) must be greater than --start ({start})")
        if has("quantiles") and get("quantiles") is not None:
            low, high = get("quantiles")
            if not 0 <= low < high <= 1:
                raise ValueError("--quantiles must satisfy 0 <= lower < upper <= 1")
        if command == "hydration" and get("confidence") and get("bootstrap_replicates") is not None \
                and get("bootstrap_replicates") < 200:
            raise ValueError("--bootstrap-replicates: hydration intervals need at least 200 replicates")
        if command == "cavity-trajectory" and get("obstacle_selection") and get("geometry_mode") != "reference-region":
            raise ValueError("--obstacle-selection requires --geometry-mode reference-region")
        if command == "hydration" and get("cavity_results") and not get("trajectory"):
            raise ValueError("--cavity-results requires --trajectory")
        if command == "trajectory":
            if get("alignment_residue") and get("alignment_residues_json"):
                raise ValueError("Choose --alignment-residue or --alignment-residues-json")
            if get("topology") and len(get("structures") or ()) != 1:
                raise ValueError("--topology requires exactly one trajectory file")
        if command == "publish" and get("cast_mode") == "rolling" and (
                get("focus") is not None or get("max_component_volume") is not None or get("cast_min_radius")):
            raise ValueError("--cast-mode rolling does not support --focus, --max-component-volume or --cast-min-radius")
        if command == "publish" and get("max_component_volume") is not None and \
                get("max_component_volume") < (get("min_component_volume") or 0):
            raise ValueError("--max-component-volume must be at least --min-component-volume")
        if command == "profile" and get("pymol") and not get("pdb"):
            raise ValueError("--pymol requires --pdb so the dummy atoms can be loaded")


def _check(command: str, dest: str, value, check: str) -> None:
    if value is None or isinstance(value, bool):
        return
    if check == "file":
        if not Path(value).exists():
            raise FileNotFoundError(f"--{dest.replace('_', '-')}: {value} does not exist")
        return
    option = "--" + ("start" if dest == "tunnel_start" else dest.replace("_", "-"))
    if not isinstance(value, (int, float)):
        return
    if not math.isfinite(value):
        raise ValueError(f"{option} must be finite (got {value})")
    rules = {"positive": (lambda v: v > 0, "positive"), "nonnegative": (lambda v: v >= 0, "non-negative"),
             "count1": (lambda v: v >= 1, "at least 1"), "count2": (lambda v: v >= 2, "at least 2"),
             "count32": (lambda v: v >= 32, "at least 32"), "fraction": (lambda v: 0 <= v <= 1, "in [0, 1]"),
             "fraction_open": (lambda v: 0 < v <= 1, "in (0, 1]"), "confidence": (lambda v: 0 <= v < 1, "in [0, 1)")}
    ok, text = rules[check]
    if not ok(value):
        raise ValueError(f"{option} must be {text} (got {value})")


def _normalise(command: str, dest: str, value):
    """Turn JSON/Python values into what argparse would have produced for ``dest``."""
    action = command_actions(command).get(dest)
    if action is None or value is None:
        return value
    from .geometry import parse_coord
    if isinstance(value, list) and action.type is parse_coord:
        return tuple(float(v) for v in value)
    if action.type is float and isinstance(value, int) and not isinstance(value, bool):
        return float(value)
    if action.type is float and isinstance(value, list):
        return [float(v) if isinstance(v, int) and not isinstance(v, bool) else v for v in value]
    return value


def flatten_settings(data: dict, command: str | None = None) -> dict:
    """Flatten a settings record (sections) into ``{argparse destination: value}``.

    Parameters
    ----------
    data : dict
        :meth:`RunSettings.to_dict` output. A flat mapping of destinations is
        also accepted.
    command : str, optional
        Command the settings are for; must match ``data["command"]`` when both
        are given.

    Returns
    -------
    dict

    Raises
    ------
    ValueError
        Wrong command, unknown section, or an option the command does not accept.
    """
    command = command or data.get("command")
    if data.get("command") not in (None, command):
        raise ValueError(f"the settings are for 'crevice {data['command']}', not 'crevice {command}'")
    if data.get("schema") not in (None, SCHEMA):
        raise ValueError(f"unsupported settings schema {data.get('schema')!r} (expected {SCHEMA})")
    actions = command_actions(command)
    values = {}
    for key, block in data.items():
        if key in ("schema", "command"):
            continue
        if key in SECTIONS and isinstance(block, dict):
            for name, value in block.items():
                values[_dest_of(command, key, name)] = value
        elif key in SECTIONS:
            raise ValueError(f"settings section '{key}' must be a mapping")
        else:
            values[key] = value_or_error(key, block, actions, command)
    unknown = sorted(d for d in values if d not in actions)
    if unknown:
        raise ValueError(f"crevice {command} has no option(s): {', '.join(unknown)}")
    return {dest: _normalise(command, dest, value) for dest, value in values.items()}


def value_or_error(key: str, value, actions: dict, command: str):
    """Return ``value`` for a flat (sectionless) option ``key``; reject unknown sections.

    Parameters
    ----------
    key : str
        Option destination.
    value : object
        Its value.
    actions : dict
        The command's actions (:func:`command_actions`).
    command : str
        Command name, for the message.

    Returns
    -------
    object

    Raises
    ------
    ValueError
        When ``key`` is neither a section nor an option of the command.
    """
    if key not in actions:
        raise ValueError(f"crevice {command}: unknown settings key '{key}' (sections: {', '.join(SECTIONS)})")
    return value


def load_settings_file(path: str | Path) -> dict:
    """Read a JSON settings file written by ``--print-config`` (or by hand).

    Parameters
    ----------
    path : str or Path
        JSON file.

    Returns
    -------
    dict

    Raises
    ------
    ValueError
        The file is not a JSON object (TOML is not supported: it would need an
        extra dependency on Python 3.10).
    """
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        hint = " (TOML is not supported; use the JSON written by --print-config)" if path.suffix == ".toml" else ""
        raise ValueError(f"--config {path}: not valid JSON{hint}: {error}") from None
    if not isinstance(data, dict):
        raise ValueError(f"--config {path}: expected a JSON object")
    return data


def setting_sources(command: str, argv: list, config_values: dict | None = None) -> dict:
    """Where each option's value came from: ``cli``, ``config`` or ``default``.

    Parameters
    ----------
    command : str
        Command name.
    argv : list of str
        The command line (without the program name).
    config_values : dict, optional
        Flattened ``--config`` values.

    Returns
    -------
    dict of str to str
    """
    tokens = [t.split("=", 1)[0] for t in argv if t.startswith("-")]
    sources = {}
    for dest, action in command_actions(command).items():
        if not action.option_strings:
            sources[dest] = "cli"
        elif any(option in tokens for option in action.option_strings):
            sources[dest] = "cli"
        elif config_values and dest in config_values:
            sources[dest] = "config"
        else:
            sources[dest] = "default"
    return sources


def file_record(path: str | Path, role: str) -> dict:
    """Provenance of one input file: resolved path, size and SHA-256.

    Parameters
    ----------
    path : str or Path
        Input file.
    role : str
        Its role (``structure``, ``trajectory``, ``volume_dx`` ...).

    Returns
    -------
    dict
    """
    path = Path(path)
    record = {"role": role, "argument": str(path)}
    if path.is_file():
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        record.update(path=str(path.resolve()), bytes=path.stat().st_size, sha256=digest.hexdigest())
    else:
        record["status"] = "not_a_local_file (an accession resolved at run time; see the input report)"
    return record


#: Input options whose files are hashed into the run record.
INPUT_FILE_OPTIONS = ("structure", "structures", "topology", "trajectory", "definition", "volume_dx",
                      "reference_volume_dx", "cavity_results", "residues_json", "residue_groups", "lining_evidence",
                      "alignment_residues_json", "scene", "radii", "left_statistics", "right_statistics")


def run_record(settings: RunSettings, argv: list | None, plan=None, log=None, source=None) -> dict:
    """The ``run`` block written into a command's manifest.

    Parameters
    ----------
    settings : RunSettings
        Effective settings.
    argv : list of str or None
        Command line (``None`` for a Python :func:`run`).
    plan : crevice.outputs.OutputPlan, optional
        Resolved output plan.
    log : crevice.outputs.OutputLog, optional
        Skipped outputs and notes.
    source : dict, optional
        ``StructureSource.to_dict()`` of a resolved structure argument.

    Returns
    -------
    dict
    """
    from . import __version__
    inputs = []
    for dest in INPUT_FILE_OPTIONS:
        if dest not in settings.applicable:
            continue
        value = settings.get(dest)
        if not value:
            continue
        if dest == "radii" and Path(str(value)).suffix.lower() not in {".json", ".csv", ".rad"}:
            continue
        for item in (value if isinstance(value, (list, tuple)) else [value]):
            path = Path(str(item))
            if dest == "structure" and source and not path.is_file() and source.get("path"):
                path = Path(source["path"])
            inputs.append(file_record(path, dest))
    record = {"crevice_version": __version__, "command": settings.command,
              "argv": list(argv) if argv is not None else None,
              "interface": "cli" if argv is not None else "python",
              "config_file": settings.config_file, "settings": settings.to_dict(),
              "non_default_settings": {d: s for d, s in sorted(settings.sources.items()) if s != "default"},
              "units": settings.units(), "inputs": inputs,
              "python": sys.version.split()[0]}
    if source:
        record["structure_source"] = source
    if plan is not None:
        outputs = plan.to_dict()
        if log is not None:
            outputs["not_written_count"] = len(log.blocked)
            if log.notes:
                outputs["notes"] = dict(log.notes)
        record["outputs"] = outputs
    return record


def settings(command: str, *args, **values) -> RunSettings:
    """Build validated settings for ``crevice <command>`` from Python.

    Positional arguments fill the command's positional inputs in order;
    keyword arguments are argparse destinations (``section_spacing=0.25``) or
    whole sections (``geometry={"samples": 101}``). Unspecified options take
    the CLI defaults.

    Parameters
    ----------
    command : str
        Command name, for example ``"publish"``.
    *args : str
        Positional inputs (for ``publish``: the structure).
    **values : object
        Options by destination or section.

    Returns
    -------
    RunSettings

    Raises
    ------
    ValueError
        Unknown options, invalid values or combinations.

    Examples
    --------
    >>> from crevice.config import settings
    >>> s = settings("publish", "1GRM", out_dir="out", samples=101, outputs="figures")
    >>> s.geometry.samples, s.outputs.outputs, s.geometry.section_spacing
    (101, 'figures', 0.5)
    """
    result = RunSettings.defaults(command)
    positional = [d for d, a in command_actions(command).items() if not a.option_strings]
    if len(args) > len(positional):
        raise ValueError(f"crevice {command} takes {len(positional)} positional input(s)")
    for dest, value in zip(positional, args):
        if dest == "structures" and isinstance(value, (str, Path)):
            value = [value]
        result.set(dest, str(value) if isinstance(value, Path) else value)
    for key, value in values.items():
        if key in SECTIONS and isinstance(value, dict):
            for name, item in value.items():
                result.set(_dest_of(command, key, name), item)
        else:
            result.set(key, value)
    missing = [d for d in positional if result.get(d) in (None, [])]
    required = [d for d, a in command_actions(command).items() if a.option_strings and a.required]
    missing += [d for d in required if result.get(d) is None]
    if missing:
        raise ValueError(f"crevice {command} needs: {', '.join(missing)}")
    return result.validate()


def run(run_settings: RunSettings) -> dict:
    """Execute ``run_settings`` exactly as the ``crevice`` command line would.

    The same handler, output plan, radius set, annotation and progress
    settings are used as by :func:`crevice.cli.main`; errors propagate as
    exceptions instead of an exit status.

    Parameters
    ----------
    run_settings : RunSettings
        Settings from :func:`settings` or :meth:`RunSettings.from_dict`.

    Returns
    -------
    dict
        ``manifest`` (path or ``None``), ``files`` (the manifest's file index,
        empty when the command writes no manifest) and ``plan`` (the resolved
        output plan record, or ``None``).
    """
    from .cli import execute
    run_settings.validate()
    return execute(run_settings.to_namespace(), settings=run_settings, argv=None)

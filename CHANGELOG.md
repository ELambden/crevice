# Changelog

All notable changes to CREVICE are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) (see the version
policy in [CONTRIBUTING.md](CONTRIBUTING.md#versioning-and-releases)).

Entries describe software behaviour. Passing tests are software and synthetic
checks, not real-system validation or functional inference. No bundled
benchmark is curated, so numbers produced by CREVICE are not validated
biological results.

## [Unreleased]

### Added

- Lateral exit legs (`--lateral-exits`) are cast: each leg gets a probe-swept
  cast of kind `lateral_exit` (`crevice.channel_exits.lateral_exit_casts`),
  bounded by the capped mouth plane and by bulk solvent, drawn in the viewer
  scenes joined to the channel cast so the widening towards each portal is
  visible. Leg volumes are reported separately (`metadata.lateral_exit_casts`
  in the cast JSON and manifest, `lateral_exit_<end>_<leg>` rows in
  `*_void_cast.csv`, `profile_status.lateral_exit_cast_volumes_A3`,
  `*_exit_casts.dx`); the axial channel volume and profile are unchanged.
- `--mouth-guides` and `--exit-centre-lines` (and
  `crevice.presentation.display_guides`) opt back in to mouth rings, profile
  mouth bands and exit-leg centre lines.

- AlphaFold DB identifiers (`AF-P02920-F1`, `AF-P02920-F1-model_v6` or the
  download file name) are accepted wherever a structure is, and by `fetch`.
  Models come from AlphaFold DB (the latest version, looked up through its API,
  unless one is given), are cached under their AlphaFold DB name with
  provenance, and are checked against `_entry.id` or the `DBREF` UniProt
  accession. Input reports flag AlphaFold files as predicted models whose
  B-factor column is pLDDT. New worked example: an AlphaFold DB model.
- `crevice.radii.write_radius_template` writes an editable CSV radius file
  starting from the default set, and `RadiusSet.has_radius`
  tells whether a set gives an atom a radius.
- Flat CSV copies of the nested summary statistics that were JSON only:
  `*_cavity_statistics_summary.csv` (`cavity-trajectory`,
  `region-trajectory`), `*_water_membership_summary.csv` and
  `*_water_membership_strata.csv` (`--water-membership`),
  `*_hydration_summary.csv` (`hydration`) and `comparison_summary.csv`
  (`region-compare`). Each has one row per interval statistic (volume, each
  residue's metric, each contact or shared-water pair, each count series):
  `quantity` (JSON path), `item` (residue or pair), frame counts, mean, SD,
  median, 2.5/97.5% quantiles, extremes, half means, statistical
  inefficiency, effective frames, interval status/method/level and the
  pointwise and simultaneous mean bounds, plus one `value` row per
  correlation. Python: `crevice.io.summary_statistics_rows` and
  `write_summary_statistics_csv`. The JSON files are unchanged.

- Atomic radius sets throughout: every command that measures geometry takes
  `--radii {default,bondi,hole,charmm_like,PATH}` and every Python analysis
  that uses atomic radii takes `radii=` (a preset, a radius file or a
  `RadiusSet`; `None` = the set in effect, by default `DEFAULT_RADII`). Radius
  files may be JSON, CSV or HOLE `.rad` (`VDWR` records read with HOLE's
  column layout and first-match wildcard rules); JSON/CSV files may set
  `ion_radii`. One context mechanism (`crevice.radii.use_radii`) carries the
  set to every nested analysis and, explicitly, to `cavity-trajectory` worker
  processes.
- **Every output records its radius set**, default runs included: `radii`
  (name, source, file SHA-256, citation, table hash, overrides, ion radii) in
  result metadata, manifests and JSON outputs, the settings of
  `cavity-trajectory` statistics and trajectory hydration, and beside every
  CSV (`residues` without `--json` writes `<stem>_provenance.json`).
  `crevice.radii.radius_set_from_provenance`, `same_radius_set` and
  `reconcile_recorded_radii` rebuild and compare recorded sets.
- **Region definitions store their radius set** (schema version 2, `radii`).
  `region-prepare` writes it; `region-trajectory` and
  `hydration --cavity-results` measure with the recorded set; a different
  `--radii`, or `region-compare` of runs with different sets, is an error
  unless `--allow-radii-mismatch` is given (recorded as `radii_check`).
  Schema version 1 definitions still load and are read as the default set.
- **Per-frame enclosure-probe fallback in `trajectory`** (Python
  `analyze_trajectory(probe_fallback=...)`, CLI
  `--probe-fallback/--no-probe-fallback`). When the probe pinned on the
  reference frame does not resolve a frame, that frame is retried with the
  automatic choice made on it (same axis and origin). On by default for the
  automatic probe, off for an explicit probe. Each frame records
  `enclosure_probe` (source `pinned`/`fallback`/`none`, probe); the frames CSV
  gains `enclosure_probe_source`; `metadata.enclosure_probe.fallback` counts
  and lists fallback frames; profile distributions pool fallback frames and
  count them.
- Force-field ion names are recognised when a file has no element column,
  by residue and atom name together: CHARMM36 `SOD` (Na), `POT` (K), `CLA`
  (Cl), `CAL` (Ca), `MG`, `LIT`, `CES`, `ZN2`/`ZN`, `BAR`, `RUB`, `CD2`/`CD`;
  GROMACS/AMBER `NA`, `K`, `CL`, `CA`, `ZN`, `Na+`, `K+`, `Cl-` and others
  (`crevice.radii.ION_NAME_ELEMENTS`, sources in the inputs guide). Explicit
  element columns still win.
- Docs: new method page "Atomic radii"; ion-name table in the inputs guide;
  `--radii` in the tool pages.

- Automatic enclosure probe for `profile` and every command that builds a
  profile (`--enclosure-radius auto`, now the default; Python
  `enclosure_radius=None`). CREVICE resolves the channel at each probe of a
  0.5–3.0 Å ladder, groups the results into channels, and uses the smallest
  probe of the channel that resolves at the most probes (at least two). The
  choice, every probe tried and the reason are recorded in
  `metadata.enclosure_probe_selection` and the publish manifest
  (`profile_status.enclosure_probe`); if no channel resolves at two probes the
  profile is unresolved and the error lists the probes and reasons. An
  explicit number behaves exactly as before. `static-suite`, `benchmark` and
  `trajectory` use it too (see Changed); `cast` keeps its 0.8 Å rolling
  probe.
- `--lateral-exits` reports every distinct exit leg of a capped end
  (`metadata.exits[end].legs`, `leg_count`), for example both C2-related
  portals of TWIK-1, instead of breaking the tie. Exits are ranked by their
  width beyond the mouth sample's inscribed sphere, then by route length, and
  are distinct when their exit points are more than twice the bulk-probe
  radius apart and their routes, away from the start, stay more than one bulk
  radius apart. The first leg can therefore differ from the single leg
  reported before (the mouth no longer ties every route). All legs are drawn
  in the scenes.
- The two profile radius plots draw each lateral exit leg as a bluish-green
  dash-dot continuation beyond its mouth, against path length along the leg.

- `profile`/`publish`/`analyze` `--lateral-exits` (Python `lateral_exits=True`):
  opt-in support for capped channels. An end whose straight axial exit is
  blocked is accepted only if a widest free path to bulk (rolling-probe
  definition, `--exit-bulk-radius`, `--exit-spacing`) exists. The exit leg is
  reported separately (`metadata.exits`, `path_type`,
  `path_bottleneck_radius_A`, manifest `profile_status.path_type`) and drawn
  as a blue centre-line tube in the viewer scenes. Without the option all
  outputs are unchanged.

### Changed

- Channel-mouth guides are off by default: no orange rings in the PyMOL, VMD
  and ChimeraX scenes and no orange mouth bands or lines on profile plots.
  Mouth positions remain in the profile JSON, scene JSON and manifests.
- Exit legs are drawn as casts by default; the thin blue centre-line tubes need
  `--exit-centre-lines`.
- Non-channel casts (cavity, rolling-probe and dominant-region casts and the
  unresolved-profile fallback) are violet (`#a855f7`) in every viewer, and the
  width-profile lines of `cavity-trajectory` use the same violet; resolved
  channel casts stay teal. `scripts/check_pymol.py` and `scripts/check_vmd.py`
  accept either cast colour.

- Atoms that get no radius from the radius set in force (no
  `RESNAME:ATOMNAME`, atom-name, ion, HOLE or element entry) are now excluded
  when a PDB or mmCIF file is loaded, with an `UnrecognisedElementWarning` that
  explains how to supply radii (`--radii FILE`, `radii=`), and are recorded in
  the input report per element and per `RESNAME:ATOMNAME`. Previously such
  atoms (for example Gd) were kept with the 1.70 Å default. MD selections
  containing them now stop with the same explanation and a selection hint;
  hydration environments and obstacle selections drop them with a warning.
  `load_structure`/`load_structure_report` take `radii=`. Structures whose
  atoms all have radii, including every bundled benchmark, are unaffected.
- Deuterium (`D`) is treated as hydrogen: the hydrogen radius (unless a set
  lists `D`), removed by the default hydrogen filter and by the alignment fit
  selection. Previously it was a 1.70 Å heavy atom.
- **Recognised ion residues get CHARMM36 radii by default.** The default
  radius set (`crevice.radii.DEFAULT_RADII`, preset `default`) is the standard
  element table plus, for recognised monatomic ion residues only (residue
  *and* atom name in the ion table, e.g. `SOD`/`SOD`, `CLA`/`CLA`, `ZN2`/`ZN`,
  GROMACS or crystallographic `NA`/`NA`, `ZN`/`ZN`, `CA`/`CA`), the CHARMM36
  Lennard-Jones Rmin/2 from `toppar_water_ions.str` (`crevice.radii.ION_RADII`:
  Li 1.2975, Na 1.41075, Mg 1.185, K 1.76375, Ca 1.367, Rb 1.90, Cs 2.10, Ba
  1.89, Zn 1.09, Cd 1.357, Cl 2.27 Å). Every other atom, including C-alpha
  `CA` and metals inside other residues, keeps its element radius; Cu, F, Br, I
  ions keep theirs. Results change only where ion residues are obstacles or
  SASA environment atoms (for example GLUT1 region trajectories with
  `obstacle_selection: not resname TIP3`). The `bondi` preset now means the
  element table for every atom (ions as neutral atoms, the previous default);
  `charmm_like` gains the ion radii. Default outputs now record `radii`, so
  default manifests and JSON files are no longer byte-identical to earlier
  versions; numerical data are unchanged for structures without ion residues.

- Network chord diagrams (`network --chord-png`, `publish`, `static-suite`)
  put their edge-class legend in its own band below the whole circle,
  residue names included, so it can no longer overlap a residue name (it
  overlapped the lowest names of larger networks such as 4PYP). The figure is
  slightly taller; the circle, nodes, edges and labels are unchanged.
- `profile` always writes its per-sample table, as `<stem>.csv` beside the
  `-o` JSON (like `cavities`, `tunnels` and `network`); `--csv PATH` now only
  chooses a different path. Before, the CSV was written only with `--csv`.
- **`static-suite`, `benchmark` and `trajectory` choose the enclosure probe
  automatically by default** (`--enclosure-radius auto`; Python
  `enclosure_radius=None` for `run_static_benchmark_suite`); before, they
  used a fixed 0.8 Å. A number, for example `--enclosure-radius 0.8`, keeps
  the previous behaviour. `static-suite` and `benchmark` choose per system and
  record `enclosure_radius_A` and `enclosure_probe_mode` in their CSVs
  (`benchmark` also `enclosure_probe_reason`; `static-suite` system summaries
  the full `enclosure_probe` record). `trajectory` chooses **once**, on the
  reference frame, and passes that probe explicitly for every frame (no
  per-frame re-selection), as `analyze_trajectory` already did; the choice and
  reason are printed and recorded in `metadata.enclosure_probe` of the output
  JSON (also new in the Python result) and in a new `enclosure_radius_A`
  column of `<stem>_frames.csv`. `cavity-trajectory` and `region-trajectory`
  build no channel profile and are unchanged. `cast --probe-radius` stays
  0.8 Å: it defines the cast volume itself rather than the channel wall, so
  it is not chosen automatically (explained in the `cast` help and page).
- `static-suite` sensitivity rows keep each system's enclosure probe
  (`enclosure_radius_A` column) and a measurement probe that does not pass
  the channel now gives a `status: unresolved` row with its reason instead of
  failing the whole system. With the default probes (0, 1.0, 1.4 Å) the 1.4 Å
  row of 1GRM is unresolved (narrowest radius about 1.33 Å), and 1GRM now
  completes.
- **The `hole` radius preset is now HOLE's own `simple.rad`** (HOLE 2.3.1;
  AMBER united-atom radii of Weiner et al. 1984: C 1.85, O 1.65, S 2.00,
  N 1.75, H 1.00, P 2.10 Å, matched by atom name, unmatched atoms an error).
  It was previously the standard table with only hydrogen set to 1.00 Å.
- **Element inference for ion names (GRO and other files without elements).**
  CHARMM `SOD` was read as S (1.80 Å) and is now Na (radius: see the CHARMM36
  ion radii above); `POT` P → K,
  `CAL` C → Ca, `CES` C → Cs, `BAR` B → Ba, `LIT`/`RUB` unknown → Li/Rb,
  `CD` in residue `CD2` C → Cd, and GROMACS `RB`/`CS` are now Rb/Cs. `CLA`,
  `NA`, `K`, `MG`, `ZN` and `CA` ions were already correct. Analyses that
  include such ions as obstacles or SASA environment atoms (for example GLUT1
  region trajectories with `obstacle_selection: not resname TIP3`) can change.
- `crevice.radii.custom_radii_set` validates keys and radii and accepts
  `base=`; `atom_vdw_radius` takes an optional `radii` argument.

- **Outputs are CSV first; no HTML, no cast preview, no tunnel figure.**
  Every tabular dataset is now written as CSV with units in the column names
  (`_A` = Å, `_A2` = Å², `_A3` = Å³, `_ps`, `_degrees`). New tables:
  `*_void_cast.csv`, `*_cavities.csv`, `*_tunnels.csv`,
  `*_tunnel_points.csv`, `*_network_nodes.csv`, `*_network_edges.csv`,
  `*_connectivity.csv`/`_groups.csv`, `*_features.csv`, trajectory
  `<stem>_frames.csv`/`<stem>_profiles.csv`, `*_hydration_frames.csv`,
  `*_hydration_region_frames.csv`, `*_residue_conformation_frames.csv`,
  `*_interaction_edge_frames.csv`, `*_cavity_section_frames.csv`,
  `*_cavity_residue_frames.csv`, `*_residue_evidence_partners.csv`,
  `--network-json` `<stem>_residues.csv`, `*_review_landmarks.csv`,
  `paired_observations.csv` and `benchmark` `<stem>.csv`. The single-analysis
  commands (`cavities`, `tunnels`, `network`, `features`, `trajectory`) write
  their CSV beside the `-o` JSON. JSON remains for manifests, provenance, full
  structured records read by other commands and viewer scene definitions (see
  "Data files" in the outputs guide). Numbers are unchanged (1GRM publish:
  minimum radius 1.3301888073 Å, cast 292.40625 Å³).
- Removed from default outputs: `*_analysis.html`, region-review
  `index.html` and `*_analysis_report_data.json` (their values are in the
  CSV/JSON files); `*_cast_preview.png` (`cast` and the `publish` fallback;
  the viewer scenes are the quality view of a cast); `*_tunnel_summary.png`
  (`publish`, `static-suite`, `analyze --png`); and duplicate JSON copies of
  CSV tables (`publish`/`static-suite` `*_residue_contacts.json`, `analyze`
  `*_residues.json` and `*_features.json`, `*_hydration_mobility.json`,
  `*_profile_sensitivity.json`).
- Figures: the residue-contact figure has a legend naming the contact roles
  its colours encode (bottleneck, bottleneck-nearby, lining, nearby, with
  their gap definitions) and an explained influence-score axis. The network
  summary now ranks residues by residue-residue contacts, stacked by contact
  class with a legend, on an axis that states the rule, for example "Residue
  contacts (degree: residues with atom-centre distance ≤ 4.5 Å)"; the region
  pseudo-node is no longer a bar. Network chord edges are now coloured by the
  interaction names the network actually writes (previously every residue
  edge fell through to grey) and the legend lists only drawn classes.
- JSON unit text uses Å (`water oxygens per Å³`; membership `units`), and the
  CLI help strings of the edited options use Å.

- Default enclosure probe: automatic instead of a fixed 0.8 Å (see Added).
  Profiles that resolved at 0.8 Å can change where the automatic rule picks a
  different probe: the 1GRM quick-start command (`crevice profile 1GRM.pdb`,
  81 samples, 0.5 Å sections) now uses 0.75 Å and reports a minimum of 1.324 Å
  instead of 1.328 Å; the canonical 1GRM settings (0.25 Å sections, 101
  samples) still choose 0.8 Å and are unchanged. Pass `--enclosure-radius 0.8`
  to reproduce earlier results.
- `metadata.exits[end]` of a lateral (capped) end now holds a `legs` list;
  per-leg fields (points, radii, bottleneck, exit point, length, vertices)
  moved into the legs. Scene JSON `exit_paths` entries gain a `leg` rank.
- The profile radius plots no longer draw a marker at every sample (line and
  fill only; the bottleneck dot stays). The grey exterior-clearance line is
  not drawn at a capped end.
- Section discovery (`SectionGeometry.sections`) is vectorised with NumPy
  using the same arithmetic; results are identical and large search grids are
  10–20 times faster.

  **Migration.** Rename columns when reading older CSVs: profile `x,y,z,
  radius, raw_clearance, axis_position` → `x_A, y_A, z_A, radius_A,
  raw_clearance_A, axis_position_A`; residue contacts `min_distance,
  mean_distance` → `min_distance_A, mean_distance_A`; `--distribution-csv`
  `coordinate, mean_radius, ...` → `coordinate_A, mean_radius_A, ...`;
  `static-suite` summary and sensitivity columns gain `_A`/`_A3`. Read
  residue contacts from `*_residue_contacts.csv` instead of the JSON copy.
  `tunnels --png`/`--dpi` and `--annotate` on `tunnels` were removed, and
  `crevice.plot_tunnel_summary` and `crevice.plot_void_cast` no longer exist;
  plot `*_tunnels.csv` or open the viewer scenes. `region-prepare` prints the
  review scene and landmark table instead of `index.html`, and its result has
  `review_landmarks_csv` instead of `review_html`.
- A profile whose resampling loses continuity now names the first axial
  position that could not be reached and mentions wall openings wider than the
  enclosure probe.


- Tunnels: each path point's `t` (JSON `axis_position`) is now the cumulative
  path length in Å from the tunnel start, like explicit-path profiles, instead
  of the step number; the step number remains the point `index`. This applies
  to `crevice tunnels`, the tunnels in `analyze`, `publish`, `static-suite` and
  `trajectory`, and `crevice.find_tunnels`/`tunnel_profile`. Every other tunnel
  field is unchanged (checked on 1GRM: identical apart from this field; PDB
  output byte-identical). `profile_distribution` still pools only axial
  profiles and now says why it refuses a path profile.
  (`crevice.analysis.find_tunnels`, the underlying grid search, still returns
  step numbers.)
- Python benchmark helpers (`benchmark_path`, `fetch_benchmark_structure`,
  `fetch_benchmark_record`, `fetch_all_benchmark_structures`) now default to
  mmCIF (`extension="cif"`), like `fetch_structure` and the `--format` default
  of `fetch`, `benchmark` and `static-suite`. Pass `extension="pdb"` for a
  PDB-format cache.


### Fixed

- `crevice.analysis.find_tunnels` (the grid engine behind
  `crevice.find_tunnels`) now also reports each tunnel point's `t` as the
  cumulative path length in Å rather than the step number.
- CLI help: `--samples` now says it is a minimum (connected profiles add
  points to keep axial steps ≤ 0.75 Å); profile options state Å units.

- `profile`: where a plane holds both the lumen and a small enclosed pocket in
  the wall that ties with it, the path now follows the widest connected route
  (largest bottleneck) instead of section discovery order. Previously the
  reported minimum radius could depend on the sample count and section
  spacing: on 1GRM some settings, including the quick-start defaults, gave
  0.84–0.87 Å from a one-point wall pocket instead of about 1.33 Å. Profiles
  that did not pass through such a pocket are unchanged.
- Element inference: an explicit element column that holds any real
  periodic-table symbol is now kept as written. Previously a two-letter element
  missing from the radius table was cut to its first letter, so mercury (`HG`)
  became hydrogen and was removed by the hydrogen filter, and `CD`/`GD` became
  carbon/`G`. Without an element column, ligand atoms such as `CAB` or heme `NA`
  are no longer read as calcium or sodium; `CA`/`NA`/`CO`/`NI` mean a metal
  only for a single-atom residue of the same name. Radii were added for further
  elements from Bondi (1964) (for example Hg, Cd, Ag, Au, Pt, Pd, Pb, Tl, Si,
  As, noble gases) and Mantina et al. (2009) (for example Al, Rb, Cs, Sr, Ba);
  existing radii are unchanged. Elements with no tabulated radius keep their
  symbol, use 1.70 Å, and are listed in the input report
  (`unsupported_radius_elements`, `radius_note`) with a warning. Results change
  only for structures containing these elements or, without element columns,
  such atom names.
- Accession identity check: a fetch of an extended accession such as
  `pdb_00001grm` was compared with the file's `_entry.id` after upper-casing
  only, so it failed against the classic `1GRM` that RCSB files state (a strict
  download raised an error). Classic and extended forms of the same entry now
  match, and a cache entry without provenance records the extended ID in its
  normalised lower-case form.
- `crevice trajectory --inspect` no longer requires `-o/--output`; it prints to
  the terminal. Analysis runs without `-o` now stop with a clear error.
- Static hydration (used by `analyze`, `cast`, `publish`, `residue-evidence`):
  a focus residue that is present but has no non-HETATM heavy atom stopped the
  run with "Hydration focus contains unknown protein residue identities" (for
  example `crevice analyze 1GRM.pdb --axis auto`, whose D-amino acids and
  modified ends are HETATM in the PDB file). Such residues are now listed in the
  hydration JSON as `settings.untyped_focus_residues` with a reason and are not
  measured; identities absent from the structure are still an error.
- Figures: length units in the profile, landmark, cavity-summary and
  tunnel-summary plots (axis labels, value call-outs, mouth labels and the PNG
  `Description` metadata) read "Å"/"Å³" instead of "A"/"A^3".
- Built-in structure readers: a residue deposited only with alternate
  locations other than `A`/`1` (for example only `B` and `C`) was dropped
  entirely. PDB and built-in mmCIF input now keep blank locations plus the
  first alternate location of each residue, the rule the Gemmi reader already
  used. The bundled PDB and cached mmCIF inputs load identically.
- Built-in mmCIF reader: `chain_ids`/`--chain-ids` was ignored and author chain
  IDs were always used. It now honours the option (default `label`, as for
  the Gemmi reader), falling back per row to the other column when one is
  null. With `--parser builtin`, ligands can therefore get their own chain ID
  (for example BNG in 4PYP moves from chain A to B, matching Gemmi). Reports
  now include `chain_id_namespace` and `zero_occupancy_policy`.
- Rolling-probe cast metadata: `enclosure_definition` described the
  dominant-mode use of `enclosure_fraction` even for `selection_mode="all"`,
  where it is a hard crop of the probe centres. The text now states the rule
  that was applied for the mode used (or that the test was not applied when the
  fraction is 0), and `enclosure_required_rays` is recorded. Cast geometry is
  unchanged.

## [0.1.0]

Initial public release.

### Added

- Structure intake: PDB and mmCIF with model and chain selection; the
  entity-aware Gemmi mmCIF loader when `crevice[structures]` is installed, with
  a built-in parser as fallback; GRO snapshots and MD trajectories through
  MDAnalysis, with atom-selection, frame-range and periodic-boundary reports.
- RCSB download by accession (`fetch`, or any structure argument), with
  format validation, atomic caching, SHA-256 provenance sidecars and
  cache verification.
- Geometry: connected HOLE-like channel radius profiles with continuous
  atom-clearance checks (`profile`); HOLLOW-like void casts and rolling-probe
  cavity casts with dominant-region selection and optional seeding (`cast`);
  enclosed-cavity detection (`cavities`); CAVER-like widest-path tunnel search
  (`tunnels`).
- Residues: profile-lining annotation (`residues`), channel/cavity-aware
  residue interaction networks (`network`), and measured-boundary residue
  attribution with nonlocal contact partners (`residue-evidence`), including
  `--stick-residues all|top:N|IDS|@FILE` to choose the residues drawn as
  sticks.
- Trajectories: aligned through-channel profile distributions with
  batch-means mean intervals (`trajectory`); all-frame cavity volumes,
  sectional widths, changing boundary residues, contacts and motion
  (`cavity-trajectory`); physical-obstacle selections and lattice-phase
  controls.
- Named regions: versioned, hashed region definitions with landmark proposals,
  review overlays and paired comparisons (`region-init`, `region-prepare`,
  `region-trajectory`, `region-compare`).
- Hydration: explicit residue-water contacts, geometric solvent-accessible
  area, occupancy, sampled persistence, aligned water density maps, typed
  interaction candidates and conformation histories (`hydration`, and standard
  in the workflow commands).
- Instantaneous-cavity water membership (`--water-membership` on
  `cavity-trajectory` and `region-trajectory`), paired per frame with the
  fixed-reference water count.
- `publish`: a complete static-structure bundle; when no unique through-profile
  resolves it writes a rolling-probe cast, records the unresolved status and
  still completes the requested cavity, tunnel and network analyses.
- `--smooth X` constrained display smoothing (default 0.4 Å) for `cast` and
  `publish` surfaces; measured grids, volumes and profiles are unaffected.
- Exports: JSON, CSV, PDB and OpenDX data, publication figures, and PyMOL, VMD
  and ChimeraX scene scripts sharing one viewer structure.
- Documentation site (Sphinx) with a page per tool and method notes.

### Known limitations

- Benchmark entries are not curated; axes and regions are chosen by the
  analysis or the user.
- Channel profiles require a passage that is monotone along a proposed
  direction with two open ends; branched or folded paths are not profiled.
- Trajectory mean intervals are approximate and assume stationary,
  representative sampling.
- The networks written by `publish` and `analyze` use the low-level network
  defaults (surface gap, heteroatoms excluded), which differ from
  `crevice network`.
- `region-compare` does not compare instantaneous-cavity water membership.

# Outputs and unresolved results

This page explains what CREVICE writes, how to read the tables, and what it
means when an analysis comes back unresolved.

## What a run writes

The low-level commands (`profile`, `residues`, `cavities`, `cavities --tunnels`, `residues --network`,
`residues --features`) write the main result file you name with `-o`. `profile`,
`cavities`, `cavities --tunnels`, `residues --network` and `residues --features` also write their tables as CSV
beside it, with the same stem (for example `-o 1GRM_profile.json` gives
`1GRM_profile.csv`, and `-o net.json` gives `net_nodes.csv` and
`net_edges.csv`); `profile --csv PATH` chooses another path for the profile
table, and `residues -o` is itself a CSV. Other files (`--pdb`, `--png`, ...)
are written only when you ask for them.

The workflow commands (`cavities --cast`, `report`, `residues --evidence`, `trajectory --water`,
`trajectory`, `trajectory --cavities`, `trajectory --region`, `report --analyze`) write a
directory (`--out-dir`) with a `*_manifest.json` that lists every file and
records the run: CREVICE version, command line, every effective setting with
where it came from (default, `--config` file or command line), units, the
SHA-256 of each input file and the output plan (under `run`).

(output-policy)=
## Choosing what is written (`--outputs`, `--all`)

Since 0.3 these commands write a **compact default set**, and everything else
is opt-in. The same rules apply on the command line and in Python
(`outputs=` on {func}`crevice.write_static_publication_bundle`,
{func}`crevice.figures.write_rolling_probe_bundle` and
{func}`crevice.static_workflow.export_static_workflow`, or
{func}`crevice.outputs.use_outputs`).

By default a bundle command writes:

- **data**, always: the manifest, the JSON records, every CSV table, NPZ arrays
  and the measured DX maps that later commands read (`*_volume.dx`);
- **key figures**: the figure that sums up the run (for example
  `*_profile_radius.png` and `*_residue_contacts.png` for `report`,
  `*_residue_evidence.png`, `*_hydration_residues.png`, or the cavity volume,
  radius and residue figures of the trajectory commands);
- **one scene set per native viewer**: PyMOL (`.pml`), VMD (`.vmd`/`.tcl`) and
  ChimeraX (`.cxc`), with the scene data they load;
- the **typed residue–residue interactions** of the lining or boundary residues
  (hydrogen bonds, salt-bridge candidates, aromatic contacts) and residue
  conformations, as before.

The **water analyses** (water contacts and occupancy, solvent-accessible area,
water bridges, water density and the hydration scenes) no longer run by default
in `cavities --cast`, `report`, `report --analyze`, `residues --evidence`, `trajectory`,
`trajectory --cavities` and `trajectory --region`. Request them with
`--outputs default,hydration` or `--all`. They still run by default where
water is the point of the command: `trajectory --water`, and `trajectory --cavities` /
`trajectory --region` with `--water-membership`. `report` also no longer runs
the separate cavity search, tunnel search and residue network unless asked
(`--outputs default,cavities,tunnels,network`). Residue contacts, lining
residues and cast segments are always written.

### What each command writes

| Command | Written by default | Added by `--all` |
|---|---|---|
| `report` | data; `*_profile_radius.png` and `*_residue_contacts.png`; one PyMOL, VMD and ChimeraX scene of the cast and its entry, lumen and exit segments; lining residues; typed residue interactions | cavity search, tunnel search and residue network (tables and figures); hydration; pore fingerprint; the annotated radius figure; input copy and dummy-atom PDBs; `*_pore_cast.*` and `*_pore_spheres.pml` scenes; HTML view |
| `cavities --cast` | data, including one map per region in `regions/`; PyMOL, VMD and ChimeraX scenes of the cast; typed residue interactions | hydration (with its scenes and figures); input copy; `*_pore_cast.*` and `*_end_view.*` scenes; HTML view |
| `report --analyze` | data only: the CSV and JSON of each `--analysis`, and typed residue interactions | figures (as `--png`); hydration with its scenes; HTML view |
| `residues --evidence` | data; `*_residue_evidence.png`; PyMOL, VMD and ChimeraX `*_residue_context.*` scenes; typed residue interactions | hydration (with its scenes and figures); HTML view |
| `trajectory --water` | data and water-density maps; `*_hydration_residues.png`; PyMOL, VMD and ChimeraX hydration scenes; typed residue interactions | the other figures (overview, legend, density projections, trajectory heatmaps and time series); HTML view |
| `trajectory --cavities`, `trajectory --region` | data; `*_cavity_volume_timeseries.png`, `*_cavity_radius_profile.png` and `*_cavity_residues.png`; typed residue interactions (with `--water-membership`, also hydration, its scenes and the membership files) | hydration and its scenes; every other figure (width and residue heatmaps, partner, section and free-sphere figures) |
| `trajectory` | data (`-o` JSON, `*_frames.csv`, `*_profiles.csv`, per-frame residue interactions) and any file you name (`--png`, `--distribution-csv`, `--network-json`) | per-frame hydration bundles with their figures and scenes |

"Data" always includes the manifest, the CSV tables, the JSON records and the
measured maps. The low-level commands (`profile`, `residues`, `cavities`,
`cavities --tunnels`, `residues --network`, `residues --features`, `profile --fingerprint`, `paths`, `trajectory --track`,
`view`) don't take `--outputs`; they write the files described on their tool
pages.

`--outputs LIST` takes comma-separated values (`--outputs help` lists them with
the command's defaults):

| Value | What it adds |
|---|---|
| `key-figures` | the command's key figures |
| `figures` | every Matplotlib figure |
| `pymol`, `vmd`, `chimerax` | that viewer's scenes; `viewers` = all three |
| `scenes` | scene data only (viewer PDB, display maps, meshes, scene JSON); added automatically with any viewer |
| `viewer-html` | a self-contained interactive HTML view (`*_view.html`, see {doc}`view <../tools/view>`) |
| `structures` | a copy of the input structure and dummy-atom PDBs |
| `extra-views` | alias and variant scenes: `*_pore_cast.*`, `*_end_view.*`, `*_pore_spheres.pml` |
| `hydration` | the water analyses (implies `interactions`) |
| `interactions` | typed residue–residue interactions and conformations |
| `cavities`, `tunnels`, `network`, `fingerprint` | the extra `report` analyses |
| `default` | the command's default selection |
| `none` | data only |

Listing file kinds (`figures`, `pymol`, ...) replaces the default file kinds;
listing analyses (`hydration`, `cavities`, ...) replaces the default optional
analyses; the other half keeps its default. So `--outputs figures` writes the
data and every figure but no scene, `--outputs pymol` the data and the PyMOL
scene only, `--outputs hydration` adds the water analyses to the default
files, and `--outputs default,vmd,hydration` is the default plus hydration.

`--all` writes every file kind and optional analysis of the command. It never
raises an allocation or frame guard; an output whose prerequisite is missing
is skipped and the reason recorded under `run.outputs.notes` in the manifest
(for example, no fingerprint when no through-channel resolves, never an
invented profile). The `--skip-*` flags still remove an analysis, also from
`--all`:

| Flag | Omits |
|---|---|
| `--skip-hydration` | the water analyses |
| `--skip-interactions` | the typed residue interactions (and, with hydration, the water bridges) |
| `--skip-water-density` | 3D water-density maps (hydration) |
| `--skip-analysis-views` | the hydration scenes, evidence table and overview figures (the numbers are still written) |
| `report --skip-cavities / --skip-tunnels / --skip-network` | the extra grid searches and network |

The plan is printed on standard error before any work when standard error is a
terminal (or with `CREVICE_PROGRESS=1`). `--dry-run` validates the settings,
prints the plan and exits without reading an input or writing a file.

`--legacy-outputs`, which reproduced the CREVICE 0.2 file set during the 0.3
series, was removed in 0.4; `--all` writes every file kind and optional
analysis. See {doc}`Upgrading to 0.4 <upgrading-0.4>`.

Files written by the synthetic test cases (`tests/inventory_cases.py`) and by
`report 1GRM|4PYP --parser gemmi --samples 101` with `--section-spacing` and
`--cast-spacing` 0.25 Å (1GRM) or 0.5 Å (4PYP), `--hydration-sasa-points 512`
and `--dpi 400`:

| Command | 0.2 default | 0.3 default | 0.3 `--all` |
|---|---|---|---|
| `report` (1GRM) | 96 | 40 | |
| `report` (4PYP, no channel resolves) | 86 | 29 | |
| `report` (synthetic channel) | 90 | 34 | 94 |
| `cavities --cast` (synthetic pocket) | 72 | 28 | 73 |
| `report --analyze` (profile, residues) | 42 | 10 | 45 |
| `residues --evidence` | 72 | 32 | 73 |
| `trajectory --water` (static / trajectory) | 43 / 53 | 40 / 48 | 44 / 54 |
| `trajectory --cavities` | 86 | 27 | 87 |
| `trajectory --region` | 80 | 28 | 80 |
| `trajectory` (2 static frames) | 79 | 15 | 79 |

The numbers are the same in every mode: the selection changes which files are
written, never a measurement.

(settings-files)=
## Settings files (`--config`, `--print-config`)

`--print-config` prints a command's effective settings as JSON, grouped into
sections (`input`, `geometry`, `radii`, `outputs`, `display`, `trajectory`,
`hydration`, `method`), and exits. Save it, edit it and pass it back with
`--config FILE`: flags on the command line override the file, the file
overrides the built-in defaults, and positional inputs (the structure) always
come from the command line. Settings files are JSON. The same settings are
available in Python as {class}`crevice.config.RunSettings`
({func}`crevice.config.settings` builds and validates them,
{func}`crevice.config.run` runs them exactly like the command line).
Every setting is checked (units, ranges, input files, combinations) before any
structure is read.

## Data files: CSV first

Every table you would analyse is a CSV: profiles, residue contacts, cast
regions, cavities, tunnels and tunnel path nodes, network nodes and edges,
connectivity paths, features, trajectory frames and profiles, hydration per
residue and per frame, typed interactions and conformations per frame, cavity
trajectories per frame, per section and per boundary residue, water membership,
residue evidence and partners, and region comparisons. The per-quantity summary
statistics of the trajectory runs (mean, SD, quantiles, effective sample size
and the interval for the mean of every volume, residue, contact and water
series) are flattened into one table per run: `*_cavity_statistics_summary.csv`,
`*_hydration_summary.csv`, `*_water_membership_summary.csv` (plus
`*_water_membership_strata.csv`) and `comparison_summary.csv`, one row per
quantity with a `quantity` (its JSON path) and `item` (residue or pair) column.
Column names are ASCII
and carry their unit: `_A` = Å, `_A2` = Å², `_A3` = Å³, `_ps` = picoseconds,
`_degrees`; counts, fractions and scores have no suffix. An empty cell is an
unresolved or unobserved value (NaN or `null` elsewhere), never zero. The tool
pages list the columns of each file.

JSON is kept only where a nested, machine-readable record is needed:

| JSON | Why it is kept |
|---|---|
| `*_manifest.json`, `*_hydration_manifest.json`, `*_cavity_manifest.json`, `*_evidence_manifest.json`, `manifest.json` | file index, settings and status of a run |
| `*_input_report.json`, `*_reader.json`, `reader.json`, `request.json`, `*_region_provenance.json`, `*_selection.json`, `*_hydration_alignment.json`, `*_cavity_frame_diagnostics.json` | input, reader, atom-selection and alignment provenance |
| `*_profile.json`, `*_void_cast.json`, `*_cavities.json`, `*_tunnels.json`, `*_network.json`, `*_connectivity.json`, `*_residue_evidence.json`, `*_interactions.json`, `*_hydration.json`, `*_water_density.json`, `*_water_membership.json`, `*_cavity_statistics.json`, `*_cavity_profile_statistics.json`, `comparison.json` | full structured records (settings, metadata, nested statistics and intervals) that other commands read, for example `hydration --cavity-results`, `region compare` and `trajectory --lining-evidence` |
| `*_scene.json`, `*_analysis_scene.json`, `*_review_scene.json`, `*_volume_metadata.json`, `*_region_annotations.json`, `*.identities.json` | scene definitions and identity maps read by the viewer scripts |
| `*_region.json`, `*_review.json` | versioned region definitions and their review diagnostics |
| `-o` results of `profile`, `cavities`, `cavities --tunnels`, `residues --network`, `residues --features`, `trajectory`, and `--network-json`, `--report-json`, `fetch --json`, `report --benchmark --output` | the file you named; its tables are also written as CSV |

NPZ and DX files hold grids and arrays (measured and display maps, density,
per-frame arrays) that viewers and later commands load; their tabular content is
also in the CSV files.

(recorded-paths)=
## Recorded paths and moving an output folder

From 0.4.0, every file path that an output records (manifest file indexes,
run records, scene and metadata JSON, region and comparison records, suite and
benchmark summaries) is written **relative to the directory of the file that
records it**, with `/` separators on every operating system. A manifest's
`"files": {"volume_dx": "p_volume.dx", ...}` means `p_volume.dx` next to the
manifest. An output folder can therefore be moved, copied to another machine
or opened on another operating system as a whole.

Some paths stay absolute (with `/` separators, for example `C:/data/prot.pdb`
on Windows), because they point outside the folder and cannot move with it:
input files elsewhere on disk (the run record keeps each input's SHA-256 next
to it, so a copied input can be checked), a trajectory on another drive, and
the reference map of a region prepared in another folder. The command line
(`run.argv`), the input arguments (`run.inputs[].argument`) and the settings
(`run.settings`) are recorded as typed, so that `--config` reads them back.

`crevice view`, `region-compare`, `hydration --cavity-results` and the Python
API resolve a relative path against the recording file's directory. They also
read 0.3.x folders, whose manifests hold absolute paths of the machine that
wrote them: a path that no longer exists is looked for inside the folder.
`crevice.config.run` returns the file index as absolute paths. If you read a
manifest yourself, join each path with the manifest's directory
(`Path(manifest).parent / value`; an absolute value is used as it is).

The viewer scripts (PyMOL, VMD, ChimeraX) contain no absolute path of the
machine that wrote them. Each script names its data files relative to its own
folder and looks for them next to itself, then in the working directory; if a
file is in neither place the script stops with an error naming that file. Keep
a script together with the other files of its folder (moving or copying the
whole folder is fine). The one exception to "next to the script" is
`profile --pymol` on a structure stored elsewhere: the script names it by a
relative path such as `../inputs/prot.pdb` (by its file name alone when the
two share no folder below the root), so keep the two folders in the same
relative position or copy the structure next to the script.

VMD's `-e` option does not tell a script where it is, so with `vmd -e` the
script finds its files only in the working directory: start VMD in the output
folder (`cd out; vmd -e 1GRM_volume.vmd`), or type `source /path/to/1GRM_volume.vmd`
in VMD's console, which works from any folder. PyMOL and ChimeraX know the
script's location and work from any folder.

No command writes HTML unless asked: every value is in the CSV and JSON files
above, and `--outputs viewer-html` (or `--all`) adds the interactive view.

`residues --evidence --stick-residues` controls which residues its viewer
overlays draw as sticks (default `top:12`; see
[`residues --evidence`](../tools/residue-evidence.md#choosing-which-residues-are-drawn-as-sticks)).

## Unresolved, missing and zero are different

CREVICE keeps three outcomes apart. They should not be merged when you summarise
results:

- **Unresolved.** The method could not produce a defensible result for this
  input. For example, `profile --axis auto` found no unique two-ended channel,
  or a trajectory frame matched no reference cavity, or matched more than one.
  The status and reason are recorded. The value is `null`, NaN or absent. It is
  never 0.
- **Unavailable / not observed.** The input does not contain what the analysis
  needs. For example, a crystal structure with no explicit water gives water
  counts with status `unavailable_no_explicit_water`, and SASA is still
  reported. A missing optional dependency is also recorded, and other analyses
  continue.
- **Measured zero.** The analysis ran and found nothing in the declared region,
  such as zero retained cavity samples or zero waters in contact.

None of these proves a biological state on its own. An unresolved through-channel
does not show that the protein is closed. A zero volume in a declared region does
not show a closed gate. Water absent from a static structure is not evidence of
a dry pocket.

### The `report` fallback

When `report` cannot resolve a unique through-profile, it warns, records the
reason and writes a 3D rolling-probe cast. It does not make up a profile. The
cavity, tunnel and network analyses you requested (`--outputs` or `--all`)
still run.
`<prefix>_profile_status.json` and the manifest's `profile_status` and
`fallback_analyses` entries record what was completed, what was skipped on
request and what could not be completed. See
[`report`](../tools/publish.md#when-no-channel-resolves).

## Figure and scene text (`--annotate`)

CREVICE's figures and viewer scenes are **clean by default**. They draw the
data with axis labels (with units), tick labels, colour bars for colour scales
that encode a value, and a short legend only where one panel shows several
series distinguished by colour or line style. They do not draw titles,
figure titles, captions, value call-outs (such as the bottleneck radius),
channel-mouth labels, residue-landmark names, interpretation notes, viewer
colour keys, 2D labels or residue labels in 3D scenes. Two kinds of text are
always kept, because a figure cannot be read without them: category labels
(residue IDs on bar charts and heatmaps, node names in the network diagrams),
legends that name the categories a colour encodes (for example the contact
roles in the residue-contact figure and the contact classes in the network
figures) and a one-line notice in a panel that has no data at all.

Add `--annotate` to any command that writes figures or scenes (Python:
`annotate=True` on any figure or scene writer, or the
{func}`crevice.presentation.figure_annotations` context manager) to draw the
explanatory text as well. It never changes a measurement or data file.

Nothing is lost when the text is not drawn:

- every PNG carries `Title`, `Description` and `CREVICE annotations` metadata:
  what the figure shows, its key values and every omitted title, note, label
  and axis detail. Read it with {func}`crevice.presentation.figure_description`
  (or any PNG metadata viewer);
- scene definitions (`*_analysis_scene.json`, `*_review_scene.json`) record
  the scene text in an `annotations` entry.

## Display guides and cast colours

Channel-mouth guides are also off by default: viewer scenes draw no mouth rings
and radius plots no orange mouth bands or lines. `--mouth-guides` (Python:
{func}`crevice.presentation.display_guides` with `mouth_guides=True`, or
`mouth_guides=True` on {func}`crevice.volume_export.write_volume_viewer_bundle`)
draws them. With `--lateral-exits`, exit legs are drawn as EXIT segment casts;
`--exit-centre-lines` adds their thin centre-line tubes. Mouth positions and
leg centre lines are always written to the JSON outputs.

A resolved channel cast is teal (`#1494a1`). Every other cast (cavity,
rolling-probe and dominant-region casts, and the rolling-probe cast written
when no channel profile resolves) is violet (`#a855f7`,
{data}`crevice.presentation.NON_CHANNEL_CAST_RGB`), in every viewer, and the
width-profile lines of non-channel casts use the same violet. The violet was
checked to stay distinguishable, including under simulated red-green colour
vision deficiency, from the channel teal and profile blue, the exit-leg green
(`#009e73`), the orange member-water and interval colours and both ends of the
red-to-blue hydration scale.

A resolved channel cast is split into separate objects (see
{doc}`report <../tools/publish>`, "Cast segments"): `crevice_entry_N` yellow
(`#f0d43a`), `crevice_lumen` teal and `crevice_exit_N` rust (`#c2410c`), each
recolourable and hideable on its own. "Entry" and "exit" are geometric labels,
not a transport direction. Residues lining each segment (atom centre within
3.3 Å of the segment's grid points, `--lining-cutoff`) are written to
`*_lining_residues.csv`, flagged in the residue and network tables
(`lumen_lining`, `lining_segments`) and included in the scenes as stick objects
`crevice_<segment>_lining` that are hidden until switched on.

The colours and representations of every figure and scene are listed in the
{doc}`tool pages <../tools/index>` and in the API documentation of each writer.

## Viewer scripts

Cast, evidence and hydration bundles include scene scripts for PyMOL (`.pml`),
VMD (`.vmd`/`.tcl`) and ChimeraX (`.cxc`) by default (`--outputs pymol` etc.
choose fewer). Keep the whole output directory
together, because scripts refer to companion meshes, maps, PDB files and
identity sidecars by relative path (see {ref}`recorded-paths`). `crevice_render` (PyMOL/VMD) and `*_render.cxc`
(ChimeraX) make high-resolution images.

Always open a scene and have a look before you use it in a figure: a script
that CREVICE wrote has not been checked in the viewer on your machine. The
source repository also has `scripts/check_pymol.py`, `check_vmd.py` and
`check_chimerax.py`, which run a scene headlessly and check its coordinates,
colours and surfaces.

## Measured maps versus display surfaces

Cast bundles contain a measured binary map (`*_volume.dx`), where each grid
sample is either occupied or empty, and separate smoothed display fields
(`*_display.dx`, PyMOL meshes). Volumes, radii and residue attribution come from
the binary map only. Display smoothing changes how the surface looks, not what
was measured. Give the binary map, not a display map, to commands that take
`--volume-dx` or `--reference-volume-dx`.

`cavities --cast` and `report` take `--smooth X` (X in Å, default 0.4; `0` shows the raw
voxel boundary; Python `smooth=X`). It applies constrained Gaussian twicing on
a supersampled display grid, so every measured sample keeps its inside/outside
state and the surface moves only within boundary grid cells. The simpler
single-Gaussian `--surface-smoothing W` is an alternative and cannot be combined
with `--smooth`. The chosen method and
width are recorded under `display_smoothing` in `*_volume_metadata.json`. See
[Display-surface smoothing](../SURFACE_SMOOTHING.md).

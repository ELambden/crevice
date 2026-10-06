# Outputs and unresolved results

## What a run writes

The low-level commands (`profile`, `residues`, `cavities`, `tunnels`, `network`,
`features`) write the main result file you name with `-o`. `profile`,
`cavities`, `tunnels`, `network` and `features` also write their tables as CSV
beside it, with the same stem (for example `-o 1GRM_profile.json` gives
`1GRM_profile.csv`, and `-o net.json` gives `net_nodes.csv` and
`net_edges.csv`); `profile --csv PATH` chooses another path for the profile
table, and `residues -o` is itself a CSV. Other files (`--pdb`, `--png`, ...)
are written only when you ask for them.

The workflow commands (`cast`, `publish`, `residue-evidence`, `hydration`,
`trajectory`, `cavity-trajectory`, `region-trajectory`, `analyze`) write a
directory (`--out-dir`) with a `*_manifest.json` that lists every file. Most of
them also add standard hydration, typed-interaction and viewer outputs. Opt-out
flags control those extras:

| Flag | Omits |
|---|---|
| `--skip-hydration` | water contacts and solvent-accessible area |
| `--skip-water-density` | 3D water-density maps |
| `--skip-interactions` | typed chemical interactions and conformation histories |
| `--skip-analysis-views` | combined viewer scenes, the evidence table and the overview figures (the numbers are still written) |
| `publish --skip-cavities / --skip-tunnels / --skip-network` | the extra grid searches and network |

A default `publish` of 1GRM writes about 85 files, and a default `cast` of 1GRM
about 70 (both include the standard hydration outputs).

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
| `*_profile.json`, `*_void_cast.json`, `*_cavities.json`, `*_tunnels.json`, `*_network.json`, `*_connectivity.json`, `*_residue_evidence.json`, `*_interactions.json`, `*_hydration.json`, `*_water_density.json`, `*_water_membership.json`, `*_cavity_statistics.json`, `*_cavity_profile_statistics.json`, `comparison.json` | full structured records (settings, metadata, nested statistics and intervals) that other commands read, for example `hydration --cavity-results`, `region-compare` and `trajectory --lining-evidence` |
| `*_scene.json`, `*_analysis_scene.json`, `*_review_scene.json`, `*_volume_metadata.json`, `*_region_annotations.json`, `*.identities.json` | scene definitions and identity maps read by the viewer scripts |
| `*_region.json`, `*_review.json` | versioned region definitions and their review diagnostics |
| `-o` results of `profile`, `cavities`, `tunnels`, `network`, `features`, `trajectory`, and `--network-json`, `--report-json`, `fetch --json`, `benchmark --output` | the file you named; its tables are also written as CSV |

NPZ and DX files hold grids and arrays (measured and display maps, density,
per-frame arrays) that viewers and later commands load; their tabular content is
also in the CSV files.

No command writes HTML: every value is in the CSV and JSON files above.

`residue-evidence --stick-residues` controls which residues its viewer
overlays draw as sticks (default `top:12`; see
[`residue-evidence`](../tools/residue-evidence.md#choosing-which-residues-are-drawn-as-sticks)).

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

### The `publish` fallback

When `publish` cannot resolve a unique through-profile, it warns, records the
reason and writes a 3D rolling-probe cast. It does not make up a profile. The
requested cavity, tunnel and network analyses still run.
`<prefix>_profile_status.json` and the manifest's `profile_status` and
`fallback_analyses` entries record what was completed, what was skipped on
request and what could not be completed. See
[`publish`](../tools/publish.md#when-no-through-profile-resolves).

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
draws them. With `--lateral-exits`, exit legs are drawn as casts;
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

The colours and representations of every figure and scene are listed in the
{doc}`tool pages <../tools/index>` and in the API documentation of each writer.

## Viewer scripts

Cast, evidence and hydration bundles include scene scripts for PyMOL (`.pml`),
VMD (`.vmd`/`.tcl`) and ChimeraX (`.cxc`). Keep the whole output directory
together, because scripts refer to companion meshes, maps, PDB files and
identity sidecars by relative path. `crevice_render` (PyMOL/VMD) and `*_render.cxc`
(ChimeraX) make high-resolution images.

A generated script has **not** been checked in a viewer just because CREVICE
wrote it. Open it in the viewer and inspect the result before relying on it.
The source repository's `scripts/check_pymol.py`, `check_vmd.py` and
`check_chimerax.py` helpers run a scene headlessly and check its coordinates,
colours and surfaces (see [Contributing](../development/contributing.md)).

## Measured maps versus display surfaces

Cast bundles contain a measured binary map (`*_volume.dx`), where each grid
sample is either occupied or empty, and separate smoothed display fields
(`*_display.dx`, PyMOL meshes). Volumes, radii and residue attribution come from
the binary map only. Display smoothing changes how the surface looks, not what
was measured. Give the binary map, not a display map, to commands that take
`--volume-dx` or `--reference-volume-dx`.

`cast` and `publish` take `--smooth X` (X in Å, default 0.4; `0` shows the raw
voxel boundary; Python `smooth=X`). It applies constrained Gaussian twicing on
a supersampled display grid, so every measured sample keeps its inside/outside
state and the surface moves only within boundary grid cells. The simpler
single-Gaussian `--surface-smoothing W` is an alternative and cannot be combined
with `--smooth`. The chosen method and
width are recorded under `display_smoothing` in `*_volume_metadata.json`. See
[Display-surface smoothing](../SURFACE_SMOOTHING.md).

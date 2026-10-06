# `region-*`: named, versioned analysis regions

The four `region-*` commands define an analysis region as a versioned JSON
record, bind it to the input files by hash, and run and compare trajectory
analyses on it. A candidate region is never promoted to a reviewed biological
assignment automatically.

| Command | Purpose |
|---|---|
| `region-init` | Create a candidate region request from residue landmarks or an existing reference map |
| `region-prepare` | Check the request, bind it to a reference frame, and write review scenes and diagnostics |
| `region-trajectory` | Run the standard geometry, hydration, density, interaction and viewer analysis on a prepared region |
| `region-compare` | Compare two separately prepared regions on identical frames |

## Scientific question

When a study has more than one candidate site (for example, an existing pocket
definition and a landmark-based sugar-site candidate), how do their
geometry and hydration compare over *the same* frames? And how can anyone later
check exactly which inputs, selections and settings defined each site?

## Method

- **Definition (schema version 2).** Region ID and label; candidate or reviewed
  status, rationale and sources; topology and trajectory paths with SHA-256
  hashes; reference frame; protein and obstacle selections; assembly,
  preparation and membrane context; alignment residues; grid spacing, probe,
  margin, axis, phase and allocation limit; landmarks; and `radii`, the
  atomic radius set (provenance record with a table hash). Paths are resolved
  relative to the definition file. `region-init` writes a version 1 request
  (no radius set chosen yet); `region-prepare` writes the prepared definition
  as version 2 with the set it used. A version 1 definition (including one
  prepared by an earlier CREVICE) still loads and is read as the
  [default radius set](../methods/atomic-radii.md#the-default-set).
- **Landmark proposal** (`landmark_sphere`, optional). Each landmark residue's
  declared atoms give a centroid, and the anchor is the mean of those centroids.
  Atom-clear probe centres are connected along atom-clear segments, and the
  component nearest the anchor is kept. It is then swept by the probe and
  clipped to a sphere (default 6 Å). If no admissible component lies within 3 Å,
  or two are about equally close (0.5 Å threshold), the proposal **fails**
  instead of forcing a result.
- **Review scenes.** Gold sticks show residue landmarks, purple spheres show
  transferred coordinate landmarks, teal shows the region and green shows
  nearby non-protein obstacle atoms. These colours describe review evidence,
  not hydration. Residue-landmark names, the status caption and the legend
  line are drawn in the scenes only with `--annotate`; they are always
  recorded in `*_review_scene.json`, and the landmark table is
  `*_review_landmarks.csv`.
- **Comparison.** `region-compare` needs identical source hashes, frames and
  times, protein identities, alignment and relevant geometry settings. Bundles
  without region provenance need an explicit prepared definition and the
  input-hash record of the run that produced them.

Landmark file entries look like:

```json
{"residue": "SYSTEM:GLN282", "atoms": ["OE1", "NE2"],
 "evidence": "Describe the source-supported landmark and its mapping to this input."}
```

## Assumptions

- The spherical crop is an analysis choice. It is not a cavity mouth or a
  biological boundary.
- The reference is fixed during the trajectory. It is not re-centred on moving
  side chains.
- A `reviewed` status needs a reviewer, date, scope and sources. It records
  provenance and does not certify anything.

## Typical run

```bash
crevice region-init system.gro --trajectory run.xtc \
  --region-id site_candidate --landmarks-json landmarks.json \
  --obstacle-selection 'not resname TIP3' --spacing .2 --output site-request.json
crevice region-prepare site-request.json --out-dir site-review
crevice region-trajectory site-review/site_candidate_region.json \
  --out-dir site-trajectory --max-frames 1100 --workers 2 --water-membership
crevice region-compare original/original_pocket_cavity_statistics.json \
  site-trajectory/site_candidate_cavity_statistics.json --out-dir comparison
```

`TIP3` is the water residue name in a CHARMM-style topology. It is not a
universal solvent selection. Output directories must be empty, so earlier
results are never overwritten.

`--radii` on `region-prepare` selects the atomic radius set (default: the
default set, the standard table plus CHARMM36 ion radii), and the prepared
definition stores it. `region-trajectory` then measures with that stored set;
giving it a different `--radii` is an error unless `--allow-radii-mismatch`
is added (the mismatch is recorded as `radii_check` in the region provenance).
`region-compare` refuses two runs that recorded different radius sets unless
`--allow-radii-mismatch` is given. See
[Atomic radii](../methods/atomic-radii.md#reusing-earlier-results).

## Outputs

`region-prepare`: resolved definition (schema version 2, with `radii`), companion binary reference DX, identity
and coordinate fingerprint, diagnostics (`*_review.json`: status, rationale,
reference volume, obstacle overlap, warnings), `*_review_landmarks.csv` (one
row per residue landmark: `residue`, `anchor_atoms`,
`nearest_reference_sample_A`, centroid `_A` columns) and PML/VMD/TCL/CXC
review scenes. No HTML page is written.
`region-trajectory`: the full [`cavity-trajectory`](cavity-trajectory.md) and
[`hydration`](hydration.md) outputs, labelled with the region status, plus the
instantaneous-cavity water-membership files when `--water-membership` is given.
`region-compare`: paired volume and water observations
(`paired_observations.csv`: `source_frame_indices`, `time_ps`,
`left_volume_A3`, `right_volume_A3`, water counts; also as NPZ), difference
statistics (`comparison.json`, flattened to `comparison_summary.csv`: one row
per region volume, fixed-neighbourhood water count and paired difference),
`residue_comparison.csv`,
separate width profiles, boundary-residue overlap and hydration tables, and
`comparison.png` (left region brown, right region teal; region IDs in the
volume legend, drawn as panel titles only with `--annotate`).

## Python equivalent

{func}`crevice.region_definition.prepare_region`,
{func}`crevice.region_definition.load_definition`,
{func}`crevice.region_definition.trajectory_arguments` and
{func}`crevice.region_comparison.compare_regions`. These are submodule functions
that are not exported from the top-level `crevice` namespace, and `region-init`
has no standalone Python function.

## Testing and validation status

- **Software / synthetic:** schema validation, hashing, failure on ambiguous
  proposals, provenance binding and comparison contracts are covered by tests.
- **Real-input use:** a landmark candidate built from literature residues
  (see the [GLUT1 example](../REGION_WORKFLOW.md#example-literature-landmarks-for-glut1))
  has been prepared, analysed over a full MD trajectory and compared with a
  separately defined pocket.
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL, VMD and ChimeraX (review overlays) for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established. Anatomical curation, membrane
  sides and substrate accessibility need separate review.

Method note: [Named biological-region workflow](../REGION_WORKFLOW.md) and
[Biological identity and region review](../BIOLOGICAL_REGIONS.md).

## Known limitations

- The mmCIF adapter and trajectory command need reference frame/model 0.
- The fixed reference map is computed once at `region-prepare`, with the radius
  set of that run. Reading an older version 1 definition as the default set
  does not recompute its reference map.
- Different crops measure different spaces. Their volume difference is not a
  change in transport state.
- The standard water counts refer to fixed neighbourhoods. Add
  `--water-membership` to `region-trajectory` to count water inside the region
  as measured in each frame
  ([definition](hydration.md#fixed-reference-counts-and-instantaneous-cavity-membership)).
  `region-compare` does not yet compare that observable.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: region-init
```

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: region-prepare
```

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: region-trajectory
```

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: region-compare
```

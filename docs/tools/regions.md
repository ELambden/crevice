# Named regions

Sometimes a study has more than one candidate site, say an existing pocket and
a site built from literature residues, and you want to compare them fairly,
over exactly the same frames, and let anyone check later what each site was.
The `region-*` commands do this. A region is a versioned JSON record tied to
your input files by their checksums; you review it, analyse it and compare it
with another.

| Command | What it does |
|---|---|
| [`region-init`](../reference/cli/region-init.md) | start a candidate region from residue landmarks or an existing reference map |
| [`region-prepare`](../reference/cli/region-prepare.md) | check it, bind it to a reference frame and write review scenes |
| [`region-trajectory`](../reference/cli/region-trajectory.md) | run the full geometry, hydration and viewer analysis on it |
| [`region-compare`](../reference/cli/region-compare.md) | compare two prepared regions on identical frames |

## A typical run

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

`TIP3` is the water name in a CHARMM-style topology; use your own. Output
directories must be empty, so earlier results are never overwritten.

Landmarks are residues with the atoms that anchor them and a note of where the
evidence comes from:

```json
{"residue": "SYSTEM:GLN282", "atoms": ["OE1", "NE2"],
 "evidence": "Describe the source-supported landmark and its mapping to this input."}
```

In Python, the building blocks are
{func}`crevice.region_definition.prepare_region`,
{func}`crevice.region_definition.load_definition`,
{func}`crevice.region_definition.trajectory_arguments` and
{func}`crevice.region_comparison.compare_regions`; `region-init` has no
standalone Python function yet.

## How it works

- **The definition** records the region's ID, label and status (candidate or
  reviewed, with rationale and sources), the topology and trajectory with
  their SHA-256 checksums, the reference frame, atom and obstacle selections,
  alignment residues, grid and probe settings, the landmarks and the atomic
  radius set. `region-init` writes a request; `region-prepare` writes the
  prepared definition (schema version 2) with the radius set it used. Older
  version 1 definitions still load and are read with the
  [default radius set](../methods/atomic-radii.md#the-default-set).
- **Proposing a region from landmarks.** The anchor is the mean centroid of
  the landmark atoms. CREVICE finds the connected free space nearest the
  anchor, sweeps it with the probe and clips it to a sphere (6 Å by default).
  If nothing suitable lies within 3 Å of the anchor, or two candidates are
  about equally close, the proposal fails rather than guessing.
- **Review scenes** show residue landmarks as gold sticks, transferred
  coordinate landmarks as purple spheres, the region in teal and nearby
  non-protein obstacles in green. Names, the status caption and the legend
  appear with `--annotate`.
- **Radii stay consistent.** `region-prepare` stores the radius set
  (`--radii`), and `region-trajectory` uses it. Asking for a different set is
  an error unless you add `--allow-radii-mismatch`, which records the
  mismatch; `region-compare` behaves the same way. See
  [Atomic radii](../methods/atomic-radii.md#reusing-earlier-results).
- **Comparing.** `region-compare` insists on identical source checksums,
  frames, times, protein identities, alignment and geometry settings. Results
  made without region provenance need their prepared definition and input
  checksums supplied explicitly.

## What you get

- **`region-prepare`**: the prepared definition, a binary reference map, an
  identity fingerprint, review diagnostics (`*_review.json`), a landmark table
  (`*_review_landmarks.csv`) and PyMOL, VMD and ChimeraX review scenes.
- **`region-trajectory`**: everything [`cavity-trajectory`](cavity-trajectory.md)
  and [`hydration`](hydration.md) write, labelled with the region's status,
  plus the water membership files with `--water-membership`.
- **`region-compare`**: paired volumes and water counts per frame
  (`paired_observations.csv`), difference statistics (`comparison.json` and
  `comparison_summary.csv`), residue and hydration comparisons, width
  profiles and `comparison.png` (left region brown, right region teal).

:::{admonition} Interpreting results
:class: crevice-interpret

A region is an analysis choice: the spherical crop is not a cavity mouth or a
biological boundary, and the reference stays fixed rather than following
moving side chains. Two different crops measure different spaces, so their
volume difference isn't a change in transport state. Marking a region
`reviewed` records who reviewed it, when and on what evidence; it doesn't
certify it. The standard water counts use fixed neighbourhoods; add
`--water-membership` to count water inside the region as measured in each
frame ([details](hydration.md#water-in-the-cavity-two-ways-to-count)), though
`region-compare` doesn't compare that count yet. The reference frame must be
frame/model 0. The full workflow, including a worked GLUT1 landmark example, is
in [Named biological-region workflow](../REGION_WORKFLOW.md) and
[Biological identity and region review](../BIOLOGICAL_REGIONS.md).
:::

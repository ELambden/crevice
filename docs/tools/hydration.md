# `hydration`: water contacts, exposure, density and typed interactions

## Scientific question

Which residues around a cavity are in contact with explicit water, how often,
and how exposed are they geometrically? Where does water sit on average once
the protein is aligned? Which residue pairs make typed chemical contacts
(explicit-H hydrogen bonds, salt-bridge and aromatic candidates, water
bridges)?

The same analysis also runs automatically inside `analyze`, `cast`, `publish`,
`residue-evidence`, `trajectory`, `cavity-trajectory` and `region-trajectory`,
unless you pass `--skip-hydration`. The standalone command is useful when you
want to reuse finished geometry.

## Method (summary)

- **Water contact.** A water oxygen within `--hydration-cutoff` (3.5 Å) of any
  heavy atom of a residue. Each water counts once per residue per frame.
  **Occupancy** is the fraction of analysed frames with at least one such
  water. This is proximity only. It is not a hydrogen bond or a hydration free
  energy, and it does not place the water inside the cavity.
- **SASA.** Shrake–Rupley with an equal-area Fibonacci sphere
  (`--hydration-sasa-points` 256 per atom) and a `--hydration-probe` of 1.4 Å.
  All nonwater heavy atoms, including lipids, ligands and ions, are obstacles.
  Reported per residue and for backbone, side-chain and N/O subsets. The subsets
  overlap, so do not add them.
- **Trajectories.** The solvent is streamed frame by frame. Water identity is
  the topology residue index. Distances use minimum-image periodic geometry
  (orthogonal or triclinic cells). The protein is aligned with matching CA atoms
  for motion and fixed-region observables.
- **Statistics.** Every selected frame is used. Frame ranges (fluctuation) and
  approximate mean confidence intervals (contiguous-batch bootstrap) are
  reported separately. Intervals are withheld when sampling is inadequate.
- **Persistence.** Saved-frame water persistence and censored contact
  episodes. These are **not** continuous residence lifetimes: exchanges between
  saved frames (for example 100 ps apart) are not resolved.
- **Density.** Aligned 3D water-oxygen density maps in declared physical units
  (`--water-density-spacing`, a display-only `--water-density-smoothing`, and an
  absolute isovalue `--water-density-level`).
- **Typed interactions and conformation.** Explicit-H directional hydrogen
  bonds, salt-bridge and aromatic candidates, water-mediated hydrogen-bond
  edges, simplified DSSP (needs MDTraj) and side-chain torsions. Residues with
  unassigned chemistry are reported and left out of the typed results.
- **Coordinated views.** PyMOL/VMD/ChimeraX scenes with residues coloured by
  absolute contact occupancy (0% red → 50% purple → 100% blue; grey = not
  observed) and a teal mean-density surface, plus an evidence table
  (`PREFIX_analysis_evidence.csv`) and overview figures.

Full definitions: [Residue hydration and solvent exposure](../HYDRATION_METHODS.md)
and [Coordinated hydration and residue evidence](../HYDRATION_NETWORK_VIEWS.md).

## Fixed-reference counts and instantaneous-cavity membership

Two regional water observables exist. They answer different questions and are
reported separately; neither replaces the other.

- **Fixed-reference count** (always computed for trajectory regions): water
  oxygens within the region margin (2 Å) of the occupied samples of the fixed
  reference map. The neighbourhood never moves, so protein moving into or out
  of it changes the count even when no water is exchanged with the measured
  space.
- **Instantaneous-cavity membership** (`--water-membership`):
  water oxygens whose centre lies in an occupied voxel of the region measured
  **in that same frame**. Voxels are cubes of edge *h* (the grid spacing)
  centred on the fixed reference lattice; a face tie goes to the higher index,
  so every oxygen is counted once. Because the reported volume is (occupied
  samples) × h³, members/volume is a consistent number density. Counts within
  0.5 and 1.0 Å of an occupied sample centre are exported only as sensitivity
  counts.

`--water-membership` is an option of [`cavity-trajectory`](cavity-trajectory.md)
and [`region-trajectory`](regions.md), not of the standalone `hydration`
command. The measured grid exists only inside the geometry worker and is never
saved, so `hydration --cavity-results` cannot add membership afterwards. The
option is off by default, and existing output files are unchanged when it is
omitted. Each frame reuses the rigid fit and minimum-image placement of the
fixed count; a region reaching half the minimum cell height is refused, and
`--pbc unwrap` is refused. Unresolved geometry (rolling mode) gives a missing
count, never zero; a measured zero volume gives zero members and an undefined
density.

New files, written to the same output directory:

| File | Contents |
|---|---|
| `PREFIX_water_membership.json` | definitions, settings, per-series statistics with mean intervals or withholding, volume-tertile strata, correlations, saved-frame membership persistence, and the cross-check against the hydration pass |
| `PREFIX_water_membership_summary.csv` | the per-series statistics and correlations of the JSON, one row per series (columns as in {func}`crevice.io.summary_statistics_rows`) |
| `PREFIX_water_membership_strata.csv` | volume-tertile strata: frames, volume range and mean (Å³), mean member and fixed-reference counts, mean and pooled number density (Å⁻³) |
| `PREFIX_water_membership_frames.csv` / `.npz` | per-frame fixed, instantaneous, both, fixed-only and instantaneous-only counts, density, member `water-resindex` IDs and per-plane counts |
| `PREFIX_water_membership_axial.csv` | 1 Å bins along the fixed axis: mean members, wet fraction, measured volume, density, mean width and mean width in wet and dry frames |
| `PREFIX_water_membership.png` | summary figure |
| `PREFIX_analysis_member_waters_NNNNNN.pdb` | member water oxygens of each displayed viewer snapshot (source frame `NNNNNN`), aligned to the reference; drawn as the `crevice_member_waters` viewer group (orange 0.7 Å spheres, shown by default). Member IDs and counts are in the `member_waters` entry of `PREFIX_analysis_scene.json`; a snapshot with no members has no PDB |

Per frame, fixed = both + fixed-only and instantaneous = both +
instantaneous-only; in reference-region mode the discretised fixed domain also
satisfies domain = members + domain non-members. These identities are checked
and recorded in every bundle.

Membership is defined by the probe-swept geometry, grid, margin, obstacles and
alignment. It is not binding, a hydrogen bond, a hydration free energy or
permeation. An oxygen in a site too narrow for the probe, or beyond the crop,
is not a member. Definition and limits:
[Instantaneous measured-region water membership](../HYDRATION_METHODS.md#instantaneous-measured-region-water-membership).

## Choosing the residue focus

| Input | Focus |
|---|---|
| `--cavity-results stats.json` (from `cavity-trajectory`) | residues that were ever on the boundary or a nonlocal partner; reuses that run's frames, times, alignment and reference |
| `--volume-dx map.dx` (static) | boundary residues and partners of that measured map |
| `--residues-json ids.json` | an explicit list of residue IDs |
| none | all selected protein residues |

Only residues with non-HETATM heavy atoms are measured. In a static structure,
a focus residue that is present but has none (for example the D-amino acids
DLE/DVA and the modified ends FVA/ETA of 1GRM read from a PDB file, where they
are HETATM records, or a ligand, ion or water) is not measured and is listed in
the JSON as `settings.untyped_focus_residues` with a reason (`hetero_record`,
`water` or `no_heavy_atoms`); its heavy atoms still occlude SASA unless it is
water. Residue IDs absent from the structure are an error. Loading the mmCIF
with the Gemmi reader treats modified polymer residues as protein.

## Assumptions

- Water is only observed if it is present in the input. A crystal structure
  with no explicit water gives water contacts with status
  `unavailable_no_explicit_water` (null/NaN, not zero). SASA is still reported.
- Water residue names come from a standard list (HOH, WAT, SOL, TIP3, ...), or
  from `--water-selection` for MD input. Each water must have exactly one
  oxygen.
- Protein and solvent must be imaged consistently. CREVICE does not move only the
  protein away from its solvent.

## Key parameters and units

| Option | Default | Unit |
|---|---:|---|
| `--hydration-cutoff` | 3.5 | Å |
| `--hydration-probe` | 1.4 | Å |
| `--hydration-sasa-points` | 256 | points per atom |
| `--water-density-spacing` | 0.5 | Å |
| `--water-density-smoothing` | 0.5 | Å (display only) |
| `--water-density-level` | 0.05 | water oxygens per Å³ |
| `--confidence` | 0.95 | nominal level; `0` = description only |
| `--start`, `--stop`, `--stride`, `--max-frames` | 0, end, 1, 500 | frames |
| `--radii` | standard table | radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file ([Atomic radii](../methods/atomic-radii.md)) |

## Outputs

Tables (CSV; column names carry units such as `_A2`, `_ps`, `_degrees`):

| File | Contents |
|---|---|
| `*_hydration_residues.csv` | per residue: water-count, SASA and occupancy means, SDs and intervals |
| `*_hydration_summary.csv` | every interval statistic of `*_hydration.json` on one row (per residue and metric, per shared-water pair, unique focus waters, fixed-region waters), with frame quantiles, effective frames and interval status, plus the correlations |
| `*_hydration_frames.csv` | per frame and residue: `frame_index`, `time_ps`, `residue`, water counts, `sasa_A2` (all, backbone, side chain, N/O) and, with a cavity trajectory, `boundary_area_A2`, `nonlining_partner`, `centroid_distance_A`, `aligned_centroid_{x,y,z}_A` |
| `*_hydration_region_frames.csv` | per frame: `unique_focus_water_count` and, with a cavity trajectory, `region_water_count`, `cavity_volume_A3` |
| `*_hydration_shared_water.csv` | residue pairs sharing a water oxygen |
| `*_hydration_observed_episodes.csv`, `*_hydration_sampled_survival.csv`, `*_hydration_mobility.csv` | MD only: contact episodes, sampled persistence, centroid motion |
| `*_interaction_edges.csv`, `*_interaction_edge_frames.csv` | typed contacts: per pair, and per frame and pair |
| `*_residue_conformations.csv`, `*_residue_conformation_frames.csv` | torsion/secondary-structure summaries, and per-frame `phi_degrees` ... `chi4_degrees`, `secondary_structure` |
| `*_analysis_evidence.csv` | coordinated per-residue evidence (hydration, contacts, footprint, exposure, conformation) |

The same per-frame values are in `*_hydration_frames.npz` and
`*_interaction_frames.npz` (with sparse water-identity contacts). JSON files
are kept only as structured records or scene data: `*_hydration.json` (report
with settings, definitions and statistics, read by `region-compare`),
`*_interactions.json` (typed-contact record read by the scenes),
`*_water_density.json` (density grid metadata), `*_hydration_alignment.json`
(MD fit provenance), `*_analysis_scene.json` (viewer scene definition) and
`*_hydration_manifest.json` (file index, limitations and definitions). Density
DX maps, typed-interaction tables and viewer scenes are written unless
`--skip-water-density`, `--skip-interactions` or `--skip-analysis-views` is
given. A pair or episode table with no rows has only its header (or is empty).
It does not mean the site is dry. No HTML report is written.

Figures: residue water-contact bars (teal) with geometric SASA bars (orange all
heavy atoms, blue-grey N/O atoms, with a legend); for trajectories, heatmaps of
contacts (`viridis`) and SASA (`cividis`) with labelled colour bars, and
unique/fixed-region water time series; the density projections use `magma`
with colour bars. Viewer scenes: residues coloured red (0%) to blue (100%)
water-contact occupancy, grey when unobserved; boundary-contributor sticks;
pale reference cavity; teal density surface; gold/teal/purple/blue typed-contact
dashes; orange 0.7 Å member-water spheres with `--water-membership` (see
[Coordinated hydration views](../HYDRATION_NETWORK_VIEWS.md)). Titles and the
viewer colour key, its caption and the region banner are drawn only with
`--annotate`; the text is always in the PNG metadata and the scene JSON.

## Python equivalent

Static structure:

```python
from crevice import load_structure
from crevice.hydration import static_hydration
from crevice.hydration_export import write_hydration_bundle

frame = load_structure(".crevice/pdb/4PYP.cif")
result = static_hydration(frame, cutoff=3.5, sasa_probe=1.4, sasa_samples=256)
write_hydration_bundle(result, "4PYP_hydration", prefix="4PYP")
```

Trajectory: {func}`crevice.hydration_trajectory.analyze_hydration_trajectory`.
The complete CLI behaviour, including focus selection, density, typed
interactions and views, is assembled in
{mod}`crevice.hydration_workflow`. Some of its helpers take CLI argument objects,
which is a known API gap (see [CLI or Python?](../getting-started/cli-vs-python.md)).

## Testing and validation status

- **Software / synthetic:** analytic sphere and occlusion controls,
  unique-water counting, periodic orthogonal and triclinic cells, rigid-motion
  invariance, source-identity checks and standard-workflow integration.
- **Instantaneous-cavity membership, software / synthetic:** voxel and tie
  rules, an independent sphere enumeration, rigid-motion and periodic-image
  invariance (orthogonal and triclinic), moving-boundary and moving-water
  controls, zero/unresolved/no-water frames and an end-to-end CLI run
  (`tests/test_water_membership.py`).
- **Real-input use:** static runs on public entries (4PYP, 1AF6, 1GRM) and runs
  on an all-atom membrane-protein MD trajectory. Counts depend strongly on the
  contact cutoff. Member counts that simply track the measured volume at a
  roughly constant number density are descriptive and are not evidence of
  wetting or dewetting transitions.
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL, VMD and ChimeraX for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established.

## Known limitations

- The standard regional water counts refer to a **fixed reference
  neighbourhood**. Water inside the cavity as measured in each frame needs
  `--water-membership` on `cavity-trajectory` or `region-trajectory`; it cannot
  be added to finished results with `hydration --cavity-results`.
- `region-compare` does not yet compare instantaneous-cavity membership between
  regions, and probe/grid/margin sensitivity of membership has not been run.
- A saved-frame spacing of 100 ps cannot resolve fast water exchange. No
  lifetimes, rates or kinetic models are fitted.
- Correlations with motion or volume are descriptive, with no p-values and no
  causal claims.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: hydration
```

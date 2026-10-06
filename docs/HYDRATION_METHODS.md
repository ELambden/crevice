# Residue hydration and solvent exposure

Implementation: {mod}`crevice.hydration`, {mod}`crevice.hydration_trajectory`,
{mod}`crevice.hydration_export` and {mod}`crevice.hydration_workflow`. The command
is described in the [`hydration` tool guide](tools/hydration.md).

## What is measured

Explicit water contacts and geometric solvent-accessible surface area (SASA)
answer different questions. A water contact is an observed water oxygen within
3.5 Å of any selected residue heavy atom. Each water is counted once per residue
per frame, even if it contacts several atoms. A water contacting two residues is
counted for both; the total number of unique waters near the whole focus is
separately deduplicated. Contact occupancy is the fraction of observed frames
with at least one such oxygen. The cutoff is configurable and affects counts.
It is a proximity definition, not a hydrogen-bond assignment or hydration free
energy. Contacts can occur on any side of a residue; they do not by themselves
assign a water to the cavity interior. Backbone, sidechain and N/O contact sets can overlap and must not be
added to obtain a total. N/O is an element subset, not a donor/acceptor model.

SASA uses an equal-area Fibonacci sphere implementation of the
[Shrake–Rupley construction](https://pubmed.ncbi.nlm.nih.gov/4760134/): for each
protein heavy atom, sample the sphere at its CREVICE van der Waals radius plus a
1.4 Å probe radius, discard points inside another expanded atom, and multiply
the unoccluded fraction by the sphere area. The default is 256 points per atom.
Areas are summed per residue and separately for backbone, sidechain and N/O
atoms. Backbone means N, CA, C, O, OXT, OT1 and OT2; all remaining selected heavy
atoms form the sidechain subset. N/O overlaps these two groups. SASA is reported
in Å²; it is not normalized by reference amino-acid maximum exposure.

All supplied **nonwater heavy atoms** are SASA obstacles, including lipids,
ligands and ions. Waters and H/D atoms are excluded from the expanded-atom
surface. Solvent is read separately from the protein geometry so adding these
measurements does not add water atoms to the pore cast. A locally unoccluded
probe surface does not prove connectivity to bulk water. No artificial water
insertion, hydration thermodynamics or energetic ranking is performed.

## Static structures

PDB/mmCIF observations use the selected model and supplied asymmetric unit;
crystal symmetry mates are not generated. GRO snapshots use their periodic
cell and original solvent/environment coordinates. Water names recognized by
default are HOH, WAT, SOL, TIP3, TIP3P, TIP4, TIP4P, TIP5, TIP5P, SPC, SPCE, H2O
and DOD. Exactly one oxygen must identify each water residue. Zero-occupancy
waters are excluded. Positive partial occupancies are counted once and reported
as partial observations, not converted into equilibrium contact probabilities.

If no explicit waters are present, water-contact data are missing (`null` in
JSON, NaN in NPZ), with status `unavailable_no_explicit_water`. SASA is still
calculated. If explicit waters are present but none contacts a residue, its
count is zero. Neither case establishes equilibrium dryness, and a single
structure receives no time-distribution or mean confidence interval.

## Trajectories, identities and periodicity

The original solvent trajectory is streamed frame by frame; water identity is
the fixed MD topology residue index (`water-resindex:<index>`, zero based).
Protein-heavy-atom frames retain the normal CREVICE peptide-link periodicity
screen. Original protein/solvent timestamps and coordinates must agree. Missing,
nonfinite or non-increasing MD times cannot support this analysis.

Water contacts and atom occlusion use minimum-image distances with orthogonal
or triclinic cells through [MDAnalysis distance routines](https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html).
Invalid cells and local search radii reaching half the smallest cell height are
rejected. `--pbc none` explicitly requests nonperiodic distances. Split or
inconsistently reimaged protein/solvent inputs must be repaired consistently;
the reader does not silently move just the protein away from its solvent.

Proper matching-CA alignment to the reference removes whole-protein translation
and rotation for residue motion and fixed-region observables. Local occlusion
vectors are rotated into reference axes, keeping finite SASA quadrature from
changing under rigid-body rotation. Distances themselves remain periodic in the
original cell. A predeclared scaffold can be supplied through the geometry
analysis. CA fitting and peptide checks do not prove biological multimer assembly
or that the selected scaffold is stationary.

With `--cavity-results`, hydration reuses the cavity analysis frame identities,
timestamps, scaffold and reference. The default focus is the union of residues
that ever contribute boundary area or occur as nonlining structural partners.
Residues can change roles across frames. Source paths, identities and times are
checked; when the cavity results record topology/trajectory SHA256 hashes, they
are verified too. The report says whether content
hashes were available. Cached source paths must match; moved input files require
an appropriately regenerated context. Subsetting retains the original reference
frame, even when that frame is outside the requested subset.

For fixed-region cavity analyses, an additional water count includes oxygen
positions within the configured neighbourhood of the occupied reference map
samples (the region margin, for example 2 Å). Periodic waters are placed around the
aligned reference anchor. This fixed spatial observable is **not** water
membership in an instantaneous rolling-probe cast, nor a ligand-accessibility
claim. Its association with instantaneous cavity volume is descriptive.

## Instantaneous measured-region water membership

Implementation: {mod}`crevice.water_membership`. With
`--water-membership` on `cavity-trajectory` or `region-trajectory`, CREVICE also
asks which aligned water oxygens lie inside the region **measured in that same
frame**. The fixed-reference count above is unchanged and computed by the same
shared function in the same pass, so the two observables are paired per frame;
the bundle cross-checks it against the hydration pass. Neither replaces the other.

- **Rule.** An oxygen is a member if its centre lies in an occupied voxel of the
  frame's measured grid: cubes of edge *h* centred on the fixed reference lattice,
  voxel index `floor(((x - anchor)·basis)/h + 1/2)` per axis. Voxel *k* owns
  [k−½, k+½), so a face tie goes to the higher index and every oxygen falls in
  exactly one voxel. Because reported volume is (occupied samples)·h³, count/volume
  is a consistent number density. Distance rules (voxel OR within 0.5/1.0 Å of an
  occupied sample centre) are exported only as sensitivity counts.
- **Alignment and periodicity.** The rigid fit is the one the frame geometry used
  (all matching CA atoms unless a declared scaffold). Oxygens are minimum-imaged
  around the aligned anchor exactly as for the fixed count. A region reaching half
  the minimum cell height is refused rather than imaged ambiguously; `--pbc unwrap`
  is refused. Solvent-reader protein coordinates must match the geometry frame.
- **Bookkeeping.** Fixed = both + fixed-only and instantaneous = both +
  instantaneous-only per frame. In reference-region mode a discretised fixed-domain
  count satisfies domain = members + domain non-members. These identities, member
  IDs and axial-plane totals are checked and recorded in every bundle.
- **Missing and zero.** Unresolved geometry (rolling mode) gives a missing count,
  never zero, while the fixed count stays measured. Measured zero volume gives
  zero members and an undefined density. Mean intervals use the existing block
  bootstrap and are withheld with missing frames or inadequate sampling.
- **Outputs.** `_water_membership.json` (definitions, statistics, volume-tertile
  strata, correlations, saved-frame membership persistence, cross-check),
  `_water_membership_frames.csv/.npz` (per-frame counts, member `water-resindex`
  IDs, per-plane counts), `_water_membership_axial.csv` (1 Å bins along the fixed
  axis: mean members, wet fraction, measured volume, pooled density, mean width and
  mean width in wet/dry frames) and a figure. None of these files is written, and
  no other output changes, when the option is omitted.
- **Viewer group.** Standard hydration scenes gain a separate `crevice_member_waters`
  object (PyMOL), molecule (VMD) or model (ChimeraX): orange 0.7 Å spheres at each
  displayed snapshot's member oxygens, **shown by default**. View-mode commands never
  toggle it; hide it with `disable crevice_member_waters`, `mol off $crevice_member_waters`,
  or `hide #<id> models` (or the viewers' object/model panels). A hidden group
  stays hidden when `crevice_frame`/`creviceFrame` changes snapshot, and it is then
  reloaded with that snapshot's members. The group is absent without
  `--water-membership`.
- **Memory.** The per-frame grid lives only inside the geometry worker and is
  never saved, so membership is computed there from a streamed one-frame water
  snapshot. It therefore needs a geometry run; `crevice hydration --cavity-results`
  cannot add it afterwards.

Limits: the measured region is a probe-swept geometric space defined by probe,
grid, margin, obstacles and alignment. An oxygen in a site too narrow for the probe,
or beyond the crop, is not a member. Membership is not binding, a hydrogen bond, a
hydration free energy or permeation. A change in the fixed count caused by protein
moving through a fixed crop is not a hydration transition. Counts conditional on
geometry, strata and correlations are descriptive. The 100 ps spacing does not
resolve exchange. The controls in `tests/test_water_membership.py` (analytic
ties, independent voxel enumeration, periodic shifts, moving boundaries and
empty frames) are synthetic software tests, not real-system validation.
`region-compare` does not compare membership.

## Temporal statistics and motion

Every selected saved frame is used. Tables include count, mean, sample SD,
empirical 2.5–97.5% frame range, missingness and first/second-half means. These
frame ranges describe fluctuations; they are not confidence intervals for a
mean. Approximate mean intervals use the existing contiguous-batch bootstrap
and residual-correlation/effective-batch safeguards described in
[Cavity trajectory methods](CAVITY_TRAJECTORY_METHODS.md). Confidence bounds are
restricted to physical count/area or occupancy domains. Nominal pointwise 95%
intervals appear only when estimable; they do not provide simultaneous coverage
across residues, establish stationarity or account for cutoff, radii, force-field,
quadrature or anatomical uncertainty. `--confidence 0` requests descriptions
without mean intervals. Insufficient sampling remains an explicit status.

Shared-water proximity means the same oxygen contacts two residues in a frame.
Its pair occupancy and uncertainty are exported. Adjacent and nonlocal pairs
are both eligible. Without donor/acceptor typing, protonation and angular
geometry these are not water-mediated hydrogen bonds; see the distinct
[MDAnalysis hydrogen-bond requirements](https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html).

Saved-frame water persistence uses the mean, over nonempty time origins, of the
fraction of initial water IDs present at every intervening saved frame, with
zero allowed missed frames. This matches the set-intersection definition in
[MDAnalysis correlations](https://docs.mdanalysis.org/stable/documentation_pages/lib/correlations.html).
Each lag reports its own eligible origins; populations can differ by lag.
Consecutive contact episodes retain start/end time, number of observations and
left/right trajectory-boundary censoring. The observed span is last minus first
sample time; a one-frame episode has span zero, not a measured zero lifetime.
Irregularly sampled trajectories receive no regular-lag persistence curve.

Typical production trajectories are saved every 10–100 ps. Waters may depart and
return between observations, so persistence and observed spans are not
continuous-time residence lifetimes. No exponential lifetime, exchange rate or
kinetic model is fitted. Such sampling remains useful for hydration populations
and slower changes, subject to correlated sampling and representativeness.

Aligned residue-heavy-atom centroid RMSF and following-frame displacement are
reported alongside hydration. Conditional wet-at-both-endpoints and dry-at-both-
endpoints step means retain observation counts. Associations with cavity volume,
SASA, centroid distance and subsequent centroid displacement are descriptive
Pearson correlations without p-values or causal claims. Endpoint wet/dry states
do not describe the intervening interval. Centroid motion does not isolate
sidechain torsions, rotamers or internal mobility; small conditional samples
cannot establish hydration-controlled movement.

## Standard workflows and explicit reuse

The numerical libraries hydration needs are required CREVICE dependencies;
`crevice[hydration]` adds MDTraj for DSSP context and `crevice[structures]` adds the
validated mmCIF loader. Standard `analyze`, `cast`, `publish`, `residue-evidence`,
`trajectory` and `cavity-trajectory` workflows append hydration for their residue
focus. Use `--skip-hydration` to omit it. If there is no focus, or a dependency
cannot be imported, standard analysis writes an explicit unavailable status and
retains its other outputs. The low-level `profile` and `network` commands do not
add hydration.

To reuse measured geometry from a completed `cavity-trajectory` run without
repeating the expensive cavity calculation:

```bash
crevice hydration system.gro --trajectory run.xtc \
  --cavity-results results/system_cavity_statistics.json \
  --max-frames 1100 --out-dir hydration-results --prefix system
```

The frame limit is a guard, not a stride. Standalone hydration can also use a
JSON residue list (`--residues-json`) or, for a static structure, a measured map
(`--volume-dx`) to derive a boundary/partner focus. Without either focus option,
the standalone command examines the selected protein residues. Contact cutoff,
SASA probe radius and quadrature are controlled by `--hydration-cutoff`,
`--hydration-probe` and `--hydration-sasa-points`. MD water names can be overridden
with `--water-selection`. A sequence of static files produces separate static
hydration reports; it does not invent MD sampling times.

The bundle includes a JSON report, residue/shared-water CSV tables, per-frame
CSV tables (`_hydration_frames.csv`, `_hydration_region_frames.csv`), compressed
all-frame arrays and sparse water-identity contacts, plots, alignment/cell
metadata, and (for MD) sampled persistence, censored episodes and centroid-motion
tables. An observed pair/episode table with no rows is not a dryness
inference. The hydration figures and tables are added to the geometry and viewer
outputs.

## Testing status and limitations

Analytic sphere/occlusion, unique-water, periodic orthogonal/triclinic,
rigid-motion, source-identity and standard-workflow controls are software tests,
not independent biological validation. Quadrature, contact cutoff and probe
radius all change the numbers; check their sensitivity on a few frames before
interpreting a system, and report the values used. Contact counts depend
strongly on the cutoff. Fixed-domain water counts do not measure
instantaneous-cavity membership; use `--water-membership` for that observable.

Binding-site and membrane curation, independent simulations, equilibration and
force-field validation are needed for functional conclusions. More finely saved
trajectories are needed to investigate rapid water exchange. The explicit-H
hydrogen-bond and typed-interaction extension is documented in
[Coordinated hydration and residue evidence](HYDRATION_NETWORK_VIEWS.md); the
proximity observables defined here do not themselves assign hydrogen bonds.
Broader chemistry and protonation validation and energetic models are not
provided.

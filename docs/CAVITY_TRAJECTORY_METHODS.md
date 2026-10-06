# Cavity trajectories: measurement and interpretation

The `cavity-trajectory` command analyzes every requested frame after a proper
matching-CA fit to the first frame. The reference DX must already share that
first frame's coordinates. It must be the measured binary map, not the smoothed
viewer map. Actual timestamps, original residue identities, fit diagnostics,
reader selection and frame coverage are retained. A frame-count guard rejects
an oversized request; it never silently subsamples it.

## Two explicit geometric observables

`--geometry-mode rolling` recomputes global rolling-probe components on a fixed
reference lattice in each frame. Candidate identity uses the harmonic mean of
two nearest-sample overlap fractions (candidate precision and reference recall),
with a distance tolerance of 0.9 times the larger spacing. Defaults require a
score of at least 0.1, and reject a second candidate scoring at least 0.8 times
the winner. These thresholds are geometric heuristics, not biological identity
validation. All candidates and scores are exported. An absent or ambiguous
match is missing, not zero volume or evidence of a closed gate. Splitting,
merging and expansion can defeat correspondence even with alignment.

`--geometry-mode reference-region` measures a predeclared local region. The
measurement domain consists of grid samples within `--region-margin` Å of the
occupied reference samples. Eligible probe centres have atom-surface clearance
at least `--probe-radius` and lie within margin + probe radius of the reference.
Centre connectivity is tested along atom-clear grid edges. Every centre
component intersecting a 0.9-grid-spacing neighbourhood of the reference is
retained. Their fixed-radius sphere sweeps are united and clipped to the
measurement domain and nonnegative atom clearance. The volume is the number of
retained samples times spacing cubed. Component counts are reported separately.
A zero measures no retained swept samples in this specified region; it does not
establish transport closure. This observable can remain defined when one
connected pocket cannot be followed.

Both use protein heavy-atom obstacles by default. Filtering waters, ions and lipids gives a
clean protein geometry, but removes membrane and ligand barriers. Neither mode
establishes substrate accessibility, substrate shape compatibility, alternating
access, or a complete permeation route. Keep that distinction when discussing a
transporter or channel. Anatomical region and scaffold curation remain
necessary. Other physical obstacles can be selected in reference-region mode;
see [Physical obstacles and numerical sensitivity](CAVITY_OBSTACLES.md).

## Fixed coordinates and sectional width

The basis, origin and grid phase are fixed after alignment. `--axis x,y,z`
sets the third basis direction explicitly; the default is the map's third
direction. The lattice anchor is the centroid of occupied reference samples.
Changing spacing, direction or phase resamples geometry; an identical protein
need not reproduce the input map's volume on a different lattice.

At every fixed axial coordinate CREVICE retains total planar area, the area of
the largest face-connected planar component, component count and two widths:

- **Area-equivalent diameter:** `2 sqrt(A/pi)` for the largest connected planar
  section. It describes cross-sectional size; it is not an inscribed diameter.
- **Local atom-clear sphere diameter:** twice the maximum exact atom-surface
  clearance among measured grid centres in the plane. The sphere can extend
  outside the measured region. This is local capacity, not a connected path.

Radius figures divide diameters by two. A successfully measured region with no
section at a coordinate contributes zero. An unresolved frame contributes NaN
at every coordinate. Profile coverage and section-presence fractions accompany
the plots. Per-frame axial extents and volume outside the profile window are
exported. These geometric extents are not assigned biological channel mouths.
The `trajectory` through-channel analysis is a separate method.

## Sampling uncertainty and convergence

Every scalar series has coverage, mean, sample SD, median, empirical 2.5–97.5%
quantiles, extrema, first/second-half means and (when estimable) statistical
inefficiency/effective frame count. All supplied frames, including initial
frames, are retained; no equilibration period is silently discarded.

Light plot shading is the observed 95% frame range. It describes fluctuations,
not confidence in the mean. Dark shading is an approximate pointwise mean
interval from contiguous batch means, with the existing autocorrelation-aware
studentized bootstrap, minimum eight batches/effective batches, and Student-t
critical-value floor. Auto batch length is the larger of ceil(N^(1/3)) and
ceil(5 times statistical inefficiency). Missing observations, irregular time
spacing, inadequate batches and observed-constant series can suppress mean
intervals. Constant observed contact occupancy does not establish certainty
about unobserved transitions.

The joint profile family uses the most correlated coordinate to size
batches and may be unavailable even when some individual coordinates can be
estimated. Additional `section` and `free_sphere` reports estimate each
coordinate separately and explicitly make **no simultaneous coverage claim**.
Bounds are conditional on representative stationary sampling; a passing
estimator is not proof of equilibration. Inspect drift and compare independent
replicas before functional conclusions. These intervals exclude discretization,
probe/region choice, force-field and anatomical errors. See
[Grossfield and Zuckerman (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2865156/)
for the distinction between correlated sampling and broader uncertainty.

## Changing residue and network evidence

The measured binary surface is triangulated afresh in every measured frame.
Nearest selected heavy-atom VDW surfaces within 1.5 Å receive triangle area;
ties are shared. Unassigned area remains in the denominator. For fixed regions,
triangles whose centres lie within 0.9 grid spacing of the domain crop boundary
are excluded from residue attribution; their area is reported separately.
This conservative exclusion avoids assigning artificial cut faces to protein.
The resulting footprint is not SASA or an interaction energy.

A direct contributor has positive assigned boundary area in that frame.
Residue outputs include its frequency, changing area/fraction, confidence status
and all-frame frequency bounds when geometry is missing. A nonlining structural
partner directly contacts a current contributor, with minimum selected atom
centre distance <=4.5 Å and either a different chain or separation greater than
two positions in observed residue order. These roles can change between frames.

All physical contact pairs are accumulated across all frames independently of
cavity resolution. Reported quantities include occupancy and its uncertainty,
distance, contact-type fractions, dynamic boundary-partner frequency, aligned
residue-centroid RMSF and whole/split-half motion correlation (DCCM). Contact
occupancy includes any contact type; changing chemical labels do not split one
physical pair. The reference-lining role is reported separately from
per-frame boundary membership.

Volume–boundary-area and volume–centroid-distance correlations are descriptive,
without p-values or causal interpretation. Boundary area and volume share
geometry, so their correlation is not independent mechanistic evidence. The
network identifies candidates for follow-up, not crucial functional
residues or causal information flow. Perturbation experiments, replicas and
predeclared alternative alignment scaffolds are needed for stronger claims.

## Numerical sensitivity

Probe radius, region margin, lattice orientation and phase, and grid spacing all
change the measured volume. A coarse grid (for example 0.5 Å) can substantially
underestimate a small regional volume relative to finer grids, and even at
0.2 Å a finer 0.15 Å grid can give volumes a few per cent higher. Check these
settings on a few frames spread across the trajectory before running all frames
(see [Physical obstacles and numerical sensitivity](CAVITY_OBSTACLES.md#checking-numerical-sensitivity)),
and report the settings with every result. Agreement across settings is not a
claim of complete convergence.

The peptide PBC screen tests consecutive peptide links in every frame; it does
not prove full molecular completeness. A reference pocket obtained from a
transferred seed or a fitted axis is an anatomical candidate, not a curated
substrate-accessible cavity. Alternative references belong in the
[named-region workflow](REGION_WORKFLOW.md), alongside the original definition
rather than replacing it.

## Residue hydration

Standard cavity trajectory analysis also writes explicit water contacts and
geometric solvent exposure for the union of observed boundary residues and
nonlining partners. Existing cavity results can be reused by `crevice hydration`
without recomputing geometry. The fixed-neighbourhood water count is distinct
from instantaneous-cavity water membership, which `--water-membership` measures
inside each geometry worker (it cannot be recovered from saved statistics). See
[Hydration methods](HYDRATION_METHODS.md).

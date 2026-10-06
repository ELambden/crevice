# Trajectory profiles and residue evidence

This note covers the `trajectory` workflow (through-channel profiles across
frames), static and trajectory residue evidence, and seeded cavity selection.
Functional residue assignment and general cavity centerlines are not provided.
The separate [cavity-trajectory workflow](CAVITY_TRAJECTORY_METHODS.md)
provides fixed-region volumes, planar widths and changing residue evidence
without requiring two open ends.

## Coordinate and sampling contract

The MD reader filters the selected atoms before analysis, preserves original
frame indices and actual timestamps, and stores coordinate snapshots with one
atom template. A proper Kabsch fit uses matching CA identities, falling back to
matching heavy atoms only when fewer than three CA atoms exist. Ambiguous atom
identities, mismatched selections and collinear fits are rejected. Rotation and
translation apply to every selected atom. Fit RMSD and the explicit scaffold
are retained. The default all-CA fit is a convenience, not an assertion of a
rigid biological core. Selecting a scaffold after examining motion can bias
subsequent interpretation; predeclared structural scaffolds and sensitivity
checks are preferable. MDAnalysis documents the importance of matching atoms
and selection order in [structural alignment](https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html).

Make molecules whole before fitting. The default PBC screen tests available
bonds for long raw distances but short minimum-image distances. Without bonded
topology it tests only consecutive peptide C–N links. A passing screen does not
validate all bonds, sidechains, separated subunits or image choices. Automatic
unwrapping is deliberately restricted to a single reliably bonded fragment.
MDAnalysis describes this ordering in its [coordinate-transformation workflow](https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/).

## Profile fluctuations and mean uncertainty

Profiles share the reference axis and seed after alignment. Distributions reject
incompatible axes or laterally displaced origins; they interpolate within each
frame's measured range without extrapolation. Missing profiles remain missing
rows. Means at partial-coverage coordinates describe only resolved observations.
The plotted frame-quantile band is a distribution of observed radii, not a
confidence interval. Canonical tables retain radii; diameter figures scale by two.

Approximate mean intervals require alignment, ordered regular sampling, validated
profile geometry and complete frame coverage at the reported coordinate. The
estimator splits the full series into balanced contiguous batches, jointly
resamples their mean vectors, and uses studentized statistics. Auto target batch
length is max(ceil(N^(1/3)), ceil(5 times the largest estimated statistical
inefficiency)). At least eight batches and eight estimated effective batches
are required. Residual positive batch autocorrelation inflates standard errors.
Pointwise intervals use a Student-t critical-value floor; the simultaneous band
also uses a conservative Bonferroni-t floor. Whole batches preserve dependence
among axial coordinates. Partial or degenerate series may receive no interval.
Diagnostics include effective frames, batch lengths, half-mean differences and
which coordinates were estimated. A user-selected block length is not exempt
from the effective-batch check.

The nominal 95% level is approximate. It assumes stationary, representative
sampling and a batch length that captures temporal dependence; these conditions
are not automatically established by fitting or by the autocorrelation estimate.
Inspect block-length sensitivity, drift, equilibration and independent replicas.
The confidence band excludes force-field error, unsampled states, probe/grid
error and wrong channel assignment. Short, strongly correlated samples remain a
calibration limitation: intervals are then often suppressed, and when reported
their coverage can fall below the nominal level.
The distinction between sampling precision and broader simulation uncertainty
follows [Grossfield et al., best practices for simulation uncertainty](https://www.nist.gov/publications/best-practices-quantification-uncertainty-and-sampling-quality-molecular-simulations).

A reference that fails the monotone two-ended channel criterion cannot
supply a defensible width profile. The run records that failure and explicitly
does not attempt later profiles without a reference. It does not infer closure.
Closed pockets, branched cavities and transporters in occluded conformations
need curated 3D paths/regions and independent geometric validation.

## Original identities

The evidence CLI restores original chain/residue numbers and insertion codes
from CREVICE viewer sidecars. For example, viewer Trp380 in 4PYP maps back to
original Trp388. Reports therefore retain the original scientific labels.
Scene colouring verifies coordinates and uses either direct PDB serials or
the reversible original-to-viewer serial map. PyMOL sorting cannot silently
reassign highlights. Share the identity sidecar with a viewer PDB.

## Direct geometric contributors

Boundary analysis accepts the **measured binary** occupancy DX, rejecting a
smoothed display field. Triangle centroids on its 0.5 isosurface are attributed
to the nearest selected heavy-atom VDW sphere within 1.5 Å by default. Tied
residues share triangle area; unassigned area is retained. The resulting area
is a footprint on the sampled cavity boundary, not solvent-accessible surface
area, interaction energy or a causal importance score. Grid/probe choices and
artificial mouth/crop faces affect it. Small apparent boundary overlap is a
sampling diagnostic, not an intentional atom-overlapping cast.

For resolved profiles, a separate exact piecewise-linear segment calculation
identifies limiting atom spheres, including bottlenecks between sample points.
Residues within 0.25 Å of the minimum are reported with atom identities and
locations. Their frequency is conditional on resolved eligible frames, not the
probability that a functional gate is open or shut.

## Nonlocal contacts and movement

CLI contacts use the minimum selected heavy-atom center distance, <=4.5 Å by default.
A nonlocal pair lies on different chains or is separated by more than two
positions in observed residue order. Missing sequence stretches can affect that
order-based definition. Distance and observed separation are exported so users
can apply another criterion. The low-level `build_residue_network` default is a
van der Waals surface-gap metric; it must not be mistaken for a center-distance
cutoff, so choose the metric explicitly in scripts.

A trajectory physical contact pools all qualifying chemical labels for each
residue pair. Occupancy is contact frames divided by all analyzed frames with
the same residue set. The default persistent threshold is 75%. These are
operational settings, inspired by the contact-network approach of
[Sethi et al. (2009)](https://doi.org/10.1073/pnas.0810961106), not universal
biophysical constants. Test the distance/occupancy sensitivity for each system.
A salt-bridge label denotes a distance-based candidate; directional hydrogen
bonds, energetic interactions and protonation-state effects are not established.

Motion uses the geometric centroid of selected atoms per residue, not a mass
center. DCCM is the normalized mean dot product of aligned displacement vectors.
It is withheld without alignment or sufficient nonzero variation. Correlations
use all frames, including frames when a retained contact is absent. Split-half
values are descriptive diagnostics; uncertainty on DCCM is not estimated.
Rigid superposition can change DCCM while leaving contacts unchanged. A
contact-filtered correlation network is not causal information flow.

In the `trajectory --network-json` workflow described here, reference lining
membership is fixed from the supplied measured reference map. The separate
`cavity-trajectory` workflow computes per-frame boundary membership, lining
occupancy and evolving boundary areas; see
[CAVITY_TRAJECTORY_METHODS.md](CAVITY_TRAJECTORY_METHODS.md).
The output distinguishes lining–lining pairs from lining–nonlining partners.
Static ranks separately show boundary area and contacted lining footprint;
there is no arbitrary composite functional score. Network centrality alone is
not evidence of a crucial residue. Manual helix/strand labels help grouping but
do not supply independent mechanistic validation.

## Evidence needed for functional claims

The output supports choosing residues for follow-up. Stronger claims
need reviewed anatomical region/path definitions and obstacle selections,
converged geometry, independent trajectories, fit/contact sensitivity and
agreement with appropriate structural/functional observations. Perturbations,
mutagenesis or other independent tests are needed to argue that a residue
controls transport rather than merely contacting or moving with the boundary.
Residue contexts include PyMOL, VMD/Tcl and ChimeraX commands. Generating a
scene is not a native viewer check: open it in the viewer (or use the
`scripts/check_*.py` helpers) and inspect it before relying on it. CREVICE does
not perform causal or pathway inference.


## Alternating access and seeded cavity selection

An alternating-access transporter such as GLUT1 should not be required to
contain a simultaneous membrane-spanning open channel.
[4PYP](https://www.rcsb.org/structure/4PYP), for example, is an inward-open
structure. A static enclosed pocket or one-sided cavity is a valid geometric
result.

`crevice cast --seed x,y,z --probe-radius R` optionally retains the interior
probe-center component connected to the specified site. Coordinates and radii
are in Å. The seed itself must accommodate the probe; the default maximum grid
snap distance is 0.75 Å, and the complete snap segment must clear every selected
atom sphere. Every lattice edge also has an analytic segment-clearance check.
Use `--cavity-selection all --enclosure-fraction 0` to inspect the fixed-probe
component without the separate dominant-core qualification. In dominant mode,
the same component restriction precedes the existing core/volume selection.
Without a seed, the normal cavity selection applies. `--seed-tolerance` sets
an explicit snap bound; the algorithm never silently relocates an overlapping
seed to another pocket.

This is connectivity **inside the chosen outer-probe envelope**, not a test of
bulk-solvent access, lipid permeability or substrate transport. A larger probe can
exclude a side space behind a narrow neck; spaces joined by a sufficiently wide
neck remain connected. A missing six-neighbor route does not prove absence of
a continuum path. Protein-only filtering removes membrane atoms as obstacles,
so it cannot by itself establish which opening faces solvent rather than lipid.
A reviewed binding-site seed, membrane/obstacle domain and probe/grid sensitivity
are needed before interpreting a region as substrate-accessible. Flexible ligand
shape, orientation, hydration and conformational changes are not modeled by a
single spherical radius.

Seeded results can depend strongly on the probe radius: in a transferred
transporter pocket, for example, probes of 0.8 and 1.4 Å can retain volumes
that differ by two orders of magnitude, and larger probes can fail because the
seed no longer fits. A transferred region's overlap with a new structure is a
diagnostic, not independent validation, and a loose global fit cannot establish
a precise binding-site correspondence. Report the seed, probe, fit and grid with
every seeded volume.

## Reading the residue-evidence figure and viewer contexts

The left panel reports 100 × (area assigned to one residue)/(total measured
boundary area). Its denominator includes unassigned faces. The right panel
shows nonlining residues that contact lining residues nonlocally; its score is
100 × the sum of those contacted lining residues' boundary fractions. It counts
each contacted lining residue once per partner, not once per atom pair. Two
partners can contact the same lining residue, so their bars overlap in meaning
and must not be summed as a partition or compared directly to left-panel bars.
Both panels show at most 12 residues, using original structure identities.

Residues shown in the plots also appear as opaque gold (direct boundary) and
magenta (nonlocal partner) sticks. Background protein uses cartoons and the
context cast has 55% opacity to expose buried sidechains. The solid volume
scene is written separately. Boundary residues with contacts remain in the
direct group; the magenta display group contains nonlining partners.

A residue-evidence run with a measured `--volume-dx` produces `.pml`, `.vmd`,
`.tcl` and `.cxc` contexts, even without `--scene`. Supplying a CREVICE base scene
reuses its exact assets and camera; without one, cubic orthogonal binary maps
are exported without redetecting the cavity. Share the whole output directory,
including PDB, identity sidecar, meshes, DX maps and residue reference files.
PyMOL/VMD validate loaded atom serials and coordinates; native checks also
verify the actual displayed sticks. ChimeraX uses corresponding viewer residue
IDs. By default the sticks are the plotted top 12 residues per role;
`--stick-residues all|top:N|IDS|@FILE` chooses another set (for example every
boundary residue and partner) without changing the plot or tables. See
[`residue-evidence`](tools/residue-evidence.md).

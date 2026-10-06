# Physical obstacles and numerical sensitivity

Cavity geometry can include selected lipids, ions, ligands and buffer molecules
as physical walls, while protein alignment, residue networks and viewer cartoons
keep their own protein selection. The `cavity-trajectory` default is protein
heavy atoms. An explicit choice is recorded in the trajectory statistics and in
each frame's diagnostics. Static `publish` has different defaults: retained
heteroatoms are included unless `--exclude-hetero` is requested, whereas
`cast --selection protein` defaults to protein. Always record the command and the
loader's actual atom selection.

```bash
crevice cavity-trajectory system.gro run.xtc \
  --reference-volume-dx reference_volume.dx \
  --geometry-mode reference-region --spacing 0.2 --axis 0,0,1 \
  --probe-radius 1.4 --region-margin 2 \
  --obstacle-selection 'not resname TIP3' \
  --max-frames 1100 --workers 2 --out-dir results --prefix system
```

`TIP3` is the water residue name in a CHARMM-style topology, not a universal
solvent selection; use the name your topology uses. This example keeps all other
heavy atoms. Hydrogen and deuterium atoms are removed from the obstacle selection
by element, and every analysed protein heavy atom must remain in it. Waters are
measured by the standard hydration workflow. Selecting water as an obstacle
instead defines the instantaneous unoccupied volume, which is a different
observable. No ligand binding, charge-dependent exclusion or membrane
accessibility is inferred from van der Waals spheres.

The obstacle option supports `reference-region` geometry. Global
rolling-component selection rejects it explicitly. The reference-region mode
keeps all eligible probe-centre components that touch a declared reference
domain; it does not demand an open path across both sides of an
alternating-access protein. A domain and a spherical probe cannot establish
accessibility for a flexible substrate.

Obstacle coordinates are read one snapshot at a time, with source frame indices
and timestamps checked against the protein trajectory. Forked workers reopen the
reader. Every obstacle undergoes the protein's proper rigid transform. Orthogonal
and triclinic lattice images are enumerated in a conservatively padded regional
box. The padding uses an upper bound on protein-only clearance and covers
eligible centre-to-centre segments. No nearest-image-only assumption is used.
`--pbc none` explicitly disables imaging; absent box data also means
nonperiodic coordinates.

Boundary triangles are attributed to the nearest selected physical spheres.
Protein contributors remain in the residue evidence table. Other obstacles have
separate boundary-area records and a separate identity namespace, so a lipid
surface cannot be attributed to a nearby protein residue. The reported assigned
area includes both categories; `protein_assigned_boundary_area_A2` and
`other_obstacle_boundary_area_A2` distinguish them. Display atoms and network
nodes remain protein residues.

## Lattice controls

`--grid-phase x,y,z` shifts the sampling lattice by fractions of its spacing;
each fraction must lie in `[0,1)`. It leaves the reference samples and the
physical region unchanged. The phase is stored in the statistics metadata.
`--max-grid-points` is an explicit allocation guard (default 2,000,000), with no
automatic coarsening. A finer or larger region can need an explicitly increased
limit.

## Checking numerical sensitivity

Cavity volumes, sectional radii and boundary-residue rankings depend on the grid
spacing, lattice phase, probe radius, region margin, alignment scaffold and
obstacle selection. Before interpreting a result, rerun `cavity-trajectory` on a
short frame range (for example `--stop 3`, or on a few snapshots written to a
separate trajectory file) with one setting changed at a time:

- `--spacing` (for example 0.15, 0.2 and 0.3 Å);
- `--grid-phase 0.5,0,0`, `0,0.5,0`, `0,0,0.5` and `0.5,0.5,0.5`;
- `--probe-radius` and `--region-margin` either side of your chosen values;
- an alternative `--alignment-residues-json` scaffold;
- with and without `--obstacle-selection`.

Compare volumes, sectional radii (interpolated onto the same axial positions),
the set of boundary residues and the rank order of their boundary areas. Decide
review thresholds before looking at the results. As a guide, grid-spacing
effects are typically a few per cent of the volume and half-grid phase effects
somewhat smaller, while probe, region and alignment choices can have larger
effects; these depend on the system. Paired snapshot differences are
diagnostics, not confidence bounds, and agreement between a few snapshots does
not show equilibrium, convergence or anatomical validity. Probe, margin and
scaffold changes have no universal tolerance.

Results for different regions should not be pooled, and a static cast that keeps
deposited nonwater heteroatoms is not a same-settings comparison with a
protein-only cast.

Identity and site review for the bundled benchmark entries is described in
[Biological identity and region review](BIOLOGICAL_REGIONS.md).

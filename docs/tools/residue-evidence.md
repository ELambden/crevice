# `residue-evidence`: boundary residues and their partners

## Scientific question

Given a measured cavity map, which residues form its boundary, and how much of
the surface does each one contribute? Which other residues contact those
boundary residues away from their sequence neighbours? `residue-evidence`
answers both for a static structure. It also writes matching PyMOL, VMD and
ChimeraX residue scenes.

## Method

1. **Boundary attribution.** The measured binary map (`--volume-dx`) is
   triangulated at its 0.5 isosurface. Each triangle's area goes to the nearest
   selected heavy-atom van der Waals sphere within `--lining-distance` (1.5 Å).
   Ties share the area. Faces further away stay unassigned but remain in the
   denominator. The result is a **boundary footprint**. It is not
   solvent-accessible area, interaction energy or an importance score.
2. **Nonlocal partners.** A residue that is not on the boundary is a partner if
   its minimum heavy-atom centre distance to a boundary residue is ≤
   `--contact-cutoff` (4.5 Å), and the two are on different chains or more than
   two positions apart in residue order.
3. **Original identities.** When the input is a CREVICE viewer PDB, its identity
   sidecar restores the original chain, residue number and insertion code. For
   example, viewer Trp380 in 4PYP maps back to Trp388.

## Assumptions

- `--volume-dx` is the **measured binary** map from `cast` or `publish`, not a
  smoothed display map. Smoothed maps are rejected.
- The map and the structure share coordinates. Use the same reference structure
  that produced the map. A similarly named crystal structure is not a
  substitute.
- Mouth and crop faces can take area from nearby residues. Grid and probe
  choices change the footprint.

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--volume-dx` | required | | measured binary map |
| `--lining-distance` | 1.5 | Å | maximum face-to-VDW-surface distance for attribution |
| `--contact-cutoff` | 4.5 | Å | partner contact distance (heavy-atom centres) |
| `--residue-groups` | none | | JSON mapping residue IDs to helix/strand groups |
| `--scene` | none | | reuse an accepted CREVICE scene's surface and camera |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

## Outputs

- `*_residue_evidence.csv`: per residue: `group`, `role`,
  `boundary_area_A2`, `boundary_fraction`, `minimum_boundary_gap_A`,
  `nonlocal_lining_contact_count`, `contacted_lining_boundary_fraction`.
- `*_residue_evidence_partners.csv`: one row per nonlocal contact with a
  boundary-lining residue: `residue`, `lining_residue`,
  `interaction_candidate`, `center_distance_A`, `surface_distance_A`,
  `closest_atoms`, `lining_boundary_fraction`.
- `*_residue_evidence.json`: the full report (also read by
  `trajectory --lining-evidence` and by the overlays).
- A two-panel figure (orange-gold and magenta bars; panel and figure titles
  only with `--annotate`). The left panel is each residue's percentage of the total
  boundary area. The right panel is, for each partner, the summed boundary
  fraction of the lining residues it contacts. Right-panel bars overlap in
  meaning. Do not add them up or compare them directly with left-panel bars.
- `*_residue_context.pml`, `.vmd`, `.tcl`, `.cxc`: cartoon protein, translucent
  cast, **gold** sticks for boundary residues and **magenta** sticks for
  nonlocal partners. These role colours are different from the hydration
  occupancy colours. The overlays contain no labels.
- Standard hydration outputs for the boundary and partner focus, unless
  `--skip-hydration`.

### Choosing which residues are drawn as sticks

`--stick-residues` (Python `stick_residues=`, see
{func}`crevice.select_stick_residues`) sets which residues appear as sticks in
the PyMOL, VMD and ChimeraX overlays:

| Value | Sticks |
|---|---|
| `top:12` (default) | the 12 highest-ranked boundary residues and 12 highest-ranked partners (the plotted residues) |
| `top:N` | the top N of each role |
| `all` | every boundary residue and every nonlocal partner |
| `A:TRP388,A:VAL391` | the listed residue IDs |
| `@residues.txt` | IDs read from a file (comma- or line-separated; `#` lines ignored) |

Listed residues keep their gold (boundary) or magenta (partner) role colour.
An ID that is neither a boundary residue nor a partner is an error. The bar
chart and the CSV always show the top 12 per role and do not change. The
evidence JSON records the selection as `stick_residue_selection`.

## Python equivalent

```python
from crevice import (boundary_residue_evidence, load_structure, read_binary_dx,
                   restore_viewer_identities, write_residue_evidence_bundle)

viewer_pdb, volume = "4PYP_cast/4PYP_viewer.pdb", "4PYP_cast/4PYP_volume.dx"
frame = load_structure(viewer_pdb)
frame, identity_provenance = restore_viewer_identities(frame, viewer_pdb)
grid, origin, deltas = read_binary_dx(volume)
report = boundary_residue_evidence(frame, grid, origin, deltas,
                                   lining_distance=1.5, contact_cutoff=4.5)
files = write_residue_evidence_bundle(report, frame, "4PYP_evidence", prefix="4PYP",
                                      volume_path=volume, stick_residues="all")
```

For a resolved through-profile, {func}`crevice.profile_constriction_evidence`
lists the atoms that limit the continuous bottleneck: every residue within
0.25 Å of the minimum clearance.

## Testing and validation status

- **Software / synthetic:** attribution, tie sharing, identity restoration,
  viewer atom and serial mapping, and every `--stick-residues` form (including
  invalid selections) are covered by tests.
- **Example observation (public entry):** for the 4PYP cast,
  `--stick-residues all` shows 48 boundary residues and 88 partners
  (1031 atoms), and the default `top:12` shows 12 + 12 residues (199 atoms).
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL, VMD and ChimeraX (memberships, colours,
  coordinates and surfaces) for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established. A large footprint or many
  partners does not show that a residue controls transport.

Method note: [Trajectory profiles and residue evidence](../TRAJECTORY_RESIDUE_METHODS.md).

## Known limitations

- The bar chart shows only the top 12 residues per role, whatever sticks are
  selected.
- The order-based nonlocal rule can be affected by missing residues.
- For several cavities, run the command on each `regions/region_*.dx` map as
  well as on the combined map.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: residue-evidence
```

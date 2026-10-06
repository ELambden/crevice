# `cast`: 3D rolling-probe cavity casts

## Scientific question

What 3D space does a pocket, cavity or channel occupy inside a protein, and
what is its measured volume? `cast` fills interior space with variable-radius
spheres that do not overlap any atom. It does **not** need a straight,
two-ended channel, so it also covers one-sided pockets, closed cavities and
irregular spaces. For example, the inward-open GLUT1 cavity in 4PYP has no
through-profile.

## Method

1. **Outer envelope.** A large probe (`--outer-radius`, default 6 Å) is rolled
   over the exterior. Wherever it can reach counts as outside, and this sets the
   finite protein envelope and its mouths.
2. **Inner probe.** Small-probe centres (`--probe-radius`, default 0.8 Å) that
   fit inside the envelope are sampled on a grid (`--spacing`, default 0.5 Å).
   Neighbouring centres connect only if the segment between them clears every
   atom sphere.
3. **Dominant cores (default).** A sample is buried if at least
   `--enclosure-fraction` (default 0.9, i.e. 24 of 26 lattice rays) is blocked.
   Burial picks cavity *cores*. Cores smaller than `--core-fraction` (default
   8%) of the strongest core are rejected. Burial is used to choose cores, not
   to cut the final shape. Rules for multimer interfaces and connected lobes are
   recorded in `metadata.dominance`.
4. **Fill.** Each qualified core is swept with spheres. Each sphere's radius is
   the smaller of the local atom clearance and the distance to the outer
   envelope. Regions below `--min-component-volume` (default 50 Å³) are removed.
   The rest are ranked by volume and capped at `--max-cavities` (default 1).
   Raising the cap allows more regions only if they qualify on their own; it
   never splits one space into several.
5. **Seeded mode.** `--seed x,y,z` keeps only the probe-connected interior
   attached to a reviewed site. A seed that overlaps atoms, or is further than
   `--seed-tolerance` (0.75 Å) from an atom-clear grid point, is rejected. It is
   never moved silently to another pocket. Add `--cavity-selection all
   --enclosure-fraction 0` for a plain fixed-probe sweep.

`--cavity-selection all` instead lists every fixed-probe void (exhaustive
enumeration). A finite cap still applies. In this mode a non-zero
`--enclosure-fraction` is a hard crop: probe centres that fail the burial test
are removed before connectivity and the swept volume are formed. The cast
metadata records which use applied (`enclosure_definition`) and the number of
blocked rays required (`enclosure_required_rays`).

## Assumptions

- `--selection protein` (the default) keeps polymer heavy atoms. With Gemmi
  mmCIF input, modified polymer residues deposited as HETATM are kept. Ligands,
  cofactors, waters and membrane are **not** walls unless you pass
  `--selection all`. Removing them can let the cast grow into space that is
  blocked in reality.
- The thresholds are explicit geometric controls. They were not trained or
  biologically calibrated.
- The requested grid spacing is never coarsened. If a grid would exceed
  `--max-grid-points`, the command fails.

## Why the probe radius is fixed

`cast --probe-radius` (default 0.8 Å) is **not** chosen automatically, unlike
the enclosure probe of [`profile`](profile.md#choosing-the-enclosure-probe---enclosure-radius-auto)
(`--enclosure-radius auto`). The two probes do different jobs:

- In a profile, the enclosure probe only decides which gaps in the channel
  wall count as closed and which sections connect. The reported radii are
  atom-surface clearances of the sampled centres and do not depend on it over
  the range of probes that resolve the channel. That is why a rule can choose
  it (the smallest probe of the channel that persists over the most probes):
  the channel either resolves or it does not, and the choice can be checked
  against neighbouring probes.
- In a cast, the rolling probe **is** the measurement. Its centres define the
  space that is kept, and the swept spheres define the reported volume,
  shape and connectivity. A different probe gives a different, equally
  well-defined cast, with no single channel whose persistence could select
  one (casts also describe closed pockets and several regions). Choosing it
  automatically would silently change volumes. Keep 0.8 Å, or pass the probe
  you want and report it.

`publish` uses the same rolling cast for its fallback when no single channel
resolves; with `--enclosure-radius auto` that cast uses 0.8 Å (recorded as
`fallback_cast_probe_A`).

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--spacing` | 0.5 | Å | grid spacing |
| `--probe-radius` | 0.8 | Å | inner rolling probe; fixed, never automatic (see above) |
| `--outer-radius` | 6.0 | Å | outer envelope probe |
| `--enclosure-fraction` | 0.9 | fraction | burial threshold for cores |
| `--core-fraction` | 0.08 | fraction | minimum core size relative to the strongest core |
| `--min-component-volume` | 50 | Å³ | minimum region volume |
| `--max-cavities` | 1 | count | maximum number of regions (alias `--max-components`) |
| `--seed`, `--seed-tolerance` | none, 0.75 | Å | restrict to the space connected to one site |
| `--smooth X` | 0.4 | Å | constrained **display** smoothing only; `0` gives raw voxels ([method](../SURFACE_SMOOTHING.md)) |
| `--surface-smoothing W` | off | Å | alternative single-Gaussian display interpolation; cannot be combined with `--smooth` |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

## Outputs

A directory containing:

- `*_volume.dx`: the **measured** binary map. Use this file for volumes and for
  every command that takes `--volume-dx` or `--reference-volume-dx`.
- `*_display.dx`, `*_pymol_mesh.npz`: smoothed display surfaces. They change
  appearance only.
- `*_pore_cast.pml/.vmd/.tcl/.cxc`, `*_volume.*`, `*_end_view.*`,
  `*_render.cxc`: viewer scenes: grey-blue cartoon (25% transparent) and
  opaque cast surfaces, one colour per region (violet `#a855f7`, the
  non-channel cast colour, then purple, blue, mauve; orange for a shared core,
  translucent grey for other voids). Teal is kept for resolved channel casts.
  No text.
- `*_void_cast.csv`: one row per retained region: `region_id`, `label`,
  `kind`, `volume_A3`, `point_count`, centre (`center_x_A`, ...), clearance
  range (`min_radius_A`, `mean_radius_A`, `max_radius_A`), `spacing_A`.
- `*_void_cast.json`, `*_selection.json`, `*_input_report.json`,
  `*_region_annotations.json`: the full cast record with candidate decisions,
  the atom selection, input provenance and the display-region labels used by
  the scenes (kept as JSON because they are nested records).
- `*_geometry_reference.npz`, `regions/`: scene reference points and one
  measured map per region.
- `*_manifest.json`: every file written.

No Matplotlib preview image is written: open the viewer scenes (or render them
with `*_render.cxc`/`crevice_render`) to look at the cast.
- Standard hydration outputs, unless `--skip-hydration`.

## Python equivalent

```python
from crevice import load_structure, rolling_probe_cast
from crevice.rolling import select_cast_atoms

frame = load_structure(".crevice/pdb/4PYP.cif")
walls, intake = select_cast_atoms(frame, "protein")   # the CLI's --selection protein
cast = rolling_probe_cast(
    walls, spacing=0.5, probe_radius=0.8, outer_radius=6.0,
    selection_mode="dominant", max_components=1, enclosure_fraction=0.9,
    min_component_volume=50.0,
)
print(cast.component_count, round(cast.total_volume, 1))
```

{func}`crevice.rolling_probe_cast` defaults to
`selection_mode="all"` with no count limit, so pass the values above to match
the CLI. The complete CLI bundle (maps, meshes, scenes, manifest) comes from
{func}`crevice.figures.write_rolling_probe_bundle`. That function is not exported
from the top-level `crevice` namespace yet.

## Testing and validation status

- **Software / synthetic:** mathematical controls check narrow-gap and
  laterally open multimers, seed rejection, segment clearance and the rule that
  the requested spacing is never silently coarsened.
- **Example observations (public entries, stated settings):** 4PYP gives one
  dominant region; with all heavy atoms including deposited heteroatoms as walls
  on a 0.5 Å grid it measured 2439.6 Å³. A protein-only default cast gives a
  different volume, because the walls differ. 1AF6 gives four regions with
  `--outer-radius 10`.
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL (triangle vertices compared with the cast samples,
  with a wrong-transform negative control), VMD and ChimeraX for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established. A cast region is not the same
  thing as a biological channel or a substrate-accessible site. No numerical
  agreement with 3V or with an anatomical benchmark has been shown.

## Known limitations

- The method can miss weakly enclosed or unusually small cavities, and it can
  keep large voids that are not of interest.
- Volumes depend on the grid spacing, probe radius and obstacle selection.
  Compare several settings before you interpret numbers
  ([Physical obstacles and numerical sensitivity](../CAVITY_OBSTACLES.md)).
- Probe connectivity inside the envelope does not show that the space is open
  to bulk solvent or reachable by a substrate.
- General centrelines, branches and per-region radius profiles are not
  implemented yet.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: cast
```

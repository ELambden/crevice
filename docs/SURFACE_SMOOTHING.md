# Display-surface smoothing (`--smooth X`)

CREVICE measures pores and cavities on a binary occupancy grid: each grid sample
is either inside the selected void or not. Volumes, radius profiles,
bottlenecks, residue attribution and connectivity all use that binary grid.
Contouring a binary grid at 0.5 gives a voxel "staircase". This is most
visible in narrow channels such as 1GRM, where a 0.25 Å grid is shown at high
magnification.

`--smooth X` (Python: `smooth=X`) changes only the **display surface** that is
written for PyMOL, VMD and ChimeraX. The measured files and numbers stay the
same, byte for byte.

## What X means

`X` is a length in ångström, `X ≥ 0`. The default is `0.4`.

1. The measured binary grid `B` (spacing `h`) is linearly resampled onto a
   display grid with spacing `h/2`. Every measured sample is also a
   display-grid node.
2. The field is smoothed by *Gaussian twicing*, `G*(2B − G*B)`, where `G` is a
   Gaussian with standard deviation `X` Å. A single Gaussian blur shrinks
   convex shapes and narrow necks. Twicing cancels that shrinkage to leading
   order while still removing the voxel steps.
3. **Constraint.** Some display nodes are clamped to the measured state:
   nodes on a measured sample, and nodes inside a measured grid cell, face or
   edge whose corner samples all agree. Occupied nodes get ≥ 0.501 and empty
   nodes get ≤ 0.499. The surface can therefore move only inside grid cells
   that already contain the unsmoothed voxel boundary. As a result:
   - every measured sample stays on its original side of the surface;
   - face-adjacent occupied samples stay connected, so a one-sample neck is
     not cut;
   - empty samples between two regions stay empty, so distinct regions are
     never bridged;
   - the surface lies within one cell diagonal (`√3·h`) of the raw voxel
     surface.
4. All three viewers contour the same field at 0.5. PyMOL loads the identical
   triangles as a standalone mesh.

`X = 0` shows the raw 0.5 voxel boundary on the measured grid.

The method, `X`, the supersampling factor, the display spacing and this
definition are recorded under `display_smoothing` in `<prefix>_volume_metadata.json`.
`unsmoothed_reference_dx` names the measured binary grid, `<prefix>_volume.dx`,
which is the unsmoothed reference.

## Alternative: `--surface-smoothing W`

`--surface-smoothing W` (Python: `surface_smoothing=W`) selects a simpler
method: a single Gaussian of width `W` on the measured grid, clamped only at the
samples. On fine grids it can pin against the sample clamps and produce
terraced surfaces, which is why `--smooth` is the default.

`--smooth` and `--surface-smoothing` cannot be combined:

- On the command line, argparse rejects the pair.
- In Python, passing both `smooth` and `surface_smoothing` raises `ValueError`.

Both options reject negative, NaN, infinite and non-numeric values.

## Choosing X

In practice, choose X between about `h` and `2h`, where `h` is the cast grid
spacing:

| Cast spacing `h` | Suggested `X` | Notes |
|---|---|---|
| 0.25 Å (narrow channels, e.g. 1GRM) | 0.3–0.5 Å (default 0.4) | Above about 0.6 Å, the smoothed surface presses against the constraint. Creases and small dimples then appear, so a larger X makes the surface *less* smooth. |
| 0.5 Å (large cavities, e.g. 4PYP) | 0.4–0.6 Å (default 0.4) | At 0.8 Å and above, deviation and atom overlap grow with little visual gain. |

Examples:

```bash
# Default display smoothing (X = 0.4 Å)
crevice publish 1GRM.cif --cast-spacing 0.25 --section-spacing 0.25 --out-dir out/1GRM
# Raw voxel surface for an unsmoothed comparison figure
crevice publish 1GRM.cif --cast-spacing 0.25 --section-spacing 0.25 --smooth 0 --out-dir out/1GRM-raw
# Slightly stronger smoothing of a large 0.5 Å cavity cast
crevice cast 4PYP.cif --spacing 0.5 --smooth 0.6 --out-dir out/4PYP
```

```python
from crevice.parser import load_structure
from crevice.figures import write_static_publication_bundle

frame = load_structure("1GRM.cif")
write_static_publication_bundle(frame, structure_path="1GRM.cif", output_dir="out",
                                cast_spacing=0.25, section_spacing=0.25, smooth=0.4)
```

To check a choice, render `--smooth 0` next to your chosen X. The constriction
should keep its apparent width, no separate regions should merge, and mouth or
crop faces should stay where they were.

**A smoother picture is not a more accurate measurement.** Smoothing does not
refine the grid, and it is not evidence of numerical convergence. To test
sensitivity to the grid, rerun with a different `--cast-spacing`.

## Display checks behind the default

The default was chosen from display checks on two public structures, 1GRM
(0.25 Å grid) and 4PYP (0.5 Å grid). These are software and display checks,
not biological validation. At `X = 0.4` Å:

| | 1GRM (h = 0.25 Å) | 4PYP (h = 0.5 Å) |
|---|---|---|
| Measured samples that change side | 0 | 0 |
| Max deviation from raw voxel surface | 0.21 Å | 0.32 Å |
| Change in enclosed display volume vs raw voxel mesh | −2.0 % | −0.4 % |
| Apparent radius at measured bottleneck (measured 1.330 Å) | 1.340 Å (raw voxel mesh 1.243 Å) | — |
| Max atom-sphere penetration of display vertices | 0.10 Å (raw 0.12 Å) | 0.24 Å (raw 0.25 Å) |
| Surface roughness, as normal deviation (raw → `--surface-smoothing 0.6` → `--smooth 0.4`) | 21.0° → 14.6° → 4.6° | 17.3° → 5.6° → 3.5° |

All the 1GRM and 4PYP measurement outputs from `publish` were byte-identical
with and without smoothing.

# Curved and branched pathways

`crevice paths` follows the free space of a structure from a starting point,
such as a buried binding site, out to the surrounding solvent, without
assuming a straight channel axis. Each route is a smooth, curved centre line
with the radius of the largest atom-free sphere measured along it, like HOLE
or CAVER, and routes that leave the same start through different openings form
a tree that splits at junctions. Use it for transporters such as GLUT1, where
no straight channel resolves, for capped channels whose exits turn sideways,
and for buried sites with several ways out.

## Quick start

```bash
# from a point, searching each side of a membrane normal (z) separately
crevice paths protein.pdb --seed 44.8,55.4,64.2 --axis z -o paths

# from the centroid of named residues (atoms optional: GLN282@OE1+NE2)
crevice paths protein.pdb --seed-residues GLN282,GLN283,ASN288,ASN317,ASN415 -o paths

# from a region definition's landmarks, on MD frames, lipids as walls
crevice paths glut1_excerpt.gro --trajectory glut1_excerpt.xtc --frames 0,10,20 \
    --md-selection 'not resname TIP3' --include-hetero --pbc none \
    --region glut1_site_candidate_region.json --axis z -o glut1_paths
```

```python
from crevice import load_structure
from crevice.paths import find_paths, write_paths_outputs, paths_view

frame = load_structure("protein.pdb")
result = find_paths(frame, seed=(44.8, 55.4, 64.2), axis="z")
for branch in result.branches:
    print(branch.id, branch.side, round(branch.length, 1), round(branch.bottleneck_radius, 2))
write_paths_outputs(result, "paths", prefix="protein", frame=frame)
paths_view(result, frame).show()        # in Jupyter (crevice[view])
```

## How it works

1. A lattice (0.5 Å by default) covers the protein. A point is free when a
   probe of the probe radius fits there. By default the probe is chosen
   automatically (below); `--probe-radius 0.8` fixes it.
2. **Bulk** is the space reached by a rolling 6 Å probe from outside the
   protein, as for [casts](cast.md) and lateral exits. A route ends at its
   **mouth**: the first point in bulk, or an *exposed* point that sees bulk
   along several straight free lines.
3. From the start, CREVICE follows the **widest** routes: those whose
   narrowest point is as wide as possible, with shorter routes winning ties.
4. Routes are kept when they leave through **different openings**: their
   mouths are more than 12 Å apart and, after they separate, they stay more
   than 6 Å apart (the same rule as the lateral exits of
   [`profile`](profile.md)). The kept routes share their first part and split
   at **junctions**.
5. Each route is smoothed into a centre line that moves toward the middle of
   the passage but never narrows below the route's original bottleneck, and
   its radius is measured exactly every 0.25 Å.

With `--axis`, each side of the plane through the start is searched
separately, as a membrane protein's two sides should be: otherwise the widest
route to the far side can go around a narrow passage through a wider crevice.
A side with no route is reported as unresolved. The full definitions are in
the [method note](../methods/paths.md).

## Choosing the probe automatically

No single probe suits every structure: a potassium-channel selectivity filter
(3UKM) is only passed by a 0.5 Å probe, while in a wide vestibule a small
probe also finds narrow gaps between atoms and reports them as extra routes.
With the default `--probe-radius auto`, CREVICE tries the probes 0.5–3.0 Å
(the ladder that [`profile`](profile.md) uses for its enclosure probe) and,
on each side, keeps the **smallest probe of the route set found at the most
probes**:

- A narrow filter is the only way out, so its routes are the ones found at
  every probe up to their width, and a small probe is chosen.
- A wide vestibule's routes are found at many probes, while narrow side
  routes appear only at the smallest ones, so a larger probe is chosen and
  the side routes are left out (they stay listed in `*_path_probes.csv`).
- With `--axis`, each side gets its own probe.
- If no probe gives a route, the result is unresolved and lists the probes
  tried.

The clearance field is computed once and shared by all probes, and probes
wider than the widest route found at the smallest probe are skipped, so the
whole ladder costs one route search per remaining probe; `--workers` runs
them in parallel (default: up to 4). On a GLUT1 membrane frame this takes
about 50 s with 4 processes instead of 20 s for one fixed probe; on 3UKM,
30 s instead of 26 s (timings in the [method note](../methods/paths.md)).
The command prints the choice, for
example `probe=auto chosen: + 0.5 Å, - 0.9 Å`. Give a number, such as
`--probe-radius 0.8`, to fix the probe; the result is then exactly what that
probe gives.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--seed`, `--seed-residues`, `--region` | one required | where the routes start |
| `--axis` | none | search each side of the start separately (`x`, `y`, `z` or `X,Y,Z`; write `--axis=-1,0,0` for a leading minus) |
| `--probe-radius` | `auto` | narrowest clearance a route may pass; `auto` chooses it per side (above), a number fixes it |
| `--workers` | automatic | processes for the automatic probe's candidate searches (results do not depend on it) |
| `--spacing` | 0.5 Å | lattice spacing; finer finds narrower passages and costs more |
| `--bulk-radius` | 6 Å | the rolling probe that defines bulk, and the distinct-route distances |
| `--include-hetero` | off | treat ligands, or lipids and ions in MD input, as walls |
| `--max-paths` | 8 | routes per search |
| `--trajectory`, `--frames` | none | analyse chosen frames of a trajectory |
| `--html` | off | also write an interactive view with one layer per branch |

All options are listed under [`crevice paths`](../reference/cli/paths.md).

## What you get

- **`*_paths.json`**: the whole tree: start, branches (points, radii, arc
  lengths, bottleneck, mouth, side, lining residues), junctions, segments,
  unresolved cases and every setting and definition.
- **`*_paths.csv`**: one row per branch: length, bottleneck radius, its
  position and arc length, start and mouth, end type, the probe it was found
  with, lining and bottleneck residues. An unresolved run writes one row with
  its status.
- **`*_path_probes.csv`** (automatic probe only): one row per candidate probe
  and side: status (`resolved`, `pruned`, unresolved), number of routes and
  their widths, end types and lengths, the route-set group and whether it was
  chosen. The JSON's `settings.probe_selection` has the same, plus the rule
  and the reason for the choice.
- **`*_path_points.csv`**, **`*_path_lining.csv`**,
  **`*_path_junctions.csv`**: samples along each branch, lining residues with
  their positions along the branch, and junctions with their arc length on
  each branch.
- **`*_paths.png`**: radius against arc length, one line per branch in its
  own colour and line style, junctions as hollow circles, bottlenecks as dots
  and the probe radius as a dotted line.
- **Viewer scenes**: `*_paths.pml` (PyMOL), `*_paths.cxc` with one
  `*_path_<n>.bild` per branch (ChimeraX) and `*_paths_vmd.tcl` (VMD), each
  with one object per branch in colour-blind-safe colours and no labels, over
  a grey cartoon. These scripts are generated; they have not been run in
  every viewer.

:::{admonition} Interpreting results
:class: crevice-interpret

A route is a geometric passage for the chosen probe, not a demonstrated path
for a substrate, ion or water. Radii, lengths and even whether a route exists
depend on the probe (the automatic choice included), the lattice spacing (a passage narrower than the probe radius
plus up to 0.87 lattice steps can be missed when no lattice point near its
axis clears the probe) and the start
point; check a finer `--spacing` before concluding that a side is closed.
Leaving out the membrane can open routes into the bilayer. Lengths near the
mouths are uncertain by a few Å.
:::

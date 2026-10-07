# Access tunnels

`crevice tunnels` starts from a point inside the protein, such as a buried
active site, and finds the routes that lead out to the surroundings. For each
route it tells you how narrow it gets at its tightest point. Among all the
routes a probe can follow, CREVICE keeps the *widest*: those whose narrowest
point is as wide as possible, with shorter routes winning ties.

## Quick start

```bash
crevice tunnels 4PYP -o 4PYP_tunnels.json --pdb 4PYP_tunnels.pdb
```

```python
from crevice import find_tunnels, load_structure, write_tunnels_json

frame = load_structure(".crevice/pdb/4PYP.cif")
tunnels = find_tunnels(frame, start=None, spacing=2.0, min_radius=0.8, max_tunnels=3)
write_tunnels_json(tunnels, "4PYP_tunnels.json")
```

Without `--start x,y,z` (or with `start=None`), the search starts from the
centre of the selected atoms. For a real question, start from the site you
care about.

## How it works

A grid covers the atoms with 3 Å of padding. A grid point can be used if a
probe of radius `--min-radius` fits there, and two points connect only if the
straight line between them also clears every atom. The start point is joined
to the nearest usable point by a checked segment, and the search follows the
widest routes to the edge of the grid, returning up to `--max-tunnels` of them.

- `bottleneck_radius` includes the narrowest points between grid points and on
  the link from the start, so it can be smaller than any radius listed along
  the path.
- `mean_radius` is a plain average over the path points, not weighted by
  length.

Heteroatoms are left out unless you pass `--include-hetero`.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--start` | centre of the atoms | where the routes start, `x,y,z` in Å |
| `--spacing` | 2.0 Å | grid spacing |
| `--min-radius` | 0.8 Å | narrowest clearance a route may have |
| `--max-tunnels` | 3 | how many routes to report |

In Python, `max_grid_points` (120 000) caps the grid; a box that would need
more points gets a coarser spacing, and the JSON records both the requested
and the spacing actually used. All options are listed under
[`crevice tunnels`](../reference/cli/tunnels.md).

## What you get

- **`-o tunnels.json`**: every route with its points, the clearance of each
  segment, the link from the start, and the requested and effective spacing.
  Each point's `axis_position` is the distance travelled from the start.
- **`tunnels.csv`**: one row per route with its start and end, `length_A`,
  `bottleneck_radius_A`, `mean_radius_A`, `point_count`, the spacing used and
  the nearest residues.
- **`tunnels_points.csv`**: one row per point along each route, with its
  radius, distance from the start and nearest atom and residue.
- **`--pdb`**: the routes as dummy atoms, ready to open in a viewer.

There is no tunnel figure; plot the CSV or open the PDB.

:::{admonition} Interpreting results
:class: crevice-interpret

Routes depend on the start point, the grid spacing and the minimum radius, so
try a couple of settings. Reaching the edge of the grid stands in for reaching
bulk solvent; the box is not a membrane or solvent model, and a coarse grid can
miss a narrow opening that exists in reality. A route is a geometric
possibility, not evidence of transport, hydration or conductance. Routes are
not clustered or tracked across trajectory frames by this command.
:::

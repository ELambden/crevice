# `tunnels`: widest grid paths to the outside

## Scientific question

Starting from a point inside the protein, which routes lead out to the
surroundings, and how narrow is each one at its tightest point? `tunnels` runs a
widest-path search on a clearance grid, in the spirit of CAVER.

## Method

A grid covers the selected atoms with `padding` (3 Å) around them. Nodes need
atom-surface clearance of at least `--min-radius`. Edges connect only if the
whole segment between two nodes clears every atom sphere. The requested start
(`--start x,y,z`, or the centroid of the selected atoms by default) is joined
to the nearest reachable node by a checked segment. The search returns up to
`--max-tunnels` paths from the start to the grid boundary. Boundary endpoints
are ordered by their widest-path bottleneck, with shorter paths first when
bottlenecks tie.

- `bottleneck_radius` includes minima between nodes and on the seed
  attachment. It can be smaller than every node radius shown.
- `mean_radius` is the plain mean over the returned nodes. It is not weighted by
  path length.

## Assumptions

- Reaching the grid boundary stands in for reaching bulk solvent. The grid box
  is not a membrane or a solvent model.
- HETATM records are dropped unless you pass `--include-hetero`.
- A path that the grid does not find may still exist in continuous space. A
  coarse or offset grid can miss a narrow opening.

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--start` | centroid | Å | start coordinate `x,y,z` |
| `--spacing` | 2.0 | Å | requested grid spacing |
| `--min-radius` | 0.8 | Å | minimum clearance to pass |
| `--max-tunnels` | 3 | count | number of paths returned |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

The Python function has a `max_grid_points` cap (default 120 000). If the box
would need more points, the spacing is **scaled up** to fit. The JSON records
both the requested and the effective spacing, so check them before you compare
runs.

## Outputs

- `-o tunnels.json`: nodes, per-segment clearances, seed attachment, requested
  and effective spacing, and the definitions used for the estimates. For each
  path point, `axis_position` is the cumulative path length in Å from the start
  (0 at the start, equal to the tunnel `length` at the end) and `index` is its
  step number along the path.
- `tunnels.csv` (beside `-o`, same stem): one row per tunnel: `tunnel_id`,
  start and end coordinates (`start_x_A`, ..., `end_z_A`), `length_A`,
  `bottleneck_radius_A`, `mean_radius_A`, `score`, `point_count`,
  `grid_spacing_A`, `nearest_residues`.
- `tunnels_points.csv`: one row per path node: `tunnel_id`, `index`, `x_A`,
  `y_A`, `z_A`, `radius_A`, `raw_clearance_A`, `axis_position_A` (cumulative
  path length), nearest atom and residue.
- `--pdb`: dummy atoms along each path.

There is no tunnel figure (the former `--png` summary was removed); plot the
CSV or open the PDB in a viewer.

## Python equivalent

```python
from crevice import find_tunnels, load_structure, write_tunnels_json

frame = load_structure(".crevice/pdb/1GRM.cif")
tunnels = find_tunnels(frame, start=None, spacing=2.0, min_radius=0.8, max_tunnels=3)
write_tunnels_json(tunnels, "tunnels.json")
```

## Testing and validation status

- **Software / synthetic:** tests cover constrictions between nodes, seed
  attachment, grid phase and diagonal collisions. A 48-case resolution and
  offset sweep on mathematical walls is available
  (`examples/connectivity_validation.py`).
- **Example observation (public entries):** in example `publish` runs, 1GRM
  and 4PYP each returned three paths.
- **Biological / functional:** not established. A grid path is not evidence of
  transport, hydration or conductance.

## Known limitations

- The result depends on the start point, spacing and minimum radius.
- Paths are not clustered or tracked across trajectory frames by this command.
- Automatic spacing scaling can make a fine request coarser without stopping
  the run.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: tunnels
```

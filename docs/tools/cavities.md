# `cavities`: enclosed grid voids

## Scientific question

Does the structure contain enclosed empty spaces that no probe can reach from
outside at the chosen resolution? `cavities` is a grid search in the spirit of
HOLLOW. It is a simpler counterpart of [`cast`](cast.md), and it uses
different definitions. Its output is not interchangeable with `cast` regions.

## Method

Grid points with atom-surface clearance of at least `--min-radius` are free.
Free points connect only if the segment between them clears every atom sphere,
including diagonal steps. An exterior flood fill starts from the grid boundary.
Free points it cannot reach are enclosed. Connected enclosed components are
reported as cavities, ranked by size and capped at `--max-cavities`. Open
channels and tunnels are left out on purpose.

## Assumptions

- The walls are protein heavy atoms. HETATM records are **dropped** unless you
  pass `--include-hetero`, and hydrogens are dropped unless you pass
  `--include-hydrogen`.
- Enclosure is judged at one resolution and one probe size. A coarse or offset
  grid can miss a narrow opening, so the space may look enclosed when it is not,
  or open when it is not.

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--spacing` | 2.0 | Å | grid spacing |
| `--min-radius` | 1.4 | Å | minimum atom-surface clearance for a free point |
| `--max-cavities` | 10 | count | maximum cavities reported |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

The Python function also has a `max_grid_points` cap (default 75 000). In
`crevice analyze`, that cap (`--cavity-max-grid-points`) **coarsens the requested
spacing** to stay under it. Raise the cap if you need a fine spacing. This is
different from `cast`, which never coarsens the grid.

## Outputs

- `-o cavities.json`: for each cavity, `center` (clearance-weighted centroid of
  its grid points), `volume` (grid points × spacing³, Å³), `radius` (largest
  clearance, Å), `grid_points`, `nearest_residues` and `score`. Here `score`
  always equals `volume`; it is a ranking value, not a probability or energy.
- `cavities.csv` (written beside `-o`, same stem): the same table, one row per
  cavity: `cavity_id`, `center_x_A`, `center_y_A`, `center_z_A`,
  `max_radius_A`, `volume_A3`, `grid_points`, `score`, `nearest_residues`
  (semicolon-separated).
- `--pdb`: dummy-atom points.
- `--png`: summary figure of the top 10 cavities by score (that is, by volume):
  volume (Å³, teal bars, left axis) and maximum radius (Å, red line, right
  axis). No title unless `--annotate`.

A {class}`~crevice.models.Cavity` is a summary without its grid points. Two
Python helpers that take one use only its centre and radius:

- {func}`crevice.classify_voids` labels a cavity by the distance from its
  centre to the nearest face of the axis-aligned bounding box of all atom
  centres (`surface_proximal` within `surface_accessibility`, default 6 Å),
  else by its radius (`large_connected_candidate` if at least that threshold,
  otherwise `buried_cavity`). The bounding box is not the molecular surface,
  and no label tests connectivity to bulk. No command uses these labels.
- {func}`crevice.build_cavity_network` with a `Cavity` links residues within
  the cutoff of the cavity **centre** only. For contacts over the whole cavity
  surface, pass the cast's {class}`~crevice.models.VoidComponent` instead.

## Python equivalent

```python
from crevice import detect_cavities, load_structure, write_cavities_json

frame = load_structure(".crevice/pdb/1GRM.cif")
cavities = detect_cavities(frame, spacing=2.0, min_radius=1.4, max_cavities=10)
write_cavities_json(cavities, "cavities.json")
```

## Testing and validation status

- **Software / synthetic:** tests cover exterior flood fill, diagonal
  collisions and runs without SciPy.
- **Example observation (public entries):** in example `publish` runs on 1GRM
  and 4PYP, the enclosed-cavity search found **zero** cavities in both. Zero
  enclosed voids at one resolution does not mean the protein has no cavity. For
  example, the 4PYP [`cast`](cast.md) region opens to one side.
- **Biological / functional:** not established.

## Known limitations

- The default 2 Å grid is coarse. Results can change with spacing and grid
  offset.
- An empty result does not prove the protein has no cavity.
- Open pockets are outside this definition. Use [`cast`](cast.md).

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: cavities
```

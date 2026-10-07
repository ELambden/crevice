# Buried cavities

`crevice cavities` is a quick search for empty spaces that are completely
enclosed by the protein, the kind a probe cannot reach from outside at the
grid resolution you choose. It works on a coarse grid, so it's fast and makes a
good first look. For pockets that open to one side, and for a finer,
better-controlled shape and volume, use a [cast](cast.md); the two use
different definitions, so their results aren't interchangeable.

## Quick start

```bash
crevice cavities 4PYP -o 4PYP_cavities.json --png 4PYP_cavities.png
```

```python
from crevice import detect_cavities, load_structure, write_cavities_json

frame = load_structure(".crevice/pdb/4PYP.cif")
cavities = detect_cavities(frame, spacing=2.0, min_radius=1.4, max_cavities=10)
write_cavities_json(cavities, "4PYP_cavities.json")
```

## How it works

A grid point is free if a probe of radius `--min-radius` fits there. Two free
points connect only if the straight line between them, diagonals included,
also clears every atom. CREVICE then floods the free space inwards from the
edges of the grid: whatever the flood can't reach is enclosed. Each connected
enclosed group is a cavity, and the largest `--max-cavities` are reported.
Open channels and tunnels are left out on purpose.

The walls are protein heavy atoms. Heteroatoms (ligands, waters, ions) are
left out unless you pass `--include-hetero`, and hydrogens unless you pass
`--include-hydrogen`.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--spacing` | 2.0 Å | grid spacing |
| `--min-radius` | 1.4 Å | smallest clearance that counts as free space |
| `--max-cavities` | 10 | how many cavities to report |
| `--include-hetero` | off | treat ligands, waters and ions as walls |

The Python function caps the grid at `max_grid_points` (75 000). In
`crevice analyze` that cap (`--cavity-max-grid-points`) coarsens the spacing to
stay under it, so raise it if you need a fine grid. All options are listed
under [`crevice cavities`](../reference/cli/cavities.md).

## What you get

- **`-o cavities.json`** and **`cavities.csv`** beside it: for each cavity,
  its centre (clearance-weighted), `volume_A3` (grid points × spacing³),
  `max_radius_A`, `grid_points`, the nearest residues and a `score`, which is
  simply the volume, used for ranking.
- **`--png`**: the ten largest cavities, with volume as teal bars and maximum
  radius as a red line.
- **`--pdb`**: the cavity grid points as dummy atoms.

Two Python helpers work from a cavity's centre and radius alone.
{func}`crevice.classify_voids` labels a cavity by its distance from the edge of
the atoms' bounding box (not the molecular surface) and by its radius; no
command uses these labels. {func}`crevice.build_cavity_network` links residues
near the cavity centre; for contacts over a whole cavity surface, pass a cast's
{class}`~crevice.models.VoidComponent` instead.

:::{admonition} Interpreting results
:class: crevice-interpret

Enclosure is judged at one grid resolution and one probe size, and the default
2 Å grid is coarse: a narrow opening can be missed, or a closed space can look
open, and results can change with the spacing or a small shift of the grid. An
empty result doesn't mean the protein has no cavity. In our `publish` runs on
1GRM and 4PYP this search found none, yet 4PYP has a large cavity that opens to
one side and casts clearly. The search is tested on synthetic structures; it
has not been benchmarked against curated cavities.
:::

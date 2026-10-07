# Cavity and channel casts

`crevice cast` shows you the 3D shape of the space inside a protein and
measures its volume. It fills the interior with spheres that don't touch any
atom, so the result follows the real shape of a pocket, a buried cavity or an
irregular channel. Unlike a [profile](profile.md), a cast doesn't need a
straight passage with two open ends: the inward-open GLUT1 cavity in 4PYP,
which has no through-channel, casts happily.

## Quick start

```bash
crevice cast 4PYP --out-dir 4PYP_cast
```

```python
from crevice import load_structure, rolling_probe_cast
from crevice.rolling import select_cast_atoms

frame = load_structure(".crevice/pdb/4PYP.cif")
walls, intake = select_cast_atoms(frame, "protein")   # same atoms as --selection protein
cast = rolling_probe_cast(
    walls, spacing=0.5, probe_radius=0.8, outer_radius=6.0,
    selection_mode="dominant", max_components=1, enclosure_fraction=0.9,
    min_component_volume=50.0,
)
print(cast.component_count, round(cast.total_volume, 1))
```

In Python, {func}`crevice.rolling_probe_cast` keeps every region by default
(`selection_mode="all"`), so pass the values above to match the command. The
complete bundle of maps, meshes and scenes comes from
{func}`crevice.figures.write_rolling_probe_bundle`.

## How it works

1. **Find the outside.** A large probe (`--outer-radius`, 6 Å) rolls over the
   protein. Wherever it reaches is outside; this defines the protein's
   envelope and its mouths.
2. **Sample the inside.** Centres where a small probe (`--probe-radius`,
   0.8 Å) fits inside the envelope are placed on a grid (`--spacing`, 0.5 Å).
   Two neighbouring centres are connected only if the straight line between
   them also clears every atom.
3. **Pick the main cavity.** A centre is *buried* if at least
   `--enclosure-fraction` (0.9, so 24 of 26 directions) of the rays from it hit
   protein. Buried centres mark cavity cores; cores smaller than
   `--core-fraction` (8 %) of the strongest are dropped. Burial is used only
   to choose cores, not to trim the final shape.
4. **Fill.** Each chosen core is swept with spheres as large as the local
   clearance and the envelope allow. Regions under `--min-component-volume`
   (50 Å³) are dropped, the rest are ranked by volume, and up to
   `--max-cavities` (1) are kept. Raising the cap only admits regions that
   qualify on their own; it never splits one space into several.

Two other modes are useful:

- **Seeded.** `--seed x,y,z` keeps only the space connected to a site you
  choose. If the seed sits inside an atom, or further than `--seed-tolerance`
  (0.75 Å) from free space, it is rejected, never moved to a different pocket.
  Add `--cavity-selection all --enclosure-fraction 0` for a plain sweep from
  the seed.
- **Everything.** `--cavity-selection all` lists every void the probe finds.
  Here a non-zero `--enclosure-fraction` trims away unburied centres before
  the regions are formed; the metadata records which rule applied.

By default the walls are the protein's heavy atoms (`--selection protein`,
including modified residues deposited as HETATM when Gemmi reads the mmCIF).
Ligands, cofactors, waters and membrane only count as walls with
`--selection all`. Leaving them out can let a cast grow into space that is
really occupied, so choose deliberately. The grid is never coarsened behind
your back: if it would exceed `--max-grid-points`, the command stops.

### Why the probe radius is fixed

Unlike the enclosure probe of a [profile](profile.md#choosing-the-enclosure-probe),
the cast probe is not chosen automatically, because here the probe *is* the
measurement. In a profile, the enclosure probe only decides which gaps count
as wall, and the radii don't depend on it as long as the channel resolves, so a
rule can pick it. In a cast, the probe centres define the space that is kept
and the swept spheres define its volume and shape: a different probe gives a
different, equally valid cast. Choosing it automatically would quietly change
volumes, so keep 0.8 Å or pass the probe you want and report it. `publish`
uses the same 0.8 Å cast as its fallback when no channel resolves.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--spacing` | 0.5 Å | grid spacing |
| `--probe-radius` | 0.8 Å | the rolling probe that defines the cast |
| `--outer-radius` | 6 Å | the probe that defines the outside |
| `--max-cavities` | 1 | how many regions to keep (alias `--max-components`) |
| `--seed`, `--seed-tolerance` | none, 0.75 Å | keep only the space connected to one site |
| `--selection` | `protein` | which atoms count as walls |
| `--enclosure-fraction`, `--core-fraction` | 0.9, 0.08 | how buried a core must be, and how large relative to the strongest |
| `--min-component-volume` | 50 Å³ | smallest region kept |

`--smooth` changes only how the surface is drawn
([Surface smoothing](../SURFACE_SMOOTHING.md)). Atomic radii, figure text,
hydration and structure input are
[shared options](../reference/cli/common-options.md); every option is listed
under [`crevice cast`](../reference/cli/cast.md).

## What you get

A directory containing:

- **`*_volume.dx`**: the measured cast as a binary map. Use this file for
  volumes and for any command that takes `--volume-dx` or
  `--reference-volume-dx`.
- **`*_void_cast.csv`**: one row per region with `volume_A3`, `point_count`,
  centre coordinates and the minimum, mean and maximum clearance.
- **Viewer scenes** for PyMOL, VMD and ChimeraX (`*_pore_cast.*`,
  `*_volume.*`, `*_end_view.*`, `*_render.cxc`): a grey-blue cartoon with
  opaque cast surfaces, one colour per region, starting with violet
  (`#a855f7`). Teal is reserved for channel casts.
- **`*_display.dx`** and **`*_pymol_mesh.npz`**: smoothed display surfaces,
  for looks only.
- The full cast record, atom selection and input report as JSON, one measured
  map per region in `regions/`, standard hydration outputs (unless
  `--skip-hydration`) and a `*_manifest.json` listing everything.

There is no Matplotlib preview: open a scene, or render one with
`*_render.cxc`. A cast of this kind has no channel axis, so it isn't split into
entry, lumen and exit segments; [`publish`](publish.md) does that for a
resolved channel.

:::{admonition} Interpreting results
:class: crevice-interpret

A cast volume depends on the grid spacing, the probe radius and which atoms
count as walls: in our 4PYP example, a cast with every heavy atom (including
deposited heteroatoms) as walls on a 0.5 Å grid measured 2439.6 Å³, and the
protein-only default gives a different number. Compare a few settings before
reading much into a volume ([Physical obstacles and numerical
sensitivity](../CAVITY_OBSTACLES.md)). The method can miss weakly enclosed or
very small cavities and can keep large voids you don't care about. Probe
connectivity doesn't show that a space is open to bulk solvent or reachable by
a substrate, and a cast region is a geometric candidate, not a biological site,
until other evidence says so.
:::

The [4PYP tutorial](../examples/4pyp-cavity-residues.md) casts the GLUT1
cavity and finds the residues that form its wall.

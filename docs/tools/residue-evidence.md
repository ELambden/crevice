# Boundary residues and their partners

Once you have a cast, `crevice residue-evidence` tells you which residues form
its wall and how much of the surface each one contributes. It also finds the
residues that touch those wall residues from further away in the sequence, the
second shell that holds the wall in place, and writes matching PyMOL, VMD and
ChimeraX scenes so you can look at both.

## Quick start

```bash
crevice cast 4PYP --out-dir 4PYP_cast
crevice residue-evidence 4PYP_cast/4PYP_viewer.pdb \
    --volume-dx 4PYP_cast/4PYP_volume.dx --out-dir 4PYP_evidence --prefix 4PYP
```

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

```{figure} ../examples/images/4pyp_residue_evidence.png
:alt: Bar charts of boundary share per residue and of partner contacts for the 4PYP cast.
:width: 640px

The GLUT1 (4PYP) cavity: each wall residue's share of the boundary (left) and,
for each partner, how much of the wall it touches (right).
```

## How it works

1. **Who forms the wall.** The measured cast map (`--volume-dx`) is turned
   into a triangulated surface. Each triangle's area goes to the nearest heavy
   atom within `--lining-distance` (1.5 Å); ties share the area, and triangles
   with no atom nearby stay unassigned but still count in the total. The
   result is a *boundary footprint*: how much of the cavity surface each
   residue makes.
2. **Who holds the wall in place.** A residue that isn't on the wall is a
   *partner* if it comes within `--contact-cutoff` (4.5 Å, heavy-atom
   centres) of a wall residue and is on a different chain or more than two
   positions away in the sequence.
3. **Original numbering.** CREVICE's viewer PDB files come with an identity
   sidecar that restores the original chain, residue number and insertion
   code; in 4PYP, viewer Trp380 maps back to Trp388.

Use the measured binary map from `cast` or `publish`; smoothed display maps
are rejected. The map and the structure must share coordinates, so use the
same structure that produced the cast.

## Choosing which residues are drawn as sticks

`--stick-residues` (Python `stick_residues=`, see
{func}`crevice.select_stick_residues`) chooses which residues appear as sticks
in the scenes:

| Value | Sticks |
|---|---|
| `top:12` (default) | the 12 highest-ranked wall residues and the 12 highest-ranked partners |
| `top:N` | the top N of each |
| `all` | every wall residue and every partner |
| `A:TRP388,A:VAL391` | just the residues you list |
| `@residues.txt` | residues read from a file (comma- or line-separated; `#` lines ignored) |

Wall residues are gold and partners magenta, whichever you choose. A residue
that is neither is an error. The bar chart and the CSV always show the top 12
of each.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--volume-dx` | required | the measured cast map |
| `--lining-distance` | 1.5 Å | how far a surface triangle may be from an atom to be assigned to it |
| `--contact-cutoff` | 4.5 Å | partner contact distance |
| `--stick-residues` | `top:12` | which residues are drawn as sticks |
| `--residue-groups` | none | a JSON file mapping residues to helices or strands |
| `--scene` | none | reuse the surface and camera of an existing CREVICE scene |

All options, including the shared ones, are listed under
[`crevice residue-evidence`](../reference/cli/residue-evidence.md).

## What you get

- **`*_residue_evidence.csv`**: one row per residue with its `role`,
  `boundary_area_A2`, `boundary_fraction`, `minimum_boundary_gap_A` and how
  many wall residues it contacts.
- **`*_residue_evidence_partners.csv`**: one row per partner contact, with the
  distances and the closest atoms.
- **`*_residue_evidence.json`**: the full report, which
  `trajectory --lining-evidence` can read.
- **A two-panel figure**: each wall residue's share of the boundary (orange-gold
  bars), and for each partner the summed share of the wall residues it touches
  (magenta bars). Partner bars overlap in meaning, so don't add them up or
  compare them with the left panel.
- **Scenes** (`*_residue_context.pml/.vmd/.tcl/.cxc`): a cartoon, a
  translucent cast, gold wall residues and magenta partners, without labels.
- Standard hydration outputs for these residues, unless `--skip-hydration`.

For a resolved channel profile, {func}`crevice.profile_constriction_evidence`
lists the atoms that set the bottleneck: every residue within 0.25 Å of the
minimum clearance.

:::{admonition} Interpreting results
:class: crevice-interpret

A footprint is geometry: it is not solvent-accessible area, an interaction
energy or an importance score, and a large footprint or many partners doesn't
show that a residue controls transport. Faces at a cavity mouth or a crop edge
can pull area towards nearby residues, and the grid and probe settings change
the footprint. Missing residues can affect the sequence-distance rule for
partners. For a cast with several regions, run the command on each
`regions/region_*.dx` map as well as on the combined map. The method is
detailed in [Trajectory profiles and residue evidence](../TRAJECTORY_RESIDUE_METHODS.md).
:::

The [4PYP tutorial](../examples/4pyp-cavity-residues.md) walks through this
on GLUT1.

# GLUT1: curved paths and cavity tracking

GLUT1 is a sugar transporter with no straight channel through it: depending on
its conformation, the substrate site opens to one side of the membrane or the
other. This tutorial follows the free space from the candidate sugar site to
the solvent on both sides in several frames of a membrane simulation, and then
tracks the cavities around the site through all 21 frames of the example.

```{admonition} Status
:class: note
Run by hand on 15 October 2026 with the commands below (CREVICE 0.3.0
development branch); not yet run in CI. The numbers are geometric measurements
under the stated probes. They are not transport pathways, a conformational
assignment or a functional claim, and the membrane sides of this simulation
are not assigned to cytoplasm or exterior.
```

## Inputs

```bash
crevice fetch-example glut1-excerpt        # prints the example directory
cd ~/.cache/crevice/examples/glut1-excerpt
```

The example holds 21 frames (every 5 ns over 100 ns) of an all-atom membrane
simulation of a human GLUT1 mutant, with lipids, water and ions, and a region
definition whose landmarks are the terminal amides of Gln282, Gln283, Asn288,
Asn317 and Asn415 (the residues that coordinate the sugar in the 4PYP
structure). In this simulation the membrane normal is *z*; the phosphate planes
of the two leaflets sit at about *z* = 40 Å and 83 Å, and the site at
*z* ≈ 62 Å, in the middle of the bilayer.

## 1. Curved paths from the site to both sides

```bash
crevice paths glut1_excerpt.gro --trajectory glut1_excerpt.xtc --frames 0,5,10,15,20 \
    --md-selection 'not resname TIP3' --include-hetero --pbc none \
    --region glut1_site_candidate_region.json --axis z -o glut1_paths --prefix glut1 --html
```

- `--md-selection 'not resname TIP3' --include-hetero` keeps the lipids and
  ions as walls, so a route cannot pass through the bilayer unless the lipids
  leave a gap, and bulk is the water on either side. `--pbc none` skips the
  bond check that wrapped lipids would fail; the protein is whole.
- `--region` starts every frame at the mean of the five landmark centroids in
  that frame; the start snaps to the widest free point within 3 Å.
- `--axis z` searches the +*z* and −*z* half-spaces separately and reports a
  side with no route as unresolved.

Each frame takes about 40 s and writes `glut1_f<N>_paths.*`. What we found
(probe 0.8 Å, 0.5 Å lattice):

| Frame (time) | −*z* widest bottleneck (Å) | −*z* length (Å) | +*z* widest bottleneck (Å) | +*z* length (Å) |
|---|---:|---:|---:|---:|
| 0 (0 ns) | 2.13 | 29.3 | 1.26 | 45.3 |
| 5 (25 ns) | 2.13 | 19.3 | no route | — |
| 10 (50 ns) | 2.63 | 18.9 | 0.92 | 35.1 |
| 15 (75 ns) | 2.89 | 18.3 | 0.94 (lined by lipids) | 83.3 |
| 20 (100 ns) | 3.05 | 17.1 | 0.98 | 46.3 |

```{figure} images/glut1_paths_bottleneck_by_frame.png
:alt: Widest-route bottleneck radius on each side of the GLUT1 site in five frames
:width: 85%

Widest-route bottleneck on each side of the site. Orange squares: −*z*; blue
circles: +*z* (hollow: the bottleneck is lined by lipids); cross: no route to
bulk on that side; dotted line: the probe radius.
```

On the −*z* side every frame has a wide route (bottleneck 2.1–3.1 Å; Met420
lines it in every frame, usually with Ile168 and Phe291) that reaches the
solvent at *z* ≈ 46–48 Å, a few Å inside the lower phosphate plane, after
17–29 Å. On the +*z* side the routes are
narrow (0.9–1.3 Å, just above the 0.8 Å probe) and long; at 25 ns there is no
route at this probe, and at 75 ns the only route squeezes between lipids
rather than through the protein. Read the +*z* routes as narrow geometric
gaps, not as an open pathway, and check them with a finer lattice
(`--spacing 0.35`) or a smaller probe before drawing conclusions.

```{figure} images/glut1_f10_paths.png
:alt: Free-sphere radius against arc length for the five GLUT1 branches at 50 ns
:width: 85%

`glut1_f10_paths.png` (50 ns): one line per branch, junctions as hollow
circles, bottlenecks as dots. Branch 1 leaves to −*z*; branches 2–5 share their
first 18 Å toward +*z* and then separate.
```

`glut1_f10_paths_view.html` opens the same branches in a browser, one layer
per branch, over the protein cartoon; `crevice view glut1_paths` does the same
for the whole directory.

## 2. Tracking the cavities around the site

```bash
crevice track-cavities glut1_excerpt.gro glut1_excerpt.xtc \
    --region glut1_site_candidate_region.json --radius 8 --workers 4 \
    -o glut1_tracking --prefix glut1_excerpt
```

Every frame is aligned to the first and cast with a 1.4 Å probe (protein
atoms only); every cavity with a voxel within 8 Å of the landmark anchor is
kept and linked to the next frame's cavities by voxel overlap. The run took
76 s with 4 workers:

```text
track-cavities frames=21 measured=21 identities=97 birth=58 death=64 split=23 merge=20
```

```{figure} images/glut1_excerpt_tracking.png
:alt: Volumes and lifetimes of tracked GLUT1 site cavities over 100 ns
:width: 90%

Top: volume of each identity (log scale); the six most persistent in colour.
Bottom: one lane per identity with births (triangles), deaths (crosses),
splits (diamonds) and merges (circles).
```

One identity (1, blue) is present in all 21 frames: the large cavity that
contains the site, 1930–3930 Å³ (mean 2740 Å³), lined in every frame by all
five landmark residues; its centroid sits 7–12 Å below the site, toward −*z*,
matching the wide −*z* routes above. Most other identities are small
(5–200 Å³) pockets that appear for one or a few frames. The many splits and
merges of identity 1 are small pieces budding off and rejoining it between
frames 5 ns apart; they summarise everything that happened in between and are
sensitive to the probe and the link thresholds.

## What to check next

- Rerun `paths` with `--spacing 0.35` or `--probe-radius 0.6` to see whether
  the narrow +*z* routes and the 25 ns gap are lattice effects.
- Rerun `track-cavities` with `--link-jaccard 0.1 --link-containment 0`: event
  counts that change a lot are threshold-driven.
- The full 1001-frame trajectory runs in about an hour with
  `--workers 4`, in constant memory.

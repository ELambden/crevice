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
- The probe is chosen automatically for each side of each frame (the default
  `--probe-radius auto`); the command prints the choice, for example
  `paths_f10 probe=auto chosen: + 0.65 Å, - 0.5 Å`.

Each frame takes about 45 s with 4 worker processes (about 25 s with a fixed
`--probe-radius 0.8`) and writes `glut1_f<N>_paths.*`, including
`glut1_f<N>_path_probes.csv` with every candidate probe. What we found
(0.5 Å lattice):

| Frame (time) | −*z* probe (Å) | −*z* widest bottleneck (Å) | −*z* length (Å) | +*z* probe (Å) | +*z* widest bottleneck (Å) | +*z* length (Å) |
|---|---:|---:|---:|---:|---:|---:|
| 0 (0 ns) | 0.9 | 2.15 | 28.1 | 0.5 | 1.18 | 42.0 |
| 5 (25 ns) | 0.55 | 2.13 | 19.3 | 0.5 | 0.77 (lined partly by lipids) | 52.6 |
| 10 (50 ns) | 0.5 | 2.62 | 19.4 | 0.65 | 0.91 | 37.0 |
| 15 (75 ns) | 0.5 | 2.89 | 17.1 | 0.5 | 0.93 (lined by lipids) | 79.4 |
| 20 (100 ns) | 0.5 | 3.05 | 17.3 | 0.5 | 0.97 | 46.4 |

```{figure} images/glut1_paths_bottleneck_by_frame.png
:alt: Widest-route bottleneck radius on each side of the GLUT1 site in all 21 frames, automatic and fixed probe
:width: 95%

Widest-route bottleneck on each side of the site in all 21 frames, with the
automatic probe (orange hollow squares) and a fixed 0.8 Å probe (blue dots;
missing at 25 ns on the +*z* side, where that probe finds no route). Dotted
grey: the probe chosen automatically.
```

On the −*z* side every frame has a wide route (bottleneck 2.0–3.1 Å; Met420
lines it in every frame, usually with Ile168 and Phe291) that reaches the
solvent at *z* ≈ 45–49 Å, a few Å inside the lower phosphate plane, after
14–28 Å. On the +*z* side the routes are narrow (0.8–1.2 Å) and long; at
25 ns they are narrower than 0.8 Å, so a fixed 0.8 Å probe finds none, and
at 25 and 75 ns the routes squeeze between lipids rather than only through
the protein. Read the +*z* routes as narrow geometric gaps, not as an open
pathway, and check them with a finer lattice (`--spacing 0.35`) before
drawing conclusions.

Compared with a fixed 0.8 Å probe on all 21 frames, the automatic probe gives
a widest −*z* route with a bottleneck within 0.1 Å of the fixed-probe one in
every frame (identical in 10 of 21), finds
the +*z* routes at 25 ns, and leaves out three narrower secondary −*z* routes
(1.0–1.9 Å, in frames 0, 8 and 12) whose route set changed from one probe to
the next; `glut1_f<N>_path_probes.csv` lists them, and `--probe-radius 0.8`
reports them.

```{figure} images/glut1_f10_paths.png
:alt: Free-sphere radius against arc length for the five GLUT1 branches at 50 ns
:width: 85%

`glut1_f10_paths.png` (50 ns): one line per branch, junctions as hollow
circles, bottlenecks as dots, dotted lines at the two chosen probes (0.5 Å
for −*z*, 0.65 Å for +*z*). Branch 1 leaves to −*z*; branches 2–5 share their
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

- Rerun `paths` with `--spacing 0.35` to see whether the narrow +*z* routes
  are lattice effects, and compare a fixed `--probe-radius` with the
  automatic choice in `*_path_probes.csv`.
- Rerun `track-cavities` with `--link-jaccard 0.1 --link-containment 0`: event
  counts that change a lot are threshold-driven.
- The full 1001-frame trajectory runs in about an hour with
  `--workers 4`, in constant memory.

# Explore results interactively

CREVICE's main 3D outputs are scenes for PyMOL, VMD and ChimeraX. To look at
a result without leaving Python, or to send it to someone who has none of
those programs, open an **interactive view** instead: the protein as a cartoon,
the cast surfaces, the lining residues, the profile centre line and the
member waters, each with its own checkbox. It works in a Jupyter notebook and
as a single HTML file that opens in any modern web browser, with no network
connection needed.

An interactive view is an extra output that you ask for. No analysis command
writes HTML by itself.

## Before you start

The standalone HTML file needs nothing beyond CREVICE. Showing a view inside
Jupyter also needs [py3Dmol](https://pypi.org/project/py3Dmol/), which the
`view` extra installs:

```bash
python -m pip install "crevice[view]"
```

The examples below use the output directories of two earlier tutorials:
`3UKM_bundle` from [3UKM: a capped channel with side exits](3ukm-lateral-exits.md)
and `glut1-run` from
[A membrane MD trajectory](glut1-trajectory.md). Any `publish`, `cast`
or `hydration` output directory works the same way, and so does a
`region-trajectory` or `cavity-trajectory` directory that has the hydration
scene (written with `--water-membership`, `--outputs default,hydration` or
`--all`; the default trajectory bundle is data and figures only).

## A standalone HTML viewer

```bash
crevice view 3UKM_bundle --open
```

```text
view file=3UKM_bundle/3UKM_view.html layers=protein,entry_1,lumen,exit_1,exit_2,lining_entry_1,lining_lumen,lining_exit_1,lining_exit_2,centre_line size_MB=2.60
```

`--open` opens the file in your default browser; `-o FILE` writes it somewhere
else. Drag to rotate, scroll to zoom and Ctrl-drag (or middle-drag) to move
the scene. Each checkbox shows or hides one layer, *Reset view* returns to the
starting orientation and *Save PNG* saves the canvas as an image.

```{figure} images/3ukm_interactive_view.png
:alt: Interactive view of 3UKM in a web browser, with the yellow entry, teal lumen and two rust-coloured exits inside a grey cartoon, and a column of layer checkboxes
:width: 90%

`3UKM_view.html` in Firefox. The channel cast is split into its ENTRY
(yellow), LUMEN (teal) and two EXIT segments (rust). The starting orientation
is the same as the native scenes: channel axis vertical.
```

The file holds everything it needs, including the
[3Dmol.js](https://3dmol.org) library that draws the scene, so you can e-mail
it or put it on a web page; it never fetches anything over the network. The
3UKM file is 2.6 MB. `--cdn` writes a smaller file that loads 3Dmol.js from
the jsDelivr CDN instead, but then the page needs an internet connection to
open.

To change which layers are visible when the page opens, use `--show` and
`--hide` with a layer name, or with `entry`, `exit` or `lining` for every
layer of that kind:

```bash
crevice view 3UKM_bundle -o 3UKM_lining.html --show lining --show centre_line --hide entry
```

## The same view in Jupyter

```python
import crevice

v = crevice.view("3UKM_bundle")
v                      # the last value in a cell is displayed
```

The notebook shows the same viewer and checkboxes as the HTML page.
`v.layer_names` lists the layers, and `toggle` changes them before the view is
displayed (each call returns the view, so calls can be chained):

```python
v.toggle("lining_lumen").toggle("entry_1", False).show()
v.save_html("3UKM_view.html")
```

`show()` loads 3Dmol.js the way py3Dmol normally does, from a CDN.
`show(offline=True)` puts the library into the notebook instead, which then
works without a network connection but makes the saved notebook about 0.7 MB
larger for every view.

## Trajectory regions and member waters

For a `region-trajectory` run with `--water-membership`, the view shows the
region cast (pale grey), the mean water density (teal) and the waters inside
the measured region in one saved frame (orange spheres):

```bash
crevice view glut1-run --show lining_region
```

```{figure} images/glut1_interactive_view.png
:alt: Interactive view of the GLUT1 candidate sugar site, with orange member waters among grey residue sticks inside the protein cartoon
:width: 90%

`glut1_site_candidate_view.html` with the region and density surfaces switched
off and the region's boundary residues switched on: the 24 member waters of
the first frame.
```

A run saves the first, middle and last analysed frames. `--snapshot 1` (or
`crevice.view("glut1-run", snapshot=1)`) shows the middle one with its own
member waters, and `--snapshot -1` the last.

## Results you have in memory

`crevice.view` also takes results straight from the Python API, so you can
look at a cast before writing any files:

```python
import crevice
from crevice.cast_segments import lining_residues, segment_casts

frame = crevice.load_structure("1GRM.cif")
profile = crevice.pore_profile(frame)
cast = crevice.void_cast(frame, profile=profile)
segments = segment_casts(frame, profile, cast)
v = crevice.view(frame=frame, cast=cast, segments=segments, profile=profile,
                 lining=lining_residues(frame, cast, segments))
```

The surfaces are computed from the cast you pass, on its own grid and without
the short display continuation that `publish` adds beyond each mouth, so they
can look coarser than those of a `publish` bundle (which uses a 0.5 Å grid).

## What you are looking at

- **Protein:** cartoon, grey-blue, 25% transparent.
- **Cast surfaces:** ENTRY yellow, LUMEN teal and EXIT rust, as in the native
  scenes. A cast that is not a channel is violet. Each surface is the native
  display surface: the 0.5 isosurface of the same display field CREVICE writes
  for PyMOL, VMD and ChimeraX (`*_display.dx`), smoothed by `--smooth`.
  Measured volumes are those of the binary `*_volume.dx` grids and are in the
  bundle's tables, not in the picture.
- **Lining residues** (hidden at first): sticks with carbons in their
  segment's colour.
- **Profile centre line** (hidden at first): blue tubes through the axial
  profile samples and any side-exit legs.
- **Member waters:** orange spheres. Membership is geometric (the oxygen
  centre lies inside that frame's measured region); it is not evidence of
  binding.
- **Region outputs:** the region cast is drawn 35% opaque and the water
  density 60% opaque, more opaque than in the native scenes (18% and 42%), so
  that they remain visible through the cartoon in a web browser.

Nothing is written into the 3D picture: no labels, titles or legends. The
layer names appear only as checkbox labels beside it, so *Save PNG* gives a
clean image.

:::{admonition} Interpreting the results
:class: crevice-interpret

An interactive view shows the same display surfaces as the native scenes, so
it is a way to look and to share, not a new measurement: read volumes and
radii from the tables. The region and density surfaces are drawn more opaque
than in the native scenes, and the member waters are one saved frame's
geometric assignment, not evidence of binding.
:::

[Interactive 3D views](../tools/view.md) summarises the options, and the
[method note](../methods/interactive-view.md) defines what is drawn. For the
same channel as a set of curves along its length, try a
[pore fingerprint](pore-fingerprint.md).

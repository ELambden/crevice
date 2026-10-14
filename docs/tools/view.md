# Interactive 3D views

`crevice view` turns an output directory into an interactive 3D scene: the
protein as a cartoon, the cast surfaces, the lining residues, the profile
centre line and, for trajectory regions, the waters inside the region, each
with its own checkbox. It is written as one HTML file that opens in any modern
web browser without a network connection, and {func}`crevice.view` shows the
same scene in a Jupyter notebook.

It doesn't measure anything. It reads the files another command already wrote,
so it takes a second or two, and no analysis command writes HTML unless you
ask for it.

## Quick start

```bash
crevice publish 3UKM --assembly 1 --exclude-hetero --lateral-exits --enclosure-radius 0.55 \
    --out-dir 3UKM_bundle --prefix 3UKM
crevice view 3UKM_bundle --open
```

```python
import crevice

v = crevice.view("3UKM_bundle")      # needs the view extra: pip install "crevice[view]"
v.toggle("lining_lumen").show()
v.save_html("3UKM_view.html")
```

The HTML file needs nothing beyond CREVICE; showing a view inside Jupyter also
needs py3Dmol, which the `view` extra installs.

```{figure} ../examples/images/3ukm_interactive_view.png
:alt: Interactive view of 3UKM in a web browser, with the yellow entry, teal lumen and two rust-coloured exits inside a grey cartoon, and a column of layer checkboxes
:width: 90%

`crevice view 3UKM_bundle` in Firefox: ENTRY (yellow), LUMEN (teal) and the two
side EXITs (rust) inside the protein cartoon.
```

## How it works

`crevice view` reads the manifest or scene JSON of a `publish`, `cast`,
`region-trajectory`, `cavity-trajectory` or `hydration` run and builds one
layer per object:

- **Cast surfaces** are the display surfaces of the native PyMOL, VMD and
  ChimeraX scenes: the 0.5 isosurface of the same `*_display.dx` field,
  including any `--smooth` you chose. Volumes are measured on the binary
  `*_volume.dx` grids and live in the tables, not in the picture.
- **Lining residues** and the **profile centre line** start hidden; switch them
  on with a checkbox, `--show` or `.toggle()`.
- **Member waters** come from one of the frames a `region-trajectory
  --water-membership` run saved (`--snapshot`: 0 the first, 1 the middle, -1
  the last).

The page includes the 3Dmol.js library (Rego and Koes 2015), so it works
offline and can be e-mailed or put on a web page; `--cdn` writes a smaller file
that loads the library from a CDN instead. The
[method note](../methods/interactive-view.md) gives the full definitions.

You can also view results you still have in memory, before writing any files:
`crevice.view(frame=..., cast=..., segments=..., profile=..., lining=...)`.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `-o/--output` | `OUTDIR/PREFIX_view.html` | where to write the HTML file |
| `--open` | off | open the file in your web browser |
| `--show LAYER`, `--hide LAYER` | | change which layers are visible at first; `entry`, `exit` and `lining` select every layer of that kind (repeatable) |
| `--snapshot N` | 0 | which saved frame's member waters to show |
| `--kind` | `auto` | show the cast scene or the hydration (region) scene of a directory that has both |
| `--cdn` | off | load 3Dmol.js from a CDN instead of inlining it |

Every option is listed under [`crevice view`](../reference/cli/view.md).

## What you get

- **`PREFIX_view.html`**: one self-contained page, about 2.6 MB for the 3UKM
  bundle (0.54 MB of it is 3Dmol.js). Drag to rotate, scroll to zoom,
  Ctrl-drag to move; *Reset view* returns to the starting orientation (channel
  axis vertical, as in the native scenes) and *Save PNG* saves the canvas.
- In Python, a {class}`~crevice.viewer.CreviceView` with `layer_names`,
  `toggle()`, `show()` and `save_html()`.

Nothing is drawn into the 3D picture itself: no labels, titles or legends, so
a saved PNG is clean.

:::{admonition} Interpreting results
:class: crevice-interpret

An interactive view is for looking, sharing and checking; it shows the same
display surfaces as the native scenes, so it inherits their choices
(isosurface level, smoothing, the short display continuation beyond each
mouth). Read numbers from the tables, not from the picture. In region views
the cast and density surfaces are drawn more opaque than in the native scenes
so they stay visible in a browser. Member waters are a geometric assignment
(the oxygen lies inside that frame's measured region), not evidence of
binding. For publication figures, the native ChimeraX or PyMOL scenes give you
more control over lighting and resolution.
:::

Try it in the [Explore results interactively](../examples/explore-interactively.md)
tutorial. Profiles and casts come from [`profile`](profile.md),
[`cast`](cast.md) and [`publish`](publish.md); a
[pore fingerprint](fingerprint.md) gives the same channel as a set of curves.

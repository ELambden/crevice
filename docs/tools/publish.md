# One-command output bundle

`crevice publish` is the quickest way to get everything for one structure. It
runs the profile, cast, lining residues, cavities, tunnels, residue network
and hydration, and writes the tables, figures and PyMOL, VMD and ChimeraX
scenes into one directory, with a manifest listing every file and setting. It
doesn't measure anything new; it runs the other tools for you, consistently.

## Quick start

```bash
crevice publish 1GRM --out-dir results/1GRM --prefix 1GRM
```

```python
from crevice import load_structure, write_static_publication_bundle

frame = load_structure(".crevice/pdb/1GRM.cif")
manifest = write_static_publication_bundle(
    frame, structure_path=".crevice/pdb/1GRM.cif", output_dir="results/1GRM",
    prefix="1GRM", samples=81)
```

The command also adds hydration outputs and an input report, which the Python
call leaves out. For finer detail on a channel, try
`--section-spacing 0.25 --cast-spacing 0.25 --cast-max-grid-points 2000000`.

```{figure} ../examples/images/3ukm_cast_exits_chimerax.png
:alt: The TWIK-1 channel cast in ChimeraX with yellow entry, teal lumen and rust exit segments.
:width: 420px

A `publish` scene of TWIK-1 (3UKM) with `--lateral-exits`, rendered in
ChimeraX: the entry vestibule in yellow, the lumen in teal and the two side
exits in rust.
```

## What it runs

1. A [profile](profile.md), with any profile options you pass.
2. A cast. With `--cast-mode auto` (the default) this is a channel cast when
   the profile resolves. If it doesn't, `publish` says so, records why and
   casts the main cavity with a [rolling probe](cast.md) instead. `channel`,
   `cavity`, `all` and `rolling` force a mode.
3. [Lining residues](residues.md), [cavities](cavities.md),
   [tunnels](tunnels.md) and a [network](network.md), unless you pass
   `--skip-cavities`, `--skip-tunnels` or `--skip-network`.
4. [Hydration](hydration.md), unless you pass `--skip-hydration`.

A channel cast fills all the free space in each slice, not just the largest
sphere. For display, it extends 2 Å beyond each mouth (`--cast-extension`);
that extension is never counted in volumes or radii.

## Entry, lumen and exit

When a channel resolves, its cast is split into separate pieces, each with its
own volume and its own object in the viewer:

| Segment | What it is | Colour |
|---|---|---|
| **LUMEN** (`crevice_lumen`) | the channel between its two mouths; this is the cast volume reported everywhere else | teal `#1494a1` |
| **ENTRY** (`crevice_entry_1`, …) | the vestibule beyond one mouth, out to bulk solvent | yellow `#f0d43a` |
| **EXIT** (`crevice_exit_1`, …) | the vestibule beyond the other mouth, or each side exit of a capped channel | rust `#c2410c` |

Each vestibule fills all the probe-accessible space attached to its centre
line, out to where bulk solvent begins. Bulk is the space a rolling probe of
`--exit-bulk-radius` (6 Å) reaches from outside, so that option sets where a
vestibule *ends*, not how wide it may be.

Which end is the entry? With `--entry-end auto`, if exactly one end is capped
with side exits, those are the exits and the open end is the entry; otherwise
the start of the axis is the entry. `--entry-end start|end` overrides this.
The names are geometric labels, not a statement about which way anything
moves.

Segments are reported in `PREFIX_cast_segments.csv` and in the cast JSON and
manifest, each with its own binary map (`PREFIX_entry_1_cast.dx`, …).
`total_volume` stays the lumen alone; `segments_total_volume_A3` is the
separate sum. A rolling-probe cast, including the fallback, has no axis or
mouths and isn't split.

### Lining residues

A residue lines a segment when any of its atoms comes within `--lining-cutoff`
(3.3 Å, a hydrogen-bond distance) of the segment's measured grid points.
They're written to `PREFIX_lining_residues.csv`, flagged in the residue and
network tables (`lumen_lining`, `lining_segments`) and added to the scenes as
stick objects `crevice_<segment>_lining`, hidden until you switch them on:

```text
crevice_lining lumen, on                      # PyMOL: on, off or toggle
crevice_lining lumen on                       # VMD
color magenta, crevice_lumen_lining_sel       # PyMOL selection
show crevice_lumen_lining atoms               # ChimeraX named selection
```

Every object can be recoloured or hidden on its own, for example
`color red, crevice_exit_1` or `disable crevice_entry_1` in PyMOL. The four
cast colours (with violet for non-channel casts) were checked to stay
distinguishable for normal and red-green colour vision.

### Capped channels

With `--lateral-exits`, a capped channel resolves through its side exits (see
[profile](profile.md#capped-channels-and-side-exits)). The lumen cast still
covers only the axial channel; each exit leg gets its own cast, reported
separately and drawn as an EXIT segment, so you can see it widen towards its
portal. `--exit-centre-lines` adds a thin blue tube along each leg.

## When no channel resolves

If the profile doesn't resolve, `publish` still runs everything else on the
fallback cast. No profile is invented: there is no profile JSON, CSV or radius
plot. The tunnel search starts from the widest point of the cast, the network
treats the whole cast as its region node, and
`<prefix>_profile_status.json` records `status: "unresolved"` with the reason.
The manifest marks each extra as `completed`, `skipped_by_request` or
`not_completed`, and the summary line ends with `profile=unresolved`. Requests
that only make sense for a channel (focus, maximum volume, minimum clearance,
hydrogens) fail instead of switching silently.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--cast-mode` | `auto` | `auto`, `channel`, `cavity`, `all` or `rolling` |
| `--cast-spacing` | 0.5 Å | cast grid spacing |
| `--cast-max-cavities` | 1 | regions in a rolling cast |
| `--focus X Y Z`, `--focus-radius` | none, 8 Å | crop the cast to spheres around points (repeatable) |
| `--entry-end` | `auto` | which end is the entry |
| `--lining-cutoff` | 3.3 Å | distance for lining residues |
| `--skip-cavities`, `--skip-tunnels`, `--skip-network`, `--skip-hydration` | off | leave parts out |

Profile settings, atomic radii, figure text, mouth guides and smoothing are
[shared options](../reference/cli/common-options.md); every option is listed
under [`crevice publish`](../reference/cli/publish.md).

## What you get

A default 1GRM run writes about 85 files:

| File | Contents |
|---|---|
| `*_profile.csv` | the profile, one row per point |
| `*_residue_contacts.csv` | residues along the channel, with roles |
| `*_lining_residues.csv` | residues lining each segment |
| `*_void_cast.csv`, `*_cast_segments.csv` | cast regions and segment volumes |
| `*_cavities.csv`, `*_tunnels.csv`, `*_tunnel_points.csv` | grid cavities and tunnels |
| `*_network_nodes.csv`, `*_network_edges.csv` | the residue network |
| hydration tables | see [Hydration](hydration.md#what-you-get) |

Plus the measured `*_volume.dx` and display maps, viewer scenes and render
scripts, the radius, residue, cavity and network figures, dummy-atom PDB files
for the cast, cavities and tunnels, and JSON records for the manifest, the
input report and the full results. Figures and scenes are clean unless you
add `--annotate`.

:::{admonition} Interpreting results
:class: crevice-interpret

A bundle is only as good as the settings behind it, and they're all in the
manifest; report them. `publish` keeps heteroatoms as walls unless you pass
`--exclude-hetero`, unlike `cast`, and its network uses the Python defaults
(surface gap, no heteroatoms) rather than those of `crevice network`. An
unresolved profile doesn't show the protein is closed: in our runs 4PYP is
unresolved at every probe, and its fallback cast still measures a large
cavity. Segment labels are geometry, not transport direction. Open the scenes
and check them before using them in a figure.
:::

Worked bundles: the [1GRM](../examples/1grm-profile-cast.md) and
[3UKM](../examples/3ukm-lateral-exits.md) tutorials.

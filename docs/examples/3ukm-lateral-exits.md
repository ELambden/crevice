# 3UKM: a capped channel with lateral exits

[PDB 3UKM](https://www.rcsb.org/structure/3UKM) is the two-pore-domain
potassium channel TWIK-1 (K2P1). Its pore is capped on the extracellular side
by a domain that blocks the straight axial exit; the structure has two side
portals, related by the dimer's two-fold axis, through which ions are thought
to leave. By default `profile`
requires both ends to continue straight out of the protein, so it refuses this
channel. `--lateral-exits` accepts a capped end when a free path leads from it
sideways to bulk solvent, places that end where the widest such path leaves
the lumen, and reports every distinct exit. The example takes a few minutes.

## The biological assembly, without ions and ligands

```bash
crevice fetch 3UKM --assembly 1 --cache-dir .crevice/pdb
```

```text
.crevice/pdb/3UKM-assembly1.cif	download	409063 bytes	unverified (no accession recorded: _entry.id holds the placeholder XXXX)
```

The deposited entry holds two channel dimers; RCSB's biological assembly 1
(chains A and B) is one dimer. Assembly files carry a placeholder entry ID, so
the accession cannot be checked against the file; the provenance sidecar
records that, with the SHA-256 of what was downloaded. Any command can read the
assembly with `--assembly 1`.

## Profile with lateral exits

```bash
crevice profile 3UKM --assembly 1 --exclude-hetero --lateral-exits \
    -o 3UKM_profile.json --png 3UKM_profile.png --dpi 150
```

```text
profile points=81 min_radius=0.641 mean_radius=3.141 bottleneck_index=22
enclosure_radius=per-plane 0.30-3.00 (auto: per-plane enclosure probe, 0.30-3.00 A over 81 samples)
```

The last status line on standard error says what kind of path was found:

```text
crevice: pore radius found in 4.6 s: minimum 0.641 A at s = -1.4 A (no lower than 0.621 A between samples); mean 3.141 A over 81 samples, s = -12.8 to +28.6 A; per-plane enclosure probes 0.30-3.00 A; capped end with a lateral exit
```

`--exclude-hetero` removes the HETATM records (5 K⁺ ions, 22 atoms of the
undecane ligand and 4 waters); the potassium ions in the filter would otherwise
block the channel. The narrowest axial sample is 0.641 Å (nearest residue
B:THR118), in the selectivity filter. Each slice chose its own enclosure
probe: 0.30 Å through the filter, and up to 3.0 Å in the wide inner cavity
and vestibule at the open (upper) end.

`metadata.path_type` is `capped_channel_lateral_exit`. In `metadata.exits`,
the upper end is `axial` (a straight exit) and the lower end is `lateral`.
Which end is "lower" depends on the sign of the axis the automatic search
chose. Below the filter the lumen widens into the chamber under the cap
(radius up to 3.5 Å near t = −12.2 Å). The chamber's sideways openings are
narrower than the chamber itself, so every slice there is enclosed, and the
slice-by-slice trace alone continues past them into the cap, which narrows to
1.2 Å and leaves through a single 1.09 Å route
(`metadata.channel_mouths.lower.traced_mouth`, t = −17.9 Å).

A capped end leaves sideways in any case, so `--lateral-exits` places it where
the widest route to bulk leaves the lumen. From the slices between the filter
and the traced mouth, CREVICE searches for a route to bulk that is clearly
wider than the rest of the axial path (1.09 Å here: the cap and its exit) and
that therefore cannot continue through the lumen. The widest, 2.20 Å, leaves
from the slice at t = −12.75 Å, in the chamber under the cap, so the lumen ends
there (mouth reason `lateral_branch`; the decision is recorded in
`metadata.exits.lower.lateral_branch`), and the side-exit search from that
mouth finds the two portals (`leg_count` 2). Each leg in
`metadata.exits.lower.legs` records its samples and free-sphere radii, its
bottleneck, length and exit point:

| Leg | Bottleneck radius (Å) | Length (Å) | Exit point offset from the mouth (Å): axial, lateral |
|---:|---:|---:|---|
| 1 | 2.20 | 16.0 | −9.5, 10.1 |
| 2 | 2.17 | 15.1 | −9.5, 9.4 |

Legs are ranked by their width beyond the mouth sphere
(`ranking_width_raw_clearance_A`, 2.161 and 2.118 Å), then by length. Leg
radii are distances to the nearest atom surface along a 3D path, like tunnel
radii, not planar section radii, so they never enter the axial `min_radius`.

```{figure} images/3ukm_profile.png
:alt: 3UKM pore radius against the axial coordinate with the default per-slice probe, with two lateral exit legs
:width: 85%

`3UKM_profile.png`. Blue line and fill: axial pore radius (Å) between the two
mouths (light grey outside them). Red dashed line and dot: the narrowest sample
(0.641 Å). Bluish-green dash-dot lines: the two lateral exit legs, continuing
from the capped (lower) mouth; beyond that mouth the horizontal axis is the
mouth's axial coordinate minus the path length along the leg, not an axial
position, as the legend says.
```

## The same portals with a fixed probe

A fixed probe smaller than the portals but still able to pass the filter
stops the lumen where the chamber first leaks, and the side-exit search from
there finds the same two portals:

```bash
crevice profile 3UKM --assembly 1 --exclude-hetero --lateral-exits --enclosure-radius 0.55 \
    -o 3UKM_profile_0.55.json --png 3UKM_profile_0.55.png --dpi 150
```

```text
profile points=81 min_radius=0.635 mean_radius=1.285 bottleneck_index=48
```

The lumen now runs from t = −12.6 to +5.8 Å: below the filter it ends where
the chamber leaks sideways at 0.55 Å (0.2 Å from where the default probe's
lumen ends), and above it where the inner cavity is too wide for the probe to
enclose. The lower end has `leg_count` 2, one leg through each portal:

| Leg | Bottleneck radius (Å) | Length (Å) | Exit point offset from the mouth (Å): axial, lateral |
|---:|---:|---:|---|
| 1 | 2.17 | 16.2 | −10.0, 10.1 |
| 2 | 2.22 | 15.8 | −10.0, 9.7 |

Here the ranking widths are 2.170 and 2.166 Å.

```{figure} images/3ukm_profile_fixed055.png
:alt: 3UKM pore radius against the axial coordinate with a fixed 0.55 Å enclosure probe, with two lateral exit legs
:width: 85%

`3UKM_profile_0.55.png`, drawn as above: the shorter lumen of the fixed
0.55 Å probe and its two exit legs (bluish-green dash-dot lines).
```

Both probes give the two portals with bottlenecks within 0.05 Å of each
other. The default keeps the wide upper vestibule in the lumen; the fixed
0.55 Å probe ends the lumen at the inner cavity.

## The cast and the exit legs in 3D

```bash
crevice report 3UKM --assembly 1 --exclude-hetero --lateral-exits \
    --out-dir 3UKM_bundle --prefix 3UKM --dpi 150
```

```text
publish files=41 out_dir=3UKM_bundle profile=resolved
```

The bundle uses the default per-slice probe, so its profile is the first one
on this page: the lumen runs from the portals' branching plane to the open
upper end, and its cast fills the wide upper vestibule as well as the pore
(with the fixed 0.55 Å probe the lumen stops at the inner cavity and the
vestibule becomes a separate entry segment). The command also warns on
standard error that the narrowest clearance (0.64 Å) is below two grid
spacings of the 0.5 Å cast grid, so narrow openings and volume may be lost
(`--cast-spacing 0.35` or finer for narrow channels). The default outputs
keep the bundle to the profile, the cast, the residue contacts and the typed
residue interactions (`--all` adds the network and the rest). The manifest records `profile_status.path_type`,
`exit_types` and `exit_leg_counts`.

The cast is split into LUMEN and EXIT segments (and an ENTRY segment when
there is space beyond an open axial end before bulk). The LUMEN is the axial
channel cast between the mouths. Each exit leg is cast: the grid points beyond
the capped mouth plane and outside bulk solvent that lie inside the free
spheres of the leg's samples, swept by the enclosure probe and connected to
the leg's centre line. The open axial end gets the same treatment along a
straight centre line continuing the axis. Bulk (a rolling 6 Å probe,
`--exit-bulk-radius`) is the outer boundary of every segment, not a width
limit: up to that boundary each segment fills the probe-accessible space. Because one end is
capped with lateral exits, those legs are the EXITs; the open upper end meets
bulk at its mouth, so it has no ENTRY segment here (`--entry-end` sets which
end is the entry; the labels are geometric, not a transport direction).
Volumes are reported separately and never added to the axial cast:

| Segment | Source | Volume (Å³, 0.5 Å grid) | Lining residues (3.3 Å) |
|---|---|---:|---:|
| `lumen` | axial channel cast (`total_volume`), pore and upper vestibule | 5801.3 | 156 |
| `exit_1` (lower end, leg 1) | `lateral_exit` | 725.4 | 20 |
| `exit_2` (lower end, leg 2) | `lateral_exit` | 636.5 | 18 |

The lumen-lining residues include the filter and pore-helix residues
S116–V124 and T225–V232 of both subunits, with N87–G89 and W100–D103 of the
cap, residues of the chamber under it (E207, W210, P233–G236, K246) and the
wall of the upper vestibule (among others L146, T150, L257–L264, C268 and
R277 of each subunit). The exit legs are lined mostly by cap residues
(G80–L93, N101) of both subunits, with H122, D209 and E235–Y237. The full
table is `3UKM_lining_residues.csv`
(segment, residue, closest atom, minimum distance), and
`3UKM_residue_contacts.csv` gains `lumen_lining` and `lining_segments`
columns. The segment volumes are in `3UKM_cast_segments.csv` and
`metadata.cast_segments`; each exit segment also has a measured map
(`3UKM_exit_1_cast.dx`, `3UKM_exit_2_cast.dx`).

They are in `metadata.lateral_exit_casts` of `3UKM_void_cast.json` and the
manifest, as `lateral_exit_lower_1`/`_2` rows of `3UKM_void_cast.csv`, in
`profile_status.lateral_exit_cast_volumes_A3`, and as the binary map
`3UKM_exit_casts.dx`. A segment volume depends on the bulk definition
(`--exit-bulk-radius`) as well as the grid, so treat it as a description of
the portal region under these settings, not as a measured permeation volume.
The scene file `3UKM_scene.json` still holds one `exit_paths` centre line per
leg; `--exit-centre-lines` draws them as thin blue tubes.

```{figure} images/3ukm_cast_exits_chimerax.png
:alt: ChimeraX render of the 3UKM channel lumen with its upper vestibule and the two exit-leg casts as separate objects
:width: 65%

`3UKM_volume.cxc` rendered by ChimeraX 1.12 (headless,
`scripts/check_chimerax.py --xvfb`). Teal: the LUMEN (`crevice_lumen`,
5801.3 Å³ on the default 0.5 Å grid), the narrow pore with the wide upper
vestibule. Rust: the two EXIT legs (`crevice_exit_1`, `crevice_exit_2`), which
widen towards the portals in the cap at the bottom. Each is a separate model
with its own colour. Grey-blue cartoon: the protein. No mouth rings, centre
lines, lining sticks or labels are drawn by default. The check confirmed the
surfaces (within 7 × 10⁻⁶ Å of the reference meshes) and the camera. The PyMOL and VMD versions of this scene were not opened for this
page.
```

:::{admonition} Interpreting the results
:class: crevice-interpret

These values come from one run of the commands above on the downloaded
assembly. The axis was found automatically. Where the capped end is placed
follows a geometric rule (the widest route to bulk, when it is clearly wider
than the rest of the axial path); the route through the cap is still there,
only narrower. The segment volumes depend on the probe, the grids and the
bulk definition (`--exit-bulk-radius 6`, `--exit-spacing 0.5`). The portals
are geometric results with these settings, not demonstrated permeation
pathways.
:::

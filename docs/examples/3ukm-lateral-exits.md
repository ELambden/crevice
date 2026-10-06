# 3UKM: a capped channel with lateral exits

[PDB 3UKM](https://www.rcsb.org/structure/3UKM) is the two-pore-domain
potassium channel TWIK-1 (K2P1). Its pore is capped on the extracellular side
by a domain that blocks the straight axial exit; ions leave sideways through
two portals related by the dimer's two-fold axis. By default `profile`
requires both ends to continue straight out of the protein, so it refuses this
channel. `--lateral-exits` accepts a capped end when a free path leads from it
sideways to bulk solvent, and reports every distinct exit. The example takes
about ten minutes, most of it for the automatic enclosure probe.

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
profile points=81 min_radius=0.618 mean_radius=1.031 bottleneck_index=44
enclosure_radius=0.550 (auto: smallest probe of channel group 1, which resolved at 2 of 21 probes tried)
```

`--exclude-hetero` removes the HETATM records (5 K⁺ ions, 22 atoms of the
undecane ligand and 4 waters); the potassium ions in the filter would otherwise
block the channel. The automatic enclosure probe found two channels that each
resolved at two probes (0.55–0.60 Å and 1.6–1.8 Å); the tie goes to the channel
found at the smaller probe, and its smallest probe, 0.55 Å, is used. The window
of probes that resolve the narrow filter is only about 0.07 Å wide. The
narrowest axial sample is 0.618 Å (nearest residue A:THR118).

`metadata.path_type` is `capped_channel_lateral_exit`. In
`metadata.exits`, the upper end is `axial` (a straight exit) and the lower end
is `lateral` with `leg_count` 2: two distinct exits, one through each portal.
Which end is "lower" depends on the sign of the axis the automatic search
chose. Each leg in `metadata.exits.lower.legs` records its samples and
free-sphere radii, its bottleneck, length and exit point:

| Leg | Bottleneck radius (Å) | Length (Å) | Exit point offset from the mouth (Å): axial, lateral |
|---:|---:|---:|---|
| 1 | 2.03 | 17.0 | −12.0, 10.1 |
| 2 | 2.03 | 16.2 | −10.0, 11.4 |

Both bottlenecks are the mouth sample's own clearance (2.03 Å), which every
route shares; the legs are therefore ranked by their width beyond the mouth
sphere (2.12 and 2.08 Å, `ranking_width_raw_clearance_A`). Leg radii are
distances to the nearest atom surface along a 3D path, like tunnel radii, not
planar section radii, so they never enter the axial `min_radius`.

```{figure} images/3ukm_profile.png
:alt: 3UKM pore radius against the axial coordinate, with two lateral exit legs
:width: 85%

`3UKM_profile.png`. Blue line and fill: axial pore radius (Å) between the two
mouths (light grey outside them). Red dashed line and dot: the narrowest sample
(0.618 Å). Bluish-green dash-dot lines: the two lateral exit legs, continuing
from the capped (lower) mouth; beyond that mouth the horizontal axis is the
mouth's axial coordinate minus the path length along the leg, not an axial
position, as the legend says.
```

## The cast and the exit legs in 3D

```bash
crevice publish 3UKM --assembly 1 --exclude-hetero --lateral-exits --enclosure-radius 0.55 \
    --skip-hydration --skip-cavities --skip-tunnels --skip-network \
    --out-dir 3UKM_bundle --prefix 3UKM --dpi 150
```

```text
publish files=32 out_dir=3UKM_bundle profile=resolved
```

Passing the chosen probe explicitly skips the ladder (a minute instead of
seven); the `--skip-*` options keep the bundle to the profile, the cast and the
residue contacts. The manifest records `profile_status.path_type`,
`exit_types` and `exit_leg_counts`.

Each exit leg is also cast: the grid points beyond the capped mouth plane and
outside bulk solvent that lie inside the free spheres of the leg's samples,
swept by the enclosure probe and connected to the leg's centre line. The leg
casts show how the space widens from the mouth towards each portal. Their
volumes are reported separately and never added to the axial cast:

| Cast | Kind | Volume (Å³, 0.5 Å grid) |
|---|---|---:|
| axial channel between the mouths | `section_channel` | 186.6 |
| lower end, leg 1 | `lateral_exit` | 445.8 |
| lower end, leg 2 | `lateral_exit` | 593.5 |

They are in `metadata.lateral_exit_casts` of `3UKM_void_cast.json` and the
manifest, as `lateral_exit_lower_1`/`_2` rows of `3UKM_void_cast.csv`, in
`profile_status.lateral_exit_cast_volumes_A3`, and as the binary map
`3UKM_exit_casts.dx`. A leg volume depends on the bulk definition
(`--exit-bulk-radius`) as well as the grid, so treat it as a description of
the portal region under these settings, not as a measured permeation volume.
The scene file `3UKM_scene.json` still holds one `exit_paths` centre line per
leg; `--exit-centre-lines` draws them as thin blue tubes.

```{figure} images/3ukm_cast_exits_chimerax.png
:alt: ChimeraX render of the 3UKM channel cast joined to two exit-leg casts
:width: 65%

`3UKM_volume.cxc` rendered by ChimeraX 1.12 (headless,
`scripts/check_chimerax.py --xvfb`). Teal: the cast of the enclosed channel
between its mouths (186.6 Å³ on the default 0.5 Å grid) joined at the capped
mouth to the casts of the two exit legs, which widen towards the portals.
Grey-blue cartoon: the protein. No mouth rings or centre lines are drawn by
default. The check confirmed the surface (within 8 × 10⁻⁶ Å of the reference
mesh) and the camera. The PyMOL version of this scene was also rendered and
checked during development; the VMD version was not opened for this page.
```

## What was checked

The commands were run with CREVICE 0.1.0 on the downloaded assembly and the
values are copied from that run. The axis was chosen automatically, not taken
from a curated reference; the two-portal result is a geometric observation
under these settings (probe, grids, `--exit-bulk-radius 6`,
`--exit-spacing 0.5`), not a validated permeation pathway.

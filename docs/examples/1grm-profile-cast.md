# 1GRM: pore radius profile and channel cast

This example measures the gramicidin A channel,
[PDB 1GRM](https://www.rcsb.org/structure/1GRM) (a solution NMR structure; the
first model is used), with [`profile`](../tools/profile.md), then writes the
complete static bundle with [`publish`](../tools/publish.md): the 3D channel
cast, residue contacts, cavities, tunnels, the residue network and viewer
scenes. It takes about a minute.

## Get the structure

```bash
crevice fetch 1GRM --cache-dir .crevice/pdb
```

```text
.crevice/pdb/1GRM.cif	download	156004 bytes	verified
```

`fetch` downloads the mmCIF file from RCSB
(`https://files.rcsb.org/download/1GRM.cif`), checks that it states the
accession 1GRM, and writes a provenance sidecar
(`1GRM.cif.provenance.json`: URL, SHA-256, size, time). Every later command can
name the accession instead of a path; it reads the cached file (`--offline`
forbids any download). `.crevice/pdb` is also the default cache.

## Radius profile

```bash
crevice profile 1GRM --axis auto -o 1GRM_profile.json \
    --png 1GRM_profile.png --annotated-png 1GRM_profile_residues.png --dpi 150
```

```text
profile points=81 min_radius=1.324 mean_radius=1.563 bottleneck_index=73
enclosure_radius=0.750 (auto: smallest probe of channel group 1, which resolved at 5 of 21 probes tried)
```

`--axis auto` tried the three principal directions and accepted the one with a
connected, laterally enclosed run of sections and two open ends. The enclosure
probe was chosen automatically: the channel resolved at five probes of the
0.5–3.0 Å ladder (0.75–1.1 Å), and the smallest of them, 0.75 Å, was used. The
narrowest sample is 1.324 Å from the nearest atom surface (sample 73, nearest
residue A:TRP13). Files written:

- `1GRM_profile.json`: every sample, the axis, mouths, the probe choice
  (`metadata.enclosure_probe_selection`) and the settings;
- `1GRM_profile.csv`: one row per sample (`index`, `x_A`, `y_A`, `z_A`,
  `radius_A`, `raw_clearance_A`, `axis_position_A`, `nearest_atom_serial`,
  `nearest_residue`), written beside the JSON by default;
- `1GRM_profile.png` and `1GRM_profile_residues.png` (below).

```text
index,x_A,y_A,z_A,radius_A,raw_clearance_A,axis_position_A,nearest_atom_serial,nearest_residue
0,-0.1885,-11.6853,3.3339,1.6850,1.6850,-11.6599,235,B:TRP13
1,-0.2821,-11.3977,3.4458,1.6625,1.6625,-11.3693,272,B:ETA16
```

(values rounded here; the file keeps full precision)

```{figure} images/1grm_profile.png
:alt: Pore radius of 1GRM against the axial coordinate
:width: 85%

`1GRM_profile.png`. Blue line and light-blue fill: pore radius (Å) at each
sample against the axial coordinate (Å). Red dashed line and dot: the
bottleneck (1.324 Å near t = 9.6 Å). Light grey: outside the two channel
mouths, where lateral enclosure is lost (`--mouth-guides` adds orange mouth
bands and lines; the positions are always in the JSON and PNG metadata). The figure
is clean by default; the values are in its PNG metadata, and `--annotate`
draws them (see the [radii and annotation example](custom-radii-annotate.md)).
```

```{figure} images/1grm_profile_residues.png
:alt: 1GRM radius profile above a lane of nearest residues
:width: 85%

`1GRM_profile_residues.png`: the same profile above a lane of the nearest
residue every fifth sample (`--residue-label-every`), coloured by chemistry
(purple charged, teal polar, blue other). Residue names are in the PNG
metadata and are drawn with `--annotate`.
```

The same profile from Python:

```python
from crevice import load_structure, pore_profile, write_profile_csv

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto")          # enclosure probe chosen automatically
print(round(profile.min_radius, 3), profile.bottleneck.nearest_residue)   # 1.324 A:TRP13
print(profile.metadata["enclosure_probe_selection"]["chosen_A"])          # 0.75
write_profile_csv(profile, "1GRM_profile.csv")
```

## Complete bundle with the channel cast

```bash
crevice publish 1GRM --out-dir 1GRM_bundle --prefix 1GRM --skip-hydration --dpi 150
```

```text
publish files=55 out_dir=1GRM_bundle profile=resolved
```

`--skip-hydration` leaves out the water and solvent-accessibility analysis
(1GRM has no explicit water) and its viewer scenes; without it a default run
writes about 85 files. The bundle holds the same 81-sample profile
(`1GRM_profile.csv`), the measured channel cast and its scenes, residue
contacts, cavities, tunnels and the residue network. Every file is listed in
`1GRM_manifest.json`, which also records the chosen enclosure probe under
`profile_status.enclosure_probe`.

The cast is the probe-swept lumen between the two mouths: grid points
(0.5 Å here, `--cast-spacing`) that clear every atom (`--cast-min-radius`,
default 0 Å), swept by the enclosure probe from the enclosed sections and
connected to the measured centre line. `1GRM_void_cast.csv` reports one region
of 263.6 Å³ (2109 samples). The cast volume depends on the cast grid: the
0.25 Å grid used for the published 1GRM figures gives 292.4 Å³. Compare grids
before interpreting a volume.

This cast is the **LUMEN**. Beyond each mouth, the space out to bulk solvent
(the boundary set by a rolling 6 Å probe, `--exit-bulk-radius`) is filled by
the same probe-swept rule and reported as a separate segment
(`1GRM_cast_segments.csv`, `metadata.cast_segments` in `1GRM_void_cast.json`):

| Segment | End | Volume (Å³, 0.5 Å grid) | Volume (Å³, 0.25 Å grid) | Lining residues (3.3 Å) |
|---|---|---:|---:|---|
| `entry_1` | axis start (lower) | 17.6 | 24.2 | B:DLE10–B:ETA16 (7) |
| `lumen` | between the mouths | 263.6 | 292.4 | all 32 residues |
| `exit_1` | axis end (upper) | 24.6 | 27.3 | A:DLE10–A:ETA16 (7) |

Both ends are open, so the default rule labels the axis start ENTRY and the
axis end EXIT; `--entry-end end` swaps them. These are geometric labels, not
a transport direction. The vestibules are small because bulk solvent reaches
almost to the mouths of this short membrane channel. Lining residues (any atom
centre within 3.3 Å of a segment's grid points, `--lining-cutoff`) are in
`1GRM_lining_residues.csv` and flagged in `1GRM_residue_contacts.csv` and
`1GRM_network_nodes.csv` (`lumen_lining`, `lining_segments`).

`1GRM_volume.dx` is the measured binary map from which the volume and
residue attribution come; `1GRM_display.dx` and the PyMOL mesh are the
smoothed display surface (`--smooth`, default 0.4 Å), which changes the picture
only. Open a scene from inside the bundle directory:

```bash
cd 1GRM_bundle
pymol 1GRM_volume.pml        # or: vmd -e 1GRM_volume.vmd, chimerax 1GRM_volume.cxc
```

```{figure} images/1grm_cast_chimerax.png
:alt: ChimeraX render of the 1GRM channel cast inside the protein cartoon
:width: 60%

`1GRM_volume.cxc` rendered by ChimeraX 1.12 (headless, through
`scripts/check_chimerax.py --xvfb`). Teal: the LUMEN (`crevice_lumen`).
Yellow, bottom: the ENTRY vestibule (`crevice_entry_1`); rust, top: the EXIT
vestibule (`crevice_exit_1`). Each is its own model, so it can be recoloured
or hidden alone; separately smoothed surfaces meet at the mouth planes.
Grey-blue, 25% transparent: the protein cartoon. No mouth guides or labels are
drawn by default; the lining-residue stick objects (`crevice_lumen_lining`,
...) are present but hidden. The check confirmed that the displayed surfaces
match the reference meshes within 2.2 × 10⁻⁶ Å and that the camera matches the
scene file. The PyMOL scene was also run and checked (`scripts/check_pymol.py`);
the VMD version was not opened for this page.
```

`1GRM_residue_contacts.csv` ranks the residues around the profile; at the
default 4.5 Å cutoff the top rows are residues whose atoms touch the sampled
region (`min_distance_A` 0, role `bottleneck`), for example A:ALA3 and A:VAL7.
`1GRM_network_chord.png` and `1GRM_network_summary.png` summarise the residue
contact network.

:::{admonition} Interpreting the results
:class: crevice-interpret

The printed values and file counts come from a CREVICE 0.1.0 run (the cast
segments and lining residues from the development version after it). They are
measurements with the settings shown; gramicidin A is a convenient, well-known
test case, not a curated reference, so treat them as such (see
[How far to trust the numbers](../tools/index.md#how-far-to-trust-the-numbers)).
:::

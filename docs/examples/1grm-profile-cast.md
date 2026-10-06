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
publish files=45 out_dir=1GRM_bundle profile=resolved
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
`scripts/check_chimerax.py --xvfb`). Teal: the cast surface. Grey-blue, 25%
transparent: the protein cartoon. No mouth guides are drawn by default
(`--mouth-guides` adds orange rings at the mouth planes). The check confirmed
that the displayed surface matches the reference mesh within 2.2 × 10⁻⁶ Å and
that the camera matches the scene file.
The PyMOL and VMD versions of this scene were not opened for this page.
```

`1GRM_residue_contacts.csv` ranks the residues around the profile; at the
default 4.5 Å cutoff the top rows are residues whose atoms touch the sampled
region (`min_distance_A` 0, role `bottleneck`), for example A:ALA3 and A:VAL7.
`1GRM_network_chord.png` and `1GRM_network_summary.png` summarise the residue
contact network.

## What was checked

These commands were run with CREVICE 0.1.0 on the files downloaded above, and
the printed values and file counts on this page are copied from that run.
They are measurements with the stated settings: no 1GRM result is a validated
biological value (see [Validation levels](../tools/index.md#validation-levels)).

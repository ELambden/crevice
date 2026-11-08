# 1OED: a wide pore with gaps in its wall

[PDB 1OED](https://www.rcsb.org/structure/1OED) is the transmembrane pore
domain of a nicotinic acetylcholine receptor: five helices around a central
pore, determined without the rest of the receptor. Between neighbouring helices
there are openings wide enough for a small enclosure probe to pass, so a
single small probe sees the lumen leak sideways. This example shows how the
per-slice enclosure probe deals with that. It takes about two minutes.

```bash
crevice fetch 1OED --cache-dir .crevice/pdb
crevice profile 1OED -o 1OED_profile.json --png 1OED_profile.png --dpi 150
```

```text
profile points=81 min_radius=3.264 mean_radius=4.282 bottleneck_index=24
enclosure_radius=per-plane 0.30-3.00 (auto: per-plane enclosure probe, 0.30-3.00 A over 81 samples)
```

Each slice chose the smallest probe at which its lumen is enclosed (see
[Choosing the enclosure probe](../tools/profile.md#choosing-the-enclosure-probe)).
In 25 of the 81 slices the wall is tight and the smallest probe, 0.30 Å,
already encloses the lumen. Where there are gaps between the helices the slice
needs a larger probe: up to 1.6 Å within 15 Å of the centre, and up to
3.0 Å towards the two ends, where the pore widens to more than 5 Å. The channel
ends at the first slice on each side whose lumen is not enclosed by any probe
smaller than its own clearance (t = −23.5 and +22.7 Å). Each slice's probe is
in `metadata.enclosure_probes_A` of `1OED_profile.json`.

```{figure} images/1oed_profile.png
:alt: Pore radius of 1OED against the axial coordinate
:width: 85%

`1OED_profile.png`: pore radius (Å) against the axial coordinate (Å) over
46.2 Å between the two mouths (light grey outside them). The narrowest sample,
3.264 Å, is at t = −9.7 Å (nearest residue A:SER248; red dashed line and
dot). The enclosure probes only decide which gaps count as wall and where the
channel ends; the plotted radius is still the distance to the nearest atom
surface.
```

## One probe for the whole pore

A number fixes the probe in every slice instead. Too small a probe and the
lumen leaks at once: with `--enclosure-radius 0.5` even the starting slice is
not enclosed, and no channel is resolved (exit status 2). A probe that closes
the gaps but still fits the lumen gives a profile with the same minimum, over
a shorter stretch:

```bash
crevice profile 1OED --enclosure-radius 2.2 -o 1OED_profile_2.2.json
```

```text
profile points=81 min_radius=3.265 mean_radius=4.063 bottleneck_index=25
```

| Enclosure probe | Extent (Å) | Minimum radius (Å) | Mean radius (Å) |
|---|---|---:|---:|
| per slice (default) | −23.5 to +22.7 | 3.264 | 4.282 |
| fixed 2.2 Å | −22.5 to +18.6 | 3.265 | 4.063 |
| fixed 1.0 Å | −20.1 to +3.9 | 3.266 | 3.751 |
| fixed 0.5 Å | not resolved | | |

With a fixed probe the channel ends where the pore opens wider than the probe
can enclose (2.2 Å) or where a gap in the wall is wider than the probe
(1.0 Å). The minimum is the same in every resolved run, because it lies in the
part of the pore that all of them follow.

## Residues lining the pore

```bash
crevice residues 1OED -o 1OED_residues.csv --png 1OED_residues.png --dpi 150
```

```text
residues contacts=125 output=1OED_residues.csv
residue,min_distance_A,mean_distance_A,atom_count,point_count,role,properties,influence_score
C:SER258,0.0,2.0014,6,34,bottleneck,polar;small,5.5333
D:THR244,0.0,2.1844,7,34,bottleneck,polar,5.5333
D:SER248,0.0,1.8985,6,34,bottleneck,polar;small,5.5333
```

(values rounded here.) `min_distance_A` is the gap between the residue's atom
surfaces and the sampled region surface; 0 means the residue bounds the
region. The influence score is a geometric ranking, not a functional
importance.

## Settings matter

The section grid and the search square decide which gaps are resolved. With
0.25 Å section grids and a 10 Å search radius (finer and wider than the
defaults used above; `--section-spacing 0.25 --search-radius 10`) the channel
has the same extent and a minimum of 3.265 Å (mean 4.277 Å). Report the
probe and the grid settings with any 1OED number, and compare settings before
interpreting small differences.

:::{admonition} Interpreting the results
:class: crevice-interpret

The values come from one run of the commands above on the downloaded entry,
including the 0.25 Å grid result. The pore axis was
found automatically. An isolated pore domain lacks the rest of the receptor and
its membrane, so read these radii as a measurement of this model with these
settings.
:::

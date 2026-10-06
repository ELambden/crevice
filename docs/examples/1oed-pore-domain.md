# 1OED: a wide pore with gaps in its wall

[PDB 1OED](https://www.rcsb.org/structure/1OED) is the transmembrane pore
domain of a nicotinic acetylcholine receptor: five helices around a central
pore, determined without the rest of the receptor. Between neighbouring helices
there are openings about 1.5 Å wide. A small enclosure probe passes through
them, so the lumen is not laterally enclosed and no channel resolves; this
example shows how the automatic enclosure probe deals with that. It takes about
two minutes.

```bash
crevice fetch 1OED --cache-dir .crevice/pdb
crevice profile 1OED -o 1OED_profile.json --png 1OED_profile.png --dpi 150
```

```text
profile points=81 min_radius=3.240 mean_radius=3.735 bottleneck_index=7
enclosure_radius=2.000 (auto: smallest probe of channel group 1, which resolved at 3 of 21 probes tried)
```

Of the 21 ladder probes (0.5–3.0 Å), only 2.0, 2.4 and 3.0 Å closed the gaps
between helices and still fitted through the pore; all three found the same
channel, so the smallest, 2.0 Å, was used. The full record of every probe tried,
with its status and reason, is in `metadata.enclosure_probe_selection` of
`1OED_profile.json`. Passing that probe explicitly gives the identical profile
in a few seconds, which is how to rerun or vary it:

```bash
crevice profile 1OED --enclosure-radius 2.0 -o 1OED_profile_2.0.json
```

```text
profile points=81 min_radius=3.240 mean_radius=3.735 bottleneck_index=7
```

(every sample radius is identical to the automatic run.)

```{figure} images/1oed_profile.png
:alt: Pore radius of 1OED against the axial coordinate
:width: 85%

`1OED_profile.png`: pore radius (Å) against the axial coordinate (Å) over
33.5 Å between the two mouths (light grey outside them). The narrowest sample,
3.240 Å, is near one end (t = −15.7 Å, nearest residue D:THR244; red dashed
line and dot). The enclosure probe only decides which gaps count as wall; the
plotted radius is still the distance to the nearest atom surface.
```

Residues lining the profile, with the same probe:

```bash
crevice residues 1OED --enclosure-radius 2.0 -o 1OED_residues.csv --png 1OED_residues.png --dpi 150
```

```text
residues contacts=107 output=1OED_residues.csv
residue,min_distance_A,mean_distance_A,atom_count,point_count,role,properties,influence_score
C:LEU265,0.0,2.3487,8,54,bottleneck,hydrophobic,5.6011
C:SER262,0.0,2.0116,6,51,bottleneck,polar;small,5.5927
A:SER248,0.0,2.1026,6,50,bottleneck,polar;small,5.5898
```

(values rounded here.) `min_distance_A` is the gap between the residue's atom
surfaces and the sampled region surface; 0 means the residue bounds the
region. The influence score is a geometric ranking, not a functional
importance.

## Settings matter

The chosen probe depends on the other settings, because they decide which
gaps are resolved. With 0.25 Å section grids and a 10 Å search radius (finer
and wider than the defaults used above) the automatic rule chooses 1.6 Å and
the narrowest radius is 3.26 Å. Report the probe and the grid settings with any
1OED number, and compare settings before interpreting small differences.

## What was checked

The commands on this page were run with CREVICE 0.1.0 on the downloaded entry,
and the printed values are copied from that run (the 0.25 Å result is from the
development record of the same version). The pore axis was found automatically,
not taken from a curated reference. These are measurements with stated
settings, not validated biological values.

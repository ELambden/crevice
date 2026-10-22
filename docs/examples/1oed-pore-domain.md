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
profile points=81 min_radius=3.264 mean_radius=3.826 bottleneck_index=23
enclosure_radius=2.200 (auto: smallest probe of channel group 2, which resolved at 4 of 21 probes tried)
```

Of the 21 ladder probes (0.5–3.0 Å), six resolved a channel, in two groups.
At 1.0 and 2.0 Å the channel ended at about t = +2 Å (group 1); at 2.2, 2.4,
2.7 and 3.0 Å the probe closed the gaps between helices and the channel ran
between both mouths (group 2). Group 2 resolved at the most probes, so its
smallest probe, 2.2 Å, was used. The full record of every probe tried,
with its status and reason, is in `metadata.enclosure_probe_selection` of
`1OED_profile.json`. Passing that probe explicitly gives the identical profile
in a few seconds, which is how to rerun or vary it:

```bash
crevice profile 1OED --enclosure-radius 2.2 -o 1OED_profile_2.2.json
```

```text
profile points=81 min_radius=3.264 mean_radius=3.826 bottleneck_index=23
```

(every sample radius is identical to the automatic run.)

```{figure} images/1oed_profile.png
:alt: Pore radius of 1OED against the axial coordinate
:width: 85%

`1OED_profile.png`: pore radius (Å) against the axial coordinate (Å) over
31.3 Å between the two mouths (light grey outside them). The narrowest sample,
3.264 Å, is at t = −9.8 Å (nearest residue D:SER248; red dashed line and
dot). The enclosure probe only decides which gaps count as wall; the
plotted radius is still the distance to the nearest atom surface.
```

Residues lining the profile, with the same probe:

```bash
crevice residues 1OED --enclosure-radius 2.2 -o 1OED_residues.csv --png 1OED_residues.png --dpi 150
```

```text
residues contacts=112 output=1OED_residues.csv
residue,min_distance_A,mean_distance_A,atom_count,point_count,role,properties,influence_score
D:SER248,0.0,1.8756,6,50,bottleneck,polar;small,5.5898
C:LEU265,0.0,2.3207,8,49,bottleneck,hydrophobic,5.5868
A:LEU251,0.0,2.4513,8,48,bottleneck,hydrophobic,5.5838
```

(values rounded here.) `min_distance_A` is the gap between the residue's atom
surfaces and the sampled region surface; 0 means the residue bounds the
region. The influence score is a geometric ranking, not a functional
importance.

## Settings matter

The chosen probe depends on the other settings, because they decide which
gaps are resolved. With 0.25 Å section grids and a 10 Å search radius (finer
and wider than the defaults used above) the automatic rule chooses 1.6 Å and
the narrowest radius is again 3.26 Å. Report the probe and the grid settings with any
1OED number, and compare settings before interpreting small differences.

:::{admonition} Interpreting the results
:class: crevice-interpret

The values come from a CREVICE 0.4.0 run on the downloaded entry, including
the 0.25 Å grid result. The pore axis was
found automatically. An isolated pore domain lacks the rest of the receptor and
its membrane, so read these radii as a measurement of this model with these
settings.
:::

# `profile`: through-channel radius profile

## Scientific question

How wide is a channel that passes right through a protein, and where is it
narrowest? `profile` measures the radius of the largest atom-clear sphere at
evenly spaced positions along the channel axis. It works in the spirit of HOLE,
with CREVICE's own connectivity checks. It applies only to passages with **two
open ends** that run roughly in one direction. Closed pockets and one-sided
cavities need [`cast`](cast.md) or [`cavities`](cavities.md).

## Method

1. **Choose a direction.** With `--axis auto` (the default), CREVICE tries the
   three principal directions of the selected atoms. A direction is accepted
   only if it gives a connected run of laterally enclosed free-space sections
   with two unobstructed axial exits. `--axis x|y|z|x,y,z` fixes the direction.
   `--origin x,y,z` gives a seed inside the intended channel. With an explicit
   origin, the channel must connect to that seed without passing through atom
   spheres.
2. **Find connected sections.** At each axial position, probe-centre space is
   sampled on a transverse grid (`--section-spacing`). A centre qualifies if a
   sphere of radius `--enclosure-radius` fits there. This is a geometric
   selection parameter, not an ion or water radius. By default it is chosen
   automatically (see "Choosing the enclosure probe" below). A section is a connected
   group of such centres that cannot reach the edge of the search region
   (`--search-radius`), i.e. it is laterally enclosed. Sections in neighbouring
   planes are joined only by a straight segment that also clears the enclosure
   probe, and they must connect from one end of the channel to the other.
3. **Choose the lumen where a plane has several sections.** Besides the lumen, a
   plane can contain small enclosed pockets in the wall, often a single grid
   point where the enclosure probe just fits between atoms. When such a pocket
   can be joined to the lumen on both sides, several equally long paths exist.
   CREVICE follows the **widest path**: the one whose narrowest straight
   segment (its bottleneck) has the largest clearance. Remaining ties go to the
   smaller sideways step, then to the path whose middle sample is nearest the
   origin, then to a fixed grid order, so the same input always gives the same
   path. Pockets are not removed by size, because a narrow lumen can itself be
   a single grid point.
4. **Measure radii.** At each sample, the reported radius is the distance from
   the local connected centre to the nearest atom van der Waals surface, minus
   `--probe-radius` (default 0). Segments between samples are also checked
   exactly, so `continuous_bottleneck_radius_A` can be smaller than every
   sampled radius.
5. **Mark the ends.** Radius plots show where lateral enclosure is lost.
   Brackets are refined to 0.025 Å or better.

If there are no candidates, or several candidates are equally good, the command
fails with an explanation. `--search-radius 0` asks for a simpler fixed-axis
clearance scan instead. Its metadata says that no channel assignment was
checked.

## Choosing the enclosure probe (`--enclosure-radius auto`)

The enclosure probe decides which openings in the channel wall count as wall.
Too small a probe lets the lumen leak sideways through gaps between atoms, so
it is not laterally enclosed; too large a probe no longer fits through the
constriction. Real structures need different probes: gramicidin A (1GRM)
resolves at 0.8–1.2 Å, the TWIK-1 filter (3UKM) only at about 0.54–0.61 Å,
and the isolated nAChR pore domain (1OED), whose helices leave gaps near
1.5 Å, only from about 1.6 Å. Within the window that resolves a channel, the
measured radii hardly change (by ≤ 0.005 Å in these cases), but the mouths
move outward as the probe grows.

With `--enclosure-radius auto` (the default; Python `enclosure_radius=None`)
CREVICE chooses the probe:

1. It resolves the channel with every probe of a fixed ladder, 0.50, 0.55,
   0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50,
   1.60, 1.80, 2.00, 2.20, 2.40, 2.70 and 3.00 Å, using all the other settings
   and the rules above unchanged. Probes smaller than 0.5 Å would pass through
   gaps that hydrogens (removed from the model) fill; probes larger than 3 Å
   would treat openings about two water molecules wide as wall. Probes at or
   above the clearance of an explicit `--origin` are skipped.
2. Results are grouped into channels: two results are the same channel when
   their axes agree (|cos| ≥ 0.99), their end types agree, their minimum radii
   agree within 0.05 Å or 2 %, and their axial ranges overlap by at least 60 %
   of their union.
3. Only a channel that resolves at **two or more** probes is eligible. If every
   channel resolved at a single ladder probe, the two half-step probes beside
   it are also tried.
4. The eligible channel that resolves at the most probes wins (ties: the one
   found at the smaller probe), and its **smallest** probe is used.

The rule prefers a channel that persists over a range of probes to one that
appears at a single threshold. Taking the smallest probe of that range keeps
the wall definition as strict as possible: only openings narrower than the
probe are closed, and the mouths are the most conservative. The largest
resolving probe was not used because it sits just under the constriction,
where the grid decides success. The smallest resolving probe alone was not
used because it can pick a different, spurious run: on 3UKM, 0.51–0.53 Å
resolves a short run that stops at the filter, and on 1OED with 0.5 Å
sections, 0.9 Å resolves a narrow run through a wall pocket.

The profile is then exactly what that probe gives when passed explicitly.
`metadata.enclosure_probe_selection` records the ladder, every probe tried
(status, short reason, minimum, continuous bottleneck, axial range, group),
the groups, the chosen probe and the reason, and the CLI prints the chosen
probe. If no channel resolves at two probes, the profile is **unresolved**:
the error lists the probes tried and why each failed, and names any channel
seen at a single probe only. It is never reported automatically, but you can
pass that probe explicitly. An explicit number (for example
`--enclosure-radius 0.8`) runs that probe only, exactly as before, and adds
no selection record. Report the chosen probe with your results.

The automatic choice costs one profile run per ladder probe (about 10 s for
1GRM with 0.25 Å sections, about a minute for 1OED or 4PYP). `static-suite`
and `benchmark` also default to `auto` (per system). `trajectory` (CLI and
Python) chooses once on the reference frame and uses that probe for every
frame. `cast` keeps its 0.8 Å rolling probe, which plays a different role
(see [`cast`](cast.md#why-the-probe-radius-is-fixed)).

## Capped channels: lateral exits (`--lateral-exits`)

Some channels are capped: a domain blocks the straight exit at one end, and
ions leave sideways through portals (for example the K2P channel TWIK-1, whose
extracellular cap sits over the pore). By default such a channel is refused.
With `--lateral-exits` (Python: `lateral_exits=True`) an end whose straight
exit is blocked is accepted **only if a lateral exit to bulk exists**:

1. The axial profile is found and sampled exactly as above, between the two
   mouths where lateral enclosure is lost. Enclosed space beyond a capped
   mouth (for example under the cap) is treated as space the exit may cross,
   not as a continuation of the channel.
2. From the last sample at the capped end, CREVICE searches for free paths
   into bulk, widest first. The path runs over a lattice (`--exit-spacing`, default
   0.5 Å). Every straight step must clear the enclosure probe. The path may
   not cross back through the plane of that last sample into the channel.
   "Widest" means the path's narrowest step has the largest clearance.
3. **Bulk** is the space swept by a large rolling probe (`--exit-bulk-radius`,
   default 6 Å, as for the rolling-probe cast) that can reach the outside.
4. If no such path exists (for example the cap closes the pore), the channel is
   not resolved.
5. **Every distinct exit is reported.** The widest-path search continues after
   the first exit. Route widths are counted beyond the inscribed sphere of the
   last sample (the mouth, already part of the axial profile); otherwise every
   route would share the mouth's width and all exits would tie. Equally wide
   routes are ranked by length. A further exit is kept when its exit point is
   more than twice the bulk-probe radius (12 Å by default) from every exit
   already kept, so that no single bulk-probe position touches both, and the
   parts of the two routes farther than one bulk radius from the start never
   come within one bulk radius of each other, so that they leave through
   different openings. Symmetric portals, such as the two C2-related portals
   of TWIK-1, therefore give one leg each. At most eight legs are reported per
   end. Each leg reports its whole-route bottleneck (including the mouth),
   `ranking_width_raw_clearance_A` (the width used to rank it) and
   `beyond_mouth_min_raw_clearance_A` (its narrowest sample outside the mouth
   sample's sphere).

The result is reported as a distinct path type. `metadata.path_type` is
`capped_channel_lateral_exit` (or `through_channel` when both ends are open).
Each end in `metadata.exits` is `axial`, or `lateral` with `leg_count`, the
largest leg `bottleneck_radius_A` and a `legs` list. Each leg (`leg` = rank)
has its points, radii, bottleneck and its position, exit point and length.
`path_bottleneck_radius_A` is the smallest of the axial bottleneck and the
best leg bottleneck of each capped end. The publish manifest
records `path_type`, `exit_types` and `exit_leg_counts`.

The leg's radii are distances to the nearest atom surface along a 3D path,
like tunnel radii. They are **not** planar cross-section radii. Therefore they
are not part of the sampled profile, its CSV or `min_radius`. The radius plots
draw each leg separately, as a bluish-green dash-dot continuation beyond its
mouth against the path length along the leg (see Outputs). In `publish`, the
measured channel cast still fills only the axial lumen; each leg gets a
separate probe-swept cast (reported with its own volume) that the scenes draw
joined to the channel cast, so the widening towards the portal is shown. The
thin blue centre-line tube of each leg (RGB 0.16, 0.47, 0.84) is optional
(`--exit-centre-lines`); it shows the route, not the portal's width.

**How this differs from HOLE.**
- HOLE keeps planes normal to one vector and lets its sphere centre slide
  sideways within each plane until the radius exceeds its end radius. Under a
  cap, HOLE's radius is therefore the largest sphere in a plane that cuts
  obliquely through a portal.
- CREVICE instead stops the planar profile where the portal opens. It then
  reports the portal separately as a free 3D path to bulk.
- CREVICE refuses a closed cap. HOLE still returns a profile.

Without `--lateral-exits` nothing changes: the same outputs are written
byte for byte. The option is not available for the `--search-radius 0` scan or
for `trajectory`.

## Wide pores with gaps between wall helices

Every plane must hold a laterally enclosed section at the enclosure probe. In
some structures (for example an isolated pore domain without lipids) gaps
between helices wider than the legacy 0.8 Å probe connect the lumen sideways
to the outside. The channel is then refused, even though its lumen is wide.
When the lumen is much wider than those gaps, a larger enclosure probe
closes the gaps without touching the lumen; the automatic choice finds it
(1OED: 1.6 Å with 0.25 Å sections). The measured radii are distances to
atom surfaces and do not depend on the probe unless the lumen itself is
narrower than it. If resampling loses continuity, the error names the first
axial position that could not be reached. Report the enclosure radius
you used.

## Assumptions

- The channel is monotone along one direction and open at both ends.
- The selected atoms are the correct walls. Heteroatoms are kept unless you
  pass `--exclude-hetero`, and hydrogens are dropped unless you pass
  `--include-hydrogen`.
- Atomic radii are CREVICE's van der Waals set (`crevice.radii`).

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--axis` | `auto` | direction | `auto`, `x`, `y`, `z` or a vector `x,y,z` |
| `--origin` | none | Å | seed inside the intended channel |
| `--samples` | 81 | count | minimum number of axial samples; connected profiles add points to keep axial steps ≤ 0.75 Å |
| `--section-spacing` | 0.5 | Å | transverse grid spacing for finding sections |
| `--enclosure-radius` | `auto` | Å | probe radius used to decide connectivity and enclosure; `auto` chooses it (see above), a number fixes it |
| `--probe-radius` | 0.0 | Å | subtracted from the reported radius |
| `--search-radius` | 6.0 | Å | transverse search half-width; `0` requests an unvalidated fixed-axis scan |
| `--lateral-exits` | off | flag | accept a capped end if a lateral exit to bulk exists (see above) |
| `--exit-bulk-radius` | 6.0 | Å | rolling probe that defines bulk for lateral exits; exits more than twice this apart are distinct |
| `--exit-spacing` | 0.5 | Å | lattice spacing of the lateral-exit search |
| `--mouth-guides` | off | flag | draw the orange mouth bands and lines on the radius plots |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

## Outputs

- `-o profile.json`: sampled centres, radii, section areas, arc lengths, chosen
  axis and seed, candidate directions, exit clearances, bottleneck and settings.
  `volume_estimate` is an integral over inscribed circles. For a volume that
  follows a non-circular shape, use a cast.
- `<stem>.csv` beside the JSON (always; `--csv PATH` chooses another path):
  one row per sample: `index`, `x_A`, `y_A`, `z_A`, `radius_A`,
  `raw_clearance_A`, `axis_position_A` (all in Å), `nearest_atom_serial`,
  `nearest_residue`.
- `--pdb`: a trace of dummy atoms. `--pymol` writes a PyMOL script for it
  (needs `--pdb`).
- `--png`: radius plot: pore radius (Å, blue line without per-sample markers,
  and light-blue fill) against axial coordinate (Å), bottleneck marked by a
  red dashed line and dot; light-grey shading outside the channel mouths (when
  present) and the display-only exterior clearance beyond an open end as a grey
  dashed line. The orange mouth bands and dotted mouth lines are drawn only with
  `--mouth-guides`; the mouth positions are always in the JSON and the PNG
  metadata. For a capped channel each lateral exit leg is a
  bluish-green (`#009e73`) dash-dot line continuing from its mouth: x is the
  mouth's axial coordinate plus (upper end) or minus (lower end) the path
  length along the leg, y the leg's free-sphere radius. Beyond the mouth the
  x axis is therefore path length, not axial position; the legend says so.
- `--annotated-png`: the same radius plot above a lane of nearest-residue
  landmarks every `--residue-label-every` samples (plus the bottleneck), shown
  as dots coloured by chemistry: purple charged, teal polar, blue other.
  Residue names are listed in the PNG metadata; they, the bottleneck value and
  the titles are drawn only with `--annotate`.
- `--pymol`: PyMOL script with the dummy atoms as translucent spheres
  coloured blue (narrow) to red (wide) by radius and a pink surface, inside a
  cyan cartoon; no text.

Figures have no titles or value labels unless `--annotate` is given (see
[Figure and scene text](../getting-started/outputs.md#figure-and-scene-text---annotate)).

## Python equivalent

```python
from crevice import load_structure, pore_profile, write_profile_json, plot_profile_radius

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto", samples=101, section_spacing=0.25)  # enclosure probe chosen automatically
print(profile.metadata["enclosure_probe_selection"]["chosen_A"])
write_profile_json(profile, "1GRM_profile.json")
plot_profile_radius(profile, "1GRM_profile_radius.png")
```

The public wrapper `crevice.pore_profile(frame, **kwargs)` passes keyword
arguments to {func}`crevice.analysis.pore_profile`. See that function for the
full parameter list.

## Testing and validation status

- **Software / synthetic:** exact segment clearance, bottlenecks between
  samples and seed attachment are tested on mathematical atom walls. A 48-case
  sweep of grid resolution and offset checks connectivity
  (`examples/connectivity_validation.py`).
- **Example observations (public entries, stated settings):** 1GRM (mmCIF,
  model 0) gives a 101-sample connected profile on 0.25 Å section and cast
  grids, with a minimum continuous-path radius of 1.330 Å; rigidly rotated
  copies of 1GRM resolve the same channel. 4PYP gives **no** unique
  through-profile, which fits an inward-open transporter structure but does not
  test it. Results depend on settings: the quick-start command (mmCIF input,
  81 samples, 0.5 Å sections, automatic enclosure probe 0.75 Å) gives
  1.324 Å (see the [1GRM example](../examples/1grm-profile-cast.md)); 3UKM
  with `--lateral-exits` and 1OED are shown in the
  [3UKM](../examples/3ukm-lateral-exits.md) and
  [1OED](../examples/1oed-pore-domain.md) examples. With a fixed axis, 61–150
  samples on 0.25 Å and 0.125 Å section grids give minima of 1.330–1.339 Å.
  The 1GRM wall contains one-point pockets that tie with the lumen at some
  sample positions; the widest-path rule (Method, step 3) keeps the profile in
  the lumen.
- **Software / synthetic (path choice):** a lumen with four one-plane wall
  pockets that tie with the lumen checks that the profile stays in the lumen
  and that the minimum does not change with the sample count (81, 101, 118,
  150) by more than the final refinement step.
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL (side and end views of the 1GRM profile scene) for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established. The 1GRM axis was not taken
  from a curated reference.

## Known limitations

- Branched, curved or folded passages, and axes far from the principal
  directions, are not handled.
- A failed profile does not show that a protein is closed.
- Grid-based section discovery can miss openings narrower than the section
  spacing. Compare spacings before you interpret small differences.
- Radius and position uncertainty from the transverse grid and enclosure probe
  are not included in the plotted end brackets.
- Each plane must hold a laterally enclosed section. If the wall has side
  openings wider than the enclosure probe (for example between helices of an
  isolated membrane domain), the lumen reaches the search edge there and the
  channel is not resolved. A larger `--enclosure-radius` closes narrower
  openings but also excludes narrower lumens.
- By default both axial exits must be straight and unobstructed. A capped
  channel is followed through a side portal only with `--lateral-exits`.
  Distinct exits are separated by a distance rule (twice the bulk probe); two
  openings closer than that are reported as one.
- The automatic enclosure probe is chosen from a fixed ladder, so a window
  narrower than the ladder step can be missed (3UKM's window is about 0.07 Å
  wide and contains two ladder probes). Its result depends on the sampling
  and grid settings, like everything else here, and the mouths depend on the
  chosen probe.
- End brackets stop early if a plane just beyond the enclosed run holds more
  than one joinable section (a branch, or a wall pocket next to the lumen).
  This was not seen for 1GRM, but it is not yet ranked like the path choice.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: profile
```

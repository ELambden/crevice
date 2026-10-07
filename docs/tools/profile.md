# Pore radius profiles

`crevice profile` tells you how wide a channel is along its length and where
it is narrowest. At evenly spaced points along the channel it measures the
radius of the largest sphere that fits without touching an atom, and it checks
the path between the points too, so a constriction can't slip through between
two samples.

Use it for passages that run right through a protein and are open at both
ends, such as ion channels and pores. For pockets, buried cavities and
one-sided tunnels, use a [cast](cast.md), [`cavities`](cavities.md) or
[`tunnels`](tunnels.md) instead.

## Quick start

```bash
crevice profile 1GRM -o 1GRM_profile.json --png 1GRM_profile.png
```

```python
from crevice import load_structure, plot_profile_radius, pore_profile, write_profile_json

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto", samples=101, section_spacing=0.25)
print(profile.min_radius, profile.metadata["enclosure_probe_selection"]["chosen_A"])
write_profile_json(profile, "1GRM_profile.json")
plot_profile_radius(profile, "1GRM_profile_radius.png")
```

{func}`crevice.pore_profile` passes its keyword arguments to
{func}`crevice.analysis.pore_profile`, which lists every parameter.

## How it works

1. **Find the direction.** With `--axis auto` (the default) CREVICE tries the
   three principal directions of the selected atoms and keeps one only if it
   holds a connected, enclosed passage with two open ends. You can fix the
   direction with `--axis x|y|z|x,y,z`, and point at the right channel with
   `--origin x,y,z`, a seed inside it.
2. **Slice and connect.** At each position along the axis, CREVICE looks for
   space where an *enclosure probe* fits (`--enclosure-radius`) and that is
   walled in sideways, so the probe can't reach the edge of the search disc
   (`--search-radius`). Neighbouring slices are joined only by straight
   segments that also clear the probe, and the joined slices must run from one
   end of the channel to the other.
3. **Follow the widest path.** A slice can contain small pockets in the
   channel wall as well as the lumen. When several equally long routes exist,
   CREVICE follows the one whose narrowest step is widest, then breaks any
   remaining tie in a fixed order, so the same input always gives the same
   path.
4. **Measure.** The radius at each point is the distance from the path to the
   nearest atom surface, minus `--probe-radius` (0 by default). The segments
   between points are checked exactly, so `continuous_bottleneck_radius_A` can
   be smaller than every sampled radius.
5. **Mark the ends.** The plot shows where the walls stop enclosing the
   channel, to within 0.025 Å.

If no channel qualifies, or two are equally good, the command stops and
explains why. `--search-radius 0` instead runs a simple fixed-axis clearance
scan, and its metadata says that no channel was checked.

### Choosing the enclosure probe

The enclosure probe decides which gaps in the wall count as wall. Too small,
and the lumen leaks sideways between atoms; too large, and it no longer fits
through the constriction. Different structures need different probes: in our
examples gramicidin A (1GRM) resolves at 0.8 to 1.2 Å, the TWIK-1 filter
(3UKM) only at about 0.54 to 0.61 Å, and the isolated nAChR pore domain (1OED),
whose helices leave gaps near 1.5 Å, only from about 1.6 Å. Within the window
that works, the radii barely change (by 0.005 Å or less in these cases), but
the channel mouths move outwards as the probe grows.

So by default (`--enclosure-radius auto`, Python `enclosure_radius=None`)
CREVICE chooses for you:

1. It tries every probe on a fixed ladder from 0.5 to 3.0 Å (0.50, 0.55, …,
   0.80, then 0.90 to 1.60 in 0.1 Å steps, 1.80, 2.00, 2.20, 2.40, 2.70, 3.00),
   with all other settings unchanged. Probes at or above the clearance of an
   explicit `--origin` are skipped.
2. It groups the results into channels: the same channel has a parallel axis
   (|cos| ≥ 0.99), the same end types, a minimum radius within 0.05 Å or 2 %,
   and axial ranges that overlap by at least 60 %.
3. A channel must appear at two or more probes. If none does, the two
   half-step probes beside each single hit are tried as well.
4. The channel found at the most probes wins (ties go to the smaller probe),
   and its smallest probe is used.

Preferring a channel that persists over a range of probes avoids one that
appears at a single threshold, and taking the smallest probe of that range
keeps the wall as strict as possible. The result is exactly what you would get
by passing that probe yourself. `metadata.enclosure_probe_selection` records
every probe tried and why it passed or failed, and the command prints the
choice. If no channel resolves at two probes, the profile is reported as
unresolved, with the reasons; you can still pass a probe explicitly. A number
such as `--enclosure-radius 0.8` runs that probe only.

The automatic choice runs one profile per probe, so it takes about 10 s for
1GRM with 0.25 Å sections and about a minute for 1OED or 4PYP. `trajectory`
chooses once, on the reference frame, and uses that probe for every frame.
`cast` uses its own fixed 0.8 Å probe, which plays a different role (see
[why the cast probe is fixed](cast.md#why-the-probe-radius-is-fixed)).

### Capped channels and side exits

Some channels are capped: a domain blocks the straight exit at one end, and
ions leave sideways through portals. The K2P channel TWIK-1 (3UKM) is an
example. By default CREVICE refuses such a channel. With `--lateral-exits`
(Python `lateral_exits=True`) it accepts a blocked end **only if a side exit to
bulk solvent exists**:

1. The axial profile is measured as usual, between the two mouths where the
   walls stop enclosing it.
2. From the last point at the capped end, CREVICE searches a lattice
   (`--exit-spacing`, 0.5 Å) for the widest free route into bulk that never
   turns back into the channel. Every step must clear the enclosure probe.
3. Bulk is the space a large rolling probe (`--exit-bulk-radius`, 6 Å) can
   reach from outside.
4. Every distinct exit is reported, up to eight per end. Two exits count as
   distinct when their exit points are more than twice the bulk radius apart
   and their routes stay apart, so symmetric portals, like the two portals of
   TWIK-1, give one leg each.

`metadata.path_type` becomes `capped_channel_lateral_exit`, and each capped
end in `metadata.exits` lists its legs with their points, radii, bottleneck,
exit point and length. A leg's radii are free-sphere clearances along a 3D
route, not slice radii, so they are kept out of the sampled profile and its
CSV. The plot draws each leg as a green dash-dot continuation beyond its mouth,
and in `publish` each leg gets its own cast. Without `--lateral-exits`, nothing
changes.

**How this compares with HOLE.** HOLE (Smart *et al.*, 1996) keeps its planes
perpendicular to one vector and lets the sphere slide sideways within each
plane, so under a cap it measures the largest sphere in a plane that cuts
obliquely through a portal, and it returns a profile even when the cap is
closed. CREVICE instead stops the planar profile where the portal opens,
reports the portal separately as a 3D route to bulk, and refuses a closed cap.

### Wide pores with gaps in the wall

Every slice must be walled in at the enclosure probe. In an isolated pore
domain without its lipids, gaps between helices can be wider than a small
probe, so the lumen leaks sideways and the channel is refused even though it
is wide. A larger probe closes those gaps without touching a much wider
lumen, and the automatic choice finds it (1OED: 1.6 Å with 0.25 Å sections).
The measured radii do not depend on the probe unless the lumen is narrower
than it.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--axis` | `auto` | channel direction: `auto`, `x`, `y`, `z` or a vector `x,y,z` |
| `--origin` | none | a seed inside the channel you want, `x,y,z` in Å |
| `--samples` | 81 | minimum number of points; extra points keep steps at 0.75 Å or less |
| `--section-spacing` | 0.5 Å | grid spacing within each slice; 0.25 Å resolves finer detail |
| `--enclosure-radius` | `auto` | the enclosure probe, chosen automatically or fixed in Å |
| `--probe-radius` | 0 Å | subtracted from every reported radius |
| `--search-radius` | 6 Å | sideways search half-width; `0` runs an unchecked fixed-axis scan |
| `--lateral-exits` | off | accept a capped end if a side exit to bulk exists |
| `--exit-bulk-radius` | 6 Å | rolling probe that defines bulk solvent for side exits |

Atomic radii (`--radii`), figure text (`--annotate`), mouth guides and
structure input are [shared options](../reference/cli/common-options.md).
Every option is listed under [`crevice profile`](../reference/cli/profile.md).

## What you get

- **`-o profile.json`**: the sampled centres, radii, slice areas, chosen axis,
  seed, bottleneck, end clearances and every setting. Its `volume_estimate`
  adds up inscribed circles; for a volume that follows the real shape, use a
  cast.
- **`<stem>.csv`** beside the JSON (or `--csv PATH`): one row per point with
  `x_A`, `y_A`, `z_A`, `radius_A`, `raw_clearance_A`, `axis_position_A`,
  `nearest_atom_serial` and `nearest_residue`.
- **`--png`**: the radius plot. The profile is a blue line with a light-blue
  fill, the bottleneck a red dashed line and dot, and the region outside the
  mouths is shaded grey. Side-exit legs appear as green (`#009e73`) dash-dot
  lines; beyond the mouth their x axis is distance along the leg, as the legend
  says.
- **`--annotated-png`**: the same plot above a lane of nearest-residue
  landmarks every `--residue-label-every` points, coloured by chemistry
  (purple charged, teal polar, blue other).
- **`--pdb`** and **`--pymol`**: the path as dummy atoms, and a PyMOL scene
  with spheres coloured blue (narrow) to red (wide) inside a cartoon.

`profile` writes no cast. [`publish`](publish.md) casts the same channel and
splits it into entry, lumen and exit segments.

:::{admonition} Interpreting results
:class: crevice-interpret

A profile is a measurement with the settings you chose: the numbers move a
little with the grid, the sample count and the input format (1GRM gives a
minimum of 1.324 Å with the quick-start settings and 1.330 to 1.339 Å on finer
grids). Report the enclosure probe with your results. The method assumes a
channel that runs roughly in one direction; branched or strongly curved
passages are not followed, openings narrower than the slice grid can be
missed, and two exits closer than twice the bulk probe are reported as one. If
no profile is found, that does not show the protein is closed: 4PYP, an
inward-open transporter, gives no through-profile, which fits its structure
but doesn't test it. The software is checked on synthetic walls with known
answers; the real-protein numbers in these pages are observations, not
validated biology.
:::

See it in action in the [1GRM](../examples/1grm-profile-cast.md),
[3UKM](../examples/3ukm-lateral-exits.md) and
[1OED](../examples/1oed-pore-domain.md) tutorials.

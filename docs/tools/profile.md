# Pore radius profiles

`crevice profile` tells you how wide a channel is along its length and where
it is narrowest. At evenly spaced points along the channel it measures the
radius of the largest sphere that fits without touching an atom, and it checks
the path between the points too, so a constriction can't slip through between
two samples.

Use it for passages that run right through a protein and are open at both
ends, such as ion channels and pores. For pockets, buried cavities and
one-sided tunnels, use a [cast](cast.md), [`cavities`](cavities.md) or
[`paths`](paths.md) instead.

## Quick start

```bash
crevice profile 1GRM -o 1GRM_profile.json --png 1GRM_profile.png
```

```python
from crevice import load_structure, plot_profile_radius, pore_profile, write_profile_json

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto", samples=101, section_spacing=0.25)
print(profile.min_radius, profile.metadata["enclosure_probe_range_A"])
write_profile_json(profile, "1GRM_profile.json")
plot_profile_radius(profile, "1GRM_profile_radius.png")
```

{func}`crevice.pore_profile` passes its keyword arguments to
{func}`crevice.analysis.pore_profile`, which lists every parameter.

## How it works

1. **Find the direction.** With `--axis auto` (the default) CREVICE tries the
   three principal directions of the selected atoms and keeps one only if it
   holds a connected, enclosed passage. You can fix the direction with
   `--axis x|y|z|x,y,z`, and point at the right channel with `--origin x,y,z`,
   a seed inside it.
2. **Follow the lumen, slice by slice.** Starting from the seed, CREVICE moves
   along the axis in steps of at most 0.75 Å. In each slice it starts from the
   previous slice's centre and climbs, within the slice plane, to the point
   farthest from any atom surface (see [Slice centres](#slice-centres)). The
   path therefore stays in the lumen and never jumps to a separate pocket in
   the wall.
3. **Check that the slice is walled in.** The slice counts as part of the
   channel only if the free space around its centre is enclosed sideways: an
   *enclosure probe* placed in it cannot reach the edge of the search square
   (`--search-radius`). Each slice uses the smallest probe at which it is
   enclosed (see [Choosing the enclosure probe](#choosing-the-enclosure-probe)).
   The search square is a window, not a wall: a slice whose lumen just touches
   its edge is accepted if it is enclosed in a square twice as wide, and a wide
   lumen gets a wider square (twice or four times `--search-radius`, so the
   square is always at least twice as wide as the lumen).
4. **Join the slices.** Neighbouring centres are joined only by a straight
   segment that clears the smaller of the two slices' probes.
5. **Find the mouths.** Going outwards from the seed, the channel ends at the
   first slice that is not enclosed at any probe, or that cannot be joined, or
   at the outermost atom surfaces along the axis. Each mouth is located to
   within 0.025 Å. An end counts as open when the straight line along the axis
   from its mouth leaves the protein without touching an atom; otherwise it is
   capped (see [Capped channels and side exits](#capped-channels-and-side-exits)).
6. **Measure.** The points are spread evenly between the mouths. The radius at
   each point is the distance from the centre to the nearest atom surface,
   minus `--probe-radius` (0 by default). The segments between points are
   checked exactly, so `continuous_bottleneck_radius_A` can be smaller than
   every sampled radius.

The same input always gives the same profile. If no channel qualifies, or two
are equally good, the command stops and explains why.

### Slice centres

Each point of the profile is the centre of its slice, and its radius is that
centre's distance to the nearest atom surface. Within the slice plane this
distance is the smallest of the distances to each atom surface, so it rises
to peaks and ridges between the atoms. A peak is always a point equidistant
from three atom surfaces whose atoms surround it (a vertex); a point on a
ridge, equidistant from only two, can still move along the ridge.

CREVICE starts from the previous slice's centre, moved onto the new slice
plane, and climbs from it. At each step length (starting at half the grid spacing) it
tries the eight grid directions and one move chosen from the local atom
geometry, which can follow a ridge at any angle, takes the best move that
increases the clearance, and halves the step when none does, down to
0.0001 Å. It then solves for the peak exactly: the highest vertex within
0.35 Å that is higher still. Every move increases the clearance, so the climb
cannot cross a wall. A wide lumen can hold more than one hill, so CREVICE also
climbs from the slice's grid point of largest clearance and keeps the higher
of the two tops, provided both lie in the same enclosed space and the
straight line between them never comes closer to an atom than 90 % of the
lower top's clearance, so the two tops belong to one lumen. A neighbouring
passage joined to the lumen only through a narrower neck is never taken.
`metadata.section_refinement` records the method and its settings, and
`--refinement-steps 0` keeps the grid points instead (`refinement_steps=0`
in Python).

### Choosing the enclosure probe

The enclosure probe decides which gaps in the wall count as wall. Too small,
and the lumen leaks sideways between atoms; too large, and it no longer fits
through the constriction. A single probe for the whole channel often cannot do
both. In KcsA (1K4C), for example, the selectivity filter is only 0.75 Å wide
(Bondi radii), while the central cavity below it leaks through gaps between
the helices at every probe up to 0.7 Å.

So by default (`--enclosure-radius auto`, Python `enclosure_radius=None`)
each slice chooses its own probe: the smallest probe of a fixed ladder (0.30,
0.35, …, 0.80, 0.90, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80, 2.00, 2.40 and
3.00 Å) at which the slice's lumen is enclosed. The probe must be smaller than
the slice's own clearance. A space enclosed at one probe stays enclosed at
every larger one, so the smallest is well defined. Each slice's probe is the
size of the widest sideways opening of the lumen there, and it is reported
for every point (`metadata.enclosure_probes_A`, with the range in
`enclosure_probe_range_A`). The command prints the range.

The probes decide only where the channel is and where it ends; the measured
radii are atom-surface distances and do not depend on them.

A slice has to be resolved on the section grid. When its clearance is within
half a grid diagonal of the smallest probe, it is re-gridded at half the
spacing, up to eight times finer. The narrowest constriction that can be
followed is therefore just above the smallest probe: 0.32 Å with 0.25 Å
sections and 0.34 Å with the default 0.5 Å. This holds wherever the sample
points fall.

A number, such as `--enclosure-radius 0.8`, runs that probe in every slice
instead: a slice whose lumen is narrower than it ends the channel there.

`trajectory` repeats the per-slice choice in every frame and keeps the axis
and seed of the first frame. `cavities --cast` uses its own fixed 0.8 Å probe,
which plays a different role (see
[why the cast probe is fixed](cast.md#why-the-probe-radius-is-fixed)).

### Capped channels and side exits

Some channels are capped: a domain blocks the straight exit at one end, and
ions leave sideways through portals. The K2P channel TWIK-1 (3UKM) is an
example. By default CREVICE refuses such a channel. With `--lateral-exits`
(Python `lateral_exits=True`) it accepts a blocked end **only if a side exit to
bulk solvent exists**:

1. The axial profile is measured as usual, between the two mouths where the
   walls stop enclosing it.
2. A capped end leaves sideways in any case, so it is placed where the widest
   route to bulk leaves the lumen. From each slice between the end and the
   channel's narrowest slice (and the seed), CREVICE looks for a route to bulk
   that is clearly wider than the rest of the axial path beyond that slice,
   its own exit included (by 0.1 Å or 10 %, whichever is larger). Such a route
   cannot continue through the lumen, so it is a side exit; if one exists, the
   capped end moves to the slice with the widest one (mouth reason
   `lateral_branch`, recorded in `exits.<end>.lateral_branch`). This finds
   side portals that leave a chamber obliquely, which no single slice shows as
   an opening wider than the lumen. Open ends never move.
3. From the last point at the capped end, CREVICE searches a lattice
   (`--exit-spacing`, 0.5 Å) for the widest free route into bulk that never
   turns back into the channel. Every step must clear a probe: the end slice's
   own probe or, if no route clears it, about half of it and then the smallest
   probe (`exit_probe_A` records the one used).
4. Bulk is the space a large rolling probe (`--exit-bulk-radius`, 6 Å) can
   reach from outside.
5. Every distinct exit is reported, up to eight per end. Two exits count as
   distinct when their exit points are more than twice the bulk radius apart
   and their routes stay apart, so symmetric portals, like the two portals of
   TWIK-1, give one leg each.

`metadata.path_type` becomes `capped_channel_lateral_exit`, and each capped
end in `metadata.exits` lists its legs with their points, radii, bottleneck,
exit point and length. A leg's radii are free-sphere clearances along a 3D
route, not slice radii, so they are kept out of the sampled profile and its
CSV. The plot draws each leg as a green dash-dot continuation beyond its mouth,
and in `report` each leg gets its own cast. Without `--lateral-exits`, nothing
changes.

An end whose last point already lies in bulk solvent (inside the rolling-probe
envelope) is an open end even when the straight line along the axis is
blocked: it is never refused, with or without `--lateral-exits`
(`exits.<end>.type` is `bulk`, or `metadata.ends_in_bulk` without
`--lateral-exits`).

**How this compares with HOLE.** HOLE (Smart *et al.*, 1996) keeps its planes
perpendicular to one vector and lets the sphere slide sideways within each
plane, so under a cap it measures the largest sphere in a plane that cuts
obliquely through a portal, and it returns a profile even when the cap is
closed. CREVICE instead stops the planar profile where the portal opens,
reports the portal separately as a 3D route to bulk, and refuses a closed cap.

### Wide pores with gaps in the wall

In an isolated pore domain without its lipids, gaps between helices can be
wider than a small probe, so the lumen leaks sideways in those slices. Each
such slice uses a larger probe, which closes the gaps without touching a much
wider lumen, while the narrow slices keep their small probes. A slice whose
lumen opens wider than any probe on the ladder ends the channel there.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--axis` | `auto` | channel direction: `auto`, `x`, `y`, `z` or a vector `x,y,z` |
| `--origin` | none | a seed inside the channel you want, `x,y,z` in Å |
| `--samples` | 81 | minimum number of points; extra points keep steps at 0.75 Å or less |
| `--section-spacing` | 0.5 Å | grid spacing within each slice; 0.25 Å resolves finer detail |
| `--enclosure-radius` | `auto` | the enclosure probe: chosen per slice (`auto`) or fixed in Å |
| `--probe-radius` | 0 Å | subtracted from every reported radius |
| `--search-radius` | 6 Å | sideways search half-width (twice or four times as wide for a wide lumen) |
| `--lateral-exits` | off | accept a capped end if a side exit to bulk exists |
| `--exit-bulk-radius` | 6 Å | rolling probe that defines bulk solvent for side exits |

Atomic radii (`--radii`), figure text (`--annotate`), mouth guides and
structure input are [shared options](../reference/cli/common-options.md).
Every option is listed under [`crevice profile`](../reference/cli/profile.md).

## What you get

- **`-o profile.json`**: the sampled centres, radii, slice areas, each slice's
  enclosure probe, search half-width and grid spacing, chosen axis, seed,
  bottleneck, end clearances and every setting. Its `volume_estimate`
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

`profile` writes no cast. [`report`](publish.md) casts the same channel and
splits it into entry, lumen and exit segments,
[`profile --fingerprint`](fingerprint.md) adds lining hydrophobicity, chemistry, water
and a conductance estimate along the same channel, and
[`crevice view`](view.md) shows the centre line and cast in 3D.

:::{admonition} Interpreting results
:class: crevice-interpret

A profile is a measurement with the settings you chose: the numbers move a
little with the grid, the sample count and the input format (1GRM gives a
minimum of 1.337 Å with the quick-start settings and 1.337 to 1.338 Å with
0.25 Å grids and 81 to 201 samples). Report the enclosure-probe setting with
your results. The method assumes a
channel that runs roughly in one direction; branched or strongly curved
passages are not followed (use [`crevice paths`](paths.md) for those), openings narrower than the slice grid can be
missed, and two exits closer than twice the bulk probe are reported as one. If
no profile is found, that does not show the protein is closed: 4PYP, an
inward-open transporter, gives no through-profile, which fits its structure
but doesn't test it. The software is checked on synthetic walls with known
answers; the real-protein numbers in these pages are observations, not
validated biology.

[Limitations](../getting-started/limitations.md) collects the caveats shared by
all tools.
:::

See it in action in the [1GRM](../examples/1grm-profile-cast.md),
[3UKM](../examples/3ukm-lateral-exits.md) and
[1OED](../examples/1oed-pore-domain.md) tutorials.

# Curved and branched pathways

This note defines what `crevice paths` (Python: {func}`crevice.paths.find_paths`)
measures. The tool page is [Curved and branched pathways](../tools/paths.md).

A pathway is a curved centre line from a start (a seed point, a set of
residues, or the landmarks of a region definition) to bulk solvent, with the
radius of the largest atom-free sphere centred on the line measured along its
arc length. No channel axis is assumed. When several distinct routes leave the
same start they form a *pathway tree* that separates at *junctions*.

```{admonition} What the numbers are
:class: note
Pathways are geometric routes for the stated probe, bulk probe and grid. They
are not transport pathways, permeation mechanisms or functional assignments,
and the absence of a route on the grid does not prove that none exists.
```

## Grid and obstacles

| Item | Definition |
|---|---|
| Obstacles | Heavy atoms, as van der Waals spheres of the [radius set](atomic-radii.md) in effect. Non-HETATM atoms by default; `--include-hetero` adds hetero atoms (for MD input with `--md-selection 'not resname TIP3'`: lipids and ions). |
| Grid | Cubic lattice with world axes, spacing `--spacing` (default 0.5 Å), over the bounding box of the protein (non-hetero) atoms padded by the largest atom radius + the bulk-probe radius + three steps. Hetero atoms whose spheres can reach this box are kept. |
| Clearance | Exact distance from each node to the nearest atom surface. |
| Free node | Clearance greater than the probe radius: chosen automatically by default (see [Automatic probe](#automatic-probe)), or `--probe-radius R` for a fixed probe. The probe is a threshold only; reported radii are not reduced by it. |
| Step | Between free nodes in the 26 lattice directions. Diagonal steps are not tested analytically during the search (the clearance along a step of at most 0.87 Å falls below its endpoints by at most a few hundredths of an Å); every reported centre line is measured exactly afterwards. |

## Bulk, exposed nodes and the mouth

* **Bulk** is the rolling-probe definition used by `crevice cast` and the
  lateral exits of `crevice profile`: nodes with clearance ≥ `--bulk-radius`
  (default 6 Å) connected along atom-clear segments to the faces of the grid
  box are bulk-probe centres, and a point lies in bulk when it is within
  `--bulk-radius` of one. With membrane obstacles included, the box faces in
  the bilayer are blocked, so only solvent on either side of the membrane is
  bulk.
* **Exposed node.** A free node from which more than 26 − ⌈26 × `--min-enclosure`⌉
  (default 0.75, so at least 7) of the 26 straight lattice rays reach a bulk
  node before meeting a blocked node (clearance ≤ probe); a ray that leaves
  the grid first does not count. Exposed nodes therefore see bulk along
  several straight free lines, and they end a route like bulk. This stops
  routes from running along surface grooves outside the protein that the
  6 Å bulk probe cannot enter; `--min-enclosure 0` uses bulk only.
* **Mouth.** A route ends at the first bulk or exposed node it reaches. Each
  branch records which (`end_type`: `bulk` or `exposed`) and the mouth
  position (`exit_point_A`).

A route that never reaches bulk or an exposed node is not reported, and
nothing is interpolated.

## Widest routes

The start snaps to the free node of largest clearance within
`--seed-tolerance` (default 3 Å) of the requested point (ties: the nearest);
the snap distance is recorded. A maximin (widest-path) Dijkstra then labels
every reachable free node with its widest route from the start, the route
whose narrowest node clearance is largest; equally wide routes prefer the
shorter one, then lattice order, so results are deterministic. This is the
same ordering as `crevice tunnels` and the lateral-exit search. The route's
**ranking width** ignores nodes inside the start's own free sphere, which
every route shares; the reported bottleneck includes them.

**Sides.** With `--axis`, the search runs once per half-space, as the lateral
exits do: the `+` search may only use nodes on the `+` side of the plane
through the start normal to the axis (plus the start's free sphere), the `-`
search only the other side, each with its own distinct routes. A widest-path
search over all space can otherwise reach the far side through a wider detour
(a crevice around a selectivity filter, for example) instead of through the
narrower direct passage. A side with no route is recorded as
`unresolved_no_bulk_path_on_side`. Without `--axis` one search covers all
free space.

## Automatic probe

With `--probe-radius auto` (the default; Python `probe_radius="auto"`), the
probe is chosen per structure, and per side with `--axis`. A number instead
runs that probe alone, exactly as a fixed-probe run always has, and adds no
selection record.

1. **Candidates.** The probes of the automatic enclosure ladder of
   [`crevice profile`](../tools/profile.md): 0.50, 0.55, 0.60, 0.65, 0.70,
   0.75, 0.80, 0.90, 1.00 … 1.60, 1.80, 2.00, 2.20, 2.40, 2.70 and 3.00 Å (steps
   of 6–12 %), those below `--bulk-radius`. Below 0.5 Å the probe passes
   between atoms where the omitted hydrogens would be.
2. **One shared field.** The clearance grid and the bulk do not depend on the
   probe and are computed once. Exposure does (a ray stops at a node with
   clearance ≤ probe), but for each ray the largest probe at which it still
   reaches bulk is the smallest clearance met before bulk; one sweep per ray
   direction therefore gives, for every node, the largest ladder probe at
   which it is exposed. Each candidate then costs one route search only.
3. **One start.** The start snaps (as above) at the smallest candidate that
   resolves it and is shared by all candidates; a candidate the start node does
   not clear is recorded as `unresolved_seed_blocked`. A fixed-probe run at
   the chosen probe snaps to the same node unless a wider node within the
   seed tolerance was exposed at the smaller probe (not seen in the cases
   below).
4. **Pruning.** A route that clears a larger probe also clears the smallest
   one and reaches a mouth of the smallest one's search, so its width there
   exceeds the larger probe. Candidates at or above the widest route width
   found at the smallest probe on that side therefore cannot have a route;
   they are recorded as `pruned` and not searched.
5. **Grouping.** On each side, the route sets of the candidates (in
   ascending probe order) are grouped: two sets are the *same* when they have
   the same number of routes and match one to one, each pair being *not
   distinct* under the distinct-route rule below and having ranking widths
   within max(0.05 Å, 2 %) (the tolerance of `profile`'s automatic probe). A
   set joins the first group whose first member it matches.
6. **Choice.** The group found at the most candidate probes wins (at least
   two; ties go to the smaller probe), and its **smallest** probe is used. A
   route set seen at only one probe is a threshold effect: a passage just
   wider than that probe, or one that appears because a smaller probe makes
   nearby nodes exposed. When no set persists over two probes, the smallest
   probe with a route is used and the record says so. When no candidate has a
   route, the result is unresolved and lists every probe tried.
7. **Record.** `settings.probe_selection` holds the ladder, the rule, every
   candidate (probe, side, status, routes, their widths, mouths, end types and
   lengths, group, search time), the groups, the chosen probe per side
   (`chosen_by_side`) and the reason; `settings.probe_radius_A` is the
   smallest chosen probe and each branch records its own `probe_radius_A`.
   `*_path_probes.csv` has one row per candidate and side.
8. **Left-out routes.** On each side with a chosen probe, every route kept
   by another searched candidate is compared with the routes of the chosen
   result under the distinct-route rule below. A route that is not distinct
   from one of them is the same opening seen at another probe and is not
   listed. The others are *left out*: they exist at some candidate probe but
   are not in the result. Going through the candidates in ascending probe
   order, a left-out route joins the first listed one it is not distinct
   from, so a route found at several probes is listed once; the first
   (smallest-probe) instance supplies its width (ranking width), mouth, end
   type, lattice length and bottleneck (the node of smallest clearance
   outside the start's free sphere). Its lining residues use the lining shell
   below, on the unsmoothed lattice nodes: at the bottleneck node, and along
   its *own* nodes, those farther than one lattice diagonal (√3 × spacing)
   from every node of the result's routes on that side, in order from the
   start. `rerun_probe_A` is the largest probe at which the route was found:
   that candidate's route set contains it, so a fixed run at that probe
   reports it as a branch (up to the start-snapping caveat of step 3). The
   list is `settings.probe_selection.left_out_routes` (with
   `left_out_definition`), `chosen_by_side` counts them per side,
   `*_path_left_out.csv` has one row each, and the command prints how many
   were left out. Listing them changes nothing in the result: no centre
   line, radius profile or scene is computed for them, and the branches,
   their numbering and every other output are those of the chosen probe.

**Why the smallest probe of the most persistent set.** Within a group the
routes are the same, so the probe only changes what counts as wall and where
the mouths are. The smallest probe treats the fewest gaps as wall, gives the
most conservative mouths and the fewest lattice fallbacks, and matches
`profile`. The largest probe of the group sits just under the bottleneck,
where the lattice decides whether a node clears it. The smallest probe with
any route is not used either: small probes add narrow side routes and move
mouths inward through exposure, and these change from one probe to the next.
Persistence separates the two situations that a fixed probe cannot: in a
narrow filter, the filter routes are the only routes and persist up to their
width, so a small probe is chosen; in a wide vestibule, the wide routes
persist over many probes while narrow side routes come and go at the
smallest ones, so a larger probe is chosen and the side routes are left out.
They are listed as left-out routes (step 8), and a fixed `--probe-radius`
reports them.

**Per side.** With `--axis` the two half-space searches are already
independent (their own distinct routes, their own unresolved status), and the
two sides of a membrane protein are often very different: in the GLUT1 example
one side is a 2–3 Å vestibule and the other a gap of about 1 Å. One probe for
both would either drop the narrow side or add small-probe side routes to the
wide one, so each side gets its own probe; the start is shared. Without
`--axis` there is one search and one probe.

**Exposure is held fixed.** `--min-enclosure` (7 of 26 rays by default) is
not varied with the probe: it decides what counts as reaching the solvent,
not which gaps are passable, and varying both would make route sets at
different probes incomparable.

**Cost.** Wall-clock times on a shared 8-core workstation (October 2026,
load about 5 from other jobs), whole `find_paths` call:

| Structure (lattice) | Fixed probe | Automatic, 1 process | Automatic, 4 processes | 21 separate fixed-probe runs |
|---|---:|---:|---:|---:|
| 1GRM, seed in the pore (0.35 Å) | 1.2 s | 2.1 s | 1.6 s | 22 s |
| 3UKM, seed in the filter (0.35 Å, 13 M points) | 26 s | 41 s | 30 s | 416 s |
| GLUT1 MD frame 0, lipids as walls (0.5 Å, 5.2 M points) | 21 s | 116 s | 49 s | 400 s |
| GLUT1 MD frame 10 | 19 s | 94 s | 43 s | — |

Sharing the clearance field makes the ladder 3.5–10 times cheaper than
separate runs; the remaining cost is one route search per candidate and side
(on GLUT1 about 1–15 s each, larger probes being cheaper). Pruning skipped
8–16 searches per case but saved 0–19 % of the time, because the probes it
skips are the cheapest. Four processes run the GLUT1 ladder 2.4 times faster
than one; the clearance field (11–19 s) is computed once, before the
searches, and does not run in parallel. Each worker process used up to
0.45 GB.

## Distinct routes, junctions and segments

Mouths are visited widest route first (then shortest). A further route is
**distinct** from every route already kept when both hold:

1. its mouth is more than `--separation` (default 2 × bulk radius = 12 Å)
   from the other's mouth, so no single bulk-probe position touches both; and
2. beyond the *junction* of the two routes (the last lattice node they share),
   the parts of both routes farther than `--route-separation` (default the
   bulk radius, 6 Å) from that junction never come within `--route-separation`
   of each other: they leave through different openings.

For two routes that share only the start, this is exactly the distinct-exit
rule of the lateral exits ({func}`crevice.channel_exits.lateral_exit_legs`).
At most `--max-paths` (default 8) routes are kept per search.

All kept routes of a search are paths of one shortest-path tree, so they
share prefixes exactly. A **junction** is a lattice node where kept routes
separate (the start itself when they separate at once; `is_start`). The tree
is cut at the start, the junctions and the mouths into **segments**; a
**branch** is the whole route from the start to one mouth, listed widest
first.

## Centre line and radius

Each segment is resampled every lattice spacing and relaxed with fixed ends,
60 iterations by default (`--relax-iterations`):

1. a Laplacian step halfway toward the mean of the two neighbours (shortens
   lattice zig-zags);
2. a centring step of half a spacing times the clearance gradient, projected
   perpendicular to the line (moves the line toward the local clearance
   maximum, like HOLE's in-plane optimisation);
3. a move is accepted only if the interpolated clearance at the new position
   is at least the smaller of the old clearance and the lattice route's
   bottleneck, so the relaxation can never narrow the route below the
   bottleneck it was found with.

If the exact clearance of the relaxed line falls to the probe radius or below,
the segment keeps the straightened lattice route (straight cuts that are never
narrower than the stretch of route they replace); `centre_line` records
`relaxed`, `lattice_straightened` or `lattice_fallback`.

The **radius** is measured every 0.25 Å along each branch as the exact
distance to the nearest atom surface (a free-sphere radius, as in HOLE and
CAVER; not a planar section radius). The **bottleneck** is the sample of
smallest radius, with its position, arc length and the residues lining it.
`below_probe` is set when the measured bottleneck is at or below the probe
radius: the lattice route passed, but the continuous line touches the probe
(a grid effect worth checking with a finer `--spacing`).

## Lining residues

The lining shell of the [pore fingerprint](fingerprint.md): an atom lines a
sample when `|x − c| − r_atom − R ≤ --lining-shell` (default 1.4 Å), with `c`
the sample centre and `R` its radius; hydrogens, water and recognised
monatomic ions are excluded. A residue lines a branch when any of its atoms
lines one of its samples; its anchor is the arc length of its smallest gap.

## Unresolved cases

| `status` | Meaning |
|---|---|
| `resolved` | at least one branch |
| `unresolved_seed_blocked` | no free node within the seed tolerance (automatic probe: at any candidate) |
| `unresolved_seed_in_bulk` | every free node near the start is in bulk or exposed |
| `unresolved_no_bulk_path` | no route of free nodes from the start reaches a mouth (automatic probe: at any candidate; the reason gives the range tried) |

## Limitations

* **Lattice resolution.** A passage narrower than the probe radius plus up to
  half a lattice diagonal (0.87 × spacing) can be missed when no lattice node near its axis
  clears the probe. In 3UKM, the 0.71 Å constriction below the selectivity
  filter is passed at 0.35 Å spacing but not at 0.5 Å with a 0.5 Å probe. The
  automatic probe never goes below 0.5 Å, so run a finer `--spacing` before
  concluding that a route is closed.
* **Automatic probe.** The choice depends on the ladder, the grouping
  tolerances and the grid; narrow side routes seen only at small probes are
  left out when a wider route set persists. They are listed in
  `*_path_left_out.csv` with their width and probes, but only as lattice
  routes; give a fixed `--probe-radius` to measure what a particular probe
  finds. A route kept by none of the searched candidates (for example one
  beyond `--max-paths` at every probe) is not listed. The rule is not
  transitive: a route that a candidate keeps as distinct from that
  candidate's own widest route can still be *not* distinct from the chosen
  probe's route (mouths and junctions move with the probe), and is then not
  listed. In the GLUT1 example, secondary −*z* routes that a fixed 0.8 Å
  probe reports had mouths 11.5–13.7 Å from the chosen route's mouth, at
  the 12 Å threshold.
* **Maximin semantics.** The widest route is preferred however long it is;
  the half-space searches (`--axis`) and exposed nodes limit detours but do
  not forbid them. Read each branch's length and lining residues.
* **Mouth position** depends on the bulk-probe radius and the exposure rule;
  lengths near mouths are uncertain by a few Å.
* **One start.** The start is one point (seed, residue centroid or region
  landmark anchor); a different start in the same cavity can give a
  different tree.
* Membranes, crystal contacts and missing residues are only modelled as far
  as the obstacle atoms include them.

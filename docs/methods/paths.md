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
| Free node | Clearance greater than `--probe-radius` (default 0.8 Å, as `crevice tunnels`). The probe is a threshold only; reported radii are not reduced by it. |
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
| `unresolved_seed_blocked` | no free node within the seed tolerance |
| `unresolved_seed_in_bulk` | every free node near the start is in bulk or exposed |
| `unresolved_no_bulk_path` | no route of free nodes from the start reaches a mouth |

## Limitations

* **Lattice resolution.** A passage narrower than the probe radius plus up to
  half a lattice diagonal (0.87 × spacing) can be missed when no lattice node near its axis
  clears the probe. In 3UKM, the 0.71 Å constriction below the selectivity
  filter is passed at 0.35 Å spacing but not at 0.5 Å with a 0.5 Å probe. Run
  a finer `--spacing` or a smaller `--probe-radius` before concluding that a
  route is closed.
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

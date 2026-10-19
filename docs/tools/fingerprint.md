# Pore fingerprints

`crevice fingerprint` puts what CREVICE measures about a channel on one path
coordinate *s* (Å): the radius, how hydrophobic the residues lining the path
are, which kinds of residue line each part, the ENTRY / LUMEN / EXIT segments,
and, for an MD trajectory, the radius spread over frames and the water found
along the path. It also gives a geometric conductance estimate from the radius
profile.

Use it once `crevice profile` resolves the channel; it takes the same channel
options (`--axis`, `--origin`, `--enclosure-radius`, `--lateral-exits` ...).

The conductance is a geometric, order-of-magnitude estimate, not a
measurement: it treats the pore as a tube of bulk electrolyte (see
[Interpreting results](#fingerprint-interpreting)).

## Quick start

```bash
crevice fingerprint 1GRM --axis auto --out-dir 1GRM_fp --prefix 1GRM
crevice fingerprint 3UKM --assembly 1 --exclude-hetero --lateral-exits --enclosure-radius 0.55 \
    --out-dir 3UKM_fp --prefix 3UKM
crevice fingerprint glut1_excerpt.gro --trajectory glut1_excerpt.xtc \
    --selection 'protein and not name H*' --axis z --lateral-exits --enclosure-radius 0.55 \
    --out-dir glut1_fp --prefix glut1
crevice publish 1GRM --fingerprint --out-dir 1GRM_bundle     # the same files inside a publish bundle
```

```python
from crevice import load_structure
from crevice.fingerprint import fingerprint_structure, write_fingerprint_bundle

frame = load_structure(".crevice/pdb/1GRM.cif")
profile, cast, fp = fingerprint_structure(frame, axis="auto", electrolyte="kcl-1m")
print(fp.conductance["g_macro_pS"], [s["segment"] for s in fp.segments])
write_fingerprint_bundle(fp, "1GRM_fp", "1GRM")
```

For trajectories use {func}`crevice.fingerprint.fingerprint_trajectory`.

## How it works

**One coordinate.** The lumen is the channel profile itself and runs from
*s* = 0 to its length *L*. Past each mouth the path follows one centre line out
to bulk solvent: the straight vestibule line at an open end, or the widest side
exit at a capped end (`--lateral-exits`). Those samples sit at negative *s*
(lower end) or beyond *L* (upper end). ENTRY and EXIT follow the same rule as
`publish`'s cast segments (`--entry-end`).

**Lining residues.** An atom lines a sample when its surface is within
`--lining-shell` (1.4 Å, one water radius) of the largest atom-free sphere
there. Each residue is counted once, where it comes closest.

**Hydrophobicity.** Each lining residue takes its value from the Wimley–White
whole-residue octanol scale (Wimley, Creamer & White 1996; White & Wimley
1999), with the sign reversed so that positive means hydrophobic (kcal/mol),
and the values are smoothed along *s* with a Gaussian of width
`--hydrophobicity-bandwidth` (σ = 3 Å). Residues are also sorted into
hydrophobic, polar, charged and aromatic classes, counted per segment.

**Water** (trajectories). Each frame is fitted to the first, the radius is
measured again on the same path, and a water oxygen is counted at its nearest
sample when it sits inside that frame's free sphere. You get the mean count,
the linear and number densities and the occupancy along *s*.

**Conductance.** The lumen's radius profile gives the HOLE geometric estimate
(Smart et al. 1993, 1997), *G* = 1 / (ρ ∫ d*s* / π*r*²), with ρ the bulk
resistivity of `--electrolyte` (150 mM KCl by default). Also reported are the
value with the access resistance of both mouths added (Hall 1975) and HOLE's
empirically corrected value.

The [method note](../methods/fingerprint.md) gives every definition, the
electrolyte table, the formulas and the references.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--trajectory` | none | add a trajectory: radius spread over frames and water along the path |
| `--lining-shell` | 1.4 Å | how close an atom's surface must come to the free sphere to count as lining |
| `--hydrophobicity-bandwidth` | 3 Å | width σ of the Gaussian smoothing along *s* |
| `--electrolyte` | `kcl-150mm` | bulk resistivity for the conductance: `kcl-150mm`, `nacl-150mm`, `kcl-1m`, `nacl-1m` |
| `--resistivity` | none | any other bulk resistivity, in Ω m |
| `--entry-end` | `auto` | which end is labelled ENTRY |
| `--no-water` | off | with `--trajectory`, skip the water count |

The channel options (`--axis`, `--origin`, `--enclosure-radius`,
`--lateral-exits` ...) are those of [`profile`](profile.md), and every option
is listed under [`crevice fingerprint`](../reference/cli/fingerprint.md).

## What you get

| File | Contents |
|---|---|
| `PREFIX_fingerprint.csv` | One row per path sample: `s_A`, `segment`, `source`, `t_A`, `x_A/y_A/z_A`, `radius_A`, `raw_clearance_A`, `hydrophobicity_kcal_mol`, `hydrophobicity_support`, lining counts per class and the nearest lining residue. Trajectories add `radius_mean_A`, `radius_sd_A`, `radius_q_low_A`, `radius_q_high_A`, `water_count_mean`, `water_linear_density_per_A`, `water_number_density_per_A3` and `water_occupancy` |
| `PREFIX_fingerprint.json` | Conductance with its assumptions and citations, segment chemistry and residue lists, every lining residue, exit legs and parameters; for trajectories, per-frame conductance statistics and water definitions |
| `PREFIX_fingerprint_frames.csv` | Trajectories: time, fit RMSD, lumen minimum radius, *G*<sub>macro</sub>, *G*<sub>access</sub> and the number of lumen waters in each frame |
| `PREFIX_fingerprint.png` | Stacked panels against *s*: radius (frame mean and quantile band for trajectories), hydrophobicity (orange above zero, blue below), water (trajectories) and a chemistry strip with one row per class. ENTRY is shaded yellow, LUMEN teal and EXIT rust. Titles and the conductance call-out appear only with `--annotate` |

```{figure} ../examples/images/1grm_fingerprint.png
:alt: Pore fingerprint of gramicidin A: radius, hydrophobicity and a chemistry strip along the path coordinate, with ENTRY, LUMEN and EXIT shading.
:width: 560px

The 1GRM fingerprint from the quick start.
```

(fingerprint-interpreting)=
:::{admonition} Interpreting results
:class: crevice-interpret

A fingerprint is a way of reading one channel, not a verdict on what it
does. The hydrophobicity is assigned per residue, so a lumen lined only by
backbone carbonyls (as in gramicidin A) still reads as hydrophobic. Water
counts use free spheres and are lower bounds in wide or oval parts of the
pore, and in trajectories the path and lining residues are those of the
reference frame. The conductance treats the pore as a tube of bulk salt water:
for gramicidin A it comes out 10–40 times above the measured value (2–7 times
with HOLE's correction), so compare it between channels or frames as an order
of magnitude, never as a prediction of function. The
[method note](../methods/fingerprint.md#limitations) lists every assumption.

[Limitations](../getting-started/limitations.md) collects the caveats shared by
all tools.
:::

Try it in the [pore fingerprint tutorial](../examples/pore-fingerprint.md),
add it to a whole bundle with [`publish --fingerprint`](publish.md) (or
`--outputs default,fingerprint`, or `--all`), or look
at the same channel in 3D with [`crevice view`](view.md).

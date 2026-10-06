# Quick start

This page runs one short analysis from the command line and the same one from
Python. The examples use the gramicidin A channel
[1GRM](https://www.rcsb.org/structure/1GRM) and the human glucose transporter
GLUT1 [4PYP](https://www.rcsb.org/structure/4PYP). CREVICE can download both
(see [Inputs and selections](inputs-and-selections.md#fetching-structures)).

:::{note}
The numbers CREVICE prints depend on the grid spacing, probe radii, sample count
and input format you choose. They are measurements made with those settings,
not validated biological properties. See
[Validation levels](../tools/index.md#validation-levels).
:::

## 1. A channel radius profile

```bash
crevice fetch 1GRM --cache-dir .crevice/pdb
crevice profile 1GRM --axis auto -o 1GRM_profile.json --png 1GRM_profile_radius.png
```

`--axis auto` tries the three principal directions of the molecule. It accepts
a direction only if it finds a connected, laterally enclosed passage with two
unobstructed ends. The enclosure probe (the sphere that decides which wall gaps
count as closed) is also chosen automatically, from a 0.5–3.0 Å ladder. The
command prints a one-line summary and the chosen probe:

```text
profile points=81 min_radius=1.324 mean_radius=1.563 bottleneck_index=73
enclosure_radius=0.750 (auto: smallest probe of channel group 1, which resolved at 5 of 21 probes tried)
```

It writes:

- `1GRM_profile.json`: every sampled centre, radius, section area and the
  settings used;
- `1GRM_profile.csv`: one row per sample (`radius_A`, `axis_position_A`, ...);
- `1GRM_profile_radius.png`: radius against the axial coordinate.

If no unique through-channel is found, the command stops with an error. It does
not make up a profile. That result does not mean the protein is closed. See
[Outputs and unresolved results](outputs.md).

## 2. A 3D cavity cast

```bash
crevice fetch 4PYP --cache-dir .crevice/pdb
crevice cast 4PYP --out-dir 4PYP_cast
```

`cast` fills the dominant interior cavity with atom-clear probe spheres. It
writes a measured binary map (`4PYP_volume.dx`) and a separate smoothed display
map. It also writes PyMOL, VMD and ChimeraX scene scripts and a manifest that
lists every file. It reports one region by default. Use `--max-cavities N` to
allow up to N regions that qualify on their own; this cap never forces extra
regions. The [cast guide](../tools/cast.md) explains the method and its limits.

## 3. The same profile from Python

```python
from crevice import annotate_residues, load_structure, pore_profile

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto", samples=81)
contacts = annotate_residues(frame, profile.points)

print(profile.min_radius, profile.bottleneck.nearest_residue)
for contact in contacts[:5]:
    print(contact.residue, contact.role, round(contact.min_distance, 2))
```

The Python functions return result objects and do not write files. Writers such
as {func}`crevice.write_profile_json` and plotting functions such as
{func}`crevice.plot_profile_radius` are called separately. See
[CLI or Python?](cli-vs-python.md).

## 4. Everything in one bundle

```bash
crevice publish 1GRM --out-dir results/1GRM --prefix 1GRM
```

`publish` runs the profile, cast, residue, cavity, tunnel, network and hydration
analyses and writes CSV tables, figures and viewer scenes (about 85 files for
1GRM). Use the `--skip-*` options to leave out analyses you do not need. See
[Outputs](outputs.md).

## Next steps

- One page per command: the [tool guide](../tools/index.md).
- Trajectory analysis: [`trajectory`](../tools/trajectory.md),
  [`cavity-trajectory`](../tools/cavity-trajectory.md) and
  [`hydration`](../tools/hydration.md).
- Worked examples: [Examples](../examples/index.md).
- Definitions and limitations: the [method notes](../methods/index.md).

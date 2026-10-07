# Basic usage

This page gets you from an installed package to your first results, first on
the command line and then in Python. We use two public structures that CREVICE
can download for you: the gramicidin A channel
[1GRM](https://www.rcsb.org/structure/1GRM) and the human glucose transporter
GLUT1 [4PYP](https://www.rcsb.org/structure/4PYP).

## On the command line

Every analysis is a sub-command of `crevice`. Give it a structure (a file, a
PDB accession or an AlphaFold DB identifier) and somewhere to write the
results. `crevice <command> --help` lists every option with its units and
default.

### Measure a channel

```bash
crevice profile 1GRM -o 1GRM_profile.json --png 1GRM_profile_radius.png
```

CREVICE downloads 1GRM into `.crevice/pdb/`, finds the channel axis, picks a
suitable enclosure probe and measures the radius of the largest atom-free
sphere at each point along the channel. It prints a one-line summary:

```text
profile points=81 min_radius=1.324 mean_radius=1.563 bottleneck_index=73
enclosure_radius=0.750 (auto: smallest probe of channel group 1, which resolved at 5 of 21 probes tried)
```

and writes three files:

- `1GRM_profile.csv`: one row per point, with the radius (`radius_A`) and the
  position along the axis (`axis_position_A`);
- `1GRM_profile.json`: the full record, including every setting that was used;
- `1GRM_profile_radius.png`: the radius plotted along the channel.

```{image} ../examples/images/1grm_profile.png
:alt: Pore radius of gramicidin A along the channel axis, with the bottleneck marked.
:width: 520px
```

If CREVICE cannot find a single clear channel that runs right through the
protein, it stops and tells you why rather than guessing. That is common for
transporters and enzymes, and it does not mean the protein is closed; a cast
(below) is usually the next thing to try. See
[Outputs and unresolved results](outputs.md).

### Cast a cavity

```bash
crevice cast 4PYP --out-dir 4PYP_cast
```

`cast` fills the largest interior space of GLUT1 with atom-clear probe
spheres. You get the cast as a 3D map (`4PYP_volume.dx`), its volume, ready-made
PyMOL, VMD and ChimeraX scenes, and a `*_manifest.json` listing every file.
Ask for more than one region with `--max-cavities N`. The
[cast guide](../tools/cast.md) explains how regions are chosen.

### Everything at once

```bash
crevice publish 1GRM --out-dir results/1GRM --prefix 1GRM
```

`publish` runs the profile, cast, lining residues, cavities, tunnels, residue
network and hydration in one go and writes the tables, figures and viewer
scenes (about 85 files for 1GRM). Use the `--skip-*` options to leave out
what you don't need.

## In Python

The same analyses are plain functions. They take a loaded structure and return
result objects; nothing is written to disk until you ask for it.

```python
from crevice import annotate_residues, load_structure, pore_profile

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto", samples=81)
contacts = annotate_residues(frame, profile.points)

print(profile.min_radius, profile.bottleneck.nearest_residue)
for contact in contacts[:5]:
    print(contact.residue, contact.role, round(contact.min_distance, 2))
```

To save results, call a writer or a plotting function:

```python
from crevice import plot_profile_radius, write_profile_csv, write_profile_json

write_profile_json(profile, "1GRM_profile.json")
write_profile_csv(profile, "1GRM_profile.csv")
plot_profile_radius(profile, "1GRM_profile_radius.png")
```

Trajectories work the same way: {func}`crevice.load_trajectory` reads a
topology and trajectory through MDAnalysis, and
{func}`crevice.analyze_trajectory` runs the analyses frame by frame.

## A few things worth knowing

- **Settings matter.** Radii and volumes depend on the grid spacing, probe
  radii, atom selection and atomic radius set. CREVICE records all of them in
  its output; report them alongside your numbers.
- **Figures are clean by default.** Add `--annotate` (or `annotate=True`) for
  titles, labels and legends.
- **CLI and Python defaults occasionally differ.** The differences are listed
  in [CLI or Python?](cli-vs-python.md#defaults-that-differ-between-the-cli-and-python).

## Where next?

```{toctree}
:maxdepth: 1

inputs-and-selections
outputs
cli-vs-python
```

- The [Analysis tools](../tools/index.md) explain each method and its options.
- The [Tutorials](../examples/index.md) work through real structures from
  start to finish.

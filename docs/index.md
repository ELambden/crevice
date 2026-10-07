# CREVICE

::::{container} crevice-hero

:::{container}
**CREVICE** finds and measures the spaces inside proteins: the pores and
channels that run through them, the cavities buried within them and the
tunnels that lead out. It tells you which residues line those spaces, which
residues they touch, and where the water is, for a single structure or across a
whole molecular-dynamics trajectory.

```bash
pip install crevice
crevice publish 1GRM --out-dir results/1GRM
```
:::

![A channel through the TWIK-1 potassium channel (PDB 3UKM), cast by CREVICE: yellow entry vestibule, teal lumen and two rust-coloured side exits.](examples/images/3ukm_cast_exits_chimerax.png)

::::

{{ full_name }} is a Python package and a command-line tool. Everything it
does can be run as a `crevice` command or called from Python, and every result
comes as tables you can analyse, figures you can publish and scenes you can
open in PyMOL, VMD or ChimeraX.

## Your first analysis

Measure how wide the gramicidin A channel is along its length. CREVICE
downloads the structure for you:

```bash
crevice profile 1GRM -o 1GRM_profile.json --png 1GRM_profile.png
```

```text
profile points=81 min_radius=1.324 mean_radius=1.563 bottleneck_index=73
enclosure_radius=0.750 (auto: smallest probe of channel group 1, which resolved at 5 of 21 probes tried)
```

You now have the radius at every point along the channel in
`1GRM_profile.csv`, the full record in `1GRM_profile.json` and a plot of the
profile. The same thing in Python:

```python
import crevice

frame = crevice.load_structure(".crevice/pdb/1GRM.cif")
profile = crevice.pore_profile(frame)
print(f"narrowest point: {profile.min_radius:.2f} Å")
```

[Basic usage](getting-started/basic-usage.md) walks through this and a few
more first steps.

## What can it do?

::::{container} crevice-cards

:::{container} crevice-card
![](examples/images/1grm_profile.png)

**[Pores and channels](tools/profile.md)**

Radius profiles along a channel and 3D casts of the space inside it.
:::

:::{container} crevice-card crevice-violet
![](examples/images/4pyp_residue_context_chimerax.png)

**[Cavities and tunnels](tools/cast.md)**

Casts of pockets and buried cavities, and the widest routes to the outside.
:::

:::{container} crevice-card crevice-entry
![](examples/images/4pyp_residue_evidence.png)

**[Residues and networks](tools/residue-evidence.md)**

Which residues form a boundary, and who they talk to.
:::

:::{container} crevice-card crevice-exit
![](examples/images/1grm_nmr_profiles.png)

**[Trajectories](tools/trajectory.md)**

Follow a channel, a cavity or its water through every frame of a simulation.
:::

::::

Browse all of them in [Analysis tools](tools/index.md), or learn by example in
the [Tutorials](examples/index.md).

## Install

```bash
pip install crevice              # the core package
pip install "crevice[science]"   # plus every optional extra
```

CREVICE needs Python 3.10 or newer. [Installation](getting-started/installation.md)
covers the extras, installing from source and the molecular viewers.

## Getting help

Questions, ideas and bug reports are all welcome as
[GitHub issues](https://github.com/ELambden/crevice/issues). If something
looks broken, please include the CREVICE version,
the command you ran and, if you can, the structure (an accession is perfect)
and the `*_manifest.json` from the run.

## Citing CREVICE

If CREVICE helps your research, please [cite it](project/citing.md) and say
which version and settings you used.

```{toctree}
:hidden:

Home <self>
getting-started/installation
getting-started/basic-usage
tools/index
examples/index
reference/cli/index
api/index
project/citing
project/changelog
```

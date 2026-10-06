# CREVICE documentation

**CREVICE** ({{ full_name }}) is a Python package
and command-line tool for measuring pores, channels, tunnels and cavities in
protein structures and molecular-dynamics (MD) trajectories. It also records the
residues that line those spaces, the residues they touch, and the water around
them.

:::{admonition} Geometric measurements, not biological validation
:class: important

These pages describe `crevice` {{ version }}, an early release. Where a page
reports a numerical result from a real protein, it describes an observation made
with the stated settings. CREVICE's benchmark systems are not curated, so none of
these numbers is a validated biological result, and a cavity, path or residue
ranking found by CREVICE is a geometric candidate until it is checked against
independent evidence. See [Validation levels](tools/index.md#validation-levels).
:::

## What CREVICE measures

- **Through-channel radius profiles** along an automatically proposed or
  user-supplied axis, with continuous atom-clearance checks
  ([`profile`](tools/profile.md)).
- **3D void casts** of pockets, cavities and channels, filled with atom-clear
  probe spheres and exported as measured binary maps ([`cast`](tools/cast.md),
  [`publish`](tools/publish.md)).
- **Cavity and tunnel searches** on clearance grids
  ([`cavities`](tools/cavities.md), [`tunnels`](tools/tunnels.md)).
- **Residue evidence**: which residues form a measured boundary, and which other
  residues contact them ([`residues`](tools/residues.md),
  [`residue-evidence`](tools/residue-evidence.md), [`network`](tools/network.md)).
- **Trajectory statistics**: aligned per-frame volumes, widths, contacts and
  motion, with approximate uncertainty on means
  ([`trajectory`](tools/trajectory.md), [`cavity-trajectory`](tools/cavity-trajectory.md)).
- **Hydration**: explicit water contacts, solvent-accessible area, water density
  and typed interactions ([`hydration`](tools/hydration.md)).
- **Named regions**: versioned, reviewable region definitions for trajectory
  analysis ([`region-*`](tools/regions.md)).

Every analysis runs from the `crevice` command line or from Python.

## Where to start

- New users: [Installation](getting-started/installation.md), then
  [Quick start](getting-started/quickstart.md).
- Choosing between the command line and Python:
  [CLI or Python?](getting-started/cli-vs-python.md)
- Preparing structures and choosing atoms:
  [Inputs and selections](getting-started/inputs-and-selections.md)
- What a run writes, and what "unresolved" means:
  [Outputs and unresolved results](getting-started/outputs.md)
- Worked examples on public structures (1GRM, 3UKM, 1OED, 4PYP, an NMR
  ensemble as a small trajectory): [Examples](examples/index.md).
- Method detail: the [tool guide](tools/index.md) and the
  [method notes](methods/index.md).

```{toctree}
:hidden:
:caption: Getting started

getting-started/installation
getting-started/quickstart
getting-started/cli-vs-python
getting-started/inputs-and-selections
getting-started/outputs
```

```{toctree}
:hidden:
:caption: Tool guide

tools/index
```

```{toctree}
:hidden:
:caption: Method notes

methods/index
```

```{toctree}
:hidden:
:caption: Examples

examples/index
```

```{toctree}
:hidden:
:caption: Reference

reference/cli
api/index
api/modules
```

```{toctree}
:hidden:
:caption: Project

development/contributing
development/documentation
development/releasing
project/changelog
project/citing
project/community
```

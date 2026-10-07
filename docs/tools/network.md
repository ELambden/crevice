# Residue networks

`crevice network` builds a map of which residues touch each other in a
structure. By default it adds the channel itself as a node, so you can see how
the residues lining the channel connect to the rest of the protein. It can also
trace the shortest contact paths outwards from the lining residues.

## Quick start

```bash
crevice network 1GRM -o 1GRM_network.json --png 1GRM_network.png --chord-png 1GRM_chord.png
```

```python
from crevice import (annotate_residues, build_cavity_network, load_structure,
                     network_metrics, pore_profile)

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto")
contacts = annotate_residues(frame, profile.points)
network = build_cavity_network(frame, profile, contacts=contacts,
                               distance_metric="center", include_hetero=True)
print(network_metrics(network)["density"])
```

The Python functions default to the surface-gap distance and leave
heteroatoms out, so pass `distance_metric="center", include_hetero=True` as
above to match the command.

## How it works

- **Nodes** are residues. With `--region profile` (the default) the channel is
  added as a node and linked to the residues that line it or form its
  bottleneck, as found by [`residues`](residues.md). This needs a resolved
  profile; `--region none` gives a residue-only network.
- **Edges** join residues within `--cutoff` (4.5 Å). The command measures the
  closest heavy-atom centres by default (`--contact-metric center`); `surface`
  measures the gap between atom surfaces instead.
- **Edge labels** sort contacts into simple classes by distance and residue
  type: `salt_bridge_candidate` (charged atoms within 4.0 Å),
  `oppositely_charged_residue_contact`, `hydrophobic_residue_contact`,
  `polar_residue_contact` and `distance_contact`. They don't assign hydrogen
  bonds or energies; for typed chemical interactions, see
  [Hydration](hydration.md).
- **Paths.** `--connectivity-json` and `--connectivity-png` find weighted
  paths outwards from the lining residues (needs the `networks` extra).
  `--remove-residue` repeats the search without a residue, and
  `--residue-groups groups.json` lets you summarise by helix or strand.

In Python, the region you pass to {func}`crevice.build_cavity_network`
decides which points residues are measured against: a profile or tunnel gives
its path points, a cast {class}`~crevice.models.VoidComponent` its sample
points, and a {class}`~crevice.models.Cavity` only its centre.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--cutoff` | 4.5 Å | contact distance |
| `--contact-metric` | `center` | `center` (heavy-atom centres) or `surface` (gap between atom surfaces) |
| `--region` | `profile` | add the channel as a node, or `none` |
| `--png`, `--chord-png` | none | contact bar chart and chord diagram |
| `--connectivity-json`, `--connectivity-png` | none | paths from the lining residues |

All options are listed under [`crevice network`](../reference/cli/network.md).

## What you get

- **`-o network.json`** with **`network_nodes.csv`** and
  **`network_edges.csv`** beside it: every residue with its position, group and
  number of contacts, and every contact with its class, distance and closest
  atoms.
- **`--png`**: the 20 most-connected residues, with bars split by contact
  class (purple salt-bridge candidate, orange oppositely charged, green
  hydrophobic, teal polar, grey other). The channel node isn't counted here.
- **`--chord-png`**: a chord diagram of the 24 most-connected nodes, with the
  channel in red and residues in blue, and edges coloured by class.
- **`--connectivity-json`/`--connectivity-png`**: the paths as tables, and a
  figure placing residues in columns by how many contacts away from the lining
  they are.

:::{admonition} Interpreting results
:class: crevice-interpret

A contact means two residues are close in this one structure; it doesn't show
that the interaction is stable or favourable (follow it over a
[trajectory](trajectory.md) for that). Network paths and centrality are
associations, not evidence of signalling or allostery. One practical catch:
networks written by `publish` and `analyze` use the Python defaults (surface
gap, no heteroatoms), not the `crevice network` defaults, so compare networks
only when they were built the same way. Each network records its settings in
`metadata`.
:::

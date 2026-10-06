# `network`: residue contact networks

## Scientific question

Which residues touch each other, and how do the residues that line a channel
connect to the rest of the protein? `network` builds a residue contact graph
from one static structure. By default it adds a pseudo-node for the channel
profile.

## Method

- **Nodes** are residues. With `--region profile` (the default), a channel
  pseudo-node is added, and `region_lining` / `region_bottleneck` edges join it
  to residues from [`residues`](residues.md). This needs a resolved profile. Use
  `--region none` for a residue-only graph.
- **Edges** join residue pairs within `--cutoff` (4.5 Å). The CLI uses the
  minimum heavy-atom **centre** distance by default (`--contact-metric
  center`). `surface` uses the van der Waals surface gap.
- **Edge labels** are simple distance and residue-type categories:
  `salt_bridge_candidate` (charged atoms ≤ 4.0 Å),
  `oppositely_charged_residue_contact`, `hydrophobic_residue_contact`,
  `polar_residue_contact` and `distance_contact`. These labels do not assign
  hydrogen bonds or energies. Directional explicit-H chemistry is part of the
  [hydration workflow](hydration.md).
- `--residue-groups groups.json` maps residue IDs to helix or strand names so
  the graph can be summarised by group. `--connectivity-json` and
  `--connectivity-png` find weighted graph paths from lining residues. These
  need NetworkX. `--remove-residue` repeats the path search with a residue
  removed.

## Assumptions

- A contact is proximity in one structure. It does not show a persistent or
  energetically favourable interaction. For persistence across frames, see
  [`trajectory`](trajectory.md) and [`cavity-trajectory`](cavity-trajectory.md).
- Graph paths and centrality are **associations**. They do not show signal
  transmission, allostery or causal information flow.

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--cutoff` | 4.5 | Å | contact distance |
| `--contact-metric` | `center` | | `center` (heavy-atom centres) or `surface` (VDW gap) |
| `--region` | `profile` | | `profile` adds the channel pseudo-node; `none` omits it |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

## Outputs

- `-o network.json`: nodes, edges with labels and distances, and simple graph
  metrics (`density`, degree summaries).
- `network_nodes.csv`, `network_edges.csv` (beside `-o`, same stem plus
  `_nodes`/`_edges`): one row per node (`node_id`, `kind`, residue identity,
  `group`, `properties`, centroid `x_A`/`y_A`/`z_A`, `degree` counting all
  edges and `residue_degree` counting residue-residue edges only) and one row
  per edge (`source`, `target`, `interaction`, `distance_A`,
  `distance_metric`, `cutoff_A`, `weight`, closest atoms, centre and surface
  distances, charged atom pairs, sequence separation).
- `--png`: the 20 residues with most residue-residue contacts, as horizontal
  bars stacked by contact class (purple salt-bridge candidate, orange
  oppositely charged, green hydrophobic, teal polar, grey other distance
  contact; legend below the axes). The x axis states the metric, for example
  "Residue contacts (degree: residues with atom-centre distance ≤ 4.5 Å)". The
  region pseudo-node and its edges are not counted here: `degree` in the node
  table includes them, and the [`residues`](residues.md) figure ranks the
  residues lining the region.
- `--chord-png`: chord diagram of the 24 highest-degree nodes: region node red,
  residues blue, node names around the circle; edges coloured by interaction
  class (red channel/cavity contact, purple salt-bridge candidate, orange
  oppositely charged, green hydrophobic, teal polar, grey other distance
  contact), with a legend of the classes drawn, width scaled by weight.
- `--connectivity-json`, `--connectivity-png`: weighted lining-residue paths
  (the tables are also written beside the JSON as `<stem>.csv` and
  `<stem>_groups.csv`);
  the figure places residues in columns by contact hops from the lining
  (squares = lining residues), coloured by `--residue-groups` group, with
  residue names beside the markers.

Figures have no titles unless `--annotate` is given.

## Python equivalent

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

The low-level defaults are `distance_metric="surface"` and
`include_hetero=False`. Pass the values above to match `crevice network`.

:::{warning}
The networks written by `publish` and `analyze` use the low-level defaults
(surface gap, heteroatoms excluded), not the `crevice network` defaults.
Compare networks only when they were produced with the same metric and
obstacle set. Unifying these defaults is an open interface-design decision;
until it is made, the difference is intentional and recorded in each network's
`metadata` (`distance_metric`, `include_hetero`).
:::

In Python, the region passed to {func}`crevice.build_cavity_network` sets the
contact points: a profile or tunnel contributes its path points with their
clearance spheres, a cast {class}`~crevice.models.VoidComponent` its sample
points, and a {class}`~crevice.models.Cavity` only its centre point (so only
residues within the cutoff of the centre are linked).

## Testing and validation status

- **Software / synthetic:** contact metrics, label rules and path removal are
  covered by unit tests.
- **Example observation (public entries):** explicit 4.5 Å centre-contact
  networks were built for 1GRM and 4PYP; 4PYP had seven typed salt-bridge
  candidates in the separate typed-interaction analysis.
- **Biological / functional:** not established.

## Known limitations

- The sequence-separation and nonlocal rules used elsewhere (for example in
  `residue-evidence`) are not applied here.
- Group labels are user-supplied. Secondary structure is not assigned
  automatically by this command.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: network
```

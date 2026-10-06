# `residues`: residues along a channel profile

## Scientific question

Which residues line a through-channel, and which sit at or near its narrowest
point? `residues` first computes a [`profile`](profile.md), then records every
residue whose atoms come within a surface-distance cutoff of the sampled
channel.

For the boundary of a measured 3D cast rather than a 1D profile, use
[`residue-evidence`](residue-evidence.md).

## Method

Each profile point defines a local sphere whose radius is its atom clearance.
For each residue, CREVICE records the smallest and mean distance from its atom
surfaces to those spheres, and the profile points it touches. Roles:

| Role | Rule |
|---|---|
| `bottleneck` | minimum distance ≤ 1.0 Å, and it touches a point within two samples of the narrowest point |
| `lining` | minimum distance ≤ 1.0 Å, not at the bottleneck |
| `bottleneck-nearby` | within the cutoff and near the bottleneck, but further than 1.0 Å |
| `nearby` | within the cutoff only |

`influence_score` combines inverse distance, the number of points touched and a
role bonus. It is a **ranking aid for geometric contact**. It does not measure
functional importance. Residue labels also get simple chemical property tags
(for example `hydrophobic`, `small`).

## Assumptions

- The profile resolves. When it does not, the command fails in the same way
  `profile` does.
- Profile options (`--axis`, `--origin`, `--samples`, ...) apply as in
  [`profile`](profile.md).

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--cutoff` | 4.5 | Å | residue atom-surface distance cutoff |
| profile options | see [`profile`](profile.md) | | |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

## Outputs

- `-o residues.csv`: one row per residue: `residue, min_distance_A,
  mean_distance_A, atom_count, point_count, role, properties,
  influence_score`. The distances are the smallest and mean atom-surface gap
  to the profile's region surface (Å).
- `--json`: the same, as JSON (optional).
- `--png`: a ranked contact figure: influence score bars for the top 20
  residues (residue IDs as tick labels), coloured by role, with a legend below
  the axes naming the roles drawn: red `bottleneck` (gap ≤ 1 Å within two
  samples of the narrowest sample), orange `bottleneck-nearby` (gap > 1 Å
  there), blue `lining` (gap ≤ 1 Å elsewhere), grey `nearby` (gap > 1 Å
  elsewhere, within `--cutoff`). The x axis is "Influence score (higher =
  closer contact along more of the path)". No title unless `--annotate`.

## Python equivalent

```python
from crevice import annotate_residues, load_structure, pore_profile, write_residue_contacts_csv

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto")
contacts = annotate_residues(frame, profile.points, cutoff=4.5)
write_residue_contacts_csv(contacts, "residues.csv")
```

## Testing and validation status

- **Software / synthetic:** role and distance rules are covered by unit tests.
- **Example observation (public entry):** used on the 1GRM profile.
- **Biological / functional:** not established. A lining or bottleneck role
  does not show that a residue controls permeation.

## Known limitations

- The 1.0 Å role threshold and the score weights are fixed design choices. They
  were not calibrated.
- This command only works along a resolved profile.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: residues
```

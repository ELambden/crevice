# Lining residues

`crevice residues` tells you which residues line a channel and which ones sit
at its narrowest point. It measures a [profile](profile.md) first, then lists
every residue whose atoms come within a cutoff of the channel, with how close
it gets and how much of the channel it touches.

Working with a 3D cast rather than a channel profile? Use
[boundary residues](residue-evidence.md) instead.

## Quick start

```bash
crevice residues 1GRM -o 1GRM_residues.csv --png 1GRM_residues.png
```

```python
from crevice import annotate_residues, load_structure, pore_profile, write_residue_contacts_csv

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto")
contacts = annotate_residues(frame, profile.points, cutoff=4.5)
write_residue_contacts_csv(contacts, "1GRM_residues.csv")
```

## How it works

Each profile point is the centre of a sphere as wide as the channel there.
For every residue, CREVICE measures the smallest and mean gap between its atom
surfaces and those spheres, and counts the points it touches. Each residue
then gets a role:

| Role | Rule |
|---|---|
| `bottleneck` | within 1.0 Å, touching a point within two samples of the narrowest point |
| `lining` | within 1.0 Å, elsewhere along the channel |
| `bottleneck-nearby` | within the cutoff and near the bottleneck, but further than 1.0 Å |
| `nearby` | within the cutoff only |

An `influence_score` combines closeness, the number of points touched and the
role, to help you rank residues by how much they contact the channel. Residues
also get simple chemical tags such as `hydrophobic` or `small`.

The profile options (`--axis`, `--origin`, `--samples` and so on) work exactly
as in [`profile`](profile.md); if the profile doesn't resolve, neither does
this command.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--cutoff` | 4.5 Å | how close a residue's atom surfaces must come to be listed |
| `--json` | none | also write the table as JSON |
| `--png` | none | a ranked contact figure |

Channel profile settings, atomic radii and structure input are
[shared options](../reference/cli/common-options.md); every option is listed
under [`crevice residues`](../reference/cli/residues.md).

## What you get

- **`-o residues.csv`**: one row per residue with `min_distance_A`,
  `mean_distance_A`, `atom_count`, `point_count`, `role`, `properties` and
  `influence_score`.
- **`--png`**: the top 20 residues as bars of influence score, coloured by
  role: red bottleneck, orange bottleneck-nearby, blue lining, grey nearby.

For a picture of where along the channel each residue sits, `profile
--annotated-png` draws the radius above a lane of nearest residues:

```{figure} ../examples/images/1grm_profile_residues.png
:alt: Gramicidin A pore radius above a lane of nearest residues.
:width: 520px

Gramicidin A (1GRM): the nearest residue every fifth point, coloured by
chemistry (purple charged, teal polar, blue other).
```

:::{admonition} Interpreting results
:class: crevice-interpret

The roles describe geometric contact with the measured channel. The 1.0 Å
threshold and the score weights are fixed design choices, not calibrated
against experiment, so a `bottleneck` residue is one that sits at the
narrowest point, not one shown to control permeation. Use the list to choose
residues worth testing, for example by mutagenesis.
:::

# Hydration

`crevice hydration` looks at the water around the residues you care about.
It tells you which residues touch explicit water and how often, how exposed
each residue is, where water sits on average once the protein is aligned, and
which residue pairs make hydrogen bonds, salt bridges, aromatic contacts or
water bridges.

You often don't need to run it separately: `analyze`, `cast`, `publish`,
`residue-evidence`, `trajectory`, `cavity-trajectory` and `region-trajectory`
all include the same analysis unless you pass `--skip-hydration`. The
standalone command is handy when you want to add hydration to geometry you
have already computed.

## Quick start

For a static structure with crystallographic waters:

```bash
crevice hydration 4PYP --out-dir 4PYP_hydration --prefix 4PYP
```

For a trajectory, reusing a finished [cavity trajectory](cavity-trajectory.md)
so the same residues, frames and alignment are used:

```bash
crevice hydration topology.gro --trajectory production.xtc \
    --cavity-results cavity/crevice_cavity_statistics.json --out-dir hydration
```

```python
from crevice import load_structure
from crevice.hydration import static_hydration
from crevice.hydration_export import write_hydration_bundle

frame = load_structure(".crevice/pdb/4PYP.cif")
result = static_hydration(frame, cutoff=3.5, sasa_probe=1.4, sasa_samples=256)
write_hydration_bundle(result, "4PYP_hydration", prefix="4PYP")
```

For trajectories in Python, use
{func}`crevice.hydration_trajectory.analyze_hydration_trajectory`; the whole
command, with residue focus, density, interactions and scenes, is assembled in
{mod}`crevice.hydration_workflow` (some of its helpers still take command-line
argument objects).

## What it measures

- **Water contacts.** A water touches a residue when its oxygen is within
  `--hydration-cutoff` (3.5 Å) of any of the residue's heavy atoms. Each water
  counts once per residue per frame, and *occupancy* is the fraction of frames
  with at least one contact.
- **Solvent-accessible area (SASA).** Shrake–Rupley, with 256 points per atom
  and a 1.4 Å probe. Every non-water heavy atom, including lipids, ligands and
  ions, blocks the probe. Areas are reported per residue and for backbone,
  side-chain and N/O atoms; these subsets overlap, so don't add them up.
- **Water density.** Aligned 3D maps of where water oxygens sit, in water
  oxygens per Å³, with an absolute isovalue for display.
- **Typed interactions.** Directional hydrogen bonds (explicit hydrogens
  needed), salt-bridge and aromatic candidates and water-mediated bridges,
  plus simplified DSSP (with the `hydration` extra) and side-chain torsions.
  Residues whose chemistry can't be assigned are listed and left out.
- **Over a trajectory**, CREVICE streams the solvent frame by frame, follows
  each water by its topology index, uses periodic minimum-image distances
  (orthogonal or triclinic boxes) and aligns the protein on matching CA atoms.
  Every selected frame is used, and the frame-to-frame range is reported
  separately from an approximate confidence interval for the mean, which is
  withheld when sampling is too thin. Water persistence is measured between
  saved frames, so it is not a continuous residence time.

The full definitions are in
[Residue hydration and solvent exposure](../HYDRATION_METHODS.md) and
[Coordinated hydration and residue evidence](../HYDRATION_NETWORK_VIEWS.md).

## Choosing which residues to look at

| Input | Residues analysed |
|---|---|
| `--cavity-results stats.json` (from `cavity-trajectory`) | every residue that was ever on the cavity wall or a partner, using that run's frames, alignment and reference |
| `--volume-dx map.dx` | the wall residues and partners of a measured cast |
| `--residues-json ids.json` | a list of residues you choose |
| none | every selected protein residue |

Only residues with non-HETATM heavy atoms are measured. A focus residue
without any, such as a ligand, ion or water, or the D-amino acids of 1GRM when
read from a PDB file, is listed in the JSON with the reason
(`settings.untyped_focus_residues`) and still blocks SASA unless it is water.
Residue IDs that aren't in the structure are an error.

## Water in the cavity: two ways to count

For trajectory regions there are two water counts, and they answer different
questions:

- **Fixed-reference count** (always on): water oxygens within 2 Å of the
  cavity as it was in the reference. The neighbourhood never moves, so protein
  moving in or out of it changes the count even if no water moves.
- **Instantaneous membership** (`--water-membership` on
  [`cavity-trajectory`](cavity-trajectory.md) or
  [`region-trajectory`](regions.md)): water oxygens that sit inside the cavity
  as measured in that same frame. Each oxygen is assigned to exactly one grid
  cell, so members divided by volume is a consistent number density.

Membership needs the per-frame cavity grid, which is never saved, so it can't
be added to finished results with `hydration --cavity-results`. It writes its
own files (`PREFIX_water_membership*.json/csv/npz/png`), adds the member
waters to the viewer scenes as orange 0.7 Å spheres, and checks in every
bundle that the fixed, instantaneous and shared counts add up. An unresolved
frame gives a missing count, never zero.

## Options you'll use most

| Option | Default | What it does |
|---|---|---|
| `--hydration-cutoff` | 3.5 Å | water-oxygen contact distance |
| `--hydration-probe` | 1.4 Å | SASA probe radius |
| `--water-density-spacing` | 0.5 Å | density map bin size |
| `--water-density-level` | 0.05 Å⁻³ | display isovalue (water oxygens per Å³) |
| `--confidence` | 0.95 | level of the mean intervals; `0` for description only |
| `--start`, `--stop`, `--stride`, `--max-frames` | 0, end, 1, 500 | which trajectory frames to use |
| `--water-selection` | standard names | which residues are water, for MD input |

The hydration settings are [shared options](../reference/cli/common-options.md#hydration-settings);
everything is listed under [`crevice hydration`](../reference/cli/hydration.md).

## What you get

Tables are CSV, with units in the column names (`_A2`, `_ps`, `_degrees`):

| File | Contents |
|---|---|
| `*_hydration_residues.csv` | per residue: water counts, SASA and occupancy, with spread and intervals |
| `*_hydration_summary.csv` | every summary statistic on one row each, with its interval |
| `*_hydration_frames.csv` | per frame and residue: water counts and SASA (plus wall area and position with a cavity trajectory) |
| `*_hydration_region_frames.csv` | per frame: unique waters around the focus, and the region count and cavity volume |
| `*_hydration_shared_water.csv` | residue pairs that share a water |
| `*_hydration_observed_episodes.csv`, `*_hydration_sampled_survival.csv`, `*_hydration_mobility.csv` | trajectories only: contact episodes, persistence and residue motion |
| `*_interaction_edges.csv`, `*_interaction_edge_frames.csv` | typed contacts per pair, and per frame |
| `*_residue_conformations.csv`, `*_residue_conformation_frames.csv` | torsions and secondary structure |
| `*_analysis_evidence.csv` | everything per residue in one table |

You also get density maps (DX), figures (water-contact and SASA bars;
contact and SASA heatmaps and water time series for trajectories; density
projections) and PyMOL, VMD and ChimeraX scenes in which residues are coloured
from red (0 % water contact) through purple to blue (100 %), with a teal
density surface and dashed typed contacts. `--skip-water-density`,
`--skip-interactions` and `--skip-analysis-views` leave those parts out. An
empty pair or episode table has just its header; it doesn't mean the site is
dry.

:::{admonition} Interpreting results
:class: crevice-interpret

Water is only seen if it is in the input: a structure without explicit water
gives water contacts marked `unavailable_no_explicit_water` (empty, not zero),
while SASA is still reported. A contact is proximity, not a hydrogen bond or a
free energy, and counts depend strongly on the cutoff. Saved frames 100 ps
apart can't resolve fast exchange, so no lifetimes or rates are fitted, and
correlations with motion or volume are descriptive. Member counts that simply
follow the cavity volume at a steady density aren't evidence of wetting or
drying. Protein and solvent must be imaged consistently in your trajectory.
:::

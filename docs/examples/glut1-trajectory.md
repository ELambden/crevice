# A membrane MD trajectory: water in the GLUT1 sugar site

This example follows one region through a real membrane simulation: a human
GLUT1 glucose transporter (mutant N45Q/E329Q/N411Q) in a lipid bilayer, with
explicit water and ions. We download a short excerpt of the run, check what is
in it, measure the candidate sugar site in every frame, count the waters
inside it, and plot both over time. It takes about two minutes on two cores.

The commands and numbers below were run with CREVICE {{ version }} on the
excerpt; this page is not executed when the documentation is built.

## Download the excerpt

The full simulation is 1001 frames over 100 ns (465 MB). The example uses
every 50th frame, 21 frames 5 ns apart, with the protein, lipids, water and
ions (about 12 MB to download):

```bash
crevice fetch-example glut1-excerpt
```

The output (paths shortened):

```text
~/.cache/crevice/examples/glut1-excerpt/glut1_excerpt.gro.gz	fetched	sha256 verified
~/.cache/crevice/examples/glut1-excerpt/glut1_excerpt.gro	fetched	sha256 verified
~/.cache/crevice/examples/glut1-excerpt/glut1_excerpt.xtc	fetched	sha256 verified
~/.cache/crevice/examples/glut1-excerpt/glut1_site_candidate_region.json	fetched	sha256 verified
~/.cache/crevice/examples/glut1-excerpt/glut1_site_candidate_reference.dx	fetched	sha256 verified
~/.cache/crevice/examples/glut1-excerpt/glut1_excerpt_provenance.json	fetched	sha256 verified
~/.cache/crevice/examples/glut1-excerpt/glut1_excerpt_README.md	fetched	sha256 verified
example directory=~/.cache/crevice/examples/glut1-excerpt
```

Every file is checked against a SHA-256 hash before it is used, and a second
call reuses the verified copies. `--cache-dir` puts the files somewhere else,
and `--source-dir DIR` copies them from a local folder instead of the network.
The same thing from Python:

```python
from crevice.examples import fetch_example

data = fetch_example("glut1-excerpt")
print(data.directory, data["glut1_excerpt.xtc"])
```

`glut1_excerpt_README.md` describes the system, and
`glut1_excerpt_provenance.json` lists which source frames were kept
(0, 50, …, 1000) and their times.

## Look before you analyse

```bash
cd ~/.cache/crevice/examples/glut1-excerpt
crevice trajectory --topology glut1_excerpt.gro --selection 'protein and not name H*' \
    --inspect glut1_excerpt.xtc
```

```text
total_atom_count	123650
selected_atom_count	3814
selected_residue_count	492
frame_count	21
time_step_ps	5000.0
total_time_ps	100000.0
box_dimensions	[95.63850402832031, 96.19970703125, 131.70550537109375, 90.0, 90.0, 90.0]
has_elements	False
has_chain_ids	False
segment_ids	['SYSTEM']
```

The topology uses CHARMM names (`TIP3` water, `SOD`/`CLA` ions, DOPC, DOPE,
DOPS and cholesterol lipids) and has no element column or chain IDs, which is
typical of a GROMACS run. Residues are therefore written as `SYSTEM:GLN282`
and so on. The full simulation also contains 29 Tris molecules (`TRI`,
551 atoms); they are left out of the excerpt. The simulation protocol
(CHARMM36, 310 K, 150 mM NaCl, GROMACS 2021.5) is in
`glut1_excerpt_README.md`.

## The region

The excerpt comes with a prepared region definition,
`glut1_site_candidate_region.json`. It is the candidate sugar site around the
side-chain amides of Gln282, Gln283, Asn288, Asn317 and Asn415, residues that
line the bound sugar in the GLUT1 crystal structure
([4PYP](https://www.rcsb.org/structure/4PYP)). The definition records the
input hashes, the grid (0.2 Å spacing, 1.4 Å probe), the atomic radius set
and the obstacle selection `not resname TIP3`: protein, lipids and ions all
bound the region, while water is free to fill it. It was made with
`region-init` and `region-prepare`, as described in
[Regions](../tools/regions.md).

## Measure the region in every frame

```bash
crevice region-trajectory glut1_site_candidate_region.json \
    --out-dir glut1-run --workers 2 --water-membership
```

```text
radii=default (definition_set_used, table ac428b65f578)
cavity trajectory frames=21 statistics=glut1-run/glut1_site_candidate_cavity_statistics.json
```

Each frame is aligned to the first one on the protein Cα atoms, the empty
space inside the fixed reference region is measured, and, with
`--water-membership`, every water oxygen inside that frame's measured space is
counted as a member. The output directory holds the full cavity-trajectory
and hydration bundle; two tables answer the main questions.

`glut1_site_candidate_cavity_frames.csv` has one row per frame with the
measured volume:

```text
frame_index,time_ps,status,volume_A3,...
0,0.0,regional_measured,857.8160000000003,...
1,5000.0,regional_measured,764.2640000000002,...
```

`glut1_site_candidate_water_membership_frames.csv` adds the water counts.
`cavity_member_count` is the number of waters inside the measured space, and
`fixed_reference_count` the number inside the fixed reference region, for
comparison:

```text
frame_index,time_ps,status,volume_A3,fixed_reference_count,domain_voxel_count,cavity_member_count,...
0,0.0,measured,857.8160000000003,24,24,24,...
1,5000.0,measured,764.2640000000002,24,24,23,...
2,10000.0,measured,819.2400000000002,26,26,25,...
```

`frame_index` counts frames of the file you analysed, so here frame 10 is
source frame 500 (50 ns); `time_ps` is the simulation time either way.

## See it in a viewer

Open the scene in PyMOL from the output directory:

```bash
cd glut1-run
pymol glut1_site_candidate_analysis.pml
```

The scene shows the protein as a pale cartoon, the region cast as a grey
surface (the space measured over the run), the water density as a teal
surface and the member waters of one frame as orange spheres. Residues that
line the region are sticks coloured from red (rarely in contact with water)
to blue (in contact with water in every frame). `crevice_frame 2` and
`crevice_frame 3` switch the protein and member waters to the frames at 50 and
100 ns.
The same scene is written for ChimeraX (`.cxc`) and VMD (`.vmd`); this page
was not checked in PyMOL, and the image below is the ChimeraX scene, zoomed
onto the site, with the cartoon made partly transparent.

```{figure} images/glut1_excerpt_region_chimerax.png
:alt: ChimeraX view of the GLUT1 candidate sugar site with the region cast, water density and member waters
:width: 75%

`glut1_site_candidate_analysis.cxc` rendered by ChimeraX 1.12 (headless,
under Xvfb), followed by `transparency #1 70 target c`, `view #2,4` and
`zoom 0.55`. Grey: region cast; teal: water density; orange: the 24 member
waters of the first frame.
```

## Plot volume and water over time

```python
import pandas as pd
import matplotlib.pyplot as plt

frames = pd.read_csv("glut1-run/glut1_site_candidate_water_membership_frames.csv")
time_ns = frames["time_ps"] / 1000

fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(6, 4.5))
top.plot(time_ns, frames["volume_A3"], marker="o", color="#5b6770")
top.set_ylabel("Region volume (Å³)")
bottom.plot(time_ns, frames["cavity_member_count"], marker="o", color="#e08a1e")
bottom.set_ylabel("Member waters")
bottom.set_xlabel("Time (ns)")
fig.tight_layout()
fig.savefig("glut1_volume_waters.png", dpi=150)

print(frames[["volume_A3", "cavity_member_count"]].describe().round(1))
```

```text
       volume_A3  cavity_member_count
count       21.0                 21.0
mean       761.1                 22.0
std         47.3                  2.1
min        703.9                 18.0
25%        719.5                 20.0
50%        746.6                 22.0
75%        806.4                 24.0
max        857.8                 25.0
```

```{figure} images/glut1_excerpt_volume_waters.png
:alt: Region volume and number of member waters in each of the 21 excerpt frames
:width: 80%

The figure written by the snippet: measured region volume (top) and member
waters (bottom) in each frame.
```

The site holds roughly 18 to 25 waters throughout, and the count tends to
rise and fall with the volume (`glut1_site_candidate_water_membership_summary.csv`
reports a correlation of 0.65 over these 21 frames).

## Interpreting results

These are observations of one candidate region, chosen from literature
landmarks, in one simulation of a mutant; they are not a validated
measurement of the GLUT1 binding site, and the region boundary is an analysis
choice. With 21 frames 5 ns apart the summary intervals are marked
`insufficient_independent_blocks`, so treat the means as descriptive. The
excerpt's per-frame volumes and member-water counts are the same as those of
the full 1001-frame run at the same frames; quantities pooled over the run
(for example which residues count as lining the region) depend on which
frames were analysed. Because the excerpt omits the Tris molecules, a few
per-residue solvent-accessible areas of surface residues that Tris touched
differ slightly from the full simulation.

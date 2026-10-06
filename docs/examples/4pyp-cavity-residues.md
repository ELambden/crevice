# 4PYP: interior cavity cast and its boundary residues

[PDB 4PYP](https://www.rcsb.org/structure/4PYP) is human GLUT1 crystallised
in an inward-open conformation. Its central cavity is open to one side only, so
it is not a through-channel: `profile` refuses it and `cast` measures it
instead. [`residue-evidence`](../tools/residue-evidence.md) then lists the
residues that form the cast boundary and the residues that contact them, and
draws them as sticks in the viewer scenes. It takes about two minutes.

## `profile` does not resolve a through-channel

```bash
crevice fetch 4PYP --cache-dir .crevice/pdb
crevice profile 4PYP -o 4PYP_profile.json
```

```text
crevice: error: No enclosure probe resolved a stable channel (automatic choice, 21 probes tried, 0.50-3.00 A). Failures: No connected channel with two open ends was resolved (0.5, 0.55, ..., 3 A). An explicit enclosure_radius (--enclosure-radius) runs a single probe
```

(exit status 2; probe list shortened here.) No probe of the ladder found a
laterally enclosed run with two open ends, so no profile is reported. This is
consistent with an inward-open structure but does not show that the protein is
closed; see [Outputs and unresolved results](../getting-started/outputs.md).

## Cast the dominant cavity

```bash
crevice cast 4PYP --out-dir 4PYP_cast --skip-hydration
```

```text
cast regions=1 volume_A3=3245.250 spacing_A=0.5 out_dir=4PYP_cast
```

`cast` rolls a 0.8 Å probe (`--probe-radius`) through the protein interior,
inside the envelope of a 6 Å outer probe, and keeps the dominant buried region
(`--max-cavities 1`). On a 0.5 Å grid that region is 3245.25 Å³
(25 962 samples; `4PYP_void_cast.csv`). The volume depends on the probe, the
outer envelope and the grid; compare settings before interpreting it.
`4PYP_volume.dx` is the measured binary map, and `4PYP_volume.pml`, `.vmd`
and `.cxc` show the cast inside the protein. `--skip-hydration` leaves out the
water analysis: this entry contains no water molecules, so water contacts would
be reported as unavailable. The default `--selection protein` uses
non-HETATM heavy atoms, so the bound detergent (BNG) is not an obstacle;
`--selection all` keeps it.

## Boundary residues and their partners

```bash
crevice residue-evidence 4PYP_cast/4PYP_viewer.pdb --volume-dx 4PYP_cast/4PYP_volume.dx \
    --scene 4PYP_cast/4PYP_volume.pml --out-dir 4PYP_evidence --prefix 4PYP --skip-hydration --dpi 150
```

```text
residue-evidence lining=62 out_dir=4PYP_evidence
```

The structure argument is the cast bundle's viewer PDB, the file the scene
shows: its `.identities.json` sidecar restores the original chain and residue
numbers, and the command checks that the scene's atoms match the analysis
(giving the original `4PYP` here makes `--scene` fail with "Scene protein
coordinates and atom identities do not match the analysis reference", because
the cast scene contains only the selected protein atoms).

Each triangle of the measured cast surface is assigned to the residue whose
atom surface is nearest (within 1.5 Å, `--lining-distance`). 62 residues
receive boundary area (`boundary_lining`) and 100 others make a nonlocal
contact (heavy-atom centres within 4.5 Å, `--contact-cutoff`; different chain
or more than two residues apart) with a lining residue
(`nonlocal_lining_partner`). 92% of the 1818 Å² surface was assigned.

```text
residue,group,role,boundary_area_A2,boundary_fraction,minimum_boundary_gap_A,nonlocal_lining_contact_count,contacted_lining_boundary_fraction
A:TRP388,unassigned,boundary_lining,94.08,0.0517,-0.165,7,0.114
A:GLN161,unassigned,boundary_lining,68.89,0.0379,-0.219,5,0.064
A:ILE404,unassigned,boundary_lining,61.22,0.0337,-0.181,4,0.117
```

(first rows of `4PYP_residue_evidence.csv`, rounded.) A negative
`minimum_boundary_gap_A` means the voxel surface passes slightly inside the
atom sphere, which is expected for a surface built from 0.5 Å voxels.
`4PYP_residue_evidence_partners.csv` lists every partner contact with its
lining residue, interaction candidate, distances and closest atoms.

```{figure} images/4pyp_residue_evidence.png
:alt: Bar charts of boundary contributors and nonlocal partners of the 4PYP cavity
:width: 95%

`4PYP_residue_evidence.png`. Left (orange-gold): the 12 residues with the
largest share of assigned boundary area (%). Right (magenta): the 12 nonlocal
partners ranked by the boundary fraction of the lining residues they contact
(%). Residue IDs are category labels and are always drawn.
```

```{figure} images/4pyp_residue_context_chimerax.png
:alt: ChimeraX render of the 4PYP cavity cast with boundary and partner residues as sticks
:width: 70%

`4PYP_residue_context.cxc` rendered by ChimeraX 1.12 (headless,
`scripts/check_chimerax.py --xvfb`). Violet, 55% opaque: the cast surface
(violet marks a cast that is not a resolved channel; channel casts are teal).
Grey-blue cartoon: the protein. Gold sticks: the top 12 boundary residues;
magenta sticks: the top 12 nonlocal partners (`--stick-residues top:12`, the
default). No labels are drawn (`--annotate` adds them). The check confirmed
the stick residues, their colours, the protein coordinates (error 0.000 Å),
the surface (within 6 × 10⁻⁵ Å of the reference mesh) and the camera. The
PyMOL and VMD versions were not opened for this page.
```

`--stick-residues all` draws every boundary residue and partner, and
`--stick-residues A:TRP388,A:GLN161` or `@ids.txt` draws a chosen list; the
tables and the bar chart do not change.

## What was checked

The commands were run with CREVICE 0.1.0 on the downloaded entry; values and
file names are copied from that run. Boundary area and contacts are geometric
evidence: they identify candidates for follow-up, not residues shown to
control transport.

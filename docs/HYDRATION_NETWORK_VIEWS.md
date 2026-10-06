# Coordinated hydration and residue evidence

CREVICE's standard hydration workflow produces aligned water-oxygen density,
chemically typed residue contacts, conformation histories, a combined evidence
table, and matching PyMOL, VMD and ChimeraX scenes. These are geometric and
observational measurements. They do not establish that a residue controls
transport, that a cavity is substrate-accessible, or that a contact is energetically
favourable. The underlying water-contact and SASA definitions are in
[Residue hydration and solvent exposure](HYDRATION_METHODS.md).

## Usage

```bash
crevice hydration system.gro --trajectory run.xtc \
  --cavity-results path/to/system_cavity_statistics.json \
  --max-frames 1100 --out-dir results --prefix system

crevice hydration structure.cif --volume-dx path/to/structure_volume.dx \
  --out-dir results --prefix structure
```

Install the optional `crevice[hydration]` extra (MDTraj) for the complete
workflow; the broader `science` extra also includes it. MDAnalysis and
scikit-image are required CREVICE dependencies. Without MDTraj, DSSP is explicitly
unavailable. If scikit-image cannot be imported, PyMOL's triangle surfaces are
unavailable; the density DX remains usable in VMD/ChimeraX.
Unassigned or ambiguous residue chemistry is reported and excluded from typed
observations, while other analyses continue.

These options are available wherever the standard hydration analysis runs:

| Option | Default | Meaning |
|---|---:|---|
| `--water-density-spacing` | 0.5 Å | Spatial bin size; no silent coarsening |
| `--water-density-smoothing` | 0.5 Å | Gaussian sigma for the display field only |
| `--water-density-level` | 0.05 oxygen/Å³ | Absolute isosurface level |
| `--skip-water-density` | off | Omit spatial density |
| `--skip-interactions` | off | Omit typed chemistry/conformation analysis |
| `--skip-analysis-views` | off | Omit coordinated scene scripts, evidence table and overview figures |

`--skip-hydration` is the standard-workflow opt-out. A supplied cavity context
keeps its original alignment reference and verifies matching source identities,
frame indices and times, and input hashes when they are recorded.
Unaligned trajectories do not produce a mean spatial density map.

## Colour and presentation

Residue colour is the absolute **contact occupancy**: the fraction of analysed
frames with at least one water oxygen within the configured contact cutoff
(default 3.5 Å) of a residue heavy atom. The linear RGB scale is `(1−p, 0, p)`:
0% red, 50% purple, 100% blue. It is not divided by the wettest residue or based
on its mean number of waters. Missing water observations are grey. A static
structure provides only 0/100% observed contacts; neither absence of deposited
waters nor a single dry observation establishes equilibrium dryness.

The same colours, residue identities, reference coordinates and camera rotation
are used in the three viewers. Direct geometric boundary contributors appear
as sticks over a neutral cartoon protein, with a subtle reference cavity and a
teal water-density surface. By default the scenes contain no text; with
`--annotate` they also show a red-to-blue colour bar with 0/50/100% labels, its
caption and, for named regions, the region status banner. The scene JSON
records this text either way. Surface smoothing changes display fields only.
Original measured cavity occupancy and unsmoothed water counts remain available.
No `crevice_grid` or `crevice_region_grids` helper selections are introduced.

Open `PREFIX_analysis.pml`, `PREFIX_analysis.vmd` (or `.tcl`), or
`PREFIX_analysis.cxc` with all companion files together. Open the
`PREFIX_analysis_network` counterpart to start with typed contact lines visible.

| View/control | PyMOL / VMD | ChimeraX |
|---|---|---|
| Boundary hydration | `crevice_view hydration` | `creviceView hydration` |
| Typed contacts and their endpoints | `crevice_view network` | `creviceView network` |
| All focused partner residues | `crevice_view partners` | `creviceView partners` |
| Cavity geometry | `crevice_view geometry` | `creviceView geometry` |
| Representative snapshot | `crevice_frame 1` | `creviceFrame 1` |
| Hide / show member waters (own group, on by default) | `disable`/`enable crevice_member_waters`; VMD `mol off`/`mol on $crevice_member_waters` | `hide`/`show #<id> models` (model `crevice_member_waters`) |

With `--water-membership`, orange spheres mark the water oxygens inside the
displayed snapshot's measured region ([definition](HYDRATION_METHODS.md#instantaneous-measured-region-water-membership)).
Unlike the all-frame colours and density, this group belongs to the displayed
snapshot and follows `crevice_frame`, keeping a user's hidden state.

Trajectories provide first, middle and last analysed protein snapshots, aligned
to the same reference. Colours and density remain all-frame summaries. Reference
contact lines are hidden outside the first snapshot, because their actual atom
positions belong to that reference. Complete per-frame typed observations remain
in the numerical export. These controls are not a full trajectory movie player.

`crevice_render` exports a 3200×2400 PyMOL PNG or 2400×1800 VMD TGA. ChimeraX has
`PREFIX_analysis_render.cxc` for a supersampled PNG. These high-resolution exports
can be expensive for transparent molecular scenes. Native renderers have different
lighting/transparency implementations; coordinate and colour agreement does not
imply pixel-identical images.

## Geometric role colours versus hydration colours

The separate `residue-evidence` contexts use gold for direct boundary residues
and magenta for nonlocal partners, both as opaque sticks. Those colours describe
geometric roles, not water occupancy. The standard context shows the top 12 per
role by default, and `residue-evidence --stick-residues all` shows every
boundary residue and partner. Neither changes this hydration view's absolute
occupancy scale. Missing water observations remain grey
in hydration scenes, even when geometric role highlights are available.

## Spatial hydration definition

Every analysed frame contributes to the denominator. Water oxygen coordinates
are periodically imaged into a fixed reference neighbourhood and transformed by
exactly the protein's proper rigid alignment. Orthogonal and triclinic lattice
images are enumerated. The domain is the union of distance neighbourhoods around
occupied reference cavity grid points (2 Å margin), or around the selected
reference residue heavy atoms when no cavity map is supplied (contact-cutoff
margin). It is a **fixed neighbourhood**, not membership in a recomputed
instantaneous cast. Nearby waters can lie outside the geometric cavity.

For voxel volume ΔV and N observed frames, `density = Σ oxygen counts / (N ΔV)`.
The DX origin is the first voxel centre, axes are in input/reference Å, and the
last grid index varies fastest. The map integral is the mean number of oxygen
observations in the domain. Voxel occupancy instead counts frames with one or
more oxygens in the voxel, divided by N; its value depends on bin size. A domain
spanning periodic cells can contain multiple images of a water identity, which
is explicitly reported. The default grid guard is eight million samples.

`_water_density.dx` is unsmoothed number density. `_water_occupancy.dx` is voxel
occupancy. `_water_density_display.dx` applies the declared Gaussian for viewing.
The NPZ retains integer counts, occupied-frame counts, first-half counts and
per-frame domain counts. The JSON records units, normalization, domain, grid,
frame coverage and integrals. No bulk-relative free energy is inferred. When no
explicit water population exists, the status is unavailable and no artificial
zero-density DX is written.

The aligned histogram and density-unit conventions follow the
[MDAnalysis density documentation](https://docs.mdanalysis.org/stable/documentation_pages/analysis/density.html).
Compare the map integral with the regional water count as a consistency check.
Fine-bin spatial peaks and half-window agreement are resolution-sensitive, and a
density map does not establish a converged equilibrium hydration landscape.

## Typed residue interactions

Standard amino-acid atom templates and named protonation variants are used.
Supported variants include HID/HIE/HIP, HSD/HSE/HSP, ASH/GLH, LYN and CYM/CYX.
Generic histidine without an explicit imidazole H has unresolved tautomerism;
its imidazole acceptors are not guessed. Protonated acids without their H do not
guess the OH oxygen. CCD DLE/DVA share LEU/VAL chemical templates with D identities and signed
torsions preserved; they are excluded from unvalidated DSSP classification.
See [the component review](BIOLOGICAL_REGIONS.md). Unrecognised modified residues and ambiguous duplicate atom
names are excluded and listed. There is no general ligand/cofactor perception,
terminal-charge model, hydrogen construction, pKa assignment or energy model.

Hydrogens use supplied bonds when present. Otherwise, attachment requires the
nearest heavy atom within the same residue, a covalent-distance bound and a
0.15 Å nearest-neighbour ambiguity separation. This fallback is reported. A
partial bond topology is not silently completed. Missing hydrogens or unknown
chemistry are not negative interaction observations.

| Evidence | Required geometry |
|---|---|
| Directional H bond | Typed D/A; explicit D–H; D–A ≤3.5 Å, H–A ≤2.5 Å, D–H–A ≥150° |
| Salt-bridge candidate | Oppositely charged named sidechain groups within 4 Å |
| Aromatic candidate | Complete planar template rings; plane RMS ≤0.15 Å; centroid distance 2.5–5.5 Å |
| Parallel aromatic candidate | Plane angle ≤30° and both centroid offsets ≤2 Å |
| Edge/face aromatic candidate | Plane angle ≥60° and at least one centroid offset ≤2 Å |
| Directional water bridge | The same water oxygen forms qualifying protein–water H bonds to two residues |

Distances and angles use periodic minimum-image vectors. Protein–protein H bonds
retain donor/acceptor direction. Aromatic and salt candidates do not assign
interaction energies. A water bridge is distinct from the shared-water
**proximity** table, which requires only two residue contacts. Water–water
networks are not included. Only edges with both endpoints in the focused residue
set are analysed; an empty focused edge set is not a whole-protein negative.

The geometry is consistent with conventional directional analyses described in
[MDAnalysis H-bond documentation](https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html).
Aromatic thresholds are informed by the primary
[PLIP configuration](https://github.com/pharmai/plip/blob/master/plip/basic/config.py);
CREVICE is a conservative implementation, not a reproduction of PLIP chemistry.
CREVICE's chemical typing has not been independently validated against another
tool.

Every edge has detected-frame occupancy, a declared all-frame denominator,
example atom identities/geometry and approximate autocorrelation-aware mean
bounds where estimable. Constant, sparse or correlated data can have no supported
interval. No family-wide inference, parameter uncertainty or replica confidence
is implied. Native network views show up to 12 persistent edges present in the
reference, ordered by occupancy; the complete edge CSV/JSON/NPZ is retained.

## Conformation and combined evidence

Backbone φ/ψ and available χ1–χ4 torsions use named atom quadruplets. Missing,
degenerate or broken geometry gives NaN. Circular means and mean resultant
lengths avoid ordinary arithmetic averaging across ±180°. Coarse χ1 sectors
(0–120°, 120–240°, 240–360°) are labelled gauche+, trans and gauche−. These are
not a residue-specific rotamer library. Wet/dry endpoint χ1-step comparisons
are descriptive and cannot resolve events between saved frames.

[MDTraj simplified DSSP](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.compute_dssp.html)
provides H/E/C backbone context for complete supported residues. Its virtual-H
backbone geometry is separate from the explicit-H interaction network. Missing
backbones or the optional dependency give explicit unavailable status. Conformation
fractions are descriptive; they have no sampling intervals.

`_analysis_evidence.csv` joins role, boundary footprint, hydration, geometric SASA,
counts/occupancy sums by contact type, torsions and secondary structure. These
remain separate evidence axes; there is no universal importance score. Contact
occupancy sums can exceed one and are not probabilities. Anatomical curation,
parameter and ranking stability, independent replicas and functional comparators
are needed before functional interpretation.

No HTML report is written. The per-frame histories are CSV tables:
`_hydration_frames.csv` (contact waters, SASA, boundary area and centroid
distance per residue and frame), `_hydration_region_frames.csv` (regional volume
and fixed-neighbourhood waters per frame), `_residue_conformation_frames.csv`
(φ, ψ, χ1–χ4 and secondary structure per residue and frame) and
`_interaction_edge_frames.csv` (typed edges observed in each frame). When a cached
cavity context is available, the original cavity radius/diameter plots and their
statistical bounds are copied as `_analysis_cavity_*` figures. Those profile
images retain the context's original analysis window, which can differ from a
requested subset.

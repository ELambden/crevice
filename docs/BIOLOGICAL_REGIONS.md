# Biological identity and region review

CREVICE's benchmark registry (`crevice fetch --list`, `crevice benchmark`) names
a small set of public PDB entries used to exercise the workflows. Their deposited
records support specific identities and assemblies, summarised below. These
records do not validate the regions CREVICE selects, and they do not establish
substrate access. Every benchmark entry therefore carries
`curation_status: uncurated` (or `set_aside`), and `curated_benchmark_ids()`
returns an empty tuple.

| Input | Deposited identity / assembly | Geometric interpretation still requiring review |
|---|---|---|
| 4PYP | Human GLUT1, inward-open; assembly 1 is a monomer | Extent and accessibility of the alternating-access cavity |
| 1AF6 | *E. coli* maltoporin–sucrose complex; assembly 1 is a trimer | Each protomer channel and the separate central interfacial space |
| 1GRM | Gramicidin A; deposited assembly 1 is a dimer | Selected channel geometry and experimental membrane context |
| 2CHB | Cholera toxin B-pentamer complex | Set aside: not a channel benchmark; no pore geometry is claimed |

Sources: [4PYP](https://www.rcsb.org/structure/4PYP),
[1AF6](https://www.rcsb.org/structure/1AF6),
[1GRM](https://www.rcsb.org/structure/1GRM),
[2CHB](https://www.rcsb.org/structure/2CHB).

## Site landmarks versus geometric candidates

A cavity found by CREVICE is a geometric candidate. Where the literature
identifies a binding site, those residues are independent landmarks, and the two
should be compared explicitly rather than assumed to coincide.

For example, the deposited 4PYP ligand is **BNG, nonyl beta-D-glucopyranoside**:
a detergent ligand with a sugar headgroup, not a deposited glucose molecule. Deng
and colleagues describe sugar coordination by Gln282, Gln283, Asn288, Asn317 and
Asn415 (Extended Data Figure 2 of the
[original structure paper](https://www.nature.com/articles/nature13306)). These
are source-backed site landmarks, separate from CREVICE's geometric residue
candidates.

When the analysed structure is a simulation model rather than the deposited
entry, compare its sequence with the reference (for example
[UniProt P11166](https://www.uniprot.org/uniprotkb/P11166/entry) for GLUT1) and
record mutations, missing residues and the fitted scaffold before transferring
landmarks or ligand coordinates. Large fitted RMSDs or scaffold-dependent
transfers mean that a region cannot be assigned to a transport state or to the
substrate-binding cavity from geometry alone.

Good practice for a region review:

- record deposited ligand coordinates and literature residues as positive site
  landmarks, and occupied protein or ligand atom cores as physical negative
  controls;
- do not label an unreviewed inter-helical space as ligand-inaccessible solely
  because of a cast, a fitted axis or proximity to membrane lipids;
- add a different reference region as a distinct analysis compared with the
  original, never as a silent replacement.

The named-region workflow that supports this is described in
[Named biological-region workflow](REGION_WORKFLOW.md).

## Modified amino-acid coverage

CCD [DLE](https://www.rcsb.org/ligand/DLE) and
[DVA](https://www.rcsb.org/ligand/DVA) identify D-leucine and D-valine. Their
heavy-atom names, elements and bond orders match the LEU/VAL templates in the CCD
records. CREVICE shares those chemical templates while retaining D residue
identities, coordinates and signed torsions. It does not mirror them into L
residues or assign a residue-specific rotamer library.

In 1GRM, twelve D-amino-acid residues use these conservative chemical templates.
Four terminal components (two FVA and two ETA) remain unsupported; they are not
flattened into ordinary amino-acid templates. DLE/DVA are excluded from DSSP
classification. No hydrogens or protonation states are constructed, so the
deposited 1GRM entry cannot supply explicit-H hydrogen-bond observations or
water hydration observations.

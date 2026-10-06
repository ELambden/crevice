# Named biological-region workflow

CREVICE can prepare and analyse separate named regional observables while
retaining their input identities, coordinate frame, landmark evidence and review
status. It never promotes a candidate to a biological assignment automatically.

## Commands

```bash
crevice region-init system.gro --trajectory run.xtc \
  --region-id site_candidate --landmarks-json landmarks.json \
  --obstacle-selection 'not resname TIP3' --spacing .2 \
  --output site-request.json
crevice region-prepare site-request.json --out-dir site-review
crevice region-trajectory site-review/site_candidate_region.json \
  --out-dir site-trajectory --max-frames 1100 --workers 2
```

`landmarks.json` is a list of exact source residue identities, explicit heavy-atom
names and evidence. For example, an entry has this form:

```json
{
  "residue": "SYSTEM:GLN282",
  "atoms": ["OE1", "NE2"],
  "evidence": "Describe the source-supported landmark and its mapping to this input."
}
```

Use the identities actually present in the input; chain names, insertion codes
and modified residue names matter. Missing or ambiguous atoms are errors. A
landmark's terminal-atom centroid is an analysis choice, not a ligand pose.

For an existing reference, supply `--reference-volume-dx` to `region-init`.
Landmarks are optional in this mode. Static PDB, mmCIF and GRO snapshots can be
prepared without `--trajectory`. Static mmCIF uses the canonical entity-aware
parser and original chain/residue identities; MDAnalysis evaluates selections on
those atoms. Static mmCIF crystal cells are not treated as MD periodic boxes, and
assembly operators are not applied automatically. The current mmCIF adapter and
trajectory command require reference frame/model 0; unsupported indices fail.
GRO and trajectory physical obstacles retain the declared periodic treatment.

Preparation and analysis require empty output directories. Alternative reference
regions use distinct IDs and directories. Existing evidence is never replaced by
these commands. All normal geometry, hydration, density, interaction and viewer
outputs are produced by `region-trajectory`; standard hydration opt-outs and
sampling/uncertainty options remain available. The prepared definition supplies
geometry, alignment and obstacle settings, with no conflicting CLI overrides.

## What the definition records

Schema version 2 records the region ID and label; candidate/reviewed status,
rationale and sources; input topology/trajectory paths and SHA256 digests;
reference frame; protein and obstacle selections; assembly, preparation and
membrane-context descriptions; alignment residue IDs and rationale; spacing,
probe, margin, axis, grid phase and allocation limit; and residue/coordinate
landmarks; and `radii`, the atomic radius set used at preparation (name,
source, table hash; see [Atomic radii](methods/atomic-radii.md)). Paths
resolve relative to the definition, independent of the working directory.
Large input files are hashed incrementally. `region-init` writes a version 1
request without `radii`; `region-prepare` writes version 2. Later steps use the
stored set (`region-trajectory`, `hydration --cavity-results`) or check it
(`region-compare`), and stop on a different set unless
`--allow-radii-mismatch` is given. A version 1 definition is read as the
default radius set; its reference map is not recomputed.

`region-prepare` writes a companion binary reference DX, an identity/coordinate
fingerprint for the selected reference frame, a resolved definition, diagnostics,
and PML/VMD/TCL/CXC review scenes. These viewer bundles are portable when their
companion assets are kept together. The simulation inputs are not copied into
viewer bundles; rerunning an analysis still requires those hashed inputs.

A reviewed status requires a reviewer, review date, scope and supporting sources.
It remains declared provenance, not software certification or evidence that the
trajectory samples a transport transition. Edit a new version of a prepared
reference when recording a review, and keep the preceding candidate version.

## Bounded landmark proposal

The optional `landmark_sphere` proposal gives each residue equal weight, computes
the mean of its declared local atom centroids, and samples a sphere around the
result. Defaults are a 6 Å sphere, a maximum 3 Å distance to an eligible seed,
and a 0.5 Å ambiguity-distance threshold. All are explicit in the request.

Eligible probe centres must clear the selected atom spheres. Neighbouring centres
are connected only when the intervening segment clears those spheres. CREVICE
selects the component nearest the landmark anchor. If no admissible centre is
close enough, or competing disconnected components are similarly close, it
fails rather than forcing a result. The selected component is swept by the probe,
clipped to the sphere, and restricted to atom-clear sample centres. No grid is
silently coarsened.

The spherical crop is not a cavity mouth or a biological boundary. Swept boundary
voxels approximate a continuous surface; atom-clear centres do not prove every
point of a voxel is clear. A missing sampled connection does not prove continuum
disconnection. The resulting reference remains fixed during trajectory analysis;
it is not re-centred on changing sidechains each frame. Subsequent measurements
use the existing reference-neighbourhood method and the declared growth margin.

## Review overlays and diagnostics

Gold sticks show declared residue landmarks, purple spheres show optional
transferred coordinate landmarks, teal shows the reference region, and green
sticks show nearby nonprotein physical-obstacle atoms. The neutral protein is a
cartoon. These colours describe review evidence, not hydration. Trajectory
hydration figures keep the absolute red-to-blue water-contact occupancy scale.
The named region's status is recorded in every output and shown in the scenes
with `--annotate`.

Diagnostics report landmark-to-reference distances, transferred-point coverage,
fit RMSDs when supplied, selected/competing probe components, artificial crop
contact, and reference samples that overlap selected obstacles. Existing maps may
contain such samples when they were calculated with different obstacles; they
remain anchor domains, and every trajectory measurement recalculates free space.
A supplied transfer RMSD above 2 Å triggers a review flag, not an automatic
biological rejection. No membrane side or substrate access is inferred from an
axis, attractive cast, or a spherical probe.

## Comparing regions

```bash
crevice region-compare original/original_pocket_cavity_statistics.json \
  site-trajectory/site_candidate_cavity_statistics.json --out-dir comparison
```

The comparison requires identical source hashes, saved frames/times, protein
identities, alignment and relevant geometric/obstacle settings. Each region keeps
its own reference and profile coordinate origin. A bundle without region
provenance (for example one written by `cavity-trajectory`) requires an explicit
prepared `--left-definition` or `--right-definition`, together with the input-hash
record of the run that produced it, supplied as `--left-input-provenance` or
`--right-input-provenance`. Source, reference and numerical settings are checked
before accepting the observations, because current file contents alone cannot
establish which inputs produced an older result. This binding does not rerun the
analysis.

Outputs include paired all-frame volume/water observations, volume-difference
statistics, separate shaded width profiles, boundary-residue overlap and tables
of boundary footprint and hydration. Empirical frame ranges and estimable
pointwise mean confidence bounds are labelled separately. Inadequate sampling
still withholds intervals. Missing residue hydration is not converted to zero.

Different crops measure different spaces. Their volume difference is not a
transport-state change, and residue boundary ranks are not functional importance
scores. Typed-edge counts depend on the region-specific residue focus; an edge
absent from a focus is not an observed zero contact occupancy. Hydration maps and
the standard regional water counts refer to fixed neighbourhoods. `region-trajectory
--water-membership` adds a separate, paired count of oxygens inside each frame's
measured region ([definition](HYDRATION_METHODS.md#instantaneous-measured-region-water-membership));
`region-compare` does not yet compare that observable.

## Example: literature landmarks for GLUT1

For human GLUT1 (4PYP), a landmark request can use the residues that coordinate
the sugar headgroup of the deposited BNG ligand, Gln282, Gln283, Asn288, Asn317
and Asn415, described in
[Deng et al., Extended Data Figure 2](https://www.nature.com/articles/nature13306),
mapped to the identities present in your input. BNG is a detergent ligand, not
glucose. Transferred deposited coordinates are best shown only as independent
review overlays, and an existing pocket definition should stay a separate region.
Because the landmarks define such a candidate, their contribution to its boundary
is not independent site validation. Starting-model provenance, anatomical
extent, membrane sides and substrate accessibility still need review; running
the workflow does not provide biological validation.

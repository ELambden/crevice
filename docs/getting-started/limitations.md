(limitations)=
# Limitations

CREVICE measures geometry, and every number it gives depends on the settings
behind it. This page collects the caveats that apply across the tools, so you
can keep them in mind when you read and report results. Each tool page has
its own "Interpreting results" note with the details for that tool.

## Geometry is not biology

- **A measurement is not a validated biological result.** A radius, volume,
  route or residue list is a well-defined measurement under the settings you
  chose. The software is tested on synthetic structures with known answers;
  the real proteins in these pages are worked examples, not curated
  references, and none of the bundled benchmark entries is curated. Whether a
  space matters for function needs independent evidence, such as mutagenesis,
  functional data or other structures (see
  {ref}`How far to trust the numbers <validation-levels>`).
- **Contacts and correlations are associations.** Residue networks, motion
  correlations and hydration statistics point to candidates for follow-up,
  not to causal pathways or residues that control transport.

## Choices you make

- **The channel axis is chosen, not discovered.** With `--axis auto`,
  [`profile`](../tools/profile.md) tries the three principal directions of the
  selected atoms and keeps one only if it holds a connected, enclosed passage
  with two open ends; you can also fix the direction and a seed. The method
  assumes a channel that runs roughly in one direction, so strongly curved or
  branched passages need [`crevice paths`](../tools/paths.md). Trajectories,
  fingerprints and cavity trajectories then keep one fixed axis for every
  frame. In the bundled benchmark structures, too, the axis is chosen by the
  analysis or set by hand, not taken from a curated reference.
- **Atomic radii.** Radii come from a radius set (by default Bondi radii, with
  CHARMM36 sizes for recognised ions), and a different set changes the
  results by design. No set is correct in general: HOLE's radii, for example,
  are larger to stand in for missing hydrogens. An atom whose element has no
  radius in the set (most transition metals and lanthanides, for example) is
  left out of a structure with a warning, and an MD selection containing one
  is refused; give such atoms radii with `--radii FILE`. Report the radius set
  with your results. See [Atomic radii](../methods/atomic-radii.md).
- **Grids, probes and margins.** Cast and cavity volumes change with the
  probe radius, grid spacing and grid phase. In a
  [cavity trajectory](../tools/cavity-trajectory.md) with
  `--geometry-mode reference-region`, the region is the free space within
  `--region-margin` (2 Å by default) of the reference, so the volume also
  depends on that margin. A coarse grid can underestimate a small volume
  substantially. Check a few settings on your system before you interpret a
  number, and report them
  ([how](../CAVITY_OBSTACLES.md#checking-numerical-sensitivity)).
- **Automatic probe choices are rules.** The enclosure probe of `profile`
  (`--enclosure-radius auto`) and the route probe of `paths`
  (`--probe-radius auto`) are chosen by a stated rule over a fixed ladder of
  probes from 0.5 to 3.0 Å. For `profile`, the same channel must appear at two
  or more probes, channels are grouped by a parallel axis (|cos| ≥ 0.99), the
  same end types, a minimum radius within 0.05 Å or 2 % and at least 60 %
  axial overlap, and the smallest probe of the channel found at the most
  probes is used. For `paths`, route sets with widths within 0.05 Å or 2 % are
  grouped, and the smallest probe of the set found at the most probes is
  used, so narrow side routes seen only at small probes are left out (they
  stay listed in `*_path_probes.csv`). The thresholds are practical choices,
  not physical constants; every probe tried is recorded, and a number fixes
  the probe. See [Choosing the enclosure probe](../tools/profile.md#choosing-the-enclosure-probe)
  and [Choosing the probe automatically](../tools/paths.md#choosing-the-probe-automatically).
  The cast probe is never chosen automatically, because there the probe *is*
  the measurement. In the [benchmarks](../methods/benchmarks.md), changing
  the grouping tolerances changed nothing and no variant moved a resolved
  minimum radius by more than 0.01 Å, but *whether* a channel resolves
  depends on the support rule: requiring the channel at three or four probes
  instead of two left 4 or 7 of 9 channels unresolved. If a channel you
  expect is reported as unresolved, try a fixed `--enclosure-radius`.
  When two channels are found at the same number of probes, the one at the
  smaller probe is used. This can pick a short capped part of a pore over the
  whole channel: for aquaporin-1 (1J4N) with `simple.rad` the reported
  profile misses the narrowest point, which a fixed `--enclosure-radius 0.75`
  recovers (see the [benchmarks](../methods/benchmarks.md)). Check the
  resolved axial range, and fix the probe when it is shorter than expected.

## What the input can and cannot show

- **Hydrogens in static structures.** Most crystal and cryo-EM structures have
  no hydrogen atoms. Directional hydrogen bonds are only counted between
  explicit hydrogens, so in such a structure the typed residue interactions
  list salt-bridge and aromatic candidates but no hydrogen bonds, and their
  absence is not evidence that none form. Protonation states are taken from
  residue names, never guessed. For hydrogen bonds, add hydrogens with a
  protonation tool you trust, or use a simulation that has them.
- **Water needs explicit waters.** The water analyses see only the waters in
  the input. A structure without them gives water contacts marked
  `unavailable_no_explicit_water` (empty, not zero), and missing water in a
  static structure is not evidence of a dry pocket. Saved frames far apart
  cannot resolve fast exchange. See [Hydration](../tools/hydration.md).
- **Polymer residues recorded as HETATM are found by a distance rule.** In a
  PDB file, a modified residue counts as protein when MODRES lists it or when
  it is peptide-bonded (C-N within 1.75 Å) into the chain
  ([details](inputs-and-selections.md#polymer-residues-recorded-as-hetatm)).
  Links other than a backbone C-N bond (nucleic-acid backbones, side-chain
  cross-links) are not followed, so such a residue stays a heteroatom unless
  MODRES lists it. In trajectory commands only the selection `protein` is
  extended this way; other MDAnalysis selections (`--obstacle-selection`,
  water and permeation selections) are used as written, so name non-standard
  residues explicitly (`protein or resname DLE`).
- **Unresolved profiles are reported, not forced.** If no unique two-ended
  channel is found, the profile is reported as unresolved with the reasons,
  and no profile is invented; `report` then measures a fallback cast. An
  unresolved profile does not show that the protein is closed. See
  [Outputs and unresolved results](outputs.md).

## Estimates to read with care

- **The conductance is an order-of-magnitude estimate.** The
  [fingerprint](../tools/fingerprint.md)'s conductance is HOLE's geometric
  estimate: the lumen treated as a tube of bulk electrolyte. It ignores
  ion–wall interactions, selectivity, dehydration and surface charge, and it
  is not a measured conductance or a prediction of function. On gramicidin A
  the uncorrected value is 10–40 times the measured single-channel
  conductance, and HOLE's empirically corrected value 2–7 times
  ([details](../methods/fingerprint.md#conductance-estimate)). Across four
  open channels with published single-channel conductances (gramicidin A,
  OmpF, GLIC, open MscS) the estimate was 1.5–17.6 times the measured value,
  and it cannot be computed for K⁺ channels whose filter is narrower than the
  smallest probe ([benchmarks](../methods/benchmarks.md)). Use it to compare
  pores with each other, not to predict a conductance.
- **Trajectory uncertainty is approximate.** Intervals for a mean assume
  stationary, representative sampling and leave out force-field, grid and
  probe errors; one simulation is not a set of replicas. Frame ranges
  describe fluctuations, not confidence in the mean.
- **Trajectory memory grows with the frame count.** A cavity or region
  trajectory keeps every frame's atoms and per-frame record in memory until
  the run ends. For a 1000-frame GLUT1 membrane simulation the run peaked at
  about 7.5 GB, roughly 7 MB per frame whatever the grid spacing
  ([benchmarks](../methods/benchmarks.md)). For long trajectories, analyse
  frame ranges with `--start`, `--stop` and `--stride`.
- **Viewer scenes should be checked.** The PyMOL, VMD and ChimeraX scenes are
  generated scripts; open them and check them before using them in a figure.

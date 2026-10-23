# Benchmarks against other tools

The repository's `benchmarks/` directory holds a reproducible comparison of
CREVICE with established programs on curated inputs. Its README lists every
input, setting, comparator version and licence; the result tables are in
`benchmarks/results/`. This page summarises the method and the main results.

All comparisons are software comparisons on static structures (and one short
trajectory). Agreement with another program checks the implementation; it
is not a biological validation. CREVICE runs with its defaults unless an
option is named, and no setting was tuned to match a comparator.

## Method

* **Matched inputs.** Every tool in a comparison reads the same atoms
  (polymer heavy atoms; waters, ions, lipids, detergents and ligands
  removed), the same atomic radii (HOLE's `simple.rad` or Bondi radii, given
  to CREVICE as `--radii`) and the same axis and start point. Axes come from
  a least-squares symmetry operator of the assembly and agree with the RCSB
  symmetry axis to within about 1°.
* **Channels** (gramicidin A 1GRM, KcsA 1K4C, MthK 1LNQ, aquaporin-1 1J4N,
  OmpF 2OMF, GLIC 4HFI, nicotinic receptor 2BG9 and 1OED, MscS 2VV5 and 2OAU,
  TWIK-1 3UKM) are profiled by `crevice profile` and HOLE 2.3.1. The compared
  quantity is the radius of the largest sphere centred in each plane normal to
  the axis. CREVICE first runs with defaults; if no channel resolves,
  `--lateral-exits` and then a 12 Å search window are tried, and each attempt
  is recorded.
* **Conductance.** The geometric estimate of the
  [pore fingerprint](fingerprint.md) is compared with measured single-channel
  conductances read from the literature, using the bulk resistivity of each
  recording solution.
* **Tunnels.** `crevice paths` is compared with CAVER 3.0.3 on haloalkane
  dehalogenases DhaA (1CQW) and LinB (1MJ5) and the acetylcholinesterase
  gorge (1EA5), from the catalytic residues.
* **Volumes.** `crevice cast` and `crevice cavities` are compared with
  pyKVFinder and POVME on T4 lysozyme L99A (1L83), FABP4 (2NNQ) and the GLUT1
  central cavity (4PYP), and the GLUT1 example trajectory's region volume
  with MDpocket.
* **Automatic probe.** The automatic enclosure-probe rule is rerun with
  its grouping tolerances, minimum support and probe ladder changed one at a
  time.

## Results in brief

* **Profiles.** Where CREVICE resolves a channel, its radii agree with HOLE's
  to 0.002–0.009 Å RMSD and the minimum radii to within 0.005 Å (8 of 11
  channels with `simple.rad`; the compared range is the part CREVICE
  resolves, which for aquaporin-1 is a short capped fragment that misses the
  narrowest point, see the README). Channels whose narrowest point is smaller than
  about 0.5 Å (the smallest automatic enclosure probe), such as potassium
  channel selectivity filters and the closed KcsA gate, are not resolved;
  HOLE still reports a plane-wise radius there. With default options only 3
  of 11 channels resolve; capped channels need `--lateral-exits`, and
  CREVICE's resolved range then stops at the cap, so a constriction beyond
  it (for example in the MscS cytoplasmic domain) is not in the profile.
* **Conductance.** The geometric estimate is 1.5 to 17.6 times the measured
  conductance on the four comparable channels CREVICE resolves (HOLE's own
  estimate: 0.2 to 12 times), and it fails for potassium channels. Treat it as
  an order-of-magnitude descriptor, not a prediction.
* **Tunnels.** The main tunnel of each enzyme is found by both programs in
  the same direction, with bottleneck radii within 0.03–0.11 Å. Narrower
  secondary tunnels that CAVER reports appear in CREVICE's list of routes left
  out by the automatic probe, or as routes with `--probe-radius 0.9`.
* **Volumes.** Swept-probe volumes (CREVICE `cast`, pyKVFinder at the same
  probe) agree within 3 % for the enclosed T4 lysozyme cavity and within
  20–25 % for larger, partly open pockets; tools that count probe centres or
  distance-cutoff points report several times smaller volumes by definition.
* **Automatic probe.** On the 9 channels resolved with the default rule,
  changing the minimum-radius tolerance changes nothing, and changing the
  axial-overlap fraction (0.4 or 0.8) changes the chosen probe for 2 channels
  without changing the minimum radius by more than 0.01 Å. The rule is
  sensitive to its support requirement: most channels resolve at only two or
  three ladder probes, so requiring 3 or 4 probes or cutting the ladder at
  2.0 Å leaves 3–7 of the 9 unresolved. No variant
  changed a reported minimum radius by more than 0.01 Å
  (`benchmarks/results/probe_sensitivity/`).
* **Run time.** On single structures, HOLE is 40–400 times faster than a
  profile with the automatic probe, which costs 3.5–8 times a fixed-probe
  profile. On a 1000-frame GLUT1 membrane simulation, the region volume took
  1665 s on one core at the definition's 0.2 Å grid and 370 s at 0.7 Å;
  MDpocket took 9326 s for the same frames. In those 0.4.0 runs CREVICE's
  peak memory grew with the frame count, to about 7.5 GB at 1000 frames,
  because every frame's contact networks were kept until the run ends; since
  0.4.0 they are released as each frame finishes, and the same 1000 frames
  (0.7 Å grid) peak at 0.57 GB. Trajectory profile timings against HOLE
  follow with a channel simulation (`benchmarks/README.md`, "Timing").

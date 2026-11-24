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

## Results in brief

* **Profiles.** CREVICE resolves all 11 channels with `simple.rad` (5 with
  default options; capped channels need `--lateral-exits`), covering 32–100 %
  of each protein's axial extent, including the KcsA (0.62 Å) and TWIK-1
  (0.51 Å) selectivity filters. Where both programs follow the same lumen the
  radii agree to 0.002–0.005 Å RMSD and the minimum radii to within 0.004 Å.
  HOLE's minimum within the protein lies inside CREVICE's range for 7 of the
  11; three of the others are in domains outside the transmembrane pore (the
  nAChR extracellular domain and the MscS cytoplasmic cage, beyond a capped
  end). Planes with two clearance maxima explain the larger differences (OmpF
  vestibules, one aquaporin-1 plane). The MthK pore domain is occluded on the
  axis in the curated structure: CREVICE's widest route leaves it sideways just
  above the occlusion, while HOLE stays on the axis and reaches 0 Å.
* **Conductance.** The geometric estimate is 1.06 to 12.7 times the measured
  conductance on the four comparable open channels (HOLE's own estimate: 0.22
  to 12.1 times), close to HOLE's where CREVICE covers the pore, and it fails
  for potassium channels. Treat it as an order-of-magnitude descriptor, not a
  prediction.
* **Tunnels.** The main tunnel of each enzyme is found by both programs in
  the same direction, with bottleneck radii within 0.03–0.11 Å. Narrower
  secondary tunnels that CAVER reports appear in CREVICE's list of routes left
  out by the automatic probe, or as routes with `--probe-radius 0.9`.
* **Volumes.** Swept-probe volumes (CREVICE `cast`, pyKVFinder at the same
  probe) agree within 3 % for the enclosed T4 lysozyme cavity and differ by
  19–24 % (CREVICE smaller) for larger, partly open pockets; tools that count probe centres or
  distance-cutoff points report several times smaller volumes by definition.
* **Run time.** On single structures (one process per run, one core, with the
  numba kernels, on an idle machine) a profile took 0.6 s for
  gramicidin, 1.3 s for the nicotinic receptor pore domain 1OED and 2.7 s for TWIK-1 with its
  lateral-exit and lateral-branch search, against 0.06–0.13 s for HOLE; about 0.3 s of each
  CREVICE run is start-up (`benchmarks/README.md`, "Timing"). Larger channels
  (KcsA, OmpF, GLIC) take 1–3 s. With the NumPy path forced (`CREVICE_DISABLE_NUMBA=1`) the same profiles take 1–31 s.
  On a 1000-frame GLUT1 membrane simulation, the region volume took
  1665 s on one core at the definition's 0.2 Å grid and 370 s at 0.7 Å;
  MDpocket took 9326 s for the same frames. In those 0.4.0 runs CREVICE's
  peak memory grew with the frame count, to about 7.5 GB at 1000 frames,
  because every frame's contact networks were kept until the run ends; since
  0.4.0 they are released as each frame finishes, and the same 1000 frames
  (0.7 Å grid) peak at 0.57 GB. On a 1001-frame simulation of the OmpF porin,
  the channel profile of one pore took 986 s on one core (933 frames
  resolved; 53 others had too little clearance at the fixed start point and
  15 a capped end, both refused by default) and HOLE,
  run frame by frame through mdahole2, 1187 s (`benchmarks/README.md`,
  "Trajectory timing").

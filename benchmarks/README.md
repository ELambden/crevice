# CREVICE benchmark suite

Reproducible comparisons of CREVICE with established tools on curated
structures: channel radius profiles against HOLE, geometric conductance
against measured single-channel conductances, enzyme tunnels against CAVER,
pocket volumes against pyKVFinder and POVME, a trajectory pocket volume
against MDpocket, the sensitivity of CREVICE's automatic enclosure probe, and
run times.

Everything here is a **software comparison on static structures (and one
short trajectory) with curated inputs**. Agreement with another program is
not a biological validation, and no setting of CREVICE was tuned to match a
comparator: CREVICE runs with its defaults unless a documented option is
listed, and every attempt (including failures) is recorded.

## Layout

| Path | Content |
| --- | --- |
| `inputs/channels.csv` | channel set: PDB entry, assembly, kept chains/residues, what was removed, axis and seed rules, why each entry was chosen |
| `inputs/tunnels.csv` | enzyme tunnel set (DhaA, LinB, acetylcholinesterase) with the start-point atoms |
| `inputs/pockets.csv` | pocket/cavity set with the ligand that marks each pocket |
| `inputs/conductance_literature.csv` | measured single-channel conductances, recording conditions, electrolyte resistivity and the source each value was read from |
| `fetch_inputs.py` | downloads the structures (`crevice fetch`) and RCSB metadata |
| `install_comparators.sh` | downloads/builds HOLE, CAVER (+ OpenJDK), fpocket/mdpocket, pyKVFinder and POVME into `<workspace>/tools`; nothing is vendored |
| `curate_channels.py` | writes the curated wall atoms, axis and seed of each channel |
| `channels_vs_hole.py` | matched CREVICE/HOLE profiles |
| `conductance.py` | geometric conductance estimates vs measurements |
| `tunnels_vs_caver.py` | CREVICE `paths` vs CAVER 3.0.3 |
| `volumes_vs_comparators.py`, `kvfinder_worker.py` | pocket volumes: CREVICE `cavities --cast`/`cavities` vs pyKVFinder and POVME |
| `mdpocket_vs_region.py` | GLUT1 trajectory pocket: `trajectory --region` vs MDpocket |
| `timing_vs_comparators.py` | run times of CREVICE and the comparators on the same inputs |
| `timing_trajectory.py`, `trajectory_worker.py` | run time, CPU time and memory over an MD trajectory: CREVICE `trajectory` (channel profile, region volume at several grid spacings) vs HOLE per frame (mdahole2) and MDpocket |
| `inputs/trajectory_glut1_excerpt.json` | system description for `timing_trajectory.py` on the public GLUT1 excerpt |
| `nongeom/` | verification of the non-geometric analyses (permeation, pf/pd, current, uncertainty, water, boundary residues) against exact answers and independent references; see `nongeom/README.md` |
| `make_tables.py` | collects the per-comparison CSVs into `results/validation_channels.csv` and `results/validation_summary.csv` |
| `results/` | the small result tables (one CSV per comparison) |

Downloaded structures and raw comparator outputs go to a work directory
outside version control: `<workspace>/benchmarks`, where `<workspace>` is the
first entry of `$CREVICE_WORKSPACE_DIR` or `.crevice` (every script takes
`--work-dir`).

## Running

```bash
bash benchmarks/install_comparators.sh            # once; see the licences listed in the script
python benchmarks/fetch_inputs.py                 # structures + RCSB metadata
python benchmarks/curate_channels.py
python benchmarks/channels_vs_hole.py --radius-set simple
python benchmarks/channels_vs_hole.py --radius-set bondi
python benchmarks/conductance.py --radius-set simple   # and --radius-set bondi
python benchmarks/tunnels_vs_caver.py
CREVICE_KVFINDER_PYTHON=... CREVICE_POVME_PYTHON=... CREVICE_POVME_SCRIPT=... \
    python benchmarks/volumes_vs_comparators.py
crevice fetch --example glut1-excerpt --cache-dir <work-dir>/examples
python benchmarks/mdpocket_vs_region.py
python benchmarks/timing_vs_comparators.py --repeat 3
python benchmarks/timing_trajectory.py benchmarks/inputs/trajectory_glut1_excerpt.json --repeat 3 --cpus 0
python benchmarks/make_tables.py
```

Comparator versions used for the committed results: HOLE 2.3.1 (Apache-2.0),
CAVER 3.0.3 (GPL-3.0) on OpenJDK 17, pyKVFinder 0.9.5 (GPL-3.0-or-later),
POVME 3.0.34 (MIT; Python 2.7 only), fpocket/mdpocket tag 4.1.5 (MIT; mdpocket aborts with a
memory error when started through a long path, so the script runs it through a
short symlink with relative arguments). CHAP was not run: it needs the legacy GROMACS
C++ API (GROMACS 2016-2021 built with the legacy headers), and the available
GROMACS 2022 installation does not provide it. MOLE was not run.

## Curation

Every entry is marked `curated`: its assembly, atom selection, axis or start
point and reason were fixed in the input lists before the comparison, with
provenance recorded by the scripts (RCSB metadata, source checksums, the
operator fit and its RMSD, the angle and offset to RCSB's own symmetry axis).
Curated means the inputs are fixed and documented, not that any computed
radius or volume is a validated biological value.

* Kept atoms: polymer-entity heavy atoms of model 1 (altloc blank or A);
  waters, ions, lipids, detergents and ligands are removed (listed per
  entry).
* Channel axes: the least-squares cyclic operator over all subunits
  (homo-oligomers), the normal of the plane of subunit centroids (the
  pseudo-symmetric acetylcholine receptor) or the backbone principal axis
  (gramicidin). Every axis agrees with RCSB's symmetry axis to within 1.1°
  (gramicidin: RCSB lists only the perpendicular dimer two-fold).
* Channel seeds: on the axis level with the protein centroid; where that point
  lies in a wide vestibule (GLIC, the full acetylcholine receptor, MscS) the
  seed is on the axis at a named pore-lining ring (M2 9′, MscS L105);
  monomer pores (aquaporin-1, OmpF) are seeded at named pore residues. A seed
  with less than 1 Å clearance is moved by a fixed, recorded rule.
* Radii: both tools in a comparison read the same radius definition (HOLE's
  `simple.rad` or Bondi for channels; Bondi, which equals CAVER's, POVME's and
  CREVICE's defaults for C/N/O/S, elsewhere).

## Results

The committed tables in `results/` are the reference. Headline numbers
(machine and dates are recorded in the CSV headers):

### Channel radius profiles vs HOLE (matched atoms, radii, axis, seed)

Per-plane maximal-sphere radius, HOLE interpolated at CREVICE's samples, over
CREVICE's resolved range trimmed by 2 Å per end. HOLE `simple.rad` in both
tools (the Bondi run gives the same picture, `channels_vs_hole_bondi.csv`).
Coverage is CREVICE's resolved range as a fraction of the atoms' axial extent;
the last two columns give HOLE's minimum within that extent and whether it
lies inside CREVICE's range.

| Entry | default options | options needed | CREVICE range (Å) | coverage | RMSD (Å) | max \|Δ\| (Å) | min radius CREVICE / HOLE (Å) | HOLE min in protein (Å) @ s | inside range |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1GRM gramicidin A | resolved | – | −14.6…14.6 | 1.00 | 0.004 | 0.014 | 1.200 / 1.199 | 1.20 @ −9.6 | yes |
| 1K4C KcsA (closed gate) | resolved | – | −37.1…23.2 | 0.95 | 0.004 | 0.034 | 0.618 / 0.617 | 0.62 @ 5.0 | yes |
| 1LNQ MthK pore | unresolved | unresolved | – | – | – | – | – / – | −0.01 @ −5.6 (global) | – |
| 1J4N aquaporin-1 | unresolved | lateral exits | −19.1…21.5 | 0.63 | 0.073 | 0.561 | 0.983 / 0.982 | 0.98 @ 7.0 | yes |
| 2OMF OmpF | resolved | – | −21.9…22.5 | 0.76 | 0.270 | 1.258 | 3.563 / 3.559 | 3.56 @ −0.2 | yes |
| 4HFI GLIC | resolved | – | −70.6…27.7 | 0.86 | 0.003 | 0.026 | 2.593 / 2.593 | 2.59 @ 11.8 | yes |
| 2BG9 nAChR (full) | unresolved | lateral exits | −23.0…27.7 | 0.32 | 0.005 | 0.041 | 2.566 / 2.564 | 2.18 @ −56.0 | no |
| 1OED nAChR pore | resolved | – | −23.7…22.7 | 0.84 | 0.004 | 0.040 | 3.114 / 3.113 | 3.11 @ −9.8 | yes |
| 2VV5 MscS (open) | unresolved | lateral exits | −29.9…41.4 | 0.59 | 0.002 | 0.017 | 6.030 / 6.030 | 2.06 @ 66.8 | no |
| 2OAU MscS (closed) | unresolved | lateral exits | −39.2…38.6 | 0.63 | 0.002 | 0.024 | 2.429 / 2.429 | 1.75 @ 75.6 | no |
| 3UKM TWIK-1 | unresolved | lateral exits | −25.1…16.5 | 0.53 | 0.004 | 0.020 | 0.505 / 0.504 | 0.50 @ 5.2 | yes |

Reading the table:

* The two programs measure the same quantity, so this checks the
  implementation, not biology. Where both follow the same lumen the per-plane
  radii agree to 0.002–0.005 Å RMSD and the minimum radii to ≤ 0.004 Å
  (Bondi: ≤ 0.005 Å except MthK).
* With default options CREVICE resolves 5 of 11 entries; with the documented
  options (`--lateral-exits` for capped ends) 10 (Bondi: 11). Every plane uses
  its own enclosure probe, so narrow filters and leaky wide regions resolve in
  the same channel: KcsA's 0.62 Å filter and TWIK-1's 0.51 Å filter are
  followed. The smallest traceable clearance is about 0.32 Å with 0.25 Å
  sections.
* HOLE's minimum within the protein lies inside CREVICE's range for 7 of the
  10 resolved entries. The other three are outside the transmembrane pore:
  the extracellular domain of the full nAChR (2BG9, s = −56 Å) and the
  cytoplasmic cage of MscS (2VV5, 2OAU), beyond the capped end where CREVICE
  leaves through lateral exits.
* Larger differences come from planes with more than one clearance maximum.
  In the OmpF vestibules (r > 5.8 Å) CREVICE's centres are exact maxima of the
  enclosed lumen that are up to 1.26 Å wider than the point HOLE's search
  settles on; in aquaporin-1 one plane at s = −12 Å holds two maxima, and
  CREVICE keeps the one continuous with its path. A higher maximum is taken
  only when it lies in the same lumen (the straight line between the two
  keeps 90 % of the lower clearance).
* MthK (1LNQ): the curated pore domain is occluded (HOLE gives 0.0 Å at
  s = −5.6 Å); with Bondi radii CREVICE resolves only the part above the
  occlusion.
* TWIK-1 (3UKM): the capped end stops where the two side portals leave the
  chamber under the cap (portal legs about 2.0 Å, against 1.0 Å for the route
  through the cap; the lateral-branch rule of `crevice profile`), so the
  resolved range ends at s = 16.5 Å rather than inside the cap.
* Every resolved entry is compared on its dense (~0.2 Å step) pass. The
  channel's extent comes from its mouths and does not depend on the sample
  count.

### Geometric conductance vs measured single-channel conductance

`results/conductance/conductance_vs_measured_simple.csv` (and `_bondi`).
`G_macro = 1/(ρ F)`, `F = ∫ ds/(π r²)` over CREVICE's resolved lumen, ρ of the
recording solution; HOLE's own `Gmacro` (whole HOLE profile) rescaled to the
same ρ. Ratios are estimate / measured, per pore (OmpF: trimer / 3).

| Entry | measured (pS) | conditions | CREVICE G_macro (pS) | ratio | lumen (Å) | HOLE G_macro (pS) | ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gramicidin (1GRM) | 10 | 150 mM NaCl | 33.4 | 3.3 | 29.3 | 31.9 | 3.2 |
| MthK (1LNQ) | 242 | 200 mM K⁺, −200 mV | unresolved | – | – | ~0 | ~0 |
| OmpF (2OMF, per pore) | 1367 | 1 M KCl | 2211 | 1.6 | 44.4 | 1897 | 1.4 |
| GLIC (4HFI) | 9.3 | 160 mM NaCl / 145 mM CsCl | 118.6 | 12.7 | 98.3 | 112.6 | 12.1 |
| MscS open (2VV5) | 1200 | 200 mM KCl + 90 mM MgCl₂ | 1266 | 1.06 | 71.3 (capped) | 263 | 0.22 |
| KcsA closed (1K4C) † | 97 | 100 mM KCl | 13.4 | 0.14 | 60.4 | 13.2 | 0.14 |
| MscS closed (2OAU) † | 1200 | as 2VV5 | 506 | 0.42 | 77.8 (capped) | 178 | 0.15 |

† closed structure vs open-channel measurement: not a like-for-like comparison.

The estimate is not a predictor of single-channel conductance: over the four
comparable open channels the ratio spans 1.06–12.7 (HOLE: 0.22–12.1). Where
CREVICE now covers the pore, its estimate is close to HOLE's (gramicidin,
OmpF, GLIC, KcsA); for MscS the two differ because HOLE integrates through the
cytoplasmic cage that CREVICE leaves through its lateral portals. For K⁺
channels the geometric picture fails outright (the filter is narrower than a
water molecule). Sources and conditions are in
`inputs/conductance_literature.csv`; every value was read in the cited
article.

On a 100 ns OmpF trimer simulation at 561 mV (`nongeom/md_systems/ompf.md`), the same
estimate on every frame's profile, at the experimental 1 M KCl resistivity, is 1.5–2.0
times the simulated per-pore conductance. Over 5 ns windows it follows the simulated
current only weakly (Spearman ρ = 0.32), and it did not change when one pore's K+ current
fell threefold.

### Enzyme tunnels vs CAVER 3.0.3

`results/tunnels/`: same atoms, start point and Bondi radii; CAVER defaults
(probe 0.9 Å), CREVICE `paths` with its automatic probe and with 0.9 Å.

| Enzyme | CAVER top cluster: bottleneck / length (Å) | CREVICE route in the same direction: bottleneck / length (Å) | angle |
| --- | --- | --- | --- |
| DhaA (1CQW) | 1.39 / 16.8 | 1.36 / 15.1 | 3° |
| LinB (1MJ5) | 1.51 / 14.5 | 1.40 / 13.7 | 1° |
| AChE gorge (1EA5) | 2.47 / 11.1 | 2.44 / 17.5 | 9° |

The main tunnel of each enzyme is found by both tools in the same direction
with bottleneck radii within 0.03–0.11 Å. CREVICE's automatic probe reports
one route per enzyme and lists the narrower ones it leaves out (DhaA's second
CAVER tunnel, 0.95 Å, appears there as a 0.93 Å left-out route; with
`--probe-radius 0.9` it is reported, 0.98 Å). CAVER splits the AChE gorge
mouth into two clusters that CREVICE merges into one route (CREVICE keeps
routes whose exits are closer than 12 Å as one). Lengths differ by
definition: CAVER ends at its molecular-surface shell, CREVICE at bulk or an
exposed node.

### Pocket and cavity volumes

`results/volumes/volumes.csv`: same atoms, ligand-centroid seed, Bondi radii;
each tool's own volume definition.

| Pocket | CREVICE cast (0.8 Å probe) | pyKVFinder probe-in 0.8 Å | pyKVFinder default (1.4 Å) | POVME (12 Å sphere) | CREVICE cavities (1.4 Å centres) |
| --- | --- | --- | --- | --- | --- |
| T4L L99A (1L83) | 385 | 375 | 196 | 65 | 31 |
| FABP4 (2NNQ) | 1042 | 1281 | 949 | 579 | 355 |
| GLUT1 central cavity (4PYP) | 3245 | 4272 | 1329 | 912 | open (none) |

Volumes in Å³. Swept-probe volumes (CREVICE `cavities --cast`, pyKVFinder) agree within
3 % for the enclosed T4L cavity at the same probe and within 20–25 % for the
larger, partly open pockets, where each tool bounds the open side
differently. POVME and `crevice cavities` count probe-centre or
distance-cutoff points and are several times smaller by construction.

### Trajectory pocket volume vs MDpocket (GLUT1 excerpt)

`results/trajectory_pocket/`: the 21-frame public GLUT1 excerpt and its
prepared candidate-site region; CREVICE `trajectory --region` (defaults) and
MDpocket (fpocket 4.1.5) with the reference region as the selected pocket,
frames aligned on the same protein CA atoms, lipids and ions as obstacles.

| | CREVICE region volume | MDpocket pocket volume |
| --- | --- | --- |
| mean ± SD over 21 frames (Å³) | 761 ± 47 | 555 ± 84 |

Ratio of means 0.73; Pearson r over frames 0.22. The two volumes are
defined differently (CREVICE: probe-accessible volume inside the region at
each frame; MDpocket: alpha-sphere volume of the pockets overlapping the
selected points), and on this excerpt their frame-to-frame changes are only
weakly correlated. Twenty-one frames 5 ns apart do not support a stronger
statement.

### Timing

`results/timing/`: one process per run, median of 3, one core for CREVICE
(`--workers 1`, single-threaded BLAS, with the optional numba kernels);
CAVER is a Java program and uses more than one thread (its CPU time is given).
Intel Core Ultra 7 268V, 8 logical CPUs, Linux, Python 3.13; measured on
3 Nov 2026 on an otherwise idle machine (1-minute load 0.7–2.1; 0.7 for the profiles). Times
include start-up (about 0.3 s to import CREVICE and its kernels); the numba
kernels add about 80 MB to CREVICE's peak memory.

| Case | Tool | median wall (s) | median CPU (s) | peak RSS (MB) |
| --- | --- | --- | --- | --- |
| 1GRM profile | HOLE 2.3.1 | 0.06 | 0.06 | 19 |
| | CREVICE profile, automatic probe | 0.61 | 0.61 | 171 |
| 1OED profile | HOLE | 0.13 | 0.13 | 19 |
| | CREVICE, automatic probe | 1.3 | 1.3 | 201 |
| 3UKM profile (lateral exits) | HOLE | 0.10 | 0.10 | 19 |
| | CREVICE, automatic probe | 2.7 | 2.7 | 248 |
| DhaA tunnels | CAVER 3.0.3 | 0.85 | 4.2 | 451 |
| | CREVICE paths, automatic probe | 7.5 | 7.4 | 246 |
| | CREVICE paths, 0.9 Å | 3.5 | 3.5 | 246 |
| LinB tunnels | CAVER | 0.90 | 4.6 | 444 |
| | CREVICE paths, automatic / 0.9 Å | 6.0 / 3.3 | 6.0 / 3.3 | 241 |
| AChE gorge | CAVER | 1.0 | 4.6 | 447 |
| | CREVICE paths, automatic / 0.9 Å | 14.3 / 6.0 | 14.3 / 6.0 | 329 |
| T4L cavity | pyKVFinder | 0.23 | 0.99 | 37 |
| | CREVICE cast | 2.1 | 2.1 | 184 |
| | POVME 3 | 18.3 | 18.2 | 12 813 |
| FABP4 pocket | pyKVFinder | 0.18 | 0.75 | 35 |
| | CREVICE cast | 2.3 | 2.3 | 208 |
| | POVME 3 | 25.1 | 25.0 | 12 800 |

HOLE is 10–27 times faster than a CREVICE profile, because CREVICE checks
that every slice is laterally enclosed at the smallest probe of a ladder and
joins neighbouring slices, and searches for lateral exits and wider lateral
branches at capped ends;
about half of CREVICE's time on 1GRM is start-up. CAVER and pyKVFinder are
faster than CREVICE's path search and cast. The 0.4 profile code took 2.4,
6.1 and 45.6 s for the same three channels (22 Oct 2026, `results/timing/` of
that release). Trajectory timings are in the next section.

### Trajectory timing

`timing_trajectory.py` times whole trajectory runs. A system is described by
a small JSON file (topology, trajectory, analysed selection, a prepared region
definition and, for a channel, its axis and a point inside it), so another
trajectory needs only a new file:

```bash
python benchmarks/timing_trajectory.py SYSTEM.json --frames 1000 --repeat 3 \
    --cpus 0 --spacings 0.35,0.5,0.7
```

| Case | Command |
| --- | --- |
| `crevice_profile` | `crevice trajectory --topology ... --analysis profile` with the system's axis and channel point (automatic per-slice probe in every frame) |
| `hole` | HOLE 2.3.1 on every frame through mdahole2: same atoms, frames fitted on the same CA atoms, HOLE `simple.rad`, same channel point and direction |
| `crevice_region` | `crevice trajectory --region` at the definition's grid spacing |
| `crevice_region_<s>` | the same region volume at grid spacing `s` Å (`trajectory --cavities --geometry-mode reference-region` with the definition's selections, probe, margin and axis) |
| `mdpocket_input` | fit of the frames on the same CA atoms and writing of the non-water atoms (MDpocket needs superimposed frames) |
| `mdpocket` | MDpocket 4.1.5 characterising the reference region (`--selected_pocket`, as in `mdpocket_vs_region.py`) on those frames |

Each command is one process, pinned with `taskset` (`--cpus`); CREVICE uses
`--workers 1` unless `--workers` is given, and single-threaded BLAS. Recorded
per run: wall time; user and system CPU time; peak resident memory of the
largest process (`os.wait4`) and of the whole process tree (summed resident
memory sampled from `/proc` every 0.2 s; pages shared between worker
processes are counted once per process, so with workers this is an upper
bound); CPUs, their maximum clock and workers; grid spacing; frames requested,
frames written and frames with a result (a resolved profile, a measured
volume); the 1-minute load before and after. `summary.md` holds the table
below; `summary.csv` adds seconds per frame and CPU/wall. HOLE is run through
mdahole2 (one PDB file and one HOLE process per frame), the usual way to run
HOLE on a trajectory; MDpocket reads the trajectory directly.

On processors that mix performance and efficiency cores, pin every run to
the same kind of core (the maximum clock of the pinned CPUs is recorded).

Results (`results/timing_trajectory/glut1-md-1001/` and
`glut1-md-1001-mdpocket/`), measured on 22 Oct 2026 with the 0.4.0 code on
an Intel Core Ultra 7 268V. Each command ran on one performance core (CPU 1
for CREVICE, CPU 2 for MDpocket; the two ran at the same time, and the
1-minute load was 0.9–2.3). The CREVICE rows are medians of 3 runs. MDpocket
was run once, because one run takes 2.6 hours.

| System (frames) | Case | Tool | grid (Å) | wall (s) | CPU (s) | frames/s | peak RSS (MB) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GLUT1, 1000 | region volume | CREVICE `trajectory --region` | 0.2 | 1665 | 1665 | 0.60 | 7570 |
| | | CREVICE, reference region | 0.35 | 607 | 606 | 1.65 | 7450 |
| | | CREVICE, reference region | 0.5 | 432 | 432 | 2.32 | 7403 |
| | | CREVICE, reference region | 0.7 | 370 | 370 | 2.70 | 7310 |
| | | MDpocket input (fit, write) | – | 9.3 | 9.2 | 108 | 187 |
| | | MDpocket 4.1.5 (1 run) | 1.0 | 9326 | 9323 | 0.11 | 115 |

The region volume is timed on a 1000-frame GLUT1 membrane simulation with a
prepared candidate-site region. The two volume definitions differ (see
"Trajectory pocket volume vs MDpocket" above); this table compares cost only.
At the definition's 0.2 Å grid, CREVICE measured the region 5.6 times as
fast as MDpocket characterised the same frames, and 25 times as fast at
0.7 Å. Both counts exclude MDpocket's input preparation.

CREVICE's peak memory hardly depends on the grid. In these 0.4.0 runs it
grows with the number of frames, because 0.4.0 keeps every frame's contact
networks and boundary record in memory until the run ends: about 7 MB per
frame for this system (0.45 GB at 20 frames, 2.4 GB at 300, 7.6 GB at 1000;
0.7 Å grid).

After 0.4.0 (unreleased), each frame's contact networks are released as soon
as the frame and all earlier frames have finished, and only a compact record
is kept. Measured on 23 Oct 2026 on the same machine with the same
reference-region command (0.7 Å grid, `--workers 1`, one performance core;
peak resident memory from `/usr/bin/time`, 1 MB = 1024² bytes; not part of
the timing table above, which was not rerun):

| Frames | 0.4.0 peak RSS (MB) | after the change (MB) |
| --- | --- | --- |
| 20 | 434 | 324 |
| 100 | 981 | 355 |
| 300 | 2402 | 395 |
| 1000 | 7307 | 568 |

Memory now grows by about 0.2 MB per frame. All output files are identical
to 0.4.0's.

The channel-profile rows (CREVICE `trajectory --analysis profile` against
HOLE 2.3.1 through mdahole2) are not measured yet. GLUT1 has no
through-channel, so they need a channel trajectory: the 100 ns aquaporin-4 (3GD8)
monomer trajectories are now available for them (`nongeom/md_systems/aqp4.md`).

## Not done

* CHAP: not built (see Comparator versions above).
* MOLE: not run (licence terms not reviewed for this suite).
* Literature conductance for TWIK-1 and the nicotinic receptor: no source was
  read in full text; the receptor structures are closed and aquaporin-1 does
  not conduct ions, so no ratio is given for them.

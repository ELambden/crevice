# Verification of CREVICE's non-geometric analyses

The main benchmark suite (`benchmarks/README.md`) compares geometry with other
programs. The scripts here check the analyses that are not geometry:
permeation counting, the water permeabilities pf and pd, the ion current,
the correlation-aware uncertainty of trajectory means, water contacts and
region membership, and boundary-residue attribution with residue contacts.

Three kinds of check are kept apart:

* **software verification against exact answers**: synthetic structures and
  trajectories built so that the right answer is known by construction or
  analytically;
* **independent cross-checks**: the same quantity from a separate
  implementation that does not call CREVICE (MDAnalysis, pymbar, emcee, a
  transcription of a published reference code, or a reference written here
  from the published definition, `nongeom_common.py`);
* **real-system comparisons**: MD trajectories compared with the same
  references, or with literature values. Agreement with a literature value is
  a comparison of one simulation with one experiment, not a validation of
  CREVICE or of the force field.

None of these results validates a biological interpretation.

## Files

| Script | Checks | Results (`benchmarks/results/nongeom/`) |
| --- | --- | --- |
| `nongeom_common.py` | synthetic-system builders, stochastic models with known answers, the independent reference permeation counter and collective coordinate, the exact MSD of n(t) for walkers on a ring | – |
| `permeation_counting.py` | scripted paths, ring walkers (exact q0 and D_n), single file, 3D pore with lumen variants | `permeation_cases.csv`, `permeation_ring.csv`, `permeation_single_file.csv`, `permeation_lumen_sensitivity.csv`, `permeation_counting_summary.json` |
| `permeation_statistics.py` | D_n estimator bias and interval coverage, pd interval coverage, ion current vs saving interval and its interval coverage | `dn_estimator.csv`, `dn_end_to_end.csv`, `pd_interval.csv`, `ion_current.csv`, `ion_current_coverage.csv`, `permeation_statistics_summary.json` |
| `uncertainty_coverage.py`, `uncertainty_reference_worker.py` | coverage of the trajectory-mean intervals on AR(1) and two-state Markov series; tau_int, g and blocking vs emcee, pymbar and Jonsson's reference code | `uncertainty_coverage.csv`, `uncertainty_reference.csv`, `uncertainty_coverage_summary.json` |
| `water_occupancy.py` | per-residue water contacts and region membership on an MD trajectory vs MDAnalysis | `hydration_contacts_glut1.csv`, `membership_bulk_glut1.csv`, `water_occupancy_summary_glut1.json` |
| `lining_residues.py` | boundary attribution on an analytic sphere, 4PYP lining vs distance and ligand definitions, residue contacts and DCCM vs MDAnalysis | `lining_synthetic.csv`, `lining_4pyp.csv`, `contacts_glut1.csv`, `lining_residues_summary_*.json` |
| `md_channel_analysis.py`, `inputs/md_systems.json`, `inputs/nongeom_literature.csv` | permeation, pf/pd and current of an MD channel trajectory, CREVICE vs the independent references, with channel-definition variants and literature values | written to `--out-dir` |

Fast checks of the helpers are in `tests/test_nongeom_benchmarks.py`.

## Running

```bash
export PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python benchmarks/nongeom/permeation_counting.py              # ~4 min
python benchmarks/nongeom/permeation_statistics.py            # ~5 min
python -m venv .venv-ref && .venv-ref/bin/pip install pymbar==4.0.3 emcee==3.1.6 scipy
python benchmarks/nongeom/uncertainty_coverage.py --reference-python .venv-ref/bin/python   # ~1 min
python benchmarks/nongeom/water_occupancy.py --topology system.gro --trajectory production.xtc --system glut1
python benchmarks/nongeom/lining_residues.py --parts ab       # 4PYP from <workspace>/pdb (crevice fetch 4PYP)
python benchmarks/nongeom/lining_residues.py --parts c --topology system.gro --trajectory production.xtc
python benchmarks/nongeom/md_channel_analysis.py aqp4_monomer_porewaters_v2 \
    --topology prod.tpr --trajectory prod.xtc --out-dir analysis/aqp4_monomer         # ~3 min, 1001 frames
python benchmarks/nongeom/md_channel_analysis.py ompf_trimer_efield \
    --topology step7_1.tpr --trajectory step7_1.xtc --out-dir analysis/ompf    # ~2 min, 1001 frames
```

Every script takes `--results DIR`; the stochastic ones take `--seed` and
`--quick` (small sizes). pymbar and emcee are used only by the reference
worker, run in a separate interpreter.

## Definitions compared

* **Permeation event**, compartment rule (CREVICE's default and the
  reference's `compartment` convention): a molecule remembers the bulk side it
  last occupied; a passage through the axial extent that ends on the other
  side is an event, attributed to the channel if at least one saved frame of
  the passage is in the lumen and listed as `outside_lumen` otherwise. This is
  the rule of GROMACS computational electrophysiology (`swap/swapcoords.cpp`;
  Kutzner et al., Biophys. J. 101, 809 (2011)) and a "complete transport from
  one reservoir to the other" (Zhu, Tajkhorshid and Schulten, Biophys. J. 86,
  50 (2004)). Strict rule (`--permeation-convention strict`, the reference's
  `strict`): the molecule must also stay in the lumen without interruption,
  and a frame outside the lumen laterally forgets the side.
* **Collective coordinate**: dn = sum of axial displacements of molecules
  inside the channel / L (Zhu, Tajkhorshid and Schulten, Phys. Rev. Lett. 93,
  224501 (2004)); pf = v_w D_n, pd = v_w q0.
* **Displacement current** (independent of any channel definition):
  I = sum_i q_i dz_i / (dt L_z) over every ion in the box (Aksimentiev and
  Schulten, Biophys. J. 88, 3745 (2005)).
* **Water contact**: water oxygen within 3.5 Å of any heavy atom of a residue
  (same definition in CREVICE and in the MDAnalysis reference).
* **Lining residue**: CREVICE assigns each triangle of the cast surface to the
  nearest atom surface within 1.5 Å. The distance-based reference counts every
  residue with any heavy-atom surface within 1.5 Å + h/2 of an occupied voxel,
  so it includes residues that are occluded by nearer ones.

## Results

Machine, versions and seeds are recorded in the summary JSON files.

### Permeation counting

* **Scripted catalogue** (21 paths: full crossings up and down, entry and
  exit through the same mouth, 10 re-entries at a mouth without crossing, an
  immediate re-crossing, two crossings through periodic images, a periodic
  loop outside the radius, censored starts and ends, lateral exit inside the
  channel, lateral approach at the mouth, a passage beside the lumen, a 16 Å
  jump between frames, exact boundary positions, a channel centred on the
  periodic face, frame strides 1/3/6, a two-frame trajectory): every count
  matches the hand-written answer and the reference counter under both
  conventions. They differ only where a molecule leaves the lumen sideways
  inside the extent or enters the extent beside the lumen.
* **Independent walkers on a periodic line** (100 walkers, 20 000 frames,
  L = 10 Å, period 20 Å, 4 replicates per step size): 52 of 54 runs give
  event lists identical to the reference counter; the 2 differing events are
  both a molecule exactly half a period from the channel centre (a float32
  tie that the two programs round to opposite images). The rate q0 matches
  N D / ((L + 2βσ) Λ) (Siegmund's discrete-monitoring correction) to
  0.97–1.00 for σ/L = 0.025–0.1 and 0.91 at σ/L = 0.01; it is 0.84–0.95
  of the continuum value as σ grows, so pd depends on the saving interval.
* **D_n** recovers the exact value implied by its fit window: for this model
  the MSD of n(t) is known in closed form (`ring_n_msd`), and the fitted D_n
  divided by the window's expected value is 0.86–1.06 (single-run SD
  0.07–0.19). When the window (lags 1 frame to a tenth of a segment) is short
  compared with the relaxation of n(t), the window's expected value exceeds
  the long-time D_n by up to 1.78×.
* **Single file** (16 point particles with hard-core reflection, same
  trajectories relabelled): D_n is identical to that of independent walkers
  in every run and saving interval, as theory requires, while q0 drops and
  pf/pd is 4.4–6.6 (independent walkers 0.34–0.68). From the same trajectory
  saved every 1, 2, 5, 10 or 20 steps, D_n changes by < 2 % and q0 by up to
  −30 %.
* **3D pore in an impenetrable membrane** (600 Brownian particles, pore
  radius 4 Å, membrane 10 Å thick; every membrane transit passes through the
  pore, 722 transits at 1 ps resolution; `permeation_lumen_sensitivity.csv`).
  With the extent equal to the membrane, the compartment rule counts all 722
  for lateral radii from 1 Å inside the pore wall to 2 Å outside it. The
  strict rule counts 722 with the exact pore radius and 1 with a radius 1 Å
  smaller. With the extent reaching 2 or 4 Å into the bulk, the compartment
  rule counts 0.81 and 0.76 of the membrane transits; the difference from 1 is
  re-crossings inside the longer extent, not lost attribution
  (`compartment_over_axial_only` = 1.00). The strict rule counts 0.02–0.54.
  Saved every 10 and 50 ps, the matching cylinder counts 0.90 and 0.79 of the
  1 ps transits under both rules (passages faster than the saving interval
  are missed). At 50 ps a radius 1 Å inside the pore lists 85 transits as
  `outside_lumen` (no saved frame in the narrow lumen), and an extent 4 Å into
  the bulk lists 222–292. Every CREVICE event list, including the
  `outside_lumen` passages, equals the reference's under the same convention
  and definition.

Counts therefore depend on the axial extent, which defines what a transit is;
under the compartment rule the lateral radius only attributes transits to a
pore and can be generous. Use an extent that ends where the pore meets the
bulk, a radius that contains the whole lumen (and stays smaller than the
distance to the next pore), check `transits_outside_lumen`, and report the
variants (`md_channel_analysis.py` does this).

### pf, pd and current estimators

* `collective_diffusion` on random walks with known D_n: mean estimate
  1.00–1.01 of the truth for 201–20 001 frames (relative SD 0.17–0.21 for one
  trajectory) and interval coverage 0.94–0.95. A bounded component of n(t)
  (SD 1, relaxation time 0.5–20 frames) biases short trajectories (+8 to +31 % at
  201 frames, 0 to +7 % at 1001 frames, < 2 % from 5001 frames; +226 % at
  201 frames when that component is large).
* Through `analyze_permeation` (40 replicates of 100 walkers, 5000 frames):
  D_n / window expectation 0.95, pf interval coverage 0.95 of that value.
* pd interval, before (exact Poisson, Garwood) and after (Garwood widened by
  the dispersion index of the counts in 10 time blocks, floor 1), 300
  replicates each (`pd_interval.csv`):

  | Events | dispersion index (across replicates / mean within-run) | coverage before | after | width after / before |
  | --- | --- | --- | --- | --- |
  | independent walkers | 0.78 / 0.98 | 0.973 | 0.977 | 1.07 |
  | single file | 2.01 / 1.93 | 0.857 | 0.953 | 1.40 |
  | gated, mean dwell T/100 | 3.07 / 3.21 | 0.743 | 0.940 | 1.78 |
  | gated, mean dwell T/20 | 12.1 / 10.2 | 0.390 | 0.913 | 3.29 |

  The within-run block estimate underestimates the dispersion when the
  correlation time approaches a block (dwell T/20 = half a block).
* Ion current (cations and anions drifting at the drift speed of a
  0.05 V/nm field, 30 Å channel in a 100 Å box, exact I = 2.40 nA): count
  current 1.010–1.011 of exact and collective-coordinate current 1.016 for
  saving intervals of 1–100 ps (one trajectory), 1.007 at 200 ps; at 500 ps
  the count current drops to 0.83 (8 % of crossings unresolved) and the
  collective current to 0.94.
* Current interval coverage over 40 replicates at 10 ps
  (`ion_current_coverage.csv`):

  | Drift | net-charge dispersion index | Poisson (before) | block bootstrap (before) | primary interval (after) |
  | --- | --- | --- | --- | --- |
  | constant | 0.19 | 1.00 | 0.925 | 1.00 |
  | gated, mean dwell T/100 | 6.9 | 0.525 | 0.875 | 0.95 |
  | gated, mean dwell T/20 | 25 | 0.25 | 0.975 | 0.975 |

  Under constant drift the net charge is less variable than Poisson counts
  (re-crossings cancel), so the interval stays at the conservative Poisson
  floor.

### Uncertainty of trajectory means

See the table in `docs/methods/uncertainty.md`. The studentized batch
interval covers 0.93–0.98 once N/τ ≥ about 50 and is withheld for most
shorter series; the moving-block bootstrap covers 0.87–0.96 and mean ±
1.96 blocking SE 0.82–0.96 at N/τ ≥ 50. tau_int equals emcee's
`integrated_time` to 1e-15 (after clipping emcee's values below 1, which
CREVICE clips), g agrees with pymbar's `statistical_inefficiency` within a
few per cent except for series much shorter than τ (pymbar's estimator
truncates differently), and the blocking SE equals Jonsson's reference code
on 160 series of 2^d frames.

### Water

On a 100 ns GLUT1 membrane trajectory (CHARMM, 27 173 waters; 21 frames
every 5 ns), the per-residue water-contact counts of all 492 residues are
identical to MDAnalysis `capped_distance` in all 10 332 residue-frame pairs
(23 547 contacts). Region membership of a 8 Å lattice sphere in bulk water
that straddles the periodic x face: member identities equal an independent
voxel count in 100 of 101 frames (the other differs by one oxygen exactly on
a voxel face, 7.25 Å = 14.5 voxels, rounded differently), and the mean
member density 0.03294 Å⁻³ is 0.996 of the bulk oxygen density in a slab at
the same height (0.03306 Å⁻³).

### Boundary residues and contacts

* Analytic sphere (radius 6 Å, wall split into 5 polar bands with exact area
  fractions; partner shell on another chain; one distant residue): roles all
  correct; per-residue boundary fraction within 0.004 (h = 0.25 Å), 0.006
  (0.5 Å) and 0.012 (0.7 Å) of the exact fraction and of the reference
  (nearest-surface attribution of points on the exact sphere). The marching-
  cubes surface of the binary sphere is about 4 % larger than the sphere
  (area is not used for the fractions).
* Partners are "nonlocal" by position in the observed residue order (another
  chain, or more than two residues apart), not by residue number: a residue
  numbered 100 written directly after residue 5 counts as a sequence
  neighbour.
* 4PYP (`cavities --cast`, 3245 Å³, 62 lining residues): all 62 are also
  lining residues under the distance-based definition (74 residues, Jaccard
  0.84; the other 12 are within the distance but occluded, all but one at
  gaps of 1.46–1.73 Å), and all 12 residues within 4 Å of the bound detergent
  BNG are lining residues.
* GLUT1 trajectory (21 frames): residue-pair contact occupancy at 4.5 Å
  identical to MDAnalysis for all 2770 pairs; centroid DCCM within 1.4e-14.

### MD systems

Each analysed simulation has a record in `md_systems/`: its build and settings (checked against the run input), the
definitions used, results with intervals, convergence and the next simulations recommended.

* `md_systems/ompf.md`: OmpF trimer, 1 M KCl, static 0.05 V/nm (561 mV), 100 ns. Trimer conductance 3.04
  (2.86–3.21) nS from the displacement current and 2.98 (2.77–3.19) nS from counted crossings; CREVICE's event lists
  equal the reference counter's for every definition. The three pores differ (0.81–1.14 nS each; K+/Cl- current ratio
  0.73–2.40), so the conductance is converged within the run but not across pore states.
* `md_systems/aqp4.md`: aquaporin-4 (3GD8), three builds: the pore-water monomer (His201 HSP, 100 ns), the first
  monomer (His201 HSD, 100 ns) and the tetramer (His201 HSD, static 0.05 V/nm, 183.8 ns of a planned 1 µs). The
  pre-registered cylinders proved wrong for these builds and are kept as the record; the `*_v2` entries of
  `inputs/md_systems.json` correct them. Water events and pf/pd per pore with intervals, zero ion events (upper
  bounds), ar/R geometry and gate state over time. pf and pd are not converged (95 % half-widths per tetramer pore
  38–53 % for pd and 63–110 % for pf); the tetramer's ar/R radius is.

## Known limitations found by these checks

* n(t) loses the part of a frame step that wraps through the periodic face
  opposite the channel and lands inside it (minimal case in
  `permeation_cases.csv`, row `n_entry_through_opposite_periodic_face`). It
  needs a frame step longer than half the bulk path, (box − L)/2, and does
  not arise at MD saving intervals in membrane boxes.
* The dispersion-widened intervals assume events correlated over less than a
  tenth of the run (pd coverage 0.91 for a gate with a dwell of half a
  block); under constant drift the current interval is conservative (Poisson
  floor). In the OmpF trimer one pore changed its K+ rate for more than 60 of
  100 ns, so its within-run interval understates the spread between pores.
* The D_n fit window is fixed at lags of 1 frame to a tenth of a segment; its
  bias for a slowly relaxing n(t) is not flagged.

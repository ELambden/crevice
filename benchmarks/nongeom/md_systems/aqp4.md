# Aquaporin-4 (3GD8): two monomer builds and a tetramer under an applied field

The definition keys are in `../inputs/md_systems.json`. The driver is `../md_channel_analysis.py`. Each system has a
pre-registered entry (`aqp4_monomer_porewaters`, `aqp4_monomer`, `aqp4_tetramer_efield`), kept unchanged as the
record, and a corrected entry with the suffix `_v2`. The headline results below use the `_v2` entries.

Each build is one simulation. Comparing it with an experiment means comparing one simulation with one experiment. It
is not a validation of CREVICE or of the force field, and no functional conclusion is drawn. Every interval below is a
95 % interval. Residue numbers are those of 3GD8.

## Systems

| | Pore-water monomer | First monomer | Tetramer |
| --- | --- | --- | --- |
| Keys | `aqp4_monomer_porewaters(_v2)` | `aqp4_monomer(_v2)` | `aqp4_tetramer_efield(_v2)` |
| Protein | 3GD8 chain A, residues 32-254, NTER/CTER, no ligands | the same | 3GD8 chains A-D (OPM tetramer), residues 32-254 each |
| Histidines | **His201 HSP**; every other His HSD | all HSD (**His201 HSD**) | all HSD (**His201 HSD**) |
| Crystal waters | HOH 278-281 (pore waters) kept as TIP3 | none | none |
| Orientation at production start | pore axis 0.0° from the membrane normal (OPM orientation) | 19.9° (CHARMM-GUI reoriented the protein); settles within ~20 ns | four-fold axis 24.9° (reoriented); settles to 2-3° within ~20 ns |
| Membrane | CHL1 30, DOPC 120, DOPE 120, DOPG 30 | the same | CHL1 44, **POPC** 176, DOPE 176, DOPG 44 |
| Solution | 19 373 TIP3; 80 Na+, 51 Cl- | 20 497 TIP3; 85 Na+, 55 Cl- | 37 008 TIP3; 145 Na+, 101 Cl- |
| Atoms | 99 791 | 103 171 | 179 978 |
| CHARMM-GUI job | 9283484867 | 9244985212 | 9245173890 |
| Field | none | none | **static** E0 = 0.05 V/nm along z (`electric-field-z = 0.05 0 0 0`; ω = 0, σ = 0), about 0.54 V across the mean box height (10.81 nm) |
| Production | 100 ns, complete | 100 ns, complete | **183.8 ns analysed of a planned 1 µs**; the job stopped at 183.9 ns after 48 h and is being continued |

The solutions are 0.15 M NaCl plus neutralising Na+, with CHARMM TIP3P water.

The following apply to all three systems:

- **Force field:** CHARMM36m protein with CHARMM36 lipids, ions and TIP3P water. The systems were built with
  CHARMM-GUI Membrane Builder, with the topology from its FF-Converter. The release is as stated by the systems'
  author.
- **Engine:** GROMACS 2021.5. Time step 2 fs, LINCS on bonds to hydrogen.
- **Electrostatics and cut-offs:** PME with rcoulomb 1.2 nm and Fourier spacing 0.12 nm. Lennard-Jones cut-off
  1.2 nm with a force switch from 1.0 nm.
- **Ensemble:** 310 K with v-rescale (τ 1 ps; groups SOLU, MEMB, SOLV). 1 bar with C-rescale, semi-isotropic (τ 5 ps).
- **Centre-of-mass removal:** `comm-grps SOLU_MEMB SOLV`.
- **Preparation:** minimisation, then the six CHARMM-GUI equilibration stages. Production runs without restraints.
- **Output:** frames every 100 ps (1001 frames per monomer; 1839 for the tetramer), all atoms, xtc precision
  0.001 nm.

The tetramer production restarted from the same equilibrated state as an earlier production that is not used. In that
earlier run, `electric-field-z = 1 0.050 0 0` meant a 1 V/nm field oscillating at 0.05 ps⁻¹, and the system collapsed.

These settings were checked against each run's input (`prod.tpr`, read with `gmx dump`), the processed `.mdp`, the
topology, the run log, the trajectory and the build notes.

## Definitions analysed

- **Pores.** Each pore is a cylinder 22 Å long, from 12 Å below to 10 Å above the centroid of the Asn97/Asn213 ND2
  atoms of the NPA motif, with a 3 Å radius. Frames are fitted on protein CA to the first analysed frame. Crossings
  are counted with the compartment convention; the strict convention and radius/extent variants are reported as
  sensitivities.
- **Two pre-registered choices proved wrong for these builds.** The pre-registered entries are kept as the record and
  their results are listed in the evidence tables. The corrections:
  1. **Cylinder position.** The NPA ND2 centroid lies 2.5-3.6 Å to one side of the line along which the pore waters
     run. The 3 Å cylinder around it held only part of the lumen. In the pore-water monomer, the mean occupancy was
     3.5 waters, the strict convention found 1 event, and pf changed 28-fold between radii of 2 and 6 Å. The `_v2`
     entries shift each cylinder (`point_offset_A`) onto the median position of the pore waters between the NPA and
     the ar/R filter. The compartment event counts hardly change (57 vs 59 events).
  2. **Tetramer axis.** The tetramer starts 24.9° tilted, and every frame is fitted to the first. In frame 0, z is
     about 24° from the pore axes, so the pre-registered slab cut the membrane obliquely. That gave 1.45 million
     transits outside the lumen and a meaningless pf of order 1e-12 per pore. `aqp4_tetramer_efield_v2` discards the
     first 20 ns, during which the protein reorients. It uses the four-fold axis at frame 200, which is the normal of
     the plane of the four NPA centroids, 5.5° from z.
  - The monomers have no initial frames discarded. The pore-water monomer is not stationary, as described below. The
    first monomer's gate closes during its tilt relaxation.
- **Cross-check.** Every complete event equals that of the independent reference counter in
  `../nongeom_common.py`, for every definition, convention and species.
- **Code version.** The results use the collective-coordinate fix made on 28 Oct 2026. A step over the whole channel
  extent now counts towards n(t) only through the lumen; before the fix, pf was 5-19 % higher in these runs.

## Results

### Water permeation

These use the `_v2` definitions and the compartment convention.

- **pd** uses the quasi-Poisson interval: the Poisson interval widened by φ, the dispersion index of event counts
  over 10 time blocks.
- **pf** uses the 5-segment t interval of D_n, clipped at 0. The fit uses lags from 1 frame to a tenth of a segment
  (0.1-2 ns for the monomers, 0.1-3.2 ns for the tetramer).
- **Strict events** are the events counted with the strict convention.
- **Occupancy** is the mean number of waters in the cylinder.

| System, pore | Time | Events up/down | φ | pd (1e-15 cm³/s) | pf (1e-14 cm³/s) | pf/pd | Strict events | Occupancy |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Pore-water monomer | 0-100 ns | 35/24 | 3.9 | 8.8 (4.9-14.5) | 2.27 (0-4.87) | 2.6 | 8 | 7.3 |
| First monomer | 0-100 ns | 10/12 | 4.7 | 3.3 (1.0-7.9) | 1.01 (0-2.86) | 3.1 | 6 | 6.9 |
| Tetramer A | 20-183.8 ns | 49/35 | 4.8 | 7.7 (4.5-12.2) | 4.34 (0-8.92) | 5.7 | 66 | 7.4 |
| Tetramer B | 20-183.8 ns | 73/48 | 4.6 | 11.0 (7.2-16.1) | 3.68 (0-7.75) | 3.3 | 84 | 7.6 |
| Tetramer C | 20-183.8 ns | 12/11 | 1.1 | 2.1 (1.3-3.2) | 0.66 (0-1.54) | 3.1 | 5 | 6.5 |
| Tetramer D | 20-183.8 ns | 21/36 | 4.1 | 5.2 (2.8-8.7) | 1.37 (0.26-2.47) | 2.6 | 25 | 6.8 |
| Tetramer, per monomer (4 pores pooled) | 655 pore-ns | 155/130 | 24.2 between pores | 6.5 (3.4-11.6) | 2.51 (0-5.33) | 3.9 | 180 | 7.1 |

**How the tetramer rows were pooled**

- **pd.** The pooled interval uses φ = 24.2, the dispersion between pores.
- **pf.** The pooled value is the mean of the four pores; its interval is the t(3) interval of those four values.

**Pores and sensitivity**

- **The four pores of the tetramer do not share one rate** (χ² p = 1e-15).
- **No net water flux** in the field direction: 155 up, 130 down.
- **Sensitivity of the event counts.**
  - Radius 2 or 4 Å changes no event count.
  - An extent 2 Å longer at each end gives 7-22 % fewer events per tetramer pore.
  - In the monomers, an extent 2 Å shorter at each end gives up to 10 % more.

**Ions: zero Na+ and Cl- events in every pore of every system.** The 95 % upper bounds per pore and species:

- **Tetramer** (20-183.8 ns, 0.54 V across the box): fewer than 0.023 events per ns, current below 3.6 pA,
  conductance below 6.7 pS.
- **Monomers** (no field): fewer than 0.037 events per ns.

These bounds are set by the length of the run only.

**Comparison with experiment and with another simulation**

- The tetramer pf per monomer, 2.5e-14 cm³/s, is 0.10 times the experimental 2.4e-13 cm³/s (Yang and Verkman, J.
  Biol. Chem. 272, 16140 (1997)). That measurement is on rat AQP4 in Xenopus oocytes, normalised to expression and
  referenced to an AQP1 pf of 6e-14 cm³/s.
- It is 0.34 times one simulation estimate, 7.4e-14 cm³/s (Hashido, Kidera and Ikeguchi, Biophys. J. 93, 373 (2007);
  5 ns × 2 copies × 2 runs, channel length 16 Å).
- The pore-water monomer gives 0.095 and 0.31 times these values.
- Our pf is not converged (see Convergence). The experimental value is a relative one.

### Pore geometry and the ar/R filter over time

**How the pores were profiled**

- **Method.** CREVICE trajectory profile with the default per-plane probe, default radii, a 10 Å search radius and
  lateral exits.
- **Frames.** Every frame of the pore-water monomer; every 2nd frame of the first monomer and of the tetramer
  (0-183.8 ns).
- **Seed.** The NPA centroid itself has no clearance in most frames. The seed is therefore the point near the
  pore-water line with the largest clearance across frames.
- **Unresolved frames.** Profiles that do not resolve are recorded as unresolved, not as zero. In these runs they are
  ends capped without a lateral exit.

| System | Resolved | Profile reaches the ar/R | ar/R minimum radius (Å) | NPA minimum (Å) | ar/R ≥ 1.0 Å |
| --- | --- | --- | --- | --- | --- |
| Pore-water monomer | 954/1001 | 75 % of frames | 1.23 ± 0.39 | 1.77 | 69 % |
| First monomer | 479/501 | 62 % | 0.92 ± 0.33 | 1.81 | 38 % |
| Tetramer, pores A-D | 86-91 % | 86-91 % | 1.16 ± 0.09 in every pore | 1.82-1.89 | 96-97 % |

The monomer ar/R values come from frames whose profile reaches the ar/R. Those are mostly open frames, so the
statistics lean towards the open state.

**Pore-water monomer (His201 HSP)**

- Water events fall from 21, 18 and 15 per 20 ns (0-60 ns) to 5 and then 0 (60-100 ns).
- Over the same time the ar/R narrows: its minimum radius is 1.32-1.47 Å in 0-60 ns, then 0.92 and 0.82 Å. The
  single file of water breaks 3-9 Å above the NPA.
- Arg216 approaches the Phe77 ring (CZ to ring centroid 6.6 → 4.6 Å). It never contacts His201: the CZ to NE2
  distance stays at 7-8 Å.
- Arg216 χ2 visits the trans rotamer in 21 % of frames, without closing the pore.

**First monomer (His201 HSD)**

- From about 14 ns, Arg216 turns to the trans rotamer (χ2 −88° → −175°). Its guanidinium contacts His201 NE2
  (N…NE2 < 3.2 Å in 71 % of frames), and the ar/R closes.
- It partly reopens at 87-95 ns. The closed fraction is 0.78.

**Tetramer (His201 HSD)**

- The ar/R never closes in 4 × 184 ns. Arg216 χ2 is −90 ± 8° and CZ to His201 NE2 is 6.0 ± 0.3 Å in all four
  chains.
- In pores A, B and D the narrowest point moves for tens of nanoseconds at a time from the ar/R (about 7.7 Å above the
  NPA) to the cytoplasmic vestibule (about 7 Å below). There it is lined by His95, Ile174, Ile193 and Val85.
- 20 ns windows with this vestibule constriction have 0-9 water events and 5.7-6.6 waters in the pore, against 14-27
  events and 8.5-9.3 waters without it.
- Pore C has a low event rate without that constriction; this is not explained here.
- These are descriptive associations in one trajectory.

**Closing rates**

- Comparing closings per open time (exact conditional binomial test), the first monomer closed twice in 17 ns of open
  time, while the tetramer did not close in 736 pore-ns: one-sided p = 0.0005.
- Against the pore-water monomer (0 closings in 100 ns), p = 0.022.
- Tilt, oligomeric state, His201 tautomer and starting waters all differ between the builds, so these differences
  cannot be assigned to any one of them.

**Agreement with HOLE 2.3.1** (same radii, `simple.rad`, 0.25 Å sections)

- **Pore-water monomer** (every 50th frame): 19 of 21 frames resolved, median profile RMSD 0.019 Å.
- **Tetramer** (every 200th frame × 4 pores): 34 of 40 resolved, median RMSD 0.010 Å. The median minimum-radius
  difference is 0.002 Å.
- **Location of HOLE's minimum.** HOLE's in-protein minimum lies inside CREVICE's profile in 18/19 and 33/34 frames.

## Convergence

**Equilibration**

- Discard the tetramer's first 20 ns: it reorients in the membrane from 24.9° to 2-3°.
- The pore-water monomer is not stationary. There is one ar/R transition at about 60 ns, and the event rate before
  20 ns differs from that after it (p = 0.005). No start time makes this run stationary.

**Event rate and pd: not converged.**

- The tetramer's relative half-widths are 38-53 % per pore and about 57 % pooled.
- The counts are correlated over tens of nanoseconds: in the pore-water monomer the dispersion index rises from 1.8 for
  2 ns blocks to 6.8 for 20 ns blocks.

**pf: not converged.**

- Relative half-widths are 63-110 % per pore.
- pf depends on the fit window. Fitting lags of 2-20 ns instead of 0.1-1 ns changes it by −82 % to +29 % per tetramer
  pore, by −35 % in the pore-water monomer and by −85 % in the first monomer. The fall with longer lags points to a
  bounded component of the collective coordinate n(t); the long-lag fits are themselves noisy.

**pf/pd: not converged.** It is 2.6-5.7 per pore, with both uncertainties combined.

**ar/R geometry**

- **Tetramer: converged within this run.** The ar/R minimum is 1.16 ± 0.09 Å with a blocking standard error of
  0.004-0.007 Å, and the four pores agree to 0.01 Å.
- **Tetramer overall minimum: not converged.** It includes the His95 constriction, which switches over tens of ns:
  correlation time 12-16 ns, 10-14 effective samples in pores A and B.
- **Monomers: not converged.** 9-17 effective samples, and one transition per run.

**Gate-state populations: not converged.**

- The pore-water monomer had no closings (closing rate below 0.037 per ns).
- The first monomer had 3 transitions.
- The tetramer had no closings (below 0.005 per ns per pore).
- The tetramer's vestibule-constriction fraction ranges from 0 to about 0.8 between pores.

**Saving interval: cannot be judged from these data.**

- Thinning the 100 ps frames to 200-1000 ps adds ambiguous jumps of bulk water across the periodic box. Unresolved
  crossings go from 9 to 35 472, and spurious events are added.
- How many short passages the 100 ps frames miss needs a segment saved more often.

**Projection**

The projection assumes the system is stationary after 20 ns, and that pores and replicas are exchangeable and
independent. The relative half-width then scales as z·a/√T, where T is the total pore-time. a is the relative standard
deviation times √(pore-ns):

- For pd, a = √(φ/rate) = 7.5, with φ = 24.2 between pores and 0.44 events per pore-ns.
- For pf, a = 9.1, from the spread between pores.

| Target, 95 % relative half-width | pd | pf |
| --- | --- | --- |
| ±20 % | 1.3 µs of tetramer | 2.0 µs |
| ±10 % | 5.3 µs | 8.0 µs |
| At the planned 1 µs (20 ns discarded) | ±23 % | ±29 % |

If the pores were as consistent as the within-pore statistics suggest (φ of about 4), ±10 % on pd would need about
1 µs of tetramer.

## Planned and recommended

1. **The tetramer is being continued to 1 µs.** Rerun `aqp4_tetramer_efield_v2` on the full trajectory, with pf fit
   lags extended to 5-10 ns.
2. **For about ±10 %, add independent 1 µs tetramer replicas.** Start them from the equilibrated state with new
   velocities: about 5 for pd and 8 for pf. About 4 would give roughly ±15 % on both. Because the pores switch
   slowly (dwell 50-100 ns), independent sampling matters more than length.
3. **Saving interval.** Save one 20 ns segment of the tetramer every 10 ps (about 1.4 GB) to measure how many
   passages 100 ps frames miss.
4. **His201 tautomer and oligomeric state.** These effects can only be separated with an HSP201 tetramer, or 3
   replicas of about 200 ns of each monomer build.

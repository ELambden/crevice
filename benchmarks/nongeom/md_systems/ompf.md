# OmpF trimer in 1 M KCl under an applied field

The definition key is `ompf_trimer_efield` in `../inputs/md_systems.json`; the driver is `../md_channel_analysis.py`.

This is one 200 ns simulation, run as two continuous 100 ns parts. Comparing it with an experiment means comparing one
simulation with one experiment. It is not a validation of CREVICE or of the force field. Every interval below is a 95 %
interval.

## System

| | |
| --- | --- |
| Protein | OmpF (PDB 2OMF), biological trimer. Chains PROA, P2PROA and P3PROA hold residues 1-340 each, with NTER/CTER termini. Each chain has one His, set as HSD. Each chain carries a net charge of -12 e. |
| Orientation | OPM; tilt 0.01°. The N and C termini (periplasmic side) lie below the bilayer midplane, so **+z is extracellular**. |
| Membrane | 354 DOPE and 117 DOPG (3:1). The upper leaflet has 180 DOPE and 58 DOPG; the lower leaflet has 174 DOPE and 59 DOPG. |
| Solution | 47 179 TIP3 waters (CHARMM TIP3P), 1001 K+ and 848 Cl-. That is 1.0 M from the Cl- count. In the bulk slab (more than 25 Å from the bilayer midplane), the run averages 1.11 M K+ and 1.06 M Cl-. The excess K+ neutralises DOPG and the protein. The system has 219 592 atoms and zero net charge. |
| Force field | CHARMM36m protein with CHARMM36 lipids, ions and TIP3P water, built with CHARMM-GUI (topology from its FF-Converter); the release is as stated by the system's author. |
| Build | CHARMM-GUI, repaired before simulation. CHARMM-GUI had translated the protein 7.2 Å too deep, so it was moved +7.18 Å along z onto the OPM membrane depth. Molecules overlapping the moved protein were removed whole: 25 lipids, 661 waters, 7 K+ and 11 Cl-; a further 11 bulk K+ were removed to restore neutrality. The three barrels had no water across the membrane, so 793 TIP3 waters taken from the system's bulk were placed inside them. The system was then minimised (5126 steepest-descent steps) and equilibrated with the CHARMM-GUI stages. |
| Engine | GROMACS 2021.5 for 0-100 ns; GROMACS 2026.1 for 100-200 ns, continued from the 100 ns `.gro` with `continuation = yes` and `gen-vel = no`, so velocities are at `.gro` precision and the thermostat and barostat integrals restart. Time step 2 fs, LINCS on bonds to hydrogen. PME with rcoulomb 1.2 nm and Fourier spacing 0.12 nm. Lennard-Jones cut-off 1.2 nm with a force switch from 1.0 nm. |
| Ensemble | 310 K with v-rescale (τ 1 ps; groups SOLU, MEMB, SOLV). 1 bar with C-rescale, semi-isotropic (τ 5 ps, compressibility 4.5e-5 bar⁻¹). |
| Field | `electric-field-z = 0.050 0 0 0`, production only: E0 = 0.05 V/nm with ω = 0 and σ = 0, so the field is **static**. V = E0 × ⟨Lz⟩ = **561 mV**, with Lz = 112.29 ± 0.92 Å (SD) over 200 ns and a drift of 0.2 Å per 100 ns. The field drives K+ towards +z (extracellular) and Cl- towards −z. |
| Centre-of-mass removal | `comm-mode linear`, `nstcomm 100`, `comm-grps SOLU_MEMB SOLV`. The solvent (water and ions) has its centre-of-mass velocity removed separately from the protein and membrane. This damps any net solvent flow through the pores (see Water). It does not affect the ion current. |
| Length and output | 2 × 50 000 000 steps = 200 ns. The two parts are joined into 2001 frames from 0 to 200 ns, saved every 100 ps (the second part's times restart at 0; its first frame equals the first part's last within 0.01 Å and is kept once). All atoms are written; xtc precision is 0.001 nm. |

These settings were checked against the run input (`step7_1.tpr`, read with `gmx dump`), the processed `.mdp`, the
topology, the run logs and the trajectory. The second part's `mdout.mdp` equals the first part's processed `.mdp` apart
from GROMACS 2026 defaults; its `.tpr` (tpx 138) cannot be read by MDAnalysis 2.x, so the first part's `.tpr`, with
the same topology, is used for both.

## Definitions analysed

The definitions were fixed before the trajectory was analysed.

- **Pores.** One cylinder per monomer, along z. Its centre is the centroid of the constriction arginines (Arg42,
  Arg82 and Arg132, atom CZ; `resindex 41 81 131`, plus 340 per chain) in the first frame. The cylinder extends
  ±15 Å along the axis and has a lateral radius of 8 Å.
- **Fit and counting rule.** Frames are fitted on the protein CA atoms. Crossings are counted with the compartment
  convention.
- **Sensitivity variants.** Radius 6 and 10 Å; extent ±5 Å; no lateral radius; the strict convention; the extent
  between the lipid phosphate planes (`--membrane-planes "name P"`).
- **The pre-registered definition held over 200 ns.**
  - With no lateral radius, 1067 K+ and 800 Cl- cross each pore's slab; the pores account for 1058 K+ and 801 Cl-
    (the three slabs lie up to 0.7 Å apart, so one transit can complete one slab and not another).
  - The radius-10 Å and extent-±5 Å variants change no per-pore count by more than 5.
  - The radius-6 Å variant loses 23-36 K+ per pore.
  - The phosphate-plane extent stays within 3 events of the base definition.
  - The strict convention keeps only 67 K+ and 435 Cl-, because ions leave the 8 Å cylinder in OmpF's wide
    vestibules. It is reported as a sensitivity only.
- **Cross-check.** CREVICE's event lists are identical to the independent reference counter's for every definition,
  convention and species. Both report 0 unresolved crossings.
- **No frame is excluded** from any result, including the frames where CREVICE and HOLE disagree (see Pore geometry).

## Results

**Ionic current and conductance**

| Estimator | Current (pA) | Trimer G (nS) |
| --- | --- | --- |
| Displacement current over every ion in the box (Aksimentiev and Schulten 2005), 0-200 ns | 1512 (1420-1604) | **2.69 (2.53-2.86)** |
| Counted crossings, sum over the three pores, 0-200 ns | 1489 (1335-1644) | 2.65 (2.38-2.93) |
| Per pore (counted) | 418 / 582 / 489 | 0.74 (0.44-1.05) / 1.04 (0.86-1.21) / 0.87 (0.77-0.97) |
| First 100 ns (displacement / counted) | 1703 / 1671 | 3.03 / 2.98 |
| Second 100 ns (displacement / counted) | 1321 / 1308 | 2.35 / 2.33 |

- Every complete crossing went in the field direction: 1058 K+ moved up and 801 Cl- moved down, with none reversed.
- **The two halves differ by 23 %**, far beyond either half's within-run interval. The run is not stationary on a
  100 ns scale, so the 0-200 ns intervals, which assume stationarity, understate the uncertainty.
- The collective-coordinate current, 1285 pA, is lower than the counted current. It sums the displacement of ions
  inside the 8 Å cylinders only, so it depends on the lateral radius.

**Selectivity**

Selectivity here is the K+/Cl- current ratio at 561 mV. It is not a permeability ratio from a reversal potential.

- The trimer ratio is 1.32 (1.08-1.56) from the counts and 1.28 (1.06-1.49) from the displacement current: 1.61 in
  the first 100 ns, 1.03 in the second.
- The per-pore ratios are 2.24 (1.58-2.90), 1.83 (1.35-2.30) and **0.59 (0.29-0.90)**.
- The three pores differ by more than counting noise allows: for K+, χ² = 84 on 2 degrees of freedom, and p = 0.01
  after allowing for over-dispersion.

**Changes of state at the constriction**

- **Pore 3**, from 33 ns to the end of the run: Glu117's carboxylate bridges Lys16 instead of Arg42 (E117-K16 bridge
  in 89-100 % of each 10 ns window from 40 ns on). Its K+ rate fell from 2.21 per ns before the switch to 0.9 per ns
  over 40-200 ns, and its Cl- rate rose from 1.21 to 2.1 per ns.
- **Pore 1**, from about 110-130 ns: the same switch. Its E117-R42 bridge holds in 83-100 % of each 10 ns window up to
  110 ns, 39 % in 120-130 ns and 0-8 % after 130 ns; E117-K16 then holds in 73-100 %. Its K+ crossings fell from
  19-37 per 10 ns to 1-9 between 100 and 180 ns (the fall began in 100-110 ns, before the switch was complete) and
  were 18-19 per 10 ns in 180-200 ns; its Cl- crossings fell too (8-16 to 1-9 per 10 ns).
- **Pore 2** kept its E117-R42 bridge until the last 10 ns (3 % in 190-200 ns, with no K16 bridge).
- Over 5 ns windows across the three pores, the E117-K16 bridge fraction is associated with fewer K+ crossings
  (Spearman ρ = −0.62) and more Cl- crossings (ρ = 0.25).
- Two pores of one trajectory now show the same change of state. It is still a descriptive association, not a
  demonstrated mechanism, and it is the main reason the two halves differ.

**Water**

- Net flux towards +z through pores 1 to 3 is 3.1 (0.8-5.5), 3.4 (1.0-5.8) and −3.8 (−7.1 to −0.5) molecules per ns.
  The trimer total is 2.7 (−1.8 to 7.2) per ns. Each pore's flux follows its majority ion.
- Because of the centre-of-mass removal groups, this is not a measurement of free electro-osmotic flow.

**Pore geometry over time**

- **Method.** CREVICE profiled each pore in every frame with the per-plane probe, with its widest-route centres. The
  axis is z and the origin is the frame-0 constriction centroid (Arg42/82/132 CZ, Asp113 CG, Glu117 CD). Settings:
  HOLE `simple.rad` radii, 0.25 Å sections, 81 samples.
- **Resolution.** 5464 of 6003 profiles resolved (pores 1-3: 1668, 1914, 1882 of 2001). The rest are recorded as
  unresolved, never as zero: 269 where the fixed origin has no clearance, 211 with a capped end, 59 not enclosed.
  Most are pore 1 after 90 ns, so its narrowest frames are under-represented in the means below.
- **Minimum radius.** 2.38 ± 0.49, 2.74 ± 0.26 and 2.47 ± 0.34 Å for pores 1 to 3; the crystal structure gives 3.56 Å.
  Over the first 100 ns 2.71, 2.71 and 2.64 Å; over the second 1.94, 2.77 and 2.30 Å.
- **Pore 1 narrowed first.** Its mean minimum radius per 10 ns was 2.5-2.9 Å up to 90 ns, 2.18 Å in 90-100 ns, 1.5-2.0 Å
  from 100 to 180 ns and 2.2-2.3 Å in 180-200 ns: the narrowing came before its current fell and before its salt
  bridge switched, and eased when its current recovered.
- **Geometric conductance.** At the experimental resistivity of 1 M KCl (0.0936 Ω m), the estimate is 1.46, 1.66 and
  1.58 nS per pore, 1.6-2.0 times the simulated per-pore conductance.
- **Geometry and current.** Over 5 ns windows the geometric conductance follows the simulated one with ρ = 0.56
  (n = 120). Within pores ρ = 0.88, 0.32 and −0.10: pore 1's change is visible in
  its profile, pore 3's at 33 ns is not (its minimum radius was 2.64 Å before and 2.67 Å after).
- **Agreement with HOLE** (every 20th frame, 303 pore profiles; 276 resolved by CREVICE's dense pass). The median
  minimum-radius difference is 0.003 Å and the median profile RMSD 0.067 Å. HOLE's in-protein minimum lies inside
  CREVICE's range in all 276. Every minimum that differs by more than 0.05 Å (68) was checked against the widest
  through-route on an exact 0.2 Å clearance grid (rule fixed before the run): both tools are on the route within the
  grid margin in 41; HOLE has left it in 25 (16 a one-point dip at its starting plane, 9 a side or narrower parallel
  branch, down to 0.55 Å); both have left it in 2 (pore 1, frames 1180 and 1680), where the wider route passes 8-10 Å
  from the pore axis, outside CREVICE's default 6 Å search square.

**Model bulk conductivity (indicator only)**

- From lateral ion diffusion in the bulk slab, D(K+) = 3.06e-5 cm²/s and D(Cl-) = 3.02e-5 cm²/s (the two halves agree to 1 %).
- The ideal Nernst-Einstein conductivity is 23.8 S/m. That is 2.2 times the 10.68 S/m of 1 M KCl at 20 °C.
- Ideal Nernst-Einstein ignores ion correlations, so it overestimates the conductivity of a 1 M electrolyte.

## Convergence

**Trimer conductance at 561 mV: not converged across conformational states.**

- At 200 ns the within-run relative half-width is 6 % (displacement current) or 10 % (counts), but the first and
  second 100 ns give 3.03 and 2.35 nS: a change of state in one pore moves the trimer value by more than that.
- Automated equilibration detection finds t0 = 0 for the whole-box displacement series and for pores 2 and 3, and
  t0 = 90 ns for pore 1's net charge (its change of state).
- The displacement series is uncorrelated from frame to frame (statistical inefficiency 1). Per-pore K+ counts are
  over-dispersed and their dispersion index grows with block length (to 15 at 20 ns in pore 1), the signature of slow
  changes of state.
- The pores differ by a coefficient of variation of 0.17. Taking that spread for independent replicas, about 7
  independent 200 ns trimer replicas would give ±10 % and about 4 would give ±20 %.

**Selectivity: not converged.** The within-run interval is ±18 %, and the halves give 1.61 and 1.03. The per-pore
ratio varies with a coefficient of variation of 0.55, so ±10 % needs about 120 independent pore samples and ±20 %
about 32.

**Water flux: not converged**, and affected by the centre-of-mass setting.

**Box and voltage: converged.** The voltage's standard error is 0.1 %.

## Comparison with experiment

The experimental value is **4.1 ± 0.4 nS per trimer**, measured in 1 M KCl at pH 7 in a planar bilayer at ±50 mV
(Wang, Prajapati, Kleinekathöfer and Winterhalter, Chem. Sci. 11, 10344 (2020), doi 10.1039/d0sc03486j). The
simulated value over 200 ns is 2.69 (2.53-2.86) nS, a ratio of **0.66** (first 100 ns: 3.03 nS, 0.74; second:
2.35 nS, 0.57). The within-run interval does not cover the change of state, so the ratio is uncertain by more than
it suggests.

The conditions differ in several ways:

- 561 mV against 50 mV. There is only one voltage, so the linearity of the I-V curve is not tested.
- 310 K, against the experiment's recording temperature.
- DOPE:DOPG, which is charged, against the experimental bilayer.
- Bulk ion concentrations of 1.06-1.11 M.
- The TIP3P/CHARMM36 electrolyte (CHARMM36 ion parameters).
- One 200 ns run, in which two of the three pores changed constriction state.

## Recommended next simulations

1. **Replicas before length.** The second 100 ns added a second change of state rather than narrowing the estimate.
   Run 4-6 independent trimer replicas of 200 ns at 0.05 V/nm with new velocities, ideally starting from different
   Glu117 states, and record that state (E117-R42 and E117-K16 distances per frame).
2. **An I-V curve.** Run E0 = ±0.025, ±0.05 and ±0.1 V/nm (±280, ±561 and ±1122 mV), including the reversed field.
   With the over-dispersion seen over 200 ns (dispersion index 3.9 for the trimer net charge), ±10 % on the trimer
   conductance needs about 320, 160 and 80 ns per voltage, and ±20 % about 80, 40 and 20 ns, if the channel is ohmic.
   Near the experimental 50 mV, ±10 % would take about 1.8 µs. Selectivity needs the replica sampling of item 1.
3. **Saving interval.** 100 ps is enough for counting ions: no ion crossed in 100 ps or less (one in 200 ps); the
   median crossing took 1.3 ns, and there were no unresolved crossings. Saving only protein, ions and water would
   reduce storage. Save every 10-20 ps only if the collective water permeability (pf) is wanted.
4. **Centre-of-mass removal.** For water or electro-osmotic flux, reconsider the `comm-grps` grouping.
5. **Continuations.** Continue a run with `gmx mdrun -cpi` and the same GROMACS version where possible; this run's
   second part restarted from the `.gro` with a newer version.

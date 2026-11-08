# OmpF trimer in 1 M KCl under an applied field

The definition key is `ompf_trimer_efield` in `../inputs/md_systems.json`; the driver is `../md_channel_analysis.py`.

This is one 100 ns simulation. Comparing it with an experiment means comparing one simulation with one experiment. It is
not a validation of CREVICE or of the force field. Every interval below is a 95 % interval.

## System

| | |
| --- | --- |
| Protein | OmpF (PDB 2OMF), biological trimer. Chains PROA, P2PROA and P3PROA hold residues 1-340 each, with NTER/CTER termini. Each chain has one His, set as HSD. Each chain carries a net charge of -12 e. |
| Orientation | OPM; tilt 0.01°. The N and C termini (periplasmic side) lie below the bilayer midplane, so **+z is extracellular**. |
| Membrane | 354 DOPE and 117 DOPG (3:1). The upper leaflet has 180 DOPE and 58 DOPG; the lower leaflet has 174 DOPE and 59 DOPG. |
| Solution | 47 179 TIP3 waters (CHARMM TIP3P), 1001 K+ and 848 Cl-. That is 1.0 M from the Cl- count. In the bulk slab (more than 25 Å from the bilayer midplane), the run averages 1.11 M K+ and 1.06 M Cl-. The excess K+ neutralises DOPG and the protein. The system has 219 592 atoms and zero net charge. |
| Force field | CHARMM36m protein with CHARMM36 lipids, ions and TIP3P water, built with CHARMM-GUI (topology from its FF-Converter); the release is as stated by the system's author. |
| Build | CHARMM-GUI, repaired before simulation. CHARMM-GUI had translated the protein 7.2 Å too deep, so it was moved +7.18 Å along z onto the OPM membrane depth. Molecules overlapping the moved protein were removed whole: 25 lipids, 661 waters, 7 K+ and 11 Cl-; a further 11 bulk K+ were removed to restore neutrality. The three barrels had no water across the membrane, so 793 TIP3 waters taken from the system's bulk were placed inside them. The system was then minimised (5126 steepest-descent steps) and equilibrated with the CHARMM-GUI stages. |
| Engine | GROMACS 2021.5. Time step 2 fs, LINCS on bonds to hydrogen. PME with rcoulomb 1.2 nm and Fourier spacing 0.12 nm. Lennard-Jones cut-off 1.2 nm with a force switch from 1.0 nm. |
| Ensemble | 310 K with v-rescale (τ 1 ps; groups SOLU, MEMB, SOLV). 1 bar with C-rescale, semi-isotropic (τ 5 ps, compressibility 4.5e-5 bar⁻¹). |
| Field | `electric-field-z = 0.050 0 0 0`, production only: E0 = 0.05 V/nm with ω = 0 and σ = 0, so the field is **static**. V = E0 × ⟨Lz⟩ = **561 mV**, with Lz = 112.15 ± 0.90 Å (SD) and no drift beyond its fluctuation. The field drives K+ towards +z (extracellular) and Cl- towards −z. |
| Centre-of-mass removal | `comm-mode linear`, `nstcomm 100`, `comm-grps SOLU_MEMB SOLV`. The solvent (water and ions) has its centre-of-mass velocity removed separately from the protein and membrane. This damps any net solvent flow through the pores (see Water). It does not affect the ion current. |
| Length and output | 50 000 000 steps = 100 ns. 1001 frames from 0 to 100 ns, saved every 100 ps. All atoms are written; xtc precision is 0.001 nm. |

These settings were checked against the run input (`step7_1.tpr`, read with `gmx dump`), the processed `.mdp`, the
topology, the run log and the trajectory.

## Definitions analysed

The definitions were fixed before the trajectory was analysed.

- **Pores.** One cylinder per monomer, along z. Its centre is the centroid of the constriction arginines (Arg42,
  Arg82 and Arg132, atom CZ; `resindex 41 81 131`, plus 340 per chain) in the first frame. The cylinder extends
  ±15 Å along the axis and has a lateral radius of 8 Å.
- **Fit and counting rule.** Frames are fitted on the protein CA atoms. Crossings are counted with the compartment
  convention.
- **Sensitivity variants.** Radius 6 and 10 Å; extent ±5 Å; no lateral radius; the strict convention. One extra
  variant, analysed afterwards, takes the extent between the lipid phosphate planes (`--membrane-planes "name P"`).
- **The pre-registered definition held.**
  - Of the transmembrane transits found with no lateral radius, the pores account for 643 of 645 K+ and all 400 Cl-.
  - The radius-10 Å and extent-±5 Å variants change no per-pore count by more than 3.
  - The radius-6 Å variant loses 13-24 K+ per pore.
  - The phosphate-plane extent stays within 2 events of the base definition.
  - The strict convention keeps only 34 K+ and 220 Cl-, because ions leave the 8 Å cylinder in OmpF's wide
    vestibules. It is reported as a sensitivity only.
- **Cross-check.** CREVICE's event lists are identical to the independent reference counter's for every definition,
  convention and species. Both report 0 unresolved crossings.

## Results

**Ionic current and conductance**

| Estimator | Current (pA) | Trimer G (nS) |
| --- | --- | --- |
| Displacement current over every ion in the box (Aksimentiev and Schulten 2005) | 1703 (1604-1802) | **3.04 (2.86-3.21)** |
| Counted crossings, sum over the three pores | 1671 (1552-1790) | 2.98 (2.77-3.19) |
| Per pore (counted) | 638 / 578 / 455 | 1.14 (1.01-1.27) / 1.03 (0.83-1.24) / 0.81 (0.65-0.97) |
| The three pores treated as replicas (trimer = 3 × mean, t(2)) | 1671 (977-2365) | 2.98 (1.74-4.22) |

- Every complete crossing went in the field direction: 643 K+ moved up and 400 Cl- moved down, with none reversed.
- The collective-coordinate current, 1462 pA, is lower than the counted current. It sums the displacement of ions
  inside the 8 Å cylinders only, so it depends on the lateral radius.

**Selectivity**

Selectivity here is the K+/Cl- current ratio at 561 mV. It is not a permeability ratio from a reversal potential.

- The trimer ratio is 1.61 (1.39-1.83) from the counts and 1.51 (1.28-1.75) from the displacement current.
- The per-pore ratios are 2.40 (1.75-3.05), 2.03 (1.51-2.56) and **0.73 (0.23-1.24)**.
- The three pores differ by more than counting noise allows. For K+, the counts give χ² = 66 on 2 degrees of freedom,
  and p = 1e-4 after allowing for over-dispersion.

**Pore 3**

- From 33 ns until 99.5 ns, Glu117's carboxylate in pore 3 bridged Lys16 instead of Arg42. The E117-R42 bridge is
  present in 36 % of frames in pore 3, compared with 92-94 % in pores 1 and 2.
- Over that interval, pore 3's K+ rate fell from 2.21 to 0.71 per ns, and its Cl- rate rose from 1.21 to 1.83 per ns.
- Over 5 ns windows across all three pores, the E117-K16 bridge fraction is associated with fewer K+ crossings
  (Spearman ρ = −0.56) and more Cl- crossings (ρ = 0.41).
- This is one change of state in one trajectory. It is a descriptive association, not a demonstrated mechanism.

**Water**

- Net flux towards +z through pores 1 to 3 is 5.2 (1.7-8.8), 4.9 (2.7-7.0) and −2.4 (−6.9 to 2.1) molecules per ns.
  The trimer total is 7.7 (3.2-12.2) per ns. Each pore's flux follows its majority ion.
- Because of the centre-of-mass removal groups, this is not a measurement of free electro-osmotic flow.

**Pore geometry over time**

- **Method.** CREVICE profiled each pore in every frame with the per-plane probe. The axis is z and the origin is the
  constriction centroid (Arg42/82/132 CZ, Asp113 CG, Glu117 CD). Settings: HOLE `simple.rad` radii and 0.25 Å
  sections.
- **Resolution.** 2873 of 3003 profiles resolved. The rest are recorded as unresolved, never as zero.
- **Minimum radius.** 2.70 ± 0.28, 2.70 ± 0.25 and 2.64 ± 0.23 Å for pores 1 to 3; the crystal structure gives 3.56 Å.
- **Geometric conductance.** At the experimental resistivity of 1 M KCl (0.0936 Ω m), the estimate is about 1.67 nS
  per pore, which is 1.5-2.0 times the simulated per-pore conductance.
- **Geometry and current.** Geometry follows the current only weakly: over 5 ns windows, ρ = 0.32.
- **Pore 3.** Its geometric conductance did not change when its current fell.
- **Agreement with HOLE.** On every 20th frame, HOLE and CREVICE give the same minimum radius to within 0.002 Å
  (median), with a median profile RMSD of 0.049 Å. In 31 of those frames, HOLE's minimum is a one-point dip at its
  starting plane.

**Model bulk conductivity (indicator only)**

- From lateral ion diffusion in the bulk slab, D(K+) = 3.05e-5 cm²/s and D(Cl-) = 3.02e-5 cm²/s.
- The ideal Nernst-Einstein conductivity is 23.8 S/m. That is 2.2 times the 10.68 S/m of 1 M KCl at 20 °C.
- Ideal Nernst-Einstein ignores ion correlations, so it overestimates the conductivity of a 1 M electrolyte.

## Convergence

**Trimer conductance at 561 mV: converged within the run, not across conformational states.**

- At 100 ns the relative half-width of the interval is 6 % (displacement current) or 7 % (counts). It falls below
  10 % after about 34-50 ns.
- No equilibration period needs discarding. Automated detection finds t0 = 0, and discarding up to 50 ns moves G only
  between 2.99 and 3.08 nS.
- The displacement series is uncorrelated from frame to frame (statistical inefficiency 1).
- The pores differ by a coefficient of variation of 0.17, and pore 3's change of state lasted more than 60 ns. Taking
  the observed between-pore spread for independent replicas, about 7 independent 100 ns trimer replicas would give
  ±10 % and about 4 would give ±20 %.

**Selectivity: not converged.** One change of state in one pore dominates it.

- The within-run interval is ±14 %, which would need about 190 ns to reach ±10 %.
- The per-pore ratio varies with a coefficient of variation of 0.51, so ±10 % needs about 100 independent pore samples
  and ±20 % about 28.

**Water flux: not converged**, and affected by the centre-of-mass setting.

**Box and voltage: converged.** The voltage's standard error is 0.1 %.

## Comparison with experiment

The experimental value is **4.1 ± 0.4 nS per trimer**, measured in 1 M KCl at pH 7 in a planar bilayer at ±50 mV
(Wang, Prajapati, Kleinekathöfer and Winterhalter, Chem. Sci. 11, 10344 (2020), doi 10.1039/d0sc03486j). The
simulated value is 3.04 (2.86-3.21) nS, a ratio of **0.74**.

The conditions differ in several ways:

- 561 mV against 50 mV. There is only one voltage, so the linearity of the I-V curve is not tested.
- 310 K, against the experiment's recording temperature.
- DOPE:DOPG, which is charged, against the experimental bilayer.
- Bulk ion concentrations of 1.06-1.11 M.
- The TIP3P/CHARMM36 electrolyte (CHARMM36 ion parameters).
- One 100 ns run.

## Recommended next simulations

1. **Replicas before length.** Run 4-6 independent trimer replicas of 200-300 ns at 0.05 V/nm. Use new velocities and
   ideally different starting Glu117 states, and record that state.
2. **An I-V curve.** Run E0 = ±0.025, ±0.05 and ±0.1 V/nm (±280, ±561 and ±1122 mV), including the reversed field.
   - For ±10 % on the trimer conductance, this needs about 76, 38 and 19 ns per voltage respectively.
   - For ±20 %, it needs about 19, 9 and 5 ns.
   - Near the experimental 50 mV (E0 = 0.0045 V/nm), ±10 % would take about 420 ns.
   - These estimates assume an ohmic channel and the dispersion observed here.
   - Selectivity needs the replica sampling of item 1 at every voltage.
3. **Saving interval.** 100 ps is enough for counting ions. No ion crossed in 100 ps or less; the median crossing took
   1.2 ns, and there were no unresolved crossings. Saving only protein, ions and water would reduce storage. Save every
   10-20 ps only if the collective water permeability (pf) is wanted.
4. **Centre-of-mass removal.** For water or electro-osmotic flux, reconsider the `comm-grps` grouping.

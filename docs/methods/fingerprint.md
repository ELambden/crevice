# Pore fingerprint definitions

This note gives the definitions behind
[`crevice profile --fingerprint`](../tools/fingerprint.md) and
{mod}`crevice.fingerprint`. All distances are in ångström, and atoms are van
der Waals spheres of the radius set in force
([Atomic radii](atomic-radii.md)).

## The path coordinate *s*

The fingerprint puts everything on one path coordinate *s* through a resolved
channel.

- **Lumen.** The lumen samples are the channel profile's own samples, at
  *s* = *t* − *t*<sub>first</sub>, where *t* is the profile's axial
  coordinate (or the cumulative path length for an explicit path). The lumen
  therefore runs from *s* = 0 to its length *L*.
- **Beyond each mouth** the path follows one centre line out to bulk solvent.
  At an open end this is the straight vestibule centre line used for the cast
  segments ({func}`crevice.cast_segments.axial_vestibule_leg`, sampled every
  0.25 Å). At a capped end (`--lateral-exits`) it is the widest lateral exit
  leg ({func}`crevice.channel_exits.lateral_exit_legs`, rank 1). These samples
  sit at *s* = −*arc* below the lower mouth and *s* = *L* + *arc* above the
  upper one, where *arc* is the path length along that centre line.
- **Other exit legs** are reported in the JSON (rank, length, bottleneck) but
  are not placed on *s*.

**Segments.** Lumen samples are labelled `lumen`; the samples beyond each mouth
take that end's ENTRY or EXIT label (`entry_1`, `exit_1`), following the same
rule as `report`'s cast segments
({func}`crevice.cast_segments.segment_roles`, `--entry-end`). The labels are
geometric and say nothing about the direction of transport.

**Radius.** For lumen samples this is the profile radius (free-sphere
clearance minus `--probe-radius`). Along a vestibule or exit leg it is the
clearance of the centre line minus the same probe radius: a free-sphere radius
along a 3D path, not the radius of a planar section.

## Lining shell

An atom lines a path sample when its van der Waals surface lies within the
lining shell *δ* (`--lining-shell`, default 1.4 Å, one water radius) of the
largest atom-free sphere at that sample:

|**x** − **c**| − *r*<sub>atom</sub> − *R* ≤ *δ*,

where **c** is the sample centre and *R* its free-sphere (raw) clearance.
Hydrogens, water and recognised monatomic ions are not counted. A residue
lines a sample when any of its atoms does, and its **anchor** *s*<sub>j</sub>
is the sample where its surface gap is smallest. Each residue is counted once,
at its anchor.

## Hydrophobicity

**Scale.** Each lining residue *j* takes the value *h*<sub>j</sub> from the
Wimley–White whole-residue octanol scale (Wimley, Creamer & White 1996; White
& Wimley 1999), with the sign reversed, *h* = −Δ*G*(water → octanol) in
kcal/mol, so that positive means hydrophobic. We use the octanol scale because
it is measured for whole residues (side chain plus peptide unit), covers both
charged and neutral forms, and describes moving into a low-dielectric
environment, which is the contrast that matters at a pore wall.

**Charge states.** Asp, Glu, Lys and Arg are taken as charged and His as
neutral (pH 7). Protonation-variant names (HSD, HSE, HIE, HID, HSP, HIP, ASH,
GLH, LYN, CYX and others) and common D-amino-acid and modified residues (DLE,
DVA, MSE and others) map to their parent and charge state
(`crevice.fingerprint.RESIDUE_PARENTS`). Ligands, lipids and caps have no
value and are counted as `other`.

**Gaussian smoothing.** The profile is a kernel-weighted average of the
residue values along *s*:

*H*(*s*) = Σ<sub>j</sub> *h*<sub>j</sub> *K*(*s* − *s*<sub>j</sub>) /
Σ<sub>j</sub> *K*(*s* − *s*<sub>j</sub>), with
*K*(*d*) = exp(−*d*² / 2*σ*²),

one weight per residue, and *σ* set by `--hydrophobicity-bandwidth` (default
3 Å, about one helical turn of rise). Where Σ *K* < 0.01, that is, no lining
residue within about 3*σ*, the value is left empty. CHAP uses the same kernel
form (Klesse et al. 2019).

**Chemistry.** Residues are also classed as hydrophobic, polar, charged or
aromatic by {func}`crevice.radii.residue_properties` of their parent residue,
and the counts are given per segment. A residue can belong to several classes
(Tyr is hydrophobic, polar and aromatic), so the class counts need not add up
to the number of residues.

## Water along the path (trajectories)

The path is fixed at the reference frame. Every frame is fitted to it with the
same rigid fit as the protein, the radius is measured again at the fixed path
samples, and water oxygens are imaged about the path centroid with the
periodic minimum image
({func}`crevice.water_membership.image_waters_around_anchor`). A water oxygen
belongs to sample *i* when *i* is its nearest path sample and the oxygen lies
inside that frame's free sphere at the sample (distance to the centre at most
the frame's raw clearance there).

With bin width Δ*s*<sub>i</sub> (half the distance to the neighbouring
samples), each sample reports:

- the mean count per frame;
- the linear density, mean count / Δ*s*<sub>i</sub> (molecules per Å);
- the number density, mean count / mean(π*R*²Δ*s*<sub>i</sub>) (Å⁻³; bulk
  water is about 0.0334 Å⁻³);
- the occupancy, the fraction of frames with at least one oxygen in the bin.

A free sphere is smaller than a wide or non-circular cross-section, so these
are lower bounds on the water that fills the pore there.

## Conductance estimate

The conductance follows the geometric estimate of HOLE (Smart et al. 1993,
1997), over the lumen only:

*G*<sub>macro</sub> = 1 / (*ρ F*), with
*F* = ∫ d*s* / π*r*(*s*)²

integrated with the trapezoid rule. *ρ* is the bulk resistivity of the
electrolyte at 25 °C, from the conductivities in the CRC Handbook of Chemistry
and Physics (150 mM values interpolated):

| `--electrolyte` | Conductivity (S/m) | *ρ* (Ω m) |
|---|---|---|
| `kcl-150mm` (default) | 1.88 | 0.532 |
| `nacl-150mm` | 1.56 | 0.641 |
| `kcl-1m` | 11.2 | 0.0894 |
| `nacl-1m` | 8.6 | 0.116 |

`--resistivity` sets any other value (HOLE itself uses 1/12 Ω m for 1 M KCl).
Two further values are reported:

- ***G*<sub>access</sub>** adds the access resistance *ρ*/(4*a*) of each
  mouth, where *a* is the radius at that end of the lumen (Hall 1975; Hille
  2001): *G*<sub>access</sub> = 1 / [*ρ* (*F* + 1/4*a*<sub>1</sub> +
  1/4*a*<sub>2</sub>)].
- ***G*<sub>HOLE-corrected</sub>** = *G*<sub>macro</sub> / 5.594, the
  constant empirical divisor of HOLE 2.3.1, which was calibrated on eight
  channels in 1 M KCl (Smart et al. 1997). It is shown for orientation only.

A zero radius anywhere in the lumen gives *G* = 0 (occluded). For a
trajectory, the conductance is computed in every frame and summarised.

**How it compares.** On gramicidin A (1GRM) with *ρ* = 1/12 Ω m, CREVICE gives
*G*<sub>macro</sub> ≈ 390 pS with the default radii and 318 pS with the HOLE
radius set; HOLE 2.3.1 on the same structure gives 280 pS (its integration
runs past the mouths to a 5 Å end radius). With the HOLE divisor, CREVICE gives
57–70 pS. The measured single-channel conductance of gramicidin A in 1 M salt
is about 10–30 pS, so the uncorrected value is 10–40 times too high and the
corrected value 2–7 times too high: the overestimate Smart et al. (1997)
describe for narrow pores.

## Limitations

- The conductance treats the pore as a tube of bulk electrolyte. It ignores
  ion–wall interactions, selectivity, dehydration and surface charge, and bulk
  resistivity does not hold where ions pass in single file. It is an
  order-of-magnitude geometric estimate, not a measurement or a prediction of
  function.
- Hydrophobicity is assigned per residue, so it cannot tell a backbone-lined
  lumen from a side-chain-lined one (gramicidin A is lined by carbonyl oxygens
  of residues the scale calls hydrophobic), and it knows protonation states
  only through residue names.
- Only one exit leg per end is placed on *s*.
- In trajectories the path and the lining residues are those of the reference
  frame; neither is resolved again in each frame.
- Water counts use free spheres, so they are lower bounds in wide or
  non-circular parts of the pore.

The software tests check these definitions on synthetic inputs with known
answers (an analytic cylinder, rings of Leu and Asp, placed water oxygens);
they are checks of the software, not validation on real channels.

## References

- Hall JE (1975) Access resistance of a small circular pore. *J Gen Physiol*
  66:531–532.
- Hille B (2001) *Ion Channels of Excitable Membranes*, 3rd ed., Sinauer,
  ch. 11.
- Klesse G, Rao S, Sansom MSP, Tucker SJ (2019) CHAP: a versatile tool for the
  structural and functional annotation of ion channel pores. *J Mol Biol*
  431:3353–3365.
- *CRC Handbook of Chemistry and Physics*, "Electrical Conductivity of Aqueous
  Solutions", CRC Press.
- Smart OS, Goodfellow JM, Wallace BA (1993) The pore dimensions of
  gramicidin A. *Biophys J* 65:2455–2460.
- Smart OS, Breed J, Smith GR, Sansom MSP (1997) A novel method for structure-based
  prediction of ion channel conductance properties. *Biophys J* 72:1109–1126.
- White SH, Wimley WC (1999) Membrane protein folding and stability: physical
  principles. *Annu Rev Biophys Biomol Struct* 28:319–365.
- Wimley WC, Creamer TP, White SH (1996) Solvation energies of amino acid side
  chains and backbone in a family of host–guest pentapeptides.
  *Biochemistry* 35:5109–5124.

# Permeation and conductance

`crevice trajectory --permeation` follows water molecules and ions through a
channel in an MD trajectory. It lists every complete crossing (which molecule,
when it entered and left, in which direction and how long it took), counts
crossings per direction, measures how many molecules are inside the channel in
each frame, and derives:

- for water, the osmotic permeability **pf** (collective-diffusion model) and
  the diffusive permeability **pd** (from the crossing rate);
- for ions, the **current** from the net charge carried across, and the
  **conductance** when the voltage across the membrane is known.

The Python interface is {mod}`crevice.permeation`.

## Quick start

With the channel taken from the profile that `crevice trajectory` measures in
each frame:

```bash
crevice trajectory production.xtc --topology system.tpr \
    --axis z --origin 52.1,48.0,80.5 \
    --permeation TIP3,SOD,CLA -o run/channel.json
```

With an explicit channel (faster, and the route when no profile resolves):

```bash
crevice trajectory production.xtc --topology system.tpr --permeation-only \
    --permeation TIP3,SOD,CLA --axis z --channel-extent 68,92 \
    --channel-radius 4 --origin 52.1,48.0,80.5 -o run/channel.json
```

Four pores of a tetramer, counted separately:

```bash
crevice trajectory production.xtc --topology system.tpr --permeation-only \
    --permeation SOL,NA,CL --axis z --channel-extent 68,92 --channel-radius 4 \
    --permeation-pore 40,40,80 --permeation-pore 60,40,80 \
    --permeation-pore 60,60,80 --permeation-pore 40,60,80 -o run/aqp.json
```

```python
from crevice.permeation import CylinderChannel, analyze_permeation, write_permeation_outputs

channel = CylinderChannel("pore", axis=(0, 0, 1), zmin=68, zmax=92,
                          point=(52.1, 48.0, 80.5), radius=4.0)
result = analyze_permeation("system.tpr", "production.xtc",
                            species=["TIP3", "SOD", "CLA"], channels=[channel])
write_permeation_outputs(result, "run", prefix="channel")
print(result.summary["channels"][0]["species"][0]["events_up"])
```

## Choosing what to follow

`--permeation` takes residue names, atom names or MDAnalysis selections:

- `--permeation TIP3,SOD,CLA` (CHARMM), `--permeation SOL,NA,CL` (GROMACS),
  `--permeation WAT,Na+,Cl-` (AMBER), `--permeation W,ION` (Martini).
- A value containing spaces is one selection; give it a label with
  `LABEL=selection`, for example `--permeation "glycerol=resname GOL and name C2"`.
- Names that match nothing are listed in the summary and skipped, so a
  generic list such as `W,ION,TIP3,SOL,SOD,NA,CLA,CL,POT,K` works on any system.

Each residue is one molecule, followed through one atom: the only atom of an
ion or bead, otherwise the water oxygen (`OH2`, `OW`, `O`, `OW1`, `OH`) or
bead `W`, otherwise the first selected atom of the residue. A molecule is
counted in one species only.

**Charges** come from a table of common names (`NA`, `SOD`, `Na+`, `K`, `POT`,
`K+` = +1; `CL`, `CLA`, `Cl-` = −1; `MG`, `CA`, `CAL` = +2, and similar), then
from the topology charges; `--charge NAME=q` overrides either. The summary
records where each charge came from.

## Defining the channel

All channel coordinates are in the frame of the first analysed frame: every
frame is first superimposed on it using the CA atoms of `--selection` (Martini
`BB` beads for coarse-grained input). Exclude flexible tails from
`--selection` if they move a lot; the summary reports the fit RMSD.

| Definition | Options | Inside the lumen when |
| --- | --- | --- |
| Measured profile (default) | the trajectory command's profile options | within the frame's profile extent and within the profile radius plus `--lumen-tolerance` (1 Å) of the profile centre line |
| Explicit cylinder | `--channel-extent ZMIN,ZMAX` along `--axis`, optional `--channel-radius`, `--origin` or `--permeation-pore` | `ZMIN ≤ x·axis ≤ ZMAX` and within the radius of the axis line |
| Membrane-referenced cylinder | `--membrane-planes SELECTION` (for example `"name P"`) along `--axis`, optional `--channel-radius`, `--origin` or `--permeation-pore` | between the frame's lower and upper membrane planes and within the radius of the axis line |
| Region map | `--permeation-region FILE` (`.dx` map or prepared region `.json`) | in a region voxel (dilated by `--lumen-tolerance`, default 0) |

The extent is measured as the projection `x·axis` of reference-frame
coordinates (Å). Without `--channel-radius` every molecule within the extent
counts as inside, which is only appropriate when nothing else (lipid, other
pores) lies within that slab. A frame without a resolved profile reuses the
latest earlier profile; such frames are counted in the summary.

**Following the membrane.** `--membrane-planes "name P"` replaces a fixed
extent by the bilayer itself: in every frame the extent runs from the mean
position of the lower headgroup atoms to that of the upper ones (the selected
atoms are split at their mean along `--axis`, after the protein fit, as nearest
images about each pore). An event is then a passage from one side of the
bilayer to the other, however the membrane moves or thickens relative to the
protein. This is the reference of membrane-based ion counters (for example the
phosphate planes of Gapsys and Kopec's ion permeation analyzer); the
compartment rule, periodic unwrapping and lateral attribution to each pore are
unchanged. The summary records the planes (`membrane_planes`: reference-frame
planes, mean and standard deviation of each, mean thickness). `n(t)` uses the
reference frame's planes as its fixed extent. On a 100 ns OmpF trimer run in
1 M KCl with an applied field (1001 frames saved every 100 ps), CREVICE and a
transcription of that analyzer's counter both found 643 K⁺ and 398 Cl⁻
crossings, 1040 of them the same ion in the same frame. The other two were
ions that moved through more than half the membrane in one saved interval: the
analyzer's jump rule counted one downward move of 24 Å as an upward crossing
and discarded one real crossing, while CREVICE follows the unwrapped
displacement.


## How crossings are counted

In each frame a molecule is **below** the channel (`x·axis < ZMIN`),
**above** it (`> ZMAX`), **inside** the lumen, or **out**: within the extent
but outside the lumen laterally, for example in the protein wall, the
membrane or another pore.

A **permeation event** is a complete transport from one side to the other
(the *compartment* rule, the default). Each molecule remembers the bulk side
it last occupied. A *passage* is the run of frames it spends within the axial
extent between two bulk frames; when the passage ends on the side opposite
to the one the molecule came from, that is an event: *up* (along `+axis`) or
*down*. The event is attributed to the channel when the molecule was inside
the lumen in at least one saved frame of the passage; otherwise it is listed
with status `outside_lumen` and not counted. This is the rule of GROMACS
computational electrophysiology (each ion remembers its compartment, and a
change of compartment is credited to the channel cylinder it last passed
through) and of a "complete transport of a water molecule through the
channel" (Zhu, Tajkhorshid and Schulten, 2004a).

- Entering and leaving on the same side is not an event, however long the
  molecule stays.
- Leaving the lumen sideways within the extent (*out*) does not forget the
  side: a molecule that steps out of the lumen and back in, or enters the
  extent beside the lumen and then moves into it, is counted when it reaches
  the other side.
- `outside_lumen` passages (listed in the events table and counted as
  `transits_outside_lumen`) are transits through the extent that never
  touched the lumen: through another pore, the membrane, or a lumen drawn
  too narrow. With several pores, each pore's count is independent; their
  cylinders must not overlap.
- Periodic boundaries are handled per molecule: positions are taken as the
  nearest image about the channel centre, and frame-to-frame displacements
  are unwrapped. A molecule that reaches the other side through the periodic
  bulk is not counted (its side simply changes), and one that crosses the
  channel several times is counted each time.
- A molecule that crosses the whole extent between two saved frames is an
  **unresolved crossing**: listed, not counted, and reported as a warning.
  Save frames more often if there are many.

`--permeation-convention strict` keeps the stricter rule of earlier
versions: the molecule must stay inside the lumen without interruption from
one side to the other, and a frame outside the lumen laterally forgets the
side. It undercounts as soon as the lumen radius or extent does not match
the pore exactly. In a synthetic test (Brownian particles crossing an
impenetrable membrane through a 4 Å pore, saved every 1 ps), a cylinder
1 Å narrower than the pore counted all 722 membrane transits under the
compartment rule and 1 under the strict rule. With an extent reaching 2 Å into
the bulk, the counts were 0.81 and 0.15 of the transits.

**Choosing the lumen.** The axial extent decides what a transit is; make it
end where the pore meets the bulk on each side. The lateral radius (or region)
only attributes transits to a pore, so make it contain the whole pore lumen.
A radius 1–2 Å wider than the pore does not change the counts. The
radius must be smaller than the distance to the next pore. With frames saved
far apart, a short passage may have no saved frame in a narrow lumen, so it
is listed as `outside_lumen` instead of counted. In the same test saved every
50 ps, a radius 1 Å narrower than the pore lost 15 % of the transits this way;
the pore radius lost none. An extent reaching 4 Å into the bulk made it worse
(292 passages outside a pore-sized lumen), because the membrane crossing
often fell between two saved frames. Check `transits_outside_lumen` and the
warning. For a single pore in an otherwise impermeable membrane, the
`outside_lumen` passages are transits too, so events plus
`transits_outside_lumen` is then the better count.

Entry and exit frames, transit times (exit time minus entry time) and dwell
times are resolved to the saved-frame interval. The entry frame is the first
frame of the passage within the extent (strict rule: the first frame in the
lumen). Lumen visits cut by the start or end of the trajectory are marked as
censored and left out of the dwell statistics, and a molecule that starts
within the extent has no known origin, so its first passage is not counted.

## Water permeability

Two permeabilities characterise a water channel (Zhu, Tajkhorshid and
Schulten, 2004a,b):

- **Diffusive permeability** `pd = v_w q0`, where `q0` is the rate of
  permeation events in one direction (at equilibrium,
  `q0 = (up + down) / 2T`) and `v_w = 18 cm³ mol⁻¹ / N_A` is the volume of one
  water molecule. Its interval is the exact (Garwood) Poisson interval of the
  event count, widened when events are correlated (see
  [Intervals for correlated events](#intervals-for-correlated-events)).
- **Osmotic permeability** `pf = v_w D_n`, from the collective coordinate
  `n(t)`: in each frame interval the axial displacements of the molecules
  inside the channel are summed and divided by the channel length `L` (only
  the part of a displacement inside the extent counts; a step over the whole
  extent counts only if its straight path crosses the middle of the extent
  inside the lumen), and `n(t)` is the running sum. At equilibrium `n`
  diffuses, `⟨Δn²⟩ = 2 D_n τ`. CREVICE fits
  the mean-square displacement over all time origins with `a + 2 D_n τ`
  (lags of 1 frame to a tenth of a segment by default) and reports an
  interval from the spread of the same fit on `--msd-blocks` (5) contiguous
  segments. Because `D_n` cannot be negative, the lower bound is clipped at
  zero (the raw value is kept).

`pf` and `pd` are per channel, in cm³ s⁻¹. Their ratio `pf/pd` is also
reported. Both assume a stationary equilibrium trajectory, and their
intervals describe sampling only.

## Ion current and conductance

For each charged species and for all of them together, the current is
`I = e Σ q (up − down) / T`, positive when net positive charge moves along
`+axis`. Its interval (`interval_pA`) is the Poisson interval (normal
approximation to the difference of two Poisson counts; with no events, the
exact upper bound on the number of crossings), widened when the net charge
varies more between time blocks than Poisson counts would (see below). The
plain Poisson interval and a block-bootstrap interval of the net charge per
block of frames are also reported for comparison; the block bootstrap
resamples only 10 blocks and was too narrow in synthetic runs with a constant
drift (0.85–0.93 coverage).

The **conductance** `G = I / V` (pS) is reported only when the voltage is
known:

- `--voltage MV`, positive when it drives cations along `+axis`;
- an applied field: `--efield 0.05` (V/nm along z, or `EX,EY,EZ`) or
  `--efield-from run.mdp|run.tpr` (reading a `.tpr` needs GROMACS `gmx` on the
  `PATH`). The voltage is the field times the mean box length along it, with
  the sign of the channel axis; time-dependent fields are refused, and the
  conductance is withheld if the channel axis is more than 45° from the
  field.

At equilibrium (no voltage), CREVICE reports the events and the net flux, and
states that the conductance is undefined.

## Intervals for correlated events

Permeation events are often not independent. Water in a single-file pore
moves as a column, so its crossings come in correlated runs; a gating channel
conducts in bursts. A Poisson interval is then too narrow. CREVICE splits the
analysed frames into 10 contiguous time blocks (`--event-blocks`), counts the
events in each, and computes the **dispersion index**: the variance of the
block counts divided by their mean (1 for a Poisson process). The rate, q0
and pd intervals are quasi-Poisson intervals: the Garwood interval of
`k / phi` scaled by `phi`, with `phi` the dispersion index floored at 1, so
independent events keep the exact Garwood interval. The current uses the block
variance of the net charge (with a Student t factor for the 10 blocks) when
it exceeds the Poisson variance.

The dispersion indices are reported (`event_dispersion` per channel and
species, `water.event_dispersion_index`, `ion_current_total.dispersion_index`),
and the plain Garwood intervals are kept (`*_garwood_interval_*`). On
synthetic trajectories with a known rate, the 95 % pd interval covered the
true value in 0.98 of runs for independent walkers (the Garwood floor), and
in 0.95 instead of 0.86 for single-file water. For a channel gated on and off
with a mean dwell of a hundredth of the run, coverage went from 0.74 to 0.94
for pd and from 0.53 to 0.95 for the current. The widening assumes that
events are correlated over less than a block (a tenth of the run). With a
gate whose mean dwell is half a block, pd coverage was 0.91 (Garwood: 0.39),
so slower processes are still somewhat underestimated.

## Where molecules go: position traces and the axial profile

**Position traces.** Every molecule with a counted event is followed along the
channel axis: through every frame it spends within the extent, from the bulk
frame before each passage to the bulk frame after it. The traces are written
to `STEM_permeation_positions.csv` (one row per molecule and frame:
`axial_A`, `state` = `below`, `inside`, `out` or `above`, and the direction of
its first event) and drawn in `STEM_permeation_positions.png`: one panel per
species, axial position against time, orange for molecules whose first event
went up and teal for down, with the extent as dashed lines. They show how
crossings happen: a steady drift under an applied field, hops between sites,
or long waits at a constriction. Each trace is continuous in time (unwrapped
frame to frame) and placed in the channel frame at a frame inside the lumen.
At most `--max-traces` molecules per channel and species are kept (200 by
default, evenly spaced in order of their first event; `0` keeps all); the
summary records how many had events and how many were traced
(`position_traces`).

**Axial density and free energy.** The molecules inside the lumen are counted
in 0.5 Å bins along the axis (`STEM_permeation_profile.csv`: `molecules_per_frame`,
`density_per_A`). For a run in equilibrium, CREVICE also gives the free energy
along the axis by Boltzmann inversion, `G(s) = −ln(n(s)/n_ref)` in kT, with
`n_ref` the mean density of the bins within 1 Å of either end of the extent,
so the mouths are the zero. Its interval comes from the same time blocks as
the event intervals (`--event-blocks`, Student t over the blocks); a bin that
is empty in some block has no interval. With `--temperature K` the profile is
also given in kcal/mol. `STEM_permeation_profile.png` has one density panel
per species (water and ion densities differ by orders of magnitude) and, when
defined, the free energies together in a last panel.

- The free energy is withheld when a voltage or an applied field is given
  (`--voltage`, `--efield`, `--efield-from`): a driven run has a steady-state
  distribution, not an equilibrium one. The density is still reported, and
  shows where ions accumulate. A warning also says when the net flux of a
  species is larger than chance, which suggests a drive the analysis was not
  told about.
- It is the free energy of the distribution inside the lumen as defined: it
  includes the effect of the lumen's cross-section and of the definition's
  radius or region (a cylinder that misses part of the lumen at one end lowers
  the density there), and it is not corrected for the accessible area.
- Narrow sites need many visits: check the interval band before reading a well
  or barrier.

## Coarse-grained (Martini) input

- Martini water beads `W` (also `WF`, `PW`) stand for four water molecules
  (`--water-mapping` changes it). Event counts, `q0` and `pd` are multiplied
  by the mapping, and `n` is expressed in water molecules (so `D_n` is
  multiplied by the mapping squared) before `pf = v_w D_n`. The summary states
  the scaling.
- Martini 2 `ION` residues are split by atom name (`ION:NA+`, `ION:CL-`), and
  the charge is read from the atom name.
- Coarse-grained input is flagged (`coarse_grained: true`). Coarse-grained
  time is not converted, so rates are in simulation time.
- Beads have no atomic radii, so the profile cannot define the channel: use
  `--permeation-only` with `--channel-extent` (or a region map).

## Outputs

Next to `-o` (with its stem as prefix):

| File | Contents |
| --- | --- |
| `STEM_permeation.json` | summary: definitions, units, species, channel definitions, counts, rates, occupancy (mean with interval, and `frame_fraction_by_count`: the fraction of frames with 0, 1, 2, ... molecules in the lumen), dwell and transit times, pf/pd and current with their intervals, voltage, warnings and limitations |
| `STEM_permeation_events.csv` | one row per event, unresolved crossing or passage outside the lumen (`status` = `complete`, `unresolved_single_interval`, `outside_lumen`): channel, species, residue, atom index, direction, entry/exit frame and time, transit time (ps) |
| `STEM_permeation_visits.csv` | one row per lumen visit: frames inside, dwell time, side entered from and left to, censoring |
| `STEM_permeation_frames.csv` | per frame and per channel/species: molecules inside, cumulative net crossings, `n(t)` |
| `STEM_permeation.png` | molecules inside, cumulative net crossings and water `n(t)` against time (written with `--outputs figures` or `key-figures`) |
| `STEM_permeation_positions.csv` | axial traces of the molecules with counted events (see above) |
| `STEM_permeation_profile.csv` | per channel, species and 0.5 Å bin: molecules per frame, density per Å and, when defined, the free energy in kT with its interval (and in kcal/mol with `--temperature`) |
| `STEM_permeation_positions.png`, `STEM_permeation_profile.png` | the traces and the axial density and free energy (written with `--outputs figures`) |

The trajectory JSON (`-o`) records the permeation files and totals under
`metadata.permeation`; with `--permeation-only` it holds only that record.
The figure has no title or labels beyond the axes; `--annotate` adds a title
with the event counts.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--permeation SPECIES` | | species to follow (repeatable) |
| `--permeation-only` | off | skip the structural analysis; needs an explicit extent or a region |
| `--channel-extent ZMIN,ZMAX` | profile | axial extent along `--axis`, Å |
| `--channel-radius R` | none | lateral radius of the explicit channel, Å |
| `--membrane-planes SELECTION` | | extent from the lower and upper planes of these atoms in every frame, instead of `--channel-extent` |
| `--permeation-pore X,Y,Z` | `--origin` | axis point of one pore (repeat for several pores) |
| `--permeation-region FILE` | | region map or prepared region definition |
| `--lumen-tolerance Å` | 1 (profile), 0 (region) | added to the profile radius / dilation of the region |
| `--charge NAME=q` | table | charge override |
| `--voltage MV` | | known voltage |
| `--efield V/nm`, `--efield-from FILE` | | applied static field |
| `--water-mapping N` | 4 | water molecules per Martini water bead |
| `--permeation-stride N` | 1 | read every Nth frame |
| `--msd-blocks N` | 5 | segments for the `D_n` interval |
| `--permeation-convention` | `compartment` | event rule: `compartment` or `strict` (see above) |
| `--event-blocks N` | 10 | time blocks for the dispersion index of the rate, pd and current intervals, and for the free-energy interval |
| `--temperature K` | | simulation temperature; adds the free energy in kcal/mol |
| `--max-traces N` | 200 | most molecules per channel and species in the position traces; 0 keeps all |

`--start`, `--stop`, `--selection`, `--align`, `--confidence`,
`--bootstrap-replicates` and `--bootstrap-seed` of `crevice trajectory` apply
as well.

## Limitations

- Counts depend on the axial extent, which defines what a transit is; the
  lateral limit only attributes transits to a pore. Check how they change when
  you vary these, and check `transits_outside_lumen`.
- Events need the molecule to be seen inside in at least one saved frame;
  times are resolved only to the frame interval.
- pf and pd assume a stationary equilibrium trajectory. Few events give wide
  intervals; a single short trajectory rarely pins down pd.
- The current counts complete crossings of the channel extent. A conductance
  needs a trajectory with an applied voltage; it is a simulation estimate, not
  a validated single-channel conductance.
- Interval widening assumes events correlated over less than a tenth of the
  run; slower processes (rare gating) are still underestimated.
- Synthetic tests check the counting, unwrapping, multi-pore and
  coarse-grained bookkeeping and the `D_n` estimator against known answers;
  they do not validate any particular channel. `benchmarks/nongeom/` holds the
  longer checks: counts against an independent counter and exact rates,
  the dependence of counts on the lumen definition and on the saving
  interval, and the coverage of the pf, pd and current intervals.

## References

- Zhu F, Tajkhorshid E, Schulten K (2004a). Theory and simulation of water
  permeation in aquaporin-1. *Biophys. J.* 86, 50–57.
- Zhu F, Tajkhorshid E, Schulten K (2004b). Collective diffusion model for
  water permeation through microscopic channels. *Phys. Rev. Lett.* 93,
  224501.
- Kutzner C, Grubmüller H, de Groot BL, Zachariae U (2011). Computational
  electrophysiology: the molecular dynamics of ion channel permeation and
  selectivity in atomistic detail. *Biophys. J.* 101, 809–817.
  doi:10.1016/j.bpj.2011.06.010. The per-channel counting rule is that of
  its GROMACS implementation (`swap/swapcoords.cpp`).
- Gapsys V, Kopec W. *ion permeation analyzer*: membrane-referenced ion
  permeation counting with MDAnalysis,
  <https://github.com/Kopec-Lab/ion-permeation-analyzer> (commit 0b8a263,
  April 2026). Its phosphate-plane reference is the model for
  `--membrane-planes`.

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
| Region map | `--permeation-region FILE` (`.dx` map or prepared region `.json`) | in a region voxel (dilated by `--lumen-tolerance`, default 0) |

The extent is measured as the projection `x·axis` of reference-frame
coordinates (Å). Without `--channel-radius` every molecule within the extent
counts as inside, which is only appropriate when nothing else (lipid, other
pores) lies within that slab. A frame without a resolved profile reuses the
latest earlier profile; such frames are counted in the summary.

## How crossings are counted

In each frame a molecule is **below** the channel (`x·axis < ZMIN`),
**above** it (`> ZMAX`), **inside** the lumen, or **out**: within the extent
but outside the lumen laterally, for example in the protein wall or the
membrane. A **permeation event** is a complete crossing: the molecule leaves
one side, stays inside without interruption and leaves to the other side.
Below → inside → above is *up* (along `+axis`), the reverse is *down*. This is
the definition of Zhu, Tajkhorshid and Schulten (2004a).

- Entering and leaving on the same side is not an event.
- Leaving the lumen laterally (*out*) ends the visit, and the molecule's
  origin is forgotten.
- Periodic boundaries are handled per molecule: positions are taken as the
  nearest image about the channel centre, and frame-to-frame displacements
  are unwrapped. A molecule that reaches the other side through the periodic
  bulk is therefore not counted, and one that crosses the channel several
  times is counted each time.
- A molecule that crosses the whole extent between two saved frames is an
  **unresolved crossing**: listed, not counted, and reported as a warning.
  Save frames more often if there are many.

Entry and exit frames, transit times (exit time minus entry time) and dwell
times are resolved to the saved-frame interval. Visits cut by the start or end
of the trajectory are marked as censored and left out of the dwell statistics.

## Water permeability

Two permeabilities characterise a water channel (Zhu, Tajkhorshid and
Schulten, 2004a,b):

- **Diffusive permeability** `pd = v_w q0`, where `q0` is the rate of
  permeation events in one direction (at equilibrium,
  `q0 = (up + down) / 2T`) and `v_w = 18 cm³ mol⁻¹ / N_A` is the volume of one
  water molecule. Its interval is the exact (Garwood) Poisson interval of the
  event count.
- **Osmotic permeability** `pf = v_w D_n`, from the collective coordinate
  `n(t)`: in each frame interval the axial displacements of the molecules
  inside the channel are summed and divided by the channel length `L` (only
  the part of a displacement inside the extent counts), and `n(t)` is the
  running sum. At equilibrium `n` diffuses, `⟨Δn²⟩ = 2 D_n τ`. CREVICE fits
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
`+axis`. Two intervals are given: a Poisson interval (normal approximation to
the difference of two Poisson counts; with no events, the exact upper bound
on the number of crossings) and a block-bootstrap interval of the net charge
carried per block of frames.

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
| `STEM_permeation.json` | summary: definitions, units, species, channel definitions, counts, rates, occupancy, dwell and transit times, pf/pd and current with their intervals, voltage, warnings and limitations |
| `STEM_permeation_events.csv` | one row per event or unresolved crossing: channel, species, residue, atom index, direction, status, entry/exit frame and time, transit time (ps) |
| `STEM_permeation_visits.csv` | one row per lumen visit: frames inside, dwell time, side entered from and left to, censoring |
| `STEM_permeation_frames.csv` | per frame and per channel/species: molecules inside, cumulative net crossings, `n(t)` |
| `STEM_permeation.png` | molecules inside, cumulative net crossings and water `n(t)` against time (written with `--outputs figures` or `key-figures`) |

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
| `--permeation-pore X,Y,Z` | `--origin` | axis point of one pore (repeat for several pores) |
| `--permeation-region FILE` | | region map or prepared region definition |
| `--lumen-tolerance Å` | 1 (profile), 0 (region) | added to the profile radius / dilation of the region |
| `--charge NAME=q` | table | charge override |
| `--voltage MV` | | known voltage |
| `--efield V/nm`, `--efield-from FILE` | | applied static field |
| `--water-mapping N` | 4 | water molecules per Martini water bead |
| `--permeation-stride N` | 1 | read every Nth frame |
| `--msd-blocks N` | 5 | segments for the `D_n` interval |

`--start`, `--stop`, `--selection`, `--align`, `--confidence`,
`--bootstrap-replicates` and `--bootstrap-seed` of `crevice trajectory` apply
as well.

## Limitations

- Counts depend on the channel definition: its extent, its lateral limit and
  the lateral tolerance. Check how they change when you vary these.
- Events need the molecule to be seen inside in at least one saved frame;
  times are resolved only to the frame interval.
- pf and pd assume a stationary equilibrium trajectory. Few events give wide
  intervals; a single short trajectory rarely pins down pd.
- The current counts complete crossings of the channel extent. A conductance
  needs a trajectory with an applied voltage; it is a simulation estimate, not
  a validated single-channel conductance.
- Synthetic tests check the counting, unwrapping, multi-pore and
  coarse-grained bookkeeping and the `D_n` estimator against known answers;
  they do not validate any particular channel.

## References

- Zhu F, Tajkhorshid E, Schulten K (2004a). Theory and simulation of water
  permeation in aquaporin-1. *Biophys. J.* 86, 50–57.
- Zhu F, Tajkhorshid E, Schulten K (2004b). Collective diffusion model for
  water permeation through microscopic channels. *Phys. Rev. Lett.* 93,
  224501.

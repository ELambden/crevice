# Tool guide

There is one page for each CREVICE command or command family. Every page has the
same layout:

1. **Scientific question**: what the tool measures, and what it does not.
2. **Method**: how the measurement is made.
3. **Assumptions**: conditions the result depends on.
4. **Key parameters and units**.
5. **Outputs**.
6. **Python equivalent**.
7. **Testing and validation status**: labelled with the levels below.
8. **Known limitations**.
9. **Command-line options**: generated from the `crevice` argument parser, so it
   always matches the code the site was built from.

(validation-levels)=
## Validation levels

CREVICE records evidence at four separate levels. A tool can be strong at one
level and have no evidence at another.

| Level | Meaning | Example |
|---|---|---|
| **Software / synthetic** | Unit and contract tests, analytic controls and mathematical atom walls with known answers | continuous-path clearance checks on synthetic walls |
| **Real-input observation** | The tool ran on a real structure or trajectory with stated settings, and the output was inspected | 1GRM profile and cast on 0.25 Å grids |
| **Native viewer check** | The generated PyMOL/VMD/ChimeraX scene was opened in that viewer, and its geometry, colours and camera were checked | representative 1GRM and 4PYP scenes |
| **Biological / functional validation** | Agreement with curated anatomy, independent data or experiments | **not established for any CREVICE result** |

CREVICE's benchmark registry has no curated entries. Channel axes and regions are
chosen by the analysis or by the user, not taken from reviewed references. A
cavity count, volume or radius from a real protein is therefore a measurement
under stated settings, not a validated biological result. Each tool page ends
with a **Testing and validation status** section that uses these levels.

## Commands

| Command | Question | Page |
|---|---|---|
| `profile` | How wide is a through-channel along its axis? | [profile](profile.md) |
| `cast` | What 3D space does a pocket, cavity or channel occupy? | [cast](cast.md) |
| `cavities` | Are there enclosed grid voids? | [cavities](cavities.md) |
| `tunnels` | Which widest grid paths lead from a point to the outside? | [tunnels](tunnels.md) |
| `residues` | Which residues line a sampled profile? | [residues](residues.md) |
| `network` | Which residues touch each other and the channel? | [network](network.md) |
| `residue-evidence` | Which residues form a measured cavity boundary, and which residues contact them? | [residue-evidence](residue-evidence.md) |
| `hydration` | Where are explicit waters, and how exposed are residues? | [hydration](hydration.md) |
| `trajectory` | How do through-channel profiles and contacts vary across frames? | [trajectory](trajectory.md) |
| `cavity-trajectory` | How do a cavity's volume, widths and boundary residues change over an MD run? | [cavity-trajectory](cavity-trajectory.md) |
| `region-init`, `region-prepare`, `region-trajectory`, `region-compare` | How do named, versioned regions behave over identical frames? | [regions](regions.md) |
| `publish` | Can I get a complete static-structure bundle in one run? | [publish](publish.md) |
| `fetch`, `analyze`, `features`, `benchmark`, `static-suite` | Utility and batch commands | [utilities](utilities.md) |

```{toctree}
:maxdepth: 1
:hidden:

profile
cast
cavities
tunnels
residues
network
residue-evidence
hydration
trajectory
cavity-trajectory
regions
publish
utilities
```

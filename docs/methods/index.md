# Method notes

These notes give the detailed definitions, algorithms, parameters, units,
assumptions and limitations behind the [analysis tools](../tools/index.md).

| Note | Covers |
|---|---|
| [Atomic radii](atomic-radii.md) | the standard radius table, `--radii`/`radii=` presets and radius files (JSON, CSV, HOLE `.rad`), precedence, how the choice reaches workers, provenance |
| [Pore fingerprint definitions](fingerprint.md) | the path coordinate *s*, lining shell, Wimley–White hydrophobicity and Gaussian smoothing, water along the path, the conductance estimate with access resistance and the HOLE divisor, limitations |
| [Curved and branched pathways](paths.md) | `paths`: grid, bulk and exposed mouths, widest routes, half-space sides, distinct routes, junctions, centre-line relaxation, free-sphere radius, lining shell, unresolved cases |
| [Cavity tracking across frames](cavity-tracking.md) | `track-cavities`: fixed aligned lattice, components, Jaccard and containment links, continuation, birth/death/split/merge, memory, thresholds |
| [Cavity trajectories](../CAVITY_TRAJECTORY_METHODS.md) | `cavity-trajectory` observables, fixed coordinates, sectional widths, sampling uncertainty, changing residue evidence |
| [Trajectory profiles and residue evidence](../TRAJECTORY_RESIDUE_METHODS.md) | alignment and PBC contract, profile uncertainty, boundary attribution, nonlocal contacts, seeded casts, reading evidence figures |
| [Residue hydration and solvent exposure](../HYDRATION_METHODS.md) | water contacts, SASA, periodic identities, instantaneous-cavity water membership, temporal statistics |
| [Coordinated hydration and residue evidence](../HYDRATION_NETWORK_VIEWS.md) | density maps, typed interactions, colour scales, coordinated viewer scenes |
| [Physical obstacles and numerical sensitivity](../CAVITY_OBSTACLES.md) | obstacle selection, periodic imaging, lattice controls, sensitivity checks |
| [Display surface smoothing](../SURFACE_SMOOTHING.md) | `--smooth X` constrained display surfaces, `--surface-smoothing`, choosing X; display only, measurements unchanged |
| [Interactive view surfaces](interactive-view.md) | what `crevice view` draws: the native display isosurfaces, layers, colours and opacities, offline rendering |
| [Named biological-region workflow](../REGION_WORKFLOW.md) | region definitions, landmark proposals, review overlays, comparisons |
| [Biological identity and region review](../BIOLOGICAL_REGIONS.md) | benchmark identities and assemblies, site landmarks, modified amino acids |

```{toctree}
:hidden:

atomic-radii
fingerprint
paths
cavity-tracking
/CAVITY_TRAJECTORY_METHODS
/TRAJECTORY_RESIDUE_METHODS
/HYDRATION_METHODS
/HYDRATION_NETWORK_VIEWS
/CAVITY_OBSTACLES
/SURFACE_SMOOTHING
interactive-view
/REGION_WORKFLOW
/BIOLOGICAL_REGIONS
```

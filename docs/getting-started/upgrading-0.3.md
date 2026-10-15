# Upgrading from 0.2 to 0.3

CREVICE 0.3 changes **which files** the bundle commands write by default. No
measurement changes: profiles, cast volumes, residue contacts, interactions
and every other number are identical to 0.2 for the same settings.

## What changed

- **Compact defaults.** `publish`, `cast`, `analyze`, `residue-evidence`,
  `hydration`, `trajectory`, `cavity-trajectory` and `region-trajectory` write
  the manifest, the data (JSON, CSV, NPZ, measured DX maps), the key figures
  and PyMOL, VMD and ChimeraX scenes. A 1GRM `publish` writes 40 files instead
  of 96 (see {ref}`output-policy`).
- **Water analyses are opt-in** (`--outputs default,hydration` or `--all`)
  except in `hydration` and in `cavity-trajectory`/`region-trajectory` with
  `--water-membership`. This covers water contacts and occupancy,
  solvent-accessible area, water bridges, water density and the hydration
  scenes.
- **Typed residue–residue interactions still run by default** (hydrogen bonds,
  salt-bridge candidates, aromatic contacts and residue conformations of the
  lining or boundary residues), now without the water analyses. Their numbers
  equal those of 0.2; the water bridges and residue–water hydrogen bonds are
  reported only with `hydration`.
- **`publish` extras are opt-in**: the separate cavity search, tunnel search
  and residue network (`--outputs default,cavities,tunnels,network`), the
  annotated radius figure (`figures`), the input copy and dummy-atom PDBs
  (`structures`), and the alias and variant scenes `*_pore_cast.*`,
  `*_pore_spheres.pml` and `cast`'s `*_end_view.*` (`extra-views`).
- **New options** on these commands: `--outputs`, `--all`, `--dry-run` and the
  deprecated `--legacy-outputs`; on every command: `--config FILE` and
  `--print-config` (JSON settings files). See {ref}`settings-files`.
- **Manifests record the run**: version, command line, every effective setting
  and its source, units, input file hashes and the output plan (`run` key).
- **Python**: the bundle writers keep their signatures and, called directly
  without `outputs=`, still write the complete 0.2 bundle (including the
  cavity, tunnel and network searches). `outputs="default"` gives the compact
  CLI set. New: {mod}`crevice.static_workflow`
  ({func}`~crevice.static_workflow.compute_static_workflow` measures without
  writing; {func}`~crevice.static_workflow.export_static_workflow` writes),
  {mod}`crevice.outputs` and {mod}`crevice.config`.

## Getting the old files back

| You used | In 0.3 |
|---|---|
| a 0.2 `publish`/`cast` bundle as is | add `--all` (everything, plus the fingerprint and HTML view), or `--legacy-outputs` for the exact 0.2 file set |
| hydration files from `publish`, `cast`, `analyze`, `residue-evidence`, `trajectory`, `cavity-trajectory`, `region-trajectory` | add `--outputs default,hydration` |
| `publish` cavities, tunnels and network | add `--outputs default,cavities,tunnels,network` |
| `*_pore_cast.pml` / `.vmd` / `.cxc` | open `*_volume.pml` / `.vmd` / `.cxc` (the same scene), or add `extra-views` |
| `*_profile_radius_annotated.png` | add `--outputs default,figures` |
| `*_structure.<ext>` (input copy), `*_pore_cast.pdb` | add `structures` |
| `--skip-hydration` in scripts | still accepted; now usually already the default |

`--legacy-outputs` reproduces the 0.2 file set exactly. It is **deprecated**:
it prints a warning, stays available during the 0.3 series and will be
removed in a later minor release. Prefer `--all` or an explicit `--outputs`
list.

## Checking a run before it starts

```bash
crevice publish 1GRM --out-dir out --all --dry-run        # prints the plan, writes nothing
crevice publish 1GRM --out-dir out --print-config > run.json
crevice publish 1GRM --config run.json --samples 101       # flags override the file
```

# `publish`: static-structure bundle

## Scientific question

`publish` is a convenience workflow rather than a separate measurement. It runs
the static analyses on one structure and writes a shareable bundle of data,
figures and viewer scenes with a manifest.

## What it runs

1. **Profile** ([`profile`](profile.md)) with the profile options.
2. **Void cast** for the scenes and measured volume. `--cast-mode auto`
   (default) uses a channel cast when a unique through-profile resolves. If it
   does not, `publish` warns, records the reason and falls back to a dominant
   [rolling-probe cast](cast.md) of all heavy atoms. `channel`, `cavity`, `all`
   and `rolling` force a mode. Explicit `channel` fails if the profile does not
   resolve.
3. **Residues** along the profile, **cavities**, **tunnels** and the
   **network**, unless skipped with `--skip-cavities`, `--skip-tunnels` or
   `--skip-network`. The bundle's network uses the low-level surface-gap
   metric with heteroatoms excluded, which differs from the `crevice network`
   defaults (see [`network`](network.md)).
4. **Hydration** extensions, unless `--skip-hydration` (see
   [`hydration`](hydration.md)).

Channel casts fill every enclosed probe centre in each section, not only the
inscribed sphere. By default the display cast extends 2 Å beyond each mouth
(`--cast-extension`, display only). That extension is excluded from volumes
and radii.

### Capped channels (`--lateral-exits`)

With `--lateral-exits`, a profile may resolve as a capped channel whose end
leaves through a lateral exit (see [`profile`](profile.md), "Capped channels").
The manifest's `profile_status` then also records `path_type`
(`capped_channel_lateral_exit` or `through_channel`), `exit_types` and
`exit_leg_counts` per end. The measured channel cast (`void_cast.total_volume`,
`PREFIX_volume.dx`, `PREFIX_pore_cast.pdb`) still covers only the axial lumen
between the mouths. Every distinct exit leg (for example both C2-related
portals of TWIK-1) also gets its own probe-swept cast, kind `lateral_exit`:
grid points beyond the capped mouth and outside bulk solvent, inside the free
spheres of the leg's samples, swept by the enclosure probe and connected to the
leg's centre line (`crevice.channel_exits.lateral_exit_casts`). Each leg's
volume is reported separately, never added to the axial volume:
`metadata.lateral_exit_casts` in `PREFIX_void_cast.json` and the manifest, one
`lateral_exit_<end>_<leg>` row in `PREFIX_void_cast.csv`,
`profile_status.lateral_exit_cast_volumes_A3`, and the binary map
`PREFIX_exit_casts.dx`. The scenes draw the leg casts as part of the cast
surface, joined to the channel at its mouth, so the widening towards each
portal is visible. `--exit-centre-lines` also draws each leg's centre line as a
thin blue tube (RGB 0.16, 0.47, 0.84; PyMOL object `crevice_exits`). The two
radius plots draw each leg as a bluish-green dash-dot continuation beyond the
mouth. No labels are drawn. Without `--lateral-exits` the bundle is unchanged.

### Display guides and colours

Channel-mouth guides are off by default: no orange rings in the PyMOL, VMD and
ChimeraX scenes and no orange mouth bands or lines in the radius plots.
`--mouth-guides` restores both. The mouth positions are always recorded
(`metadata.channel_mouths` in the profile JSON, `mouth_rings` in
`PREFIX_scene.json`, `channel_mouths` in `PREFIX_volume_metadata.json`).
A resolved channel cast is teal (RGB 0.08, 0.58, 0.63). Every other cast,
including the rolling-probe cast written when no channel resolves (for example
4PYP), is violet (`#a855f7`, RGB 0.659, 0.333, 0.969), and the width-profile
lines of non-channel casts (`cavity-trajectory`) use the same violet.

### Enclosure probe

By default (`--enclosure-radius auto`) the profile's enclosure probe is chosen
automatically, as described on the [`profile`](profile.md) page. The manifest's
`profile_status.enclosure_probe` records the mode, the chosen probe, the reason
and the probes tried; the profile JSON has the full record. When no probe
resolves a channel, `<prefix>_profile_status.json` holds the full selection
record, `profile_status.enclosure_probe` has `chosen_A: null`, and the
fallback rolling-probe cast uses 0.8 Å (`fallback_cast_probe_A`), as before.
`--cast-mode rolling` also uses 0.8 Å unless a number is given. An explicit
`--enclosure-radius` value runs exactly that probe, as before.

### When no through-profile resolves

The unresolved-profile fallback still runs the requested extras. No profile is
invented: no profile JSON, CSV or radius figure is written. Instead:

- **Cavities:** {func}`crevice.detect_cavities` defaults, with the publication
  atom selection as obstacles.
- **Tunnels:** {func}`crevice.find_tunnels` defaults. The start point is the
  measured fallback cast sample with the largest atom clearance.
- **Network:** the standard `publish` network. Its region node covers every
  measured cast sample centre.
- **Status records:** `<prefix>_profile_status.json` records
  `status: "unresolved"`, the reason and an interpretation. The manifest
  carries `profile_status` and `fallback_analyses`, where each extra is marked
  `completed`, `skipped_by_request` (a `--skip-*` flag) or `not_completed` with
  a reason. When a profile resolves, the manifest has
  `profile_status: {"status": "resolved", ...}`. The CLI summary line ends
  with `profile=resolved` or `profile=unresolved`.

An unresolved profile does not show that the protein is closed. Focus,
maximum-volume, minimum-clearance or hydrogen-inclusive requests do not
silently switch to the fallback: they fail instead.

## Assumptions

- Heteroatoms are **kept** unless you pass `--exclude-hetero`. This differs
  from `cast`, which defaults to protein only. Record the command and check the
  input report.
- Focus, maximum-volume, minimum-clearance or hydrogen-inclusive requests do
  not silently switch to the rolling fallback.

## Key parameters and units

| Option | Default | Unit | Meaning |
|---|---:|---|---|
| `--cast-mode` | `auto` | | `auto`, `channel`, `cavity`, `all`, `rolling` |
| `--cast-spacing` | 0.5 | Å | cast grid spacing |
| `--cast-min-radius` | 0.0 | Å | minimum atom-surface clearance for cast points |
| `--cast-extension` | 2.0 | Å | display-only extension beyond mouths |
| `--cast-max-cavities` | 1 | count | maximum rolling-cast regions |
| `--smooth X` | 0.4 | Å | constrained display smoothing; `0` gives raw voxels |
| `--surface-smoothing W` | off | Å | alternative single-Gaussian display interpolation; cannot be combined with `--smooth` |
| `--focus X Y Z`, `--focus-radius` | none, 8.0 | Å | crop the cast to local spheres (repeatable) |
| `--min/max-component-volume` | 0, none | Å³ | inclusive limits on final components |
| `--mouth-guides` | off | flag | draw mouth rings in the scenes and mouth bands/lines on the radius plots |
| `--exit-centre-lines` | off | flag | with `--lateral-exits`: also draw each exit leg's centre line as a thin blue tube |
| `--radii` | standard table | set | atomic radius set: `bondi`, `hole`, `charmm_like` or a JSON/CSV/HOLE `.rad` file; see [Atomic radii](../methods/atomic-radii.md) |

For finer channel refinement, use `--section-spacing 0.25 --cast-spacing 0.25
--cast-max-grid-points 2000000`.

## Outputs

Every table is a CSV whose column names carry units (`_A` = Å, `_A2` = Å²,
`_A3` = Å³, `_ps`, `_degrees`):

| File | Contents |
|---|---|
| `*_profile.csv` | one row per profile sample: `x_A`, `y_A`, `z_A`, `radius_A`, `raw_clearance_A`, `axis_position_A`, nearest atom and residue |
| `*_residue_contacts.csv` | residues along the profile: `min_distance_A`, `mean_distance_A`, `role`, `influence_score` |
| `*_void_cast.csv` | one row per cast region: `volume_A3`, centre, clearance range |
| `*_cavities.csv` | one row per cavity: centre, `max_radius_A`, `volume_A3`, `grid_points` |
| `*_tunnels.csv`, `*_tunnel_points.csv` | one row per tunnel (`length_A`, `bottleneck_radius_A`, ...) and per path node |
| `*_network_nodes.csv`, `*_network_edges.csv` | residues with `degree`/`residue_degree`; contacts with `interaction`, `distance_A`, `cutoff_A` |
| `*_connectivity.csv`, `*_connectivity_groups.csv` | with `--residue-groups`: paths to the lining and group contacts |
| hydration tables | see [`hydration`](hydration.md#outputs) |

Other files: the measured `*_volume.dx` and display maps; PyMOL/VMD/ChimeraX
scenes and render scripts; the radius, annotated radius, residue-contact,
cavity-summary, network-summary and network-chord figures; `*_pore_cast.pdb`,
`*_cavities.pdb` and `*_tunnels.pdb` dummy atoms; and these JSON records,
kept because they hold nested settings, metadata or scene definitions that
other commands and the viewer scripts read: `*_manifest.json` (file index,
cast summary, `profile_status`), `*_input_report.json` (input provenance),
`*_profile.json`, `*_void_cast.json`, `*_cavities.json`, `*_tunnels.json`,
`*_network.json` (full records), `*_scene.json` and `*_volume_metadata.json`
(viewer scene data). A default 1GRM run writes about 85 files; use the
`--skip-*` options to leave out analyses you do not need. No HTML report,
3D preview image or tunnel summary figure is written.

Figures (colours as on the [`profile`](profile.md), [`residues`](residues.md),
[`cavities`](cavities.md) and [`network`](network.md) pages) have no titles,
value call-outs, mouth labels or residue-landmark names unless `--annotate` is
given; the text is kept in each PNG's metadata. The residue-contact figure has
a legend naming the contact roles, and the network summary a legend naming the
contact classes, because their colours encode categories. The
viewer scenes show a grey-blue cartoon (25% transparent) around an opaque cast
surface (teal for a channel, violet for any other cast) and, for a capped
channel resolved with `--lateral-exits`, the exit-leg casts joined to it; they
contain no text. Mouth rings (`--mouth-guides`) and exit centre lines
(`--exit-centre-lines`) are optional. Worked bundles: [1GRM](../examples/1grm-profile-cast.md) and
[3UKM](../examples/3ukm-lateral-exits.md).

Display surfaces use `--smooth X` (default 0.4 Å; see
[Display-surface smoothing](../SURFACE_SMOOTHING.md)). Measured files are the
same for every X.

## Python equivalent

```python
from crevice import load_structure, write_static_publication_bundle

frame = load_structure(".crevice/pdb/1GRM.cif")
manifest = write_static_publication_bundle(
    frame, structure_path=".crevice/pdb/1GRM.cif", output_dir="results/1GRM",
    prefix="1GRM", samples=81)
```

The CLI also appends hydration outputs and an input report, which this call
does not do.

## Testing and validation status

- **Software / synthetic:** synthetic tests cover resolved and unresolved
  `publish`, skipped extras and the CLI summary line
  (`tests/test_publish_fallback_and_sticks.py`).
- **Example observations (public entries, 0.25 Å grids for 1GRM):** 1GRM gives
  a 101-sample profile with a minimum continuous-path radius of 1.330 Å and a
  cast volume of 292.4 Å³, with `--enclosure-radius 0.8` and with the automatic
  choice, which selects 0.8 Å for this input. 4PYP is `unresolved` at every
  candidate probe, so the fallback applies; its cast is 2439.6 Å³, with zero
  enclosed cavities and three tunnels.
- **Native viewer check:** scenes of this kind were opened and checked in
  PyMOL, VMD and ChimeraX for representative cases during development. A newly
  generated scene is not checked automatically; inspect it in the viewer.
- **Biological / functional:** not established.

## Command-line options

```{argparse}
:module: crevice.cli
:func: build_parser
:prog: crevice
:path: publish
```

# Inputs and selections

A geometric result depends on which atoms count as walls. Before trusting a
number, check which atoms a run actually used. Every command records that
choice in its input report or manifest.

## Supported inputs

| Input | Read by | Notes |
|---|---|---|
| PDB (`.pdb`) | built-in parser | ATOM/HETATM records; one model is selected (`--model`) |
| mmCIF (`.cif`) | Gemmi when the `structures` extra is installed, otherwise the built-in parser | Gemmi reads `_entity` types, so modified polymer residues deposited as HETATM stay in the polymer |
| GRO snapshot (`.gro`) | MDAnalysis | coordinates converted from nm to Å; `--md-selection` chooses atoms |
| MD trajectory plus topology | MDAnalysis | used by `trajectory`, `cavity-trajectory`, `hydration`, `region-*`; GRO + XTC is the most exercised combination |
| Several static files | built-in / Gemmi | `trajectory a.pdb b.pdb ...` treats them as frames |

CREVICE does not expand biological-assembly operators while loading. To analyse
an assembly, fetch the assembly file (`--assembly N`) or prepare one yourself.

## Fetching structures

Any command that takes a structure also accepts an RCSB accession in the same
place. If a file with that name exists, the file is used, so a local file is
never replaced by a download.

```bash
crevice fetch 1GRM 4PYP --cache-dir .crevice/pdb --json downloads.json
crevice analyze 6MVY --assembly 1 --out-dir results/6mvy-assembly
```

Downloads go to `--cache-dir` (default `.crevice/pdb`) and are reused later.
`--offline` never downloads, and `--refetch` forces a new copy. mmCIF is the
default format (`--format cif`). Before a download is cached, CREVICE checks that
it parses in the requested format. It then writes a `<file>.provenance.json`
sidecar with the URL, SHA-256 hash, size and retrieval time.
{func}`crevice.verify_cached_structure` re-checks a cached file against that
record. If the file names a different accession from the one requested, that is
an error. A classic ID and its extended form name the same entry, so
`pdb_00001grm` matches a file stating `1GRM` and vice versa. Placeholder accessions, as found in assembly files, are recorded as
unverified.

### AlphaFold DB models

AlphaFold DB identifiers are accepted in the same places: `AF-P02920-F1`
(latest version), `AF-P02920-F1-model_v6`, or the download file name
`AF-P02920-F1-model_v6.pdb` (any case). An AlphaFold file already on disk is
read like any other PDB or mmCIF file.

```bash
crevice fetch AF-P02920-F1 --format pdb
crevice cast AF-P02920-F1 --format pdb --out-dir results/lacy-model
```

Without a version, a cached copy is reused if there is one; otherwise the
latest version is looked up with the AlphaFold DB API
(`https://alphafold.ebi.ac.uk/api/prediction/<UniProt>`) and the file is
downloaded from `https://alphafold.ebi.ac.uk/files/`. It is cached under its
AlphaFold DB name with the same provenance sidecar and checked against
`_entry.id` (mmCIF) or the `DBREF` UniProt accession (PDB); a mismatch is an
error. `--assembly` does not apply to AlphaFold models; `--offline` and
`--refetch` work as for RCSB entries. See the
[AlphaFold example](../examples/alphafold-model.md).

AlphaFold models are predictions, not experimental structures. Their B-factor
column holds pLDDT (per-residue confidence, 0–100), and the input report marks
them `model_type: predicted`, `predictor: AlphaFold`. Treat geometry in
low-confidence regions with caution.

## Parser, model and chain identifiers

- `--parser auto|gemmi|builtin` chooses the mmCIF reader. The report always
  names the parser that ran.
- `--model N` picks a model (0-based) from multi-model files.
- `--chain-ids label|auth` picks the mmCIF chain namespace (default `label`)
  for either reader. Label IDs give ligands and assembly copies their own
  chain.

Every reader keeps blank alternate locations plus the first alternate location
(alphabetically) of each residue, and keeps hydrogens, waters, ions and
ligands. The two mmCIF readers still differ on purpose in two ways, and the
input report records both:

| | Gemmi (`--parser gemmi`, default when installed) | Built-in (`--parser builtin`, and all PDB files) |
|---|---|---|
| What counts as a heteroatom | non-polymer entities (`hetero_definition: nonpolymer`), so a modified residue deposited as HETATM stays in the polymer | the HETATM record type (`group_pdb`) |
| Zero-occupancy atoms | dropped (`zero_occupancy_policy: excluded`) | kept (`retained`) |

This matters for `--exclude-hetero` and for residue-level analyses. In 1GRM,
for example, the D-amino acids are HETATM records: the built-in reader treats
them as heteroatoms and the Gemmi reader as polymer.

Residues are labelled `CHAIN:RESNAMEnumber[insertion]`, for example `B:ALA3`.
GRO input has no chain field and uses `SYSTEM`, for example `SYSTEM:GLN282`.
Viewer PDB files may renumber atoms or residues. A `*_viewer.identities.json`
sidecar maps them back, and evidence tables always use the original labels.
Keep that sidecar with the viewer PDB.

## Elements and atomic radii

Every clearance is measured to van der Waals surfaces, so each atom's element
sets its radius.

- **Element column first.** If the file gives an element (PDB columns 77-78,
  mmCIF `type_symbol`) and it is a real periodic-table symbol, it is used as
  written. `HG` is mercury and `CA` in the element column is calcium.
- **Name fallback.** Without an element column (for example GRO files) the
  element is inferred from the atom name. Standard protein carbons such as `CA`
  are carbon. `CA`, `NA`, `CO` and `NI` are read as metal ions only for a
  single-atom residue of the same name, so ligand carbon `CAB` and heme
  nitrogen `NA` are not turned into calcium or sodium.
- **Force-field ion names.** Without an element column, monatomic ions of the
  common MD force fields are recognised by their residue name *and* atom name
  together (table below), so C-alpha `CA`, proline `CD` and heme `NA` are never
  affected. An explicit element column still wins. Other ion names fall back to
  the name rules and can be misread as light elements; give such files an
  element column.
- **Radii.** Values come from Bondi (1964), with Mantina et al. (2009) for
  main-group elements Bondi did not tabulate; {func}`crevice.radii.radius_source`
  names the source of each element. A recognised ion residue (table below,
  residue *and* atom name) instead gets its CHARMM36 Lennard-Jones Rmin/2 radius
  by default (Na 1.41075 Å, Cl 2.27 Å, K 1.76375 Å...; full list and source in
  [Atomic radii](../methods/atomic-radii.md#ion-radii)); Cu, F, Br and I, which
  CHARMM36 does not parameterise there, keep their element radius. An element without a tabulated radius
  (most transition metals and lanthanides, for example Gd) keeps its symbol but
  has no radius, so the atom is **excluded** (next section). `--radii` (Python
  `radii=`) selects another radius set; see [Atomic radii](../methods/atomic-radii.md).
- **Atoms without a radius.** When the radius set in force gives an atom no
  radius (no `RESNAME:ATOMNAME`, atom-name, ion, HOLE or element entry; a set's
  catch-all default does not count), PDB and mmCIF readers leave the atom out
  and print an `UnrecognisedElementWarning` naming the elements, atom counts and
  `RESNAME:ATOMNAME` keys. The input report records them per element
  (`excluded_no_radius_atoms`, `unsupported_radius_elements`) and per name
  (`excluded_no_radius_names`), with the set checked (`radius_set`). To include
  them, which matters for designed or modified residues with unusual elements,
  write a radius file and pass it with `--radii FILE` (Python: `radii=FILE`, also
  to `load_structure`, or a `use_radii(FILE)` block around loading and
  analysis):

  ```bash
  python -c "from crevice.radii import write_radius_template as w; w('radii.csv', elements=['GD'])"
  # edit radii.csv: uncomment and fill e.g.  element,GD,2.10  or  atom,XGD:GD1,2.10
  crevice cast designed.pdb --radii radii.csv --out-dir results/designed
  ```

  The template starts from the default set (`base,default,`), so only the new
  rows need values; the file's SHA-256 and the added radii are recorded in every
  output. MD trajectories cannot drop atoms silently (frames must match the
  MDAnalysis selection), so loading stops with the same explanation and a
  selection hint such as `'(protein) and not resname XGD'`.
- **Deuterium.** `D` counts as hydrogen: it gets the hydrogen radius (unless a
  set lists `D` itself) and is removed with hydrogens by default.

| Residue / atom name | Element | Force field (source) |
|---|---|---|
| `SOD` / `SOD` | Na | CHARMM36 `toppar_water_ions.str` |
| `POT` / `POT` | K | CHARMM36 |
| `CLA` / `CLA` | Cl | CHARMM36 |
| `CAL` / `CAL` | Ca | CHARMM36 |
| `MG` / `MG` | Mg | CHARMM36; GROMACS `ions.itp` |
| `LIT` / `LIT` | Li | CHARMM36 |
| `CES` / `CES` | Cs | CHARMM36 |
| `BAR` / `BAR` | Ba | CHARMM36 |
| `RUB` / `RUB` | Rb | CHARMM36 |
| `ZN2` / `ZN` | Zn | CHARMM36 |
| `CD2` / `CD` | Cd | CHARMM36 |
| `NA`, `K`, `CL`, `CA`, `ZN`, `LI`, `RB`, `CS`, `CU`, `F`, `BR`, `I` (residue = atom name) | Na, K, Cl, Ca, Zn, Li, Rb, Cs, Cu, F, Br, I | GROMACS 2022.4 `share/top/*.ff/ions.itp` (AMBER, CHARMM27, GROMOS and OPLS-AA ports) |
| `CU1` / `CU` | Cu | GROMACS GROMOS `ions.itp` |
| `NA+`, `CL-` (residue = atom name) | Na, Cl | GROMACS `residuetypes.dat` (older names); AMBER legacy |
| `Li+`, `Na+`, `K+`, `Rb+`, `Cs+`, `F-`, `Cl-`, `Br-`, `I-` (residue = atom name) | Li, Na, K, Rb, Cs, F, Cl, Br, I | AMBER legacy names (Joung and Cheatham 2008) |

The CHARMM36 names were checked against the RESI/ATOM records of a
CHARMM-GUI copy of `toppar_water_ions.str` and the GROMACS names against the
GROMACS 2022.4 topology files. No AMBER file was available when the table was
written, so the legacy AMBER names follow the Joung-Cheatham ion parameter
naming and were not checked against a file. Names are compared case-insensitively.
Earlier CREVICE versions read CHARMM `SOD` as S and `POT` as P, and later
versions gave recognised ions their neutral-atom element radius (Na 2.27 Å);
ions in GRO-based analyses that include them as obstacles (for example a region
`obstacle_selection` of `not resname TIP3`, or the full-system hydration SASA
environment) now get the CHARMM36 ion radius. `--radii bondi` restores the
neutral-atom radii.

## Choosing atoms

The atom-selection flags differ between commands:

| Command(s) | Flag | Meaning |
|---|---|---|
| structure commands with GRO input | `--md-selection` (default `protein`) | MDAnalysis selection applied when reading a GRO file |
| `profile`, `residues`, `network`, `features`, `publish`, `analyze` | `--exclude-hetero` | drop HETATM records (kept by default) |
| `cavities`, `tunnels` (standalone) | `--include-hetero` | keep HETATM records (dropped by default; `analyze` keeps them) |
| most structure commands | `--include-hydrogen` | keep H atoms (dropped by default) |
| `cast` | `--selection protein\|all` | `protein`: polymer heavy atoms; `all`: also ligands and waters |
| `trajectory`, `cavity-trajectory`, `hydration`, `region-init` | `--selection` | MDAnalysis selection string (for example `protein and not name H*`) |
| `cavity-trajectory`, `region-init` | `--obstacle-selection` | separate MDAnalysis selection of atoms that act as walls |

### Walls versus displayed atoms

Keep two decisions separate: which atoms are **physical obstacles** in the
geometry, and which atoms are **shown or analysed as residues**. Removing water,
lipids and ligands gives a clean protein geometry. It also removes membrane and
ligand barriers, so the measured space can reach into regions that are
physically blocked. `cavity-trajectory --obstacle-selection` keeps chosen
non-protein atoms as walls while residue evidence stays protein-only. See the
[obstacle note](../CAVITY_OBSTACLES.md).

Static `publish` keeps heteroatoms unless you pass `--exclude-hetero`. `cast
--selection protein` drops them. Record the command you ran and check the input
report.

## Coordinates and units

Coordinates, radii, spacings and cutoffs are in ångström (Å). Volumes are in
Å³ and areas in Å². GRO/XTC nanometres are converted on reading. Coordinate
options take `x,y,z` without spaces (for example `--origin 0,0,12.5`).
`publish --focus X Y Z` is the exception: it takes three separate numbers.

## Trajectories: frames, alignment and periodic boundaries

- **Frames.** `--start`, `--stop` and `--stride` choose frames while reading, so
  they limit memory use. `--max-frames` is a guard against loading too many
  frames by accident. It is not a stride: a request above the limit fails and
  tells you how to fix it. Use `trajectory --inspect` to see frame counts first.
- **Alignment.** By default, matching CA atoms are fitted to the first frame
  with a proper rigid transform. Fit RMSDs are kept. A predeclared scaffold can
  be given as a JSON list of residue IDs (`--alignment-residues-json`). Choosing
  the scaffold after looking at the motion can bias the result.
- **Periodic boundaries.** `--pbc check` (the default) looks for bonds split
  across the periodic box. It uses peptide C–N links when the topology has no
  bonds. `--pbc unwrap` needs a reliably bonded single fragment. `--pbc none`
  turns screening off. CREVICE does not repair molecules split across the box.
  Make molecules whole before analysis.

Method detail: [Trajectory profiles and residue evidence](../TRAJECTORY_RESIDUE_METHODS.md#coordinate-and-sampling-contract).

## Worked example: a GRO snapshot

```bash
crevice cast system.gro --md-selection protein --out-dir output/system
```

`--md-selection protein` removes water, ions and lipids before the cast.
Hydrogens are dropped. The intake report lists how many atoms were selected and
excluded, with residue names. Periodic images are not repaired.

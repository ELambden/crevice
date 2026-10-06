# Atomic radii

Every quantity CREVICE measures is a distance to atom *surfaces*: a profile
radius is the clearance to the nearest van der Waals sphere, a cast is the
space a probe reaches between spheres, contacts and residue evidence use
surface distances, and SASA rolls a probe over the spheres. The atomic radii
therefore set every reported radius, volume, area and contact. This page
describes which radii are used, how to choose another set, and how the choice
is recorded.

## The default set

With no radius set chosen, every analysis uses the set
{data}`crevice.radii.DEFAULT_RADII` (preset name `default`):

- **Every atom except recognised ions**: the standard element table
  {data}`crevice.radii.VDW_RADII`: Bondi (1964) radii, with Mantina et al.
  (2009) for main-group elements Bondi did not tabulate and three legacy values
  (Fe, Mn, Co; {func}`crevice.radii.radius_source` names the source of each).
  The element comes from the file's element column or, without one, from the
  atom and residue names (see [Elements and atomic radii](../getting-started/inputs-and-selections.md#elements-and-atomic-radii)).
  An element without a tabulated radius has no radius in this set: structure
  readers exclude such atoms with a warning (see
  [Atoms without a radius](#atoms-without-a-radius)). Deuterium uses the
  hydrogen radius.
- **Recognised monatomic ion residues**: the CHARMM36 Lennard-Jones Rmin/2 of
  that ion ({data}`crevice.radii.ION_RADII`, table below).

### Ion radii

An atom is a *recognised ion* when its residue name and atom name together are
one of the monatomic ion names of the [ion table](../getting-started/inputs-and-selections.md#elements-and-atomic-radii) and its element
is that ion's element ({func}`crevice.radii.is_recognised_ion`). Both names must
match, so C-alpha `CA`, proline `CD` or heme nitrogen `NA` are never ions, and
a metal inside a larger residue (heme iron, a zinc in a ligand residue with
other atoms, a metal bound to a modified amino acid) keeps its element radius.
A single-atom residue named like its element (`ZN`/`ZN`, `NA`/`NA`, `CA`/`CA`,
`MG`/`MG`, `K`/`K`, `CL`/`CL`) is a recognised ion whether it comes from an MD
topology or a crystal structure, so a crystallographic zinc or calcium ion
gets the CHARMM36 value too.

| Element | CHARMM36 type | Radius (Å) | Element-table radius (Å) |
|---|---|---:|---:|
| Li | `LIT` | 1.2975 | 1.82 |
| Na | `SOD` | 1.41075 | 2.27 |
| Mg | `MG` | 1.185 | 1.73 |
| K | `POT` | 1.76375 | 2.75 |
| Ca | `CAL` | 1.367 | 2.31 |
| Rb | `RUB` | 1.90 | 3.03 |
| Cs | `CES` | 2.10 | 3.43 |
| Ba | `BAR` | 1.89 | 2.68 |
| Zn | `ZN` | 1.09 | 1.39 |
| Cd | `CAD` | 1.357 | 1.58 |
| Cl | `CLA` | 2.27 | 1.75 |

Source: the `NONBONDED` section of CHARMM36 `toppar_water_ions.str` (the
column after epsilon is Rmin/2, half the distance of the Lennard-Jones minimum
between two like atoms). The copy read is the one distributed with CHARMM-GUI
v3.7 (the file itself carries no version string; SHA-256
`41c6474470fa66e927198582a314f9e390cd6bf6fbd1449475d35e3346d72bdc`); its
GROMACS `SOD.itp`/`CLA.itp` give the same values (sigma × 2^(1/6) / 2). The ion
parameters are those of Beglov and Roux (1994) and later values cited in the
file; zinc is from Stote and Karplus (1995). Ions that file does not
parameterise (Cu, F, Br, I) keep their element radius. Rmin/2 is a
force-field size, not a crystal ionic radius (Shannon Na⁺ is about 1.02 Å); it
is closer to how the ion excludes neighbours in the simulation that produced
the coordinates than the neutral-atom Bondi value. Choose `--radii bondi` to
treat ions as neutral atoms, or set your own ion radii (`ion_radii` in a radius
file).

`radius_source(element, resname, atom_name)` reports the source of the radius
an individual atom gets.

## Choosing a radius set

| Where | How |
|---|---|
| Command line | `--radii default`, `--radii bondi`, `--radii hole`, `--radii charmm_like` or `--radii PATH` on every command that measures geometry (`profile`, `residues`, `cavities`, `tunnels`, `network`, `features`, `cast`, `publish`, `analyze`, `residue-evidence`, `hydration`, `trajectory`, `cavity-trajectory`, `region-prepare`, `region-trajectory`, `benchmark`, `static-suite`) |
| Python, one call | `radii=` on every analysis that uses radii, for example `pore_profile(frame, radii="hole")`, `rolling_probe_cast(frame, radii=my_set)`, `write_static_publication_bundle(..., radii="path/to/radii.json")` |
| Python, a block of code | `with crevice.radii.use_radii("hole"): ...` |

`radii` accepts a preset name, the path of a radius file, or a
{class}`crevice.radii.RadiusSet`. `None` (the Python default) means "the set
already in effect", which is the default set unless a caller chose another.
A preset name wins over a file of the same name; write `./hole` for a file
called `hole`.

### Presets

| Name | Values | Matching | Source |
|---|---|---|---|
| `default` | the standard table, plus CHARMM36 radii for recognised ion residues | element; ion residues by residue and atom name | Bondi 1964; Mantina et al. 2009; CHARMM36 `toppar_water_ions.str` |
| `bondi` | the standard table for every atom, ions included | element | Bondi 1964; Mantina et al. 2009 |
| `hole` | C 1.85, O 1.65, S 2.00, N 1.75, H 1.00, P 2.10 Å; GLN `E2?` and ASN `D2?` 1.00 Å; `LP??` 0.00 Å | first letter of the **atom name**, as in HOLE | HOLE 2.3.1 `rad/simple.rad` (Smart et al. 1996), whose header gives the AMBER united-atom radii of Weiner et al. (1984) |
| `charmm_like` | the standard table with C 2.00, N 1.85, O 1.70, S 2.00 Å, plus the default ion radii | element | a CREVICE approximation for C, N, O, S, **not** read from a CHARMM parameter file |

`bondi` differs from the default only for recognised ion residues, which it
treats as neutral atoms (a CHARMM `SOD` sodium is 2.27 Å instead of
1.41075 Å). This was the default of earlier CREVICE versions.

`hole` reproduces HOLE's default radius file, including HOLE's rules: an atom is
matched by its name, not its element (a CHARMM sodium ion named `SOD` gets the
`S???` radius, 2.00 Å; chloride `CLA` gets 1.85 Å), and an atom that no record
matches (for example zinc, `ZN`) is an error, as HOLE stops on it. Earlier
CREVICE versions had a `hole` preset that only changed hydrogen to 1.00 Å; that
preset was not HOLE's radius file and has been replaced.

### Radius files

The format is chosen by the file extension.

**JSON** (`.json`): an object with any of

```json
{
  "name": "bondi-ionic",
  "base": "bondi",
  "default_radius": 1.70,
  "element_radii": {"C": 1.70},
  "atom_radii": {"OW": 1.50},
  "ion_radii": {"NA": 1.02, "CL": 1.81},
  "citation": "where the values come from"
}
```

`base` starts from a preset and the file's entries override it; without a base,
only the listed radii exist and every other element gets `default_radius`
(1.70 Å unless given) and ions are ordinary atoms of their element.
`element_radii` keys must be element symbols; `atom_radii` keys are `ATOMNAME`
or `RESNAME:ATOMNAME`; `ion_radii` keys are element symbols and apply only to
recognised ion residues (`"base": "default"` keeps the CHARMM36 ion radii,
`"base": "bondi"` has none). Other keys are an error.
Radii must be finite and non-negative.

**CSV** (`.csv`): a header containing `kind,key,radius`, then one row per entry.

```text
# kind: element | atom | ion | default | base
kind,key,radius
base,default,
element,C,1.70
ion,NA,1.02
atom,OW,1.50
default,,1.70
```

**HOLE radius file** (`.rad`): only `VDWR` records are used, read exactly as
HOLE 2.3.1 reads them (`tsradr.f`): each line is upper-cased; a line starting
`VDWR` has the atom-name pattern in columns 6-9, the residue-name pattern in
columns 11-13 and the radius in columns 14-23. In the patterns `?` matches any
character, including a blank; other characters, including blanks, must match
exactly. The **first** matching record wins, so specific records go before
general ones, and an atom that no record matches is an error. `BOND` records,
`remark` lines and anything else are ignored.

```text
remark: specific records first
VDWR CB   ALA 2.00
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR ???? ??? 1.70
```

CREVICE refuses a `VDWR` line containing a tab, and a radius field without a
decimal point or left blank, because HOLE's fixed-format read would silently
change such a value (`2` would become 0.002 Å). For matching, CREVICE writes
the atom name left-justified in the four-character field and truncates the
residue name to three characters. HOLE instead reads PDB columns 14-17 as the
name (adding column 13 only when it holds `H`). The two agree for every
standard atom name of up to three characters and for four-character hydrogen
names; they differ for other four-character names, names with a leading digit
(`1HB2`) and two-letter element names written from column 13 (HOLE reads a
PDB zinc named `ZN` as `N`).

## Precedence

For each atom the first rule that applies gives the radius:

1. `RESNAME:ATOMNAME` in the set's atom radii;
2. `ATOMNAME` in the set's atom radii;
3. a recognised ion residue whose element is in the set's ion radii;
4. HOLE `VDWR` records, in file order (sets read from `.rad` files and the
   `hole` preset); with such records, no match is an error;
5. the element in the set's element radii;
6. the set's default radius.

Atom and residue names are compared after stripping blanks and upper-casing.

### Atoms without a radius

Rule 6 is a fallback, not a radius: when none of rules 1–5 applies to an atom
({meth}`crevice.radii.RadiusSet.has_radius` is false), the structure readers
exclude the atom and print an {class}`~crevice.radii.UnrecognisedElementWarning`
that lists the elements, atom counts and `RESNAME:ATOMNAME` keys and explains
how to add radii. Input reports record the exclusions
(`excluded_no_radius_atoms`, `excluded_no_radius_names`, `radius_set`). MD
trajectory loading refuses such a selection instead, because frames must stay
atom-for-atom identical to the MDAnalysis selection. Atoms you build directly
in Python, without a reader, still get the default radius.

To analyse such atoms, give them radii keyed by element, atom name or
`RESNAME:ATOMNAME` (most specific wins) in a radius file, for example from
{func}`crevice.radii.write_radius_template`:

```text
kind,key,radius
base,default,
element,GD,2.10
atom,XGD:GD1,2.30
```

and pass it as `--radii radii.csv` (Python `radii="radii.csv"`, both when
loading and when analysing). The values are your choice and are recorded in
every output's `radii` provenance (`element_overrides_A`, `atom_radii_A`,
`source_sha256`).

## How the choice reaches every calculation

All analyses obtain atomic radii from one function,
{func}`crevice.radii.atom_vdw_radius`. It reads the radius set in effect, which
is held in a {mod}`contextvars` variable. `--radii`, every `radii=` argument and
{func}`crevice.radii.use_radii` all set that variable for the duration of the
command, call or `with` block, and restore it afterwards, even after an error.
Nested analyses (a publish bundle's profile, cast, tunnels, network and
hydration) therefore use the same set.

- **Threads.** The variable is private to each thread. A thread you start
  yourself does not inherit the set: pass `radii=` to the analysis inside the
  thread or run it in `contextvars.copy_context()`.
- **Worker processes.** `cavity-trajectory` and `region-trajectory` with
  `--workers N` pass the set to each worker explicitly (as an initialiser
  argument) and enter it there, so forked and spawned workers measure with the
  same radii. A test checks that two workers give the same volumes as one.
- **Radii computed once.** Some objects compute their radii when they are
  created (for example the obstacle set of a region trajectory). Create them
  inside the same `radii=` call or `use_radii` block as the analysis, as the
  commands do.

## Provenance

Every output records the set it used as `radii`, including runs with the
default set:

- in the `metadata` of results that have one (profiles, casts, tunnels,
  networks, trajectory analyses), and hence in their JSON files;
- in every manifest: `publish`, `cast`, `analyze`, `static-suite` (per system and
  in the summary), `benchmark`, `residue-evidence`, `cavity-trajectory` and
  `region-trajectory` (cavity manifest and region provenance), hydration
  manifests and `region-prepare`;
- in the JSON outputs of `cavities`, `tunnels`, `residues --json`, `network`,
  `network --connectivity-json`, `features`, `trajectory` (and its
  `--distribution-csv` JSON and `--network-json`), and in the settings of
  `cavity-trajectory` statistics and trajectory hydration results;
- beside every CSV: each CSV is written next to a JSON output or manifest that
  carries `radii`; `residues` without `--json` writes
  `<stem>_provenance.json` for that purpose. CSV files themselves carry no
  provenance column.

The record holds the set's `name`, `source` (preset description or absolute
file path), `source_sha256` (the radius file's SHA-256), `citation`,
`table_sha256` (a hash of the effective table, so two runs can be compared),
`default_radius_A`, `element_overrides_A` (radii that differ from the standard
table), `standard_elements_missing` (standard elements the set lacks, which get
the default radius), `atom_radii_A`, `ion_radii_A`, `ion_scope` (which atoms
count as ions), `hole_records` and the `precedence` above.
{func}`crevice.radii.radius_set_from_provenance` rebuilds a set from this
record and {func}`crevice.radii.same_radius_set` compares two records by
`table_sha256`.

Outputs of earlier CREVICE versions recorded nothing for the default set. For
those, the absence of `radii` means the default *of that version*, which
treated ion residues as neutral atoms (and, in still earlier versions, read
CHARMM `SOD`/`POT` names as sulfur/phosphorus).

## Reusing earlier results

Steps that build on earlier outputs measure with the same radii:

- **Region definitions** store their radius set. `region-prepare` writes the
  prepared definition as schema version 2 with a `radii` record of the set it
  used (`--radii`, or the default). `region-trajectory` measures with that set.
  A schema version 1 definition (no `radii`) is read as the default set.
- **`hydration --cavity-results`** uses the set recorded in the cavity results
  (their statistics settings or region provenance).
- **`region-compare`** checks that both runs recorded the same set.

If you pass a different `--radii` to `region-prepare` (for a definition that
already declares a set), `region-trajectory` or `hydration --cavity-results`,
or compare runs with different sets, the command stops with an error naming
both sets. `--allow-radii-mismatch` proceeds anyway and records the outcome as
`radii_check` (`mismatch_overridden`). Results that recorded no set (earlier
versions) cannot be checked; the command proceeds and records
`unverified_no_recorded_set` (hydration) or `unverified_not_recorded`
(comparison).

## Interpreting results

Radius sets change results by design. On a synthetic carbon cylinder the
measured radius changes by exactly the change in the carbon radius. On real
structures the change is not uniform, because different atoms change by
different amounts. No radius set is "correct" in general: HOLE's `simple.rad`
uses united-atom radii that are larger than Bondi's to stand in for missing
hydrogens, while Bondi radii describe neutral atoms. Neutral-atom radii are
large for monatomic ions (Bondi Na 2.27 Å, K 2.75 Å), which is why the default
uses CHARMM36 sizes for ion residues; crystal ionic radii are smaller still. If
ions act as obstacles in your analysis, report which radius they had, and
consider a sensitivity run with other ion radii (`"ion_radii"` in a radius
file) or with the ions removed from the obstacle selection. Report the radius
set with any comparison between tools.

## Limitations

- A region's fixed reference map is computed once, at `region-prepare`.
  Reading a version 1 definition as the default set does not recompute that
  map; only the per-frame measurements use the current radii.
- Ions are recognised by name only (the ion table); an ion with another
  residue or atom name is an ordinary atom of its inferred element.
- `vdw_radius(element)` and `VDW_RADII` always describe the element table;
  they know nothing about ions and are not affected by the set in effect.

## Python reference

{class}`crevice.radii.RadiusSet`, {func}`crevice.radii.radii_preset`,
{func}`crevice.radii.custom_radii_set`, {func}`crevice.radii.load_radius_file`,
{func}`crevice.radii.read_hole_radius_file`, {func}`crevice.radii.resolve_radii`,
{func}`crevice.radii.use_radii`, {func}`crevice.radii.active_radii`,
{func}`crevice.radii.radii_provenance`, {func}`crevice.radii.atom_vdw_radius`,
{data}`crevice.radii.DEFAULT_RADII`, {data}`crevice.radii.ION_RADII`,
{func}`crevice.radii.is_recognised_ion`, {func}`crevice.radii.effective_radii`,
{func}`crevice.radii.radius_set_from_provenance`,
{func}`crevice.radii.same_radius_set`,
{func}`crevice.radii.reconcile_recorded_radii`.

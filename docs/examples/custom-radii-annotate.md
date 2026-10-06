# Atomic radii and figure text: `--radii` and `--annotate`

Every distance CREVICE reports is measured to atom van der Waals surfaces, so
the radius set is part of the result. This example measures the 1GRM channel
with three radius sets, then draws the same figure with and without its
explanatory text. It takes under a minute. See
[Atomic radii](../methods/atomic-radii.md) for the presets and file formats.

The enclosure probe is fixed at 0.75 Å (the value the automatic rule chooses
for this structure, see the [1GRM example](1grm-profile-cast.md)) so that only
the radii change between runs.

```bash
crevice fetch 1GRM --cache-dir .crevice/pdb

# Standard table (Bondi 1964 element radii; the default)
crevice profile 1GRM --axis auto --enclosure-radius 0.75 -o standard.json --png standard.png --dpi 150

# HOLE's simple.rad radii (AMBER united-atom values, matched by atom name)
crevice profile 1GRM --axis auto --enclosure-radius 0.75 --radii hole -o hole.json

# Your own file: the Bondi set with a larger oxygen radius
crevice profile 1GRM --axis auto --enclosure-radius 0.75 --radii my_radii.json -o custom.json
```

with `my_radii.json`:

```json
{
  "name": "bondi-larger-oxygen",
  "base": "bondi",
  "element_radii": {"O": 1.60},
  "citation": "worked example: O 1.52 -> 1.60 A to show the effect of one radius"
}
```

| Radius set | Minimum radius (Å) | Nearest residue at the bottleneck |
|---|---:|---|
| standard (Bondi) | 1.324 | A:TRP13 |
| `hole` (C 1.85, O 1.65, N 1.75 Å, ...) | 1.180 | A:TRP13 |
| `my_radii.json` (O 1.60 Å) | 1.282 | A:TRP11 |

Larger atoms leave a narrower channel, and with the larger oxygen radius a
different residue limits the narrowest sample. The difference between the
standard and `hole` sets (0.14 Å) is larger than many effects one might want
to interpret, so always state the radius set.

When you choose a set, the output records it under `metadata.radii`: its
`name`, `source` (the file path), `source_sha256`, `citation`,
`element_overrides_A` (here `{"O": 1.6}`), `table_sha256` (a hash of the whole
effective table, so two runs can be compared) and the matching precedence. The
same `--radii` option exists on every command that measures geometry, and one
set applies to the whole command, including worker processes. In Python:

```python
from crevice import load_structure, pore_profile

frame = load_structure(".crevice/pdb/1GRM.cif")
profile = pore_profile(frame, axis="auto", enclosure_radius=0.75, radii="my_radii.json")
print(round(profile.min_radius, 3), profile.metadata["radii"]["name"])   # 1.282 bondi-larger-oxygen
```

## Clean figures and `--annotate`

Figures are clean by default: data, axis labels with units, tick labels and a
legend where colours or line styles distinguish several series. Titles, value
call-outs and labels are drawn only with `--annotate`:

```bash
crevice profile 1GRM --axis auto --enclosure-radius 0.75 --annotate -o annotated.json --png annotated.png --dpi 150
```

```{figure} images/1grm_profile_annotated.png
:alt: 1GRM radius profile with title and bottleneck value
:width: 85%

`annotated.png`: the figure of `standard.png` (shown clean in the
[1GRM example](1grm-profile-cast.md)) with the title "CREVICE pore radius
profile" and the bottleneck value "1.32 Å". The two mouth labels with their
axial positions (−11.67 Å and 11.60 Å) are drawn only when the mouth guides
are also requested (`--annotate --mouth-guides`). The data, colours and axes are the
same; the measurement files are identical.
```

Nothing is lost in a clean figure: the text is stored in the PNG metadata.

```python
from crevice.presentation import figure_description

info = figure_description("standard.png")
print(info["CREVICE annotations"])
print(info["Description"])
```

```text
omitted (recorded here; draw with annotate=True or --annotate)
Pore radius against axial coordinate for one resolved channel profile; minimum radius 1.324 Å at the bottleneck (sample 73, t = 9.56 Å, nearest residue A:TRP13); mean radius 1.563 Å over 81 samples
Label: 1.32 Å
Note: Lower mouth -11.67 Å
Note: Upper mouth 11.60 Å
Title: CREVICE pore radius profile
```

`--annotate` (Python: `annotate=True`, or
{func}`crevice.presentation.figure_annotations`) works on every command that
writes figures or viewer scenes; in the scenes it adds the colour key, captions
and residue labels.

## What was checked

All commands and the Python snippets were run with CREVICE 0.1.0 and the
values above are copied from that run. They show how settings change a
measurement; none of them is a validated biological value.

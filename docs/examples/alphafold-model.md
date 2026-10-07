# An AlphaFold DB model: fetch by identifier and cast

CREVICE accepts AlphaFold DB identifiers wherever it accepts a structure.
This example fetches the AlphaFold model of *E. coli* lactose permease (LacY,
UniProt P02920) and fills its interior with the default rolling-probe cast.
It takes under a minute, but it needs network access to `alphafold.ebi.ac.uk`,
so, like every tutorial, it isn't run when the documentation is built.

```bash
crevice fetch AF-P02920-F1 --format pdb --cache-dir .crevice/pdb
crevice cast AF-P02920-F1 --format pdb --out-dir cast
```

Printed output:

```text
.crevice/pdb/AF-P02920-F1-model_v6.pdb	download	273212 bytes	verified
cast regions=1 volume_A3=3276.500 spacing_A=0.5 out_dir=cast
```

`AF-P02920-F1` without a version fetched the latest model (version 6 at the
time of this run; the AlphaFold DB API is asked for the latest version). The
file is cached as `AF-P02920-F1-model_v6.pdb` with a provenance sidecar, and
"verified" means that its `DBREF` record names UniProt P02920. You can also
write `AF-P02920-F1-model_v6` or `AF-P02920-F1-model_v6.pdb` to pin a
version; a file of that name on disk is read like any other PDB file.

The input report (`cast/AF-P02920-F1-model_v6_input_report.json`) marks the
structure as a prediction:

```json
"model_type": "predicted",
"predictor": "AlphaFold",
"bfactor_meaning": "pLDDT per-residue confidence (0-100), not a crystallographic B-factor; treat low-confidence regions with caution"
```

`crevice profile AF-P02920-F1` finds no channel with two open ends in this
model, so the cast is the measurement to use here.

:::{admonition} Interpreting the results
:class: crevice-interpret

The 3276.5 Å³ cast is the interior of one predicted conformation at the
default settings (0.5 Å grid, 0.8 Å probe, main region only). It is not an
experimentally determined cavity, and regions with low pLDDT shouldn't be
read as structure.
:::

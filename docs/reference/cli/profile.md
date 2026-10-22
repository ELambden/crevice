# crevice profile

Radius profile of a channel that runs through the protein; `--fingerprint` adds lining chemistry, water and a geometric conductance estimate along it.

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice profile`

Measure the radius profile of a channel that runs through the protein. For the method, worked commands and outputs, see
[Pore radius profiles](../../tools/profile.md).

```{crevice-cli} profile
```

## `crevice profile --fingerprint`

*CREVICE 0.3: `crevice fingerprint`.*

Pore fingerprint along a channel: radius, lining hydrophobicity and chemistry, ENTRY/LUMEN/EXIT segments, water
along the path (MD) and a geometric conductance estimate. For the definitions, worked commands and outputs, see
[Pore fingerprints](../../tools/fingerprint.md).

```{crevice-cli} profile --fingerprint
```

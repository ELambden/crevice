# crevice residues

Residues around a channel or a measured volume: lining residues (default), contact networks (`--network`), the residues that form a measured boundary (`--evidence`) and numeric descriptors (`--features`).

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice residues`

Annotate channel-lining residues. For the method, worked commands and outputs, see
[Lining residues](../../tools/residues.md).

```{crevice-cli} residues
```

## `crevice residues --network`

*CREVICE 0.3: `crevice network`.*

Build residue and channel/cavity interaction networks. For the method, worked commands and outputs, see
[Residue networks](../../tools/network.md).

```{crevice-cli} residues --network
```

## `crevice residues --evidence`

*CREVICE 0.3: `crevice residue-evidence`.*

Explain a measured cast boundary and nonlocal residue partners. For the method, worked commands and outputs, see
[Boundary residues and their partners](../../tools/residue-evidence.md).

```{crevice-cli} residues --evidence
```

## `crevice residues --features`

*CREVICE 0.3: `crevice features`.*

Export profile and residue feature JSON. For the method, worked commands and outputs, see
[Utility and batch commands](../../tools/utilities.md).

```{crevice-cli} residues --features
```

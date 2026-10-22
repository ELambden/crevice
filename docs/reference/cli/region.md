# crevice region

The named-region workflow: `init` writes a candidate definition, `prepare` checks it and writes review scenes, `compare` compares two regions on the same frames. A prepared region is measured with `crevice trajectory --region DEFINITION`.

Each mode runs its own analysis with its own options and defaults; the first line of each section shows how to select it. The CREVICE 0.3 command name of each mode still works as a deprecated alias (see [Upgrading to 0.4](../../getting-started/upgrading-0.4.md)).

## `crevice region init`

*CREVICE 0.3: `crevice region-init`.*

Create a candidate region definition with hashed inputs and explicit landmarks. For the method, worked commands and outputs, see
[Named regions](../../tools/regions.md).

```{crevice-cli} region init
```

## `crevice region prepare`

*CREVICE 0.3: `crevice region-prepare`.*

Validate a named region, bind its reference frame and write native review overlays. For the method, worked commands and outputs, see
[Named regions](../../tools/regions.md).

```{crevice-cli} region prepare
```

## `crevice region compare`

*CREVICE 0.3: `crevice region-compare`.*

Compare separate named regions on identical trajectory frames. For the method, worked commands and outputs, see
[Named regions](../../tools/regions.md).

```{crevice-cli} region compare
```

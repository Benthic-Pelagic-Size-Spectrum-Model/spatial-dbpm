# Adapter prototypes (FishMIP / LME ↔ dbpmr)

Work-in-progress glue connecting the LME ISIMIP3a forcing to the `dbpmr` engine,
per the [design docs](../design). This is **not** part of the `dbpmr` package
(it lives outside the package dir so it doesn't affect the build).

## `stage0_prototype.R`

End-to-end Stage 0 / Stage 2 harness for one LME:

1. reads the per-LME forcing parquet (`stable-spin` constants) and the reference
   no-fishing equilibrium JSON;
2. maps the LME `sizemodel()` parameters to `dbpmr` `Setup.*` arguments via the
   crosswalk in [`design/stage0-reference.md`](../design/stage0-reference.md)
   (log10→ln grid, `lambda = slope`, `u_0 = 10^intercept`, temperature folded
   into the search constant `A`);
3. runs a matched no-fishing aspatial `dbpmr` simulation to equilibrium;
4. compares the equilibrium pelagic spectrum to the reference and writes
   `stage0_compare.png`.

### Requirements

- An installed `dbpmr`, plus `arrow` and `jsonlite`.
- The local LME data (not in this repo). Point `DBPM_DATA` at the folder holding
  `dbpm_inputs/` and `equilibrium_runs/`:
  ```sh
  DBPM_DATA=/path/to/DBPM_dev Rscript adapter/stage0_prototype.R
  ```

### Status / finding

The pipeline runs end-to-end, but the **literal crosswalk does not reproduce the
LME equilibrium**: dbpmr's pelagic spectrum decays toward extinction across the
search-rate range (final biomass ~1e-7 vs a non-trivial reference) and becomes
numerically unstable when `A` is pushed high. So the divergence is **deeper than
a single search-volume scale factor** — the resource/plankton coupling, the
reproduction boundary condition, and the assimilation/mortality balance need
reconciling between the two engines. Tracked in issue #8 (Stage 0).

This is the expected Stage 0 outcome: the harness now *quantifies* the gap, which
defines the reconciliation work before the temperature/fishing stages.

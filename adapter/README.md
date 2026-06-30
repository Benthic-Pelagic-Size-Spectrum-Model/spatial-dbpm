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

### Status / findings

The pipeline runs end-to-end. Progress so far:

- **Stable non-trivial equilibrium achieved** with `A_pel = 64`, `A_ben = 6.4`
  (0.1×64) and **weekly** steps (pelagic biomass ~0.56). Earlier configs
  collapsed to extinction or went `NaN`.
- **dbpmr is numerically unstable at monthly steps** for these LME-scale params,
  but **stable at weekly/daily** — and weekly is the FishMIP time step the LME
  workflow uses. (Relevant to issue #7.)
- **Plankton confirmed fixed** through the run, matching `sizemodel()`.

**Remaining (issue #8):** the equilibrium *shape* and *scale* don't yet match —
dbpmr's pelagic spectrum is much shallower and offset by many orders of magnitude
from the reference. Next: reconcile the density normalisation/units and the
spectral slope (growth-vs-mortality balance, recruitment boundary) by comparing
the C rate equations to `sizemodel()`.

# Sandbox: dbpmr as the FishMIP engine (build space)

A **sandbox** for wiring dbpmr into the LME/FishMIP workflow **without touching the
production repo**. It follows the plan in
[`design/dbpmr-fishmip-integration-plan.md`](../../design/dbpmr-fishmip-integration-plan.md)
and epic [#26](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/spatial-dbpm/issues/26).

**Principle: import, don't copy.** The LME pipeline is a **git submodule**
(`lme-workflow/`, pinned to `new_features`), sourced as-is — so forcing,
`sizeparam`, `run_model`, calibration and gravity are the real workflow, unchanged.
We override **only the integrator** via an `engine=` switch. When Tier-1 validates,
the plan is to fold that switch back into the LME repo as a PR — no permanent fork.

## Files
- `lme-workflow/` — the ISIMIP3a workflow (submodule). Run `git submodule update --init`.
- `dbpmr_engine.R` — **`run_model_dbpmr()`**: the dbpmr engine adapter. Same call
  signature as the workflow's `run_model()`; derives parameters via the workflow's
  own `sizeparam()` (identical biology), runs the dbpmr C engine, returns biomass +
  spectra in a `run_model`-compatible shape.
- `tier1_engine_swap.R` — Tier-1 driver: runs one FAO-LME through **both** engines
  via a single `engine = "sizemodel" | "dbpmr"` switch and reports biomass.

## Run
```sh
git submodule update --init --recursive
DBPM_DATA=/path/to/DBPM_dev Rscript adapter/sandbox/tier1_engine_swap.R 14
```
Requires an installed `dbpmr` (the script points `.libPaths` at `/tmp/dbpmrlib`)
and the submodule's R deps (arrow, dplyr, …).

## Status
- [x] Submodule import + `engine=` switch runs both engines on one LME (aspatial,
      stable-spin, no fishing). **Proven** (LME-14).
- [ ] Reconcile magnitudes: apply ln(10) fix (PR #14) to the submodule + the
      ln-vs-log10 density convention, so the two engines match like-for-like.
- [ ] Fishing (`Fmort` / gravity) + catch output in the adapter.
- [ ] Optimiser calibration of `search_vol` (#17) through the swapped engine.
- [ ] Transient obsclim forcing (needs the in-memory column driver, #5/#11).
- [ ] Fold `engine=` into the LME repo as a PR (#24).

## Scope
Tier-1 only (aspatial calibration seam). Gridded runs + gravity feedback are
Tier-3 and need the Python binding (#24) on the in-memory column driver (#5).

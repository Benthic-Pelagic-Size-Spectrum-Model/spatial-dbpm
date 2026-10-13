# Sandbox: dbpmr as the FishMIP engine (build space)

A **sandbox** for wiring dbpmr into the LME/FishMIP workflow **without touching the
production repo**. It follows the plan in
[`design/dbpmr-fishmip-integration-plan.md`](../../design/dbpmr-fishmip-integration-plan.md)
and epic [#26](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/spatial-dbpm/issues/26).

**Principle: import, don't copy.** The LME pipeline is the **separate
[`dbpm_isimip3a`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_isimip3a)
repo** — a standalone workflow, kept deliberately distinct from this package build
(it was previously embedded here as a git submodule; now decoupled). Clone it into
`lme-workflow/` to run the sandbox scripts, so forcing, `sizeparam`, `run_model`,
calibration and gravity are the real workflow, unchanged. We override **only the
integrator** via an `engine=` switch.

## Files
- `lme-workflow/` — the ISIMIP3a workflow, now the standalone
  [`dbpm_isimip3a`](https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_isimip3a)
  repo. Clone it here (not tracked by this repo):
  `git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_isimip3a.git adapter/sandbox/lme-workflow`.
- `dbpmr_engine.R` — **`run_model_dbpmr()`**: the dbpmr engine adapter. Same call
  signature as the workflow's `run_model()`; derives parameters via the workflow's
  own `sizeparam()` (identical biology), runs the dbpmr C engine, returns biomass +
  spectra in a `run_model`-compatible shape.
- `tier1_engine_swap.R` — Tier-1 driver: runs one FAO-LME through **both** engines
  via a single `engine = "sizemodel" | "dbpmr"` switch and reports biomass.

## Run
```sh
# clone the workflow repo into lme-workflow/ (one-time)
git clone https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_isimip3a.git adapter/sandbox/lme-workflow
DBPM_DATA=/path/to/DBPM_dev Rscript adapter/sandbox/tier1_engine_swap.R 14
```
Requires an installed `dbpmr` (the script points `.libPaths` at `/tmp/dbpmrlib`)
and the workflow's R deps (arrow, dplyr, …).

## Status
- [x] Submodule import + `engine=` switch runs both engines on one LME (aspatial,
      stable-spin, no fishing). **Proven** (LME-14). — `tier1_engine_swap.R`
- [x] **Fishing + single-Q transient calibration** — `F=Q·s(x)·effort_norm(t)`,
      knife-edge at 10 g, effort normalised [0,1], catch-time-series log-MSE
      objective, `Q∈[0,3]`, `A=64` fixed. **Proven** (LME-14: `optimise` 13 evals,
      `Q=0.018`, catch-time-series corr 0.95). — `tier1_fishing_calib.R`
- [ ] Time-varying **environmental** forcing (temperature/plankton) — currently
      held constant; needs the in-memory column driver (#5/#11).
- [ ] Reconcile magnitudes vs sizemodel (apply ln(10) PR #14 to submodule +
      density convention) — for like-for-like comparison (not required to run).
- [ ] Gravity effort split across groups (∝ biomass) + creep + sectors (DBPM.md).
- [ ] Fold `engine=` / the calibration into the LME repo as a PR (#24).

## Scope
Tier-1 only (aspatial calibration seam). Gridded runs + gravity feedback are
Tier-3 and need the Python binding (#24) on the in-memory column driver (#5).

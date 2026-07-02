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
      stable-spin, no fishing). **Proven** (LME-14). — `tier1_engine_swap.R`
- [x] **Fishing + single-Q transient calibration** — `F_g=Q·s(x)·(B_g/ΣB)·effort_norm(t)`,
      knife-edge at 10 g, effort normalised [0,1], catch-time-series log-MSE
      objective, `Q∈[0,3]`, `A=64` fixed. — `tier1_fishing_calib.R`
- [x] **Gravity effort split across groups** (∝ fishable biomass, DBPM.md):
      effort split pel/ben by share, `F_g=Q·(B_g/ΣB)·s(x)·effort`.
      **Proven** (LME-14: shares 0.93/0.07). — `tier1_fishing_calib.R`
- [x] **Environmental forcing PLUGIN** (#11) — NEW dbpmr C-engine forcing driver:
      auto-detected `<run>/Input/forcing_ts.txt` with a **header-named channel** line
      + one row per timestep, read each step. Channels (extensible — add a column):
      `pel_tempeff`/`ben_tempeff` scale feeding `A` + background mortality `mu_0`
      from stored base values (as sizemodel's `pel/ben_tempeffect`; senescence/fishing
      unscaled); `sinking_rate` scales surface-origin detritus inputs in `g_det`
      (detritus mortality kept un-temperatured, per sizemodel); `depth` enables the
      Dunne-2007 burial loss (opt-in — absent ⇒ original detritus dynamics). Plankton
      stays on the `ts_flag` spectrum input. All fed from obsclim (`tos`/`tob`,
      `export_ratio`, `intercept`/`slope`) + `depth`→`pref_benthos=0.8·exp(-depth/250)`.
      Verified: constant temp channel reproduces folded `A=S·pt`; `sinking_rate<1` and
      burial lower benthos; unit tests in `test-temperature.R`. **Proven** (LME-14:
      `Q=0.095`, catch corr 0.95). — `dbpmr/src/SizeSpectra.c`, `tier1_fishing_calib.R`
- [x] **Full detritus match to sizemodel** — `g_det` now adds dead benthos (background +
      senescence mortality, NO sinking — already on the bed) alongside the sinking
      surface inputs, and the pool loses the Dunne-2007 burial fraction. Both engines
      now use the same detritus forcings. — `dbpmr/src/SizeSpectra.c`
- [x] **Engine adapter uses the same forcings** — `dbpmr_engine.R` sets
      `pref_benthos` from depth and writes a constant `sinking_rate`/`depth`
      `forcing_ts.txt` (stable-spin), so the engine-swap comparison is like-for-like on
      the environment. — `dbpmr_engine.R`
- [x] **ln-vs-log10 units audited/fixed** — (a) plankton intercept is per-log10;
      dbpmr wants per-ln, so `u_0=10^intercept/ln10` (was 2.3× too high);
      (b) engine biomass integral had a spurious ×ln(10). Verified against the C
      engine's own integrals (SizeSpectra.c:2263/2293/2336). — `dbpmr_engine.R`
- [ ] Reconcile magnitudes vs sizemodel (apply ln(10) PR #14 to submodule +
      density convention) — for like-for-like comparison (not required to run).
- [ ] Gravity refinements: creep (annual Q increase) + sectors (DBPM.md).
- [ ] Fold `engine=` / the calibration into the LME repo as a PR (#24).

## Scope
Tier-1 only (aspatial calibration seam). Gridded runs + gravity feedback are
Tier-3 and need the Python binding (#24) on the in-memory column driver (#5).

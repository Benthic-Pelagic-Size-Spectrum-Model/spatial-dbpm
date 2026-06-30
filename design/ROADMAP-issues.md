# Proposed issues / roadmap checklist

A curated backlog distilled from the design docs
([integration](dbpmr-lme-integration.md), [input schema](lme-input-schema.md),
[Stage 0 reference](stage0-reference.md)) and known bugs. Intended to be copied
into GitHub issues (or worked through directly). Suggested labels in brackets.

## Bugs

- [ ] **`Average.Time()` references non-existent slots** `[bug]`
  Uses `species@grid@tout` (no such slot — should be `toutstep`), and treats
  `@run`/`@spatial.dim`/`@data` as if a data.frame; errors at runtime.
  File: `dbpmr/R/Average.Time.r`.

- [ ] **`Read.In()` / `Setup.Pelagic()` preference & gamma typos** `[bug]`
  `Read.In` swaps `pref_pla`/`pref_pel` and sets `gamma_prey <- gamma_pred`;
  `Setup.Pelagic` likewise sets `species@gamma_prey <- gamma_pred`.
  Files: `dbpmr/R/Read.In.r`, `dbpmr/R/Setup.Pelagic.r`.

## Engine / performance

- [ ] **Replace text-file I/O with an in-memory `.Call` return** `[enhancement][performance]`
  Densities/growth/mortality/budgets currently round-trip through CSV under
  `<run>/<species>/`. This is the #1 bottleneck for gridded (thousands-of-cells)
  FishMIP runs. Return arrays directly from C instead.

- [ ] **Finish `.C` → `.Call` modernization** `[enhancement]`
  Native routine is registered; move the interface to `.Call` so forcing and
  outputs can pass in memory (enables the item above).

- [ ] **Monthly time-step support + validation (`tstep = 1/12`)** `[enhancement]`
  Confirm/validate stability of the implicit scheme at monthly steps (the LME
  engine runs monthly, `timesteps_years = 0.0833`); document with `diff_method`.

## FishMIP integration roadmap

- [ ] **Stage 0 — reproducibility cross-check vs LME no-fishing equilibrium** `[integration]`
  Configure a matched no-fishing, constant-forcing `dbpmr` 0-D run per the
  parameter crosswalk; compare equilibrium pelagic/benthic spectra to the
  `equilibrium_runs/*.json` reference. Needs no engine change (temperature folds
  into `A`). See `design/stage0-reference.md`.

- [ ] **Stage 1 — benchmark current performance** `[performance]`
  Time/memory a single 0-D `dbpmr` run, its text-I/O cost, and a parquet read;
  establish the per-cell baseline that gridded scaling multiplies.

- [ ] **Stage 2 — R adapter (forcing → `dbpmr`)** `[integration]`
  Pure-R adapter reading the per-LME parquet forcing + reference JSON params and
  building `Setup.*` objects via the crosswalk; run `dbpmr`; collect outputs.
  No engine change. See `design/lme-input-schema.md`.

- [ ] **Stage 3 — environmental forcing input** `[enhancement]`
  Add a forcing input carrying, per time step: `tos`/`tob` (Boltzmann–Arrhenius
  scaling of pelagic/benthic rates), time-varying plankton `intercept`/`slope`
  (`lambda = slope`, `u_0 = 10^intercept`), and detritus export. Extend the
  `.C`/`.Call` parameters accordingly.

- [ ] **Stage 4 — effort / gravity fishing** `[integration]`
  Drive `dbpmr`'s fishing input from the LME `gravitymodel()` output; support
  size-based selectivity per functional group (`min/max_fished_weight_class`).

- [ ] **Stage 5 — MPA / accessibility masking** `[enhancement]`
  Generalise the sea-ice `simask` concept to MPA masks applied in the adapter
  (zeroing effort/accessibility per cell-time). No engine change.

- [ ] **Stage 6 — additional functional groups** `[enhancement]`
  Extend the adapter + parameter mapping to FishMIP functional groups
  (`dbpmr` already supports multiple pelagic/benthic species).

- [ ] **Stage 7 — FishMIP-compatible outputs** `[integration]`
  Adapter writes FishMIP netCDF (reuse LME writers) from `dbpmr` results;
  density ↔ tonnes via `area_m2`.

- [ ] **Stage 8 — performance optimisation (after correctness)** `[performance]`
  Profile-guided only; first target is the in-memory I/O item above, then
  consider batching cells into one C call. Do not port science into C for speed.

## Calibration

- [ ] **Wrap calibration around `dbpmr`** `[integration]`
  Keep the LME `getError()` / `LHSsearch()` optimisation in the workflow,
  calling `dbpmr` (via the adapter) as the engine; target = `catch_tonnes`
  (confirm vs `_pauly`).

## Open scientific confirmations (small)

- [ ] Confirm `dynamic_reproduction = 1` ↔ `dbpmr` `rep_method` (2 vs 3).
- [ ] Confirm net pelagic growth-efficiency mapping (`defecate_prop` +
  `growth_pred` + `energy_pred` → single `epsilon`).
- [ ] Confirm target model time step for FishMIP runs (daily vs weekly vs monthly).

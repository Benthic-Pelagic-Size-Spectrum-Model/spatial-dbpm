# LME forcing input schema and adapter mapping

Companion to [`dbpmr-lme-integration.md`](dbpmr-lme-integration.md). Documents the
real FishMIP/ISIMIP3a 0-D forcing files so the **adapter layer** (Stage 2) and the
**`dbpmr` forcing input object** (Stage 3) can be designed against actual data,
not assumptions. **Design only — no code.**

## 1. File inventory

- Location (local): `…/DBPM_mizer/DBPM_dev/dbpm_inputs/`
- **83 files**, one per FAO Large Marine Ecosystem:
  `dbpm_clim-fish-inputs_fao_lme-<N>_1641-2010.parquet`
- Each file: a single LME-aggregated (non-spatial, 0-D) time series.
- **Monthly** resolution, **4440 rows** = 370 years × 12 months (1641–2010).
- Columns: 28. One row = one LME-month.

### Scenario / spin-up structure (the `scenario` column)

| scenario | years | months | role |
|---|---|---|---|
| `stable-spin` | 1641–1840 | 2400 | constant climatology spin-up (equilibrate the spectrum) |
| `spinup` | 1841–1960 | 1440 | transient spin-up toward the historical period |
| `obsclim` | 1961–2010 | 600 | ISIMIP3a observed-climate forcing (the analysis period) |

This maps cleanly onto a `dbpmr` run: integrate the whole series, take the
`stable-spin` block to reach equilibrium initial conditions, then `spinup` +
`obsclim` as the forced transient. Catch (calibration target) is only present in
`obsclim`.

## 2. Column dictionary

| Column | Type | Model role | Notes / units |
|---|---|---|---|
| `scenario` | str | run phase | stable-spin / spinup / obsclim |
| `region`, `region_name` | str | id | e.g. "LME 1", "East Bering Sea" |
| `time`, `year`, `month` | date/num/str | time axis | monthly |
| `area_m2` | dbl | LME area | for density ↔ total conversions |
| **`intercept`** | dbl | **plankton resource spectrum intercept** | **log10** space; range ≈ −2.4…+0.1 |
| **`slope`** | dbl | **plankton resource spectrum slope** | **log10** space; ≈ −1.0 |
| **`tos`** | dbl | **sea-surface temperature** → pelagic rates | °C, ≈ 0.2…9.8 |
| **`tob`** | dbl | **sea-floor temperature** → benthic rates | °C, ≈ 0.4…4.1 |
| `export_ratio` | dbl | detritus export ratio | ≈ 0.11…0.26 |
| `expc_bot` | dbl | POC flux at bottom (detritus input) | very small |
| `input_w` | dbl | detritus/sinking input weight | ≈ 2…80 |
| `depth` | dbl | LME mean depth | constant per LME (capped at 200 m in LME runs) |
| `simask` | dbl | sea-ice / accessibility mask | 1 = accessible |
| `total_nom_active`, `total_nom_active_area_m2` | dbl | nominal fishing effort (absolute) | |
| `nom_active_relative`, `nom_active_area_m2_relative` | dbl | **relative effort (0–1)** | ramps up through obsclim |
| `catch_tonnes`, `catch_tonnes_area_m2` | dbl | **observed catch (calibration target)** | NA in spin-up |
| `catch_tonnes_pauly`, `catch_pauly_tonnes_area_m2` | dbl | alt. catch reconstruction | |
| `min_catch_density`, `max_catch_density` | dbl | catch-density bounds | |
| `min_fished_weight_class`, `max_fished_weight_class` | dbl | **fishing selectivity size range** | **log10 g**; 1.70 ≈ 50 g, 5.9–6.3 ≈ 0.8–2 t |

## 3. Mapping to `dbpmr` (keep / wrap / force)

| LME column(s) | Becomes in `dbpmr` | Mechanism | Stage |
|---|---|---|---|
| `intercept`, `slope` | time-varying **plankton resource** (replaces static `plankton.params` spectrum) | new **forcing input** (Stage 3); log10→ln conversion | 3 |
| `tos` | temperature scaling of **pelagic** rates | new forcing + Boltzmann–Arrhenius in core | 3 |
| `tob` | temperature scaling of **benthic** rates | same | 3 |
| `export_ratio`, `expc_bot`, `input_w` | **detritus input** forcing (replaces static `detritus.params@w_0`) | forcing input | 3 |
| `depth` | benthic habitat depth / coupling | run/grid parameter; cap at 200 m as in LME | 2–3 |
| `simask` | accessibility mask on effort | adapter mask (generalises to MPA, Stage 5) | 5 |
| `nom_active_relative` + `min/max_fished_weight_class` | size-selective **fishing mortality** | effort × catchability × selectivity → existing fishing input | 4 |
| `catch_tonnes` (obsclim) | **calibration target** | stays in LME calibration loop (`getError`/`LHSsearch`) | 7 |
| `area_m2` | density ↔ tonnes conversion | adapter / output formatting | 7 |

## 4. Implications for the `dbpmr` forcing input object (Stage 3 spec)

A new forcing input must carry, **per model time step**, at least:

```
time, tos, tob, plankton_intercept, plankton_slope,
export_ratio (+ expc_bot / input_w), effort_relative, accessibility_mask
```

Design notes:

1. **Time-base mismatch.** Forcing is **monthly**; `dbpmr`'s default `tstep` is
   daily (the LME repo also runs **weekly**). The adapter must expand/interpolate
   monthly forcing to the model step (step-wise hold vs linear — match the LME
   convention for continuity). Forcing should be supplied as a per-`toutstep`/
   per-`tstep` series, analogous to the existing `_ts.txt` inputs.
2. **Log base.** `intercept`/`slope` and the fished weight classes are **log10**;
   `dbpmr`'s core works in **natural log** mass. Conversion (and the exact
   plankton-spectrum parameterisation) is the single most important
   scientific-continuity detail — must be validated in the Stage 0 cross-check.
3. **Surface vs bottom temperature** drive pelagic vs benthic groups separately
   (as in the LME `sizemodel(tos, tob)`).
4. **Spin-up** is encoded in `scenario`; the adapter should use `stable-spin` to
   derive equilibrium initial conditions, then force through `spinup`+`obsclim`.

## 5. What this unlocks now (no code)

- **Stage 0 (reference + harness):** these files *are* the reference inputs for a
  per-LME 0-D run. A reference **output** (a saved `sizemodel()`/gridded result)
  or the LME R helpers (`useful_functions.R`) would let us run the reference
  engine and define the agreement metric. **Still needed from the user.**
- **Stage 1 (benchmark):** pick one LME, time a matched `dbpmr` 0-D run and its
  text-I/O cost vs the cost of reading one parquet.
- **Stage 2 (adapter):** the column→`dbpmr` mapping above is the adapter's
  contract; it can be prototyped in R reading parquet via `arrow` (already
  installed) without touching the engine.

## 6. Open questions for the user

1. Is there a **saved reference run output** (per LME) to compare against, or
   should we run the LME `sizemodel()` ourselves (need `useful_functions.R`)?
2. Confirm the **plankton spectrum parameterisation**: how `intercept`/`slope`
   (log10) define the resource that the model's smallest sizes feed on.
3. Target **model time step** for FishMIP runs: daily or weekly?
4. Which **catch series** is the calibration target — `catch_tonnes` or the
   `_pauly` reconstruction?

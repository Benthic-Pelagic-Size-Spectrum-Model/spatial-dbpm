# dbpmr × FishMIP integration — the efficient/optimal plan

How to combine **dbpmr's engine** (correct, guarded, tested, fast C) with the
**LME/FishMIP Python workflow** (GFDL-MOM6-COBALT2 forcing, calibration, gravity
fishing) so that dbpmr *replaces* `sizemodel()` and `gridded_sizemodel_rk4()` while
linking to `gravitymodel()`. Supersedes/extends
[`dbpmr-lme-integration.md`](dbpmr-lme-integration.md) and
[`lme-input-schema.md`](lme-input-schema.md), now concrete because we know the
units, the `ln(10)` bug, the senescence choice, and the calibration objective.

Core premise: the two are the **same DBPM equations**, so ~90% of the FishMIP
pipeline is *not* the integrator — it's forcing, calibration, fishing and I/O.
Keep all of that; swap only the ~200 lines that step the PDE.

## Layered architecture (keep vs replace)

| Layer | Keep / Replace | Detail |
|---|---|---|
| **1. Forcing** (scripts 00–04) | **Keep — Python** | THREDDS GFDL-MOM6-COBALT2 → per-cell time series: SST (`tos`), seafloor `tob`, plankton intercept+slope (`GetPPIntSlope` from `phyc`/`zooc`), detritus input (`expc-bot`), depth (`deptho`), mask → zarr. Unchanged. |
| **2. Engine** (`sizemodel`, `gridded_sizemodel_rk4`) | **Replace → dbpmr C** | one correct/guarded/tested engine, run per 0-D water column |
| **3a. Fishing effort** (`gravitymodel`, `effort_calculation`) | **Keep — Python** | fleet-gravity: distributes effort across cells by biomass/depth/accessibility → per-cell fishing mortality. Spatial coupling of the **fleet**, not the fish. |
| **3b. Calibration** (`LHSsearch`→optimiser, `getError`) | **Keep+simplify — Python** | minimise catch RMSE over `{search_vol, fmort_u, fmort_v}` — see "Calibration" below |
| **4a. Orchestration (parity)** | **Keep — Python** | loop regions/cells (dask), call the engine per cell, write zarr/netcdf. Cells **independent** (no fish movement) — matches FishMIP today, embarrassingly parallel. |
| **4b. Fish movement (optional extension)** | **dbpmr-enabled, later** | dbpmr natively supports advection+diffusion of fish across cells (`spatial_dim`, `Cfun`/`Dfun`, diffusion solver) — a capability FishMIP lacks. Couples neighbours → not trivially parallel → beyond parity. Relevant to migratory species (tuna range across basins). |

**Two different "spatials", kept separate:** fish population dynamics are
per-cell independent (4a) *unless* movement is switched on (4b); fishing effort is
spatially coupled via the gravity model (3a) regardless. The gravity model moves
**effort**, not fish.

## The one enabling change: dbpmr issue #5 (in-memory column driver)

dbpmr currently talks through **text files** — fine for one run, fatal for a grid.
The prerequisite is a clean `.Call` interface:

```
dbpm_column(params, forcing[t], init_state) -> {biomass[t], catch[t], spectra}
```

one water column, arrays in / arrays out, no I/O, accepting **time-varying**
forcing (temperature, plankton intercept/slope, detritus input, fishing mortality
per step). Everything else hangs off this. Bind to Python via `cffi`/`ctypes` on
the refactored C step (keeps R out of the hot loop).

## Unit / parameter mapping (all known)

| item | mapping |
|---|---|
| density convention | dbpmr ln-mass vs FishMIP log10 — the `ln(10)` factor (dbpmr is the correct side; the sizemodel PRs fix the other) |
| output units | apply FishMIP `min(depth,200)` (pelagic) / `benthic_habitat_depth` (benthos) areal multipliers on output (§6a) |
| energy budget | faithful `K/R/Ex` from `(1-def)·{K, 1-(K+AM), AM}` (#22) |
| **senescence** | **use dbpmr's hyperbolic form**, not sizemodel's power law — this is what preserves the realistic large-fish tail (tuna in oligotrophic gyres, ~uniform spectrum slope). A *feature* of the swap, not just parity. |

## Calibration — replace LHS with an optimiser

`LHSsearch` currently random-samples ~3 free params (`search_vol`, `fmort_u`,
`fmort_v`; the two `fminx` are fixed) and **picks the best draw** against a scalar
catch RMSE. That's an inefficient way to minimise a scalar objective in 3-D.

- **Use a bounded, derivative-free optimiser.** The RMSE surface has **collapse
  cliffs** (below a search-volume threshold) and — pre-guard — **NaN regions**
  (overflow), so a naïve gradient method fails. Robust choices: DE (`DEoptim`) or
  a global `nloptr` (ISRES/CRS) to find the basin, then Nelder-Mead/BOBYQA to
  refine; or Bayesian optimisation (`ParBayesianOptimization`) for a
  sample-efficient expensive-simulator fit. Bound/transform to the stable region.
- **Decompose when possible:** `search_vol` sets biomass scale, `fmort` sets catch
  given biomass → two ~1-D problems, each ~10–20 runs.
- **Synergy with the swap:** dbpmr is fast **and guarded**, so evaluations are
  cheap and the surface is smoother (fewer NaN cliffs). LHS was partly a defence
  against the fragility we've now fixed → the fix *enables* the optimiser.

**Demonstrated (Tier-1 test, LME-14, local data):** a 1-D `optimise()`/Brent
calibration of `search_vol` to a target consumer biomass converged in **8
evaluations for dbpmr** (search_vol 7.8, biomass 4.99 vs target 5) and 14 for the
fixed CMIP5 sizemodel — versus the hundreds a comparable LHS needs. Honest
caveats surfaced: the two engines have different productivity (the CMIP5 target
was unreachable in-range → it hit the ceiling), so each engine is calibrated
**separately** and targets must be reachable — both expected on an engine swap.

## Staged plan (each step verifiable)

1. **Land the correctness/robustness PRs** (ln(10): #2/#14/#15/#16; NaN guard: #3)
   so both engines agree and the surface is well-behaved. *(done/open)*
2. **dbpmr #5:** in-memory column driver accepting time-varying forcing. *(enabling)*
3. **Python binding** (`cffi`) + swap **`getError` first** (aspatial calibration,
   one LME) with an optimiser — smallest verifiable seam. *(Tier-1, partly demoed)*
4. **Transient obsclim** per-LME run (time-varying forcing) + optimiser calibration
   to observed catch time series. *(Tier-2)*
5. **Gridded runner (4a)** — map the column driver over cells (dask), with
   `gravitymodel()` computing effort between yearly chunks. Validate one region
   vs the hardened Python before global. *(Tier-3)*
6. *(optional)* **Fish movement (4b)** — turn on dbpmr's advection/diffusion for
   migratory dynamics. *(extension, beyond parity)*

Steps 2–3 are the real work; 4–5 are wiring.

## Efficiency levers
- In-memory C (#5) removes the file-I/O bottleneck.
- Compiled inner loop + **embarrassingly parallel over cells** (one cell in memory
  → no global size×cell×time array → low memory, trivial dask scaling).
- **Coarse-grained calls** (a whole year of sub-steps per C call) minimise overhead;
  region-independent kernel matrices computed once.
- dbpmr sub-steps internally where productive cells need a finer step — the
  orchestrator stays coarse.
- Optimiser calibration: tens of runs, not hundreds.

## Honest trade-off
Because `gridded_sizemodel_rk4` is already the same equations, once the ln(10)+NaN
PRs land it is *effectively* "dbpmr's correct engine in Python." So:
- **Harden-the-Python** (PRs + guards + hyperbolic senescence): least effort, no
  cross-language interface, the team already runs it — best for *correct output soon*.
- **Swap-in-dbpmr-C** (this plan): one tested/fast/guarded engine for both
  calibration and gridded runs, no duplicated math, and the speed/memory win at
  global scale — best for *one engine, long-term*, which is the stated goal.

Recommend the swap, **staged**: land the PRs → #5 → bind + swap `getError` →
gridded + gravity → (optional) movement.

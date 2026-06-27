# Design report: evolving `dbpmr` into a FishMIP-ready spatial ecosystem platform

**Status:** design only — no code changes proposed in this document.
**Scope:** integrate the strengths of `dbpmr` (ecological engine) and
`lme_scale_calibration_ISMIP3a` (FishMIP workflow) under the target architecture:

```
Climate forcing → LME preprocessing → Adapter layer → dbpmr → FishMIP outputs
```

`dbpmr` is the long-term modelling framework; the LME repository contributes
reusable workflow components (forcing prep, calibration, temperature functions,
gravity fishing, FishMIP I/O). Scientific equations are preserved unless a change
is explicitly requested, and existing results must remain reproducible.

---

## 1. `dbpmr` architecture

### 1.1 Package structure

```
dbpmr/
  DESCRIPTION, NAMESPACE        # roxygen2-generated NAMESPACE
  R/                            # S4 classes + setup/run/read/plot functions
    classes.R                   # all S4 class definitions + show() methods
    Setup.Run.r, Setup.Grid.r
    Setup.plankton.r, Setup.Pelagic.r, Setup.Benthic.r, Setup.Detritus.r
    Setup.ts.r, Setup.Rep.r, Setup.Fishing.r   # input-file writers
    SizeSpectrum.r              # marshals params → .C(C_SizeSpectrum)
    Read.In.r                   # parses text output → results objects
    Extract.Time.r, Average.Time.r, Plot.Spectrum.r, Points.Spectrum.r, Animate.r
    Calculate.du.r              # experimental, internal
    dbpmr-package.R             # @useDynLib, imports, references
  src/
    SizeSpectra.c               # ~3,700-line C ecological core
    init.c                      # R_registerRoutines (native-symbol registration)
  man/  tests/testthat/  vignettes/dbpmr.Rmd
```

The ecological core is **C (not C++)**, called through the registered `.C`
interface (`.C(C_SizeSpectrum, ...)`).

### 1.2 S4 classes

| Class | Role | Created by |
|---|---|---|
| `run.params` | run-level config (no. species, spatial dim, coupling, diff method) | `Setup.Run()` |
| `grid.params` | mass/time/space discretisation | `Setup.Grid()` |
| `plankton.params` | plankton resource spectrum | `Setup.Plankton()` |
| `pelagic.params` | pelagic species (physiology, feeding prefs, reproduction, fishing) | `Setup.Pelagic()` |
| `benthic.params` | benthic species | `Setup.Benthic()` |
| `detritus.params` | detritus pool | `Setup.Detritus()` |
| `plankton/pelagic/benthic/detritus.results` | per-component output (densities, growth, mortality, fishing, biomass budgets) | `Read.In()` |
| `singlefile.results` | a single named output file | `Read.In(filename=)` |
| `timestep.data` | one time slice of a spectrum | `Extract.Time()` / `Average.Time()` |

The S4 parameter system is the **most modular and reusable asset** in either
repository and should remain the public configuration surface.

### 1.3 Compiled code (`src/SizeSpectra.c`)

C structs mirror the S4 classes: `RUN`, `GRID`, `PLANKTON`, `PELAGIC`,
`BENTHIC`, `DETRITUS`, plus `COMMUNITY` (aggregate) and `MATRIX` (tridiagonal
workspace). Key routines:

- `SizeSpectrum(...)` — entry point; builds structs from the flat param vectors.
- `calculate_results(...)` — the main time loop.
- `calculate_g_and_mu(...)`, `g_pel/g_ben/g_det`, `mu_pel/mu_ben/mu_det`,
  `*_biomass(...)`, `phi(...)` — **growth and mortality (the rate calculations)**.
- `calculate_reproduction(...)` — reproduction methods 0–3.
- `calculate_fishing(...)` — fishing mortality (prescribed time series or computed).
- `mass_solver(...)` (`tridag`/`trimul`) — implicit tridiagonal McKendrick–von
  Foerster size solver.
- `xmove_solver/ymove_solver`, `Cfun/Dfun/Diffun` — **spatial movement**
  (diffusion/advection of the pelagic spectrum; Castle 2011).
- `free_mem`, `safe_malloc/fopen`, `print_*` — memory and text-file output.

The core is **temperature-independent**; plankton is a fixed analytic spectrum;
fishing is prescribed. These three are exactly the extension points (§1.7).

### 1.4 Run workflow

```
Setup.Run → Setup.Grid → Setup.Plankton / Setup.Pelagic / Setup.Benthic / Setup.Detritus
          → [optional] Setup.ts / Setup.Rep / Setup.fishing   (write <run>/Input/*.txt)
          → SizeSpectrum(...)                                  (C engine; writes <run>/<sp>/*.txt)
          → Read.In(...)                                       (text → results objects)
          → Extract.Time / Average.Time → Plot.Spectrum / Points.Spectrum / Animate
```

### 1.5 Input objects

| Input | Representation | C side | Used by | Required? |
|---|---|---|---|---|
| **Run** | `run.params` S4 | `RUN` struct (int vector) | all setup/run fns | **Required** |
| **Grid** | `grid.params` S4 | `GRID` struct (double vector) | `SizeSpectrum`, input writers | **Required** |
| **Species** (plankton/pelagic/benthic/detritus) | `*.params` S4 | `PLANKTON/PELAGIC/BENTHIC/DETRITUS` structs (double + flag vectors) | `SizeSpectrum` | ≥1 pelagic required; benthic/detritus optional |
| **Initial conditions** | `<run>/Input/<sp>_ts.txt` (written by `Setup.ts`) | read as text in `setup_*` | `SizeSpectrum` | Optional (`initial_flag`) |
| **Fishing** | `<run>/Input/<sp>_fishing_ts.txt` (written by `Setup.fishing`) | read in `calculate_fishing` | `SizeSpectrum` | Optional (`fishing_flag`) |
| **Reproduction** | `<run>/Input/<sp>_rep_ts.txt` (written by `Setup.Rep`) | read in `calculate_reproduction` | `SizeSpectrum` | Optional (`rep_method = 1`) |
| **Environmental forcing** (temperature, ESM plankton) | **does not exist** | **does not exist** | — | **Missing — must be added** |

The two structural gaps for FishMIP are clear: **(a) there is no forcing input
object at all**, and **(b) all I/O is via text files** rather than in-memory
arrays.

### 1.6 Output objects

- C writes per-component text files under `<run>/<species>/`: `results.txt`
  (densities), `growth.txt`, `mortality.txt`, `fishing.txt`, `summary.txt`
  (biomass budgets), etc.
- `Read.In()` parses these into the `*.results` S4 objects.
- There is **no direct in-memory return** from the C call — every quantity round-
  trips through disk as CSV.

### 1.7 Extension points (where new science hooks in)

1. **Rate calculations** (`calculate_g_and_mu`, `g_*`, `mu_*`) — the natural
   home for temperature scaling of search/metabolism/mortality.
2. **Plankton resource** (`setup_plankton`, and where the plankton spectrum is
   read in the time loop) — the hook for ESM-forced, time-varying plankton.
3. **`calculate_fishing`** — the hook for effort/gravity-based fishing mortality.
4. **The `.C` parameter vectors** — must be widened to carry forcing time series
   (or a new forcing input file/array).
5. **R layer** (`Setup.*`, a future gridded driver) — pure R, the cheapest place
   to add adapters, masking, and per-cell orchestration.

---

## 2. LME repository (`lme_scale_calibration_ISMIP3a`) architecture

Two parallel engine implementations exist — an **R** one (used for equilibrium /
no-fishing / calibration runs) and a **Python** one (used for production gridded
runs). Component locations:

| Component | Location | Notes |
|---|---|---|
| **Size-spectrum engine (R)** | `scripts/useful_functions.R` → `sizemodel()`, `sizeparam()`, `run_model()` | implicit-upwind finite difference; pelagic + benthic + detritus; `temp_effect=T` |
| **Size-spectrum engine (Python, gridded)** | `useful_functions.py` / `scripts/08_run_dbpm_gridded.py` → `gridded_sizemodel_rk4(rk4_substeps=4)` | RK4, xarray/zarr gridded, **per-cell, no inter-cell movement** |
| **Temperature functions** | inside `sizemodel()` | Boltzmann–Arrhenius on `tos` (surface) and `tob` (sea-floor); coeffs `c1`, `activation_energy`, `boltzmann`; `benthic_habitat_depth` |
| **Gravity fishing model** | `gravitymodel(effort, prop_b, depth, iter)` (R); `effort_calculation()` (Python) | Walters & Bonfil effort redistribution by biomass × depth suitability |
| **Sea-ice / access masking** | `sea_ice_mask` field + depth correction (`depth_corr`, cap 200 m) in `08_run_dbpm_gridded.py` | applied to effort/accessibility |
| **Forcing preparation** | `scripts/00–03*` (Python) | GFDL-MOM6-COBALT2 → regional/gridded inputs; plankton intercept/slope, export ratio, spinup |
| **Fishing input processing** | `scripts/04_processing_effort_fishing_inputs.R` | ISIMIP3a reconstructed effort → model inputs |
| **Calibration workflow** | `getError()`, `LHSsearch()` (Latin hypercube), `corr_calib_plots()`; `scripts/07_estimating_best_vals_fishing_params.R` | optimises search volume / catchability vs observed catch (RMSE) |
| **`gridded_sizemodel()`** | the Python gridded engine (`gridded_sizemodel_rk4`) | the component the target architecture intends `dbpmr` to **replace** as the integrator |
| **FishMIP outputs** | `scripts/08–09*` | netCDF, FishMIP-compliant directory structure |
| **Diet/metabolism lookups** | `phi_f`, `gphi_f`, `mphi_f`, `expax_f` | precomputed kernels — conceptually the same role as `phi()` in the C core |

**Key finding:** the LME `sizemodel()` is the *same DBPM science* as the C core
but a separate, less modular implementation, **plus** three things the C core
lacks (temperature forcing, ESM plankton forcing, gravity fishing) and **minus**
the C core's spatial-movement solver (the gridded runs treat cells independently).

---

## 3. Compatibility assessment

For each component: **keep** (already in dbpmr), **copy** (port the equation into
dbpmr's core), **wrap** (dbpmr accepts its output as an input; logic stays in
R/LME), or **remain** (stays a workflow step).

| Component | Decision | Why |
|---|---|---|
| S4 parameter system, results classes | **Keep** in dbpmr | Already the most modular surface; the adapter targets it. |
| Compiled size solver, growth/mortality core | **Keep** in dbpmr | The long-term engine; scientifically equivalent to `sizemodel()` but faster and modular. |
| Spatial movement (0D/1D/2D) | **Keep** in dbpmr | Unique to dbpmr; not in LME. Optional for per-cell FishMIP runs but a long-term differentiator. |
| **Temperature scaling** | **Copy** the equation into the C rate calculations, **but driven by R-supplied forcing** | The effect must be applied *inside* the integration loop, so it has to reach C. Recommended split: keep the Boltzmann–Arrhenius *formula* as a reusable R function (LME), have it produce per-timestep scaling factors, and extend the C engine to accept those as forcing. Preserves the published equation and keeps it testable in R. |
| ESM plankton forcing (intercept/slope time series) | **Wrap** (new forcing input to dbpmr) | dbpmr's plankton is currently static; add a time-varying plankton-resource forcing input. The ESM-extraction stays in LME preprocessing. |
| **Gravity fishing model** | **Remain** in LME (R), feeding dbpmr's fishing input | It redistributes *effort over space* using biomass fields — a workflow-level spatial operation. It can emit an F(size, space, time) field that dbpmr consumes through the (now-corrected) fishing-input path. No need to put it in C. |
| Sea-ice / accessibility / depth masking | **Remain** in LME / adapter | Pure masking of forcing and effort; belongs in preprocessing. MPA masking (roadmap §4.5) is the same pattern. |
| Forcing preparation (00–03) | **Remain** in LME | netCDF/ESM wrangling; not modelling. The adapter consumes its products. |
| Calibration (`getError`, `LHSsearch`) | **Remain** in LME, calling dbpmr as the engine | An optimisation loop *around* the model; engine-agnostic once the adapter exists. |
| `gridded_sizemodel()` (the integrator) | **Superseded** — decompose | The per-cell *driver* (cell loop, I/O, masking) → workflow/adapter; the *integration* → dbpmr. This is the central substitution of the project. |
| FishMIP output writers | **Remain** in LME / adapter | netCDF + directory conventions; dbpmr returns raw spectra, the adapter formats them. |
| Diet/metabolism lookups (`phi_f` etc.) | **Keep** dbpmr's `phi()` | Equivalent role already in the core; reconcile parameter definitions during validation rather than copying code. |

**Guiding principle:** put *equations that act inside the time step* in dbpmr
(temperature scaling), and keep *everything that prepares inputs or post-processes
outputs* (forcing prep, gravity effort allocation, masking, calibration, FishMIP
I/O) in the LME workflow, connected through a thin **adapter layer**.

---

## 4. Proposed development roadmap

Staged, prioritising only what the next FishMIP experiments need. **Stage 0 is
added** to honour the "existing results remain reproducible" requirement.

**Stage 0 — Reproducibility baseline & scientific cross-check.**
Pin a reference LME run (one region, one forcing file). Run the same
configuration through both `sizemodel()` and `dbpmr` and document where outputs
agree/diverge (numerics: implicit-upwind vs tridiagonal; scheme differences).
This defines "scientific continuity" concretely before anything is ported.

**Stage 1 — Benchmark current performance.**
Time and memory-profile: (a) a single `dbpmr` 0-D run; (b) a single LME
`sizemodel()` run; (c) the text-file I/O cost of a `dbpmr` run. Establish the
per-cell cost that gridded scaling will multiply.

**Stage 2 — Adapter layer (R).**
New functions converting LME gridded/regional forcing → `dbpmr` `Setup.*` objects
and forcing inputs, running `dbpmr` per cell, and collecting outputs. Begin 0-D
per-cell (matches LME's independent-cell approach). This is the backbone; it can
ship and add value before any C changes (using existing `dbpmr` features).

**Stage 3 — Temperature forcing.**
Extend the C engine to accept per-timestep temperature (or precomputed scaling)
and apply Boltzmann–Arrhenius in the rate calculations; keep the formula in an R
function. Validate one cell against `sizemodel(temp_effect=T)`.

**Stage 4 — Gravity fishing integration.**
Drive `dbpmr`'s fishing input from the LME `gravitymodel()` output. Requires the
fishing-input path to be correct (note: the `Setup.fishing`/`Setup.ts`/`Setup.Rep`
path bug was fixed this session) and possibly extended for size-based selectivity
per group.

**Stage 5 — MPA / accessibility masking.**
Generalise the sea-ice-mask concept to MPA masks applied in the adapter (zeroing
effort/accessibility per cell-time). No engine change needed.

**Stage 6 — Additional functional groups.**
`dbpmr` already supports multiple pelagic/benthic species; extend the adapter and
parameter mapping to the FishMIP functional groups.

**Stage 7 — FishMIP-compatible outputs.**
Adapter writes FishMIP netCDF (reusing the LME writers) from `dbpmr` results.

**Stage 8 — Performance optimisation (only after correctness).**
Profile-guided; see §5. The likely first target is replacing text-file I/O with
an in-memory return path, not rewriting science.

The first externally useful milestone is **Stage 2 + Stage 3** (forced per-cell
runs through dbpmr); Stages 4–7 make it a full FishMIP pipeline.

---

## 5. Performance

### Likely bottlenecks (ranked for gridded FishMIP scale)

1. **File I/O (dbpmr) — the dominant risk.** Every run writes/reads multiple
   per-species CSV files. At thousands of grid cells × decades this is enormous
   small-file churn and parsing. **This, not the maths, is what blocks gridded
   scaling.**
2. **R workflow overhead.** Per-cell loop in R assembling inputs, calling `.C`,
   parsing text outputs. Significant when multiplied by cell count.
3. **Compiled ecological core.** The tridiagonal size solver per timestep is
   already efficient C; per cell it is fast. Cost scales with `tmax/tstep` ×
   cells but is well-optimised.
4. **Movement.** Only for 1D/2D. FishMIP gridded runs are per-cell (independent),
   so movement is **not on the critical path** for the primary use case.
5. **Fishing / gravity model.** Effort redistribution over the whole grid each
   step can be costly at high resolution — but it is a vectorisable R/array
   operation, not engine code.
6. **Temperature calculations.** Cheap — a scalar (or per-cell) multiplier per
   timestep. Negligible wherever it lives.

### What should stay in R vs move to C

**Keep in R (do not move):** forcing preparation, calibration loop, gravity
fishing, masking, FishMIP output formatting, and the per-cell driver loop. These
are I/O- and orchestration-bound, benefit from the geospatial stack, and must
stay readable/testable for scientific review.

**Engine-side priorities (in order):**
1. **Replace text-file I/O with an in-memory `.Call` return** (densities/budgets
   as R arrays). Single highest-impact change for gridded scale; eliminates the
   #1 and much of the #2 bottleneck. (Aligns with finishing the `.C` → `.Call`
   modernisation deferred earlier.)
2. **Accept forcing in memory** (temperature, plankton) via the same interface,
   so per-cell runs need no disk round-trip.
3. **Only then** consider batching many cells into one C call to amortise R/`.C`
   overhead — and only if profiling (Stage 1/8) shows R-call overhead dominates
   after I/O is fixed.

**Do not** port gravity fishing, calibration, masking, or temperature *formulae*
into C for speed — none are bottlenecks, and doing so would reduce modularity and
scientific transparency. Temperature scaling enters C only because it must act
inside the time step, not for performance.

### Summary recommendation

The aim is evolution, not rewrite. The combination that delivers a modern,
modular, FishMIP-ready platform while preserving scientific continuity is:

- **dbpmr** keeps the modular S4 system, the compiled solver, spatial movement,
  and gains temperature forcing + a forcing input + an in-memory return path.
- **The LME repository** keeps forcing prep, gravity fishing, masking,
  calibration, and FishMIP I/O as reusable workflow components.
- **A thin R adapter** connects them, and **Stage 0 reproducibility checks**
  guard scientific continuity at every step.

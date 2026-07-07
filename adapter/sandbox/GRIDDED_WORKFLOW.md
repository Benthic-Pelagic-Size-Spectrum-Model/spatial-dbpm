# Gridded LME/FAO DBPM + Spatial Gravity — Workflow & Assumptions

Reproducible workflow for the per-cell, time-varying, dynamic gridded DBPM run driven by an
ideal-free-distribution ("gravity") fishing-effort model, calibrated per LME/FAO region.
All compute runs on a laptop (Apple M3 Pro, 12 cores). Scripts live in `~/dbpm_compare_scratch/`.

--------------------------------------------------------------------------------
## 0. Overview

The ocean is tiled into **83 FAO-LME regions** (66 LMEs + 17 FAO high-seas areas; mask =
`FAO-LME-corrected_1degmask_DBPM`, 1° grid). For each region:
1. **Calibrate** the 0-D size-spectrum model (`dbpmr`) to the region's observed catch → q_pel, q_ben.
2. **Grid it**: run every 1° ocean cell in the region as an independent 0-D `dbpmr` column with
   its own forcing, coupled only through an annual **spatial gravity** re-allocation of the
   region's total fishing effort.
3. **Stitch** all regions into global fields; validate against observed catch and Global Fishing
   Watch / Clawson effort.

Engine: `dbpmr` (compiled C size-spectrum solver; two spectra — U = pelagic+demersal predators,
V = benthic detritivores — plus a detritus pool). Fishing minimum size = 10 g knife-edge (v1).

--------------------------------------------------------------------------------
## 1. Data inputs

| Input | Source | Used for |
|---|---|---|
| DBPM region inputs `dbpm_clim-fish-inputs_fao_lme-<L>_*.parquet` | ISIMIP3a preprocessing (`DBPM_dev/dbpm_inputs`) | LME monthly plankton intercept/slope, tos, tob, export_ratio, effort, catch |
| Equilibrium init `init_dbpm_nonspatial_fao_lme-<L>_searchvol_*.json` | `DBPM_dev/equilibrium_runs` | biology params (global constants), depth |
| 1° FAO-LME mask CSV | portal.sf.utas.edu.au THREDDS `masks/FAO-LME_masks` | region → cells; keying |
| `phyc`, `phypico` (per-level) | ISIMIP3a THREDDS OPeNDAP | per-cell biomass-weighted plankton intercept |
| `tos`, `tob`, `expc-bot`, `intpp`, `phyc-vint` | ISIMIP3a THREDDS OPeNDAP | per-cell temperature, export, plankton-temporal |
| `siconc` (monthly) | ISIMIP3a THREDDS OPeNDAP | per-cell time-varying sea-ice gate |
| `gfw_static_spatial_measures.csv` (elevation, dist-to-shore/port) | GFW / Clawson `effort_manuscript` repo | per-cell depth + accessibility |
| Clawson gridded effort (Zenodo 19600603) | data.mapping-global-fishing.cloud.edu.au | effort validation |
| `catch_histsoc_1869_2017_EEZ_addFAO.csv` + `TaxGrps.xlsx` | FishMIP THREDDS + Julia | v2 group-split catch (NOT in v1) |

Variable definitions: FishMIP2.0 protocol Table 6 (github.com/Fish-MIP/FishMIP2.0_ISIMIP3a).

--------------------------------------------------------------------------------
## 2. Pipeline stages (scripts, in order)

### Stage A — LME calibration (0-D)  →  `calib_dintall/`, `calib_final/`
`adapter/sandbox/tier1_fishing_calib.R <L>` (env: `DINT_CSV`, `TWO_Q`, `QMAX`, `SAVE_RDS`)
- Spins the 0-D column to unfished equilibrium, applies the LME effort series, fits catchability
  to observed catch by minimising **log-MSE** (mean squared log10 catch residual).
- **single-Q**: one Q ∈ [0, QMAX=4], split q_pel=Q·s_pel, q_ben=Q·s_ben by unfished biomass share.
- **two-Q** (`TWO_Q=1`): free (q_pel, q_ben), grid seed + `nloptr` BOBYQA.
- **Biomass-weighted plankton intercept**: `DINT_CSV=lme_dint_all.csv` shifts each LME's plankton
  intercept by its measured biomass-weighted value (∫phyc²/∫phyc vs fixed 0-200 m mean;
  `all_int.R` over the FAO-LME mask). Median +0.38 dex; de-saturates Q in upwelling/polar LMEs.
  Skill vs the fixed 0-200 m intercept (`plot_200m_vs_bw.R` → `calib_200m_vs_bw.pdf`, for Denisse):
  median corr 0.46→0.51, median RMSE(log10) 0.28→0.23; 56% of LMEs improve correlation.
  - **CAVEAT — this is a "proxied" static offset, not a fully per-timestep series.** `dint` is ONE
    number per LME (int_bw − int200 from the `phyc` climatology) added to the parquet's fixed-200 m
    *time-varying* series. The temporal variation is still the fixed-200 m signal; only the mean level
    is biomass-weighted. The protocol-faithful **v2** change is to re-derive the intercept per timestep
    from `phyc` via new_features `integrating_phyto(mode="biomass_weighted")` + `GetPPIntSlope`, so both
    level AND temporal shape are biomass-weighted. Same distinction applies to the gridded driver.
- **Model selection** (`aic_select.R`): per LME, ΔAIC = n·ln(MSE₂Q/MSE₁Q)+2; pick two-Q only if
  ΔAIC<0 AND correlation not worsened (guard — plain ΔAIC over-selects on ~50-yr autocorrelated
  series). Final: 40 two-Q / 40 single-Q. Writes `calib_final/`.
- Batch: `run_all` style loop over the 80 LMEs; the 3 gaps (32,148,188) calibrated separately.

### Stage B — per-cell forcing assembly (OPeNDAP)
- `all_int.R` / `prep_lme.R <L>` → per-cell biomass-weighted intercept `int_lme<L>.csv`;
  per-cell annual `siconc_lme<L>.csv` (high-lat only). Dateline-safe (per contiguous lon-run), retried.
- `pull_tempexport_global.R` + `build_forcing_percell.R` → `tempexport_percell.csv`
  (per-cell tos, tob, export_ratio=expc-bot/intpp, phyc-vint clim) + `phycvint_annual.rds` (temporal).
- `gfw_static.csv` (depth, dist-to-shore) downloaded once.

### Stage C — gridded dynamic run  →  `grid_v2/` (or `grid_all/`)
`gridded_gravity_dynamic.R <L> --ncell=0 --spinyr=30 --cores=10 --int_csv=int_lme<L>.csv
   --qpel=.. --qben=.. [--siconc_csv=siconc_lme<L>.csv]`
Orchestrated by `run_all_v2.sh` (resumable, small-region-first, `caffeinate`, TMPDIR on a RAM disk).
- Reads q from `calib_final/lme<L>.rds`; per-cell forcing auto-loaded from `tempexport_percell.csv`.
- Loop per year 1961-2010: gravity-allocate effort → run every cell 1 model-year via **warm-restart**.

### Stage D — synthesis & validation
- `render_global.R` (or `render_global_dir.R` with `GRID_DIR`): stitches cells → global catch
  time series + world maps of effort/catch/biomass (`global_*.png`, `global_synth.rds`).
- `compare_catch.R`: per-LME observed vs 0-D-aggregated vs gridded catch (`catch_validation.pdf`),
  with log-RMSE and bias (dex). Gridded catch = area-weighted (cos-lat) mean of per-cell catch;
  cells with catch<0 or non-finite (numerical blow-ups) excluded.

--------------------------------------------------------------------------------
## 3. Model equations

**Per cell (0-D dbpmr), one model-year, warm-restart:** two size spectra advanced at dt=1/48 yr;
fishing mortality F(w) = q · effort · sel(w), sel = knife-edge at 10 g (v1). Detritus = Dunne-2007
burial closure; benthic food from sinking export + dead bodies.

**Warm-restart (validated exact, 6e-11):** C engine reads the initial spectrum from
`<run>/Input/<name>_ts.txt` when `obj@initial_flag<-TRUE`; write year-t `finaluvals` as year-(t+1)
initial state (detritus = scalar `@biomass`).

**Spatial gravity (annual, single fleet):**
- attractiveness  a_i = (q_pel·B^U_i + q_ben·B^V_i) · acc_i · ice_i
- accessibility   acc_i = exp(-depth_i/H) · exp(-shore_i/S)   (H=800 m, S=150 km)
- ice gate        ice_i = 1[siconc_i < 15%]   (hard; fisheries avoid ice)
- allocation      e_i = E_tot(year) · a_i / Σ_j a_j   (total effort conserved; mean cell effort = E_tot)
- inertia (opt)   e_i(t+1) = (1-α)e_i(t) + α·e_i^IFD  (α=1 default = instant)

**Per-cell forcing factorization** (keeps LME-mean = parquet series so calibrated q stays valid):
- plankton  intercept_i(t) = LME_intercept(t) + (int_bw_i − LME_mean) + log10(phyc-vint_i(t)/clim)_centered
- temp      tos_i(t) = LME_tos(t) + (tos_clim_i − LME_mean);  tob likewise
- export    er_i(t)  = LME_er(t)  + (er_clim_i − LME_mean)

--------------------------------------------------------------------------------
## 4. ASSUMPTIONS — resolution of every variable

| Variable | Resolution | Notes |
|---|---|---|
| Depth, accessibility (depth+shore) | **per-cell** | GFW static layers |
| Sea-ice gate | **per-cell, per-year** | real siconc, hard 15% |
| Plankton — spatial | **per-cell** | biomass-weighted intercept |
| Plankton — temporal | **per-cell** | phyc-vint annual (v2) |
| Temperature tos/tob | **per-cell** | climatology anomaly (v2) |
| Export ratio | **per-cell** | expc-bot/intpp anomaly (v2) |
| Biomass state | **per-cell** | warm-restart |
| Allocated effort | **per-cell** | gravity share |
| **Total effort** | **per-LME** | fleet total; per-cell only via gravity — *correct, not an artifact* |
| **q_pel, q_ben** | **per-LME** | calibrated to region catch; cannot be per-cell without per-cell catch → causes the residual boundary "rectangles" |
| Biology params (metabolic exponent 0.82, natural mort 0.20, growth/energy K/R/Ex, reproduction, senescence, c1/E/Boltzmann) | **global constants** | identical across LMEs → no blocks |
| Search volume A.u=64 (pelagic), A.v=6.4 (benthic) | **global** | hardcoded; JSON `searchvol` label unused |
| Fishing min size | **global 10 g** knife-edge (v1) | → per-region per-spectrum in v2 |

**Key structural assumptions:** catch = q·effort·biomass per cell; effort reallocation is annual
(effort data is annual); cells coupled only at year boundaries; warm-restart exact; all forcing
anomalies LME-centered (preserves q consistency). Biology/q per-LME; forcing per-cell.

--------------------------------------------------------------------------------
## 5. Known limitations (v1) & the v2 plan

1. **Rectangles remain even with per-cell forcing** — the between-LME steps are driven by per-LME
   **q** (a multiplicative catch factor) + LME-mean levels, which the LME-centered factorization
   preserves. Per-cell forcing adds within-block texture, not between-block smoothing
   (between-LME biomass variance 40%→41%). To smooth: spatially interpolate q at boundaries, or
   accept as genuine regional differences.
2. **Oligotrophic/reef over-prediction** — model treats 100% of benthos as fishable, but detritus
   is buffered (and on reefs consumed by corals/sponges that aren't fished). s_ben→1 at low
   productivity (cor(s_ben, plankton intercept) = −0.61) → 6-10× catch over-prediction, q pinned.
   Australian shelves, Antarctic, Pacific high-seas. NOT an effort/MPA issue (effort is coastal &
   matches obs; over-prediction is CPUE/biomass).
3. **Krill uncatchable** — 1-2 g, below the 10 g knife-edge → FAO 48 (83% krill) uncalibratable.

**v2 (staged, data in hand):** split observed catch into pelagic (fish + krill + cephalopod) vs
benthic (shrimp, lobster/crab, mollusc) using `catch_histsoc` + `TaxGrps.xlsx`; calibrate
q_pel/q_ben against the two series separately (properly identified). Resolves the oligotrophic/reef
bias, the krill gap, and the single-Q-pin/two-Q-overfit at once.

**Fished size — PER-SPECTRUM (U/V), time-varying, derived UPSTREAM (protocol).** The parquet's
single `min/max_fished_weight_class` is fish-derived and taxon-less; applied to both spectra it
excludes krill (U end) and the small-bodied benthic fishery (shrimp 3–60 g, mollusc 2–500 g; V end),
which guts benthos (California s_ben 0.24→0.05, corr→−0.16). Fix = per-spectrum U/V window derived
at the SAME workflow step that makes the fish window:
- **Upstream:** `lme-workflow/scripts/04_processing_effort_fishing_inputs.R` (`uv_summ` block) now
  derives `min/max_fished_U` and `min/max_fished_V` from the FGroup catch (`catch_histsoc`, already
  loaded there): FGroup→gram range + class MIDPOINT (matches the parquet's `log10mid_wt`; edges
  over-extend into the huge small-size tail); window = [min mid, max mid] over FGroups ≥0.5% of that
  spectrum-year catch, single-group years widen to [lo,hi]; U = fish+krill+ceph, V = shrimp/lobster/
  mollusc. These columns flow into the parquet via the existing merge. FAO 48 U: 7290–80000 g (1970,
  finfish) → 1–2 g (2010, krill); California U min ≈49 g (≈ the fish-derived 50 g), V 13–632 g.
- **Downstream:** `tier1_fishing_calib.R` reads them per (region, year) with fallback chain
  `FSIZE=knife` → parquet U/V cols → `FSIZE_UV_CSV` → single `min/max_fished_weight_class` → knife.
  Committed interim CSV for laptop runs (before the parquet is regenerated on Gadi):
  `fished_size_UV_tv.csv` (built by `build_fished_size_uv.R`, mirrors the script-04 logic).
- **GRAVITY (0-D and gridded) must use per-GROUP fishable biomass WITHIN each window**, not total
  biomass: `B_U^fish = ∫U over [minU,maxU]`, `B_V^fish = ∫V over [minV,maxV]`, then the gravity/
  effort split uses `q_pel·B_U^fish + q_ben·B_V^fish`. tier1's `catch_ts` already integrates each
  spectrum over its own window; the gridded driver's attractiveness `a_i` needs the same restriction.
- **Known tension (not a bug):** realistic size windows can *worsen* small-pelagic LMEs (California)
  where the fishery targets sizes at/below the derived min — the coarse `<30cm` FGroup can't resolve
  anchovy, and Reg's fine size data is fish-only. Two-Q + finer sizes are the next step; judge net
  effect on the full re-run, not per-LME.

--------------------------------------------------------------------------------
## 6. Reproduce from scratch (commands)

```sh
export R_LIBS=$HOME/dbpm_compare_scratch/dbpmrlib          # dbpmr installed to persistent lib
# RAM disk (dbpmr file I/O is the bottleneck; ~6 min spin vs would-not-finish on SSD):
DEV=$(hdiutil attach -nomount ram://12582912|awk '{print $1}'); diskutil erasevolume HFS+ dbpmram "$DEV"
export TMPDIR=/Volumes/dbpmram

# A. calibrate all regions (single-Q + two-Q + AIC guard -> calib_final/)
#    (batch tier1_fishing_calib.R over LMEs; then aic_select.R)
# B. assemble per-cell forcing (prep_lme.R per region; build_forcing_percell.R once)
# C. gridded run (resumable, caffeinated):
caffeinate -is sh run_all_v2.sh          # -> grid_v2/ ; writes grid_v2/V2_DONE.marker
# D. synthesis:
GRID_DIR=grid_v2 Rscript render_global_dir.R
Rscript compare_catch.R grid_v2/catch_validation.pdf
```

**Timing (M3 Pro, 10 cores):** ~0.66 s/cell wall; 43k cells ≈ 8-16 h (giant FAO/Arctic regions
dominate). Per-cell tempdirs MUST be unlinked (`on.exit(unlink)`) or they overflow the RAM disk.

**Robustness notes:** completion markers must be *files* (not `pgrep`-on-string — a stale monitor
process containing the script name will fake it out). `tier1` and the gridded driver both glob the
init JSON `searchvol_*` (Arabian Sea = 6.4, not 12.8). `runcell` returns NULL for negative/
non-finite biomass (numerical blow-ups) so they're excluded, not propagated.

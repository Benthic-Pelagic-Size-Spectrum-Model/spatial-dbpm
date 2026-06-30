# Stage 0: reference equilibrium runs and dbpmr ↔ LME parameter crosswalk

Companion to [`dbpmr-lme-integration.md`](dbpmr-lme-integration.md) and
[`lme-input-schema.md`](lme-input-schema.md). Documents the **no-fishing
equilibrium reference outputs** and the parameter mapping needed to configure a
matched `dbpmr` run for the Stage 0 scientific cross-check. **Design only.**

## 1. Reference output files

- Location: `…/DBPM_dev/equilibrium_runs/`
- **83 files**, one per LME: `init_dbpm_nonspatial_fao_lme-<N>_searchvol_12.8.json`
  (~37 MB each), plus rendered PDFs of spectra/biomass under `…/no-fishing/`.
- Produced by the LME R engine (`sizemodel()`), **no fishing**, run to
  equilibrium over the `stable-spin` block (constant climatology).

### Structure of one JSON

| Key | Shape | Meaning |
|---|---|---|
| `predators` | `[181][2400]` | pelagic predator spectrum: 181 size bins × 2400 monthly steps |
| `detritivores` | `[181][2400]` | benthic detritivore spectrum |
| `detritus` | `[2400]` | detritus pool through time |
| `growth_int_pred`, `growth_det` | `[181][2401]` | growth rates |
| `pred_mort_pred`, `pred_mort_det` | `[181][2401]` | predation mortality |
| `catch_pred`, `catch_det` | `[181][2401]` | catch (zero here — no fishing) |
| `params` | dict | full model configuration (below) |

The **equilibrium state** is the final time column of `predators` /
`detritivores` / `detritus`. (181 bins, all non-zero at equilibrium.)

## 2. Grid alignment (the encouraging part)

| | LME `sizemodel()` | `dbpmr` | Note |
|---|---|---|---|
| mass axis | **log10** g | **natural log** | factor `ln 10 ≈ 2.3026` |
| range | −12 … 6 (log10) | −28 … 14 (default, ln) | −12…6 log10 = −27.6…13.8 ln ≈ dbpmr default |
| step | 0.1 (log10) | 0.2 (default, ln) | 0.1 log10 = 0.23 ln ≈ dbpmr default |
| bins | 181 | ~211 (default) | match exactly by setting dbpmr `mmin=-12*ln10`, `mmax=6*ln10`, `mstep=0.1*ln10` |

The two grids are the same design; an exact match needs only setting `dbpmr`'s
`grid.params` to the log10×ln10 values.

## 3. Parameter crosswalk (LME `params` → `dbpmr`)

| LME param (value, LME 10) | `dbpmr` equivalent | Match? / conversion |
|---|---|---|
| `metabolic_req_pred` = 0.82 | pelagic `alpha` (default 0.82) | **exact** |
| `metabolic_req_detritivore` = 0.75 | benthic `alpha` (default 0.75) | **exact** |
| `natural_mort` = 0.2 | `mu_0` (default 0.2) | **exact** |
| `log10_pred_prey_ratio` = 2 | `q_0` = `log(100)` | **exact** (100:1) |
| `log_prey_pref` = 1 | `sig` (kernel width) | log10 vs ln — **convert** |
| `detritus_coupling` = 1 | `run@coupled_flag` | **exact** |
| `dynamic_reproduction` = 1 | `rep_method` (2/3, energy/biomass) | confirm which |
| `hr_volume_search` = 12.8 (calibrated) | pelagic `A` (search rate, default 640) | **different parameterisation — reconcile units** |
| `growth_pred` = 0.3, `energy_pred` = 0.5 | `epsilon` (assimilation/growth eff.) | reconcile (two LME coeffs vs one) |
| `defecate_prop` = 0.3, `def_low` = 0.5 | detritus routing | reconcile |
| `size_senescence` = 3, `exp_senescence_mort` = 0.3, `const_senescence_mort` = 0.2 | `mu_s` + senescence terms | reconcile |
| `c1`=25.22, `activation_energy`=0.63, `boltzmann`=8.62e-5 | **none — temperature scaling absent in dbpmr** | Stage 3 addition |
| `sea_surf_temp`, `sea_floor_temp` | **none** (forcing) | Stage 3 |
| `int_phy_zoo`, `slope_phy_zoo` (log10) | plankton `u_0`/`lambda` | log10→ln **convert**; static for now |
| `min_log10_pred`=−3 … `max_log10_pred`=6 | pelagic `mmin`/`mmax` | ×ln10 |
| `numb_size_bins`=181, `log_size_increase`=0.1 | grid `mstep` | ×ln10 |
| `fish_mort_pred`=0, `effort` (off) | fishing disabled | matches no-fishing run |

**Headline:** the core feeding/mortality/grid parameters already coincide with
`dbpmr` defaults (same lineage). The reconciliations are concentrated in:
(a) the **search-volume parameterisation** (`hr_volume_search` vs `A`),
(b) **temperature** (absent in dbpmr — but *constant* during the spin-up, so it
reduces to a fixed multiplier on rates that can be folded into the rate constants
for the Stage 0 static comparison), and
(c) **log10 ↔ ln** unit conversions for grid, plankton slope/intercept, and
kernel width.

## 4. Stage 0 comparison plan

1. Pick one LME (e.g. LME 10). Read its `params` from the equilibrium JSON.
2. Configure a `dbpmr` **no-fishing, aspatial (0-D)** run with the matched grid
   (§2) and parameters (§3); set the (constant) temperature effect as a fixed
   rate multiplier.
3. Run `dbpmr` to equilibrium; extract the final pelagic + benthic spectra.
4. Interpolate both onto a common natural-log mass axis and compare:
   - metric: max / RMS absolute difference in log10 abundance density across the
     overlapping size range;
   - target: shapes agree within numerical-scheme tolerance (LME uses
     implicit-upwind; dbpmr uses tridiagonal — small differences expected).
5. Document agreement and any systematic offset (this *defines* "scientific
   continuity" before any porting).

## 5. Blockers / asks before running Stage 0

The exact `sizemodel()` equations have now been read from source (see §7), which
**removes the main blocker**: at constant spin-up forcing, Stage 0 needs **no
engine change**. Remaining confirmations only:

1. Confirm `dynamic_reproduction = 1` ↔ `dbpmr` `rep_method` (2 vs 3).
2. Confirm the net pelagic growth efficiency mapping (LME splits it across
   `defecate_prop`, `growth_pred`, `energy_pred`; dbpmr uses a single `epsilon`).

## 7. Resolved equations (from `sizemodel()` source)

Verbatim relationships extracted from `scripts/useful_functions.R`:

**Search / encounter rate** — same functional form as `dbpmr`:
```
feed_mult_pel = hr_volume_search * 10^(log10_size * metabolic_req_pred) * pref_pelagic
              = A * w^alpha            # w = mass, alpha = metabolic_req_pred
```
→ set `dbpmr` pelagic `A = hr_volume_search * pref_pelagic` (× temp factor below),
`alpha = metabolic_req_pred` (already the default 0.82).

**Temperature (Boltzmann–Arrhenius)** — multiplies the feeding rate:
```
pel_tempeffect = exp(c1 - activation_energy/(boltzmann*(sea_surf_temp + 273)))
ben_tempeffect = exp(c1 - activation_energy/(boltzmann*(sea_floor_temp + 273)))
```
For LME 10 (T_surf = 24.74, T_floor = 1.28): pel ≈ **1.96**, ben ≈ **0.24**.
Constant during spin-up ⇒ **fold into `A`** for Stage 0:
`A_pel = hr_volume_search * pel_tempeffect`,
`A_ben = hr_volume_search * pref_benthos * ben_tempeffect`. (Becomes a true
per-timestep forcing in Stage 3.)

**Plankton resource** — clean base conversion:
```
density = 10^int_phy_zoo * w^slope_phy_zoo
```
→ `dbpmr` plankton `lambda = slope_phy_zoo`, `u_0 = 10^int_phy_zoo`
(the slope needs **no** log-base conversion — it is the exponent on mass itself).

**Assimilation / growth efficiency:**
```
growth_prop = 1 - defecate_prop            # 0.7
net pelagic growth efficiency = growth_prop * growth_pred   # 0.7 * 0.3 = 0.21
```
→ `dbpmr` `epsilon ≈ 0.21` (pelagic); benthic via `growth_detritivore`.

**Implication:** a matched **no-fishing, constant-forcing** `dbpmr` 0-D run can be
configured entirely from the crosswalk + these formulas, with temperature folded
into `A` — so Stage 0 can run today against the engine as-is.

## 6. Adapter prototype (can start now, no engine changes)

A pure-R adapter (using `arrow` for parquet, `jsonlite` for the reference JSONs)
can already: read one LME's forcing + reference params, build the matched
`dbpmr` `Setup.*` objects via this crosswalk, run `dbpmr`, and load the reference
equilibrium for comparison. This is the concrete Stage 0/Stage 2 deliverable and
needs nothing from the engine.

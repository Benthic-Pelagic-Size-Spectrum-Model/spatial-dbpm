# v2 LME calibration — two-group catchability (method note)

**One-line summary.** Per FAO-LME region, estimate two catchabilities — `q_pel` (pelagic/U spectrum)
and `q_ben` (benthic/V spectrum) — by fitting each to its **own** observed catch time series, using
each spectrum's **own** effort driver. No single-Q, no ΔAIC (those were v1).

## Model
The DBPM (`dbpmr` engine) is run per region with two coupled size spectra — **U** (pelagic
water-column predators, incl. demersal fish) and **V** (benthic detritivores) — plus a detritus pool.
Each region is spun to unfished equilibrium (40 yr at year-1 forcing), then driven by the obsclim
time series (plankton intercept/slope, surface & seafloor temperature, export ratio, effort),
1950–2010. Environmental forcing and the biomass-weighted plankton intercept are unchanged from the
current setup.

## Fishing (two catchabilities, two effort drivers)
Fishing mortality is applied **independently** to each spectrum:

    F_U(w,t) = q_pel · s_U(w,t) · E_U(t)
    F_V(w,t) = q_ben · s_V(w,t) · E_V(t)

- `q_pel`, `q_ben` — free parameters estimated per region, each acting **directly** on its own
  spectrum (the s_pel/s_ben biomass-share split used to divide a single Q is dropped).
- `s_U`, `s_V` — **min–max selectivity windows** (1 inside `[min,max]`, else 0), **region- and
  year-specific**, read from the parquet columns `min/max_fished_U`, `min/max_fished_V`. These are
  derived from the catch composition: min = lower edge of the smallest fished FGroup with a ~10 g
  gear-retention floor (small pelagics anchored at 5 cm ≈ 1.2 g; shrimp 5 g; krill the only sub-10 g
  exception at 1 g); max = Reg Watson's real `WtMax` for U (hybrid) and the invert-group max for V.
- `E_U`, `E_V` — **two effort drivers**. Nominal effort (`effort_histsoc_1841_2017_EEZ_addFAO`,
  `NomActive`) is split by the **same U/V FGroup map as the catch** (U = fish + krill + cephalopods;
  V = shrimp + lobster/crab + demersalmollusc), assigning effort by its **target** FGroup. Each series
  is normalised to its own maximum. (Effort enters via the R fishing-mortality function, so the two
  drivers need no engine change; per-gear efforts within a spectrum — e.g. bottom-trawl coupling —
  are a deferred extension.)

## Observed catch (two series)
The total observed catch `cat(t)` (parquet `catch_tonnes_area_m2 × 1e6`, **g m⁻² yr⁻¹**) is
partitioned into pelagic/benthic using the **observed group fraction** from `catch_split_lme.csv`
(catch_histsoc Reported+IUU, same U/V FGroup map):

    pel_frac(t) = pel_tonnes(t) / (pel_tonnes(t) + ben_tonnes(t))
    obs_pel(t)  = cat(t) · pel_frac(t)
    obs_ben(t)  = cat(t) · (1 − pel_frac(t))

Only the **ratio** is taken from the split file; the **level** is the parquet total. This keeps units
identical to the current single-series calibration and guarantees `obs_pel + obs_ben = cat` (avoids
absolute-level mismatch between the two catch aggregations).

## Objective
Modelled catch per group (standard `C = q·E·B`, linear in fishable biomass):

    C_pel(t) = q_pel · E_U(t) · B_U^fished(t)
    C_ben(t) = q_ben · E_V(t) · B_V^fished(t)

where `B^fished` is the spectrum biomass integrated over that year's fished window and converted to
areal density (pelagic × column depth capped at 200 m; benthic × benthic height). Estimate
`(q_pel, q_ben)` by minimising the **catch-volume-weighted two-series log-MSE**:

    J = Σ_t ω_pel(t)·(log10 C_pel − log10 obs_pel)²  +  Σ_t ω_ben(t)·(log10 C_ben − log10 obs_ben)²

with per-year weights `ω_g(t) ∝ obs_g(t)` (that group's observed catch, normalised within the group
to sum to 1) — so high-catch years dominate and noisy near-zero years are down-weighted. Optimisation:
coarse grid seed + `nloptr` BOBYQA, `q ∈ [0, QMAX=4]`.

## Assumptions
1. **One-group years** (a spectrum with no catch that year — e.g. no benthic fishery, the "grey"
   region-years): that year contributes only to the group with catch; the absent group's Q is left
   unidentified (reported NA), not forced to 0.
2. **Catchability is constant in time** per region (q_pel, q_ben are single scalars per LME).
3. **One effort per spectrum** (target-FGroup split); within-spectrum per-gear effort is deferred.
4. **No single-Q, no ΔAIC** — v2 is always the two-group fit (the old "two-Q vs one total series" was
   the overfitting case; splitting the observed series into two is the fix).

## Data / inputs
- `~/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs/*.parquet` — forcing + `min/max_fished_U/V` windows.
- `catch_split_lme.csv` — observed pelagic/benthic catch fraction per region-year.
- `effort_split_lme.csv` — **to be built** (mirror of catch split) — U/V effort per region-year.
- `lme_dint_all.csv` — biomass-weighted plankton intercept Δ (unchanged).

## Validation (separate)
Compare model **output** catch size spectra (slope/intercept) against Reg's observed `SSStats.csv`
per LME/year. Reg's slope/intercept are observed CATCH spectra, NOT model inputs.

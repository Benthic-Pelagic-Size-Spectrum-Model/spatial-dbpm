# 3-way comparison figures (dbpmr vs corrected CMIP5 vs corrected LME sizemodel)

Six representative LMEs — productive (1, 48), most-oligotrophic/marginal (64, 10),
mid (30, 44) — at `A.u = 64`, consumer-min −3, temperature on, with the `ln(10)`
fix applied to both sizemodel lineages (design/sizemodel-investigation.md §9).

- `fig_spectra*.png` — equilibrium predator size spectrum.
- `fig_timeseries*.png` — total predator biomass over the 60-yr spin-up.
- `*_abs.png` — **unit-reconciled, absolute** scale: LME-fixed ÷ `min(depth,200)`
  (its per-m² areal output, `useful_functions.R` ~l.795/813), dbpmr × ln(10)
  (per-ln → per-log10), CMIP5 raw. The non-`_abs` versions are shape-normalised.

**Findings.** After reconciliation, **CMIP5-fixed and LME-fixed overlay exactly**
(same model). **dbpmr agrees closely with the corrected sizemodel for productive
LMEs and at small sizes everywhere**, but **diverges at the large-size tail of
oligotrophic LMEs** (dbpmr sustains large fish; sizemodel's spectrum crashes) —
the residual engine differences (senescence functional form, numerics,
recruitment), not units or the ln(10). LME-fixed overflows to NaN on productive
LMEs at the weekly step (needs a finer timestep).

## All-82-LME set (paginated)

`fig_spectra_p1..4.png` / `fig_spectra_all83.pdf` — equilibrium predator size
spectra for all 82 LMEs (LME-32 has no searchvol_12.8 input), unit-reconciled
absolute scale. `fig_ts_p1..4.png` — matching biomass time series. Sizemodels run
at a 3-daily step (CMIP5-fixed) for stability; LME-fixed at weekly (parquet-locked)
so it is absent where it overflowed to NaN on productive LMEs.

**Across all 82:** CMIP5-fixed and LME-fixed overlay essentially everywhere (same
model, once depth-units reconciled). dbpmr overlays them across the recruit-to-mid
size range in most LMEs; the two sizemodels' spectra fall below dbpmr at the
**largest sizes** — modestly for most LMEs (invisible on the 100-decade axis) but
dramatically (crashing toward zero) for the deep, warm, oligotrophic subset
(e.g. LME-10, 12, 16, 19-21, 30, 31, 44, ...; median depth ~2670 m vs ~1110 m for
the closely-agreeing LMEs). That large-size residual is the senescence functional
form + recruitment boundary + numerics (issue #8), not units and not the ln(10).

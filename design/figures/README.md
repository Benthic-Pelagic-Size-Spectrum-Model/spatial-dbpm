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

# Plankton input aggregation — average-then-fit (method note)

**One-line summary.** Aggregate the plankton *fields* (phyc, phypico) horizontally across a region's
cells **first**, then fit **one** intercept & slope from the region-mean spectrum — instead of fitting
an intercept/slope per cell and then averaging the fitted values. This makes the plankton input
aggregation consistent with how catch and temperature are aggregated, and fixes the spectrum-collapse
that made the oligotrophic high-seas / heterogeneous LME calibrations anti-correlate.

## The problem (fit-then-average)
The current derivation is **fit-then-average**:
1. per cell: fit intercept `a_i` & slope from that cell's plankton spectrum;
2. region value = area-mean of the fitted intercepts, `sum(w_i·a_i)/sum(w_i)`, `w_i = cos(lat)`.

Because the intercept is `log10(biomass)`, averaging intercepts is a **geometric mean of biomass**. For a
spatially heterogeneous region (productive hotspots + oligotrophic majority) the geometric mean sits far
below the arithmetic mean, so the region looks far less productive than it is. In the 0-D model that
low mean productivity makes the size spectrum **collapse under fishing** → modelled pelagic catch falls
as effort rises → catch anti-correlates with observations (FAO 71 corr −0.92; also FAO 31/77/34/61,
Black Sea, East Brazil, Beaufort, Hudson, etc.). Confirmed: raising FAO 71's productivity ~+1 dex flips
corr −0.92 → +0.93.

## The consistency requirement (Julia)
**All spatial inputs, and the catch, must be aggregated the same way.**
- **Catch**: summed over cells, ÷ total area → `Σcatch/Σarea` = **arithmetic** area-mean density.
- **Temperature (tos, tob), export**: **arithmetic** area-mean (linear fields, averaged directly).
- **Plankton intercept**: currently area-mean of the **log** = **geometric** mean biomass — the lone
  inconsistency.

So the plankton aggregate must be an **arithmetic** (linear/biomass-space) mean too, to match. Note the
biomass-*weighted* form `Σb²/Σb` is the opposite error — it over-weights hotspots beyond their area
share, over-predicts, and is *also* inconsistent with `Σcatch/Σarea` (it numerically blew up Gulf of
Mexico at +1.6 dex). The consistent choice is the plain arithmetic mean of the plankton biomass.

## CRITICAL: intercept and slope MUST use the same averaging (fit them together)
The size spectrum is `N(m) = 10^intercept · exp(slope·m)` — intercept AND slope define it jointly. If the
intercept is average-then-fit but the slope is fit-then-average (or vice versa), the two do not
correspond to the same plankton spectrum and the spectrum is wrong. The only way to guarantee a matched
pair is to **fit both together from the same averaged fields**: average the (whole-column,
biomass-weighted) phyc/phypico across cells, then call `GetPPIntSlope` ONCE → intercept and slope come
out matched by construction. Do NOT correct only the intercept (e.g. via a DINT offset while leaving the
parquet slope) — that leaves them mismatched. Also: the vertical weighting is over the WHOLE water
column (`∫phyc²·dz / ∫phyc·dz`), not surface and not only 0-200 m (note the pipeline's
`integrating_phyto` currently uses `thresh_depth=200` — reconcile to whole-column).

## Corrected pipeline (average-then-fit)
Per region:
1. **Per cell — vertical** biomass-weighted integration of the plankton profile (unchanged; the
   `∫phyc²/∫phyc` vertical step) → each cell's representative `phyc_i`, `phypico_i`.
2. **Horizontal — arithmetic area-mean** of the plankton fields across the region's cells:
   `phyc_reg = Σ(w_i·phyc_i)/Σw_i`, `phypico_reg = Σ(w_i·phypico_i)/Σw_i`, `w_i = cos(lat)` (cell area).
3. **Fit once**: intercept & slope = `GetPPIntSlope(phyc_reg, phypico_reg)` — one region spectrum.
4. **Temperature, export**: keep the existing arithmetic area-mean (already consistent — do NOT
   biomass-weight them).

This is `average-then-fit`. It (a) gives the arithmetic-mean biomass that matches `Σcatch/Σarea`,
(b) derives the slope correctly from the region-mean spectrum (fit-then-average also mis-averaged the
slope), and (c) is exactly what the gridded model does per cell (no aggregation), so the two stay
consistent.

## Where to implement (upstream, like the fished-size hybrid)
- The region plankton intercept/slope come from the ISIMIP3a input processing (per-cell phyc/phypico →
  `GetPPIntSlope` → region parquet). Change that step from fit-then-average to average-then-fit.
- Local measurement/prototype: `all_int.R` (currently returns `int_bw = Σ(cos·int_cell)/Σcos`, geometric)
  → change to aggregate `phyc`/`phypico` arithmetically across cells, then fit. `int_lme<region>.csv`
  holds the per-cell values for a laptop redo; `all_dint.csv` holds the resulting per-region DINT.

## Also needed (calibration side)
The recovered (higher-productivity) regions need a **much lower q** to match their small catch density,
below the current optimiser grid seed (0.05). Fix the two-group optimiser to reach very low q — e.g.
optimise on `log10(q)`, or extend the seed grid down to ~1e-4 — so these regions fit their *level*
(their *shape* is already fixed by the aggregation change).

## Plan
1. Switch plankton aggregation to average-then-fit (arithmetic phyc/phypico mean → GetPPIntSlope);
   regenerate the per-region intercept/slope and DINT. Temperature/export unchanged.
2. Fix the optimiser to reach very low q.
3. Re-run the two-group calibration (all 83); check the ~17 weak regions recover and the 66 good LMEs
   are unchanged (their arithmetic mean ≈ geometric mean where cells are uniform, so little change).
4. Carry the same average-then-fit aggregation into the gridded inputs.

## Expected outcome
Recover most of the ~17 weak high-seas/heterogeneous regions (anti-correlation was a low-mean-
productivity artifact of geometric averaging), giving usable q's for the gridded model across all
regions — with a single, consistent aggregation rule for every spatial input.

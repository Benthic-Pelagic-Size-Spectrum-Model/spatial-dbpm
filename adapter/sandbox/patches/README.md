# Upstream-candidate patches (lme_scale_calibration_ISIMIP3a)

Proposed fixes to the LME workflow (the `lme-workflow/` submodule). Kept here as
tracked patches because they belong upstream and need a **preprocessing re-run on
the cluster** (the raw GFDL zarr data lives on NCI `/g/data/vf71/...`, not locally)
to regenerate the per-LME `dbpm_inputs` parquets before the effect shows in a run.

## `lme-workflow-mld-aware-phyto.patch`

`integrating_phyto()` averaged the phytoplankton concentration over a **fixed top
200 m**, which dilutes surface production in deep, productive columns (e.g. the
California Current / eastern-boundary upwelling) with phyto-poor water below the
euphotic zone. Result: the size-spectrum `intercept` is systematically ~10× lower
for deep than shallow LMEs (`corr(depth, intercept) = -0.45`), so deep LMEs look
oligotrophic and their pelagic under-produces (LME-3 California Current: Q pinned
at the [0,3] calibration bound, catch under-shooting).

The fix makes the averaging **depth** selectable via an `averaging` argument, and
reduces **both** size fractions (total `phyc`, small `phypico`) over the **same**
layer so the small/large ratio — which sets the spectrum **slope** in `GetPPIntSlope`
(`lphy = phyc − phypico`) — stays consistent:

- `fixed200` — legacy thickness-weighted mean over the top 200 m (depth-biased).
- `mld` — thickness-weighted mean over the mixed layer (`get_threshold_depth`);
  de-biases but truncates a deep chlorophyll maximum.
- `cumulative90` — mean over the layer holding 90 % of the column phyto biomass.
- **`biomass_weighted` (default, recommended)** — total-phyto-biomass-weighted mean
  concentration over the whole column, `∫phyc²/∫phyc` for total (and
  `∫phypico·phyc/∫phyc` for small) = *column biomass ÷ effective productive
  thickness*. Threshold- and light-free (no MLD field, no DCM truncation), and it's
  the density a **vertically-migrating, food-tracking pelagic forager** actually
  experiences — captures surface + DCM, ignores phyto-poor water it swims through
  but doesn't feed in. (CMIP5 used `min(depth,100)`; ISIMIP3a's fixed 200 m doubled
  the dilution.)

**Emulated effect** (intercept offset `+log10(200/MLD)` in `tier1_fishing_calib.R`
via `INT_OFFSET`, pending the real re-run): a `+0.3` shift (~MLD 100 m) on LME-3
takes Q off the bound (2.99 → 1.08), matches the observed catch (0.35 vs 0.29 g
m⁻² yr⁻¹), and restores a pelagic-dominated system (share 0.15 → 0.70).

### Input files required

The MLD path (`get_threshold_depth` + `integrating_phyto`) reads:

| variable | grid | status |
|---|---|---|
| **`mlotst-0125`** (mixed-layer depth, 0.125 kg m⁻³) | monthly | **NEW** — added to `dbpm_var` (was not downloaded before) |
| `phyc`, `phypico` (3-D depth-resolved phyto) | monthly | already pulled |
| `thkcello` (layer thickness) | fixed (grid_dir) | already pulled |
| `deptho` (bathymetry) | fixed/per-exp | already in `dbpm_var` (else = `thkcello.sum('lev')`) |

All are published on the UTAS THREDDS portal
(`.../obsclim/global/monthly/.../GFDL-MOM6-COBALT2/` has `mlotst-0125`, `phyc`,
`phypico`; `.../fixed/...` has `thkcello`). Only `mlotst-0125` is new — the patch
adds it to the download list so the workflow is self-contained.

Apply: `cd adapter/sandbox/lme-workflow && git apply ../patches/lme-workflow-mld-aware-phyto.patch`
(already applied in the working tree). Then re-run `00_processing_dbpm_global_inputs.py`
on the cluster (fetches `mlotst-0125`) to regenerate the parquets.

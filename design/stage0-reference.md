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

## 8. First prototype run (findings)

The harness `adapter/stage0_prototype.R` runs the full read → map → run → compare
pipeline for LME 10 (181 log10 bins, 200-year monthly spin-up). The forcing
parquet and reference JSON params match exactly (e.g. `intercept = int_phy_zoo =
-3.108`, `tos = sea_surf_temp = 24.74`), confirming the data path.

**Result: the literal crosswalk does not reproduce the LME equilibrium.** The
reference predator spectrum is non-zero across all 181 bins; the matched `dbpmr`
run **decays toward extinction** (final pelagic biomass ~1e-7). A search-rate
sensitivity sweep (`A × 1, 5, 25, 100`) does **not** recover it — biomass stays
~1e-7 up to ×25 and the run goes **numerically unstable (NaN)** at ×100.

**Follow-up diagnostics narrowed the cause:**

- **Plankton resource — RULED OUT.** `dbpmr` holds the plankton spectrum
  *constant* through the run (biomass 0.34217 unchanged over 40 yr), i.e. it is
  already a fixed boundary resource like `sizemodel()`. Predators collapse
  *despite* a constant resource, so starvation-by-resource-decay is not the cause.
- **Search rate `A` does not resolve it.** Low `A` → slow decay to extinction
  (pelagic biomass 0.21 → 1.6e-5 over 40 yr); raising the **benthic** `A` to 64
  (or dbpmr's defaults 640/64) makes the coupled run **numerically unstable
  (NaN)**. So `A = 64` is not the fix — the benthic search rate is in fact the
  instability driver.

**Remaining (the real Stage 0 work) — the predator energy / recruitment balance:**
- **Reproduction / recruitment boundary** — `dynamic_reproduction = 1` vs
  `rep_method = 2/3`: what sets the abundance flux at the smallest predator size
  (`xmin_consumer`) determines whether the spectrum self-sustains.
- **Growth vs mortality balance** — the `defecate_prop`/`growth_pred`/
  `energy_pred` → single `epsilon` mapping, and the senescence-mortality terms
  (`size_senescence`, `exp_senescence_mort`, `const_senescence_mort`) that
  `dbpmr` parameterises differently (`mu_s`).

**Decisive next step:** a line-by-line comparison of the growth/mortality/
reproduction rate equations between `dbpmr`'s C core (`calculate_g_and_mu`,
`g_pel`/`mu_pel`, `calculate_reproduction`) and the LME `sizemodel()`, to find
where the predator energy balance diverges. Tracked in issue #8.

### 8.1 Stable equilibrium reached (update)

Using `A_pel = 64`, `A_ben = 0.1 × 64 = 6.4` (canonical dbpmr-scale search rates,
per JB — not the literal `hr_volume_search × tempeffect`, which collapses) and a
**weekly** time step, `dbpmr` reaches a **stable, non-trivial equilibrium**
(pelagic biomass ~0.56, not extinction).

- **Time step matters:** `dbpmr` is numerically unstable at **monthly** steps for
  these LME-scale parameters (`NaN`), but stable at **weekly** and **daily**.
  Weekly is the FishMIP step the LME workflow already uses (issue #7).
- **Plankton** is confirmed held fixed at the equilibrium input through the run.

**Still open:** the equilibrium *shape* and *absolute scale* do not yet match the
reference — `dbpmr`'s pelagic spectrum is much **shallower** and offset by many
orders of magnitude (density normalisation/units + spectral-slope reconciliation:
growth-vs-mortality balance and the recruitment boundary). This is the remaining
Stage 0 work (issue #8).

### 8.2 Matched `A` + temperature correction (update)

Corrected two errors in the comparison: (a) run `dbpmr` at the **same** search
volume as the reference (`A_pel = hr_volume_search = 12.8`, `A_ben = 0.1·A_pel`),
not 64 — `sizemodel()` itself can't run at 64; and (b) apply the **temperature
effect as `sizemodel()` does** — the Boltzmann–Arrhenius factor scales both the
**feeding rate** and the **background ("other") mortality** (pelagic = surface
temp, benthic = floor temp; senescence and predation mortality are *not* scaled).
Folded into `dbpmr` as `A·tempeffect` and `mu_0·tempeffect`.

Findings (LME 10, `pel_te = 1.96`, `ben_te = 0.24`):

- **The reference is a benthic-only equilibrium.** At `searchvol = 12.8`,
  no fishing, the `sizemodel()` reference **predators are collapsed**
  (density ~1e-35 across the predator range) while the **detritivores are alive**
  (proper spectrum, slope ≈ −1/decade). The collapsed unfished pelagic looks
  like a problem on the `sizemodel()` side (the "issues" JB flagged), not a
  healthy target.
- **`dbpmr` sustains the pelagic** even with the full temperature correction
  (pelagic biomass ~0.87) — i.e. `dbpmr` does *not* reproduce the reference's
  predator collapse. So at matched `A` + temperature the engines genuinely
  disagree on pelagic persistence.
- **Benthic agrees in magnitude** (~1e3 at the recruitment boundary, both) but
  `dbpmr`'s detritivore spectrum is **steeper** (reaches a smaller maximum size)
  than the reference — a growth-vs-mortality slope difference.
- **Numerics:** `dbpmr` is `NaN` at monthly even at `A = 12.8`, stable at
  weekly/daily; `sizemodel()` is stable at monthly (#7).

**Implication:** the cleanest reconcilable target is the **benthic slope** (both
engines keep detritivores alive); the pelagic comparison is confounded by the
reference's collapsed predators and needs a trustworthy reference (a different
LME / search volume, or confirmation that the collapse is the known sizemodel
bug).

### 8.3 Canonical CMIP5 `sizemodel()` cross-check (`adapter/cmip5_reference.R`)

Ran the canonical DBPM `sizemodel()` (from `dbpm_isimip_3b`,
`dynamic_sizebased_model_functions_CMIP52019.R`) directly, with LME-10 constant
stable-spin forcing. Confirms canonical parameters: `A.u = 64`, `A.v = 0.1·A.u`,
default timestep `tstepspryr = 48` (≈ weekly), `alpha 0.82/0.75`, `mu0 = 0.2`,
depth-dependent `pref.ben = 0.8·e^(−depth/250)` (= 5.4e-8 at 4128 m, matching the
reference param exactly).

**Findings:**

- **The collapse is genuine, not a calibration artifact.** At canonical
  `A.u = 64`, the CMIP5 `sizemodel()` **also collapses the LME-10 predators**
  (max ~1e-31) while detritivores stay alive (max ~4.9). So both sizemodel
  versions agree the LME-10 unfished pelagic collapses — it is the model's real
  behaviour for this deep (4128 m), cold-floored (1.3 °C), low-export LME.
- **`A = 64` is unviable for matching:** the canonical model *crashes* at low
  search volumes for LME-10 unless the export ratio (`sinking_rate`, not
  `export_ratio`) is supplied; once supplied it runs but the predators collapse.
- **The two engines are near mirror images** at matched `A = 64` + temperature:
  - CMIP5 → **detritivore-dominated** (extended detritivore spectrum, predators
    collapsed);
  - `dbpmr` → **pelagic-dominated** (extended pelagic spectrum; benthos collapses
    once the cold-floor temperature ×0.24 is applied to `A_ben`).
  Strikingly, dbpmr's *pelagic* and CMIP5's *detritivore* have nearly the same
  shape and extent — each engine sustains the **opposite** functional group.

**Interpretation.** This is no longer a parameter-tuning gap: the engines
partition the **benthic–pelagic energy balance** in opposite directions. The
prime suspects are (a) **reproduction** — dbpmr `rep_method = 2` (R·intake,
default R = 0.2) vs the canonical egg-integral `R.u·biomass`; and (b) how each
applies the **benthic (cold-floor) temperature** and the **detritus → detritivore
→ predator coupling**. Reconciling this is a modelling decision (which partition
is correct for LME-10) for the domain experts — tracked in issue #8.

### 8.4 Resolution — dbpmr gives the correct healthy coexistence

**Domain call (JB):** the correct LME-10 behaviour is **both pelagic and benthic
thriving** — so the canonical `sizemodel()`'s collapsed pelagic is *wrong*, and
dbpmr's sustained pelagic is right.

> **Correction (see [`sizemodel-investigation.md`](sizemodel-investigation.md)):**
> the attribution to reproduction below is **superseded**. The sizemodel predator
> collapse is driven by **background mortality** at small sizes (lowering `mu0`
> restores the predators; even maximum reproduction does not), not by
> reproduction strength. dbpmr's persistence comes from a different small-size
> growth/mortality balance, still to be pinned down at the rate level.

**Reproduction (earlier hypothesis, now superseded).** The canonical reproductive
allocation is the *residual* assimilated energy:
```
R.u = (1-def.high)(1-(K.u+AM.u))·f.pel + (1-def.high)(1-(K.v+AM.v))·f.ben = 0.14·f.pel + 0.07·f.ben
R.v = (1-def.low )(1-(K.d+AM.v))·f.det                                    = 0.05·f.det
```
dbpmr `rep_method = 2` uses fixed efficiencies `R_pla/R_pel/R_ben = R_det = 0.2`
of intake — **higher**, so dbpmr invests more in recruitment and **sustains both
groups**, whereas the canonical residual scheme is leaner and collapses the
pelagic. (Naively setting dbpmr's `R` to the canonical values just destabilises
the run — the difference is structural: dbpmr keeps reproduction independent of
the growth/maintenance budget, the canonical makes it the residual of a conserved
budget.)

**Working config (healthy coexistence, no fishing, LME-10):**
`A_pel = 64·tempeffect`, `A_ben = 6.4·tempeffect`, weekly step, temperature on
(feeding + background mortality), **consumer-min = -3**, `rep_method = 2`
(default R), `epsilon = 0.21/0.14`. → pelagic biomass ~0.44, benthic ~1.71, both
with sensible spectra (`adapter/stage0_prototype.R`).

**Caveat — numerical fragility:** dbpmr `NaN`s at monthly steps, at
consumer-min = -7, and when reproduction is pushed up. consumer-min = -3 + weekly
is the stable operating point. The fragility (the implicit scheme's stability and
the text-file I/O) is the real engineering work for gridded scale — issues #5/#7.

## 6. Adapter prototype (can start now, no engine changes)

A pure-R adapter (using `arrow` for parquet, `jsonlite` for the reference JSONs)
can already: read one LME's forcing + reference params, build the matched
`dbpmr` `Setup.*` objects via this crosswalk, run `dbpmr`, and load the reference
equilibrium for comparison. This is the concrete Stage 0/Stage 2 deliverable and
needs nothing from the engine.

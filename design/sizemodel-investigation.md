# Investigation: why `sizemodel()` collapses the LME-10 predators

Companion to [`stage0-reference.md`](stage0-reference.md). Reports the diagnosis
of the pelagic-predator collapse seen in **both** `sizemodel()` lineages (the LME
calibration version that produced the reference equilibrium JSONs, and the
canonical `dbpm_isimip_3b` CMIP52019 version run directly here). Both collapse the
LME-10 unfished pelagic; the analysis below uses the CMIP5 source.

## 1. The collapse is a mortality-driven recruitment failure

Tracing the predator (consumer) biomass over a 150-yr LME-10 spin-up:

```
t=0.1 yr  0.0038      t=30 yr  1.6e-9
t=1.5 yr  0.0012      t=75 yr  1.1e-16
t=7.5 yr  3.1e-5      t=150 yr 5.2e-27   -> extinction
```

It is a smooth exponential decay. At the **recruit size** (log10 mass −3) at
equilibrium:

| rate | value (yr⁻¹) |
|---|---|
| growth `GG.u` | **0.44** |
| background mortality `tempeff·OM.u` | **2.2** |
| predation mortality `PM.u` | ~0 (predators already gone) |
| senescence | 0 |

**Recruits die ~5× faster than they grow.** They cannot grow up to maturity, so
the spectrum cannot sustain itself, and the few mature individuals' reproduction
declines with the population — a positive-feedback collapse.

## 2. It is NOT a reproduction problem

The canonical reproduction is the *residual* assimilated energy
(`R.u = (1-def)(1-(K+AM))·f`). Tripling it (set `AM.u = 0` → reproduction
fraction 0.14 → 0.49) **does not save the predators** — they still collapse.
Sensitivity (LME-10, temp on, consumer-min −3):

| change | predator |
|---|---|
| baseline (`mu0=0.2`, `AM.u=0.5`) | COLLAPSED |
| `AM.u = 0.2` (more reproduction) | COLLAPSED |
| `AM.u = 0.0` (max reproduction) | COLLAPSED |
| **`mu0 = 0.05`** (lower mortality) | **ALIVE (0.74)** |
| **`mu0 = 0.01`** | **ALIVE (5.96)** |

So the lever is **background mortality**, not reproduction. The default
`mu0 = 0.2` with `OM = mu0·w^−0.25` gives a high small-size mortality, and the
**warm surface temperature** (24.7 °C → `pel_tempeffect = 1.96`) *amplifies it*
(feeding is amplified by the same factor, but at small sizes the mortality term
dominates the balance). LME-10 is deep, warm-surfaced and low-productivity — the
combination makes the predators non-viable in this model.

## 2b. It is not a detritus → mortality cascade (causation is reversed)

Tested whether a growing detritus pool drives a detritivore boom that then spikes
mortality on the pelagic. Time trace (LME-10):

| yr | detritus `W` | detritivore bio | predator bio | predation-mort@recruit |
|---|---|---|---|---|
| 0.02 | 0.003 | 4e-5 | 0.0041 | 0.23 |
| 1.5 | 0.18 | 4e-5 | 0.0012 | 0.07 |
| 15 | 0.74 | 0.17 | 9e-7 | 3e-5 |
| 45 | 0.21 | 0.32 | 5e-12 | ~0 |

- The **predator decays from the first step**, before any detritus build-up; it
  never booms.
- Detritus accumulates and the **detritivore boom happens later** (after the
  pelagic is already gone).
- **Predation** mortality on recruits is small and *falling* (≤0.24/yr → 0);
  background mortality is a constant **2.2/yr** (~10× larger).

So causation is reversed: the pelagic collapses first (background mortality ≫
growth), then dead bodies + defecation feed the detritus, then detritivores boom
on it. The detritus/detritivore rise is a **consequence** of the collapse, not
its cause, and the preventive lever is `mu0`, not the detritus pathway.

## 3. The model is also numerically fragile

The detritus pool `W` readily goes `NaN`, crashing the run — e.g. with
`temp.effect = FALSE`, or with higher productivity (`pp = -1`). The guard that is
meant to catch this is a bug:

```r
if (W[i] == "NaN" | W[i] < 0)        # never matches a numeric NaN
```

`W[i] == "NaN"` compares a number to the *string* "NaN" (always `FALSE` / `NA`);
it should be `is.nan(W[i])`. So `NaN`s propagate instead of being caught, and the
loop later errors with "missing value where TRUE/FALSE needed". This affects both
the crash behaviour and any silent `NaN` contamination.

## 4. Why dbpmr does not collapse (correction)

An earlier note (`stage0-reference.md` §8.4) attributed dbpmr's surviving pelagic
to its **higher reproduction** (`rep_method = 2`, R≈0.2 vs the canonical 0.14).
**That was premature** — §2 above shows reproduction is not what drives the
sizemodel collapse. At the same nominal background mortality
(`mu_0·w^−0.25`, temperature-folded), dbpmr's predators nonetheless persist, so
the real difference is in the **small-size growth-vs-mortality balance** (dbpmr's
recruits survive long enough to grow), not reproduction. Pinning down exactly
which rate differs (dbpmr's effective growth higher, or its small-size mortality
lower) needs a direct rate-level comparison of `g_pel`/`mu_pel` at the recruit
size against `sizemodel()`'s `GG.u`/`Z.u` — the remaining Stage 0 task.

## 4b. Rate-level comparison at the recruit size (the 2.3× growth)

Extracted growth and mortality at the recruit (log10 −3), converting dbpmr's
`d ln m/dt` to `d log10 m/dt` (÷ln10):

| | growth (log10/yr) | total mortality (/yr) |
|---|---|---|
| **sizemodel** (early, collapsing) | 0.44 | 2.2 (mostly background) |
| **dbpmr** (sustained) | **0.996** (≈2.3×) | 13.3 (predation-dominated) |

dbpmr runs a **high-throughput** pelagic — ~2.3× faster growth, much higher
(predation) mortality, and correspondingly high reproduction — which sustains the
spectrum; sizemodel runs a **low-throughput** one (slow growth, lean reproduction)
that decays from the low initial density. So it is *not* a simple
growth/mortality ratio at the recruit (dbpmr's ratio is actually lower).

**Growth formulas (structurally identical, efficiencies near-equal):**
```
sizemodel:  GG.u = (1-def.high)*K.u*f.pel + (1-def.low)*K.v*f.ben   = 0.21*f.pel + 0.1*f.ben
dbpmr g_pel: K_pla*(1/w)*pla_bio + K_pel*(1/w)*pel_bio + K_ben*(1/w)*ben_bio   (K~0.2)
```
Both are `efficiency × mass-specific feeding`, with sizemodel effective 0.21 and
dbpmr default `K = 0.2`. The conversion efficiencies are near-equal, so the ~2.3×
sits in the **feeding integral** — see §4d, which pins it down exactly.

## 4d. The 2.3× growth IS ln(10) — but it is NOT why dbpmr survives (#23)

Reconciling the two intake calculations term-by-term, they are **structurally
identical**:

```
dbpmr   pla_bio = pref*A*w_pred^alpha * SUM_i[ exp(m_i) * phi(size-i) * u_i * mstep ]
        g_pel   = K*(1/w_pred)*pla_bio = K*pref*A*w_pred^(alpha-1) * SUM[ w_prey*phi*u*mstep ]

sizemodel f.pel = A*w_pred^alpha*pref * (U*dx) %*% gphi,   gphi(q) = 10^(-q)*phi(q) = (w_prey/w_pred)*phi
        GG.u    = (1-def)*K * f.pel  = (1-def)*K*pref*A*w_pred^(alpha-1) * SUM[ w_prey*phi*u*dx ]
```

dbpmr's explicit `exp(m_i)` prey-biomass weighting is exactly sizemodel's
`gphi = 10^(-q)·phi` factor; same `K`, `A`, `pref`, and `w_pred^(alpha-1)`
size-scaling. **The only difference is the integration step**: dbpmr sums over
**ln-mass** with `mstep = 0.1·ln10 = 0.2303`; sizemodel sums over **log10-mass**
with `dx = 0.1`. Their ratio is

```
mstep / dx = ln(10) = 2.3026
```

— exactly the measured growth ratio `0.996 / 0.44 = 2.26×` (§4b). So the 2.3× is
**`ln(10)`**, a log10-vs-ln feeding-integration convention difference: with the
plankton intercept `10^int` fed identically to both engines (Stage 0 mapping),
dbpmr's feeding integral over-counts the prey by `ln(10)` relative to sizemodel's.

**But this factor does not explain the survival difference.** Stripping it from
dbpmr the clean way — dividing `A` by `ln(10)` (it multiplies the *whole*
integral, so scaling plankton alone is insufficient) — **does not collapse the
pelagic**; it actually *rises*:

| dbpmr config (LME-10, weekly, mu0=0.2) | pelagic |
|---|---|
| `A = 64` (Stage 0) | ALIVE (32.4) |
| `A = 64/ln10 = 27.8` (ln10-corrected) | **ALIVE (50.2)** |
| plankton `10^int / ln10` | ALIVE (21.8) |

dbpmr is **robust to feeding strength** because its predation mortality (`death.u`)
scales with `A` too, so growth and death rebalance around a new equilibrium.
sizemodel is **fragile** because its collapse is driven by *background* mortality
(`mu0·w^−0.25·temp`), which is independent of `A`/feeding and so cannot be
rebalanced away.

**Conclusion (corrects §4b's implication).** The feeding kernel explains the *2.3×
growth number* (it is `ln(10)`), but **not** dbpmr's persistence. Survival is set
by the **mortality structure** (§2): dbpmr sits in a high-throughput,
predation-dominated equilibrium insensitive to the feeding scale, whereas
sizemodel's recruits are killed by A-independent background mortality. The lever
remains `mu0`, not the feeding kernel.

**Parameter-mapping note (Stage 0 setup):** the dbpmr runs set `epsilon = 0.21`,
but in dbpmr `epsilon` is the **senescence** constant — growth efficiency is
`K_pla/K_pel/K_ben`. The runs therefore used dbpmr's *default* `K = 0.2`
(≈ the canonical 0.21 by coincidence), not a deliberately-mapped value. A clean
reconciliation should set `K_*` explicitly and map `epsilon`/senescence properly.

## 4c. Should temperature force background mortality? (design question)

`sizemodel()` applies the Boltzmann factor to **feeding** (→ growth and, via
satiation, **predation mortality**) *and* explicitly to **background ("other")
mortality** (`Z.u = PM.u + pel.Tempeffect*OM.u + SM.u`); senescence is unscaled.

Test (LME-10): removing temperature from background mortality only — keeping it on
feeding/predation — raises the equilibrium predator from ~1e-26 to ~1e-16 but
**does not prevent the collapse**. So the explicit background-mortality scaling is
a *contributing* factor, not the root cause (the growth deficit of §4b dominates).

**Design recommendation (for dbpmr Stage 3 / #11).** Apply the Boltzmann factor to
**feeding only** — it propagates consistently to growth and predation mortality
(two sides of the same encounter process) — and leave **background mortality and
senescence temperature-independent**. Scaling background "other" mortality up at a
warm surface adds death with no compensating production, an asymmetry that pushes
the pelagic toward collapse; it is also harder to justify biologically (it is a
closure term, not a metabolic flux). dbpmr in Stage 0 folded temperature into both
`A` and `mu_0` and still thrived (because of the higher growth), but the cleaner
choice going forward is feeding-only.

## 5. Summary — "the problem with `sizemodel()`"

1. **Dynamical/parameter:** for warm-surfaced, low-productivity LMEs the
   temperature-amplified small-size background mortality (`mu0·w^−0.25·temp`)
   exceeds predator growth, so the pelagic is non-viable and decays to
   extinction. Reducing `mu0` (≈0.05) or otherwise rebalancing small-size
   growth/mortality restores it. Reproduction is not the cause.
2. **Numerical/code:** the detritus pool readily produces `NaN`, and the
   `W[i] == "NaN"` guard is a bug (should be `is.nan(W[i])`), so runs crash or
   silently corrupt instead of being caught.

Both `sizemodel()` lineages share this core dynamics, so both collapse the LME-10
pelagic. dbpmr gives the healthy coexistence that (per JB) is the correct
behaviour.

## 6. LME-wide sweep — it runs everywhere, collapse is warm + oligotrophic

Ran the canonical CMIP5 `sizemodel()` at **`A.u = 64`, consumer-min −3,
`temp.effect = TRUE`** across **all 82 LME** equilibrium inputs (100-yr spin).

- **No crashes.** Every run completed; *no* `NaN`/error at this configuration.
  So the `NaN` bug (§3) does **not** bite the A.u=64 LME sweep — the failures are
  clean dynamical collapses, not numerical ones.
- **70 ALIVE, 12 COLLAPSED.** The collapses are **not random**: every one is a
  **warm-surface, low-productivity** basin.

| collapsed LME | sst (°C) | pp (intercept) |
|---|---|---|
| 10, 12, 16, 30, 31, 40, 44, 131, 134, 151, 171, 177 | **all ≥ 22** | **all ≤ −2.3** |

Warm LMEs with adequate productivity survive (LME-35: sst 28.4, pp −0.39 →
alive); cold oligotrophic ones survive (LME-64: pp −3.71 → barely alive). Only
the **combination** warm + unproductive collapses — precisely the §2 mechanism
(temperature-amplified small-size background mortality outpacing growth where food
is scarce).

### 6a. dbpmr on the same sweep — all 82 alive

Running **dbpmr** across the same 82 LME inputs at the matching setup
(`A.u = 64`, consumer-min −3, weekly, temperature on feeding + `mu_0`, faithful
K/R/Ex mapping per §8) gives:

| engine | runs completed | pelagic alive | collapsed |
|---|---|---|---|
| `sizemodel()` | 82 (no crash) | 70 | **12** (warm + oligotrophic) |
| **dbpmr** | **82 (no crash, no NaN)** | **82** | **0** |

dbpmr sustains the pelagic **and** benthos in *every* LME — including all 12 that
`sizemodel()` collapses (10, 12, 16, 30, 31, 40, 44, 131, 134, 151, 171, 177) and
the most oligotrophic LME-64 (pp −3.71, pelagic 0.045, benthos 4.2). The result
is insensitive to the benthic search rate: `A.v = 6.4` (`0.1·A.u`, sizemodel
convention) and `A.v = 0.64` (`0.01·A.u`) both give all-82-alive. This is the
high-throughput, predation-dominated robustness of §4d: dbpmr does not carry the
A-independent background-mortality collapse, so warm + unproductive basins persist
rather than die — the coexistence (per JB) is the intended behaviour.

## 7. How to fix `sizemodel()` — two levers, tested on the 12

Re-ran the 12 collapsed LMEs under each candidate. Predator `max(x>−3)`:

| LME | baseline | `mu0 = 0.05` | feeding-only temp |
|---|---|---|---|
| 10  | DEAD | **alive 0.74** | DEAD 1.5e-11 |
| 12  | DEAD | alive 4.1 | alive 0.76 |
| 16  | DEAD | alive 3.9 | alive 0.08 |
| 30  | DEAD | alive 12  | alive 5.8 |
| 31  | DEAD | alive 8.4 | alive 4.3 |
| 40  | DEAD | alive 6.0 | alive 1.4 |
| 44  | DEAD | alive 8.1 | alive 0.004 |
| 131 | DEAD | **alive 0.92** | DEAD 2.6e-11 |
| 134 | DEAD | alive 11  | alive 4.0 |
| 151 | DEAD | alive 10  | alive 2.1 |
| 171 | DEAD | alive 2.7 | alive 0.04 |
| 177 | DEAD | alive 12  | alive 5.2 |

- **Lower background mortality (`mu0 ≈ 0.05`) rescues all 12** — the robust,
  reliable lever (consistent with §2). It is a *recalibration*, not a code fix.
- **Feeding-only temperature** (do not Boltzmann-amplify background mortality;
  the §4c recommendation) **rescues 10 of 12** with no recalibration; only the two
  most starved LMEs (10, 131, pp ≈ −3.1) still die. Biologically cleaner — it
  removes "death with no compensating production" — and gets most of the way.
- **`NaN` guard fix** (`W[i] == "NaN"` → `is.nan(W[i])`, §3) is *not* what the
  LME collapses need (no `NaN` at A.u=64), but it remains a genuine robustness bug
  for other forcings (`temp.effect = FALSE`, high `pp`).

**Recommendation.** For dbpmr's FishMIP path, adopt **feeding-only temperature**
(it fixes most LMEs structurally and is the cleaner closure) *and* expose `mu0`
(and/or a small-size mortality floor) for the few warm-oligotrophic LMEs that need
recalibration; fix the `is.nan` guard regardless. dbpmr already survives these
LMEs because its mortality balance is predation-dominated and robust to feeding
scale (§4d), so the same feeding-only-temperature choice is low-risk there.

## 8. Growth-efficiency / senescence mapping (#22)

The Stage 0 dbpmr crosswalk previously passed the growth efficiency into dbpmr's
`epsilon` argument — but **`epsilon` is the senescence size offset**
(`mu_s·(log10 w − log10 w_min)/((log10 w_max + epsilon) − log10 w)`), not an
efficiency. The growth efficiencies are `K_pla`/`K_pel`/`K_ben`, with companion
reproduction (`R_*`) and excretion (`Ex_*`) fractions. Each unit of intake splits
into **defecation + (K growth + R reproduction + Ex excretion)**, with
`K + R + Ex = 1 − defecation`.

Faithful mapping from the LME JSON (sizemodel keys the budget to the **prey
type** — pelagic/plankton use `def.high`, benthic/detritus use `def.low`):

| budget channel | formula | value |
|---|---|---|
| pelagic-prey `K_pel`/`R_pel`/`Ex_pel` | `(1−def.high)·{K.u, 1−(K.u+AM.u), AM.u}` | 0.21 / 0.14 / 0.35 |
| benthic-prey `K_ben`/`R_ben`/`Ex_ben` (and detritivore `K_det`/`R_det`/`Ex_det`) | `(1−def.low)·{K.v, 1−(K.v+AM.v), AM.v}` | 0.10 / 0.05 / 0.35 |

(`def.high=0.3`, `def.low=0.5`, `K.u=growth_pred=0.3`, `AM.u=energy_pred=0.5`,
`K.v=growth_detritivore=0.2`, `AM.v=energy_detritivore=0.7`.) Both channels close
to `1−def`. The dbpmr defaults (`R_pel=0.2`, `R_ben=0.2`) over-state reproduction
(faithful 0.14 / 0.05); `adapter/stage0_prototype.R` now sets `K_*`/`R_*`/`Ex_*`
explicitly and no longer mis-assigns `epsilon`. Setting the lower, faithful
reproduction does **not** collapse dbpmr (still healthy coexistence) — again
consistent with reproduction not being the survival lever (§4d).

**Senescence does not map 1:1.** sizemodel uses a power law
`SM = const_senescence·10^(exp_senescence·(x − size_senescence))`; dbpmr uses a
hyperbolic form diverging near `w_max` (`mu_s`, `epsilon`). These are different
functional forms, so senescence is left at dbpmr's default and flagged for the
engine-reconciliation work rather than force-fit through a fake parameter map.

## 9. Root cause: a missing `ln(10)` from sizemodel's log10 conversion

The §4d "2.3× growth = ln(10)" and §6a "dbpmr all-82-alive" findings have a single
root cause, and it is **not** the inputs (dbpmr runs on the *same* inputs and
survives): `sizemodel()` was rewritten from an **ln-mass** model (like dbpmr's
engine) into **log10-mass** form, and the rewrite kept the bare log10 step `dx`
(= 0.1) in the **integrals over the size spectrum**, where the ln-mass physics
needs `dx·ln(10)` (because `d(ln w) = ln(10)·d(log10 w)`). The *advection* term was
converted correctly — it carries the compensating `1/log(10)` Jacobian
(`(1/log(10))·GG·Δt/dx` in `Ai`/`Bi`) — but the **integral measures did not get
their `ln(10)`**, so every spectrum integral comes out `ln(10) = 2.3×` too small.

Two such integrals matter for the pelagic:

1. **Feeding / predation-death convolutions** `(U*dx) %*% gphi` and
   `(U*dx) %*% mphi` (cmip5 lines 234, 236, 251, 253, 269, 271; helper funcs
   95/97/100). → growth `GG`, reproduction rate `R`, and predation mortality `PM`
   are all `ln(10)` too small. Restoring it (feeding × ln10) **rescues 9 of the 12**.
2. **Reproduction egg-influx integral** `sum(R.u·10^x·U·dx)` (lines 392, 428).
   → recruitment is `ln(10)` too weak. Restoring it **rescues the remaining 3**
   (the hottest, most oligotrophic: LME-10, 131, 171).
   (A 3rd, the detritus-input sums at lines 317–329, has the same bare `dx` and
   should be corrected for benthic/detritus consistency, but does not bear on the
   pelagic collapse.)

**Demonstrated:** restoring *both* `ln(10)` factors makes **all 12** collapsed
LMEs survive at `A.u = 64` — i.e. `sizemodel()` reproduces dbpmr's all-LME
coexistence. So dbpmr "works" simply because its engine is in ln-mass and never
dropped the factor; `sizemodel()` "fails" on the warm + oligotrophic LMEs because
the missing `ln(10)` makes its growth and recruitment `2.3×` too weak — just
enough to push the marginal LMEs below the self-sustaining threshold
(survival-to-maturity `≈ exp(−∫ mortality/growth · d log w)` is exponentially
sensitive to that 2.3×). The other 70 LMEs sit far enough above threshold that
both engines agree.

### Summary of edits needed in `sizemodel()`

Express the spectrum-integration measure as the **ln step** wherever the code
integrates over size (it is a units fix, not a recalibration):

| location | now | should be |
|---|---|---|
| feeding `f.pel`/`f.ben` convolutions (l.234, 236) | `(U[,i]*dx) %*% gphi` | `(U[,i]*dx*log(10)) %*% gphi` |
| satiation denominators (l.251, 269) | `(U[,i]*dx) %*% gphi` | `… *dx*log(10)) %*% gphi` (cancels in the ratio — keep consistent) |
| predation death `PM.u`/`PM.v` (l.253, 271) | `(U[,i]*sat*dx) %*% mphi` | `(U[,i]*sat*dx*log(10)) %*% mphi` |
| helper `growth.u`/`death.u`/`death.v` (l.95–100) | `(u*dx) %*% kernel` | `(u*dx*log(10)) %*% kernel` |
| reproduction egg-sum (l.392, 428) | `sum(R.*·10^x·U·dx)` | `log(10)*sum(R.*·10^x·U·dx)` (numerator only; leave the `/(dx·10^x[ref])` bin-width normalisation as-is) |
| detritus input sums (l.317–329) | `sum(…·dx)` | `log(10)*sum(…·dx)` (consistency; benthic/detritus only) |

Equivalently, define `dx.int <- dx*log(10)` once and use it in every spectrum
integral, keeping the plain `dx` only in the advection (`Δt/dx`) and the recruit
bin-width normalisation. Leave the advection `1/log(10)` Jacobian untouched. The
maintainers should confirm the minimal patch (the feeding fix and the advection
`1/log(10)` interact — see §4d), but the diagnosis is unambiguous: the `2.3×`
deficit is a forgotten `ln(10)` on the spectrum integrals.

**Alternative that also rescues the last 3:** feeding-only temperature (§4c/§7) —
not a units fix but a modelling choice. The `ln(10)` restoration is the one that
makes `sizemodel()` *match dbpmr*.

### 9a. Both `sizemodel()` lineages share the bug

The LME calibration model
(`lme_scale_calibration_ISMIP3a`, branch `new_features`,
`scripts/useful_functions.R`) is the same engine refactored with descriptive
names (`predators`/`detritivores`/`log_size_increase`/`constant_growth`/
`constant_mortality` for `U`/`V`/`dx`/`gphi`/`mphi`). It carries the **identical**
omission — `log(10)` appears only in the advection/recruitment Jacobians
(lines 718, 720, 737, 751, 753, 777) and on **none** of the spectrum integrals:

| integral | CMIP5 line | LME `useful_functions.R` line |
|---|---|---|
| feeding → growth (pelagic / benthic) | 234 / 236 | 574 / 580 |
| predation death `PM_u` (and `PM_v`) | 253 / 271 | 612 (and sibling) |
| reproduction egg-sum (pred / det) | 392 / 428 | 732 / 774 |
| detritus-input sums | 317–329 | 691 (+ siblings) |
| advection `1/log(10)` (correct) | `Ai`/`Bi` | 718–720, 751–753 |

Both lineages descend from one ln-mass original converted to log10 with the same
dropped factor, which is why both collapse the LME-10 (and the other warm +
oligotrophic) pelagic identically. The §9 edit table applies verbatim, with the
line numbers remapped per the table above.

**Empirically closed on the LME code.** Patched `useful_functions.R` — restored
`log(10)` on the feeding/death convolutions (l.574, 580, 612) and the
reproduction egg-sums (l.732, 774) — and re-ran the repo's own
`run_model()`/`sizemodel()` with its real per-LME parquet forcing, `search_vol = 64`,
consumer-min −3, weekly, `temp_effect = TRUE`. **All 12 collapsed LMEs flip to
coexistence:**

| LME | baseline | fixed | LME | baseline | fixed |
|---|---|---|---|---|---|
| 10  | 1e-30 | **0.017** | 44  | 4e-23 | **8250** |
| 12  | 2e-33 | **3820**  | 131 | 1e-31 | **27.3** |
| 16  | 3e-32 | **3610**  | 134 | 2e-26 | **10500** |
| 30  | 2e-26 | **12500** | 151 | 5e-24 | **10000** |
| 31  | 1e-31 | **8280**  | 171 | 1e-35 | **2200** |
| 40  | 1e-29 | **5860**  | 177 | 6e-26 | **12000** |

(LME-10 and 131 stay low — they are the most oligotrophic, marginal in dbpmr too —
but alive, not extinct.) The single units fix turns `sizemodel()`'s 70/12 into
82/0, matching the ln-mass engine.

The **complete patch** is six spectrum integrals: the four convolutions
(feeding l.574, 580; predation death `PM_u` l.612 and `PM_v` l.644) and the two
reproduction egg-sums (l.732, 774). `PM_v` (the detritivore-side predation term,
on a wrapped line) does not affect the pelagic rescue but is included for
consistency; with all six fixed, **all 12 LMEs have both pelagic and benthos
alive and finite** (pelagic values unchanged from the 5-edit patch, as expected).

### 9b. Upstream PRs and the gridded model

The `ln(10)` fix is filed as PRs on both upstream repos (issues #1 / #13):

| repo | file / model | branch | PR |
|---|---|---|---|
| `dbpm_isimip_3b` | CMIP5 `sizemodel()` | `master` | #2 |
| `lme_scale_calibration_ISMIP3a` | `scripts/useful_functions.R` `sizemodel()` | `new_features` | #14 |
| `lme_scale_calibration_ISMIP3a` | `dbpm_model_functions.R` **`gridded_sizemodel()` + `sizemodel()`** | `main` | #15 |

**Gridded model / migration note.** The spatial `gridded_sizemodel()` (the
`predators[j, , i]` per-cell loop) currently exists **only on `main`
(`dbpm_model_functions.R`)** and carries the identical bug; PR #15 fixes it there
(21 lines, both functions in that file). The `new_features` refactor
(`useful_functions.R`) has the gridded *param/forcing* plumbing
(`sizeparam(gridded=)`, `run_model(gridded_forcing=)`) but its `sizemodel()` is
aspatial — the spatial model has **not yet been migrated** into the refactor.
When the team ports `gridded_sizemodel()` into `new_features/useful_functions.R`,
the same six spectrum-integral edits (feeding/death convolutions + reproduction
egg-sums, translated to the refactor names `log_size_increase`/`constant_growth`/
`constant_mortality`/`reprod_pred`/`reprod_det`) must be carried over. Until then,
PR #15 on `main` is the working gridded fix. Not empirically validated here (no
gridded forcing on hand) — recommend a baseline-vs-fixed gridded run before merge.
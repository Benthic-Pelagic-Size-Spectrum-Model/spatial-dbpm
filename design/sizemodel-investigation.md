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

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

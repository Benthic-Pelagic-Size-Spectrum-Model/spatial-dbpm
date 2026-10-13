# Reef refuge / structural-complexity — dbpmr Phase 1 design

Port the **size-structured prey-refuge (vulnerability) mechanism** from Rogers, Blanchard & Mumby
(2014, *Curr. Biol.*) and the two Rogers et al. (2018) papers into the compiled-C `dbpmr` engine as a
**localized, opt-in extension of the 2-spectrum (pelagic U + benthic V + detritus) model**. Phase 1
delivers the *refuge mechanism only*; the herbivore spectrum + algae pool are **Phase 2** (see end).

Source model: `reef_repos/Fisheries-and-degrading-reefs/sizemodel_functions.R` (`vulnerability()` L59–72;
applied L218–229, 260–269). Paper: vulnerability = "proportion of prey of a given weight vulnerable to
predation"; two complexity axes — **density** (depth of the function) and **diversity** (size-range width).

## 1. Mechanism

Per prey body-size class `w`, an **availability multiplier** `a(w)` ∈ (0,1] (0 = fully sheltered, 1 =
fully exposed). Two equivalent ways to specify it — the engine only needs `a(w)`:

**(A) Parametric — the 2014 paper (2 scalars, no crevice data):** a "vulnerability function" (a trench in
`w`) with just two parameters, fit by grid search to the observed high-complexity predator spectrum:
- `V_min` — minimum proportion of the population vulnerable to predation = **the floor** (paper grid `0.1–0.95`).
- `s` — **maximum refuge size**: the largest body size that finds shelter (paper grid `1–1500`); sets the
  width of the trench (the *diversity* axis). Depth toward `V_min` = the *density* axis.

**(B) Empirical — the 2018 papers/repos (per-size crevice data):**
```
a(w) = max( V_min , 1 − refuge(w) / competitor(w) )
```
- `refuge(w)` — measured refuge-space density at prey size `w` (count m⁻², e.g. `Refuge_data.txt`). Depth
  of the vector = refuge density; width (range of `w` with refuge>0) = refuge diversity — the same two axes.
- `competitor(w)` — abundance density of prey competing for those spaces (see §3).
- `V_min` (≡ the `min_A` param) — floor (default 0 ⇒ `refuge=0` ⇒ `a≡1` ⇒ standard dbpmr).

Design the engine around `a(w)`; support (B) first (a per-size input vector) and add the (A) generator
(`V_min`, `s`) → `refuge(w)` as a convenience for the 2014-style parametric fit.

`a(w)` multiplies the **same** feeding-kernel term on both sides so mass balance holds:
- **Predation mortality** on prey of size `w`: the predator→prey death integral is scaled by `a(w)`.
- **Predator intake/growth**: the prey contribution at size `w` to the predator's food is scaled by `a(w)`
  — sheltered prey are neither eaten nor gained as growth. This is the paper's "reduce mortality AND
  alter growth" property (emerges automatically from the shared factor).

## 2. Engine hooks (`dbpmr/src/SizeSpectra.c`)

The feeding kernel × search volume already drives both predation mortality and growth. Insert `a(w_prey)`
at the prey-size index of these two computations:

| term | function | call site | change |
|---|---|---|---|
| predation mortality | `mu_pel_pred(s,ss,i,k,l)` → `mu_pred_values[ss][i][k][l]` | L3146 | multiply the integrand by `a[prey-size index]` |
| benthic predation mortality | `mu_ben_pred(b,i,k,l)` | L3192 | (optional Phase 1) `a` on U→V predation |
| predator intake → growth | `g_pel(...)` (search `A` × kernel `phi` × prey biomass) | def ~L?, used in growth update | scale the prey-size contribution by the same `a` |

`a[]` is a per-prey-size array recomputed each timestep from the current spectrum + the refuge input
(exactly as the R loop recomputes `a.u`/`a.h` at the top of every step). It sits alongside the existing
`tempeff` scaling of `A` — same "per-timestep multiplier on the feeding term" pattern.

## 3. Design decisions for the 2-spectrum case

- **Competitor pool.** Reef model uses `U+H` (predators + herbivores share refuges). In 2-spectrum dbpmr
  the sheltering prey are the small end of the **pelagic spectrum U** (small U hide from large U within
  the same spectrum, and from U→V predation). Phase 1: `competitor(w) = U abundance density at w`
  (optionally `+V`). Documented, single choice.
- **Units.** `refuge` is a COUNT (holes m⁻²); `competitor` must be ABUNDANCE density (numbers), so convert
  the state variable to numbers per size (biomass/`w`) before the ratio. Nail exact abundance vs biomass
  convention against the engine's state definition at implementation.
- **Benthic V.** Reef model gives inverts no refuge (`a.v=1`). Phase 1 default: refuge applies to
  predation on U prey only; U→V refuge is an opt-in flag.

## 4. Input delivery + R API

Reuse the forcing-plugin architecture (the `forcing_ts.txt` pattern already added):
- New **refuge channel** = a per-size vector `refuge(w)` (static, or per-timestep for a degradation
  trajectory). Simplest: a `refuge.txt` (col 1 = log10 size grid, col 2 = refuge density) analogous to
  the repos' `Refuge_data.txt`; or a `refuge` column in the forcing table.
- New scalar param **`min_A`** (env `REFUGE_MIN_A` or a Setup arg), default 0 (off ⇒ engine unchanged).
- R: `Setup.Refuge(pe, run, grid, refuge_by_size, min_A=0)` writer, mirroring `Setup.Temperature`/`Setup.ts`.
- **Off-by-default guarantee:** no refuge input ⇒ `a≡1` ⇒ bit-identical to current dbpmr (regression-safe).

## 5. Validation (2-spectrum, before any Phase 2)

1. **Null check:** refuge off ⇒ output identical to current engine (exact).
2. **Mechanism:** a size-structured refuge (density at small/mid sizes) ⇒ reduced predation mortality on
   those sizes ⇒ **non-linear size spectrum** (local abundance bump) — reproduce the paper's Fig-2/Fig-4
   qualitative shape and the repos' `a.u` behaviour on a matched grid.
3. **Sensitivity:** shrink refuge density+diversity (high→low complexity) and check the direction/rough
   magnitude of the productivity drop (paper: >3× for full complexity loss; we expect a strong monotone
   decline in fishable-size production, not the exact Caribbean number in a 2-spectrum ocean model).
4. Cross-check `a(w)` and the mortality field against a direct R re-run of `vulnerability()` on the same
   refuge + spectrum.
5. **Optional calibration replicate (2014 method):** grid-search `(V_min ∈ 0.1–0.95, s ∈ 1–1500)` to
   minimise RSS against a target predator size spectrum — reproduces the paper's fitting step and shows
   whether the refuge term improves spectrum fit (ΔRSS, ΔR²) over the no-refuge base, as in the paper.

## 6. Effort, risk, scope

- **Effort:** small — one `a[]` array + two multiplier insertions + a refuge input channel + `min_A` +
  R writer + tests. Comparable to the temperature/forcing-plugin work already landed.
- **Risk:** low and contained — opt-in, off-by-default, reuses the plugin pattern; main care points are
  the competitor-pool/units convention and applying `a` symmetrically (mortality *and* growth) to keep
  mass balance.
- **Delivers:** the core 2014 result (complexity → size-structured predation refuge → non-linear spectra
  and productivity sensitivity) inside dbpmr, as an honest reduced form of the reef model.

## 7. Phase 2 — herbivores + algae, via the existing multi-species engine (NOT a rewrite)

**Honest caveat first:** even the 2014 base model was DBPM **+ a herbivore spectrum** (before the
vulnerability function), so Phase-1 (2-spectrum U/V + refuge) is a *reduced form* — it isolates the
refuge mechanism but with a U-only competitor pool rather than U+H.

**Key realisation (revises the earlier "structural rewrite" framing).** In the reef repos the extra
groups are literally **copies** of existing state variables: `H`/`V`/`U` are the *same* McKendrick–von
Foerster advection–reaction solve repeated (diffs are whitespace/comments), and algae `A` is a scalar
ODE structurally identical to detritus `W`. dbpmr is **already multi-species** — `setup_pelagic`/
`setup_benthic` run in a loop over species (`community.pelagic[s]`, `benthic[b]`; group-vs-group
predation `mu_pel_pred(s,ss,...)`). So the "extra spectrum" is **not** new C machinery — it is **one
more `Setup.Benthic`/`Setup.Pelagic` instance** (an R-side config) with its feeding preferences pointed
at the right food. Herbivores = a detritivore *copy* grazing algae; predators are already present and
eat them via the existing `pref_ben`/`pref_pel` coupling.

**The ONLY genuine structural change: pluralise the resource pool.** dbpmr has a *single* detritus pool
`W` (and a single, prescribed plankton). Algae is "a copy of detritus" with a different input, so the
pool machinery must go from one → N: add a second detritus-like pool with an **autotrophic logistic
input** (`alr`→capacity `AK`, from the reef code) instead of the dead-body/sinking input, and let a
consumer group graze it (a `K_alg`/`pref_alg` link to the new pool).

**Phase-2 task list (small, mostly config):**
1. **Pluralise `DETRITUS`** → support ≥2 pools; add an algae pool with logistic input + grazing loss. *(the one real C change)*
2. **Instantiate a herbivore species** — an extra `Setup.Benthic` (detritivore copy) feeding on the algae pool. *(config)*
3. **Predators eat herbivores** via the existing `pref_ben` (or `pref_pel`) coupling. *(likely free — verify §8)*
4. **Extend the refuge competitor pool to `U+H`** (Phase-1 `a(w)` uses predators + herbivores). *(tie-in)*

**Staging:**
- **Phase 2a (cheapest):** herbivore species + a *prescribed* macroalgae resource (plankton-style, fixed).
  Gets herbivores into the food web and the `U+H` competitor pool. Algae is exogenous — no degradation
  feedback, but nothing new structural beyond a second prescribed resource.
- **Phase 2b (full):** promote algae to a *dynamic* second pool (detritus-style + logistic input +
  grazing depletion) → unlocks the algae↔herbivore feedback the 2018 degradation / predator-overfishing
  results depend on.

Because every group stays on the *same* DBPM code path (just more instances), this keeps the reef
model low-divergence — it is "more `Setup` calls + one plural resource pool", not a fork of the engine.

## 8. Open items to verify at implementation

- **Cross-species predation coupling:** confirm `pref_pel`/`pref_ben` let a predator species eat *all*
  species of that type (so predators automatically eat the new herbivore group), vs a single-target link.
- **State units** (abundance vs biomass density) for the refuge `competitor(w)` ratio (§3).
- **Consumer→pool feeding**: current benthic feeds one pool (detritus); confirm the cleanest way to point
  a consumer at a *second* pool (and, for herbivores, at algae **and** detritus as in the reef code).
- **Recruitment/reproduction** of the herbivore instance reuses the existing spectrum `rep_method`.

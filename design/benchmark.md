# Benchmark: dbpmr vs the two `sizemodel()` engines

Wall-clock for one 0-D run, **matched workload: 2400 timesteps × 181 size bins,
all at a WEEKLY step (48/yr = 50 simulated years)**, LME-10 forcing (pelagic +
benthic + detritus + plankton, no fishing). 4 reps.

| Model | language | step | best (min) | mean | ~ms/timestep |
|---|---|---|---|---|---|
| **dbpmr** | C (`.C` engine) + text I/O | weekly (`tstep=1/48`) | **1.58 s** | 2.42 s | 0.7–1.2 |
| **LME `sizemodel`** | R (`lme_scale_calibration`) | weekly | 1.80 s | 1.89 s | ~0.75 |
| **CMIP5 `sizemodel`** | R (`dbpm_isimip_3b`) | weekly (`tstepspryr=48`) | 2.24 s | 2.33 s | ~0.93 |

Note: the LME `sizemodel()` takes its step from the input data resolution (the
parquet is **monthly**), so it must be fed weekly-resolution input to run weekly.
Its per-step cost is the same either way (monthly run: 1.79 s — identical), since
each timestep does the same 181-bin work regardless of the step *duration*. Step
*count* (2400) is what is matched here.

## Finding

**The three are comparable (~2 s); the compiled C engine is *not* meaningfully
faster than the vectorised R models for a single 0-D run.**

Why C doesn't win here:
- The R models are **vectorised** — each timestep is a handful of length-181
  vector ops and 181×181 kernel convolutions, which R runs through BLAS at near-C
  speed.
- dbpmr's C core uses **scalar nested loops** over size bins, and carries the
  full **spatial-movement / coupled-solver machinery and per-run text-file I/O**
  even for a 0-D cell. dbpmr's best run is fastest, but its variance is highest
  (the file I/O).

## Implications for FishMIP scaling

Raw 0-D engine speed is **not** the differentiator — all three are ~2 s/cell. At
gridded scale (thousands of cells × centuries) the cost is dominated by:

1. **How cells are iterated.** Looping cells in R and calling the engine per cell
   (the dbpmr or per-cell R approach) pays the ~2 s × N_cells plus, for dbpmr,
   the text-file I/O per cell. The production gridded path instead **vectorises
   across cells** (the Python `gridded_sizemodel_rk4`, xarray) — one set of vector
   ops covering all cells per timestep — which is the scalable design.
2. **I/O.** dbpmr's per-run text files are the main avoidable overhead and the
   first thing to fix for scale (issue #5: in-memory `.Call` return).

So dbpmr's value is its **modular S4 framework, spatial movement, and (once the
I/O is in memory) a fast compiled core** — not a raw single-cell speed advantage.
The performance work that matters is the in-memory return (#5) and a
vectorised-across-cells gridded driver, not micro-optimising the 0-D core.

*(Caveat: timings are a single machine, single 0-D cell, modest reps; treat as
order-of-magnitude. The point is the ~equal ballpark, not the exact ordering.)*

# Spatial movement solver (xmove_solver) — regression tests for issue #41.
# The fixed solver is operator-split: conservative upwind FLUX advection + Crank-Nicolson
# diffusion, both with zero-flux walls (mass-conserving, stable, symmetric).
LN10 <- log(10)

# final-time pelagic biomass per x-cell for a 1D spatial run
.move_biomass <- function(comp, dx, prey = 0, pred = 0, seed_cell = 2, ncell = 5, tmax = 3,
                          uniform = FALSE) {
  wd <- tempfile("mv"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 1, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = 1,
                     xmin = 0, xmax = (ncell - 1) * dx, xstep = dx, xoutstep = dx)
  # uniform=TRUE: also give plankton a uniform IC + ts, overriding the default spatial x_start
  # Gaussian (the real source of cell-to-cell divergence diagnosed in issue #42).
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = uniform)
  if (uniform) Setup.ts(pl, run, grid, func = function(m, t, x, y) 0.01 * exp(-m))
  pe <- Setup.Pelagic(run, filename = "fish", comp = comp, prey = prey, pred = pred, initial_flag = TRUE)
  seedx <- seed_cell * dx
  f <- if (uniform) function(m, t, x, y) (0.01/LN10) * exp(-m)
       else         function(m, t, x, y) if (abs(x - seedx) < dx/2) (0.01/LN10)*exp(-m) else 1e-20
  Setup.ts(pe, run, grid, func = f)
  # benthic also carries the x_start bump unless given a uniform IC (same reason as plankton)
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = uniform)
  if (uniform) Setup.ts(be, run, grid, func = function(m, t, x, y) 0.01 * exp(-m))
  de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  L <- readLines("R/fish/results.txt")
  hdr <- suppressWarnings(as.numeric(strsplit(L[5], ",")[[1]])); mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
  rows <- lapply(L[6:length(L)][nchar(L[6:length(L)]) > 0],
                 function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]])))
  rows <- rows[sapply(rows, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]))]
  tout <- max(sapply(rows, `[`, 1)); last <- rows[sapply(rows, `[`, 1) == tout]
  b <- sapply(last, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE)); xx <- sapply(last, `[`, 2)
  xo <- sort(unique(xx)); unname(sapply(xo, function(q) b[which(xx == q)[1]]))  # biomass per cell, ordered by x
}

test_that("diffusion spreads symmetrically from a central source (issue #41)", {
  b <- .move_biomass(comp = 1e5, dx = 1, seed_cell = 2, ncell = 5)   # strong diffusion, seed middle cell
  expect_equal(b[2], b[4], tolerance = 0.02)   # cells either side of the seed match
  expect_equal(b[1], b[5], tolerance = 0.02)   # the two ends match
  expect_gte(b[3] + 1e-12, b[2])               # peak at (>=) the seeded middle, NOT piled at a wall
  expect_false(b[5] > 1.5 * b[3])              # the old bug: biomass piling at x-max
})

test_that("with movement off and uniform inputs, identical cells stay identical (#42)", {
  # #42 closed as not-a-bug: the per-cell divergence was the default spatial plankton x_start
  # bump, not the solver. With plankton + pelagic both uniform and movement off, every cell is
  # an independent, identical 0-D run (and the connectivity floor is now per-cell too).
  b <- .move_biomass(comp = 0, dx = 100, uniform = TRUE, ncell = 5)
  expect_lt(diff(range(b)) / mean(b), 1e-6)
})

# Fishing mortality passed to the C code as an array (issue #58, step 2) instead of a
# <run>/Input/<species>_fishing_ts.txt file. Both routes must give the same run.

LO <- log(10) + 0.05                                   # fished above this size (between two grid points)
yr <- function(t) pmin(floor(t + 1e-9), 2) + 1         # year (1..3) at model time t, for a 3-year run

# rate: [x, year] matrix of fishing mortality for the pelagic species; the benthic species is fished at half of it.
# Returns the output arrays; `how` selects the input route for the same fishing.
run_fished <- function(rate, how = c("file", "array"), assign = TRUE, both = FALSE) {
  how <- match.arg(how)
  wd <- tempfile("fish"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  nx <- nrow(rate)
  run  <- Setup.Run("R", 1, 1, spatial_dim = if (nx > 1) 1 else 0, coupled_flag = TRUE, diff_method = 1)
  grid <- if (nx > 1) Setup.Grid(run, tmax = 3, tstep = 1/48, toutstep = 0.5, xmin = 0, xmax = nx - 1, xstep = 1, xoutstep = 1)
          else Setup.Grid(run, tmax = 3, tstep = 1/48, toutstep = 0.5)
  xs <- dbpmr:::grid_seq(grid@xmin, grid@xmax, grid@xstep)
  uni <- function(sp) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1) }
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE); Setup.ts(pl, run, grid, mat = uni(pl))
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE, fishing_flag = TRUE); Setup.ts(pe, run, grid, mat = uni(pe))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE, fishing_flag = TRUE); Setup.ts(be, run, grid, mat = uni(be))
  de <- Setup.Detritus(run, filename = "detritus")
  ffun <- function(scale) function(m, t, x, y) if (m >= LO) scale * rate[which.min(abs(xs - x)), yr(t)] else 0
  if (how == "file" || both) { Setup.fishing(pe, run, grid, func = ffun(1)); Setup.fishing(be, run, grid, func = ffun(0.5)) }
  if (how == "array") {
    a <- array(rate, c(nx, 1, 3))
    p2 <- Setup.fishing(pe, run, grid, rate = a, min_size = LO, time = 0:2)
    b2 <- Setup.fishing(be, run, grid, rate = 0.5 * a, min_size = LO, layer = yr(seq(0, 3, 1/48)))
    if (assign) { pe <- p2; be <- b2 }
  }
  out <- NULL
  invisible(capture.output(out <- SizeSpectrum(run, grid, pl, pe, be, de, output = "arrays")))
  attr(out, "wrote_file") <- file.exists(file.path("R", "Input", "fish_fishing_ts.txt"))
  out
}
rates <- function(nx) outer(seq(1, 1.5, length.out = nx), c(0, 0.2, 0.4))    # unfished first year, then rising

test_that("a fishing rate array gives the same 0-D run as the fishing file, without writing one", {
  f <- run_fished(rates(1), "file"); r <- run_fished(rates(1), "array")
  expect_true(attr(f, "wrote_file")); expect_false(attr(r, "wrote_file"))
  for (v in c("abundance", "fishing", "catch", "mortality"))
    expect_equal(r[[v]], f[[v]], tolerance = 1e-12, info = v)
  # fished only above the size limit, and only from the second year
  expect_equal(max(r$fishing$fish[, , r$m < LO, ]), 0)
  expect_equal(unique(as.vector(r$fishing$fish[1, 1, r$m >= LO & r$m <= 6 * log(10), ])), c(0, 0.2, 0.4))
  expect_equal(unique(as.vector(r$fishing$benthos[1, 1, r$m >= LO & r$m <= 4 * log(10), ])), c(0, 0.1, 0.2))
  expect_gt(max(r$catch$fish), 0); expect_equal(max(r$catch$fish[, , , r$t < 1]), 0)
})

test_that("a per-cell fishing rate array gives the same spatial run as the fishing file", {
  f <- run_fished(rates(3), "file"); r <- run_fished(rates(3), "array")
  for (v in c("abundance", "fishing", "catch"))
    expect_equal(r[[v]], f[[v]], tolerance = 1e-12, info = v)
  last <- length(r$t); big <- which(r$m >= LO & r$m <= 6 * log(10))[1]
  expect_equal(r$fishing$fish[, 1, big, last], c(0.4, 0.5, 0.6))                 # the cells differ as supplied
  expect_true(all(diff(apply(r$abundance$fish[, 1, r$m >= LO, last], 1, sum)) < 0))   # more fishing, less fish
})

test_that("forgetting to assign, or supplying both a rate and a file, is an error", {
  expect_error(run_fished(rates(1), "array", assign = FALSE), "assign the result")
  expect_error(run_fished(rates(1), "array", both = TRUE), "use one or the other")
})

test_that("Setup.fishing validates a rate", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1, xmin = 0, xmax = 2, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  pe <- Setup.Pelagic(run, filename = "fish", fishing_flag = TRUE)
  out <- Setup.fishing(pe, run, grid, rate = matrix(0.2, 3, 2), min_size = log(100))
  expect_s4_class(out, "pelagic.params")
  expect_equal(dim(out@arrays$fishing$rate), c(3L, 2L, 1L)); expect_equal(out@arrays$fishing$range, c(log(100), Inf))
  expect_equal(dim(Setup.fishing(pe, run, grid, rate = 0.3)@arrays$fishing$rate), c(3L, 2L, 1L))
  expect_equal(dim(Setup.fishing(pe, run, grid, rate = c(0, 0.3), time = c(0, 0.5))@arrays$fishing$rate), c(3L, 2L, 2L))
  expect_error(Setup.fishing(pe, run, grid, rate = matrix(0.2, 2, 3)), "x = 3, y = 2")
  expect_error(Setup.fishing(pe, run, grid, rate = array(0.2, c(3, 2, 4))), "4 layers but the run has 49 time points")
  expect_error(Setup.fishing(pe, run, grid, rate = -0.1), "non-negative")
  expect_error(Setup.fishing(pe, run, grid, rate = NA_real_), "NA, NaN or Inf")
  expect_error(Setup.fishing(pe, run, grid, rate = 0.2, min_size = 3, max_size = 1), "min_size <= max_size")
  expect_error(Setup.fishing(pe, run, grid, rate = 0.2, func = function(m, t, x, y) 0), "cannot be combined")
  expect_error(Setup.fishing(pe, run, grid, func = function(m, t, x, y) 0, time = 0), "only used with rate")
  expect_error(Setup.fishing(Setup.Pelagic(run, filename = "fish"), run, grid, rate = 0.2), "fishing_flag")
})

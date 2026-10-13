# Reef prey-refuge mechanism (Setup.Refuge + <name>_refuge.txt + REFUGE_MIN_A)
# -----------------------------------------------------------------------------
# a(w) = max(min_A, 1 - refuge(w)/competitor(w)) scales predation mortality (at the
# predator index, matching Rogers 2014 death.u) and predator intake (at the prey index).
# Off (no refuge / all-zero) => a(w)=1 => bit-identical to standard dbpmr. A size-structured
# refuge reduces predation on the sheltered sizes and raises their standing abundance (bump).
LN10 <- log(10)

# run a minimal unfished pelagic+benthic model; optional refuge via Setup.Refuge.
# returns the final pelagic abundance spectrum (named by log10 body mass).
run_spectrum <- function(refuge = NULL, min_A = 0) {
  on.exit(Sys.unsetenv("REFUGE_MIN_A"), add = TRUE)
  wd <- tempfile("reftest"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = 20, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton")
  pe <- Setup.Pelagic(run, filename = "fish")
  be <- Setup.Benthic(run, filename = "benthos")
  de <- Setup.Detritus(run, filename = "detritus")
  if (!is.null(refuge)) Setup.Refuge(pe, run, grid, refuge, min_A = min_A)
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  f <- Read.In("R", "fish")
  u <- as.numeric(f@finaluvals[1, -(1:3)])
  names(u) <- f@mrange / LN10   # log10 body mass
  u
}

test_that("no refuge is bit-identical to an all-zero refuge (off-by-default)", {
  mass_n <- length(run_spectrum())                       # discover grid length
  base   <- run_spectrum(refuge = NULL)
  zero   <- run_spectrum(refuge = rep(0, mass_n), min_A = 0)
  expect_equal(base, zero)                               # all-zero refuge must not perturb anything
})

test_that("a size-structured refuge raises abundance at the sheltered sizes (bump)", {
  mass_n <- length(run_spectrum())
  log10m <- as.numeric(names(run_spectrum()))
  band   <- log10m >= 0 & log10m <= 1.5                  # shelter ~1-30 g prey
  refuge <- rep(0, mass_n); refuge[band] <- 1e20          # strong refuge -> a clamps to min_A on the band
  off <- run_spectrum(refuge = NULL)
  on  <- run_spectrum(refuge = refuge, min_A = 0.2)
  expect_false(isTRUE(all.equal(off, on)))               # mechanism actually changes the spectrum
  # sheltered sizes should be MORE abundant with the refuge than without (Rogers 2014 direction)
  expect_gt(mean(on[band]) / mean(off[band]), 1)
})

test_that("Setup.Refuge validates inputs", {
  wd <- tempfile("refval"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = 1)
  pe   <- Setup.Pelagic(run, filename = "fish")
  expect_error(Setup.Refuge(pe, run, grid, refuge = c(1, 2, 3)), "one value per mass class")
  expect_error(Setup.Refuge(pe, run, grid, refuge = function(m) 0, min_A = 2), "min_A")
})

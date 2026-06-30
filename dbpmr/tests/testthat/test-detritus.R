# Regression test for the detritus-pool NaN guard (SizeSpectra.c).
# The detritus pool W is advanced with an explicit Euler step
# W <- W + dt*(g - mu); a coarse timestep can overshoot (mu*dt > W) and drive W
# negative, which without a guard propagates to NaN and crashes the coupled run
# (the analogue of the reference sizemodel's `W == "NaN"` string-compare guard
# that never matched). The engine now clamps a non-finite or negative W to 0, so
# a coupled pelagic+benthic+detritus run stays finite at a coarse step where it
# previously went NaN. The parameters below (warm-LME-like: elevated search rate
# and background mortality, fortnightly step) are a config that goes NaN without
# the guard and finite with it.

test_that("coupled run with a coarse timestep stays finite (detritus guard)", {
  old <- setwd(tempdir())
  on.exit({ unlink("det_run", recursive = TRUE); setwd(old) })
  LN10 <- log(10)

  run  <- Setup.Run("det_run", no_pelagic = 1, no_benthic = 1,
                    spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  # log10 grid (the FishMIP convention) so the consumer minimum lands at 10^-3 g.
  grid <- Setup.Grid(run, mmin = -12 * LN10, mmax = 6 * LN10,
                     mstep = 0.1 * LN10, moutstep = 0.1 * LN10,
                     tmax = 60, tstep = 1/24, toutstep = 1)        # fortnightly step
  pl <- Setup.Plankton(run, filename = "plankton",
                       mmin = -12 * LN10, mmax = -3 * LN10, u_0 = 1e-3, lambda = -1.15)
  pe <- Setup.Pelagic(run, filename = "fish",
                      mmin = -3 * LN10, mmax = 6 * LN10, A = 125, mu_0 = 0.4,
                      rep_method = 2)
  be <- Setup.Benthic(run, filename = "benthos",
                      mmin = -3 * LN10, mmax = 4 * LN10, A = 8, mu_0 = 0.26,
                      rep_method = 2)
  de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))

  fish <- Read.In("det_run", "fish")
  ben  <- Read.In("det_run", "benthos")
  expect_true(all(is.finite(as.numeric(fish@finaluvals[1, -(1:3)]))))
  expect_true(all(is.finite(as.numeric(ben@finaluvals[1, -(1:3)]))))
})

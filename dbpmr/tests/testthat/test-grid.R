# Regression tests for the grid floating-point off-by-one (issue #21):
# an irrational step (a log10 grid expressed in natural log) must not drop the
# top size bin, and Read.In's mass axis must align with the data columns.

test_that("grid_seq returns round((to-from)/by)+1 bins and restores dropped ones", {
  LN10 <- log(10)
  # default grid: nothing dropped -> identical to plain seq (preserves the exact
  # float values the rest of the package matches with `==`)
  expect_identical(dbpmr:::grid_seq(-28, 14, 0.2), seq(-28, 14, 0.2))

  # irrational step from 6-decimal values (as dbpmr writes/reads them): plain
  # seq() drops the final bin; grid_seq restores it to the correct count.
  from <- round(-12 * LN10, 6); to <- round(6 * LN10, 6); by <- round(0.1 * LN10, 6)
  expect_lt(length(seq(from, to, by)), 181L)                  # seq dropped one
  expect_equal(length(dbpmr:::grid_seq(from, to, by)), round((to - from) / by) + 1L)
  expect_equal(length(dbpmr:::grid_seq(from, to, by)), 181L)
})

test_that("Read.In mass axis aligns with the data columns on an irrational grid", {
  old <- setwd(tempdir())
  on.exit({ unlink("grid_run", recursive = TRUE); setwd(old) })
  LN10 <- log(10)

  run  <- Setup.Run("grid_run", 1, 0, 0, TRUE, 1)
  grid <- Setup.Grid(run, mmin = -12 * LN10, mmax = 6 * LN10,
                     mstep = 0.1 * LN10, moutstep = 0.1 * LN10,
                     tmax = 1, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton", mmin = -12 * LN10, mmax = -3 * LN10)
  pe <- Setup.Pelagic(run, filename = "fish", mmin = -3 * LN10, mmax = 6 * LN10, rep_method = 2)
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe)))

  fish <- Read.In("grid_run", "fish")
  expect_equal(length(fish@mrange), ncol(fish@uvals) - 3)   # axis aligns with data
  expect_equal(length(fish@mrange), 181L)
  expect_equal(round(tail(fish@mrange, 1) / LN10), 6)        # top bin is 10^6 g
})

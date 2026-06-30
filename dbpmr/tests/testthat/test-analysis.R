# Regression tests for the analysis-side bug fixes (issues #3 and #4).
# The Average.Time test runs the compiled engine, so it needs the installed
# package; it writes its output under a temporary directory.

test_that("Setup.Pelagic stores gamma_prey, not gamma_pred", {
  old <- setwd(tempdir())
  on.exit({ unlink("g_run", recursive = TRUE); setwd(old) })

  run <- Setup.Run("g_run", no_pelagic = 1, no_benthic = 0,
                   spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  pe  <- Setup.Pelagic(run, filename = "fish", gamma_prey = 0.123, gamma_pred = 0.456)

  expect_equal(pe@gamma_prey, 0.123)   # was incorrectly set to gamma_pred
  expect_equal(pe@gamma_pred, 0.456)
})

test_that("Read.In does not swap pelagic feeding preferences", {
  old <- setwd(tempdir())
  on.exit({ unlink("p_run", recursive = TRUE); setwd(old) })

  run  <- Setup.Run("p_run", 1, 0, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = 1)
  pl   <- Setup.Plankton(run, filename = "plankton")
  pe   <- Setup.Pelagic(run, filename = "fish", pref_pla = 0.11, pref_pel = 0.77)
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe)))

  fish <- Read.In("p_run", "fish")
  expect_equal(fish@species@pref_pla, 0.11)   # was read from the pref_pel column
  expect_equal(fish@species@pref_pel, 0.77)
})

test_that("Average.Time averages the whole window and matches Extract.Time shape", {
  old <- setwd(tempdir())
  on.exit({ unlink("a_run", recursive = TRUE); setwd(old) })

  run  <- Setup.Run("a_run", 1, 0, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = 2)
  pl   <- Setup.Plankton(run, filename = "plankton")
  pe   <- Setup.Pelagic(run, filename = "fish")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe)))
  fish <- Read.In("a_run", "fish")

  full <- Average.Time(fish)                      # default = whole time range
  expect_s4_class(full, "timestep.data")
  expect_gt(length(full@trange), 1)               # the old bug collapsed this to 1 step
  expect_identical(full@spatial_dim, 0L)          # was @spatial.dim via a data.frame index

  ex <- Extract.Time(fish, max(fish@trange))
  expect_identical(names(full@data), names(ex@data))
  expect_equal(dim(full@data), dim(ex@data))

  win <- Average.Time(fish, c(min(fish@trange), fish@trange[3]))
  expect_equal(length(win@trange), 3L)
})

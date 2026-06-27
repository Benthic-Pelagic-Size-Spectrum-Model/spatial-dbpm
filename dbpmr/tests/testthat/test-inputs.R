# Regression tests for the input-writers: Setup.ts, Setup.Rep and Setup.fishing
# must write to <run>/Input/<species>_*.txt (capital "Input", no run-name
# prefix), which is exactly where the C engine reads them.

test_that("Setup.ts/Setup.Rep/Setup.fishing write to <run>/Input/<species>_*.txt", {
  old <- setwd(tempdir())
  on.exit({ unlink("inp_run", recursive = TRUE); setwd(old) })

  run  <- Setup.Run("inp_run", no_pelagic = 1, no_benthic = 0,
                    spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1)
  pe   <- Setup.Pelagic(run, filename = "fish",
                        rep_method = 1, initial_flag = TRUE, fishing_flag = TRUE)

  Setup.ts(pe,      run, grid, func = function(m, t, x, y) 1e-3)
  Setup.Rep(pe,     run, grid, func = function(t, x, y) 1e-6)
  Setup.fishing(pe, run, grid, func = function(m, t, x, y) 0.2)

  # The engine looks for exactly these names:
  expect_true(file.exists(file.path("inp_run", "Input", "fish_ts.txt")))
  expect_true(file.exists(file.path("inp_run", "Input", "fish_rep_ts.txt")))
  expect_true(file.exists(file.path("inp_run", "Input", "fish_fishing_ts.txt")))

  # The old buggy run-name-prefixed names must NOT be produced.
  expect_false(file.exists(file.path("inp_run", "Input", "inp_run_fish_rep_ts.txt")))
  expect_false(file.exists(file.path("inp_run", "Input", "inp_run_fish_fishing_ts.txt")))
})

test_that("Setup.ts (ts_flag = FALSE) populates all spatial cells", {
  old <- setwd(tempdir())
  on.exit({ unlink("sp_run", recursive = TRUE); setwd(old) })

  run  <- Setup.Run("sp_run", no_pelagic = 1, no_benthic = 0,
                    spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, xmax = 50, xstep = 50, ymax = 50, ystep = 50, tmax = 1)
  pe   <- Setup.Pelagic(run, filename = "fish", initial_flag = TRUE)

  nbins <- length(seq(pe@mmin, pe@mmax, grid@mstep))
  ncell <- length(seq(grid@xmin, grid@xmax, grid@xstep)) *
           length(seq(grid@ymin, grid@ymax, grid@ystep))   # 2 x 2 = 4

  # one distinct row per cell -> every cell should be populated (the old bug
  # filled only the first cell, leaving the rest zero)
  mat <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = ncell, ncol = nbins)
  Setup.ts(pe, run, grid, mat = mat)
  f <- as.matrix(read.csv(file.path("sp_run", "Input", "fish_ts.txt"), header = FALSE))

  expect_equal(nrow(f), ncell)
  expect_true(all(rowSums(f) > 0))            # no all-zero cells
  expect_true(all(diff(rowSums(f)) > 0))      # each cell got its own row

  # a single row is broadcast uniformly to every cell
  Setup.ts(pe, run, grid, mat = matrix(0.5, nrow = 1, ncol = nbins))
  f1 <- as.matrix(read.csv(file.path("sp_run", "Input", "fish_ts.txt"), header = FALSE))
  expect_equal(length(unique(round(rowSums(f1), 6))), 1L)

  # wrong row count is rejected
  expect_error(Setup.ts(pe, run, grid, mat = matrix(0.5, nrow = 3, ncol = nbins)))
})

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

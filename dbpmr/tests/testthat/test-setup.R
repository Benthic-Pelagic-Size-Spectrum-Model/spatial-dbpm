# These tests exercise the pure-R configuration helpers. They create their
# output directories inside a temporary directory and clean up afterwards, so
# they do not require the compiled C engine.

test_that("Setup.Run builds a valid run.params object", {
  old <- setwd(tempdir()); on.exit({ unlink("run_ok", recursive = TRUE); setwd(old) })

  run <- Setup.Run("run_ok", no_pelagic = 1, no_benthic = 0,
                   spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)

  expect_s4_class(run, "run.params")
  expect_identical(run@no_pelagic, 1L)
  expect_identical(run@no_benthic, 0L)
  expect_true(run@coupled_flag)
  expect_true(dir.exists(file.path(tempdir(), "run_ok")))
})

test_that("Setup.Run rejects invalid arguments", {
  old <- setwd(tempdir()); on.exit(setwd(old))

  expect_error(
    Setup.Run("run_bad", no_pelagic = 0, no_benthic = 0,
              spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  )
  expect_error(
    Setup.Run("run_bad", no_pelagic = 1, no_benthic = 0,
              spatial_dim = 7, coupled_flag = TRUE, diff_method = 1)
  )
})

test_that("Setup.Grid collapses unused spatial dimensions", {
  old <- setwd(tempdir()); on.exit({ unlink("run_grid", recursive = TRUE); setwd(old) })

  run  <- Setup.Run("run_grid", no_pelagic = 1, no_benthic = 0,
                    spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run)

  expect_s4_class(grid, "grid.params")
  # spatial_dim == 0 collapses both x and y to a single trivial cell
  expect_equal(grid@xmax, 0)
  expect_equal(grid@ymax, 0)
})

test_that("Setup.Plankton builds a plankton.params object", {
  old <- setwd(tempdir()); on.exit({ unlink("run_pla", recursive = TRUE); setwd(old) })

  run <- Setup.Run("run_pla", no_pelagic = 1, no_benthic = 0,
                   spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  pl  <- Setup.Plankton(run, filename = "plankton")

  expect_s4_class(pl, "plankton.params")
  expect_identical(pl@speciestype, "plankton")
  expect_error(Setup.Plankton(run))  # filename is required
})

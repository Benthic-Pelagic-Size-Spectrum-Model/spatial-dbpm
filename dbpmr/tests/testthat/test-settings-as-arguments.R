# Two settings that were only reachable through environment variables are ordinary arguments:
# the recruitment floor (imm_frac in Setup.Pelagic(), was PEL_IMM_FRAC) and whether temperature scales
# background mortality (temp_on_mort in Setup.Forcing(), was TEMP_ON_MORT). The environment variables still
# work when the argument is left at its default.

# A run of nx cells, heavily fished so that the recruitment floor binds, with a warm temperature multiplier.
run_settings <- function(nx = 1, imm_frac = NA, temp_on_mort = NULL, env = character()) {
  if (length(env)) { do.call(Sys.setenv, as.list(env)); on.exit(Sys.unsetenv(names(env)), add = TRUE) }
  wd <- tempfile("set"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = if (nx > 1) 1 else 0, coupled_flag = TRUE, diff_method = 1)
  grid <- if (nx > 1) Setup.Grid(run, tmax = 6, tstep = 1/48, toutstep = 1, xmin = 0, xmax = nx - 1, xstep = 1, xoutstep = 1)
          else Setup.Grid(run, tmax = 6, tstep = 1/48, toutstep = 1)
  spec <- function(sp) sp@u_0 * exp(sp@lambda * dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep))
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE); pl <- Setup.ts(pl, run, grid, initial = spec(pl))
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE, fishing_flag = TRUE, imm_frac = imm_frac)
  pe <- Setup.ts(pe, run, grid, initial = spec(pe)); pe <- Setup.fishing(pe, run, grid, rate = 6, min_size = log(10))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); be <- Setup.ts(be, run, grid, initial = spec(be))
  de <- Setup.Detritus(run, filename = "detritus")
  fo <- Setup.Forcing(run, grid, pel_tempeff = 1.5, ben_tempeff = 1.5, temp_on_mort = temp_on_mort)
  out <- NULL
  invisible(capture.output(out <- SizeSpectrum(run, grid, pl, pe, be, de, forcing = fo, output = "arrays", variables = "mortality")))
  out
}
differ <- function(a, b) !isTRUE(all.equal(a$abundance$fish, b$abundance$fish, tolerance = 1e-6))

test_that("the recruitment floor given as imm_frac equals the PEL_IMM_FRAC route, and overrides it", {
  off <- run_settings()
  arg <- run_settings(imm_frac = 0.9)
  env <- run_settings(env = c(PEL_IMM_FRAC = "0.9"))
  expect_true(differ(arg, off))                                                # the floor binds in this run
  expect_identical(arg$abundance, env$abundance)
  expect_identical(run_settings(imm_frac = 0, env = c(PEL_IMM_FRAC = "0.9"))$abundance, off$abundance)   # the argument wins
  expect_identical(run_settings(nx = 2, imm_frac = 0.9)$abundance, run_settings(nx = 2, env = c(PEL_IMM_FRAC = "0.9"))$abundance)
  # recorded with the run, whichever route was used
  expect_equal(arg$parameters$groups$fish$imm_frac, 0.9); expect_length(arg$parameters$environment, 0)
  expect_true(is.na(env$parameters$groups$fish$imm_frac)); expect_equal(env$parameters$environment, c(PEL_IMM_FRAC = "0.9"))
})

test_that("temp_on_mort given with the forcing equals the TEMP_ON_MORT route, and overrides it", {
  for (nx in c(1, 2)) {
    on  <- run_settings(nx)                                                    # default: temperature scales mortality
    arg <- run_settings(nx, temp_on_mort = FALSE)
    env <- run_settings(nx, env = c(TEMP_ON_MORT = "0"))
    expect_true(differ(arg, on), info = nx)
    expect_identical(arg$abundance, env$abundance, info = nx)
    expect_identical(run_settings(nx, temp_on_mort = TRUE, env = c(TEMP_ON_MORT = "0"))$abundance, on$abundance, info = nx)   # the argument wins
    expect_identical(run_settings(nx)$abundance, on$abundance, info = nx)      # and a later run in the same session is not affected
    # the smallest fish die faster when temperature also scales background mortality (compared after the first
    # year: the forcing acts from the first time step, not on the starting state)
    small <- which(on$abundance$fish[1, 1, , 1] > 0)[1]
    expect_gt(on$mortality$fish[1, 1, small, 2], arg$mortality$fish[1, 1, small, 2])
  }
  expect_false(arg$parameters$forcing$temp_on_mort); expect_null(on$parameters$forcing$temp_on_mort)
})

test_that("the new arguments are validated and work from a species table", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1)
  expect_error(Setup.Pelagic(run, filename = "fish", imm_frac = -0.1), "imm_frac")
  expect_error(Setup.Pelagic(run, filename = "fish", imm_frac = c(0.1, 0.2)), "imm_frac")
  expect_true(is.na(Setup.Pelagic(run, filename = "fish")@imm_frac))
  expect_error(Setup.Forcing(run, grid, pel_tempeff = 1, temp_on_mort = "no"), "TRUE, FALSE or NULL")
  expect_output(print(Setup.Forcing(run, grid, pel_tempeff = 1, temp_on_mort = FALSE)), "does not scale background mortality")
  sp <- Setup.Species(run, data.frame(species = c("fish", "benthos"), type = c("pelagic", "benthic"), imm_frac = c(0.15, NA)))
  expect_equal(sp$pelagic[[1]]@imm_frac, 0.15)
})

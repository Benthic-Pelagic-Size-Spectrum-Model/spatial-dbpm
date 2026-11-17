# Fishing effort allocated among cells by gravity during a spatial run (issue #58, step 5c):
# Setup.Gravity() + SizeSpectrum(gravity = ). The shares the C code applies are checked against the
# same calculation done here in R from the output arrays.

LO    <- log(10) + 0.05                       # fished above this size (between two grid points)
NX    <- 4
SCALE <- c(1, 2, 3, 0.5)                      # starting abundance of the fish differs between cells
ACC   <- c(1, 0.5, 2, 1)                      # accessibility of each cell
TPEL  <- c(100, 200, 50, 100); TBEN <- 20     # layer thickness (m) of each group
QP    <- 0.3; QB <- 0.15                      # fishing mortality at evenly spread effort

# One spatial run of NX cells in a row. `pel`, `ben`: rate arguments for Setup.fishing() (anything it accepts,
# with `time` the layer start times); `gravity`: arguments for Setup.Gravity() or NULL; `file`: write the
# fishing files instead of passing arrays (constant rates only).
run_gravity <- function(pel = QP, ben = QB, time = NULL, gravity = list(), scale = SCALE, tmax = 2, toutstep = 1/48,
                        prey = 0, comp = 0, file = FALSE) {
  wd <- tempfile("grav"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 1, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = toutstep, xmin = 0, xmax = NX - 1, xstep = 1, xoutstep = 1)
  start <- function(sp, f = 1) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); outer(f, sp@u_0 * exp(sp@lambda * m)) }
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE); Setup.ts(pl, run, grid, mat = start(pl))
  pe <- Setup.Pelagic(run, filename = "fish", comp = comp, prey = prey, pred = 0, initial_flag = TRUE, fishing_flag = TRUE)
  Setup.ts(pe, run, grid, mat = start(pe, scale))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE, fishing_flag = TRUE); Setup.ts(be, run, grid, mat = start(be))
  de <- Setup.Detritus(run, filename = "detritus")
  if (file) {
    Setup.fishing(pe, run, grid, func = function(m, t, x, y) if (m >= LO) pel else 0)
    Setup.fishing(be, run, grid, func = function(m, t, x, y) if (m >= LO) ben else 0)
  } else {
    pe <- Setup.fishing(pe, run, grid, rate = pel, min_size = LO, time = time)
    be <- Setup.fishing(be, run, grid, rate = ben, min_size = LO, time = time)
  }
  gv <- if (is.null(gravity)) NULL else do.call(Setup.Gravity, c(list(run, grid), gravity))
  out <- NULL
  invisible(capture.output(out <- SizeSpectrum(run, grid, pl, pe, be, de, output = "arrays", gravity = gv)))
  out
}
GRAV <- list(access = ACC, pelagic_thickness = TPEL, benthic_thickness = TBEN, update = 0.5)

# fishing mortality of the fished sizes of a species, [x, t] (the same for every fished size of the species)
f_of <- function(out, sp) out$fishing[[sp]][, 1, which(out$m >= LO)[1], ]
f_pel <- function(out) f_of(out, "fish"); f_ben <- function(out) f_of(out, "benthos")
# multiplier ncell * share from the abundances at output index j, for rates qp, qb at evenly spread effort (per cell)
multiplier <- function(out, j, qp = QP, qb = QB, acc = ACC) {
  w <- exp(out$m) * diff(out$m)[1] * (out$m >= LO)          # abundance is zero outside a species' size range
  catch <- function(u) as.vector(u[, 1, , j] %*% w)
  a <- acc * (TPEL * qp * catch(out$abundance$fish) + TBEN * qb * catch(out$abundance$benthos))
  NX * a / sum(a)
}

test_that("fishing mortality follows the gravity shares, recalculated at each update", {
  out <- run_gravity(gravity = GRAV)
  expect_equal(max(out$fishing$fish[, , out$m < LO, ]), 0)                       # still fished only above the size limit
  j <- seq_along(out$t)
  mult <- sapply(j - (j - 1) %% 24, function(j0) multiplier(out, j0))            # from the latest update (every half year)
  expect_equal(f_pel(out), QP * mult, tolerance = 1e-10)
  expect_equal(f_ben(out), QB * mult, tolerance = 1e-10)
  expect_equal(colMeans(f_pel(out)), rep(QP, length(out$t)), tolerance = 1e-12)  # effort is moved, not added
  expect_gt(diff(range(f_pel(out)[, 1])), 0.1)                                   # and the cells do differ
  expect_false(isTRUE(all.equal(f_pel(out)[, 1], f_pel(out)[, 49])))             # shares change as the fish change
})

test_that("a forced run with the mortalities of a gravity run gives the same run", {
  g <- run_gravity(gravity = GRAV)
  upd <- seq(1, length(g$t), by = 24)                                            # the half-yearly updates
  forced <- run_gravity(pel = array(f_pel(g)[, upd], c(NX, 1, length(upd))), ben = array(f_ben(g)[, upd], c(NX, 1, length(upd))),
                        time = g$t[upd], gravity = NULL)
  for (v in c("abundance", "fishing", "catch")) expect_equal(forced[[v]], g[[v]], tolerance = 1e-10, info = v)
})

test_that("gravity changes nothing when every cell is the same, and works with fishing files", {
  same <- list(access = 1, update = 0)
  expect_equal(run_gravity(gravity = same, scale = rep(1, NX))$abundance, run_gravity(gravity = NULL, scale = rep(1, NX))$abundance,
               tolerance = 1e-12)
  expect_equal(run_gravity(gravity = GRAV, file = TRUE)$abundance, run_gravity(gravity = GRAV)$abundance, tolerance = 1e-10)
})

test_that("a closed or inaccessible cell gets no effort and the effort moves to the other cells", {
  closed <- array(c(1, 0, 1, 1), c(NX, 1, 1))
  out <- run_gravity(pel = QP * closed, ben = QB * closed, gravity = GRAV)
  expect_equal(f_pel(out)[2, ], rep(0, length(out$t)))
  expect_equal(f_pel(out)[, 1], QP * multiplier(out, 1, QP * closed[, 1, 1], QB * closed[, 1, 1]), tolerance = 1e-10)
  expect_equal(colSums(f_pel(out)), rep(QP * NX, length(out$t)), tolerance = 1e-12)

  noacc <- run_gravity(gravity = modifyList(GRAV, list(access = c(1, 0.5, 0, 1))))
  expect_equal(f_pel(noacc)[3, ], rep(0, length(noacc$t)))
  expect_equal(colSums(f_pel(noacc)), rep(QP * NX, length(noacc$t)), tolerance = 1e-12)
})

test_that("effort follows abundance at each update and a change in the rates straight away", {
  # unfished first year; abundances are updated at 0 and 2 years
  out <- run_gravity(pel = c(0, QP), ben = c(0, QB), time = c(0, 1), gravity = modifyList(GRAV, list(update = 2)), tmax = 3)
  at <- function(t) which(abs(out$t - t) < 1e-9)
  expect_equal(max(f_pel(out)[, out$t < 1 - 1e-9]), 0)
  expect_equal(f_pel(out)[, at(1)], QP * multiplier(out, at(0)), tolerance = 1e-10)       # fishing starts: abundances of year 0
  expect_equal(f_pel(out)[, at(1.75)], f_pel(out)[, at(1)], tolerance = 1e-12)             # held
  expect_equal(f_pel(out)[, at(2)], QP * multiplier(out, at(2)), tolerance = 1e-10)       # updated
})

# A marine protected area as a mask (1 = protected) written to and read back from a NetCDF file,
# then applied to the fishing rate: [x, y] if it does not change, [x, y, time] if it does.
mpa_from_netcdf <- function(mask, years = NULL) {
  f <- tempfile(fileext = ".nc"); on.exit(unlink(f))
  d <- list(ncdf4::ncdim_def("x", "", seq_len(dim(mask)[1])), ncdf4::ncdim_def("y", "", seq_len(dim(mask)[2])))
  if (!is.null(years)) d[[3]] <- ncdf4::ncdim_def("time", "years", years)
  v <- ncdf4::ncvar_def("mpa", "", d, prec = "integer")
  nc <- ncdf4::nc_create(f, v); ncdf4::ncvar_put(nc, v, mask); ncdf4::nc_close(nc)
  nc <- ncdf4::nc_open(f); on.exit(ncdf4::nc_close(nc), add = TRUE, after = FALSE)
  list(mask = ncdf4::ncvar_get(nc, "mpa", collapse_degen = FALSE), years = if (!is.null(years)) as.vector(nc$dim$time$vals))
}

test_that("a marine protected area read from NetCDF closes cells, in space and through time", {
  skip_if_not_installed("ncdf4")
  # fixed in time: cell 2 is protected throughout
  fixed <- mpa_from_netcdf(array(c(0, 1, 0, 0), c(NX, 1)))
  open1 <- 1 - fixed$mask                                                        # [x, y], 1 = open to fishing
  out <- run_gravity(pel = QP * open1, ben = QB * open1, gravity = GRAV)
  expect_equal(f_pel(out)[2, ], rep(0, length(out$t)))
  expect_equal(colSums(f_pel(out)), rep(QP * NX, length(out$t)), tolerance = 1e-12)       # its effort goes to the open cells

  # changing in time: cell 2 is protected from year 1, cell 3 as well from year 1.5; abundances updated at 0 and 2 years
  vary <- mpa_from_netcdf(array(c(0, 0, 0, 0,  0, 1, 0, 0,  0, 1, 1, 0), c(NX, 1, 3)), years = c(0, 1, 1.5))
  open3 <- 1 - vary$mask                                                         # [x, y, time]
  out <- run_gravity(pel = QP * open3, ben = QB * open3, time = vary$years, gravity = modifyList(GRAV, list(update = 2)), tmax = 3)
  at <- function(t) which(abs(out$t - t) < 1e-9)
  expect_true(all(f_pel(out)[2, out$t < 1 - 1e-9] > 0)); expect_true(all(f_pel(out)[2, out$t > 1 - 1e-9] == 0))
  expect_true(all(f_pel(out)[3, out$t < 1.5 - 1e-9] > 0)); expect_true(all(f_pel(out)[3, out$t > 1.5 - 1e-9] == 0))
  expect_equal(colSums(f_pel(out)), rep(QP * NX, length(out$t)), tolerance = 1e-12)       # effort is displaced at once, never lost
  for (tt in c(0.5, 1, 1.25, 1.5, 1.75)) {                                       # shares: current closures, abundances of year 0
    k <- findInterval(tt + 1e-9, vary$years)
    expect_equal(f_pel(out)[, at(tt)], QP * open3[, 1, k] * multiplier(out, at(0), QP * open3[, 1, k], QB * open3[, 1, k]),
                 tolerance = 1e-10, info = tt)
  }
  expect_equal(f_pel(out)[, at(2.5)], QP * open3[, 1, 3] * multiplier(out, at(2), QP * open3[, 1, 3], QB * open3[, 1, 3]), tolerance = 1e-10)
  # the protected cells end with more large fish than the same cells left open
  unprotected <- run_gravity(gravity = modifyList(GRAV, list(update = 2)), tmax = 3)
  big <- out$m >= LO; last <- length(out$t)
  expect_true(all(rowSums(out$abundance$fish[2:3, 1, big, last]) > rowSums(unprotected$abundance$fish[2:3, 1, big, last])))
})

test_that("fishing and gravity run with movement between cells switched on", {
  still <- run_gravity(gravity = GRAV, toutstep = 0.5)
  moved <- run_gravity(gravity = GRAV, toutstep = 0.5, prey = 0.5, comp = 0.5)
  for (v in c(moved$abundance, moved$catch, moved$fishing)) expect_true(all(is.finite(v)) && all(v >= 0))
  expect_equal(colMeans(f_pel(moved)), rep(QP, length(moved$t)), tolerance = 1e-12)
  expect_gt(sum(moved$catch$fish[, , , length(moved$t)]), 0)
  last <- length(moved$t); big <- moved$m >= LO
  spread <- function(o) { b <- apply(o$abundance$fish[, 1, big, last], 1, sum); max(b) / min(b) }
  expect_lt(spread(moved), spread(still))                                        # diffusion evens the cells out
})

test_that("Setup.Gravity and the gravity argument are validated", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1, xmin = 0, xmax = 2, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  gv <- Setup.Gravity(run, grid, access = matrix(1:6, 3, 2), pelagic_thickness = 200, benthic_thickness = 20)
  expect_s3_class(gv, "dbpmr_gravity")
  expect_equal(dim(gv$access), c(3L, 2L)); expect_equal(gv$nupdate, 48L)
  expect_equal(Setup.Gravity(run, grid, update = 0)$nupdate, 1L)
  expect_output(print(gv), "3 x 2 grid")
  expect_error(Setup.Gravity(run, grid, access = matrix(1, 2, 3)), "x = 3, y = 2")
  expect_error(Setup.Gravity(run, grid, access = -1), "non-negative")
  expect_error(Setup.Gravity(run, grid, access = 0), "at least one cell")
  expect_error(Setup.Gravity(run, grid, pelagic_thickness = 0), "positive")
  expect_error(Setup.Gravity(run, grid, access = NA_real_), "NA, NaN or Inf")
  expect_error(Setup.Gravity(run, grid, update = 0.03), "whole number of time steps")
  run0 <- Setup.Run("Z", 1, 1, spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  expect_error(Setup.Gravity(run0, Setup.Grid(run0, tmax = 1, tstep = 1/48, toutstep = 1)), "spatial run")

  pl <- Setup.Plankton(run, filename = "plankton"); pe <- Setup.Pelagic(run, filename = "fish")
  be <- Setup.Benthic(run, filename = "benthos"); de <- Setup.Detritus(run, filename = "detritus")
  expect_error(SizeSpectrum(run, grid, pl, pe, be, de, gravity = list(access = 1)), "Setup.Gravity")
  expect_error(SizeSpectrum(run, grid, pl, pe, be, de, gravity = gv), "no species is fished")
  grid2 <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1, xmin = 0, xmax = 3, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  pf <- Setup.fishing(Setup.Pelagic(run, filename = "fish", fishing_flag = TRUE), run, grid2, rate = 0.2)
  expect_error(SizeSpectrum(run, grid2, pl, pf, be, de, gravity = gv), "different grid")
})

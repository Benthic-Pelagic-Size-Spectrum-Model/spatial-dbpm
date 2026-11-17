# Plankton time series passed to the C code as intercept/slope arrays (issue #58) instead of a
# <run>/Input/plankton_ts.txt file. Both routes must give the same run.

mon <- function(t) (floor(t * 12 + 1e-9) %% 12) + 1          # month (1..12) at model time t

# final pelagic biomass per cell; `how` selects the input route for the same monthly plankton series
# a, b: [x, month] matrices of intercept and slope (u(m) = a * exp(b * m))
run_plankton <- function(a, b, how = c("file", "layer", "time"), assign = TRUE) {
  how <- match.arg(how)
  wd <- tempfile("arr"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  nx <- nrow(a)
  run  <- Setup.Run("R", 1, 1, spatial_dim = if (nx > 1) 1 else 0, coupled_flag = TRUE, diff_method = 1)
  grid <- if (nx > 1) Setup.Grid(run, tmax = 3, tstep = 1/48, toutstep = 1, xmin = 0, xmax = nx - 1, xstep = 1, xoutstep = 1)
          else Setup.Grid(run, tmax = 3, tstep = 1/48, toutstep = 1)
  xs <- dbpmr:::grid_seq(grid@xmin, grid@xmax, grid@xstep)
  pl <- Setup.Plankton(run, filename = "plankton", ts_flag = TRUE)
  if (how == "file") {
    Setup.ts(pl, run, grid, func = function(m, t, x, y) { k <- which.min(abs(xs - x)); a[k, mon(t)] * exp(b[k, mon(t)] * m) })
  } else {
    A <- if (nx > 1) a else a[1, ]; B <- if (nx > 1) b else b[1, ]
    r <- if (how == "layer") Setup.ts(pl, run, grid, intercept = A, slope = B, layer = mon(seq(0, grid@tmax, grid@tstep)))
         else {                                              # 37 monthly layers (3 years + the final January), addressed by start time
           rep3 <- function(v) if (is.matrix(v)) v[, c(rep(1:12, 3), 1), drop = FALSE] else c(rep(v, 3), v[1])
           Setup.ts(pl, run, grid, intercept = rep3(A), slope = rep3(B), time = (0:36) / 12)
         }
    if (assign) pl <- r
  }
  uni <- function(sp) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1) }
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE); Setup.ts(pe, run, grid, mat = uni(pe))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); Setup.ts(be, run, grid, mat = uni(be))
  de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  wrote_file <- file.exists(file.path("R", "Input", "plankton_ts.txt"))
  L <- readLines("R/fish/results.txt")
  num <- function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]]))
  hdr <- num(L[5]); mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
  rws <- lapply(L[6:length(L)][nchar(L[6:length(L)]) > 0], num)
  rws <- rws[sapply(rws, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]))]
  last <- rws[sapply(rws, `[`, 1) == max(sapply(rws, `[`, 1))]
  b0 <- sapply(last, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE)); xx <- sapply(last, `[`, 2)
  list(biomass = unname(sapply(sort(unique(xx)), function(q) b0[which(xx == q)[1]])), wrote_file = wrote_file)
}

seasonal <- function(nx) list(a = outer(seq(0.008, 0.016, length.out = nx), 1 + 0.5 * sin(2 * pi * (1:12) / 12)),
                              b = outer(seq(-1.05, -0.95, length.out = nx), rep(1, 12)) + rep(0.02 * cos(2 * pi * (1:12) / 12), each = nx))

test_that("plankton arrays give the same 0-D run as the time-series file, without writing one", {
  s <- seasonal(1)
  f <- run_plankton(s$a, s$b, "file"); r <- run_plankton(s$a, s$b, "layer")
  expect_true(f$wrote_file); expect_false(r$wrote_file)
  expect_equal(r$biomass, f$biomass, tolerance = 1e-8)
})

test_that("per-cell plankton arrays give the same spatial run as the time-series file", {
  s <- seasonal(3)
  f <- run_plankton(s$a, s$b, "file"); r <- run_plankton(s$a, s$b, "layer")
  expect_equal(r$biomass, f$biomass, tolerance = 1e-8)
  expect_true(all(diff(r$biomass) != 0))                    # the cells really differ
  expect_equal(run_plankton(s$a, s$b, "time")$biomass, r$biomass, tolerance = 1e-12)   # `time` and `layer` agree
})

test_that("forgetting to assign the result of Setup.ts(intercept=, slope=) is an error, not a silent default", {
  s <- seasonal(1)
  expect_error(run_plankton(s$a, s$b, "layer", assign = FALSE), "assign the result")
})

test_that("Setup.ts validates intercept/slope arrays", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1, xmin = 0, xmax = 2, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton", ts_flag = TRUE)
  ok <- array(0.01, c(3, 2, 12)); sl <- array(-1, c(3, 2, 12)); lay <- mon(seq(0, 1, 1/48))
  out <- Setup.ts(pl, run, grid, intercept = ok, slope = sl, layer = lay)
  expect_s4_class(out, "plankton.params"); expect_equal(dim(out@ts_arrays$intercept), c(3L, 2L, 12L))
  expect_error(Setup.ts(pl, run, grid, intercept = array(0.01, c(2, 3, 12)), slope = sl, layer = lay), "x = 3, y = 2")
  expect_error(Setup.ts(pl, run, grid, intercept = ok, slope = sl), "12 layers but the run has 49 time points")
  expect_error(Setup.ts(pl, run, grid, intercept = ok, slope = sl, layer = c(lay[-1], 13)), "between 1 and the number of layers")
  expect_error(Setup.ts(pl, run, grid, intercept = -ok, slope = sl, layer = lay), "non-negative")
  expect_error(Setup.ts(pl, run, grid, intercept = replace(ok, 1, NaN), slope = sl, layer = lay), "NA, NaN or Inf")
  expect_error(Setup.ts(pl, run, grid, intercept = ok, slope = sl, layer = lay, time = 0:11 / 12), "only one of layer or time")
  expect_error(Setup.ts(Setup.Plankton(run, filename = "plankton"), run, grid, intercept = ok, slope = sl, layer = lay), "ts_flag = TRUE")
  expect_error(Setup.ts(Setup.Pelagic(run, filename = "fish", initial_flag = TRUE), run, grid, intercept = ok, slope = sl, layer = lay), "plankton only")
})

# ---------------------------------------------------------------------------------------------
# Forcing channels (temperature, export, depth) as arrays via Setup.Forcing(), instead of
# <run>/Input/forcing_ts.txt.

# final pelagic and benthic biomass per cell. tp, tb: [x, month] temperature multipliers; dep: depth per cell.
run_forcing <- function(tp, tb, dep, how = c("file", "array"), sink = 0.3) {
  how <- match.arg(how)
  wd <- tempfile("frc"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  nx <- nrow(tp); tmax <- 3; dt <- 1/48
  run  <- Setup.Run("R", 1, 1, spatial_dim = if (nx > 1) 1 else 0, coupled_flag = TRUE, diff_method = 1)
  grid <- if (nx > 1) Setup.Grid(run, tmax = tmax, tstep = dt, toutstep = 1, xmin = 0, xmax = nx - 1, xstep = 1, xoutstep = 1)
          else Setup.Grid(run, tmax = tmax, tstep = dt, toutstep = 1)
  uni <- function(sp) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1) }
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE); Setup.ts(pl, run, grid, mat = uni(pl))
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE); Setup.ts(pe, run, grid, mat = uni(pe))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); Setup.ts(be, run, grid, mat = uni(be))
  de <- Setup.Detritus(run, filename = "detritus")
  fo <- NULL
  if (how == "file") {
    nst <- round(tmax / dt) + 2; jm <- mon(seq_len(nst) * dt)       # row r is applied while advancing to time r * dt
    rows <- unlist(lapply(jm, function(j) sprintf("%.17g,%.17g,%.17g,%.17g", tp[, j], tb[, j], sink, dep)))   # x rows per timestep
    writeLines(c("pel_tempeff,ben_tempeff,sinking_rate,depth", rows), file.path("R", "Input", "forcing_ts.txt"))
    if (nx > 1) { Sys.setenv(FORCING_PERCELL = "1"); on.exit(Sys.unsetenv("FORCING_PERCELL"), add = TRUE) }
  } else {
    a3 <- function(m) array(m, c(nx, 1, 12))
    fo <- Setup.Forcing(run, grid, pel_tempeff = a3(tp), ben_tempeff = a3(tb), sinking_rate = sink,
                        depth = matrix(dep, nx, 1), layer = mon(seq(0, tmax, dt)))
  }
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de, forcing = fo)))
  final <- function(nm) {
    L <- readLines(file.path("R", nm, "results.txt"))
    num <- function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]]))
    hdr <- num(L[5]); mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
    rws <- lapply(L[6:length(L)][nchar(L[6:length(L)]) > 0], num)
    rws <- rws[sapply(rws, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]))]
    last <- rws[sapply(rws, `[`, 1) == max(sapply(rws, `[`, 1))]
    b0 <- sapply(last, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE)); xx <- sapply(last, `[`, 2)
    unname(sapply(sort(unique(xx)), function(q) b0[which(xx == q)[1]]))
  }
  list(pel = final("fish"), ben = final("benthos"), wrote_file = file.exists(file.path("R", "Input", "forcing_ts.txt")))
}
temps <- function(nx) list(tp = outer(seq(0.7, 1.3, length.out = nx), 1 + 0.2 * sin(2 * pi * (1:12) / 12)),
                           tb = outer(seq(0.9, 1.1, length.out = nx), 1 + 0.1 * cos(2 * pi * (1:12) / 12)),
                           dep = seq(150, 200, length.out = nx))

test_that("forcing arrays give the same 0-D run as forcing_ts.txt, without writing one", {
  s <- temps(1)
  f <- run_forcing(s$tp, s$tb, s$dep, "file"); r <- run_forcing(s$tp, s$tb, s$dep, "array")
  expect_true(f$wrote_file); expect_false(r$wrote_file)
  expect_equal(r$pel, f$pel, tolerance = 1e-10); expect_equal(r$ben, f$ben, tolerance = 1e-10)
})

test_that("per-cell forcing arrays give the same spatial run as per-cell forcing_ts.txt", {
  s <- temps(3)
  f <- run_forcing(s$tp, s$tb, s$dep, "file"); r <- run_forcing(s$tp, s$tb, s$dep, "array")
  expect_equal(r$pel, f$pel, tolerance = 1e-10); expect_equal(r$ben, f$ben, tolerance = 1e-10)
  expect_true(all(diff(r$pel) != 0))                         # the cells really differ
})

test_that("Setup.Forcing validates its input and SizeSpectrum refuses arrays plus a file", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1, xmin = 0, xmax = 2, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  lay <- mon(seq(0, 1, 1/48)); a <- array(1, c(3, 2, 12))
  fo <- Setup.Forcing(run, grid, pel_tempeff = a, ben_tempeff = 1:12 / 12, depth = matrix(200, 3, 2), sinking_rate = 0.3, layer = lay)
  expect_s3_class(fo, "dbpmr_forcing")
  expect_equal(lapply(fo$channels, dim), list(pel_tempeff = c(3L, 2L, 12L), ben_tempeff = c(3L, 2L, 12L), depth = c(3L, 2L, 1L), sinking_rate = c(3L, 2L, 1L)))
  expect_equal(fo$channels$ben_tempeff[2, 1, ], 1:12 / 12)
  expect_error(Setup.Forcing(run, grid, temperature = a, layer = lay), "unknown forcing channel")
  expect_error(Setup.Forcing(run, grid, pel_tempeff = a, ben_tempeff = array(1, c(3, 2, 6)), layer = lay), "same number of layers")
  expect_error(Setup.Forcing(run, grid, pel_tempeff = array(1, c(2, 3, 12)), layer = lay), "x = 3, y = 2")
  expect_error(Setup.Forcing(run, grid, pel_tempeff = a), "12 layers but the run has 49 time points")
  expect_error(Setup.Forcing(run, grid, pel_tempeff = replace(a, 5, NA), layer = lay), "NA, NaN or Inf")
  expect_error(Setup.Forcing(run, grid, residence_time = matrix(1:6, 3, 2)), "same in every cell")
  expect_error(Setup.Forcing(run, grid), "at least one forcing channel")
  pl <- Setup.Plankton(run, filename = "plankton"); pe <- Setup.Pelagic(run, filename = "fish")
  be <- Setup.Benthic(run, filename = "benthos"); de <- Setup.Detritus(run, filename = "detritus")
  expect_error(SizeSpectrum(run, grid, pl, pe, be, de, forcing = list(pel_tempeff = 1)), "Setup.Forcing")
  dir.create(file.path("R", "Input"), showWarnings = FALSE); writeLines("pel_tempeff", file.path("R", "Input", "forcing_ts.txt"))
  expect_error(SizeSpectrum(run, grid, pl, pe, be, de, forcing = fo), "use one or the other")
})


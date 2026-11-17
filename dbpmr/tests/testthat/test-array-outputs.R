# Model output returned to R as arrays (issue #58, step 1) instead of, or as well as, the text result
# files, and written to NetCDF. The arrays must hold exactly what the text files hold.

# a small 2-D run with movement and fishing; returns whatever SizeSpectrum() returns and leaves the run directory in place
run_out <- function(output, xoutstep = 1, moutstep = 0.1 * log(10), toutstep = 0.5, move = TRUE, ...) {
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 2, tstep = 1/48, toutstep = toutstep, moutstep = moutstep,
                     xmin = 0, xmax = 4, xstep = 1, xoutstep = xoutstep, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  uni <- function(sp) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1) }
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE)
  Setup.ts(pl, run, grid, func = function(m, t, x, y) 0.01 * exp(-m) * (1 + 0.2 * x + 0.1 * y))
  pe <- Setup.Pelagic(run, filename = "fish", prey = if (move) 0.3 else 0, comp = if (move) 0.2 else 0,
                      initial_flag = TRUE, fishing_flag = TRUE)
  Setup.ts(pe, run, grid, mat = uni(pe))
  Setup.fishing(pe, run, grid, func = function(m, t, x, y) if (m > log(10)) 0.2 * (1 + 0.1 * x) else 0)
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); Setup.ts(be, run, grid, mat = uni(be))
  de <- Setup.Detritus(run, filename = "detritus")
  res <- NULL
  invisible(capture.output(res <- SizeSpectrum(run, grid, pl, pe, be, de, output = output, ...)))
  res
}
num <- function(z) suppressWarnings(as.numeric(strsplit(z, ",", fixed = TRUE)[[1]]))
# a per-timestep text result file as an [x, y, m, t] array (rows are x-outer / y-inner within each output
# time; the final time is written a second time at the end of results.txt, so only the first block is used)
read_txt <- function(group, file, nx, ny) {
  L <- readLines(file.path("R", group, paste0(file, ".txt")))
  m <- num(L[5]); m <- m[!is.na(m)]
  rows <- lapply(L[-(1:5)], num); rows <- rows[vapply(rows, function(v) length(v) == length(m) + 3 && !anyNA(v[1:3]), TRUE)]
  tt <- vapply(rows, `[`, 0, 1); ut <- unique(tt)
  a <- array(NA_real_, c(nx, ny, length(m), length(ut)))
  for (j in seq_along(ut)) { blk <- rows[tt == ut[j]][seq_len(nx * ny)]
    for (r in seq_along(blk)) a[((r - 1) %/% ny) + 1, ((r - 1) %% ny) + 1, , j] <- blk[[r]][-(1:3)] }
  list(a = a, m = m, t = ut)
}
same <- function(a, b) identical(dim(a), dim(b)) && all(abs(a - b) <= 1e-15 + 1e-12 * abs(b))   # the files hold 16 decimal places

test_that("output arrays hold exactly what the text result files hold", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  out <- run_out("both")
  expect_s3_class(out, "dbpmr_output")
  expect_named(out, c("t", "x", "y", "m", "abundance", "growth", "mortality", "predation", "fishing", "catch", "detritus"))
  expect_named(out$abundance, c("plankton", "fish", "benthos")); expect_named(out$growth, c("fish", "benthos"))
  expect_equal(dim(out$abundance$plankton), c(5L, 2L, 181L, 5L))
  expect_equal(out$t, seq(0, 2, 0.5)); expect_equal(out$x, 0:4); expect_equal(out$y, 0:1)
  expect_type(attr(out, "files"), "character")
  f <- read_txt("fish", "results", 5, 2)
  expect_equal(out$m, f$m, tolerance = 1e-6); expect_equal(out$t, f$t)
  for (g in c("plankton", "fish", "benthos")) expect_true(same(out$abundance[[g]], read_txt(g, "results", 5, 2)$a), info = g)
  for (g in c("fish", "benthos")) {
    expect_true(same(out$growth[[g]],    read_txt(g, "growth", 5, 2)$a),    info = paste(g, "growth"))
    expect_true(same(out$mortality[[g]], read_txt(g, "mortality", 5, 2)$a), info = paste(g, "mortality"))
    expect_true(same(out$fishing[[g]],   read_txt(g, "fishing", 5, 2)$a),   info = paste(g, "fishing"))
  }
  expect_named(out$predation$fish, "fish")                                     # predation on fish, by predator species
  expect_true(same(out$predation$fish$fish, read_txt("fish", "predation0", 5, 2)$a))
  expect_true(same(out$predation$benthos,   read_txt("benthos", "predation", 5, 2)$a))
  # catch = fishing mortality x abundance x body mass, at the output time
  w <- array(rep(exp(out$m), each = 5 * 2), dim(out$fishing$fish))
  expect_equal(out$catch$fish, out$fishing$fish * out$abundance$fish * w, tolerance = 1e-14)
  expect_gt(max(out$catch$fish), 0); expect_equal(max(out$catch$benthos), 0)   # only the pelagic species is fished
  # the run exercises every variable and the cell indexing
  expect_gt(max(out$fishing$fish), 0); expect_gt(max(out$predation$benthos), 0); expect_gt(max(out$growth$fish), 0)
  expect_true(any(apply(out$abundance$fish[, , , 5], 3, function(z) diff(range(z))) > 0))
  expect_true(diff(range(out$fishing$fish[, 1, 181, 5])) > 0)                  # fishing was set to vary with x
  d <- do.call(rbind, lapply(readLines(file.path("R", "detritus", "summary.txt"))[-(1:5)], num))
  expect_true(all(abs(as.vector(aperm(out$detritus, c(2, 1, 3))) - d[, 4]) <= 1e-15 + 1e-12 * abs(d[, 4])))   # file rows: time, x, y (y fastest)
})

test_that("catch summed over size is the total catch rate of summary.txt", {
  # Without movement. With movement the summary totals are computed before the movement step, so its
  # "fishing" column lags the abundance that is output; the catch array is consistent with that abundance.
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  out <- run_out("both", move = FALSE)
  sm <- do.call(rbind, lapply(readLines(file.path("R", "fish", "summary.txt"))[-(1:5)], num))   # rows: time, x, y (y fastest)
  tot <- apply(out$catch$fish, c(1, 2, 4), sum) * diff(out$m)[1]                               # [x, y, t]
  expect_gt(min(sm[, 9]), 0)
  expect_equal(as.vector(aperm(tot, c(2, 1, 3))), sm[, 9], tolerance = 1e-12)
})

test_that("output = 'arrays' returns the same arrays and writes no per-timestep text output", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  both <- run_out("both"); unlink("R", recursive = TRUE)
  arr <- run_out("arrays")
  expect_null(attr(arr, "files"))
  attr(both, "files") <- NULL
  expect_identical(arr, both)
  for (f in c("results", "growth", "mortality", "predation0", "fishing"))
    expect_length(readLines(file.path("R", "fish", paste0(f, ".txt"))), 5)     # header only
})

test_that("variables limits what is returned; abundance always is", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  full <- run_out("arrays"); unlink("R", recursive = TRUE)
  ab <- run_out("arrays", variables = "abundance")
  expect_named(ab, c("t", "x", "y", "m", "abundance", "detritus"))
  expect_identical(ab$abundance, full$abundance)
  unlink("R", recursive = TRUE)
  gf <- run_out("arrays", variables = c("growth", "fishing"))
  expect_named(gf, c("t", "x", "y", "m", "abundance", "growth", "fishing", "detritus"))
  unlink("R", recursive = TRUE)
  ca <- run_out("arrays", variables = "catch")
  expect_named(ca, c("t", "x", "y", "m", "abundance", "catch", "detritus")); expect_identical(ca$catch, full$catch)
  expect_identical(gf$fishing, full$fishing)
  unlink("R", recursive = TRUE)
  expect_error(run_out("arrays", variables = "biomass"), "should be one of")
})

test_that("output = 'files' is unchanged and the default", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  res <- run_out("files")
  expect_type(res, "character"); expect_identical(res[2:5], c("plankton", "fish", "benthos", "detritus"))
  expect_gt(length(readLines(file.path("R", "fish", "results.txt"))), 5)
})

test_that("output arrays follow the output steps in space, size and time", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  out <- run_out("both", xoutstep = 2, moutstep = 0.2 * log(10), toutstep = 1)
  expect_equal(out$x, c(0, 2, 4)); expect_equal(out$t, 0:2); expect_length(out$m, 91)
  expect_true(same(out$abundance$fish, read_txt("fish", "results", 3, 2)$a))
  expect_true(same(out$mortality$fish, read_txt("fish", "mortality", 3, 2)$a))
})

test_that("Write.NetCDF round-trips the output arrays", {
  skip_if_not_installed("ncdf4")
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  out <- run_out("arrays")
  f <- Write.NetCDF(out, "run.nc", x = 70 + out$x, y = -55 + out$y, xname = "lon", yname = "lat", xunits = "degrees_east", yunits = "degrees_north")
  nc <- ncdf4::nc_open(f); on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_setequal(names(nc$var), c("plankton_abundance", "fish_abundance", "benthos_abundance",
                                   "fish_growth", "benthos_growth", "fish_mortality", "benthos_mortality",
                                   "fish_predation_by_fish", "benthos_predation", "fish_fishing", "benthos_fishing",
                                   "fish_catch", "benthos_catch", "detritus"))
  expect_equal(vapply(nc$var$fish_abundance$dim, `[[`, "", "name"), c("lon", "lat", "size", "time"))
  expect_equal(as.vector(nc$dim$lon$vals), 70 + out$x); expect_equal(as.vector(nc$dim$time$vals), out$t)
  expect_equal(as.vector(nc$dim$size$vals), out$m)
  get <- function(v) ncdf4::ncvar_get(nc, v, collapse_degen = FALSE)
  expect_equal(get("fish_abundance"), out$abundance$fish, ignore_attr = TRUE)
  expect_equal(get("plankton_abundance"), out$abundance$plankton, ignore_attr = TRUE)
  expect_equal(get("benthos_growth"), out$growth$benthos, ignore_attr = TRUE)
  expect_equal(get("fish_predation_by_fish"), out$predation$fish$fish, ignore_attr = TRUE)
  expect_equal(get("benthos_predation"), out$predation$benthos, ignore_attr = TRUE)
  expect_equal(get("fish_fishing"), out$fishing$fish, ignore_attr = TRUE)
  expect_equal(get("fish_catch"), out$catch$fish, ignore_attr = TRUE)
  expect_equal(get("detritus"), out$detritus, ignore_attr = TRUE)
  expect_error(Write.NetCDF(list(1), "x.nc"), "SizeSpectrum")
  expect_error(Write.NetCDF(out, "x.nc", x = 1:2), "one per output cell")
})

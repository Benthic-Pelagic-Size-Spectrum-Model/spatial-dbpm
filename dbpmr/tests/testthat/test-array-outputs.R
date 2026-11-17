# Model output returned to R as arrays (issue #58, step 1) instead of, or as well as, the text result
# files, and written to NetCDF. The arrays must hold exactly what the text files hold.

# a small 2-D run with movement; returns whatever SizeSpectrum() returns and leaves the run directory in place
run_out <- function(output, xoutstep = 1, moutstep = 0.1 * log(10), toutstep = 0.5) {
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 2, tstep = 1/48, toutstep = toutstep, moutstep = moutstep,
                     xmin = 0, xmax = 4, xstep = 1, xoutstep = xoutstep, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  uni <- function(sp) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1) }
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE)
  Setup.ts(pl, run, grid, func = function(m, t, x, y) 0.01 * exp(-m) * (1 + 0.2 * x + 0.1 * y))
  pe <- Setup.Pelagic(run, filename = "fish", prey = 0.3, comp = 0.2, initial_flag = TRUE); Setup.ts(pe, run, grid, mat = uni(pe))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); Setup.ts(be, run, grid, mat = uni(be))
  de <- Setup.Detritus(run, filename = "detritus")
  res <- NULL
  invisible(capture.output(res <- SizeSpectrum(run, grid, pl, pe, be, de, output = output)))
  res
}
num <- function(z) suppressWarnings(as.numeric(strsplit(z, ",", fixed = TRUE)[[1]]))
# results.txt as an [x, y, m, t] array (rows are x-outer / y-inner within each output time; the final
# time is written a second time at the end of the file, so only the first block of each time is used)
read_results <- function(group, nx, ny) {
  L <- readLines(file.path("R", group, "results.txt"))
  m <- num(L[5]); m <- m[!is.na(m)]
  rows <- lapply(L[-(1:5)], num); rows <- rows[vapply(rows, function(v) length(v) == length(m) + 3 && !anyNA(v[1:3]), TRUE)]
  tt <- vapply(rows, `[`, 0, 1); ut <- unique(tt)
  a <- array(NA_real_, c(nx, ny, length(m), length(ut)))
  for (j in seq_along(ut)) { blk <- rows[tt == ut[j]][seq_len(nx * ny)]
    for (r in seq_along(blk)) a[((r - 1) %/% ny) + 1, ((r - 1) %% ny) + 1, , j] <- blk[[r]][-(1:3)] }
  list(a = a, m = m, t = ut)
}
same <- function(a, b) all(abs(a - b) <= 1e-15 + 1e-12 * abs(b))     # the files hold 16 decimal places

test_that("output arrays hold exactly what the text result files hold", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  out <- run_out("both")
  expect_s3_class(out, "dbpmr_output")
  expect_equal(dim(out$plankton), c(5L, 2L, 181L, 5L))
  expect_equal(out$t, seq(0, 2, 0.5)); expect_equal(out$x, 0:4); expect_equal(out$y, 0:1)
  expect_type(attr(out, "files"), "character")
  for (g in c("plankton", "fish", "benthos")) {
    f <- read_results(g, 5, 2)
    a <- switch(g, plankton = out$plankton, fish = out$pelagic$fish, benthos = out$benthic$benthos)
    expect_equal(dim(a), dim(f$a)); expect_true(same(a, f$a), info = g)
    expect_equal(out$m, f$m, tolerance = 1e-6); expect_equal(out$t, f$t)
  }
  expect_true(any(apply(out$pelagic$fish[, , , 5], 3, function(z) diff(range(z))) > 0))   # cells differ: indexing is exercised
  d <- do.call(rbind, lapply(readLines(file.path("R", "detritus", "summary.txt"))[-(1:5)], num))
  expect_true(same(as.vector(aperm(out$detritus, c(2, 1, 3))), d[, 4]))                    # file rows: time, x, y (y fastest)
})

test_that("output = 'arrays' returns the same arrays and writes no per-timestep text output", {
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  both <- run_out("both"); unlink("R", recursive = TRUE)
  arr <- run_out("arrays")
  expect_null(attr(arr, "files"))
  attr(both, "files") <- NULL
  expect_identical(arr, both)
  expect_length(readLines(file.path("R", "fish", "results.txt")), 5)        # header only
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
  f <- read_results("fish", 3, 2)
  expect_equal(dim(out$pelagic$fish), dim(f$a)); expect_true(same(out$pelagic$fish, f$a))
})

test_that("Write.NetCDF round-trips the output arrays", {
  skip_if_not_installed("ncdf4")
  wd <- tempfile("out"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  out <- run_out("arrays")
  f <- Write.NetCDF(out, "run.nc", x = 70 + out$x, y = -55 + out$y, xname = "lon", yname = "lat", xunits = "degrees_east", yunits = "degrees_north")
  nc <- ncdf4::nc_open(f); on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_setequal(names(nc$var), c("plankton", "fish", "benthos", "detritus"))
  expect_equal(vapply(nc$var$fish$dim, `[[`, "", "name"), c("lon", "lat", "size", "time"))
  expect_equal(as.vector(nc$dim$lon$vals), 70 + out$x); expect_equal(as.vector(nc$dim$time$vals), out$t)
  expect_equal(as.vector(nc$dim$size$vals), out$m)
  expect_equal(ncdf4::ncvar_get(nc, "fish", collapse_degen = FALSE), out$pelagic$fish, ignore_attr = TRUE)
  expect_equal(ncdf4::ncvar_get(nc, "plankton", collapse_degen = FALSE), out$plankton, ignore_attr = TRUE)
  expect_equal(ncdf4::ncvar_get(nc, "detritus", collapse_degen = FALSE), out$detritus, ignore_attr = TRUE)
  expect_error(Write.NetCDF(list(1), "x.nc"), "SizeSpectrum")
  expect_error(Write.NetCDF(out, "x.nc", x = 1:2), "one per output cell")
})

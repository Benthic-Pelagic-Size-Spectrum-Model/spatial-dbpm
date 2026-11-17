# Per-cell model parameters for spatial runs (issue #58).
# Setup.CellParams() writes <run>/Input/<species>_cellpar.txt so pref_ben can vary across cells.
# Not supplied (default) => the scalar pref_ben applies everywhere, identical to standard dbpmr.
LN10 <- log(10)

# final pelagic + benthic biomass of a single 0-D cell with a scalar pref_ben
biomass_0d_prefben <- function(pb) {
  wd <- tempfile("z"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1); grid <- Setup.Grid(run, tmax = 10, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton")
  pe <- Setup.Pelagic(run, filename = "fish", pref_ben = pb)
  be <- Setup.Benthic(run, filename = "benthos"); de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  f <- Read.In("R", "fish"); b <- Read.In("R", "benthos"); m <- f@mrange
  c(pel = sum(as.numeric(f@finaluvals[1, -(1:3)]) * exp(m) * diff(m)[1], na.rm = TRUE),
    ben = sum(as.numeric(b@finaluvals[1, -(1:3)]) * exp(m) * diff(m)[1], na.rm = TRUE))
}

# final per-cell biomass of a 1-D run (movement off, uniform plankton and initial conditions).
# pb_scalar goes to Setup.Pelagic(); pb_cells (if given) is supplied per cell via Setup.CellParams().
biomass_1d_prefben <- function(ncell, pb_scalar, pb_cells = NULL) {
  wd <- tempfile("pc"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 1, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 10, tstep = 1/48, toutstep = 1,
                     xmin = 0, xmax = (ncell - 1) * 100, xstep = 100, xoutstep = 100)
  uni <- function(sp) { m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep); matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1) }
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE); Setup.ts(pl, run, grid, mat = uni(pl))
  pe <- Setup.Pelagic(run, filename = "fish", pref_ben = pb_scalar, comp = 0, prey = 0, pred = 0, initial_flag = TRUE)
  Setup.ts(pe, run, grid, mat = uni(pe))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); Setup.ts(be, run, grid, mat = uni(be))
  de <- Setup.Detritus(run, filename = "detritus")
  if (!is.null(pb_cells)) Setup.CellParams(pe, run, grid, pref_ben = pb_cells)
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  final <- function(nm) {
    L <- readLines(file.path("R", nm, "results.txt"))
    num <- function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]]))
    hdr <- num(L[5]); mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
    rws <- lapply(L[6:length(L)][nchar(L[6:length(L)]) > 0], num)
    rws <- rws[sapply(rws, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]))]
    last <- rws[sapply(rws, `[`, 1) == max(sapply(rws, `[`, 1))]
    b <- sapply(last, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE)); xx <- sapply(last, `[`, 2)
    unname(sapply(sort(unique(xx)), function(q) b[which(xx == q)[1]]))
  }
  list(pel = final("fish"), ben = final("benthos"))
}

test_that("a uniform per-cell pref_ben is identical to the scalar parameter", {
  a <- biomass_1d_prefben(3, pb_scalar = 0.4)
  b <- biomass_1d_prefben(3, pb_scalar = 1, pb_cells = rep(0.4, 3))   # per-cell value overrides the scalar
  expect_identical(b, a)
})

test_that("per-cell pref_ben makes each cell behave like its own 0-D run", {
  pb  <- c(0.2, 0.4, 0.8)
  pc  <- biomass_1d_prefben(3, pb_scalar = 1, pb_cells = pb)
  b0  <- vapply(pb, biomass_0d_prefben, numeric(2))
  expect_equal(pc$pel, unname(b0["pel", ]), tolerance = 2e-3)
  expect_equal(pc$ben, unname(b0["ben", ]), tolerance = 2e-3)
  expect_true(all(diff(pc$ben) < 0))        # stronger feeding on benthos -> less benthic biomass
})

test_that("Setup.CellParams validates its input", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, xmin = 0, xmax = 2, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  pe <- Setup.Pelagic(run, filename = "fish"); be <- Setup.Benthic(run, filename = "benthos")
  expect_error(Setup.CellParams(pe, run, grid, pref_ben = matrix(0.5, 2, 3)), "3 x 2 matrix")
  expect_error(Setup.CellParams(pe, run, grid, pref_ben = matrix(-1, 3, 2)), "non-negative")
  expect_error(Setup.CellParams(be, run, grid, pref_ben = matrix(0.5, 3, 2)), "pelagic species only")
  m <- Setup.CellParams(pe, run, grid, pref_ben = function(x, y) 0.1 * (x + 1) + y)
  expect_equal(dim(m), c(3L, 2L))
  L <- readLines(file.path("R", "Input", "fish_cellpar.txt"))
  expect_identical(L[1], "pref_ben")
  expect_equal(as.numeric(L[-1]), as.vector(t(m)))   # x outer, y inner
})

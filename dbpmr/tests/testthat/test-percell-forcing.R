# Per-cell environmental forcing for spatial_dim>=1 (issue #44).
# FORCING_PERCELL=1 => forcing_ts.txt has x*y rows per timestep (one per cell), so temperature can vary
# across space. Off (default) => one global row/timestep, bit-identical to standard dbpmr.
LN10 <- log(10)

# biomass in a single 0-D cell at a constant temperature multiplier (uses the global forcing path)
biomass_0d <- function(tempeff) {
  wd <- tempfile("z"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1); grid <- Setup.Grid(run, tmax = 15, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton"); pe <- Setup.Pelagic(run, filename = "fish")
  be <- Setup.Benthic(run, filename = "benthos"); de <- Setup.Detritus(run, filename = "detritus")
  d <- file.path("R", "Input"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  nst <- round(grid@tmax / grid@tstep) + 2
  writeLines(c("pel_tempeff,ben_tempeff", rep(sprintf("%.6g,%.6g", tempeff, tempeff), nst)),
             file.path(d, "forcing_ts.txt"))
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  f <- Read.In("R", "fish"); m <- f@mrange
  sum(as.numeric(f@finaluvals[1, -(1:3)]) * exp(m) * diff(m)[1], na.rm = TRUE)
}

# per-cell biomass for a 1D spatial run with a per-cell temperature vector (uniform plankton)
biomass_percell <- function(tcell) {
  Sys.setenv(FORCING_PERCELL = "1"); on.exit(Sys.unsetenv("FORCING_PERCELL"), add = TRUE)
  wd <- tempfile("pc"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  nc <- length(tcell)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 1, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 15, tstep = 1/48, toutstep = 1,
                     xmin = 0, xmax = (nc - 1) * 100, xstep = 100, xoutstep = 100)
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE)
  Setup.ts(pl, run, grid, func = function(m, t, x, y) 0.01 * exp(-1 * m))        # uniform plankton
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE)
  Setup.ts(pe, run, grid, func = function(m, t, x, y) (0.01/LN10) * exp(-m))     # uniform pelagic IC
  be <- Setup.Benthic(run, filename = "benthos"); de <- Setup.Detritus(run, filename = "detritus")
  d <- file.path("R", "Input"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  nst <- round(grid@tmax / grid@tstep) + 2
  rows <- unlist(lapply(seq_len(nst), function(s) sprintf("%.6g,%.6g", tcell, tcell)))  # x rows/timestep
  writeLines(c("pel_tempeff,ben_tempeff", rows), file.path(d, "forcing_ts.txt"))
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  L <- readLines("R/fish/results.txt")
  hdr <- suppressWarnings(as.numeric(strsplit(L[5], ",")[[1]])); mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
  rws <- lapply(L[6:length(L)][nchar(L[6:length(L)]) > 0],
                function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]])))
  rws <- rws[sapply(rws, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]))]
  tmax <- max(sapply(rws, `[`, 1)); last <- rws[sapply(rws, `[`, 1) == tmax]
  b <- sapply(last, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE)); xx <- sapply(last, `[`, 2)
  xo <- sort(unique(xx)); unname(sapply(xo, function(q) b[which(xx == q)[1]]))
}

test_that("per-cell temperature forcing makes each cell behave like its own 0-D run", {
  tcell <- c(0.5, 1.0, 1.5)
  bpc <- biomass_percell(tcell)
  b0  <- vapply(tcell, biomass_0d, numeric(1))
  expect_equal(bpc, b0, tolerance = 2e-3)   # each cell matches the matching 0-D run at its temperature
  expect_true(all(diff(bpc) < 0))           # warmer -> less biomass (monotone across the gradient)
})

# Per-cell connectivity/immigration floor (PEL_IMM_FRAC) for spatial_dim>=1.
# The floor keeps recruitment >= imm_frac * (running-peak recruitment). That peak reference is now
# tracked PER CELL (rec0_c[k][l]); it used to be a single per-species scalar, which coupled cells in
# a spatial run: a productive cell's peak would inflate a poor neighbour's floor. With movement off,
# a cell is physically independent of its neighbours, so its biomass must NOT change when a neighbour
# gets richer. The old scalar bug would prop the poor cell up; the per-cell reference must not.
LN10 <- log(10)

# biomass in the POOR cell (cell index 1, at x=dx) of a 2-cell, movement-off run with the
# connectivity floor on. cell 0 (at x=0) gets plankton scaled by rich_mult; cell 1 is always poor.
poor_cell_biomass <- function(rich_mult) {
  # imm_frac near 1 + heavy size-selective fishing drives the fished equilibrium below the floor,
  # so the floor is genuinely active (verified: floor-on != floor-off). A milder regime leaves the
  # floor slack and the test would not discriminate the bug.
  Sys.setenv(PEL_IMM_FRAC = "0.9"); on.exit(Sys.unsetenv("PEL_IMM_FRAC"), add = TRUE)
  wd <- tempfile("cf"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  dx <- 100
  run  <- Setup.Run("R", 1, 1, spatial_dim = 1, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 20, tstep = 1/48, toutstep = 1,
                     xmin = 0, xmax = dx, xstep = dx, xoutstep = dx)   # 2 cells: x=0 (rich), x=dx (poor)
  # heterogeneous plankton: cell 0 scaled by rich_mult, cell 1 poor. Movement off, so the only
  # possible cell-0 -> cell-1 influence is via the (now per-cell) recruitment-floor reference.
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE)
  Setup.ts(pl, run, grid, func = function(m, t, x, y) (if (x < dx/2) rich_mult else 1) * 0.01 * exp(-m))
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE,
                      fishing_flag = TRUE)
  Setup.ts(pe, run, grid, func = function(m, t, x, y) (0.01/LN10) * exp(-m))   # uniform pelagic IC
  # strong size-selective fishing depresses recruitment below imm_frac*peak so the floor BINDS
  # (its real use case) - otherwise recruitment rises monotonically and the floor is never active.
  Setup.fishing(pe, run, grid, func = function(m, t, x, y) 6 * as.numeric(m >= LN10))
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE)
  Setup.ts(be, run, grid, func = function(m, t, x, y) 0.01 * exp(-m))          # uniform benthic IC
  de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  L <- readLines("R/fish/results.txt")
  hdr <- suppressWarnings(as.numeric(strsplit(L[5], ",")[[1]])); mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
  rows <- lapply(L[6:length(L)][nchar(L[6:length(L)]) > 0],
                 function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]])))
  rows <- rows[sapply(rows, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]))]
  tout <- max(sapply(rows, `[`, 1)); last <- rows[sapply(rows, `[`, 1) == tout]
  b <- sapply(last, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE)); xx <- sapply(last, `[`, 2)
  b[which(xx == max(xx))[1]]   # the poor cell (largest x)
}

test_that("connectivity floor is per-cell: a richer neighbour does not prop up a poor cell (movement off)", {
  b_poor_rich_nbr <- poor_cell_biomass(rich_mult = 5)   # poor cell next to a rich cell
  b_poor_poor_nbr <- poor_cell_biomass(rich_mult = 1)   # poor cell next to an equally-poor cell
  # With a per-cell floor reference these are identical; the old shared scalar would make the
  # first strictly larger (the rich neighbour's recruitment peak would raise the poor cell's floor).
  expect_equal(b_poor_rich_nbr, b_poor_poor_nbr, tolerance = 1e-6)
})

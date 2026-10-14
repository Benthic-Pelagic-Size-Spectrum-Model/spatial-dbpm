# 2-D coupled movement (spatial_dim == 2). Regression tests for two engine fixes:
#  (1) the spatial_dim==2 branch used to be nested inside `if(spatial_dim==1)`, so NEITHER the x nor
#      the y move solver ran for a 2-D run — movement was completely dead. These tests fail (no
#      effect) against that bug.
#  (2) the explicit upwind advection is CFL sub-stepped, so a large prey/pred coefficient saturates
#      instead of silently flushing the field to zero.
LN10 <- log(10)

# One spatial_dim=2 run on an n x n grid; returns the final per-cell depth-integrated pelagic biomass
# as an n x n matrix. Plankton/benthic/pelagic ICs are all set per cell (no default x_start bump).
field2d <- function(n = 7, dx = 1, tmax = 1, comp = 0, prey = 0, pred = 0, plfun, pefun) {
  wd <- tempfile("m2"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1); L <- (n - 1) * dx
  grid <- Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = tmax,
                     xmin = 0, xmax = L, xstep = dx, xoutstep = dx, ymin = 0, ymax = L, ystep = dx, youtstep = dx)
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE); Setup.ts(pl, run, grid, func = plfun)
  pe <- Setup.Pelagic(run, filename = "fish", comp = comp, prey = prey, pred = pred,
                      gamma_comp = 0, gamma_prey = 0, gamma_pred = 0, initial_flag = TRUE)
  Setup.ts(pe, run, grid, func = pefun)
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE); Setup.ts(be, run, grid, func = function(m, t, x, y) 0.01 * exp(-m))
  de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  L2 <- readLines("R/fish/results.txt"); hdr <- suppressWarnings(as.numeric(strsplit(L2[5], ",")[[1]]))
  mc <- which(!is.na(hdr)); m <- hdr[mc]; dd <- diff(m)[1]
  rows <- lapply(L2[6:length(L2)][nchar(L2[6:length(L2)]) > 0], function(z) suppressWarnings(as.numeric(strsplit(z, ",")[[1]])))
  rows <- rows[sapply(rows, function(v) length(v) >= max(mc) && !is.na(v[1]) && !is.na(v[2]) && !is.na(v[3]))]
  tv <- round(sapply(rows, `[`, 1), 6); xv <- sapply(rows, `[`, 2); yv <- sapply(rows, `[`, 3)
  b  <- sapply(rows, function(v) sum(v[mc] * exp(m) * dd, na.rm = TRUE))
  ux <- sort(unique(xv)); uy <- sort(unique(yv)); last <- tv == max(tv)
  M <- matrix(NA, length(ux), length(uy))
  for (j in which(last)) M[match(xv[j], ux), match(yv[j], uy)] <- b[j]
  M
}

test_that("2-D diffusion smooths a central bump (fails if the spatial_dim==2 branch never runs)", {
  flat <- function(m, t, x, y) 0.02 * exp(-m)
  bump <- function(m, t, x, y) { r2 <- (x - 3)^2 + (y - 3)^2; (0.05 * (1 + 2 * exp(-r2 / 1.5))) * (1/LN10) * exp(-m) }
  ctr <- 4   # centre index of a 7-cell axis
  off <- field2d(n = 7, tmax = 1, comp = 0,   plfun = flat, pefun = bump)
  on  <- field2d(n = 7, tmax = 1, comp = 500, plfun = flat, pefun = bump)
  expect_false(anyNA(off)); expect_false(anyNA(on)); expect_true(all(is.finite(on)))
  r_off <- off[ctr, ctr] / off[1, ctr]    # centre-to-edge contrast without diffusion
  r_on  <- on[ctr, ctr]  / on[1, ctr]      # ... with diffusion: should be closer to 1 (bump smoothed)
  expect_gt(r_off, 1)                       # a bump is present with no movement
  expect_lt(abs(r_on - 1), abs(r_off - 1)) # diffusion flattens it (dead 2-D movement => r_on == r_off)
})

test_that("2-D prey-seeking advection concentrates fish on the food-rich side", {
  grad <- function(m, t, x, y) (0.6 + 0.9 * (1 - x / 6)) * 0.02 * exp(-m)   # west (x=0) rich, east poor
  uni  <- function(m, t, x, y) 0.05 * (1/LN10) * exp(-m)
  we <- function(M) { xm <- rowMeans(M, na.rm = TRUE); mean(xm[1:2]) / mean(xm[6:7]) }
  off <- field2d(n = 7, tmax = 2, prey = 0, plfun = grad, pefun = uni)
  on  <- field2d(n = 7, tmax = 2, prey = 5, plfun = grad, pefun = uni)
  expect_false(anyNA(on)); expect_true(all(is.finite(on)))
  expect_gt(we(on), we(off) * 1.05)   # movement raises the west:east ratio (dead movement => equal)
})

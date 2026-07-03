# Environmental forcing plugin (forcing_ts.txt)
# ---------------------------------------------
# The engine auto-detects <run>/Input/forcing_ts.txt: a header line naming channels
# (pel_tempeff, ben_tempeff, sinking_rate) then one row per timestep. Each timestep it
# scales feeding (A) and background mortality (mu_0) from base values by pel/ben_tempeff
# (mirroring sizemodel's pel/ben_tempeffect) and scales surface-origin detritus inputs
# in g_det by sinking_rate. Missing channels fall back to neutral (1).

LN10 <- log(10)

# minimal pelagic+benthic run. `forcing` = named list of constant channel values written
# to forcing_ts.txt (or NULL for no file). Returns c(pel=biomass, ben=biomass).
run_biomass <- function(A_pel, A_ben, mu_pel, mu_ben, forcing = NULL) {
  wd <- tempfile("forcetest"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = 15, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton", u_0 = 10^(-0.5)/LN10, lambda = -1)
  pe <- Setup.Pelagic(run, filename = "fish", A = A_pel, mu_0 = mu_pel,
                      K_pla = 0.2, R_pla = 0.1, Ex_pla = 0.5, rep_method = 2)
  be <- Setup.Benthic(run, filename = "benthos", A = A_ben, mu_0 = mu_ben,
                      K_det = 0.1, R_det = 0.2, Ex_det = 0.4, rep_method = 2)
  de <- Setup.Detritus(run, filename = "detritus")
  if (!is.null(forcing)) {
    d <- file.path("R", "Input"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
    nst <- round(grid@tmax / grid@tstep) + 2
    hdr <- paste(names(forcing), collapse = ",")
    rows <- do.call(paste, c(lapply(forcing, function(v) sprintf("%.10g", rep_len(v, nst))), sep = ","))
    writeLines(c(hdr, rows), file.path(d, "forcing_ts.txt"))
  }
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  f <- Read.In("R", "fish"); b <- Read.In("R", "benthos")
  m <- f@mrange; dm <- diff(m)[1]
  c(pel = sum(as.numeric(f@finaluvals[1, -(1:3)]) * exp(m) * dm),
    ben = sum(as.numeric(b@finaluvals[1, -(1:3)]) * exp(m) * dm))
}

test_that("constant temperature forcing reproduces folded A/mu_0", {
  te_pel <- 0.85; te_ben <- 0.31
  b_fold   <- run_biomass(64 * te_pel, 6.4 * te_ben, 0.2 * te_pel, 0.2 * te_ben, forcing = NULL)
  b_driver <- run_biomass(64, 6.4, 0.2, 0.2,
                          forcing = list(pel_tempeff = te_pel, ben_tempeff = te_ben))
  # Near-exact: the only (intended) departure is that g_det uses the BASE, un-temperatured
  # mortality for detritus input - which the folded run cannot replicate (it bakes te into
  # mu_0_base). That leaks a ~1e-4 difference through detritus coupling; feeding+spectrum
  # mortality match to machine precision.
  expect_equal(b_driver[["pel"]], b_fold[["pel"]], tolerance = 2e-3)
})

test_that("temperature forcing changes rates (tempeff=1 differs from folded)", {
  te_pel <- 0.85; te_ben <- 0.31
  b_fold <- run_biomass(64 * te_pel, 6.4 * te_ben, 0.2 * te_pel, 0.2 * te_ben, forcing = NULL)
  b_neut <- run_biomass(64, 6.4, 0.2, 0.2, forcing = NULL)  # no file -> tempeff=1
  expect_false(isTRUE(all.equal(b_neut[["pel"]], b_fold[["pel"]], tolerance = 1e-3)))
})

test_that("sinking_rate < 1 lowers detritus-fed detritivore biomass", {
  # less surface detritus reaching the seafloor -> smaller detritus pool -> less benthos
  b_full <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(sinking_rate = 1.0))
  b_half <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(sinking_rate = 0.3))
  expect_lt(b_half[["ben"]], b_full[["ben"]])
})

test_that("sinking_rate = 1 matches no-forcing (neutral default)", {
  b_none <- run_biomass(64, 6.4, 0.2, 0.2, forcing = NULL)
  b_one  <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(sinking_rate = 1.0))
  expect_equal(b_one[["ben"]], b_none[["ben"]], tolerance = 1e-6)
})

test_that("burial is opt-in via the depth channel (absent -> no burial)", {
  # no depth channel -> burial disabled -> identical to plain no-forcing detritus
  b_none  <- run_biomass(64, 6.4, 0.2, 0.2, forcing = NULL)
  b_sink1 <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(sinking_rate = 1))  # no depth col
  expect_equal(b_sink1[["ben"]], b_none[["ben"]], tolerance = 1e-6)
})

test_that("Dunne-2007 burial (opt-in via depth channel) lowers detritus-fed benthos", {
  # supplying the depth channel enables burial -> less detritus reaches the pool -> less benthos
  b_none <- run_biomass(64, 6.4, 0.2, 0.2, forcing = NULL)             # burial off
  b_bur  <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(depth = 1))  # burial on
  expect_lt(b_bur[["ben"]], b_none[["ben"]])
})

test_that("burial matches sizemodel: applied to per-volume flux (depth value does NOT rescale)", {
  # dbpmr densities are per-VOLUME; sizemodel applies Dunne to its per-volume input_w directly,
  # so the depth value must not change the burial fraction - only its presence toggles burial on
  b_d1   <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(depth = 1))
  b_d500 <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(depth = 500))
  expect_equal(b_d500[["ben"]], b_d1[["ben"]], tolerance = 1e-9)
})

test_that("size-dependent export is opt-in (absent -> unchanged detritus dynamics)", {
  # no export channel -> flat sinking_rate, identical to a plain sinking_rate=1 run
  b_flat <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(sinking_rate = 1))
  b_none <- run_biomass(64, 6.4, 0.2, 0.2, forcing = NULL)
  expect_equal(b_flat[["ben"]], b_none[["ben"]], tolerance = 1e-6)
})

test_that("size-dependent export: stronger attenuation lowers detritus-fed benthos", {
  # export(m) = exp(-export_attn * exp(-export_gamma*m)); larger export_attn -> less reaches the bed
  b_lo <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(export_attn = 0.1, export_gamma = 0.4))
  b_hi <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(export_attn = 5.0, export_gamma = 0.4))
  expect_lt(b_hi[["ben"]], b_lo[["ben"]])
})

# fished run; base_* = per-group base fishing mortality above 10 g; gravity toggles the
# in-engine biomass-proportional split. Returns fishable areal pelagic/benthic biomass.
run_fished <- function(base_pel, base_ben, gravity = FALSE) {
  wsel <- 1 * LN10
  wd <- tempfile("grav"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = 30, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton", u_0 = 10^(-0.5)/LN10, lambda = -1)
  pe <- Setup.Pelagic(run, filename = "fish", A = 64, mu_0 = 0.2, K_pla = 0.2, R_pla = 0.1,
                      Ex_pla = 0.5, rep_method = 2, fishing_flag = TRUE)
  be <- Setup.Benthic(run, filename = "benthos", A = 6.4, mu_0 = 0.2, K_det = 0.1, R_det = 0.2,
                      Ex_det = 0.4, rep_method = 2, fishing_flag = TRUE)
  de <- Setup.Detritus(run, filename = "detritus")
  Setup.fishing(pe, run, grid, func = function(m, t, x, y) base_pel * as.numeric(m >= wsel))
  Setup.fishing(be, run, grid, func = function(m, t, x, y) base_ben * as.numeric(m >= wsel))
  if (gravity) {
    d <- file.path("R", "Input"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
    nst <- round(grid@tmax / grid@tstep) + 2
    writeLines(c("gravity,depth", rep("1,200", nst)), file.path(d, "forcing_ts.txt"))
  }
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))
  f <- Read.In("R", "fish"); b <- Read.In("R", "benthos"); m <- f@mrange; dm <- diff(m)[1]
  fi <- (m / LN10) >= 1
  c(pel = sum(as.numeric(f@finaluvals[1, -(1:3)])[fi] * exp(m[fi]) * dm) * 200,
    ben = sum(as.numeric(b@finaluvals[1, -(1:3)])[fi] * exp(m[fi]) * dm) * 20)
}

test_that("gravity fishing concentrates effort on the dominant group", {
  # pelagic dominates the fishable biomass, so the gravity split (of a 0.4 budget) puts ~all
  # effort on pelagic -> it is fished HARDER than an even split of the same budget.
  even <- run_fished(0.2, 0.2, gravity = FALSE)   # even split: F = 0.2 each
  grav <- run_fished(0.4, 0.4, gravity = TRUE)    # gravity splits 0.4 by current biomass share
  expect_lt(grav[["pel"]], even[["pel"]])
})

test_that("residence-time closure: shorter tau -> less detritus -> less benthos, and is opt-in", {
  # first-order pool loss W/tau: shorter residence time removes detritus faster -> less benthos
  b_none  <- run_biomass(64, 6.4, 0.2, 0.2, forcing = NULL)                        # closure off
  b_long  <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(residence_time = 100))  # weak loss
  b_short <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(residence_time = 0.1))  # strong loss
  expect_lt(b_short[["ben"]], b_long[["ben"]])
  # tau=0 disables it (identical to no closure)
  b_zero  <- run_biomass(64, 6.4, 0.2, 0.2, forcing = list(residence_time = 0))
  expect_equal(b_zero[["ben"]], b_none[["ben"]], tolerance = 1e-9)
})

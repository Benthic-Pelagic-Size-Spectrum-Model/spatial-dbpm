# Species set up from a table of parameters (issue #58, step 8), and the predator-avoidance warning
# in Setup.Pelagic() (issue #58, step 6).
LN10 <- log(10)

in_tmp <- function(code) { wd <- tempfile("sp"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE); force(code) }

sp_table <- function() data.frame(
  species  = c("fish", "benthos"),
  type     = c("pelagic", "benthic"),
  mmat     = c(2, 1) * LN10,
  mmax     = c(6, 4) * LN10,
  mu_0     = c(0.15, 0.1),
  pref_ben = c(0.5, NA),       # pelagic only
  K_det    = c(NA, 0.25),      # benthic only
  stringsAsFactors = FALSE)

test_that("a species table builds the same objects as direct Setup calls", {
  in_tmp({
    run <- Setup.Run("R", 1, 1, 0, TRUE, 1)
    g <- Setup.Species(run, sp_table())
    expect_named(g, c("pelagic", "benthic")); expect_named(g$pelagic, "fish"); expect_named(g$benthic, "benthos")
    expect_identical(g$pelagic$fish, Setup.Pelagic(run, filename = "fish", mmat = 2 * LN10, mmax = 6 * LN10, mu_0 = 0.15, pref_ben = 0.5))
    expect_identical(g$benthic$benthos, Setup.Benthic(run, filename = "benthos", mmat = 1 * LN10, mmax = 4 * LN10, mu_0 = 0.1, K_det = 0.25))
    expect_identical(g$pelagic$fish@A, 640)                      # a parameter not in the table keeps its default
    expect_true(dir.exists(file.path("R", "fish")))
  })
})

test_that("a run built from a species table is identical to one built with direct calls", {
  run_it <- function(from_table) in_tmp({
    run  <- Setup.Run("R", 1, 1, 0, TRUE, 1); grid <- Setup.Grid(run, tmax = 2, tstep = 1/48, toutstep = 1)
    pl <- Setup.Plankton(run, filename = "plankton"); de <- Setup.Detritus(run, filename = "detritus")
    if (from_table) { g <- Setup.Species(run, sp_table()); pe <- g$pelagic; be <- g$benthic }
    else { pe <- Setup.Pelagic(run, filename = "fish", mmat = 2 * LN10, mmax = 6 * LN10, mu_0 = 0.15, pref_ben = 0.5)
           be <- Setup.Benthic(run, filename = "benthos", mmat = 1 * LN10, mmax = 4 * LN10, mu_0 = 0.1, K_det = 0.25) }
    out <- NULL; invisible(capture.output(out <- SizeSpectrum(run, grid, pl, pe, be, de, output = "arrays")))
    out
  })
  a <- run_it(TRUE); b <- run_it(FALSE)
  expect_identical(a, b)
  expect_gt(sum(a$abundance$fish[1, 1, , 3]), 0)
})

test_that("the table can come from a CSV file, with flags as 0/1 or TRUE/FALSE", {
  in_tmp({
    run <- Setup.Run("R", 2, 1, 0, TRUE, 1)
    tab <- data.frame(species = c("small", "large", "benthos"), type = c("pelagic", "pelagic", "benthic"),
                      mmax = c(3, 6, 4) * LN10, mmat = c(1, 3, 1) * LN10, fishing_flag = c(0, 1, 0), initial_flag = c(FALSE, FALSE, TRUE))
    write.csv(tab, "species_params.csv", row.names = FALSE)
    g <- Setup.Species(run, "species_params.csv")
    expect_named(g$pelagic, c("small", "large"))
    expect_identical(g$pelagic$large@fishing_flag, TRUE); expect_identical(g$pelagic$small@fishing_flag, FALSE)
    expect_identical(g$benthic$benthos@initial_flag, TRUE)
    expect_equal(g$pelagic$large@mmax, 6 * LN10)
  })
})

test_that("the Kerguelen biology constants can be supplied as a species table", {
  p <- jsonlite::fromJSON(system.file("extdata", "kerguelen_biology_params.json", package = "dbpmr"))$params
  dh <- p$defecate_prop; dl <- p$def_low
  Kp <- (1 - dh) * p$growth_pred; Rp <- (1 - dh) * (1 - (p$growth_pred + p$energy_pred)); Ep <- (1 - dh) * p$energy_pred
  Kl <- (1 - dl) * p$growth_detritivore; Rl <- (1 - dl) * (1 - (p$growth_detritivore + p$energy_detritivore)); El <- (1 - dl) * p$energy_detritivore
  in_tmp({
    run <- Setup.Run("R", 1, 1, 0, TRUE, 1)
    tab <- data.frame(species = c("fish", "benthos"), type = c("pelagic", "benthic"),
                      mmin = -3 * LN10, mmat = c(2, 1) * LN10, mmax = c(6, 4) * LN10,
                      alpha = c(p$metabolic_req_pred, p$metabolic_req_detritivore), A = c(64, 6.4) / 3, mu_0 = p$natural_mort * 0.5,
                      pref_ben = c(0.4, NA), K_pla = c(Kp, NA), R_pla = c(Rp, NA), Ex_pla = c(Ep, NA),
                      K_pel = c(Kp, NA), R_pel = c(Rp, NA), Ex_pel = c(Ep, NA), K_ben = c(Kl, NA), R_ben = c(Rl, NA), Ex_ben = c(El, NA),
                      K_det = c(NA, Kl), R_det = c(NA, Rl), Ex_det = c(NA, El), rep_method = 2, fishing_flag = FALSE)
    g <- Setup.Species(run, tab)
    expect_identical(g$pelagic$fish, Setup.Pelagic(run, filename = "fish", mmin = -3 * LN10, mmat = 2 * LN10, mmax = 6 * LN10,
        alpha = p$metabolic_req_pred, A = 64 / 3, mu_0 = p$natural_mort * 0.5, pref_ben = 0.4, K_pla = Kp, R_pla = Rp, Ex_pla = Ep,
        K_pel = Kp, R_pel = Rp, Ex_pel = Ep, K_ben = Kl, R_ben = Rl, Ex_ben = El, rep_method = 2, fishing_flag = FALSE))
    expect_identical(g$benthic$benthos, Setup.Benthic(run, filename = "benthos", mmin = -3 * LN10, mmat = 1 * LN10, mmax = 4 * LN10,
        alpha = p$metabolic_req_detritivore, A = 6.4 / 3, mu_0 = p$natural_mort * 0.5, K_det = Kl, R_det = Rl, Ex_det = El,
        rep_method = 2, fishing_flag = FALSE))
  })
})

test_that("Setup.Species validates the table", {
  in_tmp({
    run <- Setup.Run("R", 1, 1, 0, TRUE, 1); tab <- sp_table()
    expect_error(Setup.Species(run, tab[, -1]), "'species' and 'type'")
    expect_error(Setup.Species(run, transform(tab, type = c("pelagic", "plankton"))), "pelagic.*benthic")
    expect_error(Setup.Species(run, transform(tab, species = "fish")), "unique")
    expect_error(Setup.Species(run, transform(tab, w_inf = 1)), "unknown column in species_params: w_inf")
    expect_error(Setup.Species(run, transform(tab, pref_ben = 0.5)), "'pref_ben' does not apply to a benthic species")
    expect_error(Setup.Species(run, tab[1, ]), "1 pelagic and 0 benthic species but the run was set up for 1 and 1")
    expect_error(Setup.Species(run, "no_such_file.csv"), "not found")
    expect_error(Setup.Species(run, list(a = 1)), "data frame")
    expect_warning(Setup.Species(run, transform(tab, mmat = c(7, 1) * LN10)), "no mature size classes")   # checks of the Setup functions still apply
  })
})

test_that("Setup.Pelagic warns about predator-avoidance without diffusion in a spatial run", {
  in_tmp({
    run2 <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
    expect_warning(Setup.Pelagic(run2, filename = "fish", pred = 0.2, comp = 0), "predator-avoidance movement without diffusion")
    expect_no_warning(Setup.Pelagic(run2, filename = "fish", pred = 0.2, comp = 0.2))     # diffusion present
    expect_no_warning(Setup.Pelagic(run2, filename = "fish", prey = 0.2, pred = 0, comp = 0))   # prey-seeking alone is fine
    expect_no_warning(Setup.Pelagic(run2, filename = "fish"))                              # defaults
  })
  in_tmp({
    run0 <- Setup.Run("R", 1, 1, 0, TRUE, 1)
    expect_no_warning(Setup.Pelagic(run0, filename = "fish", pred = 0.2, comp = 0))       # no movement in a 0-D run
  })
})

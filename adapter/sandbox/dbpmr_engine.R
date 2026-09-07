# dbpmr engine adapter for the LME/FishMIP workflow (sandbox, Tier-1)
# -------------------------------------------------------------------
# run_model_dbpmr() honours the same call signature and (a subset of) the return
# structure as the LME workflow's run_model(), but runs the *dbpmr* C engine as
# the integrator instead of sizemodel(). This is the "engine = dbpmr" override:
# the pipeline (sizeparam, forcing, gravity, calibration) is imported unchanged
# from the lme-workflow submodule; only the integrator is swapped.
#
# Scope (Tier-1): aspatial, stable-spin (constant) forcing, no fishing yet.
# Parameters are derived by the workflow's own sizeparam() so the two engines
# use identical biology; only the numerical engine differs.
#
# Requires: an installed dbpmr; sizeparam() from the sourced lme-workflow.

run_model_dbpmr <- function(fishing_params, dbpm_inputs,
                            xmin_consumer_u = -3, xmin_consumer_v = -3,
                            tstep = 1/48, tmax = 100, ...) {
  stopifnot(exists("sizeparam"))                 # from the lme-workflow submodule
  LN10 <- log(10)
  p <- sizeparam(dbpm_inputs, fishing_params,
                 xmin_consumer_u = xmin_consumer_u, xmin_consumer_v = xmin_consumer_v)
  f1 <- function(v) v[[1]]                        # stable-spin forcing is constant

  # --- Boltzmann-Arrhenius temperature (folded into constant A/mu_0 for the
  #     stable spin; transient forcing is Tier-2 / issue #11) ---
  te <- function(t) exp(f1(p$c1) - f1(p$activation_energy) /
                          (f1(p$boltzmann) * (t + 273)))
  pt <- te(f1(p$sea_surf_temp)); bt <- te(f1(p$sea_floor_temp))

  # --- faithful K/R/Ex energy budget (issue #22), keyed to prey type ---
  dh <- f1(p$defecate_prop); dl <- f1(p$def_low)
  Ku <- f1(p$growth_pred);        AMu <- f1(p$energy_pred)
  Kv <- f1(p$growth_detritivore); AMv <- f1(p$energy_detritivore)
  Kp <- (1-dh)*Ku; Rp <- (1-dh)*(1-(Ku+AMu)); Ep <- (1-dh)*AMu
  Kl <- (1-dl)*Kv; Rl <- (1-dl)*(1-(Kv+AMv)); El <- (1-dl)*AMv

  S  <- f1(p$hr_volume_search)                    # search volume (calibration knob)
  depth <- f1(p$depth); dcorr <- min(depth, 200); bhd <- 20  # areal output multipliers
  sink0 <- f1(p$sinking_rate)                     # export ratio (fraction reaching seafloor)
  # predator-benthos coupling on a vertical-MIGRATION length (~1500 m), not the passive
  # sinking length (250 m): the legacy 0.8*exp(-depth/250) zeroes coupling by ~1 km and
  # collapses deep-system pelagics; migration keeps them weakly coupled. See tier1_fishing_calib.R.
  pref_ben_depth <- 0.8 * exp(-depth / 1500)

  wd <- tempfile("dbpmr_eng"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  run  <- dbpmr::Setup.Run("R", 1, 1, spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  grid <- dbpmr::Setup.Grid(run, tmax = tmax, tstep = tstep, toutstep = 1)
  # LME `intercept` is a per-log10-density intercept (sizemodel convention:
  # rho = 10^intercept * w^slope, integrated with d(log10)). dbpmr's plankton
  # u_values are per-LN density, so convert: u_0 = 10^intercept / ln(10)
  # (verified by plankton-biomass equality; slope is base-invariant). Omitting
  # the /ln(10) inflates the resource seen by predators by ln(10) ~ 2.3x.
  pl <- dbpmr::Setup.Plankton(run, filename = "plankton",
          u_0 = 10^f1(p$int_phy_zoo) / LN10, lambda = f1(p$slope_phy_zoo))
  pe <- dbpmr::Setup.Pelagic(run, filename = "fish", alpha = f1(p$metabolic_req_pred),
          A = S * pt, mu_0 = f1(p$natural_mort) * pt, pref_ben = pref_ben_depth,
          K_pla = Kp, R_pla = Rp, Ex_pla = Ep, K_pel = Kp, R_pel = Rp, Ex_pel = Ep,
          K_ben = Kl, R_ben = Rl, Ex_ben = El, rep_method = 2)
  be <- dbpmr::Setup.Benthic(run, filename = "benthos", alpha = f1(p$metabolic_req_detritivore),
          A = 0.1 * S * bt, mu_0 = f1(p$natural_mort) * bt,
          K_det = Kl, R_det = Rl, Ex_det = El, rep_method = 2)
  de <- dbpmr::Setup.Detritus(run, filename = "detritus")
  # Same detritus forcings as sizemodel (stable-spin = constant): sinking_rate (export ratio)
  # scales surface detritus input, depth enables Dunne-2007 burial. Temperature stays folded
  # into A/mu_0 (constant here). Written as a 1-row-per-timestep constant forcing_ts.txt.
  nst <- round(tmax / tstep) + 2
  di_in <- file.path("R", "Input"); dir.create(di_in, showWarnings = FALSE, recursive = TRUE)
  writeLines(c("sinking_rate,depth", rep(sprintf("%.8g,%.8g", sink0, dcorr), nst)),
             file.path(di_in, "forcing_ts.txt"))
  invisible(capture.output(dbpmr::SizeSpectrum(run, grid, pl, pe, be, de)))

  fsh <- dbpmr::Read.In("R", "fish"); bn <- dbpmr::Read.In("R", "benthos")
  U <- as.numeric(fsh@finaluvals[1, -(1:3)]); V <- as.numeric(bn@finaluvals[1, -(1:3)])
  x <- fsh@mrange / LN10; w <- exp(fsh@mrange); dm <- diff(fsh@mrange)[1]
  cons <- x > xmin_consumer_u
  # dbpmr native u_values are per-LN-mass density; dm is the ln-spacing, so the
  # biomass integral is  sum(u * w * d(ln w))  with NO ln(10) factor (confirmed
  # against the C engine's own integrals, SizeSpectra.c:2263/2293/2336). We only
  # apply ln(10) to the *reported density field* to express it per-log10 (the
  # sizemodel convention) - NOT to the integral, which stays in native ln units.
  list(predators        = U * LN10,               # per-log10 density (like sizemodel U)
       detritivores      = V * LN10,
       size_log10        = x,
       total_pred_biomass = sum(U[cons] * w[cons] * dm) * dcorr,
       total_detritivore_biomass = sum(V[cons] * w[cons] * dm) * bhd,
       params            = p,
       engine            = "dbpmr")
}

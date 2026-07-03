# Tier-1 sandbox: single-Q fishing calibration (transient, one FAO-LME)
# ---------------------------------------------------------------------
# Estimate ONE catchability Q by minimising the difference between the MODELLED
# and OBSERVED catch TIME SERIES (log MSE), given the model + forcing:
#   F_g(x,t) = Q * s(x) * (B_g/sum_g B) * effort_norm(t)
# where
#   s(x)            knife-edge selectivity, 1 for w >= 10 g (log10 w = 1),
#   effort_norm(t)  nominal effort normalised to [0,1] over the LME series,
#   (B_g/sum_g B)   biomass-proportional effort split across functional groups
#                   (DBPM.md "gravity": effort ~ biomass, so catch ~ biomass^2),
#   Q               single catchability, bounded [0,3], A (search vol) fixed = 64.
# Per the FishMIP DBPM.md spec (Fish-MIP/Global_MEM_Model_Templates).
#
# ENVIRONMENTAL FORCING (#11): plankton, temperature AND detritus sinking are TIME-VARYING,
# using the exact same obsclim forcings as sizemodel(), held at year-1 values during spin.
#  - plankton: monthly (10^intercept/ln10)*w^slope via dbpmr's plankton time-series input
#    (ts_flag), re-read every timestep by the C engine.
#  - temperature + sinking: monthly tos/tob -> Boltzmann-Arrhenius tempeffect, and
#    export_ratio -> sinking_rate, written to forcing_ts.txt (the dbpmr forcing plugin:
#    header-named channels, one row per timestep). tempeffect scales feeding A + background
#    mortality mu_0 (senescence/fishing NOT scaled); sinking_rate scales surface-origin
#    detritus inputs in g_det (detritus mortality kept un-temperatured, as in sizemodel).
#  - depth: pref_benthos = 0.8*exp(-depth/250) (predator-benthos coupling) + areal conversion.
#
# UNITS (ln vs log10): the LME `intercept` is a per-LOG10-density intercept
# (sizemodel: rho = 10^intercept * w^slope, integrated with d(log10)). dbpmr's
# u_values are per-LN density, so the plankton intercept converts as
#   u_0 = 10^intercept / ln(10)      (slope/lambda is base-invariant).
# All spectrum integrals below (biomass, catch) use dbpmr-native per-ln density
# with dm = d(ln w) and carry NO ln(10) factor - matching the C engine's own
# integrals (SizeSpectra.c:2263/2293/2336).
#
#   DBPM_DATA=/path/to/DBPM_dev Rscript adapter/sandbox/tier1_fishing_calib.R [LME]

.libPaths(c("/tmp/dbpmrlib", .libPaths()))
suppressMessages({ library(jsonlite); library(dbpmr); library(arrow); library(dplyr); library(stringr) })
LN10 <- log(10)
base <- Sys.getenv("DBPM_DATA",
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev")
L <- as.integer(c(commandArgs(TRUE), "14")[1])

# --- biology params (equilibrium init json) + Boltzmann-Arrhenius temperature ---
p  <- fromJSON(file.path(base, "equilibrium_runs",
        sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_12.8.json", L)))$params
# tempeffect(T) multiplies feeding (A) and background mortality (mu_0). It is now
# TIME-VARYING via dbpmr's temperature driver (temperature_ts.txt, read each timestep)
# rather than folded into constant A/mu_0 - so A/mu_0 below are the BASE (unscaled) rates.
te <- function(T) exp(p$c1[1] - p$activation_energy[1] / (p$boltzmann[1] * (T + 273)))
# faithful K/R/Ex energy budget (issue #22), keyed to prey type
dh <- p$defecate_prop[1]; dl <- p$def_low[1]
Ku <- p$growth_pred[1];        AMu <- p$energy_pred[1]
Kv <- p$growth_detritivore[1]; AMv <- p$energy_detritivore[1]
Kp <- (1-dh)*Ku; Rp <- (1-dh)*(1-(Ku+AMu)); Ep <- (1-dh)*AMu
Kl <- (1-dl)*Kv; Rl <- (1-dl)*(1-(Kv+AMv)); El <- (1-dl)*AMv
dcorr <- min(p$depth[1], 200); bhd <- 20; wsel <- 1 * LN10   # knife-edge at 10 g (ln units)

# --- obsclim forcing: annual effort/catch series + monthly plankton intercept/slope ---
di <- read_parquet(Sys.glob(file.path(base, "dbpm_inputs",
        sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet", L)))[1]) |>
      filter(str_detect(scenario, "obsclim")) |> arrange(year, month)
yr <- di |> group_by(year) |>
      summarise(eff = mean(total_nom_active_area_m2, na.rm = TRUE),
                cat = mean(catch_tonnes_area_m2, na.rm = TRUE) * 1e6, .groups = "drop") |>
      arrange(year)
effn <- yr$eff / max(yr$eff); nyr <- nrow(yr); spin <- 40; tmax <- spin + nyr
eff_at <- function(t) { i <- floor(t) - spin + 1; if (i < 1) effn[1] else if (i > nyr) effn[nyr] else effn[i] }
# INT_OFFSET: additive log10 shift on the plankton intercept, to emulate the MLD-aware
# phyto averaging (deep columns gain ~log10(200/MLD) ~ +0.3..0.6) before the upstream
# preprocessing is re-run. Default 0 (use the parquet's fixed-200 m intercept as-is).
int_off <- as.numeric(Sys.getenv("INT_OFFSET", "0"))
intm <- di$intercept + int_off; slpm <- di$slope; nmon <- length(intm)
int1 <- mean(di$intercept[di$year == min(di$year)]) + int_off; slp1 <- mean(di$slope[di$year == min(di$year)])
pl_at <- function(t) { if (t < spin) return(c(int1, slp1))
  j <- max(1, min(nmon, floor((t - spin) * 12) + 1)); c(intm[j], slpm[j]) }
# per-LN plankton spectrum, UNIT-MATCHED to dbpmr: (10^intercept/ln10) * exp(slope*m)
plfun <- function(m, t, x, y) { ab <- pl_at(t); (10^ab[1] / LN10) * exp(ab[2] * m) }
# --- MONTHLY environmental forcing series (obsclim), spin held at year-1 mean ---
tos <- di$tos; tob <- di$tob                     # surface / seafloor temperature
er  <- di$export_ratio                           # sinking rate (fraction reaching seafloor)
tos1 <- mean(di$tos[di$year == min(di$year)]); tob1 <- mean(di$tob[di$year == min(di$year)])
er1  <- mean(di$export_ratio[di$year == min(di$year)])
midx <- function(t) max(1, min(nmon, floor((t - spin) * 12) + 1))   # obsclim month index at time t
# forcing channels at model-time t
temp_at <- function(t) { if (t < spin) return(c(te(tos1), te(tob1))); j <- midx(t); c(te(tos[j]), te(tob[j])) }
sink_at <- function(t) { if (t < spin) return(er1); er[midx(t)] }
# depth -> predator-benthos coupling. The legacy 0.8*exp(-depth/250) reuses the passive
# SINKING attenuation length (~250 m), so coupling is ~0 by 1000 m and truly zero in the
# deep - far too steep for animals that vertically migrate hundreds of metres (it also
# contradicts sizeparam's own cited intent, Trueman: 0.75 <500 m, 0.5 to 1800 m). Decay
# instead on a vertical-MIGRATION length (~1500 m, mesopelagic + deep migrators): same form
# and 0.8 surface value, but coupling persists across the DVM range and deep systems stay
# weakly coupled (0.08 at 3437 m) rather than losing the pelagic web entirely.
coupling_scale <- 1500
depth_mean <- p$depth[1]; pref_ben_depth <- 0.8 * exp(-depth_mean / coupling_scale)

# --- one simulation at per-group catchabilities (qp pelagic, qb benthic) ---
# single-Q gravity split is qp = Q*s_pel, qb = Q*s_ben; two_Q frees them.
runsim <- function(qp, qb) {
  wd <- tempfile("tv"); dir.create(wd); old <- setwd(wd)
  run  <- Setup.Run("R", 1, 1, 0, TRUE, 1)
  grid <- Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = 1)
  pl   <- Setup.Plankton(run, filename = "plankton", lambda = slp1, ts_flag = TRUE)
  Setup.ts(pl, run, grid, func = plfun)
  fishing <- (qp > 0 || qb > 0)
  # BASE (unscaled) A/mu_0; temperature applied per-timestep by the driver below
  pe <- Setup.Pelagic(run, filename = "fish", mmin = -3*LN10, mmat = 2*LN10, mmax = 6*LN10,
          alpha = p$metabolic_req_pred[1], A = 64, mu_0 = p$natural_mort[1], pref_ben = pref_ben_depth,
          K_pla = Kp, R_pla = Rp, Ex_pla = Ep, K_pel = Kp, R_pel = Rp, Ex_pel = Ep,
          K_ben = Kl, R_ben = Rl, Ex_ben = El, rep_method = 2, fishing_flag = fishing)
  be <- Setup.Benthic(run, filename = "benthos", mmin = -3*LN10, mmat = 1*LN10, mmax = 4*LN10,
          alpha = p$metabolic_req_detritivore[1], A = 6.4, mu_0 = p$natural_mort[1],
          K_det = Kl, R_det = Rl, Ex_det = El, rep_method = 2, fishing_flag = fishing)
  de <- Setup.Detritus(run, filename = "detritus")
  # environmental forcing plugin: write forcing_ts.txt (header + one row per timestep, row j ->
  # model time j*tstep), auto-detected by the C engine. Channels: pel/ben temperature effect
  # (scale feeding A + background mortality mu_0) and sinking_rate (scale surface detritus input).
  nst  <- round(tmax / (1/48)) + 2
  ts_t <- seq_len(nst) * (1/48)
  tmat <- vapply(ts_t, temp_at, numeric(2)); smat <- vapply(ts_t, sink_at, numeric(1))
  di_in <- file.path("R", "Input"); dir.create(di_in, showWarnings = FALSE, recursive = TRUE)
  # detritus closure (env DET_CLOSURE): "burial" (default) -> Dunne via the depth channel;
  # "residence:TAU" -> first-order pool loss W/tau via residence_time;
  # "export:GAMMA:ATTN0" -> SIZE-DEPENDENT export: a surface particle of ln-mass m reaches the
  #   seafloor with fraction exp(-attn(t)*exp(-GAMMA*m)) (large carcasses/pellets sink faster).
  #   attn(t) = ATTN0 * (dcorr/200) * exp(0.05*(tob(t)-mean tob))  (deeper/warmer -> more attenuation).
  cl <- Sys.getenv("DET_CLOSURE", "burial")
  if (grepl("^residence", cl)) {
    tau <- as.numeric(sub("^residence:?", "", cl)); if (is.na(tau)) tau <- 1
    hdr <- "pel_tempeff,ben_tempeff,sinking_rate,residence_time"
    rows <- sprintf("%.8g,%.8g,%.8g,%.8g", tmat[1, ], tmat[2, ], smat, tau)
  } else if (grepl("^export", cl)) {
    prm  <- as.numeric(strsplit(sub("^export:?", "", cl), ":")[[1]])
    gam  <- if (length(prm) >= 1 && !is.na(prm[1])) prm[1] else 0.4
    att0 <- if (length(prm) >= 2 && !is.na(prm[2])) prm[2] else 1
    tobt <- vapply(ts_t, function(t) tob[midx(t)], numeric(1)); tobt[ts_t < spin] <- tob1
    attn <- att0 * (dcorr / 200) * exp(0.05 * (tobt - mean(tob)))
    hdr  <- "pel_tempeff,ben_tempeff,sinking_rate,export_attn,export_gamma"
    rows <- sprintf("%.8g,%.8g,%.8g,%.8g,%.8g", tmat[1, ], tmat[2, ], smat, attn, gam)
  } else {
    hdr <- "pel_tempeff,ben_tempeff,sinking_rate,depth"
    rows <- sprintf("%.8g,%.8g,%.8g,%.8g", tmat[1, ], tmat[2, ], smat, dcorr)
  }
  writeLines(c(hdr, rows), file.path(di_in, "forcing_ts.txt"))
  if (fishing) {
    Setup.fishing(pe, run, grid, func = function(m, t, x, y) qp * as.numeric(m >= wsel) * eff_at(t))
    Setup.fishing(be, run, grid, func = function(m, t, x, y) qb * as.numeric(m >= wsel) * eff_at(t))
  }
  ok <- tryCatch({ invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de))); TRUE },
                 error = function(e) { message(conditionMessage(e)); FALSE })
  if (!ok) { setwd(old); return(NULL) }
  f <- Read.In("R", "fish"); b <- Read.In("R", "benthos"); setwd(old)
  list(Um = as.matrix(f@uvals[, -(1:3)]), Vm = as.matrix(b@uvals[, -(1:3)]), tt = f@uvals[, 1],
       x = f@mrange / LN10, w = exp(f@mrange), dm = diff(f@mrange)[1])
}

# --- 1) unfished equilibrium -> fishable-biomass shares (DBPM.md gravity split) ---
u0 <- runsim(0, 0); fi <- u0$x >= 1; nT <- length(u0$tt)
Bpel <- sum(u0$Um[nT, fi] * u0$w[fi] * u0$dm) * dcorr
Bben <- sum(u0$Vm[nT, fi] * u0$w[fi] * u0$dm) * bhd
s_pel <- Bpel / (Bpel + Bben); s_ben <- Bben / (Bpel + Bben)
two_Q <- nzchar(Sys.getenv("TWO_Q"))   # TWO_Q set -> estimate q_pel, q_ben separately
cat(sprintf("=== FAO-LME %d: %s + time-varying plankton ===\n", L,
            if (two_Q) "per-group Q (q_pel,q_ben)" else "single-Q + gravity split"))
cat(sprintf("  unfished fishable biomass  pel=%.3g  ben=%.3g  -> shares  pel=%.2f  ben=%.2f\n",
            Bpel, Bben, s_pel, s_ben))

# --- catch time series at catchabilities (qp,qb) + log-MSE objective ---
catch_ts <- function(qp, qb) {
  r <- runsim(qp, qb); if (is.null(r)) return(rep(NA, nyr))
  k <- match(spin + (1:nyr), r$tt)
  sapply(k, function(kk) { e <- eff_at(r$tt[kk])
    sum(qp * e * r$Um[kk, fi] * r$w[fi] * r$dm) * dcorr +
    sum(qb * e * r$Vm[kk, fi] * r$w[fi] * r$dm) * bhd })
}
logmse <- function(m) { if (all(is.na(m))) return(1e6)
  ok <- is.finite(m) & m > 0 & yr$cat > 0
  if (sum(ok) < 3) return(1e6); mean((log10(m[ok]) - log10(yr$cat[ok]))^2) }

# --- 2) calibrate to the observed catch time series (log-MSE) ---
n <- 0; Qf <- NA_real_
if (!two_Q) {
  # single common Q, gravity split: q_pel = Q*s_pel, q_ben = Q*s_ben (1-D Brent)
  o  <- optimise(function(Q) { n <<- n + 1; logmse(catch_ts(Q * s_pel, Q * s_ben)) },
                 lower = 0, upper = 3, tol = 0.01)
  Qf <- o$minimum; qpf <- Qf * s_pel; qbf <- Qf * s_ben; mse <- o$objective
  cat(sprintf("  optimise[0,3]: evals=%d  Q=%.3f  ", n, Qf))
} else {
  # per-group catchabilities. The objective is bimodal (yield-curve low-F / high-F
  # basins) and non-smooth (knife-edge selectivity), so use a coarse grid seed to find
  # the global basin, then a BOUNDED derivative-free polish. nloptr:
  #   TWO_Q_OPT=grid_bobyqa (default) grid seed + NLOPT_LN_BOBYQA (bounded local)
  #   TWO_Q_OPT=directL               NLOPT_GN_DIRECT_L (deterministic global)
  obj2 <- function(q) { n <<- n + 1; logmse(catch_ts(q[1], q[2])) }
  alg  <- Sys.getenv("TWO_Q_OPT", "grid_bobyqa")
  if (alg == "directL") {
    res <- nloptr::nloptr(c(1, 1), obj2, lb = c(0, 0), ub = c(3, 3),
             opts = list(algorithm = "NLOPT_GN_DIRECT_L", maxeval = 80, xtol_rel = 1e-3))
  } else {
    qs <- seq(0.05, 3, length.out = 5); G <- expand.grid(qp = qs, qb = qs)
    seed <- as.numeric(G[which.min(apply(G, 1, function(q) obj2(as.numeric(q)))), ])
    res  <- nloptr::nloptr(seed, obj2, lb = c(0, 0), ub = c(3, 3),
             opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 40, xtol_rel = 1e-3))
  }
  qpf <- res$solution[1]; qbf <- res$solution[2]; mse <- res$objective
  cat(sprintf("  %s: evals=%d  q_pel=%.3f  q_ben=%.3f  ", alg, n, qpf, qbf))
}
mc <- catch_ts(qpf, qbf)
cat(sprintf("MSE(log catch)=%.3f  corr=%.2f\n",
            mse, suppressWarnings(cor(log10(mc), log10(yr$cat), use = "complete.obs"))))
cat(sprintf("  catch mean: model=%.3f  obs=%.3f (g m-2 yr-1)\n", mean(mc, na.rm = TRUE), mean(yr$cat)))

# optional: save the calibrated catch time series + stats (for batch plotting)
if (nzchar(Sys.getenv("SAVE_RDS"))) {
  cr <- suppressWarnings(cor(log10(mc), log10(yr$cat), use = "complete.obs"))
  reg <- if ("region_name" %in% names(di)) as.character(di$region_name[1]) else NA
  saveRDS(list(L = L, region = reg, year = yr$year, obs = yr$cat, model = mc,
               Q = Qf, q_pel = qpf, q_ben = qbf, two_Q = two_Q, mse = mse, corr = cr,
               s_pel = s_pel, s_ben = s_ben, depth = depth_mean),
          Sys.getenv("SAVE_RDS"))
}

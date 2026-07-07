# Tier-1 sandbox: single-Q fishing calibration (transient, one FAO-LME)
# ---------------------------------------------------------------------
# Estimate ONE catchability Q by minimising the difference between the MODELLED
# and OBSERVED catch TIME SERIES (log MSE), given the model + forcing:
#   F_g(x,t) = Q * s(x) * (B_g/sum_g B) * effort_norm(t)
# where
#   s(x,t)          min-max WINDOW selectivity, 1 for min_fished(t) <= w <= max_fished(t),
#                   read TIME-VARYING per year from the parquet's min_fished_weight_class /
#                   max_fished_weight_class (log10 g) columns (FSIZE=knife reverts to 10 g knife-edge),
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
# glob the init json regardless of the searchvol suffix (most LMEs are 12.8; some, e.g. Arabian
# Sea LME-32, are 6.4) so every region resolves without hardcoding the coefficient.
p  <- fromJSON(Sys.glob(file.path(base, "equilibrium_runs",
        sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_*.json", L)))[1])$params
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
dcorr <- min(p$depth[1], 200); bhd <- 20; wsel <- 1 * LN10   # 10 g knife-edge (fallback default)

# --- obsclim forcing: annual effort/catch series + monthly plankton intercept/slope ---
# INPUT_PARQUET_DIR overrides just the parquet directory (e.g. a copy augmented with the
# per-spectrum min/max_fished_U/V columns); base still supplies equilibrium_runs/*.json etc.
pq_dir <- Sys.getenv("INPUT_PARQUET_DIR", file.path(base, "dbpm_inputs"))
di <- read_parquet(Sys.glob(file.path(pq_dir,
        sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet", L)))[1]) |>
      filter(str_detect(scenario, "obsclim")) |> arrange(year, month)
yr <- di |> group_by(year) |>
      summarise(eff = mean(total_nom_active_area_m2, na.rm = TRUE),
                cat = mean(catch_tonnes_area_m2, na.rm = TRUE) * 1e6, .groups = "drop") |>
      arrange(year)
nyr <- nrow(yr); spin <- 40; tmax <- spin + nyr
# --- v2 TWO-GROUP: per-spectrum effort drivers (U/V) + observed catch split (see CALIB_v2_METHODS.md)
# effort_split_lme.csv / catch_split_lme.csv: LME,Year,class(U/V),{effort|catch}, region=LME or fao+100.
EFF_SPLIT   <- Sys.getenv("EFFORT_SPLIT_CSV", "~/dbpm_compare_scratch/effort_split_lme.csv")
CATCH_SPLIT <- Sys.getenv("CATCH_SPLIT_CSV",  "~/dbpm_compare_scratch/catch_split_lme.csv")
rd_split <- function(f, val) {   # -> matrix cols U,V aligned to yr$year (0 where absent)
  x <- read.csv(path.expand(f)); x <- x[x$LME == L, ]
  key <- paste(x$Year, x$class)
  gv <- function(cl) { v <- x[[val]][match(paste(yr$year, cl), key)]; v[!is.finite(v)] <- 0; v }
  cbind(U = gv("U"), V = gv("V"))
}
es <- rd_split(EFF_SPLIT, "effort"); cs <- rd_split(CATCH_SPLIT, "catch")
# two effort drivers, each normalised to its OWN maximum (q absorbs absolute level)
effU_n <- if (max(es[, "U"]) > 0) es[, "U"] / max(es[, "U"]) else es[, "U"] * 0
effV_n <- if (max(es[, "V"]) > 0) es[, "V"] / max(es[, "V"]) else es[, "V"] * 0
eidx    <- function(t) max(1, min(nyr, floor(t) - spin + 1))
effU_at <- function(t) effU_n[eidx(t)]
effV_at <- function(t) effV_n[eidx(t)]
# observed catch: fraction from the split x parquet total (units g m-2 yr-1); series sum to yr$cat
ctot <- cs[, "U"] + cs[, "V"]; pel_frac <- ifelse(ctot > 0, cs[, "U"] / ctot, NA_real_)
obs_pel <- yr$cat * pel_frac; obs_ben <- yr$cat * (1 - pel_frac)
# --- TIME-VARYING PER-SPECTRUM fished-size WINDOW (log10 g) from the parquet ------------
# U (pelagic) window from min/max_fished_U, V (benthic) from min/max_fished_V -- the columns
# script 04 derives from the FGroup catch (fish+krill+ceph = U; shrimp/lobster/mollusc = V).
# A single fish-derived window applied to both spectra guts the small-bodied benthic fishery,
# so U and V are kept separate. Degenerate (hi<=lo) -> open the top; NA U -> 10 g knife-edge;
# NA V -> spectrum not fished that year (empty window). Fallbacks in priority order:
#   FSIZE=knife            -> both spectra fixed 10 g..Inf (reproduce v1)
#   parquet min/max_fished_U/_V columns present -> use them (protocol)
#   FSIZE_UV_CSV=<file>    -> region,year,min_fished_U,max_fished_U,min_fished_V,max_fished_V
#   single min/max_fished_weight_class present -> both spectra = that one window
#   else                   -> knife-edge
fmode <- Sys.getenv("FSIZE", "parquet"); uvcsv <- Sys.getenv("FSIZE_UV_CSV", "")
mkwin <- function(lo, hi, na_open) {                 # align to yr$year, tidy degenerate/NA
  hi <- ifelse(is.finite(hi) & hi > lo, hi, ifelse(is.finite(lo), Inf, NA))
  bad <- !is.finite(lo)
  if (na_open) { lo[bad] <- 1; hi[bad] <- Inf }      # U: NA -> 10 g knife-edge
  else         { lo[bad] <- Inf; hi[bad] <- Inf }    # V: NA -> not fished (empty)
  list(lo = lo, hi = hi) }
have_uv <- fmode != "knife" &&
           all(c("min_fished_U","max_fished_U","min_fished_V","max_fished_V") %in% names(di))
if (fmode == "knife") {
  Uw <- list(lo = rep(1, nyr), hi = rep(Inf, nyr)); Vw <- Uw; src <- "knife 10g..Inf"
} else if (have_uv) {
  fs <- di |> group_by(year) |> summarise(ul=min_fished_U[1], uh=max_fished_U[1],
        vl=min_fished_V[1], vh=max_fished_V[1], .groups="drop") |> arrange(year)
  m <- match(yr$year, fs$year)
  Uw <- mkwin(fs$ul[m], fs$uh[m], TRUE); Vw <- mkwin(fs$vl[m], fs$vh[m], FALSE); src <- "parquet U/V"
} else if (nzchar(uvcsv) && file.exists(uvcsv)) {
  cs <- read.csv(uvcsv); cs <- cs[cs$region == L, ]; m <- match(yr$year, cs$year)
  Uw <- mkwin(cs$min_fished_U[m], cs$max_fished_U[m], TRUE)
  Vw <- mkwin(cs$min_fished_V[m], cs$max_fished_V[m], FALSE); src <- paste("UV csv", basename(uvcsv))
} else if (all(c("min_fished_weight_class","max_fished_weight_class") %in% names(di))) {
  fs <- di |> group_by(year) |> summarise(l=min_fished_weight_class[1],
        h=max_fished_weight_class[1], .groups="drop") |> arrange(year)
  m <- match(yr$year, fs$year); Uw <- mkwin(fs$l[m], fs$h[m], TRUE); Vw <- Uw
  src <- "single parquet min/max_fished_weight_class (both spectra)"
} else { Uw <- list(lo=rep(1,nyr), hi=rep(Inf,nyr)); Vw <- Uw; src <- "knife (no cols)" }
winU_at <- function(t) { i <- max(1, min(nyr, floor(t)-spin+1)); c(Uw$lo[i], Uw$hi[i]) }  # log10 g
winV_at <- function(t) { i <- max(1, min(nyr, floor(t)-spin+1)); c(Vw$lo[i], Vw$hi[i]) }
wU_ref <- winU_at(spin+nyr-1); wV_ref <- winV_at(spin+nyr-1)   # ref windows for the gravity split
cat(sprintf("FSIZE [%s] U %.2f..%.2f V %.2f..%.2f (last yr, log10 g)\n",
            src, wU_ref[1], wU_ref[2], wV_ref[1], wV_ref[2]))
# INT_OFFSET: additive log10 shift on the plankton intercept. It corrects the parquet's
# fixed-0..200 m mean phyto to a biomass-weighted vertical average (deep upwelling columns
# gain ~+0.5, broad shelves ~+0.2). Preference order for this LME:
#   1) DINT_CSV=<file> with columns lme,dint -> use the measured Delta for LME `L`
#   2) INT_OFFSET=<x>  scalar fallback (emulation), default 0 (parquet fixed-200 m as-is).
# The CSV values come from the exact FAO-LME-mask area-average of ISIMIP3a phyc/phypico
# (int_biomass_weighted - int_fixed200); see ~/dbpm_compare_scratch/lme_dint.csv.
int_off <- as.numeric(Sys.getenv("INT_OFFSET", "0")); slp_off <- 0
dint_csv <- Sys.getenv("DINT_CSV", "")
if (nzchar(dint_csv) && file.exists(dint_csv)) {
  dt <- read.csv(dint_csv); if ("lme" %in% names(dt) && L %in% dt$lme) {
    i <- match(L, dt$lme)[1]
    if (is.finite(dt$dint[i])) int_off <- dt$dint[i]
    # DSLOPE: MATCHED slope correction (average-then-fit) applied with the intercept so the
    # size spectrum 10^intercept * exp(slope*m) stays internally consistent (same averaging).
    if ("dslope" %in% names(dt) && is.finite(dt$dslope[i])) slp_off <- dt$dslope[i]
    cat(sprintf("DINT_CSV: LME %d intercept shift = %+0.2f  slope shift = %+0.3f\n", L, int_off, slp_off))
  }
}
intm <- di$intercept + int_off; slpm <- di$slope + slp_off; nmon <- length(intm)
int1 <- mean(di$intercept[di$year == min(di$year)]) + int_off
slp1 <- mean(di$slope[di$year == min(di$year)]) + slp_off
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
  on.exit({ setwd(old); unlink(wd, recursive = TRUE, force = TRUE) }, add = TRUE)  # clean per-run scratch
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
  # "export:GAMMA:ATTN0:MAGG" -> SIZE-DEPENDENT export: a surface particle of ln-mass m reaches the
  #   seafloor with fraction exp(-attn(t)*exp(-GAMMA*max(m, MAGG*ln10))) (large carcasses/pellets sink
  #   faster). attn(t) = ATTN0*(dcorr/200)*exp(0.05*(tob(t)-mean tob)) (deeper/warmer -> more attenuation).
  #   MAGG (log10 g, default -3) is the aggregation floor: production below it coagulates into marine
  #   snow of that effective sinking size (so tiny plankton export via aggregates, not as single cells).
  cl <- Sys.getenv("DET_CLOSURE", "burial")
  if (grepl("^residence", cl)) {
    tau <- as.numeric(sub("^residence:?", "", cl)); if (is.na(tau)) tau <- 1
    hdr <- "pel_tempeff,ben_tempeff,sinking_rate,residence_time"
    rows <- sprintf("%.8g,%.8g,%.8g,%.8g", tmat[1, ], tmat[2, ], smat, tau)
  } else if (grepl("^export", cl)) {
    prm  <- as.numeric(strsplit(sub("^export:?", "", cl), ":")[[1]])
    gam  <- if (length(prm) >= 1 && !is.na(prm[1])) prm[1] else 0.4
    att0 <- if (length(prm) >= 2 && !is.na(prm[2])) prm[2] else 1
    magg <- (if (length(prm) >= 3 && !is.na(prm[3])) prm[3] else -3) * LN10   # aggregation floor (ln-mass)
    tobt <- vapply(ts_t, function(t) tob[midx(t)], numeric(1)); tobt[ts_t < spin] <- tob1
    attn <- att0 * (dcorr / 200) * exp(0.05 * (tobt - mean(tob)))
    hdr  <- "pel_tempeff,ben_tempeff,sinking_rate,export_attn,export_gamma,export_magg"
    rows <- sprintf("%.8g,%.8g,%.8g,%.8g,%.8g,%.8g", tmat[1, ], tmat[2, ], smat, attn, gam, magg)
  } else {
    hdr <- "pel_tempeff,ben_tempeff,sinking_rate,depth"
    rows <- sprintf("%.8g,%.8g,%.8g,%.8g", tmat[1, ], tmat[2, ], smat, dcorr)
  }
  # DYNAMIC gravity: add the gravity channel so the C engine splits F by current biomass share
  if (nzchar(Sys.getenv("GRAVITY"))) { hdr <- paste0(hdr, ",gravity"); rows <- paste0(rows, ",1") }
  writeLines(c(hdr, rows), file.path(di_in, "forcing_ts.txt"))
  if (fishing) {
    selU <- function(m, t) { w <- winU_at(t); ml <- m / LN10; as.numeric(ml >= w[1] & ml <= w[2]) }
    selV <- function(m, t) { w <- winV_at(t); ml <- m / LN10; as.numeric(ml >= w[1] & ml <= w[2]) }
    Setup.fishing(pe, run, grid, func = function(m, t, x, y) qp * selU(m, t) * effU_at(t))
    Setup.fishing(be, run, grid, func = function(m, t, x, y) qb * selV(m, t) * effV_at(t))
  }
  ok <- tryCatch({ invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de))); TRUE },
                 error = function(e) { message(conditionMessage(e)); FALSE })
  if (!ok) { setwd(old); return(NULL) }
  f <- Read.In("R", "fish"); b <- Read.In("R", "benthos"); setwd(old)
  list(Um = as.matrix(f@uvals[, -(1:3)]), Vm = as.matrix(b@uvals[, -(1:3)]), tt = f@uvals[, 1],
       x = f@mrange / LN10, w = exp(f@mrange), dm = diff(f@mrange)[1])
}

# --- 1) unfished equilibrium (fishable-biomass shares reported for info only) ---
u0 <- runsim(0, 0); nT <- length(u0$tt)
fiU <- u0$x >= wU_ref[1] & u0$x <= wU_ref[2]
fiV <- u0$x >= wV_ref[1] & u0$x <= wV_ref[2]
Bpel <- sum(u0$Um[nT, fiU] * u0$w[fiU] * u0$dm) * dcorr
Bben <- sum(u0$Vm[nT, fiV] * u0$w[fiV] * u0$dm) * bhd
s_pel <- Bpel / (Bpel + Bben); s_ben <- Bben / (Bpel + Bben)
# v2: which groups are identifiable (>=3 non-zero observed years)?
fitP <- sum(is.finite(obs_pel) & obs_pel > 0) >= 3
fitB <- sum(is.finite(obs_ben) & obs_ben > 0) >= 3
cat(sprintf("=== FAO-LME %d: TWO-GROUP (q_pel,q_ben) + time-varying plankton ===\n", L))
cat(sprintf("  fit groups: pelagic=%s benthic=%s | obs benthic frac (mean)=%.2f\n",
            fitP, fitB, mean(1 - pel_frac, na.rm = TRUE)))

# --- TWO modelled catch series (pelagic, benthic); C_g = q_g * E_g(t) * B_g^fished(t) ---
catch_ts <- function(qp, qb) {
  r <- runsim(qp, qb); if (is.null(r)) return(list(pel = rep(NA, nyr), ben = rep(NA, nyr)))
  k <- match(spin + (1:nyr), r$tt); cp <- numeric(nyr); cb <- numeric(nyr)
  for (i in seq_len(nyr)) { kk <- k[i]
    wu <- winU_at(r$tt[kk]); fu <- r$x >= wu[1] & r$x <= wu[2]
    wv <- winV_at(r$tt[kk]); fv <- r$x >= wv[1] & r$x <= wv[2]
    bU <- sum(r$Um[kk, fu] * r$w[fu] * r$dm) * dcorr   # areal fishable pelagic biomass (U window)
    bV <- sum(r$Vm[kk, fv] * r$w[fv] * r$dm) * bhd      # areal fishable benthic biomass (V window)
    cp[i] <- qp * effU_at(r$tt[kk]) * bU
    cb[i] <- qb * effV_at(r$tt[kk]) * bV }
  list(pel = cp, ben = cb)
}
# catch-volume-weighted log-MSE for one group (weight each year by that group's observed catch).
# Evaluated over observed>0 years; the model is FLOORED (not dropped) so a q->0 zero-catch
# prediction is penalised as a large under-shoot rather than escaping via NA.
wlogmse <- function(m, obs) {
  ok <- is.finite(obs) & obs > 0 & is.finite(m)
  if (sum(ok) < 3) return(NA_real_)
  mfl <- pmax(m[ok], 1e-9); w <- obs[ok] / sum(obs[ok])
  sum(w * (log10(mfl) - log10(obs[ok]))^2) }
wcor <- function(m, obs) { ok <- is.finite(m) & m > 0 & is.finite(obs) & obs > 0
  if (sum(ok) < 3) return(NA_real_); suppressWarnings(cor(log10(m[ok]), log10(obs[ok]))) }

# --- 2) estimate (q_pel, q_ben): joint two-series objective, grid seed + BOBYQA ---
qmax <- as.numeric(Sys.getenv("QMAX", "4")); n <- 0
obj2 <- function(q) { n <<- n + 1; r <- catch_ts(q[1], q[2])
  parts <- c(if (fitP) wlogmse(r$pel, obs_pel) else NA,
             if (fitB) wlogmse(r$ben, obs_ben) else NA)
  parts <- parts[is.finite(parts)]; if (length(parts) == 0) 1e6 else sum(parts) }
# optimise on LOG10(q) so q can span orders of magnitude (high-productivity regions need q~1e-4)
qmin <- as.numeric(Sys.getenv("QMIN", "1e-4"))
obj2L <- function(lq) obj2(10^lq)
lg <- seq(log10(qmin), log10(qmax), length.out = 7); GL <- expand.grid(lqp = lg, lqb = lg)
seedL <- as.numeric(GL[which.min(apply(GL, 1, function(q) obj2L(as.numeric(q)))), ])
res   <- nloptr::nloptr(seedL, obj2L, lb = c(log10(qmin), log10(qmin)),
          ub = c(log10(qmax), log10(qmax)),
          opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 80, xtol_rel = 1e-3))
qpf <- if (fitP) 10^res$solution[1] else NA_real_
qbf <- if (fitB) 10^res$solution[2] else NA_real_
J   <- res$objective

# --- per-group fit + report ---
mfit  <- catch_ts(if (is.na(qpf)) 0 else qpf, if (is.na(qbf)) 0 else qbf)
if (is.na(qpf)) mfit$pel[] <- NA
if (is.na(qbf)) mfit$ben[] <- NA
mse_p <- wlogmse(mfit$pel, obs_pel); mse_b <- wlogmse(mfit$ben, obs_ben)
cor_p <- wcor(mfit$pel, obs_pel);   cor_b <- wcor(mfit$ben, obs_ben)
cat(sprintf("  evals=%d  q_pel=%s  q_ben=%s  J=%.3f\n", n,
            if (is.na(qpf)) "NA" else sprintf("%.3f", qpf),
            if (is.na(qbf)) "NA" else sprintf("%.3f", qbf), J))
cat(sprintf("  PELAGIC  wMSE=%s corr=%s  model=%.3g obs=%.3g\n",
            if (is.na(mse_p)) "NA" else sprintf("%.3f", mse_p),
            if (is.na(cor_p)) "NA" else sprintf("%.2f", cor_p),
            mean(mfit$pel, na.rm = TRUE), mean(obs_pel, na.rm = TRUE)))
cat(sprintf("  BENTHIC  wMSE=%s corr=%s  model=%.3g obs=%.3g\n",
            if (is.na(mse_b)) "NA" else sprintf("%.3f", mse_b),
            if (is.na(cor_b)) "NA" else sprintf("%.2f", cor_b),
            mean(mfit$ben, na.rm = TRUE), mean(obs_ben, na.rm = TRUE)))

# --- save (per-group series + stats; total for backward-compatible plotting) ---
if (nzchar(Sys.getenv("SAVE_RDS"))) {
  reg <- if ("region_name" %in% names(di)) as.character(di$region_name[1]) else NA
  mtot <- rowSums(cbind(ifelse(is.na(mfit$pel), 0, mfit$pel),
                        ifelse(is.na(mfit$ben), 0, mfit$ben)), na.rm = TRUE)
  saveRDS(list(L = L, region = reg, year = yr$year,
               obs = yr$cat, model = mtot, obs_pel = obs_pel, obs_ben = obs_ben,
               model_pel = mfit$pel, model_ben = mfit$ben,
               q_pel = qpf, q_ben = qbf, two_Q = TRUE, J = J,
               mse_pel = mse_p, mse_ben = mse_b, corr_pel = cor_p, corr_ben = cor_b,
               mse = J, corr = wcor(mtot, yr$cat),
               fitP = fitP, fitB = fitB, s_pel = s_pel, s_ben = s_ben, depth = depth_mean),
          Sys.getenv("SAVE_RDS"))
}

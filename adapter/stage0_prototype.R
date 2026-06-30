# Stage 0 / Stage 2 adapter prototype
# ------------------------------------
# Reads one LME's ISIMIP3a forcing (parquet) + reference equilibrium (JSON),
# configures a matched no-fishing dbpmr 0-D run via the parameter crosswalk
# (design/stage0-reference.md), runs it, and compares the equilibrium pelagic
# spectrum to the LME `sizemodel()` reference.
#
# STATUS: read -> map -> run -> compare pipeline works end-to-end.
#   * search rate matched to sizemodel: A_pel = hr_volume_search (12.8),
#     A_ben = 0.1 * A_pel; temperature applied as in sizemodel() (feeding AND
#     background mortality * Boltzmann-Arrhenius; pelagic=surface, benthic=floor).
#   * plankton held fixed at the equilibrium input, confirmed constant.
#   * dbpmr is numerically UNSTABLE at monthly steps for these params (NaN) even
#     at A=12.8, but STABLE at weekly/daily; sizemodel() is stable at monthly --
#     a real discretisation difference despite the shared upwind scheme (#7).
# RESULT (issue #8): with A_pel=64, A_ben=6.4, weekly steps, temperature applied
#   (feeding & background mortality), and consumer-min = -3 (the LME convention),
#   dbpmr reaches a HEALTHY COEXISTING equilibrium -- BOTH pelagic (biomass ~0.4)
#   and benthic (~1.7) thriving -- which JB confirms is the correct LME-10
#   behaviour. The canonical sizemodel(), by contrast, collapses the pelagic for
#   LME-10. dbpmr's rep_method=2 (R~0.2 of intake) sustains both groups; the
#   canonical residual-energy reproduction (0.14 pel / 0.05 ben) is leaner and
#   collapses the pelagic.
#   Caveats: dbpmr is numerically fragile -- it NaNs at monthly steps, at
#   consumer-min=-7, and when reproduction is pushed; consumer-min=-3 + weekly is
#   the stable, correct operating point.
#
# REQUIRES: local LME data (not in this repo) + an installed dbpmr.
#   - arrow, jsonlite
#   - set DBPM_DATA to the folder holding dbpm_inputs/ and equilibrium_runs/

suppressMessages({library(arrow); library(jsonlite); library(dbpmr)})

base <- Sys.getenv("DBPM_DATA",
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev")
LME  <- 10
LN10 <- log(10)

## ---- 1. read forcing (stable-spin constants) + reference -------------------
fc <- as.data.frame(read_parquet(file.path(base, "dbpm_inputs",
        sprintf("dbpm_clim-fish-inputs_fao_lme-%d_1641-2010.parquet", LME))))
ss  <- fc[fc$scenario == "stable-spin", ][1, ]
ref <- fromJSON(file.path(base, "equilibrium_runs",
        sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_12.8.json", LME)))
p <- ref$params

## ---- 2. crosswalk LME params -> dbpmr ---------------------------------------
# Search rate matched to sizemodel: A_pel = hr_volume_search, A_ben = 0.1*A_pel.
# Temperature: sizemodel() scales BOTH the feeding rate AND the background
# ("other") mortality by the Boltzmann-Arrhenius factor (pelagic = surface temp,
# benthic = floor temp); senescence and predation mortality are NOT scaled.
# Forcing is constant in the spin-up, so we fold the constant factor into the
# effective A and mu_0 (a true per-timestep forcing is the Stage 3 engine work).
tempeff <- function(temp) exp(p$c1[1] - p$activation_energy[1] /
                                (p$boltzmann[1] * (temp + 273)))
pel_te  <- tempeff(p$sea_surf_temp[1])
ben_te  <- tempeff(p$sea_floor_temp[1])

A_pel   <- 64 * pel_te          # canonical search volume A.u = 64 (per JB / CMIP5)
A_ben   <- 6.4 * ben_te         # A.v = 0.1 * A.u = 6.4
mu0_pel <- p$natural_mort[1] * pel_te
mu0_ben <- p$natural_mort[1] * ben_te
# Energy budget (#22): each unit of intake -> defecation + assimilated (1-def);
# the assimilated fraction splits into growth (K), reproduction (R) and
# excretion (Ex), so K + R + Ex = 1 - def. sizemodel keys the budget to the PREY
# type: pelagic/plankton prey use def.high + growth_pred (K.u) + energy_pred
# (AM.u); benthic/detritus prey use def.low + growth_detritivore (K.v) +
# energy_detritivore (AM.v). NB dbpmr's `epsilon` is the senescence size offset,
# NOT a growth efficiency, so it must not be set from the growth coefficients.
defh <- p$defecate_prop[1]; defl <- p$def_low[1]
Ku <- p$growth_pred[1];          AMu <- p$energy_pred[1]
Kv <- p$growth_detritivore[1];   AMv <- p$energy_detritivore[1]
# pelagic predator on plankton/pelagic prey (def.high budget):
K_pel  <- (1-defh)*Ku ; R_pel  <- (1-defh)*(1-(Ku+AMu)) ; Ex_pel <- (1-defh)*AMu
# pelagic predator on benthic prey, and detritivore on detritus (def.low budget):
K_lo   <- (1-defl)*Kv ; R_lo   <- (1-defl)*(1-(Kv+AMv)) ; Ex_lo  <- (1-defl)*AMv
u0_pla  <- 10^p$int_phy_zoo[1]   # plankton held fixed at the equilibrium input
lam_pla <- p$slope_phy_zoo[1]

## ---- 3. configure & run dbpmr (no fishing, aspatial, to equilibrium) -------
wd <- tempfile("stage0"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
run  <- Setup.Run("LMErun", 1, 1, spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
grid <- Setup.Grid(run, mmin = -12*LN10, mmax = 6*LN10, mstep = 0.1*LN10,
          moutstep = 0.1*LN10, tmax = p$n_years[1], tstep = 1/52, toutstep = 1)
plankton <- Setup.Plankton(run, filename = "plankton",
              mmin = -12*LN10, mmax = -3*LN10, u_0 = u0_pla, lambda = lam_pla)
pelagic  <- Setup.Pelagic(run, filename = "fish", mmin = -3*LN10, mmax = 6*LN10,
              alpha = p$metabolic_req_pred[1], A = A_pel, mu_0 = mu0_pel,
              K_pla = K_pel, R_pla = R_pel, Ex_pla = Ex_pel,
              K_pel = K_pel, R_pel = R_pel, Ex_pel = Ex_pel,
              K_ben = K_lo,  R_ben = R_lo,  Ex_ben = Ex_lo,
              rep_method = 2, fishing_flag = FALSE)
benthic  <- Setup.Benthic(run, filename = "benthos", mmin = -3*LN10, mmax = 4*LN10,
              alpha = p$metabolic_req_detritivore[1], A = A_ben, mu_0 = mu0_ben,
              K_det = K_lo, R_det = R_lo, Ex_det = Ex_lo,
              rep_method = 2, fishing_flag = FALSE)
detritus <- Setup.Detritus(run, filename = "detritus")
invisible(capture.output(SizeSpectrum(run, grid, plankton, pelagic, benthic, detritus)))

## ---- 4. read dbpmr equilibrium + reference, compare ------------------------
fish <- Read.In("LMErun", "fish")
fu   <- as.numeric(fish@finaluvals[1, -(1:3)])
m_ln <- fish@grid@mmin + (seq_along(fu) - 1) * fish@grid@mstep
dbpmr_log10w <- m_ln / LN10
ref_pred     <- ref$predators[, ncol(ref$predators)]
ref_log10w   <- p$log10_size_bins

cat(sprintf("dbpmr nonzero bins: %d / %d ; reference nonzero: %d\n",
            sum(fu > 0), length(fu), sum(ref_pred > 0)))

png(file.path(old, "stage0_compare.png"), width = 820, height = 560)
plot(ref_log10w, log10(ref_pred), type = "l", lwd = 2, col = "black",
     xlab = "log10 body mass (g)", ylab = "log10 density",
     main = sprintf("Stage 0: dbpmr vs LME reference (no fishing) - LME %d", LME))
lines(dbpmr_log10w, log10(fu), lwd = 2, col = "firebrick")
legend("topright", c("LME reference (sizemodel)", "dbpmr"),
       col = c("black", "firebrick"), lwd = 2, bty = "n")
invisible(dev.off())

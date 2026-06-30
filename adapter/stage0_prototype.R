# Stage 0 / Stage 2 adapter prototype
# ------------------------------------
# Reads one LME's ISIMIP3a forcing (parquet) + reference equilibrium (JSON),
# configures a matched no-fishing dbpmr 0-D run via the parameter crosswalk
# (design/stage0-reference.md), runs it, and compares the equilibrium pelagic
# spectrum to the LME `sizemodel()` reference.
#
# STATUS (first run): the read -> map -> run -> compare pipeline works
# end-to-end, but the literal crosswalk does NOT reproduce the LME equilibrium
# -- dbpmr's pelagic spectrum decays toward extinction across the search-rate
# range and is numerically unstable when A is pushed high. The divergence is
# therefore deeper than a single search-volume scale factor (resource coupling /
# reproduction boundary / assimilation balance need reconciling). See issue #8.
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

## ---- 2. crosswalk LME params -> dbpmr (constant forcing; temp folds into A) -
hr <- p$hr_volume_search[1]
tempeff <- function(temp) exp(p$c1[1] - p$activation_energy[1] /
                                (p$boltzmann[1] * (temp + 273)))
A_pel   <- hr * p$pref_pelagic[1] * tempeff(p$sea_surf_temp[1])
A_ben   <- hr * p$pref_benthos[1] * tempeff(p$sea_floor_temp[1])
eps_pel <- (1 - p$defecate_prop[1]) * p$growth_pred[1]
eps_ben <- (1 - p$defecate_prop[1]) * p$growth_detritivore[1]
u0_pla  <- 10^p$int_phy_zoo[1]
lam_pla <- p$slope_phy_zoo[1]

## ---- 3. configure & run dbpmr (no fishing, aspatial, to equilibrium) -------
wd <- tempfile("stage0"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
run  <- Setup.Run("LMErun", 1, 1, spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
grid <- Setup.Grid(run, mmin = -12*LN10, mmax = 6*LN10, mstep = 0.1*LN10,
          moutstep = 0.1*LN10, tmax = p$n_years[1], tstep = 1/12, toutstep = 1)
plankton <- Setup.Plankton(run, filename = "plankton",
              mmin = -12*LN10, mmax = -3*LN10, u_0 = u0_pla, lambda = lam_pla)
pelagic  <- Setup.Pelagic(run, filename = "fish", mmin = -3*LN10, mmax = 6*LN10,
              alpha = p$metabolic_req_pred[1], A = A_pel, epsilon = eps_pel,
              mu_0 = p$natural_mort[1], rep_method = 2, fishing_flag = FALSE)
benthic  <- Setup.Benthic(run, filename = "benthos", mmin = -3*LN10, mmax = 4*LN10,
              alpha = p$metabolic_req_detritivore[1], A = A_ben, epsilon = eps_ben,
              mu_0 = p$natural_mort[1], rep_method = 2, fishing_flag = FALSE)
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

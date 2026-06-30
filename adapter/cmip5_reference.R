# Canonical CMIP5 sizemodel() reference driver
# ---------------------------------------------
# Runs the canonical DBPM `sizemodel()` (from the dbpm_isimip_3b repo) with one
# LME's constant stable-spin forcing, to produce a *trustworthy* Stage 0
# reference (the LME calibration repo's equilibrium JSONs are a different,
# possibly-problematic parameterisation). Compare its output against
# adapter/stage0_prototype.R (dbpmr) at the same A.
#
# Get the model file (not vendored here -- it lives in another repo):
#   curl -fsSL https://raw.githubusercontent.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_isimip_3b/master/size-based-models/dynamic_sizebased_model_functions_CMIP52019.R -o cmip5_sizemodel.R
#
# Requires: jsonlite; the LME data folder (DBPM_DATA).
#
# FINDING: at the canonical A.u = 64, this sizemodel() *also* collapses the
# LME-10 pelagic predators (max ~1e-31) while detritivores stay alive -- i.e.
# the collapsed pelagic is the sizemodel's genuine behaviour for this deep,
# cold-floored, low-export LME, reproduced by the canonical code (not a
# calibration artifact). dbpmr at the same config does the OPPOSITE: it sustains
# the pelagic and (with the cold-floor temperature) collapses the benthos. The
# two engines partition the benthic-pelagic balance very differently -- a
# coupling/reproduction reconciliation, not parameter tuning. See issue #8.

suppressMessages(library(jsonlite))
source("cmip5_sizemodel.R")   # the downloaded canonical model

base <- Sys.getenv("DBPM_DATA",
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev")
LME <- 10
p <- fromJSON(file.path(base, "equilibrium_runs",
        sprintf("init_dbpm_nonspatial_fao_lme-%d_searchvol_12.8.json", LME)))$params

# LME constant stable-spin forcing as length-(Neq+2) vectors (sizemodel indexes
# pp[i], sst[i], ... each step). NOTE: the export ratio lives in `sinking_rate`.
tmax <- 150; tspy <- 48; Neq <- tmax * tspy
cnst <- function(v) rep(v, Neq + 2)

param <- sizeparam(equilibrium = FALSE, dx = 0.1, xmin = -12, xmax = 6,
            xmin.consumer.u = -7, xmin.consumer.v = -7, tmax = tmax, tstepspryr = tspy,
            er = cnst(p$sinking_rate[1]), pp = cnst(p$int_phy_zoo[1]),
            slope = cnst(p$slope_phy_zoo[1]), depth = p$depth[1],
            sst = cnst(p$sea_surf_temp[1]), sft = cnst(p$sea_floor_temp[1]))
# param$A.u defaults to 64; param$A.v = 0.1 * A.u. Override here to compare at
# other search volumes if desired.

res <- sizemodel(param, U_mat = NA, V_mat = NA, W_mat = NA, temp.effect = TRUE)

x  <- param$x
Uf <- res$U[, ncol(res$U)]   # equilibrium predator spectrum
Vf <- res$V[, ncol(res$V)]   # equilibrium detritivore spectrum
cat(sprintf("CMIP5 sizemodel @A.u=%g: predator max(x>-3)=%.3g (%s); detritivore max=%.3g\n",
            param$A.u, max(Uf[x > -3]),
            ifelse(max(Uf[x > -3]) < 1e-10, "COLLAPSED", "ALIVE"), max(Vf[x > -3])))
saveRDS(list(x = x, U = Uf, V = Vf, A = param$A.u), "cmip5_ref_lme10.rds")

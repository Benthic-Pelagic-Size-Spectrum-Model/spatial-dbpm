# dbpmr all-LME sweep driver
# --------------------------
# Runs dbpmr across every LME equilibrium input at the canonical setup matched to
# the reference sizemodel() (A.u = 64, A.v = 0.1*A.u = 6.4, consumer minimum
# 10^-3 g, weekly step, temperature on feeding + background mortality, the
# faithful K/R/Ex energy budget of design/sizemodel-investigation.md S8), and
# tallies how many LMEs stay alive vs collapse vs go non-finite.
#
# FINDING (see design/sizemodel-investigation.md S6a): dbpmr completes ALL 82
# LMEs with no crash and no NaN, and sustains pelagic + benthos in EVERY one --
# including the 12 warm + oligotrophic basins that sizemodel() collapses
# (10, 12, 16, 30, 31, 40, 44, 131, 134, 151, 171, 177) and the most oligotrophic
# LME-64. sizemodel() at the identical setup gives 70 alive / 12 collapsed. The
# result is insensitive to A.v (6.4 or 0.64 -> all alive). This is the
# high-throughput, predation-dominated robustness of S4d: dbpmr does not carry
# sizemodel's A-independent background-mortality collapse.
#
# Requires: an installed dbpmr, plus jsonlite; the LME equilibrium_runs JSONs.
# Point DBPM_DATA at the folder holding equilibrium_runs/:
#   DBPM_DATA=/path/to/DBPM_dev Rscript adapter/lme_sweep.R
# Optional env: AVRATIO (A.v / A.u, default 0.1), TMAX (spin years, default 100),
#   OUT (results CSV path, default lme_sweep_results.csv).

suppressMessages({library(jsonlite); library(dbpmr)})

base <- Sys.getenv("DBPM_DATA",
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev")
eqdir   <- file.path(base, "equilibrium_runs")
Avratio <- as.numeric(Sys.getenv("AVRATIO", "0.1"))   # A.v = Avratio * A.u
tmax    <- as.numeric(Sys.getenv("TMAX", "100"))
outcsv  <- Sys.getenv("OUT", "lme_sweep_results.csv")
LN10    <- log(10)

files  <- list.files(eqdir, pattern = "init_dbpm_nonspatial_fao_lme-\\d+_searchvol_12.8.json")
if (!length(files)) stop("No LME JSONs found under ", eqdir, " (set DBPM_DATA)")
lme_of <- function(f) as.integer(sub(".*lme-(\\d+)_.*", "\\1", f))
files  <- files[order(sapply(files, lme_of))]

# Boltzmann-Arrhenius temperature factor (constant in the stable spin-up, so we
# fold it into the search constant A and the background mortality mu_0).
tempeff <- function(p, temp) exp(p$c1[1] - p$activation_energy[1] /
                                   (p$boltzmann[1] * (temp + 273)))

run_one <- function(p) {
  pt <- tempeff(p, p$sea_surf_temp[1]); bt <- tempeff(p, p$sea_floor_temp[1])
  # Faithful energy budget (K growth / R reproduction / Ex excretion), keyed to
  # prey type: pelagic/plankton use def.high, benthic/detritus use def.low.
  defh <- p$defecate_prop[1]; defl <- p$def_low[1]
  Ku <- p$growth_pred[1];        AMu <- p$energy_pred[1]
  Kv <- p$growth_detritivore[1]; AMv <- p$energy_detritivore[1]
  Kp <- (1-defh)*Ku; Rp <- (1-defh)*(1-(Ku+AMu)); Exp <- (1-defh)*AMu
  Kl <- (1-defl)*Kv; Rl <- (1-defl)*(1-(Kv+AMv)); Exl <- (1-defl)*AMv

  wd <- tempfile("lme"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  run  <- Setup.Run("LMErun", 1, 1, spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
  # Mass grid + consumer/plankton ranges default to the canonical FishMIP log10
  # grid (10^-12..10^6 g, consumer min 10^-3 g); only the time grid is set here.
  grid <- Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = 1)
  pl <- Setup.Plankton(run, filename = "plankton",
          u_0 = 10^p$int_phy_zoo[1], lambda = p$slope_phy_zoo[1])
  pe <- Setup.Pelagic(run, filename = "fish", alpha = p$metabolic_req_pred[1],
          A = 64 * pt, mu_0 = p$natural_mort[1] * pt,
          K_pla = Kp, R_pla = Rp, Ex_pla = Exp,
          K_pel = Kp, R_pel = Rp, Ex_pel = Exp,
          K_ben = Kl, R_ben = Rl, Ex_ben = Exl, rep_method = 2)
  be <- Setup.Benthic(run, filename = "benthos", alpha = p$metabolic_req_detritivore[1],
          A = 64 * Avratio * bt, mu_0 = p$natural_mort[1] * bt,
          K_det = Kl, R_det = Rl, Ex_det = Exl, rep_method = 2)
  de <- Setup.Detritus(run, filename = "detritus")
  invisible(capture.output(SizeSpectrum(run, grid, pl, pe, be, de)))

  fsh <- Read.In("LMErun", "fish"); bn <- Read.In("LMErun", "benthos")
  fu <- as.numeric(fsh@finaluvals[1, -(1:3)]); bu <- as.numeric(bn@finaluvals[1, -(1:3)])
  fx <- (fsh@grid@mmin + (seq_along(fu) - 1) * fsh@grid@mstep) / LN10
  finite <- all(is.finite(c(fu, bu)))
  pm <- suppressWarnings(max(fu[fx > -3])); bm <- suppressWarnings(max(bu))
  status <- if (!finite) "NaN" else if (pm < 1e-10) "COLLAPSED" else "ALIVE"
  list(pred = pm, det = bm, status = status)
}

res <- data.frame()
for (f in files) {
  p <- fromJSON(file.path(eqdir, f))$params
  out <- tryCatch(run_one(p),
                  error = function(e) list(pred = NA, det = NA, status = "ERROR"))
  res <- rbind(res, data.frame(
    lme = lme_of(f), sst = round(p$sea_surf_temp[1], 1),
    pp = round(p$int_phy_zoo[1], 2), pred = signif(out$pred, 3),
    det = signif(out$det, 3), status = out$status))
}

cat(sprintf("dbpmr LME sweep: A.u=64, A.v=%.2f, consumer-min -3, weekly(1/48), %g-yr spin\n\n",
            64 * Avratio, tmax))
print(res, row.names = FALSE)
cat("\n=== tally ===\n"); print(table(res$status))
write.csv(res, outcsv, row.names = FALSE)
cat("\nwrote", outcsv, "\n")

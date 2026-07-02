# Tier-1 sandbox: engine-swap on one FAO-LME (aspatial, stable-spin, no fishing)
# -----------------------------------------------------------------------------
# Demonstrates dbpmr slotting into the LME workflow as an alternative engine.
# The pipeline (sizeparam, forcing loaders, run_model) is IMPORTED from the
# lme-workflow submodule; only the integrator is switched via engine = ...
#
#   DBPM_DATA=/path/to/DBPM_dev Rscript adapter/sandbox/tier1_engine_swap.R [LME]
#
# Requires: installed dbpmr (.libPaths below), arrow + the submodule's deps.

.libPaths(c("/tmp/dbpmrlib", .libPaths()))
SB   <- "/Users/juliab6/spatial-dbpm/adapter/sandbox"
base <- Sys.getenv("DBPM_DATA",
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev")
L <- as.integer(c(commandArgs(TRUE), "14")[1])

suppressMessages({ library(dbpmr); library(arrow); library(dplyr); library(stringr) })
# import the LME pipeline (sizeparam, run_model, sizemodel, ...) - not copied
source(file.path(SB, "lme-workflow", "scripts", "useful_functions.R"))
source(file.path(SB, "dbpmr_engine.R"))

pq <- Sys.glob(file.path(base, "dbpm_inputs",
        sprintf("dbpm_clim-fish-inputs_fao_lme-%d_*.parquet", L)))[1]
dbpm_inputs <- read_parquet(pq) |> filter(str_detect(scenario, "stable"))
fishing_params <- data.frame(region = sprintf("lme-%d", L),
        fmort_u = 0, fmort_v = 0, fminx_u = 0, fminx_v = 0, search_vol = 64)

# --- engine = "dbpmr" | "sizemodel", same call, same forcing/params ---
run_engine <- function(engine) {
  if (engine == "dbpmr") {
    r <- run_model_dbpmr(fishing_params, dbpm_inputs, xmin_consumer_u = -3, xmin_consumer_v = -3)
    c(pred = r$total_pred_biomass, det = r$total_detritivore_biomass)
  } else {
    r <- run_model(fishing_params, dbpm_inputs, withinput = FALSE,
                   xmin_consumer_u = -3, xmin_consumer_v = -3)
    p <- r$params; x <- p$log10_size_bins; dx <- p$log_size_increase
    dcorr <- min(p$depth, 200); cons <- x > -3; nc <- ncol(r$predators)
    c(pred = sum(r$predators[cons, nc] * 10^x[cons] * dx) * dcorr,
      det  = sum(r$detritivores[cons, nc] * 10^x[cons] * dx) * 20)
  }
}

cat(sprintf("=== Tier-1 engine swap, FAO-LME %d (A.u=64, no fishing, stable-spin) ===\n", L))
for (eng in c("sizemodel", "dbpmr")) {
  b <- tryCatch(run_engine(eng), error = function(e) { cat("  ", eng, "ERROR:", conditionMessage(e), "\n"); c(pred=NA, det=NA) })
  cat(sprintf("  engine=%-10s  total predator biomass=%.4g  detritivore=%.4g  (g m-2)\n",
              eng, b["pred"], b["det"]))
}
cat("\nNote: magnitudes are not yet reconciled - the submodule sizemodel is the\n",
    "*pre-fix* new_features code (carries the ln(10) bug), and the ln-vs-log10\n",
    "density convention + temperature folding differ. Tier-1's purpose is to prove\n",
    "the engine= swap RUNS through the imported pipeline; matching biomass (apply\n",
    "PR #14 to the submodule + reconcile the density convention) is the next step.\n", sep="")

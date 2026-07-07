# augment_parquets_uv.R -- add the per-spectrum fished-size columns (min/max_fished_U,
# min/max_fished_V) to the LOCAL DBPM_dev parquets, so the LME calibration reads them
# straight from the parquet (the have_uv path in tier1) -- identical to what script 04
# will write on Gadi. ADDITIVE ONLY: existing columns are untouched; re-runnable.
#
# Source of the window values: fished_size_UV_tv.csv (built by build_fished_size_uv.R,
# mirrors the script-04 uv_summ logic). Joined per (region, year); constant within year
# like the existing min/max_fished_weight_class.
suppressMessages({library(arrow); library(dplyr)})

dir_in <- Sys.getenv("PARQUET_DIR",
  "/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/DBPM_mizer/DBPM_dev/dbpm_inputs")
uv <- read.csv(Sys.getenv("UV_CSV",
  "/Users/juliab6/spatial-dbpm/adapter/sandbox/fished_size_UV_tv.csv"))  # region,year,min/max_fished_U/V
dir_out <- Sys.getenv("OUT_DIR", dir_in)                          # default: in place
dry <- nzchar(Sys.getenv("DRY_RUN", ""))

files <- list.files(dir_in, pattern="dbpm_clim-fish-inputs_fao_lme-\\d+_.*\\.parquet$", full.names=TRUE)
cat(sprintf("parquets: %d  |  writing to: %s%s\n", length(files), dir_out, if(dry) " (DRY RUN)" else ""))
done <- 0
for (f in files) {
  L <- as.integer(sub(".*fao_lme-(\\d+)_.*", "\\1", basename(f)))
  d <- read_parquet(f)
  u <- uv[uv$region == L, c("year","min_fished_U","max_fished_U","min_fished_V","max_fished_V")]
  if (nrow(u) == 0) { cat(sprintf("  LME %-4d no UV rows -- skipped\n", L)); next }
  # drop any prior UV cols (idempotent), then join by year
  d <- d[, setdiff(names(d), c("min_fished_U","max_fished_U","min_fished_V","max_fished_V"))]
  d <- left_join(d, u, by = "year")
  # HYBRID max: pelagic U max <- the parquet's real (Reg WtMax-based) max_fished_weight_class;
  # benthic V max stays the invert-group max (parquet max is fish-dominated -> too large for benthos).
  # min (U and V) stays FGroup-derived. See build note.
  if ("max_fished_weight_class" %in% names(d)) d$max_fished_U <- d$max_fished_weight_class
  if (!dry) write_parquet(d, file.path(dir_out, basename(f)))
  done <- done + 1
  if (L %in% c(1,3,148)) cat(sprintf("  LME %-4d yrs=%d  U[%.2f,%.2f] V[%.2f,%.2f] (last)\n", L, nrow(u),
      tail(u$min_fished_U,1), tail(u$max_fished_U,1), tail(u$min_fished_V,1), tail(u$max_fished_V,1)))
}
cat(sprintf("augmented %d parquets\n", done))

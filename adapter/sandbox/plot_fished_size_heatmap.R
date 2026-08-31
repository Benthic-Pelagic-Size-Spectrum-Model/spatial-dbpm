# plot_fished_size_heatmap.R -- supplement figures: MIN and MAX fished size through time by
# LME/FAO region, pelagic (U) vs benthic (V). Reads the ACTUAL USED values from the augmented
# parquets (min/max_fished_U/V) -- U max is the hybrid (Reg WtMax), V max the invert-group max.
# Outputs vector PDF + PNG for each of {min, max}.
#   Rscript adapter/sandbox/plot_fished_size_heatmap.R
# env: PARQUET_DIR (augmented parquets), OUTDIR (default adapter/sandbox/figs)
suppressMessages({library(ggplot2); library(tidyr); library(arrow)})

pdir   <- Sys.getenv("PARQUET_DIR", "~/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs")
outdir <- Sys.getenv("OUTDIR", "adapter/sandbox/figs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

files <- Sys.glob(file.path(pdir, "dbpm_clim-fish-inputs_fao_lme-*_*.parquet"))
rows <- do.call(rbind, lapply(files, function(f) {
  L <- as.integer(sub(".*fao_lme-([0-9]+)_.*", "\\1", basename(f)))
  d <- read_parquet(f, col_select = c("year", "region_name",
         "min_fished_U", "max_fished_U", "min_fished_V", "max_fished_V"))
  d <- d[!duplicated(d$year), ]
  pre <- if (L < 100) paste0("LME", L) else paste0("FAO", L - 100)
  d$region <- L; d$lab <- sprintf("%-6s %s", pre, substr(d$region_name[1], 1, 20)); d }))
rows <- as.data.frame(rows)
rows <- rows[rows$year >= 1950 & rows$year <= 2010, ]
ord  <- unique(rows$lab[order(-rows$region)])
rows$lab <- factor(rows$lab, levels = ord)

mkfig <- function(bound) {                                   # bound = "min" or "max"
  L <- pivot_longer(rows, paste0(bound, c("_fished_U", "_fished_V")),
                    names_to = "spectrum", values_to = "v")
  L$spectrum <- factor(ifelse(grepl("_U$", L$spectrum), "U  (pelagic)", "V  (benthic)"))
  brk <- if (bound == "min") c(0, 1, 2, 2.43, 3.86) else c(2, 3, 4, 5, 6)
  lab <- if (bound == "min") c("1 g","10 g","100 g","270 g","7 kg") else c("100 g","1 kg","10 kg","100 kg","1 t")
  ttl <- sprintf("%s fished size through time by region -- pelagic (U) vs benthic (V)",
                 ifelse(bound == "min", "Minimum", "Maximum"))
  sub <- if (bound == "min")
    "Per region-year from catch composition; grey = no fishery (no reported catch). Colour = log10 g."
  else
    "U max = Reg's real WtMax (hybrid); V max = benthic invert-group max. grey = no fishery. Colour = log10 g."
  p <- ggplot(L, aes(year, lab, fill = v)) + geom_tile() + facet_wrap(~spectrum) +
    scale_fill_viridis_c(name = paste(bound, "fished\nsize"), option = "magma", na.value = "grey85",
      breaks = brk, labels = lab) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1950, 1970, 1990, 2010)) +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(axis.text.y = element_text(size = 5.2), panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 9),
          plot.title = element_text(face = "bold", size = 11))
  base <- file.path(outdir, sprintf("fished_size_UV_%s_heatmap", bound))
  ggsave(paste0(base, ".pdf"), p, width = 11, height = 13)
  ggsave(paste0(base, ".png"), p, width = 11, height = 13, dpi = 150)
  cat("wrote", base, ".pdf/.png\n")
}
for (b in c("min", "max")) mkfig(b)

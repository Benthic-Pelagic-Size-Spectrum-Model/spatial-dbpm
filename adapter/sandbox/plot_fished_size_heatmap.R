# plot_fished_size_heatmap.R -- supplement figure: minimum FISHED size through time by
# LME/FAO region, pelagic (U) vs benthic (V). Reads the derived per-spectrum window
# (fished_size_UV_tv.csv) + region names from the augmented parquets. Outputs vector PDF + PNG.
#   Rscript plot_fished_size_heatmap.R
# env: PARQUET_DIR (for region_name), CSV (fished_size_UV_tv.csv), OUT (basename, no ext)
suppressMessages({library(ggplot2); library(tidyr); library(arrow)})

csv <- Sys.getenv("CSV", "adapter/sandbox/fished_size_UV_tv.csv")
if (!file.exists(csv)) csv <- "fished_size_UV_tv.csv"   # if run from adapter/sandbox
pdir <- Sys.getenv("PARQUET_DIR", "~/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs")
out  <- Sys.getenv("OUT", "adapter/sandbox/figs/fished_size_UV_heatmap")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

d <- read.csv(csv); d <- d[d$year >= 1950 & d$year <= 2010, ]
files <- Sys.glob(file.path(pdir, "dbpm_clim-fish-inputs_fao_lme-*_*.parquet"))
nm <- do.call(rbind, lapply(files, function(f) {
  L  <- as.integer(sub(".*fao_lme-([0-9]+)_.*", "\\1", basename(f)))
  rn <- read_parquet(f, col_select = "region_name")$region_name[1]
  pre <- if (L < 100) paste0("LME", L) else paste0("FAO", L - 100)
  data.frame(region = L, lab = sprintf("%-6s %s", pre, substr(rn, 1, 20))) }))
d <- merge(d, nm, by = "region")

L <- pivot_longer(d, c(min_fished_U, min_fished_V), names_to = "spectrum", values_to = "minlog")
L$spectrum <- factor(ifelse(L$spectrum == "min_fished_U", "U  (pelagic)", "V  (benthic)"))
L$lab <- factor(L$lab, levels = nm$lab[order(-nm$region)])

p <- ggplot(L, aes(year, lab, fill = minlog)) +
  geom_tile() + facet_wrap(~spectrum) +
  scale_fill_viridis_c(name = "min fished\nsize", option = "magma", na.value = "grey85",
    breaks = c(0, 1, 2, 2.43, 3.86), labels = c("1 g", "10 g", "100 g", "270 g", "7 kg")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1950, 1970, 1990, 2010)) +
  labs(title = "Minimum fished size through time by region -- pelagic (U) vs benthic (V) spectrum",
       subtitle = "Derived per region-year from catch composition; grey = no fishery (no reported catch). Colour = log10 g.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(axis.text.y = element_text(size = 5.2), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 11))

ggsave(paste0(out, ".pdf"), p, width = 11, height = 13)
ggsave(paste0(out, ".png"), p, width = 11, height = 13, dpi = 150)
cat("wrote", out, ".pdf/.png\n")

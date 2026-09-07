# plot_obs_model_panels.R -- supplement figure: observed vs modelled catch per region, BOTH groups
# (pelagic U, benthic V), from the v2 two-group calibration rds (calib_uv/lme*.rds).
# Pelagic = blue, benthic = red; observed = solid, model = dashed. Paginated PDF + PNG (page 1).
#   Rscript adapter/sandbox/plot_obs_model_panels.R
# env: CALIB_DIR (calib_uv), OUT (basename), PERPAGE (panels/page)
suppressMessages({library(ggplot2); library(tidyr); library(dplyr)})

cdir <- Sys.getenv("CALIB_DIR", "~/dbpm_compare_scratch/calib_uv")
out  <- Sys.getenv("OUT", "adapter/sandbox/figs/obs_model_panels")
pp   <- as.integer(Sys.getenv("PERPAGE", "24"))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

r <- lapply(Sys.glob(file.path(path.expand(cdir), "lme*.rds")), readRDS)
mk <- function(x) {
  pre <- if (x$L < 100) paste0("LME", x$L) else paste0("FAO", x$L - 100)
  lab <- sprintf("%s %s  (rU=%.2f rV=%.2f)", pre, substr(x$region, 1, 16),
                 ifelse(is.null(x$corr_pel) || is.na(x$corr_pel), NA, x$corr_pel),
                 ifelse(is.null(x$corr_ben) || is.na(x$corr_ben), NA, x$corr_ben))
  data.frame(L = x$L, lab = lab, year = x$year,
    Uobs = x$obs_pel, Umodel = x$model_pel, Vobs = x$obs_ben, Vmodel = x$model_ben)
}
D <- bind_rows(lapply(r, mk)) |>
  pivot_longer(-c(L, lab, year), names_to = "series", values_to = "catch") |>
  mutate(group = ifelse(grepl("^U", series), "U pelagic", "V benthic"),
         type  = ifelse(grepl("obs$", series), "obs", "model")) |>
  filter(is.finite(catch), catch > 0)
D$group <- factor(D$group, c("U pelagic", "V benthic"))

labs <- unique(D$lab[order(D$L)]); npg <- ceiling(length(labs) / pp)
pdf(paste0(out, ".pdf"), width = 12, height = 9)
for (i in seq_len(npg)) {
  ll <- labs[((i - 1) * pp + 1):min(i * pp, length(labs))]
  sub <- D[D$lab %in% ll, ]; sub$lab <- factor(sub$lab, levels = ll)
  p <- ggplot(sub, aes(year, catch, colour = group, linetype = type)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~lab, scales = "free_y", ncol = 6) +
    scale_y_log10() +
    scale_colour_manual(values = c("U pelagic" = "#3366CC", "V benthic" = "#CC3333")) +
    scale_linetype_manual(values = c(obs = "solid", model = "22")) +
    labs(title = sprintf("Observed vs modelled catch by region -- two-group calibration (page %d/%d)", i, npg),
         subtitle = "Pelagic (blue) / benthic (red); observed = solid, model = dashed. y = catch (g m-2 yr-1, log).",
         x = NULL, y = NULL, colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 7) +
    theme(strip.text = element_text(size = 6.2), legend.position = "top",
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10))
  print(p)
  if (i == 1) ggsave(paste0(out, "_p1.png"), p, width = 12, height = 9, dpi = 130)
}
invisible(dev.off())
cat(sprintf("wrote %s.pdf (%d pages) + _p1.png\n", out, npg))

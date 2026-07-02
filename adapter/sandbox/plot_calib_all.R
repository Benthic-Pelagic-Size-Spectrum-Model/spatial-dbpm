# Plot modelled vs observed catch time series for every calibrated LME, with
# per-panel calibration stats (Q, correlation, log-MSE) in the corner.
# Reads the per-LME RDS written by tier1_fishing_calib.R (SAVE_RDS=...).
#
#   Rscript adapter/sandbox/plot_calib_all.R [rds_dir] [out.png]
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })
args    <- commandArgs(TRUE)
rds_dir <- ifelse(length(args) >= 1, args[1], "/tmp/calib_all")
out_png <- ifelse(length(args) >= 2, args[2], "/tmp/calib_all/catch_calibration_all_LMEs.png")

files <- Sys.glob(file.path(rds_dir, "lme*.rds"))
res   <- lapply(files, function(f) tryCatch(readRDS(f), error = function(e) NULL))
res   <- Filter(function(r) !is.null(r) && !is.null(r$year) && any(is.finite(r$model)), res)
cat(sprintf("loaded %d LME results\n", length(res)))

# per-observation long frame + per-panel stats
ts <- bind_rows(lapply(res, function(r) tibble(
  L = r$L, year = r$year, Observed = r$obs, Modelled = r$model)))
st <- bind_rows(lapply(res, function(r) tibble(
  L = r$L, region = ifelse(is.na(r$region), "", substr(r$region, 1, 20)),
  Q = r$Q, corr = r$corr, mse = r$mse, depth = r$depth,
  pel = r$s_pel)))

# ordered facet label "LME-3  California Current"
st <- st |> mutate(lab = sprintf("LME-%d  %s", L, region)) |> arrange(L)
lev <- st$lab
ts <- ts |> left_join(st |> select(L, lab), by = "L") |>
  mutate(lab = factor(lab, levels = lev))
st <- st |> mutate(lab = factor(lab, levels = lev),
                   txt = sprintf("Q=%.2f  r=%.2f\nMSE=%.3f", Q, corr, mse))

tsl <- ts |> pivot_longer(c(Observed, Modelled), names_to = "Series", values_to = "catch") |>
  filter(is.finite(catch))

ncol <- 8; nlme <- nrow(st); nrow_facet <- ceiling(nlme / ncol)
p <- ggplot(tsl, aes(year, catch, colour = Series)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~lab, scales = "free_y", ncol = ncol) +
  geom_text(data = st, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.15, size = 2.1, colour = "grey20", lineheight = 0.9) +
  scale_colour_manual(values = c(Observed = "black", Modelled = "#c0392b")) +
  labs(x = "Year", y = expression(Catch~(g~m^-2~yr^-1)),
       colour = NULL,
       title = "DBPM (dbpmr engine) single-Q calibration: modelled vs observed catch, per FAO-LME",
       subtitle = "gravity split + time-varying plankton/temperature/sinking; Q∈[0,3], log-MSE objective") +
  theme_bw(base_size = 8) +
  theme(legend.position = "top", strip.text = element_text(size = 6.2),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8.5))

ggsave(out_png, p, width = ncol * 2.15, height = nrow_facet * 1.7 + 0.6,
       dpi = 110, limitsize = FALSE)
cat("wrote", out_png, sprintf("(%d panels, %dx%d grid)\n", nlme, nrow_facet, ncol))

# also a compact summary table sorted by fit
summ <- st |> arrange(desc(corr)) |>
  mutate(across(c(Q, corr, mse, pel), ~round(.x, 3)), depth = round(depth))
saveRDS(summ, file.path(rds_dir, "calib_summary.rds"))
cat("\n=== fit summary (sorted by corr) ===\n"); print(as.data.frame(summ[, c("L","region","depth","Q","corr","mse","pel")]), row.names = FALSE)

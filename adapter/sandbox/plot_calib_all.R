# Plot modelled vs observed catch time series for every calibrated LME, with
# per-panel calibration stats (Q, correlation, log-MSE) in the corner.
# Reads the per-LME RDS written by tier1_fishing_calib.R (SAVE_RDS=...).
#
#   Rscript adapter/sandbox/plot_calib_all.R [rds_dir] [out_prefix] [per_page] [sort]
#     per_page  panels per page (default 20 -> ~4 pages for 82 LMEs; 0 = single page)
#     sort      'lme' (default) or 'fit' (by descending correlation)
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })
args     <- commandArgs(TRUE)
rds_dir  <- if (length(args) >= 1) args[1] else "/tmp/calib_all"
out_pref <- if (length(args) >= 2) args[2] else file.path(rds_dir, "catch_calibration")
per_page <- if (length(args) >= 3) as.integer(args[3]) else 20L
sort_by  <- if (length(args) >= 4) args[4] else "lme"

files <- Sys.glob(file.path(rds_dir, "lme*.rds"))
res   <- lapply(files, function(f) tryCatch(readRDS(f), error = function(e) NULL))
res   <- Filter(function(r) !is.null(r) && !is.null(r$year) && any(is.finite(r$model)), res)
cat(sprintf("loaded %d LME results\n", length(res)))

ts <- bind_rows(lapply(res, function(r) tibble(
  L = r$L, year = r$year, Observed = r$obs, Modelled = r$model)))
st <- bind_rows(lapply(res, function(r) tibble(
  L = r$L, region = ifelse(is.na(r$region), "", substr(r$region, 1, 20)),
  Q = r$Q, corr = r$corr, mse = r$mse, depth = r$depth, pel = r$s_pel)))

# order panels (LME number, or by fit) and assign page numbers
st <- st |> arrange(if (sort_by == "fit") desc(corr) else L)
if (per_page <= 0) per_page <- nrow(st)
st <- st |> mutate(lab  = sprintf("LME-%d  %s", L, region),
                   page = ((row_number() - 1L) %/% per_page) + 1L,
                   txt  = sprintf("Q=%.2f  r=%.2f\nMSE=%.3f", Q, corr, mse))
lev <- st$lab
ts  <- ts |> left_join(st |> select(L, lab, page), by = "L")
st  <- st |> mutate(lab = factor(lab, levels = lev))
ts  <- ts |> mutate(lab = factor(lab, levels = lev))

tsl <- ts |> pivot_longer(c(Observed, Modelled), names_to = "Series", values_to = "catch") |>
  filter(is.finite(catch))

npages <- max(st$page); ncol <- 5
# if out_pref ends in .pdf -> one multi-page PDF; else one PNG per page (<pref>_pN.png)
to_pdf <- grepl("\\.pdf$", out_pref, ignore.case = TRUE)
if (to_pdf) pdf(out_pref, width = ncol * 2.9,
                height = ceiling(per_page / ncol) * 2.3 + 0.9, onefile = TRUE)
for (pg in seq_len(npages)) {
  sp <- st  |> filter(page == pg) |> mutate(lab = droplevels(lab))
  tp <- tsl |> filter(page == pg) |> mutate(lab = droplevels(lab))
  nlme <- nrow(sp); nr <- ceiling(nlme / ncol)
  p <- ggplot(tp, aes(year, catch, colour = Series)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~lab, scales = "free_y", ncol = ncol) +
    geom_text(data = sp, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE,
              hjust = -0.08, vjust = 1.15, size = 3.0, colour = "grey20", lineheight = 0.9) +
    scale_colour_manual(values = c(Observed = "black", Modelled = "#c0392b")) +
    labs(x = "Year", y = expression(Catch~(g~m^-2~yr^-1)), colour = NULL,
         title = sprintf("DBPM (dbpmr) single-Q calibration: modelled (red) vs observed (black) catch  -  page %d/%d",
                         pg, npages),
         subtitle = "gravity split + time-varying plankton/temperature/sinking; Q in [0,3], log-MSE objective") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", strip.text = element_text(size = 8.5),
          panel.grid.minor = element_blank())
  if (to_pdf) {
    print(p); cat(sprintf("  page %d/%d (%d panels)\n", pg, npages, nlme))
  } else {
    out <- sprintf("%s_p%d.png", out_pref, pg)
    ggsave(out, p, width = ncol * 2.9, height = nr * 2.3 + 0.7, dpi = 120, limitsize = FALSE)
    cat(sprintf("wrote %s  (%d panels, %dx%d)\n", out, nlme, nr, ncol))
  }
}
if (to_pdf) { invisible(dev.off()); cat("wrote", out_pref, "\n") }

# build_effort_split.R -- split the ISIMIP3a nominal-effort forcing into pelagic (U) and benthic (V)
# time series per region, by the SAME U/V FGroup map as the catch. Mirrors catch_split_lme.csv.
#   U (pelagic) = all fish + krill + cephalopods; V (benthic) = shrimp, lobsterscrab, demersalmollusc.
# Effort is assigned by its TARGET FGroup. region = LME (>0) else fao_area+100 (high-seas).
# Output: effort_split_lme.csv (LME, Year, class, effort)  -- effort = summed NomActive.
suppressMessages({library(arrow); library(dplyr)})

src <- Sys.getenv("EFFORT_SRC",
  "/Users/juliab6/Downloads/effort_histsoc_1841_2017_EEZ_addFAO.csv")
out <- Sys.getenv("OUT_CSV", "~/dbpm_compare_scratch/effort_split_lme.csv")
Vgrp <- c("shrimp", "lobsterscrab", "demersalmollusc")

d <- read_csv_arrow(src, col_select = c("Year", "FGroup", "LME", "fao_area", "NomActive")) |>
  mutate(region = if_else(LME == 0, fao_area + 100, as.double(LME)),
         class  = if_else(FGroup %in% Vgrp, "V", "U")) |>
  filter(is.finite(NomActive), NomActive > 0) |>
  group_by(region, Year, class) |>
  summarise(effort = sum(NomActive, na.rm = TRUE), .groups = "drop") |>
  rename(LME = region) |>
  arrange(LME, Year, class)

write.csv(as.data.frame(d), path.expand(out), row.names = FALSE)
cat(sprintf("wrote %s: %d rows, %d regions, years %d-%d\n", out, nrow(d),
            length(unique(d$LME)), min(d$Year), max(d$Year)))

# quick sanity: V effort fraction for a few regions (recent year)
chk <- d |> filter(Year == 2010) |> tidyr::pivot_wider(names_from = class, values_from = effort)
chk$Vfrac <- round(100 * chk$V / (chk$U + chk$V), 1)
cat("\n2010 V-effort fraction (sample):\n")
print(chk[chk$LME %in% c(1, 3, 5, 148, 158), ], row.names = FALSE)

# build_catch_split.R -- split observed catch (catch_histsoc Reported+IUU) into pelagic (U) and
# benthic (V) time series per region, same U/V FGroup map + region coding as build_effort_split.R.
#   U (pelagic) = all fish + krill + cephalopods; V (benthic) = shrimp, lobsterscrab, demersalmollusc.
#   region = LME (>0) else fao_area+100 (high-seas). catch = Reported + IUU (tonnes).
# Output: catch_split_lme.csv (LME, Year, class, catch). Supersedes the earlier LME-only version
# (which dropped the high-seas FAO regions 121-188, incl. Antarctic krill/toothfish 148/158/188).
suppressMessages({library(arrow); library(dplyr)})

src <- Sys.getenv("CATCH_SRC", "~/dbpm_compare_scratch/catch_histsoc.csv")
out <- Sys.getenv("OUT_CSV", "~/dbpm_compare_scratch/catch_split_lme.csv")
Vgrp <- c("shrimp", "lobsterscrab", "demersalmollusc")

d <- read_csv_arrow(path.expand(src),
       col_select = c("Year", "fao_area", "LME", "Reported", "IUU", "FGroup")) |>
  mutate(region = if_else(LME == 0, fao_area + 100, as.double(LME)),
         class  = if_else(FGroup %in% Vgrp, "V", "U"),
         catch  = Reported + IUU) |>
  filter(is.finite(catch), catch > 0) |>
  group_by(region, Year, class) |>
  summarise(catch = sum(catch, na.rm = TRUE), .groups = "drop") |>
  rename(LME = region) |>
  arrange(LME, Year, class)

write.csv(as.data.frame(d), path.expand(out), row.names = FALSE)
cat(sprintf("wrote %s: %d rows, %d regions, years %d-%d\n", out, nrow(d),
            length(unique(d$LME)), min(d$Year), max(d$Year)))

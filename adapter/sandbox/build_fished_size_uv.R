# build_fished_size_uv.R -- LAPTOP/standalone build of the per-spectrum, time-varying
# fished-size window (U pelagic / V benthic), written to fished_size_UV_tv.csv.
#
# This MIRRORS the derivation now embedded in the workflow at
#   lme-workflow/scripts/04_processing_effort_fishing_inputs.R  (uv_summ block)
# which writes the same columns (min/max_fished_U, min/max_fished_V) into the DBPM parquet.
# Use this script only to regenerate the committed CSV on a laptop (before the parquet has
# been re-generated on Gadi); the calibration/driver otherwise read the columns from the parquet.
#
# WHY (Julia): the parquet's single min/max_fished_weight_class is fish-derived and taxon-less;
# applied to both spectra it excludes krill (U end) and the small-bodied benthic fishery
# (shrimp/mollusc, V end). This splits the FGroup catch into U vs V, includes krill/inverts,
# and is time-varying (fishery composition shifts through the years, e.g. FAO 48 finfish->krill).
#
# METHOD: FGroup -> [lo,hi] gram range + class MIDPOINT (mid). Window per (region,year,spectrum)
# = [min mid, max mid] over FGroups holding >=0.5% of that spectrum's catch (drops trace bycatch);
# midpoints (not edges) match the parquet's log10mid_wt convention. Single-group/degenerate years
# widen to the group's full [lo,hi]. U = fish+krill+cephalopods; V = shrimp,lobsterscrab,mollusc.
# CATCH_HISTSOC env or default path points at catch_histsoc_1869_2017_EEZ_addFAO.csv.
suppressMessages({library(dplyr); library(readr); library(tidyr)})

src <- Sys.getenv("CATCH_HISTSOC", "~/dbpm_compare_scratch/catch_histsoc.csv")
out <- Sys.getenv("OUT_CSV",       "fished_size_UV_tv.csv")

# invert [min FISHED, max] g. min floored ~10 g (gear retention) EXCEPT krill (fine-mesh, 1 g).
inv <- list(krill=c(1,2), shrimp=c(10,60), lobsterscrab=c(100,4000),
            cephalopods=c(20,6000), demersalmollusc=c(20,500))
Vgrp <- c("shrimp","lobsterscrab","demersalmollusc")
# cm_lo = realistic MINIMUM FISHED length (cm) -- smallest classes anchored at 10 cm (~10 g via
# 0.01*L^3) = typical small-pelagic gear onset (anchovy/sardine), NOT 4 cm (larvae, 0.6 g).
cm_lo <- function(g) if(grepl("<30cm",g))10 else if(grepl("30-90cm",g))30 else if(grepl(">=90cm",g))90 else if(grepl("<90cm",g))10 else NA
cm_hi <- function(g) if(grepl("<30cm",g))30 else if(grepl("30-90cm",g))90 else if(grepl(">=90cm",g))200 else if(grepl("<90cm",g))90 else NA
grp_lohi <- function(g){ if(!is.null(inv[[g]])) inv[[g]] else { l<-cm_lo(g); h<-cm_hi(g); if(is.na(l)) c(NA,NA) else 0.01*c(l,h)^3 } }

d <- read_csv(src, show_col_types=FALSE,
              col_select=c("Year","fao_area","LME","Reported","IUU","FGroup"))
d$catch <- d$Reported + d$IUU
d <- d[is.finite(d$catch) & d$catch>0, ]
d$region <- ifelse(d$LME>0, d$LME, d$fao_area+100)
grps <- unique(d$FGroup); lohi <- t(sapply(grps, grp_lohi)); rownames(lohi)<-grps
d$lo <- lohi[d$FGroup,1]; d$hi <- lohi[d$FGroup,2]
d$sp  <- ifelse(d$FGroup %in% Vgrp, "V", "U")
d <- d[is.finite(d$lo),]

# window = [LOWER edge of smallest fished group, UPPER edge of largest] over groups >=0.5% of
# that spectrum-year catch. min = lower edge (smallest FISHED size), not the midpoint.
w <- d |> group_by(region, Year, sp) |>
  mutate(frac = catch/sum(catch)) |> filter(frac >= 0.005) |>
  summarise(mn = log10(min(lo)), mx = log10(max(hi)), .groups="drop") |>
  select(region, year=Year, sp, mn, mx) |>
  pivot_wider(names_from=sp, values_from=c(mn,mx)) |>
  rename(min_fished_U=mn_U, max_fished_U=mx_U, min_fished_V=mn_V, max_fished_V=mx_V) |>
  arrange(region, year)

write_csv(w, out)
cat(sprintf("wrote %s: %d region-years, %d regions\n", out, nrow(w), length(unique(w$region))))

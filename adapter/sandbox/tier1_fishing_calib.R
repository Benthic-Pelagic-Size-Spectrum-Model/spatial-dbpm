# Tier-1 sandbox: single-Q fishing calibration (transient, one FAO-LME)
# ---------------------------------------------------------------------
# Estimate ONE catchability Q by minimising the difference between the MODELLED
# and OBSERVED catch TIME SERIES (log MSE), given the model + forcing:
#   F(x,t) = Q * s(x) * effort_norm(t)   -- s(x) knife-edge at log10(w)=1 (10 g),
#   effort normalised to [0,1] over the LME series, Q bounded [0,3], A fixed = 64.
# Per the FishMIP DBPM.md spec (Fish-MIP/Global_MEM_Model_Templates).
#
# CAVEAT (Tier-1): environmental forcing is held CONSTANT here (only effort varies)
# - the catch time series is driven by effort + fishing dynamics, not environmental
# variability. Time-varying temperature/plankton forcing is issue #11 (needs the
# in-memory column driver #5). So this validates the Q-calibration mechanism; a
# full fit adds the environmental time series.
#
#   DBPM_DATA=/path/to/DBPM_dev Rscript adapter/sandbox/tier1_fishing_calib.R [LME]


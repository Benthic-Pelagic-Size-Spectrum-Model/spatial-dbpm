#!/bin/sh
# Launch the v2 two-group calibration for all 83 regions, 6 in parallel (RAM disk).
# Writes calib_uv/lme<L>.rds; touches UV_CALIB_DONE.marker when finished.
OUT="$HOME/dbpm_compare_scratch/calib_uv"; mkdir -p "$OUT"
ONE=/Users/juliab6/spatial-dbpm/adapter/sandbox/calib_uv_one.sh
rm -f "$HOME/dbpm_compare_scratch/UV_CALIB_DONE.marker"
awk -F, 'NR>1{print $1}' "$HOME/dbpm_compare_scratch/all_dint.csv" | \
  xargs -P 6 -I{} sh "$ONE" {}
echo ALL_UV_CALIB_DONE
touch "$HOME/dbpm_compare_scratch/UV_CALIB_DONE.marker"

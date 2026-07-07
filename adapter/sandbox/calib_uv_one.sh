#!/bin/sh
# One-region v2 TWO-GROUP calibration (q_pel, q_ben to separate catch series). Reads the
# per-spectrum fished-size window + effort/catch splits. TMPDIR isolated on the RAM disk (dbpmr
# file-I/O bottleneck) and cleaned. Output: calib_uv/lme<L>.rds + .log
L=$1
export R_LIBS="$HOME/dbpm_compare_scratch/dbpmrlib"
export DINT_CSV="$HOME/dbpm_compare_scratch/lme_dint_hw.csv"
export INPUT_PARQUET_DIR="$HOME/dbpm_compare_scratch/DBPM_dev_uv/dbpm_inputs"
export EFFORT_SPLIT_CSV="$HOME/dbpm_compare_scratch/effort_split_lme.csv"
export CATCH_SPLIT_CSV="$HOME/dbpm_compare_scratch/catch_split_lme.csv"
TIER1=/Users/juliab6/spatial-dbpm/adapter/sandbox/tier1_fishing_calib.R
OUT="$HOME/dbpm_compare_scratch/calib_uv_hw"
RAM=/Volumes/dbpmram; [ -d "$RAM" ] || RAM=/tmp
TMPD=$(mktemp -d "$RAM/uv$L.XXXXXX"); export TMPDIR="$TMPD"
[ -f "$OUT/lme$L.rds" ] || QMAX=4 SAVE_RDS="$OUT/lme$L.rds" Rscript --vanilla "$TIER1" "$L" > "$OUT/lme$L.log" 2>&1
rm -rf "$TMPD"
echo "done LME $L $(date +%H:%M:%S)"

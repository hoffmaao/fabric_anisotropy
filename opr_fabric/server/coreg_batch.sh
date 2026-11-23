#!/bin/bash
# Build the quad-pol coregistration caches for every frame of a survey,
# K frames at a time:   coreg_batch.sh <site> [K] [mode]     (default K=3)
#
# mode is coreg_one.sh's: 'firn' for the node method (one delay per channel
# from the firn, about a minute a frame, caches under
# stages/quadpol_cf/coreg_cache), 'tiles' (default) for the least-squares
# pipeline (the OPR tile field, ~50 min a frame, stages/quadpol/coreg_cache).
#
# Sites mirror run_quadpol_survey.m. The caches (~680 MB/frame) are
# estimator-independent, so this is the one pass over a survey that reads
# the raw channels; every later rerun against them is minutes. Frames
# whose cache exists are skipped, so re-running after failures is cheap
# and safe. Frames without a
# polarimetric product (7 of the 36 Thwaites quad-pol frames) are
# reported and skipped: the pipeline takes its window and coregistration
# settings from that product.
#
#   nohup bash coreg_batch.sh ridge_a 3 firn > coreg_ridge_a.log 2>&1 &
set -u
# Workers come from this script's own directory, <work> from the walk in
# fabric_paths.sh. FABRIC_ROOT overrides <work> only. No username is baked
# in, and nothing has to be deployed anywhere.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
site=${1:?usage: coreg_batch.sh <site> [K] [mode]}
K=${2:-3}
MODE=${3:-tiles}
case $MODE in
  tiles) STAGE=quadpol;    SFX= ;;
  firn)  STAGE=quadpol_cf; SFX=_firn ;;
  *) echo "coreg_batch.sh: mode is '$MODE'; it must be tiles or firn"; exit 1 ;;
esac
CDIR=$FAB/stages/$STAGE/coreg_cache
case $site in
  ridge_a)     ROOT=/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2 ;;
  taylor_dome) ROOT=/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2 ;;
  thwaites)    ROOT=/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground ;;
  eastwind)    ROOT=/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground ;;
  *) echo "unknown site $site"; exit 1 ;;
esac
HH=$ROOT/CSARP_standardphase_HH
mkdir -p "$FAB/coreg_logs"
list=$FAB/coreg_logs/frames_$site$SFX.txt
find "$HH" -name 'Data_*.mat' ! -name 'Data_img_*' | sort | \
while read -r f; do
  seg=$(basename "$(dirname "$f")")
  n=$(basename "$f" .mat)
  n=${n##*_}
  echo "$ROOT $seg $((10#$n)) $MODE"
done > "$list"
total=$(wc -l < "$list")
echo "=== $site ($MODE): $total frames, $K concurrent, started $(date)"
xargs -P "$K" -L 1 bash "$FAB_DIR/coreg_one.sh" < "$list"
have=$(cut -d' ' -f2,3 "$list" | while read -r s n; do
  printf '%s_%03d\n' "$s" "$n"
done | while read -r tag; do
  [ -s "$CDIR/creg_$tag.mat" ] && echo x
done | wc -l)
echo "=== $site ($MODE) finished $(date): $have of $total caches present in $CDIR"
locks=$(ls -d "$CDIR"/.lock_* 2>/dev/null || true)
if [ -n "$locks" ]; then
  echo "WARNING: stale locks (remove and re-run those frames):"
  echo "$locks"
fi

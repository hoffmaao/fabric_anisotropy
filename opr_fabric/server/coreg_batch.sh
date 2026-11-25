#!/bin/bash
# Coregister every frame of a survey into the cache, K frames at a time:
#
#   coreg_batch.sh <site> [K]     (default K=3)
#
#   nohup bash coreg_batch.sh ridge_a > coreg_ridge_a.log 2>&1 &
#
# Step 1 of the node method for a whole survey: coreg_one.sh on every frame
# of the site. The caches (0.3-3.2 GB a frame at Ridge A, under
# <work>/stages/quadpol_cf/coreg_cache) are the one pass that reads the
# raw channels; everything after works from them. Frames whose cache
# exists are skipped, so re-running after failures is cheap and safe.
# Frames without a polarimetric product (11 of the 47 at Ridge A) are
# reported and skipped: the range window and the surface come from it.
set -u
# Workers come from this script's own directory, <work> from the walk in
# fabric_paths.sh. FABRIC_ROOT overrides <work> only. No username is baked
# in, and nothing has to be deployed anywhere.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
usage="usage: coreg_batch.sh <site> [K]"
site=${1:?$usage}
K=${2:-3}
if ! [[ $K =~ ^[1-9][0-9]*$ ]]; then
  echo "coreg_batch.sh: K must be a positive integer (got '$K')"; echo "$usage"; exit 1
fi
CDIR=$FAB/stages/quadpol_cf/coreg_cache
case $site in
  ridge_a)     ROOT=/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2 ;;
  taylor_dome) ROOT=/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2 ;;
  thwaites)    ROOT=/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground ;;
  eastwind)    ROOT=/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground ;;
  *) echo "unknown site $site"; echo "$usage"; exit 1 ;;
esac
HH=$ROOT/CSARP_standardphase_HH
mkdir -p "$FAB/coreg_logs"
list=$FAB/coreg_logs/frames_$site.txt
find "$HH" -name 'Data_*.mat' ! -name 'Data_img_*' | sort | \
while read -r f; do
  seg=$(basename "$(dirname "$f")")
  n=$(basename "$f" .mat)
  n=${n##*_}
  echo "$ROOT $seg $((10#$n))"
done > "$list"
total=$(wc -l < "$list")
echo "=== $site: $total frames, $K concurrent, started $(date)"
xargs -P "$K" -L 1 bash "$FAB_DIR/coreg_one.sh" < "$list"
have=$(cut -d' ' -f2,3 "$list" | while read -r s n; do
  printf '%s_%03d\n' "$s" "$n"
done | while read -r tag; do
  [ -s "$CDIR/creg_$tag.mat" ] && echo x
done | wc -l)
echo "=== $site finished $(date): $have of $total caches present in $CDIR"
locks=$(ls -d "$CDIR"/.lock_* 2>/dev/null || true)
if [ -n "$locks" ]; then
  echo "WARNING: stale locks (remove and re-run those frames):"
  echo "$locks"
fi

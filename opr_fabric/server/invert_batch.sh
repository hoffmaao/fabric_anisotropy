#!/bin/bash
# Run the quad-pol inversion over every cached frame of a survey, K at a
# time, FOLLOWING the cache builder: frames are picked up as their caches
# appear, and the loop ends when no work remains and coreg_batch.sh is no
# longer running.  invert_batch.sh <site> [K]     (default K=2: polite
# alongside the 3-way cache batch; raise it once that finishes)
#
#   nohup bash invert_batch.sh ridge_a 2 > invert_chain.log 2>&1 &
set -u
FAB=/kucresis/scratch/hoffmana_sta/fabric
site=${1:?usage: invert_batch.sh <site> [K]}
K=${2:-2}
MAX_FAIL=2
case $site in
  ridge_a)     ROOT=/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2 ;;
  taylor_dome) ROOT=/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2 ;;
  thwaites)    ROOT=/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground ;;
  eastwind)    ROOT=/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground ;;
  *) echo "unknown site $site"; exit 1 ;;
esac
HH=$ROOT/CSARP_standardphase_HH
mkdir -p "$FAB/invert_logs"
sweep=0
while :; do
  sweep=$((sweep + 1))
  list=$FAB/invert_logs/work_$site.txt
  find "$HH" -name 'Data_*.mat' ! -name 'Data_img_*' | sort | \
  while read -r f; do
    seg=$(basename "$(dirname "$f")")
    n=$(basename "$f" .mat)
    n=${n##*_}
    n=$((10#$n))
    tag=$(printf '%s_%03d' "$seg" "$n")
    if [ -s "$FAB/stages/quadpol/coreg_cache/creg_$tag.mat" ] && \
       [ ! -s "$FAB/stages/quadpol/quadpol_section_$tag.mat" ]; then
      # a frame that keeps failing is dropped after MAX_FAIL attempts
      # (invert_one.sh appends one line per failure, clears on success);
      # remove fail_<tag>.count from invert_logs to retry it by hand
      fails=$FAB/invert_logs/fail_$tag.count
      if [ -f "$fails" ] && [ "$(wc -l < "$fails")" -ge "$MAX_FAIL" ]; then
        continue
      fi
      echo "$ROOT $seg $n"
    fi
  done > "$list"
  n_work=$(wc -l < "$list")
  echo "=== sweep $sweep: $n_work frames ready, $(date)"
  if [ "$n_work" -gt 0 ]; then
    xargs -P "$K" -L 1 bash "$FAB/invert_one.sh" < "$list"
  else
    # the [c] keeps pgrep from matching this script's own command line
    if pgrep -f "[c]oreg_batch.sh" > /dev/null; then
      sleep 600
    else
      break
    fi
  fi
done
done_n=$(ls "$FAB"/stages/quadpol/quadpol_section_*.mat 2>/dev/null | wc -l)
echo "=== $site inversion sweep finished $(date): $done_n sections on disk"

#!/bin/bash
# Rebuild the EastGRIP frames whose PRODUCTION reference trajectory is null.
#
# ref_20240620_01, ref_20240621_01 and ref_20240622_01 in the group tree are
# 100% null island - every one of 581926, 869130 and 301949 samples sits at
# about 0 N, 2 E - so 21 of the 29 constant-orientation products cannot be
# placed on a map and drop out of any position-selected figure. The repaired
# trajectories already exist, from the earlier GNSS re-sync: the private
# records carry real positions for every segment and egrip_records_resync
# wrote matching reference trajectories under the private opr_data tree.
# Nothing needed re-deriving; the pipeline was simply reading the group copy.
#
# So this points ref_root at the repaired tree and re-runs those segments.
# The existing products are MOVED ASIDE, not deleted: they are a record of
# what a null trajectory produces and cost hours to make.
#
# The other two segments (20240619_01, 20240620_02) are untouched - their
# group trajectories were always fine.
set -u
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
ML=${MATLAB_BIN:-/opt/sw/matlab/2024b/bin/matlab}
ROOT=/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2
REFROOT=/cresis/users/hoffmana_sta/scratch/opr_support_egrip/opr_data/accum/2024_Greenland_Ground2
SEGS="20240620_01 20240621_01 20240622_01"
QP="$FAB/stages/quadpol"
mkdir -p "$FAB/invert_logs" "$QP/null_traj_products"

cd "$FAB_DIR" || exit 1
GUARD="p=which('run_quadpol_pipeline'); q=fullfile('$FAB_DIR','run_quadpol_pipeline.m'); \
if ~strcmp(p,q), error('pipeline shadowed by %s (expected %s)',p,q); end;"

# Build the work list WITHOUT touching anything, so DRY really is dry.
# The first version moved the products before testing DRY, which meant the
# dry run did the one irreversible thing in the script.
: > "$FAB/egrip_repos.txt"
for s in $SEGS; do
  for f in "$QP"/quadpol_section_${s}_*_ct.mat "$QP/null_traj_products"/quadpol_section_${s}_*_ct.mat; do
    [ -f "$f" ] || continue
    frm=$(basename "$f" | sed -E "s/.*_([0-9]{3})_ct\.mat/\1/")
    grep -q "^$s $frm$" "$FAB/egrip_repos.txt" 2>/dev/null || echo "$s $frm" >> "$FAB/egrip_repos.txt"
  done
done
echo "$(wc -l < "$FAB/egrip_repos.txt") frame(s) to reposition"
[ -n "${DRY:-}" ] && { echo "DRY: nothing moved, nothing run"; exit 0; }

# Only now move the old products aside, keeping them as the record of what
# a null trajectory produces.
for s in $SEGS; do
  for f in "$QP"/quadpol_section_${s}_*_ct.mat; do
    [ -f "$f" ] && mv "$f" "$QP/null_traj_products/"
  done
done
echo "old products preserved in null_traj_products/"

worker () {
  while read -r seg frm; do
    tag=$(printf '%s_%03d' "$seg" "$((10#$frm))")
    out="$QP/quadpol_section_${tag}_ct.mat"
    [ -f "$out" ] && continue
    lock="$FAB/invert_logs/lock_rp_${tag}"
    mkdir "$lock" 2>/dev/null || continue
    echo "start $tag $(date +%d/%H:%M)"
    nice -n 10 $ML -batch "addpath('$FAB_DIR'); $GUARD maxNumCompThreads(8); \
      site_root='$ROOT'; ref_root='$REFROOT'; day_seg='$seg'; frm=$((10#$frm)); \
      theta_const=true; run_quadpol_pipeline" \
      > "$FAB/invert_logs/rp_${tag}.log" 2>&1
    if [ -f "$out" ]; then echo "done  $tag ok $(date +%d/%H:%M)"
    else echo "done  $tag FAIL $(date +%d/%H:%M)"; fi
    rmdir "$lock" 2>/dev/null
  done < "$FAB/egrip_repos.txt"
}

for w in 1 2 3 4 5 6; do worker & done
wait
echo "=== reposition finished $(date): $(ls $QP/quadpol_section_2024062*_ct.mat 2>/dev/null | wc -l) products"

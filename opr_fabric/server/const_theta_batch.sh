#!/bin/bash
# Re-run the quad-pol pipeline with a CONSTANT fabric axis in depth
# (theta_const), per along-track segment, over every site where that
# assumption is defensible.
#
#   const_theta_batch.sh            # DRY=1 to only build and report the work list
#
# THWAITES IS DELIBERATELY EXCLUDED, EXCEPT THREE CONTROL FRAMES. Its
# margin has depth-rotating axes - a measured result of this project, not a
# guess - and the assumption is measured 6x worse on a rotating synthetic
# (test_quadpol_const_theta). The exclusion is by TAG, not by prefix,
# because Thwaites shares the 2024 prefix with WAIS Divide (20240120) and
# McMurdo (20240202/03): a prefix filter would either drop those two or
# include the margin. The three CONTROL frames run precisely BECAUSE the
# assumption is wrong there: the per-window spread that says "the axis
# rotates" has only been calibrated on a synthetic (2 deg constant vs 25
# rotating), and compare_const_theta.py's verdict thresholds need a real
# rotating column beside the real divides before they mean anything. The
# _ct tag keeps these from touching any Thwaites number in use.
#
# Products are written with the _ct suffix the pipeline appends, so the
# existing products - and every number published from them - stay intact
# and the two can be compared frame for frame. Idempotent: existing _ct
# products are skipped, failures get capped markers, mkdir locks arbitrate.
set -u
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
ML=${MATLAB_BIN:-/opt/sw/matlab/2024b/bin/matlab}
mkdir -p "$FAB/invert_logs" "$FAB/stages/quadpol"

# MATLAB resolves the CURRENT FOLDER ahead of the path, so a stale flat copy of
# run_quadpol_pipeline.m in the launch directory silently beats the repo mirror
# no matter what we addpath. An old pipeline does not know theta_const: it drops
# the _ct tag and overwrites the very products the tag exists to protect, while
# reporting a perfectly healthy fit. That happened on 1 Sep 2026 from a work-root
# copy dated 24 Aug. Run from the repo directory, so the current folder resolves
# to the file we mean, and make MATLAB prove which file it loaded before any
# frame is touched - a silent wrong-version run is the one failure this batch
# cannot be allowed to have.
cd "$FAB_DIR" || exit 1
GUARD="p=which('run_quadpol_pipeline'); q=fullfile('$FAB_DIR','run_quadpol_pipeline.m'); \
if ~strcmp(p,q), error('pipeline shadowed by %s (expected %s)',p,q); end;"

# season root per tag family; the season NUMBER never matches the tag year
root_for () {
  case "$1" in
    2022*|2023*) echo /cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground ;;
    202401*|20240202*|20240203*) echo /cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground ;;  # Thwaites, WAIS, McMurdo
    2025*) echo /cresis/nvme/opr_data/accum/2024_Antarctica_Ground2 ;;
    2026*) echo /cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2 ;;
    *) echo "" ;;
  esac
}

# Work list from the coreg caches, exact frame names only (the _z deep
# reruns cannot be parsed by the frame regex and would pass an empty frame
# number through to matlab).
ls "$FAB"/stages/quadpol/coreg_cache/creg_*.mat 2>/dev/null \
  | grep -E '/creg_[0-9]{8}_[0-9]{2}_[0-9]{3}\.mat$' \
  | sed -E 's|.*/creg_([0-9]{8}_[0-9]{2})_([0-9]{3})\.mat|\1 \2|' \
  | sort -u > "$FAB/ct_all.txt"

: > "$FAB/ct_work.txt"
while read -r seg frm; do
  r=$(root_for "$seg")
  [ -z "$r" ] && { echo "skip $seg $frm (no season root)"; continue; }
  case "$seg" in
    2024*)
      case "$seg $frm" in
        20240120*|20240202*|20240203*) ;;                 # WAIS / McMurdo: keep
        "20240108_01 001"|"20240103_01 001"|"20240104_01 001") ;;  # Thwaites CONTROL
        *) continue ;;                                    # Thwaites: excluded
      esac ;;
  esac
  echo "$seg $frm $r" >> "$FAB/ct_work.txt"
done < "$FAB/ct_all.txt"
# Two different exclusions, and they must not be conflated in the report:
# Thwaites is dropped on purpose by the tag filter above, while the June 2024
# EastGRIP/NEGIS segments are dropped by root_for having no Greenland season
# root. The second is a gap, not a decision - EastGRIP qlook needs its own
# adaptation (stationary-trace cull, no settings product) and its first frame
# does not yet agree with the core, so it is not ready to batch.
n_thw=$(awk '{print $1}' "$FAB/ct_all.txt" | grep -E '^202401' | grep -vcE '^20240120' || true)
n_ctl=$(awk '{print $1}' "$FAB/ct_work.txt" | grep -E '^202401' | grep -vcE '^20240120' || true)
n_egrip=$(awk '{print $1}' "$FAB/ct_all.txt" | grep -cE '^202406' || true)
echo "$(wc -l < "$FAB/ct_work.txt") frames queued; \
$((n_thw - n_ctl)) Thwaites excluded by design ($n_ctl kept as rotating-axis controls), \
$n_egrip EastGRIP skipped (no Greenland season root)"
# DRY=1 builds and reports the work list without starting MATLAB
[ -n "${DRY:-}" ] && exit 0

worker () {
  while read -r seg frm root; do
    tag=$(printf '%s_%03d' "$seg" "$((10#$frm))")
    out="$FAB/stages/quadpol/quadpol_section_${tag}_ct.mat"
    [ -f "$out" ] && continue
    fail="$FAB/invert_logs/fail_ct_${tag}.count"
    n=$(cat "$fail" 2>/dev/null || echo 0)
    [ "$n" -ge 2 ] && continue
    lock="$FAB/invert_logs/lock_ct_${tag}"
    mkdir "$lock" 2>/dev/null || continue
    echo "start $tag $(date +%d/%H:%M)"
    nice -n 10 $ML -batch "addpath('$FAB_DIR'); $GUARD maxNumCompThreads(8); site_root='$root'; day_seg='$seg'; frm=$((10#$frm)); theta_const=true; run_quadpol_pipeline" \
      > "$FAB/invert_logs/ct_${tag}.log" 2>&1
    if [ -f "$out" ]; then
      echo "done  $tag ok $(date +%d/%H:%M)"; rm -f "$fail"
    else
      echo "done  $tag FAIL $(date +%d/%H:%M)"; echo $((n+1)) > "$fail"
    fi
    rmdir "$lock" 2>/dev/null
  done < "$FAB/ct_work.txt"
}
# Fail fast, once, rather than letting every frame burn its two retries: prove
# MATLAB loads the repo pipeline AND that this pipeline understands theta_const.
echo "preflight: checking which pipeline MATLAB loads"
$ML -batch "addpath('$FAB_DIR'); $GUARD \
  s=fileread(q); assert(contains(s,'PTAG'), 'pipeline predates the _ct product tag: %s', q); \
  fprintf('preflight OK: %s\n', p);" || { echo "PREFLIGHT FAILED - batch not started"; exit 1; }

# Six workers x 8 MATLAB threads on the 112-core node: a frame is ~70 min
# with the jackknife and split-half (frame pass ~50, section ~15), so 118
# frames take about a day at this width and leave two thirds of the node
# for everyone else.
for w in 1 2 3 4 5 6; do worker & done
wait
echo "=== const-theta batch finished $(date): \
$(ls "$FAB"/stages/quadpol/quadpol_section_*_ct.mat 2>/dev/null | wc -l) products"

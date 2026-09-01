#!/bin/bash
# Reflection-ratio retrieval (run_ershadi_r.m) for every frame that has a
# coreg cache, filtered by an optional tag prefix:
#   ershadi_r_batch.sh [<prefix>]      e.g. 2025 for Ridge A
# Idempotent: frames whose ershadi_r product exists are skipped; failures
# get capped fail markers; mkdir locks arbitrate the workers - so after ANY
# interruption the same launch command resumes where it left off.
set -u
# code root and <work> from the shared walk - see fabric_paths.sh. No
# username is baked in, and nothing has to be deployed anywhere.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
PREFIX=${1:-}
ML=${MATLAB_BIN:-/opt/sw/matlab/2024b/bin/matlab}
mkdir -p "$FAB/invert_logs" "$FAB/stages/ershadi_r"

# grep -E keeps ONLY exact frame caches: the non-default-depth reruns are
# named creg_<tag>_z<N>.mat, which the sed below cannot parse, and an
# unparsed line passes the whole path through as $seg with $frm empty -
# "10#: invalid integer constant" on one frame of the first Ridge A run.
# Filtering before the sed is what the other launchers do; matching them.
ls "$FAB"/stages/quadpol/coreg_cache/creg_${PREFIX}*.mat 2>/dev/null \
  | grep -E '/creg_[0-9]{8}_[0-9]{2}_[0-9]{3}\.mat$' \
  | sed -E 's|.*/creg_([0-9]{8}_[0-9]{2})_([0-9]{3})\.mat|\1 \2|' \
  | sort -u > "$FAB/ershadi_r_work.txt"
echo "$(wc -l < "$FAB/ershadi_r_work.txt") frames queued (prefix '$PREFIX')"

worker() {
  while read -r seg frm; do
    tag=$(printf '%s_%03d' "$seg" "$((10#$frm))")
    out="$FAB/stages/ershadi_r/ershadi_r_${tag}.mat"
    [ -f "$out" ] && continue
    fail="$FAB/invert_logs/fail_er_${tag}.count"
    n=$(cat "$fail" 2>/dev/null || echo 0)
    [ "$n" -ge 2 ] && continue
    lock="$FAB/invert_logs/lock_er_${tag}"
    mkdir "$lock" 2>/dev/null || continue
    echo "start $tag $(date +%d/%H:%M)"
    nice -n 10 $ML -batch "addpath('$FAB_DIR'); maxNumCompThreads(8); day_seg='$seg'; frm=$((10#$frm)); run_ershadi_r" \
      > "$FAB/invert_logs/er_${tag}.log" 2>&1
    if [ -f "$out" ]; then
      echo "done  $tag ok $(date +%d/%H:%M)"; rm -f "$fail"
    else
      echo "done  $tag FAIL $(date +%d/%H:%M)"; echo $((n+1)) > "$fail"
    fi
    rmdir "$lock" 2>/dev/null
  done < "$FAB/ershadi_r_work.txt"
}
for w in 1 2 3; do worker & done
wait
echo "=== ershadi_r batch finished $(date): $(ls "$FAB"/stages/ershadi_r/ershadi_r_${PREFIX}*.mat 2>/dev/null | wc -l) products"

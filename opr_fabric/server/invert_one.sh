#!/bin/bash
# Run the quad-pol pipeline (both estimators + section) for ONE frame,
# from its coreg cache. Worker for invert_batch.sh; fine by hand:
#   invert_one.sh <site_root> <day_seg> <frm>
#
# Preconditions enforced here so the batch stays collision-free with the
# cache builder: the cache must already exist (frames without one belong
# to coreg_batch.sh) and the section output must not (idempotent reruns).
set -u
# <work> comes from the walk in fabric_paths.sh, so this worker lands on
# the same root whether it is run from the checkout or from a copy in
# <work>, and on the same root MATLAB derives. FABRIC_ROOT overrides.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
root=$1
seg=$2
frm=$3
tag=$(printf '%s_%03d' "$seg" "$frm")
cache=$FAB/stages/quadpol/coreg_cache/creg_$tag.mat
outmat=$FAB/stages/quadpol/quadpol_section_$tag.mat
log=$FAB/invert_logs/inv_$tag.log
if [ ! -s "$cache" ]; then
  echo "skip  $tag (no cache yet)"
  exit 0
fi
if [ -s "$outmat" ]; then
  echo "skip  $tag (section exists)"
  exit 0
fi
lock=$FAB/stages/quadpol/.invlock_$tag
if ! mkdir "$lock" 2>/dev/null; then
  echo "skip  $tag (locked)"
  exit 0
fi
trap 'rmdir "$lock" 2>/dev/null' EXIT
echo "start $tag $(date '+%H:%M')"
nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \
  "addpath('$FAB_DIR'); maxNumCompThreads(8); site_root='$root'; day_seg='$seg'; frm=$frm; run_quadpol_pipeline" \
  > "$log" 2>&1
rc=$?
if [ -s "$outmat" ]; then
  st=ok
  rm -f "$FAB/invert_logs/fail_$tag.count"
else
  st=FAILED
  # one line per failed attempt; invert_batch.sh stops re-listing the
  # frame once this reaches its retry cap
  echo "$(date '+%Y-%m-%d %H:%M') rc=$rc" >> "$FAB/invert_logs/fail_$tag.count"
fi
echo "done  $tag rc=$rc $st $(date '+%H:%M')"
exit 0

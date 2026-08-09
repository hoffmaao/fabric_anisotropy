#!/bin/bash
# Run the quad-pol pipeline (both estimators + section) for ONE frame,
# from its coreg cache. Worker for invert_batch.sh; fine by hand:
#   invert_one.sh <site_root> <day_seg> <frm>
#
# Preconditions enforced here so the batch stays collision-free with the
# cache builder: the cache must already exist (frames without one belong
# to coreg_batch.sh) and the section output must not (idempotent reruns).
set -u
FAB=/kucresis/scratch/hoffmana_sta/fabric
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
  "maxNumCompThreads(8); site_root='$root'; day_seg='$seg'; frm=$frm; run_quadpol_pipeline" \
  > "$log" 2>&1
rc=$?
if [ -s "$outmat" ]; then st=ok; else st=FAILED; fi
echo "done  $tag rc=$rc $st $(date '+%H:%M')"
exit 0

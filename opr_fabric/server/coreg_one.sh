#!/bin/bash
# Coregister ONE quad-pol frame into the cache. Worker for coreg_batch.sh;
# also fine by hand: coreg_one.sh <site_root> <day_seg> <frm>
#
# Idempotent: a frame whose cache file already exists is skipped, so the
# batch can be re-run after any failure (license denial, dropped node) and
# only the missing frames pay again. The mkdir lock guards against two
# batches racing on one frame; a stale lock from a hard kill just makes
# the frame report "locked" - remove .lock_<tag> and re-run.
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
name=$(printf 'Data_%s_%03d.mat' "$seg" "$frm")
cache=$FAB/stages/quadpol/coreg_cache/creg_$tag.mat
log=$FAB/coreg_logs/creg_$tag.log
if [ -s "$cache" ]; then
  echo "skip  $tag (cache exists)"
  exit 0
fi
# the pipeline reads its window and settings from the polarimetric
# product (plain or _unwrap); without one the frame cannot run, so say so
# here instead of burning a MATLAB launch to find out
if [ ! -s "$root/CSARP_polarimetric/$seg/$name" ] && \
   [ ! -s "$root/CSARP_polarimetric_unwrap/$seg/$name" ]; then
  echo "skip  $tag (no polarimetric product)"
  exit 0
fi
lock=$FAB/stages/quadpol/coreg_cache/.lock_$tag
if ! mkdir "$lock" 2>/dev/null; then
  echo "skip  $tag (locked)"
  exit 0
fi
trap 'rmdir "$lock" 2>/dev/null' EXIT
echo "start $tag $(date '+%H:%M')"
# nice + a thread cap: coregistration would happily take every core, and
# the node is shared. 8 threads x 3-4 concurrent frames is the polite
# footprint; raise the batch's K rather than the threads if it is idle.
nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \
  "addpath('$FAB_DIR'); maxNumCompThreads(8); site_root='$root'; day_seg='$seg'; frm=$frm; coreg_only=true; run_quadpol_pipeline" \
  > "$log" 2>&1
rc=$?
if [ -s "$cache" ]; then st=ok; else st=FAILED; fi
echo "done  $tag rc=$rc $st $(date '+%H:%M')"
exit 0

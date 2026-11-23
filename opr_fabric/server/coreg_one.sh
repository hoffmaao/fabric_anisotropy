#!/bin/bash
# Coregister ONE frame into the cache: run_coreg_frame.m under a lock, with
# a log. Worker for coreg_batch.sh; also fine by hand:
#
#   coreg_one.sh <site_root> <day_seg> <frm>
#
# One range delay per channel, fitted in the firn; 1-3 min a frame. The
# cache is <work>/stages/quadpol_cf/coreg_cache/creg_<tag>.mat and the log
# <work>/coreg_logs/creg_<tag>.log.
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
usage="usage: coreg_one.sh <site_root> <day_seg> <frm>"
root=${1:?$usage}
seg=${2:?$usage}
frm=${3:?$usage}
tag=$(printf '%s_%03d' "$seg" "$frm")
name=$(printf 'Data_%s_%03d.mat' "$seg" "$frm")
cdir=$FAB/stages/quadpol_cf/coreg_cache
cache=$cdir/creg_$tag.mat
log=$FAB/coreg_logs/creg_$tag.log
if [ -s "$cache" ]; then
  echo "skip  $tag (cache exists)"
  exit 0
fi
# the range window and the surface come from the polarimetric product
# (plain or _unwrap); without one the frame cannot run, so say so here
# instead of burning a MATLAB launch to find out
if [ ! -s "$root/CSARP_polarimetric/$seg/$name" ] && \
   [ ! -s "$root/CSARP_polarimetric_unwrap/$seg/$name" ]; then
  echo "skip  $tag (no polarimetric product)"
  exit 0
fi
# the cache and log directories, made here and not left to the batch: by
# hand on a fresh work root neither exists yet, and a lock that cannot be
# made for want of its parent would be reported as "locked"
mkdir -p "$cdir" "$FAB/coreg_logs" || exit 1
lock=$cdir/.lock_$tag
if ! mkdir "$lock" 2>/dev/null; then
  echo "skip  $tag (locked)"
  exit 0
fi
trap 'rmdir "$lock" 2>/dev/null' EXIT
echo "start $tag $(date '+%H:%M')"
# nice + a thread cap: MATLAB would happily take every core, and the node
# is shared. 8 threads x 3-4 concurrent frames is the polite footprint;
# raise the batch's K rather than the threads if it is idle.
nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \
  "addpath('$FAB_DIR'); maxNumCompThreads(8); site_root='$root'; day_seg='$seg'; frm=$frm; run_coreg_frame" \
  > "$log" 2>&1
rc=$?
if [ -s "$cache" ]; then st=ok; else st=FAILED; fi
echo "done  $tag rc=$rc $st $(date '+%H:%M')"
exit 0

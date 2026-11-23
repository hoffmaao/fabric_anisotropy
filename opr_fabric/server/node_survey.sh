#!/bin/bash
# NODE_SURVEY The node method over a survey: run_node_frame on every
# firn-coregistered frame of a site, a few at a time.
#
#   nohup bash node_survey.sh <site> [max_jobs] [threads] [pattern] [blk_m] [tagx] [extra] [cache_dir] > node_<site>.log 2>&1 &
#
#   site       ridge_a, taylor_dome, thwaites or eastwind (coreg_batch.sh's)
#   max_jobs   frames fitted at once (4), each its own MATLAB
#   threads    MATLAB threads per frame (4)
#   pattern    grep expression on the frame tag; default every frame.
#              20250111_01 is one segment, 20250111_01_003 one frame
#   blk_m      block length in metres (2000)
#   tagx       suffix for a variant pass, so it never replaces the
#              standard products ('' for the standard pass)
#   extra      MATLAB code run before run_node_frame, e.g. "psi_step_deg=2.5;"
#   cache_dir  where the creg_<tag>.mat caches are; default
#              <work>/stages/quadpol_cf/coreg_cache, which
#              `coreg_batch.sh <site> <K> firn` builds
#
# An empty argument ("") takes its default, so a later one can be given
# without the ones before it.
#
# Every creg_<tag>.mat in cache_dir that is one of the site's frames and
# matches pattern is fitted; there is no calibration stage, the channel
# gains come from each block. Products are
# <work>/stages/node/node_<tag>_B<blk_m><tagx>.mat and logs
# <work>/logs/node_<tag><run>.log, where run is tagx for 2 km blocks and
# _B<blk_m><tagx> for any other length.
#
# Idempotent and cooperative: a frame whose product exists is skipped, a
# frame being fitted holds locks/node_<tag><run>, and a second launch of
# the same pass exits at once on locks/node_survey_<site><run>. A run that
# was killed leaves its locks behind: remove <work>/locks/node_* by hand
# before relaunching. MATLAB runs from <work>/logs, a directory with no .m
# files (a script in the current directory shadows the path copy).
set -u
# Workers come from this script's own directory, <work> from the walk in
# fabric_paths.sh, the same root MATLAB's fabric_paths derives.
# FABRIC_ROOT overrides <work> only.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
usage="usage: node_survey.sh <site> [max_jobs] [threads] [pattern] [blk_m] [tagx] [extra] [cache_dir]"
site=${1:?$usage}
maxj=${2:-4}
thr=${3:-4}
PAT=${4:-.}
BLK=${5:-2000}
tagx=${6:-}
extra=${7:-}
CACHE=${8:-$FAB/stages/quadpol_cf/coreg_cache}
case $site in
  ridge_a)     ROOT=/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2 ;;
  taylor_dome) ROOT=/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2 ;;
  thwaites)    ROOT=/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground ;;
  eastwind)    ROOT=/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground ;;
  *) echo "unknown site $site"; echo "$usage"; exit 1 ;;
esac
for v in "$maxj" "$thr" "$BLK"; do
  if ! [[ $v =~ ^[1-9][0-9]*$ ]]; then
    echo "node_survey.sh: max_jobs, threads and blk_m must be positive integers (got '$v')"
    echo "$usage"
    exit 1
  fi
done
if [ ! -d "$CACHE" ]; then
  echo "no cache directory $CACHE - build the caches first:"
  echo "  bash $FAB_DIR/coreg_batch.sh $site 3 firn"
  exit 1
fi
run=$tagx
[ "$BLK" = 2000 ] || run=_B$BLK$tagx
M=/opt/sw/matlab/2024b/bin/matlab
if [ ! -x "$M" ]; then
  echo "no MATLAB at $M - this launcher is for the CReSIS machines"
  exit 1
fi
mkdir -p "$FAB/locks" "$FAB/logs" "$FAB/stages/node" || exit 1
slock=$FAB/locks/node_survey_$site$run
mkdir "$slock" 2>/dev/null || { echo "another runner holds $slock"; exit 0; }
trap 'rmdir "$slock" 2>/dev/null' EXIT
cd "$FAB/logs" || exit 1

# the frames of this pass: every cache matching pattern that is one of the
# SITE's frames. The runner takes each trace's position from the season's
# polarimetric product, so a cache of another site's frame in the same
# directory has nothing to be cut along and is left alone.
list=$FAB/logs/node_frames_$site$run.txt
ls "$CACHE" | grep "^creg_[0-9]\{8\}_[0-9][0-9]_[0-9][0-9][0-9]\.mat$" | \
  sed 's/^creg_//; s/\.mat$//' | grep -- "$PAT" | sort | \
while read -r tag; do
  seg=${tag%_*}
  if [ -s "$ROOT/CSARP_polarimetric/$seg/Data_$tag.mat" ] || \
     [ -s "$ROOT/CSARP_polarimetric_unwrap/$seg/Data_$tag.mat" ]; then
    echo "$tag"
  fi
done > "$list"
total=$(wc -l < "$list")
echo "=== $site: $total cached frames match, $maxj at a time, started $(date)"

while read -r tag; do
  seg=${tag%_*}
  fr=$((10#${tag##*_}))
  out=$FAB/stages/node/node_${tag}_B$BLK$tagx.mat
  if [ -s "$out" ]; then
    echo "skip  $tag (product exists)"
    continue
  fi
  if ! mkdir "$FAB/locks/node_$tag$run" 2>/dev/null; then
    echo "skip  $tag (locked)"
    continue
  fi
  while [ "$(jobs -rp | wc -l)" -ge "$maxj" ]; do sleep 30; done
  (
    echo "start $tag $(date '+%H:%M')"
    nice -n 10 "$M" -batch "addpath('$FAB_DIR'); site_root='$ROOT'; day_seg='$seg'; frm=$fr; threads=$thr; tagx='$tagx'; blk_m=$BLK; cache_dir='$CACHE'; $extra run_node_frame" \
      > "node_$tag$run.log" 2>&1 < /dev/null
    rc=$?
    rmdir "$FAB/locks/node_$tag$run"
    if [ -s "$out" ]; then st=ok; else st=FAILED; fi
    echo "done  $tag rc=$rc $st $(date '+%H:%M')"
  ) &
done < "$list"
wait
have=0
while read -r tag; do
  [ -s "$FAB/stages/node/node_${tag}_B$BLK$tagx.mat" ] && have=$((have + 1))
done < "$list"
echo "=== $site finished $(date): $have of $total products present in $FAB/stages/node"

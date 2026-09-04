#!/bin/bash
# Quad-pol fabric for the January 2023 segments of 2022_Antarctica_Ground:
# the WAIS Divide grid and the southwest spur that runs toward Kamb.
#
#   wais_kamb_batch.sh            # DRY=1 to only build and report the work list
#
# WHAT THESE SEGMENTS ARE. Mapped from their own coordinates
# (scripts/figures/kamb_or_wais.py): a local grid whose median distance
# from WAIS Divide camp is 19.8 km, plus an outbound leg - 20230118_11/13,
# then 20230120_04 (29 km) and 20230120_05 (49 km) - reaching 84 km
# southwest before returning on 21 January. The Kamb Ice Stream trunk is
# 518 km from camp and the closest any of this gets to it is 440 km, so
# these are WAIS Divide profiles with a spur in Kamb's direction. Named
# accordingly, and nothing here should be published as "Kamb".
#
# WHY THEY NEED THE NO-POLARIMETRIC PATH. They are SAR focused - all four
# CSARP_standardphase channels exist - but the polarimetric step never ran
# on them: CSARP_polarimetric_unwrap holds a directory per segment and no
# files. The pipeline's no-polarimetric mode supplies what that product
# would have carried (coregistration defaults, a picked surface, a range
# window from z_max), and now prefers the SAR-focused channels over qlook
# when both could exist.
#
# IMAGE 02. Each frame ships as two images, which are two range windows of
# the same traces: img_01 spans about 1100 m of ice and img_02 about
# 4600 m. Only img_02 reaches the depths this work quotes, so it is what
# runs here; the pipeline puts the image in the product and cache names so
# the two can never overwrite each other.
set -u
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
ML=${MATLAB_BIN:-/opt/sw/matlab/2024b/bin/matlab}
ROOT=/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground
IMG=${IMG:-img_02}
# CT=1 runs the CONSTANT-ORIENTATION variant: one fabric axis per ~2 km
# along-track segment, held through depth, with dlam(z) still free. Its
# products carry the '_ct' suffix after the image, so they sit beside the
# per-window ones rather than replacing them and
# scripts/compare_const_theta.py pairs them by stripping that suffix.
# Run it AFTER the default pass: that pass builds the coregistration
# caches, which turn a 50-70 min frame into a few minutes here.
CT=${CT:-}
if [ -n "$CT" ]; then CT_ARG="theta_const=true;"; SUF="_${IMG}_ct"; TAGP="ct"; else CT_ARG=""; SUF="_${IMG}"; TAGP="wk"; fi
mkdir -p "$FAB/invert_logs" "$FAB/stages/quadpol"

# Same shadowing guard as const_theta_batch.sh: MATLAB resolves the current
# folder ahead of the path, so a stale flat copy of run_quadpol_pipeline.m
# in the launch directory would silently win. An old pipeline does not know
# `img` and would look for a plain Data_<seg>_<frm>.mat that does not exist
# in this season - a loud failure rather than a silent one, but the guard
# is cheap and the failure mode is the reason it exists.
cd "$FAB_DIR" || exit 1
GUARD="p=which('run_quadpol_pipeline'); q=fullfile('$FAB_DIR','run_quadpol_pipeline.m'); \
if ~strcmp(p,q), error('pipeline shadowed by %s (expected %s)',p,q); end;"

# Work list: every January 2023 segment with all four channels present for
# this image. Frame numbers come from the files themselves rather than a
# range, because the segments carry 1 to 3 frames each.
: > "$FAB/wk_work.txt"
for s in $(ls "$ROOT/CSARP_standardphase_HH" 2>/dev/null | grep '^2023'); do
  for f in "$ROOT/CSARP_standardphase_HH/$s/Data_${IMG}_${s}_"*.mat; do
    [ -f "$f" ] || continue
    frm=$(basename "$f" | sed -E "s/.*_([0-9]{3})\.mat/\1/")
    ok=1
    for c in VV HV VH; do
      [ -f "$ROOT/CSARP_standardphase_$c/$s/Data_${IMG}_${s}_${frm}.mat" ] || ok=0
    done
    [ "$ok" -eq 1 ] && echo "$s $frm" >> "$FAB/wk_work.txt"
  done
done
echo "$(wc -l < "$FAB/wk_work.txt") frame(s) queued for image $IMG${CT:+, CONSTANT-ORIENTATION variant}"
[ -n "${DRY:-}" ] && exit 0

worker () {
  while read -r seg frm; do
    tag=$(printf '%s_%03d' "$seg" "$((10#$frm))")
    out="$FAB/stages/quadpol/quadpol_section_${tag}${SUF}.mat"
    [ -f "$out" ] && continue
    fail="$FAB/invert_logs/fail_${TAGP}_${tag}.count"
    n=$(cat "$fail" 2>/dev/null || echo 0)
    [ "$n" -ge 2 ] && continue
    lock="$FAB/invert_logs/lock_${TAGP}_${tag}"
    mkdir "$lock" 2>/dev/null || continue
    echo "start $tag $(date +%d/%H:%M)"
    nice -n 10 $ML -batch "addpath('$FAB_DIR'); $GUARD maxNumCompThreads(8); site_root='$ROOT'; day_seg='$seg'; frm=$((10#$frm)); img='$IMG'; $CT_ARG run_quadpol_pipeline" \
      > "$FAB/invert_logs/${TAGP}_${tag}.log" 2>&1
    if [ -f "$out" ]; then
      echo "done  $tag ok $(date +%d/%H:%M)"; rm -f "$fail"
    else
      echo "done  $tag FAIL $(date +%d/%H:%M)"; echo $((n+1)) > "$fail"
    fi
    rmdir "$lock" 2>/dev/null
  done < "$FAB/wk_work.txt"
}

echo "preflight: checking which pipeline MATLAB loads"
$ML -batch "addpath('$FAB_DIR'); $GUARD \
  s=fileread(q); assert(contains(s,'chan_name'), 'pipeline predates image-split support: %s', q); \
  fprintf('preflight OK: %s\n', p);" || { echo "PREFLIGHT FAILED - batch not started"; exit 1; }

# Four workers, not six: these frames are far larger than the Antarctic
# ones this repo usually runs - HH alone is 2.7 GB on 20230120_05 against
# ~700 MB at Ridge A - and coregistration holds all four channels at once.
for w in 1 2 3 4; do worker & done
wait
echo "=== WAIS/Kamb-spur batch finished $(date): \
$(ls "$FAB"/stages/quadpol/quadpol_section_*${SUF}.mat 2>/dev/null | wc -l) products"

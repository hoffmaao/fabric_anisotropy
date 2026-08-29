#!/bin/bash
# EGRIP chain: gate on the revalidated frame's section, then batch every
# complex frame of the season. Launch once with nohup and it runs
# unattended: the gate reads the section left by the revalidation rather
# than recomputing it, it stops a bad estimator from burning days of
# compute, and the batch is idempotent - sections that exist are skipped,
# failures get fail markers capped at 2 retries, mkdir locks arbitrate the
# workers - so after ANY interruption the same launch command resumes
# where it left off. Segment 20240618_01 is excluded: it is a real-only
# setup day (the pipeline would error on it deliberately). Frames missing
# a VH channel (9 in the season) fail loudly into their own logs and are
# capped like any other failure.
set -u
# <work> comes from the walk in fabric_paths.sh, so the chain lands on the
# same root whether it is run from the checkout or from a copy in <work>,
# and on the same root MATLAB derives. FABRIC_ROOT overrides. No username
# is baked in, and nothing has to be deployed anywhere.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
F=$(fabric_work_root "$FAB_DIR") || exit 1
R=/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2
ML=/opt/sw/matlab/2024b/bin/matlab
cd "$F" || { echo "cannot cd to $F" >&2; exit 1; }
mkdir -p invert_logs
echo "=== egrip chain start $(date)"

# --- step 1: gate on the ALREADY-REVALIDATED section (seg_reval.sh, 17 Aug:
# section resid 0.562 -> 0.476 with the segmented first pass, dlam 0.229 =
# the domain median, finite 0.91). The residual floor was diagnosed as the
# 2psi-dominant chan_equal-family calibration signature - the same term
# Ridge A carries at a third the amplitude - so the gate is RE-STATED per
# the author's 18 Aug decision: finite fraction + block self-consistency;
# resid is RECORDED as the documented site floor, not gated on.
# <work> arrives as argv[1] so the heredoc can stay QUOTED: unquoting it to
# interpolate $F would expand every $ in the Python body too.
gate=$(python3 - "$F" <<'PY'
import os, sys
import h5py, numpy as np
try:
    fn = os.path.join(sys.argv[1], "stages", "quadpol",
                      "quadpol_section_20240619_01_001.mat")
    with h5py.File(fn) as f:
        r = f["res"]
        z = np.array(r["z"]).ravel()
        dl = np.array(r["sec_dlam_ls"]); rs = np.array(r["sec_resid_ls"])
    if dl.shape[0] == z.size:
        dl = dl.T; rs = rs.T
    m = (z > 300) & (z < 1100)
    med = np.array([np.nanmedian(dl[b, m]) for b in range(dl.shape[0])])
    fin = float(np.mean(np.isfinite(dl[:, m])))
    resid = float(np.nanmedian(rs[:, m]))
    ok = np.isfinite(med[:-1]) & np.isfinite(med[1:])
    adj = float(np.nanmedian(np.abs(med[1:][ok] - med[:-1][ok])))
    verdict = "PASS" if (fin > 0.60 and adj < 0.05) else "FAIL"
    print("%s finite=%.2f adj_ddlam=%.3f resid_floor=%.3f dlam=%.3f" % (
        verdict, fin, adj, resid, float(np.nanmedian(med))))
except Exception as exc:
    print("FAIL error=%s" % exc)
PY
)
echo "gate: $gate"
case "$gate" in
  PASS*) ;;
  *) echo "=== egrip chain STOPPED at the gate $(date)"; exit 1;;
esac

# --- step 2: the batch. Work list from the HH channel, img files and the
# real-only setup day excluded.
ls $R/CSARP_qlook_HH/*/Data_2024*.mat 2>/dev/null | grep -v Data_img \
  | grep -v 20240618_01 \
  | sed -E 's|.*/Data_([0-9]{8}_[0-9]{2})_([0-9]{3})\.mat|\1 \2|' \
  | sort -u > egrip_work.txt
echo "$(wc -l < egrip_work.txt) frames queued"

# Stale locks from a killed run would make every worker skip those frames
# forever; the chain is launched once, so any lock present now is stale.
stale=$(find "$F/invert_logs" -maxdepth 1 -type d -name 'lock_egrip_*' 2>/dev/null)
if [ -n "$stale" ]; then
  echo "removing $(echo "$stale" | wc -l) stale lock(s)"
  echo "$stale" | xargs rmdir
fi

worker() {
  while read -r seg frm; do
    tag=$(printf '%s_%03d' "$seg" "$((10#$frm))")
    out="$F/stages/quadpol/quadpol_section_${tag}.mat"
    [ -f "$out" ] && continue
    fail="$F/invert_logs/fail_egrip_${tag}.count"
    n=$(cat "$fail" 2>/dev/null || echo 0)
    [ "$n" -ge 2 ] && continue
    lock="$F/invert_logs/lock_egrip_${tag}"
    mkdir "$lock" 2>/dev/null || continue
    echo "start $tag $(date +%d/%H:%M)"
    nice -n 10 $ML -batch "addpath('$FAB_DIR'); maxNumCompThreads(8); site_root='$R'; day_seg='$seg'; frm=$((10#$frm)); run_quadpol_pipeline" \
      > "$F/invert_logs/egrip_${tag}.log" 2>&1
    if [ -f "$out" ]; then
      echo "done  $tag ok $(date +%d/%H:%M)"; rm -f "$fail"
    else
      echo "done  $tag FAIL $(date +%d/%H:%M)"; echo $((n+1)) > "$fail"
    fi
    rmdir "$lock" 2>/dev/null
  done < "$F/egrip_work.txt"
}
for w in 1 2 3; do worker & done
wait
echo "=== egrip chain finished $(date): $(ls "$F"/stages/quadpol/quadpol_section_202406*.mat 2>/dev/null | wc -l) sections on disk"

#!/bin/bash
# EGRIP chain: revalidate one frame with the raised dlam ceiling, gate on
# its residuals, then batch every complex frame of the season. Launch once
# with nohup and it runs unattended: the revalidation reuses the existing
# coregistration cache, the gate stops a bad estimator from burning days of
# compute, and the batch is idempotent - sections that exist are skipped,
# failures get fail markers capped at 2 retries, mkdir locks arbitrate the
# workers - so after ANY interruption the same launch command resumes
# where it left off. Segment 20240618_01 is excluded: it is a real-only
# setup day (the pipeline would error on it deliberately). Frames missing
# a VH channel (9 in the season) fail loudly into their own logs and are
# capped like any other failure.
F=/kucresis/scratch/hoffmana_sta/fabric
R=/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2
ML=/opt/sw/matlab/2024b/bin/matlab
cd $F
mkdir -p invert_logs
echo "=== egrip chain start $(date)"

# --- step 1: revalidation of 20240619_01_001 (coreg cache exists, ~25 min)
$ML -batch "maxNumCompThreads(8); site_root='$R'; day_seg='20240619_01'; frm=1; run_quadpol_pipeline" > egrip_reval.log 2>&1
gate=$(python3 - <<'PY'
import h5py, numpy as np
try:
    fn = "/kucresis/scratch/hoffmana_sta/fabric/stages/quadpol/quadpol_section_20240619_01_001.mat"
    with h5py.File(fn) as f:
        r = f["res"]
        zw = np.array(r["ls_zw"]).ravel()
        rs = np.array(r["ls_resid"]).ravel()
        dl = np.array(r["ls_dlam"]).ravel()
    m = (zw > 300) & (zw < 1100)
    med = float(np.nanmedian(rs[m]))
    fin = float(np.mean(np.isfinite(dl[m])))
    dlm = float(np.nanmedian(dl[m]))
    ok = "PASS" if (med < 0.40 and fin > 0.60) else "FAIL"
    print("%s resid=%.3f finite=%.2f dlam=%.3f" % (ok, med, fin, dlm))
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

worker() {
  while read -r seg frm; do
    tag=$(printf '%s_%03d' "$seg" "$((10#$frm))")
    out=$F/stages/quadpol/quadpol_section_${tag}.mat
    [ -f "$out" ] && continue
    fail=$F/invert_logs/fail_egrip_${tag}.count
    n=$(cat "$fail" 2>/dev/null || echo 0)
    [ "$n" -ge 2 ] && continue
    lock=$F/invert_logs/lock_egrip_${tag}
    mkdir "$lock" 2>/dev/null || continue
    echo "start $tag $(date +%d/%H:%M)"
    nice -n 10 $ML -batch "maxNumCompThreads(8); site_root='$R'; day_seg='$seg'; frm=$((10#$frm)); run_quadpol_pipeline" \
      > $F/invert_logs/egrip_${tag}.log 2>&1
    if [ -f "$out" ]; then
      echo "done  $tag ok $(date +%d/%H:%M)"; rm -f "$fail"
    else
      echo "done  $tag FAIL $(date +%d/%H:%M)"; echo $((n+1)) > "$fail"
    fi
    rmdir "$lock" 2>/dev/null
  done < $F/egrip_work.txt
}
for w in 1 2 3; do worker & done
wait
echo "=== egrip chain finished $(date): $(ls $F/stages/quadpol/quadpol_section_202406*.mat 2>/dev/null | wc -l) sections on disk"

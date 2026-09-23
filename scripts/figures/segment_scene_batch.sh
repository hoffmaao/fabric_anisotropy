#!/bin/bash
# Build a per-segment scene for every mirrored section product.
#
#   segment_scene_batch.sh          # build what is missing
#   DRY=1 segment_scene_batch.sh    # list only; touches nothing
#   FORCE=1 segment_scene_batch.sh  # rebuild even where a figure exists
#
# BUILT TO SURVIVE INTERRUPTS, because this is 146 figures and several of
# them read a BedMachine window off a 848 MB file. Three mechanisms, all of
# them the ones the server batches in opr_fabric/server already use:
#
#   RESUME   a segment whose PNG exists is skipped, so re-running after a
#            kill picks up where it stopped. FORCE=1 overrides.
#   LOCKS    mkdir is atomic, so two copies of this script - or a second
#            run started before the first was noticed - cannot work the
#            same segment. Stale locks are cleared at start only when no
#            other copy is running.
#   FAIL CAP a segment that fails twice is left alone rather than retried
#            forever. The count survives a restart; a success clears it.
#
# Figures land in figs/ (gitignored) and the log in figs/.scene_batch/.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
SCAR="${SCAR_DIR:-$HOME/data/opr/scar}"
OUT="${OUT_DIR:-$REPO/figs}"
STATE="$OUT/.scene_batch"
PY="${PYBIN:-/opt/anaconda3/envs/pygmt/bin/python}"
export PROJ_DATA="${PROJ_DATA:-/opt/anaconda3/share/proj}"
mkdir -p "$STATE"
LOG="$STATE/batch.log"

# One instance at a time. A second invocation reports and exits rather than
# racing the first through the same work list.
if ! mkdir "$STATE/RUNNING" 2>/dev/null; then
  echo "another segment_scene_batch is running (remove $STATE/RUNNING if not)"
  exit 1
fi
trap 'rmdir "$STATE/RUNNING" 2>/dev/null; echo "interrupted; rerun to resume" | tee -a "$LOG"' EXIT INT TERM

# Work list: free-axis products only. The _ct products hold one axis per
# segment and the _z reruns cannot be parsed by the frame regex.
LIST=$(ls "$SCAR"/quadpol_section_*.mat 2>/dev/null \
       | grep -v '_ct\.mat' | grep -v '_z[0-9]' \
       | xargs -n1 basename 2>/dev/null \
       | sed 's/quadpol_section_//; s/\.mat$//' | sort)
TOTAL=$(echo "$LIST" | grep -c . || true)

todo=0
for f in $LIST; do
  [ -z "${FORCE:-}" ] && [ -f "$OUT/segment_scene_$f.png" ] && continue
  n=$(cat "$STATE/fail_$f" 2>/dev/null || echo 0)
  [ "$n" -ge 2 ] && continue
  todo=$((todo+1))
done
echo "$TOTAL segment(s) mirrored; $todo to build" | tee -a "$LOG"
if [ -n "${DRY:-}" ]; then
  echo "DRY: nothing built"
  trap - EXIT INT TERM; rmdir "$STATE/RUNNING" 2>/dev/null; exit 0
fi

ok=0; bad=0; skip=0
for f in $LIST; do
  png="$OUT/segment_scene_$f.png"
  if [ -z "${FORCE:-}" ] && [ -f "$png" ]; then skip=$((skip+1)); continue; fi
  n=$(cat "$STATE/fail_$f" 2>/dev/null || echo 0)
  if [ "$n" -ge 2 ]; then
    echo "skip  $f (failed $n times)" | tee -a "$LOG"; skip=$((skip+1)); continue
  fi
  mkdir "$STATE/lock_$f" 2>/dev/null || { skip=$((skip+1)); continue; }
  if "$PY" "$HERE/segment_scene.py" "$f" "$OUT" >> "$STATE/$f.log" 2>&1 && [ -f "$png" ]; then
    ok=$((ok+1)); rm -f "$STATE/fail_$f"
    printf 'ok    %-22s (%d/%d)\n' "$f" "$((ok+bad))" "$todo" | tee -a "$LOG"
  else
    bad=$((bad+1)); echo $((n+1)) > "$STATE/fail_$f"
    printf 'FAIL  %-22s attempt %d - see %s\n' "$f" "$((n+1))" "$STATE/$f.log" | tee -a "$LOG"
  fi
  rmdir "$STATE/lock_$f" 2>/dev/null
done
trap - EXIT INT TERM
rmdir "$STATE/RUNNING" 2>/dev/null
echo "=== finished $(date): $ok built, $bad failed, $skip skipped" | tee -a "$LOG"

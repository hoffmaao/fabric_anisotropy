#!/bin/bash
# Run additional survey sites end to end (cache build + inversion
# follower), one site at a time, starting only after the currently
# running batches have gone fully quiet:  site_chain.sh <site> [site...]
#
#   nohup bash site_chain.sh taylor_dome eastwind > site_chain.log 2>&1 &
set -u
# The batch scripts come from this script's own directory, <work> from the
# walk in fabric_paths.sh. FABRIC_ROOT overrides <work> only. No username
# is baked in, and nothing has to be deployed anywhere.
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
# the per-site chain logs below are written relative to here
cd "$FAB" || { echo "cannot cd to $FAB" >&2; exit 1; }
echo "waiting for current batches to finish, $(date)"
while pgrep -f "[c]oreg_batch.sh" > /dev/null || \
      pgrep -f "[i]nvert_one.sh" > /dev/null || \
      pgrep -f "[i]nvert_batch.sh" > /dev/null; do
  sleep 600
done
for site in "$@"; do
  echo "=== site $site starting $(date)"
  nohup bash "$FAB_DIR/coreg_batch.sh" "$site" 3 > "coreg_${site}.log" 2>&1 &
  sleep 30
  # the follower runs in the foreground and exits when the site's caches
  # are all swept and its coreg batch is gone; wait then reaps the coreg
  bash "$FAB_DIR/invert_batch.sh" "$site" 2 > "invert_chain_${site}.log" 2>&1
  wait
  echo "=== site $site done $(date)"
done
echo "SITE CHAIN DONE $(date)"

#!/bin/bash
# Run additional survey sites end to end (cache build + inversion
# follower), one site at a time, starting only after the currently
# running batches have gone fully quiet:  site_chain.sh <site> [site...]
#
#   nohup bash site_chain.sh taylor_dome eastwind > site_chain.log 2>&1 &
set -u
FAB=/kucresis/scratch/hoffmana_sta/fabric
cd "$FAB"
echo "waiting for current batches to finish, $(date)"
while pgrep -f "[c]oreg_batch.sh" > /dev/null || \
      pgrep -f "[i]nvert_one.sh" > /dev/null || \
      pgrep -f "[i]nvert_batch.sh" > /dev/null; do
  sleep 600
done
for site in "$@"; do
  echo "=== site $site starting $(date)"
  nohup bash coreg_batch.sh "$site" 3 > "coreg_${site}.log" 2>&1 &
  sleep 30
  # the follower runs in the foreground and exits when the site's caches
  # are all swept and its coreg batch is gone; wait then reaps the coreg
  bash invert_batch.sh "$site" 2 > "invert_chain_${site}.log" 2>&1
  wait
  echo "=== site $site done $(date)"
done
echo "SITE CHAIN DONE $(date)"

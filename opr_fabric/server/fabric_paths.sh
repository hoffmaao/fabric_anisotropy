#!/bin/bash
# Work root for the shell launchers, derived. SOURCE this, do not run it:
#
#   FAB_DIR=$(cd "$(dirname "$0")" && pwd)
#   . "$FAB_DIR/fabric_paths.sh"
#   FAB=$(fabric_work_root "$FAB_DIR") || exit 1
#
# The shell mirror of fabric_paths.m, and it keeps the SAME split, so one
# name cannot mean two things across the two languages:
#
#   CODE  location comes from the file. A launcher's workers are its
#         siblings in the repo and stay siblings wherever the tree is
#         copied, so they are always read from "$FAB_DIR", never from the
#         work root, and FABRIC_ROOT does not move them.
#   WORK  root holds stages/ and the logs, and is overridable with
#         FABRIC_ROOT for a layout whose products do not sit beside the
#         code - exactly as it overrides work_root in fabric_paths.m.
#
# Conflating the two meant exporting FABRIC_ROOT for its documented purpose
# left every launcher unable to find its own worker.
#
# The walk looks for stages/ rather than counting directories, mirroring
# the +ptt walk in fabric_paths.m, so the launchers work unmoved from the
# checkout (<work>/code/opr_fabric/server) and from a copy in <work>, and
# no deploy step has to be remembered. It FAILS rather than falling back on
# the script's own directory: the silent failure that guard exists to stop
# is a run that takes its locks and its skip-if-cached checks in a tree
# that holds no products, finds nothing to skip, and recomputes every frame
# at ~50 min each while MATLAB writes to the real work root.
fabric_work_root() {
  local d parent k
  if [ -n "${FABRIC_ROOT:-}" ]; then
    printf '%s\n' "$FABRIC_ROOT"
    return 0
  fi
  d=$1
  for k in 1 2 3 4 5; do
    if [ -d "$d/stages" ]; then
      printf '%s\n' "$d"
      return 0
    fi
    parent=$(dirname "$d")
    if [ -z "$parent" ] || [ "$parent" = "$d" ]; then
      break
    fi
    d=$parent
  done
  echo "no stages/ directory above $1 - the launchers must sit inside a" \
    "work root (mkdir -p <work>/stages/quadpol), or set FABRIC_ROOT to" \
    "one" >&2
  return 1
}

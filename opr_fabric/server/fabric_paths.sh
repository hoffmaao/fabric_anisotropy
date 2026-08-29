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
# The work root is DERIVED FROM THE CHECKOUT, on the same +ptt anchor
# fabric_paths.m uses, so the two cannot disagree about which tree a run
# belongs to: the repo root is the directory holding +ptt, and the work
# root is its parent. That yields ONE candidate, never a list of them, so
# there is nothing to fall through to. A launcher copied out of the repo
# (into <work> beside stages/, the old deploy layout) has no +ptt above it
# and is its own work root.
#
# stages/ is then CONFIRMATION on that one candidate, and its absence is
# fatal. Walking on to an ancestor that happens to have stages/ would be
# worse than not resolving at all: a second checkout beside a live one -
# which this layout exists to allow - would take its mkdir locks and write
# its fail markers in the LIVE tree while MATLAB, anchored on +ptt, wrote
# products into the test tree. So the walk stops at the checkout boundary
# and says what it looked for.
#
# This is the one place the two files deliberately differ: MATLAB creates
# stages/ under its own output root, so fabric_paths.m does not require it
# to pre-exist. The launchers do, because their skip-if-cached checks and
# mkdir locks all live under stages/, and against a missing one they find
# nothing to skip, arbitrate nothing, and recompute every frame at ~50 min
# each. Both sides agree on the rule that matters: code from the file,
# product root overridable and validated.
fabric_work_root() {
  local start repo d parent k work

  start=$1
  if [ -n "${FABRIC_ROOT:-}" ]; then
    if [ ! -d "$FABRIC_ROOT" ]; then
      echo "FABRIC_ROOT is $FABRIC_ROOT, which is not a directory - an" \
        "override that does not exist would put every lock and every" \
        "skip-if-cached check in a tree holding no products" >&2
      return 1
    fi
    # absolute, because site_chain.sh and egrip_chain.sh cd away from here
    ( cd "$FABRIC_ROOT" && pwd )
    return 0
  fi

  repo=
  d=$start
  for k in 1 2 3 4 5; do
    if [ -d "$d/+ptt" ]; then
      repo=$d
      break
    fi
    parent=$(dirname "$d")
    if [ -z "$parent" ] || [ "$parent" = "$d" ]; then
      break
    fi
    d=$parent
  done

  if [ -n "$repo" ]; then
    work=$(dirname "$repo")
  else
    work=$start
  fi

  if [ ! -d "$work/stages" ]; then
    echo "no stages/ in $work, the work root for $start - create it" \
      "(mkdir -p $work/stages/quadpol) or set FABRIC_ROOT to a work root." \
      "Refusing to search upwards: an ancestor's stages/ belongs to" \
      "another run." >&2
    return 1
  fi
  printf '%s\n' "$work"
}

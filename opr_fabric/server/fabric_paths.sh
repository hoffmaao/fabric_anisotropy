#!/bin/bash
# Work root for the shell launchers, derived. SOURCE this, do not run it:
#
#   FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
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
# A WORK ROOT IS A WORK ROOT BY THE SAME TEST HOWEVER IT WAS OBTAINED.
# FABRIC_ROOT picks the candidate; it does not exempt it. Both candidates
# then go through the one block below - exists, can be entered, is
# absolute, holds stages/. Why an override is held to the derived value's
# test rather than a weaker one is an argument about both helpers, so it is
# owned by docs/scripts.md, "Running on the CReSIS machines".
#
# stages/ is the confirmation, and its absence is fatal either way. Walking
# on to an ancestor that happens to have stages/ would be worse than not
# resolving at all: a second checkout beside a live one - which this layout
# exists to allow - would take its mkdir locks and write its fail markers
# in the LIVE tree while MATLAB, anchored on +ptt, wrote products into the
# test tree. So nothing here searches upwards for stages/, and every
# failure says what it looked for and where.
#
# Every failure path also RETURNS NON-ZERO. Callers guard with
# `FAB=$(fabric_work_root ...) || exit 1`, and a function that printed
# nothing while reporting success would set FAB to the empty string, which
# set -u does not catch: every path would then resolve against /, the
# skip-if-cached checks would miss every frame, and the survey would
# recompute at ~50 min each.
#
# Everything else agrees with fabric_paths.m: code from the file, product
# root overridable, and the override validated by the same test as the
# derived value. The two files diverge in exactly two ways - the stages/
# requirement above, which MATLAB does not make, and symlinks - and both
# are claims about the pair, so docs/scripts.md owns the reasoning for each
# and neither header restates it.
#
# What symlinks mean HERE: bash's pwd is LOGICAL by default, so a work root
# reached through one comes back as the path written, where MATLAB returns
# the resolved target. Use `pwd -P` in the resolve below if identical
# strings are ever wanted.
fabric_work_root() {
  local start repo d parent k work abs origin

  start=$1
  if [ -n "${FABRIC_ROOT:-}" ]; then
    work=$FABRIC_ROOT
    origin="FABRIC_ROOT"
  else
    repo=
    d=$start
    for k in 1 2 3 4 5; do
      if [ -d "$d/+ptt" ]; then
        repo=$d
        break
      fi
      parent=$(dirname "$d") || return 1
      if [ -z "$parent" ] || [ "$parent" = "$d" ]; then
        break
      fi
      d=$parent
    done
    if [ -n "$repo" ]; then
      work=$(dirname "$repo") || return 1
      origin="the work root for $start, the parent of the $repo checkout"
    else
      work=$start
      origin="the work root for $start, which has no +ptt checkout above it"
    fi
  fi

  if [ -z "$work" ]; then
    echo "$origin resolved to an empty path" >&2
    return 1
  fi
  if [ ! -d "$work" ]; then
    echo "$origin is $work, which is not a directory - a root that does" \
      "not exist would put every lock and every skip-if-cached check in a" \
      "tree holding no products, so every frame would recompute" >&2
    return 1
  fi
  # absolute, because site_chain.sh and egrip_chain.sh cd away from here.
  # Captured separately: assigning straight into work would blank it before
  # the failure message could name it.
  abs=$( cd "$work" && pwd ) || {
    echo "$origin is $work, which cannot be entered - check permissions" >&2
    return 1
  }
  if [ -z "$abs" ]; then
    echo "$origin is $work, which did not resolve to an absolute path" >&2
    return 1
  fi
  work=$abs
  if [ ! -d "$work/stages" ]; then
    echo "no stages/ in $work - $origin. Create it (mkdir -p" \
      "$work/stages/quadpol) or point FABRIC_ROOT at a work root that has" \
      "one. Not searching upwards: an ancestor's stages/ belongs to" \
      "another run." >&2
    return 1
  fi
  printf '%s\n' "$work"
}

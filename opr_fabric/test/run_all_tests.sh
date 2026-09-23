#!/bin/bash
# Run the opr_fabric MATLAB test suite from a shell, exiting non-zero if any
# test fails. Arguments pass through to run_all_tests.m:
#
#   bash opr_fabric/test/run_all_tests.sh             # every test
#   bash opr_fabric/test/run_all_tests.sh quick       # a few minutes
#   bash opr_fabric/test/run_all_tests.sh test_fabric_task test_quadpol
#
# MATLAB is $MATLAB_BIN if set, else `matlab` on PATH, else the newest
# install in the usual places: /Applications on a Mac, /opt/sw/matlab on the
# CReSIS machines. MATLAB_THREADS=N caps its threads. Uncapped, MATLAB takes
# every core it sees, which on a 112-core CReSIS node was ~65 cores for one
# test, so on a shared node run it as
#
#   MATLAB_THREADS=8 nice bash opr_fabric/test/run_all_tests.sh quick
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)

ml=${MATLAB_BIN:-}
if [ -z "$ml" ]; then
  ml=$(command -v matlab 2>/dev/null || true)
fi
if [ -z "$ml" ]; then
  ml=$( (ls -d /Applications/MATLAB_R*.app/bin/matlab \
    /opt/sw/matlab/*/bin/matlab 2>/dev/null || true) | sort | tail -n 1)
fi
if [ -z "$ml" ] || [ ! -x "$ml" ]; then
  echo "run_all_tests.sh: no MATLAB found; set MATLAB_BIN" >&2
  exit 2
fi

if [ $# -eq 0 ]; then
  call="run_all_tests"
elif [ $# -eq 1 ]; then
  call="run_all_tests('$1')"
else
  list=$(printf "'%s'," "$@")
  call="run_all_tests({${list%,}})"
fi
cap=
if [ -n "${MATLAB_THREADS:-}" ]; then
  if ! [[ $MATLAB_THREADS =~ ^[1-9][0-9]*$ ]]; then
    echo "run_all_tests.sh: MATLAB_THREADS must be a positive integer" >&2
    exit 2
  fi
  cap="maxNumCompThreads($MATLAB_THREADS); "
fi
echo "run_all_tests.sh: $ml -batch \"$cap$call\"" >&2
exec "$ml" -batch "${cap}addpath('$here'); $call"

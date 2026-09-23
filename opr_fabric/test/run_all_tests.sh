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
# CReSIS machines. On a shared node, run it under nice.
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
  echo "run_all_tests.sh: no MATLAB found; set MATLAB_BIN to its bin/matlab" >&2
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
echo "run_all_tests.sh: $ml -batch \"$call\"" >&2
exec "$ml" -batch "addpath('$here'); $call"

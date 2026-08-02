#!/usr/bin/env bash
# Download one quad-pol frame from the public OPR data portal into a local
# directory layout that mirrors the CReSIS servers (CSARP_<product>/<seg>/).
#
# Usage: fetch_frame.sh [day_seg] [frame] [dest_root]
# Defaults: 20250108_02, 9, ~/data/opr/accum/2024_Antarctica_Ground2

set -euo pipefail

DAY_SEG="${1:-20250108_02}"
FRM="${2:-9}"
SEASON="2024_Antarctica_Ground2"
DEST_ROOT="${3:-$HOME/data/opr/accum/$SEASON}"
BASE="https://data.cresis.ku.edu/data/accum/$SEASON"

FN=$(printf "Data_%s_%03d.mat" "$DAY_SEG" "$FRM")

for pol in HH VV HV VH; do
  dir="$DEST_ROOT/CSARP_standard_${pol}/$DAY_SEG"
  mkdir -p "$dir"
  echo "Fetching $pol: $FN"
  curl -sS --fail -o "$dir/$FN" -C - "$BASE/CSARP_standard_${pol}/$DAY_SEG/$FN"
done

echo "Done. Files under $DEST_ROOT"

#!/bin/bash
# Rebuild EVERY constant-orientation product with the currently deployed
# package, from the coregistration caches.
#
#   DRY=1 regen_const_theta.sh   # list what would move and run; touches nothing
#   regen_const_theta.sh         # move products aside, then rebuild
#
# WHY. The 182 constant-orientation products on disk were built with the
# pre-review +ptt package. The review found two independent defects in the
# orientation error bars - replicates gated out for an artefact of grid
# width, and a delete-one standard error that is unbounded on a circular
# quantity. NOTE: the original justification here - that 24% of segments
# on disk reported a standard error above 180 degrees - was a units error
# in the diagnostic (rad2deg applied to a field stored in degrees). Read
# correctly the old error bars have a median of 1.6 deg and are healthy,
# and the point estimates are identical between packages (measured on
# 20250108_02_001: dlam 0.0609 both, axis 20.9 deg both). This rebuild is
# therefore OPTIONAL - it only adds the newer diagnostic fields
# (se_theta_sat, se_theta_r) - and costs about 44 hours because the
# jackknife runs 11 replicates per segment. Not launched.
#
# The point estimates should NOT move: the review was directed to leave
# the full-grid theta_0 and dlam path unchanged. Whether it did is
# MEASURED before this is run (pkgcmp on 20250108_02_001), not assumed.
#
# NOTHING IS DELETED. Every existing product is moved to
# stages/quadpol/pre_fix_products/ first, keeping the record of what the
# old package produced. const_theta_batch.sh is idempotent and skips
# products that exist, so moving them aside is what makes it rebuild.
set -u
FAB_DIR=$(cd "$(dirname "$0")" && pwd) || exit 1
. "$FAB_DIR/fabric_paths.sh"
FAB=$(fabric_work_root "$FAB_DIR") || exit 1
QP="$FAB/stages/quadpol"
KEEP="$QP/pre_fix_products"

# Build the list WITHOUT touching anything, so DRY really is dry. A dry run
# that moves files performed the one irreversible thing it had, and that
# happened once already in egrip_reposition.sh.
LIST=$(ls "$QP"/quadpol_section_*_ct.mat 2>/dev/null | grep -vE '_z[0-9]+' )
N=$(echo "$LIST" | grep -c . || true)
echo "$N constant-orientation product(s) would be moved to pre_fix_products/ and rebuilt"
echo "by tag family:"
echo "$LIST" | sed -E 's|.*/quadpol_section_([0-9]{6}).*|\1|' | sort | uniq -c | sed 's/^/   /'
[ -n "${DRY:-}" ] && { echo "DRY: nothing moved, nothing run"; exit 0; }

mkdir -p "$KEEP"
echo "$LIST" | while read -r f; do [ -f "$f" ] && mv "$f" "$KEEP/"; done
echo "moved $(ls "$KEEP" | wc -l) product(s) aside; launching const_theta_batch.sh"
exec bash "$FAB_DIR/const_theta_batch.sh"

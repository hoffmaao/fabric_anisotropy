"""The one-pair-per-crossing invariant of crossing_pairs.pair_families.

A crossing is a PLACE, so the dedup has to hold no matter how the frame
boundaries fall across it. The two cases that broke when the dedup was
scoped to a frame-pair key are the ones worth pinning, because neither is
visible in the Ridge A totals - they are rare enough to shift n by one or
two without changing the median:

  split     one N-S line cut into consecutive frames A1 and A2 with the
            cut inside the crossing corridor of row line B; the old
            per-key dedup returned one pair under (A1, B) and another
            under (A2, B) for a single crossing
  curved    frames C and D each contribute blocks to BOTH families near
            one crossing, as a curved connector does; the old dedup put
            those candidates under the distinct keys (C, D) and (D, C)

and the case the fix must NOT break, since global clustering is only safe
while the radius stays well under the line spacing:

  distinct  two genuinely separate crossings 3 km apart stay two pairs

Synthetic geometry in EPSG:3031 metres, built from the module's own
MAX_SEP and FAMS so it tracks the defaults rather than hard-coding them.
No data files are read.

Run: python scripts/figures/test_crossing_pairs.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crossing_pairs import FAMS, MAX_SEP, pair_families  # noqa: E402

STEP = 0.25 * MAX_SEP         # block spacing well inside the pairing radius
ROW_AZ = np.radians(FAMS['row'][0])
ROW_DIR = (np.sin(ROW_AZ), np.cos(ROW_AZ))


def ns_blocks(tag, cx, cy, ts):
    """Blocks along a N-S line through (cx, cy), offsets ts."""
    return [(tag, cx, cy + t, 'ns') for t in ts]


def row_blocks(tag, cx, cy, ts):
    """Blocks along a row line through (cx, cy), offsets ts."""
    return [(tag, cx + t * ROW_DIR[0], cy + t * ROW_DIR[1], 'row')
            for t in ts]


def npairs(blocks):
    tags = np.array([b[0] for b in blocks])
    x = np.array([b[1] for b in blocks], float)
    y = np.array([b[2] for b in blocks], float)
    fam = np.array([b[3] for b in blocks], dtype=object)
    dl = np.arange(len(blocks), dtype=float) / 100.0
    pairs, diffs = pair_families(tags, x, y, fam, dl, 'ns', 'row')
    assert len(pairs) == diffs.size, 'pairs and diffs disagree'
    return len(pairs)


def check(name, blocks, want):
    got = npairs(blocks)
    print('%-9s %d pair(s), expected %d   %s'
          % (name, got, want, 'ok' if got == want else 'FAIL'))
    return got == want


def main():
    split = (ns_blocks('A1', 0.0, 0.0, [-3 * STEP, -2 * STEP, -STEP])
             + ns_blocks('A2', 0.0, 0.0, [0.0, STEP, 2 * STEP])
             + row_blocks('B', 0.0, 0.0, [-STEP, 0.0, STEP]))

    # C and D both turn through the crossing, so each holds blocks of both
    # families near it - the geometry that produced keys (C, D) and (D, C).
    curved = [('C', 0.0, 0.0, 'ns'),
              ('C', 0.6 * STEP, 0.1 * STEP, 'row'),
              ('D', 0.0, 0.8 * STEP, 'row'),
              ('D', 0.4 * STEP, 0.4 * STEP, 'ns')]

    distinct = (ns_blocks('A', 0.0, 0.0, [0.0])
                + row_blocks('B1', 0.0, 0.0, [0.0])
                + ns_blocks('A', 0.0, 3000.0, [0.0])
                + row_blocks('B2', 0.0, 3000.0, [0.0]))

    ok = [check('split', split, 1),
          check('curved', curved, 1),
          check('distinct', distinct, 2)]
    if not all(ok):
        raise SystemExit('crossing_pairs dedup invariant violated')
    print('all ok')


if __name__ == '__main__':
    main()

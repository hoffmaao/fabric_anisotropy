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
  repeat    one line re-flown on another day under a second tag, which is
            what actually accounts for the 45 -> 28 drop on the real
            sections and which no frame-keyed dedup can collapse

and the two cases the fix must NOT break, since global clustering is only
safe while the radius stays well under the line spacing:

  distinct  two genuinely separate crossings 3 km apart stay two pairs
  audit     a cluster that DOES fuse two parallel lines is reported by
            merge_audit rather than passing silently
  abstain   a cluster nothing could measure reads UNVERIFIED, not clean,
            since those two must never share a value

Synthetic geometry in EPSG:3031 metres, built from the module's own
defaults and FAMS so it tracks them rather than hard-coding them. No data
files are read, so this runs anywhere the imports resolve.

Run: python scripts/figures/test_crossing_pairs.py   (or under pytest)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crossing_pairs import (DEFAULT_MAX_SEP, FAMS,  # noqa: E402
                            analyse, line_spread, pair_families)

STEP = 0.25 * DEFAULT_MAX_SEP   # block spacing well inside the radius
ROW_AZ = np.radians(FAMS['row'][0])
ROW_DIR = (np.sin(ROW_AZ), np.cos(ROW_AZ))


def ns_blocks(tag, cx, cy, ts):
    """Blocks along a N-S line through (cx, cy), offsets ts."""
    return [(tag, cx, cy + t, 'ns') for t in ts]


def row_blocks(tag, cx, cy, ts):
    """Blocks along a row line through (cx, cy), offsets ts."""
    return [(tag, cx + t * ROW_DIR[0], cy + t * ROW_DIR[1], 'row')
            for t in ts]


def arrays(blocks):
    """Blocks are given in along-track order, so bidx is that order."""
    tags = np.array([b[0] for b in blocks])
    seen = {}
    bidx = []
    for t in tags:
        seen[t] = seen.get(t, -1) + 1
        bidx.append(seen[t])
    x = np.array([b[1] for b in blocks], float)
    y = np.array([b[2] for b in blocks], float)
    fam = np.array([b[3] for b in blocks], dtype=object)
    dl = np.arange(len(blocks), dtype=float) / 100.0
    return tags, np.array(bidx), x, y, fam, dl


def npairs(blocks):
    tags, _, x, y, fam, dl = arrays(blocks)
    pairs, diffs = pair_families(tags, x, y, fam, dl, 'ns', 'row')
    assert len(pairs) == diffs.size, 'pairs and diffs disagree'
    return len(pairs)


def test_split_frames_are_one_crossing():
    blocks = (ns_blocks('A1', 0.0, 0.0, [-3 * STEP, -2 * STEP, -STEP])
              + ns_blocks('A2', 0.0, 0.0, [0.0, STEP, 2 * STEP])
              + row_blocks('B', 0.0, 0.0, [-STEP, 0.0, STEP]))
    assert npairs(blocks) == 1


def test_curved_connector_is_one_crossing():
    # C and D both turn through the crossing, so each holds blocks of both
    # families near it - the geometry that produced keys (C, D) and (D, C).
    blocks = [('C', 0.0, 0.0, 'ns'),
              ('C', 0.6 * STEP, 0.1 * STEP, 'row'),
              ('D', 0.0, 0.8 * STEP, 'row'),
              ('D', 0.4 * STEP, 0.4 * STEP, 'ns')]
    assert npairs(blocks) == 1


def test_repeat_pass_is_one_crossing():
    # The same two lines flown again on another day, 30 m off the first
    # pass. Four frame tags, four ordered keys, one place.
    blocks = (ns_blocks('A', 0.0, 0.0, [-STEP, 0.0, STEP])
              + row_blocks('B', 0.0, 0.0, [-STEP, 0.0, STEP])
              + ns_blocks('A2', 30.0, 0.0, [-STEP, 0.0, STEP])
              + row_blocks('B2', 0.0, 30.0, [-STEP, 0.0, STEP]))
    assert npairs(blocks) == 1


def test_distinct_crossings_stay_apart():
    blocks = (ns_blocks('A', 0.0, 0.0, [0.0])
              + row_blocks('B1', 0.0, 0.0, [0.0])
              + ns_blocks('A', 0.0, 3000.0, [0.0])
              + row_blocks('B2', 0.0, 3000.0, [0.0]))
    assert npairs(blocks) == 2


def audit(blocks, fa='ns', fb='row'):
    tags, bidx, x, y, fam, dl = arrays(blocks)
    return analyse(tags, bidx, x, y, fam, dl, fa, fb)


def test_audit_passes_a_clean_repeat():
    blocks = (ns_blocks('A', 0.0, 0.0, [-STEP, 0.0, STEP])
              + row_blocks('B', 0.0, 0.0, [-STEP, 0.0, STEP])
              + ns_blocks('A2', 30.0, 0.0, [-STEP, 0.0, STEP]))
    _, _, au = audit(blocks)
    assert not au['bad'], au['bad']
    assert au['on_line'] < 150.0, au['on_line']
    assert au['checked'] >= 1 and au['unverified'] == 0, au
    assert au['radius'] > 0


def test_audit_catches_two_parallel_lines_fused():
    # Two N-S lines 400 m apart, both crossing one row line: two DISTINCT
    # crossings 407 m apart, inside the 458 m radius. One cluster swallows
    # both, which is the loss the audit exists to name - and note the
    # count alone would not reveal it, since greedy clustering leaves a
    # straggler behind that keeps the reported n at 2.
    ts = [t * STEP for t in (-3, -2, -1, 0, 1, 2, 3)]
    blocks = (ns_blocks('A', -200.0, 0.0, ts)
              + ns_blocks('C', 200.0, 0.0, ts)
              + row_blocks('B', 0.0, 0.0, [t * STEP for t in range(-4, 5)]))
    _, _, au = audit(blocks)
    assert au['bad'], 'audit missed two parallel lines fused into a cluster'
    assert ['A', 'C'] in [b[0] for b in au['bad']], au['bad']
    assert au['on_line'] > 150.0, au['on_line']
    # and the candidate-set margin collapses too, unlike a separation
    # measured among survivors, which the radius would still bound below
    assert au['gap'] < 150.0, au['gap']


def test_unmeasurable_cluster_abstains_not_passes():
    # Every block isolated in its own frame, as MIN_CELLS holes in the
    # bidx run leave them: no index-adjacent sibling anywhere, so no
    # direction can be taken. This must read UNVERIFIED, never as the
    # clean 0.0 it would share a value with.
    blocks = [('A%d' % k, 0.0, t, 'ns')
              for k, t in enumerate([-STEP, 0.0, STEP])]
    blocks += [('B%d' % k, t * ROW_DIR[0], t * ROW_DIR[1], 'row')
               for k, t in enumerate([-STEP, 0.0, STEP])]
    tags, bidx, x, y, fam, dl = arrays(blocks)
    assert np.isnan(line_spread(tags, bidx, x, y, [0, 1, 2])), 'not NaN'
    _, _, au = audit(blocks)
    assert au['unverified'] >= 1, au
    assert au['checked'] == 0, au
    assert np.isnan(au['on_line']), au['on_line']
    assert not au['bad'], au['bad']


CHECKS = [test_split_frames_are_one_crossing,
          test_curved_connector_is_one_crossing,
          test_repeat_pass_is_one_crossing,
          test_distinct_crossings_stay_apart,
          test_audit_passes_a_clean_repeat,
          test_audit_catches_two_parallel_lines_fused,
          test_unmeasurable_cluster_abstains_not_passes]


def main():
    bad = 0
    for fn in CHECKS:
        try:
            fn()
            print('%-42s ok' % fn.__name__)
        except AssertionError as e:
            bad += 1
            print('%-42s FAIL  %s' % (fn.__name__, e))
    if bad:
        raise SystemExit('%d of %d checks failed' % (bad, len(CHECKS)))
    print('all %d ok' % len(CHECKS))


if __name__ == '__main__':
    main()

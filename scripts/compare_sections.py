"""Block-for-block comparison of two generations of quadpol sections.

Written to measure what the 0/180 heading-wrap interpolation fix
(0e793d5) did to the products: the bug lived on the sub-block rotations
of the N-S family, which is exactly the family the crossing-pair test
finds a deep-dlam bias on, so it had to be excluded as the cause before
chan_equal could be built against that bias. It is general, though - any
two generations of `quadpol_section_*.mat` for the same frames can be
differenced this way.

Per block it reports the change in the deep principal contrast, the
median of `sec_dlam_ls` over 1150-1500 m with the same >= 5 finite cells
requirement `scripts/figures/crossing_pairs.py` uses, so the two agree on
what a usable block is. Per frame it also reports the change in the LS
frame axis theta0, as the doubled-angle phasor difference of the
circular mean of `ls_theta0_geo` over the 200-1200 m band the pipeline
itself summarises the frame over - an axis is never unwrapped.

Deltas are AFTER minus BEFORE, i.e. second directory minus first. Blocks
are matched by along-track index, which is what the same frame reswept
with the same block length gives; a frame whose block count differs is
named and skipped rather than silently mis-aligned, as is one present in
only one directory. z_max special runs (_z<N> tags) are skipped: they
duplicate batch frames.

The delta sample covers only blocks usable in BOTH generations: a block
that clears the finite-cell gate in one directory and not the other has
no delta, so it cannot enter n. Those usability flips are reported
beside n instead of being dropped silently - the `gain` column counts
NaN -> finite and `lost` counts finite -> NaN, per frame and pooled.
The conclusion this tool is used to draw is a null result, and a block
that gained or lost a usable deep contrast is exactly the change a null
must not be able to hide, so read n, gain and lost together.

Reproducing the recorded heading-wrap measurement needs BOTH product
generations on disk. The pre-fix copies are not in this repo and are not
regenerable from it - they survive only in a local mirror captured
before the resweep - so that measurement is a one-time record, not
something CI can re-derive.

Usage: python compare_sections.py <before_dir> <after_dir> [prefix]
"""
import glob
import os
import re
import sys

import h5py
import numpy as np

Z_DEEP = (1150.0, 1500.0)     # deep contrast band, as in crossing_pairs
MIN_CELLS = 5                 # finite section cells required in the band
Z_FRAME = (200.0, 1200.0)     # band the pipeline's frame theta0 uses
BIG = (0.001, 0.005)          # delta thresholds the fractions count

HDR = ('%-15s%5s%6s%6s%9s%9s%8s%8s%7s%6s%6s%6s'
       % ('frame', 'n', 'gain', 'lost', 'median', 'mean', 'p10', 'p90',
          'max', '>1e-3', '>5e-3', 'dth0'))
ROW = '%-15s%5d%6d%6d%+9.5f%+9.5f%+8.4f%+8.4f%7.4f%6.3f%6.3f%s'


def frame_files(d, prefix):
    """Tag -> path for every batch section product under d."""
    out = {}
    pat = os.path.join(d, 'quadpol_section_%s*.mat' % prefix)
    for fn in sorted(glob.glob(pat)):
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        if re.search(r'_z\d+$', tag):
            continue
        out[tag] = fn
    if not out:
        raise SystemExit('no quadpol_section_%s*.mat under %s (z_max '
                         'special runs are skipped)' % (prefix, d))
    return out


def deep_dlam(fn):
    """Per-block deep contrast, NaN where the band is too sparse."""
    with h5py.File(fn) as f:
        r = f['res']
        z = np.array(r['z']).ravel()
        sec = np.array(r['sec_dlam_ls']).T
        zw = np.array(r['ls_zw']).ravel()
        th = (np.array(r['ls_theta0_geo']).ravel()
              if 'ls_theta0_geo' in r else np.array([np.nan]))
    mz = (z > Z_DEEP[0]) & (z < Z_DEEP[1])
    out = np.full(sec.shape[1], np.nan)
    for b in range(sec.shape[1]):
        col = sec[mz, b]
        if np.isfinite(col).sum() >= MIN_CELLS:
            out[b] = np.nanmedian(col)
    mw = (zw > Z_FRAME[0]) & (zw < Z_FRAME[1]) & np.isfinite(th)
    ph = np.mean(np.exp(2j * np.radians(th[mw]))) if mw.any() else np.nan
    return out, ph


def axis_delta(pa, pb):
    """Axis change [deg] from two doubled-angle phasors, a minus b."""
    if not (np.isfinite(pa) and np.isfinite(pb)):
        return np.nan
    return np.degrees(np.angle(pa * np.conj(pb))) / 2.0


def flips(da, db):
    """Usability gained and lost between two per-block samples."""
    ga, gb = np.isfinite(da), np.isfinite(db)
    return int((~ga & gb).sum()), int((ga & ~gb).sum())


def report(name, d, dth, gain, lost):
    """One table row for a delta sample, flip counts beside n."""
    if d.size == 0:
        print('%-15s%5d%6d%6d  no block with a usable deep contrast in both'
              % (name, 0, gain, lost))
        return
    print(ROW % (name, d.size, gain, lost, np.median(d), d.mean(),
                 np.percentile(d, 10), np.percentile(d, 90),
                 np.abs(d).max(),
                 np.mean(np.abs(d) > BIG[0]), np.mean(np.abs(d) > BIG[1]),
                 '   n/a' if not np.isfinite(dth) else '%+6.2f' % dth))


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    before, after = sys.argv[1], sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) > 3 else ''
    fa = frame_files(before, prefix)
    fb = frame_files(after, prefix)

    both = sorted(set(fa) & set(fb))
    print(HDR)
    pooled = []
    skipped = []
    tot_gain = tot_lost = 0
    for tag in both:
        da, pa = deep_dlam(fa[tag])
        db, pb = deep_dlam(fb[tag])
        if da.size != db.size:
            skipped.append('%s (%d vs %d blocks)' % (tag, da.size, db.size))
            continue
        gain, lost = flips(da, db)
        tot_gain += gain
        tot_lost += lost
        d = db - da
        d = d[np.isfinite(d)]
        pooled.append(d)
        report(tag, d, axis_delta(pb, pa), gain, lost)
    for d, only in ((before, sorted(set(fa) - set(fb))),
                    (after, sorted(set(fb) - set(fa)))):
        if only:
            print('only in %s: %s' % (d, ', '.join(only)))
    if skipped:
        print('block count differs, not compared: %s' % ', '.join(skipped))
    if not pooled:
        raise SystemExit('no frame comparable between %s and %s'
                         % (before, after))
    print(HDR)
    report('POOLED', np.concatenate(pooled), np.nan, tot_gain, tot_lost)
    print('  gain/lost are blocks whose usability FLIPPED between the two '
          'generations\n  (NaN -> finite, finite -> NaN). They have no delta, '
          'so they are outside n.')


if __name__ == '__main__':
    main()

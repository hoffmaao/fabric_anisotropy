"""Same-ice crossing-pair test of the heading-family dlam systematic.

The decisive test for a family-dependent bias in the deep principal
contrast: global comparisons across heading families are confounded with
the real NW-SE spatial gradient at Ridge A, but where two survey lines
from different families CROSS, their nearest along-track blocks see the
same ice, so the paired difference isolates whatever the acquisition
geometry adds. Measured 10 Aug 2026 on the pre-wrap-fix sections: N-S
minus row = +0.009 median over 67 pairs, 94% positive - a ~15-20%
family-dependent bias on identical ice, consistent with residual
pedestal-fabric coupling growing with the axis-to-antenna angle.

This script makes that test durable and re-runnable: it is the
acceptance criterion for the chan_equal raw-channel calibration (the
paired same-ice difference should go to ~0 once the cross-pol pedestal
is removed at source), and the arbiter of how much of the systematic was
instead the 0/180 heading-wrap interpolation bug fixed at the 10 Aug
gate, which corrupted exactly the N-S family's sub-block rotations.

Per block: deep dlam = median LS sec_dlam_ls over 1150-1500 m (>= 5
finite cells); position = the block centre; family from the block's OWN
sec_az (curved connectors change family mid-frame), by nearest family
axis: N-S ~ 0/180, rows ~ 79, NW-SE ~ 148. A pair = the closest
cross-family block pair within max_sep, deduped to one per frame pair
(one grid crossing produces one pair, not a cluster). z_max special runs
(_z<N> tags) are excluded: they duplicate batch frames.

Usage: python crossing_pairs.py <out_dir> [max_sep_m] [tag_prefix]
"""
import glob
import os
import re
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
from pyproj import Transformer           # noqa: E402
from scipy.stats import binomtest        # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import DATA, INK, MUTED  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
MAX_SEP = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0
PREFIX = sys.argv[3] if len(sys.argv) > 3 else '2025'
Z_DEEP = (1150.0, 1500.0)     # same contrast band the map figure encodes
MIN_CELLS = 5                 # finite section cells required in the band
BLUE = '#2a78d6'

# Family axes [deg, mod 180] measured from the survey itself; a block
# joins the nearest axis within its tolerance, else stays unclassified.
FAMS = {'ns': (0.0, 20.0), 'row': (79.0, 15.0), 'nwse': (148.0, 15.0)}

T = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)


def axdist(a, b):
    """Distance between axis-valued angles [deg], mod 180."""
    return np.abs((np.asarray(a) - b + 90.0) % 180.0 - 90.0)


def classify(az):
    fam = np.full(az.shape, '', dtype=object)
    best = np.full(az.shape, np.inf)
    for name, (ctr, tol) in FAMS.items():
        d = axdist(az, ctr)
        hit = (d < tol) & (d < best)
        fam[hit] = name
        best[hit] = d[hit]
    return fam


def load_blocks():
    """One row per along-track block with a usable deep contrast."""
    rows = []
    for fn in sorted(glob.glob(os.path.join(
            DATA, 'quadpol_section_%s*.mat' % PREFIX))):
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        if re.search(r'_z\d+$', tag):
            continue
        with h5py.File(fn) as f:
            r = f['res']
            z = np.array(r['z']).ravel()
            secl = np.array(r['sec_dlam_ls']).T
            blat = np.array(r['sec_lat']).ravel()
            blon = np.array(r['sec_lon']).ravel()
            baz = np.array(r['sec_az']).ravel()
        mz = (z > Z_DEEP[0]) & (z < Z_DEEP[1])
        for b in range(secl.shape[1]):
            col = secl[mz, b]
            if np.isfinite(col).sum() < MIN_CELLS:
                continue
            rows.append((tag, blat[b], blon[b], baz[b],
                         float(np.nanmedian(col))))
    tags = np.array([r[0] for r in rows])
    lat = np.array([r[1] for r in rows])
    lon = np.array([r[2] for r in rows])
    az = np.array([r[3] for r in rows])
    dl = np.array([r[4] for r in rows])
    x, y = T.transform(lon, lat)
    return tags, x, y, classify(az), dl


def pair_families(tags, x, y, fam, dl, fa, fb):
    """Closest block pair per frame crossing, families fa vs fb."""
    ia = np.flatnonzero(fam == fa)
    ib = np.flatnonzero(fam == fb)
    best = {}
    for i in ia:
        d = np.hypot(x[ib] - x[i], y[ib] - y[i])
        for k in np.flatnonzero(d < MAX_SEP):
            j = ib[k]
            key = (tags[i], tags[j])
            if key not in best or d[k] < best[key][0]:
                best[key] = (d[k], i, j)
    pairs = [(i, j, s) for s, i, j in best.values()]
    diffs = np.array([dl[i] - dl[j] for i, j, _ in pairs])
    return pairs, diffs


def report(name, pairs, diffs):
    if len(diffs) == 0:
        print('%-12s no pairs' % name)
        return
    npos = int((diffs > 0).sum())
    p = binomtest(npos, len(diffs)).pvalue
    print('%-12s n %3d  median %+.4f  mean %+.4f  positive %3.0f%%  '
          'sign-test p %.2g' % (name, len(diffs), np.median(diffs),
                                diffs.mean(), 100.0 * npos / len(diffs), p))


def main():
    tags, x, y, fam, dl = load_blocks()
    for name, (ctr, _) in FAMS.items():
        print('%-5s %4d blocks' % (name, (fam == name).sum()))

    combos = [('ns', 'row'), ('nwse', 'row'), ('ns', 'nwse')]
    results = {}
    for fa, fb in combos:
        pairs, diffs = pair_families(tags, x, y, fam, dl, fa, fb)
        results[(fa, fb)] = (pairs, diffs)
        report('%s - %s' % (fa, fb), pairs, diffs)

    pairs, diffs = results[('ns', 'row')]
    a = np.array([dl[i] for i, _, _ in pairs])
    b = np.array([dl[j] for _, j, _ in pairs])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.4),
                                   layout='constrained')
    lim = (0.0, max(0.1, 1.05 * max(a.max(), b.max())))
    ax1.plot(lim, lim, '-', color=MUTED, lw=1.0, zorder=1)
    ax1.scatter(b, a, s=26, color=BLUE, edgecolors='white',
                linewidths=0.6, zorder=2)
    ax1.set_xlim(lim)
    ax1.set_ylim(lim)
    ax1.set_aspect('equal')
    ax1.set_xlabel('row block deep $\\Delta\\lambda$', color=INK)
    ax1.set_ylabel('N-S block deep $\\Delta\\lambda$', color=INK)

    ax2.hist(diffs, bins=21, color=BLUE, edgecolor='white', lw=0.6)
    ax2.axvline(0.0, color=MUTED, lw=1.0)
    med = np.median(diffs)
    ax2.axvline(med, color=INK, lw=1.4)
    ax2.annotate('median %+.4f' % med, xy=(med, 1.0),
                 xycoords=('data', 'axes fraction'), xytext=(4, -12),
                 textcoords='offset points', color=INK, fontsize=10)
    ax2.set_xlabel('same-ice N-S minus row $\\Delta\\lambda$', color=INK)
    ax2.set_ylabel('crossing pairs', color=INK)
    for ax in (ax1, ax2):
        ax.tick_params(colors=INK)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)

    os.makedirs(OUT, exist_ok=True)
    fig_fn = os.path.join(OUT, 'crossing_pairs.png')
    fig.savefig(fig_fn, dpi=200)
    npz_fn = os.path.join(OUT, 'crossing_pairs.npz')
    np.savez(npz_fn, **{
        '%s_%s_diffs' % (fa, fb): results[(fa, fb)][1]
        for fa, fb in combos})
    print('wrote %s and %s' % (fig_fn, npz_fn))


if __name__ == '__main__':
    main()

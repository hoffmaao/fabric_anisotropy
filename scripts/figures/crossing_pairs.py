"""Same-ice crossing-pair test of the heading-family dlam systematic.

The decisive test for a family-dependent bias in the deep principal
contrast: global comparisons across heading families are confounded with
the real NW-SE spatial gradient at Ridge A, but where two survey lines
from different families CROSS, their nearest along-track blocks see the
same ice, so the paired difference isolates whatever the acquisition
geometry adds. On the Ridge A sections as of 11 Aug 2026: N-S minus row
= +0.0068 median over 45 crossings, 91% positive - a ~15-20%
family-dependent bias on identical ice, consistent with residual
pedestal-fabric coupling growing with the axis-to-antenna angle. The
NW-SE minus row combination (-0.0271, n = 16) comes from one cluster of
end-of-row stubs and is not quotable in either direction.

This script makes that test durable and re-runnable: it is the
acceptance criterion for the chan_equal raw-channel calibration - the
paired same-ice difference should go to ~0 once the cross-pol pedestal
is removed at source. The 0/180 heading-wrap interpolation bug fixed at
the 10 Aug gate corrupted exactly the N-S family's sub-block rotations,
so it was the leading suspect for this systematic; the resweep settled
it separately (median deep-dlam change 0.00000 over 660 blocks), which
leaves instrument physics.

Per block: deep dlam = median LS sec_dlam_ls over 1150-1500 m (>= 5
finite cells); position = the block centre; family from the block's OWN
sec_az (curved connectors change family mid-frame), by nearest family
axis: N-S ~ 0/180, rows ~ 79, NW-SE ~ 148. A pair = the closest
cross-frame, cross-family block pair within max_sep, deduped to one per
CROSSING by clustering a frame pair's candidates on midpoint proximity,
at a radius derived per family combo from that combo's crossing corridor
so a shallow-angle combo groups as reliably as a near-perpendicular one.
Two frames that cross twice therefore count twice; none does in this
survey, so each of the 45 contributes a single pair.
The sign test treats pairs as independent when they are not
- one block can be the closest match at several crossings - so its p is
optimistic. z_max special runs (_z<N> tags) are excluded: they duplicate
batch frames.

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
# What sizes the cluster radius is the crossing corridor, not the block
# spacing: two lines whose axes differ by A stay within max_sep over a
# corridor of half-length max_sep / sin(A), and since candidates are taken
# closest-first the representative sits at the corridor centre, so the
# radius must cover that half-length. Derived per family combo rather than
# fixed, because max_sep is a CLI argument and the FAMS axes are
# survey-specific: a constant silently stops covering the corridor as soon
# as either moves. At the defaults this is 458 m for ns-row (corridor 306),
# 482 m for nwse-row (321) and 849 m for ns-nwse (566).
CROSS_MARGIN = 1.5
MIN_AXIS_SEP = 1.0            # deg; near-parallel families have no crossing
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
    if not rows:
        raise SystemExit('no quadpol_section_%s*.mat block with a usable '
                         'deep contrast under %s' % (PREFIX, DATA))
    tags = np.array([r[0] for r in rows])
    lat = np.array([r[1] for r in rows])
    lon = np.array([r[2] for r in rows])
    az = np.array([r[3] for r in rows])
    dl = np.array([r[4] for r in rows])
    x, y = T.transform(lon, lat)
    return tags, x, y, classify(az), dl


def pair_families(tags, x, y, fam, dl, fa, fb):
    """Closest block pair per frame crossing, families fa vs fb.

    The crossing, not the frame pair, is the unit: two frames that cross
    twice are two independent samples of the ice, while the cluster of
    blocks meeting at one crossing is a single sample. Candidates of one
    frame pair are therefore CLUSTERED by midpoint proximity rather than
    snapped to a fixed grid, which is not translation invariant: a cluster
    straddling a cell boundary would split and vote twice, the exact
    double-count this exists to prevent. Closest candidate first, so each
    cluster is represented by its own tightest pair, and the radius covers
    this combo's crossing corridor so a shallow-angle combo groups as
    reliably as a near-perpendicular one.
    """
    sep = max(float(axdist(FAMS[fa][0], FAMS[fb][0])), MIN_AXIS_SEP)
    radius = CROSS_MARGIN * MAX_SEP / np.sin(np.radians(sep))
    ia = np.flatnonzero(fam == fa)
    ib = np.flatnonzero(fam == fb)
    cand = {}
    for i in ia:
        d = np.hypot(x[ib] - x[i], y[ib] - y[i])
        for k in np.flatnonzero(d < MAX_SEP):
            j = ib[k]
            # One continuous track at its turn is not a crossing. On the
            # current Ridge A sections this removes nothing: a connector's
            # turn blocks fall in the gaps between the family tolerances
            # and go unclassified, so it is an invariant, not a filter.
            if tags[i] == tags[j]:
                continue
            cand.setdefault((str(tags[i]), str(tags[j])), []).append(
                (float(d[k]), int(i), int(j)))
    pairs = []
    for key in sorted(cand):
        seen = []
        for s, i, j in sorted(cand[key]):
            mx, my = 0.5 * (x[i] + x[j]), 0.5 * (y[i] + y[j])
            if any(np.hypot(mx - px, my - py) <= radius
                   for px, py in seen):
                continue
            seen.append((mx, my))
            pairs.append((i, j, s))
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
    for name in FAMS:
        print('%-5s %4d blocks' % (name, (fam == name).sum()))

    combos = [('ns', 'row'), ('nwse', 'row'), ('ns', 'nwse')]
    results = {}
    for fa, fb in combos:
        pairs, diffs = pair_families(tags, x, y, fam, dl, fa, fb)
        results[(fa, fb)] = (pairs, diffs)
        report('%s - %s' % (fa, fb), pairs, diffs)
    print('  sign-test p is optimistic: one block can be the closest match '
          'at several\n  crossings, so the pairs are not independent.')

    pairs, diffs = results[('ns', 'row')]
    if not pairs:
        raise SystemExit('no ns-row crossing within %.0f m in '
                         'quadpol_section_%s*.mat under %s'
                         % (MAX_SEP, PREFIX, DATA))
    a = np.array([dl[i] for i, _, _ in pairs])
    b = np.array([dl[j] for _, j, _ in pairs])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.4),
                                   layout='constrained')
    # dlam is signed, so the lower bound follows the data: clipping it at 0
    # would drop negative blocks from the scatter that the histogram, the
    # median and the reported n all still count.
    lim = (min(a.min(), b.min(), 0.0),
           max(0.1, 1.05 * max(a.max(), b.max())))
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

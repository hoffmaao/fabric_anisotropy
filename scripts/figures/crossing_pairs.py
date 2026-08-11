"""Same-ice crossing-pair test of the heading-family dlam systematic.

The decisive test for a family-dependent bias in the deep principal
contrast: global comparisons across heading families are confounded with
the real NW-SE spatial gradient at Ridge A, but where two survey lines
from different families CROSS, their nearest along-track blocks see the
same ice, so the paired difference isolates whatever the acquisition
geometry adds. On the Ridge A sections as of 11 Aug 2026: N-S minus row
= +0.0088 median over 28 crossings, 93% positive (26 of 28) - a ~17%
family-dependent bias on identical ice, consistent with residual
pedestal-fabric coupling growing with the axis-to-antenna angle. The
NW-SE minus row combination (+0.0098, n = 2) is one cluster of
end-of-row stubs and is not quotable in either direction, nor is N-S
minus NW-SE (+0.0110, n = 3).

This script makes that test durable and re-runnable: it is the
acceptance criterion for the chan_equal raw-channel calibration - the
paired same-ice difference should go to ~0 once the cross-pol pedestal
is removed at source. The 0/180 heading-wrap interpolation bug fixed at
the 10 Aug gate corrupted exactly the N-S family's sub-block rotations,
so it was the leading suspect for this systematic; the resweep settled
it separately (median deep-dlam change +0.00000 across 1574 blocks over
all 17 reswept frames), which leaves instrument physics.

Per block: deep dlam = median LS sec_dlam_ls over 1150-1500 m (>= 5
finite cells); position = the block centre; family from the block's OWN
sec_az (curved connectors change family mid-frame), by nearest family
axis: N-S ~ 0/180, rows ~ 79, NW-SE ~ 148. A pair = the closest
cross-frame, cross-family block pair within max_sep, deduped to one per
CROSSING by clustering EVERY candidate of the family combo on midpoint
proximity - globally, because a crossing is a place and a place does not
belong to a frame pair - at a radius derived per family combo from that
combo's crossing corridor, so a shallow-angle combo groups as reliably
as a near-perpendicular one. Two lines that cross twice therefore count
twice; the same place covered by several frames counts once, which at
Ridge A is common and is why n is 28 rather than the 264 raw candidates.
That the radius groups repeats without also fusing two DISTINCT
crossings was checked on the candidate set, before clustering, not on
the survivors (whose separations exceed the radius by construction and
so could not have shown the failure): no cluster is more than 101 m wide
across track, where two distinct parallel lines would be a ~1.5 km line
spacing apart, all 17 of the 45 -> 28 removals are one physical line
re-flown under a second frame tag, and the closest candidate midpoint
landing in a different cluster is 1487 m, 3.2x the 458 m ns-row radius.
merge_audit() re-derives both numbers on every run and names any cluster
that fails.
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

DEFAULT_MAX_SEP = 300.0
DEFAULT_PREFIX = '2025'
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
MIN_FAM_BLOCKS = 20           # below this main() warns the axes do not fit
# How wide across track a cluster may be and still be taken as one
# physical line covered more than once. Repeat passes over a line land
# within navigation scatter of each other (worst 101 m over the current
# sections), while two DISTINCT parallel lines are a line spacing apart
# (~1.5 km here), so anything between the two settles the question. Set
# nearer the low end: a spurious warning costs a look, a missed fusion
# costs a crossing that no count can show is gone.
LINE_TOL = 150.0
BLUE = '#2a78d6'

# Family axes [deg, mod 180] measured from the Ridge A survey itself; a
# block joins the nearest axis within its tolerance, else stays
# unclassified. These are survey-specific while the tag prefix is a CLI
# argument, so main() warns loudly when a family a requested combo needs
# comes back near-empty rather than reporting confident nonsense for
# another grid.
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


def load_blocks(prefix=DEFAULT_PREFIX):
    """One row per along-track block with a usable deep contrast."""
    rows = []
    for fn in sorted(glob.glob(os.path.join(
            DATA, 'quadpol_section_%s*.mat' % prefix))):
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
            rows.append((tag, b, blat[b], blon[b], baz[b],
                         float(np.nanmedian(col))))
    if not rows:
        raise SystemExit('no quadpol_section_%s*.mat block with a usable '
                         'deep contrast under %s' % (prefix, DATA))
    tags = np.array([r[0] for r in rows])
    bidx = np.array([r[1] for r in rows])
    lat = np.array([r[2] for r in rows])
    lon = np.array([r[3] for r in rows])
    az = np.array([r[4] for r in rows])
    dl = np.array([r[5] for r in rows])
    x, y = T.transform(lon, lat)
    return tags, bidx, x, y, classify(az), dl


def cluster_radius(fa, fb, max_sep=DEFAULT_MAX_SEP):
    """Radius [m] covering this family combo's crossing corridor."""
    sep = max(float(axdist(FAMS[fa][0], FAMS[fb][0])), MIN_AXIS_SEP)
    return CROSS_MARGIN * max_sep / np.sin(np.radians(sep))


def candidates(tags, x, y, fam, fa, fb, max_sep=DEFAULT_MAX_SEP):
    """Cross-frame, cross-family block pairs within max_sep, closest first."""
    ia = np.flatnonzero(fam == fa)
    ib = np.flatnonzero(fam == fb)
    out = []
    for i in ia:
        d = np.hypot(x[ib] - x[i], y[ib] - y[i])
        for k in np.flatnonzero(d < max_sep):
            j = ib[k]
            # One continuous track at its turn is not a crossing. On the
            # current Ridge A sections this removes nothing: a connector's
            # turn blocks fall in the gaps between the family tolerances
            # and go unclassified, so it is an invariant, not a filter.
            if tags[i] == tags[j]:
                continue
            out.append((float(d[k]), int(i), int(j)))
    return sorted(out)


def cluster(cand, x, y, radius):
    """Group candidates by midpoint proximity, closest candidate first.

    The crossing is the unit, and a crossing is a PLACE, not a pair of
    frames: the blocks meeting at one crossing are a single sample of the
    ice, wherever the frame boundaries happen to fall and however many
    times the line was flown. Candidates are therefore clustered across
    the WHOLE family combo, not within a frame-pair key. Scoping the
    dedup by key cannot express one-per-crossing and lets three cases
    vote more than once for one place: a line split across consecutive
    frames crossing another line lands under keys (A1, B) and (A2, B); a
    curved connector contributing blocks to both families near one
    crossing lands under both (C, D) and (D, C); and a line re-flown on
    another day lands under a second tag entirely. Normalising the key to
    an unordered pair would fix only the second, and nothing keyed on
    frames can fix the third.

    Clustering rather than snapping to a fixed grid, because a grid is
    not translation invariant: a cluster straddling a cell boundary would
    split and vote twice, the exact double-count this exists to prevent.
    Closest candidate first, so each cluster is headed by its own
    tightest pair, and the radius covers this combo's crossing corridor
    so a shallow-angle combo groups as reliably as a near-perpendicular
    one. Whether that radius is small enough to keep two DISTINCT
    crossings apart is not assumed - merge_audit() measures it.
    """
    reps = []
    groups = []
    for s, i, j in cand:
        mx, my = 0.5 * (x[i] + x[j]), 0.5 * (y[i] + y[j])
        hit = None
        for gi, (px, py) in enumerate(reps):
            if np.hypot(mx - px, my - py) <= radius:
                hit = gi
                break
        if hit is None:
            reps.append((mx, my))
            groups.append([])
            hit = len(groups) - 1
        groups[hit].append((s, i, j))
    return groups, reps


def track_dirs(tags, bidx, x, y, idx):
    """Unit along-track direction at each block, from its own frame.

    Taken from the block's index-adjacent siblings rather than by fitting
    a line to the set, because the set is exactly what is in question: if
    it holds two parallel lines, a fit orients ACROSS them whenever their
    separation beats the along-track extent, and then reports the fused
    pair as tight. Sibling displacements cannot do that. They are also in
    the projected grid, so they need none of the ~69 deg convergence that
    separates the geographic sec_az from an EPSG:3031 bearing here.
    """
    out = []
    for i in idx:
        m = np.flatnonzero((tags == tags[i]) & (np.abs(bidx - bidx[i]) <= 1))
        if m.size < 2:
            continue
        m = m[np.argsort(bidx[m])]
        dx, dy = x[m[-1]] - x[m[0]], y[m[-1]] - y[m[0]]
        n = np.hypot(dx, dy)
        if n > 0:
            out.append(np.arctan2(dx, dy))
    return out


def line_spread(tags, bidx, x, y, idx):
    """How many metres wide, across track, these blocks are [m].

    ~0 for one line however often it was re-flown, and the full line
    separation for two distinct parallel lines.
    """
    idx = list(idx)
    a = track_dirs(tags, bidx, x, y, idx)
    if len(idx) < 2 or not a:
        return 0.0
    th = 0.5 * np.angle(np.mean(np.exp(2j * np.array(a))))
    p = np.column_stack([x[idx], y[idx]])
    return float(np.ptp(p @ np.array([np.cos(th), -np.sin(th)])))


def merge_audit(tags, bidx, x, y, fam, fa, fb, max_sep=DEFAULT_MAX_SEP):
    """Did clustering fuse crossings of two DIFFERENT line pairs?

    The failure mode global clustering could have is over-merging: two
    genuinely distinct crossings closer together than the radius would be
    collapsed into one and one of them silently lost. Separations
    measured among the SURVIVORS cannot detect that, because greedy
    clustering leaves every survivor more than a radius from every other
    by construction, so a fused pair leaves no trace in that statistic.
    Both numbers here are therefore taken from the CANDIDATE set:

    on_line   the widest any cluster is ACROSS track, per side. Repeat
              passes over one line coincide to navigation scatter; two
              distinct parallel lines are a line spacing apart. A cluster
              over LINE_TOL is fusing line pairs, which is the failure,
              and is named.
    gap       the closest two candidate midpoints that landed in
              DIFFERENT clusters. Unlike a survivor separation this has
              no lower bound built into it, so it is the honest margin
              the radius has to sit under.
    """
    radius = cluster_radius(fa, fb, max_sep)
    cand = candidates(tags, x, y, fam, fa, fb, max_sep)
    groups, _ = cluster(cand, x, y, radius)
    on_line = 0.0
    bad = []
    for g in groups:
        if len(g) < 2:
            continue
        sa = line_spread(tags, bidx, x, y, [e[1] for e in g])
        sb = line_spread(tags, bidx, x, y, [e[2] for e in g])
        on_line = max(on_line, sa, sb)
        if max(sa, sb) > LINE_TOL:
            bad.append((sorted({str(tags[e[1]]) for e in g}),
                        sorted({str(tags[e[2]]) for e in g}), sa, sb))
    gap = np.inf
    if len(groups) > 1:
        m = np.array([[0.5 * (x[i] + x[j]), 0.5 * (y[i] + y[j])]
                      for g in groups for _, i, j in g])
        lab = np.array([gi for gi, g in enumerate(groups) for _ in g])
        d = np.hypot(m[:, None, 0] - m[None, :, 0],
                     m[:, None, 1] - m[None, :, 1])
        d[lab[:, None] == lab[None, :]] = np.inf
        gap = float(d.min())
    return radius, on_line, gap, bad


def pair_families(tags, x, y, fam, dl, fa, fb, max_sep=DEFAULT_MAX_SEP):
    """Closest block pair per crossing, families fa vs fb."""
    cand = candidates(tags, x, y, fam, fa, fb, max_sep)
    groups, _ = cluster(cand, x, y, cluster_radius(fa, fb, max_sep))
    pairs = [(i, j, s) for g in groups for s, i, j in g[:1]]
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


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    out_dir = argv[0] if len(argv) > 0 else '.'
    max_sep = float(argv[1]) if len(argv) > 1 else DEFAULT_MAX_SEP
    prefix = argv[2] if len(argv) > 2 else DEFAULT_PREFIX

    tags, bidx, x, y, fam, dl = load_blocks(prefix)
    counts = {name: int((fam == name).sum()) for name in FAMS}
    for name in FAMS:
        print('%-5s %4d blocks' % (name, counts[name]))

    combos = [('ns', 'row'), ('nwse', 'row'), ('ns', 'nwse')]
    for name in sorted({f for combo in combos for f in combo
                        if counts[f] < MIN_FAM_BLOCKS}):
        print('WARNING: family %r holds only %d blocks (< %d). The FAMS axes'
              '\n  (ns %.0f, row %.0f, nwse %.0f deg) are Ridge A survey '
              'geometry; a survey\n  flown on other headings needs its own '
              'axes, and until it has them every\n  combination using %r '
              'below is drawn from a near-empty family.'
              % (name, counts[name], MIN_FAM_BLOCKS, FAMS['ns'][0],
                 FAMS['row'][0], FAMS['nwse'][0], name))

    results = {}
    for fa, fb in combos:
        pairs, diffs = pair_families(tags, x, y, fam, dl, fa, fb, max_sep)
        results[(fa, fb)] = (pairs, diffs)
        report('%s - %s' % (fa, fb), pairs, diffs)
    print('  sign-test p is optimistic: one block can be the closest match '
          'at several\n  crossings, so the pairs are not independent.')

    print('merge audit (candidate set, before clustering):')
    for fa, fb in combos:
        radius, on_line, gap, bad = merge_audit(tags, bidx, x, y, fam,
                                                fa, fb, max_sep)
        print('  %-12s radius %4.0f m  widest cluster across track %3.0f m  '
              'nearest candidate\n               in another cluster %s'
              % ('%s - %s' % (fa, fb), radius, on_line,
                 'n/a' if not np.isfinite(gap)
                 else '%.0f m (%.1fx radius)' % (gap, gap / radius)))
        for ta, tb, sa, sb in bad:
            print('  WARNING: one cluster spans %s x %s and sits %.0f/%.0f m '
                  'off a single\n  line pair (> %.0f m): the radius is fusing '
                  'DISTINCT crossings, so n is too\n  low and this combo is '
                  'not usable until max_sep is reduced.'
                  % ('+'.join(ta), '+'.join(tb), sa, sb, LINE_TOL))

    pairs, diffs = results[('ns', 'row')]
    if not pairs:
        raise SystemExit('no ns-row crossing within %.0f m in '
                         'quadpol_section_%s*.mat under %s'
                         % (max_sep, prefix, DATA))
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

    os.makedirs(out_dir, exist_ok=True)
    fig_fn = os.path.join(out_dir, 'crossing_pairs.png')
    fig.savefig(fig_fn, dpi=200)
    npz_fn = os.path.join(out_dir, 'crossing_pairs.npz')
    np.savez(npz_fn, **{
        '%s_%s_diffs' % (fa, fb): results[(fa, fb)][1]
        for fa, fb in combos})
    print('wrote %s and %s' % (fig_fn, npz_fn))


if __name__ == '__main__':
    main()

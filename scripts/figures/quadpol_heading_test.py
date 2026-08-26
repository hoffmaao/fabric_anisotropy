"""Does quad-pol orientation follow the ice, or the antennas?

The claim made for a full scattering matrix is that fabric orientation is
recoverable from a SINGLE line, with no azimuth diversity. There is one
decisive test of that on an existing raster survey, and it needs no new
acquisition: drive the legs at different headings and ask whether the
recovered orientation, expressed GEOGRAPHICALLY, is the same on all of
them.

  theta_geo constant, independent of heading -> the method sees the ice
  theta_geo tracking the heading             -> the method sees the
                                                antennas, and has measured
                                                nothing about the fabric

Panel (a) is that test. The failure mode is the diagonal OFFSET by the
measured theta_ant, because run_quadpol_survey.m forms theta_geo as
theta_ant + track azimuth (mod 180): any point on that line has theta
fixed in the ANTENNA frame. The horizontal is the success mode.

The test only means something where the cross-polarized channels are
trustworthy in the first place, so panel (c) carries the reciprocity
screen that decides which frames are allowed into (a): a reciprocal medium
has S_hv == S_vh, and frames where the two channels are uncorrelated are
measuring system, not ice.

Usage: python quadpol_heading_test.py <out_dir> [site]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import h5py                              # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import DATA, INK, MUTED  # noqa: E402
from scar_style import FIGS  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else FIGS
SITE = sys.argv[2] if len(sys.argv) > 2 else 'ridge_a'
FN = os.path.join(DATA, 'quadpol_%s.mat' % SITE)

RECIP_MIN = 0.5       # HV/VH coherence a frame must clear to be believed
Z_BAND = (200.0, 1200.0)
# Scale bar for panel (d): the quad-pol LS contrast over the SAME 200-1200 m
# Z_BAND each frame here is summarised over, so the two numbers are directly
# comparable. Per site, because drawing one site's contrast across another's
# frames would be a comparison to nothing. Ridge A: median of sec_dlam_ls
# over the 36-frame survey (1.9e6 cells, p25/p75 0.029/0.067). This replaces
# the earlier two-azimuth solve reference, which the LS estimator superseded
# - its correlated (theta, P) errors biased it high - and it is the same bar
# ershadi_heading_test.py draws, since the two figures are shown as a pair
# and must not disagree about what the contrast is compared against.
REF_DLAM = {'ridge_a': 0.052}
C_OK, C_BAD = '#2a78d6', '#eb6834'


def circ_mean(a_deg):
    """Mean of an AXIS-valued angle: doubled, so 179 and 1 average to 0."""
    a = np.asarray(a_deg, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    m = np.mean(np.exp(2j * np.radians(a)))
    return np.degrees(0.5 * np.angle(m)) % 180


def circ_sd(a_deg):
    a = np.asarray(a_deg, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    R = np.abs(np.mean(np.exp(2j * np.radians(a))))
    return np.degrees(np.sqrt(-2 * np.log(max(R, 1e-12))) / 2)


def load():
    if not os.path.exists(FN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_quadpol_survey.m for this '
            'site and mirror the result there' % FN)
    rows = []
    with h5py.File(FN) as f:
        g = f['Q']
        for i in range(g['tag'].shape[0]):
            def a(name):
                return np.array(f[g[name][i, 0]]).ravel()
            tag = ''.join(chr(c) for c in a('tag'))
            z = a('z')
            m = (z > Z_BAND[0]) & (z < Z_BAND[1])
            if m.sum() < 10:
                continue
            with np.errstate(invalid='ignore'):
                rows.append(dict(
                    tag=tag, seg=tag[:11],
                    track=float(a('track_az')[0]),
                    theta=circ_mean(a('theta_geo')[m]),
                    dlam=np.nanmedian(a('dlam')[m]),
                    aniso=np.nanmedian(a('aniso')[m]),
                    xr=np.nanmedian(a('xr')[m]),
                    recip=np.nanmedian(a('recip')[m])))
    return rows


def main():
    rows = load()
    ok = [r for r in rows if r['recip'] > RECIP_MIN]
    bad = [r for r in rows if not r['recip'] > RECIP_MIN]
    ta = np.array([r['track'] for r in ok])
    th = np.array([r['theta'] for r in ok])
    dl = np.array([r['dlam'] for r in ok])
    ant = (th - ta) % 180
    print('%s: %d frames, %d reciprocal, %d not'
          % (SITE, len(rows), len(ok), len(bad)))
    print('  theta_geo   circ-sd %.1f deg' % circ_sd(th))
    print('  theta_ant   circ-sd %.1f deg' % circ_sd(ant))
    # nan-aware, as ershadi_heading_test.py is: a frame whose coherence
    # gate leaves no finite sample contributes a NaN, and plain
    # median/percentile would print "nan" for the whole survey over it.
    print('  dlam        median %.3f  IQR %.3f-%.3f  (%d of %d finite)'
          % (np.nanmedian(dl), *np.nanpercentile(dl, [25, 75]),
             int(np.isfinite(dl).sum()), dl.size))

    fig = plt.figure(figsize=(14.0, 4.9), layout='constrained')
    gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 1.0, 1.0, 1.0])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[0, 2])
    axd = fig.add_subplot(gs[0, 3])

    # (a) the test. run_quadpol_survey.m forms theta_geo = theta_ant +
    # track_az (mod 180), so the antenna-frame hypothesis is the diagonal
    # offset by the measured theta_ant, not the 1:1 line: drawn at 1:1 the
    # frames sit a whole theta_ant above their own reference.
    th_ant0 = circ_mean(ant)
    gr = np.linspace(0, 180, 361)
    pred = (gr + th_ant0) % 180
    # split at the mod-180 wrap so no vertical stroke crosses the panel
    cut = np.flatnonzero(np.abs(np.diff(pred)) > 90) + 1
    axa.plot(np.insert(gr, cut, np.nan), np.insert(pred, cut, np.nan),
             '-', color='#c0392b', lw=1.6, zorder=2,
             label=(r'fixed in ANTENNA frame ($\theta_{\rm ant}$ = '
                    r'%.0f$^\circ$)' % th_ant0) + '\n(measures nothing)')
    axa.axhline(circ_mean(th), color='#1baf7a', lw=1.6, ls='--', zorder=2,
                label='fixed in GEOGRAPHIC frame\n(measures the ice)')
    axa.plot(ta, th, 'o', color=C_OK, ms=8, markeredgecolor='white',
             markeredgewidth=0.8, zorder=4)
    axa.set_xlim(0, 180)
    axa.set_ylim(0, 180)
    axa.set_xticks([0, 45, 90, 135, 180])
    axa.set_yticks([0, 45, 90, 135, 180])
    axa.set_aspect('equal')
    axa.set_xlabel('track azimuth (deg E of N)', color=INK)
    axa.set_ylabel(r'recovered $\theta$ (deg E of N)', color=INK)
    axa.set_title('(a) heading test, %d reciprocal frames' % len(ok),
                  fontsize=10, color=INK)
    axa.legend(loc='best', fontsize=7, frameon=False)
    axa.grid(alpha=0.25, lw=0.6)

    # (b) the same angle in the antenna frame, where the scatter collapses
    axb.plot(ta, ant, 'o', color=C_OK, ms=7, markeredgecolor='white',
             markeredgewidth=0.8)
    axb.axhline(circ_mean(ant), color='#c0392b', lw=1.4, ls='--')
    axb.set_xlim(0, 180)
    axb.set_ylim(0, 180)
    axb.set_xticks([0, 45, 90, 135, 180])
    axb.set_xlabel('track azimuth (deg E of N)', color=INK)
    axb.set_ylabel(r'$\theta$ in the antenna frame (deg)', color=INK)
    axb.set_title('(b) antenna frame: sd %.1f$^\\circ$ vs %.1f$^\\circ$'
                  % (circ_sd(ant), circ_sd(th)), fontsize=10, color=INK)
    axb.grid(alpha=0.25, lw=0.6)

    # (c) the screen that decides which frames (a) is allowed to use
    segs = sorted({r['seg'] for r in rows})
    for r in rows:
        c = C_OK if r['recip'] > RECIP_MIN else C_BAD
        axc.plot(segs.index(r['seg']), r['recip'], 'o', color=c, ms=6,
                 markeredgecolor='white', markeredgewidth=0.6)
    axc.axhline(RECIP_MIN, color=MUTED, lw=1.0, ls=':')
    axc.set_xticks(range(len(segs)))
    axc.set_xticklabels([s[4:] for s in segs], rotation=90, fontsize=6)
    axc.set_ylim(-0.05, 1.05)
    axc.set_ylabel('HV/VH coherence', color=INK)
    axc.set_title('(c) reciprocity, by segment', fontsize=10, color=INK)
    axc.grid(alpha=0.25, lw=0.6)

    # (d) the contrast, which does not depend on the cross-pol channels
    axd.plot(ta, dl, 'o', color=C_OK, ms=7, markeredgecolor='white',
             markeredgewidth=0.8)
    if SITE in REF_DLAM:
        axd.axhline(REF_DLAM[SITE], color='black', lw=1.4, ls='--',
                    label='quad-pol LS, same band')
    axd.set_xlim(0, 180)
    axd.set_xticks([0, 45, 90, 135, 180])
    axd.set_xlabel('track azimuth (deg E of N)', color=INK)
    axd.set_ylabel(r'$\Delta\lambda$', color=INK)
    axd.set_title('(d) contrast vs heading', fontsize=10, color=INK)
    if SITE in REF_DLAM:
        axd.legend(loc='upper right', fontsize=7.5, frameon=False)
    axd.grid(alpha=0.25, lw=0.6)

    for ax in (axa, axb, axc, axd):
        ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)

    fig.suptitle('Quad-pol orientation on %s: does it follow the ice or the '
                 'antennas?' % SITE.replace('_', ' ').title(),
                 fontsize=12, color=INK)
    out = os.path.join(OUT, 'scar_quadpol_heading_%s.png' % SITE)
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', out)


if __name__ == '__main__':
    main()

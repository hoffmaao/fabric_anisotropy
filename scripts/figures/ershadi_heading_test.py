"""Ershadi et al. (2022) on a towed survey: ice-fixed or antenna-fixed?

The published method assumes a stationary sounding at one antenna
orientation, so run_ershadi_survey.m applies it in 200-trace blocks of
near-constant heading rather than a frame at a time. That turns the
vehicle's heading variation from a violated assumption into the
experiment: every block yields a (heading, theta) pair, and

  theta_geo independent of heading -> the method is measuring the ice
  theta_geo tracking the heading    -> it is measuring the antennas

Blocks make this far sharper than the frame-level version in
quadpol_heading_test.py. There are hundreds rather than tens, and because
a single frame contributes many blocks at slightly different headings the
test runs WITHIN a frame as well as between frames - same ice, same
calibration, same processing, only the heading differing. A frame-level
test cannot separate "heading changed" from "we drove somewhere else".

Panel (c) is the discriminator stated as one number: the circular spread
of theta in each frame. Whichever frame it is smaller in is the frame the
signal is fixed in.

Usage: python ershadi_heading_test.py <out_dir> [site]
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

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
SITE = sys.argv[2] if len(sys.argv) > 2 else 'ridge_a'
FN = os.path.join(DATA, 'ershadi_%s.mat' % SITE)

RECIP_MIN = 0.5     # HV/VH coherence a block must clear
SPREAD_MAX = 3.0    # deg of heading wander allowed inside one block
QUAL_MIN = 1.0      # dB of cross-pol azimuthal modulation
REF_P = 0.05        # the independent two-azimuth solve, for scale only
# Ershadi et al. gate on |C_HHVV| > 0.4. This survey never reaches it - the
# maximum over 1790 blocks is 0.359 and the median 0.221 - so the number is
# carried here to be drawn as the bar the data does not clear, rather than
# applied as a filter that would reject everything.
ERSHADI_COH_GATE = 0.4
C_OK, C_BAD = '#2a78d6', '#eb6834'


def circ_mean(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    m = np.mean(np.exp(2j * np.radians(a)))
    return np.degrees(0.5 * np.angle(m)) % 180


def circ_sd(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size < 2:
        return np.nan
    R = np.abs(np.mean(np.exp(2j * np.radians(a))))
    return np.degrees(np.sqrt(-2 * np.log(max(R, 1e-12))) / 2)


def load():
    if not os.path.exists(FN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_ershadi_survey.m for this '
            'site and mirror the result there' % FN)
    out = {}
    with h5py.File(FN) as f:
        g = f['B']
        n = g['az'].shape[0]
        for k in ('az', 'az_spread', 'theta_ant', 'theta_geo', 'dlam',
                  'cmag', 'qual', 'xr', 'recip', 'lat', 'lon'):
            out[k] = np.array([np.array(f[g[k][i, 0]]).ravel()[0]
                               for i in range(n)], float)
        out['tag'] = [''.join(chr(c) for c in
                              np.array(f[g['tag'][i, 0]]).ravel())
                      for i in range(n)]
    return out


def main():
    B = load()
    n = B['az'].size
    keep = ((B['recip'] > RECIP_MIN) & (B['az_spread'] < SPREAD_MAX)
            & (B['qual'] > QUAL_MIN) & np.isfinite(B['theta_geo']))
    print('%s: %d blocks, %d pass (recip>%.1f, wander<%.0f deg, qual>%.0f dB)'
          % (SITE, n, keep.sum(), RECIP_MIN, SPREAD_MAX, QUAL_MIN))
    if keep.sum() < 10:
        raise SystemExit('too few blocks pass the screens to test anything')

    az = B['az'][keep]
    tg = B['theta_geo'][keep]
    ta = B['theta_ant'][keep]
    dl = B['dlam'][keep]
    tags = [t for t, k in zip(B['tag'], keep) if k]

    sd_geo, sd_ant = circ_sd(tg), circ_sd(ta)
    print('  theta circular sd:  geographic %.1f deg   antenna %.1f deg'
          % (sd_geo, sd_ant))
    print('  verdict: fixed in the %s frame'
          % ('GEOGRAPHIC (ice)' if sd_geo < sd_ant
             else 'ANTENNA (instrument)'))
    # nan-aware: the coherence gate leaves some blocks without a dlam, and
    # plain median/percentile return NaN if even one is present, which
    # printed "dlam median nan" for a set that was 92% finite.
    print('  dlam median %.3f  IQR %.3f-%.3f  (%d of %d finite)'
          % (np.nanmedian(dl), *np.nanpercentile(dl, [25, 75]),
             np.isfinite(dl).sum(), dl.size))

    # Within-frame test: frames that themselves span a range of headings
    per = {}
    for t, a, g_, an in zip(tags, az, tg, ta):
        per.setdefault(t, []).append((a, g_, an))
    wf = [(t, v) for t, v in per.items()
          if len(v) >= 6 and (max(x[0] for x in v) - min(x[0] for x in v)) > 8]
    print('  %d frames span >8 deg of heading internally' % len(wf))

    fig = plt.figure(figsize=(14.0, 4.9), layout='constrained')
    gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 1.0, 1.0, 1.0])
    axa, axb, axc, axd = [fig.add_subplot(gs[0, i]) for i in range(4)]

    gr = np.linspace(0, 180, 2)
    axa.plot(gr, gr, '-', color='#c0392b', lw=1.6, zorder=2,
             label='fixed in ANTENNA frame')
    axa.axhline(circ_mean(tg), color='#1baf7a', lw=1.6, ls='--', zorder=2,
                label='fixed in GEOGRAPHIC frame')
    axa.plot(az, tg, 'o', color=C_OK, ms=4, alpha=0.55,
             markeredgecolor='none', zorder=4)
    axa.set_xlim(0, 180)
    axa.set_ylim(0, 180)
    axa.set_aspect('equal')
    axa.set_xticks([0, 45, 90, 135, 180])
    axa.set_yticks([0, 45, 90, 135, 180])
    axa.set_xlabel('block heading (deg E of N)', color=INK)
    axa.set_ylabel(r'recovered $\theta$ (deg E of N)', color=INK)
    axa.set_title('(a) heading test, %d blocks' % keep.sum(), fontsize=10,
                  color=INK)
    axa.legend(loc='upper left', fontsize=7, frameon=False)

    axb.plot(az, ta, 'o', color=C_OK, ms=4, alpha=0.55,
             markeredgecolor='none')
    axb.axhline(circ_mean(ta), color='#c0392b', lw=1.4, ls='--')
    axb.set_xlim(0, 180)
    axb.set_ylim(0, 180)
    axb.set_xticks([0, 45, 90, 135, 180])
    axb.set_xlabel('block heading (deg E of N)', color=INK)
    axb.set_ylabel(r'$\theta$ in the antenna frame (deg)', color=INK)
    axb.set_title('(b) antenna frame', fontsize=10, color=INK)

    # (c) the verdict as one number, per frame and overall
    labs = ['geographic\n(ice)', 'antenna\n(instrument)']
    vals = [sd_geo, sd_ant]
    cols = ['#1baf7a', '#c0392b']
    axc.bar(labs, vals, color=cols, width=0.6)
    for i, v in enumerate(vals):
        axc.annotate('%.1f$^\\circ$' % v, xy=(i, v), xytext=(0, 3),
                     textcoords='offset points', ha='center', fontsize=9,
                     color=INK)
    axc.set_ylabel(r'circular sd of $\theta$ (deg)', color=INK)
    axc.set_title('(c) which frame is it fixed in?', fontsize=10, color=INK)

    axd.plot(az, dl, 'o', color=C_OK, ms=4, alpha=0.55,
             markeredgecolor='none')
    axd.axhline(REF_P, color='black', lw=1.4, ls='--',
                label='two-azimuth solve')
    axd.set_xlim(0, 180)
    axd.set_xticks([0, 45, 90, 135, 180])
    axd.set_ylim(0, max(0.12, float(np.nanpercentile(dl, 97))))
    axd.set_xlabel('block heading (deg E of N)', color=INK)
    axd.set_ylabel(r'$\Delta\lambda$', color=INK)
    axd.set_title('(d) contrast vs heading', fontsize=10, color=INK)
    axd.legend(loc='upper right', fontsize=7.5, frameon=False)

    for ax in (axa, axb, axc, axd):
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)

    fig.suptitle('Ershadi et al. (2022) on %s, per heading block'
                 % SITE.replace('_', ' ').title(), fontsize=12, color=INK)
    out = os.path.join(OUT, 'scar_ershadi_heading_%s.png' % SITE)
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', out)


if __name__ == '__main__':
    main()

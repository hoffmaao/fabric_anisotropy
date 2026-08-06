"""EastGRIP core eigenvalues against the radar's horizontal contrast.

Left: the core's three orientation-tensor eigenvalues with depth
(Weikusat et al. 2022, PANGAEA.949248, CC-BY-4.0, area-weighted). The
smallest collapses to ~0.01 below 450 m - one direction with essentially
no c-axes, which under along-flow extension is the flow direction - while
the other two separate to ~0.41 and ~0.58. That is a strong vertical
girdle.

Right: what the radar should see, against what it does. At nadir the
radar is sensitive ONLY to the difference of the two HORIZONTAL diagonal
components of the orientation tensor; the vertical eigenvalue does not
enter (eps_x - eps_y = deps*(lam_x - lam_y)). The core gives eigenvalue
MAGNITUDES only - no eigenvectors, no azimuth - so which pair is
horizontal is not determined by this dataset, and all three candidate
differences are drawn:

    e2-e1   if the LARGEST eigenvalue is vertical
    e3-e1   if the MIDDLE eigenvalue is vertical
    e3-e2   if the SMALLEST is vertical (unlikely here: it would mean
            almost no c-axes near vertical)

The radar points are P = lam_max - lam_min from the azimuthal solve over
the lines within RADIUS_KM of the borehole (egrip_azimuthal.py), which
needs no core reorientation - the SAME solve panel (c) of that figure
draws, called through the shared solve_bands, so the two EastGRIP slides
cannot quote different P for one band. Vertical bars are the depth band;
horizontal bars span that solve and the two run on the individual
fringe-rate estimators, which fail in opposite directions and so bracket
the answer.

Usage: python egrip_core_compare.py <out_dir>
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import INK, MUTED         # noqa: E402
from egrip_azimuthal import (BANDS, RADIUS_KM, core_table,  # noqa: E402
                             solve_bands, stage_lines)

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
# validated categorical slots 1-3 (dataviz reference palette, light mode)
C_E = ['#2a78d6', '#eb6834', '#1baf7a']
C_RAD = '#4a3aa7'


def runmed(z, v, win=40.0):
    """Running median in depth, so the scatter does not hide the trend."""
    out = np.full(v.shape, np.nan)
    for i, zi in enumerate(z):
        m = np.abs(z - zi) <= win
        if m.sum() >= 3:
            out[i] = np.median(v[m])
    return out


def bracket(P, P_lag, P_uw):
    """Lowest and highest solved P per band, over the three solves.

    The point drawn is the headline solve, so the bracket has to CONTAIN
    it: the headline P is the norm of the mean of the two estimator
    solutions, which is not bounded below by the smaller of their norms
    when the two solutions are not colinear, and an errorbar computed as
    P - min(P_lag, P_uw) would then go negative.

    Reduced band by band over the finite entries rather than with
    np.nanmin/np.nanmax over the stack, because the estimator solves are
    rejected in most bands - the agreement gate fails below ~650 m - and
    an all-NaN column would make those emit an empty-slice RuntimeWarning
    straight into the printed table.
    """
    lo = np.full(len(BANDS), np.nan)
    hi = np.full(len(BANDS), np.nan)
    stack = np.vstack([P, P_lag, P_uw])
    for bi in np.flatnonzero(np.isfinite(P)):
        v = stack[np.isfinite(stack[:, bi]), bi]
        lo[bi], hi[bi] = v.min(), v.max()
    return lo, hi


def main():
    os.makedirs(OUT, exist_ok=True)
    Z, E = core_table()
    print('core: %d sections, %.0f-%.0f m' % (len(Z), Z.min(), Z.max()))

    lines = stage_lines()
    P, _, _, _ = solve_bands(lines)
    P_lag, _, _, _ = solve_bands(lines, 0)
    P_uw, _, _, _ = solve_bands(lines, 1)
    zc = np.array([0.5 * (a + b) for a, b in BANDS])
    Plo, Phi = bracket(P, P_lag, P_uw)

    print('\n%10s %8s %8s %8s | %8s %8s %8s'
          % ('depth_m', 'P', 'P_lag', 'P_unwrap', 'e2-e1', 'e3-e1', 'e3-e2'))
    for bi, (z0, z1) in enumerate(BANDS):
        m = (Z >= z0) & (Z < z1)
        if not np.isfinite(P[bi]) or m.sum() < 3:
            continue
        d21 = np.median(E[m, 1] - E[m, 0])
        d31 = np.median(E[m, 2] - E[m, 0])
        d32 = np.median(E[m, 2] - E[m, 1])
        print('%10s %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f'
              % ('%d-%d' % (z0, z1), P[bi], P_lag[bi], P_uw[bi],
                 d21, d31, d32))

    fig, (axe, axd) = plt.subplots(1, 2, figsize=(11.6, 6.4), sharey=True,
                                   layout='constrained')

    for k, lbl in enumerate([r'$e_1$ (smallest)', r'$e_2$',
                             r'$e_3$ (largest)']):
        axe.plot(E[:, k], Z, '.', color=C_E[k], ms=2.0, alpha=0.22)
        axe.plot(runmed(Z, E[:, k]), Z, '-', color=C_E[k], lw=2.0, label=lbl)
    axe.axvline(1 / 3, color=MUTED, lw=0.8, ls=':')
    axe.annotate('isotropic 1/3', xy=(1 / 3, 1690), xytext=(4, 0),
                 textcoords='offset points', fontsize=7.5, color=MUTED,
                 va='bottom')
    axe.set_xlabel('orientation-tensor eigenvalue', color=INK)
    axe.set_ylabel('depth (m)', color=INK)
    axe.set_title('EastGRIP core (Weikusat et al. 2022, area-weighted)',
                  fontsize=11, color=INK)
    axe.legend(loc='upper right', fontsize=8.5, frameon=False)

    d21 = runmed(Z, E[:, 1] - E[:, 0])
    d31 = runmed(Z, E[:, 2] - E[:, 0])
    d32 = runmed(Z, E[:, 2] - E[:, 1])
    axd.plot(d31, Z, '-', color='0.25', lw=1.8,
             label=r'$e_3-e_1$  (middle eigenvalue vertical)')
    axd.plot(d21, Z, color='0.45', lw=1.8, ls='--',
             label=r'$e_2-e_1$  (largest vertical)')
    axd.plot(d32, Z, color='0.65', lw=1.6, ls=':',
             label=r'$e_3-e_2$  (smallest vertical)')
    ok = np.isfinite(P)
    axd.errorbar(P[ok], zc[ok],
                 xerr=[P[ok] - Plo[ok], Phi[ok] - P[ok]],
                 yerr=[zc[ok] - np.array([b[0] for b in BANDS])[ok],
                       np.array([b[1] for b in BANDS])[ok] - zc[ok]],
                 fmt='o', color=C_RAD, ms=7, lw=1.6, capsize=3, zorder=6,
                 label='radar, azimuthal solve\n(%d lines within %.1f km)'
                       % (len(lines), RADIUS_KM))
    axd.set_xlabel('horizontal eigenvalue difference', color=INK)
    axd.set_title('what the radar is sensitive to', fontsize=11, color=INK)
    axd.set_xlim(0, 0.72)
    axd.legend(loc='lower right', fontsize=8, frameon=False)

    for ax in (axe, axd):
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    axe.set_ylim(1715, 0)

    # Upper left: the lower right is the legend, and the three candidate
    # curves run through the middle of the panel
    axd.annotate('radar solves only where its two independent\n'
                 'fringe-rate estimators agree, which fails\n'
                 'below ~650 m',
                 xy=(0.03, 0.97), xycoords='axes fraction', fontsize=7.5,
                 color=MUTED, va='top')
    fig.suptitle('EastGRIP: core eigenvalues and the radar horizontal '
                 'contrast at the same site', fontsize=12.5, color=INK)
    out = os.path.join(OUT, 'scar_egrip_core_compare.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('\nwrote', out)


if __name__ == '__main__':
    main()

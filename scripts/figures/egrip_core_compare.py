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
nine lines within 2 km of the borehole (egrip_azimuthal.py), which needs
no core reorientation. Vertical bars are the depth band; horizontal bars
are the spread between the two independent fringe-rate estimators, which
fail in opposite directions and so bracket the answer.

Usage: python egrip_core_compare.py <out_dir>
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import h5py                              # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from egrip_azimuthal import (BANDS, C_ICE, EG_LAT, EG_LON,  # noqa: E402
                             FRAMES, K_CONV, LAG, AGREE_TOL, P_BOUND,
                             DATA, FLOW_AZ, CORE, haversine_km,
                             load_line, near_mask, track_azimuth)

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
INK, MUTED = '#0b0b0b', '#52514e'
# validated categorical slots 1-3 (dataviz reference palette, light mode)
C_E = ['#2a78d6', '#eb6834', '#1baf7a']
C_RAD = '#4a3aa7'


def core_table():
    with open(CORE) as f:
        lines = f.read().split('\n')
    h = next(i for i, l in enumerate(lines) if l.startswith('Bag\t'))
    cols = lines[h].split('\t')
    iz = cols.index('Depth ice/snow [m]')
    ie = [cols.index('EVA%d (Weighted (statistic))' % k) for k in (1, 2, 3)]
    Z, E = [], []
    for l in lines[h + 1:]:
        p = l.split('\t')
        if len(p) <= max(ie):
            continue
        try:
            Z.append(float(p[iz]))
            E.append([float(p[i]) for i in ie])
        except ValueError:
            pass
    o = np.argsort(Z)
    return np.array(Z)[o], np.array(E)[o]


def runmed(z, v, win=40.0):
    """Running median in depth, so the scatter does not hide the trend."""
    out = np.full(v.shape, np.nan)
    for i, zi in enumerate(z):
        m = np.abs(z - zi) <= win
        if m.sum() >= 3:
            out[i] = np.median(v[m])
    return out


def band_rates_both(d, near):
    """(lag, unwrap) dlam and coherence per band - as in egrip_azimuthal."""
    dt = float(np.median(np.diff(d['t'])))
    tb = (d['t'] - np.nanmedian(d['surf'])) * 1e6
    w = np.nan_to_num(d['coh'][:, near])
    uw = np.unwrap(np.angle(np.nansum(d['ifg'][:, near] * w, axis=1)))
    X = d['ifg'][LAG:, near] * np.conj(d['ifg'][:-LAG, near])
    out = []
    for z0, z1 in BANDS:
        t0, t1 = 2 * z0 / C_ICE * 1e6, 2 * z1 / C_ICE * 1e6
        m = (tb >= t0) & (tb < t1)
        ml = (tb[:-LAG] >= t0) & (tb[:-LAG] < t1)
        if m.sum() < 8 or ml.sum() < 5:
            out.append((np.nan, np.nan, np.nan))
            continue
        r_uw = np.polyfit(tb[m], uw[m], 1)[0] / (2 * np.pi)
        r_lag = np.angle(np.nansum(X[ml, :])) / (2 * np.pi * LAG * dt) * 1e-6
        out.append((r_lag / K_CONV, r_uw / K_CONV,
                    float(np.nanmedian(d['coh'][m][:, near]))))
    return np.array(out)


def solve(lines, col):
    """P and theta per band using one estimator column (0 lag, 1 unwrap)."""
    az = np.array([d['az'] for d in lines])
    P = np.full(len(BANDS), np.nan)
    TH = np.full(len(BANDS), np.nan)
    for bi in range(len(BANDS)):
        y = np.array([d['rows'][bi, col] for d in lines])
        y2 = np.array([d['rows'][bi, 1 - col] for d in lines])
        w = np.array([d['rows'][bi, 2] for d in lines])
        den = np.maximum(np.abs(y), np.abs(y2))
        with np.errstate(invalid='ignore', divide='ignore'):
            agree = np.abs(y - y2) / np.where(den > 0, den, np.nan)
        ok = (np.isfinite(y) & np.isfinite(w) & (w > 0.35)
              & np.isfinite(agree) & (agree < AGREE_TOL))
        if ok.sum() < 3:
            continue
        a = np.radians(az[ok])
        A = np.c_[-np.cos(2 * a), -np.sin(2 * a)]
        W = np.sqrt(w[ok])[:, None]
        sol, *_ = np.linalg.lstsq(A * W, y[ok] * W[:, 0], rcond=None)
        Pb = np.hypot(*sol)
        if Pb > P_BOUND:
            continue
        P[bi] = Pb
        TH[bi] = np.degrees(0.5 * np.arctan2(sol[1], sol[0])) % 180
    return P, TH


def main():
    os.makedirs(OUT, exist_ok=True)
    Z, E = core_table()
    print('core: %d sections, %.0f-%.0f m' % (len(Z), Z.min(), Z.max()))

    lines = []
    for tag in FRAMES:
        d = load_line(tag)
        if d is None:
            continue
        nm = near_mask(d)
        if nm.sum() < 8:
            continue
        d['az'] = track_azimuth(d['lat'][nm], d['lon'][nm])
        d['rows'] = band_rates_both(d, nm)
        lines.append(d)
    P_lag, TH_lag = solve(lines, 0)
    P_uw, TH_uw = solve(lines, 1)
    zc = np.array([0.5 * (a + b) for a, b in BANDS])
    P = np.nanmean(np.vstack([P_lag, P_uw]), axis=0)
    Plo = np.nanmin(np.vstack([P_lag, P_uw]), axis=0)
    Phi = np.nanmax(np.vstack([P_lag, P_uw]), axis=0)

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

    for k, lbl in enumerate([r'$e_1$ (smallest)', r'$e_2$', r'$e_3$ (largest)']):
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
    axd.plot(d21, Z, '-', color='0.45', lw=1.8, ls='--',
             label=r'$e_2-e_1$  (largest vertical)')
    axd.plot(d32, Z, '-', color='0.65', lw=1.6, ls=':',
             label=r'$e_3-e_2$  (smallest vertical)')
    ok = np.isfinite(P)
    axd.errorbar(P[ok], zc[ok],
                 xerr=[P[ok] - Plo[ok], Phi[ok] - P[ok]],
                 yerr=[zc[ok] - np.array([b[0] for b in BANDS])[ok],
                       np.array([b[1] for b in BANDS])[ok] - zc[ok]],
                 fmt='o', color=C_RAD, ms=7, lw=1.6, capsize=3, zorder=6,
                 label='radar, azimuthal solve\n(9 lines within 2 km)')
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

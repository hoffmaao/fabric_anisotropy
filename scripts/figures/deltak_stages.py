"""Per-stage delta-k diagnostic figure for one SLC pair.

Plots the traveltime-difference profile at each rung of the delta-k
ladder against the two phase-based estimators, to locate where the Ridge
A amplitude is lost. deltak_vs_joint.py shows the delta-k dlam standard
deviation falling from 1.3x the joint chain near the surface to 18x
below 1500 m, but cannot say which step does the suppressing:

  stage A   adjacent narrow sub-band pairs, smoothed 5x5 on the cell grid
  stage Q   quarter-band pairs, integers resolved by A
  stage B   half-band pair, integers resolved by Q (what the chain ships)
  snaphu    ptt.blendTraveltime on the unwrapped SNAPHU phase
  coreg     ptt.blendTraveltime on the coregistration offsets alone

Read it this way: amplitude already missing at stage A is the analysis
cell (dk.cell_twtt x dk.cell_ntr) low-passing dtau in depth before the
ladder runs, so no later rung can recover it. Amplitude lost between A
and B is the ladder itself, where each rung rounds to the integer the
previous one implies.

The right panel is the RMS depth increment of each profile in sliding
traveltime bands, which is what the inversion actually reads: dlam
responds to the gradient of dtau, so a profile that tracks SNAPHU in
absolute terms can still invert to a suppressed fabric contrast if its
increments are smoothed away. Raw increments are too noisy to compare by
eye, so they are reduced to a band RMS - the amplitude retained.

Input: the deltak_stages_<day_seg>_<frm>.mat written on mem1 by
opr_fabric/server/run_deltak_stages.m.

Usage: python deltak_stages.py <deltak_stages_*.mat> [out_dir]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figs')

C_ICE = 1.68e8  # m/s, two-way traveltime to depth below the firn column

# Drawn in this order, so the reference estimators sit under the ladder
# stages the diagnostic is about rather than covering them.
STAGES = [
    ('prof_coreg', 'coregistration offsets', '#8172b3', 1.2, ':', 0.55),
    ('prof_snaphu', 'SNAPHU unwrapped phase', '#55a868', 1.8, '--', 0.9),
    ('prof_A', 'stage A (narrow sub-bands)', '#4c72b0', 1.6, '-', 0.9),
    ('prof_Q', 'stage Q (quarter-band)', '#dd8452', 1.6, '-', 0.9),
    ('prof_B', 'stage B (half-band, shipped)', '#c44e52', 2.2, '-', 1.0),
]


def col(d, name):
    """1-D float column, or None when the runner left the leg empty."""
    if name not in d:
        return None
    v = np.atleast_1d(np.asarray(d[name], float)).ravel()
    return None if v.size == 0 else v


def increment_rms(t, v, step, band):
    """RMS increment of v per `step` of traveltime, within `band` windows.

    Resampling every estimator onto the same `step` before differencing is
    what makes the comparison fair: the ladder stages live on the coarse
    analysis-cell grid and the phase estimators on the full fast-time
    grid, so raw sample-to-sample increments would differ by grid spacing
    alone rather than by retained signal.
    """
    ok = np.isfinite(v)
    if ok.sum() < 3:
        return np.array([]), np.array([])
    t, v = t[ok], v[ok]
    edges = np.arange(t[0], t[-1], step)
    if edges.size < 3:
        return np.array([]), np.array([])
    dv = np.diff(np.interp(edges, t, v))
    tc = 0.5 * (edges[:-1] + edges[1:])
    centres = np.arange(tc[0], tc[-1], band / 2.0)
    out_t, out_v = [], []
    for c in centres:
        sel = np.abs(tc - c) <= band / 2.0
        if sel.sum() >= 3:
            out_t.append(c)
            out_v.append(np.sqrt(np.mean(dv[sel] ** 2)))
    return np.asarray(out_t), np.asarray(out_v)


def main():
    fn = sys.argv[1]
    d = loadmat(fn, squeeze_me=True)
    surf = float(d['surf_mean'])
    t_full = np.asarray(d['Time'], float).ravel() - surf
    t_cell = np.asarray(d['t_cell'], float).ravel() - surf
    tag = '%s_%03d' % (str(d['day_seg']), int(d['frm']))

    # The runner puts every rung on the full fast-time grid and reduces it
    # with the same ptt.blockAverage call, so all five profiles share an
    # axis. t_cell only sets the increment step below: differencing at the
    # raw 3.3 ns sampling would measure interpolation noise, not signal,
    # since the ladder cannot resolve structure finer than its cell.
    step = float(np.median(np.diff(t_cell))) if t_cell.size > 2 else 100e-9
    axes_for = {k: t_full for k, _, _, _, _, _ in STAGES}

    band = 20 * step  # RMS window: wide enough to be readable, narrow
                      # enough to keep the near-surface/deep contrast

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 8), sharey=True)
    for key, label, c, lw, ls, alpha in STAGES:
        v = col(d, key)
        if v is None:
            print('%-14s absent from %s' % (key, os.path.basename(fn)))
            continue
        t = axes_for[key]
        n = min(t.size, v.size)
        ax0.plot(1e9 * v[:n], 1e6 * t[:n], ls, color=c, lw=lw, alpha=alpha,
                 label=label)
        tc, dv = increment_rms(t[:n], v[:n], step, band)
        if tc.size:
            ax1.plot(1e9 * dv, 1e6 * tc, ls, color=c, lw=lw, alpha=alpha,
                     label=label)

    # ptt.imgCombSeam places each combine at surface + ic(1) in absolute
    # TWTT and masks [t_c - win, t_c + 2*win], so on this below-surface
    # axis the band is [ic(1) - win, ic(1) + 2*win].
    ic = np.atleast_1d(np.asarray(d.get('img_comb', []), float)).ravel()
    for b in range(ic.size // 3):
        t_after, win = ic[3 * b], ic[3 * b + 2]
        if not np.isfinite(t_after):
            continue
        for ax in (ax0, ax1):
            ax.axhspan(1e6 * (t_after - win),
                       1e6 * (t_after + 2 * win),
                       color='0.85', zorder=0,
                       label='_' if b else 'waveform-combine seam (masked)')

    ax0.set_xlabel(r'$\Delta\tau$ (ns)')
    ax0.set_ylabel(r'TWTT below surface ($\mu$s)')
    ax0.set_title('Traveltime difference by estimator')
    ax1.set_xlabel(r'RMS $\Delta(\Delta\tau)$ per %.0f ns (ns)' % (step * 1e9))
    ax1.set_title('Retained increment (what the inversion reads)')
    ax1.set_xlim(left=0)
    ax0.invert_yaxis()  # once: the sharey axis is common to both panels
    for ax in (ax0, ax1):
        ax.grid(alpha=0.3)
    ax0.axvline(0, color='0.5', lw=0.8)
    ax0.legend(loc='lower right', fontsize=8)

    secax = ax0.secondary_yaxis(
        'right', functions=(lambda us: us * 1e-6 * C_ICE / 2,
                            lambda z: 2 * z / C_ICE * 1e6))
    secax.set_ylabel('approx. depth (m)')

    fig.suptitle('Delta-k ladder stages vs phase estimators, %s' % tag)
    fig.tight_layout()

    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'deltak_stages_%s.png' % tag)
    fig.savefig(out, dpi=160, bbox_inches='tight')
    print('wrote %s' % out)

    # The number the figure is really about: amplitude retained per stage.
    print('\nprofile std by TWTT band below surface (ns)')
    bands = [(0, 2), (2, 5), (5, 10), (10, 20)]
    hdr = '  %-10s' % 'band(us)' + ''.join(
        '%10s' % k.replace('prof_', '') for k, _, _, _, _, _ in STAGES)
    print(hdr)
    for lo, hi in bands:
        row = '  %-10s' % ('%g-%g' % (lo, hi))
        for key, _, _, _, _, _ in STAGES:
            v = col(d, key)
            if v is None:
                row += '%10s' % '-'
                continue
            t = axes_for[key]
            n = min(t.size, v.size)
            sel = (t[:n] >= lo * 1e-6) & (t[:n] < hi * 1e-6)
            vv = v[:n][sel]
            vv = vv[np.isfinite(vv)]
            row += '%10s' % ('%.3f' % (1e9 * vv.std()) if vv.size > 1 else '-')
        print(row)


if __name__ == '__main__':
    main()

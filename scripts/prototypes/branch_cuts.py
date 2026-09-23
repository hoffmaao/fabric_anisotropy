"""InSAR-style residue/branch-cut analysis of the margin interferogram.

The 'break points' in the wrapped polarimetric interferogram are phase
residues in the InSAR sense: 2x2 loops where the wrapped gradients do
not close. This prototype applies the standard InSAR chain to the
Thwaites eastern-margin extract:

  1. residue map (Goldstein & Werner 1988) + density vs coherence,
  2. Goldstein-Werner adaptive spectral filtering of the COMPLEX
     interferogram (the InSAR fix for broken fringes),
  3. residues recounted after filtering,
  4. Herraez/skimage reliability-sorted unwrapping of the filtered phase
     vs the production SNAPHU solution (disagreement in integer fringes),
  5. fast-time phase-gradient profiles (complex conjugate product) raw
     vs filtered, the quantity the banded inversion actually needs.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'figures'))
import scar_style                            # noqa: E402
from scipy.io import loadmat
from skimage.restoration import unwrap_phase

SEG = '20240108_01_001'
OUT = scar_style.out_dir(sys.argv, 1, os.path.expanduser('~/projects/fabric_anisotropy/figs'))


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def residues(ph):
    """Charge map (+1/-1) on the dual grid of a wrapped phase image."""
    d1 = wrap(ph[:-1, 1:] - ph[:-1, :-1])
    d2 = wrap(ph[1:, 1:] - ph[:-1, 1:])
    d3 = wrap(ph[1:, :-1] - ph[1:, 1:])
    d4 = wrap(ph[:-1, :-1] - ph[1:, :-1])
    s = d1 + d2 + d3 + d4
    return np.rint(s / (2 * np.pi)).astype(np.int8)


def goldstein_filter(C, alpha=0.8, win=64, step=16):
    """Goldstein-Werner adaptive filter of a complex interferogram.

    Overlapping windows; each window's spectrum S is reweighted by
    |smooth(S)|^alpha, then overlap-added with a raised-cosine taper.
    """
    from scipy.ndimage import uniform_filter
    H, W = C.shape
    out = np.zeros((H, W), np.complex64)
    wsum = np.zeros((H, W), np.float32)
    taper1 = np.hanning(win)
    taper = np.outer(taper1, taper1).astype(np.float32) + 1e-6
    for r0 in range(0, H - win + 1, step):
        for c0 in range(0, W - win + 1, step):
            tile = C[r0:r0 + win, c0:c0 + win]
            S = np.fft.fft2(tile)
            mag = uniform_filter(np.abs(S), 3)
            Sf = S * (mag / (mag.max() + 1e-30)) ** alpha
            f = np.fft.ifft2(Sf).astype(np.complex64)
            out[r0:r0 + win, c0:c0 + win] += f * taper
            wsum[r0:r0 + win, c0:c0 + win] += taper
    ok = wsum > 0
    out[ok] /= wsum[ok]
    return out


def main():
    d = loadmat(os.path.expanduser(f'~/data/opr/margin/margin_{SEG}.mat'))
    ph = d['phase_wrapped'].astype(np.float32)
    coh = np.nan_to_num(d['coherence'].astype(np.float32))
    snaphu = d['phase_unwrapped'].astype(np.float32)
    t_us = d['Time'].ravel() * 1e6
    lats = d['Latitude'].ravel()
    lons = d['Longitude'].ravel()
    R = 6371.0
    p1, p2 = np.radians(lats[:-1]), np.radians(lats[1:])
    dl = np.radians(lons[1:] - lons[:-1])
    a = np.sin((p2 - p1) / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    dist = np.concatenate([[0], np.cumsum(2 * R * np.arcsin(np.sqrt(a)))])

    C = np.exp(1j * ph).astype(np.complex64)

    # 1-2. residues before/after Goldstein filtering
    res0 = residues(ph)
    Cf = goldstein_filter(C)
    phf = np.angle(Cf)
    res1 = residues(phf)
    n0, n1 = np.abs(res0).sum(), np.abs(res1).sum()
    npix = res0.size
    print(f'residues: raw {n0} ({1e2*n0/npix:.2f}% of loops), '
          f'filtered {n1} ({1e2*n1/npix:.2f}%)  [x{n0/max(n1, 1):.1f} fewer]')

    # residue density vs coherence
    cmid = 0.25 * (coh[:-1, :-1] + coh[1:, :-1] + coh[:-1, 1:] + coh[1:, 1:])
    bins = np.linspace(0, 1, 21)
    dens0 = np.full(20, np.nan)
    dens1 = np.full(20, np.nan)
    for i in range(20):
        m = (cmid >= bins[i]) & (cmid < bins[i + 1])
        if m.sum() > 500:
            dens0[i] = np.abs(res0[m]).mean()
            dens1[i] = np.abs(res1[m]).mean()

    # 4. Herraez unwrap of the filtered phase vs production SNAPHU
    herraez = unwrap_phase(np.ma.array(phf, mask=coh < 0.1)).filled(np.nan)
    diff = (herraez - snaphu)
    # remove the global constant (both unwrappers float by a constant)
    diff -= np.nanmedian(diff)
    fringes = diff / (2 * np.pi)
    frac_agree = np.nanmean(np.abs(fringes) < 0.5)
    print(f'Herraez(filtered) vs SNAPHU: {1e2*frac_agree:.1f}% of pixels '
          f'within half a fringe (after constant removal)')

    # 5. fast-time gradient profiles, raw vs filtered (block-averaged)
    dt = np.median(np.diff(t_us)) * 1e-6
    fc = 750e6

    def grad_profile(Cin, cols, smooth=151):
        pp = Cin[1:, cols] * np.conj(Cin[:-1, cols])
        w = (coh[1:, cols] * coh[:-1, cols])**2
        acc = np.nansum(w * pp, axis=1)
        # running coherent average in fast time (still no unwrapping)
        k = np.ones(smooth) / smooth
        acc = np.convolve(acc, k, mode='same')
        g = np.angle(acc) / dt / (2 * np.pi * fc) * 1e6  # ps/us
        g[np.convolve(np.nanmean(coh[1:, cols], axis=1), k, 'same') < 0.15] = np.nan
        return g

    nb = 10
    edges = np.linspace(0, C.shape[1], nb + 1).astype(int)
    b = 4  # a mid-margin block
    cols = np.arange(edges[b], edges[b + 1])
    g_raw = grad_profile(C, cols)
    g_fil = grad_profile(Cf, cols)

    # ---- figure
    fig = plt.figure(figsize=(15, 10), layout='constrained')
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1])

    sl = slice(0, ph.shape[0], 2)
    ax = fig.add_subplot(gs[0, 0])
    ax.pcolormesh(dist, t_us[sl], ph[sl], cmap='twilight', vmin=-np.pi,
                  vmax=np.pi, shading='auto', rasterized=True)
    ry, rx = np.where(res0[sl] != 0)
    ax.plot(dist[rx], t_us[sl][ry], ',', color='lime', alpha=0.25)
    ax.invert_yaxis()
    ax.set_ylabel('TWTT (us)')
    ax.set_title(f'raw wrapped phase + residues (n={n0})', fontsize=11)
    ax.set_xlabel('distance (km)')

    ax = fig.add_subplot(gs[0, 1])
    ax.pcolormesh(dist, t_us[sl], phf[sl], cmap='twilight', vmin=-np.pi,
                  vmax=np.pi, shading='auto', rasterized=True)
    ry, rx = np.where(res1[sl] != 0)
    ax.plot(dist[rx], t_us[sl][ry], ',', color='lime', alpha=0.4)
    ax.invert_yaxis()
    ax.set_title(f'Goldstein-filtered + residues (n={n1})', fontsize=11)
    ax.set_xlabel('distance (km)')

    ax = fig.add_subplot(gs[0, 2])
    pc = ax.pcolormesh(dist, t_us[sl], fringes[sl], cmap='RdBu_r', vmin=-3,
                       vmax=3, shading='auto', rasterized=True)
    ax.invert_yaxis()
    ax.set_title('Herraez(filtered) - SNAPHU (fringes)', fontsize=11)
    ax.set_xlabel('distance (km)')
    fig.colorbar(pc, ax=ax, label='fringes', pad=0.02)

    ax = fig.add_subplot(gs[1, 0])
    mids = 0.5 * (bins[:-1] + bins[1:])
    ax.semilogy(mids, dens0, 'o-', label='raw')
    ax.semilogy(mids, dens1, 's-', label='Goldstein-filtered')
    ax.set_xlabel('coherence')
    ax.set_ylabel('residue density (per loop)')
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title('break points live where coherence dies', fontsize=11)

    ax = fig.add_subplot(gs[1, 1:])
    ax.plot(g_raw, t_us[:-1], color='0.5', lw=1.0, label='raw conj-product gradient')
    ax.plot(g_fil, t_us[:-1], color='crimson', lw=1.8,
            label='after Goldstein filter')
    ax.set_ylim(t_us[-1], 0)
    ax.set_xlim(-900, 900)
    ax.set_xlabel(r'd$(\Delta\tau)$/dt (ps/us)')
    ax.set_ylabel('TWTT (us)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    ax.set_title(f'fast-time traveltime-difference gradient, block {b} '
                 '(the observable the banded inversion needs)', fontsize=11)

    fig.suptitle(f'InSAR residue / branch-cut view of the margin '
                 f'interferogram, {SEG}', fontsize=14)
    out = os.path.join(OUT, f'branch_cut_analysis_{SEG}.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

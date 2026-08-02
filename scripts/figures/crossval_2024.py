"""Cross-validate the two independent 2024-25 polarimetric processings.

Matches co-located inversion blocks from the Lilien-product and
Paden-product CSARP_fabric_joint outputs and compares dlam sample by
sample (scatter + season-median profiles).

Usage: python crossval_2024.py [fabric_batch_root] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')

DEPTH = np.arange(0, 1900, 5.0)


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def load_blocks(pattern):
    blocks = []
    for fn in sorted(glob.glob(pattern)):
        d = loadmat(fn)
        dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
        qual = d['dlam_quality']
        clip = d.get('dlam_clipped')
        for b in range(dlam.shape[1]):
            if not np.any(np.isfinite(dlam[:, b])):
                continue
            col = np.full(DEPTH.size, np.nan)
            for k in range(dlam.shape[0]):
                pegged = (clip is not None and clip.size and clip[k, b] == 1) \
                    or abs(dlam[k, b]) > 0.6
                ok = np.isfinite(dlam[k, b]) and not pegged \
                    and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
                if ok:
                    m = (DEPTH >= top[k, b]) & (DEPTH < bot[k, b])
                    col[m] = dlam[k, b]
            blocks.append((d['Latitude'][0, b], d['Longitude'][0, b], col))
    return blocks


def main():
    L = load_blocks(f'{ROOT}/joint/2024_Antarctica_Ground2/*/Data_*.mat')
    J = load_blocks(f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat')
    print(f'blocks: Lilien-based {len(L)}, Paden-based {len(J)}')

    Jlat = np.array([b[0] for b in J])
    Jlon = np.array([b[1] for b in J])
    pairs = []
    for lat, lon, col in L:
        dists = haversine_km(lat, lon, Jlat, Jlon)
        i = int(np.argmin(dists))
        if dists[i] < 1.0:
            pairs.append((col, J[i][2]))
    print(f'matched pairs within 1 km: {len(pairs)}')

    a = np.array([p[0] for p in pairs])
    b = np.array([p[1] for p in pairs])
    both = np.isfinite(a) & np.isfinite(b)
    x, y = a[both], b[both]
    r = np.corrcoef(x, y)[0, 1]
    mad = np.median(np.abs(x - y))
    print(f'N: {x.size}, r = {r:.3f}, median |diff| = {mad:.4f}, '
          f'median bias = {np.median(x - y):+.4f}')

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    hb = ax.hexbin(x, y, gridsize=60, extent=(-0.12, 0.12, -0.12, 0.12),
                   cmap='magma', mincnt=1)
    ax.plot([-0.12, 0.12], [-0.12, 0.12], 'c--', lw=1)
    ax.set_xlabel(r'$\Delta\lambda$ (Lilien product)')
    ax.set_ylabel(r'$\Delta\lambda$ (Paden product)')
    ax.set_title(f'2024-25 season, {len(pairs)} co-located blocks\n'
                 f'r = {r:.2f}, median |diff| = {mad:.3f}')
    fig.colorbar(hb, ax=ax, label='samples')

    ax = axes[1]
    with np.errstate(invalid='ignore'):
        ax.plot(np.nanmedian(a, axis=0), DEPTH, 'b-', lw=2, label='Lilien product')
        ax.plot(np.nanmedian(b, axis=0), DEPTH, 'r--', lw=2, label='Paden product')
    ax.invert_yaxis()
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel(r'median $\Delta\lambda$')
    ax.set_ylabel('Depth (m)')
    ax.set_title('Season-median profiles')
    ax.legend()
    ax.grid(alpha=0.3)

    fig.suptitle('Cross-validation: independent polarimetric processings, same inversion',
                 fontsize=12)
    out = os.path.join(OUT, 'crossval_2024_lilien_vs_paden.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

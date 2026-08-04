"""Season fabric transects from CSARP_fabric outputs.

Usage: python transect_seasons.py [product_subdir] [fabric_batch_root] [out_dir]
where product_subdir is e.g. 'joint' (default) selecting
<fabric_batch_root>/<product_subdir>/<season>/... inputs.
"""
import glob, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

import os
ROOT = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def load_season(season, product):
    """Return per-block columns on a common depth grid + stats."""
    depth_grid = np.arange(0, 1900, 5.0)
    fns = sorted(glob.glob(f'{ROOT}/{product}/{season}/*/Data_*.mat'))
    cols, lats, lons = [], [], []
    n_peg = n_tot = 0
    for fn in fns:
        d = loadmat(fn)
        dlam, top, bot, qual = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth'], d['dlam_quality']
        clip = d.get('dlam_clipped')
        itp = d.get('dlam_interpolated')
        for b in range(dlam.shape[1]):
            if not np.any(np.isfinite(dlam[:, b])):
                continue
            n_tot += 1
            pegged = np.abs(dlam[:, b]) > 0.6
            if clip is not None and clip.size:
                pegged |= clip[:, b] == 1
            n_peg += int(np.any(pegged))
            # dtau interpolated across a masked gap (e.g. a waveform-combine
            # seam) rather than measured at these nodes
            if itp is not None and itp.size:
                pegged |= itp[:, b] == 1
            col = np.full(depth_grid.size, np.nan)
            for k in range(dlam.shape[0]):
                ok = np.isfinite(dlam[k, b]) and not pegged[k] \
                    and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
                if ok:
                    m = (depth_grid >= top[k, b]) & (depth_grid < bot[k, b])
                    col[m] = dlam[k, b]
            cols.append(col); lats.append(d['Latitude'][0, b]); lons.append(d['Longitude'][0, b])
    if not cols:
        return None
    cols = np.array(cols).T
    lats, lons = np.array(lats), np.array(lons)
    step = haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
    dist = np.concatenate([[0], np.cumsum(np.minimum(step, 2.0))])
    gap = np.concatenate([[False], step > 2.0])
    cols[:, gap] = np.nan
    return dict(cols=cols, dist=dist, depth=depth_grid, n_tot=n_tot,
                n_peg=n_peg, n_frames=len(fns))


seasons = ['2022_Antarctica_Ground', '2023_Antarctica_Ground',
           '2024_Antarctica_Ground2', '2025_Antarctica_Ground2']
product = sys.argv[1] if len(sys.argv) > 1 else 'joint'

fig, axes = plt.subplots(4, 1, figsize=(14, 13))
pc = None
for ax, season in zip(axes, seasons):
    s = load_season(season, product)
    if s is None:
        ax.set_title(f'{season}: no data'); continue
    edges = np.concatenate([[s['dist'][0] - 0.25],
                            (s['dist'][:-1] + s['dist'][1:]) / 2,
                            [s['dist'][-1] + 0.25]])
    depth_edges = np.concatenate([s['depth'], [s['depth'][-1] + 5]])
    pc = ax.pcolormesh(edges, depth_edges, s['cols'], cmap='RdBu_r',
                       vmin=-0.15, vmax=0.15)
    ax.invert_yaxis(); ax.set_ylabel('Depth (m)')
    ax.set_title(f"{season.replace('_', ' ')}: {s['n_frames']} frames, "
                 f"{s['n_tot']} blocks, {s['n_peg']} pegged "
                 f"({100 * s['n_peg'] / s['n_tot']:.0f}%)")
    print(f"{season}: frames={s['n_frames']} blocks={s['n_tot']} "
          f"pegged={s['n_peg']} ({100 * s['n_peg'] / s['n_tot']:.0f}%)")
axes[-1].set_xlabel('Along-traverse distance (km, segment gaps compressed)')
if pc is None:
    print(f'no data found under {ROOT}/{product}; nothing to plot')
    sys.exit(1)
fig.colorbar(pc, ax=axes, extend='both', label=r'$\Delta\lambda = \lambda_{cross} - \lambda_{along}$',
             shrink=0.6)
fig.suptitle('Horizontal fabric contrast, EAGER traverses 2022-2026 (joint inversion)',
             fontsize=13)
os.makedirs(OUT, exist_ok=True)
fig.savefig(f'{OUT}/transect_4season.png', dpi=140, bbox_inches='tight')
print(f'saved {OUT}/transect_4season.png')

"""Delta-k fabric sections: eigenvalue difference along each traverse.

Plots the delta-k batch's inverted fabric contrast dlam = lam_cross -
lam_along (the horizontal eigenvalue difference, survey-line frame) as
along-profile depth sections, one panel per product:
  - Thwaites 2023_Antarctica_Ground (CSARP_fabric_deltak): the chain is
    validated against SNAPHU there (fringe-fan magnitudes reproduced).
  - Ridge A 2024_Antarctica_Ground2 (CSARP_fabric_deltak_jp): shown for
    completeness but delta-k is ~5-6x amplitude-suppressed vs the joint
    chain at Ridge A (see deltak_vs_joint.py) - treat as unvalidated.

Blocks from all frames are ordered by GPS time and drawn on a common
depth grid; pegged intervals (|dlam| at the 2/3 bound) are hatched out
as NaN.

Usage: python deltak_sections.py [fabric_batch_dir] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

BATCH = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figs')

PRODUCTS = [
    # seg_prefix keeps the Thwaites panel to the January drive; the
    # 2023 season's February segments are a different site
    ('deltak', 'Thwaites 2023-24 (validated regime)', '202401'),
    ('deltak_jp', 'Ridge A 2024-25 (UNVALIDATED: delta-k amplitude '
     'suppressed vs joint chain)', ''),
]
BOUND = 2.0 / 3.0
DZ = 25.0  # depth grid (m)


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def load_blocks(product_dir, seg_prefix=''):
    """All blocks of all frames, ordered by GPS time."""
    blocks = []
    for fn in sorted(glob.glob(os.path.join(
            product_dir, seg_prefix + '*', 'Data_*.mat'))):
        d = loadmat(fn, squeeze_me=True)

        def col2d(a):
            a = np.asarray(a, float)
            return a.reshape(-1, 1) if a.ndim == 1 else a  # 1-block frame

        dlam = col2d(d['dlam'])
        top = col2d(d['dlam_top_depth'])
        bot = col2d(d['dlam_bot_depth'])
        gps = np.atleast_1d(d['GPS_time'])
        lat = np.atleast_1d(d['Latitude'])
        lon = np.atleast_1d(d['Longitude'])
        for b in range(dlam.shape[1]):
            blocks.append({'gps': gps[b], 'lat': lat[b], 'lon': lon[b],
                           'dlam': dlam[:, b], 'top': top[:, b],
                           'bot': bot[:, b]})
    blocks.sort(key=lambda k: k['gps'])
    return blocks


def section(blocks):
    """Common-depth-grid section (nz x nblocks) + distance axis (km)."""
    lat = np.array([b['lat'] for b in blocks])
    lon = np.array([b['lon'] for b in blocks])
    x = np.concatenate([[0], np.cumsum(
        haversine_km(lat[:-1], lon[:-1], lat[1:], lon[1:]))])
    zmax = np.nanmax(np.concatenate([b['bot'] for b in blocks]))
    z = np.arange(0, zmax + DZ, DZ)
    sec = np.full((len(z) - 1, len(blocks)), np.nan)
    zc = 0.5 * (z[:-1] + z[1:])
    for j, b in enumerate(blocks):
        v = b['dlam'].copy()
        v[np.abs(v) >= 0.98 * BOUND] = np.nan  # pegged at the bound
        for i in range(len(v)):
            if np.isfinite(b['top'][i]) and np.isfinite(b['bot'][i]):
                rows = (zc >= b['top'][i]) & (zc < b['bot'][i])
                sec[rows, j] = v[i]
    return x, z, sec


fig, axes = plt.subplots(len(PRODUCTS), 1, figsize=(11, 4 * len(PRODUCTS)))
axes = np.atleast_1d(axes)
for ax, (prod, title, seg_prefix) in zip(axes, PRODUCTS):
    blocks = load_blocks(os.path.join(BATCH, prod), seg_prefix)
    if not blocks:
        ax.set_visible(False)
        continue
    x, z, sec = section(blocks)
    # draw contiguous runs separately: a lone block after a traverse
    # gap must not be stretched across the whole gap
    gaps = np.where(np.diff(x) > 10.0)[0]
    starts = np.concatenate([[0], gaps + 1])
    ends = np.concatenate([gaps, [len(x) - 1]])
    pc = None
    for s, e in zip(starts, ends):
        xs = x[s:e + 1]
        if len(xs) == 1:
            xe = np.array([xs[0] - 1.0, xs[0] + 1.0])
        else:
            xe = np.concatenate([[xs[0]], 0.5 * (xs[:-1] + xs[1:]),
                                 [xs[-1]]])
        pc = ax.pcolormesh(xe, z, sec[:, s:e + 1], cmap='RdBu_r',
                           vmin=-0.15, vmax=0.15)
    ax.invert_yaxis()
    ax.set_xlabel('Distance along traverse (km)')
    ax.set_ylabel('Depth (m)')
    ax.set_title('%s - %d blocks' % (title, len(blocks)), fontsize=11)
    ax.set_facecolor('0.85')  # NaN / no-solution cells
    fig.colorbar(pc, ax=ax, pad=0.015,
                 label=r'$\Delta\lambda = \lambda_{cross} - \lambda_{along}$')

fig.suptitle('Delta-k fabric inversion: horizontal eigenvalue difference',
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(OUT, 'deltak_fabric_sections.png')
fig.savefig(out, dpi=200, bbox_inches='tight')
print('wrote', out)

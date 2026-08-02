"""Survey summary: Antarctic coverage map, fabric profiles, drive headings.

Each dataset is a survey spanning (potentially very different) flow
regimes, and even single surveys drive at varying orientations relative to
the principal fabric axes. This figure summarizes: (a) where each survey
went (South Polar Stereographic map with coastline), (b) the median
horizontal eigenvalue-contrast profile with IQR per survey, and (c) the
distribution of drive headings (mod 180, since dlam is symmetric under
line reversal) per survey.

Usage: python summary_map_profiles.py [fabric_batch_root] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')

DEPTH = np.arange(0, 1900, 5.0)
DATASETS = [
    ('2022-23 survey', f'{ROOT}/joint/2022_Antarctica_Ground/*/Data_*.mat', 'tab:gray'),
    ('2023-24 survey', f'{ROOT}/joint/2023_Antarctica_Ground/*/Data_*.mat', 'tab:blue'),
    ('2024-25 (Lilien proc.)', f'{ROOT}/joint/2024_Antarctica_Ground2/*/Data_*.mat', 'tab:orange'),
    ('2024-25 (Paden proc.)', f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat', 'tab:red'),
    ('2025-26 survey', f'{ROOT}/joint/2025_Antarctica_Ground2/*/Data_*.mat', 'tab:green'),
]


def bearing_deg(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing, degrees clockwise from north."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(x, y)) % 360.0


def load_dataset(pattern):
    """QC'd block profiles, positions, and per-frame drive headings."""
    cols, lats, lons, headings = [], [], [], []
    for fn in sorted(glob.glob(pattern)):
        d = loadmat(fn)
        dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
        qual = d['dlam_quality']
        clip = d.get('dlam_clipped')
        blat = d['Latitude'][0, :]
        blon = d['Longitude'][0, :]
        # Headings between consecutive blocks of the frame, mod 180
        if blat.size > 1:
            h = bearing_deg(blat[:-1], blon[:-1], blat[1:], blon[1:]) % 180.0
            headings.extend(h[np.isfinite(h)])
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
            cols.append(col)
            lats.append(blat[b])
            lons.append(blon[b])
    return (np.array(cols), np.array(lats), np.array(lons), np.array(headings))


def main():
    data = [(name, color) + tuple(load_dataset(pattern)[i] for i in range(4))
            for name, pattern, color in DATASETS]
    data = [(n, c, cols, la, lo, hd) for (n, c, cols, la, lo, hd) in data if cols.size]
    if not data:
        print('no data found under', ROOT)
        return

    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1.0, 0.8], height_ratios=[1, 1])

    # (a) Coverage map
    all_lat = np.concatenate([d[3] for d in data])
    all_lon = np.concatenate([d[4] for d in data])
    if HAVE_CARTOPY:
        proj = ccrs.SouthPolarStereo()
        ax = fig.add_subplot(gs[:, 0], projection=proj)
        lon_span = np.ptp(all_lon)
        if lon_span > 180:
            ax.set_extent([-180, 180, -90, -63], ccrs.PlateCarree())
        else:
            pad_lat = max(1.0, 0.2 * np.ptp(all_lat))
            pad_lon = max(2.0, 0.2 * lon_span)
            ax.set_extent([all_lon.min() - pad_lon, all_lon.max() + pad_lon,
                           all_lat.min() - pad_lat, min(-60, all_lat.max() + pad_lat)],
                          ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.6)
        ax.gridlines(draw_labels=True, lw=0.3, color='gray', alpha=0.5,
                     x_inline=False, y_inline=True)
        # Later surveys revisit earlier sites, so draw in reverse order and
        # render the co-located Lilien processing as open circles
        for name, color, cols, la, lo, hd in reversed(data):
            if 'Lilien' in name:
                ax.scatter(lo, la, s=45, facecolors='none', edgecolors=color,
                           lw=1.0, transform=ccrs.PlateCarree(), label=name, zorder=6)
            else:
                ax.scatter(lo, la, s=6, color=color, transform=ccrs.PlateCarree(),
                           label=name, zorder=5)
        ax.set_title('Survey coverage')
        ax.legend(fontsize=8, loc='lower left', markerscale=2)
    else:
        ax = fig.add_subplot(gs[:, 0])
        for name, color, cols, la, lo, hd in data:
            ax.scatter(lo, la, s=6, color=color, label=name)
        ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
        ax.set_title('Survey coverage (install cartopy for basemap)')
        ax.legend(fontsize=8)

    # (b) Median profiles with IQR
    ax = fig.add_subplot(gs[:, 1])
    for name, color, cols, la, lo, hd in data:
        with np.errstate(invalid='ignore'):
            med = np.nanmedian(cols, axis=0)
            q1 = np.nanpercentile(cols, 25, axis=0)
            q3 = np.nanpercentile(cols, 75, axis=0)
            n = np.sum(np.isfinite(cols), axis=0)
        for arr in (med, q1, q3):
            arr[n < 5] = np.nan
        ax.fill_betweenx(DEPTH, q1, q3, color=color, alpha=0.15)
        ax.plot(med, DEPTH, color=color, lw=2, label=name)
    ax.invert_yaxis()
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel(r'$\Delta\lambda = \lambda_{cross} - \lambda_{along}$')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(-0.11, 0.11)
    ax.grid(alpha=0.3)
    ax.set_title('Median profiles (IQR shaded)')
    ax.legend(fontsize=8, loc='lower right')

    # (c) Drive headings mod 180 (dlam is symmetric under line reversal)
    for row, chunk in enumerate((data[:3], data[3:])):
        ax = fig.add_subplot(gs[row, 2])
        for name, color, cols, la, lo, hd in chunk:
            if hd.size:
                ax.hist(hd, bins=np.arange(0, 181, 10), density=True,
                        histtype='step', lw=1.8, color=color, label=name)
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
        if row == 1:
            ax.set_xlabel('Drive heading mod 180 (deg from north)')
        ax.set_ylabel('density', fontsize=8)
    fig.axes[-2].set_title('Drive orientations')

    fig.suptitle('EAGER polarimetric traverses 2022-2026: coverage, fabric contrast, '
                 'and drive orientation\n(profiles are survey-frame relative: '
                 r'$\Delta\lambda$ mixes fabric strength and line-to-fabric angle)',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'summary_map_profiles.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

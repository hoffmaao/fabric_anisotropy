"""Regional map views of each survey area.

One zoomed map panel per named survey region (Thwaites, Ridge A, Taylor
Dome, McMurdo, WAIS Divide/Kamb) plus a full-Antarctica locator. Blocks
are assigned to regions by nearest reference site, so multi-region
surveys (2023-24 covers Thwaites, WAIS Divide, and McMurdo) appear in
every panel they visited. Requires cartopy.

Usage: python region_maps.py [fabric_batch_root] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')

SURVEYS = [
    ('2022-23', f'{ROOT}/joint/2022_Antarctica_Ground/*/Data_*.mat', 'tab:gray'),
    ('2023-24', f'{ROOT}/joint/2023_Antarctica_Ground/*/Data_*.mat', 'tab:blue'),
    ('2024-25 (Lilien proc.)', f'{ROOT}/joint/2024_Antarctica_Ground2/*/Data_*.mat', 'tab:orange'),
    ('2024-25 (Paden proc.)', f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat', 'tab:red'),
    ('2025-26', f'{ROOT}/joint/2025_Antarctica_Ground2/*/Data_*.mat', 'tab:green'),
]

# Region reference sites (assignment targets and star markers)
REGIONS = [
    ('Thwaites', -76.46, -105.5, ('Thwaites E shear margin', -76.46, -105.5)),
    ('WAIS Divide / Kamb', -79.42, -111.9, ('WAIS Divide core', -79.468, -112.086)),
    ('McMurdo', -77.70, 168.1, ('McMurdo Station', -77.85, 166.67)),
    ('Taylor Dome', -77.78, 158.75, ('Taylor Dome summit', -77.68, 157.72)),
    ('Ridge A', -86.60, 69.0, ('South Pole', -90.0, 0.0)),
]


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def load_positions(pattern):
    lats, lons = [], []
    for fn in sorted(glob.glob(pattern)):
        d = loadmat(fn)
        lats.extend(d['Latitude'][0, :])
        lons.extend(d['Longitude'][0, :])
    return np.array(lats), np.array(lons)


def main():
    surveys = []
    for name, pattern, color in SURVEYS:
        la, lo = load_positions(pattern)
        if la.size:
            surveys.append((name, color, la, lo))
    if not surveys:
        print('no data found under', ROOT)
        return

    # Assign every block to the nearest region reference
    ref_la = np.array([r[1] for r in REGIONS])
    ref_lo = np.array([r[2] for r in REGIONS])
    assign = {}
    for name, color, la, lo in surveys:
        d = np.stack([haversine_km(la, lo, rla, rlo)
                      for rla, rlo in zip(ref_la, ref_lo)])
        assign[name] = np.argmin(d, axis=0)

    fig = plt.figure(figsize=(17, 10))

    for ri, (rname, rla, rlo, site) in enumerate(REGIONS):
        ax = fig.add_subplot(2, 3, ri + 1,
                             projection=ccrs.SouthPolarStereo(central_longitude=rlo))
        pts_la, pts_lo = [], []
        for name, color, la, lo in surveys:
            m = assign[name] == ri
            if not np.any(m):
                continue
            if 'Lilien' in name:
                ax.scatter(lo[m], la[m], s=40, facecolors='none', edgecolors=color,
                           lw=0.8, transform=ccrs.PlateCarree(), label=name, zorder=6)
            else:
                ax.scatter(lo[m], la[m], s=4, color=color,
                           transform=ccrs.PlateCarree(), label=name, zorder=5)
            pts_la.extend(la[m]); pts_lo.extend(lo[m])
        if not pts_la:
            ax.set_title(f'{rname} (no blocks)')
            continue
        pts_la, pts_lo = np.array(pts_la), np.array(pts_lo)
        pad_la = max(0.15, 0.35 * np.ptp(pts_la))
        pad_lo = max(0.5, 0.35 * np.ptp(pts_lo))
        ax.set_extent([pts_lo.min() - pad_lo, pts_lo.max() + pad_lo,
                       pts_la.min() - pad_la, pts_la.max() + pad_la],
                      ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale('10m'), lw=0.6)
        gl = ax.gridlines(draw_labels=True, lw=0.3, color='gray', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        sname, sla, slo = site
        ax.plot(slo, sla, marker='*', ms=12, color='k',
                transform=ccrs.PlateCarree(), zorder=7)
        ax.annotate(sname, xy=(0.03, 0.03), xycoords='axes fraction', fontsize=7)
        ax.set_title(rname, fontsize=12)
        ax.legend(fontsize=7, loc='upper right', markerscale=2)

    # Locator panel
    ax = fig.add_subplot(2, 3, 6, projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180, 180, -90, -63], ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), lw=0.5)
    for rname, rla, rlo, site in REGIONS:
        ax.plot(rlo, rla, 's', ms=7, mfc='none', mec='crimson', mew=1.5,
                transform=ccrs.PlateCarree())
        ax.annotate(rname, ccrs.SouthPolarStereo().transform_point(rlo, rla,
                    ccrs.PlateCarree())[:2], fontsize=8, xytext=(4, 4),
                    textcoords='offset points')
    ax.set_title('Survey regions', fontsize=12)

    fig.suptitle('EAGER polarimetric fabric surveys: regional coverage', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'region_maps.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

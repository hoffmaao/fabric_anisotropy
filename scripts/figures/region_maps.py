"""Regional map views of each survey area.

One zoomed map panel per named survey region (Thwaites, Ridge A, Taylor
Dome, McMurdo, WAIS Divide/Kamb) plus a full-Antarctica locator. Blocks
are assigned to regions by nearest reference site, so multi-region
surveys (2023-24 covers Thwaites, WAIS Divide, and McMurdo) appear in
every panel they visited. Panels are EPSG:3031 with LIMA/MOA imagery
backgrounds (see antarctic_basemap; coastline-only when the mosaics are
absent), lat/lon graticule labels, and scale bars. Requires cartopy.

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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab
from scar_style import FIGS

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else FIGS

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
        finite = np.isfinite(la) & np.isfinite(lo)
        la, lo = la[finite], lo[finite]
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
    proj = ab.proj3031()

    for ri, (rname, rla, rlo, site) in enumerate(REGIONS):
        ax = fig.add_subplot(2, 3, ri + 1, projection=proj)
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
        # Projected extents sidestep antimeridian wrap entirely
        extent = ab.points_extent(pts_la, pts_lo)
        ax.set_extent(extent, crs=proj)
        has_img = ab.add_imagery(ax, extent)
        ax.add_feature(cfeature.COASTLINE.with_scale('10m'), lw=0.6,
                       edgecolor='k' if not has_img else 'yellow')
        # Edge (not inline) latitude labels: inline ones land anywhere in
        # the panel and can collide with the scale bar or annotations
        gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.6, y_inline=False,
                          color='gray' if not has_img else 'white')
        gl.top_labels = False
        gl.right_labels = False
        sname, sla, slo = site
        ax.plot(slo, sla, marker='*', ms=12, color='k', mec='white', mew=0.8,
                transform=ccrs.PlateCarree(), zorder=7)
        sx, sy = proj.transform_point(slo, sla, ccrs.PlateCarree())[:2]
        if not (extent[0] <= sx <= extent[1] and extent[2] <= sy <= extent[3]):
            sname = f'nearest reference: {sname}'
        # Lower-right corner, above the scale bar region (bar + label reach
        # ~0.12 axes height) and below the legend (upper right). Narrow
        # portrait panels can't fit the label on one line, so wrap on the
        # colon and shrink it there.
        narrow = (extent[1] - extent[0]) < (extent[3] - extent[2])
        if narrow:
            sname = sname.replace(': ', ':\n')
        ax.annotate(sname, xy=(0.97, 0.14), xycoords='axes fraction',
                    fontsize=6 if narrow else 7, ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none',
                              alpha=0.7))
        ab.add_scale_bar(ax, extent)
        ax.set_title(rname, fontsize=12)
        ax.legend(fontsize=7, loc='upper right', markerscale=2)

    # Locator panel
    ax = fig.add_subplot(2, 3, 6, projection=proj)
    loc_extent = ab.ANTARCTICA_EXTENT
    ax.set_extent(loc_extent, crs=proj)
    has_img = ab.add_imagery(ax, loc_extent, max_px=900)
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), lw=0.5)
    # Thwaites/WAIS Divide and McMurdo/Taylor Dome markers sit close
    # together; nudge each pair's labels apart vertically
    label_offsets = {'Thwaites': (4, 10), 'WAIS Divide / Kamb': (4, -16),
                     'Taylor Dome': (4, 8), 'McMurdo': (4, -16)}
    for rname, rla, rlo, site in REGIONS:
        ax.plot(rlo, rla, 's', ms=7, mfc='none', mec='crimson', mew=1.5,
                transform=ccrs.PlateCarree())
        ax.annotate(rname, proj.transform_point(rlo, rla,
                    ccrs.PlateCarree())[:2], fontsize=8,
                    xytext=label_offsets.get(rname, (4, 4)),
                    textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none',
                              alpha=0.7))
    ab.add_scale_bar(ax, loc_extent)
    ax.set_title('Survey regions', fontsize=12)

    fig.suptitle('EAGER polarimetric fabric surveys: regional coverage', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'region_maps.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

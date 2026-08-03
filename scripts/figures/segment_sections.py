"""Along-profile depth sections of the fabric solution per survey region.

One curtain panel per region (Thwaites, WAIS Divide/Kamb, McMurdo, Taylor
Dome, Ridge A): the depth-variable dlam solution for every along-track
block, ordered by acquisition time along the drive, with distance along
the profile on the x-axis, segment boundaries ticked, and >2 km spatial
gaps shown as breaks. The 2024-25 season uses the Paden processing (more
blocks); other regions take every survey that visited them.

Usage: python segment_sections.py [fabric_batch_root] [out_dir]
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
GAP_KM = 2.0

PRODUCTS = [
    ('2022-23', f'{ROOT}/joint/2022_Antarctica_Ground/*/Data_*.mat'),
    ('2023-24', f'{ROOT}/joint/2023_Antarctica_Ground/*/Data_*.mat'),
    ('2024-25 (Paden proc.)', f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat'),
    ('2025-26', f'{ROOT}/joint/2025_Antarctica_Ground2/*/Data_*.mat'),
]

REGIONS = [
    ('Thwaites', -76.46, -105.5),
    ('WAIS Divide / Kamb', -79.42, -111.9),
    ('McMurdo', -77.70, 168.1),
    ('Taylor Dome', -77.78, 158.75),
    ('Ridge A', -86.60, 69.0),
]


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def load_blocks(pattern, survey):
    """One dict per QC'd block with profile, position, time, and segment."""
    blocks = []
    for fn in sorted(glob.glob(pattern)):
        seg = os.path.basename(os.path.dirname(fn))
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
            lat, lon = d['Latitude'][0, b], d['Longitude'][0, b]
            gps = d['GPS_time'][0, b]
            if not (np.isfinite(lat) and np.isfinite(lon)
                    and np.isfinite(gps)):
                continue
            blocks.append(dict(col=col, lat=lat, lon=lon, gps=gps,
                               seg=seg, survey=survey))
    return blocks


def main():
    all_blocks = []
    for survey, pattern in PRODUCTS:
        all_blocks.extend(load_blocks(pattern, survey))
    if not all_blocks:
        print('no data found under', ROOT)
        return

    ref_la = np.array([r[1] for r in REGIONS])
    ref_lo = np.array([r[2] for r in REGIONS])
    for blk in all_blocks:
        d = haversine_km(blk['lat'], blk['lon'], ref_la, ref_lo)
        blk['region'] = int(np.argmin(d))

    fig, axes = plt.subplots(len(REGIONS), 1, figsize=(15, 14))
    # Room for each panel's tick labels under the next panel's title
    fig.subplots_adjust(hspace=0.45)
    pc = None
    for ri, (rname, _, _) in enumerate(REGIONS):
        ax = axes[ri]
        blks = sorted((b for b in all_blocks if b['region'] == ri),
                      key=lambda b: b['gps'])
        if not blks:
            ax.set_title(f'{rname}: no blocks')
            continue
        cols = np.array([b['col'] for b in blks]).T
        lats = np.array([b['lat'] for b in blks])
        lons = np.array([b['lon'] for b in blks])
        step = haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
        dist = np.concatenate([[0], np.cumsum(np.minimum(step, GAP_KM))])
        gap = np.concatenate([[False], step > GAP_KM])
        cols[:, gap] = np.nan

        edges = np.concatenate([[dist[0] - 0.1],
                                (dist[:-1] + dist[1:]) / 2, [dist[-1] + 0.1]])
        depth_edges = np.concatenate([DEPTH, [DEPTH[-1] + 5]])
        pc = ax.pcolormesh(edges, depth_edges, cols, cmap='RdBu_r',
                           vmin=-0.1, vmax=0.1)

        # Segment boundaries and coherent depth range
        segs = [b['seg'] for b in blks]
        n_seg = 1
        for i in range(1, len(segs)):
            if segs[i] != segs[i - 1]:
                n_seg += 1
                ax.axvline(edges[i], color='k', lw=0.4, alpha=0.5)
        finite_rows = np.where(np.any(np.isfinite(cols), axis=1))[0]
        ymax = DEPTH[finite_rows[-1]] + 100 if finite_rows.size else DEPTH[-1]
        ax.set_ylim(ymax, 0)
        surveys = sorted({b['survey'] for b in blks})
        ax.set_title(f"{rname} - {', '.join(surveys)} "
                     f"({len(blks)} blocks, {n_seg} segments)", fontsize=10)
        ax.set_ylabel('Depth (m)')
    axes[-1].set_xlabel('Distance along drive (km, gaps >2 km compressed)')
    if pc is not None:
        fig.colorbar(pc, ax=axes, shrink=0.6,
                     label=r'$\Delta\lambda = \lambda_{cross} - \lambda_{along}$')
    fig.suptitle('Depth-variable fabric solutions along survey profiles '
                 '(joint inversion; thin vertical lines mark segment boundaries)',
                 fontsize=13)
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'segment_sections.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

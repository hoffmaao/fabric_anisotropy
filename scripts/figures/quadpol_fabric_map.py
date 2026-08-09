"""Map view of the quad-pol LS fabric result over a survey grid.

Two layers per survey line. STRENGTH runs along the track: consecutive
125-trace blocks are drawn as line segments coloured by the deep
principal contrast (median LS sec_dlam_ls over 1150-1500 m), single-hue
sequential, with abstained blocks left as gaps over the grey track.
ORIENTATION is one fixed-length ink bar per frame: the LS fabric axis
theta0 expressed geographically (circular median over the 200-1200 m
windows plus the frame's track azimuth). Strength lives in colour and
orientation in the bars, so neither is double-encoded. SCAR style
otherwise: no legend box, a plain black scale bar, a projected north
arrow, and a small horizontal colorbar as the sections use.

Bar geometry is built geodetically - the endpoints are offset from the
midpoint along the azimuth in lat/lon and only then projected - so the
rendered orientation is correct at any longitude with no convergence
arithmetic.

Reads every quadpol_section_<prefix>*.mat under the SCAR staging dir
(the run_quadpol_pipeline.m outputs mirrored locally). The tag prefix
selects the survey: 2025 is Ridge A, 2024 is Thwaites.

Usage: python quadpol_fabric_map.py <out_dir> [tag_prefix] [title]
"""
import glob
import os
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
from pyproj import Transformer           # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import DATA, INK, MUTED  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
PREFIX = sys.argv[2] if len(sys.argv) > 2 else '2025'
TITLE = sys.argv[3] if len(sys.argv) > 3 else 'Ridge A grid'
BLUE = '#2a78d6'
Z_BAND = (200.0, 1200.0)      # theta0 window band
Z_DEEP = (1150.0, 1500.0)     # contrast band the colour encodes
BAR_KM = 2.2                  # fixed bar half-length [km], orientation only
VMAX = 0.08                   # colour scale top; sequential from 0

T = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)
Ti = Transformer.from_crs('EPSG:3031', 'EPSG:4326', always_xy=True)


def circmedian_axis(deg):
    """Circular mean of axis-valued angles (mod 180), in degrees."""
    ph = np.exp(2j * np.deg2rad(deg[np.isfinite(deg)]))
    return np.rad2deg(np.angle(np.mean(ph))) / 2 % 180


def main():
    frames = []
    for fn in sorted(glob.glob(os.path.join(
            DATA, 'quadpol_section_%s*.mat' % PREFIX))):
        with h5py.File(fn) as f:
            r = f['res']

            def g(k):
                return float(np.array(r[k]).ravel()[0])

            zw = np.array(r['ls_zw']).ravel()
            th0 = np.array(r['ls_theta0']).ravel()
            track = g('track_az')
            z = np.array(r['z']).ravel()
            secl = np.array(r['sec_dlam_ls']).T
            blat = np.array(r['sec_lat']).ravel()
            blon = np.array(r['sec_lon']).ravel()
            m = (zw > Z_BAND[0]) & (zw < Z_BAND[1]) & np.isfinite(th0)
            mz = (z > Z_DEEP[0]) & (z < Z_DEEP[1])
            if m.sum() < 5:
                continue
            frames.append(dict(
                th=circmedian_axis(th0[m] + track),
                blk_dl=np.nanmedian(secl[mz], axis=0),
                blat=blat, blon=blon,
                lat=g('lat'), lon=g('lon'), lat0=g('lat0'), lon0=g('lon0'),
                lat1=g('lat1'), lon1=g('lon1')))
    if not frames:
        raise SystemExit('no quadpol_section_%s*.mat with usable theta0 '
                         'under %s' % (PREFIX, DATA))
    print('%d frames' % len(frames))

    from matplotlib.collections import LineCollection
    cmap = plt.get_cmap('Blues')
    norm = plt.Normalize(0, VMAX)
    fig, ax = plt.subplots(figsize=(7.8, 7.2), layout='constrained')
    for fr in frames:
        x, y = T.transform([fr['lon0'], fr['lon1']], [fr['lat0'], fr['lat1']])
        ax.plot(x, y, '-', color=MUTED, lw=1.0, alpha=0.8, zorder=1)
        # strength along the line: segments between block centres coloured
        # by the mean of the two blocks; NaN blocks leave the grey track
        bx, by = T.transform(fr['blon'], fr['blat'])
        pts = np.column_stack([bx, by])
        segs, vals = [], []
        bd = fr['blk_dl']
        for k in range(len(bd) - 1):
            v = 0.5 * (bd[k] + bd[k + 1])
            if np.isfinite(v):
                segs.append([pts[k], pts[k + 1]])
                vals.append(v)
        if segs:
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=4.2,
                                capstyle='round', zorder=2)
            lc.set_array(np.asarray(vals))
            ax.add_collection(lc)

    def bar(lat, lon, th_deg, color=INK, lw=1.6):
        th = np.deg2rad(th_deg)
        dlat = BAR_KM * np.cos(th) / 111.2
        dlon = BAR_KM * np.sin(th) / (111.2 * np.cos(np.deg2rad(lat)))
        x, y = T.transform([lon - dlon, lon + dlon],
                           [lat - dlat, lat + dlat])
        ax.plot(x, y, '-', color=color, lw=lw, solid_capstyle='round',
                zorder=3)

    for fr in frames:
        if np.isfinite(fr['th']):
            bar(fr['lat'], fr['lon'], fr['th'])

    ax.set_aspect('equal')
    xs, ys = ax.get_xlim(), ax.get_ylim()
    pad = 0.06 * (xs[1] - xs[0])
    ax.set_xlim(xs[0] - pad, xs[1] + pad)
    ax.set_ylim(ys[0] - pad, ys[1] + pad)
    xs, ys = ax.get_xlim(), ax.get_ylim()

    mean_axis = circmedian_axis(np.array([fr['th'] for fr in frames]))
    xr = xs[0] + 0.15 * (xs[1] - xs[0])
    yr = ys[0] + 0.10 * (ys[1] - ys[0])
    rlon, rlat = Ti.transform(xr, yr)
    bar(rlat, rlon, mean_axis)
    ax.text(xr, yr - 0.045 * (ys[1] - ys[0]), 'fabric axis', ha='center',
            fontsize=8.5, color=INK)
    cax = ax.inset_axes([0.36, 0.045, 0.28, 0.022])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                 orientation='horizontal')
    cax.set_title(r'$\Delta\lambda$ (%d-%d m)' % Z_DEEP, fontsize=8.5,
                  color=INK, pad=3)
    cax.tick_params(labelsize=7.5)

    xb = xs[1] - 0.22 * (xs[1] - xs[0])
    yb = ys[0] + 0.09 * (ys[1] - ys[0])
    ax.plot([xb, xb + 5000], [yb, yb], '-', color=INK, lw=2.5)
    ax.text(xb + 2500, yb + 0.018 * (ys[1] - ys[0]), '5 km', ha='center',
            fontsize=8.5, color=INK)

    clat = np.mean([fr['lat'] for fr in frames])
    clon = np.mean([fr['lon'] for fr in frames])
    x0, y0 = T.transform(clon, clat)
    x1, y1 = T.transform(clon, clat + 0.05)
    nv = np.array([x1 - x0, y1 - y0])
    nv /= np.hypot(*nv)
    xa = xs[0] + 0.10 * (xs[1] - xs[0])
    ya = ys[1] - 0.13 * (ys[1] - ys[0])
    L = 0.055 * (xs[1] - xs[0])
    ax.annotate('', xy=(xa + nv[0] * L, ya + nv[1] * L), xytext=(xa, ya),
                arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.4))
    ax.text(xa + nv[0] * L * 1.45, ya + nv[1] * L * 1.45, 'N', fontsize=10,
            color=INK, ha='center', va='center')

    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title('%s: horizontal fabric axis and contrast from single-pass\n'
                 r'quad-pol (LS fit) - bars: axis $\theta_0$; track colour: '
                 r'$\Delta\lambda$ per 125-trace block'
                 % TITLE, fontsize=10.5, color=INK)
    out = os.path.join(OUT, 'scar_quadpol_fabric_map_%s.png' % PREFIX)
    fig.savefig(out, dpi=200, facecolor='white', bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()

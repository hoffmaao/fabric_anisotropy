"""Fabric strength through the column: one map per depth window, animated.

The map figures fix a deep band (1150-1500 m) and show how fabric varies
laterally. This sweeps that window DOWN the column instead, so the survey
map is redrawn every 20 m and the movie shows how the fabric strengthens
with depth across the whole grid at once.

Everything except the depth window is held fixed, deliberately: the same
colour scale on every frame, the same extent, the same graticule. A frame
that autoscaled its colours would make every depth look equally
anisotropic, which is the one thing this movie exists to disprove.

The depth axis is the LS estimator's own sec_dlam_ls, per 125-trace block,
so a block enters a frame only where that window holds enough finite cells
to support a median; blocks that fail leave the grey track showing rather
than being interpolated.

ON THIS SURVEY NOTHING FAILS: coverage is 2157 of 2157 blocks in EVERY
window from 200 m to 1450 m. So the pale shallow frames are MEASUREMENTS
OF A WEAK FABRIC, not gaps, and must not be read as missing data - which
is the opposite of the gap-banded behaviour the published direct chain
shows over the same frames, and one of the things the LS fit bought. The
measured profile, median over all blocks:

    200-350 m  0.013     600-750 m  0.062     1200-1350 m  0.064
    300-450 m  0.016     800-950 m  0.055     1300-1450 m  0.073
    400-550 m  0.040    1000-1150 m 0.064

so the contrast rises about fivefold down the column, with most of the
climb between 350 and 550 m and a mild reversal around 800-950 m. Note the
shallow half is where the co-polarized comparison is still unresolved at
about a factor of two, so treat the top few frames as the LS estimator's
own answer rather than as a cross-validated one.

Usage: python quadpol_depth_movie.py <out_dir> [site]
"""
import glob
import os
import re
import shutil
import subprocess
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.collections import LineCollection      # noqa: E402
from pyproj import Transformer                         # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab                         # noqa: E402
import quadpol_sites as qs                             # noqa: E402
from scar_style import DATA, INK, MUTED                # noqa: E402
from scar_style import FIGS  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else FIGS
SITE = sys.argv[2] if len(sys.argv) > 2 else 'ridge_a'

FPS, HOLD = 25, 6
BAR_FRAC = 0.035       # orientation bar half-length as a fraction of span

CFG = qs.get(SITE)
TITLE = CFG['title']
VMAX = CFG['vmax']
# Ridge A stops at 1050 m, not at the bottom of the record: below that its
# mostly N-S lines stop following the regular fabric cycle and the deeper
# windows would read as fabric development when they are something else.
# Every other site runs to the depth its own coverage supports.
Z_FIRST, Z_LAST = CFG['movie']
# window and step come from the site: see quadpol_sites for why a single
# window height cannot serve both a 1850 m column and a 300 m ice shelf
WIN_M = CFG.get('win', 150.0)
STEP_M = CFG.get('step', 20.0)
MIN_CELLS = 3 if WIN_M < 80.0 else 5

T = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)


def find_ffmpeg():
    for c in (shutil.which('ffmpeg'), '/opt/anaconda3/envs/ar-env/bin/ffmpeg',
              '/opt/homebrew/bin/ffmpeg', '/usr/local/bin/ffmpeg'):
        if c and os.path.exists(c):
            return c
    raise SystemExit('no ffmpeg found; frames are in '
                     '<out>/depth_frames_%s' % SITE)


def load():
    frames = []
    for fn in sorted(glob.glob(os.path.join(
            DATA, 'quadpol_section_*.mat'))):
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        if re.search(r'_z\d+$', tag):
            continue
        with h5py.File(fn) as f:
            r = f['res']
            z = np.array(r['z']).ravel()
            secl = np.array(r['sec_dlam_ls']).T
            blat = np.array(r['sec_lat']).ravel()
            blon = np.array(r['sec_lon']).ravel()
            zw = np.array(r['ls_zw']).ravel()
            th_geo = np.array(r['ls_theta0_geo']).ravel()
            lat = float(np.array(r['lat']).ravel()[0])
            lon = float(np.array(r['lon']).ravel()[0])
            lat0 = float(np.array(r['lat0']).ravel()[0])
            lon0 = float(np.array(r['lon0']).ravel()[0])
            lat1 = float(np.array(r['lat1']).ravel()[0])
            lon1 = float(np.array(r['lon1']).ravel()[0])
        if not qs.in_site(CFG, lat, lon, tag):
            continue
        ok = np.isfinite(blat) & np.isfinite(blon)
        if not ok.any():
            continue
        bx, by = T.transform(blon[ok], blat[ok])
        ex, ey = T.transform([lon0, lon1], [lat0, lat1])
        frames.append(dict(tag=tag, z=z, sec=secl[:, ok],
                           pts=np.column_stack([bx, by]),
                           zw=zw, th_geo=th_geo, lat=lat, lon=lon,
                           ex=ex, ey=ey))
    if not frames:
        raise SystemExit('no sections within %.0f deg of %s under %s'
                         % (qs.SITE_RADIUS_DEG, SITE, DATA))
    return frames


def main():
    frames = load()
    centres = np.arange(Z_FIRST + WIN_M / 2, Z_LAST - WIN_M / 2 + 1e-6,
                        STEP_M)
    print('%d frames, %d depth windows (%.0f m window, %.0f m step)'
          % (len(frames), centres.size, WIN_M, STEP_M))

    proj = ab.proj3031()
    cmap, norm = plt.get_cmap('Blues'), plt.Normalize(0.0, VMAX)
    allx = np.concatenate([f['pts'][:, 0] for f in frames])
    ally = np.concatenate([f['pts'][:, 1] for f in frames])
    sq, span_x, span_y = qs.square_extent(allx, ally)
    bar_km = BAR_FRAC * max(span_x, span_y) / 1000.0
    print('  survey span %.1f x %.1f km, bar %.2f km'
          % (span_x / 1e3, span_y / 1e3, bar_km))

    # per-site frame dir, cleared up front: leftover seq_/d* frames from a
    # crashed run (or another site) would be numbered after this run's and
    # ffmpeg would append them to the tail of the movie
    fdir = os.path.join(OUT, 'depth_frames_%s' % SITE)
    if os.path.isdir(fdir):
        shutil.rmtree(fdir)
    os.makedirs(fdir)
    n = 0
    for k, zc in enumerate(centres):
        lo, hi = zc - WIN_M / 2, zc + WIN_M / 2
        fig = plt.figure(figsize=(7.6, 7.6), layout='constrained')
        ax = fig.add_subplot(1, 1, 1, projection=proj)

        def arm(lat, lon, th_deg, half_km=None, color=INK, lw=1.6):
            half_km = bar_km if half_km is None else half_km
            th = np.deg2rad(th_deg)
            dlat = half_km * np.cos(th) / 111.2
            dlon = half_km * np.sin(th) / (111.2 * np.cos(np.deg2rad(lat)))
            x, y = T.transform([lon - dlon, lon + dlon],
                               [lat - dlat, lat + dlat])
            ax.plot(x, y, '-', color=color, lw=lw, solid_capstyle='round',
                    zorder=4, transform=proj)

        cov = 0
        for fr in frames:
            ax.plot(fr['ex'], fr['ey'], '-', color=MUTED, lw=1.0, alpha=0.8,
                    zorder=1, transform=proj)
            mz = (fr['z'] > lo) & (fr['z'] < hi)
            sub = fr['sec'][mz]
            bd = np.where(np.isfinite(sub).sum(axis=0) >= MIN_CELLS,
                          np.nanmedian(np.where(np.isfinite(sub), sub,
                                                np.nan), axis=0), np.nan)
            cov += int(np.isfinite(bd).sum())
            segs, vals = [], []
            for j in range(len(bd) - 1):
                v = 0.5 * (bd[j] + bd[j + 1])
                if np.isfinite(v):
                    segs.append([fr['pts'][j], fr['pts'][j + 1]])
                    vals.append(v)
            if segs:
                lc = LineCollection(segs, cmap=cmap, norm=norm,
                                    linewidths=4.2, capstyle='round',
                                    zorder=2, transform=proj)
                lc.set_array(np.asarray(vals))
                ax.add_collection(lc)

        # Orientation AT THIS DEPTH, not the frame average: theta0 comes
        # from the LS windows inside this same slice, combined as a
        # doubled-angle phasor (never an unwrapped angle - an unwrap over
        # gappy axis samples can slip a branch and hand a silent 90 deg
        # error to the bar). A frame whose windows all abstained here draws
        # no bar, so the axis appears only where it is actually resolved.
        th_here = []
        for fr in frames:
            mw = (fr['zw'] > lo) & (fr['zw'] < hi) & np.isfinite(fr['th_geo'])
            if mw.sum() < 2:
                continue
            ph = np.mean(np.exp(2j * np.deg2rad(fr['th_geo'][mw])))
            if not np.isfinite(ph) or abs(ph) < 1e-9:
                continue
            th = np.rad2deg(np.angle(ph)) / 2 % 180
            arm(fr['lat'], fr['lon'], th)
            th_here.append(th)

        ax.set_xlim(sq[0], sq[1])
        ax.set_ylim(sq[2], sq[3])
        xs, ys = ax.get_xlim(), ax.get_ylim()
        gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='gray',
                          y_inline=False, x_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 8.5, 'color': INK}
        gl.ylabel_style = {'size': 8.5, 'color': INK}

        cax = ax.inset_axes([0.36, 0.045, 0.28, 0.022])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                     orientation='horizontal')
        cax.set_title(r'$\Delta\lambda$', fontsize=8.5, color=INK, pad=3)
        cax.tick_params(labelsize=7.5)

        sb_m = ab._nice_length(0.22 * max(span_x, span_y))
        xb = xs[1] - 0.06 * (xs[1] - xs[0]) - sb_m
        yb = ys[0] + 0.07 * (ys[1] - ys[0])
        ax.plot([xb, xb + sb_m], [yb, yb], '-', color=INK, lw=2.5,
                transform=proj)
        ax.text(xb + sb_m / 2, yb - 0.028 * (ys[1] - ys[0]),
                '%g km' % (sb_m / 1000.0) if sb_m >= 1000 else '%g m' % sb_m,
                ha='center', va='top', fontsize=8.5, color=INK,
                transform=proj)

        # Depth readout only. A survey-wide circular mean of theta0 was
        # drawn here too, and a bar tracking the window down the column;
        # both are gone. The mean angle averaged 36 frames whose axes are
        # a real spatial field, not repeat measurements of one number, so
        # its spread read as uncertainty when it is mostly variation - and
        # the bars already show the axis where it is measured. The column
        # bar restated what the depth text says.
        ax.text(xs[0] + 0.04 * (xs[1] - xs[0]),
                ys[1] - 0.055 * (ys[1] - ys[0]),
                '%.0f - %.0f m' % (lo, hi), fontsize=13, color=INK,
                ha='left', va='center', transform=proj)

        ax.set_title('%s: fabric contrast through the column\n'
                     r'track colour: $\Delta\lambda$ per block; '
                     r'bars: fabric axis $\theta_0$ at this depth' % TITLE,
                     fontsize=10.5, color=INK)
        png = os.path.join(fdir, 'd%04d.png' % k)
        fig.savefig(png, dpi=150, facecolor='white')
        plt.close(fig)
        for r in range(1, HOLD):
            os.link(png, os.path.join(fdir, 'd%04d_%03d.png' % (k, r)))
        n += HOLD
        print('  [%2d/%d] %.0f-%.0f m, %d blocks, %d axes'
              % (k + 1, centres.size, lo, hi, cov, len(th_here)),
              flush=True)

    seq = os.path.join(fdir, 'seq_%05d.png')
    for i, p in enumerate(sorted(glob.glob(os.path.join(fdir, 'd*.png')))):
        os.rename(p, seq % i)
    mp4 = os.path.join(OUT, 'quadpol_depth_movie_%s.mp4' % SITE)
    subprocess.run([find_ffmpeg(), '-y', '-framerate', str(FPS), '-i', seq,
                    '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20',
                    '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', mp4],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    shutil.rmtree(fdir)
    print('wrote %s (%d frames, %.1f s)' % (mp4, n, n / FPS))


if __name__ == '__main__':
    main()

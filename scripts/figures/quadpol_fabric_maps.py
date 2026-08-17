"""Quad-pol fabric maps: principal, grid-north projection, or both axes.

One script, three modes, so the survey geometry, the axes furniture and the
colour rules cannot drift apart between them:

  principal  track coloured by dlam, the axis-resolved contrast, sequential
  gridnorth  track coloured by the eigenvalue difference PROJECTED onto GRID
             north, signed, diverging about zero
  cross      both horizontal principal axes drawn, arm-length difference
             carrying dlam, track coloured by dlam

GRID NORTH, NOT TRUE NORTH. On EPSG:3031 the meridian convergence equals
the longitude, which at Ridge A is ~68.8 deg - so a projection onto true
north and one onto grid north are nearly 69 deg apart and DISAGREE IN SIGN
there (true north -0.065, grid north +0.057), because theta0 sits ~15 deg
off grid north but ~96 deg off true north. The convergence also varies
across a survey (~3.7 deg over Ridge A's 25 km, since a degree of longitude
is only 6.7 km at 86.6 S), so it is evaluated PER BLOCK from that block's
own position rather than once at the survey centre. It is measured
numerically - step north a little and see which way that points in grid
coordinates - rather than assumed equal to the longitude, so the same code
is correct on any projection the basemap module might be switched to.

NO NORTH ARROW. The graticule carries the orientation instead, drawn and
labelled the way the two-panel site figures do it, so these maps and those
read the same way. An arrow would also have to say WHICH north it meant,
which is exactly the ambiguity that motivated the grid-north change.

The projection law, verified against the co-polarized ribbons across all 36
Ridge A frames (corr +0.989, 100% sign agreement):

    lam(a) - lam(a + 90) = dlam * cos 2(a - theta0)

Usage: python quadpol_fabric_maps.py <out_dir> <mode> [site]
"""
import glob
import os
import re
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import cartopy.crs as ccrs                             # noqa: E402
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.collections import LineCollection      # noqa: E402
from matplotlib.colors import TwoSlopeNorm             # noqa: E402
from pyproj import Transformer                         # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab                         # noqa: E402
from scar_style import DATA, INK, MUTED                # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
MODE = sys.argv[2] if len(sys.argv) > 2 else 'principal'
SITE = sys.argv[3] if len(sys.argv) > 3 else 'ridge_a'
if MODE not in ('principal', 'gridnorth', 'cross'):
    raise SystemExit("mode must be principal, gridnorth or cross")

import quadpol_sites as qs                              # noqa: E402

CFG = qs.get(SITE)
TITLE = CFG['title']
Z_BAND = CFG['z_band']
Z_DEEP = CFG['z_deep']
VMAX = CFG['vmax']
VLIM = 0.09
# Marker sizes are FRACTIONS OF THE SURVEY, not fixed kilometres. A 2.2 km
# bar is right on Ridge A's ~25 km grid and absurd on Taylor Dome's ~3 km
# wide strip, where fixed-length bars sprawl clear across the survey and
# hide the tracks they annotate.
BAR_FRAC = 0.040              # orientation bar half-length, fraction of span
CROSS_FRAC = 0.022            # MEAN cross arm half-length; carries no meaning
GAIN_FRAC = 0.36              # arm-length difference per unit dlam, of span

T = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)
Ti = Transformer.from_crs('EPSG:3031', 'EPSG:4326', always_xy=True)


def circmedian_axis(deg):
    ph = np.exp(2j * np.deg2rad(deg[np.isfinite(deg)]))
    return np.rad2deg(np.angle(np.mean(ph))) / 2 % 180


def load():
    frames = []
    for fn in sorted(glob.glob(os.path.join(
            DATA, 'quadpol_section_*.mat'))):
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        if re.search(r'_z\d+$', tag):
            continue
        with h5py.File(fn) as f:
            r = f['res']

            def g(k):
                return float(np.array(r[k]).ravel()[0])

            if not qs.in_site(CFG, g('lat'), g('lon'), tag):
                continue
            zw = np.array(r['ls_zw']).ravel()
            th_geo = np.array(r['ls_theta0_geo']).ravel()
            z = np.array(r['z']).ravel()
            secl = np.array(r['sec_dlam_ls']).T
            blat = np.array(r['sec_lat']).ravel()
            blon = np.array(r['sec_lon']).ravel()
            m = (zw > Z_BAND[0]) & (zw < Z_BAND[1]) & np.isfinite(th_geo)
            mz = (z > Z_DEEP[0]) & (z < Z_DEEP[1])
            if m.sum() < 5:
                continue
            th = circmedian_axis(th_geo[m])
            dl = np.nanmedian(secl[mz], axis=0)
            ok = np.isfinite(blat) & np.isfinite(blon)
            az_gn = np.full(blat.shape, np.nan)
            if ok.any():
                az_gn[ok] = qs.grid_north_az(blat[ok], blon[ok])
            frames.append(dict(
                tag=tag, th=th, dl=dl,
                proj=dl * np.cos(2 * np.deg2rad(az_gn - th)),
                dl_frame=float(np.nanmedian(dl)),
                blat=blat, blon=blon,
                lat=g('lat'), lon=g('lon'), lat0=g('lat0'), lon0=g('lon0'),
                lat1=g('lat1'), lon1=g('lon1')))
    if not frames:
        raise SystemExit('no sections with usable theta0 within %.0f deg of '
                         '%s (%.2f, %.2f) under %s'
                         % (qs.SITE_RADIUS_DEG, SITE, CFG['lat'],
                            CFG['lon'], DATA))
    return frames


def main():
    frames = load()
    key = 'proj' if MODE == 'gridnorth' else 'dl'
    allv = np.concatenate([f[key][np.isfinite(f[key])] for f in frames])
    print('%s / %s: %d frames, %d blocks' % (SITE, MODE, len(frames),
                                             allv.size))
    print('  value: median %+.4f  p10/p90 %+.4f/%+.4f'
          % (np.median(allv), *np.percentile(allv, [10, 90])))
    gn = qs.grid_north_az([np.mean([f['lat'] for f in frames])],
                       [np.mean([f['lon'] for f in frames])])[0]
    print('  grid north at survey centre: %.2f deg E of true north' % gn)

    proj = ab.proj3031()
    if MODE == 'gridnorth':
        cmap, norm = plt.get_cmap('RdBu_r'), TwoSlopeNorm(0.0, -VLIM, VLIM)
    else:
        cmap, norm = plt.get_cmap('Blues'), plt.Normalize(0.0, VMAX)

    # Figure shape follows the survey. A square canvas around a long thin
    # strip wastes most of the slide and shrinks the tracks to a thread.
    bx0 = np.concatenate([T.transform(f['blon'], f['blat'])[0]
                          for f in frames])
    by0 = np.concatenate([T.transform(f['blon'], f['blat'])[1]
                          for f in frames])
    sq_extent, span_x, span_y = qs.square_extent(bx0, by0)
    fw = fh = 7.6
    span = max(span_x, span_y)
    bar_km = BAR_FRAC * span / 1000.0
    cross_km = CROSS_FRAC * span / 1000.0
    gain_km = GAIN_FRAC * span / 1000.0
    print('  survey span %.1f x %.1f km, bar %.2f km' %
          (span_x / 1e3, span_y / 1e3, bar_km))

    fig = plt.figure(figsize=(fw, fh), layout='constrained')
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    for fr in frames:
        x, y = T.transform([fr['lon0'], fr['lon1']], [fr['lat0'], fr['lat1']])
        ax.plot(x, y, '-', color=MUTED, lw=1.0, alpha=0.8, zorder=1,
                transform=proj)
        bx, by = T.transform(fr['blon'], fr['blat'])
        pts = np.column_stack([bx, by])
        segs, vals = [], []
        bd = fr[key]
        for k in range(len(bd) - 1):
            v = 0.5 * (bd[k] + bd[k + 1])
            if np.isfinite(v):
                segs.append([pts[k], pts[k + 1]])
                vals.append(v)
        if segs:
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=4.2,
                                capstyle='round', zorder=2, transform=proj)
            lc.set_array(np.asarray(vals))
            ax.add_collection(lc)

    def arm(lat, lon, th_deg, half_km, color=INK, lw=1.6):
        th = np.deg2rad(th_deg)
        dlat = half_km * np.cos(th) / 111.2
        dlon = half_km * np.sin(th) / (111.2 * np.cos(np.deg2rad(lat)))
        x, y = T.transform([lon - dlon, lon + dlon],
                           [lat - dlat, lat + dlat])
        ax.plot(x, y, '-', color=color, lw=lw, solid_capstyle='round',
                zorder=3, transform=proj)

    def cross(lat, lon, th_deg, dl, color=INK, lw=1.6):
        # clamped so a large dlam cannot drive the short arm through zero,
        # which would draw a cross whose short axis lay along theta0
        d = 0.5 * gain_km * float(np.clip(dl, 0.0, 2 * cross_km / gain_km))
        arm(lat, lon, th_deg, cross_km + d, color, lw)
        arm(lat, lon, th_deg + 90.0, max(cross_km - d, 0.02), color, lw)

    for fr in frames:
        if not np.isfinite(fr['th']):
            continue
        if MODE == 'cross' and np.isfinite(fr['dl_frame']):
            cross(fr['lat'], fr['lon'], fr['th'], fr['dl_frame'])
        elif MODE != 'cross':
            arm(fr['lat'], fr['lon'], fr['th'], bar_km)

    ax.set_xlim(sq_extent[0], sq_extent[1])
    ax.set_ylim(sq_extent[2], sq_extent[3])
    xs, ys = ax.get_xlim(), ax.get_ylim()

    # The graticule replaces the north arrow: it says which way north is
    # AND where the survey sits, and matches the two-panel site figures.
    gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='gray',
                      y_inline=False, x_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8.5, 'color': INK}
    gl.ylabel_style = {'size': 8.5, 'color': INK}

    mean_axis = circmedian_axis(np.array([fr['th'] for fr in frames]))
    if MODE == 'cross':
        for i, dref in enumerate((0.03, 0.08)):
            xr = xs[0] + (0.075 + 0.105 * i) * (xs[1] - xs[0])
            yr = ys[0] + 0.165 * (ys[1] - ys[0])
            rlon, rlat = Ti.transform(xr, yr)
            cross(rlat, rlon, mean_axis, dref, MUTED, 1.4)
            ax.text(xr, yr - 0.075 * (ys[1] - ys[0]),
                    r'$\Delta\lambda$=%.2f' % dref, ha='center', fontsize=8,
                    color=INK, transform=proj)
        ax.text(xs[0] + 0.128 * (xs[1] - xs[0]),
                ys[0] + 0.045 * (ys[1] - ys[0]), 'arm-length difference',
                ha='center', fontsize=8, color=INK, transform=proj)
    else:
        xr = xs[0] + 0.13 * (xs[1] - xs[0])
        yr = ys[0] + 0.11 * (ys[1] - ys[0])
        rlon, rlat = Ti.transform(xr, yr)
        arm(rlat, rlon, mean_axis, bar_km)
        ax.text(xr, yr - 0.045 * (ys[1] - ys[0]), 'fabric axis', ha='center',
                fontsize=8.5, color=INK, transform=proj)

    lab = (r'$\lambda_{\rm gridN} - \lambda_{\rm gridE}$ (%d-%d m)' % Z_DEEP
           if MODE == 'gridnorth' else r'$\Delta\lambda$ (%d-%d m)' % Z_DEEP)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                      orientation='horizontal', location='bottom',
                      shrink=0.62, pad=0.04, aspect=32)
    cb.set_label(lab, fontsize=9, color=INK)
    # cap the tick count: matplotlib's default put 9 five-decimal labels on
    # the 0-0.02 sites and they ran into each other unreadably
    from matplotlib.ticker import MaxNLocator
    cb.locator = MaxNLocator(nbins=5)
    cb.update_ticks()
    cb.ax.tick_params(labelsize=8, colors=INK)

    # scale bar length follows the survey too, rounded to something sayable
    sb_m = ab._nice_length(0.22 * span)
    xb = xs[1] - 0.06 * (xs[1] - xs[0]) - sb_m
    yb = ys[0] + 0.06 * (ys[1] - ys[0])
    ax.plot([xb, xb + sb_m], [yb, yb], '-', color=INK, lw=2.5, transform=proj)
    ax.text(xb + sb_m / 2, yb - 0.030 * (ys[1] - ys[0]),
            '%g km' % (sb_m / 1000.0) if sb_m >= 1000 else '%g m' % sb_m,
            ha='center', va='top', fontsize=8.5, color=INK, transform=proj)

    sub = {
        'principal': (r'bars: axis $\theta_0$; track colour: '
                      r'$\Delta\lambda$ per 125-trace block'),
        'gridnorth': (r'bars: axis $\theta_0$; track colour: '
                      r'$\Delta\lambda\,\cos 2(\alpha_{\rm gridN}-\theta_0)$ '
                      r'per block'),
        'cross': (r'cross arms: long along $\theta_0$, short across it, '
                  r'their DIFFERENCE $\propto\Delta\lambda$'),
    }[MODE]
    head = {
        'principal': 'horizontal fabric axis and contrast from single-pass '
                     'quad-pol (LS fit)',
        'gridnorth': 'horizontal eigenvalue difference projected onto GRID '
                     'north',
        'cross': 'both horizontal principal axes from single-pass quad-pol '
                 '(LS fit)',
    }[MODE]
    ax.set_title('%s: %s\n%s' % (TITLE, head, sub), fontsize=10.5, color=INK)
    out = os.path.join(OUT, 'scar_quadpol_fabric_%s_%s.png' % (MODE, SITE))
    fig.savefig(out, dpi=200, facecolor='white', bbox_inches='tight')
    print('  wrote', out)


if __name__ == '__main__':
    main()

"""Per-segment scene: setting, fringes, strength, and the axis against flow.

Generalises thwaites_margin_scene.py and ridge_a_scene.py into one figure
that works at any site, and adds the piece neither had - the depth-averaged
fabric axis drawn on the map and referenced to the local flow direction.

    LEFT   map view of the segment. Background is surface SPEED where the
           ice moves (ITS_LIVE covers the track and something in the season
           domain exceeds SPEED_MIN), and surface ELEVATION where it does
           not. Both cases draw the other frames of the site as context,
           the segment itself heavy, and the depth-averaged fabric axis as
           bars along it. Where speed is available the flow direction is
           drawn too and the axis-to-flow angle is annotated.
    RIGHT  two panels on a shared along-track distance axis: the wrapped
           HH-VV interferogram phase - the fringes, in the format the
           existing scene figures use - and the horizontal eigenvalue
           difference from the inversion over the same profile.

WHY THE BACKGROUND SWITCHES. At a divide the speed is small and its
DIRECTION is poorly determined, so a flow-referenced angle there would be
an arrow drawn through noise; ITS_LIVE has no coverage at Ridge A at all.
Elevation is what carries the setting at those sites. The switch is
therefore about which field is meaningful, not which file happens to be on
disk, and the figure says which one it used.

UNITS. Every angle field in a section product is stored in DEGREES
(run_quadpol_pipeline writes theta through rad2deg). Feeding those numbers
into exp(2i*theta) as if they were radians is a mistake this project has
already made and published from, so every circular statistic here converts
first and the helper below is the only place it happens.

Usage: python segment_scene.py [frame] [out_dir]
"""
import glob
import os
import sys

import h5py
import matplotlib
import numpy as np
from scipy.io import loadmat

matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.gridspec import GridSpec              # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import quadpol_sites as qs                            # noqa: E402
import scar_style                                     # noqa: E402
from scar_style import INK, MUTED                     # noqa: E402

FRAME = sys.argv[1] if len(sys.argv) > 1 else '20250108_02_009'
OUT = scar_style.out_dir(sys.argv, 2)
SCAR = os.path.expanduser('~/data/opr/scar')
ACCUM = os.path.expanduser('~/data/opr/accum')
BASE = os.path.expanduser('~/data/opr/basemap')

SPEED_MIN = 30.0          # m/yr; above this somewhere in the domain, map speed
BAND = (200.0, 1200.0)    # depth band for the depth-averaged axis
C_AIR = 299792458.0


# ---------------------------------------------------------------- helpers
def axis_mean_deg(theta_deg):
    """Circular mean of an AXIS given in DEGREES, returned in degrees.

    The doubling is what makes it an axis mean: theta and theta+180 are the
    same fabric, and averaging them naively gives their perpendicular.
    """
    a = np.deg2rad(np.asarray(theta_deg, float))
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    return np.degrees(np.angle(np.mean(np.exp(2j * a)))) / 2.0 % 180.0


def axis_to_direction_deg(axis_deg, dir_deg):
    """Angle between an AXIS and a DIRECTION, folded to [0, 90] degrees."""
    d = np.deg2rad(np.asarray(axis_deg, float) - np.asarray(dir_deg, float))
    return np.abs(np.degrees(np.angle(np.exp(2j * d))) / 2.0)


def cumdist_km(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    dla, dlo = np.diff(la), np.diff(lo)
    a = np.sin(dla / 2) ** 2 + np.cos(la[:-1]) * np.cos(la[1:]) * np.sin(dlo / 2) ** 2
    return np.r_[0.0, np.cumsum(6371.0 * 2 * np.arcsin(np.sqrt(a)))]


def read_ct(frame):
    fn = os.path.join(SCAR, 'quadpol_section_%s_ct.mat' % frame)
    if not os.path.exists(fn):
        raise SystemExit('no constant-orientation product for %s' % frame)
    with h5py.File(fn) as f:
        r = f['res']
        g = lambda k: np.array(r[k]).ravel() if k in r else None   # noqa: E731
        out = dict(z=g('z'), lat=g('sec_lat'), lon=g('sec_lon'),
                   th_geo=g('ls_theta0_geo'), zw=g('ls_zw'))
        out['dlam'] = np.array(r['sec_dlam_ls'])       # h5py gives (nz, nblk)
    return out


def read_fringes(frame):
    """Wrapped HH-VV phase, in the format the existing scenes draw."""
    seg = frame[:11]
    fn = os.path.join(ACCUM, '2024_Antarctica_Ground2',
                      'CSARP_polarimetric_unwrap', seg, 'Data_%s.mat' % frame)
    if os.path.exists(fn):
        d = loadmat(fn, variable_names=['interferogram_mlook', 'Time',
                                        'Latitude', 'Longitude', 'Surface'])
        return dict(phase=np.angle(d['interferogram_mlook']),
                    time=d['Time'].ravel(), lat=d['Latitude'].ravel(),
                    lon=d['Longitude'].ravel(), surf=d['Surface'].ravel())
    # Thwaites keeps the same chain in the margin extracts, already wrapped.
    mg = os.path.expanduser('~/data/opr/margin/margin_%s.mat' % frame)
    if os.path.exists(mg):
        d = loadmat(mg, variable_names=['phase_wrapped', 'Time', 'Latitude',
                                        'Longitude', 'Surface'])
        return dict(phase=d['phase_wrapped'], time=d['Time'].ravel(),
                    lat=d['Latitude'].ravel(), lon=d['Longitude'].ravel(),
                    surf=d['Surface'].ravel())
    return None


def read_elevation(frame):
    seg = frame[:11]
    fn = os.path.join(ACCUM, '2024_Antarctica_Ground2', 'CSARP_standard_HH',
                      seg, 'Data_%s.mat' % frame)
    if not os.path.exists(fn):
        return None, None
    with h5py.File(fn) as f:
        elev = np.array(f['Elevation']).ravel()
        surf = np.array(f['Surface']).ravel()
    return elev - surf * C_AIR / 2.0, None


def read_speed(bx, by):
    """ITS_LIVE speed and flow azimuth at projected points, or None."""
    fn = os.path.join(BASE, 'itslive_thwaites_window.nc')
    if not os.path.exists(fn):
        return None
    import xarray as xr
    d = xr.open_dataset(fn)
    x, y = d['x'].values, d['y'].values
    if not (x.min() <= np.nanmin(bx) and np.nanmax(bx) <= x.max()
            and y.min() <= np.nanmin(by) and np.nanmax(by) <= y.max()):
        return None                     # window does not cover the track
    v = d['v'].values
    vx, vy = d['vx'].values, d['vy'].values
    ix = np.clip(np.searchsorted(x, bx), 0, x.size - 1)
    iy = np.clip(np.searchsorted(y[::-1], by), 0, y.size - 1)
    iy = y.size - 1 - iy
    sp = v[iy, ix]
    az = np.degrees(np.arctan2(vx[iy, ix], vy[iy, ix])) % 360.0
    return dict(field=(x, y, v), speed=sp, az=az)


# -------------------------------------------------------------------- main
def main():
    ct = read_ct(FRAME)
    fr = read_fringes(FRAME)
    blat, blon = ct['lat'], ct['lon']
    cfg_key = None
    for k in qs.SITES:
        c = qs.get(k)
        if qs.in_site(c, np.nanmedian(blat), np.nanmedian(blon), FRAME):
            cfg_key = k
            break
    cfg = qs.get(cfg_key) if cfg_key else qs.get('ridge_a')
    T = qs.site_transformer(cfg)
    bx, by = T.transform(blon, blat)

    # Context tracks first: the speed rule is about the SEASON DOMAIN, not
    # this segment. Testing the segment alone put two frames of one site on
    # different backgrounds - 20240108_01_001 tops out at 27.6 m/yr and fell
    # to elevation while its neighbour got speed, which is the inconsistency
    # a domain rule exists to prevent.
    ctx_tracks = []
    ctx_x, ctx_y = np.array([]), np.array([])
    for fn in sorted(glob.glob(os.path.join(SCAR, 'quadpol_section_*_ct.mat'))):
        tag = os.path.basename(fn)[16:-7]
        if tag == FRAME:
            continue
        try:
            with h5py.File(fn) as f:
                la = np.array(f['res']['sec_lat']).ravel()
                lo = np.array(f['res']['sec_lon']).ravel()
        except Exception:
            continue
        if not qs.in_site(cfg, np.nanmedian(la), np.nanmedian(lo), tag):
            continue
        cx, cy = T.transform(lo, la)
        ctx_tracks.append((cx, cy))
        ctx_x = np.r_[ctx_x, cx]; ctx_y = np.r_[ctx_y, cy]

    dom_x = np.r_[bx, ctx_x]
    dom_y = np.r_[by, ctx_y]
    sp = read_speed(bx, by)
    sp_dom = read_speed(dom_x, dom_y)
    dom_max = np.nanmax(sp_dom['speed']) if sp_dom is not None else np.nan
    use_speed = sp is not None and np.isfinite(dom_max) and dom_max > SPEED_MIN

    # depth-averaged fabric axis over the trusted band, per the units note
    zw, tg = ct['zw'], ct['th_geo']
    n = min(zw.size, tg.size)
    m = (zw[:n] >= BAND[0]) & (zw[:n] <= BAND[1]) & np.isfinite(tg[:n])
    axis_deg = axis_mean_deg(tg[:n][m])

    fig = plt.figure(figsize=(15.4, 7.6), layout='constrained')
    gsp = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.62])
    axm = fig.add_subplot(gsp[:, 0])
    axf = fig.add_subplot(gsp[0, 1])
    axd = fig.add_subplot(gsp[1, 1], sharex=axf)

    # ---- map
    if use_speed:
        x, y, v = sp['field']
        axm.pcolormesh(x / 1e3, y / 1e3, np.log10(np.maximum(v, 1.0)),
                       cmap='viridis', shading='auto', rasterized=True)
        axm.set_title('surface speed (ITS_LIVE), log$_{10}$ m yr$^{-1}$   '
                      '(domain max %.0f m yr$^{-1}$)' % dom_max,
                      fontsize=10, color=INK)
        flow_az = np.nanmedian(sp['az'])
    else:
        elev, _ = read_elevation(FRAME)
        if elev is not None and np.isfinite(elev).any():
            e = np.interp(np.linspace(0, 1, blat.size),
                          np.linspace(0, 1, elev.size), elev)
            s = axm.scatter(bx / 1e3, by / 1e3, c=e, s=46, cmap='cividis',
                            zorder=5, edgecolors='none')
            cb = fig.colorbar(s, ax=axm, orientation='horizontal',
                              shrink=0.7, pad=0.08, aspect=34)
            cb.set_label('surface elevation (m)  -  a few metres over the '
                         'whole segment, which is what a divide looks like',
                         fontsize=8.6)
        why = ('domain max %.0f m yr$^{-1}$' % dom_max
               if np.isfinite(dom_max) else 'no ITS_LIVE coverage')
        axm.set_title('surface elevation  -  slow ice, no usable flow '
                      'direction (%s)' % why, fontsize=10, color=INK)
        flow_az = np.nan

    for cx, cy in ctx_tracks:
        axm.plot(cx / 1e3, cy / 1e3, '-', color='0.72', lw=1.0, zorder=2)
    axm.plot(bx / 1e3, by / 1e3, '-', color=INK, lw=2.6, zorder=4)

    # depth-averaged axis as bars along the segment
    if np.isfinite(axis_deg):
        span = max(np.ptp(bx), np.ptp(by)) / 1e3
        half = 0.13 * span
        gnaz = qs.grid_north_az(blat, blon, cfg=cfg)
        for i in np.linspace(0, blat.size - 1, 6).astype(int):
            a = np.deg2rad(axis_deg - gnaz[i])       # true -> grid azimuth
            dx, dy = half * np.sin(a), half * np.cos(a)
            for lw_, col, zo in ((5.2, 'white', 6), (2.6, '#c0392b', 7)):
                axm.plot([bx[i] / 1e3 - dx, bx[i] / 1e3 + dx],
                         [by[i] / 1e3 - dy, by[i] / 1e3 + dy],
                         '-', color=col, lw=lw_, zorder=zo,
                         solid_capstyle='round')
    if np.isfinite(flow_az):
        a = np.deg2rad(flow_az - np.nanmedian(qs.grid_north_az(blat, blon, cfg=cfg)))
        span = max(np.ptp(bx), np.ptp(by)) / 1e3
        x0, y0 = np.nanmedian(bx) / 1e3, np.nanmedian(by) / 1e3
        axm.annotate('', xy=(x0 + 0.16 * span * np.sin(a),
                             y0 + 0.16 * span * np.cos(a)), xytext=(x0, y0),
                     arrowprops=dict(arrowstyle='-|>', color='#2a78d6', lw=2.0),
                     zorder=7)
        ang = axis_to_direction_deg(axis_deg, flow_az)
        note = ('fabric axis %.0f$^\\circ$, flow %.0f$^\\circ$\n'
                'axis to flow: %.0f$^\\circ$' % (axis_deg, flow_az, ang))
    else:
        note = ('fabric axis %.0f$^\\circ$ (true)\n'
                'no flow reference: ice too slow here' % axis_deg)
    axm.text(0.03, 0.03, note, transform=axm.transAxes, fontsize=9.4,
             color=INK, va='bottom', linespacing=1.6,
             bbox=dict(fc='white', ec='0.8', alpha=0.9, pad=3.0))
    # Zoom to the segment and its context tracks. Drawing the full speed
    # raster extent put a 5 km segment on a 280 km map.
    ex = np.r_[bx, ctx_x] / 1e3
    ey = np.r_[by, ctx_y] / 1e3
    if ex.size and np.isfinite(ex).any():
        cx, cy = np.nanmean(ex), np.nanmean(ey)
        half = 0.62 * max(np.nanmax(ex) - np.nanmin(ex),
                          np.nanmax(ey) - np.nanmin(ey), 4.0)
        axm.set_xlim(cx - half, cx + half)
        axm.set_ylim(cy - half, cy + half)
    axm.set_aspect('equal')
    axm.tick_params(labelsize=9)
    axm.set_xlabel('EPSG:%d easting (km)   -   bars drawn in grid azimuth'
                   % qs.site_epsg(cfg))
    axm.set_ylabel('northing (km)')

    # ---- fringes
    dist_b = cumdist_km(blat, blon)
    if fr is not None:
        dist_h = cumdist_km(fr['lat'], fr['lon'])
        depth = (fr['time'][:, None] - fr['surf'][None, :]) * C_AIR / 2 / 1.78
        dm = np.nanmedian(depth, axis=1)
        keep = (dm > 0) & (dm < 1600)
        axf.pcolormesh(dist_h, dm[keep], fr['phase'][keep, :], cmap='hsv',
                       vmin=-np.pi, vmax=np.pi, shading='auto', rasterized=True)
        axf.set_ylim(1600, 0)
        axf.set_ylabel('depth (m)')
        axf.set_title('wrapped HH$-$VV interferogram phase (fringes)',
                      fontsize=10, color=INK)
    else:
        axf.text(0.5, 0.5, 'no polarimetric product mirrored for this frame',
                 transform=axf.transAxes, ha='center', color=MUTED)
        axf.set_xlim(dist_b.min(), dist_b.max())

    # ---- strength
    dl = ct['dlam']
    if dl.shape[0] != ct['z'].size:
        dl = dl.T
    pm = axd.pcolormesh(dist_b, ct['z'], dl, cmap='Blues',
                        vmin=0.0, vmax=cfg.get('vmax', 0.08),
                        shading='auto', rasterized=True)
    fig.colorbar(pm, ax=axd, shrink=0.85, pad=0.02, label=r'$\Delta\lambda$')
    axd.set_ylim(1600, 0)
    axd.set_ylabel('depth (m)')
    axd.set_xlabel('along-segment distance (km)')
    axd.set_title(r'fabric strength $\Delta\lambda$, axis held through the column',
                  fontsize=10, color=INK)

    fig.suptitle('%s  -  %s' % (FRAME, cfg.get('title', cfg_key or '')),
                 fontsize=13, color=INK)
    fn = os.path.join(OUT, 'segment_scene_%s.png' % FRAME)
    fig.savefig(fn, dpi=180, facecolor='white')
    print('wrote %s' % fn)


if __name__ == '__main__':
    main()

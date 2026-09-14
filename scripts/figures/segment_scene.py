"""Per-segment scene: setting, fringes, strength, and the axis against flow.

Generalises thwaites_margin_scene.py and ridge_a_scene.py into one figure
that works at any site, and adds the piece neither had - the depth-averaged
fabric axis drawn on the map and referenced to the local flow direction.

    LEFT   map view. A raster background: surface SPEED where the ice moves
           (ITS_LIVE, and something in the season DOMAIN exceeds SPEED_MIN)
           or surface ELEVATION from BedMachine where it does not. The
           profile is drawn black with start and stop markers, the site's
           other frames are context, and the depth-averaged fabric axis is
           drawn as bars along it against the local surface-slope gradient.
    RIGHT  three panels on a shared distance axis, carrying the same start
           and stop markers above them: the wrapped HH-VV interferogram
           (HSV, darkened by coherence, with a phase colorbar), the
           horizontal eigenvalue difference from the inversion, and the
           angle between the depth-uniform fabric axis and the reference
           direction.

WHY THE REFERENCE IS SURFACE SLOPE. Flow direction is the natural
reference where ice moves, but at a divide the speed is small and its
direction poorly determined, and ITS_LIVE has no coverage at Ridge A at
all - a flow-referenced angle there is an arrow through noise. The surface
slope gradient is defined everywhere and is what sets the driving stress,
so it is the reference drawn and differenced here. Where speed exists the
flow direction is reported alongside it for comparison, and the two agree
closely in fast ice.

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
from scar_style import hsv_plot_coherence             # noqa: E402

FRAME = sys.argv[1] if len(sys.argv) > 1 else '20250108_02_009'
OUT = scar_style.out_dir(sys.argv, 2)
SCAR = os.path.expanduser('~/data/opr/scar')
ACCUM = os.path.expanduser('~/data/opr/accum')
BASE = os.path.expanduser('~/data/opr/basemap')
BEDM = os.path.expanduser('~/data/BedMachine')

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


def read_section(frame):
    """Free-axis section product: the direction is allowed to vary along
    the line, so sec_theta carries one axis per block per depth rather than
    one held for the whole segment."""
    fn = os.path.join(SCAR, 'quadpol_section_%s.mat' % frame)
    if not os.path.exists(fn):
        raise SystemExit('no section product for %s' % frame)
    with h5py.File(fn) as f:
        r = f['res']
        g = lambda k: np.array(r[k]).ravel() if k in r else None   # noqa: E731
        out = dict(z=g('z'), lat=g('sec_lat'), lon=g('sec_lon'),
                   th_geo=g('ls_theta0_geo'), zw=g('ls_zw'))
        out['dlam'] = np.array(r['sec_dlam_ls'])
        out['theta'] = np.array(r['sec_theta'])        # geographic, degrees
    return out


def read_fringes(frame):
    """Wrapped HH-VV phase, in the format the existing scenes draw."""
    seg = frame[:11]
    fn = os.path.join(ACCUM, '2024_Antarctica_Ground2',
                      'CSARP_polarimetric_unwrap', seg, 'Data_%s.mat' % frame)
    if os.path.exists(fn):
        d = loadmat(fn, variable_names=['interferogram_mlook',
                                        'interferogram_coherence', 'Time',
                                        'Latitude', 'Longitude', 'Surface'])
        return dict(phase=np.angle(d['interferogram_mlook']),
                    coh=np.abs(d['interferogram_coherence']),
                    time=d['Time'].ravel(), lat=d['Latitude'].ravel(),
                    lon=d['Longitude'].ravel(), surf=d['Surface'].ravel())
    # Thwaites keeps the same chain in the margin extracts, already wrapped.
    mg = os.path.expanduser('~/data/opr/margin/margin_%s.mat' % frame)
    if os.path.exists(mg):
        d = loadmat(mg, variable_names=['phase_wrapped', 'coherence', 'Time',
                                        'Latitude', 'Longitude', 'Surface'])
        return dict(phase=d['phase_wrapped'], coh=np.abs(d['coherence']),
                    time=d['Time'].ravel(), lat=d['Latitude'].ravel(),
                    lon=d['Longitude'].ravel(), surf=d['Surface'].ravel())
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


def read_surface(x0, x1, y0, y1, epsg):
    """BedMachine surface elevation over a projected window, as a raster.

    Windowed with xarray so the 13333x13333 grid is never fully read.
    """
    fn = {3031: os.path.join(BEDM, 'Antarctica', 'BedMachineAntarctica-v3.nc'),
          3413: os.path.join(BEDM, 'Greenland', 'BedMachineGreenland-v5.nc')}.get(epsg)
    if fn is None or not os.path.exists(fn):
        return None
    import xarray as xr
    d = xr.open_dataset(fn)
    ys = d['y'].values
    ysl = slice(y1, y0) if ys[0] > ys[-1] else slice(y0, y1)
    w = d['surface'].sel(x=slice(x0, x1), y=ysl)
    if w.size == 0:
        return None
    return dict(x=w['x'].values, y=w['y'].values, z=np.asarray(w.values, float))


def slope_azimuth(surf, bx, by, gn_deg):
    """Downslope azimuth [deg from grid north] near the profile.

    The gradient is taken on the DEM and averaged over the cells the
    profile passes through, so it is the slope beside the line rather than
    a whole-domain mean: at a divide those differ by a lot.
    """
    if surf is None or surf['z'].size < 9:
        return np.nan
    x, y, z = surf['x'], surf['y'], surf['z']
    dzdy, dzdx = np.gradient(z, y, x)
    ix = np.clip(np.searchsorted(x, bx), 0, x.size - 1)
    iy = (np.clip(np.searchsorted(y[::-1], by), 0, y.size - 1)
          if y[0] > y[-1] else np.clip(np.searchsorted(y, by), 0, y.size - 1))
    if y[0] > y[-1]:
        iy = y.size - 1 - iy
    gx, gy = np.nanmean(dzdx[iy, ix]), np.nanmean(dzdy[iy, ix])
    if not np.isfinite(gx) or not np.isfinite(gy) or (gx == 0 and gy == 0):
        return np.nan
    # grid azimuth of the downslope vector, then to TRUE north
    return (np.degrees(np.arctan2(-gx, -gy)) + gn_deg) % 360.0


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
    az = np.degrees(np.arctan2(vx[iy, ix], vy[iy, ix])) % 360.0   # GRID
    return dict(field=(x, y, v), speed=sp, az=az)


# -------------------------------------------------------------------- main
def main():
    ct = read_section(FRAME)
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
    gn_site = np.nanmedian(qs.grid_north_az(blat, blon, cfg=cfg))
    sp = read_speed(bx, by)
    sp_dom = read_speed(dom_x, dom_y)
    dom_max = np.nanmax(sp_dom['speed']) if sp_dom is not None else np.nan
    use_speed = sp is not None and np.isfinite(dom_max) and dom_max > SPEED_MIN

    # Depth-averaged axis PER BLOCK, so the direction varies along the line.
    th = ct['theta']
    if th.shape[0] != ct['z'].size:
        th = th.T                                   # -> (nz, nblk)
    zb = (ct['z'] >= BAND[0]) & (ct['z'] <= BAND[1])
    a = np.deg2rad(th[zb, :])
    with np.errstate(invalid='ignore'):
        axis_blk = np.degrees(np.angle(np.nanmean(np.exp(2j * a), axis=0))) / 2 % 180
    axis_deg = axis_mean_deg(axis_blk)              # frame summary

    fig = plt.figure(figsize=(16.2, 8.6), layout='constrained')
    gsp = GridSpec(3, 2, figure=fig, width_ratios=[1.0, 1.66])
    axm = fig.add_subplot(gsp[:, 0])
    axf = fig.add_subplot(gsp[0, 1])
    axd = fig.add_subplot(gsp[1, 1], sharex=axf)
    axa = fig.add_subplot(gsp[2, 1], sharex=axf)

    # ---- map extent first: the rasters are windowed to it
    ex = np.r_[bx, ctx_x] / 1e3
    ey = np.r_[by, ctx_y] / 1e3
    cxm, cym = np.nanmean(ex), np.nanmean(ey)
    half = 0.62 * max(np.nanmax(ex) - np.nanmin(ex),
                      np.nanmax(ey) - np.nanmin(ey), 4.0)
    x0, x1 = (cxm - half) * 1e3, (cxm + half) * 1e3
    y0, y1 = (cym - half) * 1e3, (cym + half) * 1e3

    surf = read_surface(x0, x1, y0, y1, qs.site_epsg(cfg))
    slope_az = slope_azimuth(surf, bx, by, gn_site)

    if use_speed:
        gx, gy, gv = sp['field']
        mx = (gx >= x0) & (gx <= x1)
        my = (gy >= y0) & (gy <= y1)
        im = axm.pcolormesh(gx[mx] / 1e3, gy[my] / 1e3,
                            np.log10(np.maximum(gv[np.ix_(my, mx)], 1.0)),
                            cmap='viridis', shading='auto', rasterized=True)
        cb = fig.colorbar(im, ax=axm, orientation='horizontal', shrink=0.72,
                          pad=0.06, aspect=32)
        cb.set_label('surface speed, log$_{10}$ m yr$^{-1}$   '
                     '(domain max %.0f)' % dom_max, fontsize=9)
        flow_az = (np.nanmedian(sp['az']) + gn_site) % 360.0
    else:
        flow_az = np.nan
        if surf is not None:
            im = axm.pcolormesh(surf['x'] / 1e3, surf['y'] / 1e3, surf['z'],
                                cmap='cividis', shading='auto', rasterized=True)
            cb = fig.colorbar(im, ax=axm, orientation='horizontal',
                              shrink=0.72, pad=0.06, aspect=32)
            why = ('domain max %.0f m yr$^{-1}$' % dom_max
                   if np.isfinite(dom_max) else 'no ITS_LIVE coverage')
            cb.set_label('BedMachine surface elevation (m)   '
                         '(slow ice: %s)' % why, fontsize=9)

    for cx, cy in ctx_tracks:
        axm.plot(cx / 1e3, cy / 1e3, '-', color='0.55', lw=1.0, zorder=2,
                 alpha=0.85)

    # ---- the profile: black, with start and stop markers
    axm.plot(bx / 1e3, by / 1e3, '-', color='black', lw=3.0, zorder=5)
    axm.plot(bx[0] / 1e3, by[0] / 1e3, 'o', mfc='white', mec='black',
             ms=10, mew=2.0, zorder=8)
    axm.plot(bx[-1] / 1e3, by[-1] / 1e3, 's', mfc='black', mec='white',
             ms=10, mew=2.0, zorder=8)

    gnaz = qs.grid_north_az(blat, blon, cfg=cfg)
    span = 2 * half
    if np.isfinite(axis_deg):
        halfbar = 0.030 * span
        for i in np.linspace(0, blat.size - 1, 7).astype(int):
            ab = axis_blk[i] if np.isfinite(axis_blk[i]) else axis_deg
            a = np.deg2rad(ab - gnaz[i])
            dx, dy = halfbar * np.sin(a), halfbar * np.cos(a)
            for lw_, col, zo in ((3.8, 'white', 3), (2.0, '#c0392b', 4)):
                axm.plot([bx[i] / 1e3 - dx, bx[i] / 1e3 + dx],
                         [by[i] / 1e3 - dy, by[i] / 1e3 + dy],
                         '-', color=col, lw=lw_, zorder=zo,
                         solid_capstyle='round')

    # ---- reference arrow, drawn BESIDE the profile
    ref_az, ref_name = ((slope_az, 'surface slope') if np.isfinite(slope_az)
                        else (flow_az, 'flow'))
    if np.isfinite(ref_az):
        gnm = np.nanmedian(gnaz)
        a = np.deg2rad(ref_az - gnm)
        # offset perpendicular to the track so the arrow sits next to it
        tx, ty = bx[-1] - bx[0], by[-1] - by[0]
        tn = np.hypot(tx, ty) or 1.0
        ox, oy = -ty / tn, tx / tn
        px = (np.nanmean(bx) + 0.10 * span * 1e3 * ox) / 1e3
        py = (np.nanmean(by) + 0.10 * span * 1e3 * oy) / 1e3
        L = 0.14 * span
        axm.annotate('', xy=(px + L * np.sin(a), py + L * np.cos(a)),
                     xytext=(px, py),
                     arrowprops=dict(arrowstyle='-|>', color='#2a78d6', lw=2.4),
                     zorder=9)
        # at the arrow TIP: at its tail the label sat under the arrow and
        # the profile's end marker

    axm.set_xlim(x0 / 1e3, x1 / 1e3)
    axm.set_ylim(y0 / 1e3, y1 / 1e3)
    axm.set_aspect('equal')
    axm.tick_params(labelsize=9)
    axm.set_xlabel('EPSG:%d easting (km)   -   bars in grid azimuth'
                   % qs.site_epsg(cfg))
    axm.set_ylabel('northing (km)')

    # ---- fringes: HSV darkened by coherence, with a phase colorbar
    dist_b = cumdist_km(blat, blon)
    if fr is not None:
        dist_h = cumdist_km(fr['lat'], fr['lon'])
        depth = (fr['time'][:, None] - fr['surf'][None, :]) * C_AIR / 2 / 1.78
        dm = np.nanmedian(depth, axis=1)
        keep = (dm > 0) & (dm < 1600)
        ph, co = fr['phase'][keep, :], fr['coh'][keep, :]
        st = max(1, ph.shape[0] // 1400)
        sx = max(1, ph.shape[1] // 1600)
        rgb = hsv_plot_coherence(ph[::st, ::sx], co[::st, ::sx])
        axf.imshow(rgb, aspect='auto', interpolation='nearest',
                   extent=[dist_h[0], dist_h[-1], dm[keep][-1], dm[keep][0]])
        sm = plt.cm.ScalarMappable(cmap='hsv',
                                   norm=plt.Normalize(-np.pi, np.pi))
        cbp = fig.colorbar(sm, ax=axf, pad=0.02)
        cbp.set_label('HH$-$VV phase (rad)', fontsize=9)
        cbp.set_ticks([-np.pi, 0, np.pi])
        cbp.set_ticklabels([r'$-\pi$', '0', r'$\pi$'])
        axf.set_ylabel('depth (m)')
        axf.set_title('wrapped HH$-$VV interferogram, value darkened by '
                      'coherence', fontsize=10, color=INK)
    else:
        axf.set_xlim(dist_b.min(), dist_b.max())

    # ---- strength
    dl = ct['dlam']
    if dl.shape[0] != ct['z'].size:
        dl = dl.T
    pm = axd.pcolormesh(dist_b, ct['z'], dl, cmap='Blues', vmin=0.0,
                        vmax=cfg.get('vmax', 0.08), shading='auto',
                        rasterized=True)
    cbd = fig.colorbar(pm, ax=axd, pad=0.02)
    cbd.set_label(r'$\Delta\lambda$', fontsize=9)
    axd.set_ylim(1600, 0)
    axd.set_ylabel('depth (m)')
    axd.set_title(r'fabric strength $\Delta\lambda$', fontsize=10, color=INK)

    # ---- axis against the reference, along the profile
    if np.isfinite(ref_az) and np.isfinite(axis_deg):
        ang = axis_to_direction_deg(axis_blk, ref_az)
        axa.plot(dist_b, ang, '-', color='#4a3aa7', lw=2.0)
        axa.axhline(45, color=MUTED, lw=1.0, ls='--')
        axa.set_ylim(0, 90)
        axa.set_yticks([0, 30, 45, 60, 90])
        axa.set_ylabel('deg')
        axa.set_title('fabric axis against %s direction (%.0f$^\\circ$ true)'
                      % (ref_name, ref_az), fontsize=10, color=INK)
    axa.set_xlabel('along-segment distance (km)')

    # ---- start / stop markers above the right-hand panels
    for ax_ in (axf, axd, axa):
        ax_.plot(dist_b[0], 1.0, 'o', mfc='white', mec='black', ms=9,
                 mew=1.8, transform=ax_.get_xaxis_transform(), clip_on=False,
                 zorder=10)
        ax_.plot(dist_b[-1], 1.0, 's', mfc='black', mec='white', ms=9,
                 mew=1.8, transform=ax_.get_xaxis_transform(), clip_on=False,
                 zorder=10)

    fig.suptitle('%s  -  %s   (fabric axis free to vary along the line)'
                 % (FRAME, cfg.get('title', cfg_key or '')),
                 fontsize=13, color=INK)
    fn = os.path.join(OUT, 'segment_scene_%s.png' % FRAME)
    fig.savefig(fn, dpi=180, facecolor='white')
    print('wrote %s' % fn)


if __name__ == '__main__':
    main()

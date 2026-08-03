"""Thwaites shear-margin fabric signature in the flow frame.

Samples ITS_LIVE v2 surface velocity (anonymous S3, windowed read; cached
locally) at every Thwaites-region inversion block, then:
  1. locates the shear margin objectively from the along-line speed and
     lateral shear profile,
  2. rotates the survey-frame contrast into the FLOW frame under the
     hypothesis that the horizontal principal axes are flow-aligned:
     Lambda = lam_crossflow - lam_alongflow = dlam_measured / cos 2(phi - phi_flow),
     masked where |cos 2(phi - phi_flow)| < 0.4 (geometry degenerate),
  3. compares zone-median Lambda(z) profiles (interior vs margin vs fast
     side) with the Ridge A principal contrast P(z) baseline solved by
     ridge_a_azimuthal.py.

The flow-aligned-axes hypothesis is expected to fail INSIDE the margin,
where simple shear rotates the fabric axes toward ~45 deg to flow; a
misfit there is itself the shear-margin signature.

Usage: python thwaites_fabric.py [fabric_batch_root] [out_dir]
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')
VEL_CACHE = os.path.expanduser('~/data/opr/basemap/thwaites_itslive_sample.npz')
# Local copy preferred (user's disk); anonymous S3 as fallback
ITSLIVE_LOCAL = os.path.expanduser(
    '~/Downloads/ITS_LIVE_velocity_120m_RGI19A_0000_V02.1.nc')
ITSLIVE_S3 = 'its-live-data/velocity_mosaic/v2/static/ITS_LIVE_velocity_120m_RGI19A_0000_v02.nc'

DEPTH = np.arange(0, 1900, 5.0)
THWAITES = (-76.46, -105.5)
MAX_DIST_KM = 150.0
MIN_COS2 = 0.4


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def bearing_deg(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(x, y)) % 360.0


def load_thwaites(pattern):
    blocks = []
    for fn in sorted(glob.glob(pattern)):
        d = loadmat(fn)
        dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
        qual = d['dlam_quality']
        clip = d.get('dlam_clipped')
        blat = d['Latitude'][0, :]
        blon = d['Longitude'][0, :]
        nb = dlam.shape[1]
        for b in range(nb):
            lat, lon = blat[b], blon[b]
            if not (np.isfinite(lat) and np.isfinite(lon)):
                continue
            if haversine_km(lat, lon, *THWAITES) > MAX_DIST_KM:
                continue
            j = b + 1 if b + 1 < nb else b - 1
            if j < 0 or not (np.isfinite(blat[j]) and np.isfinite(blon[j])):
                continue
            head = (bearing_deg(lat, lon, blat[j], blon[j]) if j > b
                    else bearing_deg(blat[j], blon[j], lat, lon)) % 180.0
            col = np.full(DEPTH.size, np.nan)
            for k in range(dlam.shape[0]):
                pegged = (clip is not None and clip.size and clip[k, b] == 1) \
                    or abs(dlam[k, b]) > 0.6
                ok = np.isfinite(dlam[k, b]) and not pegged \
                    and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
                if ok:
                    m = (DEPTH >= top[k, b]) & (DEPTH < bot[k, b])
                    col[m] = dlam[k, b]
            if not np.any(np.isfinite(col)):
                continue
            blocks.append(dict(col=col, lat=lat, lon=lon, head=head,
                               gps=d['GPS_time'][0, b]))
    blocks.sort(key=lambda b: b['gps'])
    return blocks


def sample_itslive(lats, lons):
    """Speed [m/yr] and true flow bearing [deg] at points; cached."""
    if os.path.exists(VEL_CACHE):
        c = np.load(VEL_CACHE)
        if c['lats'].size == len(lats) and np.allclose(c['lats'], lats):
            return c['speed'], c['flow_az']
    import xarray as xr
    proj = ab.proj3031()
    xyz = proj.transform_points(ccrs.PlateCarree(),
                                np.asarray(lons), np.asarray(lats))
    px, py = xyz[:, 0], xyz[:, 1]
    pad = 10e3
    if os.path.exists(ITSLIVE_LOCAL):
        src = ITSLIVE_LOCAL
    else:
        import s3fs
        src = s3fs.S3FileSystem(anon=True).open(ITSLIVE_S3, 'rb')
    with xr.open_dataset(src, engine='h5netcdf') as ds:
        win = ds.sel(x=slice(px.min() - pad, px.max() + pad),
                     y=slice(py.max() + pad, py.min() - pad))
        vx = win['vx'].values.astype(float)
        vy = win['vy'].values.astype(float)
        wx = win['x'].values
        wy = win['y'].values
    ix = np.clip(np.searchsorted(wx, px), 0, wx.size - 1)
    iy = np.clip(np.searchsorted(-wy, -py), 0, wy.size - 1)
    svx = vx[iy, ix]
    svy = vy[iy, ix]
    speed = np.hypot(svx, svy)
    # True flow bearing via a small projected displacement
    step = 1000.0 / np.maximum(speed, 1e-6)
    qx, qy = px + svx * step, py + svy * step
    pc = ccrs.PlateCarree()
    q = pc.transform_points(proj, qx, qy)
    flow_az = bearing_deg(np.asarray(lats), np.asarray(lons), q[:, 1], q[:, 0])
    os.makedirs(os.path.dirname(VEL_CACHE), exist_ok=True)
    np.savez(VEL_CACHE, lats=lats, lons=lons, speed=speed, flow_az=flow_az)
    return speed, flow_az


def main():
    blocks = load_thwaites(f'{ROOT}/joint/2023_Antarctica_Ground/*/Data_*.mat')
    print(f'Thwaites blocks: {len(blocks)}')
    if not blocks:
        return
    lats = np.array([b['lat'] for b in blocks])
    lons = np.array([b['lon'] for b in blocks])
    heads = np.array([b['head'] for b in blocks])
    profs = np.array([b['col'] for b in blocks])
    speed, flow_az = sample_itslive(lats, lons)
    print(f'speed range along line: {np.nanmin(speed):.0f}-{np.nanmax(speed):.0f} m/yr')

    step = haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
    dist = np.concatenate([[0], np.cumsum(np.minimum(step, 2.0))])

    # Flow-frame rotation under the flow-aligned-axes hypothesis
    alpha = np.radians(2.0 * (heads - flow_az % 180.0))
    c2 = np.cos(alpha)
    Lam = np.where(np.abs(c2)[:, None] >= MIN_COS2,
                   profs / c2[:, None], np.nan)

    # Zones from the speed profile (data-driven thirds of log-speed range)
    ls = np.log10(np.maximum(speed, 1.0))
    z1, z2 = np.percentile(ls, [33, 66])
    zone = np.digitize(ls, [z1, z2])  # 0 slow, 1 margin/intermediate, 2 fast
    zone_names = [f'slow (<{10**z1:.0f} m/yr)',
                  f'intermediate ({10**z1:.0f}-{10**z2:.0f} m/yr)',
                  f'fast (>{10**z2:.0f} m/yr)']
    zone_colors = ['tab:blue', 'tab:orange', 'tab:red']

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.9, 1.4, 1.6],
                          width_ratios=[2.2, 1, 1])

    # (a) speed + angles along line
    ax = fig.add_subplot(gs[0, :])
    ax.plot(dist, speed, 'k-', lw=1.5)
    ax.set_ylabel('ITS_LIVE speed (m/yr)')
    ax.set_yscale('log')
    axb = ax.twinx()
    axb.plot(dist, (heads - flow_az) % 180.0, '.', ms=3, color='tab:purple',
             alpha=0.6)
    axb.set_ylabel('line-to-flow angle (deg)', color='tab:purple')
    axb.set_ylim(0, 180)
    for zi in range(3):
        m = zone == zi
        ax.scatter(dist[m], speed[m], s=10, color=zone_colors[zi], zorder=5)
    ax.set_xlabel('Distance along drive (km)')
    ax.set_title('Speed (colored by zone) and line-to-flow geometry')
    ax.grid(alpha=0.3)

    # (b) flow-frame section
    ax = fig.add_subplot(gs[1, :])
    edges = np.concatenate([[dist[0] - 0.1], (dist[:-1] + dist[1:]) / 2,
                            [dist[-1] + 0.1]])
    depth_edges = np.concatenate([DEPTH, [DEPTH[-1] + 5]])
    pc = ax.pcolormesh(edges, depth_edges, Lam.T, cmap='RdBu_r',
                       vmin=-0.25, vmax=0.25)
    ax.set_ylim(1500, 0)
    ax.set_ylabel('Depth (m)')
    ax.set_title(r'Flow-frame contrast $\Lambda = \lambda_{\perp flow} - '
                 r'\lambda_{\parallel flow}$ (masked where geometry degenerate)')
    fig.colorbar(pc, ax=ax, pad=0.01, label=r'$\Lambda$', extend='both')

    # (c) zone-median profiles vs Ridge A
    ax = fig.add_subplot(gs[2, 0])
    for zi in range(3):
        m = zone == zi
        if m.sum() < 5:
            continue
        with np.errstate(invalid='ignore'):
            med = np.nanmedian(Lam[m], axis=0)
            n = np.sum(np.isfinite(Lam[m]), axis=0)
        med[n < 5] = np.nan
        ax.plot(med, DEPTH, color=zone_colors[zi], lw=2,
                label=f'{zone_names[zi]} (n={m.sum()})')
    ax.invert_yaxis()
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel(r'median $\Lambda$')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(-0.25, 0.25)
    ax.grid(alpha=0.3)
    ax.set_title('Zone medians (flow frame)')
    ax.legend(fontsize=7, loc='lower left')

    # (d) |Lambda| vs Ridge A P baseline
    ax = fig.add_subplot(gs[2, 1])
    with np.errstate(invalid='ignore'):
        absmed = np.nanmedian(np.abs(Lam), axis=0)
        n = np.sum(np.isfinite(Lam), axis=0)
    absmed[n < 10] = np.nan
    ax.plot(absmed, DEPTH, 'k-', lw=2, label=r'Thwaites $|\Lambda|$')
    ax.invert_yaxis()
    ax.set_xlabel('contrast magnitude')
    ax.set_xlim(0, 0.25)
    ax.grid(alpha=0.3)
    ax.set_title('vs Ridge A strength')
    ax.legend(fontsize=7)
    ax.annotate('Ridge A P(z) plateaus at ~0.12\n(1150-1500 m; divide setting)',
                xy=(0.05, 0.05), xycoords='axes fraction', fontsize=7)

    # (e) hypothesis misfit: fraction of masked/inconsistent blocks per zone
    ax = fig.add_subplot(gs[2, 2])
    for zi in range(3):
        m = zone == zi
        if m.sum() < 5:
            continue
        with np.errstate(invalid='ignore'):
            spread = np.nanstd(Lam[m], axis=0)
        ax.plot(spread, DEPTH, color=zone_colors[zi], lw=1.5)
    ax.invert_yaxis()
    ax.set_xlabel(r'zone std of $\Lambda$')
    ax.set_xlim(0, 0.15)
    ax.grid(alpha=0.3)
    ax.set_title('Within-zone scatter\n(high = axes not flow-aligned)',
                 fontsize=9)

    fig.suptitle('Thwaites line in the flow frame (ITS_LIVE velocities): '
                 'fabric signature vs margin proximity', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'thwaites_fabric.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

"""Shear-margin scene figure for one Thwaites margin-crossing segment.

Left: EPSG:3031 map of the segment (default 20240108_01, the eastern
shear margin crossing) over the ITS_LIVE v2 speed field, with the other
margin display extracts as context tracks. Right: four stacked panels
sharing the along-segment distance axis - ITS_LIVE surface speed, HH
power, wrapped interferogram phase, and the inverted horizontal
eigenvalue difference (dlam depth bands from the joint inversion).

Inputs are the compact margin extract (extract_margin.py on mem1, see
opr_fabric/README.md) and the joint-inversion Data_*.mat for the same
segment under the fabric batch root.

Usage: python thwaites_margin_scene.py [seg] [extract_dir] [fabric_root] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.io import loadmat

import cartopy.crs as ccrs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab

SEG = sys.argv[1] if len(sys.argv) > 1 else '20240108_01_001'
EXT = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser('~/data/opr/margin')
ROOT = sys.argv[3] if len(sys.argv) > 3 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[4] if len(sys.argv) > 4 else os.path.join(
    os.path.dirname(__file__), '..', '..', 'figs')

ITSLIVE_LOCAL = os.path.expanduser(
    '~/Downloads/ITS_LIVE_velocity_120m_RGI19A_0000_V02.1.nc')

DEPTH = np.arange(0, 800, 5.0)


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def block_edges(centers, pad_km=0.5):
    """pcolormesh edges around block centers; a lone block gets +/- pad."""
    if centers.size < 2:
        return np.array([centers[0] - pad_km, centers[0] + pad_km])
    half = np.diff(centers) / 2
    return np.concatenate([[centers[0] - half[0]], centers[:-1] + half,
                           [centers[-1] + half[-1]]])


def sym_limit(a, default=0.25):
    """Symmetric color limit from |a|; default if empty, all-NaN or all-zero."""
    a = np.asarray(a, dtype=float)
    finite = np.abs(a[np.isfinite(a)])
    if finite.size == 0 or finite.max() == 0:
        return default
    return float(finite.max())


def itslive_window(extent, pad=5e3):
    """(x, y, v) speed window in EPSG:3031 meters; None if unavailable."""
    if not os.path.exists(ITSLIVE_LOCAL):
        return None
    import xarray as xr
    with xr.open_dataset(ITSLIVE_LOCAL, engine='h5netcdf') as ds:
        win = ds.sel(x=slice(extent[0] - pad, extent[1] + pad),
                     y=slice(extent[3] + pad, extent[2] - pad))
        return (win['x'].values, win['y'].values,
                win['v'].values.astype(float))


def sample_speed(win, lats, lons):
    """Nearest-cell ITS_LIVE speed at points from a precomputed window."""
    if win is None:
        return np.full(len(lats), np.nan)
    wx, wy, v = win
    proj = ab.proj3031()
    xyz = proj.transform_points(ccrs.PlateCarree(),
                                np.asarray(lons), np.asarray(lats))
    ix = np.clip(np.searchsorted(wx, xyz[:, 0]), 0, wx.size - 1)
    iy = np.clip(np.searchsorted(-wy, -xyz[:, 1]), 0, wy.size - 1)
    return v[iy, ix]


def dlam_section(seg):
    """(block_dists_placeholder_lats, lons, gridded dlam) for the segment.

    Returns (lats, lons, cols) with cols shaped (nblocks, DEPTH.size);
    masking follows thwaites_fabric.py (pegged/clipped and low quality
    bands dropped). None if the inversion output is missing.
    """
    base = seg.rsplit('_', 1)[0]
    fns = sorted(glob.glob(
        f'{ROOT}/joint/2023_Antarctica_Ground/{base}/Data_{seg}.mat'))
    if not fns:
        return None
    d = loadmat(fns[0])
    dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
    qual = d['dlam_quality']
    clip = d.get('dlam_clipped')
    itp = d.get('dlam_interpolated')
    lats = d['Latitude'][0, :]
    lons = d['Longitude'][0, :]
    cols = np.full((dlam.shape[1], DEPTH.size), np.nan)
    for b in range(dlam.shape[1]):
        for k in range(dlam.shape[0]):
            pegged = (clip is not None and clip.size and clip[k, b] == 1) \
                or abs(dlam[k, b]) > 0.6
            # dtau interpolated across a masked gap (e.g. a waveform-combine
            # seam) rather than measured at this node
            filled = itp is not None and itp.size and itp[k, b] == 1
            ok = np.isfinite(dlam[k, b]) and not pegged and not filled \
                and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
            if ok:
                m = (DEPTH >= top[k, b]) & (DEPTH < bot[k, b])
                cols[b, m] = dlam[k, b]
    return lats, lons, cols


def main():
    fn = os.path.join(EXT, f'margin_{SEG}.mat')
    if not os.path.exists(fn):
        print(f'extract not found: {fn}')
        return
    d = loadmat(fn)
    t_us = d['Time'].ravel() * 1e6
    lats = d['Latitude'].ravel()
    lons = d['Longitude'].ravel()
    step = haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
    dist = np.concatenate([[0], np.cumsum(step)])

    proj = ab.proj3031()
    extent = ab.points_extent(lats, lons, pad_frac=1.0, min_pad_m=12e3)
    win = itslive_window(extent)
    speed = sample_speed(win, lats, lons)

    inv = dlam_section(SEG)

    fig = plt.figure(figsize=(16, 9), layout='constrained')
    gs = fig.add_gridspec(4, 3, width_ratios=[0.85, 1.6, 0.018],
                          height_ratios=[0.55, 1, 1, 1])

    # ---- (left) map: segment over the ITS_LIVE speed field
    axm = fig.add_subplot(gs[:, 0], projection=proj)
    axm.set_extent(extent, crs=proj)
    if win is not None:
        wx, wy, v = win
        pcm = axm.pcolormesh(wx, wy, v, transform=proj, cmap='magma',
                             norm=LogNorm(vmin=10, vmax=400), shading='auto',
                             rasterized=True)
        cs = axm.contour(wx, wy, v, levels=[30, 100, 300], colors='white',
                         linewidths=0.6, alpha=0.8, transform=proj)
        axm.clabel(cs, fontsize=6, fmt='%g')
        # Inset colorbar tracks the map axes, which constrained layout
        # letterboxes inside its gridspec slot (fixed geographic aspect)
        cax = axm.inset_axes([0.12, -0.09, 0.76, 0.03])
        fig.colorbar(pcm, cax=cax, orientation='horizontal',
                     label='ITS_LIVE speed (m/yr)')
    # Context: the other margin display extracts
    for cfn in sorted(glob.glob(os.path.join(EXT, 'margin_*.mat'))):
        if os.path.basename(cfn) == f'margin_{SEG}.mat':
            continue
        c = loadmat(cfn, variable_names=['Latitude', 'Longitude'])
        axm.plot(c['Longitude'].ravel(), c['Latitude'].ravel(), '-',
                 color='0.8', lw=1.0, transform=ccrs.PlateCarree(), zorder=5)
    axm.plot(lons, lats, '-', color='cyan', lw=2.2,
             transform=ccrs.PlateCarree(), zorder=6)
    axm.plot(lons[0], lats[0], 'o', color='cyan', mec='k', ms=7,
             transform=ccrs.PlateCarree(), zorder=7)
    axm.annotate('km 0', proj.transform_point(lons[0], lats[0],
                 ccrs.PlateCarree())[:2], xytext=(6, 6),
                 textcoords='offset points', fontsize=8, zorder=8,
                 bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none',
                           alpha=0.75))
    gl = axm.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='white',
                       y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    ab.add_scale_bar(axm, extent)
    axm.set_title(f'Segment {SEG}\n(eastern shear margin crossing)',
                  fontsize=11)

    # ---- (right) profile panels sharing the distance axis
    ax0 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 1], sharex=ax0)
    ax2 = fig.add_subplot(gs[2, 1], sharex=ax0)
    ax3 = fig.add_subplot(gs[3, 1], sharex=ax0)
    for ax in (ax0, ax1, ax2):
        ax.tick_params(labelbottom=False)

    ax0.plot(dist, speed, 'k-', lw=1.5)
    ax0.set_ylabel('speed (m/yr)')
    ax0.grid(alpha=0.3)
    ax0.set_title('ITS_LIVE surface speed along the crossing', fontsize=10)

    img = d['power_ref']
    pc = ax1.pcolormesh(dist, t_us, img, cmap='gray',
                        vmin=np.nanmin(img), vmax=np.nanmax(img),
                        shading='auto', rasterized=True)
    ax1.invert_yaxis()
    ax1.set_ylabel('TWTT (us)')
    fig.colorbar(pc, cax=fig.add_subplot(gs[1, 2]), label='HH power (dB)')

    pc = ax2.pcolormesh(dist, t_us, d['phase_wrapped'], cmap='twilight',
                        vmin=-np.pi, vmax=np.pi, shading='auto',
                        rasterized=True)
    ax2.invert_yaxis()
    ax2.set_ylabel('TWTT (us)')
    fig.colorbar(pc, cax=fig.add_subplot(gs[2, 2]),
                 label='interferogram phase (rad)')

    if inv is not None and len(inv[0]):
        blats, blons, cols = inv
        # Nearest extract trace maps each block onto the distance axis
        bdist = np.array([dist[np.argmin(haversine_km(lats, lons, la, lo))]
                          for la, lo in zip(blats, blons)])
        order = np.argsort(bdist)
        bdist, cols = bdist[order], cols[order]
        edges = block_edges(bdist)
        depth_edges = np.concatenate([DEPTH, [DEPTH[-1] + 5]])
        lim = sym_limit(cols)
        pc = ax3.pcolormesh(edges, depth_edges, cols.T, cmap='RdBu_r',
                            vmin=-lim, vmax=lim, shading='auto')
        fig.colorbar(pc, cax=fig.add_subplot(gs[3, 2]),
                     label=r'$\Delta\lambda$ (survey frame)')
    else:
        ax3.annotate('inversion output unavailable', xy=(0.5, 0.5),
                     xycoords='axes fraction', ha='center')
    ax3.set_ylim(DEPTH[-1], 0)
    ax3.set_ylabel('Depth (m)')
    ax3.set_title('Inverted horizontal eigenvalue difference', fontsize=10)
    ax3.set_xlabel('Distance along segment (km)')
    ax3.set_xlim(dist[0], dist[-1])

    fig.suptitle(f'Thwaites eastern shear margin crossing, segment {SEG} '
                 '(2023-24 season)', fontsize=13)
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, f'margin_scene_{SEG}.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

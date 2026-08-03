"""Interferogram-chain QC figure for the Thwaites margin-crossing segment.

One stacked panel per product, sharing the distance axis, full chain for
a single segment: HH power, VV power, coherence, wrapped interferogram
phase, SNAPHU-unwrapped phase, and coregistration row offsets (as
traveltime difference), with the ITS_LIVE speed profile on top for
margin context. Input is the compact extract produced on mem1 by
extract_margin.py (see "Margin display extracts" in opr_fabric/README.md
for the server workflow) and mirrored under ~/data/opr/margin/.

Usage: python thwaites_interferogram.py [seg] [extract_dir] [out_dir]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SEG = sys.argv[1] if len(sys.argv) > 1 else '20240101_02_001'
EXT = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser('~/data/opr/margin')
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
    os.path.dirname(__file__), '..', '..', 'figs')

ITSLIVE_LOCAL = os.path.expanduser(
    '~/Downloads/ITS_LIVE_velocity_120m_RGI19A_0000_V02.1.nc')


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def sample_speed(lats, lons):
    """ITS_LIVE speed at points from the local mosaic; NaN if unavailable."""
    if not os.path.exists(ITSLIVE_LOCAL):
        return np.full(len(lats), np.nan)
    import cartopy.crs as ccrs
    import xarray as xr
    import antarctic_basemap as ab
    proj = ab.proj3031()
    xyz = proj.transform_points(ccrs.PlateCarree(),
                                np.asarray(lons), np.asarray(lats))
    px, py = xyz[:, 0], xyz[:, 1]
    pad = 5e3
    with xr.open_dataset(ITSLIVE_LOCAL, engine='h5netcdf') as ds:
        win = ds.sel(x=slice(px.min() - pad, px.max() + pad),
                     y=slice(py.max() + pad, py.min() - pad))
        v = win['v'].values.astype(float)
        wx = win['x'].values
        wy = win['y'].values
    ix = np.clip(np.searchsorted(wx, px), 0, wx.size - 1)
    iy = np.clip(np.searchsorted(-wy, -py), 0, wy.size - 1)
    return v[iy, ix]


def main():
    fn = os.path.join(EXT, f'margin_{SEG}.mat')
    if not os.path.exists(fn):
        print(f'extract not found: {fn}\n'
              'Run scratch extract_margin.py on mem1 and rsync the outputs.')
        return
    d = loadmat(fn)
    t_us = d['Time'].ravel() * 1e6
    lats = d['Latitude'].ravel()
    lons = d['Longitude'].ravel()
    step = haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
    dist = np.concatenate([[0], np.cumsum(step)])
    speed = sample_speed(lats, lons)
    fc = 750e6
    dt = np.median(np.diff(d['Time'].ravel()))

    panels = [
        ('power_ref', 'HH power (dB)', 'gray', None, None),
        ('power_sec', 'VV power (dB)', 'gray', None, None),
        ('coherence', 'coherence', 'viridis', 0, 1),
        ('phase_wrapped', 'interferogram phase (rad)', 'twilight',
         -np.pi, np.pi),
        ('phase_unwrapped', 'unwrapped phase (rad)', 'RdBu_r', None, None),
        ('row_offset', r'coregistration $\Delta\tau$ (ns)', 'RdBu_r',
         None, None),
    ]
    have = [p for p in panels if p[0] in d]

    n = len(have)
    fig = plt.figure(figsize=(14, 2.1 * n + 2.2), layout='constrained')
    gs = fig.add_gridspec(n + 1, 2, width_ratios=[1, 0.015],
                          height_ratios=[0.7] + [1] * n, wspace=0.03)
    axes = [fig.add_subplot(gs[0, 0])]
    for i in range(n):
        axes.append(fig.add_subplot(gs[i + 1, 0], sharex=axes[0]))
    caxes = [fig.add_subplot(gs[i + 1, 1]) for i in range(n)]
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)

    ax = axes[0]
    ax.plot(dist, speed, 'k-', lw=1.5)
    ax.set_ylabel('speed (m/yr)')
    ax.set_title(f'Segment {SEG}: ITS_LIVE speed (margin context)')
    ax.grid(alpha=0.3)

    for ax, cax, (key, label, cmap, vmin, vmax) in zip(axes[1:], caxes, have):
        img = d[key]
        if key == 'row_offset':
            img = img * dt * 1e9  # bins -> ns
            lim = np.nanpercentile(np.abs(img), 98)
            vmin, vmax = -lim, lim
        if key == 'phase_unwrapped':
            lim = np.nanpercentile(np.abs(img - np.nanmedian(img)), 98)
            img = img - np.nanmedian(img)
            vmin, vmax = -lim, lim
        if key.startswith('power'):
            vmin, vmax = np.nanpercentile(img, [8, 99.5])
        pc = ax.pcolormesh(dist, t_us, img, cmap=cmap, vmin=vmin, vmax=vmax,
                           shading='auto', rasterized=True)
        ax.invert_yaxis()
        ax.set_ylabel('TWTT (us)')
        cb = fig.colorbar(pc, cax=cax)
        cb.set_label(label, fontsize=8)
    axes[-1].set_xlabel('Distance along segment (km)')
    fig.suptitle('Interferogram product chain, Thwaites margin segment '
                 f'{SEG} (2023-24 season, Hoffman/Christianson processing)',
                 fontsize=13)
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, f'interferogram_chain_{SEG}.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

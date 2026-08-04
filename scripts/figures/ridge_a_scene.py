"""Ridge A scene figure for one grid leg (default frame 20250108_02_009).

Companion to thwaites_margin_scene.py at the site where the fabric
signal is cleanest. Left: EPSG:3031 map of the Ridge A raster grid
(every 2024-25 joint-inversion block) on a plain background, with the
scene frame highlighted; MOA over the featureless dome is pure sensor
noise once stretched, so no imagery is drawn, and ITS_LIVE has no
coverage this far south, so there is no speed layer. Right: four stacked panels sharing the
along-frame distance axis - surface/bed elevation (divide context, in
place of the Thwaites speed panel), HH power, wrapped interferogram
phase, and the inverted horizontal eigenvalue difference.

Inputs are the locally mirrored CSARP products for the frame
(~/data/opr/accum/<season>/CSARP_standard_HH and
CSARP_polarimetric_unwrap) and the joint-inversion Data_*.mat under the
fabric batch root (Paden processing, as in ridge_a_azimuthal.py).

Usage: python ridge_a_scene.py [frame] [accum_root] [fabric_root] [out_dir]
"""
import glob
import os
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

import cartopy.crs as ccrs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab
from fabric_qc import interpolated_intervals

FRAME = sys.argv[1] if len(sys.argv) > 1 else '20250108_02_009'
ACCUM = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
    '~/data/opr/accum/2024_Antarctica_Ground2')
ROOT = sys.argv[3] if len(sys.argv) > 3 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[4] if len(sys.argv) > 4 else os.path.join(
    os.path.dirname(__file__), '..', '..', 'figs')

DEPTH = np.arange(0, 1900, 5.0)
C_AIR = 299792458.0
C_ICE = C_AIR / 1.78


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


def bearing_deg(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(x, y)) % 360.0


def cumdist_km(lats, lons):
    step = haversine_km(lats[:-1], lons[:-1], lats[1:], lons[1:])
    return np.concatenate([[0], np.cumsum(step)])


def load_standard_hh(frame):
    """HH echogram (v7.3): (lats, lons, t_us, power_db, surf, bot)."""
    seg = frame.rsplit('_', 1)[0]
    fn = f'{ACCUM}/CSARP_standard_HH/{seg}/Data_{frame}.mat'
    with h5py.File(fn) as f:
        # MATLAB v7.3: Data rows are traces, columns time samples
        p = f['Data'][:]
        t = f['Time'][0, :]
        lats = f['Latitude'][:, 0]
        lons = f['Longitude'][:, 0]
        surf = f['Surface'][:, 0]
        bot = f['Bottom'][:, 0]
    with np.errstate(divide='ignore'):
        pdb = 10.0 * np.log10(p)
    return lats, lons, t * 1e6, pdb.T, surf, bot


def load_unwrap(frame):
    """(lats, lons, t_us, wrapped phase) from the polarimetric product."""
    seg = frame.rsplit('_', 1)[0]
    fn = f'{ACCUM}/CSARP_polarimetric_unwrap/{seg}/Data_{frame}.mat'
    d = loadmat(fn, variable_names=['interferogram_mlook', 'Time',
                                    'Latitude', 'Longitude'])
    return (d['Latitude'].ravel(), d['Longitude'].ravel(),
            d['Time'].ravel() * 1e6, np.angle(d['interferogram_mlook']))


def dlam_section(frame):
    """(lats, lons, cols) inversion blocks gridded on DEPTH; None if absent.

    Masking follows thwaites_fabric.py (pegged/clipped and low-quality
    bands dropped).
    """
    seg = frame.rsplit('_', 1)[0]
    fns = sorted(glob.glob(
        f'{ROOT}/joint_jp/2024_Antarctica_Ground2/{seg}/Data_{frame}.mat'))
    if not fns:
        return None
    d = loadmat(fns[0])
    dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
    qual = d['dlam_quality']
    clip = d.get('dlam_clipped')
    itp = interpolated_intervals(d)
    lats = d['Latitude'][0, :]
    lons = d['Longitude'][0, :]
    cols = np.full((dlam.shape[1], DEPTH.size), np.nan)
    for b in range(dlam.shape[1]):
        for k in range(dlam.shape[0]):
            pegged = (clip is not None and clip.size and clip[k, b] == 1) \
                or abs(dlam[k, b]) > 0.6
            filled = itp is not None and itp[k, b]
            ok = np.isfinite(dlam[k, b]) and not pegged and not filled \
                and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
            if ok:
                m = (DEPTH >= top[k, b]) & (DEPTH < bot[k, b])
                cols[b, m] = dlam[k, b]
    return lats, lons, cols


def grid_blocks():
    """Lat/lon of every Ridge A joint-inversion block (map context)."""
    lats, lons = [], []
    for fn in sorted(glob.glob(
            f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat')):
        d = loadmat(fn, variable_names=['Latitude', 'Longitude'])
        lats.extend(d['Latitude'][0, :])
        lons.extend(d['Longitude'][0, :])
    la = np.array(lats)
    lo = np.array(lons)
    ok = np.isfinite(la) & np.isfinite(lo)
    return la[ok], lo[ok]


def main():
    lats_h, lons_h, t_h, pdb, surf_twtt, bot_twtt = load_standard_hh(FRAME)
    lats_p, lons_p, t_p, phase = load_unwrap(FRAME)
    dist_h = cumdist_km(lats_h, lons_h)
    dist_p = cumdist_km(lats_p, lons_p)
    inv = dlam_section(FRAME)
    head = bearing_deg(lats_h[0], lons_h[0], lats_h[-1], lons_h[-1])

    # Divide-context elevations from the radar picks (sled at the surface)
    surf_elev = np.full(surf_twtt.shape, np.nan)
    bed_elev = np.full(surf_twtt.shape, np.nan)
    ok = np.isfinite(surf_twtt)
    with h5py.File(f'{ACCUM}/CSARP_standard_HH/'
                   f'{FRAME.rsplit("_", 1)[0]}/Data_{FRAME}.mat') as f:
        gps_elev = f['Elevation'][:, 0]
    surf_elev[ok] = gps_elev[ok] - surf_twtt[ok] * C_AIR / 2
    ok = np.isfinite(surf_twtt) & np.isfinite(bot_twtt)
    bed_elev[ok] = surf_elev[ok] - (bot_twtt[ok] - surf_twtt[ok]) * C_ICE / 2

    proj = ab.proj3031()
    gla, glo = grid_blocks()
    extent = ab.points_extent(np.concatenate([gla, lats_h]),
                              np.concatenate([glo, lons_h]),
                              pad_frac=0.35, min_pad_m=5e3)

    fig = plt.figure(figsize=(16, 9), layout='constrained')
    gs = fig.add_gridspec(4, 3, width_ratios=[0.85, 1.6, 0.018],
                          height_ratios=[0.55, 1, 1, 1])

    # ---- (left) map: the raster grid with the scene frame highlighted
    axm = fig.add_subplot(gs[:, 0], projection=proj)
    axm.set_extent(extent, crs=proj)
    # No imagery underlay: MOA over the featureless dome is pure sensor
    # noise after the contrast stretch, which just buries the grid
    axm.set_facecolor('0.94')
    axm.scatter(glo, gla, s=5, color='tab:orange', alpha=0.7,
                transform=ccrs.PlateCarree(), zorder=5,
                label='grid inversion blocks')
    axm.plot(lons_h, lats_h, '-', color='cyan', lw=2.5,
             transform=ccrs.PlateCarree(), zorder=6, label=f'frame {FRAME}')
    axm.plot(lons_h[0], lats_h[0], 'o', color='cyan', mec='k', ms=7,
             transform=ccrs.PlateCarree(), zorder=7)
    axm.annotate('km 0', proj.transform_point(lons_h[0], lats_h[0],
                 ccrs.PlateCarree())[:2], xytext=(6, 6),
                 textcoords='offset points', fontsize=8, zorder=8,
                 bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none',
                           alpha=0.75))
    gl = axm.gridlines(draw_labels=True, lw=0.4, alpha=0.5, y_inline=False,
                       color='gray')
    gl.top_labels = False
    gl.right_labels = False
    ab.add_scale_bar(axm, extent)
    axm.legend(fontsize=8, loc='upper right', markerscale=2)
    axm.set_title(f'Ridge A raster grid, frame {FRAME}\n'
                  f'(drive heading {head:.0f}\N{DEGREE SIGN} E of N)',
                  fontsize=11)

    # ---- (right) profile panels sharing the distance axis
    ax0 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 1], sharex=ax0)
    ax2 = fig.add_subplot(gs[2, 1], sharex=ax0)
    ax3 = fig.add_subplot(gs[3, 1], sharex=ax0)
    for ax in (ax0, ax1, ax2):
        ax.tick_params(labelbottom=False)

    ax0.plot(dist_h, surf_elev, 'k-', lw=1.2, label='surface')
    if np.any(np.isfinite(bed_elev)):
        ax0.plot(dist_h, bed_elev, '-', color='saddlebrown', lw=1.2,
                 label='bed')
    else:
        ax0.annotate('no bed pick in this frame', xy=(0.99, 0.9),
                     xycoords='axes fraction', ha='right', va='top',
                     fontsize=7, color='gray')
    ax0.set_ylabel('elevation (m WGS84)')
    ax0.grid(alpha=0.3)
    ax0.legend(fontsize=7, loc='center left')
    ax0.set_title('Divide setting: surface elevation along the frame '
                  '(no ITS_LIVE coverage this far south)', fontsize=10)

    # No bed pick in this frame: bound TWTT by the inversion depth range
    t_max = np.nanmedian(surf_twtt) * 1e6 + 2.2 * DEPTH[-1] / C_ICE * 1e6
    pc = ax1.pcolormesh(dist_h, t_h, pdb, cmap='gray',
                        vmin=np.nanpercentile(pdb, 5),
                        vmax=np.nanmax(pdb), shading='auto', rasterized=True)
    ax1.set_ylim(t_max, 0)
    ax1.set_ylabel('TWTT (us)')
    fig.colorbar(pc, cax=fig.add_subplot(gs[1, 2]), label='HH power (dB)')

    pc = ax2.pcolormesh(dist_p, t_p, phase, cmap='twilight',
                        vmin=-np.pi, vmax=np.pi, shading='auto',
                        rasterized=True)
    ax2.set_ylim(t_max, 0)
    ax2.set_ylabel('TWTT (us)')
    fig.colorbar(pc, cax=fig.add_subplot(gs[2, 2]),
                 label='interferogram phase (rad)')

    if inv is not None and len(inv[0]):
        blats, blons, cols = inv
        bdist = np.array([dist_p[np.argmin(haversine_km(lats_p, lons_p,
                                                        la, lo))]
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
    ax3.set_xlabel('Distance along frame (km)')
    ax3.set_xlim(dist_p[0], dist_p[-1])

    fig.suptitle(f'Ridge A fabric scene, frame {FRAME} '
                 '(2024-25 season, Paden processing)', fontsize=13)
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, f'ridge_a_scene_{FRAME}.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

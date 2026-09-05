"""Two-panel SCAR site figures: survey overview + interferogram.

One PNG per site:
  left   the survey in polar stereographic over log-scale ITS_LIVE
         surface speed (where covered), all survey lines in WHITE, the
         focused profile as a BLACK line, flow arrows, a continental
         locator inset and a scale bar in the bottom-right corner
  right  the wrapped polarimetric interferogram of the focused profile

Slide styling, the end markers and the shared panel geometry live in
scar_style.py, so this figure and the NEGIS one cannot drift apart.

The black profile line is drawn over a thin white casing so it stays
legible where the speed field is dark (slow ice is exactly where these
profiles sit); the casing reads as an outline, not as a second track.

Sites rendered here: Thwaites eastern shear margin (ITS_LIVE v2.1 window
streamed from S3; the panel shows that WHOLE window rather than a box
padded round the tracks, so the trunk the survey sits beside stays in
frame) and Ridge A (south of satellite velocity coverage, hence a plain
dark background and no speed colorbar - their absence is what says there
is no satellite velocity here, in place of an in-map caption).

Usage: python scar_two_panel.py <out_dir>
"""
import glob
import os
import sys
import warnings

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib.colors import LogNorm      # noqa: E402
from scipy.io import loadmat               # noqa: E402

warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab            # noqa: E402
import cartopy.crs as ccrs                # noqa: E402
import scar_style as sty                  # noqa: E402
from scar_style import (DATA, DPI, PROFILE_C, PROFILE_FX,  # noqa: E402
                        PROFILE_LW, SPEED_CB_LABEL, TRACK_C, zoom_inset)

OUT = sty.out_dir(sys.argv)

ACCUM = os.path.expanduser('~/data/opr/accum/2024_Antarctica_Ground2')
BATCH = os.path.expanduser('~/data/opr/fabric_batch')
MARGIN = os.path.expanduser('~/data/opr/margin')

C_ICE = 299792458.0 / np.sqrt(3.171)
GREY = np.array([0.74, 0.74, 0.74])


def fabric_panel(ax, fig, tag, dist, lats, lons, t_us, surf_t_s, vmax,
                 cb_label=r'$\Delta\lambda$ (LS, principal axes)'):
    """The depth-resolved fabric product on the interferogram panel's
    exact geometry: same distance axis, same absolute-TWTT vertical axis
    and limits, same colorbar placement - so the wrapped-phase figure and
    this one overlay slide-for-slide. Cells fade toward grey by their own
    coherence (the fabric grammar of the section figures) and abstained
    cells are white; the fabric occupies the constrained upper record and
    the canvas below stays empty exactly where the phase panel shows the
    incoherent deep record."""
    # prefer a deep variant (full-record z_max rerun) when one is staged
    cand = sorted(glob.glob(os.path.join(DATA,
                                         'quadpol_section_%s_z*.mat' % tag)))
    sec_fn = cand[-1] if cand else os.path.join(
        DATA, 'quadpol_section_%s.mat' % tag)
    with h5py.File(sec_fn) as f:
        r = f['res']
        z = np.array(r['z']).ravel().astype(float)
        dl = np.atleast_2d(np.array(r['sec_dlam_ls']).T).astype(float)
        cm = np.atleast_2d(np.array(r['sec_cmag']).T).astype(float)
        blat = np.array(r['sec_lat']).ravel()
        blon = np.array(r['sec_lon']).ravel()
    # place each block on the interferogram's own distance axis by
    # nearest trace, so the two panels use one x definition
    bdist = np.empty(blat.size)
    for b in range(blat.size):
        j = np.nanargmin((lats - blat[b])**2 +
                         (np.cos(np.deg2rad(blat[b]))*(lons - blon[b]))**2)
        bdist[b] = dist[j]
    o = np.argsort(bdist)
    dl, cm, bdist = dl[:, o], cm[:, o], bdist[o]
    t_z = (surf_t_s + 2.0*z/C_ICE) * 1e6
    w = np.clip((cm - 0.10)/(0.45 - 0.10), 0, 1)**0.7
    w[cm < 0.25] *= 0.35
    cmap = plt.get_cmap('Blues')
    rgb = cmap(np.clip(dl, 0, vmax)/vmax)[..., :3]
    ww = np.clip(w, 0, 1)[..., None]
    rgb = rgb*ww + GREY[None, None, :]*(1 - ww)
    rgb[~np.isfinite(dl)] = 1.0
    ax.imshow(rgb, aspect='auto', interpolation='nearest',
              extent=[bdist[0], bdist[-1], t_z[-1], t_z[0]])
    ax.set_xlim(dist[0], dist[-1])
    ax.set_ylim(t_us[-1], t_us[0])
    ax.set_xlabel('distance along profile (km)')
    ax.set_ylabel('TWTT (\u03bcs)')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, ax=ax, pad=0.035)
    cb.set_label(cb_label)


def two_panel(map_fn, ifg_fn, out):
    fig, axm, axi = sty.two_panel_figure(ab.proj3031())
    map_fn(fig, axm)
    ifg_fn(fig, axi)
    fig.savefig(os.path.join(OUT, out), dpi=DPI)
    plt.close(fig)
    print('wrote', out)


# ---------------------------------------------------------------- Thwaites
def thwaites():
    seg = '20240108_01_001'
    d = loadmat(os.path.join(MARGIN, f'margin_{seg}.mat'),
                variable_names=['phase_wrapped', 'coherence', 'Time',
                                'Latitude', 'Longitude', 'Surface'])
    lats, lons = d['Latitude'].ravel(), d['Longitude'].ravel()
    dist = sty.cumdist_km(lats, lons)

    tracks = []
    for fn in sorted(glob.glob(
            f'{BATCH}/joint/2023_Antarctica_Ground/*/Data_*.mat')):
        g = loadmat(fn, variable_names=['Latitude', 'Longitude'])
        la, lo = g['Latitude'][0, :], g['Longitude'][0, :]
        keep = (np.isfinite(la) & np.isfinite(lo) &
                (lo > -110) & (lo < -98) & (la > -78))
        if keep.any():
            tracks.append((la[keep], lo[keep]))

    z = np.load(os.path.join(DATA, 'itslive_thwaites_win.npz'))
    wx, wy, v = z['x'], z['y'], z['v']
    vx, vy = z['vx'], z['vy']

    proj = ab.proj3031()
    # The whole ITS_LIVE window rather than a padded box round the survey:
    # at 41 x 113 km the tracks occupy a fraction of the 180 x 240 km
    # window, and padding to them cropped away most of the Thwaites trunk
    # the survey exists to sit beside. The zoom inset and the continental
    # locator carry the other two scales.
    extent = (wx.min(), wx.max(), wy.min(), wy.max())

    norm = LogNorm(vmin=5, vmax=600)

    def speed_layer(ax, dec):
        pcm = ax.pcolormesh(wx, wy, v, transform=proj, cmap='magma',
                            norm=norm, shading='auto', rasterized=True)
        Xq, Yq = np.meshgrid(wx[::dec], wy[::dec])
        U, V = vx[::dec, ::dec], vy[::dec, ::dec]
        sp = np.hypot(U, V)
        ax.quiver(Xq, Yq, U / sp, V / sp, transform=proj, color='white',
                  scale=26, width=0.005, alpha=0.6, zorder=6)
        return pcm

    def map_panel(fig, ax):
        ax.set_extent(extent, crs=proj)
        pcm = speed_layer(ax, dec=110)
        # Below the rotated longitude labels: at this latitude they slant
        # far enough off the bottom axis to sit on a colorbar placed at
        # the usual -0.07
        cax = ax.inset_axes([0.05, -0.15, 0.62, 0.03])
        fig.colorbar(pcm, cax=cax, orientation='horizontal',
                     label=SPEED_CB_LABEL)
        for la, lo in tracks:
            ax.plot(lo, la, '-', color=TRACK_C, lw=1.1, alpha=0.9,
                    transform=ccrs.PlateCarree(), zorder=7)
        ax.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
                transform=ccrs.PlateCarree(), zorder=8,
                path_effects=PROFILE_FX)
        gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5,
                          color='gray', y_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        ab.scale_bar_br(ax, extent, ab.proj3031())
        # Hard right in the top corner: the zoom window sits high and
        # slightly left of centre here, so at 0.665 the locator painted
        # over the top-right corner of its dashed box. The locator keeps
        # its square aspect and is centred in the rect, which leaves the
        # visual margin inside the panel edge.
        sty.continental_inset(ax, 'antarctica',
                              (0.695, 0.775, 0.29, 0.20), extent, proj)
        ax.set_title('Thwaites Glacier survey', fontsize=12)

        # Zoom inset: at survey scale the 10 km profile is a blob, and
        # the end markers exist to orient the interferogram against the
        # flow, so they need a window where the two ends separate
        px, py = proj.transform_point(np.mean(lons), np.mean(lats),
                                      ccrs.PlateCarree())[:2]
        half = 9e3
        zext = (px - half, px + half, py - half, py + half)
        axz = zoom_inset(ax, [0.02, 0.30, 0.55, 0.34], proj)
        axz.set_extent(zext, crs=proj)
        speed_layer(axz, dec=18)
        axz.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
                 transform=ccrs.PlateCarree(), zorder=8,
                 path_effects=PROFILE_FX)
        sty.draw_map_ends(axz, lats, lons)
        sty.inset_frame(axz)
        sty.locator_box(ax, zext, proj)

    def ifg(fig, ax):
        sty.ifg_panel(ax, fig, dist, d['Time'].ravel() * 1e6,
                      d['phase_wrapped'], np.abs(d['coherence']))
        sty.draw_ifg_ends(ax, dist)
        ax.set_title(f'eastern shear margin crossing, segment {seg}',
                     fontsize=12, pad=26)

    two_panel(map_panel, ifg, 'scar_thwaites_2panel.png')

    def fab(fig, ax):
        fabric_panel(ax, fig, seg, dist, lats, lons, d['Time'].ravel()*1e6,
                     float(np.nanmedian(d['Surface'])), vmax=0.15)
        sty.draw_ifg_ends(ax, dist)
        ax.set_title(f'eastern shear margin crossing, segment {seg} - '
                     'depth-resolved fabric', fontsize=12, pad=26)

    two_panel(map_panel, fab, 'scar_thwaites_2panel_fabric.png')


# ---------------------------------------------------------------- Ridge A
def ridge_a():
    frame = '20250108_02_009'
    seg = frame.rsplit('_', 1)[0]
    d = loadmat(f'{ACCUM}/CSARP_polarimetric_unwrap/{seg}/Data_{frame}.mat',
                variable_names=['interferogram_mlook',
                                'interferogram_coherence', 'Time',
                                'Latitude', 'Longitude', 'Surface'])
    lats, lons = d['Latitude'].ravel(), d['Longitude'].ravel()
    dist = sty.cumdist_km(lats, lons)

    # Tracks drawn as lines like every other site. The fabric product
    # carries only a few block centres per frame, so the legs have to be
    # rebuilt per SEGMENT - one frame's three points would draw as a dash,
    # and a scatter made Ridge A look like a different kind of survey
    # than Thwaites or NEGIS.
    by_seg = {}
    for fn in sorted(glob.glob(
            f'{BATCH}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat')):
        g = loadmat(fn, variable_names=['Latitude', 'Longitude'])
        la, lo = np.atleast_1d(g['Latitude'].ravel()), \
            np.atleast_1d(g['Longitude'].ravel())
        keep = np.isfinite(la) & np.isfinite(lo)
        if keep.any():
            by_seg.setdefault(os.path.basename(os.path.dirname(fn)),
                              []).append((la[keep], lo[keep]))
    tracks = [(np.concatenate([p[0] for p in v]),
               np.concatenate([p[1] for p in v]))
              for v in by_seg.values()]
    gla = np.concatenate([t[0] for t in tracks])
    glo = np.concatenate([t[1] for t in tracks])

    proj = ab.proj3031()
    extent = ab.points_extent(np.concatenate([gla, lats]),
                              np.concatenate([glo, lons]),
                              pad_frac=0.30, min_pad_m=5e3)

    def map_panel(fig, ax):
        ax.set_extent(extent, crs=proj)
        # Ridge A has no velocity field to sit on, so the survey needs a
        # background dark enough for white tracks to read against
        ax.set_facecolor('0.35')
        for la, lo in tracks:
            ax.plot(lo, la, '-', color=TRACK_C, lw=1.1, alpha=0.9,
                    transform=ccrs.PlateCarree(), zorder=5)
        ax.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
                transform=ccrs.PlateCarree(), zorder=8,
                path_effects=PROFILE_FX)
        # Zoom inset so the two ends separate (same reasoning as Thwaites)
        px, py = proj.transform_point(np.mean(lons), np.mean(lats),
                                      ccrs.PlateCarree())[:2]
        half = 4e3
        zext = (px - half, px + half, py - half, py + half)
        axz = zoom_inset(ax, [0.02, 0.30, 0.5, 0.32], proj)
        axz.set_extent(zext, crs=proj)
        axz.set_facecolor('0.28')
        for la, lo in tracks:
            axz.plot(lo, la, '-', color=TRACK_C, lw=1.6, alpha=0.9,
                     transform=ccrs.PlateCarree(), zorder=5)
        axz.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
                 transform=ccrs.PlateCarree(), zorder=8,
                 path_effects=PROFILE_FX)
        sty.draw_map_ends(axz, lats, lons)
        sty.inset_frame(axz)
        sty.locator_box(ax, zext, proj)
        gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5,
                          color='gray', y_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        ab.scale_bar_br(ax, extent, ab.proj3031())
        sty.continental_inset(ax, 'antarctica',
                              (0.02, 0.755, 0.26, 0.225), extent, proj)
        # No in-map caption, boxed or loose - deliberate, do not re-add. The
        # absence of a speed layer and its colorbar already says there is no
        # satellite velocity here, the locator inset puts the site deep in
        # the interior, and the talk says it. On a projected slide the text
        # competed with the survey lines, which are the panel's subject.
        ax.set_title('Ridge A raster survey', fontsize=12)

    def ifg(fig, ax):
        sty.ifg_panel(ax, fig, dist, d['Time'].ravel() * 1e6,
                      np.angle(d['interferogram_mlook']),
                      np.abs(d['interferogram_coherence']))
        sty.draw_ifg_ends(ax, dist)
        ax.set_title(f'divide setting, frame {frame}', fontsize=12, pad=26)

    two_panel(map_panel, ifg, 'scar_ridge_a_2panel.png')

    def fab(fig, ax):
        fabric_panel(ax, fig, frame, dist, lats, lons, d['Time'].ravel()*1e6,
                     float(np.nanmedian(d['Surface'])), vmax=0.08)
        sty.draw_ifg_ends(ax, dist)
        ax.set_title(f'divide setting, frame {frame} - depth-resolved fabric',
                     fontsize=12, pad=26)

    two_panel(map_panel, fab, 'scar_ridge_a_2panel_fabric.png')


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    thwaites()
    ridge_a()

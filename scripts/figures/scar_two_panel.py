"""Two-panel SCAR site figures: survey overview + interferogram.

One PNG per site:
  left   the ENTIRE survey in polar stereographic over log-scale ITS_LIVE
         surface speed (where covered), all survey lines in WHITE, the
         focused profile as a BLACK line, flow arrows, and a scale bar in
         the bottom-right corner
  right  the wrapped polarimetric interferogram of the focused profile

Slide styling (Andrew, 4 Aug): no legends - the titles already name the
segment, and a legend box eats map area on a projected slide. The
profile's two ENDS carry the map-to-interferogram correspondence
instead of a distance-marker ladder: a red marker at the start and a
white marker at the end, both with a black outline, drawn on the map
inset AND at the matching left/right edges of the interferogram. Two
unambiguous ends read from the back of a room where five graded viridis
dots do not.

The black profile line is drawn over a thin white casing so it stays
legible where the speed field is dark (slow ice is exactly where these
profiles sit); the casing reads as an outline, not as a second track.

Sites rendered here: Thwaites eastern shear margin (ITS_LIVE v2.1
window streamed from S3) and Ridge A (south of satellite velocity
coverage; plain background, stated on the panel).

Usage: python scar_two_panel.py <out_dir>
"""
import glob
import os
import sys
import warnings

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.patheffects as pe        # noqa: E402
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib import transforms          # noqa: E402
from matplotlib.colors import LogNorm      # noqa: E402
from scipy.io import loadmat               # noqa: E402

warnings.filterwarnings('ignore')

REPO = os.path.expanduser('~/projects/fabric_anisotropy')
sys.path.insert(0, os.path.join(REPO, 'scripts', 'figures'))
import antarctic_basemap as ab            # noqa: E402
import cartopy.crs as ccrs                # noqa: E402

OUT = sys.argv[1]
# Inputs that are data, not code: ITS_LIVE windows streamed from S3, the
# NEGIS track/interferogram extracts pulled off mem1. They live beside the
# rest of the mirrored products rather than in the repo.
DATA = os.path.expanduser(os.environ.get('SCAR_DATA', '~/data/opr/scar'))

ACCUM = os.path.expanduser('~/data/opr/accum/2024_Antarctica_Ground2')
BATCH = os.path.expanduser('~/data/opr/fabric_batch')
MARGIN = os.path.expanduser('~/data/opr/margin')

FIGSIZE = (13.0, 5.8)
DPI = 200

# Profile endpoints: start red, end white, both outlined black. The same
# pair marks the interferogram's left and right edges, which is what
# makes the correspondence readable without a legend.
END_FACES = ('#e8000b', '#ffffff')
END_SIZE = 130
END_EDGE = 1.6

TRACK_C = 'white'       # every survey line
PROFILE_C = 'black'     # the focused profile
PROFILE_LW = 2.8
# Thin white casing under the black profile line, so it survives the dark
# end of the speed ramp - slow ice is exactly where these profiles sit -
# without the casing widening into something that reads as a white line
PROFILE_FX = [pe.withStroke(linewidth=PROFILE_LW + 1.4, foreground='white')]


def locator_box(ax, zext, proj, color='white'):
    """Outline the zoom window on the main map.

    Dashed, where the inset's own frame is solid: at slide size the two
    rectangles are otherwise the same white box and it is not obvious
    which one is the magnified view.
    """
    ax.plot([zext[0], zext[1], zext[1], zext[0], zext[0]],
            [zext[2], zext[2], zext[3], zext[3], zext[2]],
            '--', color=color, lw=1.0, dashes=(4, 2), transform=proj,
            zorder=9)


def cumdist_km(lats, lons):
    R = 6371.0
    p = np.radians(np.asarray(lats))
    dl = np.radians(np.diff(np.asarray(lons)))
    dp = np.diff(p)
    a = np.sin(dp / 2)**2 + np.cos(p[:-1]) * np.cos(p[1:]) * np.sin(dl / 2)**2
    return np.concatenate([[0], np.cumsum(2 * R * np.arcsin(np.sqrt(a)))])


def scale_bar_br(ax, extent, frac=0.25):
    """Plain black scale bar in the BOTTOM-RIGHT corner.

    No white casing on the bar and no white box behind the label: every
    site places the bar over the pale end of its background, so the
    casing only added visual weight.
    """
    width = extent[1] - extent[0]
    height = extent[3] - extent[2]
    length = ab._nice_length(frac * width)
    x1 = extent[1] - 0.06 * width
    x0 = x1 - length
    y0 = extent[2] + 0.045 * height
    proj = ab.proj3031()
    ax.plot([x0, x1], [y0, y0], color='black', lw=3, transform=proj,
            zorder=9, solid_capstyle='butt')
    label = f'{length/1000:g} km' if length >= 1000 else f'{length:g} m'
    ax.text((x0 + x1) / 2, y0 + 0.02 * height, label, ha='center',
            va='bottom', fontsize=9, color='black', transform=proj,
            zorder=9)


def draw_map_ends(ax, lats, lons):
    """Red start / white end markers on the profile, black outlined."""
    ax.scatter([lons[0], lons[-1]], [lats[0], lats[-1]], s=END_SIZE,
               c=list(END_FACES), edgecolor='black', lw=END_EDGE,
               transform=ccrs.PlateCarree(), zorder=10)


def draw_ifg_ends(ax, dist):
    """The same two markers at the interferogram's left and right edges."""
    tr = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.scatter([dist[0], dist[-1]], [1.035, 1.035], s=END_SIZE,
               c=list(END_FACES), edgecolor='black', lw=END_EDGE,
               transform=tr, clip_on=False, zorder=10)


def hsv_plot_coherence(phase, coh, coherence_limits=(0.0, 1.0)):
    """Python port of OPR display/hsv_plot_coherence.m (John Paden):
    hue = normalized phase (-pi..pi -> 0..1), sat = 1, val = coherence
    clipped to coherence_limits and scaled 0..1. Low-coherence pixels
    fade to black instead of shouting random fringe colours."""
    from matplotlib.colors import hsv_to_rgb
    hue = phase / (2 * np.pi) + 0.5
    val = np.clip(np.abs(coh), *coherence_limits)
    val = (val - coherence_limits[0]) / (
        coherence_limits[1] - coherence_limits[0])
    val = np.nan_to_num(val, nan=0.0)
    hsv = np.stack([np.clip(hue, 0, 1), np.ones_like(hue), val], axis=-1)
    return hsv_to_rgb(hsv)


def ifg_panel(ax, fig, dist, t_us, phase, coh):
    st = max(1, phase.shape[0] // 2200)
    sx = max(1, phase.shape[1] // 2000)
    rgb = hsv_plot_coherence(phase[::st, ::sx], coh[::st, ::sx])
    ax.imshow(rgb, aspect='auto', interpolation='nearest',
              extent=[dist[0], dist[-1], t_us[-1], t_us[0]])
    ax.set_xlabel('distance along profile (km)')
    ax.set_ylabel('TWTT (μs)')
    sm = plt.cm.ScalarMappable(cmap='hsv',
                               norm=plt.Normalize(-np.pi, np.pi))
    # Clear of the image frame: at pad=0.01 the axes spine and the
    # colorbar edge merge into one line on a projected slide
    cb = fig.colorbar(sm, ax=ax, pad=0.035)
    cb.set_label('interferogram phase (rad); brightness = coherence')
    cb.set_ticks([-np.pi, 0, np.pi])
    cb.set_ticklabels([r'$-\pi$', '0', r'$\pi$'])


def two_panel(map_fn, ifg_fn, out):
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI, layout='constrained')
    gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.18])
    axm = fig.add_subplot(gs[0, 0], projection=ab.proj3031())
    axi = fig.add_subplot(gs[0, 1])
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
                                'Latitude', 'Longitude'])
    lats, lons = d['Latitude'].ravel(), d['Longitude'].ravel()
    dist = cumdist_km(lats, lons)

    tracks = []
    for fn in sorted(glob.glob(
            f'{BATCH}/joint/2023_Antarctica_Ground/*/Data_*.mat')):
        g = loadmat(fn, variable_names=['Latitude', 'Longitude'])
        la, lo = g['Latitude'][0, :], g['Longitude'][0, :]
        keep = (np.isfinite(la) & np.isfinite(lo) &
                (lo > -110) & (lo < -98) & (la > -78))
        if keep.any():
            tracks.append((la[keep], lo[keep]))
    tla = np.concatenate([t[0] for t in tracks])
    tlo = np.concatenate([t[1] for t in tracks])

    z = np.load(os.path.join(DATA, 'itslive_thwaites_win.npz'))
    wx, wy, v = z['x'], z['y'], z['v']
    vx, vy = z['vx'], z['vy']

    proj = ab.proj3031()
    extent = ab.points_extent(np.concatenate([tla, lats]),
                              np.concatenate([tlo, lons]),
                              pad_frac=0.45, min_pad_m=15e3)

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
                     label='ITS_LIVE surface speed (m/yr, log scale)')
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
        scale_bar_br(ax, extent)
        ax.set_title('Thwaites Glacier survey', fontsize=12)

        # Zoom inset: at survey scale the 10 km profile is a blob, and
        # the end markers exist to orient the interferogram against the
        # flow, so they need a window where the two ends separate
        px, py = proj.transform_point(np.mean(lons), np.mean(lats),
                                      ccrs.PlateCarree())[:2]
        half = 9e3
        zext = (px - half, px + half, py - half, py + half)
        axz = ax.inset_axes([0.02, 0.30, 0.55, 0.34], projection=proj)
        axz.set_extent(zext, crs=proj)
        speed_layer(axz, dec=18)
        axz.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
                 transform=ccrs.PlateCarree(), zorder=8,
                 path_effects=PROFILE_FX)
        draw_map_ends(axz, lats, lons)
        for spine in axz.spines.values():
            spine.set_edgecolor('white')
            spine.set_linewidth(1.5)
        locator_box(ax, zext, proj)

    def ifg(fig, ax):
        ifg_panel(ax, fig, dist, d['Time'].ravel() * 1e6,
                  d['phase_wrapped'], np.abs(d['coherence']))
        draw_ifg_ends(ax, dist)
        ax.set_title(f'eastern shear margin crossing, segment {seg}',
                     fontsize=12, pad=26)

    two_panel(map_panel, ifg, 'scar_thwaites_2panel.png')


# ---------------------------------------------------------------- Ridge A
def ridge_a():
    frame = '20250108_02_009'
    seg = frame.rsplit('_', 1)[0]
    d = loadmat(f'{ACCUM}/CSARP_polarimetric_unwrap/{seg}/Data_{frame}.mat',
                variable_names=['interferogram_mlook',
                                'interferogram_coherence', 'Time',
                                'Latitude', 'Longitude'])
    lats, lons = d['Latitude'].ravel(), d['Longitude'].ravel()
    dist = cumdist_km(lats, lons)

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
        axz = ax.inset_axes([0.02, 0.30, 0.5, 0.32], projection=proj)
        axz.set_extent(zext, crs=proj)
        axz.set_facecolor('0.28')
        for la, lo in tracks:
            axz.plot(lo, la, '-', color=TRACK_C, lw=1.6, alpha=0.9,
                     transform=ccrs.PlateCarree(), zorder=5)
        axz.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
                 transform=ccrs.PlateCarree(), zorder=8,
                 path_effects=PROFILE_FX)
        draw_map_ends(axz, lats, lons)
        for spine in axz.spines.values():
            spine.set_edgecolor('white')
            spine.set_linewidth(1.5)
        locator_box(ax, zext, proj)
        gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5,
                          color='gray', y_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        scale_bar_br(ax, extent)
        ax.annotate('no satellite velocity coverage south of 82.7°S\n'
                    '(interior divide site, flow < 2 m/yr)',
                    xy=(0.03, 0.03), xycoords='axes fraction', fontsize=7.5,
                    color='white')
        ax.set_title('Ridge A raster survey', fontsize=12)

    def ifg(fig, ax):
        ifg_panel(ax, fig, dist, d['Time'].ravel() * 1e6,
                  np.angle(d['interferogram_mlook']),
                  np.abs(d['interferogram_coherence']))
        draw_ifg_ends(ax, dist)
        ax.set_title(f'divide setting, frame {frame}', fontsize=12, pad=26)

    two_panel(map_panel, ifg, 'scar_ridge_a_2panel.png')


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    thwaites()
    ridge_a()

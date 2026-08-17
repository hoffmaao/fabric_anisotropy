"""Ridge A wrapped-phase movie: one segment at a time across the survey.

Steps through every frame of the Ridge A raster in acquisition order. Each
movie frame is the scar_two_panel.py slide geometry, unchanged:

  left   the whole survey in white over the REMA surface, with the CURRENT
         segment drawn as the black profile and its red/white end markers
  right  that segment's wrapped polarimetric interferogram

so the marker walks the grid while the fringes beside it change. The map
never moves, which is the point: every segment is read against the same
survey outline rather than against a re-framed view.

NO PANEL TITLES. The still figure titles each panel ("Ridge A raster
survey", "divide setting, frame ..."), but in a movie the right-hand title
would rewrite itself 36 times and the left-hand one never changes, so both
are dropped and the moving marker carries the identification.

The zoom inset of the still figure is also dropped. It exists there to
separate the two ends of one profile; here it would re-frame on every
segment and jitter through the whole movie.

Reads the compact per-frame npz written by reduce_ifg_movie.py on mem1 -
the shipped products are ~690 MB a frame, ~25 GB for the survey, so the
phase and coherence arrive already decimated to display resolution and
quantized to uint8. Geometry stays float; see that script for the trade.

Usage: python ridge_a_ifg_movie.py <ifg_movie_dir> <out_dir> [fps] [sec]
"""
import glob
import os
import subprocess
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import cartopy.crs as ccrs                 # noqa: E402
import matplotlib.pyplot as plt            # noqa: E402
from scipy.io import loadmat               # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab             # noqa: E402
import scar_style as sty                   # noqa: E402
from scar_style import (DPI, INK, PROFILE_C, PROFILE_FX,  # noqa: E402
                        PROFILE_LW, TRACK_C)

_argv = list(sys.argv[1:])
# --still <tag> renders that ONE segment as a PNG and stops. It exists so
# the slide that precedes or follows the movie can be built by this same
# code path rather than by matching scar_two_panel.py's geometry to it by
# hand: identical figure size, identical panel split, identical map extent
# and identical absence of titles, so cutting between the slide and the
# movie moves nothing on screen. Any hand-matched copy would drift the
# moment either side is touched.
STILL = None
if '--still' in _argv:
    i = _argv.index('--still')
    STILL = _argv[i + 1]
    del _argv[i:i + 2]

IFG_DIR = _argv[0] if len(_argv) > 0 else 'ifg_movie'
OUT = _argv[1] if len(_argv) > 1 else '.'
FPS = int(_argv[2]) if len(_argv) > 2 else 25
SEC_PER_SEG = float(_argv[3]) if len(_argv) > 3 else 1.07
BATCH = os.path.expanduser('~/data/opr/fabric_batch')
BASEMAP = os.path.expanduser(os.environ.get(
    'ANT_BASEMAP_DIR', '~/data/opr/basemap'))
# REMA v2.0 32 m tile 32_34 covers the Ridge A grid (x 299904..400096,
# y 99904..200096 in EPSG:3031). Fetch once with
#   base=https://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v2.0/32m
#   curl -O $base/32_34/32_34_32m_v2.0.tar.gz
# and keep only the _dem.tif. The whole-continent 100 m mosaic is 5.6 GB,
# and the 1 km one puts only ~25 pixels across this survey, so the single
# 32 m tile is the right unit here.
REMA_FN = os.path.join(BASEMAP, 'REMA_32m_32_34_dem.tif')


def find_ffmpeg():
    """ffmpeg is not on PATH in the cartopy env these figures render in."""
    import shutil
    for c in (shutil.which('ffmpeg'),
              '/opt/anaconda3/envs/ar-env/bin/ffmpeg',
              '/opt/homebrew/bin/ffmpeg', '/usr/local/bin/ffmpeg'):
        if c and os.path.exists(c):
            return c
    raise SystemExit('no ffmpeg found; frames are in <out>/frames if you '
                     'want to encode them by hand')


def read_rema(extent):
    """REMA surface elevation over `extent`, or None if unavailable.

    Degrades the way antarctic_basemap.add_imagery does: a missing tile or
    a missing rasterio leaves the caller with its plain background rather
    than failing the render.
    """
    try:
        import rasterio
        from rasterio.windows import from_bounds
    except Exception:
        print('  rasterio unavailable; keeping the plain background')
        return None
    if not os.path.exists(REMA_FN):
        print('  no REMA tile at %s; keeping the plain background' % REMA_FN)
        return None
    with rasterio.open(REMA_FN) as src:
        w = from_bounds(extent[0], extent[2], extent[1], extent[3],
                        src.transform)
        a = src.read(1, window=w, masked=True).astype(float)
        a = np.ma.masked_invalid(np.ma.filled(a, np.nan))
        # REMA marks voids with a large negative sentinel, which would
        # stretch the colour scale over an elevation that does not exist
        a = np.ma.masked_less(a, -1000.0)
    print('  REMA %s, %.0f-%.0f m' % (a.shape, a.min(), a.max()))
    return a


def load_tracks():
    """Every Ridge A survey leg, rebuilt per segment as in scar_two_panel."""
    by_seg = {}
    for fn in sorted(glob.glob(
            f'{BATCH}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat')):
        g = loadmat(fn, variable_names=['Latitude', 'Longitude'])
        la = np.atleast_1d(g['Latitude'].ravel())
        lo = np.atleast_1d(g['Longitude'].ravel())
        keep = np.isfinite(la) & np.isfinite(lo)
        if keep.any():
            by_seg.setdefault(os.path.basename(os.path.dirname(fn)),
                              []).append((la[keep], lo[keep]))
    if not by_seg:
        raise SystemExit(
            'no track products under %s/joint_jp/2024_Antarctica_Ground2; '
            'mirror the fabric batch before building the movie' % BATCH)
    return [(np.concatenate([p[0] for p in v]),
             np.concatenate([p[1] for p in v]))
            for v in by_seg.values()]


def load_frames(d):
    fns = sorted(glob.glob(os.path.join(d, 'ifg_*.npz')))
    if not fns:
        raise SystemExit(
            'no ifg_*.npz under %s; run reduce_ifg_movie.py on mem1 and '
            'fetch its output there' % d)
    out = []
    for fn in fns:
        z = np.load(fn)
        out.append(dict(
            tag=os.path.basename(fn)[len('ifg_'):-len('.npz')],
            # uint8 -> physical units, the inverse of the reduction
            phase=(z['phase'].astype(np.float32) / 255.0) * 2 * np.pi - np.pi,
            coh=z['coh'].astype(np.float32) / 255.0,
            t_us=z['t'] * 1e6, dist=z['dist'],
            lat=z['lat_full'], lon=z['lon_full']))
    return out


def main():
    tracks = load_tracks()
    frames = load_frames(IFG_DIR)
    print('%d segments, %d survey legs' % (len(frames), len(tracks)))

    # The still is one segment of the SAME sequence, so the map extent and
    # the shared TWTT window are computed over ALL segments even when only
    # one is drawn. Deriving them from the single segment instead would
    # reframe the map and rescale the depth axis, which is precisely the
    # apparent movement the still exists to avoid.
    still_tag = None
    if STILL is not None:
        want = [f for f in frames if f['tag'] == STILL]
        if not want:
            raise SystemExit('no segment %s among %d; tags look like %s'
                             % (STILL, len(frames), frames[0]['tag']))
        still_tag = want[0]['tag']

    proj = ab.proj3031()
    gla = np.concatenate([t[0] for t in tracks]
                         + [f['lat'] for f in frames])
    glo = np.concatenate([t[1] for t in tracks]
                         + [f['lon'] for f in frames])
    extent = ab.points_extent(gla, glo, pad_frac=0.10, min_pad_m=3e3)

    # One fixed TWTT window for every segment. Letting each autoscale would
    # make the depth axis breathe from frame to frame, which reads as the
    # ice changing thickness rather than as the record changing length.
    t_lo = min(float(f['t_us'][0]) for f in frames)
    t_hi = max(float(f['t_us'][-1]) for f in frames)
    print('common TWTT window %.1f-%.1f us' % (t_lo, t_hi))

    rema = read_rema(extent)
    if rema is not None:
        # fixed stretch across the whole sweep so the background never
        # changes between frames and cannot read as the surface moving
        rema_lo = float(np.percentile(rema.compressed(), 1))
        rema_hi = float(np.percentile(rema.compressed(), 99))
    else:
        rema_lo = rema_hi = None

    fdir = os.path.join(OUT, 'frames')
    os.makedirs(fdir, exist_ok=True)
    hold = max(1, int(round(FPS * SEC_PER_SEG)))
    n = 0
    draw = [f for f in frames if still_tag is None or f['tag'] == still_tag]
    for k, f in enumerate(draw):
        fig, axm, axi = sty.two_panel_figure(proj)

        axm.set_extent(extent, crs=proj)
        axm.set_facecolor('0.35')
        if rema is not None:
            # Dark, perceptually uniform and colourblind-safe, so the white
            # tracks and the black profile both keep their contrast against
            # it - a light or high-chroma terrain ramp would swallow one or
            # the other. Elevation is context here, not the subject.
            # alpha over the grey facecolor pulls the whole ramp toward
            # mid-tone: at full strength cividis runs from near-black to a
            # bright yellow, and the white tracks lose contrast at the top
            # end while the black profile loses it at the bottom.
            im = axm.imshow(rema, extent=extent, origin='upper',
                            cmap='cividis', vmin=rema_lo, vmax=rema_hi,
                            zorder=0, alpha=0.72, interpolation='bilinear',
                            transform=proj)
            # Standard colorbar BELOW the map panel rather than an inset on
            # it: on the map the caption had to be white to be legible and
            # still collided with the panel edge. Below it, centred under
            # the panel it describes, it can be ordinary black text.
            cb = fig.colorbar(im, ax=axm, orientation='horizontal',
                              location='bottom', shrink=0.62, pad=0.045,
                              aspect=32)
            cb.set_label('REMA surface elevation (m)', fontsize=9,
                         color=INK)
            cb.ax.tick_params(labelsize=8, colors=INK)
        for la, lo in tracks:
            axm.plot(lo, la, '-', color=TRACK_C, lw=1.1, alpha=0.9,
                     transform=ccrs.PlateCarree(), zorder=5)
        axm.plot(f['lon'], f['lat'], '-', color=PROFILE_C, lw=PROFILE_LW,
                 transform=ccrs.PlateCarree(), zorder=8,
                 path_effects=PROFILE_FX)
        sty.draw_map_ends(axm, f['lat'], f['lon'])
        gl = axm.gridlines(draw_labels=True, lw=0.4, alpha=0.5,
                           color='gray', y_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        sty.continental_inset(axm, 'antarctica',
                              (0.02, 0.755, 0.26, 0.225), extent, proj)
        ab.scale_bar_br(axm, extent, proj)

        sty.ifg_panel(axi, fig, f['dist'], f['t_us'], f['phase'], f['coh'])
        # The end markers sit at axes fraction 1.035 with clip_on=False, so
        # they live in whatever headroom the layout reserves ABOVE the axes.
        # In the still that headroom comes from the panel title's pad=26;
        # dropping the title took the markers off the canvas with it. A
        # blank title restores exactly that reserve without putting the
        # words back, so the markers land where they do in the still.
        axi.set_title(' ', fontsize=12, pad=26)
        sty.draw_ifg_ends(axi, f['dist'])
        axi.set_ylim(t_hi, t_lo)

        if still_tag is not None:
            os.makedirs(OUT, exist_ok=True)
            still_fn = os.path.join(
                OUT, 'ridge_a_2panel_still_%s.png' % still_tag)
            fig.savefig(still_fn, dpi=DPI, facecolor='white')
            plt.close(fig)
            print('wrote %s (movie-frame geometry, no titles)' % still_fn)
            return
        png = os.path.join(fdir, 'f%04d.png' % k)
        fig.savefig(png, dpi=DPI, facecolor='white')
        plt.close(fig)
        # hold each segment for SEC_PER_SEG by repeating its rendered frame
        # rather than re-rendering it, which is ~40x cheaper
        for r in range(1, hold):
            os.link(png, os.path.join(fdir, 'f%04d_%03d.png' % (k, r)))
        n += hold
        print('  [%2d/%d] %s' % (k + 1, len(draw), f['tag']), flush=True)

    seq = os.path.join(fdir, 'seq_%05d.png')
    for i, p in enumerate(sorted(glob.glob(os.path.join(fdir, 'f*.png')))):
        os.rename(p, seq % i)
    mp4 = os.path.join(OUT, 'ridge_a_ifg_movie.mp4')
    subprocess.run(
        [find_ffmpeg(), '-y', '-framerate', str(FPS), '-i', seq,
         '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20',
         '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', mp4],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print('wrote %s (%d frames, %.1f s)' % (mp4, n, n / FPS))


if __name__ == '__main__':
    main()

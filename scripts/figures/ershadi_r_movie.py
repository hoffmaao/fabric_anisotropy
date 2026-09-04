"""Depth-sweep movie of the Ershadi reflection ratio r over a site's map.

Companion to quadpol_depth_movie.py, for the OTHER Fujita parameter: each
frame's retrieved r_dB profile (stages/ershadi_r/ershadi_r_<tag>.mat, built
by opr_fabric/server/run_ershadi_r.m) painted onto the frame's section
blocks as the sweep descends. r is retrieved PER FRAME per 50 m interval -
the science it feeds is the heading-family comparison, which is per frame -
so every block of a frame carries the frame's value at the current depth.

Colour is DIVERGING because r is signed: r_dB > 0 means the y-axis
scattering coefficient exceeds the x-axis one and the sign flips with the
stronger axis, so a sequential ramp would hide the physically meaningful
zero crossing (Ershadi's EDML has zones of both signs). RdBu, centred on
0 dB, ceiling VMAX_DB; blocks whose frame abstained at this depth (any r13
gate, or no interval) are drawn in grey - the same drawn-but-unclaimed
convention as the dlam movie.

The key sits in its OWN STRIP BELOW the map, never on it - the rule from
the dlam movie's occlusion fix: no spot inside the axes is safe from a
track at some site.

Usage: python ershadi_r_movie.py [<out_dir>] [<site>]
Reads $SCAR_DATA (default ~/data/opr/scar): quadpol_section_<tag>.mat for
block positions, ershadi_r/ershadi_r_<tag>.mat for the profiles.
"""
import glob
import os
import re
import shutil
import subprocess
import sys

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab           # noqa: E402
import quadpol_sites as qs               # noqa: E402
from scar_style import FIGS              # noqa: E402

DATA = os.path.expanduser(os.environ.get('SCAR_DATA', '~/data/opr/scar'))
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGS
SITE = sys.argv[2] if len(sys.argv) > 2 else 'ridge_a'

CFG = qs.get(SITE)
VMAX_DB = 15.0
Z0, Z1, ZSTEP = 225.0, 1375.0, 25.0
FPS, HOLD = 25, 6
GREY = (0.62, 0.62, 0.62)


def find_ffmpeg():
    for c in (shutil.which('ffmpeg'), '/opt/anaconda3/envs/ar-env/bin/ffmpeg',
              '/opt/homebrew/bin/ffmpeg', '/usr/local/bin/ffmpeg'):
        if c and os.path.exists(c):
            return c
    raise SystemExit('no ffmpeg found; frames are in <out>/r_frames_%s' % SITE)


def load():
    frames = []
    for fn in sorted(glob.glob(os.path.join(DATA, 'quadpol_section_*.mat'))):
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        if re.search(r'_z\d+$', tag):
            continue
        rfn = os.path.join(DATA, 'ershadi_r', 'ershadi_r_%s.mat' % tag)
        if not os.path.exists(rfn):
            continue
        with h5py.File(fn) as f:
            r = f['res']
            la = np.array(r['sec_lat']).ravel()
            lo = np.array(r['sec_lon']).ravel()
            lat = float(np.array(r['lat']).ravel()[0])
            lon = float(np.array(r['lon']).ravel()[0])
        # membership is decided on the frame's own centre, as every other
        # figure does it - a site may be pinned by segment, by season or by
        # a position box, and only in_site knows which
        if not qs.in_site(CFG, lat, lon, tag):
            continue
        with h5py.File(rfn) as f:
            rr = f['res']
            edges = np.array(rr['edges']).ravel()
            r_db = np.array(rr['r_db_int']).ravel()
        ok = np.isfinite(la) & np.isfinite(lo)
        if not ok.any():
            continue
        x, y = T.transform(lo[ok], la[ok])
        frames.append(dict(tag=tag, pts=np.column_stack([x, y]),
                           edges=edges, r_db=r_db))
    return frames


def r_at(fr, zc):
    e = fr['edges']
    k = np.searchsorted(e, zc, side='right') - 1
    if k < 0 or k >= fr['r_db'].size:
        return np.nan
    return fr['r_db'][k]


from pyproj import Transformer           # noqa: E402
T = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)


def main():
    frames = load()
    if not frames:
        raise SystemExit('no ershadi_r profiles for %s under %s'
                         % (SITE, DATA))
    centres = np.arange(Z0, Z1 + 1e-6, ZSTEP)
    print('%d frames with r profiles, %d depth steps' %
          (len(frames), centres.size))

    proj = ab.proj3031()
    cmap = plt.get_cmap('RdBu_r')
    norm = plt.Normalize(-VMAX_DB, VMAX_DB)
    allx = np.concatenate([f['pts'][:, 0] for f in frames])
    ally = np.concatenate([f['pts'][:, 1] for f in frames])
    sq, span_x, span_y = qs.square_extent(allx, ally)

    fdir = os.path.join(OUT, 'r_frames_%s' % SITE)
    if os.path.isdir(fdir):
        shutil.rmtree(fdir)
    os.makedirs(fdir)

    for k, zc in enumerate(centres):
        fig = plt.figure(figsize=(6.4, 7.2))
        gs = fig.add_gridspec(2, 1, height_ratios=[6.0, 0.9], hspace=0.06)
        ax = fig.add_subplot(gs[0], projection=proj)
        ax.set_extent(sq, crs=proj)
        ab.add_imagery(ax, sq)
        for fr in frames:
            v = r_at(fr, zc)
            c = GREY if not np.isfinite(v) else cmap(norm(v))
            ax.scatter(fr['pts'][:, 0], fr['pts'][:, 1], s=7, color=c,
                       transform=proj, linewidths=0)
        ab.scale_bar(ax, sq)
        ax.set_title('%s: reflection ratio r, %.0f m' % (CFG['title'], zc),
                     fontsize=11)
        # key strip below the map, never on it
        cax = fig.add_subplot(gs[1])
        cax.set_axis_off()
        cb_ax = cax.inset_axes([0.06, 0.55, 0.6, 0.3])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cb_ax,
                     orientation='horizontal')
        cb_ax.set_title('r  [dB, 20log10 Gy/Gx]', fontsize=8, pad=2)
        cax.add_patch(plt.Rectangle((0.74, 0.55), 0.05, 0.3, color=GREY,
                                    transform=cax.transAxes))
        cax.text(0.805, 0.70, 'abstained', fontsize=8, va='center',
                 transform=cax.transAxes)
        fig.savefig(os.path.join(fdir, 'd%05d.png' % k), dpi=110)
        plt.close(fig)
        if (k + 1) % 10 == 0 or k == centres.size - 1:
            print('  [%d/%d] %.0f m' % (k + 1, centres.size, zc))

    seq = os.path.join(fdir, 'seq_%05d.png')
    pngs = sorted(glob.glob(os.path.join(fdir, 'd*.png')))
    n = 0
    for i, p in enumerate(pngs):
        reps = HOLD if (i == 0 or i == len(pngs) - 1) else 1
        for _ in range(reps):
            shutil.copy(p, seq % n)
            n += 1
    out_fn = os.path.join(OUT, 'ershadi_r_movie_%s.mp4' % SITE)
    subprocess.run([find_ffmpeg(), '-y', '-loglevel', 'error', '-framerate',
                    str(FPS), '-i', seq, '-pix_fmt', 'yuv420p', out_fn],
                   check=True)
    print('wrote %s (%d frames)' % (out_fn, n))


if __name__ == '__main__':
    main()

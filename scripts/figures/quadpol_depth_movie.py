"""Fabric strength through the column: one map per depth window, animated.

The map figures fix a deep band (1150-1500 m) and show how fabric varies
laterally. This sweeps that window DOWN the column instead, so the survey
map is redrawn every 20 m and the movie shows how the fabric strengthens
with depth across the whole grid at once.

Everything except the depth window is held fixed, deliberately: the same
colour scale on every frame, the same extent, the same graticule. A frame
that autoscaled its colours would make every depth look equally
anisotropic, which is the one thing this movie exists to disprove.

The depth axis is the LS estimator's own sec_dlam_ls, per ~125 m block
(a length, so the trace count behind it differs with the site's spacing),
so a block enters a frame only where that window holds enough finite cells
to support a median; blocks that fail leave the grey track showing rather
than being interpolated.

ON THIS SURVEY NOTHING FAILS: coverage is 2157 of 2157 blocks in EVERY
window from 200 m to 1450 m. So the pale shallow frames are MEASUREMENTS
OF A WEAK FABRIC, not gaps, and must not be read as missing data - which
is the opposite of the gap-banded behaviour the published direct chain
shows over the same frames, and one of the things the LS fit bought. The
measured profile, median over all blocks:

    200-350 m  0.013     600-750 m  0.062     1200-1350 m  0.064
    300-450 m  0.016     800-950 m  0.055     1300-1450 m  0.073
    400-550 m  0.040    1000-1150 m 0.064

so the contrast rises about fivefold down the column, with most of the
climb between 350 and 550 m and a mild reversal around 800-950 m. Note the
shallow half is where the co-polarized comparison is still unresolved at
about a factor of two, so treat the top few frames as the LS estimator's
own answer rather than as a cross-validated one.

Usage: python quadpol_depth_movie.py <out_dir> [site]
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.patheffects as mpe                   # noqa: E402
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.collections import LineCollection      # noqa: E402
from pyproj import Transformer                         # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import antarctic_basemap as ab                         # noqa: E402
import quadpol_sites as qs                             # noqa: E402
from scar_style import DATA, INK, MUTED                # noqa: E402
from scar_style import FIGS  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else FIGS
SITE = sys.argv[2] if len(sys.argv) > 2 else 'ridge_a'

FPS, HOLD = 25, 6
BAR_FRAC = 0.035       # orientation bar half-length as a fraction of span

CFG = qs.get(SITE)
TITLE = CFG['title']
VMAX = CFG['vmax']
# Ridge A stops at 1050 m, not at the bottom of the record: below that its
# mostly N-S lines stop following the regular fabric cycle and the deeper
# windows would read as fabric development when they are something else.
# Every other site runs to the depth its own coverage supports.
Z_FIRST, Z_LAST = CFG['movie']
# window and step come from the site: see quadpol_sites for why a single
# window height cannot serve both a 1850 m column and a 300 m ice shelf
WIN_M = CFG.get('win', 150.0)
STEP_M = CFG.get('step', 20.0)
MIN_CELLS = 3 if WIN_M < 80.0 else 5
# Fit-quality fade. A window at RESID_GOOD or better draws at full colour; by
# RESID_BAD it has gone to grey. The bounds are the LS residual, measured
# across these surveys: Ridge A sits at 0.20-0.26 through its whole column
# and its displayed values are stable under gating, so 0.25 is what a
# healthy fit looks like here and 0.45 is where the value stops meaning
# anything.
RESID_GOOD, RESID_BAD = 0.25, 0.45
GREY_RGB = (0.88, 0.88, 0.89)
# Halo for the annotations that must stay ON the map, where a track may run
# under them at one site and not at the next.
HALO = [mpe.withStroke(linewidth=3.0, foreground='white')]

# Bed depth per along-track block, from the CSARP_layer bottom pick, keyed by
# frame tag. Without it a frame keeps drawing fabric colour at depths where
# its ice has already ended - and the ice varies enormously WITHIN a survey,
# 566-996 m across Taylor Dome and by 200 m inside single frames, so a
# per-site depth cut cannot express it. A block DROPS OUT of the movie for
# good once the bed enters its window - not proportionally, see the rule at
# the draw call - so a frame retires block by block as its own bed comes up.
# Blocks with no pick are drawn unchanged.
#
# Built by scripts/make_bed_by_block.py; the format is documented there and
# in docs/scripts.md. A MISSING file is legitimate - it means no picks were
# staged for this survey - and only says so. A file that is present but
# unreadable is not: it would silently render the pre-bed-masking movie, so
# it raises. A file that is present but does not name this frame warns per
# frame: that frame draws unmasked, which is indistinguishable from a survey
# with no bed to mask unless it is said out loud.
BED_FILE = os.path.join(DATA, 'bed_by_block.json')
if os.path.exists(BED_FILE):
    with open(BED_FILE) as _fh:
        BEDS = json.load(_fh)
    if not isinstance(BEDS, dict):
        raise SystemExit('%s must hold an object keyed by frame tag, got %s'
                         % (BED_FILE, type(BEDS).__name__))
else:
    BEDS = {}
    print('no %s; blocks are drawn without bed masking' % BED_FILE)

T = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)


def find_ffmpeg():
    for c in (shutil.which('ffmpeg'), '/opt/anaconda3/envs/ar-env/bin/ffmpeg',
              '/opt/homebrew/bin/ffmpeg', '/usr/local/bin/ffmpeg'):
        if c and os.path.exists(c):
            return c
    raise SystemExit('no ffmpeg found; frames are in '
                     '<out>/depth_frames_%s' % SITE)


def load():
    frames = []
    for fn in sorted(glob.glob(os.path.join(
            DATA, 'quadpol_section_*.mat'))):
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        if re.search(r'_z\d+$', tag):
            continue
        with h5py.File(fn) as f:
            r = f['res']
            z = np.array(r['z']).ravel()
            secl = np.array(r['sec_dlam_ls']).T
            resl = np.array(r['sec_resid_ls']).T
            blat = np.array(r['sec_lat']).ravel()
            blon = np.array(r['sec_lon']).ravel()
            zw = np.array(r['ls_zw']).ravel()
            th_geo = np.array(r['ls_theta0_geo']).ravel()
            lat = float(np.array(r['lat']).ravel()[0])
            lon = float(np.array(r['lon']).ravel()[0])
            lat0 = float(np.array(r['lat0']).ravel()[0])
            lon0 = float(np.array(r['lon0']).ravel()[0])
            lat1 = float(np.array(r['lat1']).ravel()[0])
            lon1 = float(np.array(r['lon1']).ravel()[0])
        if not qs.in_site(CFG, lat, lon, tag):
            continue
        ok = np.isfinite(blat) & np.isfinite(blon)
        if not ok.any():
            continue
        bx, by = T.transform(blon[ok], blat[ok])
        ex, ey = T.transform([lon0, lon1], [lat0, lat1])
        # bed depth per block, filtered by the same mask as the positions so
        # index j means the same block in both
        bed = BEDS.get(tag)
        if bed is None:
            # a bed file that does not name this frame is NOT the same as no
            # bed file at all: it means the run that wrote it did not cover
            # this survey, and an unnamed frame draws unmasked, which looks
            # exactly like a survey with no bed to mask
            if BEDS:
                print('WARNING: %s carries no bed for %s; its blocks are '
                      'drawn unmasked - rerun scripts/make_bed_by_block.py '
                      "with this survey's site root"
                      % (os.path.basename(BED_FILE), tag))
            bedv = np.full(int(ok.sum()), np.nan)
        elif len(bed) != ok.size:
            # a stale bed file cut against a different block size masks the
            # wrong ice; say so rather than quietly dropping to no masking
            print('WARNING: %s lists %d bed values for %d blocks in %s; '
                  'bed masking skipped for this frame - rebuild it with '
                  'scripts/make_bed_by_block.py'
                  % (os.path.basename(BED_FILE), len(bed), ok.size, tag))
            bedv = np.full(int(ok.sum()), np.nan)
        else:
            bedv = np.array([np.nan if b is None else b for b in bed],
                            dtype=float)[ok]
        frames.append(dict(tag=tag, z=z, sec=secl[:, ok], res=resl[:, ok],
                           bed=bedv, pts=np.column_stack([bx, by]),
                           zw=zw, th_geo=th_geo, lat=lat, lon=lon,
                           ex=ex, ey=ey))
    if not frames:
        raise SystemExit('no sections within %.0f deg of %s under %s'
                         % (qs.SITE_RADIUS_DEG, SITE, DATA))
    return frames


def main():
    frames = load()
    centres = np.arange(Z_FIRST + WIN_M / 2, Z_LAST - WIN_M / 2 + 1e-6,
                        STEP_M)
    print('%d frames, %d depth windows (%.0f m window, %.0f m step)'
          % (len(frames), centres.size, WIN_M, STEP_M))

    proj = ab.proj3031()
    cmap, norm = plt.get_cmap('Blues'), plt.Normalize(0.0, VMAX)
    allx = np.concatenate([f['pts'][:, 0] for f in frames])
    ally = np.concatenate([f['pts'][:, 1] for f in frames])
    sq, span_x, span_y = qs.square_extent(allx, ally)
    bar_km = BAR_FRAC * max(span_x, span_y) / 1000.0
    print('  survey span %.1f x %.1f km, bar %.2f km'
          % (span_x / 1e3, span_y / 1e3, bar_km))

    # per-site frame dir, cleared up front: leftover seq_/d* frames from a
    # crashed run (or another site) would be numbered after this run's and
    # ffmpeg would append them to the tail of the movie
    fdir = os.path.join(OUT, 'depth_frames_%s' % SITE)
    if os.path.isdir(fdir):
        shutil.rmtree(fdir)
    os.makedirs(fdir)
    n = 0
    for k, zc in enumerate(centres):
        lo, hi = zc - WIN_M / 2, zc + WIN_M / 2
        # The key lives in a strip of its OWN below the map, never on it.
        # No spot inside the axes is safe to put it: square_extent pads the
        # block positions to 1.18x, which leaves the blocks in the middle
        # 84.7% and looks like a clear band underneath - but the orientation
        # bars overhang their block by BAR_FRAC of the span, reaching down to
        # axes fraction 0.047, and the frame end-lines are drawn from
        # endpoints that are not in the block set bounding the extent at all.
        # So a key position tuned until it clears one survey re-occludes on
        # the next; that is how the track came to sit over the Taylor Dome
        # colour-bar title. Outside the axes the geometry cannot follow.
        fig = plt.figure(figsize=(7.6, 8.0), layout='constrained')
        gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.075])
        ax = fig.add_subplot(gs[0, 0], projection=proj)
        lax = fig.add_subplot(gs[1, 0])
        lax.set_axis_off()

        def arm(lat, lon, th_deg, half_km=None, color=INK, lw=1.6):
            half_km = bar_km if half_km is None else half_km
            th = np.deg2rad(th_deg)
            dlat = half_km * np.cos(th) / 111.2
            dlon = half_km * np.sin(th) / (111.2 * np.cos(np.deg2rad(lat)))
            x, y = T.transform([lon - dlon, lon + dlon],
                               [lat - dlat, lat + dlat])
            ax.plot(x, y, '-', color=color, lw=lw, solid_capstyle='round',
                    zorder=4, transform=proj)

        cov = 0
        for fr in frames:
            ax.plot(fr['ex'], fr['ey'], '-', color=MUTED, lw=1.0, alpha=0.8,
                    zorder=1, transform=proj)
            mz = (fr['z'] > lo) & (fr['z'] < hi)
            sub = fr['sec'][mz]
            bd = np.where(np.isfinite(sub).sum(axis=0) >= MIN_CELLS,
                          np.nanmedian(np.where(np.isfinite(sub), sub,
                                                np.nan), axis=0), np.nan)
            # Fit quality per block-window, on the same footing as the value.
            subr = fr['res'][mz]
            br = np.where(np.isfinite(subr).sum(axis=0) >= MIN_CELLS,
                          np.nanmedian(np.where(np.isfinite(subr), subr,
                                                np.nan), axis=0), np.nan)
            cov += int(np.isfinite(bd).sum())
            # Grey the whole interval as soon as the bed enters it, rather
            # than fading in proportion to how much ice is left. The bed is a
            # strong specular reflector: a window that contains it has its
            # coherence dominated by the bed echo throughout, so the fabric
            # estimate is contaminated across the whole window, not partially.
            # A block therefore draws while the window sits entirely in ice
            # and drops out for good once the bed appears in it.
            bedf = fr['bed']
            with np.errstate(invalid='ignore'):
                icef = (hi <= bedf).astype(float)
            icef = np.where(np.isfinite(bedf), icef, 1.0)
            segs, vals, wts, ice = [], [], [], []
            for j in range(len(bd) - 1):
                v = 0.5 * (bd[j] + bd[j + 1])
                if np.isfinite(v):
                    segs.append([fr['pts'][j], fr['pts'][j + 1]])
                    vals.append(v)
                    wts.append(np.nanmean([br[j], br[j + 1]]))
                    ice.append(min(icef[j], icef[j + 1]))
            if segs:
                # Blend toward neutral grey by fit quality, the same rule the
                # co-polarized sections use. Without it a badly-fit window is
                # painted at full saturation and reads exactly like a measured
                # one - and at every site except Ridge A that is most of the
                # map: gating at the RESID_GOOD level used here, resid<0.25,
                # drops the displayed p95 from 0.119 to 0.108 at Taylor Dome
                # and from 0.181 to 0.126 at Thwaites, while Ridge A barely
                # moves (0.083 -> 0.079). Those are the same numbers the
                # ceilings in quadpol_sites.py are measured from, so what the
                # ceiling scales is what the map shows at full colour. The
                # bright tail was the fit failing, not the ice.
                rgb = cmap(norm(np.asarray(vals)))[:, :3]
                w = np.clip((RESID_BAD - np.asarray(wts, float))
                            / (RESID_BAD - RESID_GOOD), 0.0, 1.0)
                w[~np.isfinite(w)] = 0.0
                # below the bed there is no fabric to report, whatever the
                # fit residual says, so the ice fraction gates the colour
                w = w * np.asarray(ice, float)
                rgb = rgb * w[:, None] + np.array(GREY_RGB) * (1 - w[:, None])
                lc = LineCollection(segs, colors=rgb, linewidths=4.2,
                                    capstyle='round', zorder=2, transform=proj)
                ax.add_collection(lc)

        # Orientation AT THIS DEPTH, not the frame average: theta0 comes
        # from the LS windows inside this same slice, combined as a
        # doubled-angle phasor (never an unwrapped angle - an unwrap over
        # gappy axis samples can slip a branch and hand a silent 90 deg
        # error to the bar). A frame whose windows all abstained here draws
        # no bar, so the axis appears only where it is actually resolved.
        th_here = []
        for fr in frames:
            mw = (fr['zw'] > lo) & (fr['zw'] < hi) & np.isfinite(fr['th_geo'])
            if mw.sum() < 2:
                continue
            ph = np.mean(np.exp(2j * np.deg2rad(fr['th_geo'][mw])))
            if not np.isfinite(ph) or abs(ph) < 1e-9:
                continue
            th = np.rad2deg(np.angle(ph)) / 2 % 180
            arm(fr['lat'], fr['lon'], th)
            th_here.append(th)

        ax.set_xlim(sq[0], sq[1])
        ax.set_ylim(sq[2], sq[3])
        xs, ys = ax.get_xlim(), ax.get_ylim()
        gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='gray',
                          y_inline=False, x_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 8.5, 'color': INK}
        gl.ylabel_style = {'size': 8.5, 'color': INK}

        cax = lax.inset_axes([0.30, 0.52, 0.26, 0.30])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                     orientation='horizontal')
        # label alongside the bar rather than above it: the strip is shallow
        # and the map's rotated longitude labels come down to meet it
        cax.text(-0.035, 0.5, r'$\Delta\lambda$', transform=cax.transAxes,
                 ha='right', va='center', fontsize=9.5, color=INK)
        cax.tick_params(labelsize=7.5)

        # Grey is OFF this ramp, not at the bottom of it: it is what a
        # segment is blended toward when its LS residual is poor or when the
        # bed has entered the depth window. Against the pale end of Blues
        # the two are otherwise indistinguishable, which is the whole point
        # of the fade - hence a swatch, as fabric_sections.py carries a
        # two-dimensional key for the same reason.
        gax = lax.inset_axes([0.635, 0.52, 0.026, 0.30])
        gax.set_facecolor(GREY_RGB)
        gax.set_xticks([])
        gax.set_yticks([])
        for s in gax.spines.values():
            s.set_linewidth(0.6)
            s.set_color('0.4')
        gax.text(1.25, 0.5, 'unresolved', transform=gax.transAxes,
                 ha='left', va='center', fontsize=7.5, color=INK)

        # The scale bar has to stay ON the map - it means nothing off it -
        # so it gets the other half of the rule: drawn above every map
        # artist, over a white halo, so a track crossing it cannot swallow
        # it the way one swallowed the colour-bar title.
        sb_m = ab._nice_length(0.22 * max(span_x, span_y))
        xb = xs[1] - 0.06 * (xs[1] - xs[0]) - sb_m
        yb = ys[0] + 0.07 * (ys[1] - ys[0])
        ax.plot([xb, xb + sb_m], [yb, yb], '-', color=INK, lw=2.5,
                zorder=6, path_effects=HALO, transform=proj)
        ax.text(xb + sb_m / 2, yb - 0.028 * (ys[1] - ys[0]),
                '%g km' % (sb_m / 1000.0) if sb_m >= 1000 else '%g m' % sb_m,
                ha='center', va='top', fontsize=8.5, color=INK,
                zorder=6, path_effects=HALO, transform=proj)

        # Depth readout only. A survey-wide circular mean of theta0 was
        # drawn here too, and a bar tracking the window down the column;
        # both are gone. The mean angle averaged 36 frames whose axes are
        # a real spatial field, not repeat measurements of one number, so
        # its spread read as uncertainty when it is mostly variation - and
        # the bars already show the axis where it is measured. The column
        # bar restated what the depth text says.
        ax.text(xs[0] + 0.04 * (xs[1] - xs[0]),
                ys[1] - 0.055 * (ys[1] - ys[0]),
                '%.0f - %.0f m' % (lo, hi), fontsize=13, color=INK,
                ha='left', va='center', zorder=6, path_effects=HALO,
                transform=proj)

        ax.set_title(('%s: fabric contrast through the column\n' % TITLE)
                     + r'track colour: $\Delta\lambda$ per block; '
                     + r'bars: fabric axis $\theta_0$ at this depth'
                     + '\n'
                     + ('grey: LS residual fading %.2f to %.2f, or the bed '
                        'inside this window' % (RESID_GOOD, RESID_BAD)),
                     fontsize=10.5, color=INK)
        png = os.path.join(fdir, 'd%04d.png' % k)
        fig.savefig(png, dpi=150, facecolor='white')
        plt.close(fig)
        for r in range(1, HOLD):
            os.link(png, os.path.join(fdir, 'd%04d_%03d.png' % (k, r)))
        n += HOLD
        print('  [%2d/%d] %.0f-%.0f m, %d blocks, %d axes'
              % (k + 1, centres.size, lo, hi, cov, len(th_here)),
              flush=True)

    seq = os.path.join(fdir, 'seq_%05d.png')
    for i, p in enumerate(sorted(glob.glob(os.path.join(fdir, 'd*.png')))):
        os.rename(p, seq % i)
    mp4 = os.path.join(OUT, 'quadpol_depth_movie_%s.mp4' % SITE)
    subprocess.run([find_ffmpeg(), '-y', '-framerate', str(FPS), '-i', seq,
                    '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20',
                    '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', mp4],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    shutil.rmtree(fdir)
    print('wrote %s (%d frames, %.1f s)' % (mp4, n, n / FPS))


if __name__ == '__main__':
    main()

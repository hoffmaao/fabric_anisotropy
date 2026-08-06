"""Depth-averaged horizontal fabric and its orientation, in map view.

Modelled on Nymand et al. (2025), "Double Reflections in Polarized Radar
Data Reveal Ice Fabric in the North East Greenland Ice Stream", GRL 52,
doi:10.1029/2024GL110453, fig. 3: survey tracks drawn as ribbons coloured
by the depth-averaged horizontal asymmetry, with the two horizontal
principal axes drawn as a cross of double-headed arrows, and companion
panels for strength and orientation.

WHAT A SINGLE LEG CAN AND CANNOT MEASURE. An HH-VV pair on one line
measures only the projection of the horizontal ellipse onto that line's
own axes,

    dlam_obs(alpha) = -P cos 2(alpha - theta),

so a leg alone cannot separate a weak fabric aligned with it from a strong
fabric at 45 degrees to it. Two legs at different azimuths determine both,
so the figure keeps the distinction explicit at every site:

  ribbons   the per-leg PROJECTION, which is what that leg measured
  crosses   P and theta, solved only in grid cells that hold at least two
            azimuths separated by AZ_MIN_SEP, and drawn nowhere else

That is the honest version of Nymand et al.'s layout. Their double-
reflection method recovers orientation along a single track, so their
panels (b) and (c) are continuous; ours cannot be, and a ribbon coloured
by an orientation this survey never resolved would be an invention.

NO VELOCITY BACKGROUND, AND NO FLOW REFERENCE. Their panels sit on ITS_LIVE
surface speed and report orientation relative to flow. Ridge A is south of
satellite velocity coverage and is a divide site flowing at under ~2 m/yr,
where a flow direction is neither available nor meaningful, so orientation
is reported as an absolute azimuth (degrees east of north) over a plain
background.

WHICH SITES CAN SHOW ORIENTATION. Measured over each product set:

  ridge_a      raster, two azimuth families          strength + crosses
  taylor_dome  raster, 60-80 and 120-140 deg         strength + crosses
  negis        radiating traverse, 33 and 128-145    strength + crosses
  eastwind     17 lines over ~1 km, 20-140 deg       strength + crosses
  thwaites     single W-E transect, all 66 frames
               between 60 and 100 deg                RIBBONS ONLY

Thwaites gets one panel, not two. Its legs never differ by the ~20 deg the
two-parameter fit needs, so a strength panel there would be empty - and an
empty panel invites the reader to conclude the fabric is weak rather than
that the survey cannot see it.

Usage: python fabric_map.py <out_dir> [site ...]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import cartopy.crs as ccrs                # noqa: E402
import h5py                               # noqa: E402
import matplotlib.pyplot as plt           # noqa: E402
from matplotlib.collections import LineCollection   # noqa: E402
from matplotlib.colors import TwoSlopeNorm          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style as sty                  # noqa: E402
from scar_style import DATA, INK, MUTED   # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'

# name -> (title, staged survey file). The file is what
# opr_fabric/server/run_survey_fabric.m writes for that site.
SITES = {
    'ridge_a': ('Ridge A divide', 'ridge_a_survey.mat'),
    'taylor_dome': ('Taylor Dome', 'taylor_dome_survey.mat'),
    'negis': ('NEGIS onset (EastGRIP)', 'negis_survey.mat'),
    'thwaites': ('Thwaites eastern shear margin', 'thwaites_survey.mat'),
    'eastwind': ('Eastwind', 'eastwind_survey.mat'),
}
# Sites whose legs never differ by AZ_MIN_SEP; ribbons only, no solve.
NO_ORIENTATION = {'thwaites'}

# Depth window the ribbons average over. Below the firn (the shallowest
# interval sits inside the surface reference band and measures ~0 fringes
# by construction) and above the deepest one or two intervals, which the
# edge test showed are regulariser output rather than data - node coherence
# there falls to 0.39 and 0.09.
Z_AVG = (200.0, 1300.0)
Q_MIN = 0.30            # node coherence below which an interval is dropped
# Cell size trades sample count against legibility. At 3 km the survey
# yields ~26 cells on a ~30 km map, so the crosses sit ~3 km apart and have
# to be drawn smaller than they are readable. 5 km gives fewer, bigger
# crosses - closer to the dozen labelled points Nymand et al. draw - and
# each still holds both azimuths.
CELL_KM = 5.0           # grid cell for the two-azimuth orientation solve
AZ_MIN_SEP = 20.0       # degrees; below this two legs do not constrain theta
MIN_PER_AZ = 2          # blocks a leg must contribute to a cell
P_BOUND = 2.0 / 3       # |lam_max - lam_min| cannot exceed 1 - lam_z
# Frames whose inversion cannot fit its OWN dtau are dropped rather than
# drawn. Ridge A and Taylor Dome sit at a median 0.045 and 0.063 ns and
# never exceed 0.45; Thwaites has a median of 0.209 ns but a tail reaching
# 8.4 ns on the margin frames, and those few frames otherwise set the
# colour range for the whole map - +-0.39 against +-0.08 at Ridge A - so
# the sites would be drawn on scales that are not comparable and the
# well-fitted majority would wash out to nothing.
RMS_MAX_NS = 1.0

RIBBON_LW = 5.0
CROSS_KM = 1.9          # half-length of the longer principal-axis arrow


def load_survey(fn):
    if not os.path.exists(fn):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_survey_fabric.m for this '
            'site on the server and mirror the result there' % fn)
    out = []
    with h5py.File(fn) as f:
        g = f['F']
        for i in range(g['tag'].shape[0]):
            def a(name):
                return np.array(f[g[name][i, 0]])
            tag = ''.join(chr(c) for c in a('tag').ravel())
            out.append(dict(
                tag=tag,
                dlam=a('dlam').T, top=a('top').T, bot=a('bot').T,
                quality=a('quality').T,
                lat=a('lat').ravel(), lon=a('lon').ravel(),
                lat0=a('lat0').ravel(), lon0=a('lon0').ravel(),
                lat1=a('lat1').ravel(), lon1=a('lon1').ravel(),
                gate=float(a('gate').ravel()[0]),
                rms=float(a('rms').ravel()[0])))
    return out


def depth_average(d):
    """Coherence-weighted mean dlam over Z_AVG, per block.

    Weighted by node coherence AND by the interval's overlap with the
    window, so an interval straddling the window edge contributes the part
    that is inside it rather than all-or-nothing at its centre.
    """
    dl, q = d['dlam'], d['quality']
    top, bot = d['top'], d['bot']
    ov = np.minimum(bot, Z_AVG[1]) - np.maximum(top, Z_AVG[0])
    w = np.where(np.isfinite(q), np.clip(q, 0, 1), 0.0)
    w = np.where((w >= Q_MIN) & (ov > 0) & np.isfinite(dl), w * ov, 0.0)
    num = np.nansum(np.nan_to_num(dl) * w, axis=0)
    den = np.nansum(w, axis=0)
    out = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    return out, den


def seg_azimuth(lat0, lon0, lat1, lon1):
    """Azimuth of each block's own segment, degrees east of north, mod 180."""
    p0, p1 = np.radians(lat0), np.radians(lat1)
    dl = np.radians(lon1 - lon0)
    y = np.sin(dl) * np.cos(p1)
    x = np.cos(p0) * np.sin(p1) - np.sin(p0) * np.cos(p1) * np.cos(dl)
    return np.degrees(np.arctan2(y, x)) % 180.0


def solve_cell(az, y, w):
    """P and theta from dlam = -P cos 2(alpha - theta) over one cell.

    Returns None unless the cell holds two azimuths at least AZ_MIN_SEP
    apart. Without that the two columns of the design matrix are nearly
    parallel and the fit trades P against theta freely - the same
    ill-conditioning that put the EastGRIP free solve past the eigenvalue
    bound.
    """
    if az.size < 3:
        return None
    a = np.sort(np.unique(np.round(az)))
    spread = np.abs(a[:, None] - a[None, :])
    spread = np.minimum(spread, 180 - spread)
    if spread.max() < AZ_MIN_SEP:
        return None
    r = np.radians(az)
    A = np.c_[-np.cos(2 * r), -np.sin(2 * r)]
    W = np.sqrt(np.clip(w, 0, None))[:, None]
    sol, *_ = np.linalg.lstsq(A * W, y * W[:, 0], rcond=None)
    P = float(np.hypot(*sol))
    if not np.isfinite(P) or P > P_BOUND:
        return None
    th = float(np.degrees(0.5 * np.arctan2(sol[1], sol[0])) % 180)
    return P, th


def cells(lines, proj):
    """Per-cell (x, y, P, theta, n) on a CELL_KM grid."""
    X, Y, V, W, A = [], [], [], [], []
    for d in lines:
        v, w = depth_average(d)
        xy = proj.transform_points(ccrs.PlateCarree(), d['lon'], d['lat'])
        az = seg_azimuth(d['lat0'], d['lon0'], d['lat1'], d['lon1'])
        ok = np.isfinite(v) & (w > 0) & np.isfinite(az)
        X.append(xy[ok, 0])
        Y.append(xy[ok, 1])
        V.append(v[ok])
        W.append(w[ok])
        A.append(az[ok])
    X = np.concatenate(X)
    Y = np.concatenate(Y)
    V = np.concatenate(V)
    W = np.concatenate(W)
    A = np.concatenate(A)

    s = CELL_KM * 1e3
    ix = np.floor(X / s).astype(int)
    iy = np.floor(Y / s).astype(int)
    out = []
    for key in set(zip(ix.tolist(), iy.tolist())):
        m = (ix == key[0]) & (iy == key[1])
        if m.sum() < 2 * MIN_PER_AZ:
            continue
        fit = solve_cell(A[m], V[m], W[m])
        if fit is None:
            continue
        # Anchor on the CENTROID of the blocks that solved, not the cell
        # centre. Where a survey only clips the corner of a cell - which is
        # most of them on a track-following survey - the cell centre can be
        # kilometres off the nearest track, and a cross floating in empty
        # space reads as a measurement of ground that was never sounded.
        out.append((float(X[m].mean()), float(Y[m].mean()),
                    fit[0], fit[1], int(m.sum())))
    return out


def ribbons(ax, lines, norm, cmap, lw=RIBBON_LW, zorder=6):
    """Each block drawn as its own coloured segment along the track."""
    segs, cols = [], []
    for d in lines:
        v, w = depth_average(d)
        ok = np.isfinite(v) & (w > 0)
        for k in np.flatnonzero(ok):
            segs.append([(d['lon0'][k], d['lat0'][k]),
                         (d['lon1'][k], d['lat1'][k])])
            cols.append(v[k])
    lc = LineCollection(segs, transform=ccrs.PlateCarree(), linewidths=lw,
                        colors=cmap(norm(np.array(cols))),
                        capstyle='round', zorder=zorder)
    ax.add_collection(lc)
    return np.array(cols)


def draw_cross(ax, x, y, P, th, proj, scale, pmax):
    """The two horizontal principal axes as a cross of double-headed arrows.

    Long dark arrow = azimuth of the LARGER horizontal eigenvalue; short
    pale arrow = the smaller, orthogonal to it. Length encodes P on the
    long axis only, so a near-isotropic cell reads as two arms of similar
    length rather than as a strong fabric pointing nowhere.
    """
    # P spans only ~0-0.11 here, so scaling length by P/P_BOUND would make
    # every cross the same stub. Normalise to the observed spread instead,
    # keeping a floor so a weak cell still shows its axis.
    L = scale * (0.45 + 0.55 * min(P / max(pmax, 1e-6), 1.0))
    tr = proj._as_mpl_transform(ax)
    for ang, col, lw, frac in ((th + 90, '#f0f0f0', 4.2, 0.60),
                               (th, '#f0f0f0', 5.4, 1.0),
                               (th + 90, '#8a8a8a', 1.6, 0.60),
                               (th, '#111111', 2.6, 1.0)):
        r = np.radians(ang)
        dx, dy = np.sin(r) * L * frac, np.cos(r) * L * frac
        ax.plot([x - dx, x + dx], [y - dy, y + dy], '-', color=col, lw=lw,
                solid_capstyle='round', transform=proj,
                zorder=8 if col == '#f0f0f0' else 9)
    del tr


def render(site):
    title, fname = SITES[site]
    lines = load_survey(os.path.join(DATA, fname))
    n_all = len(lines)
    drop = [d for d in lines if not (d['rms'] <= RMS_MAX_NS)]
    lines = [d for d in lines if d['rms'] <= RMS_MAX_NS]
    if not lines:
        print('%s: every frame exceeds the %.1f ns misfit gate; skipped'
              % (site, RMS_MAX_NS))
        return
    if drop:
        print('%s: dropped %d of %d frames over %.1f ns misfit (%s)'
              % (site, len(drop), n_all, RMS_MAX_NS,
                 ', '.join('%s %.1f' % (d['tag'], d['rms']) for d in drop)))
    print('\n=== %s: %d frames, %d blocks'
          % (site, len(lines), sum(d['lat'].size for d in lines)))

    # NEGIS is the one northern-hemisphere site here, so its tracks would
    # transform to nonsense under a south polar projection.
    proj = (ccrs.NorthPolarStereo(central_longitude=-45) if site == 'negis'
            else ccrs.SouthPolarStereo())
    two = site not in NO_ORIENTATION
    cl = cells(lines, proj) if two else []
    if two:
        print('%d grid cells (%.0f km) resolved orientation'
              % (len(cl), CELL_KM))

    allv = np.concatenate([depth_average(d)[0] for d in lines])
    allv = allv[np.isfinite(allv)]
    lim = float(np.nanpercentile(np.abs(allv), 96))
    norm = TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim)
    cmap = plt.get_cmap('RdBu_r')
    print('depth-averaged dlam over %.0f-%.0f m: %.3f..%.3f (colour +-%.3f)'
          % (*Z_AVG, allv.min(), allv.max(), lim))

    lat = np.concatenate([d['lat'] for d in lines])
    lon = np.concatenate([d['lon'] for d in lines])
    xy = proj.transform_points(ccrs.PlateCarree(), lon, lat)
    pad = 4e3
    extent = (xy[:, 0].min() - pad, xy[:, 0].max() + pad,
              xy[:, 1].min() - pad, xy[:, 1].max() + pad)

    # Two panels of equal width where orientation is solvable, one where it
    # is not. The figure is sized from the survey's own aspect: cartopy
    # holds the map aspect fixed, so a default landscape canvas around a
    # 112 x 10 km transect (Thwaites, which polar stereographic rotates to
    # near-vertical) leaves a sliver of data in a field of grey.
    # Drive the layout from a FIXED map height and let width follow the
    # aspect, rather than the reverse. Sizing from a fixed width instead
    # made Thwaites 15 in tall, because cartopy holds the map aspect and
    # that survey is 112 x 10 km rotated to near-vertical by the
    # projection. The width floor keeps the title and colourbar legible on
    # the narrowest site.
    asp = (extent[3] - extent[2]) / max(extent[1] - extent[0], 1.0)
    map_h = 6.0
    pw = float(np.clip(map_h / max(asp, 1e-3), 2.9, 6.6))
    fig = plt.figure(figsize=((2 if two else 1) * pw + 0.6, map_h + 1.5),
                     layout='constrained')
    gs = fig.add_gridspec(1, 2 if two else 1)
    axm = fig.add_subplot(gs[0, 0], projection=proj)
    axp = fig.add_subplot(gs[0, 1], projection=proj) if two else None

    for ax in filter(None, (axm, axp)):
        ax.set_extent(extent, crs=proj)
        ax.set_facecolor('#f2f1ef')
        gl = ax.gridlines(draw_labels=False, lw=0.5, color='0.85')
        gl.top_labels = gl.right_labels = False

    # (a) ribbons + orientation crosses. No angle labels: the crosses carry
    # the orientation, and the numbers only repeated it while colliding
    # with neighbouring arms.
    ribbons(axm, lines, norm, cmap)
    scale = CROSS_KM * 1e3
    pmax = max([c[2] for c in cl]) if cl else 1.0
    for x, y, P, th, n in cl:
        draw_cross(axm, x, y, P, th, proj, scale, pmax)
    sty.scale_bar_br(axm, extent, proj)
    sub = ('with horizontal principal axes where two azimuths cross' if two
           else 'single-azimuth transect: orientation not recoverable')
    axm.set_title('%s: depth-averaged $\\Delta\\lambda$ (%.0f-%.0f m)\n%s'
                  % (title, Z_AVG[0], Z_AVG[1], sub),
                  fontsize=10.5, color=INK)
    cax = axm.inset_axes([0.05, -0.075, 0.90, 0.026])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation='horizontal',
                      label=r'$\Delta\lambda$ projected on each leg')
    # Three explicit ticks. Matplotlib's default spacing collided on the
    # narrow single-panel layout, printing "-0.250.00 0.25" as one run.
    cb.set_ticks([-lim, 0, lim])
    cb.set_ticklabels(['%.2f' % -lim, '0', '%.2f' % lim])

    # (b) solved strength
    if two:
        pv = np.array([c[2] for c in cl]) if cl else np.array([])
        if pv.size:
            sc = axp.scatter([c[0] for c in cl], [c[1] for c in cl], c=pv,
                             s=110, cmap='viridis', vmin=0,
                             vmax=float(np.nanpercentile(pv, 95)),
                             edgecolor='black', lw=0.6, transform=proj,
                             zorder=7)
            cax2 = axp.inset_axes([0.05, -0.075, 0.90, 0.026])
            fig.colorbar(sc, cax=cax2, orientation='horizontal',
                         label=r'$P=\lambda_{max}-\lambda_{min}$')
        # Same ribbon weight as (a). The ribbons are the same measurement in
        # both panels, so drawing them thinner here would read as a
        # different, lesser quantity rather than the same one behind other
        # symbols.
        ribbons(axp, lines, norm, cmap, zorder=4)
        axp.set_title('fabric strength, two-azimuth solve', fontsize=10.5,
                      color=INK)

    out = os.path.join(OUT, 'scar_%s_fabric_map.png' % site)
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('wrote', out)


def main():
    os.makedirs(OUT, exist_ok=True)
    want = sys.argv[2:] or list(SITES)
    bad = [s for s in want if s not in SITES]
    if bad:
        raise SystemExit('unknown site(s): %s; known: %s'
                         % (', '.join(bad), ', '.join(SITES)))
    for site in want:
        # A site whose survey has not been run yet is skipped rather than
        # fatal, so rendering the finished ones does not wait on the batch
        # still in flight.
        if not os.path.exists(os.path.join(DATA, SITES[site][1])):
            print('no survey staged for %s; skipped' % site)
            continue
        render(site)


if __name__ == '__main__':
    main()

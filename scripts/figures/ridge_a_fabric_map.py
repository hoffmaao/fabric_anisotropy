"""Ridge A: depth-averaged horizontal fabric and its orientation, in map view.

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
fabric at 45 degrees to it. Two legs at different azimuths determine both.
Ridge A is a raster survey, so this figure keeps the distinction explicit:

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

Usage: python ridge_a_fabric_map.py <out_dir>
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
FN = os.path.join(DATA, 'ridge_a_survey.mat')

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

RIBBON_LW = 5.0
CROSS_KM = 1.9          # half-length of the longer principal-axis arrow


def load_survey():
    if not os.path.exists(FN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_ridge_a_survey.m on the '
            'server and mirror the result there' % FN)
    out = []
    with h5py.File(FN) as f:
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
        out.append(((key[0] + 0.5) * s, (key[1] + 0.5) * s,
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


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = load_survey()
    print('loaded %d frames, %d blocks'
          % (len(lines), sum(d['lat'].size for d in lines)))

    proj = ccrs.SouthPolarStereo()
    cl = cells(lines, proj)
    print('%d grid cells (%.0f km) resolved orientation' % (len(cl), CELL_KM))

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

    fig = plt.figure(figsize=(15.0, 6.6), layout='constrained')
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0])
    axm = fig.add_subplot(gs[0, 0], projection=proj)
    axp = fig.add_subplot(gs[0, 1], projection=proj)
    axt = fig.add_subplot(gs[0, 2], projection=proj)

    for ax in (axm, axp, axt):
        ax.set_extent(extent, crs=proj)
        ax.set_facecolor('#f2f1ef')
        gl = ax.gridlines(draw_labels=False, lw=0.5, color='0.85')
        gl.top_labels = gl.right_labels = False

    # (a) ribbons + orientation crosses
    ribbons(axm, lines, norm, cmap)
    scale = CROSS_KM * 1e3
    pmax = max([c[2] for c in cl]) if cl else 1.0
    for x, y, P, th, n in cl:
        draw_cross(axm, x, y, P, th, proj, scale, pmax)
    # A few angle labels, as Nymand et al. annotate theirs. Not all of them:
    # at this cell size the crosses are ~5 km apart and a label on every one
    # would collide with its neighbours' arms.
    for x, y, P, th, n in cl[::3]:
        axm.annotate('%.0f$^\\circ$' % th, xy=(x, y),
                     # Offset in POINTS, a display unit. `scale` is in
                     # metres; passing it here pushed the label ~1200 pt
                     # away and collapsed constrained_layout to zero.
                     xytext=(0, 13), textcoords='offset points',
                     xycoords=proj._as_mpl_transform(axm), fontsize=7,
                     color='#111111', ha='center', va='bottom', zorder=10,
                     bbox=dict(boxstyle='round,pad=0.14', fc='white',
                               ec='none', alpha=0.72))
    sty.scale_bar_br(axm, extent, proj)
    axm.set_title('Ridge A: depth-averaged $\\Delta\\lambda$ (%.0f-%.0f m)\n'
                  'with horizontal principal axes where two azimuths cross'
                  % Z_AVG, fontsize=10.5, color=INK)
    cax = axm.inset_axes([0.06, -0.09, 0.55, 0.028])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                 orientation='horizontal',
                 label=r'$\Delta\lambda$ projected on each leg')

    # (b) solved strength
    pv = np.array([c[2] for c in cl]) if cl else np.array([])
    if pv.size:
        sc = axp.scatter([c[0] for c in cl], [c[1] for c in cl], c=pv,
                         s=110, cmap='viridis', vmin=0,
                         vmax=float(np.nanpercentile(pv, 95)),
                         edgecolor='black', lw=0.6, transform=proj, zorder=7)
        cax2 = axp.inset_axes([0.06, -0.09, 0.55, 0.028])
        fig.colorbar(sc, cax=cax2, orientation='horizontal',
                     label=r'$P=\lambda_{max}-\lambda_{min}$')
    ribbons(axp, lines, norm, cmap, lw=1.6, zorder=4)
    axp.set_title('fabric strength, two-azimuth solve', fontsize=10.5,
                  color=INK)

    # (c) solved orientation. Cyclic colormap: theta is an AXIS, not a
    # direction, so 179 deg and 1 deg are neighbours and a linear ramp
    # would draw them as opposite extremes.
    tv = np.array([c[3] for c in cl]) if cl else np.array([])
    if tv.size:
        # Drawn as oriented ticks, not dots. theta is an angle, and this
        # survey's values sit inside a ~30 deg spread, so on any cyclic
        # ramp they resolve to nearly one colour - the reader would be
        # asked to distinguish orientations by a hue difference smaller
        # than the colourbar can show. The tick states the angle directly
        # and the colour only reinforces it.
        tk = CROSS_KM * 1e3 * 1.05
        r = np.radians(tv)
        xs = np.array([c[0] for c in cl])
        ys = np.array([c[1] for c in cl])
        segs = [[(x - np.sin(a) * tk, y - np.cos(a) * tk),
                 (x + np.sin(a) * tk, y + np.cos(a) * tk)]
                for x, y, a in zip(xs, ys, r)]
        cm2 = plt.get_cmap('twilight')
        lc2 = LineCollection(segs, transform=proj, linewidths=3.4,
                             colors=cm2(tv / 180.0), capstyle='round',
                             zorder=7)
        axt.add_collection(lc2)
        sc2 = plt.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(0, 180), cmap=cm2)
        cax3 = axt.inset_axes([0.06, -0.09, 0.55, 0.028])
        cb = fig.colorbar(sc2, cax=cax3, orientation='horizontal',
                          label=r'$\theta$ of $\lambda_{max}$ (deg E of N)')
        cb.set_ticks([0, 45, 90, 135, 180])
    ribbons(axt, lines, norm, cmap, lw=1.6, zorder=4)
    axt.set_title('fabric orientation, two-azimuth solve', fontsize=10.5,
                  color=INK)

    out = os.path.join(OUT, 'scar_ridge_a_fabric_map.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', out)


if __name__ == '__main__':
    main()

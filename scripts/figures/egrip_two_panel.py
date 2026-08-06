"""EastGRIP site figure: survey map + HH-VV interferogram.

Same two-panel form as scar_two_panel.py and negis_two_panel.py - survey
map left, interferogram right, no legend, survey lines white, focused
profile black, red start / white end markers repeated on the map and at
the interferogram's edges.

Two things differ, both because this figure is anchored to a drill site
rather than to a survey:
  - the map is windowed on the borehole rather than on the whole survey,
    so no zoom inset is needed - the profile is already legible at this
    scale;
  - the borehole is marked with a star and NO text label, and the
    interferogram's vertical axis is DEPTH rather than traveltime, so the
    panel can be set directly against core measurements.

The line is 20240619_01_001, the nearest phase-preserving line that also
inverts stably; see egrip_two_panel.py for why the two that pass closer
were rejected.

Usage: python egrip_two_panel.py <out_dir>
"""
import os
import sys
import warnings

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.colors import LogNorm    # noqa: E402
from scipy.io import loadmat             # noqa: E402

warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cartopy.crs as ccrs               # noqa: E402
import scar_style as sty                 # noqa: E402
from scar_style import (DPI, PROFILE_C, PROFILE_FX,  # noqa: E402
                        PROFILE_LW, TRACK_C)
from egrip_azimuthal import (C_ICE, DATA, EG_LAT, EG_LON,  # noqa: E402
                             haversine_km, load_line)

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
FRAME = '20240619_01_001'
CAMP_SEG = '20240618_01'
BORE_FACE = '#ffd400'
HALF_M = 11e3
PROJ = ccrs.Stereographic(central_latitude=90, central_longitude=-45,
                          true_scale_latitude=70)


def tracks_and_camp():
    T = loadmat(os.path.join(DATA, 'negis2024_tracks_fixed.mat'),
                squeeze_me=True)['T']
    tr, cla, clo = [], [], []
    for r in T:
        la = np.atleast_1d(np.asarray(r['lat'], float))
        lo = np.atleast_1d(np.asarray(r['lon'], float))
        ok = np.isfinite(la) & np.isfinite(lo)
        if ok.sum() > 2:
            tr.append((la[ok], lo[ok]))
        if str(r['seg']) == CAMP_SEG:
            cla.extend(la[ok]); clo.extend(lo[ok])
    camp = ((float(np.mean(cla)), float(np.mean(clo))) if cla
            else (EG_LAT, EG_LON))
    return tr, camp


def main():
    os.makedirs(OUT, exist_ok=True)
    d = load_line(FRAME)
    if d is None:
        raise SystemExit('no extract for %s' % FRAME)
    tr, (bore_lat, bore_lon) = tracks_and_camp()

    lat, lon = d['lat'], d['lon']
    step = haversine_km(lat[1:], lon[1:], lat[:-1], lon[:-1])
    dist = np.concatenate([[0], np.cumsum(step)])
    tb = (d['t'] - np.nanmedian(d['surf'])) * 1e6
    depth = tb * C_ICE / 2 * 1e-6
    print('%s: %.2f km, borehole %.2f km from the near end'
          % (FRAME, dist[-1],
             haversine_km(lat, lon, bore_lat, bore_lon).min()))

    fig = plt.figure(figsize=sty.FIGSIZE, dpi=DPI, layout='constrained')
    gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.18])
    axm = fig.add_subplot(gs[0, 0], projection=PROJ)
    axi = fig.add_subplot(gs[0, 1])

    z = np.load(os.path.join(DATA, 'itslive_negis_win.npz'))
    wx, wy, v = z['x'], z['y'], z['v']
    vx, vy = z['vx'], z['vy']
    px, py = PROJ.transform_point(bore_lon, bore_lat, ccrs.PlateCarree())[:2]
    extent = (px - HALF_M, px + HALF_M, py - HALF_M, py + HALF_M)
    axm.set_extent(extent, crs=PROJ)
    pcm = axm.pcolormesh(wx, wy, np.maximum(v, 1e-3), transform=PROJ,
                         cmap='magma', norm=LogNorm(vmin=2, vmax=120),
                         shading='auto', rasterized=True)
    dec = 8
    Xq, Yq = np.meshgrid(wx[::dec], wy[::dec])
    U, V = vx[::dec, ::dec], vy[::dec, ::dec]
    sp = np.hypot(U, V)
    ok = sp > 1.0
    axm.quiver(Xq, Yq, np.where(ok, U / sp, np.nan),
               np.where(ok, V / sp, np.nan), transform=PROJ, color='white',
               scale=26, width=0.005, alpha=0.6, zorder=6)
    cax = axm.inset_axes([0.05, -0.10, 0.62, 0.03])
    fig.colorbar(pcm, cax=cax, orientation='horizontal',
                 label='surface speed (m/yr)')

    for la, lo in tr:
        axm.plot(lo, la, '-', color=TRACK_C, lw=1.1, alpha=0.9,
                 transform=ccrs.PlateCarree(), zorder=7)
    axm.plot(lon, lat, '-', color=PROFILE_C, lw=PROFILE_LW,
             transform=ccrs.PlateCarree(), zorder=8, path_effects=PROFILE_FX)
    sty.draw_map_ends(axm, lat, lon)
    # Star only - the borehole needs no caption on a slide that is titled
    # for it, and the label crowded the profile's start marker
    axm.scatter([bore_lon], [bore_lat], marker='*', s=360, c=BORE_FACE,
                edgecolor='black', lw=1.4, transform=ccrs.PlateCarree(),
                zorder=11)
    gl = axm.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='gray',
                       y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    sty.scale_bar_br(axm, extent, PROJ, nice_length=sty.nice_length_up)
    axm.set_title('EastGRIP drill site', fontsize=12)

    st = max(1, d['ifg'].shape[0] // 2200)
    sx = max(1, d['ifg'].shape[1] // 2000)
    rgb = sty.hsv_plot_coherence(np.angle(d['ifg'][::st, ::sx]),
                                 d['coh'][::st, ::sx])
    axi.imshow(rgb, aspect='auto', interpolation='nearest',
               extent=[dist[0], dist[-1], depth[-1], depth[0]])
    axi.set_ylim(1300, 0)
    axi.set_xlabel('distance along profile (km)')
    axi.set_ylabel('depth (m)')
    sm = plt.cm.ScalarMappable(cmap='hsv', norm=plt.Normalize(-np.pi, np.pi))
    cb = fig.colorbar(sm, ax=axi, pad=0.035)
    cb.set_label('phase change (rad)')
    cb.set_ticks([-np.pi, 0, np.pi])
    cb.set_ticklabels([r'$-\pi$', '0', r'$\pi$'])
    sty.draw_ifg_ends(axi, dist)
    axi.set_title('borehole line %s' % FRAME, fontsize=12, pad=26)

    out = os.path.join(OUT, 'scar_egrip_2panel.png')
    fig.savefig(out, dpi=DPI)
    print('wrote', out)


if __name__ == '__main__':
    main()

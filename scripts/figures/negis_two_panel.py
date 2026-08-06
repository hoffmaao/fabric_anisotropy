"""NEGIS two-panel SCAR figure: survey overview + HH-VV interferogram.

Same slide design as scar_two_panel.py - survey lines white, focused
profile black, red start / white end markers on the map inset and at the
matching edges of the interferogram, no legend, plain black scale bar -
and literally the same constants and helpers, out of scar_style.py. Only
the projection differs (polar stereographic north, not EPSG:3031).

The interferogram comes from negis_interferogram.m on mem1: the
2024_Greenland_Ground2 qlook products are COMPLEX for every segment
except 20240618_01, so HH x conj(VV) is formable directly. Positions are
GPS-file interpolations onto the product GPS_time, because the records
time sync failed for these segments (validated against CSARP_layer on
20240618_01, agreement ~1 m).

Usage: python negis_two_panel.py <out_dir> [seg_frm]
"""
import os
import sys
import warnings

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
from matplotlib.colors import LogNorm      # noqa: E402
from scipy.io import loadmat               # noqa: E402

warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cartopy.crs as ccrs                # noqa: E402
import scar_style as sty                  # noqa: E402
from scar_style import (DATA, DPI, PROFILE_C, PROFILE_FX,  # noqa: E402
                        PROFILE_LW, TRACK_C, zoom_inset)

OUT = sys.argv[1]
# 20240626_03_001: the along-flow line beside the southeastern shear
# margin. Picked over the cross-margin crossing (20240620_01_003) because
# it holds coherence over the WHOLE profile - every along-track column
# stays above 0.35 down to 10 us and 82% of them to 15 us, against 70%
# and 51% for the crossing, whose interferogram breaks into incoherent
# blocks. It is also the longest (13.7 km), straightest (0.99) and best
# flow-aligned (0.97) frame in the survey.
TAG = sys.argv[2] if len(sys.argv) > 2 else '20240626_03_001'

PROJ = ccrs.Stereographic(central_latitude=90, central_longitude=-45,
                          true_scale_latitude=70)


def block_sum(a, nr, na):
    """Sum over (nr x na) cells, trimming the ragged edge."""
    ntc, nxc = a.shape[0] // nr, a.shape[1] // na
    return a[:ntc * nr, :nxc * na].reshape(ntc, nr, nxc, na).sum(axis=(1, 3))


def drop_stationary(lat, lon, min_step_m=1.0):
    """Keep only columns the vehicle actually moved between.

    The traverse stopped three times inside this frame (~230, 300 and 160
    traces at 0.0, 5.1 and 10.2 km, about a metre of movement each). The
    radar kept writing, so those traces multilook into high COHERENCE with
    random range phase - bright horizontal rainbow stripes, not gaps, and
    high coherence is exactly why a coherence threshold does not catch
    them. Worse, imshow spreads columns evenly over [dist[0], dist[-1]],
    so a stop that covers no ground still eats ~0.7 km of the axis. This
    is the cull Knut Christianson applies in the NEGIS reposition script.

    negis_interferogram.m now applies the same criterion to the RAW traces
    before its 4x4 multilook, which is the only place it can keep a
    stopped stretch out of the azimuth cell that straddles its edge. This
    stays as a no-op safeguard for extracts written before that change.
    """
    R = 6371e3
    p = np.radians(lat)
    dl = np.radians(np.diff(lon))
    dp = np.diff(p)
    a = np.sin(dp / 2)**2 + np.cos(p[:-1]) * np.cos(p[1:]) * np.sin(dl / 2)**2
    step = 2 * R * np.arcsin(np.sqrt(a))
    keep = np.concatenate([[True], step >= min_step_m])
    return keep


def load_ifg(tag, extra_looks=(1, 1)):
    """Interferogram, coherence and geometry, optionally re-multilooked.

    negis_interferogram.m ships the un-normalised power sums alongside the
    complex interferogram precisely so the look count can be raised here:
    coherence has to be re-formed from sums over the coarser cell, never
    by averaging the finer cells' coherence (that would average a ratio of
    noisy magnitudes and sit high wherever the phase is pure noise).
    """
    fn = os.path.join(DATA, f'negis_ifg_{tag}.mat')
    with h5py.File(fn) as f:
        # MATLAB -v7.3 stores matrices transposed in HDF5
        z = f['interferogram_mlook'][()]
        ifg = (z['real'] + 1j * z['imag']).T.astype(np.complex128)
        t = f['Time'][()].ravel()
        lat = f['Latitude_ml'][()].ravel()
        lon = f['Longitude_ml'][()].ravel()
        surf = f['Surface_ml'][()].ravel()
        if 'power_hh' in f:
            ph = f['power_hh'][()].T.astype(np.float64)
            pv = f['power_vv'][()].T.astype(np.float64)
        else:
            ph = pv = None
    # Cull before multilooking, so a stopped stretch cannot contaminate the
    # cell that straddles its edge
    keep = drop_stationary(lat, lon)
    n_drop = int((~keep).sum())
    if n_drop:
        print('  dropped %d of %d columns where the vehicle was stopped'
              % (n_drop, keep.size))
        ifg, lat, lon, surf = ifg[:, keep], lat[keep], lon[keep], surf[keep]
        if ph is not None:
            ph, pv = ph[:, keep], pv[:, keep]
    nr, na = extra_looks
    if nr > 1 or na > 1:
        if ph is None:
            raise SystemExit('extract lacks power_hh/power_vv; re-run '
                             'negis_interferogram.m before re-multilooking')
        ifg = block_sum(ifg, nr, na)
        ph = block_sum(ph, nr, na)
        pv = block_sum(pv, nr, na)
        t = t[:len(t) // nr * nr].reshape(-1, nr).mean(axis=1)
        lat = lat[:len(lat) // na * na].reshape(-1, na).mean(axis=1)
        lon = lon[:len(lon) // na * na].reshape(-1, na).mean(axis=1)
        surf = surf[:len(surf) // na * na].reshape(-1, na).mean(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        coh = np.abs(ifg) / np.sqrt(ph * pv) if ph is not None else None
    return ifg, coh, t, lat, lon, surf


def main():
    os.makedirs(OUT, exist_ok=True)
    # On top of the 4x4 the extract already carries: 128 looks total. Held
    # down in azimuth on purpose - going heavier smooths the fringes into
    # 180 m columns, which is coarser than the along-track structure the
    # profile exists to show, and this frame's coherence (0.66-0.77 over
    # the top 10 us) does not need the extra looks.
    ifg, coh, t, lats, lons, surf = load_ifg(TAG, extra_looks=(2, 4))
    dist = sty.cumdist_km(lats, lons)
    tb_all = (t - np.nanmedian(surf)) * 1e6
    shallow = coh[(tb_all >= 0) & (tb_all < 15)]
    print('%s: %s ifg, %.2f km, coherence median %.3f (0-15 us below surface)'
          % (TAG, ifg.shape, dist[-1], np.nanmedian(shallow)))

    # All survey tracks, repositioned the same way
    T = loadmat(os.path.join(DATA, 'negis2024_tracks_fixed.mat'),
                squeeze_me=True)['T']
    tracks = []
    for r in T:
        la = np.atleast_1d(np.asarray(r['lat'], float))
        lo = np.atleast_1d(np.asarray(r['lon'], float))
        ok = np.isfinite(la) & np.isfinite(lo)
        if ok.sum() > 2:
            tracks.append((la[ok], lo[ok]))
    tla = np.concatenate([a[0] for a in tracks])
    tlo = np.concatenate([a[1] for a in tracks])

    z = np.load(os.path.join(DATA, 'itslive_negis_win.npz'))
    wx, wy, v = z['x'], z['y'], z['v']
    vx, vy = z['vx'], z['vy']
    norm = LogNorm(vmin=2, vmax=120)

    def speed_layer(ax, dec):
        pcm = ax.pcolormesh(wx, wy, np.maximum(v, 1e-3), transform=PROJ,
                            cmap='magma', norm=norm, shading='auto',
                            rasterized=True)
        Xq, Yq = np.meshgrid(wx[::dec], wy[::dec])
        U, V = vx[::dec, ::dec], vy[::dec, ::dec]
        sp = np.hypot(U, V)
        ok = sp > 1.0
        U = np.where(ok, U / sp, np.nan)
        V = np.where(ok, V / sp, np.nan)
        ax.quiver(Xq, Yq, U, V, transform=PROJ, color='white',
                  scale=26, width=0.005, alpha=0.6, zorder=6)
        return pcm

    pts = PROJ.transform_points(ccrs.PlateCarree(), tlo, tla)
    pad = 0.16 * max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
    # Clipped to the speed window: padding past its edge would put a band
    # of blank page inside the map frame
    extent = (max(pts[:, 0].min() - pad, wx.min()),
              min(pts[:, 0].max() + pad, wx.max()),
              max(pts[:, 1].min() - pad, wy.min()),
              min(pts[:, 1].max() + pad, wy.max()))

    fig, ax, axi = sty.two_panel_figure(PROJ)

    ax.set_extent(extent, crs=PROJ)
    pcm = speed_layer(ax, dec=40)
    cax = ax.inset_axes([0.05, -0.12, 0.62, 0.03])
    fig.colorbar(pcm, cax=cax, orientation='horizontal',
                 label='surface speed (m/yr)')
    for la, lo in tracks:
        ax.plot(lo, la, '-', color=TRACK_C, lw=1.1, alpha=0.9,
                transform=ccrs.PlateCarree(), zorder=7)
    ax.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
            transform=ccrs.PlateCarree(), zorder=8,
            path_effects=PROFILE_FX)
    gl = ax.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='gray',
                      y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False   # they rotate into the speed colorbar
    sty.scale_bar_br(ax, extent, PROJ)
    sty.continental_inset(ax, float(np.nanmean(lats)),
                          float(np.nanmean(lons)),
                          'greenland', (0.02, 0.035, 0.24, 0.235))
    ax.set_title('NEGIS onset survey (EGRIP)', fontsize=12)

    # Zoom inset so the profile's two ends separate, as on the other sites
    px, py = PROJ.transform_point(np.mean(lons), np.mean(lats),
                                  ccrs.PlateCarree())[:2]
    # Wide enough that both end markers clear the inset frame; the profile
    # very nearly fills the window at 0.5
    half = 0.62 * 1e3 * max(6.0, dist[-1])
    zext = (px - half, px + half, py - half, py + half)
    # Top-left: the profile sits in the southeast, so the locator box lands
    # bottom-right where the scale bar is - putting the inset there too
    # would stack three things in one corner
    axz = zoom_inset(ax, [0.02, 0.66, 0.40, 0.30], PROJ)
    axz.set_extent(zext, crs=PROJ)
    # sparser than the main map's density would give at this zoom, or the
    # arrows carpet it
    speed_layer(axz, dec=14)
    for la, lo in tracks:
        axz.plot(lo, la, '-', color=TRACK_C, lw=1.4, alpha=0.9,
                 transform=ccrs.PlateCarree(), zorder=7)
    axz.plot(lons, lats, '-', color=PROFILE_C, lw=PROFILE_LW,
             transform=ccrs.PlateCarree(), zorder=8,
             path_effects=PROFILE_FX)
    sty.draw_map_ends(axz, lats, lons)
    # Dashed locator, solid inset frame, so the two boxes are telling apart
    sty.inset_frame(axz)
    sty.locator_box(ax, zext, PROJ)

    # ---- interferogram
    tb = (t - np.nanmedian(surf)) * 1e6
    sty.ifg_panel(
        axi, fig, dist, tb, np.angle(ifg), coh,
        ylabel='TWTT below surface (μs)',
        cb_label='phase change (rad)')
    axi.set_ylim(min(22.0, tb[-1]), -0.5)
    sty.draw_ifg_ends(axi, dist)
    axi.set_title('along-flow profile, frame %s' % TAG, fontsize=12,
                  pad=26)

    out = os.path.join(OUT, 'scar_negis_2panel.png')
    fig.savefig(out, dpi=DPI)
    print('wrote', out)


if __name__ == '__main__':
    main()

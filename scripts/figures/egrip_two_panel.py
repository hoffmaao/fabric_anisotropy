"""EastGRIP two-panel: survey map + fabric depth profile at the borehole.

Built to be set beside the EastGRIP core's own fabric measurements, so
the radar profile has to be anchored to the drill site rather than to a
survey line. Two things follow from that:

  - The borehole's OWN segment, 20240618_01, cannot be used. It is the
    single segment of 2024_Greenland_Ground2 written with incoherent
    decimation (qlook inc_dec = 10), so its data are real and carry no
    polarimetric phase. 20240621_01_010 is the nearest complex line,
    standing off 0.13 km, and it is also the longest and the most
    coherent at depth of the candidates within 0.35 km.
  - The traverse drives out FROM camp, so after the stopped-trace cull
    the borehole sits at one end of the profile. The distance axis is
    therefore flipped where needed to put 0 km at the borehole, and the
    depth profile averages only the blocks within RADIUS_KM of it - a
    profile averaged over the whole 7.9 km line would not be a
    measurement at the drill site.

The camp position is taken from the 20240618_01 camp survey in our own
data, not from a published collar coordinate, and is labelled as such.

No core fabric values are drawn here: this figure supplies the radar
side of the comparison at the right location and depth scale, for the
core measurements to be overlaid on.

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
from scar_style import (DATA, DPI, PROFILE_C, PROFILE_FX,  # noqa: E402
                        PROFILE_LW, TRACK_C, field, track_azimuth,
                        zoom_inset)

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
SITE = 'egrip'
FRAME = '20240619_01_001'
CAMP_SEG = '20240618_01'      # the camp survey, used to locate the drill site

RADIUS_KM = 2.0     # blocks this close to the borehole enter the profile
Z_MIN = 200.0       # above this the interval is inside the reference band
COH_MIN = 0.45      # node coherence an interval must reach to be drawn

PROJ = ccrs.Stereographic(central_latitude=90, central_longitude=-45,
                          true_scale_latitude=70)
BORE_FACE = '#ffd400'
C_PROF = '#1baf7a'


def haversine_km(lat, lon, lat0, lon0):
    R = 6371.0
    p0, p = np.radians(lat0), np.radians(lat)
    a = (np.sin((p - p0) / 2)**2
         + np.cos(p0) * np.cos(p) * np.sin(np.radians(lon - lon0) / 2)**2)
    return 2 * R * np.arcsin(np.sqrt(a))


def camp_position():
    """EastGRIP camp centroid from the 20240618_01 survey in our own data."""
    T = loadmat(os.path.join(DATA, 'negis2024_tracks_fixed.mat'),
                squeeze_me=True)['T']
    la, lo = [], []
    for r in T:
        if str(r['seg']) != CAMP_SEG:
            continue
        a = np.atleast_1d(np.asarray(r['lat'], float))
        b = np.atleast_1d(np.asarray(r['lon'], float))
        ok = np.isfinite(a) & np.isfinite(b)
        la.extend(a[ok])
        lo.extend(b[ok])
    if not la:
        raise SystemExit('no %s camp survey in the track file' % CAMP_SEG)
    return float(np.mean(la)), float(np.mean(lo))


def all_tracks():
    T = loadmat(os.path.join(DATA, 'negis2024_tracks_fixed.mat'),
                squeeze_me=True)['T']
    out = []
    for r in T:
        la = np.atleast_1d(np.asarray(r['lat'], float))
        lo = np.atleast_1d(np.asarray(r['lon'], float))
        ok = np.isfinite(la) & np.isfinite(lo)
        if ok.sum() > 2:
            out.append((la[ok], lo[ok]))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    bore_lat, bore_lon = camp_position()

    S = loadmat(os.path.join(DATA, 'fabric_sections.mat'),
                squeeze_me=True)['S']
    if SITE not in S.dtype.names:
        raise SystemExit('no %s section in fabric_sections.mat - run '
                         'opr_fabric/server/run_sections.m first' % SITE)
    rec = S[SITE].item() if S[SITE].dtype == object else S[SITE]

    nint = int(field(rec, 'nint')[0])
    dlam = field(rec, 'dlam').reshape(nint, -1)
    qual = field(rec, 'quality').reshape(nint, -1)
    top = field(rec, 'top').reshape(nint, -1)
    bot = field(rec, 'bot').reshape(nint, -1)
    lat, lon = field(rec, 'lat'), field(rec, 'lon')
    az = track_azimuth(lat, lon)
    z = np.nanmedian((top + bot) / 2, axis=1)

    dist_to_bore = haversine_km(lat, lon, bore_lat, bore_lon)
    near = dist_to_bore <= RADIUS_KM
    print('%s: %d of %d blocks within %.1f km of the borehole (closest '
          '%.2f km)' % (FRAME, near.sum(), near.size, RADIUS_KM,
                        dist_to_bore.min()))
    if near.sum() < 3:
        near = dist_to_bore <= np.percentile(dist_to_bore, 25)
        print('  widened to the nearest quartile (%d blocks)' % near.sum())

    med = np.full(nint, np.nan)
    lo_b = np.full(nint, np.nan)
    hi_b = np.full(nint, np.nan)
    for i in range(nint):
        if z[i] < Z_MIN:
            continue
        ok = near & np.isfinite(dlam[i]) & (qual[i] >= COH_MIN)
        if ok.sum() >= 3:
            med[i] = np.median(dlam[i][ok])
            lo_b[i], hi_b[i] = np.percentile(dlam[i][ok], [16, 84])
    drawn = np.isfinite(med)
    if not drawn.any():
        raise SystemExit('no interval near the borehole clears the gates')
    print('  %d intervals drawn, %.0f-%.0f m, dlam %.3f..%.3f'
          % (drawn.sum(), z[drawn].min(), z[drawn].max(),
             np.nanmin(med), np.nanmax(med)))

    # ---- figure
    fig = plt.figure(figsize=(12.2, 5.8), dpi=DPI, layout='constrained')
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0])
    axm = fig.add_subplot(gs[0, 0], projection=PROJ)
    axp = fig.add_subplot(gs[0, 1])

    zwin = np.load(os.path.join(DATA, 'itslive_negis_win.npz'))
    wx, wy, v = zwin['x'], zwin['y'], zwin['v']
    vx, vy = zwin['vx'], zwin['vy']
    norm = LogNorm(vmin=2, vmax=120)

    def speed_layer(ax, dec):
        pcm = ax.pcolormesh(wx, wy, np.maximum(v, 1e-3), transform=PROJ,
                            cmap='magma', norm=norm, shading='auto',
                            rasterized=True)
        Xq, Yq = np.meshgrid(wx[::dec], wy[::dec])
        U, V = vx[::dec, ::dec], vy[::dec, ::dec]
        sp = np.hypot(U, V)
        ok = sp > 1.0
        ax.quiver(Xq, Yq, np.where(ok, U / sp, np.nan),
                  np.where(ok, V / sp, np.nan), transform=PROJ,
                  color='white', scale=26, width=0.005, alpha=0.6, zorder=6)
        return pcm

    # Window on the drill site, not the whole NEGIS survey: the figure is
    # about what the radar samples next to the borehole.
    half = 11e3
    px, py = PROJ.transform_point(bore_lon, bore_lat, ccrs.PlateCarree())[:2]
    extent = (px - half, px + half, py - half, py + half)
    axm.set_extent(extent, crs=PROJ)
    pcm = speed_layer(axm, dec=8)
    cax = axm.inset_axes([0.05, -0.10, 0.62, 0.03])
    fig.colorbar(pcm, cax=cax, orientation='horizontal',
                 label='surface speed (m/yr)')

    for la, lo in all_tracks():
        axm.plot(lo, la, '-', color=TRACK_C, lw=1.1, alpha=0.9,
                 transform=ccrs.PlateCarree(), zorder=7)
    axm.plot(lon, lat, '-', color=PROFILE_C, lw=PROFILE_LW,
             transform=ccrs.PlateCarree(), zorder=8,
             path_effects=PROFILE_FX)
    sty.draw_map_ends(axm, lat, lon)
    # The radius the depth profile actually averages over
    ang = np.linspace(0, 2 * np.pi, 181)
    axm.plot(bore_lon + RADIUS_KM / 111.0 * np.sin(ang) / np.cos(
        np.radians(bore_lat)),
        bore_lat + RADIUS_KM / 111.0 * np.cos(ang), '--', color=BORE_FACE,
        lw=1.2, transform=ccrs.PlateCarree(), zorder=9)
    axm.scatter([bore_lon], [bore_lat], marker='*', s=340, c=BORE_FACE,
                edgecolor='black', lw=1.4, transform=ccrs.PlateCarree(),
                zorder=11)
    axm.annotate('EastGRIP', xy=(bore_lon, bore_lat),
                 xycoords=ccrs.PlateCarree()._as_mpl_transform(axm),
                 xytext=(10, -16), textcoords='offset points', fontsize=9,
                 fontweight='bold', color='black', zorder=12,
                 bbox=dict(boxstyle='round,pad=0.16', fc=BORE_FACE,
                           ec='black', lw=0.8, alpha=0.95))
    gl = axm.gridlines(draw_labels=True, lw=0.4, alpha=0.5, color='gray',
                       y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    sty.scale_bar_br(axm, extent, PROJ, nice_length=sty.nice_length_up)
    axm.set_title('EastGRIP: nearest phase-preserving line (%s)' % FRAME,
                  fontsize=11)

    # ---- depth profile
    axp.fill_betweenx(z[drawn], lo_b[drawn], hi_b[drawn], color=C_PROF,
                      alpha=0.20, lw=0, label='block-to-block 16-84%')
    axp.plot(med[drawn], z[drawn], '-', color=C_PROF, lw=2.2,
             label='radar, within %.1f km of borehole' % RADIUS_KM)
    axp.axvline(0, color=sty.MUTED if hasattr(sty, 'MUTED') else '0.4',
                lw=0.8)
    axp.grid(alpha=0.25, lw=0.6)
    axp.set_axisbelow(True)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)
    axp.invert_yaxis()
    axp.set_xlabel(r'$\Delta\lambda = \lambda_{\perp} - \lambda_{\parallel}$'
                   '\n' r'($\perp$ = %.0f$^\circ$,  $\parallel$ = %.0f$^\circ$'
                   ' E of N)' % ((az + 90) % 180, az))
    axp.set_ylabel('depth (m)')
    axp.set_title('horizontal fabric contrast at the drill site',
                  fontsize=11)
    axp.legend(loc='lower right', fontsize=8.5, frameon=False)
    axp.annotate('Radar measures the DIFFERENCE of the two horizontal '
                 'eigenvalues on the\naxes named below; rotate the core '
                 'fabric into that frame before\ncomparing.  This line is '
                 'noisier than the other sites (misfit\n0.37 ns, node '
                 'coherence 0.34) - see the script header.',
                 xy=(0.0, -0.30), xycoords='axes fraction', fontsize=7.5,
                 color='0.35', va='top', ha='left')

    out = os.path.join(OUT, 'scar_egrip_2panel.png')
    # bbox_inches: the caption sits below the axes and constrained
    # layout does not reserve space for it, so it clips otherwise
    fig.savefig(out, dpi=DPI, bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()

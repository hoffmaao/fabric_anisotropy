"""EastGRIP: interferogram plus fabric strength and orientation from azimuth.

The traverse radiates from EastGRIP camp, so several phase-preserving
lines pass within a kilometre of the borehole at different headings. A
single line only measures the horizontal fabric contrast projected onto
its own axes; the set measures the whole horizontal ellipse.

Every quantity here is taken from the FRINGE RATE, not from the layer
stripping. The rate per unit two-way traveltime is fc*d(dtau)/dt, which
the forward model makes a function of dlam alone - 4.021 cycles/us per
unit dlam at 750 MHz, measured by calling ptt.twttDifference with a
uniform dlam (opr_fabric/server/fringe_check.m). That matters here: the
EastGRIP-line inversion is the weakest of the four sites (0.37 ns misfit,
and its sign flips near 560 m where the raw fringe rate does not), so
going through it would propagate that instability into the answer.

For a horizontal fabric with principal azimuth theta and contrast
P = lam_max - lam_min, a line whose TRACK azimuth is alpha measures

    dlam(alpha) = lam_perp - lam_par = -P * cos(2*(alpha - theta))

which is linear in (P cos 2theta, P sin 2theta), so the fit is a
two-parameter least squares per depth band, weighted by coherence.

The 90 degree branch of theta depends on the sign convention linking the
measured fringe rate to dlam, which the radar cannot settle on its own.
It is resolved against the EastGRIP core: the core's smallest eigenvalue
collapses to ~0.01 below 450 m, which under along-flow extension is the
flow direction, so the LARGEST horizontal eigenvalue must lie cross-flow
near 123 deg. The figure prints both branches so the choice is visible.

Usage: python egrip_azimuthal.py <out_dir>
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import h5py                              # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402
from scipy.io import loadmat             # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import DATA, hsv_plot_coherence      # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
# The EastGRIP core table is data, not code: Weikusat et al. (2022),
# PANGAEA.949248, CC-BY-4.0, mirrored beside the other staged inputs.
CORE = os.path.join(DATA, 'egrip_fabric_949248.tab')

EG_LAT, EG_LON = 75.6294, -35.9672      # camp survey centroid (20240618_01)
FLOW_AZ = 32.7                          # ITS_LIVE at EGRIP, deg E of N
K_CONV = 4.021                          # cycles/us per unit dlam at 750 MHz
C_ICE = 1.68e8
RADIUS_KM = 2.0
LAG = 8                                 # fast-time lag, in 4x-looked samples
SHOW = '20240619_01_001'                # the line drawn in the left panel

FRAMES = ['20240628_01_001', '20240626_03_005', '20240626_01_001',
          '20240619_01_001', '20240621_01_001', '20240620_01_001',
          '20240621_01_010', '20240619_01_005', '20240620_02_003']

BANDS = [(200, 350), (350, 500), (500, 650), (650, 800), (800, 950),
         (950, 1100), (1100, 1250)]

AGREE_TOL = 0.35     # fractional agreement the two rate estimators must reach
P_BOUND = 2.0 / 3    # |lam_max - lam_min| cannot exceed 1 - lam_z

INK, MUTED = '#0b0b0b', '#52514e'
C_FIT, C_PT = '#2a78d6', '#eb6834'


def haversine_km(lat, lon, lat0, lon0):
    R = 6371.0
    p0, p = np.radians(lat0), np.radians(lat)
    a = (np.sin((p - p0) / 2)**2
         + np.cos(p0) * np.cos(p) * np.sin(np.radians(lon - lon0) / 2)**2)
    return 2 * R * np.arcsin(np.sqrt(a))


def load_line(tag):
    fn = os.path.join(DATA, 'negis_ifg_%s.mat' % tag)
    if not os.path.exists(fn):
        return None
    with h5py.File(fn) as f:
        z = f['interferogram_mlook'][()]
        ifg = (z['real'] + 1j * z['imag']).T.astype(complex)
        ph = f['power_hh'][()].T.astype(float)
        pv = f['power_vv'][()].T.astype(float)
        t = f['Time'][()].ravel()
        surf = f['Surface_ml'][()].ravel()
        lat = f['Latitude_ml'][()].ravel()
        lon = f['Longitude_ml'][()].ravel()
    with np.errstate(invalid='ignore'):
        coh = np.abs(ifg) / np.sqrt(ph * pv)
    return dict(ifg=ifg, coh=coh, t=t, surf=surf, lat=lat, lon=lon, tag=tag)


def near_mask(d):
    return haversine_km(d['lat'], d['lon'], EG_LAT, EG_LON) <= RADIUS_KM


def track_azimuth(lat, lon):
    """Principal axis of the track, deg E of N, over 0-180."""
    y = lat - lat.mean()
    x = (lon - lon.mean()) * np.cos(np.radians(lat.mean()))
    pts = np.c_[x, y]
    pts = pts - pts.mean(0)
    _, _, vt = np.linalg.svd(pts, full_matrices=False)
    return np.degrees(np.arctan2(vt[0][0], vt[0][1])) % 180


def band_rates(d, near):
    """Fringe rate and coherence per depth band, over the near columns.

    Returns TWO independent rate estimates per band, plus coherence.

    The first is the angle of a single-lag coherent sum,
    angle(sum I(t+D) conj(I(t))); the second is the slope of the UNWRAPPED
    phase of the coherence-weighted stacked phasor. Neither is reliable
    alone once coherence falls, and they fail in OPPOSITE directions: the
    lag form is biased toward zero (the same mechanism that suppresses
    delta-k at depth, here under-reading the 500-950 m contrast by about a
    third against a direct cycle count), while unwrapping slips whole
    cycles and runs away from zero - below 650 m it returned P above 1,
    which is impossible for a difference of eigenvalues.

    Because they fail in opposite directions, their AGREEMENT is a
    reliability test that neither provides on its own, and the caller
    accepts a band only where they agree.
    """
    dt = float(np.median(np.diff(d['t'])))
    tb = (d['t'] - np.nanmedian(d['surf'])) * 1e6
    w = np.nan_to_num(d['coh'][:, near])
    stack = np.nansum(d['ifg'][:, near] * w, axis=1)
    uw = np.unwrap(np.angle(stack))
    X = d['ifg'][LAG:, near] * np.conj(d['ifg'][:-LAG, near])
    out = []
    for z0, z1 in BANDS:
        t0, t1 = 2 * z0 / C_ICE * 1e6, 2 * z1 / C_ICE * 1e6
        m = (tb >= t0) & (tb < t1)
        ml = (tb[:-LAG] >= t0) & (tb[:-LAG] < t1)
        if m.sum() < 8 or ml.sum() < 5:
            out.append((np.nan, np.nan, np.nan))
            continue
        r_uw = np.polyfit(tb[m], uw[m], 1)[0] / (2 * np.pi)
        r_lag = np.angle(np.nansum(X[ml, :])) / (2 * np.pi * LAG * dt) * 1e-6
        c = float(np.nanmedian(d['coh'][m][:, near]))
        out.append((r_lag / K_CONV, r_uw / K_CONV, c))
    return np.array(out)


def core_bracket():
    """Weighted eigenvalues from the EGRIP core, as a horizontal bracket."""
    with open(CORE) as f:
        lines = f.read().split('\n')
    h = next(i for i, l in enumerate(lines) if l.startswith('Bag\t'))
    cols = lines[h].split('\t')
    iz = cols.index('Depth ice/snow [m]')
    ie = [cols.index('EVA%d (Weighted (statistic))' % k) for k in (1, 2, 3)]
    Z, E = [], []
    for l in lines[h + 1:]:
        p = l.split('\t')
        if len(p) <= max(ie):
            continue
        try:
            Z.append(float(p[iz]))
            E.append([float(p[i]) for i in ie])
        except ValueError:
            pass
    Z, E = np.array(Z), np.array(E)
    zc, lo, hi = [], [], []
    for z0, z1 in BANDS:
        m = (Z >= z0) & (Z < z1)
        if m.sum() < 3:
            zc.append(np.nan); lo.append(np.nan); hi.append(np.nan); continue
        zc.append(0.5 * (z0 + z1))
        lo.append(np.median(E[m, 1] - E[m, 0]))   # e2-e1, largest vertical
        hi.append(np.median(E[m, 2] - E[m, 0]))   # e3-e1, largest horizontal
    return np.array(zc), np.array(lo), np.array(hi)


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []
    for tag in FRAMES:
        d = load_line(tag)
        if d is None:
            print('missing extract for', tag)
            continue
        nm = near_mask(d)
        if nm.sum() < 8:
            print('%s: only %d columns within %.1f km; skipped'
                  % (tag, nm.sum(), RADIUS_KM))
            continue
        d['az'] = track_azimuth(d['lat'][nm], d['lon'][nm])
        d['rows'] = band_rates(d, nm)
        d['nnear'] = int(nm.sum())
        lines.append(d)
        print('%-18s azimuth %6.1f deg, %4d near columns' % (tag, d['az'],
                                                             d['nnear']))
    if len(lines) < 3:
        raise SystemExit('need at least 3 lines for an azimuthal solve')

    az = np.array([d['az'] for d in lines])
    print('\nazimuth coverage: %s' % np.array2string(np.sort(az), precision=1))

    # ---- per-band two-parameter fit: dlam = -P cos 2(alpha - theta)
    P = np.full(len(BANDS), np.nan)
    TH = np.full(len(BANDS), np.nan)
    RES = np.full(len(BANDS), np.nan)
    NL = np.zeros(len(BANDS), int)
    for bi in range(len(BANDS)):
        y = np.array([d['rows'][bi, 0] for d in lines])
        y2 = np.array([d['rows'][bi, 1] for d in lines])
        w = np.array([d['rows'][bi, 2] for d in lines])
        denom = np.maximum(np.abs(y), np.abs(y2))
        with np.errstate(invalid='ignore', divide='ignore'):
            agree = np.abs(y - y2) / np.where(denom > 0, denom, np.nan)
        y = 0.5 * (y + y2)
        ok = (np.isfinite(y) & np.isfinite(w) & (w > 0.35)
              & np.isfinite(agree) & (agree < AGREE_TOL))
        if ok.sum() < 3:
            continue
        a = np.radians(az[ok])
        A = np.c_[-np.cos(2 * a), -np.sin(2 * a)]
        W = np.sqrt(w[ok])[:, None]
        sol, *_ = np.linalg.lstsq(A * W, y[ok] * W[:, 0], rcond=None)
        Pb = np.hypot(*sol)
        if Pb > P_BOUND:
            print('  %d-%d m: P = %.2f exceeds the eigenvalue bound; rejected'
                  % (*BANDS[bi], Pb))
            continue
        P[bi] = Pb
        TH[bi] = np.degrees(0.5 * np.arctan2(sol[1], sol[0])) % 180
        RES[bi] = np.sqrt(np.mean((A @ sol - y[ok])**2))
        NL[bi] = ok.sum()

    zc = np.array([0.5 * (a + b) for a, b in BANDS])
    print('\n%9s %5s %8s %9s %9s' % ('depth_m', 'n', 'P', 'theta_deg', 'resid'))
    for bi in range(len(BANDS)):
        print('%9s %5d %8.3f %9.1f %9.3f'
              % ('%d-%d' % BANDS[bi], NL[bi], P[bi], TH[bi], RES[bi]))

    czc, clo, chi = core_bracket()

    # ---- figure
    fig = plt.figure(figsize=(16.0, 5.4), layout='constrained')
    gs = fig.add_gridspec(1, 4, width_ratios=[1.7, 1.05, 0.85, 0.85])
    axi = fig.add_subplot(gs[0, 0])
    axa = fig.add_subplot(gs[0, 1])
    axp = fig.add_subplot(gs[0, 2])
    axt = fig.add_subplot(gs[0, 3], sharey=axp)

    # (a) the interferogram of the line we use
    d = next(x for x in lines if x['tag'] == SHOW)
    nm = near_mask(d)
    dist = np.concatenate([[0], np.cumsum(haversine_km(
        d['lat'][1:], d['lon'][1:], d['lat'][:-1], d['lon'][:-1]))])
    tb = (d['t'] - np.nanmedian(d['surf'])) * 1e6
    zz = tb * C_ICE / 2 * 1e-6
    st = max(1, d['ifg'].shape[0] // 2000)
    sx = max(1, d['ifg'].shape[1] // 1400)
    rgb = hsv_plot_coherence(np.angle(d['ifg'][::st, ::sx]),
                             d['coh'][::st, ::sx])
    axi.imshow(rgb, aspect='auto', interpolation='nearest',
               extent=[dist[0], dist[-1], zz[-1], zz[0]])
    axi.set_ylim(1300, 0)
    axi.set_xlabel('distance along profile (km)', color=INK)
    axi.set_ylabel('depth (m)', color=INK)
    axi.set_title('HH-VV interferogram, %s  (azimuth %.0f$^\\circ$)'
                  % (SHOW, d['az']), fontsize=11, color=INK)
    # where the borehole sits along this line
    kb = int(np.argmin(haversine_km(d['lat'], d['lon'], EG_LAT, EG_LON)))
    axi.axvline(dist[kb], color='#ffd400', lw=1.6, ls='--')
    # The borehole sits at one end of every line (the traverse drives out
    # from camp), so anchor the label inside the panel rather than centring
    # it on the marker, where it would hang off the axis.
    side = 'left' if dist[kb] < 0.5 * dist[-1] else 'right'
    axi.annotate('EastGRIP', xy=(dist[kb], 55),
                 xytext=(8 if side == 'left' else -8, 0),
                 textcoords='offset points', fontsize=9, color='black',
                 fontweight='bold', va='top',
                 ha='left' if side == 'left' else 'right',
                 bbox=dict(boxstyle='round,pad=0.16', fc='#ffd400',
                           ec='black', lw=0.8))
    sm = plt.cm.ScalarMappable(cmap='hsv', norm=plt.Normalize(-np.pi, np.pi))
    cb = fig.colorbar(sm, ax=axi, pad=0.02)
    cb.set_label('phase change (rad)', color=INK)
    cb.set_ticks([-np.pi, 0, np.pi])
    cb.set_ticklabels([r'$-\pi$', '0', r'$\pi$'])

    # (b) the azimuthal fit, for the bands that solved
    show_b = [bi for bi in range(len(BANDS)) if np.isfinite(P[bi])][:4]
    phi = np.linspace(0, 180, 181)
    for j, bi in enumerate(show_b):
        off = j * 0.30
        y = 0.5 * np.array([d['rows'][bi, 0] + d['rows'][bi, 1]
                            for d in lines])
        w = np.array([d['rows'][bi, 2] for d in lines])
        ok = np.isfinite(y) & (w > 0.35)
        col = plt.get_cmap('viridis')(j / max(len(show_b) - 1, 1))
        axa.plot(phi, -P[bi] * np.cos(2 * np.radians(phi - TH[bi])) + off,
                 '-', color=col, lw=1.8)
        axa.scatter(az[ok], y[ok] + off, s=34, facecolor=col,
                    edgecolor='black', lw=0.7, zorder=5)
        axa.annotate('%d-%d m' % BANDS[bi], xy=(182, off), fontsize=8,
                     color=col, va='center', annotation_clip=False)
        axa.axhline(off, color='0.85', lw=0.6, zorder=0)
    axa.axvline(FLOW_AZ, color=MUTED, lw=1.0, ls=':')
    axa.annotate('flow %.0f$^\\circ$' % FLOW_AZ, xy=(FLOW_AZ + 2, axa.get_ylim()[1]),
                 fontsize=7.5, color=MUTED, va='top')
    axa.set_xlim(0, 180)
    axa.set_xticks([0, 45, 90, 135, 180])
    axa.set_xlabel('line azimuth (deg E of N)', color=INK)
    axa.set_ylabel(r'$\Delta\lambda$ from fringe rate  (offset per band)',
                   color=INK)
    axa.set_title('azimuthal fit', fontsize=11, color=INK)
    axa.grid(alpha=0.25, lw=0.6)
    axa.set_axisbelow(True)
    for s in ('top', 'right'):
        axa.spines[s].set_visible(False)

    # (c) strength with the core bracket
    ok = np.isfinite(P)
    axp.fill_betweenx(czc, clo, chi, color='0.75', alpha=0.55, lw=0,
                      label='EGRIP core\n($e_2-e_1$ to $e_3-e_1$)')
    axp.plot(P[ok], zc[ok], 'o-', color=C_FIT, lw=2.0, ms=5,
             label='radar (this solve)')
    axp.set_xlabel('horizontal contrast $P$', color=INK)
    axp.set_ylabel('depth (m)', color=INK)
    axp.set_title('fabric strength', fontsize=11, color=INK)
    axp.set_xlim(left=0)
    axp.invert_yaxis()
    axp.grid(alpha=0.25, lw=0.6)
    axp.set_axisbelow(True)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)
    axp.legend(loc='lower right', fontsize=7.5, frameon=False)

    # (d) orientation, with both branches and the expected cross-flow axis
    axt.plot(TH[ok], zc[ok], 'o-', color=C_PT, lw=2.0, ms=5, label='solved')
    axt.plot((TH[ok] + 90) % 180, zc[ok], 'o', color=C_PT, ms=4, alpha=0.30,
             label='other 90$^\\circ$ branch')
    axt.axvline(FLOW_AZ, color=MUTED, lw=1.0, ls=':')
    axt.axvline((FLOW_AZ + 90) % 180, color='0.25', lw=1.2, ls='--')
    axt.annotate('cross-flow\n%.0f$^\\circ$' % ((FLOW_AZ + 90) % 180),
                 xy=((FLOW_AZ + 90) % 180 + 3, zc[ok].max()), fontsize=7.5,
                 color='0.25', va='bottom')
    axt.annotate('flow', xy=(FLOW_AZ + 3, zc[ok].max()), fontsize=7.5,
                 color=MUTED, va='bottom')
    axt.set_xlim(0, 180)
    axt.set_xticks([0, 45, 90, 135, 180])
    axt.set_xlabel(r'$\theta$, azimuth of $\lambda_{max}$ (deg E of N)',
                   color=INK)
    axt.set_title('fabric orientation', fontsize=11, color=INK)
    axt.grid(alpha=0.25, lw=0.6)
    axt.set_axisbelow(True)
    for s in ('top', 'right'):
        axt.spines[s].set_visible(False)
    axt.legend(loc='lower right', fontsize=7.5, frameon=False)
    plt.setp(axt.get_yticklabels(), visible=False)

    fig.suptitle('EastGRIP: horizontal fabric from %d phase-preserving lines '
                 'within %.1f km of the borehole' % (len(lines), RADIUS_KM),
                 fontsize=12.5, color=INK)
    out = os.path.join(OUT, 'scar_egrip_azimuthal.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('\nwrote', out)


if __name__ == '__main__':
    main()

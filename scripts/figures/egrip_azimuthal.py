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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style as sty                 # noqa: E402
# track_azimuth is re-exported deliberately: the sibling EastGRIP scripts
# import it from here alongside FRAMES/BANDS, and they must use the same
# NaN-guarded axis fit the rest of the deck does.
from scar_style import DATA, cumdist_km, track_azimuth   # noqa: E402,F401

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
COH_MIN = 0.35       # band coherence below which a line does not enter the fit
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


def band_estimates(lines, bi, col=None):
    """Per-line value, weight and acceptance mask for one depth band.

    The acceptance rule lives here and nowhere else, because it is applied
    in three places - the solve, the scatter drawn over the fitted curve,
    and the two bracketing solves in egrip_core_compare.py - and a copy
    that drifts would draw points the fit never saw.

    A line is accepted only where the band coherence clears COH_MIN AND
    the two independent rate estimators agree to AGREE_TOL. `col` picks
    which estimator is fitted: None averages them (the headline solve),
    0 the lag form, 1 the unwrapped slope.
    """
    est = np.array([[d['rows'][bi, 0] for d in lines],
                    [d['rows'][bi, 1] for d in lines]])
    w = np.array([d['rows'][bi, 2] for d in lines])
    denom = np.maximum(np.abs(est[0]), np.abs(est[1]))
    with np.errstate(invalid='ignore', divide='ignore'):
        agree = np.abs(est[0] - est[1]) / np.where(denom > 0, denom, np.nan)
    y = est.mean(0) if col is None else est[col]
    ok = (np.isfinite(y) & np.isfinite(w) & (w > COH_MIN)
          & np.isfinite(agree) & (agree < AGREE_TOL))
    return y, w, ok


def solve_band(az, y, w, ok):
    """Weighted fit of dlam = -P cos 2(alpha - theta) over the accepted lines.

    Returns (P, theta_deg, rms residual), or None if fewer than three
    lines survive - two free parameters need three points to be more than
    an interpolation. The P_BOUND test is left to the caller so it can
    report the rejected value.
    """
    if ok.sum() < 3:
        return None
    a = np.radians(az[ok])
    A = np.c_[-np.cos(2 * a), -np.sin(2 * a)]
    W = np.sqrt(w[ok])[:, None]
    sol, *_ = np.linalg.lstsq(A * W, y[ok] * W[:, 0], rcond=None)
    return (float(np.hypot(*sol)),
            float(np.degrees(0.5 * np.arctan2(sol[1], sol[0])) % 180),
            float(np.sqrt(np.mean((A @ sol - y[ok])**2))))


def core_table():
    """Depth and the three weighted eigenvalues from the EGRIP core table.

    A row contributes only if the depth AND all three eigenvalues parse.
    Appending the depth first would leave Z one entry longer than E on any
    row with a blank eigenvalue field, and every caller builds its mask on
    Z and applies it to E, so the two must not desynchronise.
    """
    with open(CORE) as f:
        rows = f.read().split('\n')
    h = next(i for i, l in enumerate(rows) if l.startswith('Bag\t'))
    cols = rows[h].split('\t')
    iz = cols.index('Depth ice/snow [m]')
    ie = [cols.index('EVA%d (Weighted (statistic))' % k) for k in (1, 2, 3)]
    Z, E = [], []
    for l in rows[h + 1:]:
        p = l.split('\t')
        if len(p) <= max(ie + [iz]):
            continue
        try:
            z, e = float(p[iz]), [float(p[i]) for i in ie]
        except ValueError:
            continue
        Z.append(z)
        E.append(e)
    Z, E = np.array(Z), np.array(E)
    o = np.argsort(Z)
    return Z[o], E[o]


def core_bracket():
    """Weighted eigenvalues from the EGRIP core, as a horizontal bracket."""
    Z, E = core_table()
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
        y, w, ok = band_estimates(lines, bi)
        fit = solve_band(az, y, w, ok)
        if fit is None:
            continue
        Pb, th, res = fit
        if Pb > P_BOUND:
            print('  %d-%d m: P = %.2f exceeds the eigenvalue bound; rejected'
                  % (*BANDS[bi], Pb))
            continue
        P[bi], TH[bi], RES[bi], NL[bi] = Pb, th, res, ok.sum()

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
    d = next((x for x in lines if x['tag'] == SHOW), None)
    if d is None:
        raise SystemExit('%s did not load, so the left panel has nothing to '
                         'draw; re-run the extract for it or point SHOW at '
                         'one of: %s'
                         % (SHOW, ', '.join(x['tag'] for x in lines)))
    dist = cumdist_km(d['lat'], d['lon'])
    tb = (d['t'] - np.nanmedian(d['surf'])) * 1e6
    zz = tb * C_ICE / 2 * 1e-6
    sty.ifg_panel(axi, fig, dist, zz, np.angle(d['ifg']), d['coh'],
                  ylabel='depth (m)')
    axi.set_ylim(1300, 0)
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

    # (b) the azimuthal fit, for the bands that solved. Filled markers are
    # the lines the fit used; hollow ones are lines the agreement or
    # coherence test rejected, drawn so the curve is not seen to miss
    # points it was never fitted to.
    show_b = [bi for bi in range(len(BANDS)) if np.isfinite(P[bi])][:4]
    phi = np.linspace(0, 180, 181)
    for j, bi in enumerate(show_b):
        off = j * 0.30
        y, _, ok = band_estimates(lines, bi)
        drop = np.isfinite(y) & ~ok
        col = plt.get_cmap('viridis')(j / max(len(show_b) - 1, 1))
        axa.plot(phi, -P[bi] * np.cos(2 * np.radians(phi - TH[bi])) + off,
                 '-', color=col, lw=1.8)
        axa.scatter(az[ok], y[ok] + off, s=34, facecolor=col,
                    edgecolor='black', lw=0.7, zorder=5)
        axa.scatter(az[drop], y[drop] + off, s=30, facecolor='none',
                    edgecolor=col, lw=0.9, alpha=0.55, zorder=4)
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
    # The agreement gate is allowed to reject every band, so anchor the
    # axis labels on the band grid rather than on the solved depths, which
    # can be empty.
    z_lbl = zc[ok].max() if ok.any() else zc.max()
    axt.annotate('cross-flow\n%.0f$^\\circ$' % ((FLOW_AZ + 90) % 180),
                 xy=((FLOW_AZ + 90) % 180 + 3, z_lbl), fontsize=7.5,
                 color='0.25', va='bottom')
    axt.annotate('flow', xy=(FLOW_AZ + 3, z_lbl), fontsize=7.5,
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

"""EastGRIP: our inversion against Zeising et al. (2023), on the core.

Zeising, O., Gerber, T. A., Eisen, O., Ershadi, M. R., Stoll, N.,
Weikusat, I. and Humbert, A. (2023), "Improved estimation of the bulk ice
crystal fabric asymmetry from polarimetric phase co-registration", The
Cryosphere 17, 1097-1105, reaches an RMS of 0.03 against the EastGRIP
core. This figure is the analogue of their Fig. 3b, drawn for all three
methods on the SAME lines so the numbers are comparable rather than
merely similar.

WHAT IS HELD FIXED. Every curve uses the same nine borehole-proximal
lines, the same near-borehole aperture, the same leading-edge surface
pick, the same depth bands and the same azimuthal solve. They differ only
in how the delay becomes a fabric estimate:

  Zeising      per-segment coarse lag + coherence phase, gated at
               coherence 0.65, smoothed over 100 m, differentiated over a
               200 m window - a LOCAL, pointwise estimator
               (opr_fabric/server/egrip_zeising.m)
  fringe rate  our per-band fringe rate straight through the forward
               model, with no layer stripping (egrip_azimuthal.py)
  joint        our full path: SNAPHU-unwrapped dtau into ptt.invertBlocks /
               invertHorizontalFabricJoint, a regularised solve that
               additionally fits the returned POWER. This is the one the
               comparison is about (run_negis_fabric.m, inversion =
               'joint'; gathered by egrip_collect_inversion.m)

The conversions to eigenvalues agree across methods: their
Deltaeps' = 0.034 implies 6.39e-11 s/m per unit dlam, ours 6.35e-11 from
an independent calibration of ptt.twttDifference. So a difference in the
answer is a difference of estimator, not of constants.

WHAT IT SHOWS. Over the three bands all three methods solve, RMS against
the core is 0.282 for the joint inversion, 0.280 for the bare fringe rate
and 0.356 for Zeising's estimator. The joint inversion is the one being
tested and it now edges the reference method, at parity with the far
simpler fringe rate.

That took fixing the dtau source, not the inversion. On delta-k the same
joint inversion scored 0.411 - WORSE than either - returning about a
third of the fringe rate's amplitude above 650 m and flipping sign below
it, while agreeing with the fringe rate in sign at every azimuth above
that. It was the delta-k stage-A suppression diagnosed on Ridge A (15 MHz
sub-bands SNR-starved at depth, so smoothing noise-dominated phasors
biases the recovered angle toward zero) arriving through dtau. Moving
NEGIS/EastGRIP onto Goldstein-filtered SNAPHU phase roughly tripled the
recovered amplitude over 500-950 m and removed the sign flip down to
950 m. The bare fringe rate never went through delta-k, so it is
unchanged; Zeising's estimator escapes the same way, its coarse lag
being an amplitude cross-correlation.

Two honest caveats. The per-frame dtau misfit ROSE when SNAPHU replaced
delta-k (0.21-1.16 ns against 0.19-0.43), which is the expected
direction: suppressed dtau is smooth and easy to fit and wrong, while
unwrapped phase is larger and noisier and right - so misfit alone would
have preferred the worse answer, exactly as run_sections.m warns about
rms endorsing an overfit. And all three still under-read a core sitting
at 0.35-0.55, and only three of seven bands solve at all, which is a
property of this two-azimuth acquisition rather than of any estimator.

WHICH EIGENVALUE PAIR THE CORE CONTRIBUTES. dlam is the difference of the
two HORIZONTAL eigenvalues, so which core pair to use depends on which
eigenvector is vertical, and at EastGRIP that changes with depth. This
follows their choice exactly: e2-e1 above 250 m, where the largest
eigenvalue is still vertical, and e3-e1 below it, where NEGIS shear has
rotated the largest into the horizontal plane and the middle one is
vertical. Using either pair throughout would misstate the core.

Usage: python egrip_method_compare.py <out_dir>
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import h5py                              # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import DATA, INK, MUTED  # noqa: E402
from scar_style import FIGS  # noqa: E402
from egrip_azimuthal import (BANDS, EG_LAT, EG_LON, FLOW_AZ,  # noqa: E402
                             P_BOUND, RADIUS_KM, band_estimates, core_table,
                             haversine_km, solve_band, solve_bands,
                             stage_lines)

OUT = sys.argv[1] if len(sys.argv) > 1 else FIGS
ZFN = os.path.join(DATA, 'egrip_zeising.mat')
IFN = os.path.join(DATA, 'egrip_inversion.mat')

C_OURS, C_ZEIS, C_CORE = '#2a78d6', '#eb6834', '#3d3d3d'
C_INV = '#1baf7a'
Z_SWITCH = 250.0        # where the vertical eigenvector changes, per Fig 3b
ZEIS_RMS_PUB = 0.03     # their published RMS against this core
MIN_SEG = 5             # segments a band needs before it enters the solve

# Orientation held fixed in the one-parameter solve. Their ApRES rotates
# through full azimuth at a spot; our nine lines sit at 33 deg and 128-145
# deg, so the two-parameter fit is effectively a two-point measurement and
# the single 33 deg line carries all of the orthogonal leverage. It is
# reliable to 350 m and not below, which is exactly where the free fit
# starts returning P past the eigenvalue bound. Fixing theta at cross-flow
# - independently supported by the core, whose smallest eigenvalue
# collapses to 0.01 along flow - removes the ill-conditioned direction and
# leaves a one-parameter fit the geometry can actually support.
THETA_FIX = (FLOW_AZ + 90) % 180


def solve_fixed(az, y, w, ok, theta):
    """Amplitude of -P cos 2(alpha - theta) at KNOWN theta, weighted.

    One free parameter, so it needs only one azimuth rather than three,
    and it cannot trade strength against orientation - the failure mode
    that puts the free fit past the eigenvalue bound.
    """
    if ok.sum() < 2:
        return np.nan
    c = np.cos(2 * np.radians(az[ok] - theta))
    ww = w[ok]
    den = float(np.sum(ww * c * c))
    if den <= 0:
        return np.nan
    return float(-np.sum(ww * y[ok] * c) / den)


def zeising_lines():
    """Per-line dlam(z) profiles from the ported estimator.

    Returns (segment lengths, list of per-variant {tag: profile} dicts).
    The server script runs several segment lengths because that is the one
    parameter of theirs that cannot be carried over as a length - it sets
    how far the sub-bin phase drifts across a segment, which scales with
    centre frequency. All are returned so the pick is visible.
    """
    if not os.path.exists(ZFN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/egrip_zeising.m on the '
            'server and mirror the result there' % ZFN)
    with h5py.File(ZFN) as f:
        g = f['L']
        segs = np.array(f[g['seg_list'][0, 0]]).ravel()
        out = [{} for _ in segs]
        for i in range(g['tag'].shape[0]):
            tag = ''.join(chr(c) for c in np.array(f[g['tag'][i, 0]]).ravel())
            zc = f[g['z'][i, 0]]
            dc = f[g['dlam'][i, 0]]
            cc = f[g['coh'][i, 0]]
            for vi in range(len(segs)):
                out[vi][tag] = dict(
                    z=np.array(f[zc[vi, 0]]).ravel(),
                    dlam=np.array(f[dc[vi, 0]]).ravel(),
                    coh=np.array(f[cc[vi, 0]]).ravel())
    return segs, out


def zeising_bands(zl, tag):
    """Band-average one line's dlam and coherence onto our depth bands."""
    d = zl.get(tag)
    y = np.full(len(BANDS), np.nan)
    w = np.full(len(BANDS), np.nan)
    if d is None:
        return y, w
    for bi, (z0, z1) in enumerate(BANDS):
        m = (d['z'] >= z0) & (d['z'] < z1) & np.isfinite(d['dlam'])
        if m.sum() < MIN_SEG:
            continue
        y[bi] = np.median(d['dlam'][m])
        w[bi] = np.nanmedian(d['coh'][m])
    return y, w


def inversion_lines():
    """Per-line dlam from our joint (dtau + power) inversion, by tag.

    This is the path the comparison is actually about: ptt.invertBlocks /
    invertHorizontalFabricJoint, which fits the returned power alongside
    dtau under a regulariser, rather than differentiating a pointwise
    delay. Produced by run_negis_fabric.m over the same nine frames and
    gathered by egrip_collect_inversion.m.
    """
    if not os.path.exists(IFN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_negis_fabric.m with the '
            'nine EastGRIP frames, then egrip_collect_inversion.m, and '
            'mirror the result there' % IFN)
    out = {}
    with h5py.File(IFN) as f:
        g = f['R']
        for i in range(g['tag'].shape[0]):
            tag = ''.join(chr(c) for c in np.array(f[g['tag'][i, 0]]).ravel())
            out[tag] = dict(
                dlam=np.array(f[g['dlam'][i, 0]]).T,
                top=np.array(f[g['top'][i, 0]]).T,
                bot=np.array(f[g['bot'][i, 0]]).T,
                quality=np.array(f[g['quality'][i, 0]]).T,
                lat=np.array(f[g['lat'][i, 0]]).ravel(),
                lon=np.array(f[g['lon'][i, 0]]).ravel())
    return out


def inversion_bands(il, tag):
    """Band-average one line's inverted dlam over near-borehole blocks.

    Intervals are ~146 m and the bands 150 m, so an interval is credited
    to a band by OVERLAP rather than by its centre: at this spacing a
    centre test drops roughly one interval per band, and the ones it drops
    are not random - they are whichever intervals happen to straddle a
    band edge, which biases the average toward the interval interiors.
    """
    d = il.get(tag)
    y = np.full(len(BANDS), np.nan)
    w = np.full(len(BANDS), np.nan)
    if d is None:
        return y, w
    near = haversine_km(d['lat'], d['lon'], EG_LAT, EG_LON) <= RADIUS_KM
    if near.sum() == 0:
        return y, w
    dl = d['dlam'][:, near]
    q = d['quality'][:, near]
    top = d['top'][:, near]
    bot = d['bot'][:, near]
    for bi, (z0, z1) in enumerate(BANDS):
        ov = np.minimum(bot, z1) - np.maximum(top, z0)
        m = np.isfinite(dl) & np.isfinite(q) & (ov > 0)
        if m.sum() < 2:
            continue
        ww = np.nan_to_num(q[m]) * ov[m]
        if ww.sum() <= 0:
            continue
        y[bi] = float(np.sum(dl[m] * ww) / np.sum(ww))
        w[bi] = float(np.nanmedian(q[m]))
    return y, w


def fixed_ours(lines):
    """Our per-band P at fixed theta, off the same accepted-line rule."""
    az = np.array([d['az'] for d in lines])
    P = np.full(len(BANDS), np.nan)
    for bi in range(len(BANDS)):
        y, w, ok = band_estimates(lines, bi)
        P[bi] = solve_fixed(az, y, w, ok, THETA_FIX)
    return P


def _stack(lines, bandfn, src):
    Y, W = [], []
    for d in lines:
        y, w = bandfn(src, d['tag'])
        Y.append(y)
        W.append(w)
    return np.array(Y), np.array(W)


def fixed_from(lines, bandfn, src):
    """Per-band P at fixed theta for any per-line band estimate."""
    az = np.array([d['az'] for d in lines])
    Y, W = _stack(lines, bandfn, src)
    P = np.full(len(BANDS), np.nan)
    for bi in range(len(BANDS)):
        y, w = Y[:, bi], W[:, bi]
        ok = np.isfinite(y) & np.isfinite(w)
        P[bi] = solve_fixed(az, y, w, ok, THETA_FIX)
    return P


def fixed_zeising(lines, zl):
    """Their per-band P at fixed theta, off the same band averaging."""
    return fixed_from(lines, zeising_bands, zl)


def solve_zeising(lines, zl):
    """P and theta per band from the Zeising per-line dlam.

    Deliberately the SAME solve_band as our own path: only the per-line
    dlam differs, so anything the fit does to one method it does to the
    other, and the comparison cannot be an artefact of the azimuthal step.
    """
    az = np.array([d['az'] for d in lines])
    Y, W = _stack(lines, zeising_bands, zl)
    n = len(BANDS)
    P = np.full(n, np.nan)
    TH = np.full(n, np.nan)
    NL = np.zeros(n, int)
    for bi in range(n):
        y, w = Y[:, bi], W[:, bi]
        ok = np.isfinite(y) & np.isfinite(w)
        fit = solve_band(az, y, w, ok)
        if fit is None:
            continue
        Pb, th, _ = fit
        if Pb > P_BOUND:
            print('  %d-%d m (Zeising): P = %.2f exceeds the eigenvalue '
                  'bound; rejected' % (*BANDS[bi], Pb))
            continue
        P[bi], TH[bi], NL[bi] = Pb, th, ok.sum()
    return P, TH, NL


def core_reference():
    """Core dlam per band, on Zeising's depth-dependent eigenvalue pair.

    Returns the band-median plus the 16-84 spread, which is the honest
    error bar on 'the core value in this band': the core samples every
    ~0.55 m and the bands are 150 m, so a single median hides real
    structure that neither radar method could be expected to follow.
    """
    Z, E = core_table()
    zc = np.array([0.5 * (a + b) for a, b in BANDS])
    med = np.full(len(BANDS), np.nan)
    lo = np.full(len(BANDS), np.nan)
    hi = np.full(len(BANDS), np.nan)
    for bi, (z0, z1) in enumerate(BANDS):
        m = (Z >= z0) & (Z < z1)
        if m.sum() < 3:
            continue
        pair = (E[m, 1] - E[m, 0]) if 0.5 * (z0 + z1) < Z_SWITCH \
            else (E[m, 2] - E[m, 0])
        med[bi] = np.median(pair)
        lo[bi], hi[bi] = np.percentile(pair, [16, 84])
    return zc, med, lo, hi


def rms(a, b, m):
    return float(np.sqrt(np.mean((a[m] - b[m])**2))) if m.any() else np.nan


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = stage_lines()
    segs, zvars = zeising_lines()

    print('\n--- our solve ---')
    P_ours, TH_ours, _, NL_ours = solve_bands(lines)
    zc, core, clo, chi = core_reference()

    # Run every segment length and keep their BEST. Scoring their method on
    # a parameter our radar forced on it and then reporting that number
    # would be a straw man; giving it its best variant is the only version
    # of "at least as well" worth claiming.
    F_ours = fixed_ours(lines)
    il = inversion_lines()
    F_inv = fixed_from(lines, inversion_bands, il)
    cand = []
    for vi, seg in enumerate(segs):
        print('--- Zeising solve, %.1f m segments ---' % seg)
        Pz, _, _ = solve_zeising(lines, zvars[vi])
        Fz = fixed_zeising(lines, zvars[vi])
        m = np.isfinite(F_ours) & np.isfinite(Fz) & np.isfinite(core)
        cand.append((rms(Fz, core, m), int(m.sum()), seg, Pz, Fz))
        print('    free-theta bands %d; fixed-theta bands %d, RMS %.3f'
              % (np.isfinite(Pz).sum(), m.sum(), cand[-1][0]))
    nmax = max(c[1] for c in cand)
    best = min((c for c in cand if c[1] == nmax), key=lambda c: c[0])
    _, _, seg_best, P_zeis, F_zeis = best
    print('\nusing their %.1f m variant (best of %s)'
          % (seg_best, ', '.join('%.1f' % s for s in segs)))

    # Score both on the SAME bands: a method that solved an extra band the
    # other could not would otherwise be compared on a different problem.
    both = (np.isfinite(F_ours) & np.isfinite(F_zeis) & np.isfinite(F_inv)
            & np.isfinite(core))
    r_ours = rms(F_ours, core, both)
    r_zeis = rms(F_zeis, core, both)
    r_inv = rms(F_inv, core, both)

    print('\n%9s %7s | %8s %8s %8s'
          % ('depth_m', 'core', 'joint', 'fringe', 'Zeising'))
    for bi in range(len(BANDS)):
        print('%9s %7.3f | %8.3f %8.3f %8.3f'
              % ('%d-%d' % BANDS[bi], core[bi], F_inv[bi], F_ours[bi],
                 F_zeis[bi]))
    print('\nOne-parameter solve at theta = %.0f deg (cross-flow); %d bands '
          'all three solved.' % (THETA_FIX, both.sum()))
    print('RMS against the EastGRIP core:')
    print('  ours, joint inversion (dtau + power)   %.3f' % r_inv)
    print('  ours, fringe rate only                 %.3f' % r_ours)
    print('  Zeising et al. (2023), %.1f m segments  %.3f'
          % (seg_best, r_zeis))
    print('  their published value on ApRES data    %.2f' % ZEIS_RMS_PUB)

    # ---- figure
    fig = plt.figure(figsize=(11.4, 5.6), layout='constrained')
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0])
    axp = fig.add_subplot(gs[0, 0])
    axs = fig.add_subplot(gs[0, 1])

    # (a) their Fig. 3b: dlam against depth, both methods over the core
    axp.fill_betweenx(zc, clo, chi, color='0.80', alpha=0.6, lw=0,
                      label='EGRIP core, 16-84% in band')
    axp.plot(core, zc, 'o-', color=C_CORE, lw=1.6, ms=4.5,
             label='EGRIP core (Weikusat et al. 2022)')
    mo = np.isfinite(F_ours)
    mz = np.isfinite(F_zeis)
    mi = np.isfinite(F_inv)
    axp.plot(F_zeis[mz], zc[mz], 's--', color=C_ZEIS, lw=2.0, ms=5.5,
             label='Zeising et al. (2023), %.1f m segments' % seg_best)
    axp.plot(F_ours[mo], zc[mo], '^:', color=C_OURS, lw=1.8, ms=6,
             label='ours, fringe rate only')
    axp.plot(F_inv[mi], zc[mi], 'o-', color=C_INV, lw=2.6, ms=6.5,
             label='ours, joint inversion (dtau + power)')
    axp.axhline(Z_SWITCH, color=MUTED, lw=0.9, ls=':')
    # Below the switch line and inset from the right spine: at the top of
    # the axis this label ran into the title, and the shallowest band
    # centre is only 25 m below the line, so there is no room above it.
    axp.annotate('core pair switches\n$e_2-e_1$  $\\rightarrow$  $e_3-e_1$',
                 xy=(0.98, Z_SWITCH + 14), xycoords=('axes fraction', 'data'),
                 fontsize=7.5, color=MUTED, ha='right', va='top')
    axp.set_xlabel(r'horizontal asymmetry  $\Delta\lambda$', color=INK)
    axp.set_ylabel('depth (m)', color=INK)
    axp.set_xlim(left=0)
    # Explicit, so the switch line at 250 m has clear space above it rather
    # than sitting on the top spine wherever the shallowest solved band lands
    axp.set_ylim(1300, 150)
    axp.grid(alpha=0.25, lw=0.6)
    axp.set_axisbelow(True)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)
    # Centre left: the radar curves all stop by 650 m and the core sits at
    # dlam ~ 0.55, so the deep low-dlam corner is the only block of the
    # panel no series crosses. At lower right the entries ran through the
    # core line.
    axp.legend(loc='center left', fontsize=8, frameon=False)
    axp.set_title('EastGRIP borehole, %d lines within 2 km  '
                  '($\\theta$ fixed at cross-flow)' % len(lines),
                  fontsize=11, color=INK)

    # (b) against the core one-to-one, which is where "at least as well"
    # is actually read: distance from the diagonal IS the error being
    # summarised by the two RMS numbers.
    lim = [0, max(np.nanmax(core), np.nanmax(F_ours), np.nanmax(F_inv),
                  np.nanmax(F_zeis)) * 1.15]
    axs.plot(lim, lim, '-', color='0.6', lw=1.0, zorder=1)
    axs.plot(core[both], F_zeis[both], 's', color=C_ZEIS, ms=8,
             markeredgecolor='white', markeredgewidth=0.8, zorder=3,
             label='Zeising  RMS %.3f' % r_zeis)
    axs.plot(core[both], F_ours[both], '^', color=C_OURS, ms=8,
             markeredgecolor='white', markeredgewidth=0.8, zorder=4,
             label='ours, fringe rate  RMS %.3f' % r_ours)
    axs.plot(core[both], F_inv[both], 'o', color=C_INV, ms=9,
             markeredgecolor='white', markeredgewidth=0.8, zorder=5,
             label='ours, joint  RMS %.3f' % r_inv)
    for bi in np.flatnonzero(both):
        axs.annotate('%d' % zc[bi], xy=(core[bi], F_inv[bi]),
                     xytext=(6, -3), textcoords='offset points',
                     fontsize=6.5, color=MUTED)
    axs.set_xlim(lim)
    axs.set_ylim(lim)
    axs.set_aspect('equal')
    axs.set_xlabel(r'core  $\Delta\lambda$', color=INK)
    axs.set_ylabel(r'radar  $\Delta\lambda$', color=INK)
    axs.grid(alpha=0.25, lw=0.6)
    axs.set_axisbelow(True)
    for s in ('top', 'right'):
        axs.spines[s].set_visible(False)
    axs.legend(loc='upper left', fontsize=8.5, frameon=False)
    axs.set_title('against the core (labels: band centre, m)',
                  fontsize=11, color=INK)

    out = os.path.join(OUT, 'scar_egrip_method_compare.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('\nwrote', out)


if __name__ == '__main__':
    main()

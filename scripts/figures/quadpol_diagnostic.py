"""Quad-pol azimuth sweep on a real frame, and whether to believe it.

Draws what ptt.quadpolAzimuth actually produced for one frame, because the
method's validity is visible in the sweep itself and not in the numbers it
returns. For a birefringent column the cross-polarized power goes as

    Pxc(z, psi) ~ sin^2(2(psi - theta)) * sin^2(delta(z)/2)

so panel (a) should show TWO things at once: vertical nulls at psi = theta
and theta + 90, holding at every depth, and horizontal nulls wherever
delta passes a multiple of 2*pi. Vertical structure alone means an
orientation with no measurable birefringence; horizontal structure alone
means the reverse; a flat panel means the cross-pol is leakage, since real
cross-polarized power has nowhere to hide from those nulls.

Panel (b) is the leakage test made quantitative: finite antenna isolation
adds a floor that does NOT null, so it shows up as cross/co power that
sits flat with depth instead of oscillating.

Usage: python quadpol_diagnostic.py <out_dir> [tag] [site]
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

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
TAG = sys.argv[2] if len(sys.argv) > 2 else '20250108_02_009'
SITE = sys.argv[3] if len(sys.argv) > 3 else 'ridge_a'
FN = os.path.join(DATA, 'quadpol_%s.mat' % TAG)

# Scale bar for panel (d): the median quad-pol LS contrast over the
# 200-1200 m band. Per site, because drawing one site's contrast across
# another site's frame would be a comparison to nothing, and no bar at all
# is drawn for a site with no entry. Ridge A: median of sec_dlam_ls over the
# 36-frame survey (1.9e6 cells, p25/p75 0.029/0.067). This is the same bar
# ershadi_heading_test.py and quadpol_heading_test.py draw, since the three
# figures are shown in one deck and must not disagree about what the
# contrast is compared against.
#
# It replaces the two-azimuth solve (fabric_map.py) this panel used to draw,
# which is NOT an independent check of the quad-pol path, for two separate
# reasons. First, the LS estimator superseded it: its correlated (theta, P)
# errors bias it high. Second, at Ridge A that solve is close to singular -
# its conditioning goes as |sin 2*da| in the azimuth separation, which
# vanishes at BOTH 0 and 90 deg, and Ridge A's rows at 80.5 deg against the
# N-S connectors at 179.1 deg sit 81.4 deg apart, 9 deg from the 90 deg
# singularity, for 3.4x error amplification. That is how a leg-level bias of
# ~0.007 becomes a factor-of-two disagreement in P. Do not try to fix that
# conditioning here: fabric_map.py's AZ_MIN_SEP gate has the same one-sided
# flaw and is a separate change.
REF_DLAM = {'ridge_a': 0.052}
# Orientation from that same two-azimuth solve, drawn in panel (c) labelled
# with its provenance. Context for the eye, not an independent check.
REF_THETA = 90.0      # median theta, deg E of N
Z_SHOW = (0.0, 1400.0)


def main():
    if not os.path.exists(FN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_quadpol_frame.m on the '
            'server and mirror the result there' % FN)
    with h5py.File(FN) as f:
        z = np.array(f['z']).ravel()
        A = f['A']
        Pxc = np.array(A['Pxc']).T
        psi = np.degrees(np.array(A['psi']).ravel())
        o = f['out']
        theta = np.degrees(np.array(o['theta']).ravel())
        dlam = np.array(o['dlam']).ravel()
        dlam_node = np.array(o['dlam_node']).ravel()
        aniso = np.array(o['aniso']).ravel()
        xr = np.array(f['xr']).ravel()
        cx = np.array(f['c_x'])
        track_az = float(np.array(f['track_az']).ravel()[0])
    cx = np.abs(cx['real'] + 1j * cx['imag']) if cx.dtype.names else np.abs(cx)
    cx = cx.ravel()
    theta_geo = (theta + track_az) % 180.0

    m = (z >= Z_SHOW[0]) & (z <= Z_SHOW[1])
    print('%s: %d depths, aniso median %.3f, cross/co median %.1f dB'
          % (TAG, m.sum(), np.nanmedian(aniso[m]), np.nanmedian(xr[m])))

    fig = plt.figure(figsize=(14.5, 6.2), layout='constrained')
    gs = fig.add_gridspec(1, 4, width_ratios=[1.5, 0.85, 0.85, 0.95])
    axs = fig.add_subplot(gs[0, 0])
    axx = fig.add_subplot(gs[0, 1], sharey=axs)
    axt = fig.add_subplot(gs[0, 2], sharey=axs)
    axd = fig.add_subplot(gs[0, 3], sharey=axs)

    # (a) the sweep itself, normalized per depth so the azimuthal SHAPE is
    # visible at every depth rather than only where the returned power
    # happens to be strong - the shape is the measurement, the level is not
    P = Pxc[m]
    P = P / np.maximum(np.nanmax(P, axis=1, keepdims=True), 1e-30)
    im = axs.pcolormesh(psi, z[m], P, cmap='magma', vmin=0, vmax=1,
                        shading='nearest', rasterized=True)
    # Both principal axes, in the ANTENNA frame, which is the frame the
    # sweep is drawn in. Subtracting track_az from theta_geo - already
    # wrapped to [0,180) - was the first version and put the line off the
    # panel for a track heading near 180 deg.
    th_a = theta[m] % 180.0
    axs.plot(th_a, z[m], '.', color='#67e8f9', ms=1.4, alpha=0.9,
             label=r'solved $\theta$ (both axes)')
    axs.plot((th_a + 90) % 180, z[m], '.', color='#67e8f9', ms=1.4,
             alpha=0.55)
    axs.set_xlim(0, 178)
    axs.set_xticks([0, 45, 90, 135])
    axs.set_xlabel('synthetic antenna azimuth $\\psi$ (deg, antenna frame)',
                   color=INK)
    axs.set_ylabel('depth (m)', color=INK)
    axs.set_title('(a) cross-pol power, azimuth-synthesized\n'
                  '$P_{xc}\\sim\\sin^2 2(\\psi-\\theta)\\,\\sin^2(\\delta/2)$',
                  fontsize=10, color=INK)
    axs.legend(loc='lower right', fontsize=8, frameon=False,
               labelcolor='white')
    cb = fig.colorbar(im, ax=axs, orientation='horizontal', pad=0.02,
                      fraction=0.045)
    cb.set_label('cross-pol power, normalized per depth', fontsize=8)

    # (b) the leakage test
    axx.plot(xr[m], z[m], '-', color='#eb6834', lw=1.0)
    axx.axvline(np.nanmedian(xr[m]), color=MUTED, lw=0.9, ls=':')
    axx.set_xlabel('cross / co power (dB)', color=INK)
    axx.set_title('(b) leakage test', fontsize=10, color=INK)
    axx.grid(alpha=0.25, lw=0.6)

    # (c) orientation
    axt.plot(theta_geo[m], z[m], '.', color='#2a78d6', ms=1.6, alpha=0.6)
    axt.axvline(REF_THETA, color='black', lw=1.4, ls='--')
    axt.axvline((REF_THETA + 90) % 180, color='0.55', lw=1.2, ls=':')
    axt.annotate('two-azimuth\nsolve', xy=(REF_THETA - 4, Z_SHOW[1] * 0.06),
                 fontsize=7.5, color='black', ha='right', va='top')
    axt.set_xlim(0, 180)
    axt.set_xticks([0, 45, 90, 135, 180])
    axt.set_xlabel(r'$\theta$ (deg E of N)', color=INK)
    axt.set_title('(c) orientation', fontsize=10, color=INK)
    axt.grid(alpha=0.25, lw=0.6)

    # (d) the two contrast estimators against the quad-pol LS reference
    axd.plot(dlam[m], z[m], '-', color='#2a78d6', lw=1.6,
             label='coherence phase gradient')
    axd.plot(dlam_node[m], z[m], '-', color='#1baf7a', lw=1.2, alpha=0.8,
             label='cross-pol fringe rate')
    if SITE in REF_DLAM:
        axd.axvline(REF_DLAM[SITE], color='black', lw=1.4, ls='--',
                    label='quad-pol LS, 200-1200 m')
    # Scaled to the estimator that works. The cross-pol fringe rate runs to
    # ~2 here and letting it set the axis would compress the phase-gradient
    # curve - the one that agrees with the LS reference - into the axis
    # line.
    axd.set_xlim(0, 0.35)
    axd.set_xlabel(r'$\Delta\lambda$', color=INK)
    axd.set_title('(d) contrast, two ways', fontsize=10, color=INK)
    axd.grid(alpha=0.25, lw=0.6)
    axd.legend(loc='lower right', fontsize=7.5, frameon=False)

    for ax in (axs, axx, axt, axd):
        ax.set_ylim(Z_SHOW[1], Z_SHOW[0])
    for ax in (axx, axt, axd):
        plt.setp(ax.get_yticklabels(), visible=False)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)

    # Site follows SITE rather than being fixed to Ridge A, so a frame from
    # another survey is not labelled with this one's name.
    fig.suptitle('Quad-pol scattering-matrix inversion, %s %s  '
                 '(track %.0f$^\\circ$)'
                 % (SITE.replace('_', ' ').title(), TAG, track_az),
                 fontsize=12, color=INK)
    out = os.path.join(OUT, 'scar_quadpol_%s.png' % TAG)
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', out)


if __name__ == '__main__':
    main()

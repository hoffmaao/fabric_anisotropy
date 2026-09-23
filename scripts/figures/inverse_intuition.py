"""Why the fabric inverse solver needs a prior, shown on the real estimator.

Slide 13 of the Tubingen talk gives the forward physics: the medium
INTEGRATES the horizontal contrast into a delay,

    Dt(z) = (DS^2 / S_iso) INT_H^z [lam_x(z') - lam_y(z')] dz'

The solver runs that backwards, and this figure is about what the reversal
costs. Integration is a smoother, so inversion is a differentiator, and
differentiating a noisy delay is catastrophic even when the delay itself
is measured well: MEASURED here, a delay with a signal-to-noise of 71
gives a naive-derivative contrast whose error (0.104) is LARGER than the
contrast being measured (0.02-0.09).

Panels are computed with ptt.traveltimeFabricML itself, not sketched, so
the error bars and the resolution are the ones the estimator reports.

Usage: python inverse_intuition.py [out_dir]
"""
import os
import sys

import matplotlib
import numpy as np
from scipy.io import loadmat

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style                        # noqa: E402
from scar_style import INK, MUTED        # noqa: E402

OUT = scar_style.out_dir(sys.argv)
SRC = os.environ.get('INTUITION_DIR',
                     os.path.dirname(os.path.abspath(__file__)))
a = loadmat(os.path.join(SRC, 'intuition.mat'))
b = loadmat(os.path.join(SRC, 'intuition2.mat'))

z = a['z'].ravel()
dl_true = a['dl_true'].ravel()
dt_true = a['dt_true'].ravel()
dt_obs = a['dt_obs'].ravel()
dl_naive = a['dl_naive'].ravel()
S_NS = float(np.ravel(a['S_NS'])[0])
zz = b['zz'].ravel()
dl = b['dl'].ravel()
sg = b['sg'].ravel()
rs = b['rs'].ravel()
dl2 = b['dl2'].ravel()
rs2 = b['rs2'].ravel()

TRUE, MEAS, SOLVE, BAD = '#1baf7a', '#2a78d6', '#4a3aa7', '#eb6834'

fig, ax = plt.subplots(1, 5, figsize=(17.2, 6.4), sharey=True,
                       layout='constrained',
                       gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 0.62]))

# ---- 1. what we want
ax[0].plot(dl_true, z, color=TRUE, lw=2.2)
ax[0].set_title('what we want\n$\\Delta\\lambda(z) = \\lambda_x - \\lambda_y$',
                fontsize=11.5, color=INK)
ax[0].set_xlabel(r'$\Delta\lambda$')
ax[0].set_ylabel('Depth (m)')
ax[0].set_xlim(0, 0.10)

# ---- 2. what the medium gives us: the running integral
ax[1].plot(dt_obs, z, color=MEAS, lw=0.9, alpha=0.85)
ax[1].plot(dt_true, z, color=INK, lw=1.4, ls='--')
ax[1].set_title('what the radar measures\n'
                r'$\Delta t(z) \propto \int_H^z \Delta\lambda\,dz\prime$',
                fontsize=11.5, color=INK)
ax[1].set_xlabel(r'$\Delta t$ (ns)')
ax[1].text(0.04, 0.03, 'the integral SMOOTHS:\n'
           'noise %.2f ns on a %.2f ns span,\nsignal-to-noise %.0f'
           % (S_NS, dt_true[-1], dt_true[-1]/S_NS),
           transform=ax[1].transAxes, fontsize=9.2, color=MUTED,
           va='bottom', linespacing=1.6)

# ---- 3. inverting it the obvious way
g = np.isfinite(dl_naive)
ax[2].plot(dl_naive[g], z[g], color=BAD, lw=0.7, alpha=0.9)
ax[2].plot(dl_true, z, color=TRUE, lw=1.8)
ax[2].set_title('inverted by differencing\n(the obvious way)',
                fontsize=11.5, color=INK)
ax[2].set_xlabel(r'$\Delta\lambda$')
ax[2].set_xlim(-0.35, 0.45)
err_n = np.nanstd(dl_naive[1:] - dl_true[1:])
ax[2].text(0.04, 0.03, 'error %.3f,\nLARGER than the\ncontrast itself' % err_n,
           transform=ax[2].transAxes, fontsize=9.4, color=BAD,
           va='bottom', linespacing=1.6)

# ---- 4. the solver: prior + exact posterior + resolution
ok = np.isfinite(dl)
ax[3].fill_betweenx(zz[ok], (dl - sg)[ok], (dl + sg)[ok],
                    color=SOLVE, alpha=0.22, lw=0)
ax[3].plot(dl[ok], zz[ok], color=SOLVE, lw=2.0)
ax[3].plot(dl_true, z, color=TRUE, lw=1.6, ls='--')
ax[3].set_title('inverted with a prior\n(ptt.traveltimeFabricML)',
                fontsize=11.5, color=INK)
ax[3].set_xlabel(r'$\Delta\lambda$')
ax[3].set_xlim(0, 0.10)
err_s = np.nanstd(dl[ok] - dl_true[ok])
ax[3].text(0.97, 0.985,
           'error %.4f, %.0fx better,\n'
           'and it reports its own $\\pm1\\sigma$'
           % (err_s, err_n / max(err_s, 1e-12)),
           transform=ax[3].transAxes, fontsize=9.0, color=SOLVE,
           ha='right', va='top', linespacing=1.6)
ax[3].text(0.97, 0.845, 'dotted: a stiffer prior,\n600 m instead of 200 m',
           transform=ax[3].transAxes, fontsize=8.8, color='#c0392b',
           ha='right', va='top', linespacing=1.6)

# The stiffer prior, to show what the prior actually controls. Same data,
# same solver, correlation length 600 m instead of 200 m: the wiggles the
# short prior permits are not information, they are how much structure the
# prior was willing to believe.
ok2 = np.isfinite(dl2)
ax[3].plot(dl2[ok2], zz[ok2], color='#c0392b', lw=1.6, ls=':')

# ---- 5. how much of the answer is data
# Shading panel 4 by resolution was the first attempt and it greyed out the
# whole panel, because resolution really IS below 0.5 nearly everywhere.
# That is the lesson rather than a plotting fault, so it gets its own axis:
# at 10 m sampling under a 200 m prior no single depth is data-dominated,
# and what the delay actually pins down is the trend, not the samples.
ax[4].plot(rs, zz, color=SOLVE, lw=1.8)
ax[4].plot(rs2, zz, color='#c0392b', lw=1.6, ls=':')
ax[4].axvline(0.5, color=MUTED, lw=1.0, ls='--')
ax[4].set_title('how much is data?\nresolution diagonal', fontsize=11.5,
                color=INK)
ax[4].set_xlabel('data $\\leftrightarrow$ prior')
ax[4].set_xlim(0, 1)
ax[4].text(0.5, 0.02, '1 = the datum decided it\n0 = the prior did',
           transform=ax[4].transAxes, fontsize=8.6, color=MUTED,
           ha='center', va='bottom', linespacing=1.6)

for a_ in ax:
    a_.grid(alpha=0.25, lw=0.6)
    a_.tick_params(labelsize=9.5)
ax[4].tick_params(labelsize=8.8)
ax[0].set_ylim(1400, 0)

fig.suptitle('The medium integrates, so the solver must differentiate '
             '- which is why it needs a prior', fontsize=13.2, color=INK)
fn = os.path.join(OUT, 'inverse_intuition.png')
fig.savefig(fn, dpi=200, facecolor='white')
print('wrote %s' % fn)

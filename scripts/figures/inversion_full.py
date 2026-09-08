"""The inversion written out in full, continuing the talk's own notation.

Bridges the gap in the Tubingen talk between the processing slides and the
case studies. Those slides end at VV . HH*, the wrapped phase difference,
and the statement

    d(Dphi)/dz  proportional to  lam_x(z) - lam_y(z)

with the caption "more fringes means stronger horizontal fabric". This
slide says how that becomes a number, and how the same fit also returns
the DIRECTION, which counting fringes on one channel pair cannot give.

The through-line is that the wrapped phase on the previous slide is never
unwrapped. Everything here works on the fringe RATE and on the azimuthal
SHAPE of the coherence, both of which survive wrapping.

Every equation is transcribed from ptt.quadpolFabricLS.

Usage: python inversion_full.py [out_dir]
"""
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import FancyBboxPatch         # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style                                     # noqa: E402
from scar_style import INK, MUTED                     # noqa: E402

OUT = scar_style.out_dir(sys.argv)
PREV, NEW, OUTB = '#eef2f7', '#f6f9fd', '#eef7f1'
EDGE, ACC, RED = '#b9c3cf', '#2a78d6', '#c0392b'
MF = 'dejavuserif'

fig = plt.figure(figsize=(14.2, 9.4))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')


def t(x, y, s, size=12, color=INK, ha='center', weight=None):
    ax.text(x, y, s, ha=ha, va='center', fontsize=size, color=color,
            zorder=3, math_fontfamily=MF, fontweight=weight, linespacing=1.55)


def box(cx, cy, w, h, face, ec=EDGE):
    ax.add_patch(FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                                boxstyle='round,pad=0.45,rounding_size=0.9',
                                linewidth=1.1, edgecolor=ec, facecolor=face,
                                zorder=1))


t(50, 96.8, 'Inverting the wrapped phase, without unwrapping it', 16.5, INK,
  weight='bold')

# ------------------------------------------------- what the talk established
box(50, 88.0, 92, 11.0, PREV)
t(50, 92.0, 'ESTABLISHED ON THE PREVIOUS SLIDES', 10, MUTED, weight='bold')
t(26, 87.4, r'$VV\cdot HH^{*}$', 14)
t(51, 87.4, r'$\dfrac{d\,\Delta\phi}{dz} \;\propto\; \lambda_x(z)-\lambda_y(z)$', 14)
t(78, 87.4, r'$\Delta t = \dfrac{\Delta\phi}{2\pi f_c},\quad '
  r'g \equiv \dfrac{\Delta S^{2}}{S_{iso}}$', 13)
t(50, 83.4, 'one channel pair, wrapped, and a fringe rate that gives strength but not direction',
  9.4, MUTED)

# ------------------------------------------------------- 1. every azimuth
t(6, 78.0, 'ONE', 10, ACC, ha='left', weight='bold')
t(6, 75.2, 'every azimuth\nfrom four channels', 9.4, MUTED, ha='left')
t(56, 77.6, r'$T(\psi) = R^{\prime}(\psi)\,S\,R(\psi),\qquad '
  r'S = \left[\ S_{hh}\ \ S_{hv}\ ;\ \ S_{vh}\ \ S_{vv}\ \right]$', 13.5)
t(56, 73.8, r'$C(\psi,z) \;=\; \dfrac{\langle T_{hh}\,T_{vv}^{*}\rangle}'
  r'{\sqrt{\langle |T_{hh}|^{2}\rangle\,\langle |T_{vv}|^{2}\rangle}}$', 14)
t(56, 68.2, 'every observable is QUADRATIC in $S$, so the sweep is a fixed trig '
  'combination\nof 16 moments per range bin - the images are never rotated', 9.4, MUTED)

# ------------------------------------------------------------- 2. the model
box(50, 57.2, 92, 19.0, NEW)
t(6, 64.6, 'TWO', 10, ACC, ha='left', weight='bold')
t(6, 61.8, 'the slide-13\nphysics, per window', 9.4, MUTED, ha='left')
t(56, 64.4, r'$\mu = \cos 2(\psi-\theta_0)$'
  r'$\qquad$'
  r'$\delta(z) = \delta_0 + \delta^{\prime}(z-z_c)$'
  r'$\qquad$'
  r'$\delta^{\prime} = g\,(\lambda_x-\lambda_y)$', 12.6)
t(56, 59.8, r'$\hat{C}(\psi,z) \;=\; \gamma\,'
  r'\dfrac{\frac{1-\mu^{2}}{2} + \frac{1+\mu^{2}}{2}\cos\delta + i\,\mu\sin\delta}'
  r'{1-\frac{1-\mu^{2}}{2}\left(1-\cos\delta\right)} \;+\; P(\psi)$', 14)
t(56, 53.4, r'$P(\psi) = (a_1+i\,a_2)\sin 2\psi + a_3\sin^{2}2\psi$', 12.4)
t(56, 50.2, 'the antenna pedestal: instrument leakage carried as a term to be SOLVED,\n'
  'not left to contaminate the axis', 9.4, ACC)

# --------------------------------------------------------- 3. the objective
t(6, 43.4, 'THREE', 10, ACC, ha='left', weight='bold')
t(6, 40.6, 'fit all azimuths\nand depths at once', 9.4, MUTED, ha='left')
t(56, 43.0, r'$J \;=\; \sum_{\psi}\sum_{z\in W} w(\psi,z)\,'
  r'\left|\,C(\psi,z)-\hat{C}(\psi,z)\,\right|^{2},\qquad '
  r'w = \dfrac{|C|^{2}}{\max(1-|C|^{2},\,0.02)}$', 13)
t(56, 38.6, 'no unwrapping anywhere: the fit sees the fringe RATE and the azimuthal SHAPE',
  9.4, RED)

# ------------------------------------------------------------- 4. the solve
box(50, 27.0, 92, 15.0, NEW)
t(6, 33.6, 'FOUR', 10, ACC, ha='left', weight='bold')
t(6, 29.0, 'separable:\nfour of seven\nsolved exactly', 9.4, MUTED, ha='left')
t(56, 32.0, r'$\left(\Phi^{T}\tilde{W}\Phi+\Lambda\right)\hat{x} = \Phi^{T}\tilde{W}\tilde{d},'
  r'\qquad x=(\gamma,a_1,a_2,a_3)^{T},\qquad \Lambda=\mathrm{diag}(0,\tau,\tau,\tau)$', 12.4)
t(56, 27.4, r'$J(\theta_0,\delta_0,\delta^{\prime}) = \tilde{d}^{T}\tilde{W}\tilde{d}'
  r' - \left(\Phi^{T}\tilde{W}\tilde{d}\right)^{T}'
  r'\left(\Phi^{T}\tilde{W}\Phi+\Lambda\right)^{-1}'
  r'\left(\Phi^{T}\tilde{W}\tilde{d}\right)$', 13)
t(56, 22.6, r'$(\hat{\theta}_0,\hat{\delta}_0,\hat{\delta}^{\prime}) '
  r'= \arg\min_{\ \Theta\times D_0\times D^{\prime}} J,'
  r'\qquad \delta^{\prime}\geq 0,\qquad \theta_0 \ \mathrm{mod}\ \pi$', 12.4)

# ------------------------------------------------------------- 5. the result
box(27, 10.6, 44, 11.0, OUTB)
t(27, 14.6, 'STRENGTH  (what the fringes encode)', 9.6, MUTED, weight='bold')
t(27, 10.4, r'$\lambda_x-\lambda_y \;=\; \hat{\delta}^{\prime}/g$', 14.5)
t(27, 6.8, 'the fringe rate of the previous slide,\nnow a number', 9.0, MUTED)

box(73, 10.6, 44, 11.0, OUTB)
t(73, 14.6, 'DIRECTION  (what one pair cannot give)', 9.6, MUTED, weight='bold')
t(73, 10.4, r'$\hat{\theta}_0$', 15.5)
t(73, 6.8, 'recovered from the azimuthal shape,\nnot from a cross-pol minimum', 9.0, MUTED)

fn = os.path.join(OUT, 'inversion_full.png')
fig.savefig(fn, dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.28)
print('wrote %s' % fn)

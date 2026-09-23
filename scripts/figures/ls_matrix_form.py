"""The orientation-and-strength inversion, written in matrix form.

This is ptt.quadpolFabricLS - the solve that returns the fabric DIRECTION
as well as the contrast - stated as what it is: a SEPARABLE NONLINEAR
WEIGHTED LEAST SQUARES problem, solved by VARIABLE PROJECTION
(Golub & Pereyra 1973).

Separable because the seven parameters split cleanly. Three are non-linear
- the axis theta_0, the accumulated phase delta_0 and its rate delta' -
and the model depends on them through cos and sin. Four are LINEAR - the
coherence scale gamma and the three antenna-pedestal amplitudes - and the
model is a plain weighted sum of them. Variable projection exploits that:
at each trial value of the non-linear three, the linear four are solved
exactly by a 4x4 normal system and substituted back, leaving a cost that
depends on the non-linear parameters alone.

Transcribed from the code: the grid triple is the nested loop over ths,
dd_grid and d0_grid; the 4x4 system with its ridge on the nuisance
diagonal is H_solve.

Usage: python ls_matrix_form.py [out_dir]
"""
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style                                      # noqa: E402
from scar_style import INK, MUTED                      # noqa: E402

OUT = scar_style.out_dir(sys.argv)
NL, LIN, ACC = '#eef7f1', '#fdf6ec', '#2a78d6'
EDGE, RED = '#b9c3cf', '#c0392b'
MF = 'dejavuserif'

fig = plt.figure(figsize=(13.4, 10.4))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')


def t(x, y, s, size=12, color=INK, ha='center', weight=None):
    ax.text(x, y, s, ha=ha, va='center', fontsize=size, color=color,
            zorder=3, math_fontfamily=MF, fontweight=weight, linespacing=1.6)


def box(cx, cy, w, h, face, ec=EDGE):
    ax.add_patch(FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                                boxstyle='round,pad=0.5,rounding_size=1.0',
                                linewidth=1.1, edgecolor=ec, facecolor=face,
                                zorder=1))


t(50, 96.6,
  'Separable nonlinear weighted least squares, by variable projection',
  16, INK, weight='bold')
t(50, 92.8, 'the solve that returns fabric DIRECTION as well as strength '
  r'($\mathtt{ptt.quadpolFabricLS}$)', 11, MUTED)

# ---------------------------------------------------------------- the split
box(25, 84.0, 44, 9.5, NL)
t(25, 87.0, 'NON-LINEAR   (searched on a grid)', 10, MUTED, weight='bold')
t(25, 83.0, r'$m_n = (\theta_0,\ \delta_0,\ \delta\prime)$', 14)
t(25, 79.8, 'axis, accumulated phase, phase rate', 9, MUTED)

box(74, 84.0, 44, 9.5, LIN)
t(74, 87.0, 'LINEAR   (eliminated exactly)', 10, MUTED, weight='bold')
t(74, 83.0, r'$x = (\gamma,\ a_1,\ a_2,\ a_3)^{T}$', 14)
t(74, 79.8, 'coherence scale, three pedestal amplitudes', 9, MUTED)

# ------------------------------------------------------------- the stacking
t(50, 74.6, 'Stack the window: azimuths $\\psi_k$ by depths $z_j$, real part '
  'over imaginary part', 10.5, MUTED)
t(50, 70.8, r'$\tilde{d} = \left[\,\mathrm{Re}\,C\ ;\ \mathrm{Im}\,C\,\right]'
  r'\in \mathbb{R}^{2KJ}$'
  r'$\qquad\qquad$'
  r'$\tilde{W} = \mathrm{diag}\left[\,w\ ;\ w\,\right],\quad '
  r'w = \dfrac{|C|^{2}}{\max(1-|C|^{2},\,0.02)}$', 13)

# --------------------------------------------------------- the design matrix
t(15.5, 61.6, r'$\Phi(m_n) \;=$', 14.5, ha='center')

X0, Y0, CW, RH = 29.0, 57.0, 8.4, 4.6
cols = [r'$\mathrm{Re}\,H$', r'$s$', r'$0$', r'$q$',
        r'$\mathrm{Im}\,H$', r'$0$', r'$s$', r'$0$']
for r_ in range(2):
    for c_ in range(4):
        x = X0 + c_*CW
        y = Y0 + (1-r_)*RH
        ax.add_patch(Rectangle((x, y), CW, RH, facecolor='white',
                               edgecolor=EDGE, lw=0.9, zorder=2))
        t(x + CW/2, y + RH/2, cols[r_*4 + c_], 12.5)
# brackets
for xb in (X0 - 0.9, X0 + 4*CW + 0.9):
    ax.plot([xb, xb], [Y0, Y0 + 2*RH], color=INK, lw=1.6, zorder=3)
    for yy in (Y0, Y0 + 2*RH):
        ax.plot([xb, xb + (1.1 if xb < X0 else -1.1)], [yy, yy],
                color=INK, lw=1.6, zorder=3)
for c_, lab in enumerate([r'$\gamma$', r'$a_1$', r'$a_2$', r'$a_3$']):
    t(X0 + c_*CW + CW/2, Y0 + 2*RH + 2.2, lab, 12, ACC)
t(X0 - 2.1, Y0 + 1.5*RH, 'real', 8.8, MUTED, ha='right')
t(X0 - 2.1, Y0 + 0.5*RH, 'imag', 8.8, MUTED, ha='right')
t(X0 + 2*CW, Y0 - 3.0,
  r'$H_{kj} = \dfrac{A_k + B_k\cos\delta_j + i\,\mu_k\sin\delta_j}'
  r'{\max\left[1 - A_k(1-\cos\delta_j),\ 0.05\right]}$', 12)
t(X0 + 2*CW, Y0 - 7.4,
  r'$\mu_k = \cos 2(\psi_k-\theta_0),\quad '
  r'A_k = \frac{1-\mu_k^{2}}{2},\quad B_k = \frac{1+\mu_k^{2}}{2}$', 10.5)
t(X0 + 2*CW, Y0 - 10.8,
  r'$\delta_j = \delta_0 + \delta\prime(z_j - z_c),\quad '
  r's_k = \sin 2\psi_k,\quad q_k = s_k^{2}$', 10.5)

# --------------------------------------------------------- normal equations
box(50, 34.0, 84, 12.0, LIN)
t(50, 38.4, 'At every grid node, solve the linear block exactly', 10.5, MUTED,
  weight='bold')
t(50, 34.2, r'$\left(\Phi^{T}\tilde{W}\Phi + \Lambda\right)\hat{x} '
  r'\;=\; \Phi^{T}\tilde{W}\tilde{d}$'
  r'$\qquad$'
  r'$\Lambda = \mathrm{diag}(0,\ \tau,\ \tau,\ \tau)$', 14)
t(50, 30.2, r'a 4x4 system; the ridge $\tau$ steadies the pedestal where one '
  r'window cannot separate it from the fabric, and $\gamma$ carries no ridge',
  9.2, MUTED)

# ------------------------------------------------------------ the projection
box(50, 18.6, 84, 12.6, NL)
t(50, 23.4, 'Substitute back: a cost in the non-linear three alone', 10.5,
  MUTED, weight='bold')
t(50, 19.4, r'$J(m_n) \;=\; \tilde{d}^{T}\tilde{W}\tilde{d} \;-\; '
  r'\left(\Phi^{T}\tilde{W}\tilde{d}\right)^{T}'
  r'\left(\Phi^{T}\tilde{W}\Phi + \Lambda\right)^{-1}'
  r'\left(\Phi^{T}\tilde{W}\tilde{d}\right)$', 13.5)
t(50, 15.0, r'$\hat{m}_n = '
  r'\arg\min_{\,\Theta\,\times\,D_0\,\times\,D\prime} J,'
  r'\qquad \delta\prime \geq 0,\qquad \gamma \in [0,\ 1.05]$', 12.5)

# ------------------------------------------------------------------- output
box(24, 5.4, 40, 7.6, '#f6f9fd')
t(24, 7.6, 'STRENGTH', 9.5, MUTED, weight='bold')
t(24, 4.2, r'$\Delta\lambda = \hat{\delta}\prime / g$', 13)
box(70, 5.4, 46, 7.6, '#f6f9fd')
t(70, 7.6, 'DIRECTION', 9.5, MUTED, weight='bold')
t(70, 4.2, r'$\hat{\theta}_0$ , unique modulo $\pi$', 13)

fn = os.path.join(OUT, 'ls_matrix_form.png')
fig.savefig(fn, dpi=200, facecolor='white', bbox_inches='tight',
            pad_inches=0.3)
print('wrote %s' % fn)

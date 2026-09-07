"""How the least-squares fabric fit is put together, as a flow chart.

The same equations as ls_equations.py, but with every term named and
arranged in the order the estimator uses them: the measured field in at
the top, the weight and the linear nuisances entering from the side, the
fitted strength and direction coming out at the bottom.

Transcribed from +ptt/quadpolFabricLS.m - the coherence model and the
pedestal from its header, the CRB weight from line 213, the pooled
constant-direction vote from the cost that mode minimises. If the
estimator changes this figure is stale until it is regenerated.

Usage: python ls_flowchart.py [out_dir]
"""
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style                                     # noqa: E402
from scar_style import INK, MUTED                     # noqa: E402

OUT = scar_style.out_dir(sys.argv)

MAIN, SIDE, OUTB = '#f6f9fd', '#fdf6ec', '#eef7f1'
EDGE, ACC = '#b9c3cf', '#2a78d6'
MF = 'dejavuserif'

# LAYOUT. The main path is a single column; anything entering from the
# side gets its own column that does not overlap it. Coordinates are
# explicit rather than derived, because the first version let the side
# column run under the objective box and the box titles collide with
# their own first equation row.
CX, CW = 38.0, 60.0            # main column centre and width  (8 .. 68)
SX, SW = 84.0, 27.0            # side column centre and width  (70.5 .. 97.5)

fig = plt.figure(figsize=(13.0, 15.6))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')


def panel(cx, cy, w, h, face):
    ax.add_patch(FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                                boxstyle='round,pad=0.55,rounding_size=1.1',
                                linewidth=1.1, edgecolor=EDGE,
                                facecolor=face, zorder=2))


def title(cx, y, s):
    ax.text(cx, y, s, ha='center', va='center', fontsize=10.5,
            color=MUTED, fontweight='bold', zorder=3)


def eq(cx, y, s, size=15.5, color=INK):
    ax.text(cx, y, s, ha='center', va='center', fontsize=size, color=color,
            zorder=3, math_fontfamily=MF)


def note(cx, y, s, size=9.3, color=MUTED):
    ax.text(cx, y, s, ha='center', va='center', fontsize=size, color=color,
            zorder=3, linespacing=1.55, math_fontfamily=MF)


def arrow(x0, y0, x1, y1, color=INK, lw=1.5):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>',
                                 mutation_scale=15, linewidth=lw, color=color,
                                 zorder=1, shrinkA=0, shrinkB=0))


# --------------------------------------------------------- measured field
panel(CX, 95.2, CW, 7.4, MAIN)
title(CX, 97.6, 'MEASURED FIELD')
eq(CX, 95.0, r'$C(\psi, z)$', 16)
note(CX, 92.6, r'HH$-$VV coherence at every synthetic azimuth $\psi$, every depth $z$')
arrow(CX, 91.2, CX, 88.6)

# ---------------------------------------------------------------- weight
panel(CX, 84.4, CW, 8.4, MAIN)
title(CX, 87.4, 'WEIGHT  (Cramer-Rao)')
eq(CX, 84.2, r'$w(\psi,z) \;=\; \frac{|C|^{2}}{\max(1-|C|^{2},\ 0.02)}$', 16)
note(CX, 81.2, 'down-weights decorrelated cells; the 0.02 floor caps the gain')
arrow(CX, 79.9, CX, 77.4)

# ----------------------------------------------------------------- model
panel(CX, 65.8, CW, 22.4, MAIN)
title(CX, 75.6, r'FORWARD MODEL   (one depth window $W$)')
eq(CX, 72.4, r'$\mu \;=\; \cos\left[\,2(\psi-\theta_0)\,\right]$', 16)
note(CX, 70.0, r'the fabric direction $\theta_0$ enters the fit ONLY here',
     9.6, ACC)
eq(CX, 67.2, r"$\delta(z) \;=\; \delta_0 \;+\; \delta'\,(z-z_c)$", 15)
eq(CX - 14.5, 63.4,
   r'$N = \frac{1-\mu^{2}}{2} + \frac{1+\mu^{2}}{2}\cos\delta + i\mu\sin\delta$', 13)
eq(CX + 17.0, 63.4,
   r'$D = 1 - \frac{1-\mu^{2}}{2}(1-\cos\delta)$', 13)
eq(CX, 59.4, r'$\hat{C}(\psi,z) \;=\; \gamma\,\dfrac{N}{D} \;+\; P(\psi)$', 16.5)
note(CX, 56.6,
     r'$\delta_0$  accumulated two-way phase        '
     r"$\delta'$  its rate with depth", 9.6)

# ------------------------------------------------- linear nuisances (side)
panel(SX, 65.8, SW, 15.0, SIDE)
title(SX, 71.6, 'LINEAR NUISANCES')
eq(SX, 68.6, r'$\gamma$', 15)
eq(SX, 65.6, r'$P(\psi) = (a_1 + i a_2)\sin 2\psi + a_3 \sin^{2} 2\psi$', 10.4)
note(SX, 61.6, 'coherence scale, and the antenna pedestal\n'
     'that antenna-locks power-only methods;\nboth solved in closed form', 8.8)
arrow(SX - SW/2 - 0.4, 65.8, CX + CW/2 + 0.4, 65.8, MUTED, 1.3)

arrow(CX, 54.6, CX, 52.0)

# ------------------------------------------------------------- objective
panel(CX, 47.6, CW, 8.6, MAIN)
title(CX, 50.7, 'OBJECTIVE')
eq(CX, 47.4, r'$J \;=\; \sum_{\psi}\ \sum_{z\,\in\,W}\ w(\psi,z)\ '
   r'\left|\ C(\psi,z) - \hat{C}(\psi,z)\ \right|^{2}$', 15)
note(CX, 44.3, 'all azimuths and all depths in the window, fitted at once')
arrow(CX, 43.3, CX, 40.7)

# ---------------------------------------------------------------- solve
panel(CX, 35.0, CW, 10.8, MAIN)
title(CX, 39.2, 'MINIMISE')
eq(CX, 36.2, r'$\theta_0$  on a grid$\quad$'
   r"$\delta_0,\ \delta'$  non-linear$\quad$"
   r'$\gamma,\ a$  closed form at each node', 11.4)
eq(CX, 33.2, r"$\delta' \;\geq\; 0$$\qquad$"
   r'$\theta_0$  unique modulo $\pi$', 12.6)
note(CX, 30.8, r'the constraint removes the $(\theta_0+90^\circ,\ -\delta)$ '
     r'twin, which fits identically')
arrow(CX, 29.6, CX, 27.6)

# --------------------------------------------------------------- outputs
arrow(CX, 27.6, 22.0, 25.4)
arrow(CX, 27.6, 54.0, 25.4)
panel(22.0, 20.6, 26.0, 8.6, OUTB)
title(22.0, 23.6, 'STRENGTH')
eq(22.0, 20.6, r"$\Delta\lambda \;=\; \delta' / g$", 16)
note(22.0, 17.9, r'$g = 2\pi f_c \Delta\varepsilon / (n c)$', 9.6)
panel(54.0, 20.6, 26.0, 8.6, OUTB)
title(54.0, 23.6, 'DIRECTION')
eq(54.0, 20.8, r'$\theta_0$', 17)
note(54.0, 17.9, 'one axis per window,\nmodulo $180^\\circ$', 9.3)
arrow(54.0, 16.3, 54.0, 13.6)

# ------------------------------------------------ constant-direction mode
panel(CX, 8.4, CW, 9.4, MAIN)
title(CX, 11.8, 'CONSTANT-DIRECTION MODE   (optional)')
eq(CX, 8.6, r'$\hat{\theta}_0 \;=\; \mathrm{arg\,min}_{\ \theta}\ '
   r'\sum_{W}\ \frac{J_W(\theta)}{C^{2}_{W}}$', 16)
note(CX, 5.2, 'one axis for the whole segment; each window normalised by its own\n'
     'weighted coherence power, so the surface cannot outvote the column')

fn = os.path.join(OUT, 'ls_flowchart.png')
fig.savefig(fn, dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.3)
print('wrote %s' % fn)

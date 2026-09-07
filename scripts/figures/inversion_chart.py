"""The inversion procedure, in the notation of slide 13 of the Tubingen talk.

Slide 13 gives the forward physics - the medium accumulates phase in
proportion to the integral of the horizontal contrast. This chart runs the
same symbols the other way, through the steps actually implemented here, so
it can sit beside that slide without a change of notation.

Kept in the talk's symbols throughout: lam_x - lam_y for the contrast,
Dphi and Dt for the accumulated phase and delay, DS^2 / S_iso for the
constant that links them. Internally the code calls the contrast dlam and
the phase delta; those are the same quantities.

Sized 16:9 for a slide.

Usage: python inversion_chart.py [out_dir]
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

ICE, MEAS, FIT, OPT = '#eef7f1', '#f6f9fd', '#f6f9fd', '#fdf6ec'
EDGE, ACC, RED = '#b9c3cf', '#2a78d6', '#c0392b'
MF = 'dejavuserif'

fig = plt.figure(figsize=(13.33, 7.5))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')


def panel(cx, cy, w, h, face, ec=EDGE):
    ax.add_patch(FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                                boxstyle='round,pad=0.5,rounding_size=1.0',
                                linewidth=1.2, edgecolor=ec, facecolor=face,
                                zorder=2))


def txt(cx, y, s, size, color=INK, weight=None, ha='center'):
    ax.text(cx, y, s, ha=ha, va='center', fontsize=size, color=color,
            zorder=3, math_fontfamily=MF, fontweight=weight, linespacing=1.6)


def arrow(x0, y0, x1, y1, color=INK, lw=1.8):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>',
                                 mutation_scale=17, linewidth=lw, color=color,
                                 zorder=1, shrinkA=0, shrinkB=0))


txt(50, 96.0, 'How the inversion runs, in the notation of the forward slide',
    15.5, INK, 'bold')

# ---------------------------------------------------------------- the ice
CY, H, W = 74.0, 20.0, 21.0
xs = [13.5, 37.5, 61.5, 85.5]

panel(xs[0], CY, W, H, ICE)
txt(xs[0], 82.0, 'THE ICE', 11, MUTED, 'bold')
txt(xs[0], 77.4, r'$\Delta\phi(z) \;\propto\; \int_H^z \!(\lambda_x-\lambda_y)\,dz\prime$', 14)
txt(xs[0], 72.8, r'$\Delta t = \dfrac{\Delta\phi}{2\pi f_c}$', 14)
txt(xs[0], 66.8, 'birefringence accumulates\nphase with depth', 9.0, MUTED)

# ------------------------------------------------------------- what we get
panel(xs[1], CY, W, H, MEAS)
txt(xs[1], 82.0, 'WHAT WE RECORD', 11, MUTED, 'bold')
txt(xs[1], 76.6, r'$C(\psi, z)$', 16)
txt(xs[1], 71.0, 'HH$-$VV coherence at every\nsynthesized azimuth $\\psi$,\n'
    'from the four channels', 9.4, MUTED)

# --------------------------------------------------------------- the fit
panel(xs[2], CY, W, H, FIT)
txt(xs[2], 82.0, 'FIT, PER WINDOW', 11, MUTED, 'bold')
txt(xs[2], 78.0, r'$\theta_0 \qquad \dfrac{d\Delta\phi}{dz}$', 15)
txt(xs[2], 70.4, 'one axis and one phase rate,\nfrom the SHAPE of $C(\\psi,z)$\n'
    'over all azimuths at once', 9.0, MUTED)

# --------------------------------------------------------------- convert
panel(xs[3], CY, W, H, FIT)
txt(xs[3], 82.0, 'CONVERT', 11, MUTED, 'bold')
txt(xs[3], 77.0,
    r'$\lambda_x-\lambda_y \;=\; \dfrac{S_{iso}}{\Delta S^{2}}\,'
    r'\dfrac{d\Delta t}{dz}$', 14)
txt(xs[3], 70.0, r'$\Delta S^{2}/S_{iso} = 0.0637$ ns per m' '\n'
    'per unit contrast, at 750 MHz', 9.0, MUTED)

for i in range(3):
    arrow(xs[i] + W/2 + 0.6, CY, xs[i+1] - W/2 - 0.6, CY)

# the reversal, called out
ax.annotate('', xy=(xs[0], 61.3), xytext=(xs[3], 61.3),
            arrowprops=dict(arrowstyle='-|>', color=ACC, lw=1.6,
                            connectionstyle='arc3,rad=0.07'))
txt(50, 56.6, 'the solver runs the forward physics backwards: '
    'the medium integrates, so we differentiate', 10.6, ACC)

# --------------------------------------------------------------- options
arrow(50, 53.4, 27.0, 47.0, MUTED)
arrow(50, 53.4, 73.0, 47.0, MUTED)

panel(27.0, 33.5, 40.0, 25.0, OPT)
txt(27.0, 43.0, 'IF THE AXIS DOES NOT TURN', 11, MUTED, 'bold')
txt(27.0, 38.0, r'$\theta_0$  held for the whole segment,'
    '\n' r'$\lambda_x-\lambda_y(z)$  still free', 12.5)
txt(27.0, 30.0,
    'one axis solved from every window at once,\n'
    'each normalised by its own coherence power\n'
    'so the surface cannot outvote the column', 9.4, MUTED)
txt(27.0, 24.2, 'a decorrelating column needs this:\n'
    'the deep windows carry one fewer parameter', 9.4, ACC)

panel(73.0, 33.5, 40.0, 25.0, OPT)
txt(73.0, 43.0, 'IF THE DELAY CAN BE MEASURED', 11, MUTED, 'bold')
txt(73.0, 38.6, r'$\Delta t(z)$  measured absolutely,'
    '\n' 'then the integral inverted', 12.5)
txt(73.0, 32.6, r'$\tilde{m} = (G^{T}C_d^{-1}G + C_m^{-1})^{-1}'
    r'(G^{T}C_d^{-1}\Delta t + C_m^{-1}m_{prior})$', 10.8)
txt(73.0, 28.0, r'$\tilde{C}_m = (G^{T}C_d^{-1}G + C_m^{-1})^{-1}$', 10.8)
txt(73.0, 23.6, 'uses the ACCUMULATED delay, not the local rate;\n'
    'linear, so this posterior is exact', 9.4, ACC)

# ------------------------------------------------------------- the caveat
panel(50, 9.0, 82.0, 11.0, '#fff5f5', RED)
txt(50, 12.4, 'WHAT THE PRIOR IS DOING', 10.5, RED, 'bold')
txt(50, 7.6, 'Inverting an integral differentiates the noise. Measured on a delay with '
    'signal-to-noise 71,\ndifferencing gives an error of 0.104 - larger than the '
    '0.02-0.09 contrast being measured. '
    'The prior is what makes the problem solvable, and it sets what counts as structure.',
    9.6, INK)

fn = os.path.join(OUT, 'inversion_chart.png')
fig.savefig(fn, dpi=200, facecolor='white', bbox_inches='tight', pad_inches=0.25)
print('wrote %s' % fn)

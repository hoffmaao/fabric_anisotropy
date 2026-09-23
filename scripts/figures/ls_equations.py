"""The least-squares problem solved by ptt.quadpolFabricLS, as a figure.

Every line here is transcribed from the estimator, not restated from the
literature: the model and the pedestal from the header of
+ptt/quadpolFabricLS.m, the weight from its line 213
(Wc = |C|^2 / max(1 - |C|^2, 0.02)), the pooled constant-direction vote
from the normalised cost it minimises, and the contrast conversion from
the same grad_per_dlam that +ptt/quadpolFabric.m uses. If the estimator
changes, this figure is wrong until it is regenerated.

Deliberately UNLABELLED - no equation numbers, no tags - so nothing here
can be cited by number against a version that has moved on.

Usage: python ls_equations.py [out_dir]
"""
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style                        # noqa: E402
from scar_style import INK               # noqa: E402

OUT = scar_style.out_dir(sys.argv)

# The fitted parameters per window: the fabric DIRECTION theta_0, the
# accumulated two-way eigenmode phase and its rate, the coherence scale,
# and the three antenna-frame pedestal coefficients.
LINES = [
    (r"$\hat{J}(\theta_0,\ \delta_0,\ \delta',\ \gamma,\ a) \;=\; "
     r"\sum_{\psi}\ \sum_{z\,\in\,W}\ "
     r"w(\psi,z)\ \left|\ C(\psi,z)\ -\ \hat{C}(\psi,z)\ \right|^{2}$", 1.18),
    (r"$\hat{C}(\psi,z) \;=\; \gamma\ \frac{N(\mu,\delta)}{D(\mu,\delta)}"
     r"\ +\ P(\psi)$", 1.26),
    (r"$\mu \;=\; \cos\left[\,2(\psi-\theta_0)\,\right]$", 1.06),
    (r"$N \;=\; \frac{1-\mu^{2}}{2}\ +\ \frac{1+\mu^{2}}{2}\,\cos\delta"
     r"\ +\ i\,\mu\,\sin\delta$", 1.24),
    (r"$D \;=\; 1\ -\ \frac{1-\mu^{2}}{2}\,\left(1-\cos\delta\right)$", 1.24),
    (r"$\delta(z) \;=\; \delta_0\ +\ \delta'\,(z-z_c)$", 1.06),
    (r"$P(\psi) \;=\; (a_1+i\,a_2)\,\sin 2\psi\ +\ a_3\,\sin^{2} 2\psi$", 1.06),
    (r"$w(\psi,z) \;=\; \frac{|C|^{2}}{\max\left(1-|C|^{2},\ 0.02\right)}$", 1.24),
    (r"$\Delta\lambda \;=\; \delta'\,/\,g,"
     r"\qquad g \;=\; \frac{2\pi f_c\,\Delta\varepsilon}{n\,c}$", 1.22),
    (r"$\delta' \,\geq\, 0,"
     r"\qquad \theta_0\ \ \mathrm{unique\ modulo}\ \pi$", 1.04),
    (r"$\hat{\theta}_0 \;=\; \mathrm{arg\,min}_{\ \theta}\ \sum_{W}\ "
     r"\frac{J_W(\theta)}{C^{2}_{W}}$", 1.26),
]

fig = plt.figure(figsize=(9.6, 9.4))
fig.patch.set_facecolor('white')
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')

# Vertical layout: proportional gaps so the taller displayed fractions do
# not collide with their neighbours.
gaps = [1.62, 1.52, 1.05, 1.58, 1.42, 1.05, 1.10, 1.62, 1.55, 1.05, 1.60]
total = sum(gaps)
y = 0.965
for (tex, size), gap in zip(LINES, gaps):
    ax.text(0.5, y, tex, transform=ax.transAxes, ha='center', va='top',
            fontsize=15.0 * size, color=INK, math_fontfamily='dejavuserif')
    y -= (gap / total) * 0.93

fn = os.path.join(OUT, 'ls_equations.png')
fig.savefig(fn, dpi=220, facecolor='white', bbox_inches='tight', pad_inches=0.35)
print('wrote %s' % fn)

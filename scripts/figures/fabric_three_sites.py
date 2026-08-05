"""Horizontal fabric contrast at the three SCAR sites, on one axis.

Companion to fabric_sections.py, which draws each site's 2D section; this
is the cross-site comparison. Both read the same
run_sections.m output, so the three profiles are produced at the same
retuned settings (small along-track blocks, per-frame coherence gate) and
differ only in the dtau estimator each site supports:

  Ridge A   SNAPHU + coregistration - SNAPHU there is exactly consistent
            with the wrapped phase (4.43 fringes either way of counting),
            and delta-k is suppressed at depth
  Thwaites  delta-k - the joint chain does not fit its own observations
            there (5.2 ns rms against 0.2 ns)
  NEGIS     delta-k - the only option; no SNAPHU, no coregistration

READ THE AXES BEFORE COMPARING MAGNITUDES. dlam is the difference of the
two horizontal eigenvalues on the profile axis and its perpendicular, and
those axes are different at every site (Ridge A 179 deg, Thwaites 68 deg,
NEGIS 16 deg E of N). A larger dlam therefore does not straightforwardly
mean stronger fabric - it can equally mean the profile happens to lie
closer to a principal axis. Putting all three on one dlam scale is
correct for showing what each survey measured; converting to a
flow-relative or principal-axis frame is a separate step.

Intervals whose node coherence is too low to constrain dlam are dropped
rather than drawn, which is why the profiles stop at different depths and
why none of them extends into the firn.

Usage: python fabric_three_sites.py <out_dir>
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
from scipy.io import loadmat             # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
DATA = os.path.expanduser(os.environ.get('SCAR_DATA', '~/data/opr/scar'))
INK, MUTED = '#0b0b0b', '#52514e'

# Validated categorical slots (dataviz reference palette, light mode):
# adjacent-pair CVD separation dE 24.7 / 18.9, normal-vision 33.6 / 21.1.
SITES = [
    ('ridge_a', 'Ridge A divide', '#2a78d6'),
    ('thwaites', 'Thwaites eastern shear margin', '#eb6834'),
    ('negis', 'NEGIS onset (EGRIP)', '#1baf7a'),
]
# Below this the interval is regularizer output, not a measurement: it
# lies inside the surface reference band where dtau is pinned to zero, and
# the measured fringe rate there is ~0 (see opr_fabric/server/fringe_check.m)
Z_MIN = 200.0
COH_MIN = 0.45      # node coherence an interval must reach to be drawn


def field(rec, name):
    v = rec[name]
    while isinstance(v, np.ndarray) and v.dtype == object and v.size == 1:
        v = v.item()
    return np.atleast_1d(np.asarray(v, float))


def track_azimuth(lat, lon):
    """Principal azimuth of the track, degrees E of N, as an axis (0-180)."""
    ok = np.isfinite(lat) & np.isfinite(lon)
    la, lo = lat[ok], lon[ok]
    pts = np.c_[(lo - lo.mean()) * np.cos(np.radians(la.mean())),
                la - la.mean()]
    pts = pts - pts.mean(0)
    _, _, vt = np.linalg.svd(pts, full_matrices=False)
    return np.degrees(np.arctan2(vt[0][0], vt[0][1])) % 180


def profile(rec):
    """Median dlam and 16-84% block spread per interval, coherence-gated."""
    nint = int(field(rec, 'nint')[0])
    dlam = field(rec, 'dlam').reshape(nint, -1)
    qual = field(rec, 'quality').reshape(nint, -1)
    top = field(rec, 'top').reshape(nint, -1)
    bot = field(rec, 'bot').reshape(nint, -1)
    z = np.nanmedian((top + bot) / 2, axis=1)
    med = np.full(nint, np.nan)
    lo = np.full(nint, np.nan)
    hi = np.full(nint, np.nan)
    for i in range(nint):
        if z[i] < Z_MIN:
            continue
        ok = np.isfinite(dlam[i]) & (qual[i] >= COH_MIN)
        if ok.sum() >= 3:
            med[i] = np.median(dlam[i][ok])
            lo[i], hi[i] = np.percentile(dlam[i][ok], [16, 84])
    return z, med, lo, hi


def main():
    os.makedirs(OUT, exist_ok=True)
    S = loadmat(os.path.join(DATA, 'fabric_sections.mat'),
                squeeze_me=True)['S']

    fig, ax = plt.subplots(figsize=(7.4, 7.6), layout='constrained')
    for key, label, color in SITES:
        if key not in S.dtype.names:
            print('no section for', key)
            continue
        rec = S[key].item() if S[key].dtype == object else S[key]
        z, med, lo, hi = profile(rec)
        az = track_azimuth(field(rec, 'lat'), field(rec, 'lon'))
        ok = np.isfinite(med)
        ax.fill_betweenx(z[ok], lo[ok], hi[ok], color=color, alpha=0.18, lw=0)
        ax.plot(med[ok], z[ok], '-', color=color, lw=2.2,
                label=r'%s   ($\perp$ %.0f$^\circ$, $\parallel$ %.0f$^\circ$)'
                      % (label, (az + 90) % 180, az))
        print('%-9s %2d intervals drawn, %.0f-%.0f m, dlam %.3f..%.3f'
              % (key, ok.sum(), z[ok].min(), z[ok].max(),
                 np.nanmin(med), np.nanmax(med)))

    ax.axvline(0, color=MUTED, lw=0.8)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.invert_yaxis()
    ax.set_xlabel(r'$\Delta\lambda = \lambda_{\perp} - \lambda_{\parallel}$'
                  '   (eigenvalues on each profile\'s own axes)', color=INK)
    ax.set_ylabel('depth (m)', color=INK)
    ax.set_title('Horizontal fabric contrast by layer stripping '
                 '(Rathmann 2026)', fontsize=12, color=INK)
    ax.legend(loc='lower right', fontsize=8.5, frameon=False)
    ax.annotate('bands: block-to-block 16-84%%.  Intervals above %.0f m or '
                'below node coherence %.2f are omitted -\nthey are '
                'regularizer output, not measurements.  Axes differ per '
                'site, so magnitudes are not directly comparable.'
                % (Z_MIN, COH_MIN),
                xy=(0.0, -0.085), xycoords='axes fraction', fontsize=7.5,
                color=MUTED, ha='left', va='top')

    out = os.path.join(OUT, 'fabric_three_sites.png')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', out)


if __name__ == '__main__':
    main()

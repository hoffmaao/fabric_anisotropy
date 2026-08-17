"""Cross-check our LS fabric estimate against PoRaPy's crossover method.

PoRaPy (Lilien, PoRaCo) measures horizontal anisotropy the way the SHARAD
crossover analysis does: correlate two traces sounded at different antenna
orientations in sliding depth bins, take the sub-sample delay that maximises
correlation in each bin, and read the fabric off how that delay GROWS with
depth. That is a completely different observable from our quad-pol LS fit,
which matches the complex HH-VV coherence field against the birefringence
model, so agreement between them is a real check on our implementation.

The two traces here are synthesized from one frame's full scattering matrix
at the fabric axis and across it (scripts/prototypes/extract_polpair.py on
mem1), so no crossing geometry is involved and nothing depends on
coregistering two survey lines. This is the cleanest possible version of
their measurement.

WHAT IT SHOWED on Ridge A frame 009 (11 Aug 2026), kept as the record of
the comparison: the raw differential phase between the two synthesized
polarisations accumulates 2.96 cycles = 3.94 ns over 250-1400 m and gives
dlam 0.035/0.057/0.065 by band against the LS estimator's
0.035/0.056/0.062 - an estimator-independent confirmation of the LS chain.
The sliding-bin envelope correlation recovers only ~a third of that delay
(sub-wavelength shifts against broad envelope features at 750 MHz), which
is NOT a criticism of PoRaPy: their real method correlates on the BED
reflection, which this accum radar never reaches (no bed echo at Ridge A
to 1879 m or Thwaites to 1964 m; the coincident rds/MCoRDS picks put the
Ridge A bed at ~2912 m). A bed-referenced comparison needs the rds data.

CONVERSION. A two-way delay difference dt accumulated over depth relates to
the eigenvalue contrast by the same constant the LS estimator uses:

    d(phi)/dz = grad_per_dlam * dlam   [rad/m],  phi = 2 pi fc dt

    => d(dt)/dz = grad_per_dlam * dlam / (2 pi fc)

PoRaPy returns the delay in SAMPLES, so it is scaled by the record's own
sample interval before that is applied. Both sides therefore go through the
same physical constants and the comparison is of the measurement, not of
two different calibrations.

Usage: python crossover_check.py <polpair.npz> <out_dir>
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402

# PoRaPy is Lilien's PoRaCo package (not vendored here): point PORAPY_SRC at
# its src/ directory. It imports line_profiler purely for the @profile
# decorator, so a no-op module is injected rather than requiring the
# package.
import types                             # noqa: E402
_lp = types.ModuleType('line_profiler')
_lp.profile = lambda f: f
sys.modules['line_profiler'] = _lp
sys.path.insert(0, os.environ.get('PORAPY_SRC', '.'))
from PoRaPy.alignment import binned_best_offsets, binned_corrs  # noqa: E402

NPZ = sys.argv[1] if len(sys.argv) > 1 else 'polpair_009.npz'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'

C_ICE = 1.685e8            # PoRaPy's ice velocity, m/s
FC = 750e6                 # our system centre frequency
EPS_BAR, DEPS = 3.15, 0.034
N_ICE = np.sqrt(EPS_BAR)
C0 = 299792458.0
# rad/m per unit dlam - the same constant ptt.quadpolFabricLS uses
GRAD_PER_DLAM = 2 * np.pi * FC * DEPS / (N_ICE * C0)
BIN_M = 150.0              # depth bin for the correlation
OVERSAMPLE = 3
MAX_SHIFT = 3.0            # samples
INCREMENT = 0.01


def main():
    d = np.load(NPZ)
    z = d['z'].astype(float)
    # Incoherently stacked envelopes: the two polarisations share the same
    # stratigraphy (envelope correlation 0.984 here) while their complex
    # speckle only correlates at 0.43, so the envelope is what carries the
    # travel-time offset PoRaPy is built to find. Stacking complex traces
    # along-track instead would average the speckle phase to noise.
    y0 = d['along_env'].astype(float)
    y1 = d['across_env'].astype(float)
    dz = float(np.median(np.diff(z)))
    # sample interval as TWO-WAY TIME, which is what a delay in samples means
    dt_s = 2.0 * dz / C_ICE
    win = max(8, int(round(BIN_M / dz)))
    print('%d samples, dz %.3f m, sample dt %.4f ns, bin %d samples (%.0f m)'
          % (z.size, dz, dt_s * 1e9, win, win * dz))

    corrs, inds, delays, _sc, _sd = binned_corrs(
        y0, y1, win, max_shift=MAX_SHIFT, increment=INCREMENT,
        oversample=OVERSAMPLE, surf_guesses=[0.0, 0.0], progress_bar=False)
    centre, off = binned_best_offsets(inds, corrs.real, delays, win)
    peak = np.max(corrs.real, axis=1)
    zc = z[np.clip(centre, 0, z.size - 1)]

    ok = (peak > 0.3) & (zc > 200.0) & (zc < 1450.0)
    print('%d bins, %d pass (peak corr > 0.3, 200-1450 m)'
          % (len(zc), ok.sum()))
    if ok.sum() < 4:
        raise SystemExit('too few usable bins to fit a trend')

    # delay in samples -> two-way delay in seconds, then its depth gradient
    tau = off * dt_s
    A = np.vstack([zc[ok], np.ones(ok.sum())]).T
    slope, icept = np.linalg.lstsq(A, tau[ok], rcond=None)[0]
    dlam_cross = abs(slope) * 2 * np.pi * FC / GRAD_PER_DLAM
    print('\ncrossover method (PoRaPy):')
    print('  d(tau)/dz = %.4g s/m   ->  dlam = %.4f' % (slope, dlam_cross))
    print('  total two-way delay over 200-1450 m: %.2f ns'
          % (abs(slope) * 1250 * 1e9))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 5.0), sharey=True,
                                 layout='constrained')
    a1.plot(peak, zc, '-', color='#888', lw=1.2)
    a1.axvline(0.3, color='#c0392b', ls='--', lw=1.0)
    a1.set_xlabel('peak correlation')
    a1.set_ylabel('depth (m)')
    a1.set_ylim(1500, 0)
    a2.plot(tau * 1e9, zc, 'o', ms=3.5, color='#bbb', label='all bins')
    a2.plot(tau[ok] * 1e9, zc[ok], 'o', ms=4.5, color='#2a78d6',
            label='used')
    zz = np.linspace(200, 1450, 50)
    a2.plot((slope * zz + icept) * 1e9, zz, '-', color='#c0392b', lw=2,
            label=r'fit: $\Delta\lambda$ = %.3f' % dlam_cross)
    a2.set_xlabel('two-way delay between\northogonal polarisations (ns)')
    a2.legend(fontsize=8, loc='lower left')
    for a in (a1, a2):
        a.grid(alpha=0.25)
        for s in ('top', 'right'):
            a.spines[s].set_visible(False)
    fig.suptitle('PoRaPy crossover method on our quad-pol frame '
                 '20250108_02_009', fontsize=11)
    fn = os.path.join(OUT, 'crossover_check.png')
    fig.savefig(fn, dpi=180, facecolor='white')
    print('wrote', fn)
    np.savez(os.path.join(OUT, 'crossover_check.npz'), z=zc, tau=tau,
             peak=peak, slope=slope, dlam=dlam_cross)


if __name__ == '__main__':
    main()

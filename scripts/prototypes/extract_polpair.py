"""Two orthogonal co-polarized traces from a coregistered quad-pol frame.

Run on mem1. Produces the input PoRaPy's crossover method expects - the same
place sounded at two antenna orientations - but WITHOUT needing a crossing,
because the full scattering matrix lets both orientations be synthesized at
every trace. That isolates the travel-time estimator from coregistration
and geolocation differences between two survey lines, which is what we want
when the question is whether the two METHODS agree.

Synthesis is the standard rotation of the co-polarized element,

    T_hh(psi) = c^2 S_hh + c s (S_hv + S_vh) + s^2 S_vv,   c = cos psi

evaluated along the fabric axis and across it, so the pair straddles the
birefringent delay the crossover method is designed to measure.

Writes a small npz (two traces plus geometry) so the crossover analysis
itself can run locally against the shipped PoRaPy code.

Usage: python3 extract_polpair.py <cache.mat> <theta0_ant_deg> <out.npz>
"""
import sys

import h5py
import numpy as np

CACHE = sys.argv[1]
TH_ANT = float(sys.argv[2])     # fabric axis in the ANTENNA frame [deg]
OUT = sys.argv[3]
NLOOK = 51                      # range multilook, samples
NTRACE = 400                    # along-track traces kept


def rd(f, k):
    a = np.array(f[k])
    if a.dtype.names is not None:          # MATLAB complex comes through
        a = a['real'] + 1j * a['imag']     # as a compound dtype
    return a


def main():
    with h5py.File(CACHE, 'r') as f:
        z = np.array(f['z']).ravel()
        hh, vv = rd(f, 'hh'), rd(f, 'vv')
        hv, vh = rd(f, 'hv'), rd(f, 'vh')
    # h5py returns MATLAB arrays transposed; put samples on axis 0
    if hh.shape[0] != z.size:
        hh, vv, hv, vh = hh.T, vv.T, hv.T, vh.T
    print('channels %s, z %d samples %.1f-%.1f m'
          % (hh.shape, z.size, z[0], z[-1]), flush=True)

    out = {}
    for name, psi_deg in (('along', TH_ANT), ('across', TH_ANT + 90.0)):
        c = np.cos(np.deg2rad(psi_deg))
        s = np.sin(np.deg2rad(psi_deg))
        T = c * c * hh + c * s * (hv + vh) + s * s * vv
        # incoherent stack across the frame, then a light range multilook:
        # the crossover method correlates ONE trace per orientation, and a
        # single raw trace is speckle, not stratigraphy
        # Keep a BLOCK of traces rather than one stack. Averaging complex
        # traces along-track destroys the very coherence the delay estimate
        # lives on (speckle phase is random between traces), so the choice
        # of stacking has to be made where it can be tested, not here.
        j0 = max(0, T.shape[1] // 2 - NTRACE // 2)
        Tb = T[:, j0:j0 + NTRACE]
        k = np.ones(NLOOK) / NLOOK
        env = np.apply_along_axis(
            lambda v: np.convolve(v, k, 'same'), 0, np.abs(Tb))
        out[name + '_c'] = Tb.astype(np.complex64)
        out[name + '_env'] = env.mean(axis=1).astype(np.float32)
        print('  %-6s psi_ant %7.2f deg  block %s  |T| median %.4g'
              % (name, psi_deg, Tb.shape, np.median(np.abs(Tb))), flush=True)

    np.savez_compressed(OUT, z=z, theta0_ant=TH_ANT, nlook=NLOOK, **out)
    print('wrote', OUT)


if __name__ == '__main__':
    main()

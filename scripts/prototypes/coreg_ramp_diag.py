#!/usr/bin/env python3
"""How much block coherence does along-track phase-ramp compensation recover?

Within 125-trace blocks (the pipeline's Ridge A block), compare |C| from a
plain average of the single-look HH*conj(VV) products with |C| after
demodulating each product by the block's local along-track phase ramp,
estimated from adjacent-trace phasor products of the depth-multilooked
product (never an unwrap). Per zone and depth band. Also the residual
lateral fringe rate after the linear model (curvature) as cycles per block.
"""
import os
import sys
import numpy as np
import h5py

fn = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    "~/projects/polarimetry_thwaites/data/Data_20240108_01_001.mat")
C_ICE = 3e8 / 1.78 / 2
NB_Z = 25
BLK = 125
ZONES = {"west fan  x 200-1800": (200, 1825),
         "centre    x 3500-5500": (3500, 5500),
         "east      x 7000-9500": (7000, 9500)}
BANDS = [(150, 300), (300, 600), (600, 1000), (1000, 1400)]


def cplx(a):
    return a["real"].astype(np.float64) + 1j * a["imag"].astype(np.float64)


with h5py.File(fn, "r") as f:
    T = np.array(f["Time"]).ravel()
    z = T * C_ICE
    r1 = int(np.searchsorted(z, 1450))
    nz = r1 // NB_Z
    zc = z[NB_Z//2::NB_Z][:nz]
    print("%-22s %-10s %8s %8s %8s   %s"
          % ("zone", "band", "|C|plain", "|C|ramp", "|C|1look",
             "ramp cyc/blk (median |k|)"))
    for zn, (x0, x1) in ZONES.items():
        nb = (x1 - x0) // BLK
        x1 = x0 + nb * BLK
        ref = cplx(f["ref"][x0:x1, :r1]).T
        sec = cplx(f["sec_reg"][x0:x1, :r1]).T
        ref[~np.isfinite(ref)] = 0
        sec[~np.isfinite(sec)] = 0
        prod = ref * np.conj(sec)
        p1 = np.abs(ref)**2
        p2 = np.abs(sec)**2
        P = prod[:nz*NB_Z].reshape(nz, NB_Z, -1).sum(1)
        A = p1[:nz*NB_Z].reshape(nz, NB_Z, -1).sum(1)
        B = p2[:nz*NB_Z].reshape(nz, NB_Z, -1).sum(1)
        Pb = P.reshape(nz, nb, BLK)
        Ab = A.reshape(nz, nb, BLK).sum(2)
        Bb = B.reshape(nz, nb, BLK).sum(2)
        c_plain = np.abs(Pb.sum(2)) / np.sqrt(Ab * Bb + 1e-30)
        # local ramp per block: mean adjacent-trace phasor product (smoothed
        # over 5 depth cells for stability), then demodulate
        ph = Pb / (np.abs(Pb) + 1e-30)
        adj = (ph[:, :, 1:] * np.conj(ph[:, :, :-1])).sum(2)
        k = np.angle(adj)                                   # rad per trace
        ks = np.copy(k)
        for i in range(nz):
            lo = max(0, i-2)
            hi = min(nz, i+3)
            ks[i] = np.angle(np.exp(1j*k[lo:hi]).sum(0))
        x = np.arange(BLK) - BLK/2
        demod = np.exp(-1j * ks[:, :, None] * x[None, None, :])
        c_ramp = np.abs((Pb * demod).sum(2)) / np.sqrt(Ab * Bb + 1e-30)
        c_1 = (np.abs(P) / np.sqrt(A * B + 1e-30)).reshape(nz, nb, BLK).mean(2)
        for lo_, hi_ in BANDS:
            m = (zc >= lo_) & (zc < hi_)
            print("%-22s %4d-%4d %8.3f %8.3f %8.3f   %.2f"
                  % (zn, lo_, hi_, np.nanmedian(c_plain[m]),
                     np.nanmedian(c_ramp[m]), np.nanmedian(c_1[m]),
                     np.nanmedian(np.abs(ks[m])) * BLK / (2*np.pi)))

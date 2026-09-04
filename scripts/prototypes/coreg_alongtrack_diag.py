#!/usr/bin/env python3
"""Is the margin coherence loss ALONG-TRACK phase variation inside the multilook?

For the same product as coreg_lag_diag.py: single-look HH*conj(VV) products,
then coherence with along-track windows of 1, 3, 7, 15, 51 traces (depth
window fixed at 25 bins), per zone. If |C| recovers as the window shrinks,
the fringes are dense along track and averaging kills them; if it does not,
the pair is decorrelated at the single-look level (physical or SNR). Also the
along-track phase gradient (cycles per 100 traces) from adjacent-trace
phasor products, and the trace spacing from the positions.

  python3 scripts/prototypes/coreg_alongtrack_diag.py [Data_<tag>.mat]
"""
import os, sys
import numpy as np, h5py

fn = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    "~/projects/polarimetry_thwaites/data/Data_20240108_01_001.mat")
C_ICE = 3e8 / 1.78 / 2
NB_Z = 25
WINS = [1, 3, 7, 15, 51]
ZONES = {"west fan  x 200-1800": (200, 1800), "centre    x 3500-5500": (3500, 5500), "east      x 7000-9500": (7000, 9500)}
BANDS = [(150, 300), (300, 600), (600, 1000), (1000, 1400), (1400, 1800)]

def cplx(a):
    return a["real"].astype(np.float64) + 1j * a["imag"].astype(np.float64)

with h5py.File(fn, "r") as f:
    T = np.array(f["Time"]).ravel(); z = T * C_ICE
    la = np.array(f["Latitude"]).ravel(); lo = np.array(f["Longitude"]).ravel()
    dph = np.radians(np.diff(la)); dlo = np.radians(np.diff(lo))
    aa = np.sin(dph/2)**2 + np.cos(np.radians(la[:-1]))*np.cos(np.radians(la[1:]))*np.sin(dlo/2)**2
    step = 2*6371e3*np.arcsin(np.sqrt(aa))
    print("trace spacing: median %.2f m, p10 %.2f, p90 %.2f; frame %.1f km" % (np.median(step), np.percentile(step,10), np.percentile(step,90), step.sum()/1e3))
    r1 = int(np.searchsorted(z, 1850))
    zc = z[NB_Z//2::NB_Z][: r1 // NB_Z]
    print("%-22s %-10s " % ("zone", "band") + " ".join("w=%2d" % w for w in WINS) + "   dphi/dx cyc/100tr  |C|_1look(0.5m)")
    for zn, (x0, x1) in ZONES.items():
        ref = cplx(f["ref"][x0:x1, :r1]).T; sec = cplx(f["sec_reg"][x0:x1, :r1]).T
        ref[~np.isfinite(ref)] = 0; sec[~np.isfinite(sec)] = 0
        prod = ref * np.conj(sec); p1 = np.abs(ref)**2; p2 = np.abs(sec)**2
        nz = r1 // NB_Z
        # depth-multilook first (25 bins), keep single-look along track
        P = prod[:nz*NB_Z].reshape(nz, NB_Z, -1).sum(1)
        A = p1[:nz*NB_Z].reshape(nz, NB_Z, -1).sum(1)
        B = p2[:nz*NB_Z].reshape(nz, NB_Z, -1).sum(1)
        # along-track phase gradient from adjacent-trace phasor products of the
        # depth-multilooked product (no unwrap)
        ph = P / (np.abs(P) + 1e-30)
        dphx = np.angle((ph[:, 1:] * np.conj(ph[:, :-1])))
        for lo_, hi_ in BANDS:
            m = (zc >= lo_) & (zc < hi_)
            row = "%-22s %4d-%4d " % (zn, lo_, hi_)
            for w in WINS:
                nx = P.shape[1] // w
                Pw = P[:, :nx*w].reshape(nz, nx, w).sum(2); Aw = A[:, :nx*w].reshape(nz, nx, w).sum(2); Bw = B[:, :nx*w].reshape(nz, nx, w).sum(2)
                c = np.abs(Pw) / np.sqrt(Aw * Bw + 1e-30)
                row += "%5.2f" % np.nanmedian(c[m])
            # circular-mean along-track gradient magnitude in the band
            g = np.abs(np.angle(np.mean(np.exp(1j*dphx[m]), axis=1))) / (2*np.pi) * 100
            snr_proxy = np.nanmedian(A[m]) / np.nanmedian(A[(zc >= 1700)])
            row += "   %6.2f            (HH power / deep-noise %.1f dB)" % (np.nanmedian(g), 10*np.log10(max(snr_proxy, 1e-9)))
            print(row)

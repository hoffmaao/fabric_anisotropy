#!/usr/bin/env python3
"""Diagnose HH-VV coherence loss on a Thwaites margin frame: is it the coregistration?

Reads a shipped CSARP_polarimetric product (ref = HH, sec = VV unregistered,
sec_reg = VV after the toolbox tile coregistration, row_offset = the applied
per-pixel shift) and forms, on a multilook grid:

  |C_raw|, |C_reg|   coherence before and after coregistration
  dphi/dz            depth gradient of the co-pol phase (registered), from
                     adjacent-cell phasor products - never an unwrap
  row_offset         the applied shift, with the fraction outside +-search_t
  lag scan           at a few along-track positions, |gamma(lag)| of HH vs VV
                     in 60 m windows over lags -12..12 bins: the correlation
                     the tile pick actually faces. One peak = a pure-mode pair;
                     peaks at 0 and +-dtau = mode mixing, where an argmax flips.

  python3 scripts/prototypes/coreg_lag_diag.py [Data_<tag>.mat] [out.png]
"""
import os
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'figures'))
import scar_style                            # noqa: E402

fn = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    "~/projects/polarimetry_thwaites/data/Data_20240108_01_001.mat")
tag = os.path.basename(fn)[5:-4]
out = os.path.join(scar_style.out_dir(sys.argv, 2),
                   "coreg_lag_diag_%s.png" % tag)

C_ICE = 3e8 / 1.78 / 2          # m per s of two-way time
NB_Z, NB_X = 25, 15             # multilook cell: bins x traces
Z_MAX = 1900.0
LAGS = np.arange(-12, 13)
WIN_LAG = 214                    # ~60 m of bins for the lag scan
SEARCH_T = 5


def cplx(a):
    return a["real"].astype(np.float64) + 1j * a["imag"].astype(np.float64)


with h5py.File(fn, "r") as f:
    T = np.array(f["Time"]).ravel()
    dt = T[1] - T[0]
    z = T * C_ICE
    Nx, Nt = f["ref"].shape
    r1 = int(np.searchsorted(z, Z_MAX))
    z = z[:r1]
    nz = r1 // NB_Z
    nx = Nx // NB_X
    Craw = np.zeros((nz, nx), complex); Creg = np.zeros((nz, nx), complex)
    P1 = np.zeros((nz, nx)); P2 = np.zeros((nz, nx)); P2r = np.zeros((nz, nx))
    RO = np.zeros((nz, nx)); ROmax = np.zeros((nz, nx))
    # lag scan positions: 6 along-track locations
    xs_lag = np.linspace(0.05, 0.95, 6) * Nx
    lag_pos = [int(x) for x in xs_lag]
    z_lag = np.arange(WIN_LAG // 2, r1 - WIN_LAG // 2, WIN_LAG // 2)
    G = np.full((len(lag_pos), len(z_lag), len(LAGS)), np.nan)
    CH = 300
    for x0 in range(0, nx * NB_X, CH):
        x1 = min(x0 + CH, nx * NB_X)
        ref = cplx(f["ref"][x0:x1, :r1]).T          # [Nt x nxchunk]
        sec = cplx(f["sec"][x0:x1, :r1]).T
        secr = cplx(f["sec_reg"][x0:x1, :r1]).T
        ro = np.array(f["row_offset"][x0:x1, :r1]).T
        ref[~np.isfinite(ref)] = 0; sec[~np.isfinite(sec)] = 0; secr[~np.isfinite(secr)] = 0
        c0, c1 = x0 // NB_X, x1 // NB_X
        for j in range(c0, c1):
            jj = slice(j * NB_X - x0, (j + 1) * NB_X - x0)
            a = ref[:, jj]; b = sec[:, jj]; br = secr[:, jj]; r = ro[:, jj]
            for i in range(nz):
                ii = slice(i * NB_Z, (i + 1) * NB_Z)
                Craw[i, j] = np.sum(a[ii] * np.conj(b[ii]))
                Creg[i, j] = np.sum(a[ii] * np.conj(br[ii]))
                P1[i, j] = np.sum(np.abs(a[ii]) ** 2)
                P2[i, j] = np.sum(np.abs(b[ii]) ** 2)
                P2r[i, j] = np.sum(np.abs(br[ii]) ** 2)
                RO[i, j] = np.nanmedian(r[ii])
                ROmax[i, j] = np.nanmax(np.abs(r[ii]))
        # lag scans for positions inside this chunk
        for p, xp in enumerate(lag_pos):
            if not (x0 <= xp < x1):
                continue
            jj = slice(max(xp - x0 - 7, 0), xp - x0 + 8)
            a = ref[:, jj]; b = sec[:, jj]
            for k, zc in enumerate(z_lag):
                ii = np.arange(zc - WIN_LAG // 2, zc + WIN_LAG // 2)
                pa = np.sum(np.abs(a[ii]) ** 2)
                for l, lag in enumerate(LAGS):
                    jl = ii + lag
                    if jl[0] < 0 or jl[-1] >= r1:
                        continue
                    pb = np.sum(np.abs(b[jl]) ** 2)
                    G[p, k, l] = np.abs(np.sum(a[ii] * np.conj(b[jl]))) / np.sqrt(pa * pb + 1e-30)
        print("chunk %d-%d of %d" % (x0, x1, Nx), flush=True)

coh_raw = np.abs(Craw) / np.sqrt(P1 * P2 + 1e-30)
coh_reg = np.abs(Creg) / np.sqrt(P1 * P2r + 1e-30)
zc = z[NB_Z // 2::NB_Z][:nz]
xc = np.arange(nx) * NB_X * 1.0  # trace index of cell centres
# depth phase gradient of the registered coherence: cycles per 100 m, from
# adjacent-cell phasor products (no unwrap)
ph = Creg / (np.abs(Creg) + 1e-30)
dphi = np.angle(ph[1:] * np.conj(ph[:-1]))               # rad per cell
dz_cell = NB_Z * dt * C_ICE
grad = dphi / (2 * np.pi) / dz_cell * 100.0              # cycles / 100 m
bad_off = ROmax > SEARCH_T

print("frame %s: %d traces, %d bins to %.0f m, cell %.1f m x %d traces" % (tag, Nx, r1, Z_MAX, dz_cell, NB_X))
for lo, hi in [(100, 300), (300, 600), (600, 1000), (1000, 1500), (1500, 1900)]:
    m = (zc >= lo) & (zc < hi)
    print("  %4d-%4d m: |C| raw %.3f reg %.3f | reg<0.3: %3.0f%% | |grad| median %.2f cyc/100m | offset>search_t: %4.1f%% of cells | median |offset| %.2f bins"
          % (lo, hi, np.nanmedian(coh_raw[m]), np.nanmedian(coh_reg[m]), 100 * np.mean(coh_reg[m] < 0.3),
             np.nanmedian(np.abs(grad[m[:-1]])), 100 * np.mean(bad_off[m]), np.nanmedian(np.abs(RO[m]))))
low = coh_reg < 0.3
print("  cells with reg |C|<0.3: %.1f%%; of those, offset beyond search_t in %.1f%% (vs %.1f%% overall)"
      % (100 * low.mean(), 100 * bad_off[low].mean(), 100 * bad_off.mean()))
print("  cells where raw |C| > reg |C| + 0.1 (registration made it worse): %.1f%%" % (100 * np.mean(coh_raw > coh_reg + 0.1)))

fig, ax = plt.subplots(4, 1, figsize=(14, 15), sharex=True)
ext = [xc[0], xc[-1], zc[-1], zc[0]]
im = ax[0].imshow(coh_raw, aspect="auto", extent=ext, vmin=0, vmax=1, cmap="viridis"); ax[0].set_title("|C| raw (HH vs VV unregistered)"); plt.colorbar(im, ax=ax[0])
im = ax[1].imshow(coh_reg, aspect="auto", extent=ext, vmin=0, vmax=1, cmap="viridis"); ax[1].set_title("|C| after toolbox coregistration"); plt.colorbar(im, ax=ax[1])
im = ax[2].imshow(np.clip(RO, -8, 8), aspect="auto", extent=ext, cmap="RdBu_r", vmin=-8, vmax=8); ax[2].set_title("applied row offset (bins, clipped +-8; search was +-%d)" % SEARCH_T); plt.colorbar(im, ax=ax[2])
im = ax[3].imshow(np.clip(np.abs(grad), 0, 3), aspect="auto", extent=[xc[0], xc[-1], zc[-2], zc[0]], cmap="magma", vmin=0, vmax=3); ax[3].set_title("|d phi / dz| of registered co-pol phase (cycles per 100 m)"); plt.colorbar(im, ax=ax[3])
for a in ax:
    a.set_ylabel("depth (m)")
    for xp in lag_pos:
        a.axvline(xp, color="w", lw=0.5, ls="--")
ax[3].set_xlabel("trace")
fig.suptitle("%s: coherence, coregistration offset and phase gradient" % tag)
fig.tight_layout()
fig.savefig(out, dpi=110)

fig, ax = plt.subplots(1, len(lag_pos), figsize=(3 * len(lag_pos), 7), sharey=True)
for p in range(len(lag_pos)):
    im = ax[p].imshow(G[p], aspect="auto", extent=[LAGS[0], LAGS[-1], z[z_lag[-1]], z[z_lag[0]]], vmin=0, vmax=0.8, cmap="viridis")
    ax[p].set_title("trace %d" % lag_pos[p]); ax[p].set_xlabel("lag (bins)")
    ax[p].axvline(-SEARCH_T, color="w", lw=0.5); ax[p].axvline(SEARCH_T, color="w", lw=0.5)
ax[0].set_ylabel("depth (m)")
fig.suptitle("%s: |gamma(lag)| of HH vs VV, 60 m windows (white: toolbox search bound)" % tag)
fig.tight_layout()
fig.savefig(out.replace(".png", "_lags.png"), dpi=110)
np.savez(out.replace(".png", ".npz"), zc=zc, xc=xc, coh_raw=coh_raw, coh_reg=coh_reg, RO=RO, ROmax=ROmax, grad=grad, G=G, lags=LAGS, z_lag=z[z_lag], lag_pos=lag_pos)
print("wrote", out)

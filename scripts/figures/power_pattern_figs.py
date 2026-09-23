#!/usr/bin/env python3
"""Figures for the power-extinction fabric estimator and what the 2D SAR
quad-pol context adds over a rotated point measurement.

Inputs (scratch, produced by power_figdata.m and coreg_lag_diag.py):
  power_figdata.mat                  ridge (009 frame) and thw (250 m blocks)
  coreg_lag_diag_<tag>.npz           coherence / offset cells of the product
  quadpol_section_<tag>_ct.mat       the coherence chain's blocks (server)

  python3 scripts/figures/power_pattern_figs.py <scratch_dir> [out_dir=figs]
"""
import os
import sys
import numpy as np
import h5py
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style  # noqa: E402

SP = sys.argv[1]
OUT = scar_style.out_dir(sys.argv, 2)
os.makedirs(OUT, exist_ok=True)
TAG = "20240108_01_001"
FC = 750e6
DT_NS = 3.3333          # range bin, ns
# rad/m per unit dlam (ptt.constants)
GPD = 2 * np.pi * FC * 0.034 / (np.sqrt(3.15) * 3e8)
KAPPA_MIN = 0.15

d = loadmat(os.path.join(SP, "power_figdata.mat"), squeeze_me=True,
            struct_as_record=False)
R, T = d["ridge"], d["thw"]


def deg(a):
    return np.degrees(a) % 180


# ---------------------------------------------------------------- Ridge A
fig, ax = plt.subplots(1, 3, figsize=(16, 7),
                       gridspec_kw={"width_ratios": [1.6, 1, 1]})
psi = np.degrees(R.psi)
z = R.z
im = ax[0].imshow(np.asarray(R.Pn2, float), aspect="auto",
                  extent=[psi[0], psi[-1], z[-1], z[0]], vmin=0, vmax=2,
                  cmap="magma")
th = deg(np.nanmedian(R.theta0))
for a, ls, lab in [(th, "-", "axis"), (th + 90, "-", ""),
                   ((th + 45) % 180, "--", "nodes"),
                   ((th - 45) % 180, "--", "")]:
    ax[0].axvline(a % 180, color="c", ls=ls, lw=1, label=lab or None)
ax[0].legend(loc="lower right")
ax[0].set_xlabel("synthesized azimuth (deg, geographic)")
ax[0].set_ylabel("depth (m)")
ax[0].set_title("Ridge A 009: co-pol power P(psi, z), normalised per row")
plt.colorbar(im, ax=ax[0], label="P / ((r1^2 + r2^2)/2)")
jn = np.argmin(np.abs(psi - (th + 45) % 180))
ax[1].plot(np.asarray(R.Pn2, float)[:, jn], z, color="k", lw=0.6,
           label="power at node azimuth")
ax[1].invert_yaxis()
ax[1].set_xlim(0, 1.6)
ax[1].set_xlabel("normalised power")
ax[1].set_title("node series (nulls at delta = pi mod 2pi)")
ax[1].legend(loc="lower right")
ax[2].plot(R.ls_dlam, R.ls_zw, "o-", ms=3, color="C0",
           label="coherence LS (held axis)")
ok = np.isfinite(R.dlam)
ax[2].plot(R.dlam[ok], R.zw[ok], "s-", ms=3, color="C3",
           label="power extinction |dlam|")
ax2 = ax[2].twiny()
ax2.plot(R.kappa, R.zw, color="0.5", lw=0.8, label="mode coherence kappa")
ax2.set_xlim(0, 1)
ax2.set_xlabel("kappa")
ax[2].invert_yaxis()
ax[2].set_xlim(0, 0.12)
ax[2].set_xlabel("dlam")
ax[2].set_title("dlam: power vs coherence")
ax[2].legend(loc="lower right")
ax2.legend(loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "power_pattern_ridge_a_009.png"), dpi=120)

# ---------------------------------------------------- Thwaites blocks patterns
fig, ax = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
for k, b in enumerate(np.atleast_1d(T.Pn2_which)):
    im = ax[k].imshow(np.asarray(T.Pn2_blocks[:, :, k], float), aspect="auto",
                      extent=[psi[0], psi[-1], T.z[-1], T.z[0]], vmin=0,
                      vmax=2, cmap="magma")
    thb = deg(np.nanmedian(T.theta0[:, b - 1]))
    for a in [thb, thb + 90]:
        ax[k].axvline(a % 180, color="c", lw=1)
    for a in [thb + 45, thb - 45]:
        ax[k].axvline(a % 180, color="c", lw=1, ls="--")
    ax[k].set_title("block %d, x = %.0f m: kappa %.2f"
                    % (b, T.x[b - 1], np.nanmedian(T.kappa[:, b - 1])))
    ax[k].set_xlabel("azimuth (deg)")
ax[0].set_ylabel("depth (m)")
fig.suptitle("Thwaites %s: power patterns in three 250 m blocks "
             "(solid: axes, dashed: nodes)" % TAG)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "power_pattern_thwaites_blocks.png"), dpi=120)

# -------------------------------------------------------------- margin context
n = np.load(os.path.join(SP, "coreg_lag_diag_%s.npz" % TAG))
zc, xc, coh_raw, RO = n["zc"], n["xc"], n["coh_raw"], n["RO"]
with h5py.File(os.path.join(SP, "quadpol_section_%s_ct.mat" % TAG), "r") as f:
    r = f["res"]
    zs = np.array(r["z"]).ravel()
    sd = np.array(r["sec_dlam_ls"])
    sd = sd if sd.shape[-1] == zs.size else sd.T
    la = np.array(r["sec_lat"]).ravel()
    lo = np.array(r["sec_lon"]).ravel()
    la0 = np.array(r["lat0"]).ravel()[0]
    lo0 = np.array(r["lon0"]).ravel()[0]
# block along-track positions (m from the frame start) from the block centres


def hav(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    dl = np.radians(lo2 - lo1)
    a = (np.sin((p2 - p1) / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2)
    return 2 * 6371e3 * np.arcsin(np.sqrt(a))


xb = hav(la0, lo0, la, lo)
if xb[0] > xb[-1]:
    xb = xb.max() - xb   # orientation guard: the frame start is trace 0
ok_x = np.isfinite(T.x)
X = T.x[ok_x]
K = T.kappa[:, ok_x]
DL = T.dlam[:, ok_x]
Vb = T.V[:, ok_x]


def nanmed_smooth(A, kz=3, kx=3):
    """running nan-median over kz windows x kx blocks: the per-window power
    estimates are far noisier than the coherence estimator's, and a section
    of them is only readable pooled"""
    out = np.full_like(A, np.nan, dtype=float)
    hz, hx = kz // 2, kx // 2
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            blk = A[max(0, i - hz):i + hz + 1, max(0, j - hx):j + hx + 1]
            if np.isfinite(blk).sum() >= 3:
                out[i, j] = np.nanmedian(blk)
    return out


Ks = nanmed_smooth(K, 5, 3)
DLs = nanmed_smooth(np.where(K >= KAPPA_MIN, DL, np.nan), 5, 3)
# predicted birefringent delay: from the coherence chain's own dlam where it
# reports (the pipeline product), integrated in depth; two-way phase over
# 2 pi f0 in range bins
sdm = np.where(np.isfinite(sd), sd, 0.0)
dzs = np.gradient(zs)
# [nb x Nt]
tau_ls = (np.cumsum(sdm * dzs[None, :], axis=1) * GPD / (2 * np.pi * FC)
          * 1e9 / DT_NS)

fig, ax = plt.subplots(4, 1, figsize=(14, 16), sharex=True)
ext_c = [xc[0], xc[-1], zc[-1], zc[0]]
im = ax[0].imshow(coh_raw, aspect="auto", extent=ext_c, vmin=0, vmax=1,
                  cmap="viridis")
plt.colorbar(im, ax=ax[0])
ax[0].set_title("|C| HH-VV (product, raw): where the phase estimators "
                "can work")
ext_p = [X[0], X[-1], T.zw[-1], T.zw[0]]
im = ax[1].imshow(Ks, aspect="auto", extent=ext_p, vmin=0, vmax=1,
                  cmap="viridis")
plt.colorbar(im, ax=ax[1])
ax[1].set_title("node visibility kappa from POWER alone (250 m blocks, "
                "150 m windows, 5x3 running median)")
im = ax[2].imshow(np.asarray(sd, float), aspect="auto",
                  extent=[xb.min(), xb.max(), zs[-1], zs[0]], vmin=0,
                  vmax=0.2, cmap="plasma")
ax[2].set_title("dlam: coherence chain blocks (colour; white = abstained) "
                "and power-extinction |dlam| where kappa >= %.2f "
                "(dots, same scale)" % KAPPA_MIN)
plt.colorbar(im, ax=ax[2])
XX, ZZ = np.meshgrid(X, T.zw)
m = np.isfinite(DLs)
ax[2].scatter(XX[m], ZZ[m], c=DLs[m], s=9, cmap="plasma", vmin=0, vmax=0.2,
              edgecolors="w", linewidths=0.2)
im = ax[3].imshow(np.clip(RO, -6, 6), aspect="auto", extent=ext_c,
                  cmap="RdBu_r", vmin=-6, vmax=6)
plt.colorbar(im, ax=ax[3])
XB, ZS = np.meshgrid(xb, zs)
cs = ax[3].contour(XB, ZS, tau_ls.T, levels=[0.5, 1, 1.5, 2, 3], colors="k",
                   linewidths=0.8)
ax[3].clabel(cs, fmt="%.1f", fontsize=7)
ax[3].set_title("tile coregistration offset measured by the toolbox "
                "(colour, bins) vs birefringent delay predicted from the "
                "coherence-chain dlam (contours, bins)")
for a in ax:
    a.set_ylabel("depth (m)")
    a.set_ylim(1500, 50)
ax[3].set_xlabel("along-track (m)")
fig.suptitle("Thwaites %s: the 2D SAR quad-pol context - coherence, "
             "power-only fabric, and model-predicted vs measured "
             "registration" % TAG)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "coreg_context_thwaites_%s.png" % TAG), dpi=110)
print("wrote figures to", OUT)

"""Delta-k (split-spectrum) traveltime-difference estimation, no unwrapping.

Thwaites eastern-margin frame 20240108_01: the polarimetric delay dtau
is non-dispersive, so the HHxVV* interferogram phase is linear in RF
frequency, phi(f) = s*2*pi*f*dtau. Splitting the range spectrum into
sub-bands turns the fringe-integer problem into per-pixel measurements:

  stage A: pulse-pair across M narrow sub-bands (adjacent df ~ 14 MHz,
           unambiguous |dtau| < ~36 ns)
  stage B: half-band delta-k (df ~ 80 MHz), integers resolved by A
  stage C: full-carrier interferogram (fc = 750 MHz), fringe integer
           resolved by B -> interferometric precision, absolute dtau

Validated against the SNAPHU chain, the coregistration offsets, and the
production block dtau (fabric_batch). The spectral orientation (sign of
the baseband->RF mapping) is fixed empirically by regression against
the coregistration offsets, which settles the frame's phase sign.
"""
import os

import h5py
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from scipy.io import loadmat

FN = os.path.expanduser(
    '~/projects/polarimetry_thwaites/data/Data_20240108_01_001.mat')
PROD = os.path.expanduser(
    '~/data/opr/fabric_batch/joint/2023_Antarctica_Ground/20240108_01/'
    'Data_20240108_01_001.mat')
OUT = os.path.expanduser('~/projects/fabric_anisotropy/figs')

FC = 750e6
M_SUB = 12          # narrow sub-bands for stage A
TGRP = 30           # analysis grid: fast-time samples per cell (100 ns)
XGRP = 31           # traces per cell (~31 m)


def to_complex(a):
    return a['real'].astype(np.float32) + 1j * a['imag'].astype(np.float32)


def cell_mean(arr, tg, xg):
    """Mean over (xg trace, tg time) cells; arr is (traces, time)."""
    nx = arr.shape[0] // xg * xg
    nt = arr.shape[1] // tg * tg
    a = arr[:nx, :nt].reshape(nx // xg, xg, nt // tg, tg)
    return a.mean(axis=(1, 3))


def csmooth(z, size=(5, 5)):
    """Complex boxcar smoothing (for ladder resolver fields only)."""
    from scipy.ndimage import uniform_filter
    return uniform_filter(z.real, size) + 1j * uniform_filter(z.imag, size)


def main():
    with h5py.File(FN) as f:
        t = f['Time'][0, :]
        lats = f['Latitude'][:, 0]
        lons = f['Longitude'][:, 0]
        surf = f['Surface'][:, 0]
        ref = to_complex(f['ref'][()])
        sec_raw = to_complex(f['sec'][()])
        sec = to_complex(f['sec_reg'][()])
        snaphu = f['snaphu_out_phase'][()]
        coh = np.abs(to_complex(f['interferogram_coherence'][()]))
        row_off = f['row_offset'][()]
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    Nx, Nt = ref.shape
    print(f'{Nx} traces x {Nt} samples, fs {fs/1e6:.0f} MHz')

    # ---- range spectra
    R = fft(ref, axis=1, workers=8)
    # Coarse rungs use the UNREGISTERED secondary: coregistration shifts
    # the envelope, which removes the group delay that delta-k measures.
    # Narrow sub-bands tolerate the few-bin misregistration (envelope
    # correlation width 1/B_sub >> offsets). The fine rung (full-band
    # carrier phase) is registration-invariant and uses sec_reg.
    S = fft(sec_raw, axis=1, workers=8)
    S_REG = fft(sec, axis=1, workers=8)
    del ref, sec, sec_raw
    pw = (np.abs(R[::37]) ** 2).mean(axis=0)
    pw = np.fft.fftshift(pw)
    fbb = (np.arange(Nt) - Nt // 2) * fs / Nt
    cum = np.cumsum(pw) / pw.sum()
    b0, b1 = np.searchsorted(cum, [0.02, 0.98])
    print(f'band support: baseband {fbb[b0]/1e6:.0f}..{fbb[b1]/1e6:.0f} MHz '
          f'({(fbb[b1]-fbb[b0])/1e6:.0f} MHz wide)')

    n_cell_t = Nt // TGRP
    t_common = t[0] + (np.arange(n_cell_t) + 0.5) * TGRP * dt

    def band_product(i0, i1, SPEC=None):
        """Multilooked sub-band interferogram on the analysis grid.

        Slices shifted-spectrum bins [i0:i1), inverse-transforms the
        band-limited envelopes (natural decimation), forms sec*conj(ref),
        and averages onto the (XGRP x 100 ns) grid. Also returns the
        power-weighted RF center frequency of the slice.
        """
        sl = np.arange(i0, i1)
        unshift = (sl - Nt // 2) % Nt   # back to unshifted bin indices
        if SPEC is None:
            SPEC = S
        r = ifft(R[:, unshift], axis=1, workers=8)
        s = ifft(SPEC[:, unshift], axis=1, workers=8)
        I = s * np.conj(r)
        n = I.shape[1]
        tg = max(1, round(TGRP * dt / (Nt * dt / n)))  # cells of ~100 ns
        Ic = cell_mean(I, tg, XGRP)
        # this band's cell-center time axis, then resample onto the
        # common 100 ns analysis axis so every stage shares one grid
        dt_b = Nt * dt / n
        t_b = t[0] + (np.arange(Ic.shape[1]) + 0.5) * tg * dt_b
        Ir = np.empty((Ic.shape[0], t_common.size), np.complex64)
        for i in range(Ic.shape[0]):
            Ir[i] = np.interp(t_common, t_b, Ic[i].real) + \
                1j * np.interp(t_common, t_b, Ic[i].imag)
        f_c = FC + np.average(fbb[i0:i1], weights=pw[i0:i1])
        return Ir, f_c

    # ---- stage A: M narrow sub-bands, adjacent-pair pulse-pair
    edges = np.linspace(b0, b1, M_SUB + 1).astype(int)
    subs, fcs = [], []
    for m in range(M_SUB):
        Ic, f_c = band_product(edges[m], edges[m + 1])
        subs.append(Ic)
        fcs.append(f_c)
    fcs = np.array(fcs)
    nt_c = subs[0].shape[1]
    nx_c = subs[0].shape[0]
    subs = np.array(subs)
    wsub = np.array([np.abs(s).mean() for s in subs])
    acc = np.zeros((nx_c, nt_c), np.complex64)
    dfs = np.diff(fcs)
    for m in range(M_SUB - 1):
        w = min(wsub[m], wsub[m + 1])
        acc += w * subs[m + 1] * np.conj(subs[m])
    df_A = np.average(dfs, weights=[min(wsub[m], wsub[m + 1])
                                    for m in range(M_SUB - 1)])
    tau_A = np.angle(csmooth(acc)) / (2 * np.pi * df_A)
    print(f'stage A: {M_SUB} sub-bands, mean adjacent df {df_A/1e6:.1f} MHz, '
          f'unambiguous +/-{1e9/(2*df_A):.0f} ns')

    # ---- stage Q: quarter-band rung between A and B
    qedges = np.linspace(b0, b1, 5).astype(int)
    qsubs, qfcs = [], []
    for m in range(4):
        Iq, f_q = band_product(qedges[m], qedges[m + 1])
        qsubs.append(Iq)
        qfcs.append(f_q)
    acc_q = np.zeros_like(acc)
    for m in range(3):
        acc_q += qsubs[m + 1] * np.conj(qsubs[m])
    df_Q = np.mean(np.diff(qfcs))
    tau_Qw = np.angle(csmooth(acc_q)) / (2 * np.pi * df_Q)
    print(f'stage Q: df {df_Q/1e6:.1f} MHz, unambiguous +/-{1e9/(2*df_Q):.1f} ns')

    # ---- stage B: half-band delta-k
    mid = (b0 + b1) // 2
    I_L, f_L = band_product(b0, mid)
    I_U, f_U = band_product(mid, b1)
    ntB = nt_c
    df_B = f_U - f_L
    tau_Bw = np.angle(csmooth(I_U[:, :ntB] * np.conj(I_L[:, :ntB]), (3, 3))) / (2 * np.pi * df_B)
    print(f'stage B: df {df_B/1e6:.1f} MHz, unambiguous +/-{1e9/(2*df_B):.1f} ns')

    # ---- stage C: full-band product on the same grid
    I_full, f_full = band_product(b0, b1, SPEC=S_REG)
    tau_Cw = np.angle(I_full[:, :ntB]) / (2 * np.pi * f_full)

    # analysis-grid axes
    t_cell = t_common
    x_cell = np.arange(nx_c) * XGRP + XGRP // 2
    R_E = 6371.0
    p1, p2 = np.radians(lats[:-1]), np.radians(lats[1:])
    dl = np.radians(lons[1:] - lons[:-1])
    aa = np.sin((p2-p1)/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    dist_tr = np.concatenate([[0], np.cumsum(2*R_E*np.arcsin(np.sqrt(aa)))])
    dist = dist_tr[np.clip(x_cell, 0, Nx - 1)]

    # comparison fields on the same grid
    coh_c = cell_mean(coh, TGRP, XGRP)[:, :ntB]
    snap_c = cell_mean(snaphu, TGRP, XGRP)[:, :ntB]
    coreg_c = cell_mean(row_off * dt, TGRP, XGRP)[:, :ntB]
    surf_c = np.array([np.nanmedian(surf[i*XGRP:(i+1)*XGRP])
                       for i in range(nx_c)])

    # ---- spectral orientation / sign from the coregistration
    def surface_ref(A):
        """Reference each column to a robust median over a shallow band
        (150-600 ns below the surface): single-cell references stamp
        their noise down the whole column as stripes."""
        out = A.copy()
        for i in range(nx_c):
            band = (t_cell > surf_c[i] + 150e-9) & (t_cell < surf_c[i] + 600e-9)
            out[i] -= np.nanmedian(A[i, band])
        return out

    tau_A_r = surface_ref(tau_A[:, :ntB])
    coreg_r = surface_ref(coreg_c)
    m = (coh_c > 0.35) & np.isfinite(coreg_r) & (np.abs(coreg_r) > dt / 2)
    sign = np.sign(np.nansum(tau_A_r[m] * coreg_r[m]))
    print(f'spectral orientation vs coregistration: sign {sign:+.0f} '
          f'(n={m.sum()} cells) -> frame phase_sign should be {sign:+.0f} '
          f'* (production used +1)')
    if sign < 0:
        tau_A, tau_Qw, tau_Bw, tau_Cw = -tau_A, -tau_Qw, -tau_Bw, -tau_Cw
        tau_A_r = surface_ref(tau_A[:, :ntB])

    # ---- ladder: resolve integers rung by rung (all per-pixel)
    n_Q = np.round((tau_A[:, :ntB] - tau_Qw[:, :ntB]) * df_Q)
    tau_Q = tau_Qw[:, :ntB] + n_Q / df_Q
    n_B = np.round((tau_Q - tau_Bw) * df_B)
    tau_B = tau_Bw + n_B / df_B
    # carrier-sign diagnostic: the correct sign should round cleanly in
    # high-coherence cells; a wrong sign leaves margins ~uniform (median .25)
    hi = coh_c > 0.5
    for sgn, nm in [(1, 'assumed'), (-1, 'flipped')]:
        arr = (tau_B - sgn * tau_Cw) * f_full
        fr = np.abs(arr - np.round(arr))[hi]
        print(f'  carrier sign {nm}: high-coh rounding margin median '
              f'{np.nanmedian(fr):.3f}, 90th {np.nanpercentile(fr, 90):.3f}')
    n_C = np.round((tau_B - tau_Cw) * f_full)
    tau_C = tau_Cw + n_C / f_full
    for nm, arr, q in [('Q', (tau_A[:, :ntB] - tau_Qw[:, :ntB]) * df_Q, 1),
                       ('B', (tau_Q - tau_Bw) * df_B, 1),
                       ('C', (tau_B - tau_Cw) * f_full, 1)]:
        fr = np.abs(arr - np.round(arr))
        print(f'  rung {nm}: rounding margin (want << 0.5): '
              f'median {np.nanmedian(fr):.2f}, 90th pct {np.nanpercentile(fr, 90):.2f}')

    tau_C_r = surface_ref(tau_C)
    snap_tau_r = surface_ref(snap_c / (2 * np.pi * FC))

    # agreement stats
    d_sn = (tau_C_r - snap_tau_r) * FC   # in fringes at fc
    ok = coh_c > 0.3
    print(f'delta-k vs SNAPHU chain (coh>0.3): median |diff| '
          f'{np.nanmedian(np.abs(d_sn[ok])):.2f} fringes; '
          f'{np.nanmean(np.abs(d_sn[ok]) < 0.5)*100:.0f}% within half a fringe')
    d_cg = (tau_C_r - coreg_r)
    print(f'delta-k vs coreg (strong cells): median |diff| '
          f'{np.nanmedian(np.abs(d_cg[m]))*1e9:.2f} ns')
    print(f'delta-k vs SNAPHU chain in ns: median |diff| '
          f'{np.nanmedian(np.abs(d_sn[ok]))/FC*1e9:.2f} ns; '
          f'SNAPHU chain vs coreg median |diff| '
          f'{np.nanmedian(np.abs((snap_tau_r - coreg_r)[m]))*1e9:.2f} ns')
    sA = np.nanstd((tau_A_r - tau_C_r)[ok])
    print(f'stage A scatter about final: {sA*1e9:.2f} ns '
          f'(need << {1e9/(2*df_B):.1f} ns for a clean B ladder)')

    # production block dtau for profile comparison
    prod = loadmat(PROD)
    dtau_blk = prod['dtau_blk']
    blk_n = dtau_blk.shape[1]

    # ---- figure
    fig = plt.figure(figsize=(16, 11), layout='constrained')
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(fbb / 1e6 + FC / 1e6, pw / pw.max(), 'k-', lw=1)
    for e in edges:
        ax.axvline(FC / 1e6 + fbb[e] / 1e6, color='tab:orange', lw=0.6)
    ax.axvline(FC / 1e6 + fbb[mid] / 1e6, color='crimson', lw=1.5)
    ax.set_xlabel('RF frequency (MHz)')
    ax.set_ylabel('mean power (norm)')
    ax.set_title(f'range spectrum + {M_SUB} sub-bands (red: half-band split)',
                 fontsize=11)
    ax.grid(alpha=0.3)

    tmask = np.where(~np.isfinite(coh_c), np.nan, 1.0)
    for i in range(nx_c):   # mask direct-arrival cells
        tmask[i, t_cell < surf_c[i] + 150e-9] = np.nan
    v = 8
    for k, (fld, ttl) in enumerate([
            (tau_C_r * 1e9, 'delta-k dtau (ns), no unwrapping'),
            (snap_tau_r * 1e9, 'SNAPHU-chain dtau (ns)')]):
        ax = fig.add_subplot(gs[0, 1 + k])
        pc = ax.pcolormesh(dist, t_cell * 1e6, (fld * tmask).T, cmap='RdBu_r',
                           vmin=-v, vmax=v, shading='auto', rasterized=True)
        ax.set_ylim(t_cell[-1] * 1e6, 0)
        ax.set_xlabel('distance (km)')
        if k == 0:
            ax.set_ylabel('TWTT (us)')
        ax.set_title(ttl, fontsize=11)
        fig.colorbar(pc, ax=ax, pad=0.02, label='dtau (ns)')

    ax = fig.add_subplot(gs[1, 0])
    pc = ax.pcolormesh(dist, t_cell * 1e6, (d_sn * tmask).T, cmap='RdBu_r',
                       vmin=-4, vmax=4, shading='auto', rasterized=True)
    ax.set_ylim(t_cell[-1] * 1e6, 0)
    ax.set_xlabel('distance (km)')
    ax.set_ylabel('TWTT (us)')
    ax.set_title('difference in fringes at fc\n(integer plateaus = SNAPHU '
                 'region errors)', fontsize=11)
    fig.colorbar(pc, ax=ax, pad=0.02, label='fringes')

    # block profiles
    axp = fig.add_subplot(gs[1, 1:])
    blocks = [3, 5, 8]
    colors = ['tab:blue', 'tab:red', 'tab:green']
    edges_x = np.linspace(0, nx_c, blk_n + 1).astype(int)
    for b, c in zip(blocks, colors):
        cells = np.arange(edges_x[b], edges_x[b + 1])
        w = np.where((coh_c[cells] > 0.25) & np.isfinite(tmask[cells]),
                     coh_c[cells] ** 2, 0)
        prof = np.nansum(w * tau_C_r[cells], axis=0) / np.nansum(w, axis=0)
        axp.plot(prof * 1e9, t_cell * 1e6, color=c, lw=2,
                 label=f'block {b}: delta-k')
        snp = np.nansum(w * snap_tau_r[cells], axis=0) / np.nansum(w, axis=0)
        axp.plot(snp * 1e9, t_cell * 1e6, color=c, lw=1.2, ls='--',
                 label=f'block {b}: SNAPHU chain')
        axp.plot(1e9 * dtau_blk[::TGRP, b][:len(t_cell)], t_cell * 1e6,
                 color=c, lw=1.0, ls=':', alpha=0.8,
                 label=f'block {b}: production dtau_blk')
    axp.set_ylim(t_cell[-1] * 1e6, 0)
    axp.set_xlim(-14, 14)
    axp.set_xlabel('dtau (ns)')
    axp.set_ylabel('TWTT (us)')
    axp.grid(alpha=0.3)
    axp.legend(fontsize=8, ncol=3)
    axp.set_title('block-averaged dtau profiles: ladder vs unwrap chain vs '
                  'production', fontsize=11)

    fig.suptitle('Delta-k split-spectrum dtau for Thwaites 20240108_01: '
                 'per-pixel absolute traveltime difference, no phase '
                 'unwrapping', fontsize=14)
    out = os.path.join(OUT, 'deltak_thwaites_20240108_01.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

"""Stage-A remedy prototype: deramped coherent multilook for delta-k.

PARKED - THIS REMEDY DOES NOT RESOLVE THE SUPPRESSION. Read the outcome
at the bottom of this docstring before treating anything here as working,
and do NOT port it into ptt.deltakTraveltime on the strength of these
numbers.

The per-stage diagnostic (run_deltak_stages.m, Job 2) located the Ridge A
delta-k suppression AT STAGE A: all three ladder rungs miss the smooth
~5 ns dtau ramp at 10-22 us that SNAPHU and the coregistration offsets
both see, with clean rounding margins - so the ladder is fine and the
loss happens before it, in how the per-rung cross products are formed.

Mechanism: the shipped estimator forms cross products PER PIXEL (the
ordering the Thwaites margin fringe bundle forced: sub-band products
carry the full-carrier fringe rate, so cell-averaging them first
decorrelates where fringes are dense). But per-pixel cross products
SQUARE the noise; at Ridge A depth the interferogram SNR per pixel is
low and the multilook after differencing cannot recover the signal -
the angle of the noise-dominated average sits near zero, which is
exactly the observed bias toward zero.

Remedy prototyped here: deramp each sub-band interferogram by a common
smoothed full-band interferogram phase, cell-average the deramped
products COHERENTLY (the multilook gain now happens before the noise is
squared), and only then cross-multiply the cell-level phasors:

    I_m   = ifft(S_m) . conj(ifft(R_m))        per-pixel, band m
    J_m   = I_m . exp(-i phi_ref)              phi_ref: smoothed angle of
                                               the full-band interferogram
    M_m   = cellavg(J_m)                       coherent multilook
    acc  += M_{m+1} . conj(M_m)                cross products at cell level

The deramp phase is COMMON to both factors of every cross product, so it
cancels exactly and cannot bias dtau. Its only job is to keep the
within-cell phase of J_m slow: where fringes are real, phi_ref tracks
them and the residual within-cell rate is 2*pi*(f_m - fc)*dtau (sub-band
offset scale, ~100x slower than the carrier); where the interferogram is
noise, phi_ref is a smooth random field, effectively constant within a
cell, and cancels in the cross product. Both regimes are safe - no
gating needed, which is what makes this preferable to picking between
average-then-difference and difference-then-average per cell.

Run on the local Ridge A SLC subset (Data_20250108_02_009) against
SNAPHU, the coregistration offsets, and the joint block product. The
success criterion was that the remedied stage B recover the deep ramp -
band std within ~2x of SNAPHU below 10 us - without breaking the shallow
agreement.

OUTCOME: it does not. The deramped coherent multilook does tighten the
ladder rounding margins and it does raise the 5-10 us band std, to 0.89 ns
against 0.46 ns for the shipped chain and 0.49 ns for SNAPHU. But 10-20 us
stays at 0.56 ns against SNAPHU's 1.16 ns, and - the part that actually
blocks this - BOTH chains still ANTI-CORRELATE with SNAPHU over 5-22 us:
corr -0.61 remedied, -0.75 shipped. So the stage-A suppression is not
resolved, and the anti-correlation is unexplained: a remedy that recovered
the missing amplitude would have to correlate with the estimator it is
being judged against, and neither does. Recovering more variance while
still pointing the wrong way is not evidence the mechanism above is the
right diagnosis.

The investigation is therefore PARKED pending an explanation of the
anti-correlation, not concluded. The mechanism and rationale above are
kept because they are the hypothesis under test; they are not a result.

Usage: /opt/anaconda3/bin/python deltak_remedy.py [out_dir]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt              # noqa: E402
from scipy.fft import fft, ifft              # noqa: E402
from scipy.ndimage import uniform_filter     # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'figures'))
import scar_style                            # noqa: E402

ROOT = os.path.expanduser('~/projects/fabric_anisotropy')
FRAME = os.path.expanduser('~/data/opr/stage_frame_20250108_02_001.mat')
OUT = scar_style.out_dir(sys.argv, 1, os.path.join(ROOT, 'figs'))

FC = 750e6          # run_deltak_stages.m convention
S_SIGN = -1         # forced phase_sign (settled 3 Aug)
N_SUB = 12
CELL_TWTT = 100e-9
# traces per cell; the extract is trace-decimated 2x, so this matches the
# production 31-trace ground span
CELL_NTR = 16
SMOOTH = 5          # resolver smoothing on the cell grid
BAND_FRAC = (0.02, 0.98)
REF_BAND = (150e-9, 600e-9)
DERAMP_WIN = (15, 15)   # (fast-time bins, traces) for phi_ref smoothing


def cellavg(a, tg, xg, ntc, nxc):
    """Mean over (tg x xg) cells; a is (Nt, Nx) -> (ntc, nxc)."""
    return a[:ntc * tg, :nxc * xg].reshape(ntc, tg, nxc, xg).mean(axis=(1, 3))


def csmooth(z, n):
    return uniform_filter(z.real, n) + 1j * uniform_filter(z.imag, n)


def band_ref(tau, t_cell, surf_c):
    """Subtract the shallow-band median per cell column (NaN if empty)."""
    out = np.full_like(tau, np.nan)
    for i in range(tau.shape[1]):
        band = (t_cell > surf_c[i] + REF_BAND[0]) & \
               (t_cell < surf_c[i] + REF_BAND[1])
        v = tau[band, i]
        v = v[np.isfinite(v)]
        if v.size >= 2:
            out[:, i] = tau[:, i] - np.median(v)
    return out


def ladder(tau_res, tau_w, df):
    """Resolve tau_w's integer ambiguity from tau_res; return tau, margin."""
    n = np.round((tau_res - tau_w) * df)
    return tau_w + n / df, np.abs((tau_res - tau_w) * df - n)


def profile(tau_ref):
    """Median across cell columns."""
    with np.errstate(all='ignore'):
        return np.nanmedian(tau_ref, axis=1)


def main():
    import h5py

    def to_complex(a):
        return a['real'].astype(np.float32) + 1j * a['imag'].astype(np.float32)

    # MATLAB -v7.3 stores each (Nt x Nx) matrix transposed in HDF5
    with h5py.File(FRAME) as f:
        t = f['Time'][()].ravel()
        surf = f['Surface'][()].ravel()
        R0 = to_complex(f['ref'][()]).T
        S0 = to_complex(f['sec'][()]).T
        row_off = f['row_offset'][()].astype(np.float64).T
        snaphu = f['snaphu_out_phase'][()].astype(np.float64).T
    Nt, Nx = R0.shape
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    print('%d samples x %d traces, fs %.0f MHz' % (Nt, Nx, fs / 1e6))

    bad = ~np.isfinite(R0) | ~np.isfinite(S0)
    if bad.any():
        print('zeroing %d non-finite SLC samples' % bad.sum())
        R0[bad] = 0
        S0[bad] = 0

    tg = max(1, round(CELL_TWTT / dt))
    xg = CELL_NTR
    ntc, nxc = Nt // tg, Nx // xg
    t_cell = t[0] + (np.arange(ntc) + 0.5) * tg * dt
    surf_c = np.array([np.nanmedian(surf[i * xg:(i + 1) * xg])
                       for i in range(nxc)])

    # ---- band support
    R = fft(R0, axis=0, workers=8)
    S = fft(S0, axis=0, workers=8)
    pw = np.fft.fftshift(np.mean(np.abs(R[:, ::max(1, Nx // 200)]) ** 2, 1))
    fbb = (np.arange(Nt) - Nt // 2) * fs / Nt
    cum = np.cumsum(pw) / pw.sum()
    b0 = int(np.searchsorted(cum, BAND_FRAC[0]))
    b1 = int(np.searchsorted(cum, BAND_FRAC[1]))
    print('band support %.0f..%.0f MHz baseband (%.0f MHz)'
          % (fbb[b0] / 1e6, fbb[b1] / 1e6, (fbb[b1] - fbb[b0]) / 1e6))

    # ---- deramp reference: smoothed full-band interferogram phase
    I_full = S0 * np.conj(R0)
    ph = csmooth(I_full, DERAMP_WIN)
    mag = np.abs(ph)
    unit = np.ones_like(ph)
    nz = mag > 0
    unit[nz] = ph[nz] / mag[nz]
    del I_full, ph, mag

    def stage_cross(nband):
        """Both orderings of one ladder rung.

        E1 (shipped): per-pixel adjacent cross products, then cellavg.
        E2 (remedy):  deramp, cellavg per band, then cross the cell means.
        Returns (acc1, acc2, df) on the cell grid.
        """
        edges = np.round(np.linspace(b0, b1, nband + 1)).astype(int)
        acc1 = np.zeros((ntc, nxc), np.complex128)
        acc2 = np.zeros((ntc, nxc), np.complex128)
        prev1 = prev2 = None
        fcs = []
        for m in range(nband):
            sl = slice(edges[m], edges[m + 1])
            un = (np.arange(edges[m], edges[m + 1]) - Nt // 2) % Nt
            Rb = np.zeros((Nt, Nx), np.complex64)
            Sb = np.zeros((Nt, Nx), np.complex64)
            Rb[un] = R[un]
            Sb[un] = S[un]
            Ib = (ifft(Sb, axis=0, workers=8)
                  * np.conj(ifft(Rb, axis=0, workers=8)))
            fcs.append(FC + np.sum(fbb[sl] * pw[sl]) / np.sum(pw[sl]))
            Mb = cellavg(Ib * np.conj(unit), tg, xg, ntc, nxc)
            if prev1 is not None:
                acc1 += cellavg(Ib * np.conj(prev1), tg, xg, ntc, nxc)
                acc2 += Mb * np.conj(prev2)
            prev1, prev2 = Ib, Mb
        return acc1, acc2, float(np.mean(np.diff(fcs)))

    accA1, accA2, dfA = stage_cross(N_SUB)
    accQ1, accQ2, dfQ = stage_cross(4)
    accB1, accB2, dfB = stage_cross(2)
    print('df A/Q/B = %.1f/%.1f/%.1f MHz' % (dfA / 1e6, dfQ / 1e6, dfB / 1e6))

    chains = {}
    for tag, (aA, aQ, aB) in {
            'E1 shipped': (accA1, accQ1, accB1),
            'E2 remedy': (accA2, accQ2, accB2)}.items():
        tau_A = S_SIGN * np.angle(csmooth(aA, SMOOTH)) / (2 * np.pi * dfA)
        tau_Qw = S_SIGN * np.angle(csmooth(aQ, SMOOTH)) / (2 * np.pi * dfQ)
        tau_Bw = S_SIGN * np.angle(csmooth(aB, 3)) / (2 * np.pi * dfB)
        tau_Q, mQ = ladder(tau_A, tau_Qw, dfQ)
        tau_B, mB = ladder(tau_Q, tau_Bw, dfB)
        print('%s: margins Q %.2f/%.2f B %.2f/%.2f (med/p90)' % (
            tag, np.nanmedian(mQ), np.nanpercentile(mQ, 90),
            np.nanmedian(mB), np.nanpercentile(mB, 90)))
        chains[tag] = {
            'A': band_ref(tau_A, t_cell, surf_c),
            'Q': band_ref(tau_Q, t_cell, surf_c),
            'B': band_ref(tau_B, t_cell, surf_c)}

    # ---- reference estimators on the same cell grid + referencing
    tau_sn = band_ref(cellavg(S_SIGN * snaphu / (2 * np.pi * FC),
                              tg, xg, ntc, nxc), t_cell, surf_c)
    tau_cg = band_ref(cellavg(row_off * dt, tg, xg, ntc, nxc),
                      t_cell, surf_c)

    # ---- profiles and the depth-band amplitude table
    tb = t_cell - np.nanmean(surf)
    profs = [
        ('coreg', tau_cg, '#8172b3', ':', 1.2),
        ('snaphu', tau_sn, '#55a868', '--', 1.8),
        ('E1 A', chains['E1 shipped']['A'], '#9ecae1', '-', 1.0),
        ('E1 B (shipped)', chains['E1 shipped']['B'], '#c44e52', '-', 1.6),
        ('E2 A', chains['E2 remedy']['A'], '#4c72b0', '-', 1.0),
        ('E2 B (remedy)', chains['E2 remedy']['B'], '#000000', '-', 2.0),
    ]
    print('\nprofile std by TWTT band below surface (ns)')
    bands = [(0, 2), (2, 5), (5, 10), (10, 20)]
    print('  %-10s' % 'band(us)' + ''.join('%16s' % n for n, *_ in profs))
    for lo, hi in bands:
        row = '  %-10s' % ('%g-%g' % (lo, hi))
        sel = (tb >= lo * 1e-6) & (tb < hi * 1e-6)
        for name, tau, *_ in profs:
            v = profile(tau)[sel]
            v = v[np.isfinite(v)]
            row += '%16s' % ('%.3f' % (1e9 * v.std()) if v.size > 1 else '-')
        print(row)

    # correlation of each chain's B against snaphu over the deep band
    sel = (tb >= 5e-6) & (tb < 22e-6)
    psn = profile(tau_sn)
    for tag in ('E1 shipped', 'E2 remedy'):
        pb = profile(chains[tag]['B'])
        ok = sel & np.isfinite(pb) & np.isfinite(psn)
        if ok.sum() > 10:
            c = np.corrcoef(pb[ok], psn[ok])[0, 1]
            g = np.polyfit(psn[ok], pb[ok], 1)[0]
            print('%s vs snaphu 5-22us: corr %.3f, gain %.2f' % (tag, c, g))

    # ---- figure
    fig, ax = plt.subplots(figsize=(7, 9))
    for name, tau, c, ls, lw in profs:
        ax.plot(1e9 * profile(tau), 1e6 * tb, ls, color=c, lw=lw, label=name)
    ax.invert_yaxis()
    ax.set_xlabel(r'$\Delta\tau$ (ns)')
    ax.set_ylabel(r'TWTT below surface ($\mu$s)')
    ax.grid(alpha=0.3)
    ax.axvline(0, color='0.5', lw=0.8)
    ax.legend(loc='lower left', fontsize=8)
    ax.set_title('Delta-k stage-A remedy, Ridge A 20250108_02_001')
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'deltak_remedy_20250108_02_001.png')
    fig.savefig(out, dpi=160, bbox_inches='tight')
    print('wrote %s' % out)


if __name__ == '__main__':
    main()

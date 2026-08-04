"""v2: apply production phase_sign; DC-guarded FFT peak; sign cross-check
against coregistration offsets where available."""
import os
import numpy as np
from scipy.io import loadmat

FC = 750e6

def pulse_pair(cplx, coh, rows, cols, dt):
    z = cplx[np.ix_(rows, cols)]
    w = np.nan_to_num(coh[np.ix_(rows, cols)]).astype(float) ** 2
    pp = z[1:, :] * np.conj(z[:-1, :])
    ww = w[1:, :] * w[:-1, :]
    acc = np.nansum(ww * pp)
    n = np.nansum(ww)
    if n <= 0 or not np.isfinite(acc):
        return np.nan
    return np.angle(acc) / dt / (2 * np.pi * FC)

def fft_peak(cplx, coh, rows, cols, dt, pad=8):
    z = np.nan_to_num(cplx[np.ix_(rows, cols)])
    w = np.nan_to_num(coh[np.ix_(rows, cols)]).astype(float) ** 2
    n = len(rows)
    win = np.hanning(n)[:, None]
    spec = (np.abs(np.fft.fft(z * w * win, n=pad * n, axis=0)) ** 2).sum(axis=1)
    f = np.fft.fftfreq(pad * n, d=dt)
    guard = np.abs(f) < 1.5 / (n * dt)   # skip DC/envelope leakage
    spec_g = spec.copy(); spec_g[guard] = 0
    k = np.argmax(spec_g)
    return f[k] / FC, spec_g[k] / np.median(spec)

def slope(prof, t, rows):
    y, x = prof[rows], t[rows]
    m = np.isfinite(y)
    if m.sum() < 10: return np.nan
    y, x = y[m], x[m]
    idx = np.linspace(0, len(y)-1, min(len(y), 80)).astype(int)
    yy, xx = y[idx], x[idx]
    dx = xx[None,:]-xx[:,None]; dy = yy[None,:]-yy[:,None]
    return np.median(dy[dx>0]/dx[dx>0])

def analyze(site, cplx, coh, t, surf, inv_fn, show=True):
    d = loadmat(inv_fn)
    dtau_blk, ps = d['dtau_blk'], float(d['phase_sign'].ravel()[0])
    nblk = dtau_blk.shape[1]
    tt = t.ravel(); dt = np.median(np.diff(tt))
    edges = np.linspace(0, cplx.shape[1], nblk+1).astype(int)
    t0 = np.nanmedian(surf) + 0.5e-6   # skip the surface-reference zone
    nwin = int(2.0e-6/dt)
    rows0 = np.where(tt > t0)[0]
    out = []
    for b in range(nblk):
        cols = np.arange(edges[b], edges[b+1])
        for r0 in range(rows0[0], len(tt)-nwin, nwin):
            rows = np.arange(r0, r0+nwin)
            mc = np.nanmean(coh[np.ix_(rows, cols)])
            if mc < 0.15: continue
            s_cur = slope(dtau_blk[:, b], tt, rows)
            s_pp = ps * pulse_pair(cplx, coh, rows, cols, dt)
            s_ff, snr = fft_peak(cplx, coh, rows, cols, dt)
            out.append((b, tt[r0]*1e6, s_cur, s_pp, ps*s_ff, snr, mc))
    out = np.array(out)
    ok = np.isfinite(out[:,2]) & np.isfinite(out[:,3])
    hi = ok & (out[:,6] > 0.3)
    for name, m in [('all', ok), ('coh>0.3', hi)]:
        if m.sum() > 3:
            r = np.corrcoef(out[m,2], out[m,3])[0,1]
            md = np.median(np.abs(out[m,2]-out[m,3]))*1e6
            print(f'{site} [{name}] n={m.sum()}: corr(unwrap,pp)={r:.3f}, '
                  f'med|diff|={md:.0f} ps/us, '
                  f'med|rate| unwrap={np.median(np.abs(out[m,2]))*1e6:.0f} '
                  f'pp={np.median(np.abs(out[m,3]))*1e6:.0f}')
    if show:
        print('  blk  t_us | unwrap     pp    fft  (ps/us) |  snr   coh')
        for b, tus, sc, sp, sf, snr, mc in out:
            print(f'  {int(b):3d} {tus:5.1f} | {sc*1e6:6.0f} {sp*1e6:6.0f} '
                  f'{sf*1e6:6.0f} | {snr:6.0f} {mc:5.2f}')
    return out

base = os.path.expanduser('~/data/opr/accum/2024_Antarctica_Ground2')
p = loadmat(f'{base}/CSARP_polarimetric_unwrap/20250108_02/Data_20250108_02_009.mat',
            variable_names=['interferogram_mlook','interferogram_coherence','Time','Surface'])
ra = analyze('RidgeA', p['interferogram_mlook'], np.abs(p['interferogram_coherence']),
             p['Time'], p['Surface'],
             os.path.expanduser('~/data/opr/fabric_batch/joint_jp/2024_Antarctica_Ground2/20250108_02/Data_20250108_02_009.mat'))

m = loadmat(os.path.expanduser('~/data/opr/margin/margin_20240108_01_001.mat'))
tw = analyze('Thwaites', np.exp(1j*m['phase_wrapped']), m['coherence'],
             m['Time'].T, m['Surface'],
             os.path.expanduser('~/data/opr/fabric_batch/joint/2023_Antarctica_Ground/20240108_01/Data_20240108_01_001.mat'),
             show=False)

# Independent sign check at the margin: pulse-pair-integrated dtau vs coreg
dt = np.median(np.diff(m['Time'].ravel()))
coreg = m['row_offset'] * dt
coh = m['coherence']; cplx = np.exp(1j*m['phase_wrapped'])
tt = m['Time'].ravel(); surf = np.nanmedian(m['Surface'])
rows = np.where((tt > surf+1e-6) & (tt < surf+18e-6))[0]
pp_col = np.cumsum(np.where(np.isfinite(np.angle(cplx[rows])), 
         np.angle(cplx[rows][1:]*np.conj(cplx[rows][:-1])), 0), axis=0)/(2*np.pi*FC)
cg = coreg[rows][1:] - coreg[rows][0]
w = (coh[rows][1:]**2) * np.isfinite(cg)
num = np.nansum(w*pp_col*np.nan_to_num(cg)); 
print(f'\nMargin sign regression (pp-integrated vs coreg): sign = {np.sign(num):+.0f} '
      f'(production used +1)')

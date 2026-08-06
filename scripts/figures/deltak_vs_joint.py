"""Delta-k vs joint-inversion fabric comparison.

Ridge A (2024_Antarctica_Ground2, Paden CSARP_polarimetric) is the
clean-data validation: the SNAPHU joint chain is healthy there, so the
delta-k chain (no unwrapping, no fringe blending) should reproduce
fabric_joint_jp within noise. The Thwaites margin (2023_Antarctica_
Ground) is where the chains are expected to DISAGREE: the production
blend corrections moved dtau away from both SNAPHU and delta-k, and the
auto-detected phase sign there was flipped.

Only intervals genuinely measured in both chains are compared: pegged
intervals (|dlam| at the 2/3 layer-stripping bound) are dropped, and so
is any interval a dlam_interpolated node contaminates in EITHER chain
(fabric_qc.interpolated_intervals; the chains' flags legitimately differ
because the delta-k seam band is wider, so they are OR-ed per interval).

Outputs (figs/): deltak_vs_joint_ridge_a.png plus printed stats for
both seasons.

Usage: python deltak_vs_joint.py [fabric_batch_dir] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fabric_qc import interpolated_intervals

BATCH = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    '~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', 'figs')

RIDGE_A = {
    'deltak': os.path.join(BATCH, 'deltak_jp', '*', 'Data_*.mat'),
    'joint': os.path.join(BATCH, 'joint_jp', '2024_Antarctica_Ground2',
                          '*', 'Data_*.mat'),
}
THWAITES = {
    'deltak': os.path.join(BATCH, 'deltak', '*', 'Data_*.mat'),
    'joint': os.path.join(BATCH, 'joint', '2023_Antarctica_Ground',
                          '*', 'Data_*.mat'),
}
BOUND = 2.0 / 3.0  # layer-stripping dlam bound; pegged blocks excluded


def frame_key(fn):
    return os.path.basename(fn).replace('Data_', '').replace('.mat', '')


def load_pairs(spec):
    """Match frames across chains; return per-block paired arrays."""
    dk_fns = {frame_key(f): f for f in glob.glob(spec['deltak'])}
    jt_fns = {frame_key(f): f for f in glob.glob(spec['joint'])}
    keys = sorted(set(dk_fns) & set(jt_fns))
    dk, jt, depth, filled = [], [], [], []
    rms_dk, rms_jt = [], []
    for k in keys:
        a = loadmat(dk_fns[k], squeeze_me=True)
        b = loadmat(jt_fns[k], squeeze_me=True)
        if a['dlam'].shape != b['dlam'].shape:
            continue
        dk.append(np.asarray(a['dlam'], float).ravel())
        jt.append(np.asarray(b['dlam'], float).ravel())
        depth.append(0.5 * (np.asarray(a['dlam_top_depth'], float) +
                            np.asarray(a['dlam_bot_depth'], float)).ravel())
        # An interval only counts if it was measured in BOTH chains: OR
        # each chain's fabricated-interval mask (None on outputs written
        # before dlam_interpolated existed, meaning nothing to drop).
        bad = np.zeros(dk[-1].size, bool)
        for itp in (interpolated_intervals(a), interpolated_intervals(b)):
            if itp is not None:
                bad |= itp.ravel()
        filled.append(bad)
        rms_dk.append(np.atleast_1d(a['dtau_rms']).astype(float))
        rms_jt.append(np.atleast_1d(b['dtau_rms']).astype(float))
    dk = np.concatenate(dk) if dk else np.array([])
    jt = np.concatenate(jt) if jt else np.array([])
    depth = np.concatenate(depth) if depth else np.array([])
    filled = np.concatenate(filled) if filled else np.array([], bool)
    ok = (np.isfinite(dk) & np.isfinite(jt) & ~filled &
          (np.abs(dk) < 0.98 * BOUND) & (np.abs(jt) < 0.98 * BOUND))
    return {'keys': keys, 'dk': dk[ok], 'jt': jt[ok], 'depth': depth[ok],
            'n_all': dk.size,
            'rms_dk': np.concatenate(rms_dk) if rms_dk else np.array([]),
            'rms_jt': np.concatenate(rms_jt) if rms_jt else np.array([])}


def stats(tag, p):
    if p['dk'].size == 0:
        print('%s: no paired frames' % tag)
        return
    diff = p['dk'] - p['jt']
    r = np.corrcoef(p['dk'], p['jt'])[0, 1]
    sgn = np.mean(np.sign(p['dk']) == np.sign(p['jt']))
    print('%s: %d frames, %d/%d clean paired intervals' %
          (tag, len(p['keys']), p['dk'].size, p['n_all']))
    print('  corr %.3f | median|diff| %.4f | bias %+.4f | sign agree %.0f%%'
          % (r, np.median(np.abs(diff)), np.mean(diff), 100 * sgn))
    print('  dtau_rms median: deltak %.3f ns, joint %.3f ns'
          % (np.nanmedian(p['rms_dk']), np.nanmedian(p['rms_jt'])))


ra = load_pairs(RIDGE_A)
stats('Ridge A (jp)', ra)
th = load_pairs(THWAITES)
stats('Thwaites', th)

if ra['dk'].size == 0:
    sys.exit('Ridge A: no clean paired intervals under %s; '
             'nothing to plot' % BATCH)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
ax = axes[0]
sc = ax.scatter(ra['jt'], ra['dk'], c=ra['depth'], s=8, cmap='viridis',
                alpha=0.75, linewidths=0)
lim = 1.05 * max(np.abs(np.r_[ra['jt'], ra['dk']]).max(), 0.1)
ax.plot([-lim, lim], [-lim, lim], color='0.4', lw=0.8, zorder=0)
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xlabel('dlam, joint inversion (SNAPHU chain)')
ax.set_ylabel('dlam, delta-k chain')
ax.set_title('Ridge A: per-interval fabric contrast', fontsize=11)
ax.grid(alpha=0.25, lw=0.5)
fig.colorbar(sc, ax=ax, pad=0.02, label='Interval mid-depth (m)')

ax = axes[1]
diff = ra['dk'] - ra['jt']
ax.hist(diff, bins=41, color='#0072B2', alpha=0.85)
ax.axvline(0, color='0.4', lw=0.8)
ax.set_xlabel('dlam difference (delta-k minus joint)')
ax.set_ylabel('Intervals')
ax.set_title('median |diff| = %.4f, bias = %+.4f'
             % (np.median(np.abs(diff)), np.mean(diff)), fontsize=11)
ax.grid(alpha=0.25, lw=0.5)

fig.suptitle('Delta-k vs joint inversion, Ridge A grid '
             '(%d frames, clean intervals)' % len(ra['keys']),
             fontsize=12)
fig.tight_layout()
os.makedirs(OUT, exist_ok=True)
out = os.path.join(OUT, 'deltak_vs_joint_ridge_a.png')
fig.savefig(out, dpi=200, bbox_inches='tight')
print('wrote', out)

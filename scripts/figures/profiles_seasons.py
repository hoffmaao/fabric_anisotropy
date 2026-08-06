"""Depth-profile figures of the horizontal fabric contrast per season.

Reads CSARP_fabric_joint outputs mirrored under ~/data/opr/fabric_batch
(see opr_fabric/server/run_fabric_scratch.m) and renders, per dataset,
the season-median profile with IQR and per-segment median profiles.

Usage: python profiles_seasons.py [fabric_batch_root] [out_dir]
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

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')

DEPTH = np.arange(0, 1900, 5.0)
DATASETS = [
    ('2022 Antarctica\n(thin ice)', f'{ROOT}/joint/2022_Antarctica_Ground/*/Data_*.mat'),
    ('2023 Antarctica', f'{ROOT}/joint/2023_Antarctica_Ground/*/Data_*.mat'),
    ('2024-25 (Lilien proc.)', f'{ROOT}/joint/2024_Antarctica_Ground2/*/Data_*.mat'),
    ('2024-25 (Paden proc.)', f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat'),
    ('2025-26', f'{ROOT}/joint/2025_Antarctica_Ground2/*/Data_*.mat'),
]


def block_cols(pattern):
    """Per-block QC'd profiles on the common depth grid, grouped by segment."""
    segs = {}
    for fn in sorted(glob.glob(pattern)):
        seg = os.path.basename(os.path.dirname(fn))
        d = loadmat(fn)
        dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
        qual = d['dlam_quality']
        clip = d.get('dlam_clipped')
        itp = interpolated_intervals(d)
        for b in range(dlam.shape[1]):
            if not np.any(np.isfinite(dlam[:, b])):
                continue
            col = np.full(DEPTH.size, np.nan)
            for k in range(dlam.shape[0]):
                pegged = (clip is not None and clip.size and clip[k, b] == 1) \
                    or abs(dlam[k, b]) > 0.6
                filled = itp is not None and itp[k, b]
                ok = np.isfinite(dlam[k, b]) and not pegged and not filled \
                    and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
                if ok:
                    m = (DEPTH >= top[k, b]) & (DEPTH < bot[k, b])
                    col[m] = dlam[k, b]
            segs.setdefault(seg, []).append(col)
    return segs


def main():
    fig, axes = plt.subplots(2, len(DATASETS), figsize=(16, 10), sharey=True)
    cmap = plt.get_cmap('viridis')

    for ci, (name, pattern) in enumerate(DATASETS):
        segs = block_cols(pattern)
        if not segs:
            axes[0, ci].set_title(f'{name}\n(no data)')
            continue
        all_cols = np.array([c for cols in segs.values() for c in cols])

        ax = axes[0, ci]
        with np.errstate(invalid='ignore'):
            med = np.nanmedian(all_cols, axis=0)
            q1 = np.nanpercentile(all_cols, 25, axis=0)
            q3 = np.nanpercentile(all_cols, 75, axis=0)
            nsamp = np.sum(np.isfinite(all_cols), axis=0)
        for arr in (med, q1, q3):
            arr[nsamp < 5] = np.nan
        ax.fill_betweenx(DEPTH, q1, q3, color='steelblue', alpha=0.3, label='IQR')
        ax.plot(med, DEPTH, 'b-', lw=2, label='median')
        ax.axvline(0, color='gray', lw=0.5)
        ax.set_title(f'{name}\n{all_cols.shape[0]} blocks, {len(segs)} segments', fontsize=10)
        ax.set_xlim(-0.11, 0.11)
        ax.grid(alpha=0.3)
        if ci == 0:
            ax.set_ylabel('Depth (m)')
            ax.legend(fontsize=8, loc='lower right')

        ax = axes[1, ci]
        seg_names = list(segs.keys())
        for si, seg in enumerate(seg_names):
            cols = np.array(segs[seg])
            with np.errstate(invalid='ignore'):
                smed = np.nanmedian(cols, axis=0)
                n = np.sum(np.isfinite(cols), axis=0)
            smed[n < 2] = np.nan
            ax.plot(smed, DEPTH, color=cmap(si / max(1, len(seg_names) - 1)),
                    lw=1, alpha=0.85)
        ax.axvline(0, color='gray', lw=0.5)
        ax.set_xlim(-0.11, 0.11)
        ax.grid(alpha=0.3)
        ax.set_xlabel(r'$\Delta\lambda = \lambda_{cross} - \lambda_{along}$')
        if ci == 0:
            ax.set_ylabel('Depth (m)')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(1, len(seg_names)))
        cb = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
        cb.set_label('segment order', fontsize=7)
        cb.ax.tick_params(labelsize=7)

    axes[0, 0].invert_yaxis()
    fig.suptitle(
        'Horizontal fabric contrast profiles, EAGER polarimetric traverses (joint inversion)\n'
        'top: season median with interquartile range; bottom: per-segment medians along traverse',
        fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'profiles_all_seasons.png')
    fig.savefig(out, dpi=140)
    print('saved', out)


if __name__ == '__main__':
    main()

"""Quad-pol fabric section along one profile, beside the co-polarized one.

run_quadpol_pipeline.m inverts each 125-trace block of a coregistered
quad-pol frame with ptt.ershadiFabric and ptt.quadpolFabricLS, giving dlam
and theta as functions of BOTH depth and along-track distance. 125 traces
is what run_sections.m uses for the co-polarized Ridge A section, so the
two are sampled the same way along track and can be read against each
other cell for cell. When the LS fields are present the figure adds a
third row - the LS coherence-field fit above the published direct chain,
on a shared color scale; older .mat files without them still draw the
two-row layout.

WHAT THE TWO SECTIONS ARE NOT. They do not measure the same quantity, and
drawing them side by side without saying so would invite reading the
difference as disagreement:

  co-polarized   dlam projected onto THIS LINE's axes, lam_perp - lam_par,
                 which is -P cos 2(alpha - theta) and so depends on the
                 heading the vehicle happened to drive
  quad-pol       lam_max - lam_min at the principal axes, recovered by
                 synthesizing the azimuth sweep, which in principle does
                 not

That is the whole claim being tested, so the panels are labelled with what
each one is rather than both with "dlam".

THE ORIENTATION PANEL CARRIES A WARNING. On Ridge A the recovered theta
tracks the ANTENNA frame rather than the ice - measured over 403 blocks at
2.0 deg circular spread in the antenna frame against 41.7 geographic, and
confirmed on two individually coregistered frames whose track azimuths
differ by 98 deg (theta_ant 90.0 and 86.7). So the theta panel is drawn as
a diagnostic, not as a fabric orientation, and says so.

Usage: python quadpol_section.py <out_dir> [tag]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import h5py                              # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.colors import TwoSlopeNorm   # noqa: E402
from scipy.io import loadmat             # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import DATA, INK, MUTED, cumdist_km, field  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
TAG = sys.argv[2] if len(sys.argv) > 2 else '20250108_02_009'
# run_quadpol_pipeline.m writes this one; run_quadpol_frame.m writes
# quadpol_<tag>.mat, which is a different layout read by
# quadpol_diagnostic.py. The prefixes are distinct so both can be staged.
FN = os.path.join(DATA, 'quadpol_section_%s.mat' % TAG)
CO_FN = os.path.join(DATA, 'fabric_sections.mat')

Z_SHOW = (0.0, 1500.0)
CMAG_MIN = 0.25     # |C_HHVV| below which a cell is greyed rather than drawn
GREY = np.array([0.74, 0.74, 0.74])


def shade(v, w, norm, cmap):
    """Diverging colour blended toward grey by confidence, as the
    co-polarized section does, so the two use the same visual grammar."""
    rgb = cmap(norm(v))[..., :3]
    ww = np.clip(w, 0, 1)[..., None]
    out = rgb * ww + GREY[None, None, :] * (1 - ww)
    out[~np.isfinite(v)] = 1.0
    return out


def main():
    if not os.path.exists(FN):
        raise SystemExit(
            'missing %s; run opr_fabric/server/run_quadpol_pipeline.m for '
            'this frame and mirror the result there' % FN)
    # h5py, not loadmat: the pipeline saves -v7.3, which scipy.io cannot
    # read at all. HDF5 stores MATLAB arrays transposed, so the [Nz x nb]
    # sections come back [nb x Nz].
    with h5py.File(FN) as f:
        res = f['res']

        def a(name):
            return np.array(res[name])

        z = a('z').ravel().astype(float)
        dl = np.atleast_2d(a('sec_dlam').T).astype(float)
        th = np.atleast_2d(a('sec_theta').T).astype(float)
        cm = np.atleast_2d(a('sec_cmag').T).astype(float)
        lat = a('sec_lat').ravel().astype(float)
        lon = a('sec_lon').ravel().astype(float)
        track_az = float(a('track_az').ravel()[0])
        nbt = int(a('nblk_tr').ravel()[0])
        # the LS-fit section exists once the pipeline has run with
        # ptt.quadpolFabricLS wired in; older .mat files draw the old layout
        has_ls = 'sec_dlam_ls' in res
        dl_ls = np.atleast_2d(a('sec_dlam_ls').T).astype(float) if has_ls \
            else None
    dist = cumdist_km(lat, lon)

    ok = np.isfinite(dl)
    print('%s: %d depths x %d blocks, %d finite cells' % (TAG, *dl.shape,
                                                          ok.sum()))
    print('  dlam %.3f..%.3f (median %.3f)'
          % (np.nanmin(dl), np.nanmax(dl), np.nanmedian(dl)))
    print('  |C_HHVV| median %.3f' % np.nanmedian(cm))
    if has_ls:
        print('  LS dlam %.3f..%.3f (median %.3f)'
              % (np.nanmin(dl_ls), np.nanmax(dl_ls), np.nanmedian(dl_ls)))

    sel = (z >= Z_SHOW[0]) & (z <= Z_SHOW[1])
    z = z[sel]
    dl, th, cm = dl[sel], th[sel], cm[sel]
    if has_ls:
        dl_ls = dl_ls[sel]

    # confidence from the coherence, on the same idea the co-pol section
    # uses for node quality: below CMAG_MIN a cell fades out rather than
    # being drawn in a colour it has not earned
    w = np.clip((cm - 0.10) / (0.45 - 0.10), 0, 1) ** 0.7
    w[cm < CMAG_MIN] *= 0.35

    trust = np.isfinite(dl) & (w > 0.4)
    vals = np.abs(dl[trust])
    if has_ls:
        vals = np.concatenate(
            [vals, np.abs(dl_ls[np.isfinite(dl_ls) & (w > 0.4)])])
    lim = float(np.nanpercentile(vals, 96)) if vals.size else 0.1
    norm = TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim)
    cmap = plt.get_cmap('RdBu_r')

    xe = np.concatenate([[dist[0] - (dist[1] - dist[0]) / 2],
                         (dist[:-1] + dist[1:]) / 2,
                         [dist[-1] + (dist[-1] - dist[-2]) / 2]])

    nrow = 3 if has_ls else 2
    fig = plt.figure(figsize=(13.0, 8.6 if has_ls else 6.4),
                     layout='constrained')
    gs = fig.add_gridspec(nrow, 2, width_ratios=[2.6, 1.0])
    ax_ls = fig.add_subplot(gs[0, 0]) if has_ls else None
    axs = fig.add_subplot(gs[1 if has_ls else 0, 0],
                          sharex=ax_ls if has_ls else None)
    axt = fig.add_subplot(gs[nrow - 1, 0], sharex=axs)
    axp = fig.add_subplot(gs[:, 1], sharey=axs)

    if has_ls:
        ax_ls.imshow(shade(dl_ls, w, norm, cmap), origin='upper',
                     aspect='auto', extent=[xe[0], xe[-1], z[-1], z[0]],
                     interpolation='nearest')
        ax_ls.set_ylabel('depth (m)', color=INK)
        ax_ls.set_title('LS coherence-field fit   '
                        r'$\Delta\lambda$ at the coherence-derived axes '
                        '(theta0 fixed per frame)', fontsize=11, color=INK)
        plt.setp(ax_ls.get_xticklabels(), visible=False)
        cax = ax_ls.inset_axes([0.02, -0.10, 0.30, 0.045])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                     orientation='horizontal', label=r'$\Delta\lambda$')

    axs.imshow(shade(dl, w, norm, cmap), origin='upper', aspect='auto',
               extent=[xe[0], xe[-1], z[-1], z[0]], interpolation='nearest')
    axs.set_ylabel('depth (m)', color=INK)
    if has_ls:
        axs.set_title('published direct chain (Ershadi) - evaluated at the '
                      'cross-pol minimum, antenna-locked; for comparison',
                      fontsize=10, color=MUTED)
    else:
        axs.set_title('quad-pol inversion, %s   '
                      r'$\Delta\lambda=\lambda_{max}-\lambda_{min}$ '
                      '(principal axes)' % TAG, fontsize=11, color=INK)
        cax = axs.inset_axes([0.02, -0.10, 0.30, 0.045])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                     orientation='horizontal', label=r'$\Delta\lambda$')
    plt.setp(axs.get_xticklabels(), visible=False)

    # orientation, drawn as a diagnostic
    tn = TwoSlopeNorm(vcenter=90, vmin=0, vmax=180)
    axt.imshow(shade(th, w, tn, plt.get_cmap('twilight')), origin='upper',
               aspect='auto', extent=[xe[0], xe[-1], z[-1], z[0]],
               interpolation='nearest')
    axt.set_xlabel('distance along profile (km)', color=INK)
    axt.set_ylabel('depth (m)', color=INK)
    axt.set_title(r'$\theta$ (deg E of N) - DIAGNOSTIC: on Ridge A this '
                  'tracks the antenna frame, not the ice', fontsize=9.5,
                  color=MUTED)

    # depth profile, with the co-polarized section's median for comparison
    def block_median(sec):
        m = np.full(z.size, np.nan)
        lo_b = np.full(z.size, np.nan)
        hi_b = np.full(z.size, np.nan)
        for i in range(z.size):
            v, ww = sec[i], w[i]
            good = np.isfinite(v) & (ww > 0.15)
            if good.sum() >= 3:
                m[i] = np.median(v[good])
                lo_b[i], hi_b[i] = np.percentile(v[good], [16, 84])
        return m, lo_b, hi_b

    med, lo_p, hi_p = block_median(dl)
    if has_ls:
        med_ls, lo_ls, hi_ls = block_median(dl_ls)
        axp.fill_betweenx(z, lo_ls, hi_ls, color='#2a78d6', alpha=0.20, lw=0,
                          label='LS block-to-block 16-84%')
        axp.plot(med_ls, z, '-', color='#2a78d6', lw=2.0, label='LS fit')
        axp.plot(med, z, '-', color='#9a9a9a', lw=1.4,
                 label='published chain (locked)')
    else:
        axp.fill_betweenx(z, lo_p, hi_p, color='#2a78d6', alpha=0.20, lw=0,
                          label='block-to-block 16-84%')
        axp.plot(med, z, '-', color='#2a78d6', lw=2.0,
                 label='quad-pol median')

    if os.path.exists(CO_FN):
        S = loadmat(CO_FN, squeeze_me=True)['S']
        if 'ridge_a' in (S.dtype.names or ()):
            r = S['ridge_a'].item() if S['ridge_a'].dtype == object \
                else S['ridge_a']
            nint = int(field(r, 'nint')[0])
            cdl = field(r, 'dlam').reshape(nint, -1)
            ctop = field(r, 'top').reshape(cdl.shape)
            cbot = field(r, 'bot').reshape(cdl.shape)
            zc = np.nanmedian((ctop + cbot) / 2, axis=1)
            axp.plot(np.nanmedian(cdl, axis=1), zc, '--', color='#eb6834',
                     lw=1.8, label='co-pol, projected on this line')
    axp.axvline(0, color=MUTED, lw=0.8)
    axp.set_xlabel(r'$\Delta\lambda$', color=INK)
    axp.set_title('depth profile', fontsize=11, color=INK)
    axp.grid(alpha=0.25, lw=0.6)
    axp.set_axisbelow(True)
    axp.legend(loc='lower right', fontsize=7.5, frameon=False)
    plt.setp(axp.get_yticklabels(), visible=False)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)

    for ax in ((ax_ls, axs, axt, axp) if has_ls else (axs, axt, axp)):
        ax.set_ylim(Z_SHOW[1], Z_SHOW[0])

    fig.suptitle('Ridge A %s: quad-pol section (track %.0f$^\\circ$, '
                 '%d blocks of %d traces)'
                 % (TAG, track_az, dl.shape[1], nbt),
                 fontsize=12, color=INK)
    out = os.path.join(OUT, 'scar_quadpol_section_%s.png' % TAG)
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', out)


if __name__ == '__main__':
    main()

"""2D fabric sections, shaded by node coherence, one figure per site.

Same idea as the interferogram panels: colour carries the quantity, and
how strongly it is painted carries how well it is constrained. There the
mapping was HSV value = coherence, which fades to black. That cannot be
reused here, because the quantity is on a DIVERGING ramp whose extremes
are already dark - dimming a strong positive dlam would make it look
like a weak one. Instead each cell is blended toward a neutral grey by
its node coherence, so low-confidence cells desaturate rather than
darken and cannot be mistaken for either the pale midpoint or a dark
extreme.

Nodes the inversion filled by interpolating across a gap (waveform-combine
seam, incoherent run) are down-weighted too: dlam_quality is the unmasked
coherence and does not see them.

The two-dimensional key at the lower right shows both axes of that
mapping, because a pair of one-dimensional colourbars would imply the
fade is independent of the value, and it is not.

Usage: python fabric_sections.py <out_dir> [site …]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt             # noqa: E402
from matplotlib.colors import TwoSlopeNorm   # noqa: E402
from scipy.io import loadmat                 # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scar_style import (DATA, INK, MUTED, cumdist_km,  # noqa: E402
                        field, track_azimuth)
from scar_style import FIGS  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else FIGS
GREY = np.array([0.74, 0.74, 0.74])   # what an unconstrained cell fades to
# Node coherence spans ~0.15-0.90 (median ~0.47 at Ridge A, Thwaites and
# NEGIS, ~0.34 on the noisier EastGRIP borehole line), so
# these limits put the fade across the range that actually occurs. The
# gamma keeps mid-coherence cells legible instead of washing out the
# middle of every column; only genuinely poor nodes go fully grey.
COH_LO, COH_HI, COH_GAMMA = 0.15, 0.70, 0.7

TITLES = {
    'ridge_a': ('Ridge A divide', '20250108_02_009', 150.0),
    'thwaites': ('Thwaites eastern shear margin', '20240108_01_001', 150.0),
    'negis': ('NEGIS onset (EGRIP)', '20240626_03_001', 200.0),
    'egrip': ('EastGRIP borehole line', '20240619_01_001', 200.0),
}
SRC_LABEL = {'phase': 'SNAPHU + coregistration', 'deltak': 'delta-k'}
# The eigenvalue pair is defined by the profile's own axes, so name the
# geographic directions rather than saying "cross" and "along".
LBL = (r'$\Delta\lambda = \lambda_{\perp} - \lambda_{\parallel}$'
       '\n' r'($\perp$ = %.0f$^\circ$,  $\parallel$ = %.0f$^\circ$ E of N)')


def shade(values, weight, norm, cmap):
    """Diverging colour for `values`, blended toward grey by `weight`."""
    rgb = cmap(norm(values))[..., :3]
    w = np.clip(weight, 0, 1)[..., None]
    out = rgb * w + GREY[None, None, :] * (1 - w)
    out[~np.isfinite(values)] = 1.0     # unsolved cells: page white
    return out


def confidence(q, interp):
    """Node coherence scaled to 0..1, halved where the node was filled in."""
    w = (q - COH_LO) / (COH_HI - COH_LO)
    w = np.clip(np.nan_to_num(w, nan=0.0), 0, 1) ** COH_GAMMA
    if interp is not None:
        w = w * (1 - 0.5 * np.nan_to_num(interp, nan=0.0))
    return w


def key_axes(ax, norm, cmap, vmin, vmax):
    """Two-dimensional key: dlam across, node coherence up."""
    nx, ny = 128, 48
    xs = np.linspace(vmin, vmax, nx)
    # Ramp the KEY through confidence() as well, not through the raw blend
    # weight: the axis is labelled in node coherence, and the COH_GAMMA
    # bend means a cell at coherence 0.425 is painted at weight 0.62 in the
    # section. A linear ramp here would have a reader match that saturation
    # to 0.50 on the key and read the coherence systematically too high.
    ys = confidence(np.linspace(COH_LO, COH_HI, ny), None)
    img = shade(np.tile(xs, (ny, 1)), np.tile(ys[:, None], (1, nx)),
                norm, cmap)
    ax.imshow(img, origin='lower', aspect='auto',
              extent=[vmin, vmax, COH_LO, COH_HI], interpolation='bilinear')
    ax.set_xlabel(r'$\Delta\lambda$', fontsize=7.5, color=INK, labelpad=1)
    ax.set_ylabel('node coh.', fontsize=7.5, color=INK, labelpad=1)
    ax.tick_params(labelsize=6.5, length=2, pad=1)
    # Label zero only when it is far enough from an end to not collide:
    # at Ridge A vmin is -0.011 and "-0.01" ran straight into "0.00".
    ticks = [vmin, vmax]
    if vmin < 0 < vmax and min(-vmin, vmax) > 0.22 * (vmax - vmin):
        ticks = [vmin, 0, vmax]
    ax.set_xticks(ticks)
    ax.set_xticklabels(['%.2f' % v for v in ticks])
    ax.set_yticks([COH_LO, COH_HI])
    for s in ax.spines.values():
        s.set_linewidth(0.6)
        s.set_color('0.4')


def render(name, rec, out_dir):
    title, frame, z_firn = TITLES[name]
    nint = int(field(rec, 'nint')[0])
    dlam = field(rec, 'dlam').reshape(nint, -1)
    top = field(rec, 'top').reshape(dlam.shape)
    bot = field(rec, 'bot').reshape(dlam.shape)
    qual = field(rec, 'quality').reshape(dlam.shape)
    try:
        interp = field(rec, 'interpolated').reshape(dlam.shape)
    except (KeyError, ValueError):
        interp = None
    lat, lon = field(rec, 'lat'), field(rec, 'lon')
    dist = cumdist_km(lat, lon)
    az = track_azimuth(lat, lon)
    az_cross = (az + 90) % 180
    src = str(np.atleast_1d(rec['src'])[0])
    gate = field(rec, 'gate')[0]
    rms = field(rec, 'rms')[0]
    adj = field(rec, 'adj')[0]
    nblk = dlam.shape[1]

    w = confidence(qual, interp)
    # Range from the cells the figure actually asks you to believe. Taking
    # percentiles over everything lets the greyed-out noise at the edges
    # set the limits and washes the trusted interior to nothing - at
    # Thwaites that stretched the scale to +-0.4 for a profile whose real
    # excursions are +-0.05.
    trust = dlam[np.isfinite(dlam) & (w > 0.4)]
    if trust.size < 20:
        trust = dlam[np.isfinite(dlam)]
    vmax = float(np.nanpercentile(trust, 98))
    vmin = float(min(np.nanpercentile(trust, 2), -0.15 * abs(vmax)))
    if vmax <= 0:
        vmax = abs(vmin)
    norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap('RdBu_r')

    zc = np.nanmedian(np.vstack([top[0:1, :], bot]), axis=1)
    xe = np.concatenate([[dist[0] - (dist[1] - dist[0]) / 2],
                         (dist[:-1] + dist[1:]) / 2,
                         [dist[-1] + (dist[-1] - dist[-2]) / 2]])

    fig = plt.figure(figsize=(12.2, 6.1), layout='constrained')
    gs = fig.add_gridspec(2, 2, width_ratios=[2.5, 1.0],
                          height_ratios=[1.0, 0.20])
    axs = fig.add_subplot(gs[0, 0])
    axp = fig.add_subplot(gs[0, 1], sharey=axs)
    # Key and caption live in their own row: as an inset the key covered
    # the section's lower right, which is data at every site.
    sub = gs[1, 0].subgridspec(1, 2, width_ratios=[0.26, 0.74])
    axk = fig.add_subplot(sub[0, 0])
    axt = fig.add_subplot(sub[0, 1])
    axt.axis('off')

    # pcolormesh cannot take per-cell RGB, so the shaded array is drawn as
    # an image on the non-uniform grid via NonUniformImage-style extents;
    # interval thicknesses here are equal in twtt, so a linear depth axis
    # is a fair approximation and imshow keeps the RGB.
    axs.imshow(shade(dlam, w, norm, cmap), origin='upper', aspect='auto',
               extent=[xe[0], xe[-1], zc[-1], zc[0]],
               interpolation='nearest')
    axs.set_xlabel('distance along profile (km)', color=INK)
    axs.set_ylabel('depth (m)', color=INK)
    axs.set_title('%s   -   %s' % (title, frame), fontsize=11, color=INK)
    axs.axhline(z_firn, color='0.25', lw=0.9, ls=':')
    axs.annotate('firn / few fringes', xy=(0.012, z_firn + 0.012 * zc[-1]),
                 xycoords=('axes fraction', 'data'), fontsize=7.5,
                 color='0.15', va='top')

    key_axes(axk, norm, cmap, vmin, vmax)

    # Depth profile: coherence-weighted median, so the profile and the
    # section are telling the same story about what is trusted
    med = np.full(nint, np.nan)
    lo = np.full(nint, np.nan)
    hi = np.full(nint, np.nan)
    for i in range(nint):
        v, ww = dlam[i], w[i]
        ok = np.isfinite(v) & (ww > 0.15)
        if ok.sum() >= 3:
            med[i] = np.median(v[ok])
            lo[i], hi[i] = np.percentile(v[ok], [16, 84])
    zmid = np.nanmedian((top + bot) / 2, axis=1)
    axp.fill_betweenx(zmid, lo, hi, color='#2a78d6', alpha=0.22, lw=0,
                      label='block-to-block 16-84%')
    axp.plot(med, zmid, '-', color='#2a78d6', lw=2.0, label='median profile')
    axp.axvline(0, color=MUTED, lw=0.8)
    axp.axhline(z_firn, color='0.25', lw=0.9, ls=':')
    axp.grid(alpha=0.25, lw=0.6)
    axp.set_axisbelow(True)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)
    axp.set_xlabel(LBL % (az_cross, az), color=INK)
    axp.set_title('depth profile', fontsize=11, color=INK)
    axp.legend(loc='lower right', fontsize=8, frameon=False)
    plt.setp(axp.get_yticklabels(), visible=False)

    axt.text(0.0, 0.85,
             '%s dtau   |   %d intervals (~%.0f m) x %d blocks (~%.0f m)   |'
             '   coherence gate %.2f, misfit %.3f ns, adjacent-block '
             'agreement %.2f\ncells fade to grey where the interval is '
             'weakly constrained; nodes interpolated across gaps are '
             'halved again'
             % (SRC_LABEL.get(src, src), nint,
                abs(np.nanmedian(np.diff(zc))), nblk,
                1000 * np.median(np.diff(dist)), gate, rms, adj),
             transform=axt.transAxes, fontsize=7.5, color=MUTED,
             ha='left', va='top')

    out = os.path.join(out_dir, 'fabric_section_%s.png' % name)
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('wrote %s  (%d x %d, dlam %.3f..%.3f)'
          % (out, nint, nblk, vmin, vmax))


def main():
    os.makedirs(OUT, exist_ok=True)
    S = loadmat(os.path.join(DATA, 'fabric_sections.mat'),
                squeeze_me=True)['S']
    names = sys.argv[2:] or list(S.dtype.names)
    for name in names:
        if name not in S.dtype.names:
            print('no section for', name)
            continue
        render(name, S[name].item() if S[name].dtype == object else S[name],
               OUT)


if __name__ == '__main__':
    main()

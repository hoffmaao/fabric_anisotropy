"""Ridge A azimuthal analysis: horizontal fabric orientation and strength.

The Ridge A raster grid samples the same fabric at many drive headings.
For a line at heading phi over fabric with principal horizontal
eigenvalues lam1 > lam2 (lam1-axis at azimuth theta), the survey-frame
contrast measured by the common-offset inversion is

    dlam(phi) = lam_cross - lam_along = -P cos 2(phi - theta),
    P = lam1 - lam2 >= 0,

which is linear in (A, B) = (-P cos 2theta, -P sin 2theta). At each depth
we fit (A, B) by quality-weighted least squares over all grid blocks,
giving the principal contrast P(z), the lam1 azimuth theta(z) (deg east
of north, mod 180), 1-sigma uncertainties from the residual covariance,
and the azimuthal-model R^2 (a test of the fabric-uniform-across-grid
assumption).

Uses the Paden 2024-25 processing (more blocks); the Lilien processing is
overlaid as an independent check on P(z).

Usage: python ridge_a_azimuthal.py [fabric_batch_root] [out_dir]
"""
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/data/opr/fabric_batch')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', '..', 'figs')

DEPTH = np.arange(0, 1900, 5.0)
RIDGE_A = (-86.60, 69.0)
MAX_DIST_KM = 200.0   # keep only blocks near Ridge A
MIN_BLOCKS = 12       # minimum azimuth samples per depth for a fit
MIN_SPREAD = 30.0     # minimum circular heading spread (deg, mod 180) per depth
EXAMPLE_DEPTHS = [400, 900, 1300, 1600]


def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def bearing_deg(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    x = np.sin(dl) * np.cos(p2)
    y = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(x, y)) % 360.0


def load_ridge_a(pattern):
    """QC'd profiles on the common grid + per-block heading and weight."""
    profs, heads, wts, lats, lons = [], [], [], [], []
    for fn in sorted(glob.glob(pattern)):
        d = loadmat(fn)
        dlam, top, bot = d['dlam'], d['dlam_top_depth'], d['dlam_bot_depth']
        qual = d['dlam_quality']
        clip = d.get('dlam_clipped')
        itp = d.get('dlam_interpolated')
        blat = d['Latitude'][0, :]
        blon = d['Longitude'][0, :]
        nb = dlam.shape[1]
        for b in range(nb):
            lat, lon = blat[b], blon[b]
            if not (np.isfinite(lat) and np.isfinite(lon)):
                continue
            if haversine_km(lat, lon, *RIDGE_A) > MAX_DIST_KM:
                continue
            # Heading from the neighboring block along the drive
            j = b + 1 if b + 1 < nb else b - 1
            if j < 0 or not (np.isfinite(blat[j]) and np.isfinite(blon[j])):
                continue
            if j > b:
                h = bearing_deg(lat, lon, blat[j], blon[j])
            else:
                h = bearing_deg(blat[j], blon[j], lat, lon)
            col = np.full(DEPTH.size, np.nan)
            for k in range(dlam.shape[0]):
                pegged = (clip is not None and clip.size and clip[k, b] == 1) \
                    or abs(dlam[k, b]) > 0.6
                # dtau interpolated across a masked gap (e.g. a waveform-
                # combine seam) rather than measured at this node
                filled = itp is not None and itp.size and itp[k, b] == 1
                ok = np.isfinite(dlam[k, b]) and not pegged and not filled \
                    and np.isfinite(qual[k, b]) and qual[k, b] > 0.35
                if ok:
                    m = (DEPTH >= top[k, b]) & (DEPTH < bot[k, b])
                    col[m] = dlam[k, b]
            if not np.any(np.isfinite(col)):
                continue
            profs.append(col)
            heads.append(h % 180.0)
            wts.append(np.nanmean(qual[:, b])**2)
            lats.append(lat)
            lons.append(lon)
    return (np.array(profs), np.array(heads), np.array(wts),
            np.array(lats), np.array(lons))


def fit_azimuthal(profs, heads, wts):
    """Per-depth weighted LSQ of dlam = A cos2phi + B sin2phi."""
    nz = DEPTH.size
    P = np.full(nz, np.nan)
    theta = np.full(nz, np.nan)
    P_sig = np.full(nz, np.nan)
    theta_sig = np.full(nz, np.nan)
    r2 = np.full(nz, np.nan)
    nfit = np.zeros(nz, int)
    phi2 = np.radians(2.0 * heads)
    for iz in range(nz):
        y = profs[:, iz]
        m = np.isfinite(y)
        if m.sum() < MIN_BLOCKS:
            continue
        # Require real azimuthal diversity, not one line direction:
        # circular spread mod 180 from the largest gap on doubled angles
        ang = np.sort((2.0 * heads[m]) % 360.0)
        gap = np.diff(ang, append=ang[0] + 360.0).max()
        if (360.0 - gap) / 2.0 < MIN_SPREAD:
            continue
        X = np.column_stack([np.cos(phi2[m]), np.sin(phi2[m])])
        w = wts[m]
        W = w / w.sum()
        XtW = X.T * W
        cov_inv = XtW @ X
        try:
            cov = np.linalg.inv(cov_inv)
        except np.linalg.LinAlgError:
            continue
        ab = cov @ (XtW @ y[m])
        resid = y[m] - X @ ab
        dof = max(m.sum() - 2, 1)
        s2 = np.sum(W * resid**2) * m.sum() / dof
        ab_cov = cov * s2 / m.sum()
        A, B = ab
        P[iz] = np.hypot(A, B)
        theta[iz] = 0.5 * np.degrees(np.arctan2(-B, -A)) % 180.0
        # First-order error propagation
        if P[iz] > 0:
            g_p = np.array([A, B]) / P[iz]
            P_sig[iz] = np.sqrt(g_p @ ab_cov @ g_p)
            g_t = np.array([-B, A]) / (2 * P[iz]**2)
            theta_sig[iz] = np.degrees(np.sqrt(g_t @ ab_cov @ g_t))
        ss_tot = np.sum(W * (y[m] - np.sum(W * y[m]))**2)
        r2[iz] = 1.0 - np.sum(W * resid**2) / ss_tot if ss_tot > 0 else np.nan
        nfit[iz] = m.sum()
    return P, theta, P_sig, theta_sig, r2, nfit


def main():
    profs, heads, wts, lats, lons = load_ridge_a(
        f'{ROOT}/joint_jp/2024_Antarctica_Ground2/*/Data_*.mat')
    print(f'Ridge A blocks (Paden): {len(profs)}; heading spread '
          f'{np.percentile(heads, 5):.0f}-{np.percentile(heads, 95):.0f} deg mod 180')
    P, theta, P_sig, theta_sig, r2, nfit = fit_azimuthal(profs, heads, wts)

    profs_l, heads_l, wts_l, _, _ = load_ridge_a(
        f'{ROOT}/joint/2024_Antarctica_Ground2/*/Data_*.mat')
    P_l = theta_l = None
    if len(profs_l) >= MIN_BLOCKS:
        P_l, theta_l, _, _, _, _ = fit_azimuthal(profs_l, heads_l, wts_l)
        print(f'Ridge A blocks (Lilien check): {len(profs_l)}')

    ok = np.isfinite(P)
    print(f'fitted depths: {ok.sum()} of {DEPTH.size}; '
          f'median R^2 where fitted: {np.nanmedian(r2[ok]):.2f}')

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 4, height_ratios=[2.2, 1])

    # (a) Strength P(z)
    ax = fig.add_subplot(gs[0, 0])
    ax.fill_betweenx(DEPTH, P - P_sig, P + P_sig, color='tab:red', alpha=0.25)
    ax.plot(P, DEPTH, 'tab:red', lw=2, label='Paden proc.')
    if P_l is not None:
        ax.plot(P_l, DEPTH, 'k--', lw=1.2, label='Lilien proc.')
    ax.invert_yaxis()
    ax.set_xlabel(r'$P = \lambda_1 - \lambda_2$')
    ax.set_ylabel('Depth (m)')
    ax.set_xlim(0, 0.15)
    ax.grid(alpha=0.3)
    ax.set_title('Principal contrast (strength)')
    ax.legend(fontsize=8)

    # (b) Azimuth theta(z)
    ax = fig.add_subplot(gs[0, 1], sharey=fig.axes[0])
    ax.fill_betweenx(DEPTH, theta - theta_sig, theta + theta_sig,
                     color='tab:blue', alpha=0.25)
    ax.plot(theta, DEPTH, 'tab:blue', lw=2)
    if theta_l is not None:
        ax.plot(theta_l, DEPTH, 'k--', lw=1.2)
    ax.set_xlabel(r'$\theta_{\lambda_1}$ (deg E of N, mod 180)')
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 45, 90, 135, 180])
    ax.grid(alpha=0.3)
    ax.set_title('Principal-axis azimuth (orientation)')

    # (c) Azimuthal-fit quality
    ax = fig.add_subplot(gs[0, 2], sharey=fig.axes[0])
    ax.plot(r2, DEPTH, 'tab:green', lw=1.5)
    ax.set_xlabel(r'azimuthal-model $R^2$')
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)
    ax.set_title('Uniform-fabric fit quality')

    # (d) Sample count
    ax = fig.add_subplot(gs[0, 3], sharey=fig.axes[0])
    ax.plot(nfit, DEPTH, 'tab:gray', lw=1.5)
    ax.set_xlabel('blocks in fit')
    ax.grid(alpha=0.3)
    ax.set_title('Azimuth samples')

    # (e) The money plots: dlam vs heading at example depths, with fit
    for i, dz in enumerate(EXAMPLE_DEPTHS):
        ax = fig.add_subplot(gs[1, i])
        iz = int(np.argmin(np.abs(DEPTH - dz)))
        y = profs[:, iz]
        m = np.isfinite(y)
        ax.scatter(heads[m], y[m], s=8, alpha=0.5, color='tab:purple')
        if np.isfinite(P[iz]):
            phi = np.linspace(0, 180, 361)
            ax.plot(phi, -P[iz] * np.cos(2 * np.radians(phi - theta[iz])),
                    'k-', lw=1.5)
            ax.set_title(f'{DEPTH[iz]:.0f} m: P={P[iz]:.3f}, '
                         f'$\\theta$={theta[iz]:.0f}$^\\circ$', fontsize=9)
        else:
            ax.set_title(f'{DEPTH[iz]:.0f} m: no fit', fontsize=9)
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 90, 180])
        ax.set_ylim(-0.15, 0.15)
        ax.axhline(0, color='gray', lw=0.5)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.set_ylabel(r'$\Delta\lambda$ measured')
        ax.set_xlabel('drive heading (deg)')

    fig.suptitle('Ridge A azimuthal inversion: horizontal fabric orientation and '
                 'strength vs depth\n'
                 r'per-depth fit of $\Delta\lambda(\varphi) = -P\,\cos 2(\varphi - \theta)$'
                 ' over all grid headings (grid-uniform fabric assumption)',
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, 'ridge_a_azimuthal.png')
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('saved', out)


if __name__ == '__main__':
    main()

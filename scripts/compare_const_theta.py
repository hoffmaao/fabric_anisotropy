#!/usr/bin/env python3
"""Compare the constant-axis (_ct) quad-pol products against the per-window ones.

The _ct batch answers one question: does holding the fabric axis constant in
depth, per along-track segment, describe these sites better than letting every
depth window find its own axis?

Constant mode changes the SEGMENT profiles the blocks inherit (ls_theta_seg)
and, through them, every block fit (sec_*). It never touches the frame-pooled
pass on a multi-segment frame, so ls_theta0_geo / ls_resid are the wrong
fields to compare: an earlier version of this script read them and reported
the two products identical when the block residual had moved by 0.05.

Three numbers decide it, and they are not interchangeable:

  th_spread_seg     Contrast-weighted circular SD of the per-window axis minima
                    BEFORE the axis is held, per segment. This is the test of
                    the ASSUMPTION. Small means the windows already agreed and
                    holding them costs nothing; large means the column
                    genuinely rotates and the held axis is an average of a
                    rotation.
  sec_resid_ls      Block fit residual. Holding an axis removes a free parameter
                    per window, so the residual can only rise. How MUCH it
                    rises is the price of the assumption.
  theta_const_q     Pooled contrast. Reports how sharply the pooled curve is
                    peaked, NOT whether the windows agreed - measured higher on
                    a rotating synthetic than a constant one. Never read it as
                    evidence the assumption holds; that is th_spread_seg.
  se_th             Jackknife SE of the held axis: its RANDOM error under
                    sub-block resampling only. At an isotropic site (McMurdo,
                    dlam 0.001) it is 1-5 deg while the axis itself is
                    meaningless - a pedestal-driven minimum is perfectly
                    repeatable. Read it with dlam and the residual.

Runs where the products live (mem1), since the products are HDF5 and stay
server-side:

  python3 scripts/compare_const_theta.py [stage_dir]   # default <work>/stages/quadpol
"""
import os
import sys
import glob
import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("needs h5py; run this where the products live (mem1)")

# The work root the server scripts use (opr_fabric/server/fabric_paths.m):
# FABRIC_ROOT if set, else the parent of the checkout this file sits in.
_WORK = os.environ.get("FABRIC_ROOT") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_STAGE = os.path.join(_WORK, "stages", "quadpol")
Z_BAND = (200.0, 1500.0)     # quotable fabric starts ~200 m (co-pol reference offset)

# Verdict thresholds on the contrast-weighted per-window spread. PROVISIONAL:
# the synthetic in test_quadpol_const_theta gives 2 deg constant vs 25 deg
# rotating; the real-data calibration is Ridge A (a divide, expected to hold)
# against a Thwaites margin control (known to rotate), and these numbers are
# to be reset from that comparison once both have run.
SPREAD_HOLDS = 8.0
SPREAD_ROTATES = 20.0
RESID_COST_OK = 0.03

# Season tag -> site. The season NUMBER never matches the tag year, and the
# 2024 tags split three ways, so this cannot be a prefix-to-season shortcut.
SITES = [
    ("20221206", "Taylor Dome"), ("20221207", "Taylor Dome"),
    ("20221209", "Taylor Dome"), ("20221210", "Taylor Dome"),
    ("20221211", "Taylor Dome"),
    ("20240120", "WAIS Divide"),
    ("20240202", "McMurdo"), ("20240203", "McMurdo"),
    ("202401", "Thwaites"),
    ("20250107", "Ridge A"), ("20250108", "Ridge A"), ("20250111", "Ridge A"),
    ("20250112", "Ridge A"), ("20250113", "Ridge A"), ("20250115", "Ridge A"),
    ("20250117", "Ridge A"),
    ("20260106", "Eastwind"), ("20260107", "Eastwind"),
    ("20260109", "Eastwind"), ("20260119", "Eastwind"), ("20260121", "Eastwind"),
]


def site_of(tag):
    for pre, name in SITES:
        if tag.startswith(pre):
            return name
    return "unknown"


def circ_mean_deg(th, w=None):
    """Weighted circular mean of a mod-180 axis via the doubled angle. Never unwrap."""
    th = np.asarray(th, float).ravel()
    w = np.ones_like(th) if w is None else np.asarray(w, float).ravel()
    g = np.isfinite(th) & np.isfinite(w) & (w > 0)
    if not g.any():
        return np.nan
    ph = np.sum(w[g] * np.exp(2j * np.radians(th[g])))
    return np.degrees(np.angle(ph)) / 2 % 180


def circ_sd_deg(th):
    """Circular SD (deg) of a mod-180 axis over the doubled angle."""
    th = np.asarray(th, float).ravel()
    g = np.isfinite(th)
    if g.sum() < 2:
        return np.nan
    R = abs(np.mean(np.exp(2j * np.radians(th[g]))))
    return np.degrees(np.sqrt(max(-2 * np.log(max(R, 1e-300)), 0))) / 2


def circ_diff_deg(a, b):
    """Signed mod-180 difference, in (-90, 90]."""
    if not (np.isfinite(a) and np.isfinite(b)):
        return np.nan
    return (a - b + 90) % 180 - 90


def get(res, key, default=np.nan):
    if key not in res:
        return np.array([default], float)
    return np.array(res[key]).astype(float)


def read(fn):
    """One product's numbers over the quotable band.

    Field names are not interchangeable: the LS block dlam is sec_dlam_ls (NOT
    ls_dlam, which is the frame pass, and NOT dlam_ls, which is a scalar). The
    segment profiles are ls_theta_seg [nseg x Nw] against ls_zw; the block
    fields are [nb x Nt] against z. th_spread_seg / held_seg exist only in
    _ct products.
    """
    with h5py.File(fn, "r") as f:
        r = f["res"]
        zw = get(r, "ls_zw").ravel()
        z = get(r, "z").ravel()
        th_seg = np.atleast_2d(get(r, "ls_theta_seg"))
        if th_seg.shape[-1] != zw.size:
            th_seg = th_seg.T
        q_seg = np.atleast_2d(get(r, "ls_q_seg"))
        if q_seg.shape[-1] != zw.size:
            q_seg = q_seg.T
        mw = (zw >= Z_BAND[0]) & (zw < Z_BAND[1])
        sec_dlam = np.atleast_2d(get(r, "sec_dlam_ls"))
        sec_res = np.atleast_2d(get(r, "sec_resid_ls"))
        if sec_dlam.shape[-1] != z.size:
            sec_dlam = sec_dlam.T
            sec_res = sec_res.T
        mz = (z >= Z_BAND[0]) & (z < Z_BAND[1])
        deep = (z >= 1150) & (z < Z_BAND[1])
        # the axis the blocks inherit: contrast-weighted over segments and
        # the quotable band; its depth spread per segment says how much the
        # free profile wandered before anything was held
        th_band = th_seg[:, mw]
        q_band = q_seg[:, mw]
        depth_sd = np.nanmedian([circ_sd_deg(row) for row in th_band]) if th_band.size else np.nan
        spread = get(r, "th_spread_seg").ravel()
        held = get(r, "held_seg", 0).ravel()
        # uncertainty fields (products built with the jackknife / split-half):
        # segment-level standard errors on the window grid, block-level
        # pooled split-half sigma on z
        se_th = np.atleast_2d(get(r, "ls_se_theta_seg"))
        if se_th.shape[-1] != zw.size and se_th.size > 1:
            se_th = se_th.T
        # a NaN in se_th is not one thing: the axis SE is bounded, so it
        # abstains where the replicates scattered past what it can resolve
        # (ls_se_theta_sat) as well as where there were too few of them.
        # unres is NOT comparable across frames of different length - the
        # resolution floor tightens with the sub-block count, so a longer
        # segment abstains at a smaller true scatter (17.2 deg at 6
        # sub-blocks, 12.2 at 11, 8.4 at 22). Read it against ls_se_theta_n,
        # never as a measure of ice quality on its own. See ptt.circAxisSE.
        # Reporting only the median of the finite cells would quietly drop
        # the unresolvable ones and read as a tighter error bar than the
        # data support, so the abstained fraction is carried alongside.
        sat = np.atleast_2d(get(r, "ls_se_theta_sat", 0)).astype(bool)
        if sat.shape != se_th.shape:
            sat = np.zeros_like(se_th, dtype=bool)
        # the sub-block count unres has to be read against, printed beside it
        n_rep = np.atleast_2d(get(r, "ls_se_theta_n", np.nan)).astype(float)
        if n_rep.shape != se_th.shape:
            n_rep = np.full_like(se_th, np.nan, dtype=float)
        se_dl = np.atleast_2d(get(r, "ls_se_dlam_seg"))
        if se_dl.shape[-1] != zw.size and se_dl.size > 1:
            se_dl = se_dl.T
        se_blk = np.atleast_2d(get(r, "sec_se_dlam_ls"))
        if se_blk.shape[-1] != z.size and se_blk.size > 1:
            se_blk = se_blk.T
        return dict(
            se_theta=np.nanmedian(se_th[:, mw]) if se_th.shape[-1] == zw.size else np.nan,
            se_theta_unres=(float(np.mean(sat[:, mw]))
                            if sat.shape[-1] == zw.size and sat.size else np.nan),
            se_theta_n=(np.nanmedian(n_rep[:, mw])
                        if n_rep.shape[-1] == zw.size and n_rep.size else np.nan),
            se_dlam=np.nanmedian(se_dl[:, mw]) if se_dl.shape[-1] == zw.size else np.nan,
            se_blk=np.nanmedian(se_blk[:, mz]) if se_blk.shape[-1] == z.size else np.nan,
            nseg=int(th_seg.shape[0]),
            th=circ_mean_deg(th_band, np.where(np.isfinite(q_band), q_band, 0) + 1e-6),
            depth_sd=depth_sd,
            dlam=np.nanmedian(sec_dlam[:, mz]) if sec_dlam.size > 1 else np.nan,
            dlam_deep=np.nanmedian(sec_dlam[:, deep]) if sec_dlam.size > 1 else np.nan,
            resid=np.nanmedian(sec_res[:, mz]) if sec_res.size > 1 else np.nan,
            finite=np.isfinite(sec_dlam[:, mz]).mean() if sec_dlam.size > 1 else np.nan,
            spread=np.nanmedian(spread) if np.isfinite(spread).any() else np.nan,
            n_held=int(np.nansum(held)),
        )


DLAM_ISO = 0.01        # below this the windows carry no axis (dlam_min_theta)
RESID_FAILS = 0.40     # block residual at which the single-column model is not fitting


def verdict(spr, dres, dlam, resid):
    """Read the three numbers together, not the spread alone.

    The first batch showed why: McMurdo has spreads near 50 deg with
    residuals of 0.1 and dlam of 0.001 - isotropic ice, where the per-window
    minima are noise around a meaningless axis and holding it costs nothing;
    Taylor Dome has the same spreads with residuals of 0.6 in BOTH modes -
    the model is not describing that ice either way. Neither is a rotating
    axis. The Thwaites controls are: spreads of 30-40 deg AND a residual
    cost of +0.05 from holding. So: no fabric -> no axis; no fit -> no
    verdict; otherwise the cost of holding is the discriminant and the
    spread the corroboration.
    """
    if np.isfinite(dlam) and dlam < DLAM_ISO:
        return "ISOTROPIC - no axis to hold (dlam < %.2f)" % DLAM_ISO
    if np.isfinite(resid) and resid > RESID_FAILS:
        return "FIT FAILS - model does not describe this ice in either mode"
    if not np.isfinite(spr):
        return "no spread reported"
    if np.isfinite(dres) and dres >= RESID_COST_OK:
        return "ROTATES - holding costs fit (or apparent rotation; see notes)"
    if spr < SPREAD_HOLDS:
        return "HOLDS - windows already agreed"
    if spr < SPREAD_ROTATES:
        return "marginal - inspect residual profile"
    return "undetermined - axes scatter but holding costs nothing"


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_STAGE
    ct_files = sorted(glob.glob(os.path.join(stage, "quadpol_section_*_ct.mat")))
    if not ct_files:
        sys.exit("no _ct products in %s - has the batch run?" % stage)

    rows, orphans = [], []
    for cf in ct_files:
        base = cf[:-len("_ct.mat")] + ".mat"
        tag = os.path.basename(cf)[len("quadpol_section_"):-len("_ct.mat")]
        if not os.path.exists(base):
            orphans.append(tag)
            continue
        try:
            a, b = read(base), read(cf)
        except (OSError, KeyError, ValueError) as e:
            orphans.append("%s (unreadable: %s)" % (tag, e))
            continue
        rows.append((site_of(tag), tag, a, b))

    print("band %.0f-%.0f m; th = contrast-weighted axis the blocks inherit; "
          "sd_z = depth spread of the FREE segment profiles; spread = per-window "
          "minima spread before holding (_ct); resid/dlam/finite = block fits"
          % Z_BAND)
    print("%-11s %-16s %4s %4s  %6s %6s %6s %5s  %6s  %6s %6s %+7s  %6s %6s  %5s %5s"
          % ("site", "frame", "nseg", "held", "th_def", "th_ct", "dth", "sd_z",
             "spread", "res_df", "res_ct", "dres", "dl_def", "dl_ct", "fin_d", "fin_c"))
    print("-" * 128)
    by_site = {}
    for site, tag, a, b in rows:
        dth = circ_diff_deg(b["th"], a["th"])
        dres = b["resid"] - a["resid"]
        print("%-11s %-16s %4d %4d  %6.1f %6.1f %+6.1f %5.1f  %6.1f  %6.3f %6.3f %+7.3f  %6.3f %6.3f  %5.2f %5.2f"
              % (site, tag, b["nseg"], b["n_held"], a["th"], b["th"], dth, a["depth_sd"],
                 b["spread"], a["resid"], b["resid"], dres, a["dlam"], b["dlam"],
                 a["finite"], b["finite"]))
        by_site.setdefault(site, []).append((a, b, dth, dres))
        # per-frame verdict too: sites mix regimes (WAIS Divide is one strong
        # divide frame and one isotropic frame; Ridge A holds on some lines
        # and costs 0.03-0.04 on others)
        print("%-11s %-16s   -> %s" % ("", "", verdict(b["spread"], dres, a["dlam"], a["resid"])))

    # what the numbers above are worth: the _ct products' own standard
    # errors, where the products carry them (segment jackknife, block
    # split-half); a dth well inside se_th is no change at all
    if any(np.isfinite(b["se_dlam"]) for _, _, _, b in rows):
        print("\nuncertainty (_ct products, band medians): se_th = held-axis SE (deg), "
              "n = sub-blocks the SE was formed from, "
              "unres = fraction of cells where the axis SE could not be resolved. "
              "unres is only comparable at equal n: the abstention threshold scales "
              "as ~1.343/(2*sqrt(n-1)) rad, so 17.2 deg of replicate scatter is "
              "resolvable at n = 6 but 8.4 deg is not at n = 22 - see ptt.circAxisSE. "
              "se_dl = segment dlam SE, se_blk = pooled block split-half sigma")
        print("%-11s %-16s %7s %4s %6s %7s %7s"
              % ("site", "frame", "se_th", "n", "unres", "se_dl", "se_blk"))
        for site, tag, a, b in rows:
            unres = b.get("se_theta_unres", float("nan"))
            n_rep = b.get("se_theta_n", float("nan"))
            print("%-11s %-16s %7.2f %4s %5.0f%% %7.4f %7.4f"
                  % (site, tag, b["se_theta"],
                     "--" if not np.isfinite(n_rep) else "%d" % round(n_rep),
                     100 * unres, b["se_dlam"], b["se_blk"]))

    print("\n%-11s %6s  %8s  %6s  %7s  %7s  %7s  %s"
          % ("SITE", "frames", "med|dth|", "sd_z", "spread", "d_resid", "d_dlam", "assumption"))
    print("-" * 96)
    for site in sorted(by_site):
        v = by_site[site]
        adth = np.nanmedian([abs(d) for _, _, d, _ in v])
        sdz = np.nanmedian([a["depth_sd"] for a, _, _, _ in v])
        spr = np.nanmedian([b["spread"] for _, b, _, _ in v])
        dres = np.nanmedian([d for _, _, _, d in v])
        ddl = np.nanmedian([b["dlam"] - a["dlam"] for a, b, _, _ in v])
        dl_site = np.nanmedian([a["dlam"] for a, _, _, _ in v])
        res_site = np.nanmedian([a["resid"] for a, _, _, _ in v])
        print("%-11s %6d  %8.1f  %6.1f  %7.1f  %+7.4f  %+7.4f  %s"
              % (site, len(v), adth, sdz, spr, dres, ddl, verdict(spr, dres, dl_site, res_site)))

    if orphans:
        print("\n%d _ct product(s) with no default counterpart to compare:"
              % len(orphans))
        for o in orphans:
            print("  %s" % o)


if __name__ == "__main__":
    main()

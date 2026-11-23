"""What the node products can say about fabric with depth, and whether
the differences between blocks are a field over the survey or the solver.

    python scripts/figures/node_representation.py \\
        [<mirror of stages/node>] [out_dir] [pattern]

Reads every node_<tag>_B2000.mat (or the products pattern names, e.g.
node_*_B500.mat for 500 m blocks), redoes each block's linear solve under
the depth model (node_solve.SMOOTH: a smooth dlam(z) about a free level
plus a bounded window phase, nodes counted in the band's own frequency),
and answers three questions with numbers that do not pass through the
prior:

  DEPTH    How many numbers does a column carry, and at what vertical
           scale does the answer stop depending on the model? The plain
           model (node_solve.POSTER: every phase increment taken as a
           delay, no window phase, which hands the reflectors' phase
           back as 10 m stripes), the depth model and the same model at
           half and twice its correlation length are compared after
           running means of increasing length; the prior-free interval
           means (summed phase increments, node_solve.intervals) are the
           reference.
  STRIPES  Is the 10 m structure a delay? A delay accumulates down the
           column; a phase the reflectors put between the polarisations
           does not. The increment residuals' autocorrelation and the
           growth of the residual delay with span tell them apart.
  LATERAL  Are the differences between blocks real? Three tests that the
           solver cannot fake: independent passes over the same ice, the
           cross-polarised POWER oscillation against the delay the PHASE
           gives (two data types), and predicting each frame from all
           the others.

  HORIZON  The change near 470 m is sharper than the depth model
           resolves. One block's evidence for a ramp there is a few log
           units, so it is not part of the estimator; pooled over the
           survey it is decisive, and here each block's ramp depth and
           step are found with the smooth model's own noise level.

Writes node_depth_model.png, node_lateral_field.png and
node_resolved.npz (profiles, errors, resolution, interval means, horizon
depth and step per block).

Nothing is written on the panels but their letters; the numbers are
printed, and these are the captions.

node_depth_model.png
  A  One block (the one nearest the survey's centre, named in the
     output). Grey: the plain model, every phase increment taken as a
     delay. Black with its one-sigma band: the depth model, which is
     what the product stores. Blue: the prior-free means over 125 m,
     from the summed phase increments. Orange stairs: dlam from the
     spacing of the nodes alone.
  B  How much the answer depends on the model, against the vertical
     averaging length. Blue: the spread between blocks, which is the
     signal. Grey: plain model minus depth model. Black: the depth
     model at half (squares) and twice (triangles) its correlation
     length. Orange dashed: what the window phase leaves on a mean.
  C  Black: the residual delay the depth model leaves on the phase,
     against span. Dashed: how a delay (fabric) would accumulate. The
     residuals' autocorrelation at one, two and three window steps is
     printed; a bounded phase, differenced, gives -0.31, -0.12, -0.04.
  D  Log evidence over the survey against the frequency the nodes are
     counted in. Dashed: the LO, where the carrier phase counts. Band:
     the band's centre as measured on 20250108_02_003, 734-738 MHz.

node_lateral_field.png
  A  The blocks in map view, coloured by dlam over 500-1450 m; grey
     lines join the blocks of a frame.
  B  The same ice on another pass: 125 m means of block pairs under
     0.6 km apart (colour: depth, dark shallow).
  C  Rms difference of dlam (500-1450 m) against separation. Black: all
     pairs of different frames. Orange: pairs separated east-west.
     Blue: pairs separated north-south.
  D  Each frame's dlam (500-1450 m) predicted from all the other frames.
  E  The fringe rate of the cross-polarised POWER over the delay rate of
     the PHASE, two data types, against strength. Blue: the band's
     centre over the LO, which is what they should give. Dashed: what
     the points would follow if every block had the weakest fabric.
  F  Black: the depth pattern of the one map that carries most of the
     variance between blocks. Grey: the survey's mean profile. Both
     relative to their means.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)
import morrow_style as ms  # noqa: E402
import node_export as ne  # noqa: E402
import node_solve as ns  # noqa: E402

R_E = 6371.0
SPAN = 125.0
EDGES = np.arange(200.0, 1451.0, SPAN)
ZC = (EDGES[:-1] + EDGES[1:]) / 2
DEEP = ZC > 500                    # below the horizon near 470 m
BAND = (250.0, 1400.0)             # depths the profile comparisons use
SEP_BINS = [(0, 0.6), (0.6, 1.5), (1.5, 2.5), (2.5, 4), (4, 6), (6, 9),
            (9, 13), (13, 18), (18, 25)]
C_NS, C_EW = ne.C_NS, ne.C_EW


def rms(v):
    v = np.asarray(v, float)
    return np.sqrt(np.nanmean(v ** 2))


def load(src, pattern="node_*_B2000.mat"):
    rows = []
    for fn in sorted(glob.glob(os.path.join(src, pattern))):
        rows += [r for r in ne.read_product(fn) if r.get("ok")]
    return rows


def local_xy(rows):
    """km east and north of the survey centroid"""
    lat = np.array([r["lat"] for r in rows])
    lon = np.array([r["lon"] for r in rows])
    y = np.radians(lat - lat.mean()) * R_E
    x = np.radians(lon - lon.mean()) * R_E * np.cos(np.radians(lat.mean()))
    return x, y


def running_mean(v, n):
    """mean over n cells along the last axis, the ends over what there is"""
    c = np.cumsum(np.insert(v, 0, 0.0, axis=-1), axis=-1)
    h = n // 2
    out = np.empty_like(v)
    for i in range(v.shape[-1]):
        a, b = max(0, i - h), min(v.shape[-1], i + h + 1)
        out[..., i] = (c[..., b] - c[..., a]) / (b - a)
    return out


# ------------------------------------------------------------- depth
def survey_evidence(rows):
    """log evidence summed over the blocks, relative to the plain model"""
    def total(model, fc):
        return sum(ns.resolve(r, model, node_fc=fc)["logev"] for r in rows)
    ref = total(ns.POSTER, ns.FC)
    freqs = np.arange(720.0, 766.0, 5.0)
    curve = np.array([total(ns.POSTER, f * 1e6) for f in freqs]) - ref
    return dict(freqs=freqs, curve=curve,
                poster_band=total(ns.POSTER, ns.NODE_FC) - ref,
                smooth=total(ns.SMOOTH, ns.NODE_FC) - ref)


def residual_statistics(rows, outs):
    """Autocorrelation of the increment residuals at 25..200 m and the rms
    of the residual delay over a span, pooled over the survey: a bounded
    phase is flat with span and negative at one lag, a delay grows."""
    zw = rows[0]["zw"]
    res = np.full((len(rows), zw.size - 1), np.nan)
    for b, (r, o) in enumerate(zip(rows, outs)):
        res[b, ns.window_index(r)] = o["res"]
    zc = (zw[:-1] + zw[1:]) / 2
    res = res[:, (zc > BAND[0]) & (zc < BAND[1])]
    acf = []
    for lag in range(1, 9):
        a, b = res[:, :-lag], res[:, lag:]
        ok = np.isfinite(a) & np.isfinite(b)
        acf.append(np.sum(a[ok] * b[ok])
                   / np.sqrt(np.sum(a[ok] ** 2) * np.sum(b[ok] ** 2)))
    spans = np.array([1, 2, 3, 4, 6, 8, 12, 16])
    sf = []
    filled = np.where(np.isfinite(res), res, 0.0)
    c = np.cumsum(np.insert(filled, 0, 0.0, axis=1), axis=1)
    gaps = np.cumsum(np.insert(~np.isfinite(res), 0, 0, axis=1), axis=1)
    for n in spans:
        d = c[:, n:] - c[:, :-n]
        ok = (gaps[:, n:] - gaps[:, :-n]) == 0
        sf.append(np.sqrt(np.mean(d[ok] ** 2)))
    step = np.median(np.diff(zw))
    return dict(acf=np.array(acf), lags=step * np.arange(1, 9),
                spans=step * spans, sf=np.array(sf), res=res)


def model_dependence(rows, outs, plain, zm):
    """rms difference between solutions after a running mean of each
    length: the plain model and the depth model at half and twice its
    correlation length, against the depth model; and the survey's own
    lateral spread at that averaging."""
    sel = (zm >= BAND[0]) & (zm <= BAND[1])
    base = np.array([o["m"] for o in outs])
    stored = np.array([o["m"] for o in plain])
    alt = {}
    for length in (250.0, 1000.0):
        mdl = dict(ns.SMOOTH, len_m=length)
        alt[length] = np.array([ns.resolve(r, mdl)["m"] for r in rows])
    lengths = np.array([1, 3, 5, 7, 11, 15, 21, 31])
    out = dict(lengths=lengths * ns.MODEL_DZ, plain=[], half=[], twice=[],
               lateral=[])
    for n in lengths:
        b = running_mean(base, n)
        out["plain"].append(rms((running_mean(stored, n) - b)[:, sel]))
        out["half"].append(rms((running_mean(alt[250.0], n) - b)[:, sel]))
        out["twice"].append(rms((running_mean(alt[1000.0], n) - b)[:, sel]))
        out["lateral"].append(np.sqrt(np.mean(np.var(b[:, sel], axis=0))))
    return {k: np.asarray(v) for k, v in out.items()}


def node_spacing_stairs(row):
    """dlam between successive nodes from their spacing alone"""
    z = np.r_[0.0, row["node_z"]]
    tau = np.r_[0.0, row["node_n"] / (2 * ns.NODE_FC * 1e-9)]
    return z, np.diff(tau) / (row["k"] * np.diff(z))


def horizon(rows, outs, z_range=(380.0, 600.0), width=30.0, sd=0.05):
    """One ramp of `width` added to each block's smooth prior, its depth
    where the block's evidence is largest in z_range. Returns the depth,
    the step (posterior mean of the ramp's amplitude) and the evidence
    gain with the depth marginalised over the range."""
    zb = np.arange(z_range[0], z_range[1] + 1e-9, ns.MODEL_DZ)
    depth, step, gain = [], [], []
    for r, o in zip(rows, outs):
        zm = r["zm"]
        g, d, sd_abs = ns.block_rows(r)
        cd = ns.data_cov(sd_abs, r["inc_s"], ns.window_index(r), r["zw"],
                         o["psi_sd_ns"], o["model"]["psi_len_m"])
        cm = ns.kernel(zm, o["model"])
        mp = np.zeros(zm.size)
        ramps = [np.clip((zm + ns.MODEL_DZ / 2 - (z - width / 2)) / width,
                         0.0, 1.0) for z in zb]
        le = np.array([ns.log_evidence(
            g, d, cd, cm + sd ** 2 * np.outer(q, q), mp) for q in ramps])
        j = int(np.argmax(le))
        q = ramps[j]
        alpha = ns.solve(g, d, cd, cm + sd ** 2 * np.outer(q, q),
                         mp)["alpha"]
        depth.append(zb[j])
        step.append(sd ** 2 * (g @ q) @ alpha)
        gain.append(le[j] + np.log(np.mean(np.exp(le - le[j])))
                    - ns.log_evidence(g, d, cd, cm, mp))
    return np.array(depth), np.array(step), np.array(gain)


# ----------------------------------------------------------- lateral
def interval_matrix(rows, outs):
    return np.array([ns.intervals(r, o["psi_sd_ns"], edges=EDGES)["dlam"]
                     for r, o in zip(rows, outs)])


def power_phase_scale(row):
    """The 45-deg cross/co POWER ratio oscillates as cos(2 pi f tau). With
    tau from the PHASE increments alone and f the LO, find the factor s
    on tau that best explains the power: s = f_band / f_LO when the phase
    delay is the delay the power sees, whatever the block's strength."""
    zm, k = row["zm"], row["k"]
    m_inc = ns.resolve(row, ns.POSTER, node_fc=ns.FC, use_nodes=False)["m"]
    tau_cells = np.cumsum(m_inc * k * ns.MODEL_DZ)
    z, v = row["z"], row["ratio45"]
    n = max(1, int(round(20.0 / np.median(np.diff(z)))))
    ok = np.isfinite(v)
    vs = (np.convolve(np.where(ok, v, 0.0), np.ones(n), "same")
          / np.maximum(np.convolve(ok.astype(float), np.ones(n), "same"), 1))
    zz, vv = z[::n // 2], vs[::n // 2]
    m = (zz >= 300) & (zz <= 1450) & np.isfinite(vv)
    zz, vv = zz[m], vv[m]
    tau = np.interp(zz, zm + ns.MODEL_DZ, tau_cells)
    t = (zz - zz.mean()) / 500.0
    a0 = np.c_[np.ones_like(t), t, t ** 2]
    r0 = vv - a0 @ np.linalg.lstsq(a0, vv, rcond=None)[0]
    scales = np.arange(0.5, 1.6, 0.005)
    r2 = np.empty(scales.size)
    for i, s in enumerate(scales):
        ph = 2 * np.pi * s * ns.FC * 1e-9 * tau
        a = np.c_[a0, np.cos(ph), np.sin(ph)]
        res = vv - a @ np.linalg.lstsq(a, vv, rcond=None)[0]
        r2[i] = 1 - np.sum(res ** 2) / np.sum(r0 ** 2)
    i = int(np.argmax(r2))
    return scales[i], r2[i]


def same_ice_pairs(x, y, tags, max_km=0.6):
    """blocks of different frames whose centres lie within max_km"""
    d = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    i, j = np.triu_indices(x.size, 1)
    m = (tags[i] != tags[j]) & (d[i, j] < max_km)
    return i[m], j[m]


def variogram(v, x, y, tags, sel=None, circular=False):
    """sqrt of the semivariance by separation, pairs of different frames"""
    i, j = np.triu_indices(x.size, 1)
    d = np.hypot(x[i] - x[j], y[i] - y[j])
    keep = tags[i] != tags[j]
    if sel is not None:
        keep &= sel(np.abs(x[i] - x[j]), np.abs(y[i] - y[j]))
    out, cnt = [], []
    for lo, hi in SEP_BINS:
        m = keep & (d >= lo) & (d < hi)
        dv = v[i[m]] - v[j[m]]
        if circular:
            dv = (dv + 90) % 180 - 90
        ok = np.isfinite(dv)
        cnt.append(int(ok.sum()))
        out.append(np.sqrt(np.mean(dv[ok] ** 2) / 2) if ok.sum() > 5
                   else np.nan)
    return np.array(out), np.array(cnt)


def krige(xt, yt, vt, xp, yp, rng_km=6.0, nugget=0.02):
    """universal kriging with a plane, exponential covariance"""
    def cov(xa, ya, xb, yb):
        return np.exp(-np.hypot(xa[:, None] - xb[None, :],
                                ya[:, None] - yb[None, :]) / rng_km)
    ki = np.linalg.inv(cov(xt, yt, xt, yt) + nugget * np.eye(xt.size))
    ft = np.c_[np.ones_like(xt), xt, yt]
    fp = np.c_[np.ones_like(xp), xp, yp]
    beta = np.linalg.solve(ft.T @ ki @ ft, ft.T @ ki @ vt)
    return fp @ beta + cov(xp, yp, xt, yt) @ ki @ (vt - ft @ beta)


def leave_group_out(v, x, y, groups, circular=False):
    """each group's blocks predicted from every other group: the kriged
    value and the others' mean"""
    ok = np.isfinite(v)
    pred = np.full(v.size, np.nan)
    flat = np.full(v.size, np.nan)
    for gname in np.unique(groups):
        te = (groups == gname) & ok
        tr = (groups != gname) & ok
        if te.any():
            pred[te] = krige(x[tr], y[tr], v[tr], x[te], y[te])
            flat[te] = v[tr].mean()
    e, e0 = pred - v, flat - v
    if circular:
        e, e0 = (e + 90) % 180 - 90, (e0 + 90) % 180 - 90
    return pred, rms(e), rms(e0)


def variance_budget(im):
    """depth against lateral variance of the interval means, and how much
    of the lateral part one map carries"""
    full = im[np.all(np.isfinite(im), axis=1)]
    prof = full.mean(0)
    dev = full - prof
    s = np.linalg.svd(dev, compute_uv=False)
    vt = np.linalg.svd(dev, full_matrices=False)[2]
    total = np.mean((full - full.mean()) ** 2)
    shape = vt[0] * np.sign(vt[0].sum())
    return dict(n=full.shape[0], profile=prof, depth=np.var(prof) / total,
                lateral=np.mean(dev ** 2) / total,
                modes=s ** 2 / np.sum(s ** 2), shape=shape / shape.mean())


# ----------------------------------------------------------- figures
def fig_depth(rows, outs, plain, ev, rs, md, pick, out_png):
    fig, ax = plt.subplots(2, 2, figsize=(12.8, 10.4),
                           constrained_layout=True)
    r, o = rows[pick], outs[pick]
    zm = r["zm"] + ns.MODEL_DZ / 2
    a = ax[0, 0]
    a.plot(plain[pick]["m"], zm, "-", color="0.6", lw=1.0)
    a.fill_betweenx(zm, o["m"] - o["sigma"], o["m"] + o["sigma"],
                    color="0.2", alpha=0.25, lw=0)
    a.plot(o["m"], zm, "-", color="k", lw=2.0)
    iv = ns.intervals(r, o["psi_sd_ns"], edges=EDGES)
    a.errorbar(iv["dlam"], ZC, xerr=iv["sigma"], yerr=SPAN / 2, fmt="o",
               color=C_NS, ms=5, lw=1.2, capsize=0)
    zs, ds = node_spacing_stairs(r)
    a.plot(np.repeat(ds, 2), np.repeat(zs, 2)[1:-1], "-", color=C_EW, lw=1.6)
    a.set_ylim(1500, 0)
    a.set_xlim(0, 0.12)
    a.set_xlabel(r"$\Delta\lambda$ along the slow axis")
    a.set_ylabel("depth [m]")
    ms.tag(a, "A", x=0.93)

    a = ax[0, 1]
    a.loglog(md["lengths"], md["lateral"], "-", color=C_NS, lw=2.2)
    a.loglog(md["lengths"], md["plain"], "o-", color="0.5", lw=1.6, ms=5)
    a.loglog(md["lengths"], md["half"], "s-", color="k", lw=1.6, ms=5)
    a.loglog(md["lengths"], md["twice"], "^-", color="k", lw=1.6, ms=5,
             mfc="w")
    a.loglog(rs["spans"], rs["sf"] / (rows[0]["k"] * rs["spans"]), "--",
             color=C_EW, lw=1.6)
    a.set_xlabel("vertical averaging length [m]")
    a.set_ylabel(r"rms difference in $\Delta\lambda$")
    a.set_xlim(9, 420)
    a.set_ylim(1e-4, 3e-2)
    a.set_xticks([10, 30, 100, 300])
    a.set_xticklabels(["10", "30", "100", "300"])
    a.set_xticks([], minor=True)
    ms.tag(a, "B", x=0.93)

    a = ax[1, 0]
    a.plot(rs["spans"], rs["sf"], "o-", color="k", lw=2.0, ms=6)
    a.plot(rs["spans"], rs["sf"][0] * np.sqrt(rs["spans"] / rs["spans"][0]),
           "--", color="0.5", lw=1.4)
    a.set_xlabel("span [m]")
    a.set_ylabel("rms residual delay over the span [ns]")
    a.set_ylim(0, 0.055)
    ms.tag(a, "C")

    a = ax[1, 1]
    a.plot(ev["freqs"], ev["curve"], "o-", color="k", lw=2.0, ms=6)
    a.axvline(ns.FC / 1e6, color="0.5", lw=1.2, ls="--")
    a.axvspan(734, 738, color=C_NS, alpha=0.25, lw=0)
    a.set_xlabel("frequency the nodes are counted in [MHz]")
    a.set_ylabel("log evidence, %d blocks, against 750 MHz" % len(rows))
    ms.tag(a, "D")
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def fig_lateral(rows, x, y, deep, im, pairs, vg, loo, scale, budget,
                out_png):
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 10.2),
                           constrained_layout=True)
    tags = np.array([r["tag"] for r in rows])
    a = ax[0, 0]
    for t in np.unique(tags):
        m = tags == t                      # blocks are in along-track order
        a.plot(x[m], y[m], "-", color="0.7", lw=0.8, zorder=1)
    sc = a.scatter(x, y, c=deep, s=46, cmap="viridis", vmin=0.05, vmax=0.10,
                   edgecolor="k", lw=0.3, zorder=2)
    a.set_aspect("equal")
    a.set_xlabel("east [km]")
    a.set_ylabel("north [km]")
    cb = fig.colorbar(sc, ax=a, shrink=0.8, pad=0.02)
    cb.set_label(r"$\Delta\lambda$, 500-1450 m")
    ms.tag(a, "A")

    a = ax[0, 1]
    i, j = pairs
    for k in range(ZC.size):
        a.plot(im[i, k], im[j, k], "o", ms=4.5, alpha=0.85,
               color=plt.cm.viridis(k / (ZC.size - 1)))
    a.plot([0, 0.11], [0, 0.11], "-", color="0.5", lw=0.8)
    a.set_xlim(0, 0.11)
    a.set_ylim(0, 0.11)
    a.set_xlabel(r"$\Delta\lambda$ over %g m, one pass" % SPAN)
    a.set_ylabel(r"the same ice, another pass")
    ms.tag(a, "B")

    a = ax[0, 2]
    sep = np.array([(lo + hi) / 2 for lo, hi in SEP_BINS])
    a.semilogy(sep, vg["all"], "o-", color="k", lw=2.0, ms=6)
    a.semilogy(sep, vg["ew"], "-", color=C_EW, lw=1.4)
    a.semilogy(sep, vg["ns"], "-", color=C_NS, lw=1.4)
    a.set_xlabel("separation [km]")
    a.set_ylabel(r"rms $\Delta\lambda$ difference $/\sqrt{2}$, 500-1450 m")
    a.set_ylim(3e-4, 3e-2)
    ms.tag(a, "C")

    a = ax[1, 0]
    a.plot(deep, loo["pred"], "o", color="k", ms=4.5, alpha=0.8)
    lim = [0.045, 0.108]
    a.plot(lim, lim, "-", color="0.5", lw=0.8)
    a.set_xlim(lim)
    a.set_ylim(lim)
    a.set_xlabel(r"$\Delta\lambda$ measured, 500-1450 m")
    a.set_ylabel("predicted from the other frames")
    ms.tag(a, "D")

    a = ax[1, 1]
    a.plot(deep, scale, "o", color="k", ms=4.5, alpha=0.8)
    a.axhline(ns.NODE_FC / ns.FC, color=C_NS, lw=1.4)
    lo = np.nanpercentile(deep, 12.5)
    xs = np.linspace(np.nanmin(deep), np.nanmax(deep), 20)
    a.plot(xs, np.nanmedian(scale) * lo / xs, "--", color="0.5", lw=1.4)
    a.set_ylim(0.55, 1.1)
    a.set_xlabel(r"$\Delta\lambda$ from the phase, 500-1450 m")
    a.set_ylabel("power fringe rate / phase delay rate")
    ms.tag(a, "E")

    a = ax[1, 2]
    a.plot(budget["shape"], ZC, "o-", color="k", lw=2.0, ms=6)
    a.plot(budget["profile"] / budget["profile"].mean(), ZC, "-",
           color="0.5", lw=1.4)
    a.set_ylim(1500, 0)
    a.set_xlim(0, 2.1)
    a.set_xlabel("depth pattern, relative to its mean")
    a.set_ylabel("depth [m]")
    ms.tag(a, "F")
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


# -------------------------------------------------------------- main
def main():
    ms.apply()
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.expanduser(os.environ.get("SCAR_DATA", "~/data/opr/scar")),
        "node")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figs")
    os.makedirs(out, exist_ok=True)
    rows = load(src, sys.argv[3] if len(sys.argv) > 3
                else "node_*_B2000.mat")
    tags = np.array([r["tag"] for r in rows])
    days = np.array([t[:11] for t in tags])
    x, y = local_xy(rows)
    zm = rows[0]["zm"]
    print("%d blocks from %d frames" % (len(rows), np.unique(tags).size))

    # ---- depth
    outs = [ns.resolve(r) for r in rows]
    plain = [ns.resolve(r, ns.POSTER, node_fc=ns.FC) for r in rows]
    sel = (zm >= BAND[0]) & (zm <= BAND[1])
    med = lambda k, oo: np.median([np.nanmedian(o[k][sel])  # noqa: E731
                                   for o in oo])
    print("DEPTH")
    print("  numbers per column: plain model %.1f, depth model %.1f (of %d "
          "cells, %d data rows)" % (
              np.median([o["n_eff"] for o in plain]),
              np.median([o["n_eff"] for o in outs]), zm.size,
              np.median([r["inc_d"].size + r["node_z"].size + 1
                         for r in rows])))
    print("  one value averages over: plain model %.0f m, depth model %.0f m;"
          " posterior sd %.4f and %.4f" % (
              med("resolution_m", plain), med("resolution_m", outs),
              med("sigma", plain), med("sigma", outs)))
    psi = np.array([o["psi_sd_ns"] for o in outs])
    print("  window phase by block: median %.4f ns (p10 %.4f, p90 %.4f) = "
          "%.1f deg of carrier" % (np.median(psi), np.percentile(psi, 10),
                                   np.percentile(psi, 90),
                                   np.median(psi) * ns.FC * 1e-9 * 360))
    ev = survey_evidence(rows)
    print("  log evidence against the plain model: nodes at %.0f MHz %+.0f;"
          " depth model %+.0f; best node frequency %.0f MHz" % (
              ns.NODE_FC / 1e6, ev["poster_band"], ev["smooth"],
              ev["freqs"][np.argmax(ev["curve"])]))
    md = model_dependence(rows, outs, plain, zm)
    for n, lab in ((0, "10 m"), (5, "150 m")):
        print("  at %s averaging: plain minus depth model %.4f, half / "
              "twice the correlation length %.4f / %.4f, spread between "
              "blocks %.4f" % (lab, md["plain"][n], md["half"][n],
                               md["twice"][n], md["lateral"][n]))
    im = interval_matrix(rows, outs)
    mi = np.array([[np.mean(o["m"][(zm >= a) & (zm < b)])
                    for a, b in zip(EDGES[:-1], EDGES[1:])] for o in outs])
    print("  depth model against the prior-free %g m means: rms %.4f, "
          "mean %+.4f (n %d)" % (SPAN, rms(mi - im), np.nanmean(mi - im),
                                 np.isfinite(im).sum()))
    rs = residual_statistics(rows, outs)
    acf1 = np.array([o["acf1"] for o in outs])
    print("STRIPES")
    print("  increment residuals: rms %.4f ns; autocorrelation at 25, 50, "
          "75 m %+.2f %+.2f %+.2f; by block median %+.2f, positive in %d "
          "of %d" % (rms(rs["res"]), rs["acf"][0], rs["acf"][1],
                     rs["acf"][2], np.nanmedian(acf1), np.sum(acf1 > 0),
                     acf1.size))
    print("  residual delay over 25 / 100 / 300 m: %.4f / %.4f / %.4f ns "
          "(a delay would reach %.4f / %.4f)" % (
              rs["sf"][0], rs["sf"][3], rs["sf"][6],
              rs["sf"][0] * 2, rs["sf"][0] * np.sqrt(12)))
    k = rows[0]["k"]
    print("  so a mean over 100 / 150 / 300 m is good to %.4f / %.4f / "
          "%.4f" % tuple(rs["sf"][i] / (k * rs["spans"][i])
                         for i in (3, 4, 6)))

    # ---- lateral
    print("LATERAL")
    deep = np.nanmean(im[:, DEEP], axis=1)
    theta = np.array([r["theta"] for r in rows])
    pi_, pj_ = same_ice_pairs(x, y, tags)
    cross = np.array([ne.is_ns(rows[a]["heading"])
                      != ne.is_ns(rows[b]["heading"])
                      for a, b in zip(pi_, pj_)])
    dth = (theta[pi_] - theta[pj_] + 90) % 180 - 90
    print("  same ice, another pass (%d pairs, %d crossing, %d on another "
          "day): per block %.4f over %g m, %.4f over 500-1450 m, axis %.1f "
          "deg" % (pi_.size, cross.sum(),
                   np.sum(days[pi_] != days[pj_]),
                   rms(im[pi_] - im[pj_]) / np.sqrt(2), SPAN,
                   rms(deep[pi_] - deep[pj_]) / np.sqrt(2),
                   rms(dth) / np.sqrt(2)))
    noise = np.sqrt(np.nanmean((im[pi_] - im[pj_]) ** 2, axis=0) / 2)
    print("  spread between blocks over that error, by interval: "
          + " ".join("%.0f" % v for v in np.nanstd(im, axis=0) / noise))
    vg = dict(all=variogram(deep, x, y, tags)[0])
    for key, along in (("ew", lambda dx, dy: dx > 2 * dy),
                       ("ns", lambda dx, dy: dy > 2 * dx)):
        v, n = variogram(deep, x, y, tags, along)
        vg[key] = np.where(n >= 20, v, np.nan)
    print("  rms difference/sqrt(2), 500-1450 m, by separation: "
          + " ".join("%g-%g km %.4f" % (lo, hi, v)
                     for (lo, hi), v in zip(SEP_BINS, vg["all"])))
    pred, e, e0 = leave_group_out(deep, x, y, tags)
    e_day = leave_group_out(deep, x, y, days)[1]
    loo = dict(pred=pred, rms=e, rms0=e0, rms_day=e_day)
    print("  each frame predicted from the others: rms %.4f (survey mean "
          "%.4f, skill %.2f); whole days left out %.4f" % (
              e, e0, 1 - e ** 2 / e0 ** 2, e_day))
    _, e, e0 = leave_group_out(theta, x, y, tags, circular=True)
    print("  slow axis likewise: rms %.1f deg (survey mean %.1f)" % (e, e0))
    sc = np.array([power_phase_scale(r) for r in rows])
    lo, hi = deep < np.nanpercentile(deep, 25), deep > np.nanpercentile(
        deep, 75)
    print("  power fringe rate over phase delay rate: median %.3f (sd %.3f;"
          " band centre / LO = %.3f); weakest quarter %.3f at dlam %.3f, "
          "strongest %.3f at %.3f, where one fabric for all would give "
          "%.2f" % (np.median(sc[:, 0]), sc[:, 0].std(), ns.NODE_FC / ns.FC,
                    sc[lo, 0].mean(), deep[lo].mean(), sc[hi, 0].mean(),
                    deep[hi].mean(),
                    sc[lo, 0].mean() * deep[lo].mean() / deep[hi].mean()))
    budget = variance_budget(im)
    print("  variance of the %g m means (%d blocks with every interval): "
          "%.0f %% depth, %.0f %% between blocks, of which one map %.0f %%"
          % (SPAN, budget["n"], 100 * budget["depth"],
             100 * budget["lateral"], 100 * budget["modes"][0]))

    # ---- the horizon
    hz, hs, hg = horizon(rows, outs)
    print("HORIZON")
    print("  a 30 m ramp at 380-600 m: log evidence %+.0f over the survey, "
          "%+.1f in the median block, positive in %d of %d" % (
              hg.sum(), np.median(hg), np.sum(hg > 0), hg.size))
    print("  depth median %.0f m (p10 %.0f, p90 %.0f); step %+.3f (p10 "
          "%+.3f, p90 %+.3f), positive in %d" % (
              np.median(hz), np.percentile(hz, 10), np.percentile(hz, 90),
              np.median(hs), np.percentile(hs, 10), np.percentile(hs, 90),
              np.sum(hs > 0)))
    _, e, e0 = leave_group_out(hz, x, y, tags)
    print("  same ice, another pass: depth to %.0f m per block; predicted "
          "from the other frames to %.0f m (survey spread %.0f m)" % (
              rms(hz[pi_] - hz[pj_]) / np.sqrt(2), e, e0))

    # ---- products and figures
    pick = int(np.argmin(np.hypot(x, y)))
    print("node_depth_model.png panel A is %s block %d, the block nearest "
          "the survey's centre" % (rows[pick]["tag"], rows[pick]["block"]))
    fig_depth(rows, outs, plain, ev, rs, md, pick,
              os.path.join(out, "node_depth_model.png"))
    fig_lateral(rows, x, y, deep, im, (pi_, pj_), vg, loo, sc[:, 0], budget,
                os.path.join(out, "node_lateral_field.png"))
    resolved = np.array([np.where(np.isfinite(r["dlam"]), o["m"], np.nan)
                         for r, o in zip(rows, outs)])
    np.savez(os.path.join(out, "node_resolved.npz"), tags=tags,
             block=np.array([r["block"] for r in rows]),
             lat=np.array([r["lat"] for r in rows]),
             lon=np.array([r["lon"] for r in rows]), x_km=x, y_km=y,
             heading=np.array([r["heading"] for r in rows]), theta=theta,
             zm=zm, dlam=resolved,
             sigma=np.array([o["sigma"] for o in outs]),
             resolution_m=np.array([o["resolution_m"] for o in outs]),
             n_eff=np.array([o["n_eff"] for o in outs]), psi_sd_ns=psi,
             acf1=acf1, interval_edges=EDGES, interval_dlam=im,
             loo_pred=pred, power_phase_scale=sc[:, 0], horizon_z=hz,
             horizon_step=hs, horizon_logev_gain=hg)
    print("wrote node_depth_model.png, node_lateral_field.png, "
          "node_resolved.npz in", out)


if __name__ == "__main__":
    main()

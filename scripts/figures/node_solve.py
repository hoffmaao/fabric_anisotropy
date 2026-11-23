"""nodeFabric's linear solve, from the measurements a node product stores.

Everything the depth solve needs is in the product: the wrapped phase
increments between windows, the accepted nodes and the firn reference.
So the posterior can be recomputed here under any depth model without
the radar data, which is how a model is compared over a whole survey
before a single frame is rerun. The functions mirror the helpers of
+ptt/nodeFabric.m (H_kernel, H_cd, H_solve,
H_psi_auto, H_intervals) and reproduce its numbers to rounding: the depth
model gives back the dlam a product stores to 1e-10.

    import node_solve as ns
    out = ns.resolve(row)                     # the depth model, as stored
    out = ns.resolve(row, dict(ns.SMOOTH, len_m=250.0))   # a finer one
    out = ns.resolve(row, ns.POSTER, node_fc=ns.FC)   # the plain model:
    #   every phase increment taken as a delay, no window phase

`row` is one block as node_export.read_product returns it.
"""
import numpy as np

FC = 750e6                    # the LO: what the carrier phase counts in
NODE_FC = 735e6               # the band's centre on Ridge A: what a node
#                               counts in (measured 734-738 MHz on
#                               20250108_02_003; the survey's evidence
#                               peaks at 730-735). Products that carry
#                               nodes.fc_hz use their own.
MODEL_DZ = 10.0
SMOOTH = dict(kernel="matern32", sd=0.04, len_m=500.0, mean=0.0,
              level_sd=0.2, psi_sd_ns="auto", psi_len_m=25.0)
POSTER = dict(kernel="exp", sd=0.05, len_m=150.0, mean=0.0,
              level_sd=0.0, psi_sd_ns=0.0, psi_len_m=25.0)
PSI_GRID = np.logspace(np.log10(0.002), np.log10(0.04), 14)


def a_abs(zz, zm, k, dz=MODEL_DZ):
    """k times the thickness of each model cell lying above each depth"""
    zz = np.atleast_1d(np.asarray(zz, float))
    return k * np.clip(zz[:, None] - zm[None, :], 0.0, dz)


def kernel(zm, model):
    """prior covariance of dlam on the model grid"""
    d = np.abs(zm[:, None] - zm[None, :])
    sd, length = model["sd"], model["len_m"]
    if model["kernel"] == "exp":
        c = sd ** 2 * np.exp(-d / length)
    elif model["kernel"] == "matern32":
        a = np.sqrt(3.0) * d / length
        c = sd ** 2 * (1 + a) * np.exp(-a)
    else:
        raise ValueError("unknown kernel %r" % model["kernel"])
    return c + 1e-12 * sd ** 2 * np.eye(zm.size) + model["level_sd"] ** 2


def increment_matrix(iw, n_windows):
    """D with psi(lower window) - psi(upper window) on each increment"""
    d = np.zeros((iw.size, n_windows))
    d[np.arange(iw.size), iw] = -1.0
    d[np.arange(iw.size), iw + 1] = 1.0
    return d


def data_cov(sd_abs, sd_inc, iw, zw, psi_sd, psi_len):
    """Covariance of [absolute rows; increments]: the absolute rows (firn
    anchor, nodes) independent, the increments their own noise plus the
    window phase differenced."""
    ci = np.diag(np.asarray(sd_inc, float) ** 2)
    if psi_sd > 0 and iw.size:
        d = increment_matrix(iw, zw.size)
        rp = np.exp(-np.abs(zw[:, None] - zw[None, :]) / max(psi_len, 1e-9))
        ci = ci + psi_sd ** 2 * d @ rp @ d.T
    na = np.size(sd_abs)
    cd = np.zeros((na + iw.size, na + iw.size))
    cd[:na, :na] = np.diag(np.asarray(sd_abs, float) ** 2)
    cd[na:, na:] = ci
    return cd


def log_evidence(g, d, cd, cm, m_prior):
    s = g @ cm @ g.T + cd
    lc = np.linalg.cholesky((s + s.T) / 2)
    y = np.linalg.solve(lc, d - g @ m_prior)
    return (-0.5 * y @ y - np.log(np.diag(lc)).sum()
            - 0.5 * d.size * np.log(2 * np.pi))


def solve(g, d, cd, cm, m_prior):
    """The linear Gaussian posterior in data space. Returns a dict: m,
    cov, R (the resolution matrix), res, chi2_dof, logev, alpha."""
    s = g @ cm @ g.T + cd
    lc = np.linalg.cholesky((s + s.T) / 2)
    r0 = d - g @ m_prior

    def si(x):
        return np.linalg.solve(lc.T, np.linalg.solve(lc, x))
    alpha = si(r0)
    k = cm @ g.T
    m = m_prior + k @ alpha
    rm = k @ si(g)
    cov = cm - k @ si(k.T)
    res = d - g @ m
    chi2 = res @ np.linalg.solve(cd, res)
    dof = max(1.0, d.size - np.trace(rm))
    logev = (-0.5 * r0 @ alpha - np.log(np.diag(lc)).sum()
             - 0.5 * d.size * np.log(2 * np.pi))
    return dict(m=m, cov=cov, R=rm, res=res, chi2_dof=chi2 / dof,
                logev=logev, alpha=alpha)


def psi_auto(g, d, sd_abs, sd_inc, iw, zw, psi_len, cm, m_prior):
    """The window phase level that maximises the block's evidence: a log
    grid, a parabola through the top three, zero when none beats none."""
    le = np.array([log_evidence(
        g, d, data_cov(sd_abs, sd_inc, iw, zw, p, psi_len), cm, m_prior)
        for p in PSI_GRID])
    le0 = log_evidence(
        g, d, data_cov(sd_abs, sd_inc, iw, zw, 0.0, psi_len), cm, m_prior)
    i = int(np.argmax(le))
    if iw.size < 5 or le0 >= le[i]:
        return 0.0
    if 0 < i < PSI_GRID.size - 1:
        x = np.log(PSI_GRID[i - 1:i + 2])
        p = np.polyfit(x - x[1], le[i - 1:i + 2], 2)
        if p[0] < 0:
            return float(np.exp(x[1] - p[1] / (2 * p[0])))
    return float(PSI_GRID[i])


def window_index(row):
    """index of each increment's upper window in the window grid"""
    zw = row["zw"]
    return np.array([int(np.argmin(np.abs(zw - z))) for z in row["inc_z1"]],
                    dtype=int)


def firn_reference_depth(row, ref_z=(30.0, 90.0)):
    zw = row["zw"]
    inref = (zw >= ref_z[0]) & (zw <= ref_z[1]) & np.isfinite(row["fine"])
    if not inref.any():
        inref = np.arange(zw.size) < 2
    return zw[inref].mean()


def node_delays(row, node_fc=None):
    """n half waves of the frequency the nodes count in, less the firn
    phase. node_fc None: the product's own nodes.fc_hz if it has them,
    else NODE_FC."""
    n = row["node_n"]
    if node_fc is None:
        fc_n = row.get("node_fc")
        fc_n = np.full(n.size, NODE_FC) if fc_n is None else fc_n
    else:
        fc_n = np.full(n.size, float(node_fc))
    return n / (2 * fc_n * 1e-9) - row["fine_ref"] * FC / fc_n


def block_rows(row, node_fc=None, use_nodes=True):
    """G, d and the row sigmas of one block: anchor, nodes, increments"""
    zm, k = row["zm"], row["k"]
    g = [a_abs(firn_reference_depth(row), zm, k)]
    d = [np.zeros(1)]
    sd_abs = [np.array([0.02])]
    if use_nodes:
        ok = np.isfinite(row["node_n"])
        g.append(a_abs(row["node_z"][ok], zm, k))
        d.append(node_delays(row, node_fc)[ok])
        sd_abs.append(row["node_sig"][ok])
    g.append(a_abs(row["inc_z2"], zm, k) - a_abs(row["inc_z1"], zm, k))
    d.append(row["inc_d"])
    return np.vstack(g), np.concatenate(d), np.concatenate(sd_abs)


def resolution_width(rm, zm):
    """width [m] of the central 68 % of each cell's averaging kernel"""
    out = np.full(zm.size, np.nan)
    for i in range(zm.size):
        a = np.abs(rm[i])
        if a.sum() > 0:
            c = np.cumsum(a) / a.sum()
            out[i] = np.interp(0.84, c, zm) - np.interp(0.16, c, zm)
    return out


def intervals(row, psi_sd, psi_len=25.0, span=125.0, edges=None):
    """dlam over fixed spans from the summed increments alone: no prior,
    no nodes. NaN where increments cover less than 60 % of the span."""
    zw, k = row["zw"], row["k"]
    z1, z2, inc = row["inc_z1"], row["inc_z2"], row["inc_d"]
    if edges is None:
        edges = np.arange(span * np.ceil(zw[0] / span), zw[-1] + 1e-9, span)
    edges = np.asarray(edges, float)
    zi = (z1 + z2) / 2
    iw = window_index(row)
    ci = data_cov(np.zeros(0), row["inc_s"], iw, zw, psi_sd, psi_len)
    n = edges.size - 1
    dlam = np.full(n, np.nan)
    sig = np.full(n, np.nan)
    cover = np.zeros(n)
    for j in range(n):
        m = (zi >= edges[j]) & (zi < edges[j + 1])
        th = (z2[m] - z1[m]).sum()
        cover[j] = th / (edges[j + 1] - edges[j])
        if cover[j] >= 0.6:
            dlam[j] = inc[m].sum() / (k * th)
            sig[j] = np.sqrt(ci[np.ix_(m, m)].sum()) / (k * th)
    return dict(z1=edges[:-1], z2=edges[1:], dlam=dlam, sigma=sig,
                coverage=cover)


def resolve(row, model=None, node_fc=None, use_nodes=True):
    """One block's posterior under a depth model. Returns a dict: zm, m,
    sigma, R, resolution_m, n_eff, psi_sd_ns, psi (window phase), res
    (increment residuals), acf1, chi2_dof, logev, node_res, model."""
    model = dict(SMOOTH if model is None else model)
    zm = row["zm"]
    g, d, sd_abs = block_rows(row, node_fc, use_nodes)
    iw = window_index(row)
    zw = row["zw"]
    cm = kernel(zm, model)
    mp = np.full(zm.size, float(model["mean"]))
    if isinstance(model["psi_sd_ns"], str):
        model["psi_sd_ns"] = psi_auto(
            g, d, sd_abs, row["inc_s"], iw, zw, model["psi_len_m"], cm, mp)
    psi_sd = model["psi_sd_ns"]
    cd = data_cov(sd_abs, row["inc_s"], iw, zw, psi_sd, model["psi_len_m"])
    o = solve(g, d, cd, cm, mp)
    na = sd_abs.size
    r_inc = o["res"][na:]
    psi = np.full(zw.size, np.nan)
    if psi_sd > 0 and iw.size:
        dmat = increment_matrix(iw, zw.size)
        rp = np.exp(-np.abs(zw[:, None] - zw[None, :]) / model["psi_len_m"])
        p_hat = psi_sd ** 2 * rp @ (dmat.T @ o["alpha"][na:])
        used = np.zeros(zw.size, bool)
        used[iw] = True
        used[iw + 1] = True
        psi[used] = p_hat[used]
    nb = np.where(np.diff(iw) == 1)[0]
    acf1 = np.nan
    if nb.size >= 5 and np.any(r_inc):
        acf1 = np.mean(r_inc[nb] * r_inc[nb + 1]) / np.mean(r_inc ** 2)
    return dict(zm=zm, m=o["m"], sigma=np.sqrt(np.clip(np.diag(o["cov"]),
                                                       0, None)),
                R=o["R"], resolution_m=resolution_width(o["R"], zm),
                n_eff=float(np.trace(o["R"])), psi_sd_ns=psi_sd, psi=psi,
                res=r_inc, acf1=acf1, chi2_dof=o["chi2_dof"],
                logev=o["logev"], node_res=o["res"][1:na], model=model)

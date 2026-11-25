"""Table and figures from the node-method products (run_node_frame,
stages/node).

    python scripts/figures/node_export.py \\
        [<mirror of stages/node>] [out_dir] [pattern]

Reads every node_<tag>_B2000.mat (MATLAB v7.3, read with h5py; the
directory defaults to $SCAR_DATA/node, $SCAR_DATA to ~/data/opr/scar),
writes figs/node_blocks.csv (one row per block: tag, block, lat, lon,
heading, slow axis geographic, plateau ratio, gains, node count,
chi2/dof, abstention depth, strength by 200 m band) and
figs/node_profiles.npz (dlam(z) and sigma per block), prints the survey's
numbers, and draws three figures. Nothing is written on the panels but
their letters; these are the captions. Blocks driven north-south are
blue and blocks driven east-west orange throughout.

  node_heading_test.png   A: slow axis of each block against its heading
                          (marker size: the plateau ratio, the
                          orientation's quality; dashed: each family's
                          circular mean). A fabric axis does not turn
                          with the antennas. B: dlam(z) of each block
                          (thin) and each family's median (thick); black
                          stairs: dlam from the spacing of the nodes
                          alone.
  node_ladders.png        node depths of every block, in file order
  node_crossings.png      N-S against E-W blocks within 1.5 km of each
                          other, the same ice driven in two directions.
                          A: depths of the nodes both blocks have. B:
                          strength by 200 m band (colour: depth, dark
                          shallow). C: slow axis. None of the three can
                          depend on heading.

read_product also returns what nodeFabric's linear solve was given (the
phase increments, the node delays, the firn reference), which is what
node_solve.py and node_representation.py work from.
"""
import glob
import os
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)
import morrow_style as ms  # noqa: E402
C_NS, C_EW = "#2a6fb0", "#d2691e"
# ns of two-way delay per m per unit dlam
K_NS = 0.063688
LO_HZ = 750e6
BANDS = [(200, 400), (400, 600), (600, 900), (900, 1200), (1200, 1450)]


def h5v(f, ref):
    """dereference an h5py object reference and return a numpy array"""
    d = f[ref]
    v = d[()]
    if isinstance(v, np.ndarray) and v.dtype.names and "real" in v.dtype.names:
        return v["real"] + 1j * v["imag"]
    return np.asarray(v)


def field(f, R, name, b):
    """field `name` of block b of R, whether MATLAB saved R as a scalar
    struct (fields are datasets or groups) or a struct array (fields are
    object reference arrays). Returns an h5py Group for struct-valued
    fields, else a numpy array."""
    d = R[name]
    if isinstance(d, h5py.Group):
        return d
    v = d[()]
    if v.dtype.kind == "O":
        obj = f[v.ravel()[b]]
        if isinstance(obj, h5py.Group):
            return obj
        return h5v(f, v.ravel()[b])
    if v.dtype.names and "real" in v.dtype.names:
        return v["real"] + 1j * v["imag"]
    a = np.asarray(v)
    return a if a.size == 1 or a.ndim == 0 else a


def chars(a):
    """a MATLAB char array (uint16) as str; anything else (e.g. the uint64
    [0 0] MATLAB writes for an empty value) is the empty string"""
    a = np.asarray(a)
    if a.dtype != np.uint16:
        return ""
    return "".join(chr(int(c)) for c in a.ravel())


def val(ds):
    """a dataset's values, honouring MATLAB's empty marker"""
    if ds.attrs.get("MATLAB_empty", 0) == 1:
        return np.array([])
    return np.asarray(ds[()]).ravel()


def scalar(f, g, name):
    v = g[name][()]
    if isinstance(v, np.ndarray) and v.dtype.kind == "O":
        v = h5v(f, v.ravel()[0])
    if isinstance(v, np.ndarray) and v.dtype.names and "real" in v.dtype.names:
        return complex(v["real"].ravel()[0], v["imag"].ravel()[0])
    return np.asarray(v).ravel()[0]


def read_product(fn):
    """One dict per block of a product. A block the frame runner skipped,
    or that failed, has ok False and the reason in err."""
    rows = []
    with h5py.File(fn) as f:
        tag = chars(f["tag"])
        R = f["R"]
        meta = f["meta"]
        lat = np.asarray(meta["lat"]).ravel()
        lon = np.asarray(meta["lon"]).ravel()
        az = np.asarray(meta["az"]).ravel()
        x = np.asarray(meta["x"]).ravel()
        zfull = np.asarray(f["z"]).ravel()
        opts = f["opts"] if "opts" in f else {}
        step_m = (float(np.asarray(opts["step_m"]).ravel()[0])
                  if "step_m" in opts else None)
        blk = R["block"]
        n = 1 if isinstance(blk, h5py.Group) else np.asarray(blk[()]).size
        for b in range(n):
            where = dict(tag=tag, block=b + 1, lat=lat[b], lon=lon[b],
                         x=x[b])

            def g(name, b=b):
                return field(f, R, name, b)
            skipped = bool(np.asarray(g("skipped")).ravel()[0])
            msg = chars(g("error"))
            if skipped or msg:
                rows.append(dict(where, heading=az[b], ok=False,
                                 err="skipped" if skipped else msg))
                continue
            try:
                rows.append(dict(where, step_m=step_m,
                                 **read_block(f, g, zfull)))
            except Exception as e:  # a block that did not write every field
                rows.append(dict(where, heading=az[b], ok=False,
                                 err=repr(e)))
    return rows


def read_block(f, g, zfull):
    """the fields of one block the estimator ran on"""
    def one(name):
        return float(np.asarray(g(name)).ravel()[0])

    def arr(name):
        return np.asarray(g(name)).ravel()
    zm, dlam, sig = arr("zm"), arr("dlam"), arr("sigma_dlam")
    theta = one("theta_geo_deg")
    plateau = float(scalar(f, g("theta_scan"), "plateau_ratio"))
    gains = g("gains")
    gr = complex(scalar(f, gains, "gr"))
    gt = complex(scalar(f, gains, "gt"))
    nodes = g("nodes")
    nz, nn = val(nodes["z"]), val(nodes["n"])
    nsrc = chars(nodes["n_source"][()]) if "n_source" in nodes else ""
    # the smoothed 45-deg cross-pol ratio the nodes were picked on, and the
    # candidate extrema the phase test rejected
    smooth = val(nodes["smoothed"]) if "smoothed" in nodes else np.array([])
    rej_z = np.array([])
    rej = nodes["rejected"] if "rejected" in nodes else None
    if isinstance(rej, h5py.Group) and "z" in rej:
        rej_z = val(rej["z"])
    chi2 = one("chi2_dof")
    abst = one("abstain_below_m")
    # a block the estimator ran on but could not read: no gains (the gain
    # estimator saw no valid rows) or no node accepted; its dlam is the
    # prior and its axis is undefined
    if not (np.isfinite(gr) and np.isfinite(gt)):
        err = "no gains"
    elif nz.size == 0:
        err = "no nodes"
    else:
        err = ""
    row = dict(
        heading=one("heading_deg"), ok=not err, err=err, theta=theta,
        plateau=plateau, gr_abs=abs(gr), gr_deg=np.degrees(np.angle(gr)),
        gt_abs=abs(gt), gt_deg=np.degrees(np.angle(gt)), n_nodes=nz.size,
        n_source=nsrc, chi2=chi2, abstain=abst, zm=zm, dlam=dlam, sigma=sig,
        node_z=nz, node_n=nn, node_resid=val(nodes["resid_ns"]),
        node_type=val(nodes["type"]), node_prom=val(nodes["prom_db"]),
        dlam_raw=arr("dlam_raw"), z=zfull, ratio_smoothed=smooth,
        rejected_z=rej_z)
    for lo, hi in BANDS:
        m = (zm >= lo) & (zm < hi)
        row["dlam_%d_%d" % (lo, hi)] = (np.nanmean(dlam[m]) if m.any()
                                        else np.nan)
    row.update(measurements(g, nodes, zfull))
    return row


def measurements(g, nodes, zfull):
    """What nodeFabric's linear solve was given, so node_solve can redo it
    under another depth model without the radar data: the eigen-pair
    windows, the wrapped phase increments between them, the node delays
    and the firn reference. node_fc is None for products made before the
    nodes carried their own frequency. ratio45 and ratio0 are the cross/co
    power ratios [dB] 45 deg off the slow axis and in its frame."""
    w = g("windows")
    inc = g("inc")

    def one(name):
        return float(np.asarray(g(name)).ravel()[0])
    return dict(
        k=one("k_ns_per_m"), fine_ref=one("fine_ref"),
        zw=zfull[val(w["idx"]).astype(int) - 1], fine=val(w["fine"]),
        gamma=val(w["gamma"]), inc_z1=val(inc["z1"]), inc_z2=val(inc["z2"]),
        inc_d=val(inc["dtau"]), inc_s=val(inc["sigma"]),
        node_tau=val(nodes["tau_ns"]), node_sig=val(nodes["sigma_ns"]),
        node_fc=val(nodes["fc_hz"]) if "fc_hz" in nodes else None,
        ratio45=np.asarray(g("ratio45")).ravel(),
        ratio0=np.asarray(g("ratio0")).ravel())


def is_ns(h):
    return abs(np.sin(np.radians(h))) < 0.5


def circular_mean(theta_deg):
    """mean of axes (mod 180) in degrees"""
    ph = np.exp(2j * np.radians(theta_deg))
    return np.degrees(np.angle(np.mean(ph))) / 2 % 180


def crossing_pairs(good, max_km=1.5):
    """(N-S block, E-W block) pairs whose centres lie within max_km: the
    same ice seen along two headings. Node depths and strength there cannot
    depend on heading."""
    R_E = 6371.0
    out = []
    ns = [r for r in good if is_ns(r["heading"])]
    ew = [r for r in good if not is_ns(r["heading"])]
    for a in ns:
        best = None
        for b in ew:
            dlat = np.radians(b["lat"] - a["lat"])
            dlon = np.radians(b["lon"] - a["lon"])
            d = R_E * np.sqrt(dlat ** 2
                              + (np.cos(np.radians(a["lat"])) * dlon) ** 2)
            if d <= max_km and (best is None or d < best[0]):
                best = (d, b)
        if best:
            out.append((a, best[1], best[0]))
    return out


def shared_nodes(a, b):
    """node depths of the counts both blocks have: {n: (z_a, z_b)}"""
    na = dict(zip(a["node_n"].astype(int), a["node_z"]))
    nb = dict(zip(b["node_n"].astype(int), b["node_z"]))
    return {n: (na[n], nb[n]) for n in sorted(set(na) & set(nb))}


def fig_crossings(good, out):
    pairs = crossing_pairs(good)
    if not pairs:
        print("no crossing pairs within 1.5 km")
        return
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2), constrained_layout=True)
    dz_all = []
    for a, b, d in pairs:
        for za, zb in shared_nodes(a, b).values():
            ax[0].plot(za, zb, "o", color="k", ms=4, alpha=0.7)
            dz_all.append(zb - za)
        for lo, hi in BANDS:
            k = "dlam_%d_%d" % (lo, hi)
            ax[1].plot(a[k], b[k], "o", ms=5, alpha=0.8,
                       color=plt.cm.viridis((lo - 200) / 1250))
        ax[2].plot([a["theta"]], [b["theta"]], "o", color="k", ms=5,
                   alpha=0.7)
    dz_med = np.median(dz_all) if dz_all else np.nan
    dz_rms = np.sqrt(np.mean(np.square(dz_all))) if dz_all else np.nan
    ax[0].plot([300, 1500], [300, 1500], "-", color="0.6", lw=0.8)
    ax[0].set_xlabel("node depth, N-S block [m]")
    ax[0].set_ylabel("node depth, E-W block [m]")
    ms.tag(ax[0], "A")
    ax[1].plot([0, 0.12], [0, 0.12], "-", color="0.6", lw=0.8)
    ax[1].set_xlabel(r"$\Delta\lambda$, N-S block")
    ax[1].set_ylabel(r"$\Delta\lambda$, E-W block")
    ms.tag(ax[1], "B")
    ax[2].plot([60, 160], [60, 160], "-", color="0.6", lw=0.8)
    ax[2].set_xlabel(r"slow axis, N-S block [$^\circ$]")
    ax[2].set_ylabel(r"slow axis, E-W block [$^\circ$]")
    ms.tag(ax[2], "C")
    ax[2].set_xlim(60, 160)
    ax[2].set_ylim(60, 160)
    fig.savefig(os.path.join(out, "node_crossings.png"), dpi=170)
    plt.close(fig)
    da = np.array([b["theta"] - a["theta"] for a, b, d in pairs])
    da = (da + 90) % 180 - 90
    print("crossing pairs: %d; slow-axis E-W minus N-S median %+.1f deg "
          "(rms %.1f); node depth E-W minus N-S median %+.0f m (rms %.0f)"
          % (len(pairs), np.median(da), np.sqrt(np.mean(da ** 2)), dz_med,
             dz_rms))
    for a, b, d in pairs:
        print("  %s b%d (heading %.0f, axis %.1f) with %s b%d (heading "
              "%.0f, axis %.1f): %.2f km apart, %d shared nodes" % (
                  a["tag"], a["block"], a["heading"], a["theta"], b["tag"],
                  b["block"], b["heading"], b["theta"], d,
                  len(shared_nodes(a, b))))
    for lo, hi in BANDS:
        k = "dlam_%d_%d" % (lo, hi)
        dd = np.array([b[k] - a[k] for a, b, d in pairs])
        print("  %4d-%4d m: dlam E-W minus N-S median %+.4f (n %d)" % (
            lo, hi, np.nanmedian(dd), np.isfinite(dd).sum()))


# the flat table's 25 m rows (centres)
ROW_Z = np.arange(35.0, 1500.0, 25.0)


def write_flow_table(good, fn, solver="half-wave node ladder (nodeFabric)",
                     prefix="node_flow_ct"):
    """The usable blocks as one flat table for gridding in 3D: one row
    per block and 25 m of depth, dlam the mean of the block's dlam(z)
    over that depth row (NaN where it abstains), theta the block's
    single slow axis in radians geographic,
    az the block heading in degrees. dlam_max = dlam: the node method has
    one axis per block, so the strength is already along it."""
    from scipy.io import savemat
    tags = sorted({r["tag"] for r in good})
    names = ("k", "blk", "lat", "lon", "az", "x", "zc", "dlam", "theta")
    cols = {k: [] for k in names}
    for r in good:
        zm = np.asarray(r["zm"], float)
        dl = np.asarray(r["dlam"], float)
        for zc in ROW_Z:
            m = np.abs(zm - zc) <= 12.5
            v = dl[m]
            v = (np.nanmean(v) if m.any() and np.isfinite(v).any()
                 else np.nan)
            row = (tags.index(r["tag"]) + 1, r["block"], r["lat"], r["lon"],
                   r["heading"], r["x"], zc, v, np.radians(r["theta"]))
            for key, value in zip(names, row):
                cols[key].append(value)
    T = {k: np.asarray(v, float) for k, v in cols.items()}
    T["dlam_max"] = T["dlam"]
    T["tags"] = np.array(tags, dtype=object)
    T["solver"] = solver
    T["prefix"] = prefix
    savemat(fn, T)


def half_wave_ns(row):
    """ns of delay between successive nodes of a block: half a wave of the
    frequency its nodes count in, the band's own centre as the product
    measured it (the LO for a product made before nodes carried one)"""
    fc = row.get("node_fc")
    if fc is None or not np.size(fc):
        return np.full(np.size(row["node_z"]), 0.5e9 / LO_HZ)
    return 0.5e9 / np.asarray(fc, float)


def node_spacing_profile(good):
    """Strength from the node depths alone, no solver: between successive
    nodes the two-way delay grows by one half wave, so the mean dlam over
    the interval is the half wave / (K_NS * spacing). Medians over blocks
    of each interval (surface to node 1, node 1 to 2, ...), returned as
    stair edges (len n + 1 depths) and values (the last repeated);
    intervals seen in fewer than 10 blocks are left out."""
    iv = {}
    for r in good:
        o = np.argsort(r["node_z"])
        z = np.r_[0.0, np.asarray(r["node_z"], float)[o]]
        n = np.r_[0, np.asarray(r["node_n"], int)[o]]
        half = half_wave_ns(r)[o]
        for a in range(1, len(z)):
            rate = ((n[a] - n[a - 1]) * half[a - 1]
                    / (K_NS * (z[a] - z[a - 1])))
            iv.setdefault(n[a], []).append((z[a - 1], z[a], rate))
    zs, ds, z_end = [], [], np.nan
    for n in sorted(iv):
        a = np.array(iv[n])
        if len(a) < 10:
            continue
        zs.append(np.median(a[:, 0]))
        ds.append(np.median(a[:, 2]))
        z_end = np.median(a[:, 1])
        print("  node spacing alone, interval to node %d (%4.0f-%4.0f m, "
              "n %d): dlam %.3f" % (n, np.median(a[:, 0]), np.median(a[:, 1]),
                                    len(a), np.median(a[:, 2])))
    zs.append(z_end)
    ds.append(ds[-1])
    return np.array(zs), np.array(ds)


def fig_heading_test(good, out):
    """slow axis against heading, and dlam(z) by family"""
    ns = np.array([is_ns(r["heading"]) for r in good])
    th = np.array([r["theta"] for r in good])
    hd = np.array([r["heading"] for r in good])
    pl = np.array([r["plateau"] for r in good])
    fig, ax = plt.subplots(1, 2, figsize=(12.8, 5.4),
                           constrained_layout=True)
    for fam, c in ((ns, C_NS), (~ns, C_EW)):
        ax[0].scatter(hd[fam], th[fam], s=22 + 6 * np.clip(pl[fam], 0, 8),
                      color=c, alpha=0.8, edgecolor="w", lw=0.5)
        if fam.any():
            ax[0].axhline(circular_mean(th[fam]), color=c, lw=0.8, ls="--")
    ax[0].set_xlabel(r"block heading [$^\circ$]")
    ax[0].set_ylabel(r"slow axis, geographic [$^\circ$]")
    ax[0].set_xlim(0, 180)
    ax[0].set_ylim(0, 180)
    ax[0].set_xticks([0, 45, 90, 135, 180])
    ax[0].set_yticks([0, 45, 90, 135, 180])
    ms.tag(ax[0], "A")
    zg = np.arange(50, 1500, 10.0)
    for fam, c in ((ns, C_NS), (~ns, C_EW)):
        prof = []
        for r in (r for r, f_ in zip(good, fam) if f_):
            ax[1].plot(r["dlam"], r["zm"], "-", color=c, lw=0.4, alpha=0.35)
            prof.append(np.interp(zg, r["zm"], r["dlam"], left=np.nan,
                                  right=np.nan))
        if prof:
            ax[1].plot(np.nanmedian(np.array(prof), axis=0), zg, "-",
                       color=c, lw=2.2)
    zs, ds = node_spacing_profile(good)
    # stairs in depth
    ax[1].plot(np.repeat(ds[:-1], 2), np.repeat(zs, 2)[1:-1], "-", color="k",
               lw=1.8)
    ax[1].set_ylim(1500, 0)
    ax[1].set_xlim(0, 0.12)
    ax[1].set_xlabel(r"$\Delta\lambda$ along the slow axis")
    ax[1].set_ylabel("depth [m]")
    ms.tag(ax[1], "B", x=0.93)
    fig.savefig(os.path.join(out, "node_heading_test.png"), dpi=170)
    plt.close(fig)


def print_families(good):
    """band medians and the mean slow axis of each heading family"""
    ns = np.array([is_ns(r["heading"]) for r in good])
    th = np.array([r["theta"] for r in good])
    print("band medians, N-S vs E-W:")
    for lo, hi in BANDS:
        k = "dlam_%d_%d" % (lo, hi)
        a = np.array([r[k] for r, f_ in zip(good, ns) if f_])
        b = np.array([r[k] for r, f_ in zip(good, ~ns) if f_])
        print("  %4d-%4d m: N-S %.4f (n %d)  E-W %.4f (n %d)" % (
            lo, hi, np.nanmedian(a), a.size, np.nanmedian(b), b.size))
    for fam, lab in ((ns, "N-S"), (~ns, "E-W")):
        if fam.any():
            ph = np.exp(2j * np.radians(th[fam]))
            sp = 0
            if fam.sum() > 1:
                sp = np.degrees(np.sqrt(-2 * np.log(np.abs(np.mean(ph))))
                                / 2)
            print("slow axis %s: %.1f deg (circular sd %.1f, n %d)" % (
                lab, circular_mean(th[fam]), sp, fam.sum()))


def fig_ladders(good, out):
    fig, ax = plt.subplots(figsize=ms.FRAME, constrained_layout=True)
    for i, r in enumerate(good):
        c = C_NS if is_ns(r["heading"]) else C_EW
        if r["node_z"].size:
            ax.plot(np.full(r["node_z"].size, i), r["node_z"], "o", color=c,
                    ms=3)
    ax.set_ylim(1500, 0)
    ax.set_xlabel("block (in file order)")
    ax.set_ylabel("node depth [m]")
    fig.savefig(os.path.join(out, "node_ladders.png"), dpi=170)
    plt.close(fig)


def main():
    ms.apply()
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.expanduser(os.environ.get("SCAR_DATA", "~/data/opr/scar")),
        "node")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figs")
    pattern = sys.argv[3] if len(sys.argv) > 3 else "node_*_B2000.mat"
    os.makedirs(out, exist_ok=True)
    rows = []
    for fn in sorted(glob.glob(os.path.join(src, pattern))):
        rows += read_product(fn)
    good = [r for r in rows if r.get("ok")]
    print("%d blocks read, %d usable, from %d frames" % (
        len(rows), len(good), len({r["tag"] for r in rows})))
    for r in rows:
        if not r.get("ok"):
            print("  not usable: %s block %d (heading %.1f): %s" % (
                r["tag"], r["block"], r["heading"], r.get("err", "")))
    # csv: every block, usable or not, with the reason
    keys = ["tag", "block", "ok", "err", "lat", "lon", "heading", "x",
            "theta", "plateau", "gr_abs", "gr_deg", "gt_abs", "gt_deg",
            "n_nodes", "n_source", "chi2", "abstain"]
    keys += ["dlam_%d_%d" % b for b in BANDS]
    with open(os.path.join(out, "node_blocks.csv"), "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
    np.savez(os.path.join(out, "node_profiles.npz"),
             tags=np.array([r["tag"] for r in good]),
             block=np.array([r["block"] for r in good]),
             heading=np.array([r["heading"] for r in good]),
             theta=np.array([r["theta"] for r in good]),
             zm=np.array([r["zm"] for r in good], dtype=object),
             dlam=np.array([r["dlam"] for r in good], dtype=object),
             sigma=np.array([r["sigma"] for r in good], dtype=object),
             allow_pickle=True)
    write_flow_table(good, os.path.join(out, "node_table.mat"))
    fig_heading_test(good, out)
    print_families(good)
    fig_ladders(good, out)
    fig_crossings(good, out)
    print("wrote node_blocks.csv, node_profiles.npz, node_heading_test.png, "
          "node_ladders.png, node_crossings.png in", out)


if __name__ == "__main__":
    main()

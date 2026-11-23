"""Node ladders along track: every Ridge A frame as a section, the accepted
half-wave nodes as dots on the data they were picked from.

    python scripts/figures/node_sections.py \\
        [src ~/data/opr/scar/node] [out figs]

Two figures, one row per segment, one panel per frame, panel width
proportional to the frame's block count, x along track in km:

  node_sections_ratio.png  background: the 45-deg cross-pol ratio (20 m
                           smoothed, each block's 200-1500 m median
                           removed). Odd nodes are its maxima (black dots),
                           even nodes its minima (white dots); candidate
                           extrema the phase test rejected are grey crosses.
                           The slow axis of each block is printed at the top.
  node_sections_dlam.png   background: the dlam(z) each product stores, the
                           same dots. A node at depth z says the two-way
                           delay integrated to z is n half waves of the
                           band's centre (0.68 ns each), so the dots should
                           ride up where dlam is high and sink where it is
                           low.

node_representation.py draws the same sections under the depth model
(node_sections_resolved.png) by calling draw() with re-solved profiles.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpecFromSubplotSpec  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts", "figures"))
from node_export import read_product, is_ns, C_NS, C_EW  # noqa: E402
import morrow_style as ms  # noqa: E402

BLK_KM = 2.0
Z_TOP, Z_BOT = 0, 1500


def frames_by_segment(rows):
    frames = {}
    for r in rows:
        frames.setdefault(r["tag"], []).append(r)
    segs = {}
    for tag in sorted(frames):
        blocks = sorted(frames[tag], key=lambda r: r["block"])
        segs.setdefault(tag[:11], []).append((tag, blocks))
    return segs


def block_edges(blocks):
    """Along-track edges [km] of each block of one frame. A frame is cut
    into equal blocks a little longer than the nominal 2 km, so a fixed
    width at each centre leaves a white sliver between neighbours: meet
    them halfway instead. A real gap (a block skipped or cut short, the
    centres more than 1.5 blocks apart) stays a gap."""
    xc = np.array([r["x"] for r in blocks]) / 1000.0
    lo, hi = xc - BLK_KM / 2, xc + BLK_KM / 2
    for i in range(xc.size - 1):
        if xc[i + 1] - xc[i] < 1.5 * BLK_KM:
            lo[i + 1] = hi[i] = (xc[i] + xc[i + 1]) / 2
    return lo, hi


def draw_block(ax, r, x_edges, background):
    """One block's column of the section and its nodes. Returns the
    mappable, or None when the block has nothing to draw."""
    xc = r["x"] / 1000.0
    if background == "ratio":
        z, v = r["z"], r["ratio_smoothed"]
        if v.size != z.size:
            return None
        base = np.nanmedian(v[(z >= 200) & (z <= 1500)])
        ze = np.concatenate([[z[0] - (z[1] - z[0]) / 2],
                             (z[:-1] + z[1:]) / 2,
                             [z[-1] + (z[-1] - z[-2]) / 2]])
        mappable = ax.pcolormesh(x_edges, ze, (v - base)[:, None],
                                 cmap="RdBu_r", vmin=-3, vmax=3,
                                 rasterized=True)
    else:
        zm, v = r["zm"], r["dlam"]
        dz = np.median(np.diff(zm))
        ze = np.concatenate([zm - dz / 2, [zm[-1] + dz / 2]])
        mappable = ax.pcolormesh(x_edges, ze, v[:, None], cmap="viridis",
                                 vmin=0, vmax=0.10, rasterized=True)
    rej = r["rejected_z"]
    if rej.size:
        ax.plot(np.full(rej.size, xc), rej, "x", color="0.45", ms=3,
                mew=0.8)
    nz, nt = r["node_z"], r["node_type"]
    if nz.size:
        mx = nt > 0
        ax.plot(np.full(mx.sum(), xc), nz[mx], "o", color="k", ms=3.2,
                mec="w", mew=0.5)
        ax.plot(np.full((~mx).sum(), xc), nz[~mx], "o", color="w", ms=3.2,
                mec="k", mew=0.7)
    ax.text(xc, 45, "%.0f" % r["theta"], ha="center", va="top", fontsize=7,
            color="k")
    return mappable


def draw(segs, out_png, background, title):
    n_rows = len(segs)
    fig = plt.figure(figsize=(18, 2.9 * n_rows + 0.9))
    outer = fig.add_gridspec(n_rows, 1, left=0.06, right=0.93, top=0.93,
                             bottom=0.05, hspace=0.6)
    mappable = None
    # one along-track scale for every row: each frame's width is its length
    # in blocks, and a row shorter than the longest is padded with an empty
    # column on the right, so a segment of one short frame is not stretched
    # across the page
    widths_all = {seg: [max(1, max(r["block"] for r in blocks))
                        for _, blocks in frames]
                  for seg, frames in segs.items()}
    row_len = {seg: sum(w) for seg, w in widths_all.items()}
    full = max(row_len.values())
    for i, (seg, frames) in enumerate(segs.items()):
        widths = widths_all[seg]
        pad = full - row_len[seg]
        ratios = widths + ([pad] if pad > 0 else [])
        inner = GridSpecFromSubplotSpec(1, len(ratios),
                                        subplot_spec=outer[i],
                                        width_ratios=ratios, wspace=0.06)
        for j, (tag, blocks) in enumerate(frames):
            ax = fig.add_subplot(inner[j])
            fam_ns = is_ns(np.median([r["heading"] for r in blocks]))
            lo, hi = block_edges(blocks)
            for r, x_lo, x_hi in zip(blocks, lo, hi):
                drawn = draw_block(ax, r, np.array([x_lo, x_hi]), background)
                mappable = drawn if drawn is not None else mappable
            x_all = np.array([r["x"] for r in blocks]) / 1000.0
            ax.set_xlim(lo.min(), hi.max())
            ax.set_ylim(Z_BOT, Z_TOP)
            ax.set_title("f%s %s" % (tag[-3:], "N-S" if fam_ns else "E-W"),
                         fontsize=10, color=C_NS if fam_ns else C_EW,
                         loc="left", pad=2)
            ax.tick_params(labelsize=9, length=2)
            if x_all.max() - x_all.min() > 6:
                ax.set_xticks(np.arange(0, x_all.max() + 1, 5))
            elif x_all.size > 1:
                ax.set_xticks([round(x_all.min(), 0), round(x_all.max(), 0)])
            else:
                ax.set_xticks([round(x_all[0], 0)])
            if j == 0:
                ax.set_ylabel("%s\ndepth [m]" % seg, fontsize=10)
                ax.set_xlabel("along track [km]", fontsize=9, labelpad=1)
                ms.tag(None, "ABCDEFGH"[i], fig=fig, x=0.018,
                       y=ax.get_position().y1 + 0.004, fontsize=13)
            else:
                ax.set_yticklabels([])
    if mappable is not None:
        cax = fig.add_axes([0.945, 0.3, 0.008, 0.4])
        cb = fig.colorbar(mappable, cax=cax)
        cb.ax.tick_params(labelsize=10)
        cb.set_label("cross-pol ratio anomaly [dB]" if background == "ratio"
                     else r"$\Delta\lambda$", fontsize=12)
    ms.label(None, title, fig=fig, x=0.06, y=0.975, va="top", fontsize=12)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main():
    ms.apply()
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.expanduser(os.environ.get("SCAR_DATA", "~/data/opr/scar")),
        "node")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figs")
    os.makedirs(out, exist_ok=True)
    rows = []
    for fn in sorted(glob.glob(os.path.join(src, "node_*_B2000.mat"))):
        rows += [r for r in read_product(fn) if r.get("ok")]
    segs = frames_by_segment(rows)
    print("%d blocks in %d frames, %d segments" % (
        len(rows), sum(len(v) for v in segs.values()), len(segs)))
    draw(segs, os.path.join(out, "node_sections_ratio.png"), "ratio",
         r"half-wave nodes along track on the 45$^\circ$ cross-pol ratio: "
         r"black = maxima (odd nodes), white = minima (even), "
         r"grey $\times$ = rejected candidates; "
         r"number = slow axis [$^\circ$ geographic]")
    draw(segs, os.path.join(out, "node_sections_dlam.png"), "dlam",
         r"the same nodes on the stored $\Delta\lambda$(z) of each block")
    print("wrote node_sections_ratio.png, node_sections_dlam.png in", out)


if __name__ == "__main__":
    main()

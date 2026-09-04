#!/usr/bin/env python3
"""Where is the "Kamb" data, really? The January 2023 tracks against both places.

The 2022_Antarctica_Ground season carries segments 20230120_* that a
processing script in the group scratch calls Kamb. This asks the
coordinates rather than the filename: it draws every moving segment of
January 2023 with WAIS Divide camp and the Kamb Ice Stream trunk marked,
at a scale that holds both, and again zoomed on the survey itself.

  python3 scripts/figures/kamb_or_wais.py <tracks.mat> [out_dir=figs]
"""
import os
import sys
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer

SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "figs"
os.makedirs(OUT, exist_ok=True)

WAIS = (-79.468, -112.086)      # WAIS Divide camp
KAMB = (-82.40, -135.50)        # Kamb Ice Stream trunk, approximate
T3031 = Transformer.from_crs("EPSG:4326", "EPSG:3031", always_xy=True)
R = 6371e3


def gc(a1, o1, a2, o2):
    p1, p2 = np.radians(a1), np.radians(a2)
    dl = np.radians(o2 - o1)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


d = loadmat(SRC, squeeze_me=True, struct_as_record=False)["T"]
segs = [(str(s.seg), np.atleast_1d(s.lat), np.atleast_1d(s.lon), float(s.len_km))
        for s in np.atleast_1d(d)]
segs = [s for s in segs if s[3] >= 2.0]          # drop static/short segments

# the outbound leg toward Kamb: segments whose far end is furthest from camp
far = {}
for name, la, lo, L in segs:
    far[name] = gc(la, lo, *WAIS).max() / 1e3
outbound = sorted(far, key=far.get, reverse=True)[:8]

wx, wy = T3031.transform(WAIS[1], WAIS[0])
kx, ky = T3031.transform(KAMB[1], KAMB[0])

# two rows, not two columns: the regional view is 520 km wide and 160 km
# tall, so side-by-side panels at equal aspect leave most of the canvas empty
fig, ax = plt.subplots(2, 1, figsize=(13, 11),
                       gridspec_kw={'height_ratios': [1, 1.5]})

# ---- panel 1: both places in one frame
a = ax[0]
for name, la, lo, L in segs:
    x, y = T3031.transform(lo, la)
    a.plot(np.asarray(x) / 1e3, np.asarray(y) / 1e3, "-",
           color="#c1272d" if name in outbound else "#1f4e79", lw=1.4, zorder=3)
a.plot(wx / 1e3, wy / 1e3, "o", ms=11, mfc="w", mec="k", mew=1.8, zorder=5)
a.annotate("WAIS Divide camp", (wx / 1e3, wy / 1e3), textcoords="offset points",
           xytext=(14, 6), fontsize=11, weight="bold")
a.plot(kx / 1e3, ky / 1e3, "*", ms=20, mfc="#f0a202", mec="k", mew=1.2, zorder=5)
a.annotate("Kamb Ice Stream trunk", (kx / 1e3, ky / 1e3), textcoords="offset points",
           xytext=(-12, -22), fontsize=11, weight="bold", ha="right")
a.plot([wx / 1e3, kx / 1e3], [wy / 1e3, ky / 1e3], "--", color="0.45", lw=1.4, zorder=2)
mid = ((wx + kx) / 2e3, (wy + ky) / 2e3)
a.annotate("%.0f km" % (gc(WAIS[0], WAIS[1], KAMB[0], KAMB[1]) / 1e3), mid,
           textcoords="offset points", xytext=(6, 8), fontsize=11, color="0.3")
a.margins(x=0.06, y=0.12)
a.set_title("The whole January 2023 season, and where Kamb actually is")

# ---- panel 2: the survey itself
b = ax[1]
for name, la, lo, L in segs:
    x, y = T3031.transform(lo, la)
    b.plot(np.asarray(x) / 1e3, np.asarray(y) / 1e3, "-",
           color="#c1272d" if name in outbound else "#1f4e79", lw=1.6, zorder=3)
b.plot(wx / 1e3, wy / 1e3, "o", ms=11, mfc="w", mec="k", mew=1.8, zorder=5)
b.annotate("WAIS Divide camp", (wx / 1e3, wy / 1e3), textcoords="offset points",
           xytext=(12, 6), fontsize=11, weight="bold")
for name in ("20230120_05", "20230120_04"):
    m = [s for s in segs if s[0] == name]
    if m:
        _, la, lo, L = m[0]
        x, y = T3031.transform(lo[-1], la[-1])
        b.annotate("%s\n%.0f km" % (name, L), (x / 1e3, y / 1e3),
                   textcoords="offset points", xytext=(8, -4), fontsize=9, color="#c1272d")
b.set_title("Zoom: a local grid plus an 84 km leg toward Kamb (red)")

for a_ in ax:
    a_.set_aspect("equal")
    a_.set_xlabel("EPSG:3031 easting (km)")
    a_.set_ylabel("EPSG:3031 northing (km)")
    a_.grid(alpha=0.3)
fig.suptitle("2022_Antarctica_Ground, January 2023: the segments called 'Kamb' are a "
             "WAIS Divide survey with a southwest spur, not a line down Kamb",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fn = os.path.join(OUT, "kamb_or_wais.png")
fig.savefig(fn, dpi=115)
print("wrote", fn)

allla = np.concatenate([s[1] for s in segs])
alllo = np.concatenate([s[2] for s in segs])
dw = gc(allla, alllo, *WAIS) / 1e3
dk = gc(allla, alllo, *KAMB) / 1e3
print("%d moving segments, %.0f km of track" % (len(segs), sum(s[3] for s in segs)))
print("distance from WAIS Divide: median %.1f km, max %.1f km" % (np.median(dw), dw.max()))
print("distance from Kamb trunk : min %.0f km, max %.0f km" % (dk.min(), dk.max()))

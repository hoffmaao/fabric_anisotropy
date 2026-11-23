"""The figure style of Morrow, Ghelichkhan & Kramer (2026): serif Times
text, STIX math, circled light-grey panel letters, units in brackets, no
figure titles and no keys drawn on the panels (the caption carries them;
label() is for a word or two that names a panel's data).

    import morrow_style as ms
    ms.apply()                      # before any figure is made
    ms.tag(ax, "A")                 # circled panel letter, upper left
    ms.label(ax, "text")            # boxed label, lower left by default
"""
import logging

import matplotlib.pyplot as plt

RC = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Liberation Serif",
                   "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 12,
    "axes.labelsize": 13,
    "xtick.labelsize": 11.5,
    "ytick.labelsize": 11.5,
    "legend.fontsize": 12,
    "axes.linewidth": 1.0,
    "axes.edgecolor": "black",
    "figure.facecolor": "white",
}
TAG = dict(ha="center", va="center", fontsize=15, zorder=20,
           bbox=dict(boxstyle="circle,pad=0.3", facecolor="lightgrey",
                     edgecolor="black"))
LABEL = dict(zorder=20, bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor="lightyellow",
                                  edgecolor="black"))
FRAME = (12.8, 7.2)          # 1920 x 1080 at 150 dpi, the Morrow movie frame


def apply():
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams.update(RC)


def tag(ax, letter, x=0.045, y=0.94, **kw):
    """circled panel letter in axes coordinates (2D axes) or figure
    coordinates (pass fig=...)"""
    fig = kw.pop("fig", None)
    if fig is not None:
        return fig.text(x, y, letter, **{**TAG, **kw})
    return ax.text(x, y, letter, transform=ax.transAxes, **{**TAG, **kw})


def label(ax, text, x=0.02, y=0.03, ha="left", va="bottom", fontsize=12, **kw):
    """boxed light-yellow label in axes coordinates (or figure coordinates
    with fig=...)"""
    fig = kw.pop("fig", None)
    d = dict(ha=ha, va=va, fontsize=fontsize, **LABEL)
    d.update(kw)
    if fig is not None:
        return fig.text(x, y, text, **d)
    return ax.text(x, y, text, transform=ax.transAxes, **d)

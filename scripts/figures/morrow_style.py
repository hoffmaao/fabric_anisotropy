"""The figure style of Morrow, Ghelichkhan & Kramer (2026): serif Times
text, STIX math, circled light-grey panel letters, units in brackets, no
figure titles and no keys drawn on the panels (the caption carries them).

    import morrow_style as ms
    ms.apply()                      # before any figure is made
    ms.tag(ax, "A")                 # circled panel letter, upper left
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

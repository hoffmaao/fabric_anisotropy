"""Shared slide styling and geometry for the SCAR site figures.

scar_two_panel.py (Thwaites, Ridge A), negis_two_panel.py (NEGIS),
egrip_two_panel.py (the EastGRIP borehole), fabric_sections.py and
fabric_three_sites.py are one deck. The site
slides only read as a set while their end markers, track colours, profile
weight, scale bars and panel proportions stay identical, so those
constants and the helpers that draw them live here rather than being
copied per script - a copy is exactly the thing that desynchronises the
deck the next time one of them is tweaked.

Everything that varies per site is a parameter. Above all the projection:
EPSG:3031 for the Antarctic sites, polar stereographic north for NEGIS.
Cartopy is imported lazily inside the two helpers that need it, so the
section figures (which draw no map) keep working without it.

Slide styling (Andrew, 4 Aug): no legends - the titles already name the
segment, and a legend box eats map area on a projected slide. The
profile's two ENDS carry the map-to-interferogram correspondence instead
of a distance-marker ladder: a red marker at the start and a white marker
at the end, both with a black outline, drawn on the map inset AND at the
matching left/right edges of the interferogram. Two unambiguous ends read
from the back of a room where five graded viridis dots do not.
"""
import os

# pyplot is imported INSIDE the two functions that need it, never at module
# scope. Importing pyplot RESOLVES AND BINDS A BACKEND, so a module-level
# import here would bind one the moment any script imported this file - and
# every consumer is a batch figure script whose first act is
# matplotlib.use('Agg'). One did land its `import scar_style` above that
# guard, and the interactive macosx backend was imported before the guard
# ran. Keeping pyplot out of module scope removes the ordering contract
# rather than restating it, so import position here cannot matter again.
import matplotlib.patheffects as pe
import numpy as np
from matplotlib import transforms
from matplotlib.colors import hsv_to_rgb

# Inputs that are data, not code: ITS_LIVE windows streamed from S3, the
# NEGIS track/interferogram extracts pulled off mem1, the fabric section
# stages copied back from the server. They live beside the rest of the
# mirrored products rather than in the repo.
DATA = os.path.expanduser(os.environ.get('SCAR_DATA', '~/data/opr/scar'))

# Outputs are data too, and are gitignored for the same reason - but they
# belong beside the code that made them, so they land in the repository's
# own figs/ rather than in a directory outside the tree. Resolved from this
# file's location so a script run from anywhere writes to the same place,
# and created on import so savefig never fails on a missing directory.
# Override per run with a script's <out_dir> argument, or globally with
# FABRIC_FIGS.
FIGS = os.environ.get('FABRIC_FIGS') or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..',
                 'figs'))
os.makedirs(FIGS, exist_ok=True)


def out_dir(argv, pos=1, default=None):
    """Output directory from argv[pos], falling back to `default` or FIGS.

    An EMPTY argument means "use the default", not "use the cwd". Scripts here
    take the output directory first and their selectors after it, so choosing a
    site or a variant forces argv[1] to be supplied and '' is the natural way
    to leave the output alone. Resolving that to the cwd drops figures and
    movies into the repository root - the one path that is not the gitignored
    figs/ - with nothing warning that it happened.

    `pos` indexes whatever list is passed, so a script that has already
    sliced or spliced sys.argv passes its own list and its own position.
    `default` is for a script whose output genuinely belongs somewhere other
    than FIGS - a prototype rooted in its own scratch tree, or a cache of
    intermediate frames - so it can take the empty-argument guard without
    having its destination moved.

    The directory returned is CREATED. FIGS exists because this module makes
    it at import, but a `default` root or a caller-supplied path had no such
    guarantee, and three prototypes raised FileNotFoundError from savefig
    whenever they ran anywhere their hand-rolled relative default did not
    already exist. Creating it here is what makes the returned path usable
    rather than merely resolved, so no caller has to remember the makedirs.
    """
    val = argv[pos].strip() if len(argv) > pos else ''
    d = val or (FIGS if default is None else default)
    os.makedirs(d, exist_ok=True)
    return d


FIGSIZE = (13.0, 5.8)
DPI = 200

INK, MUTED = '#0b0b0b', '#52514e'

# Profile endpoints: start red, end white, both outlined black. The same
# pair marks the interferogram's left and right edges, which is what
# makes the correspondence readable without a legend.
END_FACES = ('#e8000b', '#ffffff')
END_SIZE = 130
END_EDGE = 1.6

# Both colorbar captions are deliberately terse, and are NOT to be
# re-expanded with the provenance and scaling they omit. These are
# conference-talk slides, narrated live: the speed layer's source
# (ITS_LIVE v2.1) and its log scaling, and the coherence on the HSV value
# channel, are all said out loud and written in the abstract, so on the
# slide they are words the audience reads instead of looking at the
# figure. One constant each, because three panels draw the speed bar and
# four draw the interferogram bar.
SPEED_CB_LABEL = 'surface speed (m/yr)'
IFG_CB_LABEL = 'phase change (rad)'

TRACK_C = 'white'       # every survey line
PROFILE_C = 'black'     # the focused profile
PROFILE_LW = 2.8
# Thin white casing under the black profile line, so it survives the dark
# end of the speed ramp - slow ice is exactly where these profiles sit -
# without the casing widening into something that reads as a white line
PROFILE_FX = [pe.withStroke(linewidth=PROFILE_LW + 1.4, foreground='white')]


def nice_length_up(target_m):
    """Smallest 1/2/5 x power of ten at or above `target_m`."""
    for length in [1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5]:
        if length >= target_m:
            return length
    return 5e5


def cumdist_km(lats, lons):
    """Cumulative great-circle distance along a track, in km."""
    R = 6371.0
    p = np.radians(np.asarray(lats, float))
    dl = np.radians(np.diff(np.asarray(lons, float)))
    dp = np.diff(p)
    a = np.sin(dp / 2)**2 + np.cos(p[:-1]) * np.cos(p[1:]) * np.sin(dl / 2)**2
    return np.concatenate([[0], np.cumsum(2 * R * np.arcsin(np.sqrt(a)))])


def hsv_plot_coherence(phase, coh, coherence_limits=(0.0, 1.0)):
    """Python port of OPR display/hsv_plot_coherence.m (John Paden):
    hue = normalized phase (-pi..pi -> 0..1), sat = 1, val = coherence
    clipped to coherence_limits and scaled 0..1. Low-coherence pixels
    fade to black instead of shouting random fringe colours."""
    hue = phase / (2 * np.pi) + 0.5
    val = np.clip(np.abs(coh), *coherence_limits)
    val = (val - coherence_limits[0]) / (
        coherence_limits[1] - coherence_limits[0])
    val = np.nan_to_num(val, nan=0.0)
    hsv = np.stack([np.clip(hue, 0, 1), np.ones_like(hue), val], axis=-1)
    return hsv_to_rgb(hsv)


def scale_bar_br(ax, extent, proj, frac=0.25, nice_length=nice_length_up):
    """Plain black scale bar in the BOTTOM-RIGHT corner.

    No white casing on the bar and no white box behind the label: every
    site places the bar over the pale end of its background, so the
    casing only added visual weight.

    `extent` is in `proj` metres. `nice_length` is the rounding rule for
    the bar length; the Antarctic panels pass antarctic_basemap's so the
    bar matches the other figures in that family.
    """
    width = extent[1] - extent[0]
    height = extent[3] - extent[2]
    length = nice_length(frac * width)
    x1 = extent[1] - 0.06 * width
    x0 = x1 - length
    y0 = extent[2] + 0.045 * height
    ax.plot([x0, x1], [y0, y0], color='black', lw=3, transform=proj,
            zorder=9, solid_capstyle='butt')
    label = f'{length/1000:g} km' if length >= 1000 else f'{length:g} m'
    ax.text((x0 + x1) / 2, y0 + 0.02 * height, label, ha='center',
            va='bottom', fontsize=9, color='black', transform=proj,
            zorder=9)


def locator_box(ax, zext, proj, color='white'):
    """Outline the zoom window on the main map.

    Dashed, where the inset's own frame is solid: at slide size the two
    rectangles are otherwise the same white box and it is not obvious
    which one is the magnified view.
    """
    ax.plot([zext[0], zext[1], zext[1], zext[0], zext[0]],
            [zext[2], zext[2], zext[3], zext[3], zext[2]],
            '--', color=color, lw=1.0, dashes=(4, 2), transform=proj,
            zorder=9)


def inset_frame(axz, color='white', lw=1.5):
    """Solid frame on a zoom inset, against the dashed locator box."""
    for spine in axz.spines.values():
        spine.set_edgecolor(color)
        spine.set_linewidth(lw)


# Continental locators. Natural Earth 50m land is already cached under
# ~/.local/share/cartopy, so these need no network.
CONTINENT = {
    'antarctica': ('SouthPolarStereo', {}, [-180, 180, -90, -63]),
    'greenland': ('NorthPolarStereo', {'central_longitude': -42},
                  [-58, -8, 58.5, 84]),
}
LOCATOR_LAND = '0.86'
LOCATOR_EDGE = '0.45'
# Floor on the drawn AOI box, as a fraction of the locator's width.
# True footprints here are 0.7-4% of the ice sheet, sub-pixel at
# three of the four sites.
AOI_MIN_FRAC = 0.055


def continental_inset(ax, region, rect, extent, proj):
    """Small locator putting the study area on its ice sheet.

    These slides get shown to people who do not know where Ridge A or the
    EastGRIP borehole are, and the survey-scale panel cannot tell them: at
    that zoom every site is an anonymous patch of white. One box on the
    continent fixes it.

    The area of interest is drawn as a black square from the map panel's
    own `extent`, so it is the study area rather than a symbol placed near
    it. It is FLOORED at AOI_MIN_FRAC of the locator's width, because the
    true footprint is 0.7-4% of the ice sheet across these four sites and
    would be one pixel or less at three of them - a box that cannot be
    seen locates nothing. Read the square as "here", not as a scale bar.

    Lifted above the parent's layers for the same reason as zoom_inset -
    ax.inset_axes() registers the child inside the PARENT's artist
    ordering at zorder 5, below the map's own quiver, tracks and profile.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from matplotlib.patches import Rectangle
    name, kw, cext = CONTINENT[region]
    iproj = getattr(ccrs, name)(**kw)
    axc = ax.inset_axes(rect, projection=iproj)
    axc.set_zorder(ax.get_zorder() + 21)
    axc.set_extent(cext, crs=ccrs.PlateCarree())
    axc.add_feature(cfeature.LAND.with_scale('50m'), facecolor=LOCATOR_LAND,
                    edgecolor=LOCATOR_EDGE, lw=0.4)
    axc.set_facecolor('white')

    x0, x1, y0, y1 = extent
    lon, lat = ccrs.PlateCarree().transform_point(
        (x0 + x1) / 2, (y0 + y1) / 2, proj)[:2]
    ix, iy = iproj.transform_point(lon, lat, ccrs.PlateCarree())[:2]
    cx0, cx1 = axc.get_xlim()
    half = max(abs(x1 - x0) / 2, abs(y1 - y0) / 2,
               AOI_MIN_FRAC * abs(cx1 - cx0) / 2)
    axc.add_patch(Rectangle((ix - half, iy - half), 2 * half, 2 * half,
                            transform=iproj, facecolor='none',
                            edgecolor='black', lw=1.5, zorder=6))
    for sp in axc.spines.values():
        sp.set_edgecolor('0.25')
        sp.set_linewidth(0.9)
    return axc


def zoom_inset(ax, rect, proj):
    """Zoom inset that composites ABOVE the parent map's own layers.

    `ax.inset_axes()` gives the child zorder 5 (not 0) AND registers it
    with `add_child_axes`, so it lands in `ax.child_axes` rather than
    `fig.axes` and is drawn inside the PARENT's artist ordering. These maps
    draw their quiver, survey tracks, profile and scale bar at explicit
    zorders 6-9, all above 5, so the parent's layers paint straight over
    the inset: the result reads as two sets of arrows at two different
    scales inside one box.

    Measured, rather than reasoned: in a controlled render the parent's
    arrows inside the inset drop from 469 to 137 px with the zorder bump
    and are unchanged without it. Making the patch opaque does NOT help
    and is not done here - an inset axes' patch is already visible with an
    opaque white facecolor by default, so those calls are no-ops.
    """
    axz = ax.inset_axes(rect, projection=proj)
    axz.set_zorder(ax.get_zorder() + 20)
    return axz


def draw_map_ends(ax, lats, lons):
    """Red start / white end markers on the profile, black outlined."""
    import cartopy.crs as ccrs
    ax.scatter([lons[0], lons[-1]], [lats[0], lats[-1]], s=END_SIZE,
               c=list(END_FACES), edgecolor='black', lw=END_EDGE,
               transform=ccrs.PlateCarree(), zorder=10)


def draw_ifg_ends(ax, dist):
    """The same two markers at the interferogram's left and right edges."""
    tr = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.scatter([dist[0], dist[-1]], [1.035, 1.035], s=END_SIZE,
               c=list(END_FACES), edgecolor='black', lw=END_EDGE,
               transform=tr, clip_on=False, zorder=10)


def two_panel_figure(proj):
    """The shared slide geometry: map panel left, interferogram right."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI, layout='constrained')
    gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.18])
    axm = fig.add_subplot(gs[0, 0], projection=proj)
    axi = fig.add_subplot(gs[0, 1])
    return fig, axm, axi


def ifg_panel(ax, fig, dist, t_us, phase, coh, ylabel='TWTT (μs)',
              cb_label=IFG_CB_LABEL):
    """Wrapped interferogram with the shared HSV mapping and colorbar."""
    st = max(1, phase.shape[0] // 2200)
    sx = max(1, phase.shape[1] // 2000)
    rgb = hsv_plot_coherence(phase[::st, ::sx], coh[::st, ::sx])
    ax.imshow(rgb, aspect='auto', interpolation='nearest',
              extent=[dist[0], dist[-1], t_us[-1], t_us[0]])
    ax.set_xlabel('distance along profile (km)')
    ax.set_ylabel(ylabel)
    import matplotlib.pyplot as plt
    sm = plt.cm.ScalarMappable(cmap='hsv',
                               norm=plt.Normalize(-np.pi, np.pi))
    # Clear of the image frame: at pad=0.01 the axes spine and the
    # colorbar edge merge into one line on a projected slide
    cb = fig.colorbar(sm, ax=ax, pad=0.035)
    cb.set_label(cb_label)
    cb.set_ticks([-np.pi, 0, np.pi])
    cb.set_ticklabels([r'$-\pi$', '0', r'$\pi$'])


def field(rec, name):
    """One field of a squeeze_me'd MATLAB struct, as a 1-D float array."""
    v = rec[name]
    while isinstance(v, np.ndarray) and v.dtype == object and v.size == 1:
        v = v.item()
    return np.atleast_1d(np.asarray(v, float))


def track_azimuth(lat, lon):
    """Principal azimuth of the track, degrees E of N, as an AXIS (0-180).

    The two eigenvalues being differenced sit on the profile axis and its
    perpendicular, so the figures have to name those directions - dlam is
    meaningless without them. An axis, not a heading: eigenvalues are
    invariant under a 180 deg flip, so travelling the line either way
    gives the same pair. Fitted over the whole track rather than
    endpoint-to-endpoint so a wandering line still reports its trend.
    """
    ok = np.isfinite(lat) & np.isfinite(lon)
    la, lo = lat[ok], lon[ok]
    y = la - la.mean()
    x = (lo - lo.mean()) * np.cos(np.radians(la.mean()))
    pts = np.c_[x, y]
    pts = pts - pts.mean(0)
    _, _, vt = np.linalg.svd(pts, full_matrices=False)
    dx, dy = vt[0]
    return np.degrees(np.arctan2(dx, dy)) % 180

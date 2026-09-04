"""Survey sites for the quad-pol figures: where they are and how to draw them.

One definition shared by the map and the depth-movie scripts, so a site's
depth band, colour ceiling and identity cannot drift between the still and
the animation of the same data.

SITES ARE SELECTED BY POSITION, NOT BY TAG PREFIX. Prefix selection is
actively wrong across this survey set: the 2024 tags alone cover Thwaites
(76.5 S, 107.7 W), WAIS Divide (79.2 S, 111.6 W) and the McMurdo Ice Shelf
(77.7 S, 168.1 E), so a prefix map would draw three sites hundreds of km
apart on one axis and compute a single "survey" statistic across all of
them.

Each entry carries its OWN bands and ceiling because the sites are not
comparable. Ridge A reaches dlam 0.066 over 1450 m of ice; the McMurdo
record stops near 840 m and reads 0.0003. Drawing the latter on Ridge A's
0-0.08 ramp renders a whole survey blank white and implies it was measured
as isotropic when it was mostly not resolved at all.

CAUTION on the three weak sites. Thwaites, WAIS Divide and McMurdo all sit
at or below dlam 0.005, under the estimator's own 0.01 abstention threshold
for reporting an axis, so their maps largely show "not resolved" rather
than "measured isotropic". Do not quote a site value from them without
checking coverage first - the local mirror holds only 2 WAIS Divide and 3
McMurdo frames, and the Thwaites median here (0.004) disagrees with the
0.13-0.19 margin contrasts recorded from the earlier margin work.
"""
import numpy as np
from pyproj import Transformer

T3031 = Transformer.from_crs('EPSG:4326', 'EPSG:3031', always_xy=True)

# MOVIE DEPTH RANGES follow the Ridge A protocol: stop where the data stops
# behaving, not where the record ends. At Ridge A that cut was 1050 m,
# where its mostly N-S lines leave the regular fabric cycle. The measurable
# analogue applied at the other sites is the estimator's OWN abstention
# threshold - it declines to report an axis below dlam 0.01 - so each movie
# stops in the window where the site's median dlam crosses that line.
# Measured per 150 m window over every frame at each site:
#
#   taylor_dome  0.025 0.032 0.027 0.018 | 0.006 0.007 0.010 ...  cut 650 m
#                                          SUPERSEDED: cut 850 m, from the bed
#   thwaites     0.033 .. 0.012 at 900-1050 | 0.008 0.005 0.004   cut 1050 m
#   wais_divide  0.011 0.007 0.009 0.022 | 0.005 0.001 ...        cut 700 m
#   mcmurdo      0.012 | 0.004 0.002 0.0004 0.0003                cut 450 m
#
# So three of the four run out of resolvable fabric in the upper column,
# and their movies are correspondingly short. That is the honest result:
# animating the full record would sweep mostly through values the estimator
# has already declared unresolvable, which reads as measured isotropy.
# WINDOW AND STEP ARE PER SITE, because a window is only meaningful
# relative to the ice it sits in. The 150 m window that suits Ridge A's
# 1850 m column is half the ice at the McMurdo Ice Shelf, where the sounding
# reaches only 200-300 m: it would smear the whole profile into two
# independent windows and animate almost nothing. Thin-ice sites get 40 m
# windows at a 5 m step, which the 1.12 m sample spacing easily supports.
# COLOUR CEILINGS ARE MEASURED, NOT CHOSEN. vmax is the 95th percentile of
# the quantity the movie actually draws - the per-block median dlam over one
# depth window - taken over WELL-FIT windows only (LS residual < 0.25, the
# level Ridge A holds through its whole column).
#
# They were previously set from each site's per-window MEDIAN dlam, which
# clips half the map by construction: 49% of cells at Taylor Dome, 57% at
# Thwaites, 47% at McMurdo. Ridge A escaped only because its value happened
# to sit at its p95, and the rule below reproduces that 0.08 exactly, which
# is why it is trusted for the others.
#
# The well-fit restriction matters as much as the percentile. Over ALL
# windows the tail is the fit failing rather than the ice: gating from
# resid<1.0 to resid<0.25 drops the displayed p95 from 0.119 to 0.108 at
# Taylor Dome and from 0.181 to 0.126 at Thwaites, while Ridge A holds at
# 0.083 -> 0.079. quadpol_depth_movie.py fades toward grey on the same
# residual, so what the ceiling scales is what the map shows at full colour.
SITES = {
    'ridge_a': dict(lat=-86.65, lon=69.35, title='Ridge A grid',
                    z_band=(200.0, 1200.0), z_deep=(1150.0, 1500.0),
                    vmax=0.08, movie=(200.0, 1050.0), win=150.0, step=20.0),
    # TAYLOR DOME's movie floor is set from the BED, not from where dlam
    # crosses the abstention threshold. CSARP_layer picks the bottom on 51
    # of these frames and puts the ice at 566-996 m, median 830 - so the
    # old 650 m cut, chosen where median dlam fell under 0.01, was reading
    # a FABRIC limit off real ice and discarding the most interesting part
    # of the column: the contrast peaks near z/H 0.5 and decays back
    # toward zero by z/H 0.85, at residuals of 0.10-0.32 throughout, which
    # is what a dome should do as its fabric becomes axisymmetric.
    'taylor_dome': dict(lat=-77.78, lon=158.68, title='Taylor Dome',
                        z_band=(200.0, 1200.0), z_deep=(1150.0, 1500.0),
                        vmax=0.11, movie=(200.0, 850.0), win=150.0, step=20.0),
    'thwaites': dict(lat=-76.45, lon=-107.67, title='Thwaites margin',
                     z_band=(200.0, 1200.0), z_deep=(1000.0, 1450.0),
                     vmax=0.13, movie=(200.0, 1050.0), win=150.0, step=20.0),
    # THE KAMB LINE. The January 2023 traverse leg that leaves WAIS Divide
    # camp on a bearing of about 200 deg and runs 84 km southwest, into the
    # Ross-side catchment that drains through Kamb Ice Stream, then returns
    # on 21 January. It is matched by SEGMENT, not by position, and it has
    # to be: the leg sits 0.9 deg from the WAIS Divide centre below, well
    # inside the 2 deg radius, so on position alone every frame of it is
    # absorbed into the camp grid and labelled WAIS Divide - which is what
    # happened the first time a figure was made from 20230120_05.
    #
    # NAMING, stated so nobody is misled by it later: this is the line
    # TOWARD Kamb, in its upper catchment. The Kamb Ice Stream trunk is at
    # about -82.4, -135.5, which is 518 km from WAIS Divide camp, and the
    # farthest point of this leg is still 440 km from it
    # (scripts/figures/kamb_or_wais.py, figs/kamb_or_wais.png). Nothing
    # here is a measurement of the trunk, and a section from it must not be
    # presented as one.
    'kamb': dict(lat=-79.90, lon=-113.20, title='Kamb line',
                 segs=('20230118_11', '20230118_13', '20230120_04',
                       '20230120_05', '20230121_02', '20230121_03',
                       '20230121_04'),
                 z_band=(200.0, 1200.0), z_deep=(1000.0, 1450.0),
                 vmax=0.25, movie=(200.0, 1400.0), win=150.0, step=20.0),
    'wais_divide': dict(lat=-79.22, lon=-111.59, title='WAIS Divide',
                        z_band=(200.0, 1200.0), z_deep=(1000.0, 1450.0),
                        vmax=0.06, movie=(200.0, 700.0), win=150.0, step=20.0),
    # Eastwind and McMurdo sit within 200 m of each other but are DIFFERENT
    # EXPERIMENTS, so they are separated by season as well as position: the
    # 2022 frames are 13 soundings all at one spot (a rotation experiment),
    # while the 2024 frames are a 15 km transect across the shelf. A purely
    # spatial box merges them and averages a strong single-point fabric into
    # a weak transect. Coherence puts the usable ice at ~300 m for both;
    # below that Eastwind's |C| collapses to 0.02 by 400 m while its
    # apparent dlam climbs to 0.10, which is noise, not fabric, and must
    # never be animated as though it were.
    'eastwind': dict(lat=-77.6715, lon=168.135, title='Eastwind',
                     years=('2022',),
                     z_band=(60.0, 280.0), z_deep=(200.0, 300.0),
                     vmax=0.06, movie=(50.0, 300.0), win=40.0, step=5.0),
    # McMurdo: the HH-VV coherence measured per 30 m window is 0.48-0.59
    # through 50-130 m, falls to ~0.10 by 150 m, recovers to 0.27-0.35 at
    # 250-305 m - most likely the ice base - and collapses to a 0.015 noise
    # floor below 375 m. So the sounding carries ice to roughly 300 m and
    # the earlier 150-450 m movie was animating largely sub-bed noise; the
    # "signal ends at 450 m" cut was detecting the bed, not a fabric limit.
    'mcmurdo': dict(lat=-77.74, lon=167.84, title='McMurdo Ice Shelf',
                    years=('2024',),
                    z_band=(60.0, 280.0), z_deep=(200.0, 300.0),
                    vmax=0.04, movie=(50.0, 300.0), win=40.0, step=5.0),
}
SITE_RADIUS_DEG = 2.0


def get(site):
    if site not in SITES:
        raise SystemExit('site must be one of: %s' % ', '.join(sorted(SITES)))
    return SITES[site]


def in_site(cfg, la, lo, tag=None):
    """Is this frame inside the named site?

    A site may pin `segs`, a tuple of day_seg strings, and then membership
    is decided by the segment ALONE - position is not consulted at all.
    That is for a site which is a chosen subset of a larger survey rather
    than a separate place: the Kamb line leaves WAIS Divide camp and stays
    within its radius the whole way, so no box can separate the two.
    Sites pinning `segs` must be declared BEFORE the survey they sit
    inside, since callers take the first match.

    A site may pin `years`, in which case the tag's season must match too.
    Position alone cannot separate Eastwind from the McMurdo Ice Shelf
    transect: they overlap spatially and are told apart by season.

    Longitude tolerance is widened because a degree of longitude is only a
    few km at these latitudes - at Ridge A's 86.6 S it is 6.7 km, so the
    same physical radius spans far more degrees than it does at 76 S.
    """
    if cfg.get('segs'):
        # a segment-pinned site needs the tag; without one it cannot match,
        # rather than silently falling through to its position box
        if tag is None:
            return False
        return any(tag.startswith(sg) for sg in cfg['segs'])
    if cfg.get('years') and tag is not None:
        if tag[:4] not in cfg['years']:
            return False
    d_lon = abs((lo - cfg['lon'] + 180.0) % 360.0 - 180.0)
    return (abs(la - cfg['lat']) < SITE_RADIUS_DEG
            and d_lon < SITE_RADIUS_DEG * 3.0)


def square_extent(bx, by, pad=1.18):
    """A SQUARE extent centred on the data, in projected metres.

    The surveys range from Ridge A's 27 x 23 km block to WAIS Divide's
    30 x 1.2 km line. Fitting the axes to the data would give a figure a
    slide cannot use and would make two sites impossible to compare by eye;
    padding the short axis instead keeps the geometry undistorted and every
    site on the same square canvas. `span` is also what sizes the markers,
    so bars stay proportionate to the survey rather than to the page.
    """
    cx = 0.5 * (float(np.min(bx)) + float(np.max(bx)))
    cy = 0.5 * (float(np.min(by)) + float(np.max(by)))
    span_x = float(np.max(bx) - np.min(bx))
    span_y = float(np.max(by) - np.min(by))
    half = 0.5 * pad * max(span_x, span_y, 1.0)
    return (cx - half, cx + half, cy - half, cy + half), span_x, span_y


def grid_north_az(lat, lon, step=0.02):
    """Azimuth of GRID north [deg E of true north] at each lat/lon.

    Measured, not assumed: step a little way true north, see where that
    lands in projected coordinates, and read off the angle between true
    north and the grid's +y axis. On EPSG:3031 this equals the longitude,
    but measuring it keeps the code correct if the projection changes.
    """
    lat = np.atleast_1d(np.asarray(lat, float))
    lon = np.atleast_1d(np.asarray(lon, float))
    x0, y0 = T3031.transform(lon, lat)
    x1, y1 = T3031.transform(lon, lat + step)
    return (-np.degrees(np.arctan2(x1 - x0, y1 - y0))) % 360.0

"""Shared Antarctic imagery basemap and scale-bar helpers.

Provides LIMA (Landsat Image Mosaic of Antarctica, RGB, coverage north of
~82.5 S) and MODIS MOA (grayscale, full continent) backgrounds for the
figure scripts, plus a projected scale bar. Mosaics come from the PGC
Quantarctica mirror; download once into ~/data/opr/basemap (override with
ANT_BASEMAP_DIR):

  base=https://data.pgc.umn.edu/gis/packages/quantarctica/Quantarctica3/SatelliteImagery
  curl -o ~/data/opr/basemap/LIMA_Mosaic.jp2  $base/LIMA/LIMA_Mosaic.jp2
  curl -o ~/data/opr/basemap/MODIS_Mosaic.tif $base/MODIS/MODIS_Mosaic.tif

Both mosaics are in EPSG:3031 (Antarctic Polar Stereographic, true scale
71 S), so callers should build their GeoAxes with `proj3031()` and pass
extents in projected meters. Everything degrades gracefully: if rasterio
or the mosaic files are missing, `add_imagery` returns False and callers
keep their coastline-only maps.
"""
import os

import numpy as np

try:
    import rasterio
    from rasterio.windows import from_bounds
    HAVE_RASTERIO = True
except Exception:
    HAVE_RASTERIO = False

import cartopy.crs as ccrs

BASEMAP_DIR = os.environ.get('ANT_BASEMAP_DIR',
                             os.path.expanduser('~/data/opr/basemap'))
LIMA_FN = os.path.join(BASEMAP_DIR, 'LIMA_Mosaic.jp2')
MODIS_FN = os.path.join(BASEMAP_DIR, 'MODIS_Mosaic.tif')

# Whole-continent EPSG:3031 extent (xmin, xmax, ymin, ymax) in meters
ANTARCTICA_EXTENT = (-3.1e6, 3.1e6, -2.7e6, 3.3e6)

# LIMA has no Landsat coverage poleward of this latitude; the mosaic
# fills the pole hole with flat MODIS-derived color
LIMA_SOUTH_LIMIT = -82.5


def proj3031():
    """Cartopy CRS matching the mosaics (EPSG:3031)."""
    return ccrs.SouthPolarStereo(true_scale_latitude=-71.0)


def points_extent(lats, lons, pad_frac=0.3, min_pad_m=5e3):
    """Padded EPSG:3031 extent (xmin, xmax, ymin, ymax) around points."""
    proj = proj3031()
    xyz = proj.transform_points(ccrs.PlateCarree(),
                                np.asarray(lons, dtype=float),
                                np.asarray(lats, dtype=float))
    x, y = xyz[:, 0], xyz[:, 1]
    ok = np.isfinite(x) & np.isfinite(y)
    if not np.any(ok):
        raise ValueError('points_extent: no finite lat/lon positions')
    x, y = x[ok], y[ok]
    pad_x = max(min_pad_m, pad_frac * np.ptp(x))
    pad_y = max(min_pad_m, pad_frac * np.ptp(y))
    return (x.min() - pad_x, x.max() + pad_x, y.min() - pad_y, y.max() + pad_y)


def _read_window(fn, extent, max_px):
    with rasterio.open(fn) as src:
        b = src.bounds
        xmin, xmax = max(extent[0], b.left), min(extent[1], b.right)
        ymin, ymax = max(extent[2], b.bottom), min(extent[3], b.top)
        if xmin >= xmax or ymin >= ymax:
            return None
        win = from_bounds(xmin, ymin, xmax, ymax, src.transform)
        scale = max(win.width / max_px, win.height / max_px, 1.0)
        shape = (src.count, max(1, round(win.height / scale)),
                 max(1, round(win.width / scale)))
        img = src.read(window=win, out_shape=shape, boundless=False)
    if img.size == 0 or not np.any(img):
        return None
    return np.moveaxis(img, 0, -1), (xmin, xmax, ymin, ymax)


def add_imagery(ax, extent, max_px=1400):
    """Draw LIMA (preferred) or MODIS MOA under a proj3031 GeoAxes.

    extent is (xmin, xmax, ymin, ymax) in EPSG:3031 meters. Returns True
    if imagery was drawn.
    """
    if not HAVE_RASTERIO:
        return False
    proj = proj3031()
    for fn, is_rgb in ((LIMA_FN, True), (MODIS_FN, False)):
        if not os.path.exists(fn):
            continue
        try:
            r = _read_window(fn, extent, max_px)
        except Exception:
            continue
        if r is None:
            continue
        img, ext = r
        mpl_extent = (ext[0], ext[1], ext[2], ext[3])
        if is_rgb and img.shape[-1] >= 3:
            rgb = img[..., :3]
            # Fall through to MOA inside LIMA's pole hole. The Quantarctica
            # JP2 fills it with a flat color statistically indistinguishable
            # from featureless Landsat snow (its per-channel std is JP2
            # compression noise, larger than a real WAIS Divide scene's), so
            # detect it geometrically: skip LIMA when every corner of the
            # window is poleward of the Landsat coverage limit. Also keep
            # the black-nodata check for mosaic builds without the fill.
            cx = np.array([ext[0], ext[1], ext[0], ext[1]])
            cy = np.array([ext[2], ext[2], ext[3], ext[3]])
            corner_lat = ccrs.PlateCarree().transform_points(proj, cx, cy)[:, 1]
            if np.max(corner_lat) < LIMA_SOUTH_LIMIT:
                continue
            if np.mean(np.all(rgb == 0, axis=-1)) > 0.5:
                continue
            ax.imshow(rgb, extent=mpl_extent, transform=proj, origin='upper',
                      zorder=0, interpolation='nearest')
        else:
            g = img[..., 0].astype(float)
            valid = g > 0
            if not np.any(valid):
                continue
            lo, hi = np.percentile(g[valid], [2, 98])
            g[~valid] = np.nan
            ax.imshow(g, extent=mpl_extent, transform=proj, origin='upper',
                      cmap='gray', vmin=lo, vmax=max(hi, lo + 1), zorder=0,
                      interpolation='nearest')
        return True
    return False


def _nice_length(target_m):
    """Round to 1, 2, or 5 times a power of ten."""
    exp = np.floor(np.log10(target_m))
    frac = target_m / 10**exp
    nice = 1 if frac < 1.5 else (2 if frac < 3.5 else 5)
    return nice * 10**exp


def add_scale_bar(ax, extent, frac=0.25):
    """Draw a horizontal scale bar on a proj3031 GeoAxes.

    Stereographic scale is exact at 71 S and distorted by <~3% over the
    survey latitudes; the bar length is nominal projected distance.
    """
    width = extent[1] - extent[0]
    height = extent[3] - extent[2]
    length = _nice_length(frac * width)
    x0 = extent[0] + 0.06 * width
    y0 = extent[2] + 0.07 * height
    proj = proj3031()
    ax.plot([x0, x0 + length], [y0, y0], color='white', lw=5,
            transform=proj, zorder=8, solid_capstyle='butt')
    ax.plot([x0, x0 + length], [y0, y0], color='black', lw=2.5,
            transform=proj, zorder=9, solid_capstyle='butt')
    label = f'{length/1000:g} km' if length >= 1000 else f'{length:g} m'
    ax.text(x0 + length / 2, y0 + 0.03 * height, label, ha='center',
            va='bottom', fontsize=8, transform=proj, zorder=9,
            path_effects=None,
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none',
                      alpha=0.7))

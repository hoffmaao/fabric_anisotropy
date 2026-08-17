"""Reduce every Ridge A polarimetric interferogram to movie-sized arrays.

Run on mem1, not locally. The shipped CSARP_polarimetric frames are ~690 MB
each and there are 36 of them, so the whole survey is ~25 GB - far too much
to pull over this link for something that only ever becomes movie frames.

This is the same trade extract_swath.py makes for the swath cubes: decimate
to display resolution and quantize to uint8, deliberately, because the
output is pixels. Phase is an angle on a cyclic colormap, so 8 bits over
[-pi, pi] is 1.4 deg a step - far finer than the colormap renders or the
eye resolves. Coherence only drives the value channel, where 8 bits is
likewise beyond what shows.

What is NOT quantized: the geometry. Latitude, longitude, the time axis and
the along-track distance stay float64, because the map marker and the
distance axis are read quantitatively.

Writes one npz per frame into out/, each ~1.5 MB, plus a manifest listing
the frames in survey order.

Usage: python3 reduce_ifg_movie.py <site_root> <out_dir>
"""
import glob
import os
import sys

import h5py
import numpy as np

SITE = sys.argv[1] if len(sys.argv) > 1 else (
    '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2')
OUT = sys.argv[2] if len(sys.argv) > 2 else 'ifg_movie'
MAX_COL = 1100          # display columns; the panel is ~1200 px wide
MAX_ROW = 750           # display rows
C_ICE = 299792458.0 / np.sqrt(3.171)


def geodesic_km(lat, lon):
    """Cumulative along-track distance [km] from lat/lon, local-flat."""
    la, lo = np.radians(lat), np.radians(lon)
    dx = np.diff(lo) * np.cos(0.5 * (la[1:] + la[:-1])) * 6371.0
    dy = np.diff(la) * 6371.0
    return np.concatenate([[0.0], np.cumsum(np.hypot(dx, dy))])


def reduce_frame(fn):
    with h5py.File(fn) as f:
        def rd(k):
            return np.array(f[k])
        ifg = rd('interferogram_mlook')
        coh = rd('interferogram_coherence')
        t = rd('Time').ravel()
        lat = rd('Latitude').ravel()
        lon = rd('Longitude').ravel()
        srf = rd('Surface').ravel()
    # h5py hands back MATLAB's column-major arrays transposed, and complex
    # comes through as a compound dtype rather than a native complex.
    if ifg.dtype.names is not None:
        ifg = ifg['real'] + 1j * ifg['imag']
    if coh.dtype.names is not None:
        coh = coh['real'] + 1j * coh['imag']
    if ifg.shape[0] == lat.size:
        ifg, coh = ifg.T, coh.T
    nt, nx = ifg.shape
    if t.size != nt:
        raise ValueError('%s: Time is %d for %d rows' % (fn, t.size, nt))

    # decimate by striding: these are already multilooked, so striding is
    # the honest reduction. Averaging complex phase across a stride would
    # smooth fringes that the figure exists to show.
    sr = max(1, int(np.ceil(nt / MAX_ROW)))
    sc = max(1, int(np.ceil(nx / MAX_COL)))
    ph = np.angle(ifg[::sr, ::sc])
    cm = np.abs(coh[::sr, ::sc])

    # NaN must be handled BEFORE the cast: numpy casts NaN to uint8 as an
    # arbitrary value (0 on some builds, 255 on others), so a gap would
    # render as either a hard black block or a saturated one rather than as
    # the absent data it is. Coherence NaN -> 0 puts the value channel at
    # black, which is how the wrapped-phase panel already shows unusable
    # data; phase NaN is then irrelevant because that pixel is black, but
    # it is zeroed anyway so the cast is defined.
    bad = ~np.isfinite(cm) | ~np.isfinite(ph)
    ph = np.where(bad, 0.0, ph)
    cm = np.where(bad, 0.0, cm)

    ph_u8 = np.clip((ph + np.pi) / (2 * np.pi) * 255.0, 0,
                    255).astype(np.uint8)
    cm_u8 = np.clip(cm * 255.0, 0, 255).astype(np.uint8)

    dist = geodesic_km(lat, lon)
    return dict(
        phase=ph_u8, coh=cm_u8,
        t=t[::sr].astype(np.float64),
        dist=dist[::sc].astype(np.float64),
        lat=lat[::sc].astype(np.float64), lon=lon[::sc].astype(np.float64),
        lat_full=lat.astype(np.float64), lon_full=lon.astype(np.float64),
        surface=np.float64(np.nanmedian(srf)),
        stride=np.array([sr, sc]))


def main():
    os.makedirs(OUT, exist_ok=True)
    fns = sorted(glob.glob(os.path.join(
        SITE, 'CSARP_polarimetric', '*', 'Data_*.mat')))
    if not fns:
        raise SystemExit('no CSARP_polarimetric frames under %s' % SITE)
    print('%d frames' % len(fns), flush=True)
    tags = []
    for i, fn in enumerate(fns):
        tag = os.path.basename(fn)[len('Data_'):-len('.mat')]
        out_fn = os.path.join(OUT, 'ifg_%s.npz' % tag)
        if os.path.exists(out_fn):
            print('  [%2d/%d] %s cached' % (i + 1, len(fns), tag), flush=True)
            tags.append(tag)
            continue
        try:
            d = reduce_frame(fn)
        except Exception as exc:                       # noqa: BLE001
            print('  [%2d/%d] %s FAILED %s' % (i + 1, len(fns), tag, exc),
                  flush=True)
            continue
        tmp = out_fn + '.tmp.npz'
        np.savez_compressed(tmp, **d)
        os.replace(tmp, out_fn)
        tags.append(tag)
        print('  [%2d/%d] %s %s -> %.2f MB' % (
            i + 1, len(fns), tag, 'x'.join(map(str, d['phase'].shape)),
            os.path.getsize(out_fn) / 1e6), flush=True)
    with open(os.path.join(OUT, 'manifest.txt'), 'w') as fh:
        fh.write('\n'.join(tags) + '\n')
    print('wrote %d frames to %s' % (len(tags), OUT))


if __name__ == '__main__':
    main()

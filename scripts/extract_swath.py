"""Compact swath-movie extract from a CSARP_music tomographic frame.

Streams Tomo.img (Nt x Nsv x Nx in MATLAB order), multilooks in
fast-time, quantizes dB to uint8 for the movie cube, and computes
angle-energy reductions (full record + depth bands) at full precision
BEFORE quantization. Writes a compressed npz next to nothing else.

Steering angles are stored in the npz already converted to DEGREES,
under the keys 'theta_deg' and 'theta_cal_deg'. Extracts written before
that contract carry the .mat vectors verbatim (i.e. RADIANS) under the
bare keys 'theta'/'theta_cal'; the key name, not the magnitude, is what
tells a reader which one it has.

Usage: python3 extract_swath.py <music_frame.mat> <out.npz> [ml]
"""
import sys

import h5py
import numpy as np

C_ICE = 1.68e8
BANDS = [(0., 500.), (500., 1000.), (1000., 1500.), (1500., None)]

fn, out = sys.argv[1], sys.argv[2]
ML = int(sys.argv[3]) if len(sys.argv) > 3 else 3


def angle_deg(f, name):
    """Steering angles in degrees, or None if absent.

    Products differ on whether the angle vectors sit inside Tomo or at
    the top level; they are stored in radians either way, so the
    conversion to degrees happens here, once, on the .mat side. The npz
    written below therefore carries degrees under the '*_deg' key names,
    and readers of those keys must not convert again.
    """
    node = f['Tomo'] if 'Tomo' in f and name in f['Tomo'] else f
    if name not in node:
        return None
    return np.degrees(np.squeeze(np.asarray(node[name])))


def linear_power(a):
    """Normalize a raw Tomo.img block to real linear power.

    h5py surfaces MATLAB complex as a compound (real, imag) dtype; some
    frames store img already complex, others already as power.
    """
    if a.dtype.names is not None:
        a = a['real'] + 1j * a['imag']
    if np.iscomplexobj(a):
        a = np.abs(a) ** 2
    return np.asarray(a, np.float64)


f = h5py.File(fn, 'r')
img = f['Tomo']['img']            # h5py order: (Nx, Nsv, Nt)
Nx, Nsv, Nt = img.shape
time = np.asarray(f['Time'])[0, :]
theta = angle_deg(f, 'theta')
if theta is None:
    sys.exit('%s has no theta in Tomo or at top level' % fn)
theta_cal = angle_deg(f, 'theta_cal')
if theta_cal is None:
    theta_cal = theta
surf = np.asarray(f['Surface'])[:, 0]
bot = np.asarray(f['Bottom'])[:, 0]
lat = np.asarray(f['Latitude'])[:, 0]
lon = np.asarray(f['Longitude'])[:, 0]
elev = np.asarray(f['Elevation'])[:, 0]
data2d = np.asarray(f['Data'])    # (Nx, Nt) combined radargram
print('cube Nx=%d Nsv=%d Nt=%d, ml=%d' % (Nx, Nsv, Nt, ML), flush=True)

Ntm = Nt // ML
ml_cube = np.empty((Nx, Nsv, Ntm), np.float32)
full_sum = np.zeros(Nsv)
full_n = 0
band_sum = np.zeros((len(BANDS), Nsv))
band_n = np.zeros(len(BANDS))

BLK = 64
for x0 in range(0, Nx, BLK):
    a = linear_power(img[x0:x0 + BLK])         # (blk, Nsv, Nt)
    ml_cube[x0:x0 + BLK] = \
        a[:, :, :Ntm * ML].reshape(a.shape[0], Nsv, Ntm, ML).mean(axis=3)
    full_sum += a.sum(axis=(0, 2))
    full_n += a.shape[0] * Nt
    for x in range(a.shape[0]):
        depth = (time - surf[x0 + x]) * C_ICE / 2
        for b, (d0, d1) in enumerate(BANDS):
            rows = depth >= d0 if d1 is None else \
                (depth >= d0) & (depth < d1)
            if rows.any():
                band_sum[b] += a[x, :, rows].sum(axis=0)
                band_n[b] += rows.sum()
    print('  %d/%d' % (min(x0 + BLK, Nx), Nx), flush=True)

db = 10 * np.log10(np.maximum(ml_cube, np.finfo(np.float32).tiny))
del ml_cube
vmin, vmax = np.percentile(db, [1, 99.9])
dbq = np.clip((db - vmin) / (vmax - vmin) * 255, 0, 255).astype(np.uint8)
del db

np.savez_compressed(
    out, img_dbq=dbq, db_min=vmin, db_max=vmax, ml=ML,
    time=time[:Ntm * ML].reshape(Ntm, ML).mean(axis=1),
    theta_deg=theta, theta_cal_deg=theta_cal,
    lat=lat, lon=lon, elev=elev, surface=surf, bottom=bot,
    data2d_db=10 * np.log10(np.maximum(
        data2d, np.finfo(np.float32).tiny)).astype(np.float32),
    full_mean=full_sum / full_n,
    band_mean=band_sum / np.maximum(band_n, 1)[:, None],
    band_edges=np.array([(d0, -1 if d1 is None else d1)
                         for d0, d1 in BANDS], float))
print('wrote', out, flush=True)

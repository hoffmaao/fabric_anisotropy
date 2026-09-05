"""Across-track swath movie and look-angle energy profile for GHOST2.

The rds GHOST2 swath radar rode the same 2024-25 traverse train as the
accum polarimetric radar: the tracks agree to 1-2 m, and segments such
as 20250112_01 exist in both rds/2024_Antarctica_GroundGHOST2 and
accum/2024_Antarctica_Ground2. The public portal carries only the 2D
radargrams and the swath-derived bottom DEMs; the tomographic image
volumes (Tomo.img: fast-time x steering angle x along-track) live in
the server-side CSARP_music3D product. Pull a frame from mem1, e.g.

  rsync mem1:/cresis/snfs1/dataproducts/ct_data/rds/\
2024_Antarctica_GroundGHOST2/CSARP_music3D/20250112_01/\
Data_20250112_01_002.mat ~/data/opr/rds/2024_Antarctica_GroundGHOST2/\
CSARP_music3D/20250112_01/

The movie scans through the LOOK ANGLES: each frame is a full
along-track radargram (distance x elevation, dark = strong return) at
one steering angle, with the angle read out beside a rotating beam
icon - the style of the englacial-scattering animation on slide 20 of
Christianson_RS_radar_swath_may2025 (Holschuh et al., in prep.).

Outputs (under out_dir, default ../../figs):
  <frame>_swath.mp4 (or .gif)   look-angle sweep movie
  <frame>_look_angle_energy.png mean energy vs look angle, full record
                                and per depth band

Rather than pulling a multi-GB frame, run extract_swath.py on mem1
(scratch swath/ dir) to reduce it to a compressed npz (multilooked
uint8 dB cube + full-precision angle-energy reductions) and rsync
that; this script accepts either the raw .mat or the .npz extract.
Current extracts store look angles in degrees under 'theta_deg';
extracts predating that contract hold radians under a bare 'theta',
which this script converts on load.

Usage: python ghost2_swath_movie.py <tomo_mat_or_npz> [out_dir]
       python ghost2_swath_movie.py --selftest [out_dir]
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scar_style  # noqa: E402

# ar-env carries the ffmpeg binary; apres carries the science stack
FFMPEG = '/opt/anaconda3/envs/ar-env/bin/ffmpeg'
if os.path.exists(FFMPEG):
    matplotlib.rcParams['animation.ffmpeg_path'] = FFMPEG

C_ICE = 1.68e8  # m/s, as in the other figure scripts

# Okabe-Ito, fixed assignment order for the depth-band series
BAND_COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']


def angle_deg(f, name):
    """Steering angles in degrees, or None if absent.

    Mirrors extract_swath.angle_deg (the two scripts run on different
    machines, so the helper is duplicated rather than imported):
    products differ on whether the angle vectors sit inside Tomo or at
    the top level, and they are stored in radians either way, so the
    conversion to degrees happens here, once, on the .mat side.
    """
    node = f['Tomo'] if 'Tomo' in f and name in f['Tomo'] else f
    if name not in node:
        return None
    return np.degrees(np.squeeze(np.asarray(node[name])))


def load_tomo(fn):
    """Load an OPR tomographic frame (v7.3 mat).

    Returns dict with img (Nt x Nsv x Nx, linear power), theta_deg
    (Nsv), time (Nt, s), surface (Nx, s twtt), lat/lon/elev (Nx),
    frame str. MATLAB v7.3 stores arrays with dimensions reversed, so
    a MATLAB (Nt x Nsv x Nx) cube reads as (Nx x Nsv x Nt) via h5py.
    """
    import h5py
    out = {}
    with h5py.File(fn, 'r') as f:
        if 'Tomo' in f:
            img = np.asarray(f['Tomo']['img'])
            theta_deg = angle_deg(f, 'theta')
        elif 'Data' in f and np.asarray(f['Data']).ndim == 3:
            img = np.asarray(f['Data'])
            theta_deg = angle_deg(f, 'Theta')
        else:
            raise ValueError(
                '%s has no Tomo/img or 3D Data variable: not a '
                'tomographic frame (music3D product)' % fn)
        if theta_deg is None:
            raise ValueError('%s has no steering-angle vector' % fn)
        img = np.transpose(img, (2, 1, 0))  # -> Nt x Nsv x Nx
        if img.dtype.names is not None:
            img = img['real'] + 1j * img['imag']
        if np.iscomplexobj(img):
            img = np.abs(img) ** 2
        out['img'] = img.astype(np.float64)
        out['theta_deg'] = theta_deg
        out['time'] = np.squeeze(np.asarray(f['Time']))
        out['lat'] = np.squeeze(np.asarray(f['Latitude']))
        out['lon'] = np.squeeze(np.asarray(f['Longitude']))
        out['elev'] = np.squeeze(np.asarray(f['Elevation']))
        out['surface'] = (np.squeeze(np.asarray(f['Surface']))
                          if 'Surface' in f else None)
    out['frame'] = os.path.basename(fn).replace('Data_', '').replace(
        '.mat', '')
    return out


def load_extract(fn):
    """Load the compact npz written by extract_swath.py on mem1.

    The look-angle key name carries the units: 'theta_deg' is the
    current contract (degrees). A bare 'theta' means an extract written
    before that contract, which copied the .mat vector verbatim and is
    therefore in RADIANS.
    """
    z = np.load(fn)
    dbq = z['img_dbq'].astype(np.float32)  # (Nx, Nsv, Ntm) uint8
    img_db = (dbq * ((z['db_max'] - z['db_min']) / 255.0) +
              z['db_min']).transpose(2, 1, 0)  # -> Nt x Nsv x Nx
    if 'theta_deg' in z:
        theta_deg = np.asarray(z['theta_deg'], float)
    else:
        theta_deg = np.degrees(np.asarray(z['theta'], float))
        print('note: %s is a legacy extract (bare "theta" key); '
              'converting look angles from radians to degrees' % fn)
    out = {'img_db': img_db, 'theta_deg': theta_deg,
           'time': z['time'], 'lat': z['lat'], 'lon': z['lon'],
           'elev': z['elev'], 'surface': z['surface'],
           'full_mean': z['full_mean'], 'band_mean': z['band_mean'],
           'band_edges': z['band_edges'],
           'frame': os.path.basename(fn).replace('swath_', '').replace(
               '.npz', '')}
    return out


def synthetic_tomo(nt=500, nsv=61, nx=240, seed=0):
    """Synthetic swath cube for the machinery self-test: specular
    surface and internal layers (narrow angular response), rough bed
    (broad response) with an off-nadir scatterer, noise floor."""
    rng = np.random.default_rng(seed)
    theta = np.radians(np.linspace(-45, 45, nsv))
    time = np.arange(nt) * 20e-9
    img = rng.exponential(1.0, (nt, nsv, nx))  # noise floor, linear
    surf_bin = 40 + (5 * np.sin(np.linspace(0, 3 * np.pi, nx))).astype(int)
    ang = np.degrees(theta)
    for x in range(nx):
        img[surf_bin[x], :, x] += 3e3 * np.exp(-(ang / 6.0) ** 2)
        for lay in range(80, 380, 45):  # internal layers, specular
            img[surf_bin[x] + lay, :, x] += (
                300 * np.exp(-(ang / 4.0) ** 2) *
                (1 + 0.3 * np.sin(lay + x / 30)))
        bed = surf_bin[x] + 400 + int(20 * np.cos(x / 40))
        if bed < nt:
            img[bed, :, x] += 100 * np.exp(-(ang / 25.0) ** 2)
    # an off-nadir englacial scatterer: bright only at theta ~ +20 deg
    img[250:280, np.abs(ang - 20) < 6, nx // 3:nx // 2] += 400
    return {'img': img, 'theta_deg': ang, 'time': time,
            'lat': np.linspace(-86.69, -86.5, nx),
            'lon': np.linspace(68.9, 70.1, nx),
            'elev': np.full(nx, 3185.0),
            'surface': time[surf_bin], 'frame': 'selftest'}


def along_track_km(lat, lon):
    R = 6371.0
    dlat = np.radians(np.diff(lat))
    dlon = np.radians(np.diff(lon)) * np.cos(np.radians(lat[:-1]))
    d = R * np.hypot(dlat, dlon)
    return np.concatenate([[0], np.cumsum(d)])


def cube_db(t):
    if 'img_db' in t:
        return t['img_db']
    return 10 * np.log10(np.maximum(t['img'], 1e-30))


def swath_movie(t, out_dir, fps=8, max_angle=60.0):
    img_db = cube_db(t)
    keep = np.abs(t['theta_deg']) <= max_angle
    img_db = img_db[:, keep, :]
    theta = t['theta_deg'][keep]
    nt, nsv, nx = img_db.shape
    # most of the record is empty range below the bed, so the black
    # level must sit well above the median or the ice column saturates
    vmin, vmax = np.percentile(img_db, [80, 99.9])
    x_km = along_track_km(t['lat'], t['lon'])
    us = t['time'] * 1e6

    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    axb = fig.add_subplot(gs[1])

    im = ax.imshow(img_db[:, 0, :], aspect='auto', cmap='gray_r',
                   vmin=vmin, vmax=vmax,
                   extent=[x_km[0], x_km[-1], us[-1], us[0]],
                   interpolation='nearest')
    ax.set_xlabel('Distance along profile (km)')
    ax.set_ylabel('TWTT (us)')
    ax.set_title(t['frame'], fontsize=10)

    # beam icon: antenna circle with a rotating boresight line
    axb.set_xlim(-1.6, 1.6)
    axb.set_ylim(-2.6, 0.9)
    axb.set_aspect('equal')
    axb.axis('off')
    label = axb.text(0, 0.7, '', ha='center', fontsize=12)
    circ = plt.Circle((0, 0), 0.28, fill=False, color='black', lw=1.6)
    axb.add_patch(circ)
    beam, = axb.plot([], [], color='#D55E00', lw=2.2)

    def update(k):
        th = theta[k]
        im.set_data(img_db[:, k, :])
        label.set_text('Look angle %.1f$^o$' % th)
        r = np.radians(th)
        beam.set_data([0.28 * np.sin(r), 2.2 * np.sin(r)],
                      [-0.28 * np.cos(r), -2.2 * np.cos(r)])
        return im, label, beam

    ani = animation.FuncAnimation(fig, update, frames=nsv, blit=False)
    base = os.path.join(out_dir, '%s_swath' % t['frame'])
    if animation.FFMpegWriter.isAvailable():
        out = base + '.mp4'
        ani.save(out, writer=animation.FFMpegWriter(fps=fps, bitrate=3000),
                 dpi=110)
    else:
        out = base + '.gif'
        ani.save(out, writer=animation.PillowWriter(fps=fps), dpi=80)
    plt.close(fig)
    return out


def look_angle_energy(t, out_dir, bands=((0, 500), (500, 1000),
                                         (1000, 1500), (1500, None))):
    """Mean energy vs look angle: full record plus depth bands below
    the surface pick (depth from TWTT at the ice speed C_ICE). An npz
    extract carries these reductions precomputed at full precision;
    otherwise they are computed here from the linear-power cube."""
    th = t['theta_deg']
    fig, ax = plt.subplots(figsize=(7, 5))

    if 'full_mean' in t:
        ax.plot(th, 10 * np.log10(np.maximum(t['full_mean'], 1e-30)),
                color='0.15', lw=2, label='Full record')
        for b, (d0, d1) in enumerate(t['band_edges']):
            bm = t['band_mean'][b]
            if not np.any(bm > 0):
                continue
            lab = '%d m +' % d0 if d1 < 0 else '%d-%d m' % (d0, d1)
            ax.plot(th, 10 * np.log10(np.maximum(bm, 1e-30)),
                    color=BAND_COLORS[b % len(BAND_COLORS)], lw=1.4,
                    label=lab)
        return _finish_angle_plot(fig, ax, t, out_dir)

    img = t['img']
    nt, nsv, nx = img.shape
    full = 10 * np.log10(np.maximum(img.mean(axis=(0, 2)), 1e-30))
    ax.plot(th, full, color='0.15', lw=2, label='Full record')

    if t['surface'] is not None:
        surf = np.atleast_1d(t['surface'])
        if surf.size == 1:
            surf = np.full(nx, float(surf))
        depth = (t['time'][:, None] - surf[None, :]) * C_ICE / 2
        for b, (d0, d1) in enumerate(bands):
            msk = depth >= d0 if d1 is None else \
                (depth >= d0) & (depth < d1)
            acc = np.zeros(nsv)
            n = 0
            for x in range(nx):
                rows = msk[:, x]
                if rows.any():
                    acc += img[rows, :, x].sum(axis=0)
                    n += rows.sum()
            if n == 0:
                continue
            lab = '%d m +' % d0 if d1 is None else '%d-%d m' % (d0, d1)
            ax.plot(th, 10 * np.log10(np.maximum(acc / n, 1e-30)),
                    color=BAND_COLORS[b % len(BAND_COLORS)], lw=1.4,
                    label=lab)

    return _finish_angle_plot(fig, ax, t, out_dir)


def _finish_angle_plot(fig, ax, t, out_dir):
    ax.set_xlabel('Look angle (deg, + right of track)')
    ax.set_ylabel('Mean power (dB, uncal)')
    ax.set_title('Across-track energy vs look angle (%s)' % t['frame'],
                 fontsize=11)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=9)
    out = os.path.join(out_dir, '%s_look_angle_energy.png' % t['frame'])
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return out


def main():
    args = [a for a in sys.argv[1:] if a != '--selftest']
    selftest = '--selftest' in sys.argv
    extra = args if selftest else args[1:]
    out_dir = scar_style.out_dir(extra, 0)
    os.makedirs(out_dir, exist_ok=True)

    if selftest:
        t = synthetic_tomo()
    elif not args:
        sys.exit(__doc__)
    elif args[0].endswith('.npz'):
        t = load_extract(args[0])
    else:
        t = load_tomo(args[0])

    shape = t['img_db'].shape if 'img_db' in t else t['img'].shape
    print('cube: %d fast-time x %d angles x %d rlines, theta %.1f..%.1f deg'
          % (shape + (t['theta_deg'][0], t['theta_deg'][-1])))
    print('wrote', look_angle_energy(t, out_dir))
    print('wrote', swath_movie(t, out_dir))


if __name__ == '__main__':
    main()

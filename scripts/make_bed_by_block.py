"""Bed depth per along-track section block, from the CSARP_layer picks.

Writes the `bed_by_block.json` that `scripts/figures/quadpol_depth_movie.py`
reads to retire a block from the movie once the bed enters its depth window.
The movie needs the bed PER BLOCK, not per site: the ice varies by 566-996 m
across Taylor Dome and by 200 m inside single frames, so a per-site depth cut
cannot express it, and without a bed a frame keeps drawing fabric colour at
depths where its ice has already ended.

The bed is derived, not read: CSARP_layer records two-way travel times, so

    depth = (twtt_bottom - twtt_surface) * c / (2 * sqrt(EPS_ICE))

with the SAME permittivity the pipeline's own depth axis uses (3.171, which
is ptt.constants eps_bar), because a bed on a different ice model than the
fabric it masks would cut the column in the wrong place. No firn correction,
for the same reason: the quad-pol depth axis carries none either, so the two
are consistent to within the ~10 m documented in docs/method.md.

THE TWO LAYERS ARE CHOSEN BY NAME, never by position. These files routinely
carry internal reflectors beyond the surface/bottom pair, and differencing
whatever happens to be layer 2 would give a systematically shallow bed - a
600 m "bed" in EastGRIP's 2700 m column - that no plausibility filter can
catch and that retires blocks from the movie hundreds of metres above the
real ice base. A file that exposes no identifiable bottom is skipped with a
line in the log instead: not masking is recoverable, masking wrongly is not.

MATCHING IS BY POSITION, PER FRAME. A section block's centre is matched to
the nearest pick within the SAME frame's layer file and within MAX_MATCH_M;
the trace axes differ (the qlook frames are culled before they are blocked),
so an index-for-index match would be wrong. A block with no pick near it,
or whose pick gives an implausible thickness, is written as null - the movie
draws those blocks unmasked rather than guessing.

FORMAT: one JSON object keyed by frame tag (`20250108_02_009`, matching
`quadpol_section_<tag>.mat`), each value a list of bed depths in metres below
the surface, one per section block IN BLOCK ORDER, with null where there is
no usable pick. The list length must equal that frame's block count; the
movie warns and skips masking for any frame where it does not, which is the
signal to rerun this script after a block size change.

Usage: python scripts/make_bed_by_block.py <site_root> [<site_root> ...]

`site_root` is a season directory holding `CSARP_layer/<day_seg>/`. Sections
are read from, and the JSON written to, `$SCAR_DATA` (default
`~/data/opr/scar`), beside the rest of the mirrored products.
"""
import glob
import json
import os
import sys

import h5py
import numpy as np
from scipy.io import loadmat

DATA = os.path.expanduser(os.environ.get('SCAR_DATA', '~/data/opr/scar'))
OUT_FN = os.path.join(DATA, 'bed_by_block.json')

C0 = 299792458.0
EPS_ICE = 3.171
C_ICE = C0 / np.sqrt(EPS_ICE)

# A pick further than this from a block centre is describing different ice.
# Blocks are 125 m long, so half a block is the natural scale. The radius
# costs nothing measurable: over 56 Taylor Dome blocks the nearest pick sits
# a median 0.7 m from the block centre, p95 1.5 m, max 1.7 m, so 62.5 m and
# 250 m both match 100% of blocks - and the tighter one removes the case
# where a block borrows its bed from two blocks away, which at Taylor Dome's
# 200 m of within-frame thickness variation would be the wrong ice.
MAX_MATCH_M = 62.5
# Thicknesses outside this range are a mispick, not ice: the shallowest
# sounding here is the ~200 m McMurdo shelf and the deepest column is
# EastGRIP's ~2700 m.
MIN_BED_M, MAX_BED_M = 40.0, 4500.0
R_E = 6371000.0


def _deref(f, v):
    if not isinstance(v, h5py.Dataset):
        return None
    try:
        is_ref = h5py.check_dtype(ref=v.dtype) is not None
    except (TypeError, AttributeError):
        is_ref = False
    if is_ref:
        return [_deref(f, f[r]) for r in np.asarray(v).ravel()]
    return np.array(v).T


def load_mat(fn):
    """One .mat as a plain dict, whichever of the two formats it is saved in."""
    try:
        return loadmat(fn, squeeze_me=True, struct_as_record=False)
    except NotImplementedError:
        with h5py.File(fn, 'r') as f:
            return {k: _deref(f, v) for k, v in f.items()
                    if not k.startswith('#')}


def _flat(v):
    return np.atleast_1d(np.asarray(v, float)).ravel()


def _as_name(v):
    """One layer name as a lowercase string, from either .mat flavour."""
    if isinstance(v, np.ndarray) and v.dtype.kind in 'iuf':
        v = ''.join(chr(int(c)) for c in v.ravel() if int(c) > 0)
    return str(v).strip().lower()


def _find_layer(names, ids, wanted_name, wanted_id):
    """Index of the layer identified as `wanted_name`, or None.

    By NAME first, exact then substring; by the OPR layer id (1 surface,
    2 bottom) only where the file carries one. Never by position - a file
    whose second layer is an internal reflector would otherwise yield a
    systematically shallow bed that the plausibility filter cannot catch
    and that masks the movie hundreds of metres above the real ice base.
    """
    for match in (lambda nm: nm == wanted_name,
                  lambda nm: wanted_name in nm):
        for i, nm in enumerate(names):
            if match(nm):
                return i
    for i, v in enumerate(ids):
        if v == wanted_id:
            return i
    return None


def layer_picks(d):
    """(lat, lon, bed depth) arrays from one CSARP_layer file.

    Handles both layer formats in use: the OPR layerdata one (`lat`, `lon`,
    `twtt`, `lyr_name`/`lyr_id`) and the legacy CReSIS one (`Latitude`,
    `Longitude`, `layerData` cells, each carrying its own `name`).
    """
    if 'twtt' in d and 'lat' in d:
        lat, lon = _flat(d['lat']), _flat(d['lon'])
        tw = d['twtt']
        if isinstance(tw, (list, tuple)):
            layers = [_flat(t) for t in tw]
        else:
            tw = np.atleast_2d(np.asarray(tw, float))
            if tw.shape[0] > tw.shape[1]:
                tw = tw.T
            layers = [tw[i] for i in range(tw.shape[0])]
        raw_names = d.get('lyr_name', d.get('name'))
        names = ([] if raw_names is None
                 else [_as_name(v) for v in np.atleast_1d(raw_names)])
        raw_ids = d.get('lyr_id', d.get('id'))
        ids = [] if raw_ids is None else list(_flat(raw_ids))
    elif 'layerData' in d and 'Latitude' in d:
        lat, lon = _flat(d['Latitude']), _flat(d['Longitude'])
        layers, names, ids = [], [], []
        for lay in np.atleast_1d(d['layerData']):
            val = np.atleast_1d(getattr(lay, 'value', None))
            if val.size == 0:
                continue
            # value{2}.data carries the pick; value{1} is the manual layer
            entry = val[-1] if val.size > 1 else val[0]
            data = getattr(entry, 'data', None)
            if data is None:
                continue
            layers.append(_flat(data))
            names.append(_as_name(getattr(lay, 'name', '')))
            ids.append(float(np.atleast_1d(
                np.asarray(getattr(lay, 'id', np.nan), float)).ravel()[0]))
    else:
        raise ValueError('no recognised layer fields')
    if len(layers) < 2:
        raise ValueError('fewer than two layers (need surface and bottom)')
    names = names[:len(layers)]
    ids = ids[:len(layers)]
    i_s = _find_layer(names, ids, 'surface', 1)
    i_b = _find_layer(names, ids, 'bottom', 2)
    if i_s is None or i_b is None:
        raise ValueError(
            'cannot identify %s among %d layers named %s; refusing to '
            'difference layers by position'
            % (' and '.join(w for w, i in (('surface', i_s), ('bottom', i_b))
                            if i is None),
               len(layers), names if names else '<unnamed>'))
    n = min(lat.size, lon.size, layers[i_s].size, layers[i_b].size)
    bed = (layers[i_b][:n] - layers[i_s][:n]) * C_ICE / 2.0
    bed[~np.isfinite(bed)] = np.nan
    bed[(bed < MIN_BED_M) | (bed > MAX_BED_M)] = np.nan
    return lat[:n], lon[:n], bed


def find_layer_file(roots, day_seg, frm):
    for root in roots:
        seg_dir = os.path.join(root, 'CSARP_layer', day_seg)
        for pat in ('Data_%s_%03d.mat', 'layer_%s_%03d.mat'):
            fn = os.path.join(seg_dir, pat % (day_seg, frm))
            if os.path.exists(fn):
                return fn
    return None


def nearest_bed(blat, blon, plat, plon, bed):
    """Bed at the pick nearest each block centre, or NaN beyond MAX_MATCH_M."""
    ok = np.isfinite(plat) & np.isfinite(plon) & np.isfinite(bed)
    if not ok.any():
        return np.full(blat.size, np.nan)
    plat, plon, bed = plat[ok], plon[ok], bed[ok]
    out = np.full(blat.size, np.nan)
    pph, plm = np.radians(plat), np.radians(plon)
    for j in range(blat.size):
        if not (np.isfinite(blat[j]) and np.isfinite(blon[j])):
            continue
        bph, blm = np.radians(blat[j]), np.radians(blon[j])
        a = (np.sin((pph - bph) / 2)**2
             + np.cos(bph) * np.cos(pph) * np.sin((plm - blm) / 2)**2)
        dist = 2 * R_E * np.arcsin(np.minimum(1.0, np.sqrt(a)))
        i = int(np.argmin(dist))
        if dist[i] <= MAX_MATCH_M:
            out[j] = bed[i]
    return out


def section_blocks(fn):
    """(block lat, block lon) for one quadpol_section_*.mat."""
    with h5py.File(fn) as f:
        r = f['res']
        blat = np.array(r['sec_lat']).ravel()
        blon = np.array(r['sec_lon']).ravel()
    return blat, blon


def main():
    roots = [os.path.expanduser(a) for a in sys.argv[1:]]
    if not roots:
        raise SystemExit('usage: python scripts/make_bed_by_block.py '
                         '<site_root> [<site_root> ...]')
    missing = [r for r in roots if not os.path.isdir(r)]
    if missing:
        raise SystemExit('no such site root: %s' % ', '.join(missing))

    sections = sorted(glob.glob(os.path.join(DATA, 'quadpol_section_*.mat')))
    if not sections:
        raise SystemExit('no quadpol_section_*.mat under %s' % DATA)

    beds, n_picked, n_frames = {}, 0, 0
    for fn in sections:
        tag = os.path.basename(fn)[len('quadpol_section_'):-len('.mat')]
        day_seg, frm_s = tag.rsplit('_', 1)
        if not frm_s.isdigit():
            continue                      # a non-default depth rerun, _z900
        blat, blon = section_blocks(fn)
        lfn = find_layer_file(roots, day_seg, int(frm_s))
        if lfn is None:
            print('  %s: no CSARP_layer file' % tag)
            continue
        try:
            plat, plon, bed = layer_picks(load_mat(lfn))
        except (ValueError, KeyError, OSError) as err:
            print('  %s: %s unusable (%s)' % (tag, os.path.basename(lfn), err))
            continue
        vals = nearest_bed(blat, blon, plat, plon, bed)
        beds[tag] = [None if not np.isfinite(v) else round(float(v), 1)
                     for v in vals]
        n_frames += 1
        n_picked += int(np.isfinite(vals).sum())
        print('  %s: %d of %d blocks picked (%.0f-%.0f m)'
              % (tag, int(np.isfinite(vals).sum()), vals.size,
                 np.nanmin(vals) if np.isfinite(vals).any() else np.nan,
                 np.nanmax(vals) if np.isfinite(vals).any() else np.nan))

    if not beds:
        raise SystemExit('no frame matched a CSARP_layer file; nothing written')
    # Write to a temp name in the same directory and rename into place, the
    # pattern run_quadpol_pipeline.m uses for its coreg cache: the movie is
    # fatal on a present-but-unreadable file, so a run killed mid-write would
    # otherwise leave a truncated JSON that breaks every later movie run.
    tmp_fn = OUT_FN + '.tmp'
    with open(tmp_fn, 'w') as fh:
        json.dump(beds, fh, indent=1, sort_keys=True)
    os.replace(tmp_fn, OUT_FN)
    print('wrote %s: %d frames, %d blocks with a bed'
          % (OUT_FN, n_frames, n_picked))


if __name__ == '__main__':
    main()

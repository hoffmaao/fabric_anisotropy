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
600 m "bed" in EastGRIP's 2700 m column - that retires blocks from the movie
hundreds of metres above the real ice base. Anything ambiguous or unnamed is
skipped with a line in the log instead: not masking is recoverable, masking
wrongly is not. Which layers a frame's bed was actually differenced from is
recorded per frame in the output, so a wrong binding is auditable after the
run rather than invisible.

BED-PICK PREFERENCE, most trustworthy first:

  1. the per-trace median of whichever of bottom_HH, bottom_VV, bottom_HV
     and bottom_VH the file carries,
  2. `bottom`,
  3. `bottom_mc`,
  4. a uniquely bottom-named layer, then OPR layer id 2.

bottom_mc RANKS LAST because it is measurably biased shallow, not because it
is a fallback in name only. Over all 9 frames of 2022_Antarctica_Ground and
2023_Antarctica_Ground carrying both, bottom_mc sits a median 44.6 m ABOVE
the polarimetric bed (per frame -5.3 to -76.4 m) while bottom_mc_bot sits
44.4 m below it: the pair brackets a basal zone rather than picking it, and
the polarimetric pick lands almost exactly midway. The four polarimetric
layers agree with each other to under a metre, which is what a real
interface looks like and what bottom_mc does not. Preferring bottom_mc would
grey out ~45 m of real ice at Eastwind and McMurdo - precisely the two
seasons where it is the only non-polarimetric option.

SURFACE PREFERENCE puts the EXACT `surface` first, which is where the two
sides differ. All 107 layer files across 2022_Antarctica_Ground,
2023_Antarctica_Ground and 2025_Antarctica_Ground2 carry a plain `surface`;
69 of them carry a `surface_dem` beside it, and that is a digital elevation
model rather than a radar pick, so it must never win. The per-channel tier
(surface_HH/VV/HV/VH) sits second as a defence for a future season that
names its surface the way it already names its bed - no season does today.
A `_dem` layer is barred from the fuzzy tier outright for the same reason.

THE PICKS THEMSELVES PASS THROUGH UNTOUCHED. There is no thickness
plausibility filter and no substitution of frame medians for suspect picks:
thin ice is real at these sites and a small pick is not evidence of a bad
one. At Eastwind the thin blocks correlate -0.62/-0.64 with lat/lon toward
the terminus, and at McMurdo the thin blocks are non-overlapping in
longitude with the thick ones (+0.80 lat, +0.91 lon). The one validity check
is on the SIGN: a bottom at or above the surface is not thin ice but a
malformed record (an older layerData storing 0 for unpicked traces), and it
becomes null with a warning naming the frame and the count. It has to,
because a negative depth is finite and would leave `hi <= bed` false in
every window - greying that block for the whole movie, where null leaves it
drawn. Measured across the three seasons this fires on nothing: 0 negative
and 0 zero values in 66010 finite picks.

MATCHING IS BY POSITION, PER FRAME. A section block's centre is matched to
the nearest pick within the SAME frame's layer file and within MAX_MATCH_M;
the trace axes differ (the qlook frames are culled before they are blocked),
so an index-for-index match would be wrong. A block with no pick near it is
written as null - the movie draws those blocks unmasked rather than guessing.

FORMAT: one JSON object keyed by frame tag (`20250108_02_009`, matching
`quadpol_section_<tag>.mat`), each value a list of bed depths in metres below
the surface, one per section block IN BLOCK ORDER, with null where there is
no usable pick. The list length must equal that frame's block count; the
movie warns and skips masking for any frame where it does not, which is the
signal to rerun this script after a block size change. One reserved key,
`_layers`, maps each frame tag to the surface and bottom layer its bed came
from (the bottom reads `median(bottom_hh,bottom_vv,...)` where tier 1 won);
it cannot collide with a tag, and the movie never looks it up.

Usage: python scripts/make_bed_by_block.py <site_root> [<site_root> ...]

`site_root` is a season directory holding `CSARP_layer/<day_seg>/`. Sections
are read from, and the JSON written to, `$SCAR_DATA` (default
`~/data/opr/scar`), beside the rest of the mirrored products.
"""
import glob
import json
import os
import sys
import warnings

import h5py
import numpy as np
from scipy.io import loadmat
from scipy.io.matlab import MatReadError

DATA = os.path.expanduser(os.environ.get('SCAR_DATA', '~/data/opr/scar'))
OUT_FN = os.path.join(DATA, 'bed_by_block.json')
# Reserved key in the output, holding the layer names bound per frame. It
# cannot collide with a frame tag, which is always <day_seg>_<frm>.
LAYERS_KEY = '_layers'

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
R_E = 6371000.0
# The bed layers that agree with each other to under a metre, combined per
# trace rather than ranked among themselves. See the preference order in the
# module docstring for why bottom_mc is not one of them.
POLAR_BOTTOMS = ('bottom_hh', 'bottom_vv', 'bottom_hv', 'bottom_vh')
POLAR_SURFACES = ('surface_hh', 'surface_vv', 'surface_hv', 'surface_vh')
# Substrings that mark a layer as something other than a radar pick. A
# `surface_dem` is a digital elevation model: it loses to an exact `surface`
# on tier order anyway, but it must not be admitted to the fuzzy tier either,
# where in a file lacking a plain `surface` it would be the sole candidate
# and would shift every bed in the frame by the DEM's own offset.
NOT_A_PICK = ('dem',)


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
    """One .mat as a plain dict, whichever format it is saved in."""
    try:
        return loadmat(fn, squeeze_me=True, struct_as_record=False)
    except NotImplementedError:
        with h5py.File(fn, 'r') as f:
            return {k: _deref(f, v) for k, v in f.items()
                    if not k.startswith('#')}


def _flat(v):
    return np.atleast_1d(np.asarray(v, float)).ravel()


class LayerFormatError(Exception):
    """This file was read, but its layers cannot be identified."""


def _as_name(v):
    """One layer name as a lowercase string, from either .mat flavour."""
    if isinstance(v, np.ndarray) and v.dtype.kind in 'iuf':
        v = ''.join(chr(int(c)) for c in v.ravel() if int(c) > 0)
    return str(v).strip().lower()


def _name_list(raw):
    """Layer names as lowercase strings, decoded ONE ELEMENT AT A TIME.

    A MATLAB cell of char reaches here from the v7.3 path as a LIST of
    per-name uint16 code arrays, ragged because 'surface' is 7 codes and
    'bottom' is 6. Stacking that - np.atleast_1d, np.asarray - raises on any
    file whose layer names differ in length, i.e. every real one, so each
    element is handed to _as_name on its own instead.
    """
    if raw is None:
        return []
    if isinstance(raw, (list, tuple)):
        return [_as_name(v) for v in raw]
    if isinstance(raw, np.ndarray):
        if raw.dtype.kind in 'iuf':
            # char codes: a flat run is ONE name, a 2-D block one per row
            return ([_as_name(raw)] if raw.ndim <= 1
                    else [_as_name(row) for row in raw])
        return [_as_name(v) for v in raw.ravel()]
    return [_as_name(raw)]


def _scalar(v, default=np.nan):
    """First element of `v` as a float, or `default` where there is none."""
    a = np.atleast_1d(np.asarray(v, float)).ravel()
    return float(a[0]) if a.size else default


def _label(names, i, fallback):
    """Name of layer `i`, or `fallback` where the file did not name it.

    `i` may have come from the id tier, which is bounded by len(ids) and has
    no length relation to `names` - a partially named file resolves a bottom
    off an id it never named.
    """
    return names[i] if i < len(names) and names[i] else fallback


def _exact(names, want):
    return [i for i, nm in enumerate(names) if nm == want]


def _channels(names, wanted):
    return [i for i, nm in enumerate(names) if nm in wanted]


def _contains(names, want):
    return [i for i, nm in enumerate(names)
            if want in nm and not any(x in nm for x in NOT_A_PICK)]


def _by_id(ids, want):
    return [i for i, v in enumerate(ids) if v == want]


def _resolve(layers, names, tiers):
    """(twtt, label) for one layer role, or (None, None).

    `tiers` is (candidate indices, combine, fallback label) in preference
    order, and encodes the one rule this module rests on: the FIRST tier
    offering any candidate decides. Exactly one binds; more than one is
    unidentifiable and stops the walk, because letting a later tier rescue an
    ambiguous earlier one is position dressed up as preference, and taking
    the first of several is the same thing. `combine` takes the per-trace
    median of the tier instead - the per-channel picks agree to under a
    metre, so their median is any one of them, robust to whichever the file
    omits.
    """
    for cand, combine, fallback in tiers:
        if not cand:
            continue
        if combine:
            n = min(layers[i].size for i in cand)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                tw = np.nanmedian(np.vstack([layers[i][:n] for i in cand]),
                                  axis=0)
            return tw, 'median(%s)' % ','.join(names[i] for i in cand)
        if len(cand) > 1:
            return None, None
        return layers[cand[0]], _label(names, cand[0], fallback)
    return None, None


def _surface_layer(layers, names, ids):
    """(surface twtt, label) - EXACT `surface` outranks the per-channel tier.

    The reverse of the bed order, deliberately: every layer file in these
    seasons carries a plain `surface`, and 69 of 107 carry a `surface_dem`
    beside it that is a DEM rather than a radar pick. See the module
    docstring.
    """
    return _resolve(layers, names, [
        (_exact(names, 'surface'), False, 'surface'),
        (_channels(names, POLAR_SURFACES), True, None),
        (_contains(names, 'surface'), False, 'surface'),
        (_by_id(ids, 1), False, 'id 1')])


def _bottom_layer(layers, names, ids):
    """(bed twtt, label) under the preference order in the module docstring."""
    return _resolve(layers, names, [
        (_channels(names, POLAR_BOTTOMS), True, None),
        (_exact(names, 'bottom'), False, 'bottom'),
        (_exact(names, 'bottom_mc'), False, 'bottom_mc'),
        (_contains(names, 'bottom'), False, 'bottom'),
        (_by_id(ids, 2), False, 'id 2')])


def layer_picks(d):
    """(lat, lon, bed depth, bound layer names, malformed count) from one file.

    Handles both layer formats in use: the OPR layerdata one (`lat`, `lon`,
    `twtt`, `lyr_name`/`lyr_id`) and the legacy CReSIS one (`Latitude`,
    `Longitude`, `layerData` cells, each carrying its own `name`).

    Raises LayerFormatError - and nothing else - when the file is readable
    but its surface and bottom cannot be named. Any other exception is a
    defect in this reader and must escape rather than be logged as a frame
    with no picks.
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
        names = _name_list(d.get('lyr_name', d.get('name')))
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
            ids.append(_scalar(getattr(lay, 'id', np.nan)))
    else:
        raise LayerFormatError('no recognised layer fields')
    if len(layers) < 2:
        raise LayerFormatError(
            'holds %d layer(s); needs a surface and a bottom' % len(layers))
    names = names[:len(layers)]
    ids = ids[:len(layers)]
    tw_s, lbl_s = _surface_layer(layers, names, ids)
    tw_b, lbl_b = _bottom_layer(layers, names, ids)
    if tw_s is None or tw_b is None:
        raise LayerFormatError(
            'cannot unambiguously identify %s among %d layers named %s; '
            'refusing to difference layers by position'
            % (' and '.join(w for w, t in (('surface', tw_s), ('bottom', tw_b))
                            if t is None),
               len(layers), names if names else '<unnamed>'))
    n = min(lat.size, lon.size, tw_s.size, tw_b.size)
    bed = (tw_b[:n] - tw_s[:n]) * C_ICE / 2.0
    bed[~np.isfinite(bed)] = np.nan
    bad = bed <= 0
    n_bad = int(np.count_nonzero(bad))
    bed[bad] = np.nan
    bound = {'surface': lbl_s, 'bottom': lbl_b}
    return lat[:n], lon[:n], bed, bound, n_bad


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

    beds, bound_by_tag, n_picked, n_frames = {}, {}, 0, 0
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
            d = load_mat(lfn)
        except (OSError, ValueError, MatReadError) as err:
            print('  %s: %s unreadable (%s)'
                  % (tag, os.path.basename(lfn), err))
            continue
        try:
            plat, plon, bed, bound, n_bad = layer_picks(d)
        except LayerFormatError as err:
            print('  %s: %s unusable (%s)' % (tag, os.path.basename(lfn), err))
            continue
        if n_bad:
            print('  WARNING: %s: %d pick(s) put the bottom at or above the '
                  'surface; dropped as malformed, not as thin ice'
                  % (tag, n_bad))
        vals = nearest_bed(blat, blon, plat, plon, bed)
        beds[tag] = [None if not np.isfinite(v) else round(float(v), 1)
                     for v in vals]
        bound_by_tag[tag] = bound
        n_frames += 1
        n_picked += int(np.isfinite(vals).sum())
        print('  %s: %d of %d blocks picked (%.0f-%.0f m) [surface=%s '
              'bottom=%s]'
              % (tag, int(np.isfinite(vals).sum()), vals.size,
                 np.nanmin(vals) if np.isfinite(vals).any() else np.nan,
                 np.nanmax(vals) if np.isfinite(vals).any() else np.nan,
                 bound['surface'], bound['bottom']))

    if not beds:
        raise SystemExit(
            'no frame matched a CSARP_layer file; nothing written')
    # The layers each frame's bed was differenced from, so a shallow bed off
    # an internal reflector is auditable in the product rather than only in
    # a run log that nobody kept.
    out = dict(beds)
    out[LAYERS_KEY] = bound_by_tag
    # Write to a temp name in the same directory and rename into place, the
    # pattern run_quadpol_pipeline.m uses for its coreg cache: the movie is
    # fatal on a present-but-unreadable file, so a run killed mid-write would
    # otherwise leave a truncated JSON that breaks every later movie run.
    tmp_fn = OUT_FN + '.tmp'
    with open(tmp_fn, 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    os.replace(tmp_fn, OUT_FN)
    print('wrote %s: %d frames, %d blocks with a bed'
          % (OUT_FN, n_frames, n_picked))


if __name__ == '__main__':
    main()

function [z_bed, info] = bed_from_layer(site_root, day_seg, frm, la, lo, surf_t)
%BED_FROM_LAYER Per-trace bed depth of a frame from the season's CSARP_layer picks.
%
% [z_bed, info] = bed_from_layer(site_root, day_seg, frm, la, lo, surf_t)
%
% z_bed is [1 x Nx] metres on the PIPELINE'S depth axis - the bed pick's
% two-way time less the frame surface time surf_t the pipeline zeroes its
% depth at, at solid-ice velocity (eps 3.171, the same ice model as the
% depth axis and as scripts/make_bed_by_block.py) - and NaN for a trace
% with no pick within MAX_MATCH_M (62.5 m) of it, unless it lies in a gap of
% at most MAX_GAP_M (1000 m) along track between matched traces. Traces are
% matched to picks by POSITION, nearest pick wins, because the layer product
% is on its own trace axis (a few hundred picks along a frame of thousands
% of traces).
%
% A BED IS CONTINUOUS; A PICKER'S GAP IS NOT A HOLE IN IT. A run of
% unmatched traces between two matched ones takes the bed interpolated
% linearly in along-track distance from them, when the gap spans at most
% MAX_GAP_M. Left NaN, those traces stay unmasked, and where they are more
% than half a pass the pass's median bed is infinite (ptt.quadpolFrameTheta
% counts an unpicked trace as ice), so a frame whose picks stop for a
% stretch would let that stretch's coherent sub-bed return vote on the
% held axis again. The cap bounds how far a straight line is trusted
% across unseen bed: longer gaps, and the ends of the frame, stay NaN, and
% info.n_filled counts what was filled.
%
% THE READER TAKES THE OPR LAYERDATA FORMAT ONLY - the frame file's twtt
% rows, named by the segment's layer_<seg>.mat - which every season in use
% ships. THE LAYERS ARE BOUND BY NAME, never by row. Names live once per
% segment in that catalogue (lyr_name; row i of the frame file's twtt is
% entry i), and layer ORDER varies by season, so a positional bind
% differences a different quantity per site: it once read a 224 m shelf as
% 1 m. With the catalogue in hand the ROWS of twtt are its layers: a file
% whose row count disagrees with the catalogue is refused, never transposed
% until the counts agree (that is a positional bind wearing a name). The bed
% preference, most trustworthy first, is the one make_bed_by_block.py
% established (bottom_mc runs a median 44.6 m shallow of the polarimetric
% picks and ranks last): the per-trace median of the bottom_HH/VV/HV/VH
% picks the file carries; then `bottom`; then `bottom_mc`; then a uniquely
% bottom-named layer that is not a DEM (a model, not a radar pick - the
% producer's NOT_A_PICK). A frame whose layers cannot be named, whose file is
% in another layout, or with no layer file at all, returns all NaN and says
% so in info - not masking is recoverable, masking ice away is not.
%
% WHY THE PIPELINE'S SURFACE AND NOT THE LAYER FILE'S. The bed has to land
% on the axis the fabric is reported on, and that axis is zeroed at the
% frame's median product Surface (or, in qlook mode, at the picked leading
% edge). On a ground survey the per-trace surface pick differs from that
% median by well under the 20 m margin the callers apply, so the two
% conventions agree to within it; the json producer keeps the per-pick
% difference because it has no pipeline surface to refer to.
%
% Inputs
%   site_root  season product root holding CSARP_layer/<day_seg>/
%   day_seg, frm
%   la, lo     [Nx] per-trace position of the frame being inverted, after
%              any cull (the axis z_bed comes back on)
%   surf_t     the pipeline's surface two-way time [s]
%
% Output info: source ('' when no bed, including when no trace lies within
%   MAX_MATCH_M of any pick), file, n_picks, n_matched, n_filled (traces
%   given an interpolated bed across a gap), n_bad (non-positive picks
%   dropped), reason (why there is no bed; '' when there is one)
%
% See also ptt.maskBelowBed, ptt.quadpolFrameTheta.

MAX_MATCH_M = 62.5;      % the bed producer's radius, scripts/make_bed_by_block.py
MAX_GAP_M = 1000;        % longest gap along track bridged by interpolation
EPS_ICE = 3.171;         % the pipeline's depth axis, ptt.constants eps_bar
C_ICE = 299792458 / sqrt(EPS_ICE);
la = la(:).'; lo = lo(:).';
Nx = numel(la);
z_bed = nan(1, Nx);
info = struct('source', '', 'file', '', 'n_picks', 0, 'n_matched', 0, ...
  'n_filled', 0, 'n_bad', 0, 'reason', '');

seg_dir = fullfile(site_root, 'CSARP_layer', day_seg);
fn = '';
for pat = {'Data_%s_%03d.mat', 'layer_%s_%03d.mat'}
  cand = fullfile(seg_dir, sprintf(pat{1}, day_seg, frm));
  if exist(cand, 'file') == 2, fn = cand; break; end
end
if isempty(fn)
  info.reason = sprintf('no CSARP_layer file for %s_%03d under %s', day_seg, frm, site_root);
  return
end
info.file = fn;
D = load(fn);

if ~(isfield(D, 'twtt') && isfield(D, 'lat') && isfield(D, 'lon'))
  info.reason = sprintf('%s is not OPR layerdata (no twtt/lat/lon; fields: %s)', ...
    fn, strjoin(fieldnames(D).', ', '));
  return
end
if ~(isnumeric(D.twtt) && isnumeric(D.lat) && isnumeric(D.lon))
  info.reason = sprintf('%s is not OPR layerdata (twtt/lat/lon are %s/%s/%s, not numeric)', ...
    fn, class(D.twtt), class(D.lat), class(D.lon));
  return
end
plat = double(D.lat(:).'); plon = double(D.lon(:).');
tw = double(D.twtt);
names = {};
cat_fn = fullfile(seg_dir, sprintf('layer_%s.mat', day_seg));
if exist(cat_fn, 'file') == 2
  L = load(cat_fn);
  if isfield(L, 'lyr_name'), names = cellstr(L.lyr_name(:).'); end
end
if isempty(names)
  info.reason = sprintf('%s carries no layer names (no layer_%s.mat catalogue)', fn, day_seg);
  return
end
if size(tw, 1) ~= numel(names)
  info.reason = sprintf('%d layer names for %d pick rows in %s', numel(names), size(tw, 1), fn);
  return
end
names = lower(strtrim(names));

[tw_b, src] = H_bottom(tw, names);
if isempty(tw_b)
  info.reason = sprintf('no bottom layer among {%s} in %s', strjoin(names, ', '), fn);
  return
end
info.source = src;
n = min([numel(plat), numel(plon), numel(tw_b)]);
plat = plat(1:n); plon = plon(1:n); tw_b = tw_b(1:n);
bed = (tw_b - surf_t) * C_ICE / 2;
bad = ~(bed > 0);                        % non-positive or NaN picks do not vote
info.n_bad = nnz(isfinite(bed) & bed <= 0);
ok = isfinite(plat) & isfinite(plon) & ~bad;
info.n_picks = nnz(ok);
if ~any(ok)
  info.reason = sprintf('%s: every %s pick is unpositioned or malformed', fn, src);
  return
end
plat = plat(ok); plon = plon(ok); bed = bed(ok);

% nearest pick per trace, haversine, within MAX_MATCH_M
R_E = 6371000;
tl = deg2rad(la); to = deg2rad(lo);
pl = deg2rad(plat); po = deg2rad(plon);
best = inf(1, Nx); ibest = zeros(1, Nx);
for i = 1:numel(pl)
  a = sin((tl - pl(i))/2).^2 + cos(tl) .* cos(pl(i)) .* sin((to - po(i))/2).^2;
  d = 2 * R_E * asin(min(1, sqrt(a)));
  hit = d < best;
  best(hit) = d(hit); ibest(hit) = i;
end
m = isfinite(best) & best <= MAX_MATCH_M & ibest > 0;
z_bed(m) = bed(ibest(m));
info.n_matched = nnz(m);
if info.n_matched == 0
  info.reason = sprintf('%d %s picks, none within %.1f m of a trace', info.n_picks, src, MAX_MATCH_M);
  info.source = '';
  return
end

% gaps between matched traces, bridged linearly in along-track distance
% (see A BED IS CONTINUOUS); a gap through an unpositioned trace has no
% length and is left alone
a = sin(diff(tl)/2).^2 + cos(tl(1:end-1)) .* cos(tl(2:end)) .* sin(diff(to)/2).^2;
step = 2 * R_E * asin(min(1, sqrt(a)));
im = find(m);
for q = 1:numel(im) - 1
  i0 = im(q); i1 = im(q+1);
  if i1 - i0 < 2, continue; end
  sg = [0, cumsum(step(i0:i1-1))];
  if ~(sg(end) > 0 && sg(end) <= MAX_GAP_M), continue; end
  z_bed(i0+1:i1-1) = z_bed(i0) + (z_bed(i1) - z_bed(i0)) * sg(2:end-1) / sg(end);
  info.n_filled = info.n_filled + (i1 - i0 - 1);
end
end

function [tw, src] = H_bottom(tw_all, names)
% the bed preference of make_bed_by_block.py, in its order
tw = []; src = '';
pol = {'bottom_hh', 'bottom_vv', 'bottom_hv', 'bottom_vh'};
ip = find(ismember(names, pol));
if ~isempty(ip)
  tw = median(tw_all(ip, :), 1, 'omitnan');
  src = sprintf('median(%s)', strjoin(names(ip), ','));
  return
end
for want = {'bottom', 'bottom_mc'}
  i = find(strcmp(names, want{1}));
  if isscalar(i), tw = tw_all(i, :); src = want{1}; return; end
end
i = find(contains(names, 'bottom') & ~contains(names, 'dem'));
if isscalar(i), tw = tw_all(i, :); src = names{i}; end
end

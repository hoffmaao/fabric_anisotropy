%RUN_COREG_FRAME Put one frame's four channels on one range axis: one delay per channel, fitted in the firn.
%
%   matlab -batch "addpath('<code>/opr_fabric/server'); day_seg='20250111_01'; frm=3; run_coreg_frame"
%
% (coreg_one.sh runs this for one frame with a lock and a log, and
% coreg_batch.sh for a whole survey.) Step 1 of the node method:
%
%   1. read the range window from the season's polarimetric product
%      (CSARP_polarimetric, or CSARP_polarimetric_unwrap where a season
%      shipped only that) and load the four SAR-focused channels
%      CSARP_standardphase_{HH,VV,HV,VH} over it
%   2. check that the HH channel is the image the product was built from
%   3. depth axis from the product's surface pick, kept down to Z_MAX
%      (1500 m, the one depth every cache reaches)
%   4. fit ONE range delay per channel over the firn (coreg_firn_m,
%      20-150 m) and shift the whole channel by it (ptt.coregisterChannels)
%   5. write the coregistered channels to
%      <work>/stages/quadpol_cf/coreg_cache/creg_<tag>.mat
%
% WHY ONE DELAY AND NOT A FIELD. The channels differ by a sub-sample
% instrument delay, which is one number per channel. A coregistration that
% is allowed to vary with depth also follows the delay between the slow and
% fast echoes of birefringent ice and removes it, and that delay is what
% run_node_frame measures. So the delay is fitted where the ice is nearly
% isotropic, the firn, and applied everywhere.
%
% Variables set before the call:
%   day_seg, frm   the frame (required)
%   site_root      the season's data root (default Ridge A,
%                  /cresis/nvme/opr_data/accum/2024_Antarctica_Ground2)
%   coreg_firn_m   the firn window the delays are fitted over ([20 150])
%   out_dir        (<work>/stages/quadpol_cf)
%
% The cache holds hh, vv, hv, vh (single, [Nt x Nx], rows down in range),
% z [Nt x 1] depth in metres, and what made it: CO (the coregistration
% settings), r0 and r1 (the range window in samples) and Nx. A cache that
% already exists with the same settings is left alone; one made with other
% settings is rebuilt. It is written to a temporary name and renamed, so a
% job killed while writing never leaves a truncated cache.
%
% The season must have a polarimetric product for the frame: the range
% window and the surface come from it, and quadpol_frame_blocks takes the
% trace positions from it.
if ~exist('day_seg', 'var') || isempty(day_seg), error('run_coreg_frame:frame', 'set day_seg'); end
if ~exist('frm', 'var') || isempty(frm), error('run_coreg_frame:frame', 'set frm'); end
% char, not string: shell escaping around matlab -batch often forces
% day_seg="...", and a string scalar does not concatenate into a file name
day_seg = char(day_seg);
if ~exist('site_root', 'var') || isempty(site_root)
  site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
end
site_root = char(site_root);
if ~exist('coreg_firn_m', 'var') || isempty(coreg_firn_m), coreg_firn_m = [20 150]; end
% Work root and repo root are DERIVED from this file's location, see
% fabric_paths. The /cresis site_root stays absolute: it is a real mount
% point, not part of the work tree.
[fabric_code, fabric_work] = fabric_paths();
addpath(fabric_code);
if ~exist('out_dir', 'var') || isempty(out_dir)
  out_dir = fullfile(fabric_work, 'stages', 'quadpol_cf');
end
CHAN = {'hh', 'vv', 'hv', 'vh'};
NRW = 101;                 % range samples the reported coherence is averaged over
Z_BAND = [200 1200];       % depths the reported coherence is read over [m]
Z_MAX = 1500;              % deepest sample kept [m]; run_node_frame reads no deeper
tag = sprintf('%s_%03d', day_seg, frm);
name = sprintf('Data_%s.mat', tag);
t_all = tic;

%% 1. the range window, from the polarimetric product
pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
if exist(pol_fn, 'file') ~= 2
  % the 2022 and 2023 seasons shipped the product only in its
  % SNAPHU-unwrapped form, which carries the same window and surface
  pol_fn = fullfile(site_root, 'CSARP_polarimetric_unwrap', day_seg, name);
end
if exist(pol_fn, 'file') ~= 2
  error('run_coreg_frame:noProduct', ...
    ['%s has no polarimetric product (looked for %s in CSARP_polarimetric and ' ...
     'CSARP_polarimetric_unwrap under %s); the range window and the surface come from it'], ...
    tag, name, site_root);
end
w = {whos('-file', pol_fn).name};
sel = intersect({'param_polarimetric', 'Surface', 'ref'}, w);
P = load(pol_fn, sel{:});
pp = P.param_polarimetric;
if isfield(pp, 'polarimetric'), pp = pp.polarimetric; end
r0 = pp.min_rbin; r1 = pp.max_rbin;
fprintf('=== %s ===\n', tag);
% FIT THE WINDOW TO THE FRAME THAT ACTUALLY EXISTS. The thin-ice seasons
% (Eastwind, McMurdo: 200-300 m of ice) record a max_rbin from a deeper
% configuration than their own data - e.g. 20240202_01_001 asks for 11000
% samples and holds 4151. Clamping keeps every sample that exists instead
% of discarding a 19 km line over a bookkeeping mismatch; only a window
% starting past the end is fatal. The channel sizes are taken from the file
% headers rather than by loading, so this costs nothing where it changes
% nothing.
navail = inf;
for k = 1:4
  fnk = fullfile(site_root, sprintf('CSARP_standardphase_%s', upper(CHAN{k})), day_seg, name);
  if exist(fnk, 'file') ~= 2
    error('run_coreg_frame:missing', 'missing %s channel: %s', upper(CHAN{k}), fnk);
  end
  wi = whos('-file', fnk, 'Data');
  if ~isempty(wi) && numel(wi(1).size) >= 2
    navail = min(navail, wi(1).size(1));
  end
end
if isfinite(navail) && r0 >= navail
  error('run_coreg_frame:window', 'window starts at %d but only %d samples exist', r0, navail);
end
if isfinite(navail) && r1 > navail
  warning('run_coreg_frame:windowClamped', ...
    'recorded window %d:%d exceeds the %d samples present; clamped to %d:%d', r0, r1, navail, r0, navail);
  r1 = navail;
end
fprintf('window rbin %d:%d\n', r0, r1);

%% 2. the four channels over that window, and the source check
S = struct();
for k = 1:4
  fn = fullfile(site_root, sprintf('CSARP_standardphase_%s', upper(CHAN{k})), day_seg, name);
  q = load(fn, 'Data', 'Time');
  if r1 > size(q.Data, 1)
    error('run_coreg_frame:window', '%s too short for the window', upper(CHAN{k}));
  end
  S.(CHAN{k}) = q.Data(r0:r1, :);
  if k == 1, Tv = q.Time(r0:r1); end
end
clear q
[Nt, Nx] = size(S.hh);
fprintf('%d samples x %d traces\n', Nt, Nx);
% the product carries the HH image it was built from: if the channel just
% loaded is not that image, the window above belongs to other data
if isfield(P, 'ref')
  A = P.ref; if isstruct(A), A = A.Data; end
  if isequal(size(A), size(S.hh))
    rel = max(abs(A(:) - S.hh(:))) / max(abs(A(:)));
    fprintf('source check: max rel diff %.3g%s\n', rel, ...
      H_tag(rel < 1e-6, ' (IDENTICAL)', ' <-- NOT the product source'));
  end
  clear A
end

%% 3. depth axis
C0 = 299792458; C_ICE = C0 / sqrt(3.171);
st = median(P.Surface(:), 'omitnan');
clear P
z = (Tv(:) - st) * C_ICE / 2;
keep = z >= 0 & z <= Z_MAX;
for k = 1:4, S.(CHAN{k}) = S.(CHAN{k})(keep, :); end
z = z(keep);
fprintf('depth window %.0f..%.0f m (%d samples)\n', z(1), z(end), numel(z));
% the band the before/after coherences are read over, clamped at both ends
% to the record that exists (thin ice ends above 1200 m)
if z(end) < Z_BAND(2), Z_BAND(2) = z(end); end
if Z_BAND(1) >= Z_BAND(2), Z_BAND(1) = z(1); end
band = z > Z_BAND(1) & z < Z_BAND(2);
kr = ones(NRW, 1) / NRW;

%% 4. one delay per channel from the firn, or the cache of a previous run
CO = struct('mode', 'constant', 'firn_m', coreg_firn_m, 'search_t', 5, 'oversample', 16);
cache_dir = fullfile(out_dir, 'coreg_cache');
cache_fn = fullfile(cache_dir, sprintf('creg_%s.mat', tag));
if exist(cache_fn, 'file') == 2
  % the trace count is read from the cached channel itself, the one record
  % of it that cannot drift from the data
  same = false;
  try
    cq = load(cache_fn, 'z', 'CO', 'r0', 'r1');
    wi = whos('-file', cache_fn, 'hh');
    same = isequal(cq.CO, CO) && cq.r0 == r0 && cq.r1 == r1 ...
      && ~isempty(wi) && isequal(double(wi(1).size(2)), double(Nx)) ...
      && numel(cq.z) == numel(z) && max(abs(cq.z(:) - z(:))) < 1e-6;
  catch err
    fprintf('cache %s is unreadable (%s)\n', cache_fn, err.message);
  end
  if same
    fprintf('coregistered channels already cached: %s\n', cache_fn);
    return
  end
  fprintf('cache %s does not match these settings; coregistering again\n', cache_fn);
end
fprintf('coregistration: one delay per channel from the firn, %g-%g m\n', coreg_firn_m(1), coreg_firn_m(2));
pairs = {'hh', 'vv'; 'hh', 'hv'; 'hh', 'vh'; 'hv', 'vh'};
coh_before = zeros(size(pairs, 1), 1);
for p = 1:size(pairs, 1)
  coh_before(p) = H_coh(S.(pairs{p, 1}), S.(pairs{p, 2}), kr, band);
end
copts = CO;
copts.rows = z >= CO.firn_m(1) & z <= CO.firn_m(2);
[T, cinfo] = ptt.coregisterChannels(S, 'hh', copts);
fprintf('%-8s %8s %8s\n', 'pair', 'before', 'after');
for p = 1:size(pairs, 1)
  fprintf('%-8s %8.3f %8.3f\n', [upper(pairs{p, 1}) '-' upper(pairs{p, 2})], ...
    coh_before(p), H_coh(T.(pairs{p, 1}), T.(pairs{p, 2}), kr, band));
end
fprintf('row offsets: VV %+.3f  HV %+.3f  VH %+.3f\n', ...
  cinfo.row_med.vv, cinfo.row_med.hv, cinfo.row_med.vh);

%% 5. cache the coregistered channels
if exist(cache_dir, 'dir') ~= 7, mkdir(cache_dir); end
hh = single(T.hh); vv = single(T.vv); hv = single(T.hv); vh = single(T.vh);
tmp_fn = [cache_fn '.tmp'];
save(tmp_fn, '-v7.3', 'hh', 'vv', 'hv', 'vh', 'z', 'CO', 'r0', 'r1', 'Nx');
movefile(tmp_fn, cache_fn);
fprintf('cached coregistered channels -> %s (%.1f min)\n', cache_fn, toc(t_all) / 60);

function c = H_coh(a, b, kr, band)
% median coherence magnitude of a pair over the band, NRW samples in range
num = conv2(a .* conj(b), kr, 'same');
d1 = conv2(abs(a) .^ 2, kr, 'same');
d2 = conv2(abs(b) .^ 2, kr, 'same');
cc = abs(num) ./ sqrt(max(d1 .* d2, realmin));
c = median(cc(band, :), 'all', 'omitnan');
end

function s = H_tag(ok, yes, no)
if ok, s = yes; else, s = no; end
end

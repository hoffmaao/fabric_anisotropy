%QUADPOL_COREG_FRAME Coregister every polarimetric pair, as production does.
%
% Aligns VV, HV and VH onto HH with the OPR toolbox `coregistration`, using
% the SAME window and the SAME settings the shipped CSARP_polarimetric
% product used for its HH/VV pair, so the cross-polarized channels are
% treated identically to the co-polarized one.
%
% Both are read from the product rather than assumed:
%   window   param_polarimetric.min_rbin : max_rbin  (100:6700 on Ridge A)
%   tiling   param_polarimetric.coregistration       (Tt 101, Tx 301,
%                                                     overlap 50/150,
%                                                     search 5/5, 1-D)
%
% The product also records which standardphase products it was built from
% (HH_path etc). Those names carry a processing-time suffix that no longer
% exists on disk, so the first thing this does is CHECK that the product's
% own `ref` matches the standardphase HH over the recorded window - if the
% source has been renamed but is otherwise the same data, that comparison
% is exact, and if it is not the same data every number downstream is
% being computed on a different version than the products were.
%
% Saves the offset fields DECIMATED. Writing them at full resolution made a
% 499 MB file per frame - they are 6601 x 4346 doubles, i.e. bigger than
% the images they describe - and they carry nothing at that resolution:
% the field is a near-constant instrument delay, +0.777 samples for VV with
% a 0.19-1.75 5-95 spread, and col_offset is identically zero. A 32x
% decimation preserves every feature the tiling could resolve (tiles are
% 101 x 301 samples) and lands under a megabyte.
if ~exist('site_root', 'var') || isempty(site_root)
  site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
end
if ~exist('day_seg', 'var') || isempty(day_seg), day_seg = '20250108_02'; end
if ~exist('frm', 'var') || isempty(frm), frm = 9; end
% repo root and <work> derived from this script's own location - see
% fabric_paths
[fabric_code, fabric_work] = fabric_paths();
if ~exist('out_dir', 'var') || isempty(out_dir)
  out_dir = fullfile(fabric_work, 'stages', 'coreg');
end
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
addpath(fabric_code);

CHAN = {'hh','vv','hv','vh'};
NRW = 101;             % range bins in the per-trace coherence window
Z_BAND = [200 1200];
name = sprintf('Data_%s_%03d.mat', day_seg, frm);

% --- the product: window, settings, and the reference to check against
pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
if exist(pol_fn, 'file') ~= 2
  error('quadpol_coreg_frame:noProduct', 'no polarimetric product %s', pol_fn);
end
w = {whos('-file', pol_fn).name};
sel = intersect({'param_polarimetric','Time','Surface','ref','sec_reg'}, w);
P = load(pol_fn, sel{:});
pp = P.param_polarimetric;
if isfield(pp, 'polarimetric'), pp = pp.polarimetric; end
r0 = pp.min_rbin; r1 = pp.max_rbin;
co = pp.coregistration;
CO = struct('Tt', co.Tt, 'Tx', co.Tx, 'overlap_t', co.overlap_t, ...
  'overlap_x', co.overlap_x, 'search_t', co.search_t, ...
  'search_x', co.search_x, 'one_dim_search_en', co.one_dim_search_en);
fprintf('product window rbin %d:%d (%d samples); tiling Tt %d Tx %d ov %d/%d\n', ...
  r0, r1, r1-r0+1, CO.Tt, CO.Tx, CO.overlap_t, CO.overlap_x);

% --- the four channels, cropped to the product's own window
S = struct();
for k = 1:4
  fn = fullfile(site_root, ['CSARP_standardphase_' upper(CHAN{k})], ...
    day_seg, name);
  if exist(fn, 'file') ~= 2
    error('quadpol_coreg_frame:missing', 'missing %s', fn);
  end
  q = load(fn, 'Data', 'Time');
  if r1 > size(q.Data, 1)
    error('quadpol_coreg_frame:window', ...
      '%s has %d samples, product window needs %d', CHAN{k}, ...
      size(q.Data,1), r1);
  end
  S.(CHAN{k}) = q.Data(r0:r1, :);
  if k == 1, Tv = q.Time(r0:r1); end
end
fprintf('cropped to %d x %d\n', size(S.hh, 1), size(S.hh, 2));

% --- is this the same data the product was built from?
if isfield(P, 'ref')
  A = P.ref; if isstruct(A), A = A.Data; end
  if isequal(size(A), size(S.hh))
    d = abs(A - S.hh);
    rel = max(d(:)) / max(abs(A(:)));
    fprintf('ref vs standardphase HH over the window: max rel diff %.3g %s\n', ...
      rel, H_verdict(rel));
  else
    fprintf('ref is %s but the cropped HH is %s - NOT the same window\n', ...
      mat2str(size(A)), mat2str(size(S.hh)));
  end
end

% --- depth axis, for the coherence band only
C0 = 299792458; C_ICE = C0/sqrt(3.171);
st = median(P.Surface(:), 'omitnan');
z = (Tv(:) - st) * C_ICE / 2;
band = z > Z_BAND(1) & z < Z_BAND(2);
kr = ones(NRW,1) / NRW;

pairs = {'hh','vv'; 'hh','hv'; 'hh','vh'; 'hv','vh'};
before = zeros(size(pairs,1), 1);
for p = 1:size(pairs,1)
  before(p) = H_coh(S.(pairs{p,1}), S.(pairs{p,2}), kr, band);
end

t0 = tic;
[T, info] = ptt.coregisterChannels(S, 'hh', CO);
fprintf('\ncoregistration of 3 pairs: %.1f min\n', toc(t0)/60);

% Each H_coh is three full-grid convolutions, so the values printed here
% are the ones saved below rather than being measured a second time.
after = zeros(size(pairs,1), 1);
fprintf('\n%-10s %10s %10s %8s\n', 'pair', 'before', 'after', 'change');
for p = 1:size(pairs,1)
  after(p) = H_coh(T.(pairs{p,1}), T.(pairs{p,2}), kr, band);
  fprintf('%-10s %10.3f %10.3f %+8.3f\n', ...
    [upper(pairs{p,1}) '-' upper(pairs{p,2})], before(p), after(p), ...
    after(p)-before(p));
end
% The shipped pair, as the target to match
if isfield(P, 'ref') && isfield(P, 'sec_reg')
  A = P.ref; if isstruct(A), A = A.Data; end
  fprintf('%-10s %10s %10.3f  (shipped)\n', 'HH-VV', '', ...
    H_coh(A, P.sec_reg, kr, band));
end

fprintf('\nmedian row offsets (samples relative to HH):\n');
for k = 2:4
  fprintf('  %-3s %+.3f\n', upper(CHAN{k}), info.row_med.(CHAN{k}));
end

DEC = 32;
row_offset = struct(); col_offset = struct();
for k = 2:4
  ch = CHAN{k};
  row_offset.(ch) = single(info.row_offset.(ch)(1:DEC:end, 1:DEC:end));
  col_offset.(ch) = single(info.col_offset.(ch)(1:DEC:end, 1:DEC:end));
end
row_med = info.row_med;
col_med = struct();
for k = 2:4
  c = info.col_offset.(CHAN{k});
  col_med.(CHAN{k}) = median(c(isfinite(c)));
end
coh_before = before; coh_after = after;
pair_names = pairs;
out_fn = fullfile(out_dir, sprintf('coreg_%s_%03d.mat', day_seg, frm));
save(out_fn, '-v7.3', 'row_offset', 'col_offset', 'row_med', 'col_med', ...
  'z', 'day_seg', 'frm', 'CO', 'r0', 'r1', 'DEC', 'coh_before', ...
  'coh_after', 'pair_names');
fprintf('\nwrote %s\n', out_fn);

function c = H_coh(a, b, kr, band)
num = conv2(a .* conj(b), kr, 'same');
d1 = conv2(abs(a).^2, kr, 'same');
d2 = conv2(abs(b).^2, kr, 'same');
cc = abs(num) ./ sqrt(max(d1.*d2, realmin));
c = median(cc(band,:), 'all', 'omitnan');
end

function s = H_verdict(rel)
if rel < 1e-6, s = '-> IDENTICAL';
elseif rel < 1e-3, s = '-> same to rounding';
else, s = '-> DIFFERENT DATA';
end
end

%RUN_QUADPOL_PIPELINE Coregister and invert one profile, end to end.
%
% The complete quad-pol workflow for a single frame, in one pass:
%
%   1. load all four standardphase channels on the product's own window
%   2. verify the source matches what the shipped product was built from
%   3. coregister VV, HV and VH onto HH with the OPR toolbox, using the
%      settings the product recorded
%   4. invert with ptt.ershadiFabric on the COREGISTERED channels
%   5. report against the independent two-azimuth solve, and save
%
% Coregistration and inversion are deliberately in the SAME script. The
% coregistered images are ~370 MB a frame, so writing them out to be read
% back by a separate inversion step would cost more in I/O than the
% inversion itself, and would invite the two halves drifting apart on which
% window or which settings they used - which is exactly the mistake that
% produced the retracted result in 03d292e, where the inversion ran on
% channels the coregistration step had never touched.
%
% Runtime is dominated by step 3: ~71 min on a 6601 x 4346 frame, less in
% proportion for shorter ones. Run per profile.
%
%   matlab -batch "day_seg='20250108_02'; frm=9; run_quadpol_pipeline"
if ~exist('site_root', 'var') || isempty(site_root)
  site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
end
if ~exist('day_seg', 'var') || isempty(day_seg), day_seg = '20250108_02'; end
if ~exist('frm', 'var') || isempty(frm), frm = 9; end
if ~exist('out_dir', 'var') || isempty(out_dir)
  out_dir = '/kucresis/scratch/hoffmana_sta/fabric/stages/quadpol';
end
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
addpath('/kucresis/scratch/hoffmana_sta/fabric/code');

CHAN = {'hh','vv','hv','vh'};
NRW = 101;
Z_BAND = [200 1200];
Z_MAX = 1500;
FC = 750e6;
PSI_STEP_DEG = 1;
name = sprintf('Data_%s_%03d.mat', day_seg, frm);
t_all = tic;

%% 1-2. window, settings and source check, from the product
pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
if exist(pol_fn, 'file') ~= 2
  error('run_quadpol_pipeline:noProduct', 'no polarimetric product %s', pol_fn);
end
w = {whos('-file', pol_fn).name};
sel = intersect({'param_polarimetric','Time','Surface','Latitude', ...
  'Longitude','ref','sec_reg'}, w);
P = load(pol_fn, sel{:});
pp = P.param_polarimetric;
if isfield(pp, 'polarimetric'), pp = pp.polarimetric; end
r0 = pp.min_rbin; r1 = pp.max_rbin;
co = pp.coregistration;
CO = struct('Tt', co.Tt, 'Tx', co.Tx, 'overlap_t', co.overlap_t, ...
  'overlap_x', co.overlap_x, 'search_t', co.search_t, ...
  'search_x', co.search_x, 'one_dim_search_en', co.one_dim_search_en);
fprintf('=== %s_%03d ===\n', day_seg, frm);
fprintf('window rbin %d:%d, tiling Tt %d Tx %d ov %d/%d\n', r0, r1, ...
  CO.Tt, CO.Tx, CO.overlap_t, CO.overlap_x);

S = struct();
for k = 1:4
  fn = fullfile(site_root, ['CSARP_standardphase_' upper(CHAN{k})], ...
    day_seg, name);
  if exist(fn, 'file') ~= 2
    error('run_quadpol_pipeline:missing', 'missing %s channel', CHAN{k});
  end
  q = load(fn, 'Data', 'Time');
  if r1 > size(q.Data, 1)
    error('run_quadpol_pipeline:window', '%s too short for the window', ...
      CHAN{k});
  end
  S.(CHAN{k}) = q.Data(r0:r1, :);
  if k == 1, Tv = q.Time(r0:r1); end
end
[Nt, Nx] = size(S.hh);
fprintf('%d samples x %d traces\n', Nt, Nx);
if isfield(P, 'ref')
  A = P.ref; if isstruct(A), A = A.Data; end
  if isequal(size(A), size(S.hh))
    rel = max(abs(A(:) - S.hh(:))) / max(abs(A(:)));
    fprintf('source check: max rel diff %.3g%s\n', rel, ...
      H_tag(rel < 1e-6, ' (IDENTICAL)', ' <-- NOT the product source'));
  end
end

%% depth axis
C0 = 299792458; C_ICE = C0/sqrt(3.171);
st = median(P.Surface(:), 'omitnan');
z = (Tv(:) - st) * C_ICE / 2;
keep = z >= 0 & z <= Z_MAX;
for k = 1:4, S.(CHAN{k}) = S.(CHAN{k})(keep, :); end
z = z(keep);
band = z > Z_BAND(1) & z < Z_BAND(2);
kr = ones(NRW,1)/NRW;
fprintf('depth window %.0f..%.0f m (%d samples)\n', z(1), z(end), numel(z));

%% 3. coregistration
pairs = {'hh','vv'; 'hh','hv'; 'hh','vh'; 'hv','vh'};
coh_before = zeros(size(pairs,1),1);
for p = 1:size(pairs,1)
  coh_before(p) = H_coh(S.(pairs{p,1}), S.(pairs{p,2}), kr, band);
end
t0 = tic;
[T, cinfo] = ptt.coregisterChannels(S, 'hh', CO);
t_coreg = toc(t0)/60;
coh_after = zeros(size(pairs,1),1);
for p = 1:size(pairs,1)
  coh_after(p) = H_coh(T.(pairs{p,1}), T.(pairs{p,2}), kr, band);
end
fprintf('\ncoregistration %.1f min\n', t_coreg);
fprintf('%-8s %8s %8s\n', 'pair', 'before', 'after');
for p = 1:size(pairs,1)
  fprintf('%-8s %8.3f %8.3f\n', ...
    [upper(pairs{p,1}) '-' upper(pairs{p,2})], coh_before(p), coh_after(p));
end
fprintf('row offsets: VV %+.3f  HV %+.3f  VH %+.3f\n', ...
  cinfo.row_med.vv, cinfo.row_med.hv, cinfo.row_med.vh);

%% 4. inversion, on the coregistered channels
t0 = tic;
out = ptt.ershadiFabric(T, z, struct('fc', FC, 'psi_step_deg', PSI_STEP_DEG, ...
  'win_m', 30, 'grad_win_m', 25, 'coh_min', 0.4, 'deramped', true));
% and on the UNCOREGISTERED channels, so the effect of step 3 on the
% ANSWER - not just on the coherence - is measured rather than assumed
out_raw = ptt.ershadiFabric(S, z, struct('fc', FC, ...
  'psi_step_deg', PSI_STEP_DEG, 'win_m', 30, 'grad_win_m', 25, ...
  'coh_min', 0.4, 'deramped', true));
fprintf('\ninversion %.1f min\n', toc(t0)/60);

la = P.Latitude(:); lo = P.Longitude(:);
p0 = deg2rad(la(1)); p1 = deg2rad(la(end));
dl = deg2rad(lo(end) - lo(1));
track_az = mod(rad2deg(atan2(sin(dl)*cos(p1), ...
  cos(p0)*sin(p1) - sin(p0)*cos(p1)*cos(dl))), 180);

zb = band;
fprintf('\n%-14s %10s %10s\n', '', 'coreg', 'raw');
fprintf('%-14s %10.1f %10.1f\n', 'theta_ant', ...
  H_cmed(rad2deg(out.theta(zb))), H_cmed(rad2deg(out_raw.theta(zb))));
fprintf('%-14s %10.1f %10.1f\n', 'theta_geo', ...
  mod(H_cmed(rad2deg(out.theta(zb))) + track_az, 180), ...
  mod(H_cmed(rad2deg(out_raw.theta(zb))) + track_az, 180));
fprintf('%-14s %10.3f %10.3f\n', 'dlam', ...
  median(out.dlam(zb), 'omitnan'), median(out_raw.dlam(zb), 'omitnan'));
fprintf('%-14s %10.3f %10.3f\n', '|C_hhvv|', ...
  median(out.Cmag(zb,1), 'omitnan'), median(out_raw.Cmag(zb,1), 'omitnan'));
fprintf('%-14s %10.1f\n', 'track_az', track_az);
fprintf('\n(two-azimuth solve for Ridge A: dlam ~0.05, theta ~90 deg)\n');

%% 5. save, small
ZS = 4;
s = 1:ZS:numel(z);
res = struct('tag', sprintf('%s_%03d', day_seg, frm), ...
  'day_seg', day_seg, 'frm', frm, 'track_az', track_az, ...
  'z', single(z(s)), ...
  'theta', single(rad2deg(out.theta(s))), ...
  'dlam', single(out.dlam(s)), ...
  'Cmag', single(out.Cmag(s,1)), ...
  'theta_raw', single(rad2deg(out_raw.theta(s))), ...
  'dlam_raw', single(out_raw.dlam(s)), ...
  'row_med', cinfo.row_med, 'coh_before', coh_before, ...
  'coh_after', coh_after, 'pairs', {pairs}, ...
  'lat', median(la), 'lon', median(lo), ...
  'lat0', la(1), 'lon0', lo(1), 'lat1', la(end), 'lon1', lo(end), ...
  't_coreg_min', t_coreg);
out_fn = fullfile(out_dir, sprintf('quadpol_%s_%03d.mat', day_seg, frm));
save(out_fn, '-v7.3', 'res');
fprintf('\nwrote %s (total %.1f min)\n', out_fn, toc(t_all)/60);

function c = H_coh(a, b, kr, band)
num = conv2(a .* conj(b), kr, 'same');
d1 = conv2(abs(a).^2, kr, 'same');
d2 = conv2(abs(b).^2, kr, 'same');
cc = abs(num) ./ sqrt(max(d1.*d2, realmin));
c = median(cc(band,:), 'all', 'omitnan');
end

function m = H_cmed(a)
a = a(isfinite(a));
if isempty(a), m = NaN; return; end
m = mod(rad2deg(angle(mean(exp(2i*deg2rad(a)))))/2, 180);
end

function s = H_tag(ok, yes, no)
if ok, s = yes; else, s = no; end
end

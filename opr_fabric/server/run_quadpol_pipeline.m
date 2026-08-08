%RUN_QUADPOL_PIPELINE Coregister and invert one profile, end to end.
%
% The complete quad-pol workflow for a single frame, in one pass:
%
%   1. load all four standardphase channels on the product's own window
%   2. verify the source matches what the shipped product was built from
%   3. coregister VV, HV and VH onto HH with the OPR toolbox, using the
%      settings the product recorded - or load the coreg_cache from a
%      previous run, which turns a rerun from ~52 min into a few
%   4. invert the COREGISTERED channels twice: ptt.ershadiFabric (the
%      published chain, evaluated at the cross-pol minimum, which on this
%      system is antenna-locked) and ptt.quadpolFabricLS (the model fit
%      that takes the axis from the coherence field; frame-level theta0,
%      then per-block dlam with theta0 held)
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
% Traces per along-track block for the SECTION. 125 is what run_sections.m
% uses for the co-polarized Ridge A section, so the two sections have the
% same along-track sampling and can be read against each other cell for
% cell rather than approximately.
NBLK_TR = 125;
CACHE_COREG = true;   % keep the coregistered channels; see below
name = sprintf('Data_%s_%03d.mat', day_seg, frm);
t_all = tic;

%% 1-2. window, settings and source check, from the product
pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
if exist(pol_fn, 'file') ~= 2
  % The 2022/2023 seasons shipped the polarimetric product only in its
  % SNAPHU-unwrapped form. It carries the same param_polarimetric (window
  % and coregistration settings), Surface and ref, so it serves the same
  % role here.
  alt = fullfile(site_root, 'CSARP_polarimetric_unwrap', day_seg, name);
  if exist(alt, 'file') == 2
    pol_fn = alt;
  end
end
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
% Positions, needed by the section loop below as well as by the reporting,
% so they are defined once here rather than after the first use.
la = P.Latitude(:); lo = P.Longitude(:);
fprintf('depth window %.0f..%.0f m (%d samples)\n', z(1), z(end), numel(z));

%% 3. coregistration, or the cache of a previous run
% Coregistration is ~50 min and the inversions are minutes, so a rerun for
% a different estimator or setting should never pay it again. The cache is
% keyed by the settings that produced it; any mismatch falls through to a
% fresh coregistration.
pairs = {'hh','vv'; 'hh','hv'; 'hh','vh'; 'hv','vh'};
cache_fn = fullfile(out_dir, 'coreg_cache', ...
  sprintf('creg_%s_%03d.mat', day_seg, frm));
from_cache = false;
if exist(cache_fn, 'file') == 2
  cq = load(cache_fn, 'z', 'CO', 'r0', 'r1');
  if isequal(cq.CO, CO) && cq.r0 == r0 && cq.r1 == r1 ...
      && numel(cq.z) == numel(z) && max(abs(cq.z(:) - z(:))) < 1e-6
    cq = load(cache_fn, 'hh', 'vv', 'hv', 'vh');
    T = struct('hh', double(cq.hh), 'vv', double(cq.vv), ...
      'hv', double(cq.hv), 'vh', double(cq.vh));
    clear cq;
    from_cache = true;
    fprintf('\ncoregistered channels loaded from %s\n', cache_fn);
  else
    fprintf('\ncache %s does not match the settings; re-coregistering\n', ...
      cache_fn);
    clear cq;
  end
end
if from_cache
  t_coreg = 0;
  coh_before = nan(size(pairs,1),1);
  coh_after = nan(size(pairs,1),1);
  cinfo = struct('row_med', struct('vv', NaN, 'hv', NaN, 'vh', NaN));
else
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
end

%% 3b. cache the coregistered channels
% ~730 MB a frame, which is worth it: coregistration is 24-56 min and the
% inversion is 6 seconds, so anything that needs re-inverting - a different
% window, a different block size, a bug - would otherwise pay the whole
% coregistration again. The cache is keyed by frame and by the settings
% that produced it.
if CACHE_COREG && ~from_cache
  cdir = fullfile(out_dir, 'coreg_cache');
  if exist(cdir, 'dir') ~= 7, mkdir(cdir); end
  hh = single(T.hh); vv = single(T.vv); hv = single(T.hv); vh = single(T.vh);
  save(cache_fn, '-v7.3', 'hh', 'vv', 'hv', 'vh', 'z', 'CO', 'r0', 'r1');
  clear hh vv hv vh;
  fprintf('cached coregistered channels -> %s\n', cache_fn);
end

%% 3c. coreg-only mode, for building the caches in bulk
% The caches are estimator-independent, so the expensive step can run for
% the whole survey before the estimator settings are final: with
% coreg_only=true the script stops here and leaves the inversions to a
% later rerun, which the cache makes cheap.
if exist('coreg_only', 'var') && isequal(coreg_only, true)
  fprintf('coreg_only: stopping after the cache (total %.1f min)\n', ...
    toc(t_all)/60);
  return;
end

%% 4. inversion, on the coregistered channels
t0 = tic;
out = ptt.ershadiFabric(T, z, struct('fc', FC, 'psi_step_deg', PSI_STEP_DEG, ...
  'win_m', 30, 'grad_win_m', 25, 'coh_min', 0.4, 'deramped', true));
% and on the UNCOREGISTERED channels, so the effect of step 3 on the
% ANSWER - not just on the coherence - is measured rather than assumed
out_raw = ptt.ershadiFabric(S, z, struct('fc', FC, ...
  'psi_step_deg', PSI_STEP_DEG, 'win_m', 30, 'grad_win_m', 25, ...
  'coh_min', 0.4, 'deramped', true));
fprintf('\nershadi inversion %.1f min\n', toc(t0)/60);

%% 4a. the LS inversion: frame-level theta0, then dlam with theta0 fixed
% ershadiFabric evaluates Psi at the cross-polarized minimum, which on
% this system is antenna-locked (89.6 +- 1.9 deg over 1790 blocks) by the
% flat -3.6 dB cross-pol pedestal, i.e. ~25 deg off the true axis. That is
% what put the odd-pi coherence-null bands and the fringe-periodic
% roughness in the section, and scaled dlam by ~cos 2*25 deg. The LS fit
% takes the axis from the coherence field itself and models the nulls
% instead of gating on them; see ptt.quadpolFabricLS and test_quadpol_ls.
t0 = tic;
lsq = ptt.quadpolFabricLS(T, z, struct('fc', FC, 'deramped', true));
okt = isfinite(lsq.theta0);
if nnz(okt) >= 2
  % Hand the blocks a SMOOTHED axis, built on the doubled-angle phasor
  % weighted by each window's own theta0 contrast - never an unwrapped
  % angle. A weak window can fit the conjugate branch (a 90 deg flip);
  % smoothing the phasor votes it down, where an unwrap would have
  % propagated it to every window below as a silent axis slip.
  qw = lsq.q_theta;
  qw(~isfinite(qw) | qw < 0) = 0;
  ph = qw .* exp(2i * lsq.theta0);
  ph(~okt) = 0;
  ks = ones(5, 1);
  num_p = conv(ph, ks, 'same');
  den_p = conv(qw .* double(okt), ks, 'same');
  oks = den_p > 0.25 & abs(num_p) > 0;
  if nnz(oks) >= 2
    th_prof = struct('z', lsq.zw(oks), 'theta', 0.5 * angle(num_p(oks)));
  else
    th_prof = struct('z', lsq.zw(okt), 'theta', lsq.theta0(okt));
  end
else
  th_prof = [];   % nothing usable; let the blocks estimate their own
end
% The blocks inherit the frame-level pedestal along with theta0: both are
% instrument-or-site constants at block scale, and pinning them keeps the
% per-block problem two-parameter and immune to the per-window
% pedestal/fabric confusion at axes near 45 deg to the antennas.
if all(isfinite(lsq.pedestal))
  blk_ped = lsq.pedestal;
else
  blk_ped = 'frame';
end
fprintf(['LS frame pass %.1f min: theta0 constrained on %d of %d windows, ' ...
  'pedestal [%.3f %+.3fi %.3f]\n'], toc(t0)/60, nnz(okt), ...
  numel(lsq.theta0), lsq.pedestal(1), lsq.pedestal(2), lsq.pedestal(3));

%% 4b. the SECTION: both estimators, per along-track block
% The frame-average profile above answers "what is the fabric here"; this
% answers "how does it vary along the line", which is the quantity the
% co-polarized section already shows and the one worth comparing against.
% The LS blocks inherit the frame theta0 profile, so their two remaining
% parameters are purely local and the section carries no block-to-block
% axis jitter.
t0 = tic;
nb = max(1, floor(Nx / NBLK_TR));
sec_dlam = nan(numel(z), nb);
sec_theta = nan(numel(z), nb);
sec_cmag = nan(numel(z), nb);
sec_dlam_ls = nan(numel(z), nb);
sec_resid_ls = nan(numel(z), nb);
sec_lat = nan(1, nb); sec_lon = nan(1, nb); sec_az = nan(1, nb);
for b = 1:nb
  j0 = (b-1)*NBLK_TR + 1;
  j1 = min(b*NBLK_TR, Nx);
  if j1 - j0 < 16, continue; end
  Tb = struct();
  for k = 1:4, Tb.(CHAN{k}) = T.(CHAN{k})(:, j0:j1); end
  ob = ptt.ershadiFabric(Tb, z, struct('fc', FC, ...
    'psi_step_deg', PSI_STEP_DEG, 'win_m', 30, 'grad_win_m', 25, ...
    'coh_min', 0.4, 'deramped', true));
  sec_dlam(:, b) = ob.dlam;
  sec_theta(:, b) = rad2deg(ob.theta);
  sec_cmag(:, b) = ob.Cmag(:, 1);
  ob_ls = ptt.quadpolFabricLS(Tb, z, struct('fc', FC, 'deramped', true, ...
    'theta0', th_prof, 'pedestal', blk_ped));
  sec_dlam_ls(:, b) = ob_ls.dlam_z;
  sec_resid_ls(:, b) = interp1(ob_ls.zw, ob_ls.resid, z, 'linear');
  sec_lat(b) = mean(la(j0:j1));
  sec_lon(b) = mean(lo(j0:j1));
  p0b = deg2rad(la(j0)); p1b = deg2rad(la(j1));
  dlb = deg2rad(lo(j1) - lo(j0));
  sec_az(b) = mod(rad2deg(atan2(sin(dlb)*cos(p1b), ...
    cos(p0b)*sin(p1b) - sin(p0b)*cos(p1b)*cos(dlb))), 180);
  clear Tb ob ob_ls;
end
fprintf('section: %d blocks of %d traces, %.1f min\n', nb, NBLK_TR, ...
  toc(t0)/60);
fprintf('section dlam %.3f..%.3f (median %.3f) | LS median %.3f\n', ...
  min(sec_dlam(:)), max(sec_dlam(:)), median(sec_dlam(:), 'omitnan'), ...
  median(sec_dlam_ls(:), 'omitnan'));

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
zw_band = lsq.zw > Z_BAND(1) & lsq.zw < Z_BAND(2);
fprintf('\nLS fit (frame): theta0_ant %.1f  theta0_geo %.1f  dlam %.3f  ' ...
  , H_cmed(rad2deg(lsq.theta0(zw_band))), ...
  mod(H_cmed(rad2deg(lsq.theta0(zw_band))) + track_az, 180), ...
  median(lsq.dlam(zw_band), 'omitnan'));
fprintf('resid %.3f  q_theta %.2f\n', ...
  median(lsq.resid(zw_band), 'omitnan'), ...
  median(lsq.q_theta(zw_band), 'omitnan'));
fprintf(['(two-azimuth solve for Ridge A: principal contrast ~0.12 at the ' ...
  'plateau, axis ~110-120 deg geo;\n the ershadi rows above are the ' ...
  'antenna-locked projection and should sit LOW by ~cos 2*offset)\n']);

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
  'ls_zw', single(lsq.zw), 'ls_theta0', single(rad2deg(lsq.theta0)), ...
  'ls_dlam', single(lsq.dlam), 'ls_q_theta', single(lsq.q_theta), ...
  'ls_resid', single(lsq.resid), 'ls_gamma', single(lsq.gamma), ...
  'ls_leak', single(lsq.leak), 'ls_pedestal', lsq.pedestal, ...
  'theta0_ls', single(rad2deg(lsq.theta0_z(s))), ...
  'dlam_ls', single(lsq.dlam_z(s)), ...
  'sec_dlam_ls', single(sec_dlam_ls(s,:)), ...
  'sec_resid_ls', single(sec_resid_ls(s,:)), ...
  'sec_dlam', single(sec_dlam(s,:)), 'sec_theta', single(sec_theta(s,:)), ...
  'sec_cmag', single(sec_cmag(s,:)), 'sec_lat', sec_lat, ...
  'sec_lon', sec_lon, 'sec_az', sec_az, 'nblk_tr', NBLK_TR, ...
  'row_med', cinfo.row_med, 'coh_before', coh_before, ...
  'coh_after', coh_after, 'pairs', {pairs}, ...
  'lat', median(la), 'lon', median(lo), ...
  'lat0', la(1), 'lon0', lo(1), 'lat1', la(end), 'lon1', lo(end), ...
  't_coreg_min', t_coreg);
% `quadpol_section_` and not `quadpol_`: run_quadpol_frame.m already writes
% quadpol_<day_seg>_<frm>.mat with a flat out/z/psi/A layout, and both
% mirror to the same figure staging directory. Sharing the name would let
% whichever ran last break the other figure's reader, since this file is a
% single `res` struct instead.
out_fn = fullfile(out_dir, sprintf('quadpol_section_%s_%03d.mat', day_seg, frm));
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

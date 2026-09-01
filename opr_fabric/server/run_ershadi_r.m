%RUN_ERSHADI_R Reflection-ratio retrieval for one frame, from its cache.
%
% The Ershadi et al. (2022) Sect.-3.5 chain (ptt.ershadiFabric ->
% ptt.ershadiInverse) applied to the COREGISTERED channels of one frame,
% read from the coreg cache the quad-pol pipeline already built - so a
% frame costs minutes, not a coregistration. Products land in
% <work>/stages/ershadi_r/ershadi_r_<tag>.mat and never touch the
% quadpol_section products.
%
% WHY FRAME-LEVEL. The science question this feeds (see the crossing-pair
% record in scripts/figures/crossing_pairs.py) is whether the retrieved r
% differs between HEADING FAMILIES on the same ice: COF scattering predicts
% it does, a pure antenna pedestal predicts it does not. That comparison is
% per frame. Per-block r can come later if the frame-level answer warrants
% it.
%
% KNOWN CAVEAT, measured on synthetics (test_ershadi_r): an amplitude floor
% at -3 dB compresses retrieved |r_dB| from 10 to 5-8 by filling the co-pol
% nodes. The floor is common to every heading family on the same ice, so
% the FAMILY COMPARISON survives compression even where the absolute r is
% conservative; treat absolute values accordingly.
%
%   matlab -batch "addpath('<code>/opr_fabric/server'); \
%                  day_seg='20250108_02'; frm=9; run_ershadi_r"
if ~exist('day_seg', 'var') || ~exist('frm', 'var')
  error('run_ershadi_r:noFrame', 'set day_seg and frm at the call site');
end
day_seg = char(day_seg);
% Work root and repo root derived from this script's own location - see
% fabric_paths. No username is baked in.
[fabric_code, fabric_work] = fabric_paths();
addpath(fabric_code);

tag = sprintf('%s_%03d', day_seg, frm);
cache_fn = fullfile(fabric_work, 'stages', 'quadpol', 'coreg_cache', ...
  sprintf('creg_%s.mat', tag));
out_dir = fullfile(fabric_work, 'stages', 'ershadi_r');
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
out_fn = fullfile(out_dir, sprintf('ershadi_r_%s.mat', tag));
if exist(out_fn, 'file') == 2
  fprintf('%s exists; skipping (idempotent)\n', out_fn);
  return
end
if exist(cache_fn, 'file') ~= 2
  error('run_ershadi_r:noCache', ...
    'no coreg cache %s - build it with the quad-pol pipeline first', ...
    cache_fn);
end

t0 = tic;
cq = load(cache_fn, 'hh', 'vv', 'hv', 'vh', 'z');
S = struct('hh', double(cq.hh), 'vv', double(cq.vv), ...
  'hv', double(cq.hv), 'vh', double(cq.vh));
z = double(cq.z(:));
clear cq
fprintf('=== %s: %d x %d coregistered samples\n', tag, size(S.hh));

% Real CReSIS data is deramped (ershadiFabric E2 applies the conjugate).
fr = ptt.ershadiFabric(S, z, struct('fc', 750e6, 'psi_step_deg', 1, ...
  'deramped', true));

% Dome-C weights: coherence phase for theta (Ridge A's moderate anisotropy
% keeps the phase misfit usable), dP_HH for r - Table 3.
%
% INTERVALS ARE 200 m, NOT the paper's 50 m, and the reason is measured
% (see test_ershadi_r): interval length is bounded below by LEVERAGE - the
% share of the phase reaching an interval's rows that the interval itself
% contributes - and below ~20% the theta stage absorbs upstream error
% instead of fitting its own layer, an error MORE data makes worse. Ridge
% A's deep dlam ~0.06 over a 1400 m column puts a 50 m interval far under
% that bound; 200 m carries roughly a quarter of the accumulated phase at
% mid-column. r tolerates shorter intervals than theta (it is local to the
% reflector), so a finer r profile is available later by re-running with a
% shorter interval and trusting only r_db_int - the theta from such a run
% must not be quoted.
inv = ptt.ershadiInverse(fr, z, struct('interval_m', 200, ...
  'z_fit', [200, min(1400, z(end))], 'w_theta', [1 0 0], 'w_r', [0 1 0]));

% Carry the frame's geometry so the family analysis needs only this file:
% heading and the LS frame axis from the section product, when it exists.
track_az = NaN; ls_theta0_geo = NaN; ls_pedestal = nan(1, 3);
sec_fn = fullfile(fabric_work, 'stages', 'quadpol', ...
  sprintf('quadpol_section_%s.mat', tag));
if exist(sec_fn, 'file') == 2
  sec = load(sec_fn, 'res');
  if isfield(sec.res, 'track_az'), track_az = sec.res.track_az; end
  if isfield(sec.res, 'ls_theta0_geo')
    th = sec.res.ls_theta0_geo(isfinite(sec.res.ls_theta0_geo));
    if ~isempty(th)
      ls_theta0_geo = mod(rad2deg(angle(mean(exp(2i*deg2rad(th)))))/2, 180);
    end
  end
  if isfield(sec.res, 'ls_pedestal'), ls_pedestal = sec.res.ls_pedestal; end
  clear sec
end

res = struct('tag', tag, 'day_seg', day_seg, 'frm', frm, ...
  'z_int', inv.z_int, 'edges', inv.edges, ...
  'theta_int', inv.theta_int, 'r_db_int', inv.r_db_int, ...
  'dlam_int', inv.dlam_int, 'theta0_int', inv.theta0_int, ...
  'r13_db', inv.r13_db, 'r13_reason', inv.r13_reason, ...
  'J0', inv.J0, 'J', inv.J, 'exitflag', inv.exitflag, ...
  'track_az', track_az, 'ls_theta0_geo', ls_theta0_geo, ...
  'ls_pedestal', ls_pedestal, ...
  'coh_med', median(fr.Cmag(isfinite(fr.Cmag))), ...
  'runtime_s', toc(t0));
tmp_fn = [out_fn, '.tmp'];
save(tmp_fn, '-v7.3', 'res');
movefile(tmp_fn, out_fn);
fprintf('wrote %s (%.1f min)\n', out_fn, toc(t0)/60);

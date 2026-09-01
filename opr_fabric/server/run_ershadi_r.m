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

% The LS axis profile for this frame, in the ANTENNA frame the model works
% in - the section product stores theta0_ls in the same frame as the sweep,
% so no heading rotation is applied here. Read before the inversion so the
% intervals can be built from it; absent, the fit falls back to the
% data-derived initial guess per interval.
ls_theta_z = [];
sec_fn_pre = fullfile(fabric_work, 'stages', 'quadpol', ...
  sprintf('quadpol_section_%s.mat', tag));
if exist(sec_fn_pre, 'file') == 2
  sp = load(sec_fn_pre, 'res');
  if isfield(sp.res, 'theta0_ls') && isfield(sp.res, 'z')
    zz = double(sp.res.z(:));
    tt = double(sp.res.theta0_ls(:));    % degrees, antenna frame
    n = min(numel(zz), numel(tt));
    ls_theta_z = [zz(1:n), tt(1:n)];
  end
  clear sp
end

% THETA IS SUPPLIED FROM THE LS ESTIMATOR, NOT RE-FITTED. Two measured
% reasons, both on this site:
%   - the axis Ershadi's direct chain returns here is antenna-locked
%     (88.8 +- 2.0 deg over 36 Ridge A frames, i.e. pinned to 90 deg
%     regardless of heading or ice), which is why ptt.quadpolFabricLS
%     exists and why its axis is the one to trust;
%   - the Sect.-3.5 theta stage is ill-conditioned at this site whatever
%     the interval length. It needs an interval to carry a large share of
%     the phase reaching its rows; Ridge A's dlam ~0.06 gives a 200 m
%     interval ~3.6 of the ~19.5 rad accumulated by 1300 m (~18%), inside
%     the regime measured to fail. On 20250108_02_009 it threw stable
%     86-98 deg initial guesses to 13, 27 and 139 deg.
% r is NOT subject to that bound - it is local to the reflector - so it
% stays fitted, which is what this product is for.
%
% Intervals are 200 m rather than the paper's 50 m so each still holds
% enough rows for the r stage after decimation.
Z_FIT = [200, min(1400, z(end))];
INT_M = 200;
edges_pre = (0:INT_M:Z_FIT(2)).';
if edges_pre(end) < Z_FIT(2), edges_pre(end+1) = Z_FIT(2); end
th_fix = nan(numel(edges_pre) - 1, 1);
if exist('ls_theta_z', 'var') && ~isempty(ls_theta_z)
  for k = 1:numel(th_fix)
    m = ls_theta_z(:, 1) >= edges_pre(k) & ls_theta_z(:, 1) < edges_pre(k+1);
    v = ls_theta_z(m, 2); v = v(isfinite(v));
    if ~isempty(v)
      % circular mean of a mod-180 axis: never a plain average
      th_fix(k) = mod(angle(mean(exp(2i*deg2rad(v))))/2, pi);
    end
  end
end
inv = ptt.ershadiInverse(fr, z, struct('interval_m', INT_M, ...
  'z_fit', Z_FIT, 'w_theta', [1 0 0], 'w_r', [0 1 0], ...
  'theta_int_fixed', th_fix));

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
  'theta_fitted', inv.theta_fitted, ...
  'theta_phase_rad', inv.theta_phase_rad, ...
  'theta_was_supplied', inv.theta_was_supplied, ...
  'J0', inv.J0, 'J', inv.J, 'exitflag', inv.exitflag, ...
  'track_az', track_az, 'ls_theta0_geo', ls_theta0_geo, ...
  'ls_pedestal', ls_pedestal, ...
  'coh_med', median(fr.Cmag(isfinite(fr.Cmag))), ...
  'runtime_s', toc(t0));
tmp_fn = [out_fn, '.tmp'];
save(tmp_fn, '-v7.3', 'res');
movefile(tmp_fn, out_fn);
fprintf('wrote %s (%.1f min)\n', out_fn, toc(t0)/60);

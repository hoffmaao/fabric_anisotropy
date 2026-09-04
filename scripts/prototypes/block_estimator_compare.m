%BLOCK_ESTIMATOR_COMPARE Coherence LS vs power extinction on the same blocks.
%
%   blkmom_fn = '.../blkmom_20240120_06_001_B1000.mat';
%   out = block_estimator_compare(blkmom_fn)        % or run as a script with blkmom_fn set
%
% For every along-track block of a dump_block_moments.m file, fits the
% GEOGRAPHIC-frame moments with
%   - ptt.quadpolFabricLS  (coherence, per-window free axis, frame pedestal)
%   - ptt.quadpolFabricLS  in constant-orientation mode (held axis)
%   - ptt.quadpolFabricPower (power extinction: axis mod 90/180, |dlam|,
%     reflection ratio, node visibility kappa)
% and returns them side by side on the LS window grid, with the axis
% comparison done on the doubled angle (never unwrapped) and the power
% axis compared modulo 90 where its reflection ratio left the pair
% unresolved. What the comparison is for: at a site with an ice core the
% three should agree on the components the core measures - the horizontal
% eigenvalue difference dlam(z) and the orientation of the horizontal
% principal axes - and where they disagree, the disagreement is the
% measurement of each method's systematics.
function out = block_estimator_compare(blkmom_fn, opts)
if nargin < 2, opts = struct(); end
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));
fc = 750e6;
D = load(blkmom_fn);
z = double(D.z(:)); nb = numel(D.blk.x);
dlam_max = H_opt(opts, 'dlam_max', 0.25);
LSO = struct('fc', fc, 'deramped', true, 'dlam_max', dlam_max, 'psi_step_deg', 4, 'theta_step_deg', 4);
PWO = struct('fc', fc, 'psi_step_deg', 4, 'dlam_max', dlam_max, 'theta_step_deg', 4);
out = struct('tag', D.tag, 'x', D.blk.x, 'az', D.blk.az, 'lat', D.blk.lat, 'lon', D.blk.lon);
first = true;
equalise = H_opt(opts, 'equalise', true);
for b = 1:nb
  if D.blk.n(b) < 32, continue; end
  if equalise
    % range-dependent per-channel gains from reciprocity and the firn,
    % applied in the antenna frame, then the block's heading rotation
    [Ma, G] = ptt.equaliseChannels(double(squeeze(D.M_ant(:, :, :, b))), z);
    M = ptt.rotateMoments(Ma, -deg2rad(D.blk.az(b)));
    out.eq_firn_db(b) = G.firn_db; out.eq_seam_db(b) = G.seam_db;
  else
    M = double(squeeze(D.M_geo(:, :, :, b)));
  end
  t0 = tic;
  of = ptt.quadpolFabricLS(struct('M', M), z, LSO);
  oc = ptt.quadpolFabricLS(struct('M', M), z, setfield(LSO, 'theta_const', true)); %#ok<SFLD>
  op = ptt.quadpolFabricPower(struct('M', M), z, setfield(PWO, 'ped_az', deg2rad(D.blk.az(b)))); %#ok<SFLD>
  if first
    Nw = numel(of.zw);
    out.zw = of.zw;
    for f = {'th_free', 'dl_free', 'q_free', 'res_free', 'dl_held', 'res_held', ...
             'th_pow', 'dl_pow', 'kappa', 'rdb_pow', 'V_pow', 'res_pow'}
      out.(f{1}) = nan(Nw, nb);
    end
    out.th_held = nan(1, nb); out.spread = nan(1, nb); out.rdb_frame = nan(1, nb);
    out.mod90 = false(Nw, nb);
    first = false;
  end
  out.th_free(:, b) = of.theta0; out.dl_free(:, b) = of.dlam; out.q_free(:, b) = of.q_theta; out.res_free(:, b) = of.resid;
  out.dl_held(:, b) = oc.dlam; out.res_held(:, b) = oc.resid;
  out.th_held(b) = oc.theta_const; out.spread(b) = oc.theta_spread_deg;
  % power outputs on their own (150 m) window grid, interpolated to the LS grid
  out.th_pow(:, b) = angle(interp1(op.zw, exp(2i*op.theta0), of.zw, 'linear')) / 2;
  out.mod90(:, b) = interp1(op.zw, double(op.theta_mod90), of.zw, 'nearest') > 0.5;
  out.dl_pow(:, b) = interp1(op.zw, op.dlam, of.zw, 'linear');
  out.kappa(:, b) = interp1(op.zw, op.kappa, of.zw, 'linear');
  out.rdb_pow(:, b) = interp1(op.zw, op.r_db, of.zw, 'linear');
  out.V_pow(:, b) = interp1(op.zw, op.V, of.zw, 'linear');
  out.res_pow(:, b) = interp1(op.zw, op.resid, of.zw, 'linear');
  out.rdb_frame(b) = op.r_db_frame;
  % axis agreement (deg): held LS axis vs power axis, mod 180 or mod 90
  thp = angle(sum(exp(2i * op.theta0(isfinite(op.theta0))))) / 2;
  d180 = rad2deg(angle(exp(2i * (thp - oc.theta_const))) / 2);
  d90 = rad2deg(angle(exp(4i * (thp - oc.theta_const))) / 4);
  m = of.zw >= 200;
  fprintf(['%s blk %2d x %5.0f m (%.1f min): LS held %6.1f (spread %4.1f) | power %6.1f (r %+4.1f dB) ' ...
    'daxis mod180 %+5.1f mod90 %+5.1f | dlam LS %.3f held %.3f power %.3f | kappa %.2f\n'], ...
    D.tag, b, D.blk.x(b), toc(t0)/60, mod(rad2deg(oc.theta_const), 180), oc.theta_spread_deg, ...
    mod(rad2deg(thp), 180), op.r_db_frame, d180, d90, ...
    median(of.dlam(m), 'omitnan'), median(oc.dlam(m), 'omitnan'), median(out.dl_pow(m, b), 'omitnan'), ...
    median(op.kappa, 'omitnan'));
end
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

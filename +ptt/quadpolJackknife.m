function J = quadpolJackknife(Msub, nsub, z, os, ref, opts)
%QUADPOLJACKKNIFE Standard errors of a pooled quad-pol fit by sub-block jackknife.
%
% J = ptt.quadpolJackknife(Msub, nsub, z, os, ref, opts)
%
% The frame pass pools a segment's traces into one moment matrix as a
% trace-weighted sum of SUB-BLOCK moments (ptt.quadpolFrameTheta's heading
% sub-blocks, ~200 traces each). Those sub-blocks are the natural
% resampling unit: far enough apart that speckle, coregistration residual
% and small real lateral variation are all independent between them, so a
% delete-one jackknife over them measures the estimate's actual
% repeatability WITHOUT a noise model - which matters here, because the
% fitted data (a 90-azimuth coherence field synthesised from 16 moments,
% 60 m windows on a 15 m step, CRB weights) are far too correlated for a
% Jacobian covariance to be honest.
%
% Every replicate re-runs the full estimator on the pooled moments with
% one sub-block left out: it RE-VOTES the axis (constant mode) or refits
% every window's theta0 (free mode), then dlam, so the reported dlam error
% includes the axis's own uncertainty. The theta search is restricted to
% +-theta_half_deg around the full-data axis, which is what makes ~11
% replicate fits per segment affordable; a replicate whose axis wanted to
% leave that range would sit at its edge, and that is reported in
% J.n_edge so it is not silent.
%
% A RESTRICTED GRID IS NOT A WEAKER FIT. ptt.quadpolFabricLS scores an
% axis by q_theta, the cost contrast across whatever grid it was given
% normalised by the data power, and over a +-15 deg span that contrast is
% a small fraction of the same window's full-range value - not because
% the axis is worse determined but because the grid was narrowed to sit
% on the minimum. Gated on the same q_min as a full 0-180 search, free
% mode throws replicate values away for having been narrow-searched, and
% silently: n_edge does not cover it, and a window that loses enough
% replicates returns NaN from H_circ_se with nothing saying why.
% MEASURED in test_quadpol_uncertainty verdict D: 88% of replicate theta0
% values survive the raw gate on a clean synthetic, 100% once the scale
% is corrected, and the loss grows as the fabric weakens and q_theta
% approaches q_min. quadpolFabricLS therefore NORMALISES both contrasts -
% q_theta and the constant-orientation vote's curve range - for the span
% the grid actually reaches, instead of switching the gates off: whether
% a window has an axis at all was decided by the full fit passed in as
% `ref` and the replicate only measures how far it moves, but a replicate
% whose curve has genuinely gone flat must still be able to abstain.
%
% Inputs
%   Msub   cell of [Nt x 4 x 4] sub-block moment matrices (same frame,
%          same rotation convention, as summed into the full fit)
%   nsub   [1 x n] traces per sub-block (the pooling weights)
%   z      [Nt x 1] depth
%   os     the estimator options the FULL fit used (pedestal field,
%          theta_const, grids, ...); theta_grid is overridden here
%   ref    the full-data ptt.quadpolFabricLS output (for the axis to
%          search around and the window grid)
%   opts   theta_half_deg (15), theta_step_deg (3), min_rep (3; fewer
%          finite replicates at a window returns NaN there)
%
% Output struct J (per window centre, as ref.zw)
%   se_theta   [Nw x 1] standard error of theta0, radians, circular over
%              the doubled angle (constant mode: one value, repeated)
%   se_dlam    [Nw x 1] standard error of dlam
%   se_theta_c scalar, constant mode: SE of the held axis (NaN otherwise)
%   theta_rep  [Nw x n], dlam_rep [Nw x n] the replicate values
%   n          replicates run; n_edge  replicates whose axis sat at the
%              search edge (constant mode: of the held axis; free mode:
%              windows x replicates)
%
% Jackknife variance: (n-1)/n * sum_i (x_(i) - mean)^2, with the mean
% and the differences taken on the doubled-angle phasor for theta.
%
% See also ptt.quadpolFabricLS, ptt.quadpolFrameTheta.

if nargin < 6, opts = struct(); end
HALF = deg2rad(H_opt(opts, 'theta_half_deg', 15));
STEP = deg2rad(H_opt(opts, 'theta_step_deg', 3));
MIN_REP = H_opt(opts, 'min_rep', 3);

n = numel(Msub);
Nw = numel(ref.zw);
rel = (-HALF:STEP:HALF);
held = isfield(ref, 'theta_const') && isfinite(ref.theta_const) ...
  && isfield(os, 'theta_const') && os.theta_const;
if held
  grid = ref.theta_const + rel;               % one shared grid
else
  % per-window grids around each window's own axis. A window whose full
  % fit abstained on theta searches around the nearest reporting window's
  % axis (phasor-interpolated in depth): its replicates will mostly abstain
  % too, and a full-range search there would cost 5x for nothing.
  c = ref.theta0(:);
  okw = isfinite(c);
  if nnz(okw) < 2
    error('ptt:quadpolJackknife:noAxis', ...
      'the full fit reports no axis anywhere; nothing to resample around');
  end
  ph = interp1(ref.zw(okw), exp(2i * c(okw)), ...
    min(max(ref.zw(:), min(ref.zw(okw))), max(ref.zw(okw))), 'linear');
  c(~okw) = angle(ph(~okw)) / 2;
  grid = c + rel;                             % [Nw x Ng]
end
os.theta_grid = grid;

Nt = size(Msub{1}, 1);
Mtot = zeros(Nt, 4, 4);
wtot = 0;
for i = 1:n
  Mtot = Mtot + Msub{i} * nsub(i);
  wtot = wtot + nsub(i);
end

theta_rep = nan(Nw, n);
dlam_rep = nan(Nw, n);
th_c_rep = nan(1, n);
n_edge = 0;
for i = 1:n
  Mi = (Mtot - Msub{i} * nsub(i)) / max(wtot - nsub(i), eps);
  o = ptt.quadpolFabricLS(struct('M', Mi), z, os);
  if numel(o.zw) ~= Nw
    error('ptt:quadpolJackknife:windows', ...
      'replicate %d produced %d windows, full fit has %d', i, numel(o.zw), Nw);
  end
  theta_rep(:, i) = o.theta0;
  dlam_rep(:, i) = o.dlam;
  if held
    th_c_rep(i) = o.theta_const;
    if isfinite(o.theta_const) && abs(H_wrap(o.theta_const - ref.theta_const)) >= HALF - STEP/2
      n_edge = n_edge + 1;
    end
  else
    d = H_wrap(o.theta0 - ref.theta0(:));
    n_edge = n_edge + nnz(isfinite(d) & abs(d) >= HALF - STEP/2);
  end
end

se_theta = nan(Nw, 1);
se_dlam = nan(Nw, 1);
for w = 1:Nw
  se_theta(w) = H_circ_se(theta_rep(w, :), MIN_REP);
  se_dlam(w) = H_se(dlam_rep(w, :), MIN_REP);
end
se_theta_c = NaN;
if held
  se_theta_c = H_circ_se(th_c_rep, MIN_REP);
  se_theta(:) = se_theta_c;
end

J = struct('se_theta', se_theta, 'se_dlam', se_dlam, ...
  'se_theta_c', se_theta_c, 'theta_rep', theta_rep, 'dlam_rep', dlam_rep, ...
  'theta_c_rep', th_c_rep, 'n', n, 'n_edge', n_edge);
end

% -------------------------------------------------------------------------
function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function d = H_wrap(d)
% axis difference folded to (-pi/2, pi/2]
d = angle(exp(2i * d)) / 2;
end

function s = H_se(x, min_rep)
ok = isfinite(x);
m = nnz(ok);
if m < min_rep, s = NaN; return; end
x = x(ok);
s = sqrt((m - 1) / m * sum((x - mean(x)).^2));
end

function s = H_circ_se(th, min_rep)
ok = isfinite(th);
m = nnz(ok);
if m < min_rep, s = NaN; return; end
th = th(ok);
mu = angle(sum(exp(2i * th))) / 2;
d = angle(exp(2i * (th - mu))) / 2;
s = sqrt((m - 1) / m * sum(d.^2));
end

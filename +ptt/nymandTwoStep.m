function out = nymandTwoStep(obs, z, opts)
%NYMANDTWOSTEP Two-step fabric inversion: orientation then strength.
%
% out = ptt.nymandTwoStep(obs, z, opts)
%
% Implements the two-step scheme presented by Nymand for NEGIS, using the
% Fujita et al. (2006) depolarization-matrix model as the forward problem
% and simplifying it by assuming a CONSTANT EIGENFRAME with depth:
%
%   STEP 1  orientation and scattering ratio, from the co- and
%           cross-polarized POWER anomalies. The forward problem is Taylor
%           expanded and the first-order system solved by iterative least
%           squares - ptt.polarimetricInverse, which carries the analytic
%           reduction and is pinned against ptt.fujitaModel by test.
%           The scattering ratio is r = S22 / S11.
%
%   STEP 2  the horizontal contrast dlam(z), from the TRAVEL-TIME
%           anomalies, by linear maximum likelihood -
%           ptt.traveltimeFabricML. Exactly linear, so its posterior is
%           exact rather than linearised.
%
% WHY THE STEPS ARE COUPLED, AND WHY THIS ITERATES. The split is not
% clean in one direction: step 1's power anomalies depend on dlam through
% the accumulated birefringent phase, and step 2's travel-time anomalies
% must be measured between the EIGEN-polarizations, which needs theta0
% from step 1. Run once, each step uses a guess for what the other
% produces. So this runs them alternately until theta0 and dlam stop
% moving, which is the natural completion of the scheme rather than
% something the slides state. Set max_outer = 1 for the strict two-pass
% reading; the output records both, and n_outer says which you got.
%
% WHAT THIS BUYS OVER THE ESTIMATORS ALREADY HERE. ptt.quadpolFabricLS
% fits the complex coherence field and is the production path. This one
% never forms a coherence, so it survives where coherence does not - the
% Thwaites zones measured to have lost their fringes - and its strength
% step returns an exact posterior. It pays for that by ASSUMING a
% constant eigenframe, which is an assumption this repo has measured to
% be consequential. (An earlier claim here that the held axis sits 42.6
% degrees from the free-window axis was a units error - it is 1.7 degrees
% median, 8.8 at p90 - so the assumption is mild at Ridge A.) The held axis sits a median of 1.7 degrees from the
% free-window axis across 118 products, and that disagreement accounts
% for most of a one-third reduction in dlam. Do not read this estimator's
% strength against a free-axis one without accounting for that.
%
% obs fields:
%   .psi     [1 x Np] azimuths [rad]
%   .dP_hh   [Nz x Np] co-pol power anomaly [dB]
%   .dP_hv   [Nz x Np] cross-pol power anomaly [dB], optional
%   .dtau    [Nz x 1] two-way travel-time difference [ns] between the
%            eigen-polarizations. Optional if .dtau_fn is given.
%   .dtau_fn @(theta0) -> [Nz x 1], to re-measure the delay in the
%            eigenframe each outer iteration. Preferred: a dtau measured
%            once in the antenna frame is biased by cos(2*dtheta) and the
%            iteration cannot repair it.
%
% opts: passed through to both steps, plus
%   .max_outer (4)     outer iterations; 1 gives the strict two-pass form
%   .tol_theta (0.5)   convergence on theta0 [deg]
%   .tol_dlam  (0.002) convergence on median dlam
%
% out.theta0, out.r_z, out.dlam, out.sigma_dlam, out.resolution,
% out.chi2_dof, out.n_outer, out.converged, out.step1, out.step2, out.hist
%
% See also ptt.polarimetricInverse, ptt.traveltimeFabricML,
% ptt.orientationToFlow, ptt.quadpolFabricLS.
if nargin < 3, opts = struct(); end
z = z(:);
max_outer = H_opt(opts, 'max_outer', 4);
tol_th    = H_opt(opts, 'tol_theta', 0.5);
tol_dl    = H_opt(opts, 'tol_dlam', 0.002);

has_fn = isfield(obs, 'dtau_fn') && ~isempty(obs.dtau_fn);
if ~has_fn && (~isfield(obs, 'dtau') || isempty(obs.dtau))
  error('ptt:nymandTwoStep:noDelay', ...
    'supply obs.dtau (eigenframe travel-time difference) or obs.dtau_fn');
end

dl_scalar = H_opt(opts, 'dlam0', 0.03);   % step 1's starting contrast
th_prev = NaN; dl_prev = NaN; conv = false;
hist = struct('theta0', {}, 'dlam_med', {});
s1 = []; s2 = [];
for it = 1:max(1, max_outer)
  % --- step 1: orientation and scattering ratio from power anomalies
  o1 = opts;
  o1.dlam = dl_scalar;                    % scalar contrast for the phase term
  s1 = ptt.polarimetricInverse(obs, z, o1);
  th = s1.theta0;

  % --- step 2: contrast profile from travel-time anomalies
  if has_fn
    dtau = obs.dtau_fn(th);
  else
    dtau = obs.dtau;
  end
  s2 = ptt.traveltimeFabricML(dtau, z, opts);
  dl_med = median(s2.dlam, 'omitnan');
  if isfinite(dl_med), dl_scalar = dl_med; end

  hist(end+1) = struct('theta0', th, 'dlam_med', dl_med); %#ok<AGROW>
  % convergence on the doubled angle, never on a raw angle difference
  if isfinite(th_prev)
    dth = abs(rad2deg(angle(exp(2i*(th - th_prev)))/2));
    ddl = abs(dl_med - dl_prev);
    if dth < tol_th && (~isfinite(ddl) || ddl < tol_dl), conv = true; end
  end
  th_prev = th; dl_prev = dl_med;
  if conv, break; end
end

% dlam is returned on the depths step 2 actually solved; map back onto z
dl_z = nan(numel(z), 1); sg_z = nan(numel(z), 1); rs_z = zeros(numel(z), 1);
if ~isempty(s2) && ~isempty(s2.z)
  [tf, loc] = ismember(z, s2.z);
  dl_z(tf) = s2.dlam(loc(tf));
  sg_z(tf) = s2.sigma_dlam(loc(tf));
  rs_z(tf) = s2.resolution(loc(tf));
end

out = struct();
out.theta0 = s1.theta0;
out.r_z = s1.r_z;
out.dlam = dl_z;
out.sigma_dlam = sg_z;
out.resolution = rs_z;
out.chi2_dof = s2.chi2_dof;
out.n_outer = numel(hist);
out.converged = conv;
out.hist = hist;
out.step1 = s1;
out.step2 = s2;
end

function v = H_opt(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

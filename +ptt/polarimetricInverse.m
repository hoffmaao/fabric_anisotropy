function out = polarimetricInverse(obs, z, opts)
%POLARIMETRICINVERSE Joint inversion for fabric orientation and scattering ratio.
%
% out = ptt.polarimetricInverse(obs, z, opts)
%
% Solves for ONE fabric orientation theta0 and a depth profile of the
% anisotropic scattering ratio r(z) from the AZIMUTHAL POWER ANOMALIES,
% by iterative linearisation of an analytic forward model. Follows the
% unpublished Chapter 6 of Nymand (2024, PhD thesis, Niels Bohr Institute,
% "Radar investigation of the Northeast Greenland Ice Stream"), which
% formulates Fujita et al. (2006) under a depth-constant orientation and
% inverts it in the sense of Ershadi et al. (2022) but WITHOUT relying on
% the coherence method.
%
% WHY ADD THIS WHEN WE ALREADY HAVE TWO INVERSIONS. ptt.quadpolFabricLS
% fits the complex coherence field and returns theta0 and dlam; it says
% nothing about anisotropic SCATTERING, and it needs coherence to exist.
% ptt.ershadiInverse fits geometry per interval. Neither delivers a
% depth-resolved scattering ratio jointly with one column orientation,
% and r(z) is not a nuisance: Gerber et al. (2025) find anisotropic
% scattering dominates the azimuthal power response at 80% of their
% northeast Greenland points, so a fit that has no term for it is
% attributing that power to fabric.
%
% THE MODEL (Nymand eqs 6.8-6.10; derived under theta(z) = theta0 and zero
% conductivity, so the accumulated phases are real). With d = theta0 - psi
% the angle from the fabric axis to the synthesis azimuth, r_n the layer's
% scattering ratio S_y/S_x, and dpsi_n the ACCUMULATED two-way phase
% difference at layer n:
%
%   dP_hh = 20 log10 ( [cos^4 d + r^2 sin^4 d + (r/2) sin^2(2d) cos dpsi]
%                      / ([3 + 3 r^2 + 2 r cos dpsi] / 8) )
%   dP_hv = 20 log10 ( 8 cos^2 d sin^2 d )
%   phi   = atan2( r sin(dpsi) [1 - tan^4 d],
%                  r cos(dpsi) [1 + tan^4 d] + tan^2 d [1 + r^2] )
%
% Each power anomaly is normalised by its OWN azimuthal mean, which is
% what makes it independent of the signal magnitude - the property Nymand
% relies on instead of calibrating the channels.
%
% READ THIS BEFORE TRUSTING r FROM SYNTHESISED AZIMUTHS. That
% independence holds for a gain constant across the azimuths being
% compared. It is exact when azimuth is sampled PHYSICALLY, by driving
% different headings, because then each sample is one channel measured
% once and its gain divides out of its own mean. It is NOT exact when
% azimuth is SYNTHESISED from a quad-pol basis, because every synthesised
% azimuth is a different linear combination of four channels that may sit
% on different gains. And the failure is not benign here: r != 1 makes
% dP_hh 180-degree periodic while r = 1 makes it 90-degree periodic
% (Nymand, p73), and an H-versus-V gain imbalance under synthesis
% produces a 180-degree term of its own. So the scattering ratio and the
% channel imbalance are DEGENERATE in synthesised co-pol power. On this
% system the imbalance is 6-29 dB away from Ridge A (see
% ptt.equaliseChannels), so r from synthesised azimuths is an upper bound
% contaminated by instrument, not a measurement of the ice. Use
% opts.azimuth_source to record which case you are in; 'physical' is the
% one to believe.
%
% THE 90-DEGREE ALIAS APPLIES HERE TOO. The co-pol power anomaly cannot
% separate (theta0, r) from (theta0 + 90, 1/r) - the swap exchanges the
% HH and VV roles and the anomaly is normalised, so it absorbs the rest.
% Nymand notes the same trap (p76) and suggests resolving it from the sign
% of the phase gradient or from prior knowledge. ptt.polarimetricRatioInverse
% resolves it by including the coherence phase, which changes sign under
% the swap; prefer that function where the phase is available.
%
% WHAT THIS DOES NOT INVERT, AND WHY. dlam is NOT a free parameter.
% Nymand tried it and reported that linearised iteration "yields poor
% results", the forward problem being too non-linear in the eigenvalues
% for any realistic initial guess - which is the same conditioning we
% recorded independently in ptt.ershadiStrength, where theta is global in
% the observables and dlam is local. dlam enters here only through the
% accumulated phase dpsi(z), supplied by the caller:
%   - where the coherence survives, pass our own dlam(z) from
%     ptt.quadpolFabricLS. This is the extension over Nymand: at NEGIS he
%     had no usable coherence and had to set dlam = 0, and his eq 6.10
%     for the coherence phase is written but unused. At Ridge A, Taylor
%     Dome and the shallow half of most sites we DO have it.
%   - where it is gone, pass dlam = 0, which is Nymand's own device for
%     representing decoherence: uncorrelated HH and VV phases do not
%     interact, so the cos(dpsi) cross term averages away and the power
%     anomalies depend on r and theta0 alone.
%
%
% SIGN CONVENTION, PINNED BY TEST. Here d = theta0 - gamma, the SYNTHESIS
% sense, the same one ptt.quadpolAzimuth and ptt.quadpolFabricLS use, in
% which a fabric at azimuth theta0 puts its features at sweep index
% +theta0. ptt.fujitaModel deliberately uses the paper's R S R' sense, in
% which they appear at -theta0. The two agree only with theta negated:
% measured, matching the sign gives correlation 0.967 and median 0.23 dB
% against ptt.fujitaModel, mismatching it gives 0.11 and 1.96 dB. Both are
% asserted in test_polarimetric_inverse, because a silent flip here would
% negate every axis this code reports.
% Inputs
%   obs   struct of observables on a shared azimuth grid
%           .psi     [1 x Np] synthesis azimuths (radians)
%           .dP_hh   [Nz x Np] co-pol power anomaly (dB), or []
%           .dP_hv   [Nz x Np] cross-pol power anomaly (dB), or []
%           .phi     [Nz x Np] coherence phase (radians), or []
%           .w_hh/.w_hv/.w_phi  optional [Nz x Np] weights (default 1)
%         Rows that are all non-finite are dropped from the fit.
%   z     [Nz x 1] depth of each row (m)
%   opts  .dlam      [Nz x 1] eigenvalue difference profile, or a scalar
%                    (default 0 = Nymand's decohered case)
%         .fc        750e6
%         .theta0    initial orientation guess (rad; default from dP_hv
%                    if present, else 0)
%         .r0        initial scattering ratio (default 1)
%         .n_r       number of r layers (default 12; r is piecewise
%                    linear in depth between them, which is the
%                    regularisation that keeps it from chasing nodes -
%                    Nymand reports his r(z) overfits the co-pol nodes)
%         .eta       regularisation weight on the second difference of r
%                    (default 1e-2)
%         .max_iter  (30), .step (1.0 initial step length)
%         .clip_db   (-25) floor applied to both model and data power
%                    anomalies, so true nulls do not dominate the misfit
%         .theta0_scan (true) coarse 1-D scan over theta0 before
%                    linearising, which removes the pi/2 trap
%         .use       cellstr subset of {'hh','hv','phi'} (default: all
%                    present)
%         .azimuth_source 'synthetic' (default) or 'physical' - recorded
%                    in the output, and gates the r warning above
%
% Output
%   out.theta0, out.r_z [Nz x 1], out.r_nodes, out.z_nodes
%   out.loss (per iteration), out.n_iter, out.converged
%   out.pred  struct of predicted observables at the solution
%   out.resid struct of residuals, out.rms per observable
%   out.azimuth_source, out.dlam_used, out.warn (cellstr)
%
% See also ptt.quadpolFabricLS, ptt.ershadiInverse, ptt.fujitaModel,
%   ptt.equaliseChannels, ptt.quadpolFabricPower.

if nargin < 3, opts = struct(); end
psi = obs.psi(:).';
Np = numel(psi);
z = z(:);
Nz = numel(z);
fc = H_opt(opts, 'fc', 750e6);
n_r = H_opt(opts, 'n_r', 12);
eta = H_opt(opts, 'eta', 1e-2);
max_iter = H_opt(opts, 'max_iter', 30);
step0 = H_opt(opts, 'step', 1.0);
% NULL FLOOR. The co-pol anomaly is a log, and at r = 1 its nodes are true
% zeros: the synthetic column below reaches -56 dB. Least squares in dB
% then spends almost all its misfit on the few cells nearest those nodes,
% which are the least informative (their exact depth is set by where a
% ratio crosses zero) and the most non-linear. Real data never go there -
% they stop at the noise floor - so both model and data are floored here
% before differencing, which is a statement that we do not claim to know
% how deep a null is once it is far below anything measurable.
clip_db = H_opt(opts, 'clip_db', -25);
% theta0 is periodic and the loss in it is not convex, so a linearisation
% started far away converges to the axis 90 deg across (Nymand p76 warns
% of exactly this). A coarse 1-D scan costs a few forward evaluations and
% removes the trap; set false to reproduce the bare linearisation.
% Default: scan UNLESS the caller pinned theta0, in which case pinning it
% and then overriding it would make the argument meaningless.
theta0_scan = H_opt(opts, 'theta0_scan', ...
  ~(isfield(opts, 'theta0') && ~isempty(opts.theta0)));
az_src = H_opt(opts, 'azimuth_source', 'synthetic');
warn = {};

% --- accumulated two-way phase difference from the supplied dlam
C = ptt.constants();
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c * 1e9);   % rad/m per unit dlam
dlam = H_opt(opts, 'dlam', 0);
if isscalar(dlam), dlam = dlam * ones(Nz, 1); else, dlam = dlam(:); end
if numel(dlam) ~= Nz
  error('ptt:polarimetricInverse:dlam', ...
    'opts.dlam has %d entries for %d depths', numel(dlam), Nz);
end
dz = [0; diff(z)];
dpsi = cumsum(max(dlam, 0) .* dz) * gpd;    % [Nz x 1], zero at the top

% --- which observables
have = {};
if isfield(obs, 'dP_hh') && ~isempty(obs.dP_hh), have{end+1} = 'hh'; end
if isfield(obs, 'dP_hv') && ~isempty(obs.dP_hv), have{end+1} = 'hv'; end
if isfield(obs, 'phi')   && ~isempty(obs.phi),   have{end+1} = 'phi'; end
use = H_opt(opts, 'use', have);
use = intersect(use, have, 'stable');
if isempty(use)
  error('ptt:polarimetricInverse:noObs', 'no usable observable supplied');
end
if any(strcmp(use, 'hh')) && strcmpi(az_src, 'synthetic')
  warn{end+1} = ['r from SYNTHESISED azimuths is degenerate with the ' ...
    'H/V channel gain imbalance (both are 180-deg periodic); treat it ' ...
    'as an upper bound unless the channels were equalised first'];
end

% --- pack the data vector and weights
[d_obs, w_obs, idx] = H_pack(obs, use, Nz, Np, clip_db);
if isempty(d_obs)
  error('ptt:polarimetricInverse:empty', 'every observable row is non-finite');
end

% --- model: m = [theta0; r_nodes]
z_nodes = linspace(z(1), z(end), n_r).';
th0 = H_opt(opts, 'theta0', []);
if isempty(th0)
  th0 = H_theta_from_hv(obs, psi, use);
end
r0 = H_opt(opts, 'r0', 1);
m = [th0; r0 * ones(n_r, 1)];
if theta0_scan
  % coarse scan at the initial r, on the doubled angle's natural period
  th_try = th0 + (-90:5:85) * pi/180;
  best = inf; th_best = th0;
  for q = 1:numel(th_try)
    mt = m; mt(1) = th_try(q);
    mt_pred = H_pack_pred(H_forward(mt, psi, dpsi, z, z_nodes, Nz, Np), use, idx, clip_db);
    Lq = sum(w_obs .* (d_obs - mt_pred).^2);
    if Lq < best, best = Lq; th_best = th_try(q); end
  end
  m(1) = th_best;
end

% --- second-difference smoothness on the r nodes only
Gam = zeros(max(n_r-2, 0), numel(m));
for k = 1:n_r-2
  Gam(k, 1 + (k:k+2)) = [1 -2 1];
end
GtG_reg = eta^2 * (Gam.' * Gam);

fwd = @(mm) H_pack_pred(H_forward(mm, psi, dpsi, z, z_nodes, Nz, Np), use, idx, clip_db);
loss = @(mm) sum(w_obs .* (d_obs - fwd(mm)).^2) + eta^2 * sum((Gam*mm).^2);

L = nan(max_iter+1, 1);
L(1) = loss(m);
converged = false;
it = 0;
for it = 1:max_iter
  G = H_jac(fwd, m, numel(d_obs));
  W = spdiags(w_obs, 0, numel(w_obs), numel(w_obs));
  A = G.' * W * G + GtG_reg;
  b = G.' * W * (d_obs - fwd(m)) - GtG_reg * m;
  dm = A \ b;
  if ~all(isfinite(dm)), warn{end+1} = 'singular normal equations; stopped'; break; end %#ok<AGROW>
  % step halving: never accept an update that raises the loss
  a = step0; ok = false;
  for h = 1:20
    mt = m + a * dm;
    mt(2:end) = max(mt(2:end), 1e-3);     % scattering ratio stays positive
    Lt = loss(mt);
    if Lt < L(it), m = mt; L(it+1) = Lt; ok = true; break; end
    a = a / 2;
  end
  if ~ok, converged = true; L(it+1) = L(it); break; end
  if abs(L(it) - L(it+1)) / max(L(it+1), realmin) < 1e-4
    converged = true; break;
  end
end
L = L(1:min(it+1, numel(L)));

r_nodes = m(2:end);
r_z = interp1(z_nodes, r_nodes, z, 'linear', 'extrap');
P = H_forward(m, psi, dpsi, z, z_nodes, Nz, Np);

out = struct('theta0', mod(m(1), pi), 'r_nodes', r_nodes, 'z_nodes', z_nodes, ...
  'r_z', r_z, 'loss', L, 'n_iter', it, 'converged', converged, ...
  'pred', P, 'dlam_used', dlam, 'dpsi', dpsi, 'psi', psi, ...
  'azimuth_source', az_src, 'used', {use}, 'warn', {warn}, 'clip_db', clip_db);
out.resid = struct(); out.rms = struct();
for k = 1:numel(use)
  f = use{k};
  switch f
    case 'hh', O = max(obs.dP_hh, clip_db); M = max(P.dP_hh, clip_db);
    case 'hv', O = max(obs.dP_hv, clip_db); M = max(P.dP_hv, clip_db);
    case 'phi', O = obs.phi;  M = P.phi;
  end
  R = O - M;
  if strcmp(f, 'phi'), R = angle(exp(1i*R)); end
  out.resid.(f) = R;
  out.rms.(f) = sqrt(mean(R(isfinite(R)).^2));
end
for k = 1:numel(warn), fprintf('polarimetricInverse WARNING: %s\n', warn{k}); end
end

% -------------------------------------------------------------------------
function P = H_forward(m, psi, dpsi, z, z_nodes, Nz, Np)
%H_FORWARD Nymand eqs 6.8-6.10 on the [Nz x Np] grid.
th0 = m(1);
r = interp1(z_nodes, m(2:end), z, 'linear', 'extrap');   % [Nz x 1]
r = max(r, 1e-6);
d = th0 - psi;                       % [1 x Np]
c2 = cos(d).^2; s2 = sin(d).^2;
c4 = c2.^2;     s4 = s2.^2;
s2d2 = (2*cos(d).*sin(d)).^2;        % sin^2(2d)
cdp = cos(dpsi); sdp = sin(dpsi);    % [Nz x 1]
num = c4 + (r.^2) .* s4 + 0.5 * r .* s2d2 .* cdp;
den = (3 + 3*r.^2 + 2*r .* cdp) / 8;
P.dP_hh = 20 * log10(max(num, realmin) ./ max(den, realmin));
P.dP_hv = repmat(20 * log10(max(8 * c2 .* s2, realmin)), Nz, 1);
t2 = tan(d).^2; t4 = t2.^2;
P.phi = atan2(r .* sdp .* (1 - t4), ...
              r .* cdp .* (1 + t4) + t2 .* (1 + r.^2));
end

function [d, w, idx] = H_pack(obs, use, Nz, Np, clip_db)
d = []; w = []; idx = struct();
for k = 1:numel(use)
  f = use{k};
  switch f
    case 'hh', O = max(obs.dP_hh, clip_db); wf = H_w(obs, 'w_hh', Nz, Np);
    case 'hv', O = max(obs.dP_hv, clip_db); wf = H_w(obs, 'w_hv', Nz, Np);
    case 'phi', O = obs.phi;  wf = H_w(obs, 'w_phi', Nz, Np);
  end
  O = O(:); wf = wf(:);
  good = isfinite(O) & isfinite(wf) & wf > 0;
  idx.(f) = good;
  d = [d; O(good)]; %#ok<AGROW>
  w = [w; wf(good)]; %#ok<AGROW>
end
end

function v = H_pack_pred(P, use, idx, clip_db)
v = [];
for k = 1:numel(use)
  f = use{k};
  switch f
    case 'hh', M = max(P.dP_hh, clip_db);
    case 'hv', M = max(P.dP_hv, clip_db);
    case 'phi', M = P.phi;
  end
  M = M(:);
  v = [v; M(idx.(f))]; %#ok<AGROW>
end
end

function w = H_w(obs, name, Nz, Np)
if isfield(obs, name) && ~isempty(obs.(name)), w = obs.(name); else, w = ones(Nz, Np); end
end

function G = H_jac(fwd, m, nd)
% Numeric Jacobian. Nymand derives these analytically (his appendix E);
% finite differences are used here because the forward model is cheap on
% our grids and an analytic Jacobian that disagrees with the forward is a
% silent bias, whereas this cannot disagree by construction.
np = numel(m);
G = zeros(nd, np);
f0 = fwd(m);
for j = 1:np
  h = max(1e-6, 1e-4 * abs(m(j)));
  mp = m; mp(j) = mp(j) + h;
  G(:, j) = (fwd(mp) - f0) / h;
end
G(~isfinite(G)) = 0;
end

function th = H_theta_from_hv(obs, psi, use)
% The cross-pol power anomaly depends on theta0 ALONE (Nymand p73), so
% where it is supplied it gives the initial guess for free: its minima sit
% on the fabric axes. Guarded, because on THIS system the cross-pol sits
% on a -3.6 dB instrument pedestal that antenna-locks exactly this
% minimum, which is the whole reason ptt.quadpolFabricLS exists - so the
% guess is a starting point only, never a result.
th = 0;
if ~any(strcmp(use, 'hv')) || ~isfield(obs, 'dP_hv') || isempty(obs.dP_hv)
  return;
end
p = mean(obs.dP_hv, 1, 'omitnan');
if all(~isfinite(p)), return; end
% axes are where the cross-pol vanishes: pick the deepest minimum on the
% doubled angle so the two axes 90 deg apart do not fight
[~, i] = min(p);
th = psi(i);
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

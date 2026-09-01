function out = ershadiInverse(fr, z, opts)
%ERSHADIINVERSE Ershadi et al. (2022) Sect. 3.5 inversion for theta and r.
%
% out = ptt.ershadiInverse(fr, z, opts)
%
% The constrained non-linear least-squares step of Ershadi et al. (2022,
% The Cryosphere 16, 1719-1739): given the observables ptt.ershadiFabric
% extracts - the power anomalies dP_HH and dP_HV (their eq. 12), the
% coherence phase phi_HHVV (eqs. 7-8) - fit the layered Fujita forward
% model (ptt.fujitaModel, their eq. 5) for a piecewise-constant profile of
% the fabric orientation theta_i and the reflection ratio r_i, while the
% horizontal anisotropy profile is ACCEPTED from the phase-gradient
% estimate and not re-fit ("we optimize theta and r for all depth
% intervals, while at this stage we accept the estimated dlam0", 3.5.4).
%
% FAITHFUL CHOICES, with their sources:
%   - Piecewise-constant intervals (their EDML parameterization, 3.5.2;
%     50 m default). The Legendre alternative is not implemented.
%   - Initial guess (3.5.3): theta0 from ptt.ershadiFabric's theta output,
%     which is exactly their recipe (dP_HV minima, disambiguated by the
%     phase polarity and the sign of Psi); r0 = 0 dB.
%   - Cost (3.5.4, eqs. 15-18): J = l1*J_phi + l2*J_dPHH + l3*J_dPHV with
%     each field STANDARDIZED - here by the observation's own mean and
%     standard deviation over finite entries, the same affine map applied
%     to the model, so the three terms are commensurate and a misfit
%     cannot hide in an overall scale.
%   - Weights (Table 3) are 0/1 selectors set PER PARAMETER: theta is fit
%     with w_theta (Dome C row [1 0 0]: coherence phase only; EDML row
%     [0 1 0]: dP_HH only, for strong anisotropy where the phase misfit
%     "was not applicable"), r with w_r (both sites [0 1 0]).
%   - Bounds (3.5.4): 0 < theta_i < pi, -30 dB < r_i < 30 dB, enforced by
%     fmincon bound constraints (they used log-barriers inside the cost
%     with fmincon; interior-point bounds are the same mechanism, owned by
%     the solver instead of hand-rolled).
%   - r = Gamma_y/Gamma_x with r_dB = 20*log10(r) - the amplitude
%     convention their eq. (13) fixes (see ptt.fujitaModel).
%
% Also returned: the eq. (13) ANALYTIC estimate r = 1/tan^2(AD/2) at every
% depth where two co-polarization nodes are resolvable in dP_HH, as an
% optimization-free cross-check of the fitted profile.
%
% Inputs
%   fr    output struct of ptt.ershadiFabric (needs psi, dP_hh, dP_hv,
%         phi, Cmag, theta, dlam)
%   z     depth axis (m), same rows as the fields in fr
%   opts  .interval_m (50)   piecewise interval length
%         .z_fit ([200 z(end)]) band the misfits are evaluated over
%         .w_theta ([1 0 0])  Table-3 row for the theta stage
%         .w_r     ([0 1 0])  Table-3 row for the r stage
%         .fit_decim (4)      depth-row decimation inside the cost
%         .fit_psi_decim (4)  azimuth-column decimation inside the cost
%         .fc, .eps_perp, .deps, .win_m   passed to ptt.fujitaModel
%         .max_iter (150)     fmincon iteration cap per stage
%
% Output
%   .z_int, .theta_int, .r_db_int   per-interval fitted values
%   .theta, .r_db                   the same, expanded to the z rows
%   .dlam_int                       the accepted dlam per interval
%   .r13_db                         eq.-(13) analytic r per depth row (NaN
%                                   where nodes are not resolvable)
%   .J0, .J                         cost before/after, per stage
%   .exitflag                       fmincon exit flags [theta, r]
%
% See also ptt.fujitaModel, ptt.ershadiFabric.

if nargin < 3, opts = struct(); end
int_m = H_opt(opts, 'interval_m', 50);
zfit = H_opt(opts, 'z_fit', [200, z(end)]);
w_th = H_opt(opts, 'w_theta', [1 0 0]);
w_r = H_opt(opts, 'w_r', [0 1 0]);
dec_z = H_opt(opts, 'fit_decim', 4);
dec_p = H_opt(opts, 'fit_psi_decim', 4);
max_iter = H_opt(opts, 'max_iter', 150);

fwd = struct('fc', H_opt(opts, 'fc', 750e6), ...
  'eps_perp', H_opt(opts, 'eps_perp', 3.15), ...
  'deps', H_opt(opts, 'deps', 0.034), ...
  'win_m', H_opt(opts, 'win_m', 30));

z = z(:);
band = z >= zfit(1) & z <= zfit(2);
if ~any(band)
  error('ptt:ershadiInverse:band', 'z_fit selects no rows');
end

% --- intervals over the fit band; the layer stack starts at the SURFACE
% (isotropic-ish shallow rows still propagate phase), so intervals above
% the band exist too, carrying the accepted dlam and the initial theta.
z0 = 0; z1 = min(zfit(2), z(end));
edges = (z0:int_m:z1).';
if edges(end) < z1, edges(end+1) = z1; end
Nint = numel(edges) - 1;
fitset = find(edges(1:end-1) + int_m/2 >= zfit(1));  % intervals we optimize
zc_int = 0.5 * (edges(1:end-1) + edges(2:end));

% --- accepted dlam and initial theta per interval (3.5.3/3.5.4)
dlam_int = zeros(Nint, 1);
th0_int = zeros(Nint, 1);
for k = 1:Nint
  m = z >= edges(k) & z < edges(k+1);
  d = fr.dlam(m); d = d(isfinite(d));
  if ~isempty(d), dlam_int(k) = median(d); end
  t = fr.theta(m); t = t(isfinite(t));
  if ~isempty(t)
    th0_int(k) = mod(angle(mean(exp(2i * t))) / 2, pi);
  elseif k > 1
    th0_int(k) = th0_int(k-1);
  end
end
r0_int = zeros(Nint, 1);                 % "The initial guess for r0dB is zero"

% --- observation fields, decimated and standardized once
iz = find(band); iz = iz(1:dec_z:end);
ip = 1:dec_p:numel(fr.psi);
psi_fit = fr.psi(ip);
obs.phi = fr.phi(iz, ip);
obs.hh = fr.dP_hh(iz, ip);
obs.hv = fr.dP_hv(iz, ip);
zsub = z(iz);
std_of = @(A) deal_std(A);
[mu1, s1] = std_of(obs.phi); [mu2, s2] = std_of(obs.hh); [mu3, s3] = std_of(obs.hv);
nrm = {@(A) (A - mu1)/s1, @(A) (A - mu2)/s2, @(A) (A - mu3)/s3};
obs_n = {nrm{1}(obs.phi), nrm{2}(obs.hh), nrm{3}(obs.hv)};

cost = @(th_all, r_all, w) H_cost(th_all, r_all, dlam_int, edges, ...
  zsub, psi_fit, fwd, obs_n, nrm, w);

oopt = optimoptions('fmincon', 'Display', 'off', 'Algorithm', ...
  'interior-point', 'MaxIterations', max_iter, ...
  'MaxFunctionEvaluations', 200 * numel(fitset));

th = th0_int; rdb = r0_int;
J0 = [cost(th, rdb, w_th), cost(th, rdb, w_r)];
ex = [0 0];

% --- stage 1: theta intervals against w_theta (Table 3, theta row)
f1 = @(p) cost(H_place(th, fitset, p), rdb, w_th);
[p1, J1, ex(1)] = fmincon(f1, th(fitset), [], [], [], [], ...
  zeros(numel(fitset), 1), pi * ones(numel(fitset), 1), [], oopt);
th = H_place(th, fitset, p1);

% --- stage 2: r intervals against w_r (Table 3, r row), theta held
f2 = @(p) cost(th, H_place(rdb, fitset, p), w_r);
[p2, J2, ex(2)] = fmincon(f2, rdb(fitset), [], [], [], [], ...
  -30 * ones(numel(fitset), 1), 30 * ones(numel(fitset), 1), [], oopt);
rdb = H_place(rdb, fitset, p2);

% --- eq. (13) analytic r from the co-pol node angular distance, where two
% nodes are resolvable in a depth row of dP_HH; optimization-free check.
% The node pair defines two arcs (AD and pi - AD) and eq. (13) applied to
% the wrong one returns 1/r - the axis says which arc is which, since the
% nodes straddle it. The FITTED theta (sweep frame: -theta) disambiguates.
r13 = nan(numel(z), 1);
for i = find(band).'
  k = find(z(i) >= edges(1:end-1) & z(i) < edges(2:end), 1);
  if isempty(k), continue; end
  r13(i) = H_r13(fr.dP_hh(i, :), fr.psi, mod(-th(k), pi));
end

% expand to rows
th_row = nan(numel(z), 1); r_row = nan(numel(z), 1);
for k = 1:Nint
  m = z >= edges(k) & z < edges(k+1);
  th_row(m) = th(k); r_row(m) = rdb(k);
end

out = struct('z_int', zc_int, 'theta_int', th, 'r_db_int', rdb, ...
  'dlam_int', dlam_int, 'theta', th_row, 'r_db', r_row, 'r13_db', r13, ...
  'J0', J0, 'J', [J1, J2], 'exitflag', ex, 'edges', edges, ...
  'fitset', fitset);

end

% ---------------------------------------------------------------- helpers
function J = H_cost(th_all, r_all, dlam_int, edges, zsub, psi, fwd, obs_n, nrm, w)
NL = numel(th_all);
layers = struct('top_m', num2cell(edges(1:NL).'), ...
  'dlam', num2cell(dlam_int.'), 'theta', num2cell(th_all.'), ...
  'r_db', num2cell(r_all.'));
mod_ = ptt.fujitaModel(layers, zsub, psi, fwd);
J = 0;
flds = {mod_.phi, mod_.dP_hh, mod_.dP_hv};
for t = 1:3
  if w(t) == 0, continue; end
  d = nrm{t}(flds{t}) - obs_n{t};
  ok = isfinite(d);
  J = J + w(t) * sum(d(ok).^2) / max(nnz(ok), 1);
end
end

function v = H_place(v, idx, p)
v(idx) = p;
end

function [mu, s] = deal_std(A)
a = A(isfinite(A));
mu = mean(a); s = std(a); if s <= 0, s = 1; end
end

function rdb = H_r13(row, psi, ax)
%H_R13 Eq. (13): r = 1/tan^2(AD/2) from the two co-pol node azimuths.
% Nodes = local minima (on the periodic row) at least 3 dB below the row
% median; exactly two are required, else NaN. The pair splits the periodic
% azimuth into two arcs; AD is the one CONTAINING THE AXIS `ax` (sweep
% frame), because the nodes straddle the axis - eq. (13) on the other arc
% would return 1/r.
rdb = NaN;
ok = isfinite(row);
if nnz(ok) < 8, return; end
med = median(row(ok));
prv = circshift(row(:), 1).'; nxt = circshift(row(:), -1).';
ismin = row < prv & row < nxt & row < med - 3;
idx = find(ismin & ok);
if numel(idx) ~= 2, return; end
p1 = psi(idx(1)); p2 = psi(idx(2));           % p1 < p2 on [0, pi)
in_arc = ax >= p1 && ax < p2;                 % axis inside [p1, p2)?
ad = p2 - p1;
if ~in_arc, ad = pi - ad; end                 % take the arc holding the axis
ad = max(min(ad, pi - 1e-6), 1e-6);
rdb = 20 * log10(1 / tan(ad / 2)^2);
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function out = fabricGLS(obs, z, opts)
%FABRICGLS Generalized least squares inversion of the Fujita model, with uncertainty.
%
% out = ptt.fabricGLS(obs, z, opts)
%
% Estimates the fabric orientation theta(z), the horizontal eigenvalue
% difference dlam(z) and the anisotropic scattering ratio r(z) as DEPTH
% PROFILES, by generalized least squares against the full layered
% Fujita et al. (2006) forward model (ptt.fujitaModel), and returns a
% posterior covariance rather than a point estimate alone.
%
% WHAT THIS CHANGES RELATIVE TO ptt.polarimetricRatioInverse AND
% ptt.polarimetricInverse. Both of those impose theta(z) = theta0 as a
% HARD constraint, which is what buys their analytic forward model. That
% is a strong assumption dressed as a method: where the axis does rotate -
% a shear margin, the deep part of an ice stream - they cannot say so,
% they can only fit badly. Here theta(z) is a free profile and the
% constant-orientation case is recovered by making its PRIOR tight. The
% assumption becomes a number the user sets and the data can argue with,
% and out.resolution says how much of each estimate came from the data
% rather than from that prior.
%
% THE FORMULATION. With d the observations, g(m) the forward model, C_d
% the data covariance and (m_prior, C_m) the prior:
%
%   L(m) = (d - g(m))' inv(C_d) (d - g(m)) + (m - m_p)' inv(C_m) (m - m_p)
%
% minimised by iterative linearisation (Tarantola). At each step, with G
% the Jacobian of g at the current m,
%
%   dm   = (G' inv(C_d) G + inv(C_m))^-1 [G' inv(C_d) (d - g(m)) - inv(C_m)(m - m_p)]
%   C_M  = (G' inv(C_d) G + inv(C_m))^-1                    posterior covariance
%   R    = C_M G' inv(C_d) G                                resolution
%
% C_M is the linearised (Gaussian) posterior at the solution. Its diagonal
% gives the one-sigma uncertainty on every parameter, propagated to depth
% profiles in out.sigma_*. R is the resolution matrix: a diagonal entry
% near 1 means the data determined that parameter, near 0 means the prior
% did, and reporting an estimate whose resolution is near zero as a
% measurement is the main way this kind of inversion misleads.
%
% WHY GENERALIZED AND NOT ORDINARY LEAST SQUARES. The observables are not
% on one scale and are not equally certain. A dB power anomaly from N
% looks has variance (10/ln10)^2 psi'(N); a coherence phase from the same
% looks has variance (1-|C|^2)/(2 N |C|^2), which varies by orders of
% magnitude down a column as the coherence decays. Ordinary least squares
% would weight a decohered deep phase sample as heavily as a clean shallow
% one. C_d is where that belongs - and it also removes the arbitrary
% observable weight that ptt.polarimetricRatioInverse has to tune by hand
% (its w_phi, fixed at 2 dB per radian by experiment rather than by
% measurement uncertainty).
%
% CORRELATED ERRORS, HANDLED HONESTLY. Neighbouring depth samples are not
% independent: range multilooking correlates them, and so does the
% coherence window. This code does NOT build a full banded C_d; it takes
% an EFFECTIVE look count (opts.n_looks, the number of INDEPENDENT looks)
% and forms a diagonal C_d from it. That is the right first-order fix,
% because the failure mode of ignoring correlation is a posterior that is
% too tight by roughly the square root of the correlation length in
% samples, and inflating the variance corrects the size even though it
% does not capture the shape. The posterior is therefore CHECKED against
% Monte Carlo scatter in test_fabric_gls rather than trusted; a linearised
% posterior on a problem this non-linear can be optimistic no matter how
% carefully C_d is built.
%
% Inputs
%   obs  .psi    [1 x Np] azimuths of the observations (radians)
%        any of:
%        .dP_hh, .dP_hv  [Nz x Np] power anomalies (dB)
%        .phi            [Nz x Np] co-pol coherence phase (radians)
%        .Cmag           [Nz x Np] coherence magnitude, used to build the
%                        phase variance when .sigma_phi is not given
%        .sigma_*        optional explicit one-sigma for any observable
%   z    [Nz x 1] depth (m)
%   opts .fc (750e6), .win_m (30)
%        .win_power_m (0), .dz_model ([])   how the observables' powers
%                        were windowed and at what native sampling, passed
%                        to ptt.fujitaModel so the model windows the same
%                        way on a decimated grid (issue #29)
%        .n_layer (16)   layers the profiles are parametrised on
%        .n_looks (20)   INDEPENDENT looks behind each sample
%        .n_indep_psi    independent azimuths behind the Np columns. Use 4
%                        for channel-synthesized azimuths, Np (default)
%                        for physically rotated antennas. See REDUNDANCY.
%        .n_indep_z      independent depth samples behind the Nz rows,
%                        default Nz. Roughly range-window / sample-spacing.
%        .use            cellstr of observables (default: all supplied)
%        .theta_prior    [mean_rad sigma_rad corr_len_m], default
%                        [0 pi/2 1e4]: uninformative in value, and a
%                        10 km correlation length, i.e. "nearly constant
%                        with depth unless the data say otherwise".
%                        Shorten corr_len to allow rotation; set sigma
%                        small to force a constant axis.
%        .dlam_prior     [mean sigma corr_len_m], default [0.05 0.08 300]
%        .rdb_prior      [mean sigma corr_len_m], default [0 4 300]
%        .clip_db (-30)  floor on the dB anomalies, whose nodes are
%                        true zeros
%        .max_iter (25), .n_start (5)   solves run from the best starts
%        .n_theta_scan (6), .dlam_scan  the start grid; the loss is
%                        multi-modal in BOTH and starting from the prior
%                        mean lands in a local minimum
%
% Output
%   .theta_z, .dlam_z, .rdb_z          profiles on z
%   .sigma_theta_z, .sigma_dlam_z, .sigma_rdb_z   one-sigma, same grid
%   .C_M            posterior covariance over the parameter vector
%   .resolution     diagonal of the resolution matrix, per parameter
%   .chi2_dof       reduced chi-square; far from 1 means C_d is wrong
%   .loss, .converged, .n_iter, .z_layer, .pred, .m, .ip
%
% See also ptt.fujitaModel, ptt.polarimetricRatioInverse,
%   ptt.quadpolJackknife (the resampling uncertainty used elsewhere here).

if nargin < 3, opts = struct(); end
z = z(:); Nz = numel(z);
psi = obs.psi(:).';
Np = numel(psi);
fc = H_opt(opts, 'fc', 750e6);
win_m = H_opt(opts, 'win_m', 30);
% how the OBSERVABLES were formed, handed to the forward model so it forms
% its own the same way (ptt.fujitaModel: win_power_m, dz_model)
fwd_opts = struct('fc', fc, 'win_m', win_m, ...
  'win_power_m', H_opt(opts, 'win_power_m', 0), ...
  'dz_model', H_opt(opts, 'dz_model', []));
nL = H_opt(opts, 'n_layer', 16);
nlook = H_opt(opts, 'n_looks', 20);
max_iter = H_opt(opts, 'max_iter', 25);
n_start = H_opt(opts, 'n_start', 5);

th_pr = H_opt(opts, 'theta_prior', [0, pi/2, 1e4]);
dl_pr = H_opt(opts, 'dlam_prior', [0.05, 0.08, 300]);
rd_pr = H_opt(opts, 'rdb_prior',  [0, 4, 300]);

cand = {'dP_hh', 'dP_hv', 'phi'};
have = cand(cellfun(@(f) isfield(obs, f) && ~isempty(obs.(f)), cand));
use = H_opt(opts, 'use', have);
use = intersect(use, have, 'stable');
if isempty(use), error('ptt:fabricGLS:noObs', 'no usable observable supplied'); end

% --- layer grid the profiles live on
z_layer = linspace(0, max(z), nL+1).';
z_layer = z_layer(1:nL);                      % layer tops, first at 0
z_mid = z_layer + [diff(z_layer); max(z)-z_layer(end)]/2;

% --- parameter vector and its prior
ip = struct('th', 1:nL, 'dl', nL+(1:nL), 'rd', 2*nL+(1:nL));
np = 3*nL;
m_p = [th_pr(1)*ones(nL,1); dl_pr(1)*ones(nL,1); rd_pr(1)*ones(nL,1)];
Cm = blkdiag(H_gp(z_mid, th_pr(2), th_pr(3)), ...
             H_gp(z_mid, dl_pr(2), dl_pr(3)), ...
             H_gp(z_mid, rd_pr(2), rd_pr(3)));
Cmi = H_inv(Cm);
% finite-difference steps on the prior scale, one per parameter
hstep = 1e-3 * [th_pr(2)*ones(nL,1); dl_pr(2)*ones(nL,1); rd_pr(2)*ones(nL,1)];
hstep = max(hstep, 1e-9);

% --- data vector and its covariance (diagonal, from the look count)
% REDUNDANCY. A diagonal C_d is right only if the samples are independent,
% and the expectation here was that they are not: synthesized azimuths are
% linear combinations of the SAME four measured channels, so Np azimuth
% columns should carry about four independent numbers rather than Np, and a
% diagonal C_d counting each repeat as fresh evidence should come out too
% tight. opts.n_indep_psi was added to correct that.
%
% VERDICT 5 OF test_fabric_gls MEASURED THE OPPOSITE, and the option
% therefore stays OFF. Under noise confined to a rank-4 azimuthal subspace
% the diagonal posterior is not optimistic at all - the reported/empirical
% ratio comes out near 1.15, i.e. mildly CONSERVATIVE - while n_indep_psi =
% 4 multiplies every sigma by sqrt(Np/4) (2.1x at the default 18 azimuths)
% and overshoots to about 2.4. Azimuth redundancy is simply not what
% threatens these error bars. Do not turn n_indep_psi on for
% channel-synthesized azimuths; verdict 5 exists to record the measurement
% that stops exactly that reflex, and it prints both ratios every run.
%
% The option is kept, and documented, for data that genuinely does repeat
% samples. Rather than build a dense C_d for a matrix this large, the
% variances are inflated by the redundancy ratio, which is exactly
% equivalent for the posterior scale: set n_indep_psi to the number of
% INDEPENDENT azimuths and n_indep_z to the number of independent depth
% samples. Defaults assume independence, so an unset option never silently
% inflates anything.
%
% This is a SEPARATE question from how far the posterior can be trusted at
% all; see the LOWER BOUNDS note at the solution below, which measures a
% different thing and still stands.
clip_db = H_opt(opts, 'clip_db', -30);
n_ipsi = H_opt(opts, 'n_indep_psi', Np);
n_iz   = H_opt(opts, 'n_indep_z', Nz);
redund = sqrt(max(Np/max(n_ipsi,1), 1) * max(Nz/max(n_iz,1), 1));
[d_obs, sd, idx, is_ph] = H_pack(obs, use, nlook, Nz, Np, clip_db);
sd = sd * redund;
Cdi = 1 ./ (sd.^2);                           % diagonal inverse covariance

fwd = @(mm) H_pack_pred(H_forward(mm, ip, z_layer, z, psi, fwd_opts), use, idx, clip_db);
% RESIDUAL, not a plain difference. Two things make a subtraction wrong
% here and both were measured, not anticipated:
%   phi is an ANGLE. An unwrapped difference near +-pi returns ~2pi where
%   the true misfit is ~0, and since those entries carry the tightest
%   sigma where coherence is high they then dominate everything. Wrapping
%   took chi2/dof from ~2200 to order 1.
%   dP_hv has TRUE ZEROS on the fabric axes, so its dB anomaly is -Inf
%   there and the residual is set by how close two near-nulls land rather
%   than by whether the model is right. Both data and model are floored at
%   clip_db, which says we do not claim to know how deep a null is once it
%   is below anything measurable.
resid = @(mm) H_resid(d_obs, fwd(mm), is_ph);
loss = @(mm) sum(Cdi .* resid(mm).^2) + (mm-m_p).' * Cmi * (mm-m_p);

% --- MULTI-START OVER theta AND dlam. Both directions are multi-modal and
% for different reasons: theta because the response is periodic in it and
% the quarter-turn swap is a near-alias, dlam because it enters only
% through the ACCUMULATED phase, so the misfit oscillates in it. Measured
% in test_polarimetric_ratio, the converged loss varies four-fold across
% starting dlam alone. Starting from the prior mean in both - which is
% what a naive GLS does - reliably lands in a local minimum: on noise-free
% synthetic data that the model reproduces exactly, it returned a REVERSED
% dlam profile at chi2/dof of 2900.
%
% The grid is scored before any iteration and only the best few are run to
% convergence, so the cost is a handful of forward evaluations plus
% n_start solves.
th_try = H_opt(opts, 'theta_scan', []);
if isempty(th_try)
  th_try = linspace(0, pi, H_opt(opts, 'n_theta_scan', 6) + 1); th_try(end) = [];
end
dl_try = H_opt(opts, 'dlam_scan', [0.01 0.02 0.035 0.05 0.07 0.10]);
cand = []; cl = [];
for th = th_try
  for d0 = dl_try
    mt = m_p; mt(ip.th) = th; mt(ip.dl) = d0; mt(ip.rd) = rd_pr(1);
    cand(end+1, :) = [th d0]; %#ok<AGROW>
    cl(end+1) = loss(mt); %#ok<AGROW>
  end
end
[~, ord] = sort(cl);
best = struct('L', inf);
for q = 1:min(n_start, numel(ord))
  m = m_p; m(ip.th) = cand(ord(q),1); m(ip.dl) = cand(ord(q),2); m(ip.rd) = rd_pr(1);
  [m, L, conv, it] = H_solve(m, fwd, d_obs, Cdi, Cmi, m_p, loss, max_iter, ip, hstep, is_ph);
  if L(end) < best.L
    best = struct('L', L(end), 'm', m, 'hist', L, 'conv', conv, 'it', it);
  end
end
m = best.m;

% --- posterior at the solution
% HOW FAR TO TRUST THIS. test_fabric_gls compares the reported sigmas with
% the scatter of 16 independent noisy inversions and MEASURES both theta
% and dlam scattering about twice as far as this posterior claims. That is
% the ordinary cost of linearising a non-linear, multi-modal problem.
% sigma_theta_z and sigma_dlam_z are therefore LOWER BOUNDS: useful for
% ranking depths and comparing frames, not for a quoted confidence
% interval. For an interval, resample - ptt.quadpolJackknife.
G = H_jac(fwd, m, numel(d_obs), hstep);
A = G.' * (Cdi .* G) + Cmi;
C_M = H_inv(A);
Rres = C_M * (G.' * (Cdi .* G));
res = real(diag(Rres));
sig = sqrt(max(real(diag(C_M)), 0));

r_fin = H_resid(d_obs, fwd(m), is_ph);
dof = max(numel(d_obs) - sum(res), 1);
chi2_dof = sum(Cdi .* r_fin.^2) / dof;

% CHI-SQUARE SCALING. If the fit is worse than C_d claims it should be,
% something the covariance does not describe is present - unmodelled
% structure, leakage, residual correlation - and the honest response is to
% widen the error bars by that factor rather than to keep quoting a
% posterior the residuals contradict. Only ever widens: a fit better than
% assumed does not license tighter bounds, because on a non-linear problem
% that usually means a local minimum, not a good day. sigma_m_raw keeps
% the unscaled numbers so the two can be compared.
sig_raw = sig;
sig = sig * sqrt(max(chi2_dof, 1));

M = H_forward(m, ip, z_layer, z, psi, fwd_opts);
th_z = interp1(z_mid, m(ip.th), z, 'linear', 'extrap');
dl_z = interp1(z_mid, m(ip.dl), z, 'linear', 'extrap');
rd_z = interp1(z_mid, m(ip.rd), z, 'linear', 'extrap');
out = struct('theta_z', th_z, 'dlam_z', dl_z, 'rdb_z', rd_z, ...
  'sigma_theta_z', interp1(z_mid, sig(ip.th), z, 'linear', 'extrap'), ...
  'sigma_dlam_z',  interp1(z_mid, sig(ip.dl), z, 'linear', 'extrap'), ...
  'sigma_rdb_z',   interp1(z_mid, sig(ip.rd), z, 'linear', 'extrap'), ...
  'C_M', C_M, 'resolution', res, 'chi2_dof', chi2_dof, ...
  'loss', best.hist, 'converged', best.conv, 'n_iter', best.it, ...
  'z_layer', z_layer, 'z_mid', z_mid, 'pred', M, 'm', m, 'ip', ip, ...
  'used', {use}, 'n_looks', nlook, 'sigma_m', sig, ...
  'sigma_m_raw', sig_raw, 'redundancy', redund);
end

% -------------------------------------------------------------------------
function [m, L, converged, it] = H_solve(m, fwd, d_obs, Cdi, Cmi, m_p, loss, max_iter, ip, hstep, is_ph)
L = nan(max_iter+1,1); L(1) = loss(m); converged = false; it = 0;
for it = 1:max_iter
  G = H_jac(fwd, m, numel(d_obs), hstep);
  A = G.'*(Cdi .* G) + Cmi;
  b = G.'*(Cdi .* H_resid(d_obs, fwd(m), is_ph)) - Cmi*(m - m_p);
  dm = A \ b;
  if ~all(isfinite(dm)), break; end
  a = 1; ok = false;
  for h = 1:20
    mt = H_clamp(m + a*dm, ip);
    Lt = loss(mt);
    if Lt < L(it), m = mt; L(it+1) = Lt; ok = true; break; end
    a = a/2;
  end
  if ~ok, converged = true; L(it+1) = L(it); break; end
  if abs(L(it)-L(it+1))/max(abs(L(it+1)),realmin) < 1e-5, converged = true; break; end
end
L = L(1:min(it+1, numel(L)));
end

function M = H_forward(m, ip, z_layer, z, psi, fwd_opts)
nL = numel(z_layer);
lay = struct('top_m', num2cell(z_layer), ...
  'dlam', num2cell(max(m(ip.dl), 0)), ...
  'theta', num2cell(m(ip.th)), ...
  'r_db', num2cell(m(ip.rd)), ...
  'gx_db', num2cell(zeros(nL,1)));
% NOTE the sign: ptt.fujitaModel follows the paper's R S R' sense, where a
% fabric at physical azimuth theta appears at sweep index -theta. theta is
% carried here in the SAME sense as that model, so out.theta_z is directly
% comparable with ptt.fujitaModel and with ptt.ershadiFabric, and is the
% NEGATIVE of the synthesis-sense angle used by ptt.quadpolFabricLS.
M = ptt.fujitaModel(lay, z, psi, fwd_opts);
end

function [d, sd, idx, is_ph] = H_pack(obs, use, nlook, Nz, Np, clip_db)
d = []; sd = []; idx = struct(); is_ph = [];
% dB variance of a log-power from N independent looks: the log of a
% Gamma(N) variable has variance psi'(N), and 10/ln10 converts to dB.
s_db = (10/log(10)) * sqrt(psi(1, nlook));
for k = 1:numel(use)
  f = use{k};
  O = obs.(f);
  sname = ['sigma_' f];
  if isfield(obs, sname) && ~isempty(obs.(sname))
    S = obs.(sname); if isscalar(S), S = S*ones(size(O)); end
  elseif strcmp(f, 'phi')
    % CRB of a coherence phase: (1-|C|^2)/(2 N |C|^2). Needs |C|; without
    % it we cannot pretend to know the phase precision, so it is an error
    % rather than a guess - a wrong phase weight silently rotates the axis.
    if ~isfield(obs, 'Cmag') || isempty(obs.Cmag)
      error('ptt:fabricGLS:phiSigma', ...
        ['phi is used but neither obs.sigma_phi nor obs.Cmag was given; ' ...
         'the phase variance cannot be assumed']);
    end
    c = min(max(obs.Cmag, 1e-3), 0.999);
    S = sqrt((1 - c.^2) ./ (2 * nlook * c.^2));
  else
    S = s_db * ones(size(O));
  end
  if ~strcmp(f, 'phi'), O = max(O, clip_db); end
  O = O(:); S = S(:);
  g = isfinite(O) & isfinite(S) & S > 0;
  idx.(f) = g;
  d = [d; O(g)]; %#ok<AGROW>
  sd = [sd; S(g)]; %#ok<AGROW>
  is_ph = [is_ph; strcmp(f, 'phi') * ones(nnz(g), 1)]; %#ok<AGROW>
end
is_ph = logical(is_ph);
end

function v = H_pack_pred(M, use, idx, clip_db)
v = [];
for k = 1:numel(use)
  f = use{k};
  A = M.(f);
  if ~strcmp(f, 'phi'), A = max(A, clip_db); end
  A = A(:);
  v = [v; A(idx.(f))]; %#ok<AGROW>
end
end

function r = H_resid(d, g, is_ph)
r = d - g;
r(is_ph) = angle(exp(1i * r(is_ph)));   % wrap the phase entries
end

function C = H_gp(zc, sigma, L)
% Exponential-correlation prior: C_ij = sigma^2 exp(-|zi-zj|/L). This is
% what makes "the axis is nearly constant with depth" a PRIOR with a
% length scale rather than a hard constraint - lengthen L to insist on it,
% shorten L to let the data rotate the axis.
D = abs(zc - zc.');
C = sigma^2 * exp(-D / max(L, eps));
C = C + 1e-10 * sigma^2 * eye(numel(zc));     % keep it invertible
end

function Ai = H_inv(A)
[R, p] = chol(A);
if p == 0
  Ai = R \ (R.' \ eye(size(A)));
else
  Ai = pinv(A);
end
Ai = (Ai + Ai.')/2;
end

function G = H_jac(fwd, m, nd, hstep)
%H_JAC Finite-difference Jacobian, stepped on the PRIOR SCALE.
% Stepping by a fraction of each parameter's own VALUE is wrong when the
% parameters have different units and some are legitimately zero: r_db
% starts at 0 dB, so a value-scaled step collapsed to the 1e-6 floor and
% that whole block of the Jacobian was rounding noise, which poisoned
% every Gauss-Newton step. The prior sigma is the natural scale of each
% parameter and is always finite, so the step is a fixed small fraction of
% it: theta ~1.6e-3 rad, dlam ~8e-5, r_db ~4e-3 dB at the defaults.
np = numel(m); G = zeros(nd, np); f0 = fwd(m);
for j = 1:np
  h = hstep(j);
  mp = m; mp(j) = mp(j) + h;
  G(:, j) = (fwd(mp) - f0) / h;
end
G(~isfinite(G)) = 0;
end

function m = H_clamp(m, ip)
m(ip.dl) = max(m(ip.dl), 0);
m(ip.rd) = min(max(m(ip.rd), -30), 30);
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

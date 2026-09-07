function out = traveltimeFabricML(dtau, z, opts)
%TRAVELTIMEFABRICML Horizontal fabric contrast from travel-time anomalies.
%
% out = ptt.traveltimeFabricML(dtau, z, opts)
%
% Step 2 of the Nymand two-step inversion: solves the depth profile of the
% horizontal fabric contrast dlam(z) = lam_x - lam_y from the two-way
% travel-time difference between the two eigen-polarizations, by LINEAR
% maximum likelihood.
%
% WHY THIS IS LINEAR, AND WHY THAT MATTERS. The birefringent delay is the
% running integral of the contrast,
%
%     dtau(z) = k * INT_0^z dlam(z') dz',   k = deps / (n_ice * c),
%
% so dlam enters the data through a lower-triangular integration matrix
% and nothing else. The problem is therefore exactly linear, and the
% Gaussian maximum-likelihood solution is the weighted least-squares one
% with a posterior that is EXACT rather than linearised. That is the
% substantive difference from ptt.fabricGLS, whose posterior is a
% linearisation of a strongly non-linear fit and which test_fabric_gls
% MEASURES to be about 2x optimistic over 16 Monte Carlo realisations.
% Here the covariance below is the covariance, given the noise model.
%
% THE CATCH, stated plainly. Integration smooths, so inversion
% differentiates, and differentiating noise is how a travel-time profile
% turns into a spiky dlam that looks like layering. This is an
% ill-conditioned problem and it is regularised by a prior, not by luck:
% C_m is the same exponential-correlation Gaussian process ptt.fabricGLS
% uses, so "dlam varies smoothly over L metres" is a stated assumption
% with a length scale rather than an implicit one. Depths where the data
% did not determine the answer are reported through the RESOLUTION
% diagonal and abstained on when it falls below res_min - the estimator
% declines rather than handing back the prior dressed as a measurement.
%
% dtau is the travel-time difference in NANOSECONDS at each z, following
% this repo's metre/nanosecond convention (ptt.constants). It must be
% measured between the two eigen-polarizations, i.e. after synthesis to
% the theta0 frame that step 1 returns - a dtau measured in the antenna
% frame answers a different question and will bias dlam low by cos(2*dth).
% NaNs are dropped; the solve needs at least three finite samples.
%
% opts (all optional):
%   .fc        (750e6)  centre frequency [Hz], for the phase convention
%   .prior     [mean sigma corr_len_m], default [0.03 0.05 200]
%   .sigma_ns  (0.05)   1-sigma travel-time noise [ns], scalar or per-z
%   .res_min   (0.15)   resolution below which a depth abstains (NaN)
%   .nonneg    (false)  clip the posterior mean at zero. OFF by default:
%                       dlam is a signed difference once its axes are
%                       named, and clipping hides a sign error rather
%                       than fixing it.
%
% out.dlam        posterior mean contrast at each z (NaN where abstained)
% out.sigma_dlam  posterior standard deviation, EXACT for this model
% out.resolution  diagonal of the resolution matrix, 1 = data, 0 = prior
% out.chi2_dof    reduced chi-square of the travel-time fit
% out.pred        predicted dtau at the solution [ns]
% out.C_M         full posterior covariance
% out.k_ns_per_m  the forward constant actually used
%
% See also ptt.nymandTwoStep, ptt.fabricGLS, ptt.invertHorizontalFabric.
if nargin < 3, opts = struct(); end
z = z(:); dtau = dtau(:);
if numel(dtau) ~= numel(z)
  error('ptt:traveltimeFabricML:size', 'dtau and z must be the same length');
end
fc      = H_opt(opts, 'fc', 750e6);
pr      = H_opt(opts, 'prior', [0.03, 0.05, 200]);
sig_ns  = H_opt(opts, 'sigma_ns', 0.05);
res_min = H_opt(opts, 'res_min', 0.15);
nonneg  = H_opt(opts, 'nonneg', false);

C = ptt.constants();
n_ice = sqrt(C.eps_bar);
% Same phase convention as ptt.quadpolFabric's grad_per_dlam, converted
% from radians per metre to nanoseconds per metre so the two estimators
% cannot drift apart: grad_per_dlam = 2*pi*fc*2*(deps/(2*n))/(c*1e9).
grad_per_dlam = 2*pi*fc * 2 * (C.deps / (2*n_ice)) / (C.c * 1e9);   % rad/m
k = grad_per_dlam / (2*pi*fc) * 1e9;                                % ns/m

g = isfinite(dtau) & isfinite(z);
if nnz(g) < 3
  out = H_empty(z, k);
  return
end
zg = z(g); dg = dtau(g);
[zg, ord] = sort(zg); dg = dg(ord);

% --- forward operator: dtau_i = k * sum_j w_ij * dlam_j
% Trapezoidal weights on the observed depth grid, so the integral does not
% depend on how evenly the samples happen to be spaced.
nz = numel(zg);
dz = diff(zg);
W = zeros(nz, nz);
for i = 2:nz
  W(i, :) = W(i-1, :);
  W(i, i-1) = W(i, i-1) + dz(i-1)/2;
  W(i, i)   = W(i, i)   + dz(i-1)/2;
end
G = k * W;

% --- data covariance
if isscalar(sig_ns), sd = sig_ns * ones(nz, 1); else, sd = sig_ns(g); sd = sd(ord); sd = sd(:); end
sd = max(sd, eps);
Cdi = 1 ./ (sd.^2);

% --- prior
m_p = pr(1) * ones(nz, 1);
Cm  = H_gp(zg, pr(2), pr(3));
Cmi = H_inv(Cm);

% --- exact linear ML / Bayesian least squares
A   = G.' * (Cdi .* G) + Cmi;
C_M = H_inv(A);
m   = C_M * (G.' * (Cdi .* dg) + Cmi * m_p);
Rres = C_M * (G.' * (Cdi .* G));
res  = real(diag(Rres));
sig  = sqrt(max(real(diag(C_M)), 0));

pred = G * m;
r    = dg - pred;
dof  = max(nz - sum(res), 1);
chi2 = sum(Cdi .* r.^2) / dof;

if nonneg, m = max(m, 0); end
% Abstain where the data did not determine the answer. Reporting the
% prior mean there would look like a measurement of 0.03.
m(res < res_min) = NaN;

out = struct();
out.z = zg;
out.dlam = m;
out.sigma_dlam = sig;
out.resolution = res;
out.chi2_dof = chi2;
out.pred = pred;
out.resid = r;
out.C_M = C_M;
out.k_ns_per_m = k;
out.n_used = nz;
out.n_abstained = nnz(res < res_min);
end

% -------------------------------------------------------------------------
function out = H_empty(z, k)
out = struct('z', z, 'dlam', nan(numel(z),1), 'sigma_dlam', nan(numel(z),1), ...
  'resolution', zeros(numel(z),1), 'chi2_dof', NaN, 'pred', nan(numel(z),1), ...
  'resid', nan(numel(z),1), 'C_M', [], 'k_ns_per_m', k, 'n_used', 0, ...
  'n_abstained', numel(z));
end

function C = H_gp(zc, sigma, L)
% Exponential-correlation prior, identical in form to ptt.fabricGLS so the
% two estimators regularise the same way.
D = abs(zc - zc.');
C = sigma^2 * exp(-D / max(L, eps));
C = C + 1e-10 * sigma^2 * eye(numel(zc));
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

function v = H_opt(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

%TEST_NYMAND_TWO_STEP The two-step fabric inversion, and its exact posterior.
%
% ptt.nymandTwoStep implements the scheme presented for NEGIS: the Fujita
% depolarization-matrix model as the forward problem, simplified by a
% CONSTANT EIGENFRAME with depth, inverted in two steps - orientation and
% scattering ratio from the power anomalies, then the horizontal contrast
% from the travel-time anomalies.
%
% The verdicts, in the order that matters:
%
%   1. THE FORWARD CONSTANT IS THE ONE THE REST OF THE REPO USES. The
%      delay per metre per unit dlam is derived from ptt.constants by the
%      same expression as ptt.quadpolFabric's grad_per_dlam, only
%      converted to ns/m. If these two ever disagree, one estimator's
%      dlam is on a different scale from the other's and every
%      cross-comparison in the project is silently wrong.
%   2. STEP 2 RECOVERS dlam(z). Travel times generated from a known
%      profile invert back to it.
%   3. THE POSTERIOR IS EXACT - THE HEADLINE CLAIM. Because the delay is
%      the running integral of dlam, the problem is exactly linear, so
%      the Gaussian ML covariance is the covariance and not a
%      linearisation. Checked the same way ptt.fabricGLS's was and held
%      to a much tighter band: fabricGLS MEASURED ~2x optimistic, this
%      one has no linearisation error to be optimistic about. A ratio far
%      from 1 here means the noise model or the operator is wrong, not
%      that non-linearity bit.
%   4. IT ABSTAINS RATHER THAN RETURNING THE PRIOR. Integration smooths,
%      so inversion differentiates and the shallowest depths are poorly
%      resolved. Those must come back NaN, not as the prior mean wearing
%      a measurement's clothes.
%   5. END TO END, AND THE COUPLING IS REAL. The two steps are chained on
%      data generated from one known (theta0, dlam(z)). The verdict also
%      pins WHY this iterates: a dtau measured in the antenna frame
%      rather than the eigenframe is biased by cos(2*dtheta), so a run
%      given the wrong frame must do measurably worse.
%
% Run: matlab -batch "run('opr_fabric/test/test_nymand_two_step.m')"
clear; t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

fc = 750e6; z = (10:10:1400).';
dl_true = 0.03 + 0.04 * (z / 1400);
TH = deg2rad(35);
C = ptt.constants();

%% ---- verdict 1: the forward constant matches the rest of the repo
n_ice = sqrt(C.eps_bar);
gpd = 2*pi*fc * 2 * (C.deps / (2*n_ice)) / (C.c * 1e9);   % ptt.quadpolFabric
k_expect = gpd / (2*pi*fc) * 1e9;                          % ns per m per dlam
o0 = ptt.traveltimeFabricML(zeros(size(z)), z, struct('fc', fc));
fprintf('1. forward constant: estimator %.6f ns/m, quadpolFabric %.6f ns/m\n', ...
  o0.k_ns_per_m, k_expect);
ok1 = abs(o0.k_ns_per_m - k_expect) < 1e-12;
fprintf('   one delay scale across the project:                   %s\n', H_tick(ok1));

%% ---- verdict 2: step 2 recovers dlam(z)
dtau = k_expect * cumtrapz(z, dl_true);
o2 = ptt.traveltimeFabricML(dtau, z, struct('fc', fc));
e = abs(o2.dlam - dl_true);
fprintf('\n2. noise-free step 2: median |err| %.5f, max %.5f, dtau(end) %.2f ns\n', ...
  median(e, 'omitnan'), max(e, [], 'omitnan'), dtau(end));
ok2 = median(e, 'omitnan') < 0.002 && max(e, [], 'omitnan') < 0.01;
fprintf('   recovers the contrast profile:                        %s\n', H_tick(ok2));

%% ---- verdict 3: the posterior is EXACT, checked by Monte Carlo
S_NS = 0.05;
% CALIBRATION IS A STATEMENT ABOUT TRUTHS DRAWN FROM THE PRIOR, so that is
% how it is tested here. A posterior says: if the world is distributed the
% way my prior claims, my error bars cover at the stated rate. Testing it
% against ONE fixed smooth profile sitting near the prior mean measures
% something else - how lucky that profile is - and MEASURED, it is lucky:
% the dlam ramp used above comes back with a reported sigma about 2x its
% actual error, because a smooth profile near the prior mean is easier
% than the prior admits. That is the posterior being conservative for a
% convenient truth, not the posterior being wrong.
%
% So each realisation draws its own truth from the prior covariance,
% generates travel times through the SAME operator the estimator inverts,
% and the normalised error (m_hat - m_true)/sigma is pooled. If the linear
% ML solution is the poster's
%     m~  = (G' Cd^-1 G + Cm^-1)^-1 (G' Cd^-1 d + Cm^-1 m_prior)
%     Cm~ = (G' Cd^-1 G + Cm^-1)^-1
% and the operator G is right, that pooled quantity is standard normal and
% its spread is 1. Anything else means the covariance, the prior or the
% integration operator is wrong - which is exactly what this verdict is
% for.
PR = [0.03, 0.05, 200];
D = abs(z - z.');
Cm_t = PR(2)^2 * exp(-D / PR(3)) + 1e-10 * PR(2)^2 * eye(numel(z));
Lc = chol(Cm_t + 1e-12*eye(numel(z)), 'lower');
R = 24;
zsc = [];
for kk = 1:R
  rng(300+kk);
  m_true = PR(1) + Lc * randn(numel(z), 1);
  d_k = k_expect * cumtrapz(z, m_true) + S_NS * randn(numel(z), 1);
  ok = ptt.traveltimeFabricML(d_k, z, ...
    struct('fc', fc, 'sigma_ns', S_NS, 'prior', PR, 'res_min', 0));
  gk = isfinite(ok.dlam) & isfinite(ok.sigma_dlam) & ok.sigma_dlam > 0;
  zsc = [zsc; (ok.dlam(gk) - m_true(gk)) ./ ok.sigma_dlam(gk)]; %#ok<AGROW>
end
spread = std(zsc);
cov68 = mean(abs(zsc) <= 1);
fprintf('\n3. Bayesian calibration, %d realisations with truths drawn from the prior:\n', R);
fprintf('   normalised error: std %.3f (want 1.00), |z|<=1 covers %.3f (want 0.68)\n', ...
  spread, cov68);
ok3 = spread > 0.8 && spread < 1.25 && cov68 > 0.60 && cov68 < 0.76;
fprintf('   the posterior is exact for this linear model:          %s\n', H_tick(ok3));

%% ---- verdict 4: abstains where the data did not determine the answer
rng(7);
o4 = ptt.traveltimeFabricML(dtau + S_NS*randn(size(dtau)), z, ...
  struct('fc', fc, 'sigma_ns', S_NS, 'res_min', 0.15));
nab = nnz(isnan(o4.dlam));
pri = 0.03;
near_prior = nnz(abs(o4.dlam - pri) < 1e-6);
fprintf('\n4. abstained at %d of %d depths; cells sitting exactly at the prior: %d\n', ...
  nab, numel(z), near_prior);
ok4 = nab >= 1 && near_prior == 0 && all(o4.resolution(isfinite(o4.dlam)) >= 0.15);
fprintf('   declines rather than reporting the prior:              %s\n', H_tick(ok4));

%% ---- verdict 5: end to end, and the eigenframe matters
nL = 8; top = (0:nL-1).' * (1400/nL);
dl_lay = interp1(z, dl_true, top + (1400/nL)/2, 'linear', 'extrap');
lay = struct('top_m', num2cell(top), 'dlam', num2cell(dl_lay(:)), ...
  'theta', num2cell(TH*ones(nL,1)), 'r_db', num2cell(zeros(nL,1)), ...
  'gx_db', num2cell(zeros(nL,1)));
psi = (0:5:175) * pi/180;
fm = ptt.fujitaModel(lay, z, psi, struct('fc', fc, 'win_m', 30));
obs = struct('psi', psi, 'dP_hh', fm.dP_hh, 'dP_hv', fm.dP_hv, 'dtau', dtau);
base = struct('fc', fc, 'n_r', 6, 'eta', 1e-2, 'azimuth_source', 'physical', ...
  'max_outer', 4, 'theta0', deg2rad(20));
oe = ptt.nymandTwoStep(obs, z, base);
% SCORED AGAINST -TH, NOT TH. ptt.fujitaModel deliberately follows the
% paper's R S R' sense, in which a fabric at azimuth theta puts its
% features at sweep index -theta, while ptt.polarimetricInverse uses the
% synthesis sense. test_polarimetric_inverse pins that difference and
% MEASURES it (correlation 0.967 matched against 0.11 mismatched); this
% verdict inherits it because the data here is generated by fujitaModel
% and inverted by polarimetricInverse. Scoring against +TH reported an
% 18.4 deg error that was entirely the convention.
TH_EXPECT = -TH;
% theta is an AXIS: scored circularly on the doubled angle, mod 90 because
% the co-pol power anomaly cannot separate (theta0, r) from (theta0+90, 1/r)
eth = rad2deg(abs(angle(exp(4i*(oe.theta0 - TH_EXPECT)))/4));
edl = median(abs(oe.dlam - dl_true), 'omitnan');
fprintf('\n5. end to end: theta0 err %.2f deg (mod 90), median dlam err %.5f, %d outer iters\n', ...
  eth, edl, oe.n_outer);

% a delay measured in the ANTENNA frame instead of the eigenframe
dth = deg2rad(30);
oe_bad = ptt.nymandTwoStep(setfield(obs, 'dtau', dtau * cos(2*dth)), z, base); %#ok<SFLD>
edl_bad = median(abs(oe_bad.dlam - dl_true), 'omitnan');
fprintf('   with dtau in the wrong frame: median dlam err %.5f (%.1fx worse)\n', ...
  edl_bad, edl_bad / max(edl, eps));
ok5 = eth < 6 && edl < 0.006 && edl_bad > 2*edl;
fprintf('   chains correctly, and the eigenframe demonstrably matters: %s\n', H_tick(ok5));

fails = ~ok1 + ~ok2 + ~ok3 + ~ok4 + ~ok5;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_nymand_two_step:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

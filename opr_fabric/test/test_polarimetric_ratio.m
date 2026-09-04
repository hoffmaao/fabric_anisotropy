%TEST_POLARIMETRIC_RATIO Orientation and strength from channel power ratios.
%
% ptt.polarimetricRatioInverse solves the Fujita model, at a depth-constant
% orientation, for theta0, dlam(z) and the scattering ratio r(z) from the
% RATIOS between measured channels, with each ratio carrying a free
% additive offset so that unknown H-versus-V channel gains drop out.
%
% The synthetic injects gains the solver is never told: |t|^2 = +6 dB on
% the transmit-V chain and |r|^2 = -4 dB on receive-V, which is the scale
% of the real imbalance on this system (ptt.equaliseChannels).
%
%   1. RECOVERY. With four headings and power + phase, theta0 within 1 deg,
%      dlam within 0.005, r within 0.05, pedestal within 0.05.
%   2. THE 90-DEGREE ALIAS IS REAL. Power ratios alone cannot separate
%      (theta0, r) from (theta0 + 90, 1/r): the swap is absorbed by the
%      free offsets. The test asserts the alias DOES capture a power-only
%      fit, so nobody deletes the phase term believing it optional, and
%      that out.alias_unresolved is flagged when phi is absent.
%   3. THE PHASE BREAKS IT. Adding phi recovers the true axis.
%   4. NOISE. 0.5 dB on every power and 0.05 rad on the phase leaves the
%      answers inside the case-1 tolerances at 3x.
%   5. ONE HEADING IS NOT ENOUGH FOR THE AXIS. A counting argument says it
%      should be; it is not, because with a single d the fit trades
%      orientation against r(z), which is free at every depth. The test
%      asserts BOTH halves of the measured behaviour - strength still
%      comes out, orientation does not - so that nobody quotes an axis
%      from a single-heading frame on the strength of the counting
%      argument.
%
% Run: matlab -batch "run('opr_fabric/test/test_polarimetric_ratio.m')"
clear; t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

fc = 750e6; C = ptt.constants();
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar)*C.c*1e9);
z = (10:10:1500).'; dz = [0; diff(z)];
TH = deg2rad(35);
DL = 0.03 + 0.04*(z/z(end));
RR = 1 + 0.4*(z/z(end));
PED = 0.30;
GT = 10^(6/10);    % transmit V gain, +6 dB, unknown to the solver
GR = 10^(-4/10);   % receive V gain, -4 dB
delta = cumsum(DL.*dz)*gpd;
OPT = struct('fc', fc, 'n_r', 6, 'n_d', 10, 'eta_r', 3e-2, 'eta_d', 3e-2, ...
  'max_iter', 80, 'n_start', 6);

o1 = H_fit(H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0, 0), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv','phi'});
[e1, d1, r1] = H_err(o1, TH, DL, RR);
fprintf('1. four headings, power + phase: theta0 err %.2f deg, dlam err %.4f, r err %.3f, ped %.2f\n', ...
  e1, d1, r1, o1.ped);
ok1 = e1 < 1 && d1 < 0.005 && r1 < 0.05 && abs(o1.ped - PED) < 0.05;
fprintf('   recovers orientation AND strength through unknown gains: %s\n', H_tick(ok1));

o2 = H_fit(H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0, 0), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv'});
[e2, ~, ~] = H_err(o2, TH, DL, RR);
fprintf('\n2. power only: theta0 err %.1f deg, alias_unresolved %d\n', e2, o2.alias_unresolved);
ok2 = e2 > 45 && o2.alias_unresolved;
fprintf('   the 90-deg alias is REAL and is flagged:              %s\n', H_tick(ok2));
fprintf('3. adding the coherence phase resolves it:               %s\n', H_tick(e1 < 1));
ok3 = e1 < 1;

rng(7);
o4 = H_fit(H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0.5, 0.05), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv','phi'});
[e4, d4, r4] = H_err(o4, TH, DL, RR);
fprintf('\n4. with noise: theta0 err %.2f deg, dlam err %.4f, r err %.3f\n', e4, d4, r4);
ok4 = e4 < 3 && d4 < 0.015 && r4 < 0.15;
fprintf('   degrades gracefully:                                  %s\n', H_tick(ok4));

o5 = H_fit(H_obs(deg2rad(15), TH, RR, delta, PED, GT, GR, 0, 0), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv','phi'});
[e5, d5, ~] = H_err(o5, TH, DL, RR);
fprintf('\n5. single heading: theta0 err %.2f deg, dlam err %.4f (four headings gave %.2f, %.4f)\n', ...
  e5, d5, e1, d1);
ok5 = d5 < 0.03 && e5 > 5 && e5 > 10*e1;
fprintf('   strength survives one heading, orientation does NOT:   %s\n', H_tick(ok5));

fails = ~ok1 + ~ok2 + ~ok3 + ~ok4 + ~ok5;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_polarimetric_ratio:failed', '%d verdict(s) failed', fails);
end

function obs = H_obs(gam, th, r, delta, ped, gt, gr, noise_db, noise_rad)
d = th - gam; c2 = cos(d).^2; s2 = sin(d).^2; cd = cos(delta);
Phh = c2.^2 + (r.^2).*s2.^2 + 2*r.*c2.*s2.*cd;
Pvv = (r.^2).*c2.^2 + s2.^2 + 2*r.*c2.*s2.*cd;
Pxx = (1 + r.^2 - 2*r.*cd).*(c2.*s2) + ped*(Phh+Pvv)/2;
t2 = s2 ./ c2; t4 = t2.^2;
phi = atan2(r.*sin(delta).*(1-t4), r.*cos(delta).*(1+t4) + t2.*(1+r.^2));
n = @(A) A .* 10.^(noise_db*randn(size(A))/10);
obs = struct('gamma', gam, 'P_hh', n(Phh), 'P_vv', n(Pvv*gr*gt), ...
  'P_hv', n(Pxx*gt), 'P_vh', n(Pxx*gr), 'phi', phi + noise_rad*randn(size(phi)));
end

function o = H_fit(obs, z, OPT, use)
OPT.use = use;
o = ptt.polarimetricRatioInverse(obs, z, OPT);
end

function [e, d, r] = H_err(o, TH, DL, RR)
e = rad2deg(abs(angle(exp(2i*(o.theta0 - TH)))/2));
d = median(abs(o.dlam_z - DL));
r = median(abs(o.r_z - RR));
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

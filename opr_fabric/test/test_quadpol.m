%TEST_QUADPOL Round-trip test of the quad-pol scattering-matrix framework.
%
% Builds a synthetic quad-pol column from a KNOWN horizontal fabric - a
% fixed principal azimuth and a contrast that grows with depth - propagates
% it through the analytic birefringence model, and checks that
% ptt.quadpolMoments -> quadpolAzimuth -> quadpolFabric recovers both the
% orientation and the contrast.
%
% The point of the test is the two properties the co-polarized pair does
% not have, so both are asserted explicitly:
%   1. orientation is recovered from a SINGLE antenna azimuth, at every
%      azimuth the synthetic column is measured from;
%   2. the recovered contrast is the TRUE lam_max - lam_min, not the
%      projected one, so it does not vary with the antenna azimuth.
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol.m')"
clear;
rng(7);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
projRoot = fullfile(thisDir, '..', '..');
addpath(projRoot);

C = ptt.constants();
fc = 750e6;
n_ice = sqrt(C.eps_bar);
grad_per_dlam = 2*pi*fc * 2 * (C.deps/(2*n_ice)) / (C.c*1e9);   % rad/m

% --- truth
THETA_TRUE = deg2rad(35);      % principal axis, 35 deg
DLAM_TRUE = 0.30;              % uniform contrast, so d(delta)/dz is constant
Nz = 3000;
z = (0:Nz-1).' * 0.5;          % 0.5 m sampling to 1500 m
delta = grad_per_dlam * DLAM_TRUE * z;   % differential phase vs depth
Nx = 200;                      % traces to average over

fprintf('truth: theta %.1f deg, dlam %.3f, %.2f fringes over %.0f m\n', ...
  rad2deg(THETA_TRUE), DLAM_TRUE, delta(end)/(2*pi), z(end));

fails = 0;
for alpha_deg = [0 20 35 55 80 125]
  % --- forward: a birefringent column measured with antennas at alpha.
  % In the PRINCIPAL frame the scattering matrix of the column is
  % diag(exp(i delta), 1) times a common (irrelevant) layer response; the
  % antenna frame is that rotated by -(theta - alpha).
  d = THETA_TRUE - deg2rad(alpha_deg);
  cd_ = cos(d); sd = sin(d);
  ex = exp(1i * delta);
  ey = ones(Nz, 1);
  % S_ant = R(d) diag(ex, ey) R(d)'
  hh = cd_^2 * ex + sd^2 * ey;
  vv = sd^2 * ex + cd_^2 * ey;
  hv = cd_*sd * (ex - ey);
  vh = hv;                       % reciprocal medium

  % Random layer reflectivity per trace, common to all four channels -
  % this is what makes the channels coherent with each other while each
  % on its own looks like noise, exactly as a real englacial return does.
  r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  nz_amp = 0.02;                 % additive receiver noise
  S = struct();
  S.hh = hh .* r + nz_amp*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vv = vv .* r + nz_amp*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.hv = hv .* r + nz_amp*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vh = vh .* r + nz_amp*(randn(Nz,Nx)+1i*randn(Nz,Nx));

  % --- inverse
  M = ptt.quadpolMoments(S, [9 Nx]);
  psi = (0:2:178) * pi/180;
  A = ptt.quadpolAzimuth(M, psi);
  out = ptt.quadpolFabric(A, z, struct('fc', fc, 'win_m', 50, ...
    'grad_win_m', 200));

  % theta is recovered in the ANTENNA frame, so the truth for this pass is
  % the fabric axis relative to the antennas, modulo 90 deg.
  want = mod(rad2deg(d), 90);
  got = mod(rad2deg(median(out.theta, 'omitnan')), 90);
  err = abs(mod(got - want + 45, 90) - 45);

  mid = z > 300 & z < 1200;
  dl = median(out.dlam(mid), 'omitnan');
  dln = median(out.dlam_node(mid), 'omitnan');

  ok_t = err < 3;
  ok_d = abs(dl - DLAM_TRUE) < 0.03;
  ok_n = ~isfinite(dln) || abs(dln - DLAM_TRUE) < 0.06;
  fails = fails + ~(ok_t && ok_d && ok_n);
  fprintf(['alpha %3d: theta %5.1f (want %5.1f, err %4.1f) %s | ' ...
           'dlam %.3f %s | node %.3f %s\n'], alpha_deg, got, want, err, ...
    H_tick(ok_t), dl, H_tick(ok_d), dln, H_tick(ok_n));
end

% --- the projection property: the co-polarized phase alone would have
% given dlam*cos(2(alpha-theta)), which vanishes at 45 deg to the axes.
% Assert the quad-pol answer does NOT do that, since that is the entire
% claim being made for it.
fprintf('\n(a co-pol-only estimator would read %.3f at alpha=%.0f deg)\n', ...
  DLAM_TRUE*cos(2*(THETA_TRUE - deg2rad(80))), 80);

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_quadpol:failed', '%d azimuth(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

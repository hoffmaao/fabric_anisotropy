%TEST_ERSHADI Round-trip test of the Ershadi et al. (2022) implementation.
%
% Same synthetic column as test_quadpol - a known principal azimuth and a
% known contrast, measured from several antenna azimuths - so the two
% implementations are held to the same standard on the same truth and any
% difference between them is theirs, not the test's.
%
% The synthetic is built WITHOUT the deramp conjugation, because it is a
% model rather than radar data; the paper makes exactly that distinction
% ("we use Eq. 7 for the models and the conjugate of Eq. 7 for the radar
% data"), so the test passes deramped = false.
%
% Run: matlab -batch "run('opr_fabric/test/test_ershadi.m')"
clear;
rng(11);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
eps_perp = 3.15;
% Forward constant, inverted from the same relation the implementation
% uses, so the test states the physics once and both directions share it.
grad_per_dlam = 2*pi*fc*C.deps / (C.c*1e9 * sqrt(eps_perp));   % rad/m

THETA_TRUE = deg2rad(35);
DLAM_TRUE = 0.30;
Nz = 3000;
z = (0:Nz-1).' * 0.5;
delta = grad_per_dlam * DLAM_TRUE * z;
Nx = 200;

fprintf('truth: theta %.1f deg, dlam %.3f, %.2f fringes over %.0f m\n', ...
  rad2deg(THETA_TRUE), DLAM_TRUE, delta(end)/(2*pi), z(end));

fails = 0;
for alpha_deg = [0 20 35 55 80 125]
  d = THETA_TRUE - deg2rad(alpha_deg);
  cd_ = cos(d); sd = sin(d);
  ex = exp(1i * delta); ey = ones(Nz, 1);
  hh = cd_^2 * ex + sd^2 * ey;
  vv = sd^2 * ex + cd_^2 * ey;
  hv = cd_*sd * (ex - ey);
  r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  na = 0.02;
  S = struct();
  S.hh = hh .* r + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vv = vv .* r + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.hv = hv .* r + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vh = S.hv;

  out = ptt.ershadiFabric(S, z, struct('fc', fc, 'psi_step_deg', 1, ...
    'win_m', 30, 'grad_win_m', 25, 'coh_min', 0.4, 'deramped', false));

  mid = z > 300 & z < 1200;
  % theta is recovered modulo 90 deg in the ANTENNA frame
  want = mod(rad2deg(d), 90);
  th = rad2deg(out.theta(mid));
  got = mod(rad2deg(0.5*angle(mean(exp(2i*deg2rad(2*th))))) / 2, 90);
  err = abs(mod(got - want + 45, 90) - 45);
  dl = median(out.dlam(mid), 'omitnan');

  ok_t = err < 5;
  ok_d = abs(dl - DLAM_TRUE) < 0.03;
  fails = fails + ~(ok_t && ok_d);
  fprintf('alpha %3d: theta %5.1f (want %5.1f, err %4.1f) %s | dlam %.3f %s\n', ...
    alpha_deg, got, want, err, H_tick(ok_t), dl, H_tick(ok_d));
end

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_ershadi:failed', '%d azimuth(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

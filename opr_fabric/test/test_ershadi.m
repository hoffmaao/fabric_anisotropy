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
% A second check blanks a 100 m band of every channel (as ptt.maskBelowBed
% leaves rows below every trace's bed) and requires the chain to ABSTAIN on
% exactly those rows - theta, dlam, Cmag, dP_hh, dP_hv all NaN, never an
% exact 0 deg axis from an all-zero anomaly row - while every row beyond
% the kernels' reach (75 m gradient + 15 m coherence window) is unchanged
% bit for bit and the rows within reach still carry a finite axis on the
% truth (circular mean, as above - single rows at the fringe nulls are
% noise) and a finite dlam (biased low within ~2 sigma of the band edge
% by the one-sided kernel, as at the record's own top; the median over
% the 100 m on either side is held to 0.05). It runs at alpha 20, off the
% symmetry points where the cross-pol channel vanishes altogether.
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
  if alpha_deg == 20, S_20 = S; out_20 = out; want_20 = want; end
  fprintf('alpha %3d: theta %5.1f (want %5.1f, err %4.1f) %s | dlam %.3f %s\n', ...
    alpha_deg, got, want, err, H_tick(ok_t), dl, H_tick(ok_d));
end

% --- abstention on a band of blanked rows (the alpha 20 column)
EOPT = struct('fc', fc, 'psi_step_deg', 1, 'win_m', 30, 'grad_win_m', 25, ...
  'coh_min', 0.4, 'deramped', false);
band = z >= 600 & z <= 700;
S = S_20; out = out_20; want = want_20;
S_nan = S;
for k = {'hh', 'vv', 'hv', 'vh'}
  S_nan.(k{1})(band, :) = NaN;
end
out_n = ptt.ershadiFabric(S_nan, z, EOPT);
reach = z >= 600 - 100 & z <= 700 + 100;
near = reach & ~band;
ok_abs = all(isnan(out_n.theta(band))) && all(isnan(out_n.dlam(band))) ...
  && all(isnan(out_n.Cmag(band, :)), 'all') && all(isnan(out_n.dP_hh(band, :)), 'all') ...
  && all(isnan(out_n.dP_hv(band, :)), 'all');
ok_same = isequaln(out_n.theta(~reach), out.theta(~reach)) ...
  && isequaln(out_n.dlam(~reach), out.dlam(~reach)) ...
  && isequaln(out_n.Cmag(~reach, :), out.Cmag(~reach, :)) ...
  && isequaln(out_n.psi_grad(~reach, :), out.psi_grad(~reach, :)) ...
  && isequaln(out_n.dP_hh(~band, :), out.dP_hh(~band, :)) ...
  && isequaln(out_n.dP_hv(~band, :), out.dP_hv(~band, :));
th_n = rad2deg(out_n.theta(near));
got_n = mod(rad2deg(0.5*angle(mean(exp(2i*deg2rad(2*th_n))))) / 2, 90);
err_n = abs(mod(got_n - want + 45, 90) - 45);
dl_n = median(out_n.dlam(near), 'omitnan');
ok_near = all(isfinite(out_n.theta(near))) && all(isfinite(out_n.dlam(near))) ...
  && err_n < 5 && abs(dl_n - DLAM_TRUE) < 0.05;
ok_band = ok_abs && ok_same && ok_near;
fails = fails + ~ok_band;
fprintf(['blanked 600-700 m: abstains on the band %s | bit-identical beyond reach %s | ' ...
  'within reach theta err %.1f, dlam %.3f %s\n'], H_tick(ok_abs), H_tick(ok_same), ...
  err_n, dl_n, H_tick(ok_near));

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_ershadi:failed', '%d check(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

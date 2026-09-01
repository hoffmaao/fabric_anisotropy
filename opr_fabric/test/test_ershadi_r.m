%TEST_ERSHADI_R The Sect.-3.5 reflection-ratio retrieval recovers known r.
%
% Validates ptt.fujitaModel + ptt.ershadiInverse - the faithful Ershadi et
% al. (2022) chain - on speckle synthetics with a known three-zone profile
% shaped like their EDML result (two anisotropic-scattering zones of
% opposite sign under an isotropic lid):
%
%   L1    0-400 m   dlam 0.05   theta  30 deg   r   0 dB
%   L2  400-800 m   dlam 0.20   theta  30 deg   r +10 dB
%   L3 800-1200 m   dlam 0.20   theta 120 deg   r -10 dB
%
% The truth generator here builds the two-way Jones sandwich DIRECTLY in
% the physical convention - explicit R(theta) diag R(theta)' per layer,
% antennas fixed, channels read off the matrix elements - and the sweep is
% then performed by ptt.ershadiFabric from the quad-pol set exactly as on
% real data. So a sign-convention mismatch anywhere in fujitaModel's R S R'
% sweep or ershadiInverse's use of it FAILS the fit rather than cancelling:
% the repo's rule that the synthetic settles conventions the printed
% equations cannot.
%
% Verdicts:
%   1. forward-model conventions: the deterministic dP_HH field from
%      ptt.fujitaModel matches the noise-free synthetic sweep (corr > 0.99)
%      BEFORE any fitting, so a failure here is a convention bug, not an
%      optimizer one.
%   2. r recovery: fitted r_dB medians land in the right zone windows
%      (+10 within +-3 dB, -10 within +-3 dB, lid within +-2.5 dB).
%   3. theta recovery: circular error < 5 deg per zone.
%   4. null control: with r = 0 everywhere (same seeds), |fitted r_dB| < 2
%      in every zone - speckle and noise must not manufacture anisotropy.
%   5. eq.-(13) analytic r agrees with the fit where nodes resolve
%      (median |r13 - fit| < 4 dB in the strong zones).
%
% Run: matlab -batch "run('opr_fabric/test/test_ershadi_r.m')"
clear;
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

rng(11);
C = ptt.constants();
fc = 750e6;
eps_perp = 3.15;                        % the paper's value, as ershadiFabric
c0 = C.c * 1e9;
gpd1 = pi * fc * C.deps / (sqrt(eps_perp) * c0);   % ONE-way rad/m per dlam

Nz = 2401; dz = 0.5;
z = (0:Nz-1).' * dz;                    % 0..1200 m
Nx = 360;
% Additive noise level. The co-pol NODES carry the r information, and the
% moment-averaged amplitude floor sits near sqrt(2)*NA: at 0.35 that is a
% -3 dB clip which fills every node the +-10 dB truth needs at -20 dB, and
% the first run measured exactly that failure mode (fitted +8.5/-5.4).
% 0.08 puts the floor near -16 dB, below the discretized node depths, so
% the verdicts test the ESTIMATOR rather than the clip. Noise robustness
% is a property of the data regime (Ershadi gate |C| > 0.4 for the same
% reason), not of this recovery test.
NA = 0.08;

TRUTH = struct( ...
  'top',   {0, 400, 800}, ...
  'dlam',  {0.05, 0.20, 0.20}, ...
  'theta', {deg2rad(30), deg2rad(30), deg2rad(120)}, ...
  'r_db',  {0, 10, -10});

ZONES = [200 380; 450 750; 850 1150];   % assert bands, clear of boundaries

% ---- truth generator: physical-convention Jones sandwich, layer by layer
function M = H_truth_M(z, TRUTH, gpd1)
  Nz = numel(z);
  M = complex(zeros(2, 2, Nz));
  tops = [TRUTH.top];
  NL = numel(TRUTH);
  Pacc = eye(2);
  li = 1;
  for iz = 1:Nz
    while li < NL && z(iz) >= tops(li+1)
      d = tops(li+1) - tops(li);
      Pacc = H_layer(TRUTH(li), d, gpd1) * Pacc;
      li = li + 1;
    end
    P = H_layer(TRUTH(li), z(iz) - tops(li), gpd1) * Pacc;
    ra = 10^(TRUTH(li).r_db / 20);
    Rm = [cos(TRUTH(li).theta), -sin(TRUTH(li).theta); ...
          sin(TRUTH(li).theta),  cos(TRUTH(li).theta)];
    G = Rm * diag([1, ra]) * Rm.';
    M(:, :, iz) = P.' * G * P;          % reciprocal two-way, antennas at 0
  end
end
function A = H_layer(L, d, gpd1)
  del = gpd1 * L.dlam * d;
  Rm = [cos(L.theta), -sin(L.theta); sin(L.theta), cos(L.theta)];
  A = Rm * diag([exp(+0.5i*del), exp(-0.5i*del)]) * Rm.';
end

function S = H_channels(M, Nz, Nx, NA)
  % common complex speckle per depth/trace (co-located layer reflectors,
  % the Fujita coherent-layer picture), independent additive noise
  g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  S = struct();
  S.hh = squeeze(M(1,1,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vv = squeeze(M(2,2,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.hv = squeeze(M(1,2,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vh = squeeze(M(2,1,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
end

M = H_truth_M(z, TRUTH, gpd1);
S = H_channels(M, Nz, Nx, NA);

EOPT = struct('fc', fc, 'psi_step_deg', 1, 'deramped', false);
fr = ptt.ershadiFabric(S, z, EOPT);

% ---- verdict 1: forward-model conventions against the noise-free sweep
layersT = struct('top_m', {TRUTH.top}, 'dlam', {TRUTH.dlam}, ...
  'theta', {TRUTH.theta}, 'r_db', {TRUTH.r_db});
fm = ptt.fujitaModel(layersT, z, fr.psi, struct('fc', fc, ...
  'eps_perp', eps_perp));
% noise-free synthetic sweep of dP_HH straight from M, same E1 sense as
% ershadiFabric applies to data
cc = cos(-fr.psi); ss = sin(-fr.psi);
shh0 = squeeze(M(1,1,:)) * (cc.*cc) + squeeze(M(2,2,:)) * (ss.*ss) + ...
  (squeeze(M(1,2,:)) + squeeze(M(2,1,:))) * (cc.*ss);
A0 = abs(shh0);
dP0 = 20*log10(max(A0, realmin) ./ max(mean(A0, 2), realmin));
band = z >= 200 & z <= 1150;
v = corrcoef(fm.dP_hh(band, :), dP0(band, :));
ok_fwd = v(1, 2) > 0.99;
fprintf('1. fujitaModel matches the physical-convention sweep: corr %.4f  %s\n', ...
  v(1, 2), H_tick(ok_fwd));

% ---- the inversion
IOPT = struct('interval_m', 100, 'z_fit', [200 1150], ...
  'w_theta', [1 0 0], 'w_r', [0 1 0], 'fc', fc, 'eps_perp', eps_perp);
inv = ptt.ershadiInverse(fr, z, IOPT);

function m = H_zone_med(v, z, zn)
  m = median(v(z >= zn(1) & z <= zn(2)), 'omitnan');
end

r_med = zeros(3, 1); th_err = zeros(3, 1);
for k = 1:3
  r_med(k) = H_zone_med(inv.r_db, z, ZONES(k, :));
  tt = H_zone_med(inv.theta, z, ZONES(k, :));
  th_err(k) = abs(angle(exp(2i*(tt - TRUTH(k).theta)))) / 2;
end
fprintf(['   fitted r_db per zone: %+5.1f %+5.1f %+5.1f ' ...
  '(truth 0 / +10 / -10)\n'], r_med);
fprintf('   theta error per zone: %4.1f %4.1f %4.1f deg\n', rad2deg(th_err));
ok_r = abs(r_med(1)) < 2.5 && abs(r_med(2) - 10) < 3 && abs(r_med(3) + 10) < 3;
ok_th = all(th_err < deg2rad(5));
fprintf('2. r recovered in all three zones:                      %s\n', H_tick(ok_r));
fprintf('3. theta recovered in all three zones:                  %s\n', H_tick(ok_th));

% ---- verdict 4: null control, same seeds, r = 0 everywhere
rng(11);
TRUTH0 = TRUTH;
for k = 1:3, TRUTH0(k).r_db = 0; end
M0 = H_truth_M(z, TRUTH0, gpd1);
S0 = H_channels(M0, Nz, Nx, NA);
fr0 = ptt.ershadiFabric(S0, z, EOPT);
inv0 = ptt.ershadiInverse(fr0, z, IOPT);
r0_med = arrayfun(@(k) H_zone_med(inv0.r_db, z, ZONES(k, :)), 1:3);
fprintf('   null-control r_db per zone: %+5.1f %+5.1f %+5.1f\n', r0_med);
ok_null = all(abs(r0_med) < 2);
fprintf('4. no spurious anisotropy from speckle and noise:       %s\n', H_tick(ok_null));

% ---- verdict 5: eq.-(13) analytic agrees with the fit where it resolves
d13 = abs(inv.r13_db - inv.r_db);
m2 = z >= ZONES(2,1) & z <= ZONES(2,2) & isfinite(d13);
m3 = z >= ZONES(3,1) & z <= ZONES(3,2) & isfinite(d13);
ok_13 = true; msg = '';
for zi = {m2, m3}
  m = zi{1};
  if nnz(m) >= 20
    ok_13 = ok_13 && median(d13(m)) < 4;
    msg = sprintf('%s %.1f dB over %d rows;', msg, median(d13(m)), nnz(m));
  end
end
if isempty(msg), msg = ' (no rows resolved - not asserted)'; end
fprintf('5. eq.-(13) analytic vs fit:%s                %s\n', msg, H_tick(ok_13));

fails = ~ok_fwd + ~ok_r + ~ok_th + ~ok_null + ~ok_13;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_ershadi_r:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

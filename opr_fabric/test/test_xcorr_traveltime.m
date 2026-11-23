%TEST_XCORR_TRAVELTIME The cross-correlation travel-time difference, step 2's measurement.
%
% ptt.xcorrTraveltime is the poster's Figure 2: windowed cross-correlation
% of the two eigen-polarisation signals, the lag of best correlation
% tracked down the column, the envelope as the fringe check. It is tested
% on a synthetic that has what the accumulation radar has - a 160 MHz
% band at 300 MHz sampling (a fringe of 1.3 ns under a 6 ns envelope),
% speckle correlated over a few traces, a depolarised component that
% decorrelates the pair, power falling with depth - and a delay that is
% the running integral of a known dlam(z).
%
% The verdicts, in the order that matters:
%
%   1. THE DELAY IS RECOVERED WITHOUT A FRINGE SLIP. The tracked dtau(z)
%      follows the truth within a small fraction of a fringe all the way
%      down; the largest error is far below half a fringe, which is what
%      a slip would cost. The envelope check raises no flag.
%   2. IT SURVIVES DEPOLARISATION. With an uncorrelated component at the
%      power Ridge A shows (coherence ~0.5), the same holds.
%   3. THE UNCERTAINTY IS CALIBRATED. Over noise realisations the z-score
%      of the delay error against out.sigma has rms near one.
%   4. CHAINED INTO STEP 2, dlam(z) COMES BACK. traveltimeFabricML on the
%      measured dtau, with its posterior sigma from verdict 3, recovers
%      the known profile within its stated uncertainty.
%   5. THE WRONG FRAME DEGRADES THE PAIR, AND AT 45 DEG LOSES IT. A pair
%      rotated dth from the eigenframe mixes each signal with the other's
%      delayed copy, so the cross-correlation splits into peaks at 0 and
%      at +-tau. The tracked delay shortens by LESS than cos(2*dth) - the
%      small-phase rule that phase methods obey, and which the step-2
%      docstrings quote, does not hold for a fringe-scale correlation
%      lag - and at 45 deg the two rotated signals are the same mixture
%      and carry no delay at all. The coherence does NOT fall: the
%      rotated pair correlates as well as the aligned one (perfectly at
%      45 deg), so a wrong frame cannot be seen in gamma, only in the
%      delay itself - the delay against azimuth is what locates the
%      eigenframe. This is the poster scheme's reason for fixing the
%      orientation first, measured rather than assumed.
%
% Run: matlab -batch "run('opr_fabric/test/test_xcorr_traveltime.m')"
clear; t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
ptt.estimators.use('nymand');
rng(7);

fs = 300e6; dt = 1e9 / fs;                 % ns
fc = 750e6; B = 160e6; f_off = -13e6;      % LO, bandwidth, and the band's offset from the LO as stored
C = ptt.constants(); n_ice = sqrt(C.eps_bar);
dz = dt * C.c / (2 * n_ice);              % m per sample, two-way
Nt = 5000; Nx = 125;
z = (0:Nt-1).' * dz;                       % ~1400 m
k_ns = C.deps / (n_ice * C.c);             % ns per m per unit dlam (traveltimeFabricML's k)
dl_true = 0.012 + 0.055 * (z / 1400) .^ 1.5 + 0.015 * exp(-((z - 650) / 60) .^ 2);
tau_true = k_ns * cumtrapz(z, dl_true);    % ns, positive: a (slow) arrives later

fprintf('synthetic: %d samples x %d traces, delay at the bottom %.2f ns (%.1f fringes)\n', ...
  Nt, Nx, tau_true(end), tau_true(end) * fc / 1e9);

% --- verdict 1: clean pair
[a, b] = H_pair(z, tau_true, dt, fc, f_off, B, Nx, 0, 15, 3);
o = ptt.xcorrTraveltime(a, b, dt, struct('fc', fc));
e1 = o.dtau - interp1(z, tau_true, (o.idx - 1) * dz);
fprintf('\n1. clean pair: rms err %.3f ns, max |err| %.3f ns (fringe %.2f), env flags %d of %d, gamma median %.2f\n', ...
  rms(e1(~isnan(e1))), max(abs(e1)), o.fringe_ns, nnz(o.env_disagree), numel(o.dtau), median(o.gamma));
ok1 = all(isfinite(o.dtau)) && max(abs(e1)) < 0.25 * o.fringe_ns && ~any(o.env_disagree);
fprintf('   recovered without a slip:                                %s\n', H_tick(ok1));

% --- verdict 2: depolarised pair (coherence ~0.5)
[a, b] = H_pair(z, tau_true, dt, fc, f_off, B, Nx, 1.0, 15, 3);
o2 = ptt.xcorrTraveltime(a, b, dt, struct('fc', fc));
e2 = o2.dtau - interp1(z, tau_true, (o2.idx - 1) * dz);
fprintf('\n2. depolarised pair: rms err %.3f ns, max |err| %.3f ns, env flags %d, gamma median %.2f\n', ...
  rms(e2(~isnan(e2))), max(abs(e2)), nnz(o2.env_disagree), median(o2.gamma));
ok2 = all(isfinite(o2.dtau)) && max(abs(e2)) < 0.35 * o2.fringe_ns && nnz(o2.env_disagree) <= 2;
fprintf('   survives depolarisation:                                 %s\n', H_tick(ok2));

% --- verdict 3: sigma calibration over realisations
NR = 8; zs = [];
for r = 1:NR
  [a, b] = H_pair(z, tau_true, dt, fc, f_off, B, Nx, 1.0, 15, 3);
  orr = ptt.xcorrTraveltime(a, b, dt, struct('fc', fc));
  er = orr.dtau - interp1(z, tau_true, (orr.idx - 1) * dz);
  zs = [zs; er ./ orr.sigma]; %#ok<AGROW>
end
zr = rms(zs(isfinite(zs)));
fprintf('\n3. z-score rms over %d realisations: %.2f (sigma median %.3f ns)\n', NR, zr, median(orr.sigma));
ok3 = zr > 0.6 && zr < 1.6;
fprintf('   uncertainty calibrated:                                  %s\n', H_tick(ok3));

% --- verdict 4: chained into the linear ML
zw = (o2.idx - 1) * dz;
use = isfinite(o2.dtau) & zw > 30;
ml = ptt.estimators.run('traveltimeFabricML', o2.dtau(use), zw(use), ...
  struct('fc', 750e6, 'sigma_ns', o2.sigma(use), 'prior', [0 0.05 150], 'res_min', 0));
dl_t = interp1(z, dl_true, zw(use));
zsc = (ml.dlam - dl_t) ./ ml.sigma_dlam;
deep = zw(use) > 300;
fprintf('\n4. step 2 on the measured delays: median |dlam err| %.4f (deep %.4f), z-score rms %.2f, chi2/dof %.2f\n', ...
  median(abs(ml.dlam - dl_t)), median(abs(ml.dlam(deep) - dl_t(deep))), rms(zsc(deep)), ml.chi2_dof);
ok4 = median(abs(ml.dlam(deep) - dl_t(deep))) < 0.006 && rms(zsc(deep)) < 2.0;
fprintf('   dlam(z) recovered within its posterior:                  %s\n', H_tick(ok4));

% --- verdict 5: the wrong frame degrades the pair, and at 45 deg loses it
[a, b] = H_pair(z, tau_true, dt, fc, f_off, B, Nx, 0, 15, 3);
S = struct('hh', a, 'vv', b, 'hv', zeros(size(a)), 'vh', zeros(size(a)));
T25 = ptt.rotatePolarization(S, deg2rad(25));
T45 = ptt.rotatePolarization(S, deg2rad(45));
o25 = ptt.xcorrTraveltime(T25.hh, T25.vv, dt, struct('fc', fc));
o45 = ptt.xcorrTraveltime(T45.hh, T45.vv, dt, struct('fc', fc));
r25 = o25.dtau(end) / o.dtau(end);
r45 = abs(o45.dtau(end)) / o.dtau(end);
fprintf('\n5. pair rotated 25 deg: delay ratio %.3f (cos 2 dth %.3f), gamma median %.2f vs %.2f aligned\n', ...
  r25, cos(deg2rad(50)), median(o25.gamma), median(o.gamma));
fprintf('   pair rotated 45 deg: |delay| ratio %.3f, gamma median %.2f\n', r45, median(o45.gamma));
ok5 = r25 < 0.98 && r25 > cos(deg2rad(50)) && (isnan(r45) || r45 < 0.1) && ...
  median(o45.gamma) > 0.9 * median(o.gamma);
fprintf('   the wrong frame reads short, 45 deg reads nothing, coherence blind to it:  %s\n', H_tick(ok5));

fprintf('\n%s  (%.0f s)\n', H_tick(ok1 && ok2 && ok3 && ok4 && ok5), toc(t0));
if ~(ok1 && ok2 && ok3 && ok4 && ok5), error('test_xcorr_traveltime: FAILED'); end

% ---------------------------------------------------------------------
function [a, b] = H_pair(z, tau, dt, fc, f_off, B, Nx, depol, snr_db, tau_x)
% Two baseband signals of the same reflectivity in the ENGINEERING sign
% (as OPR stores them): b undelayed and a delayed by tau(z), built in
% chunks over which tau is constant, each chunk shifted in the frequency
% domain and carrying the phase exp(-2 pi i fc tau) that mixing by the LO
% fc leaves on a delayed echo. The band sits f_off from the LO, as the
% accumulation radar's does. depol adds an uncorrelated component of that
% relative power to each; noise at snr_db at the top falling 15 dB to the
% bottom; speckle correlated over tau_x traces.
Nt = numel(z);
Nf = 2 ^ nextpow2(2 * Nt);
f = ifftshift(((0:Nf-1) - floor(Nf / 2)) / (Nf * dt)).';     % GHz
fo = f - f_off / 1e9;
pulse = double(abs(fo) <= B / 2e9) .* (0.5 + 0.5 * cos(pi * fo / (B / 2e9)));
amp = 10 .^ (-(z / z(end)) * 15 / 20);                         % power falling 15 dB
r = H_speckle(Nt, Nx, tau_x) .* amp;
d = H_speckle(Nt, Nx, tau_x) .* amp * sqrt(depol);
d2 = H_speckle(Nt, Nx, tau_x) .* amp * sqrt(depol);
chunk = 40;
a = zeros(Nt, Nx); b = zeros(Nt, Nx);
for s = 1:chunk:Nt
  i = s:min(Nt, s + chunk - 1);
  ti = mean(tau(i));
  ri = zeros(Nt, Nx); ri(i, :) = r(i, :);                      % the chunk at its own depth
  R = fft(ri, Nf, 1) .* pulse;
  b = b + H_crop(ifft(R, [], 1), Nt);
  a = a + H_crop(ifft(R .* exp(-2i * pi * f * ti), [], 1), Nt) * exp(-2i * pi * fc * 1e-9 * ti);
end
a = a + H_crop(ifft(fft(d, Nf, 1) .* pulse, [], 1), Nt);
b = b + H_crop(ifft(fft(d2, Nf, 1) .* pulse, [], 1), Nt);
nz = sqrt(mean(abs(b(1:200, :)) .^ 2, 'all')) * 10 ^ (-snr_db / 20);
a = a + nz * (randn(Nt, Nx) + 1i * randn(Nt, Nx)) / sqrt(2);
b = b + nz * (randn(Nt, Nx) + 1i * randn(Nt, Nx)) / sqrt(2);
end

function x = H_crop(x, Nt)
x = x(1:Nt, :);
end

function s = H_speckle(Nt, Nx, tau_x)
s = (randn(Nt, Nx) + 1i * randn(Nt, Nx)) / sqrt(2);
if tau_x > 1
  k = exp(-(-3*tau_x:3*tau_x) .^ 2 / (2 * (tau_x / 2) ^ 2)); k = k / norm(k);
  s = conv2(1, k, s, 'same');
end
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

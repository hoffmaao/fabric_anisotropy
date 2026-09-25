%TEST_CALIBRATE_CHANNELS Complex channel gains recovered from a frame's own moments.
%
% ptt.calibrateChannels removes the V channel's transmit and receive gains,
% amplitude and phase, from an antenna-frame moment matrix using only what
% the ice supplies: reciprocity for the cross-pol pair, the firn for the
% co-pol amplitude, and the firn co-pol phase's intercept at zero depth
% for the co-pol phase (issue #23).
%
% The column is the one test_quadpol_ls uses - a birefringent single-axis
% column with fully correlated co-pol speckle, a reciprocal antenna
% leakage with a correlated and a decorrelated part, and additive noise -
% seen from the two Ridge A heading families, with the gains measured on
% the real system injected: |b| -3.2 dB, arg a - arg b = +111 deg, and a
% co-pol phase of -45 deg (Ridge A's firn intercept).
%
%   1. RECOVERY. |a| and |b| within 0.3 dB, arg a and arg b within 3 deg,
%      on both headings.
%   2. THE LS ON CALIBRATED MOMENTS. ptt.quadpolFabricLS on the calibrated
%      moments gives dlam within 0.002 and the axis within 1.5 deg of what
%      it gives on the TRUE moments, on both headings - the calibration
%      hands the estimator the ice it would have seen without the gains.
%   3. THE FAMILY GAP. Uncalibrated, the two headings disagree on the axis
%      by several degrees (the mechanism behind the +4 deg / +0.006 dlam
%      Ridge A family split); calibrated, they agree within 1.5 deg, and
%      the LS pedestal's a3 falls to within 0.05 of the true-moment value.
%   4. A PINNED PHASE is honoured: opts.phase_ab overrides the firn
%      intercept and is echoed back with the flag set; and the same phase
%      given a turn away (315 deg for -45) gives the same gains to 1e-9,
%      on the branch in phase with H, not the (-a, -b) one that would
%      mirror every axis.
%
% Run: matlab -batch "run('opr_fabric/test/test_calibrate_channels.m')"
clear;
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
rng(21);

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);
Nz = 2801; z = (0:Nz-1).' * 0.5;                 % 0..1400 m
Nx = 400;
DL = 0.06;
LEAK_C = 0.30 * exp(0.7i); LEAK_D = 0.35; NA = 0.45;
% the gains, as measured on Ridge A (see ptt.calibrateChannels)
A_TRUE = 10^(-0.5/20) * exp(1i * deg2rad(33));
B_TRUE = 10^(-3.2/20) * exp(1i * deg2rad(-78));   % arg a - arg b = 111, arg ab = -45
geoms = {'row', 12.5; 'ns', 81};
OPTS = struct('fc', fc, 'deramped', false, 'dlam_max', 0.25, ...
  'psi_step_deg', 4, 'theta_step_deg', 3, 'step_m', 30);
delta = gpd * DL * z;

fails = 0;
ax_raw = nan(1, 2); ax_cal = nan(1, 2); ax_true = nan(1, 2);
a3_cal = nan(1, 2); a3_true = nan(1, 2);
for gi = 1:2
  d = deg2rad(geoms{gi, 2});
  cd_ = cos(d); sd = sin(d);
  ex = exp(1i * delta); ey = ones(Nz, 1);
  hh = cd_^2 * ex + sd^2 * ey; vv = sd^2 * ex + cd_^2 * ey; hv = cd_*sd * (ex - ey);
  r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  xc = hv .* r + LEAK_C * ((hh + vv)/2) .* r + LEAK_D * g;
  S0 = struct('hh', hh .* r, 'vv', vv .* r, 'hv', xc, 'vh', xc);   % the ice, gain-free
  for k = {'hh', 'vv', 'hv', 'vh'}
    S0.(k{1}) = S0.(k{1}) + NA * (randn(Nz, Nx) + 1i*randn(Nz, Nx));
  end
  % the instrument: VV by ab, HV (tx V, rx H) by a, VH by b
  S = S0;
  S.vv = (A_TRUE * B_TRUE) * S0.vv; S.hv = A_TRUE * S0.hv; S.vh = B_TRUE * S0.vh;
  M0 = ptt.quadpolMoments(S0, [1 Nx]);
  M = ptt.quadpolMoments(S, [1 Nx]);
  [Mc, gc, info] = ptt.calibrateChannels(M, z);

  ea = 20*log10(abs(gc(3)) / abs(A_TRUE)); eb = 20*log10(abs(gc(4)) / abs(B_TRUE));
  pa = rad2deg(angle(gc(3) / A_TRUE)); pb = rad2deg(angle(gc(4) / B_TRUE));
  ok1 = abs(ea) < 0.3 && abs(eb) < 0.3 && abs(pa) < 3 && abs(pb) < 3;
  fprintf('1 (%s). gains: |a| err %+.2f dB, |b| err %+.2f dB, arg a err %+.1f, arg b err %+.1f deg (firn intercept %.1f, rms %.1f deg): %s\n', ...
    geoms{gi, 1}, ea, eb, pa, pb, info.phase_fit.intercept_deg, info.phase_fit.rms_deg, H_tick(ok1));
  fails = fails + ~ok1;

  o_true = ptt.quadpolFabricLS(struct('M', M0), z, OPTS);
  o_cal = ptt.quadpolFabricLS(struct('M', Mc), z, OPTS);
  o_raw = ptt.quadpolFabricLS(struct('M', M), z, OPTS);
  mz = o_true.zw > 800 & o_true.zw < 1400; mt = o_true.zw > 400;
  axis_of = @(o) mod(rad2deg(angle(mean(exp(2i * o.theta0(mt & isfinite(o.theta0))))))/2, 180);
  dl_true = median(o_true.dlam(mz), 'omitnan'); dl_cal = median(o_cal.dlam(mz), 'omitnan');
  ax_true(gi) = axis_of(o_true); ax_cal(gi) = axis_of(o_cal); ax_raw(gi) = axis_of(o_raw);
  a3_true(gi) = o_true.pedestal(3); a3_cal(gi) = o_cal.pedestal(3);
  dax = abs(mod(ax_cal(gi) - ax_true(gi) + 90, 180) - 90);
  ok2 = abs(dl_cal - dl_true) < 0.002 && dax < 1.5;
  fprintf('2 (%s). LS on calibrated vs true moments: dlam %.4f vs %.4f, axis %.1f vs %.1f (raw %.1f): %s\n', ...
    geoms{gi, 1}, dl_cal, dl_true, ax_cal(gi), ax_true(gi), ax_raw(gi), H_tick(ok2));
  fails = fails + ~ok2;

  if gi == 1
    [~, gp, ip] = ptt.calibrateChannels(M, z, struct('phase_ab', deg2rad(-45)));
    [~, gw] = ptt.calibrateChannels(M, z, struct('phase_ab', deg2rad(315)));
    ok4 = ip.phase_ab_pinned && abs(ip.phase_ab_deg + 45) < 1e-9 ...
      && abs(rad2deg(angle(gp(2))) + 45) < 1e-6 && max(abs(gw - gp)) < 1e-9 ...
      && real(gp(3) / abs(gp(3)) + gp(4) / abs(gp(4))) > 0;
    fprintf('4. a pinned co-pol phase is applied and echoed, a turn away gives the same gains (%.1e): %s\n', ...
      max(abs(gw - gp)), H_tick(ok4));
    fails = fails + ~ok4;
  end
end
% the two headings look at the same ice: their axes are 12.5 and 81 deg
% from H, i.e. the same geographic axis; compare in the ICE frame
gap_raw = abs(mod((ax_raw(2) - 81) - (ax_raw(1) - 12.5) + 90, 180) - 90);
gap_cal = abs(mod((ax_cal(2) - 81) - (ax_cal(1) - 12.5) + 90, 180) - 90);
ok3 = gap_cal < 1.5 && max(abs(a3_cal - a3_true)) < 0.05;
fprintf('3. family axis gap raw %.1f deg -> calibrated %.1f deg; pedestal a3 calibrated [%+.3f %+.3f] vs true [%+.3f %+.3f]: %s\n', ...
  gap_raw, gap_cal, a3_cal, a3_true, H_tick(ok3));
fails = fails + ~ok3;

fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_calibrate_channels:failed', '%d check(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

%TEST_QUADPOL_BED_VOTE Sub-bed windows must not vote on the held axis.
%
% The record runs past the bed, and nothing below it is ice. Measured on the
% constant-orientation products
% (23 Sep 2026), of the windows that cleared the vote gates 95% at Eastwind
% and McMurdo and 35% at Taylor Dome lay below the bed, and the held axis
% landed a median 51 / 41 / 9 deg off the axis the ice above the bed gives.
% A sub-bed window is not noise either. Measured on the block moments of the
% captured Taylor Dome frame 20260106_02_003, 50-150 m below the picked bed
% the cross-pol pair is still coherent at 0.97 and the HH-VV field reaches
% |C| 0.5-0.6 with strong azimuthal structure (the bed return's tail and
% off-nadir energy); at the Eastwind shelf the base is 17 dB above the ice
% and a basal multiple follows near twice the thickness. Coherent,
% polarised, and not fabric - the model reads it as fabric anyway, and
% with its own axis. Pure receiver noise below the bed does NOT reproduce
% the capture (tried first: held axes stayed within 1.3 deg), because it
% votes at random and averages out.
%
% Synthetic: a two-segment frame whose ice ends at 380 m under the first
% segment and runs 450 -> 600 m along the second (the bed varies inside a
% segment at every thin-ice site), a fabric axis of 30 deg from the H
% antenna above the bed, and below it a weaker coherent return with its
% own polarisation structure (an axis of 75 deg), standing in for the bed
% tail and multiples, the record running to 900 m as the pipeline's fixed
% z_max makes it.
%
%   1. REPRODUCE: without a bed the held axis of at least one segment is
%      pulled more than 5 deg off the ice. (If this ever stops failing the
%      synthetic no longer exercises the defect, and the test says so.)
%   2. FIX: given the per-trace bed, every held axis is within 3 deg of
%      the ice, and no window reaching below its segment's MEDIAN bed
%      carries a contrast, an axis or a weight - th_seg NaN and q_seg 0,
%      never filled from the neighbouring segment or the frame profile
%      (above it, the traces already past their own bed are blanked,
%      ptt.maskBelowBed).
%   3. SUPERSET: an empty bed changes nothing, bit for bit.
%   4. THE HANDOFF: ptt.thetaProfileAt midway between the two held
%      segments, whose beds differ, is ONE axis at every depth, including
%      the depths below the shallower segment's bed - a block there
%      interpolates two constants and does not switch to the deeper
%      segment's axis alone at the shallower one's bed.
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_bed_vote.m')"
clear;
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
rng(5);

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);
dz = 0.5; z = (0:dz:900).'; Nz = numel(z);
Nx = 800; dx = 5;                          % 4 km: two 2 km segments
x_along = (0:Nx-1) * dx;
az_tr = zeros(1, Nx);                      % straight line, heading 0
TH = deg2rad(30);
DL = 0.03;
BED = [380 450 600];                       % seg 1 flat; seg 2 ramps 450 -> 600
LEAK_C = 0.28 * exp(0.6i); LEAK_D = 0.30; NA = 0.35;

delta = gpd * DL * z;
c = cos(TH); s = sin(TH);
ex = exp(1i * delta);
hh0 = c^2 * ex + s^2; vv0 = s^2 * ex + c^2; hv0 = c*s * (ex - 1);
S = struct('hh', zeros(Nz, Nx), 'vv', zeros(Nz, Nx), 'hv', zeros(Nz, Nx), 'vh', zeros(Nz, Nx));
z_bed = zeros(1, Nx);
% the coherent sub-bed return: its own axis, 10 dB below the ice
TH_SUB = deg2rad(75); AMP_SUB = 10^(-10/20);
cb = cos(TH_SUB); sb = sin(TH_SUB);
hh1 = cb^2 * ex + sb^2; vv1 = sb^2 * ex + cb^2; hv1 = cb*sb * (ex - 1);
for j = 1:Nx
  if x_along(j) < 2000
    zb = BED(1);
  else
    zb = BED(2) + (BED(3) - BED(2)) * (x_along(j) - 2000) / (x_along(end) - 2000);
  end
  z_bed(j) = zb;
  ice = z < zb;
  r = (randn(Nz, 1) + 1i*randn(Nz, 1)) / sqrt(2);
  g = (randn(Nz, 1) + 1i*randn(Nz, 1)) / sqrt(2);
  a = ice + AMP_SUB * ~ice;
  hh = hh0 .* ice + hh1 .* ~ice; vv = vv0 .* ice + vv1 .* ~ice; hv = hv0 .* ice + hv1 .* ~ice;
  xc = (hv .* r + LEAK_C * ((hh + vv)/2) .* r + LEAK_D * g) .* a;
  S.hh(:, j) = hh .* r .* a; S.vv(:, j) = vv .* r .* a; S.hv(:, j) = xc; S.vh(:, j) = xc;
end
for k = {'hh', 'vv', 'hv', 'vh'}
  S.(k{1}) = S.(k{1}) + NA * (randn(Nz, Nx) + 1i*randn(Nz, Nx));
end

OPTS = struct('fc', fc, 'deramped', false, 'dlam_max', 0.25, ...
  'theta_const', true, 'psi_step_deg', 4, 'psi_step_seg_deg', 4, ...
  'theta_step_deg', 3, 'step_m', 20);
axis_err = @(th) abs(mod(rad2deg(th) - rad2deg(TH) + 90, 180) - 90);
held = @(fp) arrayfun(@(k) mod(median(fp.th_seg(isfinite(fp.th_seg(:, k)), k)), pi), 1:fp.nseg);

fails = 0;
% 1. without the bed
fp0 = ptt.quadpolFrameTheta(S, z, az_tr, x_along, OPTS);
e0 = axis_err(held(fp0));
fprintf('no bed:   held axis error per segment [%s] deg\n', sprintf(' %.1f', e0));
ok1 = any(e0 > 5);
fails = fails + ~ok1;
fprintf('  1. the synthetic reproduces the capture (some segment > 5 deg): %s\n', H_tick(ok1));

% 2. with the per-trace bed
fp1 = ptt.quadpolFrameTheta(S, z, az_tr, x_along, setfield(OPTS, 'z_bed', z_bed)); %#ok<SFLD>
e1 = axis_err(held(fp1));
fprintf('with bed: held axis error per segment [%s] deg\n', sprintf(' %.1f', e1));
ok2 = all(e1 < 3);
half = 30;                                 % win_fit_m / 2
leak = false;
seg_med = [BED(1), median(z_bed(x_along >= 2000))];
for k = 1:fp1.nseg
  deep = fp1.zw + half > seg_med(k) - 20;
  leak = leak || any(isfinite(fp1.dlam_seg(deep, k))) ...
    || any(isfinite(fp1.th_seg(deep, k))) || any(fp1.q_seg(deep, k) > 0);
end
ok2b = ~leak;
fails = fails + ~ok2 + ~ok2b;
fprintf('  2a. every held axis within 3 deg of the ice: %s\n', H_tick(ok2));
fprintf('  2b. no sub-bed window carries a contrast, an axis or a weight: %s\n', H_tick(ok2b));

% 3. superset
fp2 = ptt.quadpolFrameTheta(S, z, az_tr, x_along, setfield(OPTS, 'z_bed', [])); %#ok<SFLD>
ok3 = isequaln(fp2.th_seg, fp0.th_seg) && isequaln(fp2.dlam_seg, fp0.dlam_seg) ...
  && isequaln(fp2.ped_ant, fp0.ped_ant);
fails = fails + ~ok3;
fprintf('  3. empty bed is bit-identical to no bed: %s\n', H_tick(ok3));

% 4. the handoff between the two held segments
pg = ptt.thetaProfileAt(fp1, 2000);
dev = max(abs(angle(exp(2i * (pg.theta - pg.theta(1)))))) / 2;
ok4 = max(pg.z) > BED(1) && dev < 1e-9;
fails = fails + ~ok4;
fprintf('  4. handoff midway: one axis to %.0f m (the shallower bed %.0f m), spread %.1e rad: %s\n', ...
  max(pg.z), BED(1), dev, H_tick(ok4));

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_quadpol_bed_vote:failed', '%d check(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

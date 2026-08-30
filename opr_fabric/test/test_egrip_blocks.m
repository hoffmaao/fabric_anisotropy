%TEST_EGRIP_BLOCKS Lateral dlam variation at EastGRIP: what breaks, what holds.
%
% The first EastGRIP validation frame (20240619_01_001) inverted end to end
% but failed internally: block-to-block dlam spread 7.2x against Ridge A's
% 1.1x, resid 0.60 against 0.21, |C_hhvv| 0.17 against 0.42. The suspect was
% the along-track block, since 125 traces is ~125 m at Ridge A's ~1 m trace
% spacing but ~354 m at EastGRIP's 2.83 m, and EastGRIP's fabric is strong
% (EGRIP core: horizontal eigenvalue difference ~0.2-0.35).
%
% THE MECHANISM. Pooling traces whose birefringent phase ramps differ by
% ddlam builds a phase spread  sigma = gpd * ddlam * z  across the pooled
% aperture, and coherence falls as exp(-sigma^2/2). It applies to a BLOCK
% (ddlam across its own traces) and to the FRAME pass (ddlam across the
% whole frame), which is why one ramp can damage both.
%
% THIS TEST RUNS THE SEASON'S REAL GEOMETRY: 2.83 m trace spacing, from the
% CSARP_reference_trajectory rebuild. It previously ran the shipped qlook
% coordinates, which place the traverse in the Gulf of Guinea and inflate
% along-track distance 3.3x to ~9 m per trace. That mattered more than a
% wrong axis label: at 9 m the frame pass returned a CONFIDENTLY WRONG axis
% (~97 deg off) and the test asserted that no block size could rescue it.
% At 2.83 m that does not happen at any ramp rate - see verdict 1 - so the
% old assertions were pinning an artifact of the inflated geometry. The
% ramp rate is quoted per km, which is spacing-independent, so the numbers
% below survive any future re-measurement of the spacing.
%
% WHAT THIS GUARDS, all three measured at the real geometry:
%
%  1. THE FRAME PASS ABSTAINS RATHER THAN LIES. As lateral variation grows
%     the pooled fit gives up - finite window fraction falls - instead of
%     returning a confident wrong axis. Measured 0.010 -> 0.020 dlam/km:
%     finite 35% -> 5%, axis error never worse than 6.8 deg. This is the
%     project's abstain-rather-than-report-the-prior rule holding where it
%     matters most, because the pipeline hands the frame axis to EVERY
%     block: a lie propagates to a whole section, abstention does not.
%     The verdict asserts the PAIR, not coverage alone: no rate may
%     combine a profile the pipeline would hand on with an inaccurate
%     axis. "Would hand on" is the consumer's own rule, not a threshold
%     tuned to these numbers - ptt.quadpolFrameTheta keeps a segment's
%     profile only when at least min_seg_windows (5) theta windows are
%     live, and ptt.thetaProfileAt then interpolates that profile into
%     every block of the segment; below it the segment is marked dead and
%     the blocks fall back. So a row with >= 5 live windows in the fitted
%     band is a section-wide handoff and its axis must be right, while a
%     row below that is abstention and a wrong axis there costs nothing.
%     Measured, the axis term is live rather than short-circuited: 0.010
%     and 0.014 dlam/km both hand on (7 live windows each) at 3.2 and
%     5.5 deg, so the 10 deg bound binds on two of the three rows; only
%     0.020, at one live window, is excused as abstention.
%     Nothing else guards this, and a change that made the fit confident
%     under lateral variation would be a silent, section-wide regression.
%
%  2. THE POOLED HANDOFF IS LOAD-BEARING. It needs no dead frame to hurt.
%     At a MODEST 0.010 dlam/km the frame axis is only ~3 deg off, yet
%     blocks handed that profile recover dlam ~12x worse than blocks handed
%     the true axis (med|err| 0.025 vs 0.002). Frame-axis quality dominates
%     block dlam accuracy by an order of magnitude; the laterally-segmented
%     frame pass (ptt.quadpolFrameTheta) is the fix, and
%     test_quadpol_segmented.m asserts its rescue.
%
%  3. BLOCK LENGTH HAS MARGIN, AND 125 m IS THE RIGHT TARGET. With the axis
%     held true so length is the only variable, sigma alone governs the
%     collapse - it sets in between sigma ~6 and ~12 at z = 1100 m, and two
%     different (length, rate) pairs at sigma = 11.7 fail identically. So
%     L_max ~ 6/(gpd*rate*z_max), which scales as 1/L: a 125 m block is
%     safe to ~0.145 dlam/km while the 354 m block that 125 traces spans
%     at 2.83 m is safe only to ~0.051. The re-cut therefore buys exactly
%     the length ratio, 354/125 = 2.83x, and what that 2.83x is worth here
%     is the headroom over EastGRIP's expected ~0.01 dlam/km: 14.5x
%     instead of 5.1x. That is what run_quadpol_pipeline's
%     BLK_TARGET_M = 125 m re-cut buys, and this verdict is what
%     justifies it.
%
% The dlam CAP - the other EastGRIP culprit - is guarded by
% test_egrip_cap.m; the cap is raised here so it cannot mask these.
%
% Run: matlab -batch "run('opr_fabric/test/test_egrip_blocks.m')"
clear;
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);   % rad/m per unit dlam

THETA_TRUE = deg2rad(35);      % axis, fixed along track (antennas at 0)
Nz = 2800; dz = 0.5;
z = (0:Nz-1).' * dz;           % 0..1400 m
DX = 2.83;                     % EastGRIP trace spacing, reference-trajectory
NX = 1000;                     % 2.83 km for the frame passes; the real
x = (0:NX-1) * DX;             % frame is ~5.8 km
DL0 = 0.25;                    % EGRIP-core-like contrast

% ptt.quadpolFrameTheta keeps a segment profile - and ptt.thetaProfileAt
% then hands it to every block in the segment - only when at least this
% many theta windows are live. It is the consumer's definition of a
% usable handoff, and verdict 1 borrows it rather than inventing one.
MIN_SEG_W = 5;
% Above this the axis is wrong enough to matter: verdict 2 measures a
% 3 deg error costing 12x in block dlam, so 10 deg is already a bad
% handoff - and it is far below the ~97 deg lie the 9 m geometry produced.
MAX_AXIS_DEG = 10;
% Smallest number of block apertures a verdict-3 median may rest on.
NMIN_BLK = 12;

LEAK_C = 0.30 * exp(0.7i);
LEAK_D = 0.35;
NA = 0.5;

% dlam_max raised so the cap (test_egrip_cap.m's subject) cannot mask the
% mechanisms this test exists to pin down
OPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.45, 'deramped', false, ...
  'theta_step_deg', 4);

fprintf(['geometry: %.2f m spacing, %d traces (%.2f km), 125 m block = ' ...
  '%d traces\n\n'], DX, NX, NX*DX/1000, round(125/DX));

%% verdict 1: the frame pass abstains rather than lies
% Three rates spanning the transition. Reported together because they are
% only meaningful together: abstention shows as a FALLING finite fraction,
% a lie as a handed-on profile with a large angular error. Only the second
% is dangerous, and it is what the 9 m geometry used to produce.
RATES = [1.0 1.4 2.0] * 1e-5;          % dlam per metre
fin = nan(size(RATES)); ang = nan(size(RATES));
nlive = zeros(size(RATES));
TH_DEG = rad2deg(THETA_TRUE);
fprintf('%-9s %10s %9s %8s %7s %11s\n', ...
  'rate/km', 'frame_ddl', 'sig@800', 'finite', 'handed', 'axis error');
for k = 1:numel(RATES)
  [S, DL_X] = H_synth(RATES(k), NX, DX, Nz, z, THETA_TRUE, DL0, ...
    gpd, LEAK_C, LEAK_D, NA);
  fr = ptt.quadpolFabricLS(S, z, OPTS);
  m = fr.zw > 300 & fr.zw < 1100;
  ok = isfinite(fr.theta0);
  nlive(k) = nnz(ok & m);
  fin(k) = mean(ok(m));
  if nlive(k) > 0
    th = mod(rad2deg(angle(mean(exp(2i*fr.theta0(ok & m)))))/2, 180);
    ang(k) = min(abs(th - TH_DEG), 180 - abs(th - TH_DEG));
    as = sprintf('%.1f deg', ang(k));
  else
    ang(k) = NaN; as = '(all NaN)';
  end
  fprintf('%-9.3f %10.4f %9.1f %7.0f%% %7s %11s\n', RATES(k)*1000, ...
    RATES(k)*x(end), gpd*RATES(k)*x(end)*800, 100*fin(k), ...
    H_yn(nlive(k) >= MIN_SEG_W), as);
  % only the first rate's frame is reused (verdict 2); holding all three
  % would pin ~180 MB of synthetic per rate for nothing
  if k == 1, S1 = S; DLX1 = DL_X; fr1 = fr; end
end
clear S DL_X fr
% Coverage must FALL as the ramp grows, and no rate may pair a handoff the
% pipeline would propagate with an inaccurate axis. The first rate must
% itself be a live handoff, which is what keeps the axis term binding: a
% dead estimator cannot satisfy the verdict by abstaining everywhere.
handed = nlive >= MIN_SEG_W;
lies = handed & ~(isfinite(ang) & ang < MAX_AXIS_DEG);
ok_abstain = handed(1) && fin(end) < fin(1) && ~any(lies);
fprintf(['\nhanded on (>= %d live windows) at %d of %d rates; axis bound ' ...
  '%.0f deg binds on those\n'], MIN_SEG_W, nnz(handed), numel(RATES), ...
  MAX_AXIS_DEG);
fprintf('1. frame pass abstains rather than lies:            %s\n', ...
  H_tick(ok_abstain));

%% verdict 2: the pooled handoff is load-bearing even when the frame is OK
% The modest rate, whose frame axis is only ~3 deg off - no dead frame
% needed for the coupling to cost an order of magnitude.
okw = isfinite(fr1.theta0);
th_frame = struct('z', fr1.zw(okw), 'theta', fr1.theta0(okw));
% the true axis goes in as the SCALAR the solver documents, not a
% hand-built profile: a hand-built window grid has to match the solver's
% own (z(1)+half : step_m : z(end)-half) to mean what it looks like, and
% a mismatch is invisible only while the value happens to be constant
e_frame = H_blocks(S1, z, DLX1, NX, DX, OPTS, th_frame, fr1.pedestal, ...
  125, NMIN_BLK);
e_true = H_blocks(S1, z, DLX1, NX, DX, OPTS, THETA_TRUE, fr1.pedestal, ...
  125, NMIN_BLK);
clear S1 DLX1
fprintf(['\n125 m blocks at %.3f dlam/km: med|err| %.3f with the FRAME ' ...
  'axis, %.3f with the TRUE axis\n'], RATES(1)*1000, e_frame, e_true);
ok_handoff = e_frame > 5 * e_true;
fprintf('2. frame-axis quality dominates block dlam:         %s\n', ...
  H_tick(ok_handoff));

%% verdict 3: block length margin, with the axis held true
% Length is the only variable here, so this isolates the block mechanism
% from the handoff that verdict 2 covers. No pedestal is handed in either,
% so both lengths are treated identically and the only difference between
% rows is the aperture.
%
% The high rate is 0.050 dlam/km, NOT one of verdict 1's: the collapse is
% governed by sigma, and only sigma >~ 12 at z = 1100 m produces it. A
% 708 m block reaches sigma 4.7 at 0.020 - marginal by design, not failed -
% so asserting failure there would be asserting against the mechanism. At
% 0.050 the same block reaches sigma 11.7 and fails, while 125 m stays at
% 2.1 and does not. That contrast IS the margin BLK_TARGET_M buys.
%
% This verdict gets its OWN, longer section: at NX the 708 m row would
% rest on a median over four apertures, too thin a basis for the verdict
% that justifies BLK_TARGET_M. Its cap is lifted with it - the high rate
% climbs the ramp 0.21 over this length, and at 0.45 the far blocks would
% rail and abstain (quadpolFabricLS returns NaN above 0.98*dlam_max), so
% the row would be measuring the cap instead of the aperture. Both lengths
% see the same cap, so the contrast between them is untouched.
NX_LEN = 1500;                 % 4.24 km
OPTS_L = OPTS; OPTS_L.dlam_max = 0.60;
LEN_RATES = [RATES(1), 5e-5];
fprintf('\n%-9s %9s %10s %7s %10s\n', ...
  'rate/km', 'blk_len', 'sig@1100', 'nblk', 'med|err|');
short_ok = true; long_fails = false;
for k = 1:numel(LEN_RATES)
  rk = LEN_RATES(k);
  [Sk, Dk] = H_synth(rk, NX_LEN, DX, Nz, z, THETA_TRUE, DL0, ...
    gpd, LEAK_C, LEAK_D, NA);
  for L = [125 708]
    [e, nb] = H_blocks(Sk, z, Dk, NX_LEN, DX, OPTS_L, THETA_TRUE, [], ...
      L, NMIN_BLK);
    NB = max(4, round(L/DX));
    fprintf('%-9.3f %8.0fm %10.2f %7d %10.3f\n', rk*1000, NB*DX, ...
      gpd*rk*NB*DX*1100, nb, e);
    if L == 125, short_ok = short_ok && e < 0.010; end
    if L == 708 && k == 2, long_fails = e > 0.05; end
  end
  clear Sk Dk
end
ok_len = short_ok && long_fails;
fprintf('3. 125 m holds where 708 m fails:                   %s\n', ...
  H_tick(ok_len));

fails = ~ok_abstain + ~ok_handoff + ~ok_len;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_egrip_blocks:failed', '%d verdict(s) failed', fails);
end

% ---------------------------------------------------------------- helpers
function [S, DL_X] = H_synth(rate, NX, DX, Nz, z, TH, DL0, gpd, ...
    LEAK_C, LEAK_D, NA)
%H_SYNTH One EastGRIP-like frame with a lateral dlam ramp of `rate` per m.
% Seeded per call so every rate sees the same speckle and the comparison
% across rates is a comparison of the ramp alone.
rng(7);
x = (0:NX-1) * DX;
DL_X = DL0 + rate * x;
cd_ = cos(TH); sd = sin(TH);
S = struct('hh', zeros(Nz, NX), 'vv', zeros(Nz, NX), ...
  'hv', zeros(Nz, NX), 'vh', zeros(Nz, NX));
r = (randn(Nz, NX) + 1i*randn(Nz, NX)) / sqrt(2);
g = (randn(Nz, NX) + 1i*randn(Nz, NX)) / sqrt(2);
for j = 1:NX
  ex = exp(1i * gpd * DL_X(j) * z);
  hh = cd_^2 * ex + sd^2;
  vv = sd^2 * ex + cd_^2;
  hv = cd_*sd * (ex - 1);
  S.hh(:, j) = hh .* r(:, j);
  S.vv(:, j) = vv .* r(:, j);
  xc = hv .* r(:, j) + LEAK_C * ((hh + vv)/2) .* r(:, j) + LEAK_D * g(:, j);
  S.hv(:, j) = xc;
  S.vh(:, j) = xc;
end
for f = {'hh','vv','hv','vh'}
  S.(f{1}) = S.(f{1}) + NA*(randn(Nz,NX)+1i*randn(Nz,NX));
end
end

function [medae, nblk] = H_blocks(S, z, DL_X, NX, DX, OPTS, th, ped, L, nmin)
%H_BLOCKS Median |dlam error| over blocks of length L metres, axis handed in.
% Blocks tile the section disjointly, which is what the pipeline cuts.
% When a length is long enough that the disjoint tiling leaves fewer than
% nmin apertures, the stride shrinks until it does not: each block is
% still an independent L-metre fit - the thing the median is over - and
% the alternative, lengthening the section instead, walks the ramp into
% the dlam cap and measures that instead of the aperture.
if nargin < 10 || isempty(nmin), nmin = 1; end
NB = min(max(4, round(L / DX)), NX);
stride = NB;
if floor((NX - NB) / stride) + 1 < nmin
  stride = max(1, floor((NX - NB) / max(nmin - 1, 1)));
end
j0s = 1:stride:(NX - NB + 1);
nblk = numel(j0s);
dl = nan(1, nblk); tr = nan(1, nblk);
o = OPTS; o.theta0 = th;
if ~isempty(ped) && all(isfinite(ped)), o.pedestal = ped; end
for b = 1:nblk
  j0 = j0s(b); j1 = j0 + NB - 1;
  Sb = struct('hh', S.hh(:, j0:j1), 'vv', S.vv(:, j0:j1), ...
    'hv', S.hv(:, j0:j1), 'vh', S.vh(:, j0:j1));
  ob = ptt.quadpolFabricLS(Sb, z, o);
  m = ob.zw > 300 & ob.zw < 1100;
  dl(b) = median(ob.dlam(m), 'omitnan');
  tr(b) = mean(DL_X(j0:j1));
end
ok = isfinite(dl) & isfinite(tr);
medae = median(abs(dl(ok) - tr(ok)));
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function s = H_yn(ok)
if ok, s = 'yes'; else, s = 'no'; end
end

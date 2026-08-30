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
% old assertions were pinning an artifact of the inflated geometry. To
% keep that from recurring, every geometry below is stated in METRES and
% converted to traces through DX - frame aperture, block lengths, section
% length - and ramp rates are quoted per km. So a future re-measurement of
% the spacing changes the trace counts and nothing else: the same lengths,
% the same ramp per km, the same block counts.
%
% WHAT THIS GUARDS, all three measured at the real geometry:
%
%  1. THE FRAME PASS ABSTAINS RATHER THAN LIES. As lateral variation grows
%     the pooled fit gives up - finite window fraction falls - instead of
%     returning a confident wrong axis. Measured 0.010 -> 0.020 dlam/km:
%     coverage in the 300-1100 m band falls 35% -> 5% (7 live band windows
%     -> 1) and live windows over the whole profile fall 16 -> 8 of 34,
%     while the BAND axis error never exceeds 6.8 deg. This is the
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
%     the blocks fall back. Borrowing that rule means counting the same
%     POPULATION it gates on, so the LIVENESS COUNT is taken over the
%     solver's whole window grid (34 windows, 30:40:1350 m). The AXIS
%     ERROR is not: it is taken over the 300-1100 m band, because
%     sigma = gpd*ddlam*z grows with depth, so the shallow windows are
%     trivially right and averaging them in would let a deep-only lie pass
%     under the bound. Two questions, two populations, on purpose - see
%     MIN_SEG_W and MAX_AXIS_DEG below. The axis bound applies only where
%     the band holds enough live windows to average - fewer is the
%     abstention this verdict rewards, not a lie to catch - so the first
%     rate is REQUIRED to be both handed on and scored, which is what
%     stops the bound going hollow. Measured: all three rates keep 16, 14
%     and 8 live windows against the 5 the consumer requires; 0.010 and
%     0.014 keep 7 live band windows each and are held to the 10 deg bound
%     at 3.2 and 5.5 deg, while 0.020 keeps 1 and is exempt.
%     Nothing else guards this, and a change that made the fit confident
%     under lateral variation would be a silent, section-wide regression.
%
%  2. THE POOLED HANDOFF IS LOAD-BEARING. It needs no dead frame to hurt.
%     At a MODEST 0.010 dlam/km the frame axis is only ~3 deg off over the
%     300-1100 m band - the same figure verdict 1 reports for that rate,
%     over the same population, and the one MAX_AXIS_DEG is derived from -
%     yet
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
%     justifies it. Each row's median rests on DISJOINT apertures - 45 at
%     125 m, 8 at 708 m - and 8 is a cap-imposed ceiling, derived where
%     NX_LEN is set.
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
% Every geometry below is fixed in METRES and converted to traces through
% DX, never the other way round. A ramp rate per km only stays meaningful
% under a re-measured spacing if the length it acts over is held fixed
% too: 1000 traces would be 2.83 km today and something else tomorrow,
% which is how the original test came to be pinning 9 m artefacts.
FRAME_M = 2830;                % frame-pass aperture; the real frame ~5.8 km
NX = round(FRAME_M / DX);
x = (0:NX-1) * DX;
DL0 = 0.25;                    % EGRIP-core-like contrast
L_SHORT = 125;                 % BLK_TARGET_M, the pipeline's re-cut
L_LONG = 708;                  % the length verdict 3 drives to collapse
% The depth band every score in this file is taken over: verdict 1's axis
% error, verdict 1's coverage fraction, and the block dlam H_blocks
% medians. One owner, passed to H_blocks, so "the band the blocks consume"
% stays one range rather than two that happen to agree.
Z_BAND = [300 1100];

% Verdict 1 asks two questions of one fit, and they take DIFFERENT window
% populations on purpose. Do not unify them: doing so silently loosens the
% guard, whichever way it is unified.
%
%   IS THIS HANDOFF USABLE? - the consumer's question, so the consumer's
%   population. ptt.quadpolFrameTheta keeps a segment profile, and
%   ptt.thetaProfileAt then interpolates it into every block, only when at
%   least min_seg_windows theta windows are live, gated as
%   nnz(isfinite(theta0)) < MIN_SEG_W over the solver's FULL window grid.
%   Counting the band instead would call a profile dead that the pipeline
%   would happily propagate.
MIN_SEG_W = 5;
%
%   IS IT ACCURATE WHERE IT MATTERS? - the DEEP BAND (Z_BAND), the same
%   window range H_blocks scores block dlam over and verdict 2 measures its
%   12x in. sigma = gpd*ddlam*z grows with depth, so shallow windows are
%   trivially right; averaging them into the axis error pulls it toward
%   truth and would let a deep-only lie pass. The bound is set from
%   verdict 2, which measures a ~3 deg band error already costing 12x in
%   block dlam: 10 deg is unambiguously a bad handoff, and it is far below
%   the ~97 deg lie the 9 m geometry produced.
MAX_AXIS_DEG = 10;
%   ...but only where there IS a band axis to judge. The claim is that no
%   rate pairs USABLE coverage with an INACCURATE axis, so a band holding
%   one or two live windows makes no claim at all: that is the abstention
%   this verdict exists to reward, and failing it on an unaveraged angle -
%   or on a NaN - would punish the estimator for doing the right thing.
%   The minimum is what makes a circular mean mean something rather than
%   what keeps a row passing: five samples is the smallest set whose mean
%   is not swung past the bound by one outlier at the ~6 deg scale these
%   windows scatter over, and it is the same count ptt.quadpolFrameTheta
%   requires of a whole profile before it will treat it as usable.
MIN_BAND_W = 5;
% Smallest number of DISJOINT block apertures a verdict-3 median may rest
% on. Eight, not more, because the dlam cap bounds it - see verdict 3.
NMIN_BLK = 8;

LEAK_C = 0.30 * exp(0.7i);
LEAK_D = 0.35;
NA = 0.5;

% dlam_max raised so the cap (test_egrip_cap.m's subject) cannot mask the
% mechanisms this test exists to pin down
OPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.45, 'deramped', false, ...
  'theta_step_deg', 4);

fprintf(['geometry: %.2f m spacing, %d traces (%.2f km), %d m block = ' ...
  '%d traces\n\n'], DX, NX, NX*DX/1000, L_SHORT, H_ntr(L_SHORT, DX));

%% verdict 1: the frame pass abstains rather than lies
% Three rates spanning the transition. Reported together because they are
% only meaningful together: abstention shows as a FALLING finite fraction,
% a lie as a handed-on profile with a large angular error. Only the second
% is dangerous, and it is what the 9 m geometry used to produce.
RATES = [1.0 1.4 2.0] * 1e-5;          % dlam per metre
fin = nan(size(RATES)); ang = nan(size(RATES));
nlive = zeros(size(RATES)); nband = zeros(size(RATES));
TH_DEG = rad2deg(THETA_TRUE);
fprintf('%-9s %10s %9s %6s %6s %7s %7s %11s\n', 'rate/km', 'frame_ddl', ...
  'sig@800', 'band', 'live', 'handed', 'scored', 'axis error');
for k = 1:numel(RATES)
  [S, DL_X] = H_synth(RATES(k), NX, DX, Nz, z, THETA_TRUE, DL0, ...
    gpd, LEAK_C, LEAK_D, NA);
  fr = ptt.quadpolFabricLS(S, z, OPTS);
  ok = isfinite(fr.theta0);
  % nlive is the CONSUMER's population - every window of the profile
  % thetaProfileAt would interpolate - while fin, nband and ang are the
  % deep band the blocks are scored over; see MIN_SEG_W / MAX_AXIS_DEG
  % above for why the two questions do not share a population
  m = fr.zw > Z_BAND(1) & fr.zw < Z_BAND(2);
  nlive(k) = nnz(ok);
  nband(k) = nnz(ok & m);
  fin(k) = mean(ok(m));
  if nband(k) >= MIN_BAND_W
    th = mod(rad2deg(angle(mean(exp(2i*fr.theta0(ok & m)))))/2, 180);
    ang(k) = min(abs(th - TH_DEG), 180 - abs(th - TH_DEG));
    as = sprintf('%.1f deg', ang(k));
  else
    ang(k) = NaN; as = '(abstains)';
  end
  fprintf('%-9.3f %10.4f %9.1f %5.0f%% %6d %7s %7s %11s\n', RATES(k)*1000, ...
    RATES(k)*x(end), gpd*RATES(k)*x(end)*800, 100*fin(k), nlive(k), ...
    H_yn(nlive(k) >= MIN_SEG_W), H_yn(nband(k) >= MIN_BAND_W), as);
  % only the first rate's frame is reused (verdict 2); holding all three
  % would pin ~180 MB of synthetic per rate for nothing
  if k == 1, S1 = S; DLX1 = DL_X; fr1 = fr; end
end
clear S DL_X fr
% Coverage must FALL as the ramp grows - on the consumer's population and
% in the quoted band alike - and no rate may pair a handoff the pipeline
% would propagate with an inaccurate axis in the band the blocks consume.
% A row whose band has fewer than MIN_BAND_W live windows is exempt from
% the angle: it has no usable band coverage, so there is nothing to lie
% about, and that is the abstention the verdict rewards.
%
% The exemption is what would make the angle term hollow, so the FIRST
% rate - the mildest lateral variation, where a working estimator must
% still deliver - is required to be both handed on and scored. That is not
% the inverted incentive the deep rates get to enjoy: abstaining at
% 0.010 dlam/km would break this test's premise and verdict 2 with it,
% while abstaining at 0.020 is the behaviour being asserted. So at least
% one row is always subject to the bound, by construction rather than by
% today's numbers.
handed = nlive >= MIN_SEG_W;
scored = nband >= MIN_BAND_W;
lies = handed & scored & ~(ang < MAX_AXIS_DEG);
ok_abstain = handed(1) && scored(1) && nlive(end) < nlive(1) && ...
  fin(end) < fin(1) && ~any(lies);
fprintf(['\nhanded on (>= %d live windows of %d) at %d of %d rates, of ' ...
  'which %d have >= %d live band windows and are held to the %.0f deg ' ...
  'axis bound\n'], MIN_SEG_W, numel(fr1.zw), nnz(handed), numel(RATES), ...
  nnz(handed & scored), MIN_BAND_W, MAX_AXIS_DEG);
fprintf('%-51s %s\n', '1. frame pass abstains rather than lies:', ...
  H_tick(ok_abstain));

%% verdict 2: the pooled handoff is load-bearing even when the frame is OK
% The modest rate, whose frame axis is only ~3 deg off over the 300-1100 m
% band the blocks are scored in (verdict 1's axis-error column for that
% rate) - no dead frame needed for the coupling to cost an order of
% magnitude.
okw = isfinite(fr1.theta0);
th_frame = struct('z', fr1.zw(okw), 'theta', fr1.theta0(okw));
% the true axis goes in as the SCALAR the solver documents, not a
% hand-built profile: a hand-built window grid has to match the solver's
% own (z(1)+half : step_m : z(end)-half) to mean what it looks like, and
% a mismatch is invisible only while the value happens to be constant
[e_frame, ~, lb] = H_blocks(S1, z, DLX1, NX, DX, OPTS, th_frame, ...
  fr1.pedestal, L_SHORT, Z_BAND);
e_true = H_blocks(S1, z, DLX1, NX, DX, OPTS, THETA_TRUE, fr1.pedestal, ...
  L_SHORT, Z_BAND);
clear S1 DLX1
fprintf(['\n%.0f m blocks at %.3f dlam/km: med|err| %.3f with the FRAME ' ...
  'axis, %.3f with the TRUE axis\n'], lb, RATES(1)*1000, e_frame, e_true);
ok_handoff = e_frame > 5 * e_true;
fprintf('%-51s %s\n', '2. frame-axis quality dominates block dlam:', ...
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
% that justifies BLK_TARGET_M. The apertures are DISJOINT, which is the
% tiling run_quadpol_pipeline cuts (j0 = (b-1)*NBLK_TR + 1), so the nblk
% column below counts independent fits. Overlapping them to reach a bigger
% number would not make the median any more robust: neighbours sharing
% half their traces share their speckle and their stretch of ramp.
%
% EIGHT is the ceiling the dlam cap allows, and the arithmetic is worth
% recording because it does not depend on the block length. Holding a
% block at the failure sigma needs rate*L = sigma/(gpd*z_max), so the ramp
% climbs sigma/(gpd*z_max) = 11.7/(0.300*1100) = 0.035 dlam per DISJOINT
% block whatever L is. From DL0 = 0.25 up to the rail at 0.98*0.60 that
% leaves room for (0.588 - 0.25)/0.035 = 9.5 blocks, so 8 clears the rail
% with margin and 12 would not. Buying more by lifting the cap further is
% not free - see the grid coupling below.
%
% The section is sized from the LENGTH the verdict needs, never a trace
% count: NMIN_BLK whole long blocks, converted through DX. That makes the
% disjoint count exactly NMIN_BLK at any spacing, so the enough_blk check
% below cannot be broken by a re-measured DX - and the cap arithmetic
% above is already spacing-independent, since NMIN_BLK long blocks climb
% NMIN_BLK*L_LONG*rate metres of ramp however many traces span them.
NX_LEN = NMIN_BLK * H_ntr(L_LONG, DX);   % 5.66 km at 2.83 m
% The cap is lifted from OPTS's 0.45 because this longer section climbs
% the ramp 0.283 at the high rate, reaching 0.533, which would rail
% against 0.45 (quadpolFabricLS returns NaN above 0.98*dlam_max) and make
% the row measure the cap instead of the aperture. But dlam_max is not
% only a cap: quadpolFabricLS builds its initial search as
% dd_grid = linspace(0, dlam_max*gpd, 26), so 0.45 -> 0.60 also coarsens
% that grid from ~0.018 to ~0.024 dlam per node. The e < 0.010 asserted on
% the 125 m rows therefore sits BELOW one grid node and is met by the
% fminsearch refinement, not by the grid - and for the same reason these
% 125 m numbers are not strictly comparable with verdict 2's, which are
% fitted at 0.45. Both lengths here see the same cap, so the contrast
% between them is untouched.
OPTS_L = OPTS; OPTS_L.dlam_max = 0.60;
LEN_RATES = [RATES(1), 5e-5];
fprintf('\n%-9s %9s %10s %7s %10s\n', ...
  'rate/km', 'blk_len', 'sig@1100', 'nblk', 'med|err|');
short_ok = true; long_fails = false; enough_blk = true;
for k = 1:numel(LEN_RATES)
  rk = LEN_RATES(k);
  [Sk, Dk] = H_synth(rk, NX_LEN, DX, Nz, z, THETA_TRUE, DL0, ...
    gpd, LEAK_C, LEAK_D, NA);
  for L = [L_SHORT L_LONG]
    % lb is the aperture H_blocks actually fitted, not a parallel
    % recomputation of it, so the reported length and sigma describe the
    % rows below them
    [e, nb, lb] = H_blocks(Sk, z, Dk, NX_LEN, DX, OPTS_L, THETA_TRUE, ...
      [], L, Z_BAND);
    fprintf('%-9.3f %8.0fm %10.2f %7d %10.3f\n', rk*1000, lb, ...
      gpd*rk*lb*Z_BAND(2), nb, e);
    % the median is only worth its verdict if it rests on NMIN_BLK
    % disjoint apertures; a section too short for the length silently
    % yields fewer, and H_blocks reports 0 rather than one clamped block
    enough_blk = enough_blk && nb >= NMIN_BLK;
    if L == L_SHORT, short_ok = short_ok && e < 0.010; end
    if L == L_LONG && k == 2, long_fails = e > 0.05; end
  end
  clear Sk Dk
end
ok_len = short_ok && long_fails && enough_blk;
if ~enough_blk
  fprintf('   (a row rested on fewer than %d disjoint blocks)\n', NMIN_BLK);
end
fprintf('%-51s %s\n', ...
  sprintf('3. %d m holds where %d m fails:', L_SHORT, L_LONG), ...
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

function nb = H_ntr(L, DX)
%H_NTR Traces spanning L metres at spacing DX - the one place this is done.
% The block aperture, the section length and the printed diagnostics all
% come through here, so a report can never describe a different aperture
% than the one fitted.
nb = max(4, round(L / DX));
end

function [medae, nblk, len_m] = H_blocks(S, z, DL_X, NX, DX, OPTS, th, ...
    ped, L, zband)
%H_BLOCKS Median |dlam error| over blocks of length L metres, axis handed in.
% Blocks tile the section DISJOINTLY, which is what run_quadpol_pipeline
% cuts, so nblk is a count of independent apertures and the median means
% what it looks like. A section too short to hold one whole block returns
% nblk = 0 and medae = NaN rather than one clamped, shorter block, so the
% caller's block-count check fails loudly instead of the row silently
% measuring a different aperture than it reports. len_m is the aperture
% actually fitted, in metres, for the caller to report.
NB = H_ntr(L, DX);
len_m = NB * DX;
nblk = floor(NX / NB);
if nblk < 1, medae = NaN; nblk = 0; return; end
dl = nan(1, nblk); tr = nan(1, nblk);
o = OPTS; o.theta0 = th;
if ~isempty(ped) && all(isfinite(ped)), o.pedestal = ped; end
for b = 1:nblk
  j0 = (b-1)*NB + 1; j1 = j0 + NB - 1;
  Sb = struct('hh', S.hh(:, j0:j1), 'vv', S.vv(:, j0:j1), ...
    'hv', S.hv(:, j0:j1), 'vh', S.vh(:, j0:j1));
  ob = ptt.quadpolFabricLS(Sb, z, o);
  m = ob.zw > zband(1) & ob.zw < zband(2);
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

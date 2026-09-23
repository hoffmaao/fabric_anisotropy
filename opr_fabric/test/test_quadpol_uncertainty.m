%TEST_QUADPOL_UNCERTAINTY Coverage of the quad-pol standard errors.
%
% An uncertainty estimate is only worth reporting if it is CALIBRATED: the
% errors the estimator actually makes must scatter by about the sigma it
% reports. Nothing in the fit itself checks that, so this test measures
% it on synthetic columns whose truth is known, for both resampling
% estimators the pipeline uses:
%
%   A. SEGMENT LEVEL - ptt.quadpolJackknife, delete-one over the sub-blocks
%      whose weighted mean is the segment's pooled moment matrix, in
%      constant-orientation mode (the held axis and the dlam(z) profile).
%      R independent realisations of the same column give the empirical
%      scatter (SD across realisations) of the held axis and of dlam at
%      every window; each realisation's jackknife must reproduce it.
%      Verdicts: median over windows of se_dlam / SD_dlam in [0.6, 1.6];
%      se_theta / SD_theta in [0.5, 2.0] (R axis draws make SD_theta itself
%      ~30% uncertain). BIAS is reported separately and is NOT what the
%      standard error covers: a biased estimator with an honest sigma is
%      a different problem from a dishonest sigma, and this test must not
%      confuse the two. MEASURED (2 Sep 2026, scratch bias probe on this
%      column's variants): the dlam bias is NOISE-INDUCED and positive -
%      +0.0035 at NA = 0.35 with no pedestal, exactly 0 without noise,
%      and +0.0006-0.0010 when the leakage pedestal is present, because
%      the pedestal nuisance terms absorb most of the noise distortion.
%      Real data carry the pedestal, so they sit nearest the last case;
%      a per-SNR characterisation is still owed before dlam is quoted to
%      better than ~0.003.
%
%   B. BLOCK LEVEL - the split-half used for the section blocks: a
%      125-trace block refit on its two halves with the axis held, half
%      the difference taken as one draw of the block's error. Over R
%      block realisations, RMS(hdiff)/2 per window must match the SD of
%      the full-block dlam error: ratio in [0.7, 1.4].
%
%   C. The jackknife must SCALE: pooling twice the sub-blocks should cut
%      se_dlam by about sqrt(2) (ratio of medians in [1.15, 1.75]).
%
%   D. FREE MODE MUST REPORT AT ALL. A and C hold the axis, so they never
%      exercise the per-window theta0 abstention. In free mode the
%      replicates search only +-theta_half_deg around the full-data axis,
%      and q_theta - the cost contrast over whatever grid was given - is a
%      small fraction of its full-range value there, because the contrast
%      is measured across the grid and the grid has been narrowed to sit
%      on the minimum. Gating THAT on the same q_min as a full 0-180
%      search is a category error: it throws away replicate values for
%      being narrow-searched, not for being bad, and it does so silently -
%      n_edge does not cover it, and a window that loses enough replicates
%      returns NaN from the standard error with nothing saying why. See the
%      note in ptt.quadpolJackknife.
%      MEASURED on this column: gated on the raw contrast, 88% of
%      replicate theta0 values survive; with the contrast put back on the
%      full-range scale, 100%. The loss grows as the fabric weakens and
%      q_theta approaches q_min, so 88% here is the mild end of it.
%      Verdict: essentially every window that the full fit gave an axis
%      must carry a finite se_theta, of a plausible size rather than a
%      degenerate zero.
%
%   E. THE CURE IS NOT A CORRECTION, IT IS A HANDOFF. D says a narrow grid
%      must stop costing replicates their axis; it does not say a narrow
%      grid should stop the estimator abstaining at all, and switching the
%      gates off would trade a false abstention for a false detection.
%      Four rounds of rescaling q_theta for the grid span all failed the
%      same way, because max(cost_th) - min(cost_th) is a SAMPLED EXTREME
%      over whatever nodes the grid holds and no multiplicative factor can
%      make that independent of where the minimum sits. The quantity
%      should never have been re-estimated: whether a window carries
%      enough contrast to own an axis is a property of the WINDOW, settled
%      once by the full fit on the full grid with all the data. So
%      ptt.quadpolFabricLS exports that verdict as out.window_ok and
%      ptt.quadpolJackknife hands it back in as opts.window_ok; a
%      replicate reports for exactly the windows its parent accepted and
%      abstains for exactly the ones it rejected, whatever grid it
%      searched. Asserted: a full grid handed in explicitly behaves like
%      no grid at all; the full fit still abstains on a flat column; an
%      accepted window reports on a +-15 deg grid and on the full grid
%      alike; a rejected window reports on neither, so the handoff does
%      not become a false detection; and the resulting se_theta is finite,
%      never exceeds 180 deg (impossible for an axis) and sits inside the
%      52 deg circular-uniform SD on a column whose axis is well
%      determined.
%
%   F. THE STANDARD ERROR ITSELF MUST BE BOUNDED. E removes the CAUSE of
%      scattered replicates, but the statistic that summarises them was
%      still the linear delete-one formula sqrt((m-1)/m * sum(d.^2)), whose
%      ceiling is (pi/2)*sqrt(m-1) - 180 deg at m = 5, 285 at m = 11, 412
%      at m = 22 - so it can report more than a half turn of uncertainty
%      for a quantity defined modulo a half turn whenever replicates
%      genuinely disagree, which no removal of a cause can prevent.
%      ptt.circAxisSE derives the dispersion from a RESULTANT instead,
%      bounded in [0, 1] by construction rather than by a clamp, and
%      abstains where a Rayleigh test cannot separate the replicates from
%      a uniform spread of axes. Asserted directly on the helper, since
%      the adversarial inputs that matter are ones no synthetic column
%      produces on demand: m swept 3..21 over uniform, two-ended and
%      random draws never exceeds the bound; a uniform spread abstains
%      rather than returning a large number; clustered replicates still
%      reproduce the linear delete-one value, so nothing is made
%      optimistic; and the three states a NaN can mean stay
%      distinguishable through m, Rbar and the saturation flag.
%      A bounded statistic also stops DISCRIMINATING before it reaches its
%      bound: measured, m = 22 reported ~36.7 deg for every true scatter
%      from 16 to 40 deg. A value from inside that band is not a
%      measurement, so circAxisSE abstains once the stretched resultant
%      reaches its own sampling noise floor, 1/sqrt(m) - which moves with
%      m as the plateau does, and is the binding one from m = 7 up - and
%      F2b asserts the band is abstained on rather than reported flat,
%      that nothing comes back above the 31.2 deg resolution cap, and that
%      the surviving reports still tell 6 deg of scatter from 12.
%      ACCEPTED LIMITATION, asserted nowhere because it is intended: the
%      resolution floor sits on the sqrt(m-1)-inflated resultant, so the
%      true scatter at which the SE declines to report tightens as m grows
%      - 17.2 deg at m = 6, 12.2 deg at m = 11, 8.4 deg at m = 22. Better
%      sampled segments abstain sooner, so abstention rates are NOT
%      comparable between segments of different length. See ptt.circAxisSE.
%
%   G. NOTHING IN THIS BATCH MAY MOVE THE DEFAULT PATH. Every fix from the
%      span-correction saga onward was scoped to narrow-grid resampling and
%      to the axis SE statistic. The free-mode, full-grid per-window
%      theta0 and dlam path - the one the Ridge A depth movie draws - must
%      be untouched, so it is refitted here with the jackknife off and
%      compared cell by cell against the same fit run without any of the
%      new options set.
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_uncertainty.m')"
clear;
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

rng(17);
C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);

Nz = 2400; dz = 0.5;
z = (0:Nz-1).' * dz;                     % 0..1200 m
NA = 0.35;
LEAK_C = 0.28 * exp(0.6i);
LEAK_D = 0.30;
TH = deg2rad(35);
DL_Z = 0.02 + 0.06 * (z / z(end));

OPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.30, 'deramped', false, ...
  'theta_step_deg', 4, 'theta_const', true);
JOPTS = struct('theta_half_deg', 15, 'theta_step_deg', 3);

function S = H_col(z, th, DL_Z, Nx, gpd, LEAK_C, LEAK_D, NA)
  Nz = numel(z);
  dz = median(diff(z));
  del = 2 * cumsum(DL_Z) * dz * gpd / 2;
  r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  c = cos(th); s = sin(th);
  ex = exp(1i * del);
  hh = c^2 * ex + s^2; vv = s^2 * ex + c^2; hv = c*s * (ex - 1);
  S = struct();
  S.hh = hh .* r; S.vv = vv .* r;
  xc = hv .* r + LEAK_C * ((hh + vv)/2) .* r + LEAK_D * g;
  S.hv = xc; S.vh = xc;
  for f = {'hh','vv','hv','vh'}
    S.(f{1}) = S.(f{1}) + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  end
end

function [Mg, Msub, nsub] = H_subblocks(S, nb, nsb)
  Msub = cell(1, nb); nsub = nsb * ones(1, nb);
  Mg = 0;
  for i = 1:nb
    jj = (i-1)*nsb + (1:nsb);
    Sb = struct('hh', S.hh(:, jj), 'vv', S.vv(:, jj), 'hv', S.hv(:, jj), 'vh', S.vh(:, jj));
    Msub{i} = ptt.quadpolMoments(Sb, [1 nsb]);
    Mg = Mg + Msub{i} * nsb;
  end
  Mg = Mg / sum(nsub);
end

fails = 0;

%% A. segment jackknife, constant mode
R = 6; NB = 11; NSB = 150;
th_hat = nan(1, R); se_th = nan(1, R);
dl_hat = []; se_dl = []; edge = zeros(1, R);
for r = 1:R
  S = H_col(z, TH, DL_Z, NB*NSB, gpd, LEAK_C, LEAK_D, NA);
  [Mg, Msub, nsub] = H_subblocks(S, NB, NSB);
  o = ptt.quadpolFabricLS(struct('M', Mg), z, OPTS);
  J = ptt.quadpolJackknife(Msub, nsub, z, OPTS, o, JOPTS);
  th_hat(r) = o.theta_const; se_th(r) = J.se_theta_c;
  if isempty(dl_hat), dl_hat = nan(numel(o.zw), R); se_dl = dl_hat; zw = o.zw; end
  dl_hat(:, r) = o.dlam; se_dl(:, r) = J.se_dlam; edge(r) = J.n_edge;
  fprintf('  A%d: axis %.2f deg (se %.2f), se_dlam med %.4f, edge %d, %.1f min\n', r, ...
    rad2deg(o.theta_const), rad2deg(J.se_theta_c), median(J.se_dlam, 'omitnan'), J.n_edge, toc(t0)/60);
end
mu_th = angle(sum(exp(2i*th_hat))) / 2;
d_th = angle(exp(2i*(th_hat - mu_th))) / 2;
sd_th = sqrt(sum(d_th.^2) / (R - 1));
bias_th = angle(exp(2i*(mu_th - TH))) / 2;
r_th = median(se_th) / max(sd_th, eps);
dl_true = interp1(z, DL_Z, zw, 'linear', 'extrap');
sd_dl = std(dl_hat, 0, 2, 'omitnan');
bias_dl = mean(dl_hat, 2, 'omitnan') - dl_true;
ok_w = isfinite(sd_dl) & sum(isfinite(dl_hat), 2) >= R - 1;
r_dl = median(median(se_dl(ok_w, :), 2, 'omitnan') ./ sd_dl(ok_w));
fprintf('\nA. SEGMENT JACKKNIFE (R=%d, %d sub-blocks of %d traces)\n', R, NB, NSB);
fprintf('   axis: bias %+.2f deg, empirical SD %.2f deg, jackknife se median %.2f deg -> ratio %.2f\n', ...
  rad2deg(bias_th), rad2deg(sd_th), rad2deg(median(se_th)), r_th);
fprintf('   dlam: |bias| median %.4f, empirical SD median %.4f, se median %.4f -> ratio %.2f (%d windows)\n', ...
  median(abs(bias_dl(ok_w))), median(sd_dl(ok_w)), median(se_dl(ok_w, :), 'all', 'omitnan'), r_dl, nnz(ok_w));
okA_th = r_th >= 0.5 && r_th <= 2.0;
okA_dl = r_dl >= 0.6 && r_dl <= 1.6;
fprintf('A1. axis se within [0.5, 2.0] of empirical scatter:  %s\n', H_tick(okA_th));
fprintf('A2. dlam se within [0.6, 1.6] of empirical scatter:  %s\n', H_tick(okA_dl));
fprintf('    (bias, reported not asserted: axis %+.2f deg, dlam %.4f; search-edge hits %d)\n', ...
  rad2deg(bias_th), median(abs(bias_dl(ok_w))), sum(edge));
fails = fails + ~okA_th + ~okA_dl;

%% C. scaling with twice the sub-blocks
se_big = [];
for r = 1:2
  S = H_col(z, TH, DL_Z, 2*NB*NSB, gpd, LEAK_C, LEAK_D, NA);
  [Mg, Msub, nsub] = H_subblocks(S, 2*NB, NSB);
  o = ptt.quadpolFabricLS(struct('M', Mg), z, OPTS);
  J = ptt.quadpolJackknife(Msub, nsub, z, OPTS, o, JOPTS);
  se_big = [se_big; J.se_dlam(:)]; %#ok<AGROW>
end
r_scale = median(se_dl(:), 'omitnan') / median(se_big, 'omitnan');
okC = r_scale >= 1.15 && r_scale <= 1.75;
fprintf('\nC. SCALING: se_dlam(11 sub-blocks) / se_dlam(22) = %.2f (sqrt 2 = 1.41)\n', r_scale);
fprintf('C1. jackknife se scales with the data:                %s\n', H_tick(okC));
fails = fails + ~okC;

%% D. free-mode jackknife: the replicates must not all abstain
S = H_col(z, TH, DL_Z, NB*NSB, gpd, LEAK_C, LEAK_D, NA);
[Mg, Msub, nsub] = H_subblocks(S, NB, NSB);
FOPTS = OPTS; FOPTS.theta_const = false;
of = ptt.quadpolFabricLS(struct('M', Mg), z, FOPTS);
Jf = ptt.quadpolJackknife(Msub, nsub, z, FOPTS, of, JOPTS);
okf = isfinite(of.theta0);
rep_fin = mean(isfinite(Jf.theta_rep(okf, :)), 'all');
se_fin = mean(isfinite(Jf.se_theta(okf)));
se_med = median(Jf.se_theta(okf), 'omitnan');
fprintf('\nD. FREE-MODE JACKKNIFE (per-window axis, +-%g deg replicate grid)\n', ...
  H_optd(JOPTS, 'theta_half_deg', 15));
fprintf(['   %d of %d windows report an axis; replicate theta finite on ' ...
  '%.0f%% of them\n'], nnz(okf), numel(okf), 100*rep_fin);
fprintf('   se_theta finite on %.0f%% of reporting windows, median %.2f deg\n', ...
  100*se_fin, rad2deg(se_med));
okD = nnz(okf) >= 3 && rep_fin >= 0.98 && se_fin >= 0.98 && ...
  isfinite(se_med) && se_med > 0 && rad2deg(se_med) < 45;
fprintf('D1. free-mode replicates report a standard error:      %s\n', H_tick(okD));
fails = fails + ~okD;

%% E. a replicate inherits its parent's verdict; grid width gates nothing
GFREE = OPTS; GFREE.theta_const = false;
o_def = ptt.quadpolFabricLS(struct('M', Mg), z, GFREE);
GEXP = GFREE; GEXP.theta_grid = (0:4:176) * pi/180;
o_exp = ptt.quadpolFabricLS(struct('M', Mg), z, GEXP);
same_abst = isequal(isfinite(o_def.theta0), isfinite(o_exp.theta0));
finb = isfinite(o_def.theta0) & isfinite(o_exp.theta0);
dth = max(abs(angle(exp(2i*(o_def.theta0(finb) - o_exp.theta0(finb))))/2));
okE1 = same_abst && (isempty(dth) || rad2deg(dth) < 1e-6);

Siso = H_col(z, TH, 1e-4 * ones(size(z)), NB*NSB, gpd, LEAK_C, LEAK_D, NA);
Miso = H_subblocks(Siso, NB, NSB);
o_iso = ptt.quadpolFabricLS(struct('M', Miso), z, GFREE);
f_iso = mean(isfinite(o_iso.theta0));
f_ani = mean(isfinite(o_def.theta0));
okE2 = f_iso <= 0.25 && f_ani >= 0.75;

% A replicate is given its parent's per-window verdict, so the SAME
% accepted window reports on a +-15 deg grid and on the full grid alike,
% and a window the parent rejected reports on neither. Grid width decides
% nothing either way.
par = o_def.window_ok;
NAR = GFREE; NAR.theta_grid = TH + (-15:3:15) * pi/180; NAR.window_ok = par;
FUL = GFREE; FUL.window_ok = par;
r_nar = ptt.quadpolFabricLS(struct('M', Mg), z, NAR);
r_ful = ptt.quadpolFabricLS(struct('M', Mg), z, FUL);
acc = par(:);
okE3 = nnz(acc) >= 3 && all(isfinite(r_nar.theta0(acc))) && ...
  all(isfinite(r_ful.theta0(acc))) && ...
  all(isnan(r_nar.theta0(~acc))) && all(isnan(r_ful.theta0(~acc)));

% ... and the parent's rejections are still real rejections: the flat
% column must not be handed an axis by the narrow-grid pass either.
ISO = GFREE; ISO.theta_grid = TH + (-15:3:15) * pi/180;
ISO.window_ok = o_iso.window_ok;
r_iso = ptt.quadpolFabricLS(struct('M', Miso), z, ISO);
okE4 = mean(isfinite(r_iso.theta0)) <= 0.25 && ...
  isequal(isfinite(r_iso.theta0), logical(o_iso.window_ok(:)));

% the standard error the handoff exists to make meaningful: finite, and
% inside the 52 deg circular-uniform SD on a column whose axis is well
% determined. theta is modulo 180, so nothing may exceed 180 either.
UNIF = rad2deg(pi / sqrt(12));
se_all = rad2deg([se_th(:); Jf.se_theta(:)]);
se_all = se_all(isfinite(se_all));
okE5 = ~isempty(se_all) && all(se_all < 180) && median(se_all) < UNIF;

fprintf('\nE. PARENT VERDICT INHERITED (grid width gates nothing)\n');
fprintf(['   explicit full 0-180 grid vs no grid: abstention identical %d, ' ...
  'max |dtheta| %.1e deg\n'], same_abst, rad2deg(H_or(dth, 0)));
fprintf(['   parent accepts %d of %d windows; on its grid %d report, on a ' ...
  '+-15 deg grid %d\n'], nnz(acc), numel(acc), ...
  nnz(isfinite(r_ful.theta0)), nnz(isfinite(r_nar.theta0)));
fprintf(['   fabric column reports %.0f%% of windows, isotropic column ' ...
  '%.0f%% (narrow-grid replicate %.0f%%)\n'], ...
  100*f_ani, 100*f_iso, 100*mean(isfinite(r_iso.theta0)));
fprintf('   se_theta: median %.1f deg, max %.1f deg (uniform SD %.1f)\n', ...
  median(se_all), max(se_all), UNIF);
fprintf('E1. a full grid supplied explicitly is not "restricted": %s\n', H_tick(okE1));
fprintf('E2. the full fit abstains on a flat column:            %s\n', H_tick(okE2));
fprintf('E3. accepted windows report on narrow AND full grids:  %s\n', H_tick(okE3));
fprintf('E4. rejected windows report on neither:                %s\n', H_tick(okE4));
fprintf('E5. se_theta finite, under uniform SD, never > 180:    %s\n', H_tick(okE5));
fails = fails + ~okE1 + ~okE2 + ~okE3 + ~okE4 + ~okE5;

%% F. the axis standard error is bounded and abstains when uninformative
UNIFSD = rad2deg(pi / sqrt(12));      % 51.96 deg, a uniform spread of axes
CEIL = rad2deg(sqrt(2)/2);            % 40.51 deg, the statistic's own bound
rng(101);
worst = 0; okF1 = true; okF2 = true;
for m = 3:21
  % uniform over the axis range, two clusters a QUARTER turn apart, and a
  % run of random uniform draws - the sqrt(m-1) growth would show up here
  % if it returned. 0 and 180 deg are the SAME axis, so a 0/180 split is
  % zero scatter, not maximal; the genuinely opposed pair is 0 and 90 deg,
  % which are antipodal once doubled.
  cases = {linspace(0, pi, m+1), ...
    [zeros(1, ceil(m/2)), (pi/2)*ones(1, floor(m/2))]};
  cases{1} = cases{1}(1:m);
  for t = 1:8, cases{end+1} = pi * rand(1, m); end %#ok<SAGROW>
  for c = 1:numel(cases)
    [sv, mv, rv, satv] = ptt.circAxisSE(cases{c}, 3);
    okF1 = okF1 && mv == m && (isnan(sv) || (sv >= 0 && rad2deg(sv) <= CEIL + 1e-9));
    % every NaN must say which kind it is, and a finite value must not
    okF1 = okF1 && (isnan(sv) == (satv || mv < 3)) && (isnan(sv) || ~satv);
    if isfinite(sv), worst = max(worst, rad2deg(sv)); end
    if isnan(sv) && satv, okF1 = okF1 && isfinite(rv); end
  end
end
% the quarter-turn split is maximal dispersion for an axis: it must abstain
[sq, ~, rq, satq] = ptt.circAxisSE([zeros(1,6), (pi/2)*ones(1,5)], 3);
okF1 = okF1 && isnan(sq) && satq && rq < 0.2;
% a genuinely uniform spread must ABSTAIN, not return a large number
for m = 8:21
  [su, ~, ~, satu] = ptt.circAxisSE(linspace(0, pi - pi/m, m), 3);
  okF2 = okF2 && isnan(su) && satu;
end
% ... and so must the SATURATED band, where the bound has flattened the
% statistic: measured, m = 22 reported ~36.7 deg for every true scatter
% from 16 to 40 deg, so a value there is not a measurement. Nothing
% reported may come from inside it, and the plateau must be abstained on.
% the resolution floor caps what can be reported at 31.2 deg (30.1 at
% m = 5, where sampling binds instead), so nothing may come back above it
RESCAP = 31.3;
% NDRAW and PLAT_MAX are set from the estimator's own sampling spread, not
% from one cell's mean. MEASURED over 60 independent runs: the worst of the
% six cells (m = 11, 16 deg, which sits on the abstention knee by
% construction) averages 0.304 with sd 0.031 at 200 draws - range
% 0.250-0.435, so a 0.30 bar is decided by the stream and not by the code -
% and sd 0.012 at 2000 draws, range 0.278-0.327. A 0.45 bar is ~12 sd above
% that mean while still far under the ~100% reporting the flat plateau this
% guards against produced, so the verdict turns on behaviour, not on RNG.
NDRAW = 2000; PLAT_MAX = 0.45;
plateau = zeros(1, 0); res_max = 0;
for m = [11 22]
  for sc = [16 24 40]
    fin = 0;
    for t = 1:NDRAW
      sv = ptt.circAxisSE(TH + deg2rad(sc) * randn(1, m), 3);
      if isfinite(sv), fin = fin + 1; res_max = max(res_max, rad2deg(sv)); end
    end
    plateau(end+1) = fin / NDRAW; %#ok<SAGROW>
  end
end
% and the band must still DISCRIMINATE where it does report: at m = 22 the
% unbounded-cut version read 36.7 deg for 16, 24 and 40 deg alike
sc_lo = []; sc_hi = [];
for t = 1:400
  a = ptt.circAxisSE(TH + deg2rad(6) * randn(1, 22), 3);
  b = ptt.circAxisSE(TH + deg2rad(12) * randn(1, 22), 3);
  if isfinite(a), sc_lo(end+1) = rad2deg(a); end %#ok<SAGROW>
  if isfinite(b), sc_hi(end+1) = rad2deg(b); end %#ok<SAGROW>
end
disc = median(sc_hi) - median(sc_lo);
okF2b = max(plateau) <= PLAT_MAX && res_max <= RESCAP && ...
  numel(sc_lo) >= 50 && numel(sc_hi) >= 20 && disc >= 3;
% tightly clustered replicates keep the delete-one inflation: the bounded
% statistic must still match sqrt(m-1)*rms(deviation), not undercut it
mt = 11; sd_t = deg2rad(1.5);
th_t = TH + sd_t * randn(1, mt);
st = ptt.circAxisSE(th_t, 3);
dt = angle(exp(2i*(th_t - angle(sum(exp(2i*th_t)))/2)))/2;
lin_t = sqrt((mt - 1) / mt * sum(dt.^2));
okF3 = isfinite(st) && abs(st - lin_t) / lin_t < 0.02;
% the three states a caller must be able to tell apart: too few
% replicates, replicates that scattered past the resolvable range, and a
% genuine measurement. They must not collapse into one bare NaN.
[sf, mf, rf, satf] = ptt.circAxisSE([TH, TH + 0.01], 3);
[sc2, mc2, rc2, satc2] = ptt.circAxisSE(TH + deg2rad(1.5)*randn(1, 11), 3);
[ss2, ms2, rs2, sats2] = ptt.circAxisSE(pi * rand(1, 11), 3);
okF4 = isnan(sf) && mf == 2 && ~satf && isnan(rf) ...
  && isfinite(sc2) && mc2 == 11 && ~satc2 && isfinite(rc2) ...
  && isnan(ss2) && ms2 == 11 && sats2 && isfinite(rs2);

fprintf('\nF. AXIS SE IS BOUNDED (theta is modulo 180 deg)\n');
fprintf(['   m = 3..21, uniform / two-ended / random draws: worst finite ' ...
  'SE %.1f deg (bound %.1f, uniform SD %.1f)\n'], worst, CEIL, UNIFSD);
fprintf('   clustered %.3f deg vs linear jackknife %.3f deg (%d replicates)\n', ...
  rad2deg(st), rad2deg(lin_t), mt);
fprintf(['   saturated band (m 11/22, scatter 16-40 deg): at most %.0f%% report ' ...
  '(bar %.0f%%), worst %.1f deg (cap %.1f)\n'], ...
  100*max(plateau), 100*PLAT_MAX, res_max, RESCAP);
fprintf('   still discriminating at m = 22: 6 deg -> %.1f, 12 deg -> %.1f\n', ...
  median(sc_lo), median(sc_hi));
fprintf('F1. no input exceeds the circular bound:               %s\n', H_tick(okF1));
fprintf('F2. a uniform spread abstains rather than reporting:   %s\n', H_tick(okF2));
fprintf('F2b. the saturated band abstains, it does not plateau: %s\n', H_tick(okF2b));
fprintf('F3. clustered replicates keep the delete-one scale:    %s\n', H_tick(okF3));
fprintf('F4. few / scattered / resolved are distinguishable:    %s\n', H_tick(okF4));
fails = fails + ~okF1 + ~okF2 + ~okF2b + ~okF3 + ~okF4;

%% G. the default free-mode, full-grid path must be untouched by all of this
% Everything from the span-correction saga onward was scoped to narrow-grid
% resampling and to the axis SE. The per-window theta0 / dlam profile the
% depth movie draws is produced with no theta_grid, no window_ok and no
% jackknife, so setting none of the new options must give the same answer
% as the pipeline's own default call - bit for bit, not just close.
GDEF = OPTS; GDEF.theta_const = false;
g_ref = ptt.quadpolFabricLS(struct('M', Mg), z, GDEF);
GEXTRA = GDEF; GEXTRA.window_ok = [];        % documented "decide as usual"
g_new = ptt.quadpolFabricLS(struct('M', Mg), z, GEXTRA);
same_th = isequaln(g_ref.theta0, g_new.theta0);
same_dl = isequaln(g_ref.dlam, g_new.dlam);
same_q = isequaln(g_ref.q_theta, g_new.q_theta);
same_ped = isequaln(g_ref.pedestal, g_new.pedestal);
% q_theta must be the RAW contrast again: no span factor divides it
raw_q = (max(g_ref.theta_cost, [], 2) - min(g_ref.theta_cost, [], 2)) ./ g_ref.theta_c2;
okq = isfinite(raw_q) & isfinite(g_ref.q_theta);
d_raw = H_or(max(abs(raw_q(okq) - g_ref.q_theta(okq))), inf);
okG = same_th && same_dl && same_q && same_ped && nnz(okq) >= 3 && d_raw < 1e-9;
fprintf('\nG. DEFAULT FULL-GRID PATH UNCHANGED\n');
fprintf(['   theta0 %d, dlam %d, q_theta %d, pedestal %d identical; ' ...
  'q_theta vs raw contrast %.1e\n'], same_th, same_dl, same_q, same_ped, d_raw);
fprintf('G1. the free-mode full-grid fit is bit-identical:      %s\n', H_tick(okG));
fails = fails + ~okG;

%% B. block split-half, axis held at truth
RB = 24; NBLK = 124;
BOPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.30, 'deramped', false, ...
  'theta0', TH, 'pedestal', 'window');
dl_b = []; hd_b = [];
for r = 1:RB
  S = H_col(z, TH, DL_Z, NBLK, gpd, LEAK_C, LEAK_D, NA);
  ob = ptt.quadpolFabricLS(S, z, BOPTS);
  jm = NBLK / 2;
  S1 = struct('hh', S.hh(:, 1:jm), 'vv', S.vv(:, 1:jm), 'hv', S.hv(:, 1:jm), 'vh', S.vh(:, 1:jm));
  S2 = struct('hh', S.hh(:, jm+1:end), 'vv', S.vv(:, jm+1:end), 'hv', S.hv(:, jm+1:end), 'vh', S.vh(:, jm+1:end));
  o1 = ptt.quadpolFabricLS(S1, z, BOPTS);
  o2 = ptt.quadpolFabricLS(S2, z, BOPTS);
  if isempty(dl_b), dl_b = nan(numel(ob.zw), RB); hd_b = dl_b; zb = ob.zw; end
  dl_b(:, r) = ob.dlam; hd_b(:, r) = o1.dlam - o2.dlam;
end
dl_true_b = interp1(z, DL_Z, zb, 'linear', 'extrap');
sd_b = std(dl_b, 0, 2, 'omitnan');
se_b = 0.5 * sqrt(mean(hd_b.^2, 2, 'omitnan'));
okb = isfinite(sd_b) & isfinite(se_b) & sum(isfinite(dl_b), 2) >= RB - 4 & sum(isfinite(hd_b), 2) >= RB - 8;
r_b = median(se_b(okb) ./ sd_b(okb));
fin_h = mean(isfinite(hd_b(:)));
fprintf('\nB. BLOCK SPLIT-HALF (R=%d blocks of %d traces, axis held)\n', RB, NBLK);
fprintf('   dlam: empirical SD median %.4f, split-half sigma median %.4f -> ratio %.2f (%d windows; halves finite %.0f%%)\n', ...
  median(sd_b(okb)), median(se_b(okb)), r_b, nnz(okb), 100*fin_h);
okB = r_b >= 0.7 && r_b <= 1.4;
fprintf('B1. split-half sigma within [0.7, 1.4] of empirical:  %s\n', H_tick(okB));
fails = fails + ~okB;

fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_quadpol_uncertainty:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function v = H_optd(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function v = H_or(x, d)
if isempty(x), v = d; else, v = x; end
end

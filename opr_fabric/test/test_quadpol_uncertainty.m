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
%      returns NaN from H_circ_se with nothing saying why. See the note in
%      ptt.quadpolJackknife.
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

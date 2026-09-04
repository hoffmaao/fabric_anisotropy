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

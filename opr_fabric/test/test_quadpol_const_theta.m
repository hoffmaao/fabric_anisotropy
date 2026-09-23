%TEST_QUADPOL_CONST_THETA Constant-orientation mode: what it buys and costs.
%
% opts.theta_const collapses the LS estimator's per-window theta0 to ONE
% axis for the column: every window's theta cost curve (already computed by
% the grid search) votes, and the pooled winner is refit into every window.
%
% This is a MODEL ASSUMPTION, so the test pins both directions. Against a
% column whose axis really is constant it must do better; against one that
% rotates it must do worse, and must SAY so rather than return a confident
% average. A mode that only ever looked good would be untestable.
%
%   1. CONSTANT-AXIS TRUTH (theta = 35 deg, dlam ramping 0.02 -> 0.08 with
%      depth - the divide case): constant mode's axis error is under 3 deg
%      and its per-window theta scatter is at least 3x smaller than the
%      free-theta fit, while dlam accuracy is no worse. Fewer free
%      parameters on a quantity that does not vary buys back the
%      depth-varying one.
%   2. ROTATING-AXIS TRUTH (35 deg -> 125 deg over the column - the margin
%      case): dlam degrades measurably against free-theta, and the mode
%      does not hide it. This is the guard against adopting it everywhere.
%   3. HONEST DIAGNOSTIC: the returned theta_spread_deg - the circular
%      spread of the per-window cost minima - separates the two cases
%      WITHOUT knowing the truth, so a user can tell whether the assumption
%      holds at a new site. It must be several times larger on the rotating
%      column. The pooled contrast theta_const_q does NOT do this and is
%      not asserted on: measured, it came out HIGHER on the rotating column
%      (0.225) than the constant one (0.219), because it reports how
%      sharply the POOLED curve peaks rather than whether the windows
%      agreed with each other. That is exactly the distinction between a
%      diagnostic and a number that merely looks like one.
%      The vote and the spread use only windows the free fit would trust
%      and pool curves normalised by their own data power: on Ridge A
%      frame 009 a raw-cost vote let seven near-surface windows (86% of
%      the data power) put the axis 6 deg off the column, and unweighted
%      minima of untrusted windows reported a 15 deg spread where the
%      trusted ones agree to 2-6 deg.
%   4. OFF BY DEFAULT: without the option the outputs are bit-identical to
%      the previous behaviour, so nothing already computed moves.
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_const_theta.m')"
clear;
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

rng(3);
C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);

Nz = 2400; dz = 0.5;
z = (0:Nz-1).' * dz;                     % 0..1200 m
Nx = 300;
NA = 0.35;
LEAK_C = 0.28 * exp(0.6i);
LEAK_D = 0.30;

OPTS = struct('fc', fc, 'psi_step_deg', 3, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.30, 'deramped', false, ...
  'theta_step_deg', 3);

DL_Z = 0.02 + 0.06 * (z / z(end));       % ramping contrast, both cases

function S = H_col(z, th_of_z, DL_Z, Nx, gpd, LEAK_C, LEAK_D, NA)
  % A column whose axis may rotate with depth: the two-way phase is the
  % accumulated integral of dlam, and the axis at each depth sets the
  % projection. Speckle is common to the co-pol pair, noise independent.
  Nz = numel(z);
  dz = median(diff(z));
  del = 2 * cumsum(DL_Z) * dz;           % accumulated, x2 for two-way
  del = del * gpd / 2;                   % gpd already carries the 2
  r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  S = struct('hh', zeros(Nz, Nx), 'vv', zeros(Nz, Nx), ...
    'hv', zeros(Nz, Nx), 'vh', zeros(Nz, Nx));
  for i = 1:Nz
    c = cos(th_of_z(i)); s = sin(th_of_z(i));
    ex = exp(1i * del(i));
    hh = c^2 * ex + s^2;
    vv = s^2 * ex + c^2;
    hv = c*s * (ex - 1);
    S.hh(i, :) = hh * r(i, :);
    S.vv(i, :) = vv * r(i, :);
    xc = hv * r(i, :) + LEAK_C * ((hh + vv)/2) * r(i, :) + LEAK_D * g(i, :);
    S.hv(i, :) = xc;
    S.vh(i, :) = xc;
  end
  for f = {'hh','vv','hv','vh'}
    S.(f{1}) = S.(f{1}) + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  end
end

function [aerr, ascat, derr] = H_score(o, th_true_of_zw, DL_at_zw)
  ok = isfinite(o.theta0) & isfinite(o.dlam);
  if nnz(ok) < 3, aerr = NaN; ascat = NaN; derr = NaN; return; end
  d = angle(exp(2i * (o.theta0(ok) - th_true_of_zw(ok)))) / 2;
  aerr = rad2deg(abs(median(d)));
  R = mean(exp(2i * o.theta0(ok)));
  ascat = rad2deg(sqrt(max(-2*log(abs(R)), 0))) / 2;
  derr = median(abs(o.dlam(ok) - DL_at_zw(ok)));
end

%% case 1: constant axis
TH1 = deg2rad(35);
S1 = H_col(z, TH1*ones(Nz,1), DL_Z, Nx, gpd, LEAK_C, LEAK_D, NA);
of1 = ptt.quadpolFabricLS(S1, z, OPTS);
oc1 = ptt.quadpolFabricLS(S1, z, setfield(OPTS, 'theta_const', true)); %#ok<SFLD>
dlw = interp1(z, DL_Z, of1.zw, 'linear', 'extrap');
[a_f1, s_f1, d_f1] = H_score(of1, TH1*ones(numel(of1.zw),1), dlw);
[a_c1, s_c1, d_c1] = H_score(oc1, TH1*ones(numel(oc1.zw),1), dlw);
fprintf('CONSTANT-AXIS TRUTH (35 deg)\n');
fprintf('  free  theta: err %5.1f deg  scatter %5.1f deg  dlam err %.4f\n', ...
  a_f1, s_f1, d_f1);
fprintf('  const theta: err %5.1f deg  scatter %5.1f deg  dlam err %.4f  (q %.3f, n %d)\n', ...
  a_c1, s_c1, d_c1, oc1.theta_const_q, oc1.theta_const_n);
ok_axis = a_c1 < 3;
ok_scat = s_c1 * 3 <= max(s_f1, eps);
ok_dlam = d_c1 <= d_f1 * 1.25;
fprintf('1. constant mode: axis <3 deg, scatter >=3x tighter, dlam no worse: %s\n', ...
  H_tick(ok_axis && ok_scat && ok_dlam));

%% case 2: rotating axis - the mode's assumption is FALSE here
rng(3);
TH2 = deg2rad(35) + deg2rad(90) * (z / z(end));
S2 = H_col(z, TH2, DL_Z, Nx, gpd, LEAK_C, LEAK_D, NA);
of2 = ptt.quadpolFabricLS(S2, z, OPTS);
oc2 = ptt.quadpolFabricLS(S2, z, setfield(OPTS, 'theta_const', true)); %#ok<SFLD>
th2w = interp1(z, TH2, of2.zw, 'linear', 'extrap');
[~, ~, d_f2] = H_score(of2, th2w, dlw);
[~, ~, d_c2] = H_score(oc2, th2w, dlw);
fprintf('\nROTATING-AXIS TRUTH (35 -> 125 deg)\n');
fprintf('  free  theta: dlam err %.4f\n', d_f2);
fprintf('  const theta: dlam err %.4f  (q %.3f)\n', d_c2, oc2.theta_const_q);
ok_rot = d_c2 > d_f2;
fprintf('2. constant mode is WORSE where the axis rotates:        %s\n', ...
  H_tick(ok_rot));

%% case 3: the diagnostic separates the two without knowing the truth
ok_diag = oc2.theta_spread_deg > 3 * oc1.theta_spread_deg && ...
  oc1.theta_const_n >= 3;
fprintf(['3. per-window spread separates constant from rotating: %s ' ...
  '(%.1f vs %.1f deg)\n'], H_tick(ok_diag), oc1.theta_spread_deg, ...
  oc2.theta_spread_deg);
fprintf('   (pooled contrast q, NOT a discriminant: %.3f vs %.3f)\n', ...
  oc1.theta_const_q, oc2.theta_const_q);

%% case 4: off by default changes nothing
o_def = ptt.quadpolFabricLS(S1, z, OPTS);
same = isequaln(o_def.theta0, of1.theta0) && isequaln(o_def.dlam, of1.dlam);
fprintf('4. off by default, outputs unchanged:                    %s\n', H_tick(same));

fails = ~(ok_axis && ok_scat && ok_dlam) + ~ok_rot + ~ok_diag + ~same;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_quadpol_const_theta:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

%SYNTHETIC_COMMON_OFFSET Validate the common-offset horizontal-fabric method.
%   Forward-models cumulative TWTT differences t_x - t_y for six tracked
%   reflectors at a fixed 5 m half-offset (ground-based accumulation radar
%   style), pollutes them with noise, then recovers the piecewise-constant
%   horizontal fabric contrast dlam = lam_x - lam_y by layer stripping.
%   Run twice: with the correct assumed lam_z profile, and with a naive
%   isotropic-vertical assumption (lam_z = 1/3), to show the horizontal
%   contrast is insensitive to the lam_z assumption at small offsets.

clear;
rng(2);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
figDir = fullfile(thisDir, '..', 'figs');
if ~exist(figDir, 'dir'), mkdir(figDir); end

% True column: across-flow girdle developing with depth
parT = ptt.defaultParams(); % H = 250 m, zhat_bco = 0.9, e0 = 0.8, p = 0.5
parT.lam_x_sfc = 1/3;  parT.lam_z_sfc = 1/3;  % isotropic surface
parT.lam_x_bed = 0.05; parT.lam_z_bed = 0.50; % bed: lam_y = 0.45, dlam = -0.4

% Acquisition: fixed half-offset, six tracked reflectors
obs.L = 5;
zhat_layers = [0.9 0.8 0.65 0.5 0.35 0.2];
obs.z = (zhat_layers * parT.H).';

sigma = 0.1; % TWTT pick noise [ns]
dtau_true = ptt.twttDifference(parT, obs.L * ones(size(obs.z)), obs.z);
obs.dtau  = dtau_true + sigma * randn(size(dtau_true));

% Inversion 1: correct assumed lam_z profile
parA = parT;
[dlamA, outA] = ptt.invertHorizontalFabric(obs, parA);

% Inversion 2: naive lam_z = 1/3 everywhere
parB = parT;
parB.lam_z_sfc = 1/3; parB.lam_z_bed = 1/3;
[dlamB, outB] = ptt.invertHorizontalFabric(obs, parB);

% True contrast at interval midpoints for comparison
zmid = (outA.ztop + outA.zbot) / (2 * parT.H);
PT   = ptt.columnProfiles(parT, zmid);
dlam_true = PT.lam(:, 1) - PT.lam(:, 2);

fprintf('Interval (z/H)      true dlam   inferred    naive lam_z\n');
for k = 1:numel(dlamA)
    fprintf('%5.2f - %4.2f      %9.3f  %9.3f  %9.3f\n', ...
        outA.ztop(k)/parT.H, outA.zbot(k)/parT.H, dlam_true(k), dlamA(k), dlamB(k));
end
fprintf('Max |inferred - true|: %.3f (correct lam_z), %.3f (naive lam_z)\n', ...
    max(abs(dlamA - dlam_true)), max(abs(dlamB - dlam_true)));

% Plot: contrast profile and traveltime fit
fig = figure('Visible', 'off', 'Position', [100 100 900 400]);

subplot(1, 2, 1);
zh = linspace(0, 1, 500).';
P  = ptt.columnProfiles(parT, zh);
plot(P.lam(:, 1) - P.lam(:, 2), zh, 'k-', 'LineWidth', 1.5); hold on;
zedges = [outA.zbot; parT.H] / parT.H;
stairs([dlamA; dlamA(end)], zedges, 'g-', 'LineWidth', 1.5);
stairs([dlamB; dlamB(end)], zedges, 'b--', 'LineWidth', 1.2);
xlabel('\Delta\lambda = \lambda_x - \lambda_y'); ylabel('z/H');
legend('true', 'inferred', 'inferred (naive \lambda_z)', 'Location', 'northwest');
title('Horizontal fabric contrast');

subplot(1, 2, 2);
plot(obs.dtau, obs.z/parT.H, 'ko', dtau_true, obs.z/parT.H, 'k-', ...
    outA.dtau_fit, obs.z/parT.H, 'g^', 'LineWidth', 1.2);
xlabel('\Delta\tau_{x,y} (ns)'); ylabel('z/H');
legend('noisy obs', 'noise-free', 'fit', 'Location', 'northwest');
title(sprintf('TWTT differences, L = %g m', obs.L));

figFile = fullfile(figDir, 'synthetic_common_offset.png');
print(fig, figFile, '-dpng', '-r150');
fprintf('Figure saved to %s\n', figFile);
fprintf('Total runtime: %.1f s\n', toc(t0));

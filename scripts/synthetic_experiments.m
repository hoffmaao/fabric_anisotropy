%SYNTHETIC_EXPERIMENTS Validate the model against Rathmann (2026) sec. 4a.
%   Reproduces the three synthetic inversion experiments: forward-model TWTT
%   differences for a known fabric, pollute them with 0.25 ns Gaussian
%   noise, then recover the parameters from an uninformed initial guess
%   (isotropic column, BCO at half thickness). CMP sampling follows the
%   recommended minimal stencil: half-offsets L/H = {0, 0.25, 0.5} and four
%   well-spaced reflector levels.

clear;
rng(1);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..')); % make +ptt visible
figDir = fullfile(thisDir, '..', 'figs');
if ~exist(figDir, 'dir'), mkdir(figDir); end

par = ptt.defaultParams(); % H = 250 m, e0 = 0.8, p = 0.5
par.zhat_bco = 0.9;

% True fabric endpoints (surface -> bed), paper's experiments 1-3
exps = struct( ...
    'name', {'Vertical single maximum', 'Vertical girdle', 'Mixed'}, ...
    'sfc',  {[1/3 1/3], [1/3 1/3], [0.2 0.5]}, ... % [lam_x lam_z] at surface
    'bed',  {[0 1],     [0 1/2],   [0 0.75]});     % [lam_x lam_z] at bed

% Minimal CMP stencil
Lhat = [0 0.25 0.5];
zhat = [0.1 0.3 0.5 0.7];
[LL, ZZ] = ndgrid(Lhat * par.H, zhat * par.H);

sigma = 0.25; % TWTT noise [ns]
free  = {'zhat_bco', 'lam_x_sfc', 'lam_x_bed', 'lam_z_sfc', 'lam_z_bed'};

fig = figure('Visible', 'off', 'Position', [100 100 1200 900]);

for iexp = 1:numel(exps)
    parT = par;
    parT.lam_x_sfc = exps(iexp).sfc(1); parT.lam_z_sfc = exps(iexp).sfc(2);
    parT.lam_x_bed = exps(iexp).bed(1); parT.lam_z_bed = exps(iexp).bed(2);

    % Synthetic "observed" data (drop shadow-zone and >70 deg geometries)
    [dtau, theta0] = ptt.twttDifference(parT, LL(:), ZZ(:));
    keep = ~isnan(dtau) & theta0 <= 70;
    obs.L    = LL(keep);
    obs.z    = ZZ(keep);
    obs.dtau = dtau(keep) + sigma * randn(nnz(keep), 1);

    % Uninformed initial guess: isotropic column, BCO at half thickness
    par0 = par;
    par0.zhat_bco = 0.5;
    [par0.lam_x_sfc, par0.lam_x_bed, par0.lam_z_sfc, par0.lam_z_bed] = deal(1/3);

    [parF, res] = ptt.invertFabric(obs, par0);

    fprintf('\n=== Experiment %d: %s ===\n', iexp, exps(iexp).name);
    fprintf('  %d CMP points, misfit J: %.3g -> %.3g ns^2 (noise^2 = %.3g)\n', ...
        numel(obs.L), res.J0, res.J, sigma^2);
    fprintf('  %-10s %8s %8s\n', 'param', 'true', 'inferred');
    for k = 1:numel(free)
        fprintf('  %-10s %8.3f %8.3f\n', free{k}, parT.(free{k}), parF.(free{k}));
    end

    % Profiles: true (solid) vs inferred (dashed)
    zh = linspace(0, 1, 500).';
    PT = ptt.columnProfiles(parT, zh);
    PF = ptt.columnProfiles(parF, zh);

    subplot(numel(exps), 3, (iexp-1)*3 + 1);
    plot(PT.rho, zh, 'k-', PF.rho, zh, 'g--', 'LineWidth', 1.5);
    xlabel('relative density'); ylabel('z/H');
    title(sprintf('Exp %d: %s', iexp, exps(iexp).name));
    legend('true', 'inferred', 'Location', 'southwest');

    subplot(numel(exps), 3, (iexp-1)*3 + 2);
    plot(PT.e, zh, 'k-', PF.e, zh, 'g--', 'LineWidth', 1.5);
    xlabel('e'); xlim([0 1]);

    subplot(numel(exps), 3, (iexp-1)*3 + 3);
    h = plot(PT.lam, zh, '-', PF.lam, zh, '--', 'LineWidth', 1.5);
    set(h(1:3), {'Color'}, {'r'; 'b'; 'k'});
    set(h(4:6), {'Color'}, {'r'; 'b'; 'k'});
    xlabel('\lambda_i'); xlim([0 1]);
    legend(h(1:3), {'\lambda_x', '\lambda_y', '\lambda_z'}, 'Location', 'east');
end

figFile = fullfile(figDir, 'synthetic_experiments.png');
print(fig, figFile, '-dpng', '-r150');
fprintf('\nFigure saved to %s\n', figFile);
fprintf('Total runtime: %.1f s\n', toc(t0));

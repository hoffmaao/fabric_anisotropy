%MECH_EXPERIMENT_009 Split-sample tests of the depth-wiggle mechanism.
%
% The windowed-LS dlam profile carries 50-150 m wiggles (rms ~0.008) that
% are laterally coherent (adjacent-block correlation 0.75 of the
% DETRENDED profiles) and shared with the Ershadi chain (r = 0.68), so
% they live in the coherence field, not in either estimator's internals.
% Candidate mechanisms, and the split that discriminates each:
%
%   azimuth halves  interleaved psi grids (0,4,8,... vs 2,6,10,...): the
%                   same depths through different azimuth shapes. Wiggles
%                   that differ between halves are azimuth-structured
%                   (residual pedestal); identical wiggles enter through
%                   the depth/phase structure common to all azimuths.
%   trace halves    interleaved traces: two speckle-independent looks at
%                   the SAME ground. Wiggles repeating here are data-level
%                   (real fabric structure or stratigraphy-coupled bias),
%                   not sampling noise.
%   reflectivity    saved power profile, correlated offline against the
%                   detrended wiggle: the phase-centroid-wander mechanism
%                   (bright layers pulling each window's effective depth)
%                   predicts wiggle ~ d(log P)/dz smoothed at window
%                   scale.
%
% Runs from the frame-009 coreg cache; ~45 min sequential.
%   matlab -batch "maxNumCompThreads(8); mech_experiment_009"
addpath('/kucresis/scratch/hoffmana_sta/fabric/code');
cq = load(['/kucresis/scratch/hoffmana_sta/fabric/stages/quadpol/' ...
  'coreg_cache/creg_20250108_02_009.mat']);
T = struct('hh', double(cq.hh), 'vv', double(cq.vv), ...
  'hv', double(cq.hv), 'vh', double(cq.vh));
z = cq.z(:);
clear cq;

BASE = struct('fc', 750e6, 'deramped', true);
CH = {'hh','vv','hv','vh'};

runs = struct('name', {}, 'opts', {}, 'tr', {});
runs(1) = struct('name', 'full', 'opts', BASE, 'tr', 1);
oA = BASE; oA.psi_step_deg = 4; oA.psi_offset_deg = 0;
runs(2) = struct('name', 'psiA', 'opts', oA, 'tr', 1);
oB = BASE; oB.psi_step_deg = 4; oB.psi_offset_deg = 2;
runs(3) = struct('name', 'psiB', 'opts', oB, 'tr', 1);
runs(4) = struct('name', 'trA', 'opts', BASE, 'tr', 2);   % traces 1:2:end
runs(5) = struct('name', 'trB', 'opts', BASE, 'tr', 3);   % traces 2:2:end

R = struct();
for k = 1:numel(runs)
  Sk = T;
  if runs(k).tr == 2
    for c = 1:4, Sk.(CH{c}) = T.(CH{c})(:, 1:2:end); end
  elseif runs(k).tr == 3
    for c = 1:4, Sk.(CH{c}) = T.(CH{c})(:, 2:2:end); end
  end
  t0 = tic;
  o = ptt.quadpolFabricLS(Sk, z, runs(k).opts);
  fprintf('%s: %.1f min, dlam med %.3f, pedestal [%.3f %+.3f %.3f]\n', ...
    runs(k).name, toc(t0)/60, median(o.dlam, 'omitnan'), ...
    o.pedestal(1), o.pedestal(2), o.pedestal(3));
  R.(runs(k).name) = struct('zw', o.zw, 'dlam', o.dlam, ...
    'theta0', o.theta0, 'resid', o.resid, 'gamma', o.gamma, ...
    'pedestal', o.pedestal, 'q_theta', o.q_theta);
  clear Sk o;
end

% reflectivity profile for the offline stratigraphy correlation
R.pmean = mean(abs(T.hh).^2, 2);
R.z = z;
save('/kucresis/scratch/hoffmana_sta/fabric/stages/mech_009.mat', ...
  '-v7.3', 'R');
fprintf('wrote stages/mech_009.mat\n');

%GEN_INVERSE_INTUITION Numbers behind scripts/figures/inverse_intuition.py.
%
% Runs ptt.traveltimeFabricML on a synthetic column and writes the two
% .mat files that figure reads, so the panels show what the estimator
% actually does rather than a sketch of it. Run this first, then point
% the figure at the output with INTUITION_DIR.
%
% MATLAB here is the CReSIS one - the local licence has been failing with
% error 5001 - so this normally runs on mem1 and the .mat files come back.
% ptt_root defaults to the checkout this file sits in; out_dir to the cwd
if ~exist('ptt_root','var'), ptt_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); end % this checkout's root, holding +ptt
addpath(ptt_root);
rng(11);
fc = 750e6;
z = (10:10:1400).';
% a truth with structure the eye can follow: a ramp plus one localised band
dl_true = 0.02 + 0.035*(z/1400) + 0.030*exp(-((z-700)/70).^2);

C = ptt.constants(); n = sqrt(C.eps_bar);
k = 2*pi*fc*2*(C.deps/(2*n))/(C.c*1e9)/(2*pi*fc)*1e9;   % ns per m per unit dlam
dt_true = k * cumtrapz(z, dl_true);                      % the slide's Dt(z)
S_NS = 0.05;                                             % delay noise [ns]
dt_obs = dt_true + S_NS*randn(size(z));

% ---- route 1: differentiate the delay directly (what NOT to do)
dl_naive = [NaN; diff(dt_obs)./diff(z)/k];

% ---- route 2: the regularised linear ML solve, ptt.traveltimeFabricML
o = ptt.traveltimeFabricML(dt_obs, z, struct('fc', fc, 'sigma_ns', S_NS, ...
  'prior', [0.03 0.05 200], 'res_min', 0.15));

% ---- route 3: same solver, a prior that insists on much smoother ice
o_stiff = ptt.traveltimeFabricML(dt_obs, z, struct('fc', fc, 'sigma_ns', S_NS, ...
  'prior', [0.03 0.05 600], 'res_min', 0.15));

if ~exist('out_dir','var'), out_dir = pwd; end
save(fullfile(out_dir,'intuition.mat'), ...
  'z','dl_true','dt_true','dt_obs','dl_naive','S_NS','k', ...
  '-v7');
zz=o.z; dl=o.dlam; sg=o.sigma_dlam; rs=o.resolution;
dl2=o_stiff.dlam; sg2=o_stiff.sigma_dlam; rs2=o_stiff.resolution;
save(fullfile(out_dir,'intuition2.mat'), ...
  'zz','dl','sg','rs','dl2','sg2','rs2','-v7');
fprintf('naive derivative: std %.4f about truth\n', std(dl_naive(2:end)-dl_true(2:end)));
fprintf('regularised ML  : std %.4f, median sigma %.4f, %d abstained\n', ...
  std(dl(isfinite(dl))-dl_true(isfinite(dl))), median(sg,'omitnan'), nnz(isnan(dl)));
fprintf('delay span %.2f ns, noise %.2f ns -> SNR %.0f\n', ...
  dt_true(end), S_NS, dt_true(end)/S_NS);

%TEST_COPOL_SURFACE The near-surface mechanism of the co-polarized chain.
%
% The co-pol joint inversion reported dlam ~ 0.06 in its shallowest
% interval where the measured fringe rate says <= 0.015 (fringe_check,
% 57 m on Ridge A). This test reproduces the mechanism and validates the
% fix, at the solver level where it bites:
%
%   - observations are referenced to zero at a shallow depth z_ref, so
%     any reference error is a CONSTANT added to every node;
%   - the surface-anchored forward model can absorb a constant only in
%     its shallowest interval's dlam (0.3 ns over a ~74 m interval is
%     dlam ~ 0.06);
%   - with obs.zref the forward is differenced against the same
%     reference, the constant is carried by an explicit offset nuisance,
%     and the first interval is flagged reference-degenerate.
%
% Truth: isotropic cap (dlam = 0 above 200 m) ramping to 0.06 by 700 m.
% The OLD call (no obs.zref) must REPRODUCE the spurious near-surface
% fabric; the NEW call must recover the offset, honor the isotropic cap
% from the first quotable interval down, and stay unbiased at depth.
%
% Run: matlab -batch "run('opr_fabric/test/test_copol_surface.m')"
clear;
rng(7);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

par = ptt.defaultParams();
par.H = 2000;
par.lam_z_sfc = 1/3;
par.lam_z_bed = 1/3;

Z_REF = 60;                          % reference depth [m below surface]
node_depth = (134:74:1850).';        % node (interval bottom) depths [m]
m = numel(node_depth);
zs = par.H - node_depth;             % heights above bed, shallowest first
zref = par.H - Z_REF;

% truth per interval (piecewise-constant, same discretization the solver
% uses, so the test isolates the referencing mechanism rather than mixing
% in discretization error)
mid_depth = ([Z_REF; node_depth(1:end-1)] + node_depth) / 2;
dlam_true = min(max((mid_depth - 200) / 500, 0), 1) * 0.06;

% forward the truth through the SAME model family, referenced at z_ref
zhat_nodes = [zs(end:-1:1) / par.H; 1];
lam_z = ones(m+1, 1) / 3;
h = 1 - 1/3;
par_true = par;
par_true.fabric = struct('method', 'previous', 'zhat', zhat_nodes, ...
  'lam_x', ((h + [dlam_true(end:-1:1); dlam_true(1)])) / 2, 'lam_z', lam_z);
F = ptt.twttDifference(par_true, zeros(m+1, 1), [zref; zs]);
dtau_true = F(2:end) - F(1);

C_TRUE = 0.30;                       % reference error [ns]
NOISE = 0.03;                        % node noise [ns]
obs = struct('L', 0, 'z', zs, ...
  'dtau', dtau_true + C_TRUE + NOISE * randn(m, 1), ...
  'w', ones(m, 1));

fails = 0;

% --- OLD path: surface-anchored forward, no offset nuisance. Must
% REPRODUCE the bug, or the mechanism named here is not the mechanism.
[dl_old, out_old] = ptt.invertHorizontalFabricJoint(obs, par, ...
  struct('reg', 0.05));
ok = dl_old(1) > 0.03;
fails = fails + ~ok;
fprintf('old: dlam(1) %.3f (truth 0.000) - spurious near-surface fabric %s\n', ...
  dl_old(1), H_tick(ok));

% --- NEW path: reference-consistent forward + offset nuisance
obs.zref = zref;
[dl_new, out_new] = ptt.invertHorizontalFabricJoint(obs, par, ...
  struct('reg', 0.05));
ok = abs(out_new.ref_offset - C_TRUE) < 0.12;
fails = fails + ~ok;
fprintf('new: ref_offset %.3f ns (truth %.2f) %s\n', out_new.ref_offset, ...
  C_TRUE, H_tick(ok));
ok = out_new.ref_degenerate(1) && ~any(out_new.ref_degenerate(2:end));
fails = fails + ~ok;
fprintf('new: interval 1 flagged reference-degenerate %s\n', H_tick(ok));
e2 = dl_new(2:end) - dlam_true(2:end);
ok = sqrt(mean(e2.^2)) < 0.012;
fails = fails + ~ok;
fprintf('new: dlam(2:end) rms error %.4f vs truth %s\n', ...
  sqrt(mean(e2.^2)), H_tick(ok));
ok = dl_new(2) < 0.02;
fails = fails + ~ok;
fprintf('new: first quotable interval %.3f (isotropic cap honored) %s\n', ...
  dl_new(2), H_tick(ok));
ok = abs(out_new.ztop(1) - zref) < 1e-6;
fails = fails + ~ok;
fprintf('new: shallowest interval top at z_ref, not the surface %s\n', ...
  H_tick(ok));

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_copol_surface:failed', '%d check(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

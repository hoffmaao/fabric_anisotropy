%TEST_QUADPOL_SIGNED_RATE Where the axis is named, the LS rate is signed.
%
% ptt.quadpolFabricLS fixes its axis convention by keeping the rate
% ddelta >= 0: (theta0 + 90, -ddelta) is the same field, so a free search
% over the whole half turn always has the mirror to hand. Where the axis
% is NAMED - held by constant mode, pinned by the caller (the blocks and
% the split-half refits), or searched over a grid narrower than a quarter
% turn (a jackknife replicate's) - the mirror is out of reach, and the
% clamp at 0 piled weak-depth windows onto exactly 0: 17-41% of the
% finite held block cells at Taylor Dome, 29% at Eastwind, 25-30% at
% McMurdo, pushing every mean upward (issue #27).
%
% Synthetic: one column, axis 30 deg from H, dlam 0.05 down to 400 m and
% isotropic below, the record running to 900 m, with the reciprocal
% leakage and noise the other LS tests use.
%
%   1. HELD (constant mode), in the isotropic part: clamped, a large
%      fraction of the windows sit at exactly 0; signed, none do, some are
%      negative, and the clamped rates are exactly the signed ones with
%      their negative half folded onto 0 - which is the upward bias of
%      every clamped mean. In the fabric above, both read the same
%      contrast: the sign changes nothing where the rate is positive.
%   2. A WRONG-BRANCH AXIS - pinned a quarter turn off the fabric: signed,
%      it reads the contrast with its sign flipped (the named axis is the
%      slower eigen-direction); clamped, it read "no fabric".
%   3. A NARROW GRID (a replicate's +-15 deg) is named too: signed there,
%      with negatives in the isotropic part.
%   4. A FREE full-range search is untouched: signed_rate on and off give
%      bit-identical theta0 and dlam.
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_signed_rate.m')"
clear;
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
rng(27);

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);
dz = 0.5; z = (0:dz:900).';
TH = deg2rad(30); DL = 0.05; Z_ISO = 400;
S = H_col(z, TH, DL * (z < Z_ISO), 200, gpd);
M = struct('M', ptt.quadpolMoments(S, [1 200]));

BASE = struct('fc', fc, 'deramped', false, 'dlam_max', 0.25, 'psi_step_deg', 6, ...
  'theta_step_deg', 6, 'win_fit_m', 60, 'step_m', 20);
held = setfield(BASE, 'theta_const', true); %#ok<SFLD>
fails = 0;

% 1. held, isotropic part and fabric part
oc = ptt.quadpolFabricLS(M, z, setfield(held, 'signed_rate', false)); %#ok<SFLD>
os = ptt.quadpolFabricLS(M, z, held);
iso = oc.zw - oc.win_half > Z_ISO + 20 & isfinite(oc.dlam) & isfinite(os.dlam);
fab = oc.zw + oc.win_half < Z_ISO - 20 & oc.zw > 100;
z0_c = mean(oc.dlam(iso) == 0); z0_s = mean(os.dlam(iso) == 0);
m_c = mean(oc.dlam(iso)); m_s = mean(os.dlam(iso));
% the fold: every clearly negative signed window reads exactly 0 clamped,
% and elsewhere the two agree to a small fraction of the rates' spread
% (an isotropic window's cost is flat, so its optimum is only as precise
% as that flatness allows)
sp = max(abs(os.dlam(iso)));
fold = max(abs(oc.dlam(iso) - max(os.dlam(iso), 0))) / sp;
neg0 = all(oc.dlam(iso & os.dlam < -0.1 * sp) == 0);
f_c = median(oc.dlam(fab), 'omitnan'); f_s = median(os.dlam(fab), 'omitnan');
ok1 = z0_c > 0.2 && z0_s == 0 && any(os.dlam(iso) < 0) && neg0 && fold < 0.25 ...
  && abs(f_c - DL) < 0.005 && abs(f_s - f_c) < 1e-3;
fprintf(['1. held axis %.1f deg, %d isotropic windows: exactly 0 clamped %.0f%%, signed %.0f%%; ' ...
  'clamped = signed folded at 0 to %.2f of their spread; mean clamped %+.5f, signed %+.5f; ' ...
  'fabric median %.4f / %.4f (truth %.2f): %s\n'], ...
  rad2deg(os.theta_const), nnz(iso), 100 * z0_c, 100 * z0_s, fold, m_c, m_s, f_c, f_s, DL, H_tick(ok1));
fails = fails + ~ok1;

% 2. a wrong-branch axis, pinned a quarter turn off
pin = setfield(BASE, 'theta0', TH + pi/2); %#ok<SFLD>
wc = ptt.quadpolFabricLS(M, z, setfield(pin, 'signed_rate', false)); %#ok<SFLD>
ws = ptt.quadpolFabricLS(M, z, pin);
w_c = median(wc.dlam(fab), 'omitnan'); w_s = median(ws.dlam(fab), 'omitnan');
ok2 = abs(w_s + DL) < 0.005 && abs(w_c) < 0.005;
fprintf('2. axis pinned a quarter turn off: signed median %+.4f (truth %+.2f), clamped %+.4f: %s\n', ...
  w_s, -DL, w_c, H_tick(ok2));
fails = fails + ~ok2;

% 3. a narrow grid, as a replicate searches, is named
of = ptt.quadpolFabricLS(M, z, BASE);
nar = setfield(setfield(BASE, 'theta_grid', TH + deg2rad(-15:3:15)), 'window_ok', true(size(of.zw))); %#ok<SFLD>
on = ptt.quadpolFabricLS(M, z, nar);
nc = ptt.quadpolFabricLS(M, z, setfield(nar, 'signed_rate', false)); %#ok<SFLD>
isn = on.zw - on.win_half > Z_ISO + 20 & isfinite(on.dlam) & isfinite(nc.dlam);
ok3 = any(on.dlam(isn) < 0) && ~any(on.dlam(isn) == 0) && mean(nc.dlam(isn) == 0) > 0.2;
fprintf('3. +-15 deg grid: %d of %d isotropic windows negative signed, %.0f%% exactly 0 clamped: %s\n', ...
  nnz(on.dlam(isn) < 0), nnz(isn), 100 * mean(nc.dlam(isn) == 0), H_tick(ok3));
fails = fails + ~ok3;

% 4. a free full-range search is untouched
ofc = ptt.quadpolFabricLS(M, z, setfield(BASE, 'signed_rate', false)); %#ok<SFLD>
ok4 = isequaln(of.theta0, ofc.theta0) && isequaln(of.dlam, ofc.dlam) && all(of.dlam(isfinite(of.dlam)) >= 0);
fprintf('4. free search: signed_rate on and off bit-identical, every rate >= 0: %s\n', H_tick(ok4));
fails = fails + ~ok4;

fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_quadpol_signed_rate:failed', '%d check(s) failed', fails);
end

function S = H_col(z, th, dl_z, Nx, gpd)
% single-axis column (axis th, contrast dl_z along depth), co-pol speckle
% common to the pair, a reciprocal leakage with a correlated and a
% decorrelated part, and independent noise - the other LS tests' column
Nz = numel(z);
del = gpd * cumsum(dl_z) * median(diff(z));
c = cos(th); s = sin(th);
ex = exp(1i * del);
hh = c^2 * ex + s^2; vv = s^2 * ex + c^2; hv = c*s * (ex - 1);
r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
xc = hv .* r + 0.28 * exp(0.6i) * ((hh + vv)/2) .* r + 0.30 * g;
S = struct('hh', hh .* r, 'vv', vv .* r, 'hv', xc, 'vh', xc);
for k = {'hh', 'vv', 'hv', 'vh'}
  S.(k{1}) = S.(k{1}) + 0.35 * (randn(Nz, Nx) + 1i*randn(Nz, Nx));
end
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

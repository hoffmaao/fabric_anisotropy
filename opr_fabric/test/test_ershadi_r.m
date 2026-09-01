%TEST_ERSHADI_R The Sect.-3.5 reflection-ratio retrieval recovers known r.
%
% Validates ptt.fujitaModel + ptt.ershadiInverse - the faithful Ershadi et
% al. (2022) chain - on speckle synthetics with a known three-zone profile
% shaped like their EDML result (two anisotropic-scattering zones of
% opposite sign under an isotropic lid):
%
%   L1    0-400 m   dlam 0.05   theta  30 deg   r   0 dB
%   L2  400-800 m   dlam 0.20   theta  30 deg   r +10 dB
%   L3 800-1200 m   dlam 0.20   theta 120 deg   r -10 dB
%
% WHAT ANCHORS THE CONVENTIONS. The truth generator below (H_truth_M) is a
% second writing of the SAME Jones-sandwich algebra fujitaModel uses, so on
% its own it cannot settle a sign convention - a flip made in both places
% cancels and the correlation stays 1.0000. Two things do anchor them:
%   - verdict 1, which checks fujitaModel for a single uniform layer
%     against the closed-form channel expressions validated in
%     test_quadpol_ls, a derivation that never forms a matrix product; and
%   - the recovery verdicts, which reach the model only through
%     ptt.ershadiFabric's quadpolMoments sweep - the real data path, built
%     from weight vectors rather than from a rotation matrix - so a sweep
%     sense that disagrees with fujitaModel's misplaces the fitted axis.
% The H_truth_M-vs-fujitaModel correlation is kept as an internal
% consistency guard (1b) and labelled as one.
%
% Verdicts:
%   1. forward-model algebra against the INDEPENDENT closed form for one
%      uniform layer with r nonzero, on all three channels, BEFORE any
%      fitting - so a failure here is a convention bug, not an optimizer
%      one. The layer axis is 25 deg on purpose: at 0 or 45 deg the sweep
%      is its own mirror and a sweep-sense flip would be invisible.
%   1b. internal consistency: fujitaModel vs the multi-layer truth
%      generator (corr > 0.99). Shares its algebra with the model, so it
%      is a refactor guard, not a convention anchor.
%   2. r recovery: fitted r_dB medians land in the right zone windows
%      (+10 within +-3 dB, -10 within +-3 dB, lid within +-2.5 dB).
%   3. theta recovery: circular error < 5 deg per zone.
%   4. null control: with r = 0 everywhere (same seeds), |fitted r_dB| < 2
%      in every zone - speckle and noise must not manufacture anisotropy.
%   5. eq.-(13) analytic r agrees with the fit on the rows ershadiInverse
%      admits, which are the near-anti-phase ones where the relation holds
%      (median |r13 - fit| < 2 dB in BOTH strong zones). A zone that
%      resolves too few rows to sample FAILS this verdict rather than
%      skipping it: losing the nodes altogether is the regression it
%      exists to catch.
%
% Run: matlab -batch "run('opr_fabric/test/test_ershadi_r.m')"
clear;
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

rng(11);
C = ptt.constants();
fc = 750e6;
eps_perp = 3.15;                        % the paper's value, as ershadiFabric
% ONE-way rad/m per unit dlam, from the single owner the model and the
% inversion's anti-phase gate both read - writing it out here again is the
% desync this test would then be unable to see.
gpd1 = ptt.birefringentPhaseRate(fc, eps_perp, C.deps);

Nz = 2401; dz = 0.5;
z = (0:Nz-1).' * dz;                    % 0..1200 m
Nx = 360;
% Additive noise level. The co-pol NODES carry the r information, and the
% moment-averaged amplitude floor sits near sqrt(2)*NA: at 0.35 that is a
% -3 dB clip which fills every node the +-10 dB truth needs at -20 dB, and
% the first run measured exactly that failure mode (fitted +8.5/-5.4).
% 0.08 puts the floor near -16 dB, below the discretized node depths, so
% the verdicts test the ESTIMATOR rather than the clip. Noise robustness
% is a property of the data regime (Ershadi gate |C| > 0.4 for the same
% reason), not of this recovery test.
NA = 0.08;

TRUTH = struct( ...
  'top',   {0, 400, 800}, ...
  'dlam',  {0.05, 0.20, 0.20}, ...
  'theta', {deg2rad(30), deg2rad(30), deg2rad(120)}, ...
  'r_db',  {0, 10, -10});

ZONES = [200 380; 450 750; 850 1150];   % assert bands, clear of boundaries

% ---- truth generator: physical-convention Jones sandwich, layer by layer
function M = H_truth_M(z, TRUTH, gpd1)
  Nz = numel(z);
  M = complex(zeros(2, 2, Nz));
  tops = [TRUTH.top];
  NL = numel(TRUTH);
  Pacc = eye(2);
  li = 1;
  for iz = 1:Nz
    while li < NL && z(iz) >= tops(li+1)
      d = tops(li+1) - tops(li);
      Pacc = H_layer(TRUTH(li), d, gpd1) * Pacc;
      li = li + 1;
    end
    P = H_layer(TRUTH(li), z(iz) - tops(li), gpd1) * Pacc;
    ra = 10^(TRUTH(li).r_db / 20);
    Rm = [cos(TRUTH(li).theta), -sin(TRUTH(li).theta); ...
          sin(TRUTH(li).theta),  cos(TRUTH(li).theta)];
    G = Rm * diag([1, ra]) * Rm.';
    M(:, :, iz) = P.' * G * P;          % reciprocal two-way, antennas at 0
  end
end
function A = H_layer(L, d, gpd1)
  del = gpd1 * L.dlam * d;
  Rm = [cos(L.theta), -sin(L.theta); sin(L.theta), cos(L.theta)];
  A = Rm * diag([exp(+0.5i*del), exp(-0.5i*del)]) * Rm.';
end

function S = H_channels(M, Nz, Nx, NA)
  % common complex speckle per depth/trace (co-located layer reflectors,
  % the Fujita coherent-layer picture), independent additive noise
  g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  S = struct();
  S.hh = squeeze(M(1,1,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vv = squeeze(M(2,2,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.hv = squeeze(M(1,2,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vh = squeeze(M(2,1,:)) .* g + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
end

M = H_truth_M(z, TRUTH, gpd1);
S = H_channels(M, Nz, Nx, NA);

EOPT = struct('fc', fc, 'psi_step_deg', 1, 'deramped', false);
fr = ptt.ershadiFabric(S, z, EOPT);

% ---- verdict 1: fujitaModel against the INDEPENDENT closed form.
% For a SINGLE uniform layer the channels have a closed form reached by
% projecting the antennas onto the fabric eigenbasis, with no Jones product
% anywhere. With d the angle between the layer axis and the antenna h-axis
% and ex = exp(i*delta) the TWO-WAY birefringent phase factor:
%   hh = Gx cos^2(d) ex + Gy sin^2(d)
%   vv = Gx sin^2(d) ex + Gy cos^2(d)
%   hv = cos(d) sin(d) (Gy - Gx ex)
% test_quadpol_ls already validates this form at Gx = Gy = 1; the
% reflection ratio enters as the Gx/Gy weighting of the two eigen-returns,
% which is the one new thing fujitaModel adds. The forms above have the
% common one-way factor exp(-i*delta/2) divided out, so fujitaModel's
% amplitudes are compared against them times that factor.
% This DOES discriminate: the sweep sense enters through d = a - theta with
% a = -psi (the paper's R S R', ershadiFabric E1), so building the sweep at
% +psi shifts the whole field; and the layer sign enters ex, which is not
% symmetric between hh and vv.
L1 = struct('top_m', 0, 'dlam', 0.20, 'theta', deg2rad(25), 'r_db', 6);
fm1 = ptt.fujitaModel(L1, z, fr.psi, struct('fc', fc, 'eps_perp', eps_perp));
a_ant = -fr.psi;                       % antennas sit at -psi under R S R'
d1 = a_ant - L1.theta;
cd1 = cos(d1); sd1 = sin(d1);
del1 = gpd1 * L1.dlam * z;             % ONE-way; the two-way phase is 2*del1
Gx1 = 1; Gy1 = 10^(L1.r_db / 20);
ex1 = exp(2i * del1);
common = exp(-1i * del1);
cf_hh = common .* (Gx1*(cd1.^2).*ex1 + Gy1*(sd1.^2));
cf_vv = common .* (Gx1*(sd1.^2).*ex1 + Gy1*(cd1.^2));
cf_hv = common .* ((cd1.*sd1) .* (Gy1 - Gx1*ex1));
e_hh = max(abs(fm1.s_hh - cf_hh), [], 'all') / max(abs(cf_hh), [], 'all');
e_vv = max(abs(fm1.s_vv - cf_vv), [], 'all') / max(abs(cf_vv), [], 'all');
e_hv = max(abs(fm1.s_hv - cf_hv), [], 'all') / max(abs(cf_hv), [], 'all');
ok_fwd = max([e_hh, e_vv, e_hv]) < 1e-10;
fprintf(['1. fujitaModel matches the closed form (hh/vv/hv rel. err ' ...
  '%.1e %.1e %.1e):  %s\n'], e_hh, e_vv, e_hv, H_tick(ok_fwd));

% ---- verdict 1b: INTERNAL CONSISTENCY, not a convention anchor. The sweep
% below is the same R S R' algebra as fujitaModel.m and H_truth_M is the
% same depth loop, so a sign error present in both cancels here and the
% correlation is 1.0000 either way. Kept because it does catch a change
% made to one path and not the other across the full three-layer stack.
layersT = struct('top_m', {TRUTH.top}, 'dlam', {TRUTH.dlam}, ...
  'theta', {TRUTH.theta}, 'r_db', {TRUTH.r_db});
fm = ptt.fujitaModel(layersT, z, fr.psi, struct('fc', fc, ...
  'eps_perp', eps_perp));
cc = cos(-fr.psi); ss = sin(-fr.psi);
shh0 = squeeze(M(1,1,:)) * (cc.*cc) + squeeze(M(2,2,:)) * (ss.*ss) + ...
  (squeeze(M(1,2,:)) + squeeze(M(2,1,:))) * (cc.*ss);
A0 = abs(shh0);
dP0 = 20*log10(max(A0, realmin) ./ max(mean(A0, 2), realmin));
band = z >= 200 & z <= 1150;
v = corrcoef(fm.dP_hh(band, :), dP0(band, :));
ok_cons = v(1, 2) > 0.99;
fprintf(['1b. multi-layer internal consistency (shared algebra): corr ' ...
  '%.4f  %s\n'], v(1, 2), H_tick(ok_cons));

% ---- the inversion
IOPT = struct('interval_m', 100, 'z_fit', [200 1150], ...
  'w_theta', [1 0 0], 'w_r', [0 1 0], 'fc', fc, 'eps_perp', eps_perp);
inv = ptt.ershadiInverse(fr, z, IOPT);

function m = H_zone_med(v, z, zn)
  m = median(v(z >= zn(1) & z <= zn(2)), 'omitnan');
end

r_med = zeros(3, 1); th_err = zeros(3, 1);
for k = 1:3
  r_med(k) = H_zone_med(inv.r_db, z, ZONES(k, :));
  tt = H_zone_med(inv.theta, z, ZONES(k, :));
  th_err(k) = abs(angle(exp(2i*(tt - TRUTH(k).theta)))) / 2;
end
fprintf(['   staging cycles used: %d of %d, best cycle %d, ' ...
  'any cycle worsened: %d\n'], inv.n_cycles_used, 3, inv.cycle_best, ...
  inv.cycle_worsened);
fprintf(['   fitted r_db per zone: %+5.1f %+5.1f %+5.1f ' ...
  '(truth 0 / +10 / -10)\n'], r_med);
fprintf('   theta error per zone: %4.1f %4.1f %4.1f deg\n', rad2deg(th_err));
ok_r = abs(r_med(1)) < 2.5 && abs(r_med(2) - 10) < 3 && abs(r_med(3) + 10) < 3;
ok_th = all(th_err < deg2rad(5));
fprintf('2. r recovered in all three zones:                      %s\n', H_tick(ok_r));
fprintf('3. theta recovered in all three zones:                  %s\n', H_tick(ok_th));

% ---- verdict 4: null control, same seeds, r = 0 everywhere
rng(11);
TRUTH0 = TRUTH;
for k = 1:3, TRUTH0(k).r_db = 0; end
M0 = H_truth_M(z, TRUTH0, gpd1);
S0 = H_channels(M0, Nz, Nx, NA);
fr0 = ptt.ershadiFabric(S0, z, EOPT);
inv0 = ptt.ershadiInverse(fr0, z, IOPT);
r0_med = arrayfun(@(k) H_zone_med(inv0.r_db, z, ZONES(k, :)), 1:3);
fprintf('   null-control r_db per zone: %+5.1f %+5.1f %+5.1f\n', r0_med);
ok_null = all(abs(r0_med) < 2);
fprintf('4. no spurious anisotropy from speckle and noise:       %s\n', H_tick(ok_null));

% ---- verdict 5: eq.-(13) analytic agrees with the fit where it resolves.
% ershadiInverse returns r13 only on the near-anti-phase rows, the only
% place eq. (13) is valid (off anti-phase it reads +20 dB for a +10 dB
% truth at quadrature), so the finite rows ARE the fair sample. A zone that
% resolves fewer than MINROWS of them fails: a gate or node-finder
% regression that stops resolving nodes must not silence this verdict.
% BOTH strong zones must resolve MINROWS and agree, not just one. Requiring
% only one is what let a real bug through while this was being written: the
% deep zone sits at theta + 90, an unsigned phase accumulation put its
% anti-phase gate at the wrong depths, and it dropped to 16 rows at 18.1 dB
% while the shallow zone still read 0.0 dB over 109 and carried the verdict
% to PASS. The anti-phase window admits ~18% of depths by construction
% (2*acos(0.85)/(2*pi)), which still leaves 108 rows per zone.
% The tolerance is 2 dB, not the 4 dB an ungated r13 needed: off anti-phase
% eq. (13) was biased high by up to 10 dB and 4 dB was accommodating that
% bias rather than catching it. Measured 0.0 and 0.8 dB.
MINROWS = 20;
d13 = abs(inv.r13_db - inv.r_db);
m2 = z >= ZONES(2,1) & z <= ZONES(2,2) & isfinite(d13);
m3 = z >= ZONES(3,1) & z <= ZONES(3,2) & isfinite(d13);
ok_13 = true; msg = '';
for zi = {m2, m3}
  m = zi{1};
  msg = sprintf('%s %.1f dB over %d rows;', msg, median(d13(m)), nnz(m));
  ok_13 = ok_13 && nnz(m) >= MINROWS && median(d13(m)) < 2;
end
% per-zone accounting of WHY the other rows abstained. Each gate answers a
% different question - 1-3 are properties of the column above, 4 of the
% depth, 5 of the data - so a shift between them is a real change even when
% the surviving row count is unmoved.
for k = 2:3
  zm = z >= ZONES(k,1) & z <= ZONES(k,2);
  fprintf('   zone %d r13 rows by gate:', k);
  for c = 0:5
    fprintf(' %s=%d', inv.r13_reason_key{c+1}, nnz(inv.r13_reason(zm) == c));
  end
  fprintf('\n');
end
fprintf('5. eq.-(13) analytic vs fit (anti-phase rows):%s   %s\n', msg, ...
  H_tick(ok_13));

fails = ~ok_fwd + ~ok_cons + ~ok_r + ~ok_th + ~ok_null + ~ok_13;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_ershadi_r:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

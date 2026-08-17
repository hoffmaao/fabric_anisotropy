%TEST_EGRIP_BLOCKS Why 1.1 km blocks fail on EastGRIP, and what size works.
%
% The first EastGRIP validation frame (20240619_01_001) inverted end to end
% but failed internally: block-to-block dlam spread 7.2x against Ridge A's
% 1.1x, resid 0.60 against 0.21, |C_hhvv| 0.17 against 0.42. The suspect is
% the along-track block: 125 traces is ~125 m at Ridge A's ~1 m trace
% spacing but ~1.1 km at EastGRIP's ~9 m, and EastGRIP's fabric is strong
% (EGRIP core: horizontal eigenvalue difference ~0.2-0.35).
%
% THE MECHANISM BEING TESTED, which is sharper than "lateral averaging":
% pooling traces whose birefringent phase ramps differ by ddlam builds a
% phase spread  sigma_phi = gpd * ddlam * z  across the block, and the
% pooled coherence falls as exp(-sigma_phi^2/2). At a strong-fabric site a
% modest 3% relative along-track variation (ddlam ~ 0.01 within a km) costs
% 2+ radians by 800 m - the deep coherence collapse, the high residuals and
% the oscillating dlam the real frame shows. The SAME absolute variation
% inside a 127 m block is ~10x smaller and harmless. Block length must
% scale with 1/(dlam variation), which in practice means: strong-fabric
% sites need Ridge A LENGTH, not Ridge A trace count.
%
% Synthetic: EastGRIP-like column (dlam ramping 0.25 -> 0.35 along 4.5 km,
% axis fixed, 9 m trace spacing so traces are independent speckle), the
% test_quadpol_ls pedestal (correlated + decorrelated parts) and noise at
% the level that gives a realistic coherence floor. The full pipeline
% two-pass is mimicked exactly: frame-level theta0 + pedestal, then
% per-block dlam with theta0 handed in, at NBLK = 125, 28, 14.
%
% WHAT THIS GUARDS. Running it showed something sharper than the block
% question it was written for: a 0.10 along-frame ramp kills the FRAME
% pass itself (pooled coherence dies by ~100 m depth; theta0 lands ~90 deg
% off), and because the two-pass design hands the frame theta0 to every
% block, NO block size rescues the section. That single point of failure
% is the property asserted here, so a future change to the handoff that
% silently re-couples blocks to a dead frame fit fails this test. (The
% dlam CAP - the actual EastGRIP culprit - is guarded separately by
% test_egrip_cap.m; the cap is raised here so it cannot mask the ramp.)
%
% PASS =  125-trace blocks fail on the ramp AND 14-trace blocks fail the
%         same way while the frame pass is dead.
%
% Run: matlab -batch "run('opr_fabric/test/test_egrip_blocks.m')"
clear;
rng(7);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);   % rad/m per unit dlam

THETA_TRUE = deg2rad(35);      % axis, fixed along track (antennas at 0)
Nz = 2800; dz = 0.5;
z = (0:Nz-1).' * dz;           % 0..1400 m
Nx = 504; dx = 9.0;            % 4.5 km of independent 9 m traces
x = (0:Nx-1) * dx;
DL_X = 0.25 + 0.10 * (x / x(end));   % the along-track ramp

LEAK_C = 0.30 * exp(0.7i);
LEAK_D = 0.35;
NA = 0.5;

% dlam_max raised so the cap (test_egrip_cap.m's subject) cannot mask the
% ramp mechanism this test exists to pin down
OPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.45, 'deramped', false, ...
  'theta_step_deg', 4);

% --- forward model, per trace so dlam can vary along track
cd_ = cos(THETA_TRUE); sd = sin(THETA_TRUE);
S = struct('hh', zeros(Nz, Nx), 'vv', zeros(Nz, Nx), ...
  'hv', zeros(Nz, Nx), 'vh', zeros(Nz, Nx));
r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
for j = 1:Nx
  ex = exp(1i * gpd * DL_X(j) * z);
  hh = cd_^2 * ex + sd^2;
  vv = sd^2 * ex + cd_^2;
  hv = cd_*sd * (ex - 1);
  S.hh(:, j) = hh .* r(:, j);
  S.vv(:, j) = vv .* r(:, j);
  xc = hv .* r(:, j) + LEAK_C * ((hh + vv)/2) .* r(:, j) + LEAK_D * g(:, j);
  S.hv(:, j) = xc;
  S.vh(:, j) = xc;
end
for f = {'hh','vv','hv','vh'}
  S.(f{1}) = S.(f{1}) + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
end
clear r g;

% --- pass 1: frame-level theta0 and pedestal, as the pipeline runs it
fr = ptt.quadpolFabricLS(S, z, OPTS);
okw = isfinite(fr.theta0);
th_prof = struct('z', fr.zw(okw), 'theta', fr.theta0(okw));
ped = fr.pedestal;
mfr = fr.zw > 300 & fr.zw < 1100;
fprintf(['frame pass: theta0 %5.1f deg (true %.0f), pedestal ' ...
  '[%.2f %+.2fi %.2f], frame dlam %.3f (true mean %.3f)\n'], ...
  mod(rad2deg(angle(mean(exp(2i*fr.theta0(okw & mfr)))))/2, 180), ...
  rad2deg(THETA_TRUE), ped(1), ped(2), ped(3), ...
  median(fr.dlam(mfr), 'omitnan'), mean(DL_X));

% --- pass 2: per-block dlam at three block sizes
BLKS = [125, 28, 14];
fails = 0;
res = struct();
for bi = 1:numel(BLKS)
  NB = BLKS(bi);
  nblk = floor(Nx / NB);
  dl_blk = nan(1, nblk); tr_blk = nan(1, nblk);
  rs_blk = nan(1, nblk); gm_blk = nan(1, nblk);
  o = OPTS; o.theta0 = th_prof;
  if all(isfinite(ped)), o.pedestal = ped; end
  for b = 1:nblk
    j0 = (b-1)*NB + 1; j1 = min(b*NB, Nx);
    Sb = struct('hh', S.hh(:, j0:j1), 'vv', S.vv(:, j0:j1), ...
      'hv', S.hv(:, j0:j1), 'vh', S.vh(:, j0:j1));
    ob = ptt.quadpolFabricLS(Sb, z, o);
    m = ob.zw > 300 & ob.zw < 1100;
    dl_blk(b) = median(ob.dlam(m), 'omitnan');
    rs_blk(b) = median(ob.resid(m), 'omitnan');
    gm_blk(b) = median(ob.gamma(ob.zw > 800 & ob.zw < 1100), 'omitnan');
    tr_blk(b) = mean(DL_X(j0:j1));
  end
  ok = isfinite(dl_blk) & isfinite(tr_blk);
  medae = median(abs(dl_blk(ok) - tr_blk(ok)));
  if nnz(ok) >= 3
    cc = corrcoef(dl_blk(ok), tr_blk(ok)); cc = cc(1, 2);
  else
    cc = NaN;
  end
  fprintf(['NBLK %3d (%4.0f m): %2d blocks | med|err| %.3f  corr %5.2f  ' ...
    'resid %.3f  deep gamma %.2f\n'], NB, NB*dx, nnz(ok), medae, cc, ...
    median(rs_blk(ok)), median(gm_blk(ok), 'omitnan'));
  res.(sprintf('b%d', NB)) = struct('medae', medae, 'corr', cc, ...
    'resid', median(rs_blk(ok)), 'gamma', median(gm_blk(ok), 'omitnan'));
end

% --- verdicts
big = res.b125; small = res.b14;
% the failure must REPRODUCE at Ridge A's trace count...
ok_fail = big.resid > 0.35 || big.medae > 0.04;
fprintf('\n125-trace blocks reproduce the EastGRIP failure: %s\n', ...
  H_tick(ok_fail));
fails = fails + ~ok_fail;
% ...and the failed FRAME pass must poison the small blocks identically:
% if this ever starts passing, the handoff architecture changed and the
% coupling documented above needs re-examining
ok_arch = small.medae > 0.10;
fprintf('small blocks cannot rescue a dead frame pass:      %s\n', ...
  H_tick(ok_arch));
fails = fails + ~ok_arch;

fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_egrip_blocks:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

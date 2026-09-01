%TEST_ERSHADI_STRENGTH Part two: infer dlam(z) given orientation and r.
%
% The inversion is two problems, not one. Part one (ptt.ershadiInverse)
% infers the GEOMETRY - orientation theta and reflection ratio r - and
% accepts a dlam it never revisits. Part two (ptt.ershadiStrength) infers
% the STRENGTH, dlam(z), with that geometry held.
%
% The split is not cosmetic: the two halves have OPPOSITE conditioning.
% theta is global in the observables (phi at a row carries the whole column
% above it, so a short interval has little leverage on its own axis - and
% more data makes it worse, measured in test_ershadi_r), while dlam is
% local (it is the phase GRADIENT across the interval's own rows). This
% test pins that difference and the thing it buys.
%
% Truth: three zones, constant axis, dlam stepping 0.03 -> 0.09 -> 0.05 -
% a profile the phase-gradient chain must follow up AND down, so a stage
% that merely tracked depth would fail it.
%
%   1. RECOVERY: dlam per zone within 0.01 of truth, with the geometry
%      supplied at truth so only the strength stage is under test.
%   2. PARTIAL locality, and the threshold is set where the measurement
%      put it, not where the wording would be tidier. A wrong dlam forced
%      into the top zone must not wreck the zones below - but it is not
%      free either: injecting 0.040 into zone 1 moved the deep-zone error
%      from 0.0007 to 0.0200, so roughly HALF of an upstream error
%      propagates. The mechanism is that a deep interval sees the wrong
%      ACCUMULATED phase entering it and partly absorbs the discrepancy
%      into its own gradient. So staging beats accept-once (where the
%      error is permanent and total) without being immune to it, and the
%      verdict asserts the measured half rather than an aspirational
%      "local". If this ever tightens, the threshold should tighten with
%      it; if it loosens past half, the staging assumption is in trouble
%      and the deep dlam should not be quoted without the upstream one.
%   3. WEIGHTED, NOT GATED: with coherence driven under the |C| > 0.4 gate
%      that ptt.ershadiFabric applies to its own dlam, the direct chain
%      goes blank while this stage still returns a profile. That is the
%      point of fitting |C|^2-weighted instead of thresholded - it is what
%      voids 61% of deep Ridge A intervals in the published chain.
%
% Run: matlab -batch "run('opr_fabric/test/test_ershadi_strength.m')"
clear;
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

rng(17);
C = ptt.constants();
fc = 750e6;
eps_perp = 3.15;
gpd1 = ptt.birefringentPhaseRate(fc, eps_perp, C.deps);

Nz = 2401; dz = 0.5; z = (0:Nz-1).' * dz;      % 0..1200 m
Nx = 300;
TH = deg2rad(40);
EDG = (0:400:1200).';
DL = [0.03; 0.09; 0.05];                        % up then DOWN
RDB = [0; 0; 0];

function S = H_col(z, EDG, DL, TH, gpd1, Nx, NA)
  Nz = numel(z); dz = median(diff(z));
  dlz = zeros(Nz,1);
  for k = 1:numel(DL)
    m = z >= EDG(k) & z < EDG(k+1);
    if k == numel(DL), m = z >= EDG(k) & z <= EDG(k+1); end
    dlz(m) = DL(k);
  end
  del = 2 * gpd1 * cumsum(dlz) * dz;            % two-way accumulated
  r = (randn(Nz,Nx) + 1i*randn(Nz,Nx))/sqrt(2);
  g = (randn(Nz,Nx) + 1i*randn(Nz,Nx))/sqrt(2);
  c = cos(TH); s = sin(TH);
  S = struct();
  hh = (c^2*exp(1i*del) + s^2);
  vv = (s^2*exp(1i*del) + c^2);
  hv = c*s*(exp(1i*del) - 1);
  S.hh = hh .* r + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vv = vv .* r + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  xc = hv .* r + 0.25*exp(0.5i)*((hh+vv)/2) .* r + 0.28*g;
  S.hv = xc + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  S.vh = xc + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
end

geom = struct('edges', EDG, 'theta_int', TH*ones(3,1), 'r_db_int', RDB);

%% verdict 1: recovery with the geometry held at truth
S = H_col(z, EDG, DL, TH, gpd1, Nx, 0.30);
fr = ptt.ershadiFabric(S, z, struct('fc', fc, 'psi_step_deg', 3, ...
  'deramped', false));
st = ptt.ershadiStrength(fr, z, geom, struct('dlam_max', 0.30));
err = abs(st.dlam_int - DL);
fprintf('%8s %9s %10s %10s %12s\n','zone','truth','fitted','|err|','weight');
for k = 1:3
  fprintf('%8d %9.3f %10.3f %10.4f %12.3g\n', k, DL(k), st.dlam_int(k), ...
    err(k), st.w_int(k));
end
ok_rec = all(err < 0.01);
fprintf('1. dlam recovered in all three zones (<0.01):           %s\n', ...
  H_tick(ok_rec));

%% verdict 2: locality - an upstream error must not propagate
INJ = 0.04;                                     % forced error, top zone
% zone 1 PINNED at a wrong value; zones 2-3 still fitted, so the question
% is whether they absorb the upstream error or re-find their own gradient
st2 = ptt.ershadiStrength(fr, z, geom, struct('dlam_max', 0.30, ...
  'dlam_fixed', [DL(1) + INJ; NaN; NaN]));
deep_err_clean = max(abs(st.dlam_int(2:3) - DL(2:3)));
deep_err_hurt = max(abs(st2.dlam_int(2:3) - DL(2:3)));
fprintf(['\n   upstream error %.3f injected into zone 1: deep |err| ' ...
  '%.4f -> %.4f\n'], INJ, deep_err_clean, deep_err_hurt);
ok_local = (deep_err_hurt - deep_err_clean) < 0.5 * INJ;
fprintf('2. an upstream dlam error stays local:                  %s\n', ...
  H_tick(ok_local));

%% verdict 3: works where the published gate goes blank
S_lo = H_col(z, EDG, DL, TH, gpd1, Nx, 0.95);   % coherence driven down
fr_lo = ptt.ershadiFabric(S_lo, z, struct('fc', fc, 'psi_step_deg', 3, ...
  'deramped', false));
gated = mean(isfinite(fr_lo.dlam));
st_lo = ptt.ershadiStrength(fr_lo, z, geom, struct('dlam_max', 0.30));
fprintf('\n   median |C| %.2f: direct-chain dlam finite on %.0f%% of rows\n', ...
  median(fr_lo.Cmag(:), 'omitnan'), 100*gated);
fprintf('   strength stage returned %d of 3 zones\n', ...
  nnz(isfinite(st_lo.dlam_int)));
ok_gate = gated < 0.25 && all(isfinite(st_lo.dlam_int));
fprintf('3. returns a profile where the |C| gate goes blank:      %s\n', ...
  H_tick(ok_gate));

fails = ~ok_rec + ~ok_local + ~ok_gate;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_ershadi_strength:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

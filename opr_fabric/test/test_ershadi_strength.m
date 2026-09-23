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
%      free either: injecting 0.040 into zone 1 moves the deep-zone error
%      from 0.0007 to 0.0199, so roughly HALF of an upstream error
%      propagates. The mechanism is that a deep interval sees the wrong
%      ACCUMULATED phase entering it and partly absorbs the discrepancy
%      into its own gradient. So staging beats accept-once (where the
%      error is permanent and total) without being immune to it, and the
%      verdict asserts the measured half rather than an aspirational
%      "local". The measurement survived ptt.ershadiStrength moving to a
%      PADDED per-interval evaluation (the modelled coherence formed with
%      the interval's real neighbourhood instead of a zero-padded window,
%      verified to reproduce a full-column evaluation to 1e-16 against
%      0.035 in |C| unpadded): 48% of the injection, stable to under a
%      point over five noise seeds.
%      THE THRESHOLD CARRIES A MARGIN over that rather than sitting on it.
%      Drawn at exactly 0.5*INJ it passed by 0.5% of the bar, so any
%      change at all - including a change that makes the fit MORE correct
%      - flips it, and a verdict that cannot survive its own subject being
%      improved is measuring the wrong thing. What it is really for is
%      separating half from TOTAL: accept-once propagates 100%, so a leak
%      drifting toward that is the failure, and the deep dlam should then
%      not be quoted without the upstream one.
%   3. WEIGHTED, NOT GATED: with coherence driven under the |C| > 0.4 gate
%      that ptt.ershadiFabric applies to its own dlam, the direct chain
%      goes blank while this stage still returns a profile. That is the
%      point of fitting |C|^2-weighted instead of thresholded - it is what
%      voids 61% of deep Ridge A intervals in the published chain.
%   4. AN UNRESOLVABLE INTERVAL MUST NOT ASSERT ZERO. With the middle
%      zone's coherence weight zeroed, that zone has nothing to infer
%      from. Two things must then be true, and they pull in opposite
%      directions. It must ABSTAIN - dlam_int NaN, not a number - because
%      zero is the specific physical claim that the ice there is
%      isotropic. And the zone BELOW it must still come out, which means
%      the phase propagated past the gap cannot be zero either: at
%      dlam 0.09 over 400 m the missing two-way phase is 2*gpd1*0.09*400,
%      several radians, and zone 3 fitted against an origin short by that
%      lands on a different branch entirely. So the propagated value is
%      carried separately in .dlam_fwd, interpolated from the intervals
%      that did resolve, and what is REPORTED stays NaN.
%      The gap's own dlam is NOT recoverable and the verdict does not
%      pretend otherwise: interpolating 0.03 and 0.05 across a zone whose
%      truth is 0.09 gives 0.046, and zone 3 inherits about half of that
%      shortfall, exactly as verdict 2 says it must. What is asserted is
%      the comparison the choice actually decides. MEASURED: carrying the
%      interpolated value puts zone 3 at |err| 0.0122 against a truth of
%      0.05; propagating zero instead puts it at 0.0204, two thirds
%      worse. Both are consistent with verdict 2's half-leak - zero is
%      simply a bigger upstream error to leak.
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
  '%.4f -> %.4f (%.0f%% of the injection)\n'], INJ, deep_err_clean, ...
  deep_err_hurt, 100*(deep_err_hurt - deep_err_clean)/INJ);
ok_local = (deep_err_hurt - deep_err_clean) < 0.6 * INJ;
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

%% verdict 4: an unresolvable interval abstains without asserting zero
fr_gap = fr;
gap = z >= EDG(2) & z < EDG(3);            % zone 2 loses all its weight
fr_gap.Cmag(gap, :) = 0;
st4 = ptt.ershadiStrength(fr_gap, z, geom, struct('dlam_max', 0.30));
% the same column with the gap PINNED at zero: what propagating a
% fabricated isotropic layer costs the interval below it
st4z = ptt.ershadiStrength(fr_gap, z, geom, struct('dlam_max', 0.30, ...
  'dlam_fixed', [NaN; 0; NaN]));
e4 = abs(st4.dlam_int(3) - DL(3));
e4z = abs(st4z.dlam_int(3) - DL(3));
fprintf('\n   zone 2 weight zeroed: w_int %.3g, dlam_int %s, dlam_fwd %.3f\n', ...
  st4.w_int(2), mat2str(st4.dlam_int(2)), st4.dlam_fwd(2));
fprintf('   zone 3 below the gap: %.4f carried vs %.4f with zero propagated (truth %.3f)\n', ...
  st4.dlam_int(3), st4z.dlam_int(3), DL(3));
fprintf('   |err| %.4f carried vs %.4f zero\n', e4, e4z);
ok_abst = isnan(st4.dlam_int(2)) && st4.dlam_fwd(2) > 0 && ...
  st4.w_int(2) == 0 && e4 < 0.8 * e4z;
fprintf('4. an unresolvable interval abstains, not asserts zero:  %s\n', ...
  H_tick(ok_abst));

fails = ~ok_rec + ~ok_local + ~ok_gate + ~ok_abst;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_ershadi_strength:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

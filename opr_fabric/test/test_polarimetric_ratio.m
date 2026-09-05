%TEST_POLARIMETRIC_RATIO Orientation and strength from channel power ratios.
%
% ptt.polarimetricRatioInverse solves the Fujita model, at a depth-constant
% orientation, for theta0, dlam(z) and the scattering ratio r(z) from the
% RATIOS between measured channels, with each ratio carrying a free
% additive offset so that unknown H-versus-V channel gains drop out.
%
% The synthetic injects gains the solver is never told: |t|^2 = +6 dB on
% the transmit-V chain and |r|^2 = -4 dB on receive-V, which is the scale
% of the real imbalance on this system (ptt.equaliseChannels).
%
%   1. RECOVERY. With four headings and power + phase, theta0 within 1 deg,
%      dlam within 0.005, r within 0.05, pedestal within 0.05.
%   2. THE 90-DEGREE ALIAS IS REAL. Power ratios alone cannot separate
%      (theta0, r) from (theta0 + 90, 1/r). It is an EXACT symmetry of the
%      forward model, not an approximate one: under d -> d + 90 the roles
%      of cos^2 and sin^2 swap, and with r -> 1/r every one of Phh, Pvv
%      and Pxx is divided by r^2, so every power RATIO is unchanged and
%      the free per-block offsets have nothing left to absorb.
%      The verdict is therefore stated on the DATA, not on where a
%      particular multi-start happens to land: a power-only fit pinned to
%      the aliased axis must stay there and reproduce every ratio to far
%      inside the 0.5 dB noise that case 4 calls realistic, and
%      out.alias_unresolved must be flagged whenever phi is absent.
%      An earlier form of this verdict asserted that the fit LANDS on the
%      alias (theta0 error > 45 deg). That held only while the model side
%      of the ratios was unfloored: the residuals the model's own power
%      nulls then generated (0.8 dB rms on noise-free data the model
%      reproduces exactly) swamped everything and the winning branch was
%      arbitrary. With model and data floored symmetrically the true
%      branch does win - but it wins on the SMOOTHNESS PRIOR on r, since
%      1/r is curved where r is straight, and not on any evidence in the
%      data. That is not axis information, which is exactly why the flag
%      must stay and the phase term must not be deleted.
%   3. THE PHASE BREAKS IT. Adding phi recovers the true axis.
%   4. NOISE. 0.5 dB on every power and 0.05 rad on the phase leaves the
%      answers inside the case-1 tolerances at 3x.
%   5. ONE HEADING IS NOT ENOUGH FOR THE AXIS. A counting argument says it
%      should be; it is not, because with a single d the fit trades
%      orientation against r(z), which is free at every depth. The test
%      asserts BOTH halves of the measured behaviour - strength still
%      comes out, orientation does not - so that nobody quotes an axis
%      from a single-heading frame on the strength of the counting
%      argument.
%   6. A PHASE OFFSET AT THE BRANCH CUT. The header calls the instrument
%      phase offset "one more free constant", so an offset near +-pi must
%      cost nothing. It does only if the objective WRAPS the phi residual.
%      Differenced instead, a datum at +3.1 rad and a model at -3.1 scores
%      ~2*pi where the true misfit is ~0, and because those entries sit
%      where the phase is near the cut rather than where the model is
%      wrong, they drag theta0 - the failure ptt.fabricGLS measured and
%      wrapped away. Case 1 injects no offset, so nothing else here covers
%      the wrap. Verdict: an offset of 3.0 rad leaves case 1's answers
%      inside case 1's own tolerances.
%   7. NULLING THE PHASE WEIGHT. That wrap divides by w_phi to undo the
%      scaling the phi observable is carried with, so the documented
%      w_phi = 0 - the natural way to drop the phase term while leaving
%      'phi' in `use` - turns every phi residual into 0/0. The NaN then
%      travels: every start scores NaN, none beats the incumbent, and the
%      failure surfaces as an index error rather than as anything a
%      caller could act on. Zeroing the residual alone still leaves the
%      block's free offset with an identically zero Jacobian column and
%      so a singular normal matrix; the observable has to be dropped.
%      Verdict: w_phi = 0 reproduces the power-only fit exactly, reports
%      the alias as unresolved, and a negative weight is rejected by name.
%
% Run: matlab -batch "run('opr_fabric/test/test_polarimetric_ratio.m')"
clear; t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

fc = 750e6; C = ptt.constants();
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar)*C.c*1e9);
z = (10:10:1500).'; dz = [0; diff(z)];
TH = deg2rad(35);
DL = 0.03 + 0.04*(z/z(end));
RR = 1 + 0.4*(z/z(end));
PED = 0.30;
GT = 10^(6/10);    % transmit V gain, +6 dB, unknown to the solver
GR = 10^(-4/10);   % receive V gain, -4 dB
delta = cumsum(DL.*dz)*gpd;
OPT = struct('fc', fc, 'n_r', 6, 'n_d', 10, 'eta_r', 3e-2, 'eta_d', 3e-2, ...
  'max_iter', 80, 'n_start', 6);

o1 = H_fit(H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0, 0), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv','phi'});
[e1, d1, r1] = H_err(o1, TH, DL, RR);
fprintf('1. four headings, power + phase: theta0 err %.2f deg, dlam err %.4f, r err %.3f, ped %.2f\n', ...
  e1, d1, r1, o1.ped);
ok1 = e1 < 1 && d1 < 0.005 && r1 < 0.05 && abs(o1.ped - PED) < 0.05;
fprintf('   recovers orientation AND strength through unknown gains: %s\n', H_tick(ok1));

obs2 = H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0, 0);
OPT2 = OPT; OPT2.theta0 = TH + pi/2;      % pinned to the ALIASED axis
o2 = H_fit(obs2, z, OPT2, {'vv_vh','hh_hv','hh_vv'});
[e2, ~, ~] = H_err(o2, TH, DL, RR);
rms2 = max([o2.rms.vv_vh, o2.rms.hh_hv, o2.rms.hh_vv]);
fprintf('\n2. power only, started on the alias: theta0 %.1f deg from truth, worst ratio rms %.4f dB, alias_unresolved %d\n', ...
  e2, rms2, o2.alias_unresolved);
% it must STAY on the alias (the data gives it no reason to leave) and it
% must fit the ratios to far inside the 0.5 dB case 4 calls realistic
ok2 = e2 > 45 && rms2 < 0.05 && o2.alias_unresolved;
fprintf('   the 90-deg alias is REAL and is flagged:              %s\n', H_tick(ok2));
fprintf('3. adding the coherence phase resolves it:               %s\n', H_tick(e1 < 1));
ok3 = e1 < 1;

rng(7);
o4 = H_fit(H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0.5, 0.05), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv','phi'});
[e4, d4, r4] = H_err(o4, TH, DL, RR);
fprintf('\n4. with noise: theta0 err %.2f deg, dlam err %.4f, r err %.3f\n', e4, d4, r4);
ok4 = e4 < 3 && d4 < 0.015 && r4 < 0.15;
fprintf('   degrades gracefully:                                  %s\n', H_tick(ok4));

o5 = H_fit(H_obs(deg2rad(15), TH, RR, delta, PED, GT, GR, 0, 0), z, OPT, ...
  {'vv_vh','hh_hv','hh_vv','phi'});
[e5, d5, ~] = H_err(o5, TH, DL, RR);
fprintf('\n5. single heading: theta0 err %.2f deg, dlam err %.4f (four headings gave %.2f, %.4f)\n', ...
  e5, d5, e1, d1);
ok5 = d5 < 0.03 && e5 > 5 && e5 > 10*e1;
fprintf('   strength survives one heading, orientation does NOT:   %s\n', H_tick(ok5));

% offset applied where the observation is FORMED, then wrapped as any real
% instrument's phase would be, so the solver sees it exactly as it would
% on real data - straddling the cut, not as a tidy additive constant
PH_OFF = 3.0;
obs6 = H_obs(deg2rad([0 40 80 120]), TH, RR, delta, PED, GT, GR, 0, 0);
obs6.phi = angle(exp(1i*(obs6.phi + PH_OFF)));
o6 = H_fit(obs6, z, OPT, {'vv_vh','hh_hv','hh_vv','phi'});
[e6, d6, r6] = H_err(o6, TH, DL, RR);
fprintf('\n6. %.1f rad instrument phase offset: theta0 err %.2f deg, dlam err %.4f, r err %.3f\n', ...
  PH_OFF, e6, d6, r6);
fprintf('   (%.0f%% of phi cells land within 0.2 rad of the +-pi cut)\n', ...
  100*mean(abs(abs(obs6.phi(:)) - pi) < 0.2));
ok6 = e6 < 1 && d6 < 0.005 && r6 < 0.05;
fprintf('   a free phase constant at the branch cut costs nothing: %s\n', H_tick(ok6));

% ---- verdict 7: w_phi = 0 nulls the phase term, it does not poison it
% The phi observable is carried pre-scaled by w_phi, so the branch-cut
% wrap has to divide by w_phi to recover the raw angle. At w_phi = 0 both
% sides are identically zero and that division is 0/0: every phi residual
% goes NaN, every start scores NaN, no start ever beats the incumbent and
% the solver falls out at an index error instead of a diagnosable one.
% Zeroing the residual instead is not enough either: a zero-weighted phi
% block keeps its free additive offset, whose Jacobian column is then
% identically zero, so the normal matrix is exactly singular and the
% Gauss-Newton step comes back Inf on the first iteration - the caller
% gets a coarse multi-start grid node back, not a fit. The observable has
% to be DROPPED, which is the power-only problem exactly, and the 90-deg
% alias must then be reported as unresolved.
OPT0 = OPT; OPT0.w_phi = 0;
ws = warning('off', 'ptt:polarimetricRatioInverse:phiNulled');
oc = onCleanup(@() warning(ws));
o7 = H_fit(obs6, z, OPT0, {'vv_vh','hh_hv','hh_vv','phi'});
clear oc;
o7p = H_fit(obs6, z, OPT, {'vv_vh','hh_hv','hh_vv'});
[e7, d7, r7] = H_err(o7, TH, DL, RR);
dth7 = rad2deg(abs(angle(exp(2i*(o7.theta0 - o7p.theta0)))/2));
ok7 = isfinite(e7) && isfinite(d7) && isfinite(r7) && dth7 < 1e-6 && ...
  ~any(strcmp(o7.used, 'phi')) && o7.alias_unresolved;
fprintf(['\n7. w_phi = 0: theta0 err %.2f deg, dlam err %.4f, r err %.3f; ' ...
  'differs from power-only by %.1e deg\n'], e7, d7, r7, dth7);
OPTN = OPT; OPTN.w_phi = -1;
try
  H_fit(obs6, z, OPTN, {'hh_vv','phi'});
  ok7 = false;
catch ME
  ok7 = ok7 && strcmp(ME.identifier, 'ptt:polarimetricRatioInverse:wPhi');
end
fprintf('   a nulled phase weight drops the term, negative errors:  %s\n', H_tick(ok7));

fails = ~ok1 + ~ok2 + ~ok3 + ~ok4 + ~ok5 + ~ok6 + ~ok7;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_polarimetric_ratio:failed', '%d verdict(s) failed', fails);
end

function obs = H_obs(gam, th, r, delta, ped, gt, gr, noise_db, noise_rad)
d = th - gam; c2 = cos(d).^2; s2 = sin(d).^2; cd = cos(delta);
Phh = c2.^2 + (r.^2).*s2.^2 + 2*r.*c2.*s2.*cd;
Pvv = (r.^2).*c2.^2 + s2.^2 + 2*r.*c2.*s2.*cd;
Pxx = (1 + r.^2 - 2*r.*cd).*(c2.*s2) + ped*(Phh+Pvv)/2;
t2 = s2 ./ c2; t4 = t2.^2;
phi = atan2(r.*sin(delta).*(1-t4), r.*cos(delta).*(1+t4) + t2.*(1+r.^2));
n = @(A) A .* 10.^(noise_db*randn(size(A))/10);
obs = struct('gamma', gam, 'P_hh', n(Phh), 'P_vv', n(Pvv*gr*gt), ...
  'P_hv', n(Pxx*gt), 'P_vh', n(Pxx*gr), 'phi', phi + noise_rad*randn(size(phi)));
end

function o = H_fit(obs, z, OPT, use)
OPT.use = use;
o = ptt.polarimetricRatioInverse(obs, z, OPT);
end

function [e, d, r] = H_err(o, TH, DL, RR)
e = rad2deg(abs(angle(exp(2i*(o.theta0 - TH)))/2));
d = median(abs(o.dlam_z - DL));
r = median(abs(o.r_z - RR));
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

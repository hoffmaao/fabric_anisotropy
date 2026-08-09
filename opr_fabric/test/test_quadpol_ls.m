%TEST_QUADPOL_LS Round-trip and leakage-immunity test of ptt.quadpolFabricLS.
%
% Same synthetic column as test_ershadi (known axis, known contrast, seen
% from several antenna azimuths), plus the failure mode that motivated the
% estimator: an antenna-fixed leakage term added to the cross-polarized
% channels at the level the real system shows (cross/co ~ -4 dB, flat
% with depth, reciprocal), built as a MIX of a component correlated with
% the co-polarized speckle and one decorrelated from it. The decorrelated
% part matters: it enters the synthesized T_hh and T_vv with opposite
% signs, so it collapses the measured |C| azimuth-selectively at the
% birefringent crossings - which is what produced the "regular jumps"
% (crossing-depth dlam overshoots, gamma ~0.2, resid ~0.6) in the Ridge A
% validation before the CRB weighting. The spike assertion below guards
% that failure mode. ptt.ershadiFabric takes its axis from the
% cross-polarized minimum, which the leakage owns, so it locks; the LS fit
% never touches cross-polarized power and must not.
%
% Cases:
%   1. clean, dlam 0.30              - parity with test_ershadi
%   2. leakage, dlam 0.30            - the lock case; ershadi error reported
%   3. leakage, dlam 0.05            - Ridge A scale, crossing spikes checked
%   4. leakage+noise, dlam 0.05      - coherence floor near the real ~0.5
%   5. leakage, dlam 0 (isotropic)   - no folded floor, theta0 abstains
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_ls.m')"
clear;
rng(11);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);   % rad/m per unit dlam

THETA_TRUE = deg2rad(35);
Nz = 3000;
z = (0:Nz-1).' * 0.5;
Nx = 60;
LEAK_C = 0.30 * exp(0.7i);   % pedestal part correlated with co-pol speckle
LEAK_D = 0.35;               % pedestal part with its OWN speckle - this is
                             % the component that collapses |C| at crossings

% lean grids so the whole file runs in a few minutes; the polish restores
% the precision the coarse grid gives up
OPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.35, 'deramped', false, ...
  'theta_step_deg', 4);

cases = { ...
  'clean  dlam 0.30', 0.30, false, [0 20 55], 0.02,  5,   0.02; ...
  'leak   dlam 0.30', 0.30, true,  [0 20],    0.02,  5,   0.02; ...
  'leak   dlam 0.05', 0.05, true,  20,        0.010, 8,   0.02; ...
  'noisy  dlam 0.05', 0.05, true,  20,        0.015, 12,  0.50; ...
  'ped45  dlam 0.05', 0.05, true,  -10,       0.012, 10,  0.02; ...
  'leak   isotropic', 0.00, true,  20,        0.010, NaN, 0.02};
% ped45: the fabric axis sits at EXACTLY 45 deg to the antenna axes
% (d = 35 - (-10) = 45), where mu = cos 2(psi - theta0) = sin 2psi and the
% pedestal's azimuthal shape coincides with the fabric phase term. This is
% the geometry where a PER-WINDOW pedestal estimate is confounded with the
% fabric within each window; the frame-level pedestal (constant with depth
% while sin delta0 sweeps fringes across windows) is what separates them.

fails = 0;
for ci = 1:size(cases, 1)
  [name, DL, leak, alphas, tol_d, tol_t, na] = cases{ci, :};
  delta = gpd * DL * z;
  for alpha_deg = alphas
    d = THETA_TRUE - deg2rad(alpha_deg);
    cd_ = cos(d); sd = sin(d);
    ex = exp(1i * delta); ey = ones(Nz, 1);
    hh = cd_^2 * ex + sd^2 * ey;
    vv = sd^2 * ex + cd_^2 * ey;
    hv = cd_*sd * (ex - ey);
    r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
    S = struct();
    S.hh = hh .* r + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));
    S.vv = vv .* r + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));
    xc = hv .* r;
    if leak
      g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
      xc = xc + LEAK_C * ((hh + vv)/2) .* r + LEAK_D * g;
    end
    S.hv = xc + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));
    S.vh = xc + na*(randn(Nz,Nx)+1i*randn(Nz,Nx));

    out = ptt.quadpolFabricLS(S, z, OPTS);

    mid = out.zw > 300 & out.zw < 1300;
    dl_got = median(out.dlam(mid), 'omitnan');
    ok_d = abs(dl_got - DL) < tol_d;

    % crossing-spike guard, on the case whose fringes the windows resolve:
    % no window may sit far above truth, and abstention must not gut the
    % coverage either. This is the exact regression seen on real Ridge A
    % (dlam 0.13 windows at the crossings) before the CRB weighting.
    if abs(DL - 0.05) < 1e-9
      mm = out.zw > 250 & out.zw < 1350;
      pk = max(out.dlam(mm), [], 'omitnan');
      fin = mean(isfinite(out.dlam(mm)));
      ok_s = pk < DL + 0.04 && fin > 0.6;
      fails = fails + ~ok_s;
      fprintf('   [crossing spikes: max window dlam %.3f, finite %2.0f%% %s]\n', ...
        pk, 100*fin, H_tick(ok_s));
    end

    want = mod(rad2deg(d), 180);
    th = out.theta0(mid);
    th = th(isfinite(th));
    if isfinite(tol_t)
      got = mod(rad2deg(angle(mean(exp(2i*th))))/2, 180);
      terr = abs(mod(got - want + 90, 180) - 90);
      ok_t = terr < tol_t;
      tmsg = sprintf('theta %6.1f (want %6.1f, err %4.1f) %s', got, want, ...
        terr, H_tick(ok_t));
    else
      % isotropic: theta0 must ABSTAIN (q gate), not report confidently
      frac = numel(th) / max(1, nnz(mid));
      ok_t = frac < 0.5;
      tmsg = sprintf('theta reported on %2.0f%% of windows %s', 100*frac, ...
        H_tick(ok_t));
    end
    fails = fails + ~(ok_t && ok_d);
    fprintf('%s a=%3d: %s | dlam %.3f (want %.2f) %s\n', name, alpha_deg, ...
      tmsg, dl_got, DL, H_tick(ok_d));

    % the two-pass path the pipeline section uses: theta0 handed back in,
    % only (delta0, ddelta) free. Must reproduce the free fit's contrast.
    if DL > 0 && isfinite(tol_t)
      okw = isfinite(out.theta0);
      if nnz(okw) >= 2
        o2 = ptt.quadpolFabricLS(S, z, setfield(OPTS, 'theta0', ...
          struct('z', out.zw(okw), 'theta', out.theta0(okw)))); %#ok<SFLD>
        dl2 = median(o2.dlam(mid), 'omitnan');
        ok2 = abs(dl2 - DL) < tol_d;
        fails = fails + ~ok2;
        fprintf('   [theta0-fixed pass: dlam %.3f %s]\n', dl2, H_tick(ok2));
      end
    end

    % the motivating comparison, reported not asserted: where does the
    % published chain put the axis under the same leakage?
    if leak && DL > 0 && alpha_deg == alphas(1)
      oe = ptt.ershadiFabric(S, z, struct('fc', fc, 'psi_step_deg', 2, ...
        'win_m', 30, 'grad_win_m', 25, 'coh_min', 0.4, 'deramped', false));
      me = z > 300 & z < 1300;
      the = rad2deg(oe.theta(me)); the = the(isfinite(the));
      ge = mod(rad2deg(angle(mean(exp(2i*deg2rad(2*the)))))/4, 90);
      ee = abs(mod(ge - mod(want,90) + 45, 90) - 45);
      de = median(oe.dlam(me), 'omitnan');
      fprintf(['   [ershadiFabric on the same data: theta err %4.1f deg, ' ...
        'dlam %.3f - the lock this estimator removes]\n'], ee, de);
    end
  end
end

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_quadpol_ls:failed', '%d case(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

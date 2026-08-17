%TEST_EGRIP_CAP Does the default dlam_max = 0.25 cause the EastGRIP failure?
%
% The EastGRIP validation frame failed with a signature the block-size
% synthetic did NOT reproduce: HIGH residuals (0.30-0.64, worst mid-column)
% with dlam OSCILLATING between ~0.05 and ~0.19 by depth band. The EGRIP
% core puts the true contrast near 0.2-0.35 - and ptt.quadpolFabricLS caps
% its ddelta search at dlam_max, DEFAULT 0.25, which the pipeline never
% overrides. Ridge A (0.07) never came near the cap; EastGRIP is the first
% site above it. A window whose true rate exceeds the cap either rails
% there and abstains, or locks an aliased branch below it: oscillation and
% high resid, exactly as observed.
%
% Test: EastGRIP-like synthetic with CONSTANT dlam = 0.32 (no lateral
% structure at all, so the block question is out of the picture), the same
% two-pass the pipeline runs, at the PIPELINE'S defaults (dlam_max 0.25)
% and with the cap raised to 0.45.
%
% PASS = defaults REPRODUCE the failure signature (biased-low or
%        oscillating dlam / elevated resid / abstention), and the raised
%        cap RECOVERS 0.32 cleanly at BOTH 125- and 14-trace blocks.
%
% Run: matlab -batch "run('opr_fabric/test/test_egrip_cap.m')"
clear;
rng(5);
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);

THETA_TRUE = deg2rad(35);
DL_TRUE = 0.32;                % EGRIP-core scale, ABOVE the 0.25 default cap
Nz = 2800; dz = 0.5;
z = (0:Nz-1).' * dz;
Nx = 250;                      % 2.25 km at 9 m; lateral structure NONE
LEAK_C = 0.30 * exp(0.7i);
LEAK_D = 0.35;
NA = 0.5;

cd_ = cos(THETA_TRUE); sd = sin(THETA_TRUE);
ex = exp(1i * gpd * DL_TRUE * z);
hh0 = cd_^2 * ex + sd^2;
vv0 = sd^2 * ex + cd_^2;
hv0 = cd_*sd * (ex - 1);
r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
S = struct();
S.hh = hh0 .* r + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
S.vv = vv0 .* r + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
xc = hv0 .* r + LEAK_C * ((hh0 + vv0)/2) .* r + LEAK_D * g;
S.hv = xc + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
S.vh = xc + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
clear r g xc;

% grids lean enough to run in minutes but with the PIPELINE's step and cap
BASE = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 15, 'deramped', false, 'theta_step_deg', 4);

fails = 0;
for pass = 1:2
  o = BASE;
  if pass == 2, o.dlam_max = 0.45; end
  capv = 0.25; if pass == 2, capv = 0.45; end

  frq = ptt.quadpolFabricLS(S, z, o);
  m = frq.zw > 300 & frq.zw < 1100;
  fin = mean(isfinite(frq.dlam(m)));
  dlm = median(frq.dlam(m), 'omitnan');
  rsd = median(frq.resid(m), 'omitnan');
  thf = frq.theta0(m); thf = thf(isfinite(thf));
  terr = abs(mod(rad2deg(angle(mean(exp(2i*thf))))/2 - ...
    rad2deg(THETA_TRUE) + 90, 180) - 90);
  % band-to-band oscillation, the observed signature: medians per 200 m
  bands = 300:200:1100;
  bm = nan(1, numel(bands)-1);
  for k = 1:numel(bands)-1
    mm = frq.zw > bands(k) & frq.zw < bands(k+1);
    bm(k) = median(frq.dlam(mm), 'omitnan');
  end
  osc = max(bm) - min(bm);
  fprintf(['cap %.2f: frame dlam %5.3f (true %.2f)  resid %.3f  ' ...
    'finite %3.0f%%  theta err %4.1f  band swing %.3f\n'], ...
    capv, dlm, DL_TRUE, rsd, 100*fin, terr, osc);

  if pass == 1
    % the failure must reproduce at the default cap
    ok = (dlm < 0.27) || (osc > 0.06) || (fin < 0.7) || (rsd > 0.35);
    fprintf('  -> default cap reproduces the EastGRIP signature: %s\n', ...
      H_tick(ok));
    fails = fails + ~ok;
  else
    okf = abs(dlm - DL_TRUE) < 0.03 && rsd < 0.35 && fin > 0.8;
    fprintf('  -> raised cap recovers the frame fit:            %s\n', ...
      H_tick(okf));
    fails = fails + ~okf;
    % and per-block, both sizes, with the frame handoff as the pipeline
    okw = isfinite(frq.theta0);
    th_prof = struct('z', frq.zw(okw), 'theta', frq.theta0(okw));
    for NB = [125, 14]
      nblk = floor(Nx / NB);
      dlb = nan(1, nblk); rsb = nan(1, nblk);
      ob_o = o; ob_o.theta0 = th_prof;
      if all(isfinite(frq.pedestal)), ob_o.pedestal = frq.pedestal; end
      for b = 1:nblk
        j0 = (b-1)*NB + 1; j1 = min(b*NB, Nx);
        Sb = struct('hh', S.hh(:, j0:j1), 'vv', S.vv(:, j0:j1), ...
          'hv', S.hv(:, j0:j1), 'vh', S.vh(:, j0:j1));
        ob = ptt.quadpolFabricLS(Sb, z, ob_o);
        mm = ob.zw > 300 & ob.zw < 1100;
        dlb(b) = median(ob.dlam(mm), 'omitnan');
        rsb(b) = median(ob.resid(mm), 'omitnan');
      end
      me = median(abs(dlb - DL_TRUE), 'omitnan');
      fprintf(['     blocks of %3d: med|err| %.3f  resid %.3f  ' ...
        'spread %.3f-%.3f\n'], NB, me, median(rsb, 'omitnan'), ...
        min(dlb), max(dlb));
      okb = me < 0.04;
      fails = fails + ~okb;
    end
  end
end

fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_egrip_cap:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

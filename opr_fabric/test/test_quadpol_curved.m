%TEST_QUADPOL_CURVED The curving-line adaptation of the quad-pol chain.
%
% A straight-line frame can be fitted in the antenna frame because that
% frame is constant. On a CURVING line the fabric's azimuthal signature
% rotates through the antennas along the drive, so antenna-frame moment
% averaging smears it away - fitted dlam collapses toward zero and theta0
% is meaningless (the symptom observed on the turning Ridge A profiles).
% The adaptation fits in the GEOGRAPHIC frame: per-sub-block moments are
% rotated north-referenced (ptt.rotateMoments) before accumulation, and
% the antenna-fixed pedestal - which that rotation carries from sin2(phi)
% to sin2(psi - h) - is supplied as a precomputed field mixed over the
% measured heading distribution, with the survey-calibrated antenna-frame
% coefficients pinned (no new unknowns).
%
% Synthetic: a 90-degree arc (headings 30 -> 120 deg across the traces),
% fabric fixed GEOGRAPHICALLY at 35 deg with dlam 0.05, and the usual
% antenna-fixed reciprocal pedestal. Cases:
%   control  straight line, the standard path (also yields the pedestal
%            calibration the curved fit pins, as production pins the
%            survey value)
%   naive    curved line through the standard antenna-frame path - must
%            REPRODUCE the collapse (reported, dlam suppression asserted)
%   curveA   curved, geographic moments, NO pedestal - isolates the
%            rotation correctness
%   curveB   curved, geographic moments, pedestal field pinned - the
%            production path for curving lines
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_curved.m')"
clear;
rng(11);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);

TH_GEO = 35;                 % fabric axis, deg E of N, fixed in the ice
DL = 0.05;
Nz = 3000;
z = (0:Nz-1).' * 0.5;
Ntr = 1200;
NA = 0.02;
LEAK_C = 0.30 * exp(0.7i);
LEAK_D = 0.35;
NSUB = 100;                  % traces per rotation sub-block

OPTS = struct('fc', fc, 'psi_step_deg', 4, 'win_short_m', 10, ...
  'win_fit_m', 60, 'step_m', 40, 'dlam_max', 0.35, 'deramped', false, ...
  'theta_step_deg', 4);
psi = (0:OPTS.psi_step_deg:180-OPTS.psi_step_deg) * pi/180;

delta = gpd * DL * z;
ex = exp(1i * delta); ey = ones(Nz, 1);

function S = H_synth(h_deg, th_geo, ex, ey, Nz, Ntr, NA, LEAK_C, LEAK_D)
% Antenna-frame channels for per-trace headings h_deg (1 x Ntr): the
% fabric angle in the antenna frame is th_geo - h at each trace, while
% the pedestal is constant in the antenna frame by construction.
S = struct('hh', zeros(Nz, Ntr), 'vv', zeros(Nz, Ntr), ...
  'hv', zeros(Nz, Ntr), 'vh', zeros(Nz, Ntr));
r = (randn(Nz, Ntr) + 1i*randn(Nz, Ntr)) / sqrt(2);
g = (randn(Nz, Ntr) + 1i*randn(Nz, Ntr)) / sqrt(2);
for j = 1:Ntr
  d = deg2rad(th_geo - h_deg(j));
  cd_ = cos(d); sd = sin(d);
  hh = cd_^2*ex + sd^2*ey;
  vv = sd^2*ex + cd_^2*ey;
  hv = cd_*sd*(ex - ey);
  xc = hv.*r(:,j) + LEAK_C*((hh + vv)/2).*r(:,j) + LEAK_D*g(:,j);
  S.hh(:,j) = hh.*r(:,j);
  S.vv(:,j) = vv.*r(:,j);
  S.hv(:,j) = xc;
  S.vh(:,j) = xc;
end
S.hh = S.hh + NA*(randn(Nz,Ntr)+1i*randn(Nz,Ntr));
S.vv = S.vv + NA*(randn(Nz,Ntr)+1i*randn(Nz,Ntr));
S.hv = S.hv + NA*(randn(Nz,Ntr)+1i*randn(Nz,Ntr));
S.vh = S.vh + NA*(randn(Nz,Ntr)+1i*randn(Nz,Ntr));
end

function Mg = H_geo_moments(S, h_deg, NSUB)
% Sub-block moments rotated into the common north-referenced frame.
Ntr = size(S.hh, 2);
Mg = [];
wsum = 0;
for j0 = 1:NSUB:Ntr
  j1 = min(j0 + NSUB - 1, Ntr);
  Sb = struct('hh', S.hh(:,j0:j1), 'vv', S.vv(:,j0:j1), ...
    'hv', S.hv(:,j0:j1), 'vh', S.vh(:,j0:j1));
  Mb = ptt.quadpolMoments(Sb, [1 (j1-j0+1)]);
  hb = mod(rad2deg(angle(mean(exp(2i*deg2rad(h_deg(j0:j1))))))/2, 180);
  Mr = ptt.rotateMoments(Mb, -deg2rad(hb));
  if isempty(Mg), Mg = zeros(size(Mr)); end
  Mg = Mg + Mr * (j1-j0+1);
  wsum = wsum + (j1-j0+1);
end
Mg = Mg / wsum;
end

fails = 0;
mid = @(o) o.zw > 300 & o.zw < 1300;
axerr = @(got, want) abs(mod(got - want + 90, 180) - 90);
cmed = @(th) mod(rad2deg(angle(mean(exp(2i*th(isfinite(th))))))/2, 180);

% --- control: straight line, standard path; calibrates the pedestal
h_straight = 20 * ones(1, Ntr);
S = H_synth(h_straight, TH_GEO, ex, ey, Nz, Ntr, NA, LEAK_C, LEAK_D);
oc = ptt.quadpolFabricLS(S, z, OPTS);
m = mid(oc);
th_ant = cmed(oc.theta0(m));
dlc = median(oc.dlam(m), 'omitnan');
ok = axerr(th_ant + 20, TH_GEO) < 4 && abs(dlc - DL) < 0.012;
fails = fails + ~ok;
fprintf('control (straight): theta_geo %5.1f dlam %.3f  pedestal [%.2f %+.2f %.2f]  %s\n', ...
  mod(th_ant + 20, 180), dlc, oc.pedestal(1), oc.pedestal(2), ...
  oc.pedestal(3), H_tick(ok));
ped_cal = oc.pedestal;

% --- curved line
h_curve = linspace(30, 120, Ntr);
S = H_synth(h_curve, TH_GEO, ex, ey, Nz, Ntr, NA, LEAK_C, LEAK_D);

% naive: standard antenna-frame path on the smeared moments. The fabric
% signature must largely cancel; assert the suppression so this test
% keeps documenting WHY the adaptation exists.
on = ptt.quadpolFabricLS(S, z, OPTS);
dln = median(on.dlam(mid(on)), 'omitnan');
% Frame-level smearing suppresses partially (the hard per-block zeros in
% production come from pinning the frame theta0 to blocks whose heading
% deviates - cos 2*dev - on top of this); assert the suppression exists.
ok = dln < 0.85 * DL;
fails = fails + ~ok;
fprintf('naive  (curved, antenna frame): dlam %.3f (truth %.2f, suppressed) %s\n', ...
  dln, DL, H_tick(ok));

% adapted, no pedestal in the synthetic: rotation correctness alone
Snl = H_synth(h_curve, TH_GEO, ex, ey, Nz, Ntr, NA, 0, 0);
Mg = H_geo_moments(Snl, h_curve, NSUB);
oa = ptt.quadpolFabricLS(struct('M', Mg), z, ...
  setfield(OPTS, 'pedestal', 'window')); %#ok<SFLD>
m = mid(oa);
tha = cmed(oa.theta0(m));
dla = median(oa.dlam(m), 'omitnan');
ok = axerr(tha, TH_GEO) < 4 && abs(dla - DL) < 0.012;
fails = fails + ~ok;
fprintf('curveA (geo moments, no leak): theta_geo %5.1f dlam %.3f  %s\n', ...
  tha, dla, H_tick(ok));

% adapted, with the pedestal: pinned field from the control calibration
% mixed over the measured heading distribution
Mg = H_geo_moments(S, h_curve, NSUB);
c2 = mean(cosd(2*h_curve)); s2 = mean(sind(2*h_curve));
c4 = mean(cosd(4*h_curve)); s4 = mean(sind(4*h_curve));
pf = (ped_cal(1) + 1i*ped_cal(2)) * (c2*sin(2*psi(:)) - s2*cos(2*psi(:))) ...
  + ped_cal(3) * (0.5 - 0.5*(c4*cos(4*psi(:)) + s4*sin(4*psi(:))));
ob = ptt.quadpolFabricLS(struct('M', Mg), z, ...
  setfield(OPTS, 'pedestal', pf)); %#ok<SFLD>
m = mid(ob);
thb = cmed(ob.theta0(m));
dlb = median(ob.dlam(m), 'omitnan');
ok = axerr(thb, TH_GEO) < 5 && abs(dlb - DL) < 0.012;
fails = fails + ~ok;
fprintf('curveB (geo moments + pinned pedestal field): theta_geo %5.1f dlam %.3f  %s\n', ...
  thb, dlb, H_tick(ok));

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_quadpol_curved:failed', '%d case(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

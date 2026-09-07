%TEST_STRAIN_REFERENCE Fabric orientation referenced to principal strain axes.
%
% ptt.strainRateAxes reduces a surface velocity field to its horizontal
% principal strain-rate axes, and ptt.orientationToStrain reports a fabric
% axis against them. Strain is the reference because strain is what
% orients fabric; flow direction only coincides with it where flow is fast
% and straight.
%
% The verdicts:
%
%   1. ANALYTIC CASES, WITH THE AZIMUTH CONVENTION PINNED. Uniaxial
%      extension along grid east must return an east-west extension axis
%      (azimuth 90 deg from grid north); simple shear must return axes at
%      45 deg. These fix the sign and the north-vs-east convention, which
%      is the one thing here that is easy to get wrong by 90 degrees and
%      impossible to notice afterwards.
%   2. THE 45-DEGREE MARGIN CLAIM, TESTED NOT ASSERTED. The headers of
%      both functions justify preferring strain over flow by claiming
%      that in a shear margin the principal axes sit 45 deg off the flow
%      direction. That claim is checked here against the flow reference
%      computed by ptt.orientationToFlow, so the documentation cannot
%      drift away from what the code does.
%   3. THE GRADIENT PATH AGREES WITH THE ANALYTIC ONE. A velocity raster
%      differentiated on the grid must reproduce the components fed in
%      directly, or every field result is on a different footing from
%      every synthetic one.
%   4. EXTENSION AND COMPRESSION ARE EXACT COMPLEMENTS, and the choice is
%      reported. Reading a perpendicular result as a parallel one is the
%      failure this guards.
%   5. IT ABSTAINS, AND REFUSES MIXED FRAMES. Principal directions from a
%      near-zero strain rate are noise; grid north is not true north.
%
% Run: matlab -batch "run('opr_fabric/test/test_strain_reference.m')"
clear; t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

%% ---- verdict 1: analytic cases and the azimuth convention
a = ptt.strainRateAxes(struct('exx', 1, 'eyy', -1, 'exy', 0));
b = ptt.strainRateAxes(struct('exx', 0, 'eyy', 0, 'exy', 1));
c = ptt.strainRateAxes(struct('exx', -1, 'eyy', 1, 'exy', 0));
fprintf('1. uniaxial extension along grid EAST -> az_ext %.1f deg (want 90)\n', rad2deg(a.az_ext));
fprintf('   simple shear exy>0             -> az_ext %.1f deg (want 45)\n', rad2deg(b.az_ext));
fprintf('   uniaxial extension along NORTH  -> az_ext %.1f deg (want 0)\n', rad2deg(c.az_ext));
ok1 = abs(rad2deg(a.az_ext) - 90) < 1e-6 && abs(rad2deg(b.az_ext) - 45) < 1e-6 && ...
      abs(mod(rad2deg(c.az_ext) + 90, 180) - 90) < 1e-6;
fprintf('   azimuth convention pinned (north-clockwise):           %s\n', H_tick(ok1));

%% ---- verdict 2: the 45-degree shear-margin claim
% flow along grid east, deformation pure simple shear dvx/dy = gamma
gam = 1e-3;
s2 = ptt.strainRateAxes(struct('exx', 0, 'eyy', 0, 'exy', gam/2));
flow_az = atan2(1, 0);                       % due east
th_probe = deg2rad(0);                       % a north-south fabric axis
o_str = ptt.orientationToStrain(th_probe, s2, struct('theta_frame', 'grid'));
o_flw = ptt.orientationToFlow(th_probe, flow_az, ...
  struct('theta_frame', 'grid', 'flow_frame', 'grid'));
sep = abs(rad2deg(angle(exp(2i*(s2.az_ext - flow_az)))/2));
fprintf('\n2. shear margin: principal axis %.1f deg from flow\n', sep);
fprintf('   same fabric axis reads %.1f deg vs strain, %.1f deg vs flow\n', ...
  o_str.angle_deg, o_flw.angle_deg);
ok2 = abs(sep - 45) < 1e-6 && abs(o_str.angle_deg - o_flw.angle_deg) > 30;
fprintf('   the documented 45 deg offset is real, not asserted:    %s\n', H_tick(ok2));

%% ---- verdict 3: the gradient path reproduces the analytic components
[X, Y] = meshgrid(0:100:10000, 0:100:10000);
EXX = 2e-4; EYY = -1e-4; EXY = 5e-5;
vx = EXX * X + EXY * Y;
vy = EXY * X + EYY * Y;
g = ptt.strainRateAxes(struct('vx', vx, 'vy', vy, 'dx', 100));
i0 = round(size(X,1)/2); j0 = round(size(X,2)/2);
d_ex = abs(g.exx(i0,j0) - EXX); d_ey = abs(g.eyy(i0,j0) - EYY); d_xy = abs(g.exy(i0,j0) - EXY);
ref = ptt.strainRateAxes(struct('exx', EXX, 'eyy', EYY, 'exy', EXY));
d_az = abs(rad2deg(angle(exp(2i*(g.az_ext(i0,j0) - ref.az_ext)))/2));
fprintf('\n3. gradient path vs analytic: d_exx %.2e  d_eyy %.2e  d_exy %.2e  d_az %.2e deg\n', ...
  d_ex, d_ey, d_xy, d_az);
ok3 = d_ex < 1e-12 && d_ey < 1e-12 && d_xy < 1e-12 && d_az < 1e-9;
fprintf('   rasters and components are on one footing:             %s\n', H_tick(ok3));

%% ---- verdict 4: extension and compression are exact complements
th4 = deg2rad(20);
oe = ptt.orientationToStrain(th4, ref, struct('theta_frame','grid','reference','extension'));
oc = ptt.orientationToStrain(th4, ref, struct('theta_frame','grid','reference','compression'));
fprintf('\n4. same axis: %.2f deg from extension, %.2f deg from compression (sum %.2f)\n', ...
  oe.angle_deg, oc.angle_deg, oe.angle_deg + oc.angle_deg);
ok4 = abs(oe.angle_deg + oc.angle_deg - 90) < 1e-9 && strcmp(oe.reference,'extension') ...
      && strcmp(oc.reference,'compression');
fprintf('   complements, and the choice is reported back:          %s\n', H_tick(ok4));

%% ---- verdict 5: abstention and frame safety
quiet = ptt.strainRateAxes(struct('exx', 1e-12, 'eyy', -1e-12, 'exy', 0));
o5 = ptt.orientationToStrain(deg2rad(10), quiet, ...
  struct('theta_frame','grid','e_min',1e-6));
refused = false;
try
  ptt.orientationToStrain(0, ref, struct('theta_frame','true','strain_frame','grid'));
catch e
  refused = strcmp(e.identifier, 'ptt:orientationToStrain:convergence');
end
fprintf('\n5. near-zero strain rate -> %s;  mixed frames refused: %d\n', ...
  mat2str(o5.angle_deg), refused);
ok5 = isnan(o5.angle_deg) && refused;
fprintf('   abstains on noise, refuses to mix grid and true north:  %s\n', H_tick(ok5));

fails = ~ok1 + ~ok2 + ~ok3 + ~ok4 + ~ok5;
fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_strain_reference:failed', '%d verdict(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

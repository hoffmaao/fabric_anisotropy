%TEST_ROTATIONS ptt.rotatePolarization and ptt.rotateMoments agree.
%
% ptt.rotateMoments exists so a block of traces can be rotated into a
% common geographic frame for the cost of a 4x4 multiply instead of
% re-synthesizing the full grid, which is what run_quadpol_survey.m does
% before it averages. That shortcut is only sound if its W matrix is the
% SAME linear map ptt.rotatePolarization applies to the channels - the two
% are written out independently, so a transcription slip in either would
% silently rotate the survey's geographic average by the wrong angle and
% show up only as an inflated circular spread, which is exactly the
% quantity the antenna-frame test keys on.
%
% Three properties, all of them cheap and none of them true by accident:
%   1. moments(rotate(S, psi)) == rotateMoments(moments(S), psi)
%   2. rotating by psi then -psi returns the original channels
%   3. rotating by psi + 180 deg returns the original channels, since the
%      antenna axes are axes and not directions
%
% Run: matlab -batch "run('opr_fabric/test/test_rotations.m')"
clear;
rng(13);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

Nt = 64; Nx = 40;
S = struct();
for f = {'hh','vv','hv','vh'}
  S.(f{1}) = (randn(Nt,Nx) + 1i*randn(Nt,Nx))/sqrt(2);
end
scale = max(abs([S.hh(:); S.vv(:); S.hv(:); S.vh(:)]));

fails = 0;
M = ptt.quadpolMoments(S, [1 Nx]);
for psi_deg = [0 17 45 90 123 180]
  psi = deg2rad(psi_deg);

  % 1. the two rotations describe the same map
  Md = ptt.quadpolMoments(ptt.rotatePolarization(S, psi), [1 Nx]);
  Mr = ptt.rotateMoments(M, psi);
  e_mom = max(abs(Md(:) - Mr(:))) / max(abs(Md(:)));

  % 2. round trip
  R = ptt.rotatePolarization(ptt.rotatePolarization(S, psi), -psi);
  e_rt = max([max(abs(R.hh(:) - S.hh(:))), max(abs(R.vv(:) - S.vv(:))), ...
              max(abs(R.hv(:) - S.hv(:))), max(abs(R.vh(:) - S.vh(:)))]) / scale;

  % 3. 180 deg is the identity on an axis
  H = ptt.rotatePolarization(S, psi + pi);
  G = ptt.rotatePolarization(S, psi);
  e_ax = max([max(abs(H.hh(:) - G.hh(:))), max(abs(H.vv(:) - G.vv(:))), ...
              max(abs(H.hv(:) - G.hv(:))), max(abs(H.vh(:) - G.vh(:)))]) / scale;

  ok = e_mom < 1e-10 && e_rt < 1e-12 && e_ax < 1e-12;
  fails = fails + ~ok;
  fprintf('psi %3d: moments %.2e | round trip %.2e | 180 deg %.2e  %s\n', ...
    psi_deg, e_mom, e_rt, e_ax, H_tick(ok));
end

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_rotations:failed', '%d angle(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

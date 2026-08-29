%EXTRACT_SWEEP_EGRIP Measured C(psi,z) and P(psi,z) for the EastGRIP frame.
% Same idea as extract_sweep_band.m for Ridge A frame 009, but the whole
% usable column rather than two bands: with the dlam cap fixed the frame
% still misfits (resid 0.485 pooled, 0.476 with the segmented frame
% pass), and the suspect - anisotropic reflectivity
% from the EGRIP girdle - is an EVEN cos-2psi structure in the CO-POL POWER
% aligned with theta0, which the LS nuisance basis cannot absorb. That
% hypothesis is testable by looking at Phh(psi) directly.
% repo root and <work> derived from this script's own location - see
% fabric_paths
[fabric_code, scratch] = fabric_paths();
addpath(fabric_code);
cq = load(fullfile(scratch, 'stages', 'quadpol', 'coreg_cache', ...
  'creg_20240619_01_001.mat'));
T = struct('hh', double(cq.hh), 'vv', double(cq.vv), ...
  'hv', double(cq.hv), 'vh', double(cq.vh));
z = cq.z(:);
clear cq;
dz = median(diff(z));
nr = max(3, round(10 / dz));
M = ptt.quadpolMoments(T, [nr size(T.hh, 2)]);
psi = (0:2:178) * pi/180;
A = ptt.quadpolAzimuth(M, psi);
keep = find(z >= 100 & z <= 1400);
keep = keep(1:4:end);
C = single(conj(A.chhvv(keep, :)));   % deramped, as the LS sees it
Phh = single(A.Phh(keep, :)); Pxc = single(A.Pxc(keep, :));
zk = z(keep);
save(fullfile(scratch, 'stages', 'sweep_egrip_619_001.mat'), ...
  '-v7.3', 'C', 'Phh', 'Pxc', 'zk', 'psi');
fprintf('wrote sweep_egrip_619_001.mat: %d depths x %d azimuths\n', ...
  numel(zk), numel(psi));

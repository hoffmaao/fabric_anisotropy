function T = rotatePolarization(S, psi)
%ROTATEPOLARIZATION Synthesize the quad-pol response at a rotated azimuth.
%
% T = ptt.rotatePolarization(S, psi)
%
% Rotates a measured 2x2 scattering matrix into an antenna frame turned by
% `psi` (radians, clockwise positive in the same sense as track azimuth),
%
%   T(psi) = R(psi)' * S * R(psi),    R = [cos psi, -sin psi;
%                                          sin psi,  cos psi]
%
% which is just a change of basis: the rotated antenna vectors are the
% columns of R, and T_ij is the response between them.
%
% THIS IS THE WHOLE POINT OF HAVING FOUR CHANNELS. A co-polarized pair
% (HH, VV) alone measures one azimuth and cannot be rotated - the
% off-diagonal terms it would need are exactly the ones it did not record.
% With HV and VH present the response at EVERY azimuth is recoverable from
% a single pass, which is what a quad-pol ApRES achieves mechanically by
% turning the antennas on the snow. It removes the two-azimuth requirement
% that limits ptt-based mapping to survey crossings, and it is what makes
% orientation recoverable along a single line.
%
% Inputs
%   S    struct with fields hh, vv, hv, vh, each [Nt x Nx] complex (or any
%        common size; they are combined elementwise)
%   psi  rotation angle(s) in RADIANS. Scalar, or a vector - in which case
%        each output field is [Nt x Nx x numel(psi)].
%
% Output
%   T    struct with the same four fields at the rotated azimuth.
%
% Reciprocity (S_hv == S_vh) is NOT assumed. A real system has finite
% antenna isolation and channel imbalance, so the two cross terms differ;
% keeping both and letting the algebra carry them means their difference
% stays visible downstream as a calibration diagnostic instead of being
% averaged away here.

fn = {'hh','vv','hv','vh'};
for k = 1:numel(fn)
  if ~isfield(S, fn{k})
    error('ptt:rotatePolarization:field', ...
      'S must have fields hh, vv, hv, vh; missing %s', fn{k});
  end
end
sz = size(S.hh);
for k = 2:numel(fn)
  if ~isequal(size(S.(fn{k})), sz)
    error('ptt:rotatePolarization:size', ...
      'S.%s is %s but S.hh is %s', fn{k}, mat2str(size(S.(fn{k}))), ...
      mat2str(sz));
  end
end

psi = psi(:);
np = numel(psi);
T = struct();
for k = 1:numel(fn)
  T.(fn{k}) = zeros([sz np], 'like', complex(S.hh(1)));
end

for j = 1:np
  c = cos(psi(j));
  s = sin(psi(j));
  % T = R' S R, written out so the four outputs are explicit and each can
  % be checked against the analytic single-layer result independently.
  T.hh(:,:,j) = c*c*S.hh + c*s*(S.vh + S.hv) + s*s*S.vv;
  T.vv(:,:,j) = s*s*S.hh - c*s*(S.vh + S.hv) + c*c*S.vv;
  T.hv(:,:,j) = -c*s*S.hh + c*c*S.hv - s*s*S.vh + c*s*S.vv;
  T.vh(:,:,j) = -c*s*S.hh - s*s*S.hv + c*c*S.vh + c*s*S.vv;
end

if np == 1
  for k = 1:numel(fn)
    T.(fn{k}) = reshape(T.(fn{k}), sz);
  end
end

end

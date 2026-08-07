function A = quadpolAzimuth(M, psi)
%QUADPOLAZIMUTH Synthetic azimuth sweep of the polarimetric observables.
%
% A = ptt.quadpolAzimuth(M, psi)
%
% Given the moment matrix from ptt.quadpolMoments and a vector of synthetic
% antenna azimuths `psi` [radians], returns the response that a quad-pol
% system rotated to each of those azimuths would have recorded:
%
%   A.psi   [1 x Np]        the azimuths, echoed back
%   A.Phh   [Nt x Np] real  co-polarized power  <|S_hh(psi)|^2>
%   A.Pvv   [Nt x Np] real  co-polarized power  <|S_vv(psi)|^2>
%   A.Pxc   [Nt x Np] real  cross-polarized power, (|S_hv|^2+|S_vh|^2)/2
%   A.Chhvv [Nt x Np] cplx  <S_hh(psi) conj(S_vv(psi))>, UNNORMALIZED
%   A.chhvv [Nt x Np] cplx  the same, normalized to a coherence in [0,1]
%
% The rotation acts on the scattering matrix as T = R' S R, so each element
% of T is a fixed real trig combination of (hh, vv, hv, vh):
%
%   T_hh = [ c^2,  s^2,  cs,  cs] . S
%   T_vv = [ s^2,  c^2, -cs, -cs] . S
%   T_hv = [-cs,   cs,   c^2, -s^2] . S
%   T_vh = [-cs,   cs,  -s^2,  c^2] . S
%
% and any second-order observable is then w1' * M * w2 for the appropriate
% coefficient vectors, which is why only M is needed and never the images.
%
% WHAT THE SWEEP LOOKS LIKE FOR A SINGLE BIREFRINGENT COLUMN. With
% principal axes at theta and accumulated differential phase delta(z),
%
%   Pxc(z, psi)  ~  sin^2(2(psi - theta)) * sin^2(delta(z)/2)
%
% so the cross-polarized power vanishes along the principal axes at EVERY
% depth (that is what fixes theta, modulo 90 deg) and, along any other
% azimuth, nulls in depth wherever delta = 2*pi*m (that is what fixes the
% birefringence, with no reference to theta at all). The two dependences
% factor exactly, which is the property the co-polarized pair lacks: there,
% arg<S_hh conj(S_vv)> reduces to delta*cos(2(psi-theta)) only in the
% small-delta limit.

if ndims(M) ~= 3 || size(M,2) ~= 4 || size(M,3) ~= 4
  error('ptt:quadpolAzimuth:shape', 'M must be [Nt x 4 x 4]');
end
psi = psi(:).';
Nt = size(M, 1);
Np = numel(psi);

A = struct('psi', psi);
A.Phh = zeros(Nt, Np);
A.Pvv = zeros(Nt, Np);
A.Pxc = zeros(Nt, Np);
A.Chhvv = complex(zeros(Nt, Np));
A.chhvv = complex(zeros(Nt, Np));

for j = 1:Np
  c = cos(psi(j));
  s = sin(psi(j));
  cs = c*s;
  w_hh = [c*c, s*s,  cs,  cs];
  w_vv = [s*s, c*c, -cs, -cs];
  w_hv = [-cs,  cs, c*c, -s*s];
  w_vh = [-cs,  cs, -s*s, c*c];

  A.Phh(:, j) = real(H_quad(M, w_hh, w_hh));
  A.Pvv(:, j) = real(H_quad(M, w_vv, w_vv));
  % Both cross terms are averaged rather than one chosen: for a reciprocal
  % medium they are equal, and where they are not the difference is system
  % isolation, which should not be allowed to land preferentially in
  % whichever channel happened to be picked.
  A.Pxc(:, j) = 0.5 * real(H_quad(M, w_hv, w_hv) + H_quad(M, w_vh, w_vh));
  A.Chhvv(:, j) = H_quad(M, w_hh, w_vv);
  den = sqrt(max(A.Phh(:, j), 0) .* max(A.Pvv(:, j), 0));
  A.chhvv(:, j) = A.Chhvv(:, j) ./ max(den, realmin);
end

end

function v = H_quad(M, w1, w2)
% < (w1 . S) conj(w2 . S) > = sum_kl w1_k w2_l M_kl, with real weights.
v = complex(zeros(size(M,1), 1));
for k = 1:4
  if w1(k) == 0, continue; end
  for l = 1:4
    if w2(l) == 0, continue; end
    v = v + (w1(k) * w2(l)) * M(:, k, l);
  end
end
end

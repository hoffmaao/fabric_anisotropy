function Mr = rotateMoments(M, psi)
%ROTATEMOMENTS Rotate a quad-pol moment matrix into another antenna frame.
%
% Mr = ptt.rotateMoments(M, psi)
%
% M is [Nt x 4 x 4] from ptt.quadpolMoments and psi a rotation in radians.
% The scattering matrix rotates as T = R' S R, which is a REAL linear map
% W on the channel vector (hh, vv, hv, vh); the moments are quadratic, so
% they rotate as
%
%   Mr = W * M * W'
%
% WHY THIS EXISTS. Averaging the moments over traces is only valid in a
% frame the signal of interest is STATIONARY in. Antennas fixed to a
% vehicle define a frame that turns with the heading, and a ground survey
% wanders: on Ridge A 45 of 47 frames vary by more than 5 degrees, with
% 5-95 percentile spreads reaching 29. A signal fixed in the ICE rotates
% trace to trace in that frame and averages down - and where the fabric
% axis sits near 45 degrees to the heading the cross-polarized term
% changes SIGN across the wander and cancels outright - while a signal
% fixed in the ANTENNAS adds coherently no matter what the heading does.
% Averaging in the antenna frame therefore biases the answer toward an
% antenna-fixed one all by itself, which is indistinguishable from
% instrument leakage unless the two frames are compared.
%
% Rotating each trace block into a common GEOGRAPHIC frame before
% averaging inverts that preference: the ice adds coherently and the
% instrument smears. Running both and comparing is what separates a real
% leakage term from an artifact of the averaging.
%
% Doing it on the MOMENTS rather than on the images is what makes this
% affordable: the moments are 16 numbers per range bin, so a per-block
% rotation costs a 4x4 multiply instead of touching the full grid again.

if ndims(M) ~= 3 || size(M,2) ~= 4 || size(M,3) ~= 4
  error('ptt:rotateMoments:shape', 'M must be [Nt x 4 x 4]');
end
c = cos(psi);
s = sin(psi);
cs = c*s;
W = [c*c, s*s,  cs,  cs;
     s*s, c*c, -cs, -cs;
     -cs,  cs, c*c, -s*s;
     -cs,  cs, -s*s, c*c];

Nt = size(M, 1);
Mr = complex(zeros(Nt, 4, 4));
for k = 1:4
  for l = 1:4
    acc = complex(zeros(Nt, 1));
    for a = 1:4
      if W(k,a) == 0, continue; end
      for b = 1:4
        if W(l,b) == 0, continue; end
        acc = acc + (W(k,a) * W(l,b)) * M(:, a, b);
      end
    end
    Mr(:, k, l) = acc;
  end
end

end

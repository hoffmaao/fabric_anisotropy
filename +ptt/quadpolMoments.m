function M = quadpolMoments(S, looks)
%QUADPOLMOMENTS Trace-averaged second-order moments of the scattering matrix.
%
% M = ptt.quadpolMoments(S, looks)
%
% Returns M [Nt x 4 x 4], the Hermitian moment matrix
%
%   M_kl(t) = < S_k(t,:) conj(S_l(t,:)) >     k,l over (hh, vv, hv, vh)
%
% averaged over traces and, if `looks` has a range term, over range too.
%
% WHY THIS EXISTS RATHER THAN ROTATING THE IMAGES. Every polarimetric
% observable used here - co-pol power, cross-pol power, the HH-VV
% coherence - is QUADRATIC in the scattering matrix, and rotation is
% linear. So the response at any synthetic azimuth is a fixed trig
% combination of these 16 numbers per range bin, and the azimuth sweep
% costs nothing once they are formed. Rotating the images themselves would
% be ~5e9 complex operations per frame for a 180-point sweep on a
% 13500 x 2000 grid, and would recompute the same averages 180 times.
%
% The full 4x4 is kept, not just the (hh, vv) block: the cross terms
% <hh conj(hv)> and the hv/vh difference are what carry the orientation and
% the system-isolation diagnostic respectively.
%
% Inputs
%   S      struct with complex fields hh, vv, hv, vh, each [Nt x Nx]
%   looks  [nr na] boxcar looks in (range, trace). Default [1 Nx]. Only
%          the range term is used: the trace direction is always averaged
%          whole, because a single Nt x 4 x 4 moment matrix is the output
%          and a per-look trace axis would have nowhere to go. Callers pass
%          the trace count in looks(2) for documentation; pass a smaller
%          block by slicing S instead, which is what run_quadpol_survey.m
%          does for its along-track blocks.
%
% Output
%   M      [Nt x 4 x 4] complex; M(:,k,l) = conj(M(:,l,k))

fn = {'hh','vv','hv','vh'};
for k = 1:numel(fn)
  if ~isfield(S, fn{k})
    error('ptt:quadpolMoments:field', 'S is missing field %s', fn{k});
  end
end
[Nt, Nx] = size(S.hh);
if nargin < 2 || isempty(looks), looks = [1 Nx]; end
nr = max(1, round(looks(1)));

M = complex(zeros(Nt, 4, 4));
kr = ones(nr, 1) / nr;
for k = 1:4
  for l = k:4
    p = S.(fn{k}) .* conj(S.(fn{l}));
    % Average over traces first (that is the look direction with the most
    % samples), then smooth in range if asked.
    v = mean(p, 2, 'omitnan');
    if nr > 1
      v = conv(v, kr, 'same');
    end
    M(:, k, l) = v;
    if l ~= k
      M(:, l, k) = conj(v);
    end
  end
end

end

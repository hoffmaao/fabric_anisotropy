function out = goldsteinFilter(C, alpha, win, step)
%GOLDSTEINFILTER Goldstein-Werner adaptive spectral filter of an interferogram.
%
% out = ptt.goldsteinFilter(C, alpha, win, step)
%
% Goldstein, R. M. and Werner, C. L. (1998), "Radar interferogram filtering
% for geophysical applications", Geophys. Res. Lett. 25, 4035-4038.
%
% Overlapping windows; each window's spectrum S is reweighted by its own
% smoothed magnitude raised to `alpha`, then overlap-added under a
% separable Hann taper. Where the fringe pattern is locally coherent its
% energy concentrates in a few spectral bins, so the reweighting amplifies
% those and suppresses the rest; where it is noise-dominated the spectrum
% is flat and the reweighting is close to a no-op. The point is not
% cosmetic: unwrapping cost is set by RESIDUES, and residues are created by
% noise-driven phase excursions between adjacent pixels, so filtering the
% complex interferogram before unwrapping is what lets SNAPHU find a
% solution that is not dominated by branch cuts.
%
% Applied to the COMPLEX interferogram, never to the phase: smoothing a
% wrapped phase averages across the +-pi branch cut and invents fringes
% that were never in the data.
%
% Inputs
%   C     complex interferogram, [Nt x Nx]
%   alpha filter strength in [0,1] (default 0.8). 0 is a no-op; 1 is the
%         strongest. Above ~0.8 the filter starts inventing fringe
%         continuity across genuinely incoherent gaps.
%   win   window side in samples (default 64)
%   step  window stride in samples (default 16, i.e. 4x overlap)
%
% Output
%   out   filtered complex interferogram, same size and class as C
%
% Defaults match scripts/prototypes/branch_cuts.py, which is where the
% alpha/win/step choice was tested against residue counts on the Thwaites
% margin frame, so the production path and that diagnostic cannot drift.

if nargin < 2 || isempty(alpha), alpha = 0.8; end
if nargin < 3 || isempty(win),   win = 64;    end
if nargin < 4 || isempty(step),  step = 16;   end

validateattributes(C, {'single','double'}, {'2d'}, mfilename, 'C');
if alpha < 0 || alpha > 1
  error('ptt:goldsteinFilter:alpha', 'alpha must be in [0,1], got %g', alpha);
end

[Nt, Nx] = size(C);
if Nt < win || Nx < win
  % Nothing to do: a grid smaller than one window has no overlap structure
  % to add over, and silently returning zeros here would look like a
  % coherence collapse several stages downstream.
  warning('ptt:goldsteinFilter:tooSmall', ...
    ['grid is %dx%d, smaller than the %d-sample window; returning the ' ...
     'input unfiltered'], Nt, Nx, win);
  out = C;
  return;
end

was_single = isa(C, 'single');
C = double(C);
C(~isfinite(C)) = 0;    % SNAPHU gets a finite grid either way

t1 = hann(win);
taper = (t1 * t1.') + 1e-6;   % +eps so the overlap weight is never exactly 0
k3 = ones(3) / 9;

out = complex(zeros(Nt, Nx));
wsum = zeros(Nt, Nx);

r0s = 1:step:(Nt - win + 1);
c0s = 1:step:(Nx - win + 1);
% Always include the last full window: with a stride that does not divide
% the grid, the final (win-step) rows/columns would otherwise be covered
% by fewer windows than the interior and come out under-filtered.
if r0s(end) ~= Nt - win + 1, r0s(end+1) = Nt - win + 1; end
if c0s(end) ~= Nx - win + 1, c0s(end+1) = Nx - win + 1; end

for r0 = r0s
  ri = r0:r0+win-1;
  for c0 = c0s
    ci = c0:c0+win-1;
    S = fft2(C(ri, ci));
    mag = conv2(abs(S), k3, 'same');
    mx = max(mag(:));
    if mx <= 0, continue; end
    f = ifft2(S .* (mag / mx).^alpha);
    out(ri, ci) = out(ri, ci) + f .* taper;
    wsum(ri, ci) = wsum(ri, ci) + taper;
  end
end

ok = wsum > 0;
out(ok) = out(ok) ./ wsum(ok);
if was_single
  out = single(out);
end

end

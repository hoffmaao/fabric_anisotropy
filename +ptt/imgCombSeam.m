function seam = imgCombSeam(map, opts)
%IMGCOMBSEAM Fast-time mask of the waveform image-combination seams.
%   seam = IMGCOMBSEAM(map, opts) returns an (Nt x Nx) logical mask that is
%   true on the fast-time samples straddling each boundary where OPR
%   stitched two waveform images into the single echogram the polarimetric
%   product was formed from. Across the blend the trace is a weighted mix
%   of a short and a long pulse with different bandwidths, gains and
%   system delays, so the interferometric phase (and the delta-k group
%   delay derived from it) is not a property of the ice there. Masking
%   those samples keeps the seam out of the block averages.
%
%   The boundaries are read from the product's own metadata rather than
%   hardcoded, because they differ per frame: the accum3 grids combine at
%   0.9 us with a 0.1 us blend, the deeper settings at 8 us with a 1 us
%   blend, and single-image frames declare no boundary at all.
%
%   map fields used: Time (Nt x 1, uniform), Surface (1 x Nx), and
%   img_comb - the param.array.img_comb vector carried through by
%   fabric_task, laid out as [t_comb mult window] per boundary
%   (3*(Nimg-1) elements). An empty or absent img_comb masks nothing.
%
%   opts fields (optional): seam_mask_win (1), guard width as a multiple
%   of the blend window each boundary declares (see below).
%
%   The masked band is ASYMMETRIC. OPR crossfades over [t_comb, t_comb +
%   window], so the samples before the boundary are clean and the
%   contamination runs forward. On 20250108_02_009 (boundary 0.9 us,
%   window 0.1 us) the mean coherence holds its 0.968 baseline right up to
%   0.89 us, steps down to 0.857 by 0.92 us, and only recovers by ~1.1 us,
%   while the mean power kink reverses sign at 0.99 us. The band is
%   therefore [t_comb - k*window, t_comb + (1+k)*window] for
%   k = seam_mask_win, which at the default k = 1 spans 0.8-1.1 us there
%   and 7-10 us for the 8 us / 1 us settings.
%
%   OPR places the boundary at max(Surface*mult, t_comb), so a finite mult
%   makes it track the surface and the mask becomes column-dependent; the
%   common mult = -Inf pins it to the fixed traveltime t_comb.

Nt = numel(map.Time);
Nx = numel(map.Surface);
seam = false(Nt, Nx);

if ~isfield(map,'img_comb') || isempty(map.img_comb)
  return;
end
ic = map.img_comb(:).';
if mod(numel(ic),3) ~= 0
  warning('ptt:imgCombSeam:layout', ...
    'img_comb has %d elements; expected a multiple of 3 ([t_comb mult window] per boundary). No seam masked.', ...
    numel(ic));
  return;
end

if ~isfield(opts,'seam_mask_win') || isempty(opts.seam_mask_win)
  opts.seam_mask_win = 1;
end
if ~isfinite(opts.seam_mask_win) || opts.seam_mask_win <= 0
  return;
end

t = map.Time(:);
dt = t(2) - t(1);
surf = map.Surface(:).';

for b = 1:numel(ic)/3
  t_comb = ic(3*b-2);
  mult   = ic(3*b-1);
  win    = ic(3*b);
  if ~isfinite(t_comb)
    continue;
  end
  if isfinite(mult)
    t_c = max(surf*mult, t_comb);
    % Columns with no surface pick fall back to the fixed traveltime
    t_c(~isfinite(t_c)) = t_comb;
  else
    t_c = repmat(t_comb, 1, Nx);
  end
  % A zero-width blend is still a discontinuity: keep half a sample either
  % side so the boundary sample itself is always removed.
  guard = max(opts.seam_mask_win*abs(win), dt/2);
  seam = seam | (t >= t_c - guard & t <= t_c + abs(win) + guard);
end

end

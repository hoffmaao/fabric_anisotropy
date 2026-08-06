function seam = imgCombSeam(map, opts)
%IMGCOMBSEAM Fast-time mask of the waveform image-combination seams.
%   seam = IMGCOMBSEAM(map, opts) returns a logical mask that is true on
%   the fast-time samples straddling each boundary where OPR stitched two
%   waveform images into the single echogram the polarimetric product was
%   formed from. Across the blend the trace is a weighted mix of a short
%   and a long pulse with different bandwidths, gains and system delays, so
%   the interferometric phase (and the delta-k group delay derived from it)
%   is not a property of the ice there. Masking those samples keeps the
%   seam out of the block averages.
%
%   The boundaries are read from the product's own metadata rather than
%   hardcoded, because they differ per frame: the accum3 grids combine
%   0.9 us after the surface return with a 0.1 us image-1 guard, the deeper
%   settings 8 us after with a 1 us guard, and single-image frames declare
%   no boundary at all.
%
%   Each boundary is placed where OPR's img_combine.m actually combines:
%     t_c = max(min(max(0,Surface)*img_comb_mult, Surface + ic(1)), ic(2))
%   with ic = img_comb(3b-2:3b) for boundary b, where ic(1) is the combine
%   time AFTER the surface return, ic(2) the minimum absolute combine time
%   (-Inf on every in-scope frame), and img_comb_mult the separate
%   param.array surface multiplier, defaulted to Inf when the product does
%   not carry it so the min() term drops out. ic(3) is the guard time at
%   the end of image 1 - roughly the image-1 pulse duration - and scales
%   the masked band (see below).
%
%   map fields used: Time (Nt x 1, uniform), Surface (1 x Nx), img_comb
%   ([t_after t_min guard] per boundary, 3*(Nimg-1) elements) and
%   optionally img_comb_mult, both carried through from param.array by
%   fabric_task. An empty or absent img_comb masks nothing.
%
%   Because ic(1) is surface-relative the boundary tracks the surface pick
%   per column, so the mask is (Nt x Nx) in general; it collapses to
%   (Nt x 1) when the surface picks are all identical (or all missing) and
%   every column shares the same band. Callers must handle both, e.g.
%   coh_mask(seam,:) vs coh_mask(seam). A column with no finite surface
%   pick resolves to t_c = ic(2), which is what OPR's min/max arithmetic
%   computes there; with the usual ic(2) = -Inf that column masks nothing,
%   matching OPR clamping the combine to the top of the record where no
%   in-ice samples are blended.
%
%   opts fields (optional): seam_mask_win (1), mask half-width as a
%   multiple of the ic(3) guard each boundary declares (see below); 0
%   keeps only the half-sample guard, and seam_mask_en
%   (ptt.surfaceReference) is what turns the mask off. seam_mask_deltak
%   (false), set by ptt.deltakTraveltime to widen the band by the
%   analysis-cell reach.
%
%   The masked band is ASYMMETRIC. OPR blends forward from the boundary,
%   so the samples before it are clean and the contamination runs deeper.
%   On 20250108_02_009 (boundary 0.9 us, guard 0.1 us) the mean coherence
%   holds its 0.968 baseline right up to 0.89 us, steps down to 0.857 by
%   0.92 us, and only recovers by ~1.1 us, while the mean power kink
%   reverses sign at 0.99 us. The band is therefore
%   [t_c - k*ic(3), t_c + (1+k)*ic(3)] for k = seam_mask_win.
%
%   OPR's own crossfade is img_comb_bins wide (default 1 bin), far
%   narrower than ic(3). The ic(3)-scaled band here is NOT a reading of
%   the crossfade width: it is deliberate masking of the matched-filter
%   settling zone around the splice, which the coherence data above show
%   is contaminated over roughly the image-1 pulse duration either side,
%   not over one bin.
%
%   Delta-k needs a wider band. ptt.deltakTraveltime resolves its ladder
%   integers from a tau_A smoothed with a dk.smooth boxcar on the analysis
%   cell grid and then interpolates the cell values back to the full grid,
%   so a cell up to (dk.smooth-1)/2 + 1 cells away from the seam inherits
%   seam-contaminated tau_A - and a contaminated tau_A does not bias the
%   result slightly, it moves it by a whole 1/dfQ step. The band grows by
%   that reach (300 ns at the defaults) so that every cell whose smoothing
%   window touched the seam is dropped, i.e. [t_c - k*ic(3) - reach,
%   t_c + (1+k)*ic(3) + reach]. The phase and coregistration estimators
%   do no cell smoothing and keep the un-widened band.
%
%   The two bands cost very different amounts of record. Measured on the
%   6601-bin accum3 frame 20250108_02_009 at the default k = 1, per
%   dtau_source:
%     phase / coreg   0.9 us / 0.1 us:  0.8-1.1 us,    90 bins  (1.4%)
%                     8 us / 1 us:      7-10 us,      900 bins (13.6%)
%     deltak          0.9 us / 0.1 us:  0.5-1.4 us,   270 bins  (4.1%)
%                     8 us / 1 us:      6.7-10.3 us, 1080 bins (16.4%)
%   Budget delta-k coverage against the delta-k rows, not the phase ones:
%   the reach triples the cost on the 0.9 us settings.

Nt = numel(map.Time);

if ~isfield(map,'img_comb') || isempty(map.img_comb)
  seam = false(Nt,1);
  return;
end
ic = map.img_comb(:).';
if mod(numel(ic),3) ~= 0
  warning('ptt:imgCombSeam:layout', ...
    'img_comb has %d elements; expected a multiple of 3 ([t_after t_min guard] per boundary). No seam masked.', ...
    numel(ic));
  seam = false(Nt,1);
  return;
end

if ~isfield(opts,'seam_mask_win') || isempty(opts.seam_mask_win)
  opts.seam_mask_win = 1;
end
if ~isscalar(opts.seam_mask_win) || ~isfinite(opts.seam_mask_win) ...
    || opts.seam_mask_win < 0
  warning('ptt:imgCombSeam:guard', ...
    'seam_mask_win must be a finite non-negative multiple of the blend window; no seam masked. Set seam_mask_en = false to disable the mask deliberately.');
  seam = false(Nt,1);
  return;
end

% The separate surface multiplier img_combine.m applies; its input check
% defaults it to Inf, which makes the min() term drop out below.
mult = Inf;
if isfield(map,'img_comb_mult') && ~isempty(map.img_comb_mult)
  mult = map.img_comb_mult;
end

t = map.Time(:);
dt = t(2) - t(1);

% Reach of the delta-k cell smoothing plus the cell-to-grid interpolation,
% in the cell size ptt.deltakTraveltime actually uses (round to samples).
reach = 0;
if isfield(opts,'seam_mask_deltak') && ~isempty(opts.seam_mask_deltak) ...
    && opts.seam_mask_deltak
  dk = ptt.deltakDefaults(opts);
  cell_dt = max(1, round(dk.cell_twtt/dt))*dt;
  reach = (max(0,(dk.smooth-1)/2) + 1)*cell_dt;
end

% The boundary is surface-relative, so it varies across columns only
% through the surface pick: when every pick is identical (or every pick is
% missing) all columns share one band and the mask stays (Nt x 1);
% otherwise it needs the full (Nt x Nx).
surf = map.Surface(:).';
finite_surf = surf(isfinite(surf));
col_dep = numel(unique(finite_surf)) > 1 ...
  || (~isempty(finite_surf) && numel(finite_surf) < numel(surf));
if col_dep
  seam = false(Nt, numel(surf));
else
  if isempty(finite_surf)
    surf = NaN;   % every column resolves to t_c = ic(2), as in OPR
  else
    surf = finite_surf(1);
  end
  seam = false(Nt, 1);
end

for b = 1:numel(ic)/3
  t_after = ic(3*b-2);
  t_min   = ic(3*b-1);
  win     = ic(3*b);
  if ~isfinite(win)
    % A malformed guard must still remove the boundary itself, not fall
    % through to a comparison against NaN that masks nothing
    warning('ptt:imgCombSeam:window', ...
      'Boundary %d declares a non-finite image-1 guard time; masking the boundary sample only.', b);
    win = 0;
  end
  % OPR's boundary, elementwise over columns. MATLAB's min/max drop a NaN
  % side, so a column with no surface pick resolves to t_min, and a
  % non-finite t_c masks nothing through the comparisons below.
  t_c = max(min(max(0,surf)*mult, surf + t_after), t_min);
  % A zero-width blend is still a discontinuity: keep half a sample either
  % side so the boundary sample itself is always removed.
  guard = max(opts.seam_mask_win*abs(win), dt/2) + reach;
  seam = seam | (t >= t_c - guard & t <= t_c + abs(win) + guard);
end

end

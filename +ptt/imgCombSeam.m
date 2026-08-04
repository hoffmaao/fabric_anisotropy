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
%   The mask is (Nt x 1) when every boundary sits at a fixed traveltime -
%   the common mult = -Inf case, where the band is the same in every column
%   - and (Nt x Nx) only when a boundary tracks the surface (finite mult).
%   Callers must handle both, e.g. coh_mask(seam,:) vs coh_mask(seam).
%
%   The boundaries are read from the product's own metadata rather than
%   hardcoded, because they differ per frame: the accum3 grids combine at
%   0.9 us with a 0.1 us blend, the deeper settings at 8 us with a 1 us
%   blend, and single-image frames declare no boundary at all.
%
%   map fields used: Time (Nt x 1, uniform), Surface (1 x Nx, read only
%   when a boundary tracks the surface), and img_comb - the
%   param.array.img_comb vector carried through by fabric_task, laid out as
%   [t_comb mult window] per boundary (3*(Nimg-1) elements). An empty or
%   absent img_comb masks nothing.
%
%   opts fields (optional): seam_mask_win (1), guard width as a multiple
%   of the blend window each boundary declares (see below); 0 keeps only
%   the half-sample guard, and seam_mask_en (ptt.surfaceReference) is what
%   turns the mask off. seam_mask_deltak (false), set by
%   ptt.deltakTraveltime to widen the band by the analysis-cell reach.
%
%   The masked band is ASYMMETRIC. OPR crossfades over [t_comb, t_comb +
%   window], so the samples before the boundary are clean and the
%   contamination runs forward. On 20250108_02_009 (boundary 0.9 us,
%   window 0.1 us) the mean coherence holds its 0.968 baseline right up to
%   0.89 us, steps down to 0.857 by 0.92 us, and only recovers by ~1.1 us,
%   while the mean power kink reverses sign at 0.99 us. The band is
%   therefore [t_comb - k*window, t_comb + (1+k)*window] for
%   k = seam_mask_win.
%
%   Delta-k needs a wider band. ptt.deltakTraveltime resolves its ladder
%   integers from a tau_A smoothed with a dk.smooth boxcar on the analysis
%   cell grid and then interpolates the cell values back to the full grid,
%   so a cell up to (dk.smooth-1)/2 + 1 cells away from the seam inherits
%   seam-contaminated tau_A - and a contaminated tau_A does not bias the
%   result slightly, it moves it by a whole 1/dfQ step. The band grows by
%   that reach (300 ns at the defaults) so that every cell whose smoothing
%   window touched the seam is dropped, i.e. [t_comb - k*window - reach,
%   t_comb + (1+k)*window + reach]. The phase and coregistration estimators
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
%
%   OPR places the boundary at max(Surface*mult, t_comb), so a finite mult
%   makes it track the surface and the mask becomes column-dependent; the
%   common mult = -Inf pins it to the fixed traveltime t_comb.

Nt = numel(map.Time);

if ~isfield(map,'img_comb') || isempty(map.img_comb)
  seam = false(Nt,1);
  return;
end
ic = map.img_comb(:).';
if mod(numel(ic),3) ~= 0
  warning('ptt:imgCombSeam:layout', ...
    'img_comb has %d elements; expected a multiple of 3 ([t_comb mult window] per boundary). No seam masked.', ...
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

% A boundary pinned to a fixed traveltime masks the same rows in every
% column; only one that tracks the surface (finite mult) needs the full
% (Nt x Nx) mask.
t_combs = ic(1:3:end);
mults = ic(2:3:end);
col_dep = any(isfinite(t_combs) & isfinite(mults));
if col_dep
  Nx = numel(map.Surface);
  seam = false(Nt, Nx);
  surf = map.Surface(:).';
else
  seam = false(Nt, 1);
end

for b = 1:numel(ic)/3
  t_comb = t_combs(b);
  mult   = mults(b);
  win    = ic(3*b);
  if ~isfinite(t_comb)
    continue;
  end
  if ~isfinite(win)
    % A malformed window must still remove the boundary itself, not fall
    % through to a comparison against NaN that masks nothing
    warning('ptt:imgCombSeam:window', ...
      'Boundary %d declares a non-finite blend window; masking the boundary sample only.', b);
    win = 0;
  end
  if isfinite(mult)
    t_c = max(surf*mult, t_comb);
    % Columns with no surface pick fall back to the fixed traveltime
    t_c(~isfinite(t_c)) = t_comb;
  else
    t_c = t_comb;
  end
  % A zero-width blend is still a discontinuity: keep half a sample either
  % side so the boundary sample itself is always removed.
  guard = max(opts.seam_mask_win*abs(win), dt/2) + reach;
  seam = seam | (t >= t_c - guard & t <= t_c + abs(win) + guard);
end

end

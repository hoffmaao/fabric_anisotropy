function [coh_mask, ref_bin, surf_valid] = surfaceReference(map, opts)
%SURFACEREFERENCE Coherence mask and surface-referencing bins.
%   [coh_mask, ref_bin, surf_valid] = SURFACEREFERENCE(map, opts) is the
%   referencing convention shared by every dtau estimator in +ptt, so they
%   cannot drift apart on it:
%     surf_valid  columns with a finite surface pick (1 x Nx)
%     coh_mask    map.coherence >= opts.coherence_threshold, with the
%                 columns lacking a surface pick excluded and the
%                 waveform-combine seams removed (Nt x Nx)
%     ref_bin     fast-time bin opts.ref_twtt_offset below the surface,
%                 clamped into [1 Nt]; 1 where the pick is invalid (1 x Nx)
%
%   map fields used: Time (Nt x 1, uniform), Surface (1 x Nx), coherence,
%   and optionally img_comb (see ptt.imgCombSeam).
%   opts fields (optional): coherence_threshold (0.5),
%   ref_twtt_offset (50e-9 s), seam_mask_en (true), seam_mask_win (1).

if ~isfield(opts,'coherence_threshold') || isempty(opts.coherence_threshold)
  opts.coherence_threshold = 0.5;
end
if ~isfield(opts,'ref_twtt_offset') || isempty(opts.ref_twtt_offset)
  opts.ref_twtt_offset = 50e-9;
end
if ~isfield(opts,'seam_mask_en') || isempty(opts.seam_mask_en)
  opts.seam_mask_en = true;
end

Nt = numel(map.Time);
dt = map.Time(2) - map.Time(1);

surf_valid = isfinite(map.Surface(:).');
coh_mask = map.coherence >= opts.coherence_threshold;
coh_mask(:,~surf_valid) = false;

% The waveform-combine seams carry a blend of two pulses rather than an
% ice property; drop them before anything averages over fast time.
if opts.seam_mask_en
  seam = ptt.imgCombSeam(map, opts);
  n_seam = nnz(seam(:,surf_valid));
  if n_seam > 0
    coh_mask(seam) = false;
    fprintf('Seam mask: %d of %d bins removed at the waveform-combine boundaries\n', ...
      nnz(any(seam,2)), Nt);
  end
end

ref_bin = round((map.Surface(:).' + opts.ref_twtt_offset - map.Time(1))/dt) + 1;
ref_bin(~surf_valid) = 1;
ref_bin = min(max(ref_bin,1),Nt);

end

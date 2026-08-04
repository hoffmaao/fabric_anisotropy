function [coh_mask, ref_bin, surf_valid] = surfaceReference(map, opts)
%SURFACEREFERENCE Coherence mask and surface-referencing bins.
%   [coh_mask, ref_bin, surf_valid] = SURFACEREFERENCE(map, opts) is the
%   referencing convention shared by every dtau estimator in +ptt, so they
%   cannot drift apart on it:
%     surf_valid  columns with a finite surface pick (1 x Nx)
%     coh_mask    map.coherence >= opts.coherence_threshold, with the
%                 columns lacking a surface pick excluded (Nt x Nx)
%     ref_bin     fast-time bin opts.ref_twtt_offset below the surface,
%                 clamped into [1 Nt]; 1 where the pick is invalid (1 x Nx)
%
%   map fields used: Time (Nt x 1, uniform), Surface (1 x Nx), coherence.
%   opts fields (optional): coherence_threshold (0.5),
%   ref_twtt_offset (50e-9 s).

if ~isfield(opts,'coherence_threshold') || isempty(opts.coherence_threshold)
  opts.coherence_threshold = 0.5;
end
if ~isfield(opts,'ref_twtt_offset') || isempty(opts.ref_twtt_offset)
  opts.ref_twtt_offset = 50e-9;
end

Nt = numel(map.Time);
dt = map.Time(2) - map.Time(1);

surf_valid = isfinite(map.Surface(:).');
coh_mask = map.coherence >= opts.coherence_threshold;
coh_mask(:,~surf_valid) = false;

ref_bin = round((map.Surface(:).' + opts.ref_twtt_offset - map.Time(1))/dt) + 1;
ref_bin(~surf_valid) = 1;
ref_bin = min(max(ref_bin,1),Nt);

end

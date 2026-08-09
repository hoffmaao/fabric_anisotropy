function [dtau, info] = blendTraveltime(map, opts)
%BLENDTRAVELTIME Polarimetric traveltime differences from phase + offsets.
%   [dtau, info] = BLENDTRAVELTIME(map, opts) forms the two-way traveltime
%   difference dtau(fast-time, along-track) [s] between two polarization
%   channels by blending an interferogram phase (precise, ambiguous by
%   1/fc fringes) with image coregistration offsets (coarse, unambiguous):
%     1. both estimators are referenced to zero just below the surface
%        return (no birefringence above the surface; removes channel
%        timing/phase biases and any phase-unwrapping constant),
%     2. the sign of the phase convention is estimated by coherence-
%        weighted regression against the coregistration offsets (unless
%        forced via opts.phase_sign),
%     3. for wrapped phase, integer fringe offsets are resolved per pixel
%        from the coregistration. (For unwrapped phase, resolve one fringe
%        constant per along-track block later: see ptt.blockAverage.)
%
%   map is a struct describing the input images (Nt x Nx):
%     Time                fast-time axis [s], uniform spacing (Nt x 1)
%     Surface             surface return twtt per column [s] (1 x Nx)
%     fc                  center frequency [Hz]
%     phase               interferogram phase [rad], sec relative to ref
%     phase_is_unwrapped  true if phase is unwrapped (e.g. SNAPHU)
%     coherence           interferometric coherence magnitude in [0, 1]
%     row_offset          coregistration row offsets [bins], positive when
%                         the sec channel arrives later; [] if unavailable
%
%   opts fields (all optional; same names as the OPR fabric worksheet):
%     coherence_threshold (0.5), ref_twtt_offset (50e-9 s),
%     phase_sign (0 = auto), blend_coreg_en (true),
%     dtau_source ('phase' default, or 'coreg' to take dtau from the
%     coregistration offsets alone - for detected-power products where no
%     meaningful interferogram phase exists; sub-bin precision only)
%
%   info: phase_sign, coh_mask, dtau_coreg (referenced; [] if none),
%   phase_is_unwrapped, ref_bin.

% coherence_threshold and ref_twtt_offset are defaulted by
% ptt.surfaceReference, which owns the referencing convention
if ~isfield(opts,'phase_sign') || isempty(opts.phase_sign)
  opts.phase_sign = 0;
end
if ~isfield(opts,'blend_coreg_en') || isempty(opts.blend_coreg_en)
  opts.blend_coreg_en = true;
end
if ~isfield(opts,'dtau_source') || isempty(opts.dtau_source)
  opts.dtau_source = 'phase';
end
if ~isfield(opts,'ref_band_twtt') || isempty(opts.ref_band_twtt)
  opts.ref_band_twtt = 100e-9;
end

Nt = numel(map.Time);
Nx = numel(map.Surface);
dt = map.Time(2) - map.Time(1);

[coh_mask, ref_bin] = ptt.surfaceReference(map, opts);

has_coreg = isfield(map,'row_offset') && ~isempty(map.row_offset);
if has_coreg
  dtau_coreg = map.row_offset * dt;
else
  dtau_coreg = [];
end

% Surface referencing (must precede the sign regression, which the channel
% biases would otherwise contaminate). The reference is a ROBUST BAND MEAN
% over [ref_twtt_offset, ref_twtt_offset + ref_band_twtt] below the
% surface rather than the single bin at ref_twtt_offset: any error in a
% single-bin reference (surface sidelobes, speckle, a coregistration
% step) shifts every node below by a constant, and the joint inversion
% can only absorb a constant into its shallowest interval's dlam - 0.3 ns
% there manufactures dlam ~ 0.06, which is the spurious near-surface
% fabric the fringe check caught. Averaging over the coherent band shrinks
% that error by ~1/sqrt(bins); the residual is carried explicitly by the
% inversion's reference-offset nuisance. Falls back to the single bin
% where the band holds fewer than 3 coherent samples.
ref_idx = ref_bin + (0:Nx-1)*Nt;
t_rel_mat = map.Time(:) - map.Surface(:).';
band = t_rel_mat >= opts.ref_twtt_offset & ...
  t_rel_mat <= opts.ref_twtt_offset + opts.ref_band_twtt & coh_mask;
if has_coreg
  nbc = sum(band & isfinite(dtau_coreg), 1);
  rv = sum(dtau_coreg .* band, 1, 'omitnan') ./ max(nbc, 1);
  rv1 = dtau_coreg(ref_idx);
  bad = nbc < 3 | ~isfinite(rv);
  rv(bad) = rv1(bad);
  dtau_coreg = dtau_coreg - repmat(rv, [Nt 1]);
end

% Coregistration-only mode: dtau is the referenced row offsets; no phase
% sign or fringe logic applies, and map.phase may be empty (detected-power
% products have none, and the per-stage runner's coreg leg passes none)
if strcmp(opts.dtau_source,'coreg')
  assert(has_coreg, 'dtau_source=''coreg'' requires map.row_offset.');
  dtau = dtau_coreg;
  info.phase_sign = 1;
  info.coh_mask = coh_mask;
  info.dtau_coreg = []; % suppress fringe correction in ptt.blockAverage
  info.phase_is_unwrapped = true;
  info.ref_bin = ref_bin;
  return;
end

dtau_phase = map.phase / (2*pi*map.fc);
nbp = sum(band & isfinite(map.phase), 1);
if map.phase_is_unwrapped
  % arithmetic band mean on the unwrapped field
  pv = sum(dtau_phase .* band, 1, 'omitnan') ./ max(nbp, 1);
else
  % circular band mean on the wrapped field: the band is shallow, so the
  % true phase is nearly constant across it and the phasor mean is safe
  % where a plain mean would shred across a wrap
  pv = angle(sum(exp(1i*map.phase) .* band, 1, 'omitnan')) / (2*pi*map.fc);
end
pv1 = dtau_phase(ref_idx);
badp = nbp < 3 | ~isfinite(pv);
pv(badp) = pv1(badp);
dtau_phase = dtau_phase - repmat(pv, [Nt 1]);

% Phase sign: dtau = s*phase/(2*pi*fc)
if opts.phase_sign ~= 0
  phase_sign = sign(opts.phase_sign);
elseif has_coreg
  if map.phase_is_unwrapped
    dtau_coreg_cmp = dtau_coreg;
  else
    % Wrapped phase only resolves offsets within +/-1/(2*fc); rewrap the
    % coreg prediction to the same ambiguous interval before regressing
    dtau_coreg_cmp = dtau_coreg - round(dtau_coreg*map.fc)/map.fc;
  end
  sig = coh_mask & isfinite(dtau_coreg) & abs(dtau_coreg) > dt/2;
  if nnz(sig) > 100
    phase_sign = sign(sum(map.coherence(sig).^2 .* dtau_phase(sig) .* dtau_coreg_cmp(sig)));
    if phase_sign == 0, phase_sign = -1; end
  else
    warning('ptt:blendTraveltime:sign', ...
      'Too few significant coregistration offsets to estimate the phase sign; using -1 (matched-filter convention). Set opts.phase_sign to override.');
    phase_sign = -1;
  end
else
  warning('ptt:blendTraveltime:sign', ...
    'No coregistration offsets to estimate the phase sign from; using -1 (matched-filter convention). Set opts.phase_sign to override.');
  phase_sign = -1;
end
dtau_phase = phase_sign * dtau_phase;

% Wrapped phase: resolve integer fringes per pixel from the coregistration
% (coregistration noise then limits accuracy). Unwrapped phase: leave the
% per-block fringe constant to ptt.blockAverage.
dtau = dtau_phase;
if opts.blend_coreg_en && has_coreg && ~map.phase_is_unwrapped
  fringe = round((dtau_coreg - dtau_phase)*map.fc);
  fringe(~isfinite(fringe)) = 0;
  dtau = dtau_phase + fringe/map.fc;
end

info.phase_sign = phase_sign;
info.coh_mask = coh_mask;
info.dtau_coreg = dtau_coreg;
info.phase_is_unwrapped = map.phase_is_unwrapped;
info.ref_bin = ref_bin;

end

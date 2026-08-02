function inv = invertBlocks(blk, map, par, opts)
%INVERTBLOCKS Fabric inversion of block-averaged dtau.
%   inv = INVERTBLOCKS(blk, map, par, opts) inverts each along-track block
%   (from ptt.blockAverage) for the piecewise-constant horizontal fabric
%   contrast dlam = lam_x - lam_y over opts.num_intervals depth intervals,
%   using the solver selected by opts.inversion (exact per-interval layer
%   stripping via ptt.invertHorizontalFabric, or the smoothness-regularized
%   joint solve via ptt.invertHorizontalFabricJoint) with the firn-ice
%   column model par (see ptt.defaultParams; par.H must cover the observed
%   depth range). Interval bottom edges are spaced equally in twtt over the
%   contiguous coherent span below the surface reference depth.
%
%   opts fields (optional): num_intervals (10), half_offset (0 m),
%   min_coverage (0.3), ref_twtt_offset (50e-9 s), inversion ('stripping'
%   for the exact per-interval layer stripping, or 'joint' for the
%   smoothness-regularized joint solve of ptt.invertHorizontalFabricJoint,
%   recommended for noisy data), reg (0.05; joint mode only).
%
%   inv fields (num_intervals x Nblk): dlam, top_depth, bot_depth [m below
%   surface], dtau_obs, dtau_fit [ns], quality (mean coherence at nodes),
%   clipped (joint mode: 1 where the interval dlam is pegged at the
%   eigenvalue bound, 0 where it fitted interior; NaN where the block was
%   skipped or the stripping path runs). Per-block fields (1 x Nblk):
%   rms [ns] (coherence-weighted misfit rms) and alpha (regularization
%   weight), both NaN where the block was skipped or the stripping path
%   runs.

if ~isfield(opts,'num_intervals') || isempty(opts.num_intervals)
  opts.num_intervals = 10;
end
if ~isfield(opts,'half_offset') || isempty(opts.half_offset)
  opts.half_offset = 0;
end
if ~isfield(opts,'min_coverage') || isempty(opts.min_coverage)
  opts.min_coverage = 0.3;
end
if ~isfield(opts,'ref_twtt_offset') || isempty(opts.ref_twtt_offset)
  opts.ref_twtt_offset = 50e-9;
end
if ~isfield(opts,'inversion') || isempty(opts.inversion)
  opts.inversion = 'stripping';
end
if ~ischar(opts.inversion) || ~any(strcmp(opts.inversion, {'stripping','joint'}))
  % Octave's mat2str rejects char arrays, so format the two cases separately
  if ischar(opts.inversion)
    got = ['''' opts.inversion ''''];
  else
    got = sprintf('a %s', class(opts.inversion));
  end
  error('ptt:invertBlocks:inversion', ...
    'opts.inversion must be ''stripping'' or ''joint'' (got %s).', got);
end

Nt = numel(map.Time);
Nblk = numel(blk.starts);
Nint = opts.num_intervals;

[depth_fine, twtt_fine] = ptt.twttDepthMap(par);

inv.dlam = nan(Nint,Nblk);
inv.top_depth = nan(Nint,Nblk);
inv.bot_depth = nan(Nint,Nblk);
inv.dtau_obs = nan(Nint,Nblk);
inv.dtau_fit = nan(Nint,Nblk);
inv.quality = nan(Nint,Nblk);
inv.clipped = nan(Nint,Nblk);
inv.rms = nan(1,Nblk);
inv.alpha = nan(1,Nblk);

for b = 1:Nblk
  % Coherent twtt span below the surface reference depth. Real coherence
  % is patchy, so smooth the good-coverage indicator and keep the range
  % where it stays mostly good, tolerating local gaps (dtau is
  % interpolated across them below).
  t_rel = map.Time - blk.Surface(b);
  good = blk.coverage(:,b) >= opts.min_coverage & isfinite(blk.dtau(:,b));
  smooth_n = max(3, round(Nt/100));
  good_frac = conv(double(good), ones(smooth_n,1)/smooth_n, 'same');
  usable = good_frac >= 0.5 & t_rel > opts.ref_twtt_offset;
  if nnz(usable) < 4*Nint
    continue;
  end
  first_bin = find(usable,1);
  last_bin = find(usable,1,'last');
  t_max = t_rel(last_bin);

  % Interval bottom edges, equally spaced in twtt
  node_twtt = blk.Surface(b) + linspace(t_rel(first_bin) + ...
    (t_max - t_rel(first_bin))/Nint, t_max, Nint).';
  node_depth = interp1(twtt_fine, depth_fine, node_twtt - blk.Surface(b));
  if any(~isfinite(node_depth)) || node_depth(end) >= par.H
    warning('ptt:invertBlocks:depthRange', ...
      'Block %d: coherent depth range exceeds the model column (H = %g m); skipping. Increase par.H.', b, par.H);
    continue;
  end

  obs = [];
  obs.L = opts.half_offset;
  obs.z = par.H - node_depth;
  % Interpolate over masked/incoherent bins using the finite samples only
  fin = isfinite(blk.dtau(:,b));
  obs.dtau = 1e9*interp1(map.Time(fin), blk.dtau(fin,b), node_twtt); % [ns]
  if any(~isfinite(obs.dtau))
    continue;
  end

  node_coh = interp1(map.Time, blk.coh(:,b), node_twtt);
  try
    if strcmp(opts.inversion, 'joint')
      obs.w = node_coh.^2;
      [dlam_prof, inv_out] = ptt.invertHorizontalFabricJoint(obs, par, opts);
    else
      [dlam_prof, inv_out] = ptt.invertHorizontalFabric(obs, par);
    end
  catch ME
    warning('ptt:invertBlocks:failed', ...
      'Block %d: inversion failed (%s); skipping.', b, ME.message);
    continue;
  end

  inv.dlam(:,b) = dlam_prof;
  inv.top_depth(:,b) = par.H - inv_out.ztop;
  inv.bot_depth(:,b) = par.H - inv_out.zbot;
  inv.dtau_obs(:,b) = obs.dtau;
  inv.dtau_fit(:,b) = inv_out.dtau_fit;
  inv.quality(:,b) = node_coh;
  if strcmp(opts.inversion, 'joint')
    inv.clipped(:,b) = double(inv_out.clipped);
    inv.rms(b) = inv_out.rms;
    inv.alpha(b) = inv_out.alpha;
  end
end

end

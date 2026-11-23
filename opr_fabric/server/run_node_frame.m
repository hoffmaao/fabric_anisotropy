%RUN_NODE_FRAME The node method on one frame: a slow axis and dlam(z) for every ~2 km block.
%
%   bash <code>/opr_fabric/server/coreg_one.sh <site_root> 20250108_02 3
%   matlab -batch "addpath('<code>/opr_fabric/server'); day_seg='20250108_02'; frm=3; run_node_frame"
%
% The first line coregisters the frame's channels once (run_coreg_frame,
% 1-3 min); node_survey.sh runs the second over a whole survey.
%
% Reads the frame's channels from cache_dir (default
% stages/quadpol_cf/coreg_cache, where run_coreg_frame writes them: one
% constant delay per channel, which is what the node method needs - a
% depth-varying coregistration removes the delay it measures), cuts the
% frame into blocks of blk_m
% (2000 m) with quadpol_frame_blocks, and runs nodeFabric on each block:
% channel gains from the block itself (reciprocity + firn), the slow axis
% from the azimuth scan of the strength read, dlam(z) from the node
% ladder and the eigen-pair carrier phase by one linear solve under
% nodeFabric's depth model (a smooth dlam(z) plus a bounded window phase;
% the product carries its resolution and the prior-free interval means
% beside it). The synthesis
% assumes one antenna frame per block, so traces heading more than
% spread_max_deg (10) from the block's mean heading are left out (a
% raster corner inside the block), and a block keeping fewer than
% min_traces (500 for 2 km blocks, in proportion for other lengths)
% traces or half its own is skipped.
%
% Writes stages/node/node_<tag>_B<blk_m><tagx>.mat with
%   R      struct array, one per block: the nodeFabric output less the
%          correlation image, plus theta_geo_deg (slow axis, geographic,
%          clockwise from north) and the block's heading
%   meta   quadpol_frame_blocks' block table (x, lat, lon, az, n, used)
%   tag, blk_m, opts, sec
%
% Optional variables: blk_m (2000), z_max (1500), tagx (''), site_root,
% cache_dir, threads (4), spread_max_deg (10), min_traces (500 blk_m /
% 2000), psi_step_deg (5), fc
% (750e6), depth_model (nodeFabric's default; a struct of its fields, or
% 'poster'), node_fc (measured from the data), copol_db ([]: the gains as
% nodeFabric estimates them; a value sets the co-pol magnitude |gr gt|^2 to
% it in dB, keeping the firn's co-pol phase and the reciprocity ratio
% gr/gt - the instrument value from the heading flip of the deep co-pol
% ratio is -3.5 dB at Ridge A, the firn's -4.8; give it a tagx). To remake
% a product as it was before the depth model existed: depth_model =
% 'poster'; node_fc = 750e6; and give it a tagx.
[fabric_code, fabric_work] = fabric_paths();
addpath(fabric_code);
if ~exist('threads', 'var'), threads = 4; end
maxNumCompThreads(threads);
if ~exist('day_seg', 'var'), error('set day_seg'); end
if ~exist('frm', 'var'), error('set frm'); end
day_seg = char(day_seg);
if ~exist('site_root', 'var'), site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2'; end
if ~exist('cache_dir', 'var'), cache_dir = fullfile(fabric_work, 'stages', 'quadpol_cf', 'coreg_cache'); end
if ~exist('blk_m', 'var'), blk_m = 2000; end
if ~exist('z_max', 'var'), z_max = 1500; end
if ~exist('tagx', 'var'), tagx = ''; end
if ~exist('spread_max_deg', 'var'), spread_max_deg = 10; end
if ~exist('min_traces', 'var'), min_traces = 500 * blk_m / 2000; end
if ~exist('psi_step_deg', 'var'), psi_step_deg = 5; end
if ~exist('fc', 'var'), fc = 750e6; end
if ~exist('copol_db', 'var'), copol_db = []; end
dt = 1e9 / 300e6;                                   % accumulation radar, 300 MHz sampling
tag = sprintf('%s_%03d', day_seg, frm);
out_dir = fullfile(fabric_work, 'stages', 'node');
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
out_fn = fullfile(out_dir, sprintf('node_%s_B%d%s.mat', tag, round(blk_m), tagx));

cache_fn = fullfile(cache_dir, ['creg_' tag '.mat']);
if exist(cache_fn, 'file') ~= 2
  error('run_node_frame:noCache', ['no coregistered channels for %s (%s). Build them first:\n' ...
    '  bash %s %s %s %d'], tag, cache_fn, fullfile(fabric_code, 'opr_fabric', 'server', 'coreg_one.sh'), ...
    site_root, day_seg, frm);
end
[Dc, meta, z] = quadpol_frame_blocks(day_seg, frm, cache_dir, site_root, blk_m, z_max, 100);
opts = struct('fc', fc, 'calibrate', true, 'psi_step_deg', psi_step_deg);
if exist('depth_model', 'var'), opts.depth_model = depth_model; end
if exist('node_fc', 'var'), opts.node_fc = node_fc; end
fprintf('%s: %d blocks of ~%d m from %s\n', tag, numel(Dc), round(blk_m), cache_dir);
t0 = tic;
R = struct([]);
for b = 1:numel(Dc)
  S = Dc{b};
  hd = S.heading(:);
  % one antenna frame per block: keep the traces within spread_max_deg of
  % the block's circular-mean heading (a raster corner inside the block is
  % cut away, the straight part kept), skip only if too little is left
  dev = abs(angle(exp(2i * (hd - deg2rad(meta.az(b)))))) / 2;
  keep = dev <= deg2rad(spread_max_deg);
  spread = rad2deg(prctile(dev, 95));
  rb = struct('block', b, 'heading_deg', meta.az(b), 'heading_spread_deg', spread, ...
    'n_traces', numel(hd), 'n_kept', nnz(keep), 'skipped', false, 'error', '');
  if nnz(keep) < max(min_traces, 0.5 * numel(hd))
    fprintf('  block %d: heading spread %.1f deg, %d of %d traces within %g deg, skipped\n', ...
      b, spread, nnz(keep), numel(hd), spread_max_deg);
    rb.skipped = true;
  else
    try
      for ch = {'hh', 'vv', 'hv', 'vh'}, S.(ch{1}) = S.(ch{1})(:, keep); end
      hk = mod(rad2deg(angle(mean(exp(2i * hd(keep))))) / 2, 180);
      rb.heading_deg = hk;
      meta.az(b) = hk;
      ob = opts;
      if ~isempty(copol_db)
        ob.gains = H_copol_gains(rmfield(S, 'heading'), z, copol_db);
      end
      o = ptt.nodeFabric(rmfield(S, 'heading'), z, dt, ob);
      o.windows = rmfield(o.windows, 'corr');
      o.theta_geo_deg = mod(meta.az(b) + rad2deg(o.psi0), 180);
      N = o.nodes;
      fprintf(['  block %d: heading %.1f, gains |gr| %.2f %+.0f deg |gt| %.2f %+.0f deg, slow axis %.1f geo ' ...
        '(plateau ratio %.2f, contrast 45/0 %.2f/%.2f dB), %d nodes (%s) at %.1f MHz, chi2/dof %.2f, ' ...
        'abstain below %g m, dlam 900-1400 %.3f\n'], b, meta.az(b), abs(o.gains.gr), rad2deg(angle(o.gains.gr)), ...
        abs(o.gains.gt), rad2deg(angle(o.gains.gt)), o.theta_geo_deg, o.theta_scan.plateau_ratio, ...
        o.theta_scan.contrast_45_db, o.theta_scan.contrast_0_db, numel(N.z), N.n_source, median(N.fc_hz) / 1e6, ...
        o.chi2_dof, o.abstain_below_m, mean(o.dlam(o.zm >= 900 & o.zm <= 1400), 'omitnan'));
      fprintf('           window phase %.4f ns, %.1f numbers resolved, one value a mean over %.0f m, residual lag-one %+.2f\n', ...
        o.depth_model.psi_sd_ns, o.n_eff, median(o.resolution_m(o.zm >= 200), 'omitnan'), o.fit.resid_acf1);
      fn = fieldnames(o);
      for k = 1:numel(fn), rb.(fn{k}) = o.(fn{k}); end
    catch err
      fprintf('  block %d: FAILED: %s\n', b, err.message);
      rb.error = err.message;
    end
  end
  if isempty(R), R = rb; else, R(end+1) = H_align(rb, R(1)); end %#ok<SAGROW>
end
sec = toc(t0);
fprintf('%s: %d blocks in %.1f min\n', tag, numel(Dc), sec / 60);
save(out_fn, '-v7.3', 'R', 'meta', 'tag', 'blk_m', 'opts', 'sec', 'cache_dir', 'z', 'copol_db');
fprintf('wrote %s\n', out_fn);

function g = H_copol_gains(S, z, copol_db)
% ptt.channelGains on the block (reciprocity + firn, the window nodeFabric
% uses), then the co-pol magnitude |gr gt| replaced by the instrument value
% copol_db (|VV/HH|^2, dB); the co-pol phase and gr/gt are kept, and the
% root taken with Re gr > 0 as channelGains does
g0 = ptt.channelGains(S, z, struct('firn_m', [30 90]));
G = g0.gr * g0.gt;
G = G / abs(G) * 10 ^ (copol_db / 20);
gr = sqrt(G * g0.gr / g0.gt);
if real(gr) < 0, gr = -gr; end
g = struct('gr', gr, 'gt', G / gr);
end

function rb = H_align(rb, ref)
% struct arrays need one field set: fill what this block lacks with []
fn = fieldnames(ref);
for k = 1:numel(fn)
  if ~isfield(rb, fn{k}), rb.(fn{k}) = []; end
end
fr = fieldnames(rb);
extra = setdiff(fr, fn);
rb = rmfield(rb, extra);
rb = orderfields(rb, ref);
end

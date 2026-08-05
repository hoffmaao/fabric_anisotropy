%RUN_SECTIONS High-resolution 2D fabric sections for all three SCAR profiles.
%
% Re-runs each profile at the settings the Ridge A sweep established
% (ridge_a_section.m): small along-track blocks, more depth intervals, and
% a coherence gate chosen per frame rather than left at 0.5. The gate is
% swept and picked on adjacent-block agreement, not misfit rms - rms is a
% post-regularization residual and always improves with more free
% parameters, so it would endorse an overfit.
%
% dtau source per site is NOT uniform, deliberately:
%   Ridge A   phase   SNAPHU there is exactly consistent with the wrapped
%                     phase, and delta-k is suppressed at depth
%   Thwaites  delta-k the joint chain does not fit its own observations
%                     there (5.2 ns rms) - blend corruption
%   NEGIS     delta-k the only option; no SNAPHU, no coregistration
%
% Saves dlam, the node coherence (dlam_quality) used for shading, the
% interpolated-node flag, and block positions for each site.

scratch = '/kucresis/scratch/hoffmana_sta/fabric';
out_fn = fullfile(scratch, 'stages', 'fabric_sections.mat');

code = fullfile(scratch, 'code');
addpath(code); addpath(fullfile(code,'opr_fabric'));
addpath(fullfile(code,'opr_fabric','test','stubs'));

% name, file, dtau source
sites = { ...
  'ridge_a',  ['/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2/' ...
               'CSARP_polarimetric/20250108_02/Data_20250108_02_009.mat'], 'phase'; ...
  'thwaites', ['/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground/' ...
               'CSARP_polarimetric_unwrap/20240108_01/Data_20240108_01_001.mat'], 'deltak'; ...
  'negis',    [scratch '/2024_Greenland_Ground2/CSARP_polarimetric_negis/' ...
               '20240626_03/Data_20240626_03_001.mat'], 'deltak'};

BLOCK = 125;
NINT = 25;
REG = 0.05;
THR = [0.5 0.35 0.2 0.1];

par = ptt.defaultParams();
par.H = 2000; par.lam_z_sfc = 1/3; par.lam_z_bed = 1/3;
par.zhat_bco = 1 - 60/par.H;

S = struct();
for si = 1:size(sites,1)
  [name, in_fn, src] = sites{si,:};
  fprintf('\n=== %s (%s) ===\n%s\n', name, src, in_fn);
  if exist(in_fn,'file') ~= 2
    warning('missing %s', in_fn); continue;
  end

  deltak_en = strcmp(src, 'deltak');
  if deltak_en
    want = {'interferogram_coherence','row_offset','Time','GPS_time', ...
      'Latitude','Longitude','Elevation','Surface','param_records', ...
      'param_polarimetric','ref','sec'};
  else
    want = {'interferogram_mlook','interferogram_coherence', ...
      'snaphu_out_phase','row_offset','Time','GPS_time','Latitude', ...
      'Longitude','Elevation','Surface','param_records','param_polarimetric'};
  end
  have = whos('-file', in_fn);
  sel = intersect(want, {have.name});
  pol = load(in_fn, sel{:});

  map = [];
  map.Time = pol.Time;
  map.Surface = pol.Surface;
  map.fc = 750e6;
  map.coherence = abs(pol.interferogram_coherence);
  map.img_comb = [];
  for pname = {'param_polarimetric','param_records'}
    p = pname{1};
    if isfield(pol,p) && isstruct(pol.(p)) && isfield(pol.(p),'array') ...
        && isfield(pol.(p).array,'img_comb') && ~isempty(pol.(p).array.img_comb)
      map.img_comb = pol.(p).array.img_comb;
      if isfield(pol.(p).array,'img_comb_mult')
        map.img_comb_mult = pol.(p).array.img_comb_mult;
      end
      break;
    end
  end
  if isfield(pol,'row_offset'), map.row_offset = pol.row_offset;
  else, map.row_offset = []; end

  cm = map.coherence(isfinite(map.coherence));
  fprintf('  coherence percentiles 10/50/90: %.2f %.2f %.2f\n', ...
    prctile_local(cm,10), prctile_local(cm,50), prctile_local(cm,90));

  o = [];
  o.fc = 750e6; o.phase_sign = -1; o.min_coverage = 0.3;
  o.ref_twtt_offset = 50e-9; o.half_offset = 0; o.blend_coreg_en = true;
  o.block_size = BLOCK; o.inversion = 'joint'; o.reg = REG;
  o.num_intervals = NINT;
  o.coherence_threshold = min(THR);   % loosest for the one dtau pass

  if deltak_en
    map.phase = []; map.phase_is_unwrapped = true;
    slc = struct('ref', pol.ref, 'sec', pol.sec);
    pol = rmfield(pol, {'ref','sec'});
    [dtau, info] = ptt.deltakTraveltime(slc, map, o);
    clear slc;
  else
    map.phase = pol.snaphu_out_phase;
    map.phase_is_unwrapped = true;
    [dtau, info] = ptt.blendTraveltime(map, o);
  end
  fprintf('  phase_sign %+d\n', info.phase_sign);

  % Every gate is scored first and the winner picked afterwards, because
  % the rule below is relative to the BEST agreement, not to whichever
  % gate happened to be scored before it - comparing against a running
  % best lets each successive gate give up another 0.05 and the losses
  % compound into a gate that agrees far worse than the best one.
  cand = {};
  fprintf('  %-6s %8s %9s %9s %6s\n','gate','rms_ns','adj_corr','deep_m','nblk');
  for ti = 1:numel(THR)
    o2 = o; o2.coherence_threshold = THR(ti);
    % dtau does not depend on the gate; only the mask it is averaged under
    info2 = info;
    [info2.coh_mask, info2.ref_bin] = ptt.surfaceReference(map, o2);
    bad = ~isfinite(dtau);
    info2.coh_mask(bad) = false;
    blk = ptt.blockAverage(dtau, map, info2, o2);
    iv = ptt.invertBlocks(blk, map, par, o2);
    dl = iv.dlam;
    cs = [];
    for b = 1:size(dl,2)-1
      a = dl(:,b); c = dl(:,b+1);
      ok = isfinite(a) & isfinite(c);
      if nnz(ok) > 3 && std(a(ok)) > 0 && std(c(ok)) > 0
        cc = corrcoef(a(ok), c(ok)); cs(end+1) = cc(1,2); %#ok<SAGROW>
      end
    end
    adj = median(cs);
    rms = median(iv.rms(isfinite(iv.rms)));
    deep = max(iv.bot_depth(:));
    nblk = numel(blk.starts);
    fprintf('  %-6.2f %8.3f %9.2f %9.0f %6d\n', THR(ti), rms, adj, deep, nblk);
    % adj is NaN when no adjacent block pair had enough jointly finite
    % intervals to correlate (a single-block frame, or a gate that leaves
    % the record too sparse). Such a gate cannot be scored, so it cannot
    % be chosen.
    if ~isfinite(adj), continue; end
    cand{end+1} = struct('adj', adj, 'gate', THR(ti), 'rms', rms, ...
      'deep', deep, ...
      'nblk', nblk, 'dlam', single(iv.dlam), 'top', single(iv.top_depth), ...
      'bot', single(iv.bot_depth), 'quality', single(iv.quality), ...
      'interpolated', single(iv.interpolated), ...
      'clipped', single(iv.clipped), ...
      'lat', cellfun(@(c) mean(pol.Latitude(c)), blk.cols), ...
      'lon', cellfun(@(c) mean(pol.Longitude(c)), blk.cols), ...
      'src', src, 'nint', NINT, 'block', BLOCK, ...
      'phase_sign', info.phase_sign); %#ok<SAGROW>
  end
  if isempty(cand)
    warning('%s: no coherence gate produced a scorable inversion; skipping', ...
      name);
    clear dtau pol map;
    continue;
  end

  % Among gates that agree about as well as the BEST one does, take the
  % one that reaches DEEPEST. A looser gate admits weakly-constrained
  % cells rather than dropping them, and the figure shades by node
  % coherence, so those cells arrive visibly faded instead of silently
  % absent - which is strictly more informative than a blank lower half.
  % Thwaites is the case in point: gate 0.10 reaches 1977 m against
  % 1533 m at 0.35, for an agreement difference of 0.02.
  adjs = cellfun(@(c) c.adj, cand);
  deeps = cellfun(@(c) c.deep, cand);
  deeps(adjs < max(adjs) - 0.05) = -Inf;
  [~, bi] = max(deeps);
  best = cand{bi};
  fprintf('  chosen gate %.2f: rms %.3f ns, adj %.2f, %d blocks, %.0f m\n', ...
    best.gate, best.rms, best.adj, best.nblk, best.deep);
  S.(name) = best;
  clear cand dtau pol map;
end

save(out_fn, '-v7', 'S');
fprintf('\nwrote %s\n', out_fn);

function p = prctile_local(v, q)
v = sort(v(:));
p = v(max(1, min(numel(v), round(q/100*numel(v)))));
end

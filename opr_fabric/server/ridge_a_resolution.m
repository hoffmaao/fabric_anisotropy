%RIDGE_A_RESOLUTION How finely can the Ridge A fabric actually be resolved?
%
% The shipped Ridge A inversion runs 10 depth intervals over 5 blocks of a
% 4.3 km frame, which is far coarser than the data deserve: SNAPHU tracks
% the wrapped phase exactly there (4.43 fringes either way you count them)
% and the joint chain fits to 0.16 ns against a 5.9 ns signal.
%
% This sweeps depth resolution, along-track block size and regularization
% strength, and scores each setting three ways:
%   rms       coherence-weighted misfit (ns) - rises when the model can no
%             longer fit, i.e. too FEW intervals
%   clipped   fraction of intervals pegged at the eigenvalue bound - the
%             signature of exact-solve noise amplification, i.e. too MANY
%   adj_corr  median correlation between the dlam profiles of ADJACENT
%             blocks. Ridge A is a divide, so the fabric varies slowly
%             along track; neighbouring blocks that disagree are reporting
%             noise, not structure. This is the resolution test that
%             matters - rms alone always improves with more free
%             parameters and would happily pick an overfit.
%
% Prints a table and saves every sweep result for the 2D section figure.

% <work> is derived from this script's own location - see fabric_paths.
[~, scratch] = fabric_paths();
in_fn = ['/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2/' ...
  'CSARP_polarimetric/20250108_02/Data_20250108_02_009.mat'];
out_fn = fullfile(scratch, 'stages', 'ridge_a_resolution.mat');

code = fullfile(scratch, 'code');
addpath(code); addpath(fullfile(code,'opr_fabric'));
addpath(fullfile(code,'opr_fabric','test','stubs'));

fprintf('Loading %s\n', in_fn);
want = {'interferogram_mlook','interferogram_coherence','snaphu_out_phase', ...
  'row_offset','Time','GPS_time','Latitude','Longitude','Elevation', ...
  'Surface','param_records','param_polarimetric'};
have = whos('-file', in_fn);
sel = intersect(want, {have.name});
pol = load(in_fn, sel{:});

map = [];
map.Time = pol.Time;
map.Surface = pol.Surface;
map.fc = 750e6;
map.coherence = abs(pol.interferogram_coherence);
map.phase = pol.snaphu_out_phase;
map.phase_is_unwrapped = true;
map.row_offset = pol.row_offset;
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
Nx = size(map.coherence, 2);
fprintf('  %d samples x %d traces\n', numel(map.Time), Nx);

opts0 = [];
opts0.fc = 750e6;
opts0.phase_sign = -1;
opts0.coherence_threshold = 0.5;
opts0.min_coverage = 0.3;
opts0.ref_twtt_offset = 50e-9;
opts0.half_offset = 0;
opts0.blend_coreg_en = true;

% dtau is per pixel: form it once, it does not depend on the sweep
[dtau, info] = ptt.blendTraveltime(map, opts0);
fprintf('  phase_sign %+d\n', info.phase_sign);

par = ptt.defaultParams();
pp = struct('H', 2000, 'bco_depth', 60, 'lam_z_sfc', 1/3, 'lam_z_bed', 1/3);
for fld = fieldnames(pp).'
  par.(fld{1}) = pp.(fld{1});
end
par.zhat_bco = 1 - par.bco_depth/par.H;
par = rmfield(par, 'bco_depth');

BLOCKS = [1000 500 250];
NINT   = [10 15 20 30];
REG    = [0.05 0.01 0.002];

res = struct('block_size',{},'num_intervals',{},'reg',{},'rms',{}, ...
  'clipped',{},'adj_corr',{},'dlam',{},'bot',{},'top',{},'lat',{}, ...
  'lon',{},'nblk',{});

fprintf('\n%-7s %-5s %-7s %8s %9s %9s %6s\n', 'block', 'nint', 'reg', ...
  'rms_ns', 'clipped%', 'adj_corr', 'nblk');
for bi = 1:numel(BLOCKS)
  o = opts0; o.block_size = BLOCKS(bi);
  blk = ptt.blockAverage(dtau, map, info, o);
  nblk = numel(blk.starts);
  lat = cellfun(@(c) mean(pol.Latitude(c)), blk.cols);
  lon = cellfun(@(c) mean(pol.Longitude(c)), blk.cols);
  for ni = 1:numel(NINT)
    for ri = 1:numel(REG)
      o2 = o;
      o2.num_intervals = NINT(ni);
      o2.inversion = 'joint';
      o2.reg = REG(ri);
      try
        iv = ptt.invertBlocks(blk, map, par, o2);
      catch ME
        fprintf('%-7d %-5d %-7g   FAILED: %s\n', BLOCKS(bi), NINT(ni), ...
          REG(ri), ME.message);
        continue;
      end
      dl = iv.dlam;
      % adjacent-block agreement, over intervals both blocks resolved
      cs = [];
      for b = 1:size(dl,2)-1
        a = dl(:,b); c = dl(:,b+1);
        ok = isfinite(a) & isfinite(c);
        if nnz(ok) > 3 && std(a(ok)) > 0 && std(c(ok)) > 0
          cc = corrcoef(a(ok), c(ok));
          cs(end+1) = cc(1,2); %#ok<SAGROW>
        end
      end
      k = numel(res) + 1;
      res(k).block_size = BLOCKS(bi);
      res(k).num_intervals = NINT(ni);
      res(k).reg = REG(ri);
      res(k).rms = median(iv.rms(isfinite(iv.rms)));
      res(k).clipped = mean(iv.clipped(isfinite(iv.clipped)));
      res(k).adj_corr = median(cs);
      res(k).dlam = single(dl);
      res(k).top = single(iv.top_depth);
      res(k).bot = single(iv.bot_depth);
      res(k).lat = lat; res(k).lon = lon; res(k).nblk = nblk;
      fprintf('%-7d %-5d %-7g %8.3f %9.1f %9.2f %6d\n', BLOCKS(bi), ...
        NINT(ni), REG(ri), res(k).rms, 100*res(k).clipped, ...
        res(k).adj_corr, nblk);
    end
  end
end

save(out_fn, '-v7', 'res');
fprintf('\nwrote %s\n', out_fn);

%RIDGE_A_EDGE_TEST Is the deep Ridge A downturn physical or a domain edge?
%
% The inverted profile rises to dlam ~0.073 near 1500 m and falls to ~0.054
% by 1850 m. A model-free gradient of the block-averaged dtau shows the
% same fall, so it is not manufactured by the inversion - but the deepest
% gradient estimate swings with the smoothing window (0.043 at a 100 m
% half-window against 0.060 at 200 m), which is what an edge does. Three
% tests separate the possibilities:
%
%   nint    if the last interval turns down wherever it happens to land,
%           the downturn tracks the DISCRETIZATION, not the ice
%   H       the column model's bed sits at par.H = 2000 m while the ice at
%           Ridge A is far thicker; if moving the bed to 2800 m changes
%           the deep profile, the domain bottom is leaking in
%   trunc   invert with the record artificially cut at 1900/1700/1500 m
%           and compare over the depths all three share. An edge effect
%           propagates inward: values well above the cut should NOT move
%           when the cut moves. This is the decisive one.
%
% Also reports node coherence against depth, since a coherence-driven bias
% toward zero would mimic a real weakening.

% Repo root and <work> are derived from this script's own location - see
% fabric_paths.
[code, scratch] = fabric_paths();
in_fn = ['/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2/' ...
  'CSARP_polarimetric/20250108_02/Data_20250108_02_009.mat'];
out_fn = fullfile(scratch, 'stages', 'ridge_a_edge_test.mat');

addpath(code); addpath(fullfile(code,'opr_fabric'));
addpath(fullfile(code,'opr_fabric','test','stubs'));

want = {'interferogram_mlook','interferogram_coherence','snaphu_out_phase', ...
  'row_offset','Time','Surface','Latitude','Longitude','param_records', ...
  'param_polarimetric'};
have = whos('-file', in_fn);
sel = intersect(want, {have.name});
pol = load(in_fn, sel{:});

map0 = [];
map0.Time = pol.Time;
map0.Surface = pol.Surface;
map0.fc = 750e6;
map0.coherence = abs(pol.interferogram_coherence);
map0.phase = pol.snaphu_out_phase;
map0.phase_is_unwrapped = true;
map0.row_offset = pol.row_offset;
map0.img_comb = [];
for pname = {'param_polarimetric','param_records'}
  p = pname{1};
  if isfield(pol,p) && isstruct(pol.(p)) && isfield(pol.(p),'array') ...
      && isfield(pol.(p).array,'img_comb') && ~isempty(pol.(p).array.img_comb)
    map0.img_comb = pol.(p).array.img_comb;
    if isfield(pol.(p).array,'img_comb_mult')
      map0.img_comb_mult = pol.(p).array.img_comb_mult;
    end
    break;
  end
end

o0 = [];
o0.fc = 750e6; o0.phase_sign = -1; o0.min_coverage = 0.3;
o0.ref_twtt_offset = 50e-9; o0.half_offset = 0; o0.blend_coreg_en = true;
o0.block_size = 125; o0.inversion = 'joint'; o0.reg = 0.05;
o0.coherence_threshold = 0.10;

[dtau, info] = ptt.blendTraveltime(map0, o0);

C_ICE = 1.68e8;   % for the truncation cutoffs only
surf_med = median(map0.Surface(isfinite(map0.Surface)));

R = struct('tag',{},'nint',{},'H',{},'cut',{},'zc',{},'dlam',{},'qual',{});

fprintf('\n--- depth discretization ---\n');
for nint = [15 25 40]
  R = addcase(R, sprintf('nint%d', nint), nint, 2000, Inf, map0, dtau, ...
    info, o0, surf_med, C_ICE);
end

fprintf('\n--- column-model bed depth ---\n');
for H = [2000 2400 2800]
  R = addcase(R, sprintf('H%d', H), 25, H, Inf, map0, dtau, info, o0, ...
    surf_med, C_ICE);
end

fprintf('\n--- truncated record (the decisive test) ---\n');
for cut = [1900 1700 1500]
  R = addcase(R, sprintf('cut%d', cut), 25, 2000, cut, map0, dtau, info, ...
    o0, surf_med, C_ICE);
end

save(out_fn, '-v7', 'R');
fprintf('\nwrote %s\n', out_fn);

%% ---- local helpers (MATLAB requires them after all script code) ----
function par = mkpar(H)
  par = ptt.defaultParams();
  par.H = H; par.lam_z_sfc = 1/3; par.lam_z_bed = 1/3;
  par.zhat_bco = 1 - 60/H;
end

function R = addcase(R, tag, nint, H, cut, map, dtau, info, o0, surf_med, C_ICE)
  o = o0; o.num_intervals = nint;
  m = map;
  if isfinite(cut)
    % force the coherent span to end at `cut` metres by blanking the
    % coherence below it; blockAverage and invertBlocks then see a record
    % that simply stops there
    tcut = surf_med + 2*cut/C_ICE;
    m.coherence(m.Time > tcut, :) = 0;
  end
  [ci, ri] = ptt.surfaceReference(m, o);
  inf2 = info; inf2.coh_mask = ci & isfinite(dtau); inf2.ref_bin = ri;
  blk = ptt.blockAverage(dtau, m, inf2, o);
  iv = ptt.invertBlocks(blk, m, mkpar(H), o);
  k = numel(R) + 1;
  R(k).tag = tag; R(k).nint = nint; R(k).H = H; R(k).cut = cut;
  R(k).zc = median((iv.top_depth + iv.bot_depth)/2, 2, 'omitnan').';
  R(k).dlam = median(iv.dlam, 2, 'omitnan').';
  R(k).qual = median(iv.quality, 2, 'omitnan').';
  fprintf('  %-18s nint %2d  H %4d  cut %6s  deepest %.0f m, dlam %.3f\n', ...
    tag, nint, H, num2str(cut), max(R(k).zc), R(k).dlam(end));
end

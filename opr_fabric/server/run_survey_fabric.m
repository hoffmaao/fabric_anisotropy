%RUN_SURVEY_FABRIC Fabric inversion over every frame of one survey.
%
% run_sections.m inverts ONE leg per site at high vertical resolution to
% draw a section. This runs the same chain over every frame of a survey so
% the result can be MAPPED rather than sectioned - depth-averaged
% horizontal asymmetry along every track, plus, where legs cross, the
% fabric orientation.
%
% Pick the survey with `site`, default ridge_a:
%   matlab -batch "site='taylor_dome'; run_survey_fabric"
%
% AZIMUTH DIVERSITY IS A PROPERTY OF THE SURVEY, NOT OF THE METHOD. The
% orientation solve needs two legs at least ~20 deg apart over the same
% ground. Measured over each product set:
%   ridge_a      raster, both azimuth families        orientation solvable
%   taylor_dome  raster, 60-80 and 120-140 deg        orientation solvable
%   negis        radiating traverse, 33 and 128-145   orientation solvable
%   thwaites     single W-E transect, ALL 66 frames
%                between 60 and 100 deg               ribbons ONLY
% Thwaites is listed so it can be inverted and mapped as ribbons, not
% because a two-azimuth solve will produce anything there; the figure
% drops the strength panel for it rather than showing an empty one.
%
% Settings are run_sections.m's verbatim (BLOCK/NINT/REG/THR and the
% adjacent-block gate rule), because the map and the section have to be the
% same measurement: a map drawn at different blocks or a different gate
% would disagree with the published section along the very leg they share.
%
% dtau comes from the SNAPHU phase at every site, and that is not a
% preference. delta-k is suppressed at depth (the stage-A diagnosis on
% Ridge A), which flattens the very depth range the map averages over; the
% EastGRIP core comparison then measured the cost, delta-k scoring RMS
% 0.411 against SNAPHU's 0.282 through the same inversion.
%
% WHY ORIENTATION NEEDS CROSSING LEGS. One line measures only the
% projection of the horizontal ellipse onto its own axes, dlam_obs =
% -P cos 2(alpha - theta). Two lines at different azimuths over the same
% ground determine both P and theta; along a single leg only the projection
% is recoverable. That distinction is carried into the figure - ribbons are
% per-leg projections, orientation crosses are drawn only where two
% azimuths meet - which is why the table above matters before running a
% site rather than after.
%
% Saves per-frame dlam(interval, block) with block positions, so the depth
% average and the crossing search both happen downstream and can be retuned
% without another server pass.
% <work> is derived from this script's own location - see fabric_paths.
[~, scratch] = fabric_paths();
if ~exist('site', 'var') || isempty(site)
  site = 'ridge_a';
end
% dir / frame-glob per survey. The glob matters for Thwaites, whose season
% also holds two unrelated sites (a -79.2 deg pair and a Ross Island group
% at 167-168 deg E) that would otherwise be inverted and mapped as if they
% were part of the transect.
switch site
  case 'ridge_a'
    pol_dir = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2/CSARP_polarimetric';
    glob = 'Data_*.mat';
  case 'taylor_dome'
    pol_dir = '/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2/CSARP_polarimetric';
    glob = 'Data_*.mat';
  case 'thwaites'
    pol_dir = '/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground/CSARP_polarimetric_unwrap';
    glob = 'Data_2024010*.mat';
  case 'eastwind'
    % 17 frames, all carrying snaphu_out_phase, over a ~1 km footprint at
    % 77.67 S 168.14 E with azimuths spread 20-140 deg. The tight footprint
    % and wide azimuth spread is the best geometry of any site here - it is
    % close to a rotation about one point, which is what a quad-pol ApRES
    % does mechanically.
    pol_dir = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground/CSARP_polarimetric_unwrap';
    glob = 'Data_*.mat';
  case 'negis'
    % Written by run_negis_fabric.m, which repackages qlook HH/VV and adds
    % the Goldstein-filtered SNAPHU phase this season has no product for.
    pol_dir = fullfile(scratch, '2024_Greenland_Ground2', ...
      'CSARP_polarimetric_negis');
    glob = 'Data_*.mat';
  otherwise
    error('run_survey_fabric:site', 'unknown site %s', site);
end
out_fn = fullfile(scratch, 'stages', sprintf('%s_survey.mat', site));

code = fullfile(scratch, 'code');
addpath(code); addpath(fullfile(code,'opr_fabric'));
addpath(fullfile(code,'opr_fabric','test','stubs'));

BLOCK = 125;            % traces per along-track block   (run_sections.m)
NINT = 25;              % depth intervals                (run_sections.m)
REG = 0.05;             % regularisation                 (run_sections.m)
THR = [0.5 0.35 0.2 0.1];   % coherence gates swept      (run_sections.m)

% From the defaults, not from an empty struct: the firn model needs
% rho_sfc/rho_bco/e0/p/nz as well, and only the four overridden here are
% site choices. run_sections.m builds par the same way, which is what
% keeps this map and the published section the same measurement.
par = ptt.defaultParams();
par.H = 2000; par.lam_z_sfc = 1/3; par.lam_z_bed = 1/3;
par.zhat_bco = 1 - 60/par.H;

d = dir(fullfile(pol_dir, '*', glob));
if isempty(d)
  error('run_survey_fabric:noFrames', 'no frames matching %s under %s', ...
    glob, pol_dir);
end
fprintf('site %s: found %d frames under %s\n', site, numel(d), pol_dir);

F = struct([]);
n = 0;
for fi = 1:numel(d)
  in_fn = fullfile(d(fi).folder, d(fi).name);
  tag = regexprep(d(fi).name, '^Data_|\.mat$', '');
  fprintf('\n=== %s (%d of %d) ===\n', tag, fi, numel(d));
  try
    want = {'interferogram_mlook','interferogram_coherence', ...
      'snaphu_out_phase','row_offset','Time','GPS_time','Latitude', ...
      'Longitude','Elevation','Surface','param_records','param_polarimetric'};
    have = whos('-file', in_fn);
    sel = intersect(want, {have.name});
    pol = load(in_fn, sel{:});
    if ~isfield(pol, 'snaphu_out_phase') || isempty(pol.snaphu_out_phase)
      error('run_survey_fabric:noSnaphu', ...
        'no snaphu_out_phase in this product');
    end

    map = [];
    map.Time = pol.Time;
    map.Surface = pol.Surface;
    map.fc = 750e6;
    map.coherence = abs(pol.interferogram_coherence);
    map.img_comb = [];
    for pname = {'param_polarimetric','param_records'}
      p = pname{1};
      if isfield(pol,p) && isstruct(pol.(p)) && isfield(pol.(p),'array') ...
          && isfield(pol.(p).array,'img_comb') ...
          && ~isempty(pol.(p).array.img_comb)
        map.img_comb = pol.(p).array.img_comb;
        if isfield(pol.(p).array,'img_comb_mult')
          map.img_comb_mult = pol.(p).array.img_comb_mult;
        end
        break;
      end
    end
    if isfield(pol,'row_offset'), map.row_offset = pol.row_offset;
    else, map.row_offset = []; end

    o = [];
    o.fc = 750e6; o.phase_sign = -1; o.min_coverage = 0.3;
    o.ref_twtt_offset = 50e-9; o.half_offset = 0; o.blend_coreg_en = true;
    o.block_size = BLOCK; o.inversion = 'joint'; o.reg = REG;
    o.num_intervals = NINT;
    o.coherence_threshold = min(THR);

    map.phase = pol.snaphu_out_phase;
    map.phase_is_unwrapped = true;
    [dtau, info] = ptt.blendTraveltime(map, o);

    % Gate sweep, scored on adjacent-block agreement rather than misfit
    % rms - rms is a post-regularisation residual and always improves with
    % more free parameters, so it would endorse an overfit. Every gate is
    % scored before the winner is picked, because the rule is relative to
    % the BEST agreement, not a running best.
    best = []; best_adj = -Inf;
    for ti = 1:numel(THR)
      o2 = o; o2.coherence_threshold = THR(ti);
      info2 = info;
      [info2.coh_mask, info2.ref_bin] = ptt.surfaceReference(map, o2);
      info2.coh_mask(~isfinite(dtau)) = false;
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
      fprintf('  gate %.2f: %d blocks, rms %.3f ns, adj %.2f\n', THR(ti), ...
        numel(blk.starts), median(iv.rms(isfinite(iv.rms))), adj);
      if ~isfinite(adj), continue; end
      if adj > best_adj
        best_adj = adj;
        best = struct('gate', THR(ti), 'adj', adj, ...
          'rms', median(iv.rms(isfinite(iv.rms))), ...
          'dlam', single(iv.dlam), 'top', single(iv.top_depth), ...
          'bot', single(iv.bot_depth), 'quality', single(iv.quality), ...
          'interpolated', single(iv.interpolated), ...
          'clipped', single(iv.clipped), ...
          'lat', cellfun(@(c) mean(pol.Latitude(c)), blk.cols), ...
          'lon', cellfun(@(c) mean(pol.Longitude(c)), blk.cols), ...
          'lat0', cellfun(@(c) pol.Latitude(c(1)), blk.cols), ...
          'lon0', cellfun(@(c) pol.Longitude(c(1)), blk.cols), ...
          'lat1', cellfun(@(c) pol.Latitude(c(end)), blk.cols), ...
          'lon1', cellfun(@(c) pol.Longitude(c(end)), blk.cols));
      end
    end
    if isempty(best)
      error('run_survey_fabric:noGate', ...
        'no coherence gate produced a scorable block set');
    end
    fprintf('  chose gate %.2f (adj %.2f, rms %.3f ns, %d blocks)\n', ...
      best.gate, best.adj, best.rms, numel(best.lat));

    n = n + 1;
    best.tag = tag;
    F = H_append(F, best, n);
    % Save after EVERY frame, not once at the end. At ~3 min a frame this
    % is a two-hour batch, and losing all of it to a kill, a full disk or a
    % single bad product late in the list is a worse outcome than the
    % second it costs to rewrite a file this small.
    save(out_fn, '-v7.3', 'F');
    clear pol map dtau info;
  catch ME
    % Warn and continue: 36 frames is a batch, and one product missing its
    % unwrapped phase must not cost the other 35.
    warning('run_survey_fabric:frameFailed', '%s failed (%s): %s', ...
      tag, ME.identifier, ME.message);
    clear pol map dtau info;
    continue;
  end
end

if n < 2
  error('run_survey_fabric:tooFew', ...
    'only %d of %d frames inverted', n, numel(d));
end
save(out_fn, '-v7.3', 'F');
fprintf('\nwrote %s (%d of %d frames)\n', out_fn, n, numel(d));

function F = H_append(F, s, n)
fn = fieldnames(s);
for k = 1:numel(fn)
  F(n).(fn{k}) = s.(fn{k});
end
end

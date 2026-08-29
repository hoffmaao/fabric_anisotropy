%RUN_NEGIS_FABRIC Rathmann eigenvalue inversion for the NEGIS profiles.
%
% The NEGIS data are not in CSARP_polarimetric form, so for each frame in
% `targets` this first repackages its qlook HH/VV pair into that layout
% and then runs the ordinary fabric_task delta-k chain on it. Repackaging
% rather than reimplementing keeps the inversion identical to the one that
% produced the Thwaites and Ridge A results - the only thing that differs
% between the profiles is where dtau came from. A target that fails is
% warned about and skipped, so one bad frame does not cost the batch.
%
% Three season-specific repairs happen here and nowhere else:
%   1. Latitude/Longitude are rebuilt from the GPS file (the records time
%      sync failed for every segment except 20240618_01).
%   2. Traces the vehicle was stopped for are culled. They multilook to
%      HIGH coherence with random range phase, so they survive every
%      coherence test downstream and would drag whole blocks with them.
%   3. Surface is picked from HH power - the product ships Surface = NaN,
%      and ptt.surfaceReference needs it to place the reference band and
%      the waveform-combine seam.
%
% dtau comes from the SNAPHU-unwrapped phase, after Goldstein-Werner
% filtering of the complex interferogram. Both steps happen HERE, in the
% repackaging, because this season never went through polarimetric.m,
% which is where the Antarctic products get theirs.
%
% It used to run on delta-k, on the grounds that nothing else was
% available. That was wrong twice over: SNAPHU only needed running, and
% delta-k was actively costing accuracy. The EastGRIP comparison against
% Zeising et al. (2023) measured it - the joint inversion fed by delta-k
% scored RMS 0.411 against the core, WORSE than the bare fringe rate
% (0.280) and worse than the reference method (0.356), returning about a
% third of the fringe rate's amplitude above 650 m and flipping sign
% below. That is the stage-A suppression already diagnosed on Ridge A,
% where 15 MHz sub-bands are SNR-starved at depth and smoothing
% noise-dominated phasors biases the recovered angle toward zero,
% arriving at the inversion through dtau.
%
% There is still no coregistration for this season, so with unwrapped
% phase the per-block fringe constant is left to ptt.blockAverage and the
% level is set by the surface reference alone. phase_sign is forced to -1,
% the matched-filter convention settled by the Ridge A cell regression,
% because the sign regression needs row_offset to run.
%
% Launch on mem1 with (from the work root, whose code/ is on the path):
%   /opt/sw/matlab/2024b/bin/matlab -batch "run_negis_fabric"

season_dir = '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2';
gps_dir = '/cresis/dataproducts/opr_data/opr_support/gps/2024_Greenland_Ground2';
% <work> is derived from this script's own location - see fabric_paths.
[~, scratch] = fabric_paths();
season = '2024_Greenland_Ground2';
% Frames to repackage and invert. The first is the along-flow line beside
% the southeastern shear margin; the second stands 0.32 km off the
% EastGRIP borehole, so its fabric can be set against the core's own
% measurements. The borehole's own segment 20240618_01 carries no phase
% (incoherent decimation). Two lines pass NEARER - 20240621_01_010 at
% 0.13 km and 20240626_01_001 at 0.15 km - but both drive 0.8-1.6% of
% intervals onto the eigenvalue bound at |dlam| = 2/3 with 0.54-0.74 ns
% misfit, so this is the nearest line that inverts stably.
% Overridable so a caller can drive a different set through the SAME
% inversion without forking this script - egrip_zeising.m's nine
% borehole-proximal lines are run this way for the method comparison, and
% a forked copy would be free to drift in the repackaging, the cull or the
% surface pick, which is exactly what the comparison must hold fixed:
%   matlab -batch "targets = {'20240628_01',1; ...}; run_negis_fabric"
if ~exist('targets', 'var') || isempty(targets)
  targets = { '20240626_03', 1; '20240619_01', 1 };
end

MIN_STEP_M = 1.0;      % below this the vehicle was stopped
CROP_PRE = 0.5e-6;     % keep this much above the surface
CROP_POST = 22e-6;     % ...and this much below it
COH_WIN = [9 9];       % boxcar for the full-grid coherence

% Goldstein filtering + SNAPHU unwrapping, so the inversion can run off
% dtau_source = 'phase'. The binary is the one run_polarimetric.m uses for
% the Antarctic seasons, so this season is unwrapped by the same SNAPHU
% and the same call as the products it is compared against.
snaphu_en = true;
SNAPHU_BIN = '/kucresis/scratch/software/snaphu/bin/snaphu';
GOLD_ALPHA = 0.8;      % ptt.goldsteinFilter default; see branch_cuts.py
kmb = ones(COH_WIN, 'single') / prod(COH_WIN);   % shared multilook boxcar

code = fullfile(scratch, 'code');
addpath(code);
addpath(fullfile(code, 'opr_fabric'));
addpath(fullfile(code, 'opr_fabric', 'test', 'stubs'));

for ti = 1:size(targets, 1)
day_seg = targets{ti, 1};
frm = targets{ti, 2};

season_root = fullfile(scratch, season);
in_name = 'polarimetric_negis';
in_dir = fullfile(season_root, ['CSARP_' in_name], day_seg);
if ~exist(in_dir, 'dir'), mkdir(in_dir); end
in_fn = fullfile(in_dir, sprintf('Data_%s_%03d.mat', day_seg, frm));

%% ---- repackage ---------------------------------------------------------
% Reuse only a product that already carries what this run needs. Products
% written before the SNAPHU step have no snaphu_out_phase, and a bare
% exist() test would reuse one, let fabric_task fall back to the wrapped
% phase, and report a 'phase' run that never saw an unwrapped phase - the
% expensive failure being that it looks like a successful comparison.
reuse = exist(in_fn, 'file') == 2;
if reuse && snaphu_en
  % An existing product only counts if it carries a NON-EMPTY unwrapped
  % phase: H_snaphu's failure path stores [], and that must re-run rather
  % than be reused forever.
  reuse = false;
  if ismember('snaphu_out_phase', who('-file', in_fn))
    q = load(in_fn, 'snaphu_out_phase');
    reuse = ~isempty(q.snaphu_out_phase);
  end
end
if exist(in_fn, 'file') == 2
  if reuse
    fprintf('Reusing existing repackaged product %s\n', in_fn);
  else
    fprintf('Rebuilding %s (no unwrapped phase in the existing product)\n', ...
      in_fn);
  end
end
if ~reuse
try
  name = sprintf('Data_%s_%03d.mat', day_seg, frm);
  hh_fn = fullfile(season_dir, 'CSARP_qlook_HH', day_seg, name);
  vv_fn = fullfile(season_dir, 'CSARP_qlook_VV', day_seg, name);
  fprintf('Loading %s\n', hh_fn);
  H = load(hh_fn, 'Data', 'Time', 'GPS_time', 'Elevation', ...
    'param_qlook', 'param_records');
  V = load(vv_fn, 'Data', 'Time', 'GPS_time');
  assert(isequal(H.Time, V.Time) && isequal(H.GPS_time, V.GPS_time), ...
    'HH/VV axes differ');
  assert(~isreal(H.Data) && ~isreal(V.Data), 'Data is real (inc_dec > 0)');
  [Nt, Nx] = size(H.Data);
  fprintf('  %d samples x %d traces\n', Nt, Nx);

  % 1. positions from the GPS file
  g = load(fullfile(gps_dir, sprintf('gps_%s.mat', day_seg(1:8))), ...
    'sync_gps_time', 'sync_lat', 'sync_lon', 'sync_elev');
  [tu, iu] = unique(g.sync_gps_time(:));
  gt = H.GPS_time(:);
  Latitude = interp1(tu, g.sync_lat(iu), gt, 'linear', NaN).';
  Longitude = interp1(tu, g.sync_lon(iu), gt, 'linear', NaN).';
  fprintf('  repositioned lat %.4f..%.4f lon %.4f..%.4f\n', ...
    min(Latitude), max(Latitude), min(Longitude), max(Longitude));

  % 2. cull the stops. Traces the GPS file does not cover are dropped by
  % the same step but are a different failure, so they are counted apart -
  % a day file that misses the frame would otherwise read as a traverse
  % that never moved. A trace whose predecessor has no position has no
  % measurable step, so it is kept rather than called stopped.
  R = 6371e3;
  p = deg2rad(Latitude); dl = diff(deg2rad(Longitude)); dp = diff(p);
  a = sin(dp/2).^2 + cos(p(1:end-1)).*cos(p(2:end)).*sin(dl/2).^2;
  step = [Inf, 2*R*asin(sqrt(a))];
  located = isfinite(Latitude) & isfinite(Longitude);
  stopped = located & [false, located(1:end-1)] & step < MIN_STEP_M;
  keep = located & ~stopped;
  fprintf('  culling %d of %d traces: %d stopped, %d with no GPS coverage\n', ...
    nnz(~keep), Nx, nnz(stopped), nnz(~located));
  if any(~located)
    warning('run_negis_fabric:gpsCoverage', ...
      '%s frame %d: %d of %d traces (%.1f%%) fall outside the day GPS file', ...
      day_seg, frm, nnz(~located), Nx, 100*nnz(~located)/Nx);
  end
  if nnz(keep) < COH_WIN(2)
    % Warn and move on rather than abort: `targets` is a batch, and under
    % matlab -batch a bad day GPS file on the first target would otherwise
    % cost every later target too.
    warning('run_negis_fabric:tooFewTraces', ...
      ['%s frame %d: only %d traces survive the cull, fewer than the ' ...
       '%d-trace coherence boxcar. Check the day GPS file covers this ' ...
       'frame; skipping this target.'], day_seg, frm, nnz(keep), COH_WIN(2));
    clear H V;
    continue;
  end
  H.Data = H.Data(:, keep); V.Data = V.Data(:, keep);
  Latitude = Latitude(keep); Longitude = Longitude(keep);
  GPS_time = H.GPS_time(keep); Elevation = H.Elevation(keep);
  Nx = nnz(keep);

  % 3. surface from HH power: first sample above 5% of the trace peak,
  % which tracks the leading edge rather than the (broader) power maximum
  pw = abs(H.Data).^2;
  thr = 0.05 * max(pw, [], 1);
  Surface = nan(1, Nx);
  for k = 1:Nx
    j = find(pw(:,k) > thr(k), 1);
    if ~isempty(j), Surface(k) = H.Time(j); end
  end
  clear pw thr;
  fprintf('  surface %.3f..%.3f us (median %.3f)\n', min(Surface)*1e6, ...
    max(Surface)*1e6, median(Surface,'omitnan')*1e6);

  % crop the range window to the band the inversion can use; the deep half
  % of a 53 us record is noise and would only slow the ladder down
  s0 = median(Surface, 'omitnan');
  rows = H.Time >= s0 - CROP_PRE & H.Time <= s0 + CROP_POST;
  fprintf('  cropping range to %.2f..%.2f us (%d of %d samples)\n', ...
    H.Time(find(rows,1))*1e6, H.Time(find(rows,1,'last'))*1e6, nnz(rows), Nt);
  ref = single(H.Data(rows, :));
  sec = single(V.Data(rows, :));
  Time = H.Time(rows);
  clear H V;

  % full-grid coherence, boxcar over range and azimuth. Formed the same way
  % polarimetric.m does: multilook the cross product and both powers, then
  % divide - never smooth a per-pixel coherence.
  num = conv2(real(sec .* conj(ref)), kmb, 'same') ...
    + 1i*conv2(imag(sec .* conj(ref)), kmb, 'same');
  den = sqrt(conv2(abs(ref).^2, kmb, 'same') ...
    .* conv2(abs(sec).^2, kmb, 'same'));
  interferogram_coherence = single(abs(num) ./ den);
  interferogram_mlook = single(num);
  clear num den;
  fprintf('  coherence median %.3f\n', ...
    median(interferogram_coherence(isfinite(interferogram_coherence))));

  % Goldstein-filtered, SNAPHU-unwrapped phase. This is what lets the
  % inversion run off dtau_source = 'phase' instead of delta-k: the
  % EastGRIP method comparison showed delta-k's stage-A suppression
  % reaching the joint inversion through dtau and costing it a third of
  % the fringe-rate amplitude above 650 m, with a sign flip below.
  %
  % Filter BEFORE unwrapping, and filter the complex interferogram rather
  % than the phase. SNAPHU's cost is driven by residues, which are made by
  % noise-driven excursions between adjacent pixels; the Goldstein
  % reweighting removes most of them (measured on the Thwaites margin
  % frame in scripts/prototypes/branch_cuts.py) so the solution is not
  % dominated by branch cuts. The coherence handed to SNAPHU stays the
  % UNFILTERED one - it is the weight telling SNAPHU where to trust the
  % phase, and recomputing it from the filtered interferogram would report
  % the filter's own smoothing as data quality.
  snaphu_out_phase = [];
  if snaphu_en
    t2 = tic;
    ifg_filt = ptt.goldsteinFilter(interferogram_mlook, GOLD_ALPHA);
    fprintf('  Goldstein filter (alpha %.2f) %.1f min\n', GOLD_ALPHA, ...
      toc(t2)/60);
    % The try covers the SOLVE and nothing else. A tiled unwrap of this
    % grid is the longest step in the repackaging, so letting the catch see
    % a failure in the reporting below would throw away an hour of work and
    % - because the reuse test above keys on a non-empty phase - repeat it
    % on every later run.
    try
      snaphu_out_phase = H_snaphu(ifg_filt, interferogram_coherence, ...
        SNAPHU_BIN, COH_WIN);
    catch snerr
      % A missing binary or a SNAPHU that will not converge must not cost
      % the repackaging - the product is still usable via delta-k, and the
      % inversion warns for itself when snaphu_out_phase is absent.
      warning('run_negis_fabric:snaphuFailed', ...
        '%s frame %d: SNAPHU failed (%s: %s); saving without unwrapped phase', ...
        day_seg, frm, snerr.identifier, snerr.message);
      snaphu_out_phase = [];
    end
    if ~isempty(snaphu_out_phase)
      % max-min rather than range(): this script must run on a MATLAB
      % without the Statistics Toolbox, and a report is not worth a
      % dependency.
      lo_ph = min(snaphu_out_phase(:)); hi_ph = max(snaphu_out_phase(:));
      fprintf('  snaphu unwrapped: %.1f..%.1f rad (%.1f fringes)\n', ...
        lo_ph, hi_ph, (hi_ph - lo_ph)/(2*pi));
    end
    clear ifg_filt;
  end

  % param_records carries array.img_comb through to ptt.imgCombSeam, so the
  % 8 us waveform-combine boundary is masked exactly as on the other sites
  param_records = H_param_records_from(hh_fn);
  param_polarimetric = [];
  file_type = 'polarimetric';
  file_version = '1';

  fprintf('Saving %s\n', in_fn);
  save(in_fn, '-v7.3', 'ref', 'sec', 'interferogram_coherence', 'Time', ...
    'Surface', 'Latitude', 'Longitude', 'Elevation', 'GPS_time', ...
    'interferogram_mlook', 'snaphu_out_phase', ...
    'param_records', 'param_polarimetric', 'file_type', 'file_version');
  clear ref sec interferogram_coherence interferogram_mlook snaphu_out_phase;
catch err
  % Everything the repackaging can throw - a qlook product that is absent
  % or real-valued, mismatched HH/VV axes, a day GPS file that will not
  % load - is a property of THIS target, so it is warned about and skipped
  % like the cull guard above and the inversion guard below. A save that
  % died part way leaves a file the `exist` test would happily reuse next
  % run, so it is removed.
  warning('run_negis_fabric:repackageFailed', ...
    '%s frame %d: repackaging failed (%s: %s); skipping this target', ...
    day_seg, frm, err.identifier, err.message);
  clear H V ref sec interferogram_coherence interferogram_mlook;
  if exist(in_fn, 'file'), delete(in_fn); end
  continue;
end
end

%% ---- invert ------------------------------------------------------------
param = [];
param.day_seg = day_seg;
param.season_name = season;
param.radar_name = 'accum3';
param.load.frm = frm;
param.opr_file_lock = false;
param.stub_out_root = season_root;

pf = [];
pf.in_path = in_name;
pf.out_path = 'fabric_snaphu_negis';
pf.dtau_source = 'phase';
pf.use_snaphu_phase = true;
pf.inversion = 'joint';
pf.reg = 0.05;
pf.img = 0;
pf.out_file_exts = {'.png'};
pf.fc = 750e6;
pf.phase_sign = -1;          % no row_offset here; matched-filter convention
pf.coherence_threshold = 0.5;
pf.min_coverage = 0.3;
pf.block_size = 1000;
pf.num_intervals = 10;
pf.ref_twtt_offset = 50e-9;
pf.half_offset = 0;
pf.ptt = struct('H', 2000, 'bco_depth', 60, 'lam_z_sfc', 1/3, 'lam_z_bed', 1/3);
param.fabric = pf;

t1 = tic;
try
  ok = fabric_task(param);
  fprintf('\nfabric_task returned %d (%.1f min)\n', ok, toc(t1)/60);
catch err
  % Same reasoning as the cull guard above: one target that cannot invert
  % must not take the rest of the batch down with it.
  warning('run_negis_fabric:inversionFailed', ...
    '%s frame %d: fabric_task failed after %.1f min (%s: %s); continuing', ...
    day_seg, frm, toc(t1)/60, err.identifier, err.message);
end

end   % targets

function ph = H_snaphu(ifg, coh, bin, mlook_window)
% Unwrap a complex interferogram with SNAPHU, returning phase on the SAME
% grid. Byte layouts and the command line are taken from
% polarimetric_task.m verbatim, so the Antarctic products and this season
% are unwrapped by identical means: COMPLEX_DATA in (interleaved
% real/imag), FLOAT_DATA correlation, ALT_LINE_DATA out (magnitude and
% phase on alternating range lines).
if exist(bin, 'file') ~= 2
  error('run_negis_fabric:snaphuMissing', ...
    'SNAPHU binary not found at %s', bin);
end
fn = tempname();
fn_out = [fn '.out'];
fn_corr = [fn '.corr'];
c = onCleanup(@() H_rm({fn, fn_out, fn_corr}));

% SNAPHU wants a bounded interferogram; normalising by the global maximum
% is what polarimetric_task.m does and keeps the correlation file the only
% statement about quality.
x = double(ifg);
x(~isfinite(x)) = 0;
mx = max(abs(x(:)));
if mx <= 0
  error('run_negis_fabric:snaphuEmpty', 'interferogram is all zero');
end
x = x / mx;
fid = fopen(fn, 'wb');
if fid < 0, error('run_negis_fabric:snaphuIO', 'cannot write %s', fn); end
fwrite(fid, reshape([real(x(:)) imag(x(:))].', ...
  [size(x,1)*2 size(x,2)]), 'float32');
fclose(fid);

cc = double(abs(coh));
cc(cc > 1) = 1;
cc(cc < 0 | ~isfinite(cc)) = 0;
fid = fopen(fn_corr, 'wb');
if fid < 0, error('run_negis_fabric:snaphuIO', 'cannot write %s', fn_corr); end
fwrite(fid, cc, 'float32');
fclose(fid);

% Tiled, unlike polarimetric_task.m. Their Antarctic frames unwrap as a
% single tile; ours is 13501 x ~2000 after the 22 us crop at 1.667 ns
% sampling, and a single-tile solve on that 27M-pixel grid ran over an
% hour without finishing. Tiles are sized to ~2500 range bins x ~1200
% traces, which is the regime SNAPHU solves in seconds, and the overlap
% lets its secondary optimisation stitch them without leaving
% tile-boundary fringes.
%
% SNAPHU's grid is the TRANSPOSE of ours. The file above is written column
% by column with linelength = size(x,1), so a SNAPHU line is one of our
% traces and a SNAPHU sample is one of our range bins. `--tile` takes
% nrow ncol rowovrlp colovrlp in that frame, so the trace axis sets nrow
% and the range axis sets ncol - deriving them the other way round
% partitions 2000 traces into 5 and 13501 range bins into 2, and makes the
% row overlap half a tile rather than the intended third.
n_tile_row = max(1, round(size(x,2) / 1200));   % over our traces
n_tile_col = max(1, round(size(x,1) / 2500));   % over our range bins
ovr = min(200, floor(size(x,2) / n_tile_row / 3));
ovc = min(200, floor(size(x,1) / n_tile_col / 3));
tile = '';
if n_tile_row > 1 || n_tile_col > 1
  tile = sprintf(' --tile %d %d %d %d --nproc %d', n_tile_row, n_tile_col, ...
    ovr, ovc, min(8, n_tile_row*n_tile_col));
end
cmd = sprintf(['%s %s %d -c %s -s -C "NLOOKSRANGE %d" -C "NLOOKSAZ %d" ' ...
  '-C "CORRFILEFORMAT FLOAT_DATA"%s -v -o %s'], bin, fn, size(x,1), ...
  fn_corr, mlook_window(1), mlook_window(2), tile, fn_out);
fprintf('  running: %s\n', cmd);
[st, out] = system(cmd);
if st ~= 0
  error('run_negis_fabric:snaphuExit', ...
    'SNAPHU exited %d: %s', st, strtrim(out(max(1,end-400):end)));
end
if exist(fn_out, 'file') ~= 2
  error('run_negis_fabric:snaphuNoOutput', ...
    'SNAPHU wrote no output file');
end
fid = fopen(fn_out, 'rb');
raw = fread(fid, [size(x,1) 2*size(x,2)], 'float32');
fclose(fid);
if size(raw,2) < 2*size(x,2)
  error('run_negis_fabric:snaphuShort', ...
    'SNAPHU output has %d columns, expected %d', size(raw,2), 2*size(x,2));
end
ph = single(raw(:, 2:2:end));
end

function H_rm(files)
for k = 1:numel(files)
  if exist(files{k}, 'file') == 2, delete(files{k}); end
end
end

function pr = H_param_records_from(fn)
% param_records straight from the qlook product, so array.img_comb and
% img_comb_mult reach ptt.imgCombSeam unmodified.
%
% Both fields are optional downstream - imgCombSeam masks nothing without
% img_comb and defaults img_comb_mult to Inf - and where a qlook product
% keeps them differs by season (some carry a param_qlook.array copy, some
% only param_qlook.qlook). Every lookup is therefore guarded: a product
% that supplies neither has to report that, not abort the repackaging
% after the full HH/VV load, crop and coherence convolution have run.
d = load(fn, 'param_records', 'param_qlook');
pr = struct();
if H_has(d, {'param_records'}) && isstruct(d.param_records) ...
    && isscalar(d.param_records)
  pr = d.param_records;
end
if ~isfield(pr, 'array') || ~isstruct(pr.array)
  pr.array = struct();
end
if ~isfield(pr.array, 'img_comb') || isempty(pr.array.img_comb)
  if H_has(d, {'param_qlook','array'})
    pr.array = d.param_qlook.array;
  end
end
if ~isfield(pr.array, 'img_comb') || isempty(pr.array.img_comb)
  pr.array.img_comb = H_get(d, {'param_qlook','qlook','img_comb'});
end
if ~isfield(pr.array, 'img_comb_mult') || isempty(pr.array.img_comb_mult)
  pr.array.img_comb_mult = H_get(d, {'param_qlook','qlook','img_comb_mult'});
end
fprintf('  img_comb = %s, mult = %s\n', H_show(pr.array.img_comb), ...
  H_show(pr.array.img_comb_mult));
end

function tf = H_has(s, path)
tf = true;
for k = 1:numel(path)
  if ~isstruct(s) || ~isscalar(s) || ~isfield(s, path{k}), tf = false; return; end
  s = s.(path{k});
end
end

function v = H_get(s, path)
v = [];
if ~H_has(s, path), return; end
for k = 1:numel(path)
  s = s.(path{k});
end
v = s;
end

function str = H_show(v)
if isempty(v), str = '<absent>'; else, str = mat2str(v); end
end

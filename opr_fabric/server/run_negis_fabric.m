%RUN_NEGIS_FABRIC Rathmann eigenvalue inversion for the NEGIS profile.
%
% The NEGIS data are not in CSARP_polarimetric form, so this first
% repackages one qlook HH/VV pair into that layout and then runs the
% ordinary fabric_task delta-k chain on it. Repackaging rather than
% reimplementing keeps the inversion identical to the one that produced
% the Thwaites and Ridge A results - the only thing that differs between
% the three profiles is where dtau came from.
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
% delta-k is the only dtau source available here: there is no SNAPHU
% unwrapping and no coregistration for this season. With no row_offset the
% orientation regression cannot run, so phase_sign is forced to -1, the
% matched-filter convention settled by the Ridge A cell regression.
%
% Launch on mem1 with:
%   /opt/sw/matlab/2024b/bin/matlab -batch "run('/kucresis/scratch/hoffmana_sta/fabric/run_negis_fabric.m')"

season_dir = '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2';
gps_dir = '/cresis/dataproducts/opr_data/opr_support/gps/2024_Greenland_Ground2';
scratch = '/kucresis/scratch/hoffmana_sta/fabric';
season = '2024_Greenland_Ground2';
% Frames to repackage and invert. The first is the along-flow line beside
% the southeastern shear margin; the second stands 0.32 km off the
% EastGRIP borehole, so its fabric can be set against the core's own
% measurements. The borehole's own segment 20240618_01 carries no phase
% (incoherent decimation). Two lines pass NEARER - 20240621_01_010 at
% 0.13 km and 20240626_01_001 at 0.15 km - but both drive 0.8-1.6% of
% intervals onto the eigenvalue bound at |dlam| = 2/3 with 0.54-0.74 ns
% misfit, so this is the nearest line that inverts stably.
targets = { '20240626_03', 1; '20240619_01', 1 };

MIN_STEP_M = 1.0;      % below this the vehicle was stopped
CROP_PRE = 0.5e-6;     % keep this much above the surface
CROP_POST = 22e-6;     % ...and this much below it
COH_WIN = [9 9];       % boxcar for the full-grid coherence

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
if exist(in_fn, 'file')
  fprintf('Reusing existing repackaged product %s\n', in_fn);
else
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
    error('run_negis_fabric:tooFewTraces', ...
      ['%s frame %d: only %d traces survive the cull, fewer than the ' ...
       '%d-trace coherence boxcar. Check the day GPS file covers this ' ...
       'frame.'], day_seg, frm, nnz(keep), COH_WIN(2));
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
  k = ones(COH_WIN, 'single') / prod(COH_WIN);
  num = conv2(real(sec .* conj(ref)), k, 'same') ...
    + 1i*conv2(imag(sec .* conj(ref)), k, 'same');
  den = sqrt(conv2(abs(ref).^2, k, 'same') .* conv2(abs(sec).^2, k, 'same'));
  interferogram_coherence = single(abs(num) ./ den);
  clear num den k;
  fprintf('  coherence median %.3f\n', ...
    median(interferogram_coherence(isfinite(interferogram_coherence))));

  % param_records carries array.img_comb through to ptt.imgCombSeam, so the
  % 8 us waveform-combine boundary is masked exactly as on the other sites
  param_records = H_param_records_from(hh_fn);
  param_polarimetric = [];
  file_type = 'polarimetric';
  file_version = '1';

  fprintf('Saving %s\n', in_fn);
  save(in_fn, '-v7.3', 'ref', 'sec', 'interferogram_coherence', 'Time', ...
    'Surface', 'Latitude', 'Longitude', 'Elevation', 'GPS_time', ...
    'param_records', 'param_polarimetric', 'file_type', 'file_version');
  clear ref sec interferogram_coherence;
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
pf.out_path = 'fabric_deltak_negis';
pf.dtau_source = 'deltak';
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
ok = fabric_task(param);
fprintf('\nfabric_task returned %d (%.1f min)\n', ok, toc(t1)/60);

end   % targets

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

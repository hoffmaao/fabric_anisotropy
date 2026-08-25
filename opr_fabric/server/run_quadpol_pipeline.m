%RUN_QUADPOL_PIPELINE Coregister and invert one profile, end to end.
%
% The complete quad-pol workflow for a single frame, in one pass:
%
%   1. load all four standardphase channels on the product's own window
%   2. verify the source matches what the shipped product was built from
%   3. coregister VV, HV and VH onto HH with the OPR toolbox, using the
%      settings the product recorded - or load the coreg_cache from a
%      previous run, which turns a rerun from ~52 min into a few
%   4. invert the COREGISTERED channels twice: ptt.ershadiFabric (the
%      published chain, evaluated at the cross-pol minimum, which on this
%      system is antenna-locked) and ptt.quadpolFabricLS (the model fit
%      that takes the axis from the coherence field; a laterally-segmented
%      frame theta0 pass via ptt.quadpolFrameTheta, then per-block dlam
%      with the block's segment theta0 held via ptt.thetaProfileAt)
%   5. report against the LS frame fit - the two-azimuth solve is printed
%      as context only, not as a standard - and save
%
% Coregistration and inversion are deliberately in the SAME script. The
% coregistered images are ~370 MB a frame, so writing them out to be read
% back by a separate inversion step would cost more in I/O than the
% inversion itself, and would invite the two halves drifting apart on which
% window or which settings they used - which is exactly the mistake that
% produced the retracted result in 03d292e, where the inversion ran on
% channels the coregistration step had never touched.
%
% Runtime is dominated by step 3: ~71 min on a 6601 x 4346 frame, less in
% proportion for shorter ones. Run per profile.
%
%   matlab -batch "day_seg='20250108_02'; frm=9; run_quadpol_pipeline"
if ~exist('site_root', 'var') || isempty(site_root)
  site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
end
if ~exist('day_seg', 'var') || isempty(day_seg), day_seg = '20250108_02'; end
if ~exist('frm', 'var') || isempty(frm), frm = 9; end
if ~exist('out_dir', 'var') || isempty(out_dir)
  out_dir = '/kucresis/scratch/hoffmana_sta/fabric/stages/quadpol';
end
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
addpath('/kucresis/scratch/hoffmana_sta/fabric/code');

CHAN = {'hh','vv','hv','vh'};
NRW = 101;
Z_BAND = [200 1200];
% z_max is overridable (like day_seg/frm) for special runs that need the
% full record - e.g. the SCAR figures whose fabric panel must span the
% same depth range as the wrapped-phase panel. Non-default depths get
% their own cache and output names so they can NEVER clobber the
% standard 1500 m products the batches build and consume.
if ~exist('z_max', 'var') || isempty(z_max), z_max = 1500; end
Z_MAX = z_max;
if z_max ~= 1500
  ZTAG = sprintf('_z%d', round(z_max));
else
  ZTAG = '';
end
FC = 750e6;
PSI_STEP_DEG = 1;
% dlam search ceiling for the LS estimator. Its default 0.25 predates any
% site whose contrast approaches it: the EGRIP core sits at ~0.2-0.35,
% ABOVE the default, and a window whose true rate exceeds the cap either
% rails there and abstains or locks an aliased branch below it - the
% oscillating dlam, 0.3-0.6 residuals and 0.233 ceiling of the first
% EastGRIP validation frame. Synthetic (test_egrip_cap.m): at cap 0.25 a
% 0.32 fabric abstains on 100% of windows; at 0.45 it is recovered exactly
% at every block size. Overridable per run; qlook mode defaults high, the
% Antarctic sites keep the estimator default.
if exist('dlam_max', 'var') && ~isempty(dlam_max)
  DLAM_MAX = dlam_max; dlam_ov = true;
else
  DLAM_MAX = 0.25; dlam_ov = false;
end
% Traces per along-track block for the SECTION. 125 is what run_sections.m
% uses for the co-polarized Ridge A section, so the two sections have the
% same along-track sampling and can be read against each other cell for
% cell rather than approximately.
%
% OVERRIDABLE, because 125 traces is not a fixed LENGTH. Ridge A traces sit
% ~1 m apart so a block is ~125 m; the EastGRIP qlook traces are ~9 m apart
% so the same count is a 1.1 km block, and a whole frame becomes 16 of them.
% At a shear margin, where fabric varies over hundreds of metres, that
% averages genuinely different ice into one block - which depresses dlam and
% destabilises theta0 exactly as the first EastGRIP validation frame did.
% Set nblk_tr at the call site to match the block LENGTH, not the count.
if ~exist('nblk_tr', 'var') || isempty(nblk_tr)
  nblk_tr = 125;            % Ridge A's count; ~125 m at ~1 m trace spacing
  nblk_auto = true;
else
  nblk_auto = false;
end
NBLK_TR = nblk_tr;
% Along-track SEGMENT length for the first-pass theta0(z) fits, in metres.
% The frame pass re-fits theta0 per ~SEG_LEN_M segment so lateral fabric
% variation is resolved instead of pooled away - the pooled frame pass is
% the two-pass design's single point of failure on laterally-varying
% frames (test_egrip_blocks.m). Frames shorter than two segments keep the
% single frame fit unchanged. See ptt.quadpolFrameTheta.
if exist('seg_len_m', 'var') && ~isempty(seg_len_m)
  SEG_LEN_M = seg_len_m;
else
  SEG_LEN_M = 2000;
end
CACHE_COREG = true;   % keep the coregistered channels; see below
name = sprintf('Data_%s_%03d.mat', day_seg, frm);
t_all = tic;

%% 1-2. window, settings and source check, from the product
% Two acquisition families are handled here. The Antarctic ground seasons
% ship CSARP_polarimetric, which carries the coregistration settings, the
% range window, Surface and a `ref` copy of the source - everything the run
% needs. The 2024 Greenland (EastGRIP) season ships none of that: the
% channels are CSARP_qlook_{HH,VV,HV,VH}, there is no polarimetric product
% at all, and Surface is NaN on every frame. QLOOK MODE below supplies each
% missing piece explicitly rather than silently defaulting, so a frame that
% took the fallback path says so in its own log.
pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
if exist(pol_fn, 'file') ~= 2
  % The 2022/2023 seasons shipped the polarimetric product only in its
  % SNAPHU-unwrapped form. It carries the same param_polarimetric (window
  % and coregistration settings), Surface and ref, so it serves the same
  % role here.
  alt = fullfile(site_root, 'CSARP_polarimetric_unwrap', day_seg, name);
  if exist(alt, 'file') == 2
    pol_fn = alt;
  end
end
qlook_mode = exist(pol_fn, 'file') ~= 2;
if qlook_mode
  CHAN_DIR = 'CSARP_qlook_%s';
  qfn = fullfile(site_root, sprintf(CHAN_DIR, 'HH'), day_seg, name);
  if exist(qfn, 'file') ~= 2
    error('run_quadpol_pipeline:noProduct', ...
      'neither a polarimetric product nor %s', qfn);
  end
  P = load(qfn, 'Time', 'Surface', 'Latitude', 'Longitude', 'GPS_time');
  % QLOOK COORDINATES THAT ARE NOT ON THE ICE SHEET are rebuilt from the
  % season's reference trajectory. The 2024 Greenland qlook products carry
  % a Latitude/Longitude pair that is not the geographic position of the
  % traverse: frame 20240619_01_001 reads 0.196..0.229 N, 1.738..1.915 E -
  % the Gulf of Guinea - where CSARP_reference_trajectory/ref_<seg>.mat
  % puts the same GPS times at 75.507..75.628 N, -35.967..-35.603 E on
  % 2700 m of ice, which is EastGRIP. Everything geometric is computed
  % from these numbers, so using them as shipped is not a cosmetic error:
  % measured on that frame, the along-track length inflates 3.3x (20.0 km
  % against a true 6.1 km, so the trace spacing reads 9.07 m instead of
  % 2.83 m and the auto block length, the stationary cull and the segment
  % boundaries all scale with it), and because a degree of longitude
  % subtends 111 km at the equator against 27.6 km at 75.6 N, the track
  % azimuth comes out 98.8 deg instead of 127.1 - a 28 deg error carried
  % straight into every theta0_geo this season reports.
  %
  % THE TRIGGER IS THE DEFECT, NOT THE PRODUCT FAMILY. Every survey this
  % code serves is polar, so shipped coordinates that put the frame below
  % 60 deg of latitude - or that carry no position at all - are the
  % failure described above and nothing else. Gating on qlook_mode instead
  % would block a future qlook season whose coordinates are correct on a
  % reference trajectory it has no need of.
  %
  % Where a rebuild IS required it is an interpolation of the reference
  % trajectory onto the frame's own GPS_time, and it is NOT optional and
  % NOT silently skipped: a missing reference file, or a frame whose GPS
  % times fall outside it, errors out, because the alternative is a
  % product that looks fine and is wrong by tens of degrees.
  la_qk = P.Latitude(:).'; lo_qk = P.Longitude(:).';
  POLAR_LAT_MIN = 60;
  n_polar = nnz(isfinite(la_qk) & abs(la_qk) >= POLAR_LAT_MIN);
  if n_polar >= 0.5 * numel(la_qk)
    fprintf(['shipped qlook coordinates are polar (%.3f..%.3f N ' ...
      '%.3f..%.3f E); kept as shipped\n'], min(la_qk), max(la_qk), ...
      min(lo_qk), max(lo_qk));
  else
    ref_fn = fullfile(site_root, 'CSARP_reference_trajectory', ...
      sprintf('ref_%s.mat', day_seg));
    if exist(ref_fn, 'file') ~= 2
      error('run_quadpol_pipeline:noRefTraj', ...
        ['qlook coordinates place %d of %d traces off any ice sheet and ' ...
        '%s is missing; the frame cannot be positioned'], ...
        numel(la_qk) - n_polar, numel(la_qk), ref_fn);
    end
    RT = load(ref_fn, 'gps_time', 'lat', 'lon');
    [rt_gps, rt_ord] = sort(RT.gps_time(:));
    rt_lat = RT.lat(:); rt_lat = rt_lat(rt_ord);
    rt_lon = RT.lon(:); rt_lon = rt_lon(rt_ord);
    % interp1 rejects repeated sample points outright, so a reference
    % trajectory that logged one GPS time twice would abort the frame with
    % a generic grid-vector message instead of one of the checks here.
    [rt_gps, rt_uniq] = unique(rt_gps, 'stable');
    rt_lat = rt_lat(rt_uniq); rt_lon = rt_lon(rt_uniq);
    if numel(rt_gps) < 2
      error('run_quadpol_pipeline:refTrajShort', ...
        '%s holds %d distinct GPS times; it cannot position a frame', ...
        ref_fn, numel(rt_gps));
    end
    q_gps = P.GPS_time(:);
    if min(q_gps) < rt_gps(1) - 1 || max(q_gps) > rt_gps(end) + 1
      error('run_quadpol_pipeline:refTrajRange', ...
        ['frame GPS times %.1f..%.1f fall outside the reference trajectory ' ...
        '%.1f..%.1f'], min(q_gps), max(q_gps), rt_gps(1), rt_gps(end));
    end
    P.Latitude = interp1(rt_gps, rt_lat, q_gps, 'linear').';
    P.Longitude = interp1(rt_gps, rt_lon, q_gps, 'linear').';
    if any(~isfinite(P.Latitude)) || any(~isfinite(P.Longitude))
      error('run_quadpol_pipeline:refTrajGap', ...
        'reference trajectory left %d traces unpositioned', ...
        nnz(~isfinite(P.Latitude) | ~isfinite(P.Longitude)));
    end
    fprintf(['trajectory rebuilt from %s: shipped %.3f..%.3f N %.3f..%.3f E ' ...
      '-> %.3f..%.3f N %.3f..%.3f E\n'], ...
      sprintf('ref_%s.mat', day_seg), min(la_qk), max(la_qk), min(lo_qk), ...
      max(lo_qk), min(P.Latitude), max(P.Latitude), min(P.Longitude), ...
      max(P.Longitude));
    clear RT rt_gps rt_lat rt_lon rt_ord rt_uniq q_gps;
  end
  clear la_qk lo_qk;
  % Coregistration defaults, copied from what the Antarctic products
  % recorded (Ridge A frame 20250108_02_009), so both families are aligned
  % by the same tiling and search rather than by whatever a toolbox default
  % happens to be this release.
  CO = struct('Tt', 101, 'Tx', 301, 'overlap_t', 50, 'overlap_x', 150, ...
    'search_t', 5, 'search_x', 5, 'one_dim_search_en', 1);
  fprintf('=== %s_%03d === QLOOK MODE (no polarimetric product)\n', ...
    day_seg, frm);
  if ~dlam_ov, DLAM_MAX = 0.45; end
  fprintf('dlam search ceiling %.2f\n', DLAM_MAX);
  % Surface is NaN in these products, so the depth zero is picked here -
  % first sample above 5%% of the trace peak, the leading edge rather than
  % the broader power maximum. Same rule as negis_interferogram.m and
  % run_quadpol_frame.m, so every EastGRIP product shares one depth zero.
  % Picked from a SUBSET of traces read through matfile: the full record is
  % 32073 samples by several thousand traces and reading it whole, four
  % times, would cost gigabytes before the window is even known.
  surf_t = median(P.Surface(:), 'omitnan');
  if ~isfinite(surf_t)
    mf = matfile(qfn);
    [ntq, nxq] = size(mf, 'Data');
    % a matfile range must be evenly spaced, so stride rather than
    % linspace: an irregular index list is rejected outright
    stride = max(1, floor(nxq / 200));
    js = 1:stride:nxq;
    Dsub = mf.Data(:, js);
    pw = abs(double(Dsub)).^2;
    thr = 0.05 * max(pw, [], 1);
    st = nan(1, numel(js));
    for j = 1:numel(js)
      i0 = find(pw(:, j) > thr(j), 1);
      if ~isempty(i0), st(j) = P.Time(i0); end
    end
    surf_t = median(st, 'omitnan');
    fprintf('Surface all NaN; leading edge picked at %.3f us from %d traces\n', ...
      surf_t*1e6, numel(js));
    clear Dsub pw thr mf;
  end
  if ~isfinite(surf_t)
    error('run_quadpol_pipeline:surface', 'could not pick a surface');
  end
  % Range window from the surface and the depth ceiling, not from a
  % recorded rbin pair: with no polarimetric product there is none, and the
  % 53 us record runs far past any depth being inverted.
  C0q = 299792458; C_ICEq = C0q/sqrt(3.171);
  t_max = surf_t + 2*(z_max + 50)/C_ICEq;
  r0 = max(1, find(P.Time >= surf_t - 20e-9, 1));
  r1 = min(numel(P.Time), find(P.Time <= t_max, 1, 'last'));
  fprintf('window rbin %d:%d of %d (surface %.3f us, z_max %.0f m)\n', ...
    r0, r1, numel(P.Time), surf_t*1e6, z_max);
else
  CHAN_DIR = 'CSARP_standardphase_%s';
  w = {whos('-file', pol_fn).name};
  sel = intersect({'param_polarimetric','Time','Surface','Latitude', ...
    'Longitude','ref','sec_reg'}, w);
  P = load(pol_fn, sel{:});
  pp = P.param_polarimetric;
  if isfield(pp, 'polarimetric'), pp = pp.polarimetric; end
  r0 = pp.min_rbin; r1 = pp.max_rbin;
  co = pp.coregistration;
  CO = struct('Tt', co.Tt, 'Tx', co.Tx, 'overlap_t', co.overlap_t, ...
    'overlap_x', co.overlap_x, 'search_t', co.search_t, ...
    'search_x', co.search_x, 'one_dim_search_en', co.one_dim_search_en);
  fprintf('=== %s_%03d ===\n', day_seg, frm);
  surf_t = [];
end
% --- FIT THE WINDOW TO THE FRAME THAT ACTUALLY EXISTS.
% The window and the coregistration tiling are both read from the product
% and on short frames both can exceed it, which cost 21 of the 61 frames
% in the 20 Aug Antarctic run. The window is settled here because r0/r1
% are needed to read Data at all; the tiling waits for section 2c, where
% the trace axis is final.
%
% The channel sizes are taken from the file headers rather than by
% loading, so this costs nothing on the frames where it changes nothing.
navail = inf;
for k = 1:4
  fnk = fullfile(site_root, sprintf(CHAN_DIR, upper(CHAN{k})), day_seg, name);
  if exist(fnk, 'file') ~= 2
    error('run_quadpol_pipeline:missing', 'missing %s channel', CHAN{k});
  end
  wi = whos('-file', fnk, 'Data');
  if ~isempty(wi) && numel(wi(1).size) >= 2
    navail = min(navail, wi(1).size(1));
  end
end
% (1) RANGE. The thin-ice seasons (Eastwind, McMurdo: 200-300 m of ice)
% record a max_rbin from a deeper configuration than their own data - e.g.
% 20240202_01_001 asks for 11000 samples and holds 4151. Clamping keeps
% every sample that exists instead of discarding a 19 km line over a
% bookkeeping mismatch; only a window starting past the end is fatal.
if isfinite(navail) && r0 >= navail
  error('run_quadpol_pipeline:window', ...
    'window starts at %d but only %d samples exist', r0, navail);
end
if isfinite(navail) && r1 > navail
  warning('run_quadpol_pipeline:windowClamped', ...
    'recorded window %d:%d exceeds the %d samples present; clamped to %d:%d', ...
    r0, r1, navail, r0, navail);
  r1 = navail;
end
% (2) ALONG TRACK is fitted further down, once the trace axis is final:
% the tiling has to match the traces coregistration is actually handed,
% and in qlook mode the stationary cull below removes about a third of
% them. Only the range clamp can be settled here, because r0/r1 are
% needed to read Data at all.
fprintf('window rbin %d:%d\n', r0, r1);

S = struct();
for k = 1:4
  fn = fullfile(site_root, sprintf(CHAN_DIR, upper(CHAN{k})), day_seg, name);
  if exist(fn, 'file') ~= 2
    error('run_quadpol_pipeline:missing', 'missing %s channel', CHAN{k});
  end
  q = load(fn, 'Data', 'Time');
  if r1 > size(q.Data, 1)
    error('run_quadpol_pipeline:window', '%s too short for the window', ...
      CHAN{k});
  end
  S.(CHAN{k}) = q.Data(r0:r1, :);
  if k == 1, Tv = q.Time(r0:r1); end
  if qlook_mode && isreal(S.(CHAN{k}))
    error('run_quadpol_pipeline:notComplex', ...
      ['%s channel is REAL, not complex - the coherence fit needs the ' ...
       'full scattering matrix. Segment 20240618_01 of the EastGRIP ' ...
       'season is a real-only setup day and must be skipped.'], CHAN{k});
  end
end
[Nt, Nx] = size(S.hh);
fprintf('%d samples x %d traces\n', Nt, Nx);
if ~qlook_mode && isfield(P, 'ref')
  A = P.ref; if isstruct(A), A = A.Data; end
  if isequal(size(A), size(S.hh))
    rel = max(abs(A(:) - S.hh(:))) / max(abs(A(:)));
    fprintf('source check: max rel diff %.3g%s\n', rel, ...
      H_tag(rel < 1e-6, ' (IDENTICAL)', ' <-- NOT the product source'));
  end
end

%% 2b. cull the traces the traverse stopped for, BEFORE anything else
% A third of the EastGRIP traces have no along-track step. Coregistration
% tiles and every along-track block would otherwise average a stationary
% dwell as though it were distance, which both inflates the apparent
% coherence and puts blocks in the wrong place. Culling first keeps Data,
% the coordinates and the block geometry on one trace axis. Same rule as
% negis_interferogram.m: a trace whose predecessor has no position has no
% measurable step and is kept rather than called stopped.
if qlook_mode
  MIN_STEP_M = 0.5;
  la_all = P.Latitude(:).'; lo_all = P.Longitude(:).';
  Re = 6371e3;
  ph = deg2rad(la_all); dlo = diff(deg2rad(lo_all)); dph = diff(ph);
  aa = sin(dph/2).^2 + cos(ph(1:end-1)).*cos(ph(2:end)).*sin(dlo/2).^2;
  step = [Inf, 2*Re*asin(sqrt(aa))];
  located = isfinite(la_all) & isfinite(lo_all);
  stopped = located & [false, located(1:end-1)] & step < MIN_STEP_M;
  keep_tr = located & ~stopped;
  fprintf('stationary cull: %d of %d traces dropped (%d stopped, %d unlocated)\n', ...
    nnz(~keep_tr), Nx, nnz(stopped), nnz(~located));
  if nnz(keep_tr) < 64
    error('run_quadpol_pipeline:tooFewTraces', ...
      'only %d traces survive the cull', nnz(keep_tr));
  end
  for k = 1:4, S.(CHAN{k}) = S.(CHAN{k})(:, keep_tr); end
  P.Latitude = la_all(keep_tr); P.Longitude = lo_all(keep_tr);
  Nx = nnz(keep_tr);
end

%% 2c. FIT THE TILING TO THE TRACE AXIS COREGISTRATION WILL ACTUALLY SEE.
% coregistration() interpolates its per-tile offsets and needs at least
% FOUR tiles along track, i.e. 3*(Tx - overlap_x) + Tx <= Nx. That
% threshold is measured, not assumed: on this survey 773-trace frames
% coregister and 744-trace ones die indexing row_offset_full, and
% 3*151 + 301 = 754 falls exactly between them. Short frames therefore get
% a proportionally finer tiling - more, smaller tiles - rather than
% failing. Frames that already satisfy the condition are untouched, so
% every previously-processed frame is bit-identical.
%
% THIS RUNS AFTER THE CULL, and must. The count that matters is the width
% of the S handed to ptt.coregisterChannels, not the width recorded in the
% channel file headers: the cull above drops about a third of the EastGRIP
% traces, so a raw 900-trace qlook frame clears 754 on its header and then
% arrives at coregistration with ~600 - precisely the failure this guard
% exists to prevent. Outside qlook mode there is no cull and the two counts
% are the same, which is why the Antarctic run never exposed the gap.
need_x = 3*(CO.Tx - CO.overlap_x) + CO.Tx;
if Nx < need_x
  Tx0 = CO.Tx; ov0 = CO.overlap_x;
  while (3*(CO.Tx - CO.overlap_x) + CO.Tx) > Nx && CO.Tx > 16
    CO.Tx = floor(CO.Tx * 0.9);
    CO.overlap_x = floor(CO.Tx / 2);
  end
  if (3*(CO.Tx - CO.overlap_x) + CO.Tx) > Nx
    error('run_quadpol_pipeline:tooShort', ...
      ['%d traces cannot carry four coregistration tiles even at the ' ...
      'minimum tile width; frame is too short to coregister'], Nx);
  end
  warning('run_quadpol_pipeline:tilingAdapted', ...
    ['%d traces is short for tiling Tx %d ov %d (needs %d); adapted to ' ...
    'Tx %d ov %d'], Nx, Tx0, ov0, need_x, CO.Tx, CO.overlap_x);
end
fprintf('tiling Tt %d Tx %d ov %d/%d over %d traces\n', CO.Tt, CO.Tx, ...
  CO.overlap_t, CO.overlap_x, Nx);

%% depth axis
C0 = 299792458; C_ICE = C0/sqrt(3.171);
if qlook_mode
  st = surf_t;              % picked above; Surface is NaN in these products
else
  st = median(P.Surface(:), 'omitnan');
end
z = (Tv(:) - st) * C_ICE / 2;
keep = z >= 0 & z <= Z_MAX;
for k = 1:4, S.(CHAN{k}) = S.(CHAN{k})(keep, :); end
z = z(keep);
% The reporting band is clamped to the record that actually exists. The
% fixed 200-1200 m was written for the deep Antarctic sites; at the thin-ice
% seasons it reaches past the end of the sounding - Eastwind holds 561 m and
% McMurdo 998 - so every console figure quoted from it was averaging in
% depths the frame does not have. The saved arrays and the figures carry
% their own per-site bands, but this one is not report-only: `band` also
% selects the samples the before/after coregistration coherences are
% measured over, and those are saved into res. So it is clamped at BOTH
% ends - a record shallower than the 200 m floor would otherwise leave the
% band empty and write NaN coherences without a word - and a record that
% cannot carry any band at all is an error rather than a page of NaN.
z_band_asked = Z_BAND;
if z(end) < Z_BAND(2), Z_BAND(2) = z(end); end
if Z_BAND(1) >= Z_BAND(2), Z_BAND(1) = z(1); end
if ~isequal(Z_BAND, z_band_asked)
  fprintf('report band %.0f-%.0f m clamped to the record: %.0f-%.0f m\n', ...
    z_band_asked(1), z_band_asked(2), Z_BAND(1), Z_BAND(2));
end
band = z > Z_BAND(1) & z < Z_BAND(2);
if ~any(band)
  error('run_quadpol_pipeline:reportBand', ...
    ['the %.0f..%.0f m record holds no sample inside the reporting band ' ...
    '%.0f-%.0f m; every coherence and reported figure would be NaN'], ...
    z(1), z(end), Z_BAND(1), Z_BAND(2));
end
kr = ones(NRW,1)/NRW;
% Positions, needed by the section loop below as well as by the reporting,
% so they are defined once here rather than after the first use.
la = P.Latitude(:); lo = P.Longitude(:);

% --- BLOCK SIZE IS A LENGTH, AND IS SET FOR EVERY SEASON.
% A section block must cover the same ice everywhere or block-level
% quantities are not comparable between surveys. 125 traces is ~125 m only
% where traces sit ~1 m apart; this sizing used to run in qlook mode alone,
% so EastGRIP was corrected to 122 m while EASTWIND - an Antarctic season
% at 2.43 m spacing - silently kept 125 traces and produced 304 m blocks,
% 2.4x every other site, which is exactly the lateral averaging the
% segmented pass exists to avoid. Deriving it from the measured spacing
% covers both paths with one rule; an explicit nblk_tr at the call site
% still wins.
BLK_TARGET_M = 125;      % Ridge A's block, the length every site matches
BLK_TOL = 0.25;          % how far off before it is worth re-cutting
if nblk_auto
  R_E = 6371000;
  dph_s = deg2rad(diff(la(:)));
  dlo_s = deg2rad(diff(lo(:)));
  aa_s = sin(dph_s/2).^2 ...
    + cos(deg2rad(la(1:end-1))) .* cos(deg2rad(la(2:end))) .* sin(dlo_s/2).^2;
  sp = 2 * R_E * asin(min(1, sqrt(aa_s)));
  sp = median(sp(isfinite(sp) & sp > 0 & sp < 100));
  % Only RE-CUT when the default count spans materially the wrong length.
  % Without the tolerance this rounds 125 to 126 at every ~1 m site - the
  % measured spacing is 0.996 m, not exactly 1 - which would rewrite the
  % block boundaries of 107 of 130 finished frames to move a block by one
  % metre. That is churn, not a correction: it would invalidate the Ridge A
  % control and the Thwaites dose-response to no purpose. A quarter is wide
  % enough to leave every ~1 m season exactly as it was and still catch
  % Eastwind's 304 m and EastGRIP's 346 m.
  if isfinite(sp) && sp > 0 && abs(nblk_tr*sp - BLK_TARGET_M) > BLK_TOL*BLK_TARGET_M
    NBLK_TR = max(8, round(BLK_TARGET_M / max(sp, 0.5)));
    fprintf(['block size RE-CUT: %d traces (~%.0f m at %.2f m spacing); ' ...
      'the default %d would span %.0f m\n'], NBLK_TR, NBLK_TR*sp, sp, ...
      nblk_tr, nblk_tr*sp);
  elseif isfinite(sp) && sp > 0
    fprintf('block size: %d traces (~%.0f m at %.2f m spacing)\n', ...
      NBLK_TR, NBLK_TR * sp, sp);
  end
end
fprintf('depth window %.0f..%.0f m (%d samples)\n', z(1), z(end), numel(z));

%% 3. coregistration, or the cache of a previous run
% Coregistration is ~50 min and the inversions are minutes, so a rerun for
% a different estimator or setting should never pay it again. The cache is
% keyed by the settings that produced it; any mismatch falls through to a
% fresh coregistration.
pairs = {'hh','vv'; 'hh','hv'; 'hh','vh'; 'hv','vh'};
cache_fn = fullfile(out_dir, 'coreg_cache', ...
  sprintf('creg_%s_%03d%s.mat', day_seg, frm, ZTAG));
from_cache = false;
if exist(cache_fn, 'file') == 2
  try
    % The post-cull trace count is part of the key, so a changed cull rule
    % can never reuse a cache whose trace axis no longer matches the
    % freshly culled coordinates. It is read from the CACHED CHANNEL
    % ITSELF rather than from a stored Nx field: the array's own second
    % dimension is the trace axis the cache actually has, it exists in
    % every cache ever written, and it cannot drift from the data the way
    % a separately-saved scalar can. The earlier stored-field version had
    % a compatibility escape hatch - a missing Nx counted as valid - and
    % that hatch is exactly what let a 2058-trace EastGRIP cache load
    % against a corrected 2055-trace cull, indexing past the coordinate
    % arrays inside the frame pass. Nx is still written for provenance.
    cq = load(cache_fn, 'z', 'CO', 'r0', 'r1');
    cache_nx = NaN;
    try
      wi = whos('-file', cache_fn, 'hh');
      if ~isempty(wi) && numel(wi(1).size) >= 2
        cache_nx = wi(1).size(2);
      end
    catch
      cache_nx = NaN;
    end
    nx_ok = isfinite(cache_nx) && isequal(double(cache_nx), double(Nx));
    if isequal(cq.CO, CO) && cq.r0 == r0 && cq.r1 == r1 && nx_ok ...
        && numel(cq.z) == numel(z) && max(abs(cq.z(:) - z(:))) < 1e-6
      cq = load(cache_fn, 'hh', 'vv', 'hv', 'vh');
      T = struct('hh', double(cq.hh), 'vv', double(cq.vv), ...
        'hv', double(cq.hv), 'vh', double(cq.vh));
      clear cq;
      from_cache = true;
      fprintf('\ncoregistered channels loaded from %s\n', cache_fn);
    else
      fprintf('\ncache %s does not match the settings; re-coregistering\n', ...
        cache_fn);
      clear cq;
    end
  catch err
    fprintf('\ncache %s is unreadable (%s); re-coregistering\n', ...
      cache_fn, err.message);
    clear cq;
    from_cache = false;
  end
end
if from_cache
  t_coreg = 0;
  coh_before = nan(size(pairs,1),1);
  coh_after = nan(size(pairs,1),1);
  cinfo = struct('row_med', struct('vv', NaN, 'hv', NaN, 'vh', NaN));
else
  coh_before = zeros(size(pairs,1),1);
  for p = 1:size(pairs,1)
    coh_before(p) = H_coh(S.(pairs{p,1}), S.(pairs{p,2}), kr, band);
  end
  t0 = tic;
  [T, cinfo] = ptt.coregisterChannels(S, 'hh', CO);
  t_coreg = toc(t0)/60;
  coh_after = zeros(size(pairs,1),1);
  for p = 1:size(pairs,1)
    coh_after(p) = H_coh(T.(pairs{p,1}), T.(pairs{p,2}), kr, band);
  end
  fprintf('\ncoregistration %.1f min\n', t_coreg);
  fprintf('%-8s %8s %8s\n', 'pair', 'before', 'after');
  for p = 1:size(pairs,1)
    fprintf('%-8s %8.3f %8.3f\n', ...
      [upper(pairs{p,1}) '-' upper(pairs{p,2})], coh_before(p), coh_after(p));
  end
  fprintf('row offsets: VV %+.3f  HV %+.3f  VH %+.3f\n', ...
    cinfo.row_med.vv, cinfo.row_med.hv, cinfo.row_med.vh);
end

%% 3b. cache the coregistered channels
% ~730 MB a frame, which is worth it: coregistration is 24-56 min and the
% inversion is 6 seconds, so anything that needs re-inverting - a different
% window, a different block size, a bug - would otherwise pay the whole
% coregistration again. The cache is keyed by frame and by the settings
% that produced it.
if CACHE_COREG && ~from_cache
  cdir = fullfile(out_dir, 'coreg_cache');
  if exist(cdir, 'dir') ~= 7, mkdir(cdir); end
  hh = single(T.hh); vv = single(T.vv); hv = single(T.hv); vh = single(T.vh);
  % write to a temp name in the same directory and rename into place, so a
  % worker killed mid-write can never leave a truncated cache behind
  tmp_fn = [cache_fn '.tmp'];
  save(tmp_fn, '-v7.3', 'hh', 'vv', 'hv', 'vh', 'z', 'CO', 'r0', 'r1', 'Nx');
  movefile(tmp_fn, cache_fn);
  clear hh vv hv vh;
  fprintf('cached coregistered channels -> %s\n', cache_fn);
end

%% 3c. coreg-only mode, for building the caches in bulk
% The caches are estimator-independent, so the expensive step can run for
% the whole survey before the estimator settings are final: with
% coreg_only=true the script stops here and leaves the inversions to a
% later rerun, which the cache makes cheap.
if exist('coreg_only', 'var') && isequal(coreg_only, true)
  fprintf('coreg_only: stopping after the cache (total %.1f min)\n', ...
    toc(t_all)/60);
  return;
end

%% 4. inversion, on the coregistered channels
t0 = tic;
out = ptt.ershadiFabric(T, z, struct('fc', FC, 'psi_step_deg', PSI_STEP_DEG, ...
  'win_m', 30, 'grad_win_m', 25, 'coh_min', 0.4, 'deramped', true));
% and on the UNCOREGISTERED channels, so the effect of step 3 on the
% ANSWER - not just on the coherence - is measured rather than assumed
out_raw = ptt.ershadiFabric(S, z, struct('fc', FC, ...
  'psi_step_deg', PSI_STEP_DEG, 'win_m', 30, 'grad_win_m', 25, ...
  'coh_min', 0.4, 'deramped', true));
fprintf('\nershadi inversion %.1f min\n', toc(t0)/60);

%% 4a. the LS inversion: segmented frame theta0, then dlam with theta0 fixed
% ershadiFabric evaluates Psi at the cross-polarized minimum, which on
% this system is antenna-locked (89.6 +- 1.9 deg over 1790 blocks) by the
% flat -3.6 dB cross-pol pedestal, i.e. ~25 deg off the true axis. That is
% what put the odd-pi coherence-null bands and the fringe-periodic
% roughness in the section, and scaled dlam by ~cos 2*25 deg. The LS fit
% takes the axis from the coherence field itself and models the nulls
% instead of gating on them; see ptt.quadpolFabricLS and test_quadpol_ls.
t0 = tic;
p0 = deg2rad(la(1)); p1 = deg2rad(la(end));
dl = deg2rad(lo(end) - lo(1));
track_az = mod(rad2deg(atan2(sin(dl)*cos(p1), ...
  cos(p0)*sin(p1) - sin(p0)*cos(p1)*cos(dl))), 180);

% Per-trace heading, for the curving-line path and the per-block theta0
% handoff. On a CURVE the fabric rotates through the antennas along the
% drive: antenna-frame moment averaging smears it away (dlam collapses,
% theta0 is meaningless, and track_az itself stops meaning anything),
% while the antenna-fixed pedestal adds coherently. So: pedestal is
% calibrated from the antenna-frame pass either way; theta0/dlam on a
% curved frame are fitted in the GEOGRAPHIC frame from per-sub-block
% moments rotated north-referenced, with the pedestal entering as a
% precomputed field mixed over the measured heading distribution
% (rotation carries sin2(phi) to sin2(psi - h)); and every block inherits
% theta0 converted through its OWN heading rather than the frame's.
% test_quadpol_curved.m validates the geographic path on a 90-deg arc.
SMH = 51;
ih0 = 1:(Nx-SMH); ih1 = (1+SMH):Nx;
ph0h = deg2rad(la(ih0)); ph1h = deg2rad(la(ih1));
dlhh = deg2rad(lo(ih1) - lo(ih0));
az_tr = mod(rad2deg(atan2(sin(dlhh).*cos(ph1h), ...
  cos(ph0h).*sin(ph1h) - sin(ph0h).*cos(ph1h).*cos(dlhh))), 180);
% Interpolate the DOUBLED-ANGLE PHASOR, never the mod-180 angle: linear
% interpolation across the 0/180 wrap sweeps through ~90 deg on lines
% heading near north. Query points are clamped to the sample range (the
% same rule the theta0 handoff follows) because a linear extrapolation
% of a phasor can pass near zero.
xi = (1:numel(az_tr)).' + SMH/2;
ph2h = interp1(xi, exp(2i*deg2rad(az_tr)), ...
  min(max((1:Nx).', xi(1)), xi(end)), 'linear');
az_tr = mod(rad2deg(angle(ph2h))/2, 180);
% Along-track distance, so segment boundaries are cut in metres and
% "2 km" means 2 km at every site regardless of trace spacing.
R_E = 6371000;
dph_x = deg2rad(diff(la(:)));
dlo_x = deg2rad(diff(lo(:)));
aa_x = sin(dph_x/2).^2 ...
  + cos(deg2rad(la(1:end-1))) .* cos(deg2rad(la(2:end))) .* sin(dlo_x/2).^2;
x_along = [0; cumsum(2 * R_E * asin(min(1, sqrt(aa_x))))].';
NBLK_ROT = 200;

% The frame pass - the antenna-frame pedestal, the frame-pooled theta0(z)
% profile, and the per-SEGMENT geographic theta0(z) profiles the blocks
% inherit - lives in ptt.quadpolFrameTheta, so the synthetic gates
% (test_quadpol_segmented.m) exercise exactly the code that runs here.
% Segments re-fit theta0 laterally every ~SEG_LEN_M because the pooled
% frame pass is the two-pass design's single point of failure on
% laterally-varying frames; the pedestal stays frame-level because it is
% an instrument constant. Curved frames keep the validated geographic
% path (test_quadpol_curved.m) inside the helper.
fp = ptt.quadpolFrameTheta(T, z, az_tr, x_along, struct('fc', FC, ...
  'deramped', true, 'dlam_max', DLAM_MAX, 'seg_len_m', SEG_LEN_M, ...
  'nblk_rot', NBLK_ROT, 'track_az', track_az));
lsq = fp.lsq;
curved = fp.curved;
hspread = fp.hspread;
ped_ant = fp.ped_ant;
if ~curved
  th_geo_raw = lsq.theta0 + deg2rad(track_az);
else
  th_geo_raw = lsq.theta0;   % the geographic fit reports geographically
end
okt = isfinite(th_geo_raw);
% Blocks inherit the ANTENNA-frame pedestal (an instrument constant, so it
% is the same in every block's own frame) and their segment's geographic
% axis through ptt.thetaProfileAt below.
if all(isfinite(ped_ant))
  blk_ped = ped_ant;
else
  blk_ped = 'frame';
end
fprintf(['LS frame pass %.1f min (%s, heading p95 spread %.1f deg): ' ...
  'theta0 constrained on %d of %d windows, pedestal [%.3f %+.3fi %.3f], ' ...
  '%d segment(s) of ~%.1f km\n'], ...
  toc(t0)/60, H_tag(curved, 'GEOGRAPHIC frame', 'antenna frame'), hspread, ...
  nnz(okt), numel(th_geo_raw), ped_ant(1), ped_ant(2), ped_ant(3), ...
  fp.nseg, (x_along(end) - x_along(1)) / max(fp.nseg, 1) / 1000);

%% 4b. the SECTION: both estimators, per along-track block
% The frame-average profile above answers "what is the fabric here"; this
% answers "how does it vary along the line", which is the quantity the
% co-polarized section already shows and the one worth comparing against.
% The LS blocks inherit the frame theta0 profile, so their two remaining
% parameters are purely local and the section carries no block-to-block
% axis jitter.
t0 = tic;
nb = max(1, floor(Nx / NBLK_TR));
sec_dlam = nan(numel(z), nb);
sec_theta = nan(numel(z), nb);
sec_cmag = nan(numel(z), nb);
sec_dlam_ls = nan(numel(z), nb);
sec_resid_ls = nan(numel(z), nb);
sec_lat = nan(1, nb); sec_lon = nan(1, nb); sec_az = nan(1, nb);
for b = 1:nb
  j0 = (b-1)*NBLK_TR + 1;
  j1 = min(b*NBLK_TR, Nx);
  % Minimum traces per block, scaled: the fixed 16 predates small blocks
  % and would silently skip EVERY block of a 14-trace section, writing an
  % all-NaN product that looks like abstention. For NBLK_TR >= 17 this is
  % exactly the old guard; below that it keeps full blocks and drops only
  % sub-full tails.
  if j1 - j0 < min(16, NBLK_TR - 1), continue; end
  Tb = struct();
  for k = 1:4, Tb.(CHAN{k}) = T.(CHAN{k})(:, j0:j1); end
  ob = ptt.ershadiFabric(Tb, z, struct('fc', FC, ...
    'psi_step_deg', PSI_STEP_DEG, 'win_m', 30, 'grad_win_m', 25, ...
    'coh_min', 0.4, 'deramped', true));
  sec_dlam(:, b) = ob.dlam;
  sec_theta(:, b) = rad2deg(ob.theta);
  sec_cmag(:, b) = ob.Cmag(:, 1);
  hb_blk = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr(j0:j1))))))/2, 180);
  % the block's SEGMENT profile, phasor-interpolated across segment
  % centres at the block centre, converted into the block's antenna frame
  pg = ptt.thetaProfileAt(fp, mean(x_along(j0:j1)));
  if isstruct(pg)
    th_b = struct('z', pg.z, 'theta', pg.theta - deg2rad(hb_blk));
  else
    th_b = [];
  end
  ob_ls = ptt.quadpolFabricLS(Tb, z, struct('fc', FC, 'deramped', true, ...
    'theta0', th_b, 'pedestal', blk_ped, 'dlam_max', DLAM_MAX));
  sec_dlam_ls(:, b) = ob_ls.dlam_z;
  sec_resid_ls(:, b) = interp1(ob_ls.zw, ob_ls.resid, z, 'linear');
  sec_lat(b) = mean(la(j0:j1));
  sec_lon(b) = mean(lo(j0:j1));
  p0b = deg2rad(la(j0)); p1b = deg2rad(la(j1));
  dlb = deg2rad(lo(j1) - lo(j0));
  sec_az(b) = mod(rad2deg(atan2(sin(dlb)*cos(p1b), ...
    cos(p0b)*sin(p1b) - sin(p0b)*cos(p1b)*cos(dlb))), 180);
  clear Tb ob ob_ls;
end
fprintf('section: %d blocks of %d traces, %.1f min\n', nb, NBLK_TR, ...
  toc(t0)/60);
fprintf('section dlam %.3f..%.3f (median %.3f) | LS median %.3f\n', ...
  min(sec_dlam(:)), max(sec_dlam(:)), median(sec_dlam(:), 'omitnan'), ...
  median(sec_dlam_ls(:), 'omitnan'));

zb = band;
fprintf('\n%-14s %10s %10s\n', '', 'coreg', 'raw');
fprintf('%-14s %10.1f %10.1f\n', 'theta_ant', ...
  H_cmed(rad2deg(out.theta(zb))), H_cmed(rad2deg(out_raw.theta(zb))));
fprintf('%-14s %10.1f %10.1f\n', 'theta_geo', ...
  mod(H_cmed(rad2deg(out.theta(zb))) + track_az, 180), ...
  mod(H_cmed(rad2deg(out_raw.theta(zb))) + track_az, 180));
fprintf('%-14s %10.3f %10.3f\n', 'dlam', ...
  median(out.dlam(zb), 'omitnan'), median(out_raw.dlam(zb), 'omitnan'));
fprintf('%-14s %10.3f %10.3f\n', '|C_hhvv|', ...
  median(out.Cmag(zb,1), 'omitnan'), median(out_raw.Cmag(zb,1), 'omitnan'));
fprintf('%-14s %10.1f\n', 'track_az', track_az);
zw_band = lsq.zw > Z_BAND(1) & lsq.zw < Z_BAND(2);
fprintf('\nLS fit (frame): theta0_geo %.1f  dlam %.3f  ' ...
  , H_cmed(rad2deg(th_geo_raw(zw_band))), ...
  median(lsq.dlam(zw_band), 'omitnan'));
fprintf('resid %.3f  q_theta %.2f\n', ...
  median(lsq.resid(zw_band), 'omitnan'), ...
  median(lsq.q_theta(zw_band), 'omitnan'));
% The reference for the ershadi rows is the LS frame fit printed just
% above, which is the comparable single-line estimate. The two-azimuth
% solve is context only: the LS estimator superseded it (its correlated
% (theta, P) errors bias it high) and it is separately near-singular at
% Ridge A's geometry, so nothing here is measured against it.
fprintf(['(reference: the LS frame fit above. The ershadi rows are the ' ...
  'antenna-locked\n projection and sit LOW by ~cos 2*offset against ' ...
  'it.)\n']);
fprintf(['(two-azimuth solve for Ridge A, orientation only: principal ' ...
  'contrast ~0.12 at\n the plateau, axis ~110-120 deg geo. Not a ' ...
  'standard to measure against - it\n is a different quantity from the ' ...
  'LS dlam above, and its conditioning goes\n as |sin 2*da| in the ' ...
  'azimuth separation, which vanishes at BOTH 0 and 90\n deg; Ridge ' ...
  'A''s 81.4 deg leg separation sits 9 deg from the 90 deg\n ' ...
  'singularity, amplifying any per-leg bias ~3.4x.)\n']);

%% 5. save, small
ZS = 4;
s = 1:ZS:numel(z);
res = struct('tag', sprintf('%s_%03d', day_seg, frm), ...
  'day_seg', day_seg, 'frm', frm, 'track_az', track_az, ...
  'z', single(z(s)), ...
  'theta', single(rad2deg(out.theta(s))), ...
  'dlam', single(out.dlam(s)), ...
  'Cmag', single(out.Cmag(s,1)), ...
  'theta_raw', single(rad2deg(out_raw.theta(s))), ...
  'dlam_raw', single(out_raw.dlam(s)), ...
  'ls_zw', single(lsq.zw), 'ls_theta0', single(rad2deg(lsq.theta0)), ...
  'ls_theta0_geo', single(rad2deg(th_geo_raw)), ...
  'ls_curved', curved, 'ls_hspread', hspread, ...
  'ls_dlam', single(lsq.dlam), 'ls_q_theta', single(lsq.q_theta), ...
  'ls_resid', single(lsq.resid), 'ls_gamma', single(lsq.gamma), ...
  'ls_leak', single(lsq.leak), 'ls_pedestal', ped_ant, ...
  'ls_theta_seg', single(rad2deg(fp.th_seg)), ...
  'ls_q_seg', single(fp.q_seg), ...
  'ls_dlam_seg', single(fp.dlam_seg), ...
  'ls_resid_seg', single(fp.resid_seg), ...
  'ls_seg_x', fp.seg_x, 'ls_nseg', fp.nseg, 'seg_len_m', SEG_LEN_M, ...
  'theta0_ls', single(rad2deg(lsq.theta0_z(s))), ...
  'dlam_ls', single(lsq.dlam_z(s)), ...
  'sec_dlam_ls', single(sec_dlam_ls(s,:)), ...
  'sec_resid_ls', single(sec_resid_ls(s,:)), ...
  'sec_dlam', single(sec_dlam(s,:)), 'sec_theta', single(sec_theta(s,:)), ...
  'sec_cmag', single(sec_cmag(s,:)), 'sec_lat', sec_lat, ...
  'sec_lon', sec_lon, 'sec_az', sec_az, 'nblk_tr', NBLK_TR, ...
  'row_med', cinfo.row_med, 'coh_before', coh_before, ...
  'coh_after', coh_after, 'pairs', {pairs}, ...
  'lat', median(la), 'lon', median(lo), ...
  'lat0', la(1), 'lon0', lo(1), 'lat1', la(end), 'lon1', lo(end), ...
  't_coreg_min', t_coreg);
% `quadpol_section_` and not `quadpol_`: run_quadpol_frame.m already writes
% quadpol_<day_seg>_<frm>.mat with a flat out/z/psi/A layout, and both
% mirror to the same figure staging directory. Sharing the name would let
% whichever ran last break the other figure's reader, since this file is a
% single `res` struct instead.
out_fn = fullfile(out_dir, ...
  sprintf('quadpol_section_%s_%03d%s.mat', day_seg, frm, ZTAG));
save(out_fn, '-v7.3', 'res');
fprintf('\nwrote %s (total %.1f min)\n', out_fn, toc(t_all)/60);

function c = H_coh(a, b, kr, band)
num = conv2(a .* conj(b), kr, 'same');
d1 = conv2(abs(a).^2, kr, 'same');
d2 = conv2(abs(b).^2, kr, 'same');
cc = abs(num) ./ sqrt(max(d1.*d2, realmin));
c = median(cc(band,:), 'all', 'omitnan');
end

function m = H_cmed(a)
a = a(isfinite(a));
if isempty(a), m = NaN; return; end
m = mod(rad2deg(angle(mean(exp(2i*deg2rad(a)))))/2, 180);
end

function s = H_tag(ok, yes, no)
if ok, s = yes; else, s = no; end
end

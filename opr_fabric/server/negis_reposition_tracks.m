%NEGIS_REPOSITION_TRACKS Recover NEGIS frame trajectories from the GPS files.
%
% The 2024_Greenland_Ground2 records-level radar-time-to-GPS-time sync
% failed for every segment except 20240618_01, so the qlook products carry
% placeholder Latitude/Longitude (~0.2 N, ~1.8 E). The GNSS itself is
% fine: gps_<date>.mat holds sync_gps_time/sync_lat/sync_lon covering the
% survey. Re-interpolating position onto the product's GPS_time recovers
% the trajectory without touching the radar data - the same repair Knut
% Christianson applied by hand to 20240618_01.
%
% Writes negis2024_tracks_fixed.mat (per-frame decimated lat/lon plus the
% GPS coverage fraction, so a frame whose GPS_time falls outside its day
% file is visible rather than silently extrapolated).
season = '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2';
gps_dir = '/cresis/dataproducts/opr_data/opr_support/gps/2024_Greenland_Ground2';
out_fn = '/kucresis/scratch/hoffmana_sta/fabric/stages/negis2024_tracks_fixed.mat';

segs = dir(fullfile(season, 'CSARP_qlook_HH'));
segs = segs([segs.isdir] & ~startsWith({segs.name}, '.'));

T = struct('seg', {}, 'frm', {}, 'lat', {}, 'lon', {}, 'nx', {}, ...
  'cover', {}, 'len_km', {}, 'e2e_km', {}, 'cpx', {});
gps_cache = struct();

fprintf('%-14s %-30s %6s %6s %7s %7s %5s\n', 'segment', 'frame', 'ntr', ...
  'cover', 'len_km', 'e2e_km', 'cpx');
for i = 1:numel(segs)
  day = segs(i).name(1:8);
  gkey = ['g' day];
  if ~isfield(gps_cache, gkey)
    gfn = fullfile(gps_dir, sprintf('gps_%s.mat', day));
    if exist(gfn, 'file') ~= 2
      warning('No GPS file for %s', day);
      continue;
    end
    g = load(gfn, 'sync_gps_time', 'sync_lat', 'sync_lon');
    % Duplicate timestamps break interp1; keep the first of each
    [tu, iu] = unique(g.sync_gps_time(:));
    gps_cache.(gkey) = struct('t', tu, 'lat', g.sync_lat(iu), ...
      'lon', g.sync_lon(iu));
  end
  G = gps_cache.(gkey);

  fs = dir(fullfile(season, 'CSARP_qlook_HH', segs(i).name, 'Data_2*.mat'));
  for j = 1:numel(fs)
    fn = fullfile(fs(j).folder, fs(j).name);
    d = load(fn, 'GPS_time');
    w = whos('-file', fn, 'Data');
    gt = d.GPS_time(:);
    cover = mean(gt >= G.t(1) & gt <= G.t(end));
    la = interp1(G.t, G.lat, gt, 'linear', NaN);
    lo = interp1(G.t, G.lon, gt, 'linear', NaN);
    ok = isfinite(la) & isfinite(lo);
    len_km = NaN; e2e_km = NaN;
    if nnz(ok) > 2
      lag = la(ok); log_ = lo(ok);
      R = 6371.0; p = deg2rad(lag); dl = diff(deg2rad(log_)); dp = diff(p);
      a = sin(dp/2).^2 + cos(p(1:end-1)).*cos(p(2:end)).*sin(dl/2).^2;
      len_km = sum(2*R*asin(sqrt(a)));
      p0 = deg2rad(lag(1)); p1 = deg2rad(lag(end));
      dl2 = deg2rad(log_(end) - log_(1));
      a2 = sin((p1-p0)/2).^2 + cos(p0)*cos(p1)*sin(dl2/2)^2;
      e2e_km = 2*R*asin(sqrt(a2));
    end
    k = numel(T) + 1;
    T(k).seg = segs(i).name;
    T(k).frm = fs(j).name;
    T(k).lat = single(la(1:10:end));
    T(k).lon = single(lo(1:10:end));
    T(k).nx = numel(gt);
    T(k).cover = cover;
    T(k).len_km = len_km;
    T(k).e2e_km = e2e_km;
    T(k).cpx = w.complex;
    fprintf('%-14s %-30s %6d %6.2f %7.2f %7.2f %5d\n', segs(i).name, ...
      fs(j).name, numel(gt), cover, len_km, e2e_km, w.complex);
  end
end

save(out_fn, '-v7', 'T');
fprintf('\nwrote %s (%d frames)\n', out_fn, numel(T));

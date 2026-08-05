%NEGIS_INTERFEROGRAM Form HH-VV polarimetric interferograms for NEGIS.
%
% The 2024_Greenland_Ground2 qlook products were written with inc_dec = 0
% for every segment except 20240618_01, so Data is COMPLEX and the
% HH/VV phase difference survives - the same quantity CSARP_polarimetric
% carries for the Antarctic seasons, just without SAR focusing (this is a
% ground-based radar with coincident phase centres, so the qlook product
% is already at nadir).
%
% Positions come from the GPS file, not from the product: the records
% time sync failed for these segments (see negis_reposition_tracks.m).
% CSARP_layer would be the other source, but it exists for 20240618_01
% only - the one segment that is already georeferenced AND was written
% with incoherent decimation, so it has no phase to interferogram. On
% that segment the two sources agree to ~1 m (median 1.1 m north,
% 1.1 m east; see negis_layer_check.m), so interpolating the GPS file
% onto GPS_time is the layer positions by another route, and it is the
% only route that covers the complex segments.
%
% Traces the traverse stopped for are culled BEFORE the multilook. They
% multilook to HIGH coherence with random range phase, so no coherence
% test downstream removes them, and an azimuth cell that straddles a stop
% boundary mixes stopped and moving traces permanently - culling after the
% looks (as negis_two_panel.py used to be the only place to do) cannot
% undo that. This is the same criterion run_negis_fabric.m applies on the
% inversion path.
%
% RE-RUN REQUIRED: the pre-cull extracts already on mem1 were written
% without it and their boundary columns are contaminated. Delete
% stages/negis_ifg_*.mat and re-run this script before regenerating the
% NEGIS figures.
%
% For each requested frame this writes a multilooked interferogram,
% coherence, repositioned trajectory and a power-based surface pick,
% small enough to copy back for figures.
season = '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2';
gps_dir = '/cresis/dataproducts/opr_data/opr_support/gps/2024_Greenland_Ground2';
out_dir = '/kucresis/scratch/hoffmana_sta/fabric/stages';

% seg / frame candidates, ranked by along-track speed contrast over
% ITS_LIVE with straightness >= 0.93
frames = { '20240620_01', 3; '20240619_01', 4; '20240620_02', 2; ...
           '20240621_01', 4; '20240626_03', 1 };

NR = 4;    % range looks
NA = 4;    % azimuth looks
MIN_STEP_M = 1.0;   % below this along-track step the vehicle was stopped

for fi = 1:size(frames, 1)
  day_seg = frames{fi, 1};
  frm = frames{fi, 2};
  name = sprintf('Data_%s_%03d.mat', day_seg, frm);
  hh_fn = fullfile(season, 'CSARP_qlook_HH', day_seg, name);
  vv_fn = fullfile(season, 'CSARP_qlook_VV', day_seg, name);
  fprintf('\n=== %s frame %d ===\n', day_seg, frm);
  if exist(hh_fn,'file') ~= 2 || exist(vv_fn,'file') ~= 2
    warning('missing HH or VV for %s', name); continue;
  end

  H = load(hh_fn, 'Data', 'Time', 'GPS_time', 'Elevation');
  V = load(vv_fn, 'Data', 'Time', 'GPS_time');
  assert(isequal(H.Time, V.Time), 'HH/VV range axes differ');
  assert(isequal(H.GPS_time, V.GPS_time), 'HH/VV GPS times differ');
  assert(~isreal(H.Data) && ~isreal(V.Data), ...
    'Data is real - this frame was written with incoherent decimation');

  [Nt, Nx] = size(H.Data);
  fprintf('  %d samples x %d traces, Time %.2f..%.2f us\n', Nt, Nx, ...
    H.Time(1)*1e6, H.Time(end)*1e6);

  % --- reposition from GPS
  g = load(fullfile(gps_dir, sprintf('gps_%s.mat', day_seg(1:8))), ...
    'sync_gps_time', 'sync_lat', 'sync_lon', 'sync_elev');
  [tu, iu] = unique(g.sync_gps_time(:));
  gt = H.GPS_time(:);
  Latitude = interp1(tu, g.sync_lat(iu), gt, 'linear', NaN).';
  Longitude = interp1(tu, g.sync_lon(iu), gt, 'linear', NaN).';
  fprintf('  repositioned: lat %.4f..%.4f lon %.4f..%.4f\n', ...
    min(Latitude), max(Latitude), min(Longitude), max(Longitude));

  % --- cull the stops, before any looks are taken. Data, GPS_time,
  % Elevation and the repositioned coordinates go together so every axis
  % stays aligned with the traces that survive.
  Rearth = 6371e3;
  p = deg2rad(Latitude); dl = diff(deg2rad(Longitude)); dp = diff(p);
  a = sin(dp/2).^2 + cos(p(1:end-1)).*cos(p(2:end)).*sin(dl/2).^2;
  step = 2*Rearth*asin(sqrt(a));
  keep = [true, step >= MIN_STEP_M];
  fprintf('  culling %d of %d traces where the vehicle was stopped\n', ...
    nnz(~keep), Nx);
  H.Data = H.Data(:, keep); V.Data = V.Data(:, keep);
  Latitude = Latitude(keep); Longitude = Longitude(keep);
  H.GPS_time = H.GPS_time(keep); H.Elevation = H.Elevation(keep);
  Nx = nnz(keep);

  % --- surface pick from HH power (Surface is NaN in these products)
  pw = abs(H.Data).^2;
  [~, si] = max(pw, [], 1);
  Surface = H.Time(si).';
  fprintf('  surface pick %.3f..%.3f us (median %.3f)\n', ...
    min(Surface)*1e6, max(Surface)*1e6, median(Surface)*1e6);

  % --- multilook: coherent sum of HH*conj(VV), power sums for coherence
  ntc = floor(Nt/NR); nxc = floor(Nx/NA);
  cut = @(A) reshape(A(1:ntc*NR, 1:nxc*NA), NR, ntc, NA, nxc);
  I = H.Data(1:ntc*NR, 1:nxc*NA) .* conj(V.Data(1:ntc*NR, 1:nxc*NA));
  interferogram_mlook = squeeze(sum(sum(cut(I), 1), 3));
  power_hh = squeeze(sum(sum(cut(abs(H.Data).^2), 1), 3));
  power_vv = squeeze(sum(sum(cut(abs(V.Data).^2), 1), 3));
  interferogram_coherence = abs(interferogram_mlook) ...
    ./ sqrt(power_hh .* power_vv);
  % The power sums travel with the interferogram so the look count can be
  % raised later without coming back to the server: coherence at any
  % coarser cell is |sum I| / sqrt(sum|HH|^2 * sum|VV|^2) over the same
  % cell, which needs the un-normalised sums, not this ratio.
  interferogram_mlook = single(interferogram_mlook);
  interferogram_coherence = single(interferogram_coherence);
  power_hh = single(power_hh);
  power_vv = single(power_vv);
  clear I pw;

  % Label each cell with its CENTRE: interferogram_mlook is the sum over
  % the whole cell, so taking the first sample/trace would put both axes
  % half a cell shallower/earlier than the data they annotate.
  Time = mean(reshape(H.Time(1:ntc*NR), NR, ntc), 1).';
  Latitude_ml = mean(reshape(Latitude(1:nxc*NA), NA, nxc), 1);
  Longitude_ml = mean(reshape(Longitude(1:nxc*NA), NA, nxc), 1);
  Surface_ml = mean(reshape(Surface(1:nxc*NA), NA, nxc), 1);
  Elevation_ml = mean(reshape(H.Elevation(1:nxc*NA), NA, nxc), 1);
  GPS_time_ml = mean(reshape(H.GPS_time(1:nxc*NA), NA, nxc), 1);
  looks = [NR NA];

  fprintf('  multilooked to %d x %d, median coherence %.3f\n', ...
    size(interferogram_mlook,1), size(interferogram_mlook,2), ...
    median(interferogram_coherence(isfinite(interferogram_coherence))));

  out_fn = fullfile(out_dir, sprintf('negis_ifg_%s_%03d.mat', day_seg, frm));
  save(out_fn, '-v7.3', 'interferogram_mlook', 'interferogram_coherence', ...
    'power_hh', 'power_vv', 'Time', 'Latitude_ml', 'Longitude_ml', ...
    'Surface_ml', 'Elevation_ml', 'GPS_time_ml', 'looks', 'day_seg', 'frm');
  d = dir(out_fn);
  fprintf('  wrote %s (%.0f MB)\n', out_fn, d.bytes/1e6);
  clear H V interferogram_mlook interferogram_coherence;
end
fprintf('\nDone.\n');

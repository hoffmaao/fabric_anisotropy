%EGRIP_RECORDS_RESYNC Re-sync the 2024_Greenland_Ground2 records to the final GPS.
%
%   matlab -batch "egrip_records_resync"        % on mem1, from this directory
%
% WHAT AND WHY. Every EastGRIP day now has a post-processed GPS solution
% (gps_source cresis-final_2025, millions of samples) in the group OPR
% support tree, but the RECORDS of every segment from 20240620_01 onward
% are still synced to the Arena field GPS (gps_source arena-field): the
% records were built before the final solutions landed and the qlook
% channel products were built on those records. That is the whole reason
% only 5 of 55 EastGRIP frames are positioned well enough for the fabric
% pipeline's trajectory rebuild. Re-syncing a records file to the current
% GPS is exactly what the toolbox's records_update does; nothing else in
% the records changes.
%
% WHERE. Never in the group tree. The records (and frames) were copied to
% a PRIVATE support tree, ~/scratch/opr_support_egrip, whose gps/ is a
% symlink to the group gps directory (read-only use), and gRadar's
% support_path is overridden to that tree for the duration of the call.
% The production records stay as they are until someone with authority
% over them decides to replace them; the reprocessing that follows
% (qlook or sar+array per channel) reads from the private tree by the
% same override.
%
% Each segment's parameters come from the param_records struct stored in
% its own records file - the settings it was actually built with - since
% the season sheets still on disk no longer list the later segments.
%
% Output: the private records rewritten in place, the reference
% trajectories under <private>/opr_data/accum/<season>/CSARP_reference_
% trajectory, and a summary table
% records_resync_summary.txt beside them: per segment the old and new
% gps_source, the median and 95th-percentile horizontal shift of the
% trajectory (m), and the fraction of records left without a position.
global gRadar;
if isempty(gRadar), startup; end
priv = fullfile(getenv('HOME'), 'scratch', 'opr_support_egrip');
season = '2024_Greenland_Ground2';
rdir = fullfile(priv, 'records', 'accum', season);
assert(exist(rdir, 'dir') == 7, 'private records dir missing: %s', rdir);
assert(exist(fullfile(priv, 'gps', season), 'dir') == 7, 'gps link missing under %s', priv);
if ~exist('only', 'var'), only = ''; end   % optional day_seg filter
d = dir(fullfile(rdir, 'records_2024*.mat'));
fid = fopen(fullfile(rdir, 'records_resync_summary.txt'), 'w');
fprintf(fid, '%-12s %-16s %-18s %8s %8s %8s %6s\n', 'day_seg', 'old_source', 'new_source', 'med_m', 'p95_m', 'max_m', 'nan');
for k = 1:numel(d)
  fn = fullfile(d(k).folder, d(k).name);
  old = load(fn, 'gps_source', 'lat', 'lon', 'gps_time', 'param_records');
  day_seg = old.param_records.day_seg;
  if ~isempty(only) && ~strcmp(day_seg, only), continue; end
  param = old.param_records;
  param.records.gps.en = 1;
  param.records.gps.fn = '';            % default: gps_<day>.mat of the season
  ov = gRadar;
  ov.support_path = priv;               % records in, records out
  % records_update also regenerates the segment's REFERENCE TRAJECTORY
  % through opr_filename_out, i.e. under gRadar.out_path - the production
  % data tree. The first run of this script (2 Sep 2026, 20240620_02) did
  % exactly that before this line existed and rewrote
  % CSARP_reference_trajectory/ref_20240620_02.mat there. Everything this
  % script produces now lands under the private data tree; the
  % reprocessing that follows reads it through the same override.
  ov.out_path = fullfile(priv, 'opr_data');
  fprintf('\n=== %s (was %s) ===\n', day_seg, old.gps_source);
  try
    records_update(param, ov);
  catch e
    fprintf(fid, '%-12s %-16s %-18s FAILED: %s\n', day_seg, old.gps_source, '-', e.message);
    fprintf('FAILED %s: %s\n', day_seg, e.message);
    continue;
  end
  new = load(fn, 'gps_source', 'lat', 'lon');
  R_E = 6371000;
  dph = deg2rad(new.lat - old.lat); dlo = deg2rad(new.lon - old.lon);
  aa = sin(dph/2).^2 + cos(deg2rad(old.lat)) .* cos(deg2rad(new.lat)) .* sin(dlo/2).^2;
  sh = 2 * R_E * asin(min(1, sqrt(aa)));
  ok = isfinite(sh);
  fprintf(fid, '%-12s %-16s %-18s %8.2f %8.2f %8.2f %6.3f\n', day_seg, old.gps_source, new.gps_source, ...
    median(sh(ok)), prctile(sh(ok), 95), max(sh(ok)), mean(~isfinite(new.lat)));
  fprintf('%s: %s -> %s; shift median %.2f m, p95 %.2f m, max %.2f m; unpositioned %.3f\n', day_seg, ...
    old.gps_source, new.gps_source, median(sh(ok)), prctile(sh(ok), 95), max(sh(ok)), mean(~isfinite(new.lat)));
end
fclose(fid);
fprintf('\nsummary: %s\n', fullfile(rdir, 'records_resync_summary.txt'));

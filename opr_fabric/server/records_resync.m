%RECORDS_RESYNC Re-sync a season's records to a chosen GPS, in a private tree.
%
%   matlab -batch "season='2022_Antarctica_Ground'; only={'20230120_04','20230120_05'}; records_resync"
%   matlab -batch "season='...'; dry=true; records_resync"
%
% The general form of egrip_records_resync.m, which did this for EastGRIP
% alone. It re-runs the toolbox's records_update for the chosen segments so
% they pick up whatever GPS now sits in the private tree - a re-created
% solution, or one repaired by gnss_despike.m - and writes every product
% into that private tree rather than over the group's.
%
% WHY A PRIVATE TREE, AND WHAT BURNED ME ONCE. records_update writes the
% records back through gRadar.support_path AND regenerates the segment's
% reference trajectory through gRadar.out_path, which is the DATA tree.
% Overriding only the first still rewrites a production file: doing exactly
% that on 2 Sep 2026 replaced
% CSARP_reference_trajectory/ref_20240620_02.mat in the live
% 2024_Greenland_Ground2 season. Both paths are overridden here, and the
% script refuses to run if either still points inside the group tree.
%
% Each segment's parameters come from the param_records struct stored in
% its own records file, because the season sheets on disk no longer list
% every segment that was collected.
%
% SETUP the private tree once per season (records and frames copied, gps
% either symlinked to the group's or holding repaired files):
%   priv=~/scratch/opr_support_<season>
%   <priv>/records/accum/<season>/records_*.mat
%   <priv>/frames/accum/<season>/frames_*.mat
%   <priv>/gps/<season>/gps_*.mat
% ensure_gps below fills any gps day the private tree lacks by symlinking
% the group's, so a partially repaired season still resolves every segment.
%
% Inputs (workspace variables)
%   season  OPR season name
%   priv    private support root (default ~/scratch/opr_support_<season>)
%   only    cellstr of day_seg to do (default {} = all in the private tree)
%   dry     true to report what would run and stop (default false)
%
% Output: the private records rewritten, reference trajectories under
% <priv>/opr_data, and records_resync_summary.txt beside the records with
% the old and new gps_source and the horizontal shift each segment moved.
global gRadar;
if isempty(gRadar), startup; end
if ~exist('season', 'var'), error('set season'); end
season = char(season);
if ~exist('priv', 'var') || isempty(priv)
  priv = fullfile(getenv('HOME'), 'scratch', ['opr_support_' season]);
end
if ~exist('only', 'var'), only = {}; end
if ischar(only) || isstring(only), only = cellstr(only); end
if ~exist('dry', 'var') || isempty(dry), dry = false; end
rdir = fullfile(priv, 'records', 'accum', season);
assert(exist(rdir, 'dir') == 7, 'private records dir missing: %s', rdir);
% refuse to touch the group tree
grp = gRadar.support_path;
assert(~strncmp(priv, grp, numel(grp)), ...
  'priv (%s) is inside the group support tree (%s)', priv, grp);

% --- make sure every gps day the season needs resolves in the private tree
gdir_priv = fullfile(priv, 'gps', season);
gdir_grp = fullfile(grp, 'gps', season);
if exist(gdir_priv, 'dir') ~= 7, mkdir(gdir_priv); end
gg = dir(fullfile(gdir_grp, 'gps_*.mat'));
for k = 1:numel(gg)
  tgt = fullfile(gdir_priv, gg(k).name);
  if exist(tgt, 'file') ~= 2
    % symlink, so a repaired file already present is never overwritten and
    % the group's original is never copied twice
    system(sprintf('ln -sfn %s %s', fullfile(gdir_grp, gg(k).name), tgt));
  end
end
dspk = dir(fullfile(gdir_priv, 'gps_*.mat'));
n_real = 0;
for k = 1:numel(dspk)
  if ~H_islink(fullfile(gdir_priv, dspk(k).name)), n_real = n_real + 1; end
end
fprintf('gps in private tree: %d day(s), %d repaired in place, %d linked to the group\n', ...
  numel(dspk), n_real, numel(dspk) - n_real);

d = dir(fullfile(rdir, 'records_*.mat'));
fid = -1;
if ~dry
  fid = fopen(fullfile(rdir, 'records_resync_summary.txt'), 'a');
  fprintf(fid, '# %s  season %s\n', datestr(now, 'yyyy-mm-dd HH:MM'), season);
  fprintf(fid, '%-14s %-20s %-22s %9s %9s %9s %7s\n', 'day_seg', 'old_source', 'new_source', 'med_m', 'p95_m', 'max_m', 'nan');
end
R = 6371000;
for k = 1:numel(d)
  fn = fullfile(d(k).folder, d(k).name);
  old = load(fn, 'gps_source', 'lat', 'lon', 'param_records');
  day_seg = old.param_records.day_seg;
  if ~isempty(only) && ~any(strcmp(day_seg, only)), continue; end
  if dry
    fprintf('would resync %s (currently %s)\n', day_seg, old.gps_source);
    continue;
  end
  param = old.param_records;
  param.records.gps.en = 1;
  % PRESERVE the stored gps.fn. It is empty at EastGRIP, where the default
  % gps_<day>.mat rule resolves, but the Kamb segments carry an explicit
  % 'gps/<season>/gps_20230121.mat' because that season's day files are
  % offset from the segment dates - gps_20230121.mat spans 20 Jan 08:20 to
  % 21 Jan 03:42 UTC and there is no gps_20230120.mat at all. Blanking it
  % would send records_update looking for a file that does not exist, or
  % worse, silently onto a different day. It is a RELATIVE path, so it
  % resolves under the private support_path and picks up a repaired file
  % there automatically.
  if ~isfield(param.records.gps, 'fn'), param.records.gps.fn = ''; end
  ov = gRadar;
  ov.support_path = priv;
  ov.out_path = fullfile(priv, 'opr_data');
  fprintf('\n=== %s (was %s) ===\n', day_seg, old.gps_source);
  try
    records_update(param, ov);
  catch e
    fprintf(fid, '%-14s %-20s FAILED: %s\n', day_seg, old.gps_source, e.message);
    fprintf('FAILED %s: %s\n', day_seg, e.message);
    continue;
  end
  new = load(fn, 'gps_source', 'lat', 'lon');
  dph = deg2rad(new.lat - old.lat); dlo = deg2rad(new.lon - old.lon);
  aa = sin(dph/2).^2 + cos(deg2rad(old.lat)) .* cos(deg2rad(new.lat)) .* sin(dlo/2).^2;
  sh = 2 * R * asin(min(1, sqrt(aa)));
  ok = isfinite(sh);
  fprintf(fid, '%-14s %-20s %-22s %9.3f %9.3f %9.3f %7.3f\n', day_seg, old.gps_source, ...
    new.gps_source, median(sh(ok)), prctile(sh(ok), 95), max(sh(ok)), mean(~isfinite(new.lat)));
  fprintf('%s: %s -> %s; moved median %.3f m, p95 %.3f m, max %.3f m\n', day_seg, ...
    old.gps_source, new.gps_source, median(sh(ok)), prctile(sh(ok), 95), max(sh(ok)));
end
if fid > 0
  fclose(fid);
  fprintf('\nsummary: %s\n', fullfile(rdir, 'records_resync_summary.txt'));
end

function tf = H_islink(p)
[st, out] = system(sprintf('test -L %s && echo yes || echo no', p));
tf = st == 0 && strncmp(strtrim(out), 'yes', 3);
end

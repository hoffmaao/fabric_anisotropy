%GNSS_DESPIKE Repair impossible jumps in a ground-traverse GNSS solution.
%
%   matlab -batch "season='2022_Antarctica_Ground'; day='20230121'; gnss_despike"
%   matlab -batch "season='...'; day='...'; dev_max=3; dry=true; gnss_despike"
%
% WHY THIS EXISTS, AND WHY IT IS NOT THE SAME FIX AS THE EASTGRIP RE-SYNC.
% Two different faults produce a bad trajectory, and only one of them has a
% better solution sitting unused:
%
%   EASTGRIP (2024_Greenland_Ground2). A post-processed solution
%   (cresis-final_2025) exists for every day, but the records of
%   20240620_01 onward were still synced to the Arena field GPS, which had
%   labelled them near 0 N 2 E. The fix is to RE-SYNC - see
%   egrip_records_resync.m - and it recovers the whole trajectory.
%
%   KAMB (2022_Antarctica_Ground, segments 20230120_*). There is NO
%   post-processed solution for January 2023: the season's ie_* files stop
%   on 12 December 2022 and gps_create_2022_Antarctica_Ground.m hard-codes
%   gps_source_to_use = 'arena'. The positions are real - the traverse runs
%   from WAIS Divide southwest toward Kamb - but the Arena solution carries
%   isolated spikes: segment 20230120_05 reaches 43 m/s and 20230120_04
%   reaches 89 m/s between consecutive records, against a 5.9 m/s median.
%   Nothing better exists to re-sync to, so the trajectory has to be
%   repaired in place.
%
% WHAT A SPIKE COSTS. SAR focusing integrates phase along the synthetic
% aperture from the trajectory, so a position error enters the matched
% filter directly: at 750 MHz in ice one wavelength is ~0.22 m, and an
% aperture that crosses a metre-scale jump defocuses. A qlook product,
% which sums power incoherently, barely notices the same spike. That is
% why this matters for the SAR path and did not for the products we have
% been inverting.
%
% HOW A SPIKE IS DETECTED, AND WHY NOT BY SPEED. The obvious test - the
% speed implied by the step from the previous fix - is wrong at high
% sample rates, and silently so. EastGRIP's solution is 20 Hz, so 3 m of
% position jitter between consecutive fixes reads as 63 m/s; Kamb's is
% 1 Hz, where the same jitter reads as 3 m/s. A single speed ceiling
% therefore flags ordinary noise at one site and nothing at all at the
% other, which is exactly what a first version of this script did: it
% repaired Kamb (67.7 -> 14.7 m/s) and left EastGRIP's maximum untouched
% at 62.9 m/s while claiming 14229 spikes.
%
% So a sample is a SPIKE when it departs from its own LOCAL TREND by more
% than dev_max metres, the trend being a running median over trend_s
% seconds (a median, so the spikes themselves cannot drag it). That test
% is independent of sample rate and of how fast the platform is really
% moving: it asks whether the fix left the path and came back, which is
% what a bad fix does and what real motion does not. Spikes are replaced
% by linear interpolation in time from the nearest good fix on each side;
% runs longer than max_run samples are left alone and reported, because a
% long excursion is more likely real (a transit, or an aircraft leg inside
% a file that also holds ground survey) than a glitch, and should be
% looked at rather than smoothed away.
%
% Elevation gets the same treatment against elev_dev (default 5 m).
% Heading, roll and pitch are left untouched: they are separate
% measurements and a position spike does not imply an attitude spike.
%
% WHERE IT WRITES. Never over the group tree. The cleaned file goes to
% <priv>/gps/<season>/gps_<day>.mat with gps_source suffixed '+despike',
% so any records synced from it say so, and a summary line is appended to
% <priv>/gps/<season>/despike_summary.txt. Set dry=true to report without
% writing.
%
% Inputs (workspace variables)
%   season   OPR season name, e.g. '2022_Antarctica_Ground'
%   day      'YYYYMMDD' of the gps file to clean
%   priv     private support root (default ~/scratch/opr_support_<season>)
%   dev_max  horizontal departure from the local trend, m (default 3)
%   elev_dev vertical departure from the local trend, m (default 5)
%   trend_s  running-median length for the trend, s (default 2)
%   max_run  longest run of consecutive bad samples to repair (default 50)
%   dry      true to report only (default false)
global gRadar;
if isempty(gRadar), startup; end
if ~exist('season', 'var'), error('set season'); end
if ~exist('day', 'var'), error('set day'); end
season = char(season); day = char(day);
if ~exist('priv', 'var') || isempty(priv)
  priv = fullfile(getenv('HOME'), 'scratch', ['opr_support_' season]);
end
if ~exist('dev_max', 'var') || isempty(dev_max), dev_max = 3; end
if ~exist('elev_dev', 'var') || isempty(elev_dev), elev_dev = 5; end
if ~exist('trend_s', 'var') || isempty(trend_s), trend_s = 2; end
if ~exist('max_run', 'var') || isempty(max_run), max_run = 50; end
if ~exist('dry', 'var') || isempty(dry), dry = false; end

src = fullfile(gRadar.support_path, 'gps', season, ['gps_' day '.mat']);
if exist(src, 'file') ~= 2
  error('gnss_despike:missing', 'no gps file %s', src);
end
g = load(src);
t = g.gps_time(:);
[t, is] = sort(t);
fn_vec = {'lat', 'lon', 'elev', 'roll', 'pitch', 'heading'};
for k = 1:numel(fn_vec)
  if isfield(g, fn_vec{k}) && numel(g.(fn_vec{k})) == numel(is)
    g.(fn_vec{k}) = reshape(g.(fn_vec{k})(is), size(t));
  end
end
g.gps_time = t;
la = g.lat(:); lo = g.lon(:); el = g.elev(:);
n = numel(t);

R = 6371000;
step = @(a1, o1, a2, o2) 2*R*asin(min(1, sqrt(sin(deg2rad(a2-a1)/2).^2 ...
  + cos(deg2rad(a1)).*cos(deg2rad(a2)).*sin(deg2rad(o2-o1)/2).^2)));
dt = [Inf; diff(t)];
d = [0; step(la(1:end-1), lo(1:end-1), la(2:end), lo(2:end))];
v = d ./ max(dt, eps);
dt_typ = median(diff(t));
% A GAP is not a spike: across a long dead time the platform really has
% moved, and the trend either side is a different piece of track. Only
% samples surrounded by normal cadence are judged.
judgeable = dt <= 5 * dt_typ;
judgeable(1) = false;
% departure from the local trend, in metres. The window is odd and at
% least 3 samples; a median over trend_s seconds is short enough to follow
% a turn and long enough that a few bad fixes cannot move it.
w = max(3, 2*floor(trend_s / max(dt_typ, eps) / 2) + 1);
% movmedian, not medfilt1: it shrinks the window at the ends instead of
% zero-padding, so the first and last samples are compared against real
% neighbours rather than against zeros
la_t = movmedian(la, w, 'omitnan');
lo_t = movmedian(lo, w, 'omitnan');
el_t = movmedian(el, w, 'omitnan');
dev = step(la_t, lo_t, la, lo);
dev_el = abs(el - el_t);
bad = judgeable & (dev > dev_max);
bad_el = judgeable & (dev_el > elev_dev);

% group into runs and drop the ones too long to be glitches
[runs, kept, skipped] = H_runs(bad, max_run);
[runs_el, kept_el, skipped_el] = H_runs(bad_el, max_run);

fprintf('%s / gps_%s.mat: %d samples at %.3f s, source %s\n', season, day, n, dt_typ, g.gps_source);
fprintf('  trend window %d samples (%.1f s)\n', w, w * dt_typ);
fprintf('  horizontal: %d off-trend by > %.1f m, in %d run(s); repairing %d, leaving %d run(s) longer than %d samples\n', ...
  nnz(bad), dev_max, numel(runs), kept, skipped, max_run);
fprintf('  vertical:   %d off-trend by > %.1f m, in %d run(s); repairing %d, leaving %d\n', ...
  nnz(bad_el), elev_dev, numel(runs_el), kept_el, skipped_el);
fprintf('  departure before: median %.2f, p99.9 %.2f, max %.1f m\n', ...
  median(dev(judgeable)), prctile(dev(judgeable), 99.9), max(dev(judgeable)));
fprintf('  speed before:     median %.2f, p99.9 %.2f, max %.1f m/s\n', ...
  median(v(judgeable)), prctile(v(judgeable), 99.9), max(v(judgeable)));

fix_h = H_mask(runs, max_run, n);
fix_v = H_mask(runs_el, max_run, n);
if any(fix_h)
  la = H_interp(t, la, fix_h);
  lo = H_interp(t, lo, fix_h);
end
if any(fix_v)
  el = H_interp(t, el, fix_v);
end
d2 = [0; step(la(1:end-1), lo(1:end-1), la(2:end), lo(2:end))];
v2 = d2 ./ max(dt, eps);
dev2 = step(movmedian(la, w, 'omitnan'), movmedian(lo, w, 'omitnan'), la, lo);
fprintf('  departure after:  median %.2f, p99.9 %.2f, max %.1f m\n', ...
  median(dev2(judgeable)), prctile(dev2(judgeable), 99.9), max(dev2(judgeable)));
fprintf('  speed after:      median %.2f, p99.9 %.2f, max %.1f m/s\n', ...
  median(v2(judgeable)), prctile(v2(judgeable), 99.9), max(v2(judgeable)));

if dry
  fprintf('  dry run: nothing written\n');
  return;
end
g.lat = reshape(la, size(g.lat)); g.lon = reshape(lo, size(g.lon));
g.elev = reshape(el, size(g.elev));
if ~contains(g.gps_source, '+despike')
  g.gps_source = [g.gps_source '+despike'];
end
odir = fullfile(priv, 'gps', season);
if exist(odir, 'dir') ~= 7, mkdir(odir); end
out = fullfile(odir, ['gps_' day '.mat']);
opr_save(out, '-struct', 'g');
fid = fopen(fullfile(odir, 'despike_summary.txt'), 'a');
fprintf(fid, '%s gps_%s n=%d off_trend=%d(rep %d,skip %d) vert=%d dev %.1f->%.1f m src=%s\n', ...
  datestr(now, 'yyyy-mm-dd HH:MM'), day, n, nnz(bad), kept, skipped, nnz(bad_el), ...
  max(dev(judgeable)), max(dev2(judgeable)), g.gps_source);
fclose(fid);
fprintf('  wrote %s\n', out);

% -------------------------------------------------------------------------
function [runs, kept, skipped] = H_runs(bad, max_run)
runs = {};
i = 1; n = numel(bad); kept = 0; skipped = 0;
while i <= n
  if bad(i)
    j = i;
    while j < n && bad(j+1), j = j + 1; end
    runs{end+1} = [i j]; %#ok<AGROW>
    if (j - i + 1) <= max_run, kept = kept + (j - i + 1);
    else, skipped = skipped + 1; end
    i = j + 1;
  else
    i = i + 1;
  end
end
end

function m = H_mask(runs, max_run, n)
m = false(n, 1);
for k = 1:numel(runs)
  r = runs{k};
  if (r(2) - r(1) + 1) <= max_run, m(r(1):r(2)) = true; end
end
end

function y = H_interp(t, y, m)
good = ~m & isfinite(y);
if nnz(good) < 2, return; end
y(m) = interp1(t(good), y(good), t(m), 'linear', 'extrap');
end

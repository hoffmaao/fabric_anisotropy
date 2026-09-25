%TEST_BED_FROM_LAYER The per-trace bed reader binds layers by name and matches by position.
%
% bed_from_layer.m (opr_fabric/server) gives the pipeline the bed it masks
% the record below. What it must get right, each checked on a layer file
% written here to reproduce the traps the bed producer met on real files:
%
%   1. NAME, NOT ROW. The file carries `surface_dem` in row 1 and
%      `bottom_mc` in row 2, the Eastwind/McMurdo ordering that once turned
%      a 224 m shelf into a 1 m bed under a positional bind; `bottom` sits
%      in row 3 and must be the layer chosen (bottom_mc ranks below it, a
%      median 44.6 m shallow on real files).
%   2. POSITION MATCHING. Traces on the pick line get the pick's bed to
%      within a metre; traces farther than max_match_m (62.5 m, the bed
%      producer's radius) from every pick get NaN, never a far pick.
%   3. THE PIPELINE'S SURFACE. Depth is (twtt - surf_t) at solid-ice
%      velocity on eps 3.171, the pipeline's own axis, so the bed lands
%      where the fabric is reported.
%   4. MALFORMED PICKS DO NOT VOTE: a non-positive pick is dropped and
%      counted, and its traces take the nearest surviving pick.
%   5. The polarimetric tier: a file with bottom_HH/VV picks and no plain
%      bottom returns their per-trace median and names both.
%   6. REFUSAL: a frame whose twtt row count disagrees with the catalogue
%      (here the [Np x 3] transpose of a good file) is refused with the
%      reason in info, never transposed until the counts agree - that is a
%      positional bind wearing a name; no layer file, or no name catalogue,
%      likewise gives all-NaN with the reason - never a bed bound by
%      position. Traces that all lie beyond max_match_m of every pick get
%      all-NaN with the reason AND an empty source, so the pipeline's log
%      and its saved bed_source agree that no bed was applied. A file whose
%      only bottom-named layer is a DEM (surface_dem / bottom_dem) is
%      refused too: a DEM is a model, not a radar pick, as the bed
%      producer's NOT_A_PICK rules.
%   7. ANOTHER LAYOUT: a file whose twtt is a cell (one entry per layer,
%      as some v7 layer files store it) returns all-NaN with the reason,
%      never an error that would stop the frame.
%   8. GAPS. With the picks of a 400 m stretch removed, the traces more
%      than 62.5 m from any pick (a 300 m run) take the bed on the straight
%      line between the matched traces either side, and are counted in
%      n_filled; with 1200 m removed the 1100 m run stays NaN - a straight
%      line is not trusted that far across unseen bed.
%
% Run: matlab -batch "run('opr_fabric/test/test_bed_from_layer.m')"
clear;
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'server'));

C_ICE = 299792458 / sqrt(3.171);
root = fullfile(tempname, 'season');
seg = '20991231_01';
segdir = fullfile(root, 'CSARP_layer', seg);
mkdir(segdir);
cleanup = onCleanup(@() rmdir(fileparts(root), 's'));

% a 2 km E-W line of picks every 20 m at 77 S, bed 500 -> 700 m along it
Np = 101;
plat = -77.0 * ones(1, Np);
plon = 160.0 + (0:Np-1) * 20 / (111320 * cosd(-77.0));
bed_true = linspace(500, 700, Np);
surf_t = 3e-9;
tw_bot = surf_t + 2 * bed_true / C_ICE;
tw_mc = surf_t + 2 * (bed_true - 45) / C_ICE;    % the shallow multi-channel pick
tw_dem = (surf_t + 2 * 0.6 / C_ICE) * ones(1, Np);   % a DEM, not a radar pick
tw_bot(50) = surf_t - 1e-9;                      % one malformed (negative-depth) pick

fails = 0;

% --- the frame file and its segment catalogue, bottom in the LAST row
twtt = [tw_dem; tw_mc; tw_bot]; lat = plat; lon = plon; gps_time = 1:Np;
save(fullfile(segdir, sprintf('Data_%s_001.mat', seg)), 'twtt', 'lat', 'lon', 'gps_time');
lyr_name = {'surface_dem', 'bottom_mc', 'bottom'}; lyr_id = 1:3;
save(fullfile(segdir, sprintf('layer_%s.mat', seg)), 'lyr_name', 'lyr_id');

% traces: 4 per pick spacing on the line, offset 1 m so no trace sits
% exactly between two picks (a tie has two correct answers), plus a
% parallel line 200 m south
Nx = 4 * Np;
tx = 1 + (0:Nx-1) * 5;                           % along-line metres
tlon = 160.0 + tx / (111320 * cosd(-77.0));
tlat = -77.0 * ones(1, Nx);
tlat_off = tlat - 200 / 111320;
[zb, info] = bed_from_layer(root, seg, 1, [tlat tlat_off], [tlon tlon], surf_t);
on = 1:Nx; off = Nx + (1:Nx);

ok1 = strcmp(info.source, 'bottom');
fprintf('1. bound by name: source "%s" (rows were surface_dem, bottom_mc, bottom): %s\n', ...
  info.source, H_tick(ok1));
fails = fails + ~ok1;

% expected bed at each on-line trace: the nearest pick's, with pick 50 gone
want = interp1(1:Np, bed_true, 1 + tx / 20, 'nearest');
near = round(1 + tx / 20);
lost = near == 50;                      % nearest pick was the malformed one
err = abs(zb(on) - want);
ok2 = all(isfinite(zb(on))) && max(err(~lost)) < 1.0 && all(isnan(zb(off)));
fprintf('2. on-line traces within %.2f m of their pick, off-line (200 m) all NaN: %s\n', ...
  max(err(~lost)), H_tick(ok2));
fails = fails + ~ok2;

ok3 = abs(zb(1) - 500) < 1.0 && abs(zb(on(end)) - 700) < 1.0;
fprintf('3. depth on the pipeline axis: first %.1f m (500), last %.1f m (700): %s\n', ...
  zb(1), zb(on(end)), H_tick(ok3));
fails = fails + ~ok3;

ok4 = info.n_bad == 1 && info.n_picks == Np - 1 && all(isfinite(zb(on(lost)))) ...
  && all(abs(zb(on(lost)) - bed_true(49)) < 3 | abs(zb(on(lost)) - bed_true(51)) < 3);
fprintf('4. the malformed pick is dropped (n_bad %d) and its traces take a neighbour: %s\n', ...
  info.n_bad, H_tick(ok4));
fails = fails + ~ok4;

% --- the polarimetric tier: bottom_HH and bottom_VV, no plain bottom
twtt = [tw_dem; tw_bot + 2 * 2 / C_ICE; tw_bot - 2 * 2 / C_ICE];
save(fullfile(segdir, sprintf('Data_%s_002.mat', seg)), 'twtt', 'lat', 'lon', 'gps_time');
lyr_name = {'surface_dem', 'bottom_HH', 'bottom_VV'};
save(fullfile(segdir, sprintf('layer_%s.mat', seg)), 'lyr_name', 'lyr_id');
[zb5, info5] = bed_from_layer(root, seg, 2, tlat, tlon, surf_t);
ok5 = strcmp(info5.source, 'median(bottom_hh,bottom_vv)') && abs(zb5(1) - 500) < 1.0;
fprintf('5. polarimetric picks: source "%s", first %.1f m: %s\n', info5.source, zb5(1), H_tick(ok5));
fails = fails + ~ok5;

% --- refusals
twtt = twtt.';                       % [Np x 3]: rows are no longer layers
save(fullfile(segdir, sprintf('Data_%s_003.mat', seg)), 'twtt', 'lat', 'lon', 'gps_time');
[zb6, info6] = bed_from_layer(root, seg, 3, tlat, tlon, surf_t);
[zb7, info7] = bed_from_layer(root, seg, 4, tlat, tlon, surf_t);
[zb9, info9] = bed_from_layer(root, seg, 2, tlat_off, tlon, surf_t);
twtt = [tw_dem; tw_bot];             % a DEM under a DEM name, no pick
save(fullfile(segdir, sprintf('Data_%s_005.mat', seg)), 'twtt', 'lat', 'lon', 'gps_time');
lyr_name = {'surface_dem', 'bottom_dem'}; lyr_id = 1:2;
save(fullfile(segdir, sprintf('layer_%s.mat', seg)), 'lyr_name', 'lyr_id');
[zb10, info10] = bed_from_layer(root, seg, 5, tlat, tlon, surf_t);
delete(fullfile(segdir, sprintf('layer_%s.mat', seg)));
[zb8, info8] = bed_from_layer(root, seg, 2, tlat, tlon, surf_t);
ok6 = all(isnan(zb6)) && isempty(info6.source) && ~isempty(info6.reason) ...
  && all(isnan(zb7)) && ~isempty(info7.reason) && all(isnan(zb8)) && ~isempty(info8.reason) ...
  && all(isnan(zb9)) && isempty(info9.source) && ~isempty(info9.reason) ...
  && all(isnan(zb10)) && isempty(info10.source) && ~isempty(info10.reason);
fprintf(['6. transposed twtt -> "%s";\n   no layer file -> "%s";\n' ...
  '   no catalogue -> "%s";\n   no trace near a pick -> "%s" (source "%s");\n' ...
  '   only a DEM bottom -> "%s" (source "%s"): %s\n'], ...
  info6.reason, info7.reason, info8.reason, info9.reason, info9.source, ...
  info10.reason, info10.source, H_tick(ok6));
fails = fails + ~ok6;

% --- another layout: twtt as a cell
lyr_name = {'bottom'}; lyr_id = 1;
save(fullfile(segdir, sprintf('layer_%s.mat', seg)), 'lyr_name', 'lyr_id');
twtt = {tw_bot};
save(fullfile(segdir, sprintf('Data_%s_008.mat', seg)), 'twtt', 'lat', 'lon', 'gps_time');
[zb11, info11] = bed_from_layer(root, seg, 8, tlat, tlon, surf_t);
ok7 = all(isnan(zb11)) && isempty(info11.source) && ~isempty(info11.reason);
fprintf('7. twtt as a cell -> "%s": %s\n', info11.reason, H_tick(ok7));
fails = fails + ~ok7;

% --- gaps: picks 41-60 removed (400 m), then picks 21-80 (1200 m)
tw_good = surf_t + 2 * bed_true / C_ICE;
gaps = {41:60, 21:80}; zg = cell(1, 2); ig = cell(1, 2);
for gi = 1:2
  twtt = tw_good; twtt(gaps{gi}) = NaN;
  save(fullfile(segdir, sprintf('Data_%s_%03d.mat', seg, 8 + gi)), 'twtt', 'lat', 'lon', 'gps_time');
  [zg{gi}, ig{gi}] = bed_from_layer(root, seg, 8 + gi, tlat, tlon, surf_t);
end
far = @(g) tx > (g(1) - 2) * 20 + 62.5 & tx < g(end) * 20 - 62.5;   % beyond 62.5 m of every pick
f1 = far(gaps{1}); i0 = find(f1, 1) - 1; i1 = find(f1, 1, 'last') + 1;
line1 = zg{1}(i0) + (zg{1}(i1) - zg{1}(i0)) * (tx(f1) - tx(i0)) / (tx(i1) - tx(i0));
e_line = max(abs(zg{1}(f1) - line1));
e_true = max(abs(zg{1}(f1) - (500 + tx(f1) / 10)));
f2 = far(gaps{2});
ok8 = ig{1}.n_filled == nnz(f1) && e_line < 0.05 && e_true < 7 ...
  && all(isnan(zg{2}(f2))) && ig{2}.n_filled == 0 && all(isfinite(zg{2}(~f2)));
fprintf(['8. a %.0f m run beyond every pick filled (%d traces, n_filled %d) on the line to %.3f m, ' ...
  'within %.1f m of the bed; a %.0f m run left NaN (n_filled %d): %s\n'], ...
  tx(i1) - tx(i0), nnz(f1), ig{1}.n_filled, e_line, e_true, ...
  tx(find(f2, 1, 'last') + 1) - tx(find(f2, 1) - 1), ig{2}.n_filled, H_tick(ok8));
fails = fails + ~ok8;

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_bed_from_layer:failed', '%d check(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

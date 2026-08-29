%RUN_SEASON_TAYLOR_DOME 2025_Antarctica_Ground2: the Taylor Dome strip.
%
% Tags 2026xxxx (January 2026). ~40 frames along a narrow ~5 x 15 km
% strip. Fabric is weak-to-moderate (deep dlam ~0.015, non-monotonic with
% a possible recovery below 1000 m) and the axis scatter is large, so
% treat single-frame orientations with caution. Standard products.
%
%   matlab -batch "addpath('<code>/opr_fabric/server'); \
%                  day_seg='20260109_03'; frm=17; run_season_taylor_dome"
season = struct( ...
  'name', 'taylor_dome', ...
  'site_root', '/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2');
run_season_frame;

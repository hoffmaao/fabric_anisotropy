%RUN_SEASON_RIDGE_A 2024_Antarctica_Ground2: the Ridge A raster.
%
% Tags 202501xx (January 2025). 36 frames over a ~27 x 23 km grid at the
% Dome A flank divide, ~1 m trace spacing, ice ~2912 m thick (rds bed
% pick) with the accum record reaching ~1879 m. The reference season:
% every pipeline default was established here, so the declaration is the
% site root and nothing else. Full CSARP_polarimetric products (window,
% coregistration settings, Surface, ref). Pedestal family a3 = +0.27,
% frame-stable. Fabric: theta0_geo 96.5 +- 7.4, deep dlam ~0.06.
%
%   matlab -batch "day_seg='20250108_02'; frm=9; run_season_ridge_a"
season = struct( ...
  'name', 'ridge_a', ...
  'site_root', '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2');
run_season_frame;

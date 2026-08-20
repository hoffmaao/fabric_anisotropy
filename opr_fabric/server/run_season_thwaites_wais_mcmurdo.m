%RUN_SEASON_THWAITES_WAIS_MCMURDO 2023_Antarctica_Ground: three sites.
%
% Tags 2024xxxx (January-February 2024): the Thwaites eastern shear margin
% transect (2024010x), WAIS Divide (20240120) and the McMurdo Ice Shelf
% transect (20240203). One season, one instrument install - the pedestal
% family is a3 = -0.195, SIGN-FLIPPED from Ridge A, which is why chan_equal
% calibration must be fitted per season and never assumed global.
% Polarimetric products ship only as CSARP_polarimetric_unwrap; the
% pipeline falls back to it automatically. Thwaites: depth-rotating axes,
% honest abstention in the margin core - never depth-average theta0 there.
%
%   matlab -batch "day_seg='20240108_01'; frm=1; run_season_thwaites_wais_mcmurdo"
season = struct( ...
  'name', 'thwaites_wais_mcmurdo', ...
  'site_root', '/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground');
run_season_frame;

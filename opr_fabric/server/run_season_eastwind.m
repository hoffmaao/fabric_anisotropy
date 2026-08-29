%RUN_SEASON_EASTWIND 2022_Antarctica_Ground: the Eastwind rotation site.
%
% Tags 202212xx (December 2022). 13 soundings within ~200 m of one spot on
% the McMurdo Ice Shelf - a rotation experiment, not a survey - 15 km from
% the 2024 transect and separated from it BY SEASON in quadpol_sites.py.
% Thin ice: the record reaches ~561 m and the HH-VV coherence collapses to
% ~0.02 below ~340 m (apparent deep dlam there is noise; never animate
% it). Polarimetric products in _unwrap form; automatic fallback.
%
%   matlab -batch "addpath('<code>/opr_fabric/server'); \
%                  day_seg='20221210_02'; frm=1; run_season_eastwind"
season = struct( ...
  'name', 'eastwind', ...
  'site_root', '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground');
run_season_frame;

%RUN_SEASON_EASTGRIP 2024_Greenland_Ground2: the EastGRIP/NEGIS traverse.
%
% Tags 202406xx (June 2024). The season that motivated the season-driver
% split: it ships NONE of what the Antarctic seasons carry, and the
% pipeline's qlook mode supplies each missing piece (leading-edge surface
% pick - Surface is all NaN; stationary cull - a third of traces are
% dwells at ~9 m spacing; coregistration defaults; range window from the
% picked surface). Those engage automatically when CSARP_polarimetric is
% absent, so this declaration only carries what auto-detection cannot
% know:
%   - dlam_max 0.45: the EGRIP core contrast (~0.3) sits ABOVE the
%     estimator's 0.25 default; test_egrip_cap.m pins the failure mode.
%     (qlook mode also defaults to 0.45; declared here so the season file
%     is the readable record, not the mode's side effect.)
%   - SKIP segment 20240618_01: a real-only setup day; the pipeline
%     errors on it deliberately. ~9 frames also lack a VH channel and
%     fail loudly into capped markers.
% Pedestal family a3 = -0.17 (sign-flipped like Thwaites). DO NOT batch
% ahead of the egrip_chain.sh acceptance gate: the first validation frame
% still fails it at resid 0.476 (0.485 before the segmented frame pass
% recovered the lateral domain structure; the remaining floor is
% un-modeled site physics/calibration - the chan_equal family - a
% separate decision). Batch via egrip_chain.sh, never by hand.
%
%   matlab -batch "day_seg='20240619_01'; frm=1; run_season_eastgrip"
season = struct( ...
  'name', 'eastgrip', ...
  'site_root', '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2', ...
  'dlam_max', 0.45);
run_season_frame;

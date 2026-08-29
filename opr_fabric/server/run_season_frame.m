%RUN_SEASON_FRAME Apply a season declaration and invert one frame.
%
% The shared tail of every run_season_* driver. A season script declares
% WHAT is different about its acquisition in a `season` struct - paths,
% search ceilings, block policy - and this applies that declaration to
% run_quadpol_pipeline's variable interface and runs it. The pipeline
% itself stays season-agnostic: everything it needs to know arrives either
% from the products (preferred) or from this declaration (when the season
% shipped without the product that would carry it).
%
% Variables already set at the call site WIN over the season defaults, so
% special runs stay one-liners on top of the season driver:
%
%   matlab -batch "addpath('<code>/opr_fabric/server'); \
%                  day_seg='20250108_02'; frm=9; run_season_ridge_a"
%   matlab -batch "addpath('<code>/opr_fabric/server'); \
%                  day_seg='20250108_02'; frm=9; z_max=4000; run_season_ridge_a"
%
% This file must stay BEHAVIOR-FREE: no physics, no thresholds, no season
% facts. Those belong in the season scripts (facts) or the pipeline
% (physics), or the split stops meaning anything.
if ~exist('season', 'var')
  error('run_season_frame:noSeason', ...
    'run via a run_season_* driver, which defines the season struct');
end
if ~exist('day_seg', 'var') || ~exist('frm', 'var')
  error('run_season_frame:noFrame', 'set day_seg and frm at the call site');
end
if ~exist('site_root', 'var')
  site_root = season.site_root;
end
if isfield(season, 'z_max') && ~exist('z_max', 'var')
  z_max = season.z_max;
end
if isfield(season, 'nblk_tr') && ~exist('nblk_tr', 'var')
  nblk_tr = season.nblk_tr;
end
if isfield(season, 'dlam_max') && ~exist('dlam_max', 'var')
  dlam_max = season.dlam_max;
end
run_quadpol_pipeline;
